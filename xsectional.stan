// See the associated R file for explanations and definitions of abbreviations.

functions {
  real PL4(real x_log, real f_1, real f_2, real f_3, real f_4) {
    return f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (x_log - f_4)));
  }
}


data {
  // Actual data
  int<lower = 1> num_bat;
  int<lower = num_bat> num_con_tot;
  int<lower = 0> num_sam_id;
  int<lower = num_sam_id> num_sam_tot;
  array[num_con_tot] int<lower = 1, upper = num_bat> which_bat_con;
  array[num_sam_tot] int<lower = 1, upper = num_bat> which_bat_sam;
  array[num_sam_tot] int<lower = 1, upper = num_sam_id> which_id_sam;
  vector[num_con_tot] y_con;
  vector[num_con_tot] x_con;
  vector[num_sam_tot] y_sam;

  // Other things to keep fixed over a complete round of sampling: a binary
  // switch to control whether we sample from the prior or the posterior
  // (important to compare the difference), switches for using normal or student 
  // t distributions, and upper and lower bounds for the priors.
  int<lower = 0, upper = 1> sample_posterior_not_prior;
  int<lower = 0, upper = 1> use_student_for_x;
  real<lower = 0> student_df_x;
  int<lower = 0, upper = 1> use_student_for_obs;
  real<lower = 0> student_df_obs;
  vector[4] f_lower;
  vector[4] f_upper;
  vector[4] sigma_bat_f_lower;
  vector[4] sigma_bat_f_upper;
  real y_obs_sd_lower;
  real y_obs_sd_upper;
  real x_sam_neg_mu_lower;
  real x_sam_neg_mu_upper;
  real x_sam_pos_mu_lower;
  real x_sam_pos_mu_upper;
  real x_sam_neg_sd_lower;
  real x_sam_neg_sd_upper;
  real x_sam_pos_sd_lower;
  real x_sam_pos_sd_upper;
  real p_sam_pos_lower;
  real p_sam_pos_upper;
  real Rho_bat_prior_eta;
  
}

transformed data {
  array[num_bat] vector[4] zeros;
  for (i in 1:num_bat) zeros[i] = rep_vector(0, 4);
}

parameters {
  
  // Those constrained only by lower and upper
  vector<lower = f_lower, upper = f_upper>[4] f;
  vector<lower = sigma_bat_f_lower, upper = sigma_bat_f_upper>[4] sigma_bat_f;
  
  real<lower = y_obs_sd_lower, upper = y_obs_sd_upper> y_obs_sd;
  real<lower = x_sam_neg_mu_lower, upper = x_sam_neg_mu_upper> x_sam_neg_mu;
  real<lower = x_sam_neg_sd_lower, upper = x_sam_neg_sd_upper> x_sam_neg_sd;
  real<lower = x_sam_pos_sd_lower, upper = x_sam_pos_sd_upper> x_sam_pos_sd;
  real<lower = p_sam_pos_lower,    upper = p_sam_pos_upper>    p_sam_pos;

  // Enforce that x_sam_pos_mu > x_sam_neg_mu
  real<lower = max([x_sam_pos_mu_lower, x_sam_neg_mu]), upper = x_sam_pos_mu_upper> x_sam_pos_mu;
  
  // Those with explicit priors declared
  corr_matrix[4] Rho_bat;
  array[num_bat] vector[4] f_bat_effects_unscaled;
  vector[num_sam_id] x_sam;
}

transformed parameters{
  
  real p_sam_pos_log = log(  p_sam_pos);
  real p_sam_neg_log = log1m(p_sam_pos);
  
  array[num_bat] vector[4] f_per_bat;
  for (bat in 1:num_bat) {
    f_per_bat[bat] = f +
    [
      f_bat_effects_unscaled[bat][1] * sigma_bat_f[1],
      f_bat_effects_unscaled[bat][2] * sigma_bat_f[2],
      f_bat_effects_unscaled[bat][3] * sigma_bat_f[3],
      f_bat_effects_unscaled[bat][4] * sigma_bat_f[4]
    ]';
  }

  vector[num_con_tot] y_con_mean_per_obs;
  for (con_rep in 1:num_con_tot) {
    int bat = which_bat_con[con_rep];
    y_con_mean_per_obs[con_rep] =  PL4(x_con[con_rep],
    f_per_bat[bat][1], f_per_bat[bat][2], f_per_bat[bat][3], f_per_bat[bat][4]);
  }
  
  vector[num_sam_tot] y_sam_mean_per_obs;
  for (sam_rep in 1:num_sam_tot) {
    int bat = which_bat_sam[sam_rep];
    y_sam_mean_per_obs[sam_rep] =  PL4(x_sam[which_id_sam[sam_rep]],
    f_per_bat[bat][1], f_per_bat[bat][2], f_per_bat[bat][3], f_per_bat[bat][4]);
  }

}


model {
  
  // Priors
  Rho_bat ~ lkj_corr(Rho_bat_prior_eta);
  f_bat_effects_unscaled ~ multi_normal(zeros, Rho_bat);
  x_sam_pos_mu ~ uniform(max([x_sam_pos_mu_lower, x_sam_neg_mu]), x_sam_pos_mu_upper);
  if (use_student_for_x) {
    for (sam_id in 1:num_sam_id) {
      target += log_sum_exp(
        p_sam_pos_log + student_t_lpdf(x_sam[sam_id] | student_df_x, x_sam_pos_mu, x_sam_pos_sd),
        p_sam_neg_log + student_t_lpdf(x_sam[sam_id] | student_df_x, x_sam_neg_mu, x_sam_neg_sd));
    }
  } else {
    for (sam_id in 1:num_sam_id) {
      target += log_sum_exp(
        p_sam_pos_log + normal_lpdf(x_sam[sam_id] | x_sam_pos_mu, x_sam_pos_sd),
        p_sam_neg_log + normal_lpdf(x_sam[sam_id] | x_sam_neg_mu, x_sam_neg_sd));
    }
  }
  
  // Likelihood
  if (sample_posterior_not_prior) {
    if (use_student_for_obs) {
      y_con ~ student_t(student_df_obs, y_con_mean_per_obs, y_obs_sd);
      y_sam ~ student_t(student_df_obs, y_sam_mean_per_obs, y_obs_sd);
    } else {
      y_con ~ normal(y_con_mean_per_obs, y_obs_sd);
      y_sam ~ normal(y_sam_mean_per_obs, y_obs_sd);
    }
  }
} 

generated quantities {
  
  array[num_con_tot] real y_con_sim;
  if (use_student_for_obs) {
    y_con_sim = student_t_rng(student_df_obs, y_con_mean_per_obs, y_obs_sd);
  } else {
    y_con_sim = normal_rng(y_con_mean_per_obs, y_obs_sd);
  }
  
  vector[num_sam_id] p_sam_is_pos;
  if (use_student_for_x) {
    for (sam_id in 1:num_sam_id) {
      real p_log = p_sam_pos_log +
        student_t_lpdf(x_sam[sam_id] | student_df_x, x_sam_pos_mu, x_sam_pos_sd);
      p_sam_is_pos[sam_id] = exp(p_log - log_sum_exp(p_log, p_sam_neg_log +
        student_t_lpdf(x_sam[sam_id] | student_df_x, x_sam_neg_mu, x_sam_neg_sd)));
    }
  } else {
    for (sam_id in 1:num_sam_id) {
      real p_log = p_sam_pos_log +
        normal_lpdf(x_sam[sam_id] | x_sam_pos_mu, x_sam_pos_sd);
      p_sam_is_pos[sam_id] = exp(p_log - log_sum_exp(p_log, p_sam_neg_log +
        normal_lpdf(x_sam[sam_id] | x_sam_neg_mu, x_sam_neg_sd)));
    }
  }

}
