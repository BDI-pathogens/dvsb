// See the associated R file for explanations and definitions of abbreviations.

functions {
  real PL4(real x, real f_1, real f_2, real f_3, real exp_f_1_mult_f_4) {
    //return f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (xlog - f_4)));
    return f_2 + (f_3 - f_2) / (1 + x^(-f_1) * exp_f_1_mult_f_4);
  }
}


data {
  // Actual data
  int<lower = 1> num_plate;
  int<lower = num_plate> num_cal_tot;
  int<lower = 0> num_sam_id;
  int<lower = num_sam_id> num_sam_tot;
  array[num_cal_tot] int<lower = 1, upper = num_plate> which_plate_cal;
  array[num_sam_tot] int<lower = 1, upper = num_plate> which_plate_sam;
  array[num_sam_tot] int<lower = 1, upper = num_sam_id> which_id_sam;
  vector[num_cal_tot] y_cal;
  vector<lower = 0>[num_cal_tot] x_cal;
  vector[num_sam_tot] y_sam;
  
  // Other things to keep fixed over a complete round of sampling: a binary
  // switch to control whether we sample from the prior or the posterior
  // (important to compare the difference), switches for using normal or student 
  // t distributions, and upper and lower bounds for the priors.
  int<lower = 0, upper = 1> sample_posterior_not_prior;
  int<lower = 0, upper = 1> use_student_for_obs;
  real<lower = 0> student_df_obs;
  vector[4] f_lower;
  vector[4] f_upper;
  vector[4] sigma_f_plate_lower;
  vector[4] sigma_f_plate_upper;
  real y_obs_sd_min_lower;
  real y_obs_sd_min_upper;
  real y_obs_sd_jump_lower;
  real y_obs_sd_jump_upper;
  
  real<lower = 0> x_sam_neg_alpha_lower;
  real<lower = x_sam_neg_alpha_lower> x_sam_neg_alpha_upper;
  real<lower = 0> x_sam_pos_alpha_jump_lower;
  real<lower = x_sam_pos_alpha_jump_lower> x_sam_pos_alpha_jump_upper;
  real<lower = 0> x_sam_pos_beta_lower;
  real<lower = x_sam_pos_beta_lower> x_sam_pos_beta_upper;
  real<lower = 0> x_sam_neg_beta_jump_lower;
  real<lower = x_sam_neg_beta_jump_lower> x_sam_neg_beta_jump_upper;
  
  real p_sam_pos_lower;
  real p_sam_pos_upper;
  real Rho_plate_prior_eta;
  
}

transformed data {
  array[num_plate] vector[4] zeros;
  for (i in 1:num_plate) zeros[i] = rep_vector(0, 4);
}

parameters {
  
  // Those constrained only by lower and upper
  vector<lower = f_lower, upper = f_upper>[4] f;
  vector<lower = sigma_f_plate_lower, upper = sigma_f_plate_upper>[4] sigma_f_plate;
  
  real<lower = x_sam_neg_alpha_lower, upper = x_sam_neg_alpha_upper> x_sam_neg_alpha;
  real<lower = x_sam_pos_beta_lower, upper = x_sam_pos_beta_upper> x_sam_pos_beta;
  real<lower = x_sam_neg_beta_jump_lower, upper = x_sam_neg_beta_jump_upper> x_sam_neg_beta_jump;
  real<lower = p_sam_pos_lower,    upper = p_sam_pos_upper>    p_sam_pos;
  real<lower = y_obs_sd_min_lower, upper = y_obs_sd_min_upper> y_obs_sd_min;
  real<lower = y_obs_sd_jump_lower, upper = y_obs_sd_jump_upper> y_obs_sd_jump;
  real<lower = x_sam_pos_alpha_jump_lower, upper = x_sam_pos_alpha_jump_upper> x_sam_pos_alpha_jump;
  
  // Those with explicit priors declared
  corr_matrix[4] Rho_plate;
  array[num_plate] vector[4] f_plate_effects_unscaled;
  vector<lower = 0>[num_sam_id] x_sam;
}

transformed parameters{
  
  real p_sam_pos_log = log(  p_sam_pos);
  real p_sam_neg_log = log1m(p_sam_pos);
  real y_obs_sd_max = y_obs_sd_min + y_obs_sd_jump;
  
  real x_sam_pos_alpha = x_sam_neg_alpha + x_sam_pos_alpha_jump;
  real x_sam_neg_beta  = x_sam_pos_beta  + x_sam_neg_beta_jump;

  array[num_plate] vector[4] f_per_plate;
  for (plate in 1:num_plate) {
    f_per_plate[plate] = f +
    [
      f_plate_effects_unscaled[plate][1] * sigma_f_plate[1],
      f_plate_effects_unscaled[plate][2] * sigma_f_plate[2],
      f_plate_effects_unscaled[plate][3] * sigma_f_plate[3],
      f_plate_effects_unscaled[plate][4] * sigma_f_plate[4]
      ]';
  }
  array[num_plate] real exp_f_1_mult_f_4_per_plate;
  for (plate in 1:num_plate) {
    exp_f_1_mult_f_4_per_plate[plate] =
      exp(f_per_plate[plate][1] * f_per_plate[plate][4]);
  }

  
  vector[num_cal_tot] y_cal_mean_per_obs;
  vector[num_sam_tot] y_sam_mean_per_obs;
  profile("y_means") {
  for (cal_rep in 1:num_cal_tot) {
    int plate = which_plate_cal[cal_rep];
    y_cal_mean_per_obs[cal_rep] =  PL4(x_cal[cal_rep],
    f_per_plate[plate][1], f_per_plate[plate][2], f_per_plate[plate][3], 
    exp_f_1_mult_f_4_per_plate[plate]);
  }
  for (sam_rep in 1:num_sam_tot) {
    int plate = which_plate_sam[sam_rep];
    y_sam_mean_per_obs[sam_rep] =  PL4(x_sam[which_id_sam[sam_rep]],
    f_per_plate[plate][1], f_per_plate[plate][2], f_per_plate[plate][3], 
    exp_f_1_mult_f_4_per_plate[plate]);
  }
  }
  
  vector[num_cal_tot] y_obs_sd_cal;
  vector[num_sam_tot] y_obs_sd_sam;
  profile("y_sds") {
  for (cal_rep in 1:num_cal_tot) {
    int plate = which_plate_cal[cal_rep];
    y_obs_sd_cal[cal_rep] = PL4(x_cal[cal_rep], f_per_plate[plate][1],
    y_obs_sd_min, y_obs_sd_max, exp_f_1_mult_f_4_per_plate[plate]);
  }
  for (sam_rep in 1:num_sam_tot) {
    int plate = which_plate_sam[sam_rep];
    y_obs_sd_sam[sam_rep] = PL4(x_sam[which_id_sam[sam_rep]], f_per_plate[plate][1],
    y_obs_sd_min, y_obs_sd_max, exp_f_1_mult_f_4_per_plate[plate]);
  }
  }
  
}


model {
  
  // Priors
  Rho_plate ~ lkj_corr(Rho_plate_prior_eta);
  f_plate_effects_unscaled ~ multi_normal(zeros, Rho_plate);
  profile("mixture_model") {
  for (sam_id in 1:num_sam_id) {
    target += log_sum_exp(
      p_sam_pos_log + gamma_lpdf(x_sam[sam_id] | x_sam_pos_alpha, x_sam_pos_beta),
      p_sam_neg_log + gamma_lpdf(x_sam[sam_id] | x_sam_neg_alpha, x_sam_neg_beta));
  }
  }
  
  
  // Likelihood
  profile("likelihood") {
  if (sample_posterior_not_prior) {
    if (use_student_for_obs) {
      y_cal ~ student_t(rep_vector(student_df_obs, num_cal_tot), y_cal_mean_per_obs, y_obs_sd_cal);
      y_sam ~ student_t(rep_vector(student_df_obs, num_sam_tot), y_sam_mean_per_obs, y_obs_sd_sam);
    } else {
      y_cal ~ normal(y_cal_mean_per_obs, y_obs_sd_cal);
      y_sam ~ normal(y_sam_mean_per_obs, y_obs_sd_sam);
    }
  }
  }
} 

generated quantities {
  
  array[num_cal_tot] real y_cal_sim;
  if (use_student_for_obs) {
    y_cal_sim = student_t_rng(rep_vector(student_df_obs, num_cal_tot), y_cal_mean_per_obs, y_obs_sd_cal);
  } else {
    y_cal_sim = normal_rng(y_cal_mean_per_obs, y_obs_sd_cal);
  }
  
  array[num_sam_tot] real y_sam_sim;
  if (use_student_for_obs) {
    y_sam_sim = student_t_rng(rep_vector(student_df_obs, num_sam_tot), y_sam_mean_per_obs, y_obs_sd_sam);
  } else {
    y_sam_sim = normal_rng(y_sam_mean_per_obs, y_obs_sd_sam);
  }
  
  vector[num_sam_id] p_sam_is_pos;
  for (sam_id in 1:num_sam_id) {
    real p_log = p_sam_pos_log +
    gamma_lpdf(x_sam[sam_id] | x_sam_pos_alpha, x_sam_pos_beta);
    p_sam_is_pos[sam_id] = exp(p_log - log_sum_exp(p_log, p_sam_neg_log +
    gamma_lpdf(x_sam[sam_id] | x_sam_neg_alpha, x_sam_neg_beta)));
  }
  
  real x_sam_pos_mu = x_sam_pos_alpha / x_sam_pos_beta;
  real x_sam_neg_mu = x_sam_neg_alpha / x_sam_neg_beta;
  
  
}
