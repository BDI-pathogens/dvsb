// See the associated R file for explanations and definitions of abbreviations.

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
  int<lower = 0> num_f_pred_vars;
  array[num_f_pred_vars] int<lower = 2> num_cat_per_f_pred_var;
  matrix<lower = 0, upper = 1>[num_plate, sum(num_cat_per_f_pred_var)] design_matrix_f;
  int<lower = 0> num_p_pos_pred_vars;
  array[num_p_pos_pred_vars] int<lower = 2> num_cat_per_p_pos_pred_var;
  matrix<lower = 0, upper = 1>[num_sam_id, sum(num_cat_per_p_pos_pred_var)] design_matrix_p_pos;

  // Other things to keep fixed over a complete round of sampling: a binary
  // switch to control whether we sample from the prior or the posterior
  // (important to compare the difference), and upper and lower bounds for the priors.
  int<lower = 0, upper = 1> sample_posterior_not_prior;
  row_vector[4] f_lower;
  row_vector[4] f_upper;
  row_vector[4] sigma_f_plate_lower;
  row_vector[4] sigma_f_plate_upper;
  row_vector[4] sigma_f_pred_vars_lower;
  row_vector[4] sigma_f_pred_vars_upper;
  real sigma_p_pos_pred_vars_lower;
  real sigma_p_pos_pred_vars_upper;
  real y_obs_sd_cal_min_lower;
  real y_obs_sd_cal_min_upper;
  real y_obs_sd_cal_jump_lower;
  real y_obs_sd_cal_jump_upper;
  real y_obs_sd_sam_min_lower;
  real y_obs_sd_sam_min_upper;
  real y_obs_sd_sam_jump_lower;
  real y_obs_sd_sam_jump_upper;
  
  real x_sam_neg_mu_lower;
  real<lower = x_sam_neg_mu_lower> x_sam_neg_mu_upper;
  real<lower = x_sam_neg_mu_lower> x_sam_pos_mu_lower;
  real<lower = x_sam_pos_mu_lower> x_sam_pos_mu_upper;
  real<lower = 0> x_sam_neg_sd_lower;
  real<lower = x_sam_neg_sd_lower> x_sam_neg_sd_upper;
  real<lower = 0> x_sam_pos_sd_lower;
  real<lower = x_sam_pos_sd_lower> x_sam_pos_sd_upper;

  real p_pos_lower;
  real p_pos_upper;
  real rho_prior_eta;
  
}

transformed data {
  
  array[num_plate] vector[4] zeros;
  for (i in 1:num_plate) zeros[i] = rep_vector(0, 4);
  vector[4] zeros_4 = rep_vector(0, 4);
  
  int tot_cat_per_f_pred_var = sum(num_cat_per_f_pred_var);
  array[tot_cat_per_f_pred_var] vector[4] zeros_for_f_pred_vars;
  for (i in 1:tot_cat_per_f_pred_var) zeros_for_f_pred_vars[i] = rep_vector(0, 4);
  int predict_f = 1 ? num_f_pred_vars > 0 : 0;
  array[num_f_pred_vars] row_vector[4] sigma_f_pred_vars_lower_array;
  array[num_f_pred_vars] row_vector[4] sigma_f_pred_vars_upper_array;
  for (f_pred_var in 1:num_f_pred_vars) {
    sigma_f_pred_vars_lower_array[f_pred_var] = sigma_f_pred_vars_lower;
    sigma_f_pred_vars_upper_array[f_pred_var] = sigma_f_pred_vars_upper;
  }
  
  int tot_cat_per_p_pos_pred_var = sum(num_cat_per_p_pos_pred_var);
  int predict_p_pos = 1 ? num_p_pos_pred_vars > 0 : 0;
  
  vector[num_cal_tot] xlog_cal = log(x_cal);
  array[num_cal_tot] int x_cal_is_zero;
  for (cal_rep in 1:num_cal_tot) {
    x_cal_is_zero[cal_rep] = 1 ? x_cal[cal_rep] == 0 : 0;
  }

}

parameters {
  
  // Tensor params constrained only by lower and upper
  row_vector<lower = f_lower, upper = f_upper>[4] f;
  row_vector<lower = sigma_f_plate_lower, upper = sigma_f_plate_upper>[4] sigma_f_plate;
  array[num_f_pred_vars] row_vector<lower = sigma_f_pred_vars_lower_array,
  upper = sigma_f_pred_vars_upper_array>[4] sigma_f_pred_vars;
  array[num_p_pos_pred_vars] real<lower = sigma_p_pos_pred_vars_lower,
  upper = sigma_p_pos_pred_vars_upper> sigma_p_pos_pred_vars;
  
  // Scalar params constrained only by lower and upper
  real<lower = x_sam_neg_mu_lower, upper = x_sam_neg_mu_upper> x_sam_neg_mu;
  real<lower = x_sam_neg_sd_lower, upper = x_sam_neg_sd_upper> x_sam_neg_sd;
  real<lower = x_sam_pos_sd_lower, upper = x_sam_pos_sd_upper> x_sam_pos_sd;
  real<lower = p_pos_lower,    upper = p_pos_upper>    p_pos;
  real<lower = y_obs_sd_cal_min_lower,  upper = y_obs_sd_cal_min_upper>  y_obs_sd_cal_min;
  real<lower = y_obs_sd_cal_jump_lower, upper = y_obs_sd_cal_jump_upper> y_obs_sd_cal_jump;
  real<lower = y_obs_sd_sam_min_lower,  upper = y_obs_sd_sam_min_upper>  y_obs_sd_sam_min;
  real<lower = y_obs_sd_sam_jump_lower, upper = y_obs_sd_sam_jump_upper> y_obs_sd_sam_jump;
  
  // Enforce that x_sam_pos_mu > x_sam_neg_mu
  real<lower = max([x_sam_pos_mu_lower, x_sam_neg_mu]), upper = x_sam_pos_mu_upper> x_sam_pos_mu;

  
  // Those with explicit priors declared
  corr_matrix[4] rho;
  array[num_plate] row_vector[4] f_plate_effects_unscaled;
  vector[num_sam_id] xlog_sam;
  array[tot_cat_per_f_pred_var] row_vector[4] f_effects_by_pred_var_cat_unscaled;
  array[tot_cat_per_p_pos_pred_var] real p_pos_effects_by_pred_var_cat_unscaled;
}

transformed parameters{
  
  matrix[tot_cat_per_f_pred_var, 4] f_effects_by_pred_var_cat;
  {
    int cat_current = 1;
    for (f_pred_var in 1:num_f_pred_vars) {
      int num_cat_this_f_pred_var = num_cat_per_f_pred_var[f_pred_var];
      for (cat in cat_current:(cat_current + num_cat_this_f_pred_var - 1)) {
        f_effects_by_pred_var_cat[cat, ] = f_effects_by_pred_var_cat_unscaled[cat] .* 
        sigma_f_pred_vars[f_pred_var]; // element-wise multiplication of 4-vectors
      }
      cat_current += num_cat_this_f_pred_var;
    }
  }

  matrix[num_plate, 4] f_per_plate = rep_matrix(f, num_plate);
  for (plate in 1:num_plate) {
    f_per_plate[plate, ] += f_plate_effects_unscaled[plate] .* sigma_f_plate;
  }
  if (predict_f) {
    f_per_plate += design_matrix_f * f_effects_by_pred_var_cat;
  }
  
  vector[num_plate] f_3_min_f_2_per_plate = f_per_plate[, 3] - f_per_plate[, 2];

  vector[tot_cat_per_p_pos_pred_var] p_pos_effects_by_pred_var_cat;
  {
    int cat_current = 1;
    for (p_pos_pred_var in 1:num_p_pos_pred_vars) {
      int num_cat_this_p_pos_pred_var = num_cat_per_p_pos_pred_var[p_pos_pred_var];
      for (cat in cat_current:(cat_current + num_cat_this_p_pos_pred_var - 1)) {
        p_pos_effects_by_pred_var_cat[cat] = p_pos_effects_by_pred_var_cat_unscaled[cat] * 
        sigma_p_pos_pred_vars[p_pos_pred_var]; 
      }
      cat_current += num_cat_this_p_pos_pred_var;
    }
  }
  
  vector[num_sam_id] p_pos_log_per_sam_id;
  vector[num_sam_id] p_neg_log_per_sam_id;
  if (predict_p_pos) {
    vector[num_sam_id] p_pos_per_sam_id = inv_logit(
    rep_vector(logit(p_pos), num_sam_id) +
    design_matrix_p_pos * p_pos_effects_by_pred_var_cat);
    p_pos_log_per_sam_id = log(  p_pos_per_sam_id);
    p_neg_log_per_sam_id = log1m(p_pos_per_sam_id);
  } else {
    p_pos_log_per_sam_id = rep_vector(log(  p_pos), num_sam_id);
    p_neg_log_per_sam_id = rep_vector(log1m(p_pos), num_sam_id);
  }

  vector[num_cal_tot] y_cal_mean_per_obs;
  vector[num_sam_tot] y_sam_mean_per_obs;
  vector[num_cal_tot] y_obs_sd_cal;
  vector[num_sam_tot] y_obs_sd_sam;
  profile("y_means_and_sds") {
    
  for (cal_rep in 1:num_cal_tot) {
    if (x_cal_is_zero[cal_rep]) {
      y_cal_mean_per_obs[cal_rep] = f_per_plate[which_plate_cal[cal_rep], 2];
      y_obs_sd_cal[cal_rep] = y_obs_sd_cal_min;
    } else {
      int plate = which_plate_cal[cal_rep];
      real denominator = 
      (1 + exp(-f_per_plate[plate, 1] * (xlog_cal[cal_rep] - f_per_plate[plate, 4])));
      y_cal_mean_per_obs[cal_rep] =
      f_per_plate[plate, 2] + f_3_min_f_2_per_plate[plate] / denominator;
      y_obs_sd_cal[cal_rep] = y_obs_sd_cal_min + y_obs_sd_cal_jump / denominator;
    }
  }
  
  for (sam_rep in 1:num_sam_tot) {
    int plate = which_plate_sam[sam_rep];
    real denominator = (1 + exp(-f_per_plate[plate, 1] * (xlog_sam[which_id_sam[sam_rep]] - f_per_plate[plate, 4])));
    y_sam_mean_per_obs[sam_rep] =
    f_per_plate[plate, 2] + f_3_min_f_2_per_plate[plate] / denominator;
    y_obs_sd_sam[sam_rep] =
    y_obs_sd_sam_min + y_obs_sd_sam_jump / denominator;
  }
  }
  
  vector[num_sam_tot] y_sam_loglik_per_obs;
  profile("likelihood_sam") {
  for (sam_rep in 1:num_sam_tot) {
    y_sam_loglik_per_obs[sam_rep] = normal_lpdf(
    y_sam[sam_rep] | y_sam_mean_per_obs[sam_rep], y_obs_sd_sam[sam_rep]);
  }
  }
  
  vector[num_plate] loglik_per_plate;
  for (plate in 1:num_plate) {
    loglik_per_plate[plate] = multi_normal_lpdf(f_plate_effects_unscaled[plate] | zeros_4, rho);
  }
  for (cal_rep in 1:num_cal_tot) {
    loglik_per_plate[which_plate_cal[cal_rep]] +=
    normal_lpdf(y_cal[cal_rep] | y_cal_mean_per_obs[cal_rep], y_obs_sd_cal[cal_rep]);
  }
  
  
}


model {
  
  // Priors
  x_sam_pos_mu ~ uniform(max([x_sam_pos_mu_lower, x_sam_neg_mu]), x_sam_pos_mu_upper);
  rho ~ lkj_corr(rho_prior_eta);
  profile("mixture_model") {
  for (sam_id in 1:num_sam_id) {
    target += log_sum_exp(
      p_pos_log_per_sam_id[sam_id] + normal_lpdf(xlog_sam[sam_id] | x_sam_pos_mu, x_sam_pos_sd),
      p_neg_log_per_sam_id[sam_id] + normal_lpdf(xlog_sam[sam_id] | x_sam_neg_mu, x_sam_neg_sd));
  }
  }
  f_effects_by_pred_var_cat_unscaled ~ multi_normal(zeros_for_f_pred_vars, rho);
  p_pos_effects_by_pred_var_cat_unscaled ~ std_normal();
  
  // Mixed prior and likelihood term, breaking the separation:
  target += sum(loglik_per_plate);
  
  // Likelihood
  //profile("likelihood_cal") {
  //if (sample_posterior_not_prior) {
    target += sum(y_sam_loglik_per_obs);
  //}
  //}
} 

generated quantities {
  
  real y_obs_sd_cal_max = y_obs_sd_cal_min + y_obs_sd_cal_jump;
  real y_obs_sd_sam_max = y_obs_sd_sam_min + y_obs_sd_sam_jump;
  
  array[num_cal_tot] real y_cal_sim = normal_rng(y_cal_mean_per_obs, y_obs_sd_cal);
  array[num_sam_tot] real y_sam_sim = normal_rng(y_sam_mean_per_obs, y_obs_sd_sam);

  vector[num_sam_id] p_sam_is_pos;
  for (sam_id in 1:num_sam_id) {
    real p_log = p_pos_log_per_sam_id[sam_id] +
    normal_lpdf(xlog_sam[sam_id] | x_sam_pos_mu, x_sam_pos_sd);
    p_sam_is_pos[sam_id] = exp(p_log - log_sum_exp(p_log, p_neg_log_per_sam_id[sam_id] +
    normal_lpdf(xlog_sam[sam_id] | x_sam_neg_mu, x_sam_neg_sd)));
  }
  
}
