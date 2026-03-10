// See the associated R file for explanations and definitions of abbreviations.

data {
  // Actual data
  int<lower = 0> num_plate;
  int<lower = 0> num_cal_tot;
  int<lower = 0> num_sam_id;
  int<lower = num_sam_id> num_sam_rep;
  array[num_cal_tot] int<lower = 1, upper = num_plate> which_plate_cal;
  array[num_sam_rep] int<lower = 1, upper = num_plate> which_plate_sam;
  array[num_sam_rep] int<lower = 1, upper = num_sam_id> which_id_sam;
  vector[num_cal_tot] y_cal;
  vector<lower = 0>[num_cal_tot] x_cal;
  vector[num_sam_rep] y_sam;
  int<lower = 0> num_f_pred_vars;
  array[num_f_pred_vars] int<lower = 2> num_cat_per_f_pred_var;
  matrix<lower = 0, upper = 1>[num_plate, sum(num_cat_per_f_pred_var)] design_matrix_f;
  int<lower = 0> num_p_pos_pred_vars;
  array[num_p_pos_pred_vars] int<lower = 2> num_cat_per_p_pos_pred_var;
  matrix<lower = 0, upper = 1>[num_sam_id, sum(num_cat_per_p_pos_pred_var)] design_matrix_p_pos;
  int<lower = 0, upper = 1> num_mu_pos_pred_vars; // upper = 1 for now, for non-centered parameterisation
  array[num_mu_pos_pred_vars] int<lower = 2> num_cat_per_mu_pos_pred_var;
  matrix<lower = 0, upper = 1>[num_sam_id, sum(num_cat_per_mu_pos_pred_var)] design_matrix_mu_pos;
  int<lower = 0> num_sd_pos_pred_vars;
  array[num_sd_pos_pred_vars] int<lower = 2> num_cat_per_sd_pos_pred_var;
  matrix<lower = 0, upper = 1>[num_sam_id, sum(num_cat_per_sd_pos_pred_var)] design_matrix_sd_pos;
  int<lower = 0, upper = 1> num_mu_neg_pred_vars; // upper = 1 for now, for non-centered parameterisation
  array[num_mu_neg_pred_vars] int<lower = 2> num_cat_per_mu_neg_pred_var;
  matrix<lower = 0, upper = 1>[num_sam_id, sum(num_cat_per_mu_neg_pred_var)] design_matrix_mu_neg;
  int<lower = 0> num_sd_neg_pred_vars;
  array[num_sd_neg_pred_vars] int<lower = 2> num_cat_per_sd_neg_pred_var;
  matrix<lower = 0, upper = 1>[num_sam_id, sum(num_cat_per_sd_neg_pred_var)] design_matrix_sd_neg;
  int<lower = 0> num_p_pos_binary_pred_vars;
  matrix<lower = 0, upper = 1>[num_sam_id, num_p_pos_binary_pred_vars] design_matrix_p_pos_binary;
  

  // Other things to keep fixed over a complete round of sampling: a binary
  // switch to control whether we sample from the prior or the posterior
  // (important to compare the difference), and upper and lower bounds for the priors.
  int<lower = 0, upper = 1> sample_posterior_not_prior;
  row_vector[4] f_lower;
  row_vector<lower = f_lower>[4] f_upper;
  row_vector[4] sigma_f_plate_lower;
  row_vector<lower = sigma_f_plate_lower>[4] sigma_f_plate_upper;
  row_vector[4] sigma_f_pred_vars_lower;
  row_vector<lower = sigma_f_pred_vars_lower>[4] sigma_f_pred_vars_upper;
  real sigma_p_pos_pred_vars_lower;
  real<lower = sigma_p_pos_pred_vars_lower> sigma_p_pos_pred_vars_upper;
  real sigma_mu_pos_pred_vars_lower;
  real<lower = sigma_mu_pos_pred_vars_lower> sigma_mu_pos_pred_vars_upper;
  real sigma_sd_pos_pred_vars_lower;
  real<lower = sigma_sd_pos_pred_vars_lower> sigma_sd_pos_pred_vars_upper;
  real sigma_mu_neg_pred_vars_lower;
  real<lower = sigma_mu_neg_pred_vars_lower> sigma_mu_neg_pred_vars_upper;
  real sigma_sd_neg_pred_vars_lower;
  real<lower = sigma_sd_neg_pred_vars_lower> sigma_sd_neg_pred_vars_upper;
  real y_obs_sd_cal_min_lower;
  real sigma_p_pos_binary_pred_vars_lower;
  real<lower = sigma_p_pos_binary_pred_vars_lower> sigma_p_pos_binary_pred_vars_upper;
  real<lower = y_obs_sd_cal_min_lower> y_obs_sd_cal_min_upper;
  real y_obs_sd_cal_jump_lower;
  real<lower = y_obs_sd_cal_jump_lower> y_obs_sd_cal_jump_upper;
  real y_obs_sd_sam_min_lower;
  real<lower = y_obs_sd_sam_min_lower> y_obs_sd_sam_min_upper;
  real y_obs_sd_sam_jump_lower;
  real<lower = y_obs_sd_sam_jump_lower> y_obs_sd_sam_jump_upper;
  
  real mu_neg_lower;
  real<lower = mu_neg_lower> mu_neg_upper;
  real<lower = mu_neg_lower> mu_pos_lower;
  real<lower = mu_pos_lower> mu_pos_upper;
  real<lower = 0> sd_neg_lower;
  real<lower = sd_neg_lower> sd_neg_upper;
  real<lower = 0> sd_pos_lower;
  real<lower = sd_pos_lower> sd_pos_upper;

  real p_pos_lower;
  real<lower = p_pos_lower> p_pos_upper;
  real p_blank_lower;
  real<lower = p_blank_lower> p_blank_upper;
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
  int tot_cat_per_mu_pos_pred_var = sum(num_cat_per_mu_pos_pred_var);
  int predict_mu_pos = 1 ? num_mu_pos_pred_vars > 0 : 0;
  int tot_cat_per_sd_pos_pred_var = sum(num_cat_per_sd_pos_pred_var);
  int predict_sd_pos = 1 ? num_sd_pos_pred_vars > 0 : 0;
  int tot_cat_per_mu_neg_pred_var = sum(num_cat_per_mu_neg_pred_var);
  int predict_mu_neg = 1 ? num_mu_neg_pred_vars > 0 : 0;
  int tot_cat_per_sd_neg_pred_var = sum(num_cat_per_sd_neg_pred_var);
  int predict_sd_neg = 1 ? num_sd_neg_pred_vars > 0 : 0;
  int predict_p_pos_binary = 1 ? num_p_pos_binary_pred_vars > 0 : 0;
  
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
  array[num_mu_pos_pred_vars] real<lower = sigma_mu_pos_pred_vars_lower,
  upper = sigma_mu_pos_pred_vars_upper> sigma_mu_pos_pred_vars;
  array[num_sd_pos_pred_vars] real<lower = sigma_sd_pos_pred_vars_lower,
  upper = sigma_sd_pos_pred_vars_upper> sigma_sd_pos_pred_vars;
  array[num_mu_neg_pred_vars] real<lower = sigma_mu_neg_pred_vars_lower,
  upper = sigma_mu_neg_pred_vars_upper> sigma_mu_neg_pred_vars;
  array[num_sd_neg_pred_vars] real<lower = sigma_sd_neg_pred_vars_lower,
  upper = sigma_sd_neg_pred_vars_upper> sigma_sd_neg_pred_vars;
  array[predict_p_pos_binary] real<lower = sigma_p_pos_binary_pred_vars_lower,
  upper = sigma_p_pos_binary_pred_vars_upper> sigma_p_pos_binary_pred_vars;
  
  // Scalar params constrained only by lower and upper
  real<lower = mu_neg_lower, upper = mu_neg_upper> mu_neg;
  real<lower = sd_neg_lower, upper = sd_neg_upper> sd_neg;
  real<lower = sd_pos_lower, upper = sd_pos_upper> sd_pos;
  real<lower = p_pos_lower,    upper = p_pos_upper>    p_pos;
  real<lower = p_blank_lower,  upper = p_blank_upper>  p_blank;
  real<lower = y_obs_sd_cal_min_lower,  upper = y_obs_sd_cal_min_upper>  y_obs_sd_cal_min;
  real<lower = y_obs_sd_cal_jump_lower, upper = y_obs_sd_cal_jump_upper> y_obs_sd_cal_jump;
  real<lower = y_obs_sd_sam_min_lower,  upper = y_obs_sd_sam_min_upper>  y_obs_sd_sam_min;
  real<lower = y_obs_sd_sam_jump_lower, upper = y_obs_sd_sam_jump_upper> y_obs_sd_sam_jump;
  
  // Enforce that mu_pos > mu_neg
  real<lower = max([mu_pos_lower, mu_neg]), upper = mu_pos_upper> mu_pos;

  
  // Those with explicit priors declared
  corr_matrix[4] rho;
  array[num_plate] row_vector[4] f_plate_effects_unscaled;
  vector[num_sam_id] xlog_sam;
  array[tot_cat_per_f_pred_var] row_vector[4] f_effects_by_pred_var_cat_unscaled;
  array[tot_cat_per_p_pos_pred_var]  real  p_pos_effects_by_pred_var_cat_unscaled;
  array[tot_cat_per_sd_pos_pred_var] real sd_pos_effects_by_pred_var_cat_unscaled;
  array[tot_cat_per_sd_neg_pred_var] real sd_neg_effects_by_pred_var_cat_unscaled;
  vector[num_p_pos_binary_pred_vars] p_pos_binary_effects_by_pred_var_unscaled;
  vector<lower = (mu_neg - mu_pos) / (2 * sigma_mu_pos_pred_vars[1])>[tot_cat_per_mu_pos_pred_var]
  mu_pos_effects_by_pred_var_cat_unscaled;
  vector<upper = (mu_pos - mu_neg) / (2 * sigma_mu_neg_pred_vars[1])>[tot_cat_per_mu_neg_pred_var]
  mu_neg_effects_by_pred_var_cat_unscaled;
}

transformed parameters{
  
  real p_blank_log   = log(  p_blank);
  real p_blank_log1m = log1m(p_blank);
  
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

  vector[tot_cat_per_p_pos_pred_var]  p_pos_effects_by_pred_var_cat;
  vector[tot_cat_per_sd_pos_pred_var] sd_pos_effects_by_pred_var_cat;
  vector[tot_cat_per_sd_neg_pred_var] sd_neg_effects_by_pred_var_cat;
  vector[tot_cat_per_mu_pos_pred_var] mu_pos_effects_by_pred_var_cat;
  vector[tot_cat_per_mu_neg_pred_var] mu_neg_effects_by_pred_var_cat;
  vector[num_p_pos_binary_pred_vars]  p_pos_binary_effects_by_pred_var;

  vector[num_sam_id] p_pos_log_per_sam_id;
  vector[num_sam_id] p_neg_log_per_sam_id;
  vector[num_sam_id] p_pos_per_sam_id;
  vector[num_sam_id] mu_pos_per_sam_id = rep_vector(mu_pos, num_sam_id);
  vector[num_sam_id] mu_neg_per_sam_id = rep_vector(mu_neg, num_sam_id);
  vector[num_sam_id] sd_pos_per_sam_id;
  vector[num_sam_id] sd_neg_per_sam_id;

  profile("define_x_mix_pred_vars") {
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

    cat_current = 1;
    for (sd_pos_pred_var in 1:num_sd_pos_pred_vars) {
      int num_cat = num_cat_per_sd_pos_pred_var[sd_pos_pred_var];
      for (cat in cat_current:(cat_current + num_cat - 1)) {
        sd_pos_effects_by_pred_var_cat[cat] = sd_pos_effects_by_pred_var_cat_unscaled[cat] * 
        sigma_sd_pos_pred_vars[sd_pos_pred_var]; 
      }
      cat_current += num_cat;
    }
    cat_current = 1;
    for (sd_neg_pred_var in 1:num_sd_neg_pred_vars) {
      int num_cat = num_cat_per_sd_neg_pred_var[sd_neg_pred_var];
      for (cat in cat_current:(cat_current + num_cat - 1)) {
        sd_neg_effects_by_pred_var_cat[cat] = sd_neg_effects_by_pred_var_cat_unscaled[cat] * 
        sigma_sd_neg_pred_vars[sd_neg_pred_var]; 
      }
      cat_current += num_cat;
    }
  }
  
  if (predict_p_pos || predict_p_pos_binary) {
    p_pos_per_sam_id = rep_vector(logit(p_pos), num_sam_id);
    if (predict_p_pos) {
      p_pos_per_sam_id += design_matrix_p_pos * p_pos_effects_by_pred_var_cat;
    }
    if (predict_p_pos_binary) {
      p_pos_binary_effects_by_pred_var =
      p_pos_binary_effects_by_pred_var_unscaled * sigma_p_pos_binary_pred_vars[1];
      p_pos_per_sam_id +=
      design_matrix_p_pos_binary * p_pos_binary_effects_by_pred_var;
    }
    p_pos_per_sam_id = inv_logit(p_pos_per_sam_id);
    p_pos_log_per_sam_id = log(  p_pos_per_sam_id);
    p_neg_log_per_sam_id = log1m(p_pos_per_sam_id);
  } else {
    p_pos_per_sam_id     = rep_vector(p_pos,        num_sam_id);
    p_pos_log_per_sam_id = rep_vector(log(  p_pos), num_sam_id);
    p_neg_log_per_sam_id = rep_vector(log1m(p_pos), num_sam_id);
  }
  if (predict_mu_pos) {
    mu_pos_effects_by_pred_var_cat =
    mu_pos_effects_by_pred_var_cat_unscaled * sigma_mu_pos_pred_vars[1];
    mu_pos_per_sam_id += design_matrix_mu_pos * mu_pos_effects_by_pred_var_cat;
  } 
  if (predict_mu_neg) {
    mu_neg_effects_by_pred_var_cat =
    mu_neg_effects_by_pred_var_cat_unscaled * sigma_mu_neg_pred_vars[1];
    mu_neg_per_sam_id += design_matrix_mu_neg * mu_neg_effects_by_pred_var_cat;
  } 
  if (predict_sd_pos) {
    sd_pos_per_sam_id = exp(rep_vector(log(sd_pos), num_sam_id) +
    design_matrix_sd_pos * sd_pos_effects_by_pred_var_cat);
  } else {
    sd_pos_per_sam_id = rep_vector(sd_pos, num_sam_id);
  }
  if (predict_sd_neg) {
    sd_neg_per_sam_id = exp(rep_vector(log(sd_neg), num_sam_id) +
    design_matrix_sd_neg * sd_neg_effects_by_pred_var_cat);
  } else {
    sd_neg_per_sam_id = rep_vector(sd_neg, num_sam_id);
  }
  }

  vector[num_cal_tot] y_cal_mean_per_obs;
  vector[num_sam_rep] y_sam_mean_per_obs;
  vector[num_cal_tot] y_obs_sd_cal;
  vector[num_sam_rep] y_obs_sd_sam;
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
  for (sam_rep in 1:num_sam_rep) {
    int plate = which_plate_sam[sam_rep];
    real denominator = (1 + exp(-f_per_plate[plate, 1] * (xlog_sam[which_id_sam[sam_rep]] - f_per_plate[plate, 4])));
    y_sam_mean_per_obs[sam_rep] =
    f_per_plate[plate, 2] + f_3_min_f_2_per_plate[plate] / denominator;
    y_obs_sd_sam[sam_rep] =
    y_obs_sd_sam_min + y_obs_sd_sam_jump / denominator;
  }
  }
  
  //vector[num_sam_rep] y_sam_loglik_per_obs_from_blank;
  //vector[num_sam_rep] y_sam_loglik_per_obs_from_notblank;
  vector[num_sam_rep] y_sam_loglik_per_obs;
  profile("likelihood_sam") {
  for (sam_rep in 1:num_sam_rep) {
    real y_sam_loglik_per_obs_from_blank = p_blank_log + normal_lpdf(
    y_sam[sam_rep] | f_per_plate[which_plate_sam[sam_rep], 2], y_obs_sd_sam_min);
    real y_sam_loglik_per_obs_from_notblank = p_blank_log1m + normal_lpdf(
    y_sam[sam_rep] | y_sam_mean_per_obs[sam_rep], y_obs_sd_sam[sam_rep]);
    y_sam_loglik_per_obs[sam_rep] =
    log_sum_exp(y_sam_loglik_per_obs_from_blank,
    y_sam_loglik_per_obs_from_notblank);
  }
  }
  
  vector[num_plate] loglik_per_plate_notblanks = rep_vector(0, num_plate);
  vector[num_plate] loglik_per_plate_blanks = rep_vector(0, num_plate);
  vector[num_plate] logprob_f_effects_per_plate;
  for (plate in 1:num_plate) {
    logprob_f_effects_per_plate[plate] =
    multi_normal_lpdf(f_plate_effects_unscaled[plate] | zeros_4, rho);
  }
  for (cal_rep in 1:num_cal_tot) {
    if (x_cal_is_zero[cal_rep]) {
      loglik_per_plate_blanks[which_plate_cal[cal_rep]] +=
      normal_lpdf(y_cal[cal_rep] | y_cal_mean_per_obs[cal_rep], y_obs_sd_cal[cal_rep]);
    } else {
      loglik_per_plate_notblanks[which_plate_cal[cal_rep]] +=
      normal_lpdf(y_cal[cal_rep] | y_cal_mean_per_obs[cal_rep], y_obs_sd_cal[cal_rep]);  
    }
  }
  
}


model {
  
  // Priors
  mu_pos ~ uniform(max([mu_pos_lower, mu_neg]), mu_pos_upper);
  rho ~ lkj_corr(rho_prior_eta);
  profile("mixture_model") {
  for (sam_id in 1:num_sam_id) {
    target += log_sum_exp(
      p_pos_log_per_sam_id[sam_id] + normal_lpdf(xlog_sam[sam_id] |
      mu_pos_per_sam_id[sam_id], sd_pos_per_sam_id[sam_id]),
      p_neg_log_per_sam_id[sam_id] + normal_lpdf(xlog_sam[sam_id] |
      mu_neg_per_sam_id[sam_id], sd_neg_per_sam_id[sam_id]));
  }
  }
  f_effects_by_pred_var_cat_unscaled ~ multi_normal(zeros_for_f_pred_vars, rho);
  p_pos_effects_by_pred_var_cat_unscaled  ~ std_normal();
  sd_pos_effects_by_pred_var_cat_unscaled ~ std_normal();
  sd_neg_effects_by_pred_var_cat_unscaled ~ std_normal();
  p_pos_binary_effects_by_pred_var_unscaled ~ std_normal();
  if (predict_mu_pos) {
    mu_pos_effects_by_pred_var_cat_unscaled ~
    std_normal() T[(mu_neg - mu_pos) / (2 * sigma_mu_pos_pred_vars[1]), ];
  }
  if (predict_mu_neg) {
    mu_neg_effects_by_pred_var_cat_unscaled ~ 
    std_normal() T[, (mu_pos - mu_neg) / (2 * sigma_mu_neg_pred_vars[1])]; 
  }
  target += sum(logprob_f_effects_per_plate);
  
  // Likelihood
  if (sample_posterior_not_prior) {
    target += sum(loglik_per_plate_blanks);
    target += sum(loglik_per_plate_notblanks);
    target += sum(y_sam_loglik_per_obs);
  }

} 

// 'sim' is short for simulated, with a fresh draw of stochastic uncertainty 
generated quantities {
  
  real y_obs_sd_cal_max = y_obs_sd_cal_min + y_obs_sd_cal_jump;
  real y_obs_sd_sam_max = y_obs_sd_sam_min + y_obs_sd_sam_jump;
  
  array[num_cal_tot] real y_cal_sim = normal_rng(y_cal_mean_per_obs, y_obs_sd_cal);
  
  // (un)conditional refers to that sam's observed y values.
  // We always condition on population-level parameters and any x mix pred vars.
  array[num_sam_rep] real y_sam_sim_conditional;
  vector[num_sam_id] xlog_sam_sim_unconditional;
  vector[num_sam_rep]   y_sam_sim_unconditional;
  profile("simulation") {
      y_sam_sim_conditional = normal_rng(y_sam_mean_per_obs, y_obs_sd_sam);
  //array[num_sam_rep] real p_sam_rep_is_blank;
  //for (sam_rep in 1:num_sam_rep) {
    //real p_sam_rep_is_blank = exp(y_sam_loglik_per_obs_from_blank[sam_rep] - 
    //log_sum_exp(y_sam_loglik_per_obs_from_blank[sam_rep],
    //y_sam_loglik_per_obs_from_notblank[sam_rep]));
    //if (bernoulli_rng(p_sam_rep_is_blank)) {
    //  y_sam_sim[sam_rep] = normal_rng(f_per_plate[which_plate_sam[sam_rep], 2], y_obs_sd_sam_min);
    //} else {
    //  y_sam_sim[sam_rep] = normal_rng(y_sam_mean_per_obs[sam_rep], y_obs_sd_sam[sam_rep]);
    //}
  //}
  
  for (sam_id in 1:num_sam_id) {
    if (bernoulli_rng(p_pos_per_sam_id[sam_id])) {
      xlog_sam_sim_unconditional[sam_id] = 
      normal_rng(mu_pos_per_sam_id[sam_id], sd_pos_per_sam_id[sam_id]);
    } else {
      xlog_sam_sim_unconditional[sam_id] = 
      normal_rng(mu_neg_per_sam_id[sam_id], sd_neg_per_sam_id[sam_id]);
    }
  }
  for (sam_rep in 1:num_sam_rep) {
    int plate = which_plate_sam[sam_rep];
    if (bernoulli_rng(p_blank)) {
      y_sam_sim_unconditional[sam_rep] = normal_rng(
      f_per_plate[which_plate_sam[sam_rep], 2], y_obs_sd_sam_min);
    } else {
      real denominator = (1 + exp(-f_per_plate[plate, 1] *
      (xlog_sam_sim_unconditional[which_id_sam[sam_rep]] - f_per_plate[plate, 4])));
      y_sam_sim_unconditional[sam_rep] = normal_rng(
        f_per_plate[plate, 2] + f_3_min_f_2_per_plate[plate] / denominator,
        y_obs_sd_sam_min + y_obs_sd_sam_jump / denominator);
    }
  }
  }

  vector[num_sam_id] p_sam_is_pos;
  for (sam_id in 1:num_sam_id) {
    real p_log = p_pos_log_per_sam_id[sam_id] +
    normal_lpdf(xlog_sam[sam_id] | mu_pos_per_sam_id[sam_id], sd_pos_per_sam_id[sam_id]);
    p_sam_is_pos[sam_id] = exp(p_log - log_sum_exp(p_log, p_neg_log_per_sam_id[sam_id] +
    normal_lpdf(xlog_sam[sam_id] | mu_neg_per_sam_id[sam_id], sd_neg_per_sam_id[sam_id])));
  }
  
}
