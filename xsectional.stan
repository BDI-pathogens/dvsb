functions {
  real PL4(real x_log, real f_1, real f_2, real f_3, real f_4) {
    return f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (x_log - f_4)));
  }
}


data {
  // Actual data
  int<lower = 1> num_bat;
  int<lower = num_bat> num_obs_con;
  int<lower = 1, upper = num_bat> which_bat_con[num_obs_con];
  vector[num_obs_con] y_log_con;
  real concentration_log_con[num_obs_con];
  
  // Other things to keep fixed over a complete round of sampling: a binary
  // switch to control whether we sample from the prior or the posterior
  // (important to compare the difference), and upper and lower bounds for the
  // priors.
  int<lower = 0, upper = 1> sample_posterior_not_prior;
  vector[4] f_lower;
  vector[4] f_upper;
  vector[4] sigma_bat_f_lower;
  vector[4] sigma_bat_f_upper;
  real y_log_obs_sd_lower;
  real y_log_obs_sd_upper;
}

transformed data {
  array[num_bat] vector[4] zeros;
  for (i in 1:num_bat) zeros[i] = rep_vector(0, 4);
}

parameters {
  
  // Those constrained only by lower and upper
  vector<lower = f_lower, upper = f_upper>[4] f;
  vector<lower = sigma_bat_f_lower, upper = sigma_bat_f_upper>[4] sigma_bat_f;
  real<lower = y_log_obs_sd_lower, upper = y_log_obs_sd_upper> y_log_obs_sd;
  
  // Those with explicit priors declared
  corr_matrix[4] Rho_bat;
  array[num_bat] vector[4] f_bat_effects_unscaled;
}

transformed parameters{
  
  vector[num_obs_con] y_log_con_mean_per_obs;
  
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
  for (con in 1:num_obs_con) {
    y_log_con_mean_per_obs[con] = 
    PL4(concentration_log_con[con], f_per_bat[which_bat_con[con]][1],
    f_per_bat[which_bat_con[con]][2], f_per_bat[which_bat_con[con]][3],
    f_per_bat[which_bat_con[con]][4]);
  }
  
}


model {
  
  // Priors
  Rho_bat ~ lkj_corr(3);
  f_bat_effects_unscaled ~ multi_normal(zeros, Rho_bat);
  
  // Likelihood
  if (sample_posterior_not_prior) {
    y_log_con ~ normal(y_log_con_mean_per_obs, y_log_obs_sd);
  }
} 

generated quantities {
  real y_log_con_sim[num_obs_con] = normal_rng(y_log_con_mean_per_obs, y_log_obs_sd);
}
