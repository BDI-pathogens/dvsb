test_that("wrangling true params works on simulated data", {
  
  priors_dir <- system.file("input_priors", package = "dvsb")
  priors_path <- file.path(priors_dir, "priors.csv")
  priors_list <- read_priors(priors_path)
  df_priors_scalars <- priors_list$df_priors_scalars
  df_priors_vectors <- priors_list$df_priors_vectors
  rho_prior_eta <- priors_list$rho_prior_eta
  
  # Simulate lots of sams to ensure no categories are stochastically absent
  data <- simulate_data(num_plate = 100, num_sam_per_plate = 100)
  df_sam <- data$df_sam
  df_cal <- data$df_cal
  df_plate <- data$df_plate
  param_true_values_list <- data$params
  f_pred_vars_names <- data$f_pred_vars_names
  x_mix_pred_vars_names <- data$x_mix_pred_vars_names
  p_pos_binary_pred_vars <- data$p_pos_binary_pred_vars
  
  data_wrangled <- prepare_data_for_stan(
  df_sam = df_sam, 
  df_cal = df_cal, 
  df_priors_scalars = df_priors_scalars,
  df_priors_vectors = df_priors_vectors, 
  rho_prior_eta = rho_prior_eta, 
  x_mix_pred_vars_names = x_mix_pred_vars_names, 
  p_pos_binary_pred_vars = p_pos_binary_pred_vars,
  f_pred_vars_names = f_pred_vars_names
  )  
  
  wrangled_names <- names(
  wrangle_true_params(param_true_values_list = data$params,
                      data_descriptors = data_wrangled$data_descriptors))
  expect_true(identical(sort(wrangled_names),
                         sort(c(
  "f[1]", "f[2]", "f[3]", "f[4]", "sigma_f[1]_plate", "sigma_f[2]_plate", "sigma_f[3]_plate", "sigma_f[4]_plate", "rho[1,2]", "rho[1,3]", "rho[1,4]", "rho[2,3]", "rho[2,4]", "rho[3,4]", "sd_neg", "sd_pos", "mu_neg", "mu_pos", "p_pos",  "p_blank",  "y_obs_sd_cal_min", "y_obs_sd_cal_jump", "y_obs_sd_sam_min", "y_obs_sd_sam_jump", "y_obs_sd_min_log_shift_sd", "y_obs_sd_jump_log_shift_sd", "p_pos_effect_lettera", "p_pos_effect_letterb", "p_pos_effect_letterc", "p_pos_effect_letterd", "sigma_p_pos_pred_vars_letter", "mu_neg_effect_lettera", "mu_neg_effect_letterb", "mu_neg_effect_letterc", "mu_neg_effect_letterd", "sigma_mu_neg_pred_vars_letter", "mu_pos_effect_lettera", "mu_pos_effect_letterb", "mu_pos_effect_letterc", "mu_pos_effect_letterd", "sigma_mu_pos_pred_vars_letter", "sd_neg_effect_lettera", "sd_neg_effect_letterb", "sd_neg_effect_letterc", "sd_neg_effect_letterd", "sigma_sd_neg_pred_vars_letter", "sd_pos_effect_lettera", "sd_pos_effect_letterb", "sd_pos_effect_letterc", "sd_pos_effect_letterd", "sigma_sd_pos_pred_vars_letter", "p_pos_effect_boolA", "p_pos_effect_boolB", "p_pos_effect_boolC"
  ))))
})
