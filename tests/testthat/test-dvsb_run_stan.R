test_that("cmdstan works", {
  
  priors_path <- file.path(system.file("input_priors", package = "dvsb"), "priors.csv")
  priors_list <- read_priors(priors_path)
  
  data <- simulate_data(num_plate = 2, 
                        num_sam_per_plate = 20, 
                        num_rep_per_sam = 1,
                        num_rep_per_cal = 1, 
                        x_cals = c(0, 1))
  
  data_wrangled <- prepare_data_for_stan(
    df_sam = data$df_sam, 
    df_cal = data$df_cal, 
    df_priors_scalars = priors_list$df_priors_scalars,
    df_priors_vectors = priors_list$df_priors_vectors, 
    rho_prior_eta = priors_list$rho_prior_eta, 
    x_mix_pred_vars_names = data$x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = data$p_pos_binary_pred_vars,
    f_pred_vars_names = data$f_pred_vars_names
  )
  
  df_sam <- data_wrangled$df_sam
  df_cal <- data_wrangled$df_cal
  df_plate <- dplyr::left_join(data$df_plate, data_wrangled$df_plate, by = "plate") 
  
  iterations <- 1 
  cmdstan_output_basename <- "/Users/cwymant/temp_dvsb_out"
  
  df_posterior <- pkgcond::suppress_warnings(run_stan(
    input_to_stan = data_wrangled$stan_input_posterior,
    data_descriptors = data_wrangled$data_descriptors,
    interface = "cmdstan",
    iter_warmup = iterations,
    iter_sampling = iterations,
    cmdstan_path_to_output = cmdstan_output_basename),
    mastiff:::stan_safe_warnings())
  
  df_prior <- pkgcond::suppress_warnings(run_stan(
    input_to_stan = data_wrangled$stan_input_prior,
    data_descriptors = data_wrangled$data_descriptors,
    interface = "cmdstan",
    iter_warmup = iterations,
    iter_sampling = iterations,
    cmdstan_path_to_output = cmdstan_output_basename),
    mastiff:::stan_safe_warnings())
  
  expect_true(is.data.frame(df_posterior))
  expect_true(is.data.frame(df_prior))
  expect_all_true(stan_params_expected_in_testing_posterior %in% names(df_posterior))
  expect_all_true(stan_params_expected_in_testing_prior %in% names(df_prior))
  
})



test_that("cmdstanr works", {
  
  priors_path <- file.path(system.file("input_priors", package = "dvsb"), "priors.csv")
  priors_list <- read_priors(priors_path)
  
  data <- simulate_data(num_plate = 2, 
                        num_sam_per_plate = 20, 
                        num_rep_per_sam = 1,
                        num_rep_per_cal = 1, 
                        x_cals = c(0, 1))
  
  data_wrangled <- prepare_data_for_stan(
    df_sam = data$df_sam, 
    df_cal = data$df_cal, 
    df_priors_scalars = priors_list$df_priors_scalars,
    df_priors_vectors = priors_list$df_priors_vectors, 
    rho_prior_eta = priors_list$rho_prior_eta, 
    x_mix_pred_vars_names = data$x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = data$p_pos_binary_pred_vars,
    f_pred_vars_names = data$f_pred_vars_names
  )
  
  df_sam <- data_wrangled$df_sam
  df_cal <- data_wrangled$df_cal
  df_plate <- dplyr::left_join(data$df_plate, data_wrangled$df_plate, by = "plate") 
  
  iterations <- 1 
  
  df_posterior <- pkgcond::suppress_warnings(run_stan(
    input_to_stan = data_wrangled$stan_input_posterior,
    data_descriptors = data_wrangled$data_descriptors,
    interface = "cmdstanr",
    iter_warmup = iterations,
    iter_sampling = iterations),
    mastiff:::stan_safe_warnings())
  
  df_prior <- pkgcond::suppress_warnings(run_stan(
    input_to_stan = data_wrangled$stan_input_prior,
    data_descriptors = data_wrangled$data_descriptors,
    interface = "cmdstanr",
    iter_warmup = iterations,
    iter_sampling = iterations),
    mastiff:::stan_safe_warnings())
  
  expect_true(is.data.frame(df_posterior))
  expect_true(is.data.frame(df_prior))
  expect_all_true(stan_params_expected_in_testing_posterior %in% names(df_posterior))
  expect_all_true(stan_params_expected_in_testing_prior %in% names(df_prior))
  
})


test_that("rstan works", {
  
  priors_path <- file.path(system.file("input_priors", package = "dvsb"), "priors.csv")
  priors_list <- read_priors(priors_path)
  
  data <- simulate_data(num_plate = 2, 
                        num_sam_per_plate = 20, 
                        num_rep_per_sam = 1,
                        num_rep_per_cal = 1, 
                        x_cals = c(0, 1))
  
  data_wrangled <- prepare_data_for_stan(
    df_sam = data$df_sam, 
    df_cal = data$df_cal, 
    df_priors_scalars = priors_list$df_priors_scalars,
    df_priors_vectors = priors_list$df_priors_vectors, 
    rho_prior_eta = priors_list$rho_prior_eta, 
    x_mix_pred_vars_names = data$x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = data$p_pos_binary_pred_vars,
    f_pred_vars_names = data$f_pred_vars_names
  )
  
  df_sam <- data_wrangled$df_sam
  df_cal <- data_wrangled$df_cal
  df_plate <- dplyr::left_join(data$df_plate, data_wrangled$df_plate, by = "plate") 
  
  iterations <- 1 
  
  df_posterior <- pkgcond::suppress_warnings(run_stan(
    input_to_stan = data_wrangled$stan_input_posterior,
    data_descriptors = data_wrangled$data_descriptors,
    interface = "rstan",
    iter_warmup = iterations,
    iter_sampling = iterations,
    cores = 1),
    mastiff:::stan_safe_warnings())
  
  df_prior <- pkgcond::suppress_warnings(run_stan(
    input_to_stan = data_wrangled$stan_input_prior,
    data_descriptors = data_wrangled$data_descriptors,
    interface = "rstan",
    iter_warmup = iterations,
    iter_sampling = iterations,
    cores = 1),
    mastiff:::stan_safe_warnings())
  
  expect_true(is.data.frame(df_posterior))
  expect_true(is.data.frame(df_prior))
  expect_all_true(stan_params_expected_in_testing_posterior %in% names(df_posterior))
  expect_all_true(stan_params_expected_in_testing_prior %in% names(df_prior))
  
})

stan_params_expected_in_testing_posterior <- c(
  "f[1]", "f[2]", "f[3]", "f[4]", "sigma_f[1]_plate", "sigma_f[2]_plate", 
  "sigma_f[3]_plate", "sigma_f[4]_plate", "sigma_p_pos_pred_vars_letter", 
  "sigma_mu_pos_pred_vars_letter", "sigma_sd_pos_pred_vars_letter", 
  "sigma_mu_neg_pred_vars_letter", "sigma_sd_neg_pred_vars_letter", "mu_neg",
  "sd_neg", "sd_pos", "p_pos", "y_obs_sd_cal_min", "y_obs_sd_cal_jump",
  "y_obs_sd_sam_min", "y_obs_sd_sam_jump", "mu_pos", "rho[1,2]", "rho[1,3]",
  "rho[2,3]", "rho[1,4]", "rho[2,4]", "rho[3,4]", "xlog_sam[1]", "xlog_sam[2]",
  "xlog_sam[3]", "xlog_sam[4]", "p_pos_effect_boolA", "p_pos_effect_boolB", 
  "p_pos_effect_boolC", "f_per_plate[1,1]", "f_per_plate[2,1]", 
  "f_per_plate[1,2]", "f_per_plate[2,2]", "f_per_plate[1,3]",
  "f_per_plate[2,3]", "f_per_plate[1,4]", "f_per_plate[2,4]", 
  "p_pos_effect_letterd", "p_pos_effect_letterb", "p_pos_effect_letterc",
  "sd_pos_effect_letterd", "sd_pos_effect_letterb", "sd_pos_effect_letterc",
  "sd_neg_effect_letterd", "sd_neg_effect_letterb", "sd_neg_effect_letterc",
  "mu_pos_effect_letterd", "mu_pos_effect_letterb", "mu_pos_effect_letterc",
  "mu_neg_effect_letterd", "mu_neg_effect_letterb", "mu_neg_effect_letterc",
  "y_sam_loglik_per_obs[1]", "y_sam_loglik_per_obs[2]",
  "y_sam_loglik_per_obs[3]", "y_sam_loglik_per_obs[4]", 
  "loglik_per_plate_notblanks[1]", "loglik_per_plate_notblanks[2]",
  "loglik_per_plate_blanks[1]", "loglik_per_plate_blanks[2]", 
  "logprob_f_effects_per_plate[1]", "logprob_f_effects_per_plate[2]", 
  "y_cal_sim[1]", "y_cal_sim[2]", "y_cal_sim[3]", "y_cal_sim[4]", 
  "y_sam_sim_conditional[1]", "y_sam_sim_conditional[2]", 
  "y_sam_sim_conditional[3]", "y_sam_sim_conditional[4]", 
  "pos_sam_sim_unconditional[1]", "pos_sam_sim_unconditional[2]",
  "pos_sam_sim_unconditional[3]", "pos_sam_sim_unconditional[4]", 
  "xlog_sam_sim_unconditional[1]", "xlog_sam_sim_unconditional[2]", 
  "xlog_sam_sim_unconditional[3]", "xlog_sam_sim_unconditional[4]", 
  "y_sam_sim_unconditional[1]", "y_sam_sim_unconditional[2]", 
  "y_sam_sim_unconditional[3]", "y_sam_sim_unconditional[4]", "p_sam_is_pos[1]",
  "p_sam_is_pos[2]", "p_sam_is_pos[3]", "p_sam_is_pos[4]")

stan_params_expected_in_testing_prior <- 
  stan_params_expected_in_testing_posterior[
    !grepl(paste0("y_sam|p_sam_is_pos|pos_sam_sim_unconditional|xlog_sam",
                  "|y_cal_sim|f_per_plate|loglik_per_plate|logprob_f_effects_per_plate"),
    stan_params_expected_in_testing_posterior)]
