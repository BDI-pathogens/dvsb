test_that("cmdstan works", {
  
  stan_path <- file.path(system.file("stan", package = "dvsb"), "dvsb.stan")
  priors_path <- file.path(system.file("input_priors", package = "dvsb"), "priors.csv")
  priors_list <- read_priors(priors_path)
  
  data <- simulate_data(num_plate = 2, num_sam_per_plate = 2)
  
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
  
  iterations <- 10 
  cmdstan_installation_dir <- "/Users/cwymant/.cmdstan/cmdstan-2.37.0/"
  cmdstan_temp_file <- "/Users/cwymant/temp_dvsb.json" 
  cmdstan_output_basename <- "/Users/cwymant/temp_dvsb_out"
  
  df_posterior <- run_stan_interfaces(
    input_to_stan = data_wrangled$stan_input_posterior,
    path_to_stan_code = stan_path,
    interface = "cmdstan",
    iter_warmup = iterations,
    iter_sampling = iterations,
    cmdstan_path_to_installation = cmdstan_installation_dir, 
    cmdstan_path_to_json = cmdstan_temp_file, 
    cmdstan_overwrite_json = TRUE,
    cmdstan_path_to_output = cmdstan_output_basename)
  
  df_prior <- run_stan_interfaces(
    input_to_stan = data_wrangled$stan_input_prior,
    path_to_stan_code = stan_path,
    interface = "cmdstan",
    iter_warmup = iterations,
    iter_sampling = iterations,
    cmdstan_path_to_installation = cmdstan_installation_dir, 
    cmdstan_path_to_json = cmdstan_temp_file, 
    cmdstan_overwrite_json = TRUE,
    cmdstan_path_to_output = cmdstan_output_basename)
  
  expect_true(is.data.frame(df_posterior))
  expect_true(is.data.frame(df_prior))
  
})



test_that("cmdstanr works", {
  
  stan_path <- file.path(system.file("stan", package = "dvsb"), "dvsb.stan")
  priors_path <- file.path(system.file("input_priors", package = "dvsb"), "priors.csv")
  priors_list <- read_priors(priors_path)
  
  data <- simulate_data(num_plate = 2, num_sam_per_plate = 2)
  
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
  
  iterations <- 10 
  
  df_posterior <- run_stan_interfaces(
    input_to_stan = data_wrangled$stan_input_posterior,
    path_to_stan_code = stan_path,
    interface = "cmdstanr",
    iter_warmup = iterations,
    iter_sampling = iterations)
  
  df_prior <- run_stan_interfaces(
    input_to_stan = data_wrangled$stan_input_prior,
    path_to_stan_code = stan_path,
    interface = "cmdstanr",
    iter_warmup = iterations,
    iter_sampling = iterations)
  
  expect_true(is.data.frame(df_posterior))
  expect_true(is.data.frame(df_prior))
  
})


test_that("rstan works", {
  
  stan_path <- file.path(system.file("stan", package = "dvsb"), "dvsb.stan")
  priors_path <- file.path(system.file("input_priors", package = "dvsb"), "priors.csv")
  priors_list <- read_priors(priors_path)
  
  data <- simulate_data(num_plate = 2, num_sam_per_plate = 2)
  
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
  
  iterations <- 10 
  
  df_posterior <- run_stan_interfaces(
    input_to_stan = data_wrangled$stan_input_posterior,
    path_to_stan_code = stan_path,
    interface = "rstan",
    iter_warmup = iterations,
    iter_sampling = iterations,
    cores = 1)
  
  df_prior <- run_stan_interfaces(
    input_to_stan = data_wrangled$stan_input_prior,
    path_to_stan_code = stan_path,
    interface = "rstan",
    iter_warmup = iterations,
    iter_sampling = iterations,
    cores = 1)
  
  expect_true(is.data.frame(df_posterior))
  expect_true(is.data.frame(df_prior))
  
})

