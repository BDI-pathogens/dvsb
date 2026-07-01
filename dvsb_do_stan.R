library(data.table)
library(tidyverse)

data_was_simulated <- FALSE

# INPUT ABOUT STAN ----

path_here <- "~/repos/dvsb/"
file_input_code_wrangle_real_data <- "~/PathogenDynamics Dropbox/Vaccine Work/Lassa/code_serology_model/Xsectional_PrepareEnable_v9.R"
dir_stan <- "~/.cmdstan/cmdstan-2.37.0/"
num_mc_chains <- 4
num_mc_iterations_posterior <- 2000 # per chain, half of them warmup
num_mc_iterations_prior <- 5000
# one of: "rstan", "cmdstanr", "cmdstan". cmdstan uses cmdstanr for prior sampling.
stan_interface <- "cmdstan" 
file_stan_temp <- "/Users/cwymant/foo.json" # for writing the data for cmdstan
file_out_stan_basename <- "/Users/cwymant/enable/samples_full_run_"

# Upper and lower bounds for priors
df_priors_scalars <- tribble(
  ~param, ~lower, ~upper,
  #"mu_neg", -3.5, -1.75,
  #"sd_neg", 0.25, 1.8,
  #"sd_pos", 0, 2.5, 
  #"mu_pos", -1.25, 2.25,
  "mu_neg", -4, -2.1, # BEN
  "sd_neg", 0.2, 1, # BEN
  "sd_pos", 0.8, 1.4, # BEN
  "mu_pos", 0.3, 1.2, # BEN
  "p_pos", 0, 1,
  "p_blank", 0, 0.02,
  "y_obs_sd_cal_min", 0, 0.015,
  "y_obs_sd_cal_jump", 0, 1.5,
  "y_obs_sd_sam_min", 0, 0.06,
  "y_obs_sd_sam_jump", 0, 2,
  "sigma_p_pos_pred_vars", 0, 4,
  "sigma_mu_pos_pred_vars", 0, 2,
  "sigma_sd_pos_pred_vars", 0, 1,
  "sigma_mu_neg_pred_vars", 0, 2,
  "sigma_sd_neg_pred_vars", 0, 1,
  "p_pos_binary_effects", -4, 4
)
df_priors_vectors <- tribble(
  ~param, ~lower, ~upper,
  "sigma_f_plate", c(0, 0, 0, 0), c(0.25, 0.015, 2.5, 1),
  "f", c(0.6, -0.05, 0.5, 1.5), c(1.1, 0.05, 7, 3.8),
  "sigma_f_pred_vars", c(0, 0, 0, 0), c(0.6, 0.1, 4, 3)
)
  
# The eta parameter of the LKJ prior for rho
rho_prior_eta <- 1

# SOURCE CODE IN OTHER FILES ----

file_input_stan <- file.path(path_here, "dvsb.stan")
file_input_simulate_code <- file.path(path_here, "R", "dvsb_simulate.R")
file_input_code_prepare <- file.path(path_here, "R", "dvsb_prepare_data_for_stan.R")
file_input_code_rename <- file.path(path_here, "R", "dvsb_rename_params_from_stan.R")
file_input_code_wrangle_true <- file.path(path_here, "R", "dvsb_wrangle_true_params.R")
file_input_code_run_stan <- file.path(path_here, "R", "dvsb_run_stan_interfaces.R")
file_input_code_read_cmdstan <- file.path(path_here, "R", "dvsb_read_cmdstan_out_files.R")
stopifnot(dir.exists(path_here))
stopifnot(file.exists(file_input_stan))
stopifnot(file.exists(file_input_simulate_code))
stopifnot(file.exists(file_input_code_prepare))
stopifnot(file.exists(file_input_code_wrangle_true))
stopifnot(file.exists(file_input_code_run_stan))
stopifnot(file.exists(file_input_code_read_cmdstan))
source(file_input_simulate_code)
source(file_input_code_prepare)
source(file_input_code_rename)
source(file_input_code_wrangle_true)
source(file_input_code_run_stan)
source(file_input_code_read_cmdstan)


# GET DATA ----

if (data_was_simulated) {
  data <- simulate_data(num_plate = 100, num_sam_per_plate = 20)
  df_sam <- data$df_sam
  df_plate <- data$df_plate
  df_cal <- data$df_cal
  param_true_values_list <- data$params
  f_pred_vars_names <- data$f_pred_vars_names
  x_mix_pred_vars_names <- data$x_mix_pred_vars_names
  p_pos_binary_pred_vars <- data$p_pos_binary_pred_vars
} else {
  source(file_input_code_wrangle_real_data)
  data <- prepare_real_enable_data()
  df_sam <- data$df_sam
  df_cal <- data$df_cal
}

# FINISH PREPARING FOR STAN ----

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
# TODO: something cleaner? Adding in cols from wrangling
df_sam <- data_wrangled$df_sam
df_cal <- data_wrangled$df_cal
if (data_was_simulated) {
  df_plate <- left_join(df_plate, data_wrangled$df_plate, by = "plate") 
} else {
  df_plate <- data_wrangled$df_plate
}


params_to_ignore <- c(
  "accept_stat__",
  "stepsize__",
  "treedepth__",
  "n_leapfrog__",
  "divergent__",
  "energy__",
  "lp__",
  "rho.1.1",
  "rho.2.2",
  "rho.3.3",
  "rho.4.4",
  "rho.2.1",
  "rho.3.1",
  "rho.4.1",
  "rho.3.2",
  "rho.4.2",
  "rho.4.3",
  "f_plate_effects_unscaled",
  "y_cal_mean_per_obs",
  "y_cal_mean_per_obs",
  "y_sam_mean_per_obs",
  "y_obs_sd_cal",
  "y_obs_sd_sam",
  "f_effects_by_pred_var_cat_unscaled",
  "p_pos_effects_by_pred_var_cat_unscaled",
  "mu_pos_effects_by_pred_var_cat_unscaled",
  "sd_pos_effects_by_pred_var_cat_unscaled",
  "mu_neg_effects_by_pred_var_cat_unscaled",
  "sd_neg_effects_by_pred_var_cat_unscaled",
  "exp_f_1_mult_f_4_per_plate",
  "p_pos_log_per_sam_id",
  "p_neg_log_per_sam_id",
  "p_pos_per_sam_id",
  "mu_pos_per_sam_id",
  "mu_neg_per_sam_id",
  "sd_pos_per_sam_id",
  "sd_neg_per_sam_id",
  "f_3_min_f_2_per_plate",
  "rho[1,1]",
  "rho[2,2]",
  "rho[3,3]",
  "rho[4,4]",
  "rho[2,1]",
  "rho[3,1]",
  "rho[4,1]",
  "rho[3,2]",
  "rho[4,2]",
  "rho[4,3]",
  "p_blank_log",
  "p_blank_log1m",
  "y_obs_sd_min_log_shift_unscaled",
  "y_obs_sd_jump_log_shift_unscaled",
  "y_obs_sd_sam_min_per_plate",
  "y_obs_sd_cal_min_per_plate",
  "y_obs_sd_sam_jump_per_plate",
  "y_obs_sd_cal_jump_per_plate",
  ".chain",
  ".iteration",
  ".draw",
  # actually interesting, but memory is limited:
  #"xlog_sam",
  #"y_sam_sim_conditional", 
  #"p_sam_is_pos", 
  #"y_sam_loglik_per_obs",
  "xlog_sam_sim_unconditional"
)

# SAMPLE THE POSTERIOR AND PRIOR WITH STAN ----

# Manually write the cmdstan json input files and save image now if desired 
# (if skipping the running of stan here, to do it elsewhere)
#sites_string <- "allsites"
sites_string <- paste(unique(df_sam$site_), collapse = "_")
cmdstanr::write_stan_json(data_wrangled$stan_input_posterior, file = paste0("~/enable_input_posterior_", sites_string, ".json"))
cmdstanr::write_stan_json(data_wrangled$stan_input_prior, file = paste0("~/enable_input_prior_", sites_string, ".json"))
save.image(paste0("~/enable_", sites_string, ".RData"))
# Then restart your R session to ensure the image has only what's wanted

  time <- format(Sys.time(), "%Y-%m-%d-%Hh%Mm%S")
  
  df_fit_wide_postonly <- run_stan_interfaces(
    input_to_stan = data_wrangled$stan_input_posterior,
    path_to_stan_code = file_input_stan,
    interface = stan_interface,
    iterations = num_mc_iterations_posterior,
    chains = num_mc_chains,
    cores = parallel::detectCores(),
    params_to_ignore = params_to_ignore,
    cmdstan_path_to_installation = dir_stan,
    cmdstan_path_to_json = file_stan_temp, 
    cmdstan_overwrite_json = TRUE,
    cmdstan_read_output_into_df = FALSE,
    cmdstan_path_to_output = paste0(file_out_stan_basename, time, "_posterior"))
  
  df_ps <- run_stan_interfaces(
    input_to_stan = data_wrangled$stan_input_prior,
    path_to_stan_code = file_input_stan,
    interface = stan_interface,
    iterations = num_mc_iterations_prior,
    chains = num_mc_chains,
    cores = parallel::detectCores(),
    params_to_ignore = params_to_ignore,
    cmdstan_path_to_installation = dir_stan,
    cmdstan_path_to_json = file_stan_temp, 
    cmdstan_overwrite_json = TRUE,
    cmdstan_read_output_into_df = FALSE,
    cmdstan_path_to_output = paste0(file_out_stan_basename, time, "_prior"))

save.image(paste0(file_out_stan_basename, time, ".RData"))
