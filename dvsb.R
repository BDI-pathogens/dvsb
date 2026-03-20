library(data.table)
library(tidyverse)
theme_set(theme_classic())

data_was_simulated <- FALSE
read_posterior_from_file <- FALSE
#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2025-11-13-21h38m27_chain*.csv") # v19 on all data with 1500 iter, with sex as p_pos predictor, slow but OK convergence 
#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2025-11-25-18h44m40_chain*.csv") # v20 with age as a predictor for all 5 params, with cluster and site and job for p_pos
#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2025-12-03-10h28m49_chain*.csv") # first restriction to baseline only in some time(!), v20 with age as a predictor for all 5 params, with site and job for p_pos
#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2026-03-11-17h13m09_chain*.csv") # first run using p_pos_binary_pred_vars
files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2026-03-16-14h55m08_chain*.csv")

# INPUT ABOUT STAN ----

path_here <- "~/repos/dvsb/"
file_input_code_wrangle_real_data <- "~/PathogenDynamics Dropbox/Vaccine Work/Lassa/code_serology_model/Xsectional_PrepareEnable_v7.R"
dir_stan <- "~/.cmdstan/cmdstan-2.37.0/"
num_mc_chains <- 5
num_mc_iterations_posterior <- 500 # per chain, half of them warmup
num_mc_iterations_prior <- 2000
# one of: "rstan", "cmdstanr", "cmdstan". cmdstan uses cmdstanr for prior sampling.
stan_interface <- "rstan" 
file_stan_temp <- "/Users/cwymant/foo.json" # for writing the data for cmdstan
file_out_stan_basename <- "/Users/cwymant/enable/samples_"

# Upper and lower bounds for priors
df_priors_scalars <- tribble(
  ~param, ~lower, ~upper,
  "mu_neg", -5, -1,
  "sd_neg", 0.5, 2,
  "sd_pos", 0, 2.5, 
  "mu_pos", -1.5, 3,
  "p_pos", 0, 1,
  "p_blank", 0, 0.05,
  "y_obs_sd_cal_min", 0, 0.015,
  "y_obs_sd_cal_jump", 0.1, 1,
  "y_obs_sd_sam_min", 0, 0.03,
  "y_obs_sd_sam_jump", 0.1, 1,
  "sigma_p_pos_pred_vars", 0, 3,
  "sigma_mu_pos_pred_vars", 0, 2,
  "sigma_sd_pos_pred_vars", 0, 1,
  "sigma_mu_neg_pred_vars", 0, 2,
  "sigma_sd_neg_pred_vars", 0, 1,
  "p_pos_binary_effects", -4, 4
)
df_priors_vectors <- tribble(
  ~param, ~lower, ~upper,
  "sigma_f_plate", c(0, 0, 0, 0), c(0.25, 0.025, 2.5, 2),
  "f", c(0.8, -0.05, 3, 2.2), c(1.1, 0.05, 5.5, 3.8),
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
source(file_input_simulate_code)
source(file_input_code_prepare)
source(file_input_code_rename)
source(file_input_code_wrangle_true)
source(file_input_code_run_stan)
source(file_input_code_read_cmdstan)

# GET DATA ----

if (data_was_simulated) {
  data <- simulate_data(num_plate = 2, num_sam_per_plate = 3)
  df_sam <- data$df_sam
  df_plate <- data$df_plate
  df_cal <- data$df_cal
  param_true_values_list <- data$params
  f_pred_vars_names <- data$f_pred_vars_names
  x_mix_pred_vars_names <- data$x_mix_pred_vars_names
  p_pos_binary_pred_vars <- data$p_pos_binary_pred_vars
  x_mix_params <- c("p_pos", "mu_neg", "mu_pos", "sd_neg", "sd_pos")
} else if (! read_posterior_from_file) {
  source(file_input_code_wrangle_real_data)
  data <- prepare_real_enable_data()
  df_sam <- data$df_sam
  df_cal <- data$df_cal
}

# FINISH PREPARING FOR STAN ----

if (! read_posterior_from_file) {
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

# Set up Stan


# SAMPLE THE POSTERIOR AND PRIOR WITH STAN IF DESIRED... ----

if (! read_posterior_from_file) {
  
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
    cmdstan_path_to_output = paste0(file_out_stan_basename, time, "_prior"))

  }



#save.image(paste0(file_out_stan_basename, time, ".RData"))

# ...OR READ THE POSTERIOR AND PRIOR FROM FILE IF DESIRED ----

# TODO

data.table::setnames(df_fit_wide_postonly, mastiff::rename_params_cmdstanfile_to_rstan)
data.table::setnames(df_ps, mastiff::rename_params_cmdstanfile_to_rstan)


# WRANGLE STAN OUTPUT ----

x_mix_params <- c("p_pos", "mu_neg", "mu_pos", "sd_neg", "sd_pos")

df_fit_wide_postonly[, density_type := "posterior"]
df_fit_wide_postonly[, sample := 1:nrow(df_fit_wide_postonly)]
df_ps[, density_type := "prior"]
df_ps[, sample := 1:nrow(df_ps)]

# Merge prior and posterior samples
desired_cols <- names(df_ps)
df_fit_wide_postandprior <- rbind(df_fit_wide_postonly[,..desired_cols],
                                  df_ps) 

# Record true values of params
if (data_was_simulated) {
  df_true_pop_params <- wrangle_true_params(
    param_true_values_list = param_true_values_list,
    data_descriptors = data_wrangled$data_descriptors)
}

setnames(df_fit_wide_postonly, function(names) {
  rename_params_from_stan(names,
                          data_descriptors = data_wrangled$data_descriptors)})
setnames(df_fit_wide_postandprior, function(names) {
  rename_params_from_stan(names,
                          data_descriptors = data_wrangled$data_descriptors)})
setnames(df_ps, function(names) {
  rename_params_from_stan(names,
                          data_descriptors = data_wrangled$data_descriptors)})
if (data_was_simulated) {
  df_true_pop_params$param <- rename_params_from_stan(
    df_true_pop_params$param, data_descriptors = data_wrangled$data_descriptors)
}

# PLOT STAN OUTPUT ----

param_names_pop <- colnames(df_fit_wide_postandprior)
param_names_pop <- param_names_pop[param_names_pop != "sample"]

regex_for_params_to_plot <- "" #"house|risk"
params_desired <- param_names_pop[grepl(regex_for_params_to_plot, param_names_pop)]
params_for_log_transform <- c() # params_desired[grepl("sd_", params_desired)]
list_for_log_transform <- list()
for (param in params_for_log_transform) list_for_log_transform[[param]] <- log
if (data_was_simulated) {
  true_params_to_plot <- df_true_pop_params[grepl(regex_for_params_to_plot,
                                                  df_true_pop_params$param),]$value
  names(true_params_to_plot) <- df_true_pop_params[grepl(regex_for_params_to_plot,
                                                         df_true_pop_params$param),]$param
  p <- mastiff::plot_posterior(
    posterior_samples = df_fit_wide_postandprior %>% select(-sample) %>% filter(density_type == "posterior"),
    prior_samples     = df_fit_wide_postandprior %>% select(-sample) %>% filter(density_type == "prior"),
    params_desired = params_desired,
    true_param_values = true_params_to_plot,
    skip_stanfit_to_dt = TRUE)
} else {
  p <- mastiff::plot_posterior(
    posterior_samples = df_fit_wide_postandprior %>% select(-sample) %>% filter(density_type == "posterior"),
    prior_samples     = df_fit_wide_postandprior %>% select(-sample) %>% filter(density_type == "prior"),
    params_desired = params_desired,
    skip_stanfit_to_dt = TRUE)
}
p
 
# Compare true and estimated sample x 
quantiles <- c(0.025, 0.5, 0.975)
df_sam_x <- df_fit_wide_postonly %>%
  select(sample, starts_with("xlog_sam[")) %>%
  pivot_longer(-sample, names_to = "param") %>%
  tidyr::extract(param, 
                 into = "id_sam_int", 
                 regex = "xlog_sam\\[([0-9]+)\\]") %>%
  mutate(id_sam_int = as.integer(id_sam_int)) %>%
  group_by(id_sam_int) %>%
  reframe(value = quantile(value, probs = quantiles),
          quantile = quantiles) %>%
  pivot_wider(names_from = quantile, names_prefix = "x_q_")
if (data_was_simulated) {
  df_sam_x %>%
    left_join(df_sam %>%
                select(id_sam_int, xlog) %>%
                distinct(),
              by = "id_sam_int") %>%
    ggplot() +
    geom_errorbar(aes(xlog, ymin = x_q_0.025, ymax = x_q_0.975)) +
    geom_point(aes(xlog, x_q_0.5)) +
    geom_abline() +
    labs(x = "True Ab",
         y = "Estimated Ab")
}

# Classification plot
quantiles <- c(0.025, 0.5, 0.975)
df_prob_pos <- df_fit_wide_postonly %>%
  select(sample, starts_with("p_sam_is_pos[")) %>%
  pivot_longer(-sample, names_to = "param") %>%
  tidyr::extract(param, 
                 into = "id_sam_int", 
                 regex = "p_sam_is_pos\\[([0-9]+)\\]") %>%
  mutate(id_sam_int = as.integer(id_sam_int)) %>%
  group_by(id_sam_int) %>%
  reframe(value = quantile(value, probs = quantiles),
          quantile = quantiles) %>%
  pivot_wider(names_from = quantile, names_prefix = "prob_pos_q_")
if (data_was_simulated) {
  inner_join(df_prob_pos,
             df_sam %>% summarise(.by = id_sam_int, pos = unique(pos)), 
             by = "id_sam_int") %>%
    mutate(pos = if_else(pos, "pos", "neg")) %>%
    ggplot() +
    geom_histogram(aes(prob_pos_q_0.5, fill = pos), #y = after_stat(density))
                   position = "identity",
                   alpha = 0.6,
                   bins = 30) +
    labs(y = "Number of samples",
         x = "Probability sample is positive",
         fill = "Truth:") +
    coord_cartesian(expand = FALSE) +
    scale_x_continuous(limits = c(NA, NA))
}

# Compare prob positivity vs OD
inner_join(df_prob_pos,
           df_sam %>% summarise(.by = id_sam_int, y = mean(y)), 
           by = "id_sam_int") %>%
  ggplot() +
  geom_point(aes(y, prob_pos_q_0.5)) +
  labs(x = "Mean observed OD",
       y = "Estimated probability of being positive")  

# Compare prob positivity vs x
inner_join(df_sam_x, df_prob_pos, by = "id_sam_int") %>%
  ggplot() +
  geom_point(aes(x_q_0.5, prob_pos_q_0.5)) +
  labs(x = "Estimated concentration",
       y = "Estimated probability of being positive") 

# Plot cal data by plate  
ggplot() +
  geom_point(data = df_cal,
             aes(x, y)) +
  facet_wrap(~plate) +
  scale_x_log10() +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 

# The posteriors for the 4PL function by plate
xlog_range <- seq(log(0.002), log(40), length.out = 50)
df_4pl <- df_fit_wide_postonly %>%
  select(sample, starts_with("f_per_plate[")) %>%
  filter(sample %% 10 == 0) %>%
  pivot_longer(-sample, names_to = "param") %>%
  tidyr::extract(param, 
                 into = c("plate_int", "f_index"), 
                 regex = "f_per_plate\\[([0-9]+),([0-9]+)\\]") %>%
  mutate(plate_int = as.integer(plate_int),
         f_index = as.integer(f_index)) %>%
  left_join(df_plate %>% select(plate, plate_int), by = "plate_int") %>%
  pivot_wider(names_from = f_index, names_prefix = "f_") %>%
  expand_grid(xlog = xlog_range) %>%
  mutate(x = exp(xlog)) %>%
  mutate(y = PL4(xlog, f_1, f_2, f_3, f_4))
p <- ggplot(df_4pl %>%
              filter(sample %% 10 == 0)) +
  geom_line(aes(x, y, group = sample), alpha = 0.1) +
  geom_point(data = df_cal,
             aes(x, y), col = "blue") +
  facet_wrap(~plate, nrow = 3) +
  scale_x_log10() +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
if (data_was_simulated) {
  p <- p + 
    geom_line(data = df_plate %>%
                expand_grid(xlog = xlog_range) %>%
                mutate(x = exp(xlog)) %>%
                mutate(y = PL4(xlog, f_1, f_2, f_3, f_4)),
              aes(x, y), col = "blue", linewidth = 1)
}
p

# Plot P(x | pos), P(x | neg), P(x), P(pos | x) again but now with logx
xlogs_plot <- log(10) * -90:60 / 30
df_xlog_distributions <- df_fit_wide_postonly %>%
  filter(sample %% 10 == 0) %>%
  select("sample", "mu_pos", "sd_pos",
         "mu_neg", "sd_neg", "p_pos") %>%
  full_join(tibble(xlog = xlogs_plot,
                   x = exp(xlog)),
            by = character()) %>%
  mutate(`P(xlog | pos)` = dnorm(xlog, mean = mu_pos, sd = sd_pos),
         `P(xlog | neg)` = dnorm(xlog, mean = mu_neg, sd = sd_neg),
         `P(xlog)` = p_pos * `P(xlog | pos)` + (1 - p_pos) * `P(xlog | neg)`,
         `P(pos | xlog)` = p_pos * `P(xlog | pos)` / `P(xlog)`) %>%
  select(sample, xlog, `P(xlog | pos)`, 
         `P(xlog | neg)`, `P(xlog)`, `P(pos | xlog)`) %>%
  pivot_longer(c("P(xlog | pos)", "P(xlog | neg)", "P(xlog)", "P(pos | xlog)"))
p <- ggplot() +
  geom_line(data = df_xlog_distributions,
            aes(x = xlog, y = value, group = sample), alpha = 0.15) +
  facet_wrap(vars(name), scales = "free_y", ncol = 1) +
  labs(x = "log_e(Ab concentration)",
       y = "") +
  scale_x_continuous(expand = c(0, 0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, NA))
if (data_was_simulated) {
  df_x_distributions_truth <-
    tibble(xlog = xlogs_plot,
           x = exp(xlog),
           `P(xlog | pos)` = dnorm(xlog, mean = param_true_values_list$mu_pos, sd = param_true_values_list$sd_pos),
           `P(xlog | neg)` = dnorm(xlog, mean = param_true_values_list$mu_neg, sd = param_true_values_list$sd_neg),
           `P(xlog)` = param_true_values_list$p_pos * `P(xlog | pos)` + (1 - param_true_values_list$p_pos) * `P(xlog | neg)`,
           `P(pos | xlog)` = param_true_values_list$p_pos * `P(xlog | pos)` / `P(xlog)`) %>%
    pivot_longer(-c("x", "xlog"))
  p <- p +
    geom_line(data = df_x_distributions_truth,
              aes(x = xlog, y = value), colour = "blue") 
}
p

# Plot the posterior distribution of the population level distribution of point
# estimates of x_sam and the posterior distribution of the parametric
# pop-level distribution of x_sam
df_fit_wide_postonly %>%
  filter(sample %% 20 == 0) %>%
  select(sample, starts_with("xlog_sam[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  mutate(value = exp(value)) %>%
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = log(1.8)) +
  geom_line(data = df_xlog_distributions %>% filter(name == "P(xlog)"),
            aes(x = xlog, y = value, group = sample), alpha = 0.15, col = "blue") +
  geom_density(aes(log(value) , group = sample), alpha = 0.01) +
  coord_cartesian(expand = F) +
  labs(x = "log_e(Ab concentration)",
       y = "probability density") +
  NULL 

# Plot the posterior distribution of the population level distribution of 
# stochastically redrawn y_sam_sim_conditional 
df_fit_wide_postonly %>%
  filter(sample %% 30 == 0) %>%
  select(sample, starts_with("y_sam_sim_conditional[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  mutate(value = log10(value)) %>%
  ggplot() +
  geom_histogram(data = df_sam, 
                 aes(x = log10(y), y = after_stat(density)), 
                 bins = 60) +
  geom_density(aes(value, group = sample), alpha = 0.01) +
  coord_cartesian(expand = F) +
  scale_x_continuous(limits = c(-3, 1)) +
  labs(x = "log10(OD value)",
       y = "probability density") +
  NULL 

# Plot the posterior distribution of the population level distribution of 
# stochastically redrawn y_sam_sim_unconditional 
df_fit_wide_postonly %>%
  filter(sample %% 100 == 0) %>%
  select(sample, starts_with("y_sam_sim_unconditional[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  tidyr::extract(param, 
                 into = "which_sam_rep", 
                 regex = "y_sam_sim_unconditional\\[([0-9]+)\\]") %>%
  mutate(which_sam_rep = as.integer(which_sam_rep)) %>%
  left_join(df_sam, by = "which_sam_rep") %>%
  mutate(value = log10(value)) %>%
  ggplot() +
  geom_density(aes(value, group = sample), alpha = 0.01) +
  geom_histogram(data = df_sam, 
                 aes(x = log10(y), y = after_stat(density)),
                 col = "blue", fill = NA,
                 bins = 60) +
  coord_cartesian(expand = F) +
  scale_x_continuous(limits = c(-3, 1)) +
  labs(x = "log10(OD value)",
       y = "probability density") +
  NULL 

