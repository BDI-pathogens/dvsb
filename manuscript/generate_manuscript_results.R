library(dvsb)
library(data.table)
library(dr4pl)
library(tidyverse)
theme_set(theme_classic())

outdir <- "~/repos/dvsb/manuscript/"

# DEFINE A FUNC FOR SIMULATING AND RECORDING RELEVANT OUTPUT ----

dir_stan <- "~/.cmdstan/cmdstan-2.37.0/"
file_stan_temp <- "/Users/cwymant/foo.json" # for writing the data for cmdstan
file_out_stan_basename <- "/Users/cwymant/foo_dvsb"
file_out_stan_summary  <- "/Users/cwymant/foo_dvsb_summary.csv"

stan_dir <- system.file("stan", package = "dvsb")
priors_dir <- system.file("input_priors", package = "dvsb")
priors_path <- file.path(priors_dir, "priors.csv")
x_mix_params <- c("p_pos", "mu_neg", "mu_pos", "sd_neg", "sd_pos")

priors_list <- read_priors(priors_path)
df_priors_scalars <- priors_list$df_priors_scalars
df_priors_vectors <- priors_list$df_priors_vectors
rho_prior_eta <- priors_list$rho_prior_eta

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
  "xlog_sam_sim_unconditional",
  "y_obs_sd_cal",
  "y_sam_mean_per_obs",
  "sd_neg_per_sam_id",
  "sd_pos_per_sam_id",
  "mu_neg_per_sam_id",
  "mu_pos_per_sam_id",
  "p_pos_per_sam_id",
  "p_neg_per_sam_id",
  "p_neg_log_per_sam_id",
  "p_pos_log_per_sam_id",
  #"y_sam_sim_unconditional",
  "y_sam_sim_conditional",
  "xlog_sam_sim_unconditional",
  "pos_sam_sim_unconditional",
  "p_sam_is_pos",
  "y_sam_loglik_per_obs",
  "y_obs_sd_sam",
  "y_cal_sim",
  "y_cal_mean_per_obs",
  "p_pos_effects_by_pred_var_cat_unscaled",
  "sd_pos_effects_by_pred_var_cat_unscaled",
  "sd_neg_effects_by_pred_var_cat_unscaled",
  "mu_pos_effects_by_pred_var_cat_unscaled",
  "mu_neg_effects_by_pred_var_cat_unscaled",
  "f_plate_effects_unscaled",
  "loglik",
  "logprob",
  "xlog_sam",
  "mu_pos_per_sam_id",
  "f_per_plate",
  "f_3_min",
  "loglik_per_plate_notblanks",
  "loglik_per_plate_blanks",
  "logprob_f_effects_per_plate"
)


simulate_and_capture <-
  function(dvsb_bin = c("dvsb", "dvsb_no_mu_x_predictors"), ...) {
    dvsb_bin <- match.arg(dvsb_bin)
    stan_path <- file.path(stan_dir, paste0(dvsb_bin, ".stan"))
    stopifnot(file.exists(stan_path))
  data <- simulate_data(...)
  df_sam <- data$df_sam
  df_cal <- data$df_cal
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
  df_sam <- data_wrangled$df_sam
  df_cal <- data_wrangled$df_cal
  
  # Run Stan
  start_time <- Sys.time()
  df_fit_wide_postonly <- run_stan_interfaces(
    input_to_stan = data_wrangled$stan_input_posterior,
    path_to_stan_code = stan_path,
    interface = "cmdstan",
    iter_warmup = num_mc_iterations_posterior,
    iter_sampling = num_mc_iterations_posterior,
    chains = num_mc_chains,
    cores = parallel::detectCores(),
    params_to_ignore = params_to_ignore,
    cmdstan_path_to_installation = dir_stan,
    cmdstan_path_to_json = file_stan_temp, 
    cmdstan_overwrite_json = TRUE,
    cmdstan_read_output_into_df = TRUE,
    cmdstan_path_to_output = paste0(file_out_stan_basename, "_posterior"))
  end_time <- Sys.time()
  
  data.table::setnames(df_fit_wide_postonly, mastiff::rename_params_cmdstanfile_to_rstan)
  setnames(df_fit_wide_postonly, function(names) {
    rename_params_from_stan(names,
                            data_descriptors = data_wrangled$data_descriptors)})
  
  df_y_sam_sim <- df_fit_wide_postonly %>%
    select(starts_with("y_sam_sim_unconditional[")) %>%
    mutate(sample = row_number()) %>%
    filter(sample %% 5 == 0) %>%
    pivot_longer(-c("sample"), names_to = "param") %>%
    tidyr::extract(param, 
                   into = "which_sam_rep", 
                   regex = "y_sam_sim_unconditional\\[([0-9]+)\\]") %>%
    mutate(which_sam_rep = as.integer(which_sam_rep))
  if ("letter" %in% names(df_sam)) {
    df_y_sam_sim <- df_y_sam_sim %>%
      left_join(df_sam %>% select(which_sam_rep, letter), by = "which_sam_rep")   
  }
   

  # Add the x_mix_params regression coefficients onto the base values
  for (param in x_mix_params) {
    effect_params <- names(df_fit_wide_postonly)[grepl(paste0(param, "_effect_"), names(df_fit_wide_postonly))]
    for (effect_param in effect_params) {
      group_name <- str_remove(effect_param, paste0(param, "_effect_"))
      if (param == "p_pos") {
      df_fit_wide_postonly[[paste0(param, "_for_", group_name)]] <-
        mastiff::logistic(mastiff::logit(df_fit_wide_postonly[[param]]) + 
                            df_fit_wide_postonly[[effect_param]])
      } else if (param %in% c("sd_pos", "sd_neg")) {
        df_fit_wide_postonly[[paste0(param, "_for_", group_name)]] <-
          exp(log(df_fit_wide_postonly[[param]]) +
              df_fit_wide_postonly[[effect_param]])
      } else {
        df_fit_wide_postonly[[paste0(param, "_for_", group_name)]] <-
          df_fit_wide_postonly[[param]] + df_fit_wide_postonly[[effect_param]]
      }
    }
  }
  quantiles <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  df_quantiles <-
    df_fit_wide_postonly[, lapply(.SD, quantile, probs = quantiles)]
  df_quantiles[, quantile := quantiles]
  
  # Get Stan summary
  command <- paste0(dir_stan, "/bin/stansummary ",
                    paste0(file_out_stan_basename, "_posterior_chain*.csv"),
                    " -c ", file_out_stan_summary,
                    ' -s 4 -p 2.5,25,50,75,97.5 > /dev/null')
  if (file.exists(file_out_stan_summary)) file.remove(file_out_stan_summary)
  system(command)
  df_summary <- fread(cmd = paste("grep -v '^#'", file_out_stan_summary))
  #df_summary[, name := rename_params_from_stan(
  #  name, data_descriptors = data_wrangled$data_descriptors)]
  R_hat <- max(df_summary$R_hat, na.rm = TRUE)
  df_quantiles[, R_hat := R_hat ]
  df_quantiles[, time := difftime(end_time, start_time, units = "secs") %>% as.numeric() ]
  
  if (FALSE) {
    if ("letter" %in% names(df_sam)) {
      df_sam <- df_sam %>%
        select(letter, y)  
    } else {
      df_sam <- df_sam %>%
        select(y)
    }
  }
  df_quantiles <- df_quantiles %>%
    select(-starts_with("y_sam_sim_unconditional"))
  
  list(df_y_sam_sim = df_y_sam_sim,
       df_quantiles = df_quantiles,
       df_sam = df_sam,
       df_cal = df_cal)
}

# EXPLORE THE EFFECT OF INCREASING DATA SIZE ----

num_plates_range <- 2^(2:3)
num_sam_per_plate <- 10

# Keep constant the effects of different subpopulations (indexed by letter) on
# the x mixture distribution between stochastic resimulated datasets:
x_mix_effects <- list(
  p_pos = list(letter = c(a=0.3, b=-0.3, c=0.6, d=-0.6)),
  mu_pos = list(letter = c(a=0.1, b=-0.1, c=0.6, d=-0.6)),
  mu_neg = list(letter = c(a=0.3, b=-0.3, c=0.2, d=-0.2)),
  sd_pos = list(letter = c(a=0.2, b=-0.2, c=0.6, d=-0.6)),
  sd_neg = list(letter = c(a=0.3, b=-0.3, c=0.4, d=-0.4))
)

num_mc_chains <- 4
num_mc_iterations_posterior <- 100 # per chain, half of them warmup

list_df_sam <- list()
list_df_y_sam_sim <- list()
list_df_quantiles <- list()
for (num_plate in num_plates_range) {
  num_plate_chr <- as.character(num_plate)
  result <- simulate_and_capture(dvsb_bin = "dvsb",
                                 num_plate = num_plate,
                                 x_mix_effects = x_mix_effects)
  list_df_y_sam_sim[[num_plate_chr]] <- result$df_y_sam_sim
  list_df_quantiles[[num_plate_chr]] <- result$df_quantiles
  list_df_sam[[num_plate_chr]] <- result$df_sam
  list_df_y_sam_sim[[num_plate_chr]]$num_plate <- num_plate
  list_df_quantiles[[num_plate_chr]]$num_plate <- num_plate
  list_df_sam[[num_plate_chr]]$num_plate <- num_plate
}

df_y_sam_sim <- rbindlist(list_df_y_sam_sim, use.names = TRUE)
df_quantiles <- rbindlist(list_df_quantiles, use.names = TRUE)
df_sam <- rbindlist(list_df_sam, use.names = TRUE)

# Get the true params in wrangled form. Simulate again for convenience for this.
data <- simulate_data(x_mix_effects = x_mix_effects)
param_true_values_list <- data$params
data_wrangled <- prepare_data_for_stan(
  df_sam = data$df_sam, 
  df_cal = data$df_cal, 
  df_priors_scalars = df_priors_scalars,
  df_priors_vectors = df_priors_vectors, 
  rho_prior_eta = rho_prior_eta, 
  x_mix_pred_vars_names = data$x_mix_pred_vars_names, 
  p_pos_binary_pred_vars = data$p_pos_binary_pred_vars,
  f_pred_vars_names = data$f_pred_vars_names
)
param_true_values_vec <- wrangle_true_params(
  param_true_values_list = data$params,
  data_descriptors =  data_wrangled$data_descriptors)
df_true_pop_params <- tibble(param = names(param_true_values_vec),
                             value = param_true_values_vec)

# Add the x_mix_params regression coefficients onto the base values, for the true params
df_true_pop_params <- df_true_pop_params %>%
  pivot_wider(names_from = param)
for (param in x_mix_params) {
  effect_params <- names(df_true_pop_params)[grepl(paste0(param, "_effect_"), names(df_true_pop_params))]
  for (effect_param in effect_params) {
    group_name <- str_remove(effect_param, paste0(param, "_effect_"))
    if (param == "p_pos") {
      df_true_pop_params[[paste0(param, "_for_", group_name)]] <-
        mastiff::logistic(mastiff::logit(df_true_pop_params[[param]]) + 
                            df_true_pop_params[[effect_param]])
    } else if (param %in% c("sd_pos", "sd_neg")) {
      df_true_pop_params[[paste0(param, "_for_", group_name)]] <-
        exp(log(df_true_pop_params[[param]]) +
              df_true_pop_params[[effect_param]])
    } else {
      df_true_pop_params[[paste0(param, "_for_", group_name)]] <-
        df_true_pop_params[[param]] + df_true_pop_params[[effect_param]]
    }
  }
}
df_true_pop_params <- df_true_pop_params %>%
  pivot_longer(everything(), names_to = "param")

df_run_stats <- df_quantiles %>%
  select(num_plate, time, R_hat) %>%
  distinct()

df_quantiles <- df_quantiles %>%
  select(-c("time", "R_hat")) %>%
  pivot_longer(-c("quantile", "num_plate")) %>%
  pivot_wider(names_from = quantile)
  
df_quantiles %>%
  filter(! name %in% c("sigma_mu_neg_pred_vars_letter",
                       "sigma_mu_pos_pred_vars_letter",
                       "sigma_p_pos_pred_vars_letter",
                       "sigma_sd_neg_pred_vars_letter",
                       "sigma_sd_pos_pred_vars_letter")) %>%
  ggplot() +
  geom_hline(data = df_true_pop_params %>%
               rename(name = param) %>%
               filter(! name %in% c("sigma_mu_neg_pred_vars_letter",
                                    "sigma_mu_pos_pred_vars_letter",
                                    "sigma_p_pos_pred_vars_letter",
                                    "sigma_sd_neg_pred_vars_letter",
                                    "sigma_sd_pos_pred_vars_letter")) %>%
               filter(name %in% df_quantiles$name),
             aes(yintercept = value), col = "blue") +
  geom_boxplot(
    stat = "identity",
    aes(x = factor(num_plate),
        lower  = `0.25`,
        upper  = `0.75`,
        middle = `0.5`,
        ymin   = `0.025`,
        ymax   = `0.975`),
    width = 0.8
  ) +
  scale_x_discrete(expand = c(0.5, 0)) +
  facet_wrap(~name, scale = "free_y") +
  labs(x = "Number of plates (each with 20 samples in duplicate)",
       y = "Posterior estimate")
if (FALSE) ggsave(paste0(outdir, "vary_dataset_size_estimates.pdf"), height = 12, width = 16)

df_run_stats %>%
  ggplot() +
  geom_point(aes(factor(num_plate), time)) +
  geom_text(aes(factor(num_plate), time, label = R_hat),
            nudge_y = -0.02) +
  labs(x = "Number of plates\n(each with 20 samples in duplicate)",
       y = "Run time (seconds)") +
  scale_y_log10()
if (FALSE) ggsave(paste0(outdir, "vary_dataset_size_runtime.pdf"), height = 4, width = 4)

ggplot() +
  geom_density(data = df_y_sam_sim %>%
                 filter(sample %% 10 == 0),
               aes(x = log10(value),
                   group = sample),
               color = "black") +
  geom_density(data = df_sam,
               aes(x = log10(y)),
               color="blue") +
  coord_cartesian(expand = F) +
  facet_grid(letter ~ num_plate, scales = "free_y") +
  labs(x = "log10(OD value)",
       y = "probability density")
if (FALSE) ggsave(paste0(outdir, "vary_dataset_size_posterior_retrodictive.pdf"), height = 12, width = 16)


# EXPLORE THE EFFECT OF OVERLAP BETWEEN - AND + DISTRIBUTIONS ----

# TODO: unhardcode 50% 1/3 + 50% 2/3 as the seroprevalences & their mixture.
p_pos_low <- 1/3
p_pos_high <- 2/3
frac_high <- 0.5
mu_pos_range <- -2.75 + 2:5

x_mix_pred_vars <- list(
  p_pos = list(),
  mu_pos = list(),
  mu_neg = list(),
  sd_pos = list(),
  sd_neg = list()
)
x_mix_pred_vars_sds <- list(
  p_pos = numeric(),
  mu_pos = numeric(),
  mu_neg = numeric(),
  sd_pos = numeric(),
  sd_neg = numeric()
)

num_mc_chains <- 4
num_mc_iterations_posterior <- 400 



list_df_sam <- list()
list_df_cal <- list()
list_df_y_sam_sim <- list()
list_df_quantiles <- list()
for (mu_pos in mu_pos_range) {
  mu_pos_chr <- as.character(mu_pos)
  result = simulate_and_capture(dvsb_bin = "dvsb_no_mu_x_predictors",
                                num_plate = 40,
                                num_sam_per_plate = 25,
                                x_mix_pred_vars = x_mix_pred_vars,
                                x_mix_pred_vars_sds = x_mix_pred_vars_sds,
                                p_pos = p_pos_low,
                                p_pos_binary_effects = c(
                                  risk_factor = mastiff::logit(p_pos_high) - mastiff::logit(p_pos_low)),
                                mu_pos = mu_pos)
  list_df_y_sam_sim[[mu_pos_chr]] <- result$df_y_sam_sim
  list_df_quantiles[[mu_pos_chr]] <- result$df_quantiles
  list_df_sam[[mu_pos_chr]] <- result$df_sam
  list_df_cal[[mu_pos_chr]] <- result$df_cal
  list_df_y_sam_sim[[mu_pos_chr]]$mu_pos <- mu_pos
  list_df_quantiles[[mu_pos_chr]]$mu_pos_true <- mu_pos
  list_df_sam[[mu_pos_chr]]$mu_pos <- mu_pos
  list_df_cal[[mu_pos_chr]]$mu_pos <- mu_pos
}

df_y_sam_sim <- rbindlist(list_df_y_sam_sim, use.names = TRUE)
df_quantiles <- rbindlist(list_df_quantiles, use.names = TRUE)
df_sam <- rbindlist(list_df_sam, use.names = TRUE)
df_cal <- rbindlist(list_df_cal, use.names = TRUE)

df_run_stats <- df_quantiles %>%
  select(mu_pos, time, R_hat) %>%
  distinct()

df_quantiles <- df_quantiles %>%
  select(-c("time", "R_hat")) %>%
  pivot_longer(-c("quantile", "mu_pos_true")) %>%
  pivot_wider(names_from = quantile)

df_true_distributions <- expand_grid(mu_pos = mu_pos_range,
                                     x = seq(-7, 7, by = 0.1)) %>%
  mutate(p = 0.5 * dnorm(x, mean = -2.75, sd = 1) + 
           0.5 * dnorm(x, mean = mu_pos, sd = 1)) %>%
  mutate(label = paste(mu_pos + 2.75, "sigma separation of sero +/-")) 
ggplot(df_true_distributions) +
  geom_line(aes(x, p)) +
  facet_wrap(~label, nrow = 1) +
  labs(y = "probability density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, NA)) 
if (FALSE) ggsave(paste0(outdir, "vary_serostatus_separation_true_distribution.pdf"), height = 3, width = 9)

ggplot() +
  geom_density(data = df_y_sam_sim %>%
                 mutate(label = paste(mu_pos + 2.75, "sigma separation of sero +/-")) %>%
                 filter(sample %% 10 == 0),
               aes(x = log10(value),
                   group = sample),
               color = "black") +
  geom_density(data = df_sam %>%
                 mutate(label = paste(mu_pos + 2.75, "sigma separation of sero +/-")),
               aes(x = log10(y)),
               color="blue") +
  coord_cartesian(expand = F) +
  facet_wrap(~label, scales = "free_y", nrow = 1) +
  xlim(NA, 1) +
  labs(x = "log10(OD value)",
       y = "probability density")
if (FALSE) ggsave(paste0(outdir, "vary_serostatus_separation_posterior_retrodictive.pdf"), height = 3, width = 9)


# Get the true params in wrangled form. Simulate again for convenience for this.
# mu_pos is the only one that will be wrong.
data <- simulate_data(x_mix_pred_vars = x_mix_pred_vars,
                      x_mix_pred_vars_sds = x_mix_pred_vars_sds,
                      p_pos = p_pos_low,
                      p_pos_binary_effects = c(
                        risk_factor = mastiff::logit(p_pos_high) - mastiff::logit(p_pos_low)))
param_true_values_list <- data$params
data_wrangled <- prepare_data_for_stan(
  df_sam = data$df_sam, 
  df_cal = data$df_cal, 
  df_priors_scalars = df_priors_scalars,
  df_priors_vectors = df_priors_vectors, 
  rho_prior_eta = rho_prior_eta, 
  x_mix_pred_vars_names = data$x_mix_pred_vars_names, 
  p_pos_binary_pred_vars = data$p_pos_binary_pred_vars,
  f_pred_vars_names = data$f_pred_vars_names
)
param_true_values_vec <- wrangle_true_params(
  param_true_values_list = data$params,
  data_descriptors =  data_wrangled$data_descriptors)
df_true_pop_params <- tibble(param = names(param_true_values_vec),
                             value = param_true_values_vec)

df_quantiles %>%
  filter(! name %in% c("sigma_mu_neg_pred_vars_letter",
                       "sigma_mu_pos_pred_vars_letter",
                       "sigma_p_pos_pred_vars_letter",
                       "sigma_sd_neg_pred_vars_letter",
                       "sigma_sd_pos_pred_vars_letter")) %>%
  ggplot() +
  geom_hline(data = df_true_pop_params %>%
               rename(name = param) %>%
               filter(! name %in% c("sigma_mu_neg_pred_vars_letter",
                                    "sigma_mu_pos_pred_vars_letter",
                                    "sigma_p_pos_pred_vars_letter",
                                    "sigma_sd_neg_pred_vars_letter",
                                    "sigma_sd_pos_pred_vars_letter")) %>%
               filter(name %in% df_quantiles$name),
             aes(yintercept = value), col = "blue") +
  geom_boxplot(
    stat = "identity",
    aes(x = factor(mu_pos_true),
        lower  = `0.25`,
        upper  = `0.75`,
        middle = `0.5`,
        ymin   = `0.025`,
        ymax   = `0.975`),
    width = 0.8
  ) +
  scale_x_discrete(expand = c(0.5, 0)) +
  facet_wrap(~name, scale = "free_y") +
  labs(x = "Number of plates (each with 20 samples in duplicate)",
       y = "Posterior estimate")

# Fit ML y<->x relationships to each plate
df_ml <- df_cal %>%
  split(list(.$plate, .$mu_pos)) %>%
  map(function(df_) {
    fit <- dr4pl(data = df_, dose = x, response = y)$parameters
    tibble(param = names(fit), value = as.numeric(fit))
  }) %>%
  bind_rows(.id = "plate") %>%
  pivot_wider(names_from = param, values_from = value) %>%
  tidyr::extract(plate, regex = "([0-9]+)\\.(.*)", into = c("plate", "mu_pos")) %>%
  mutate(plate = as.integer(plate),
         mu_pos = as.numeric(mu_pos))

# Find ML x for each sample
invert_y_to_x <- function(y, theta_1, theta_2, theta_3, theta_4) {
  case_when(
    y <= theta_4 ~ 0,
    y >= theta_1 ~ Inf,
    TRUE ~ theta_2 * ( (theta_4 - theta_1) / (y - theta_1) - 1)^(1/theta_3)
  )
}
df_sam_ml <- df_sam %>%
  summarise(.by = c("id_sam", "plate", "mu_pos", "x", "risk_factor"), y = mean(y)) %>%
  left_join(df_ml, by = c("plate", "mu_pos")) %>%
  mutate(x_est = invert_y_to_x(y, theta_1, theta_2, theta_3, theta_4))
df_sam_ml %>%
  ggplot() +
  geom_point(aes(log(x), log(x_est))) 

# Classify each sam as pos or neg according to different thresholds (number of sigmas above the mean for seronegatives) and count 
df_prev_ml <- df_sam_ml %>%
  mutate(pos_est_threshold_2sigma = log(x_est) > -2.75 + 2,
         pos_est_threshold_3sigma = log(x_est) > -2.75 + 3,
         pos_est_threshold_4sigma = log(x_est) > -2.75 + 4) %>%
  summarise(.by = c("mu_pos", "risk_factor"),
            pos_2 = sum(pos_est_threshold_2sigma),
            neg_2 = n() - pos_2,
            pos_3 = sum(pos_est_threshold_3sigma),
            neg_3 = n() - pos_3,
            pos_4 = sum(pos_est_threshold_4sigma),
            neg_4 = n() - pos_4)

# Calculate binomial confints for prevalence in each group.
# Clumsy manual vectorisation because prop.test doesn't vectorise.
df_prev_ml_binomial <- df_prev_ml %>%
  mutate(prev_est_2   = pos_2 / (pos_2 + neg_2),
         prev_est_3   = pos_3 / (pos_3 + neg_3),
         prev_est_4   = pos_4 / (pos_4 + neg_4))
df_prev_ml_binomial$prev_lower_2 <- NA
df_prev_ml_binomial$prev_upper_2 <- NA
df_prev_ml_binomial$prev_lower_3 <- NA
df_prev_ml_binomial$prev_upper_3 <- NA
df_prev_ml_binomial$prev_lower_4 <- NA
df_prev_ml_binomial$prev_upper_4 <- NA
for (row in 1:nrow(df_prev_ml_binomial)) {
  confints_2 <- prop.test(df_prev_ml$pos_2[[row]], df_prev_ml$pos_2[[row]] + df_prev_ml$neg_2[[row]])$conf.int
  df_prev_ml_binomial$prev_lower_2[[row]] <- confints_2[[1]]
  df_prev_ml_binomial$prev_upper_2[[row]] <- confints_2[[2]]
  confints_3 <- prop.test(df_prev_ml$pos_3[[row]], df_prev_ml$pos_3[[row]] + df_prev_ml$neg_3[[row]])$conf.int
  df_prev_ml_binomial$prev_lower_3[[row]] <- confints_3[[1]]
  df_prev_ml_binomial$prev_upper_3[[row]] <- confints_3[[2]]
  confints_4 <- prop.test(df_prev_ml$pos_4[[row]], df_prev_ml$pos_4[[row]] + df_prev_ml$neg_4[[row]])$conf.int
  df_prev_ml_binomial$prev_lower_4[[row]] <- confints_4[[1]]
  df_prev_ml_binomial$prev_upper_4[[row]] <- confints_4[[2]]
}
df_prev_ml_binomial <- df_prev_ml_binomial %>%
  select(mu_pos, risk_factor, starts_with("prev")) %>%
  pivot_longer(starts_with("prev")) %>%
  tidyr::extract(name, regex = "[a-z]+_([a-z]+)_([0-9]+)", into = c("est", "threshold")) %>%
  pivot_wider(names_from = est, values_from = value)

df_prev_ml_binomial %>%
  mutate(mu_pos = as.numeric(mu_pos),
         method = paste(threshold, "sigma\nFrequentist\nstepwise")) %>%
  bind_rows(df_quantiles %>%
              filter(name %in% c("p_pos", "p_pos_for_risk_factor")) %>%
              rename(mu_pos = mu_pos_true,
                     lower = `0.025`,
                     upper = `0.975`,
                     est = `0.5`) %>%
              mutate(method = "dvsb",
                     risk_factor = name == "p_pos_for_risk_factor")) %>%
  mutate(subpop = if_else(risk_factor, "high risk group", "low risk group"),
         label = paste(mu_pos + 2.75, "sigma separation of sero +/-")) %>%
  ggplot() +
  geom_errorbar(aes(method, ymin = lower, ymax = upper)) +
  geom_point(aes(method, est)) +
  geom_hline(data = tibble(subpop = c("low risk group", "high risk group"), prev = c(p_pos_low, p_pos_high)),
             aes(yintercept = prev),
             col = "blue") +
  facet_grid(subpop ~ label) +
  labs(y = "seroprevalence")
if (FALSE) ggsave(paste0(outdir, "vary_serostatus_separation_absolute_seroprevalences.pdf"), height = 7, width = 13)

  
# Do a logistic regression for the difference in seroprevalence log odds between
# the two groups 
df_prev_ml_logistic <- df_prev_ml %>%
  pivot_longer(-c("mu_pos", "risk_factor"), values_to = "n") %>%
  tidyr::extract(name, regex = "([a-z]+)_([0-9])", into = c("serostatus", "threshold")) %>%
  mutate(serostatus = serostatus == "pos") %>%
  split(list(.$mu_pos, .$threshold)) %>%
  map(function(df_) {    
    glm(serostatus ~ risk_factor,
        weights = n, 
        family = "binomial",
        data = df_)}) %>%
  map(function(glm_) {
    inner_join(summary(glm_)$coefficients %>% as_tibble(rownames = "param"),
               confint(glm_) %>% as_tibble(rownames = "param"),
               by = "param")}) %>%
  bind_rows(.id = "mu_pos.threshold") %>%
  tidyr::extract(mu_pos.threshold,
                 regex = "(.*)\\.([0-9])$",
                 into = c("mu_pos", "threshold")) %>%
  mutate(threshold = as.numeric(threshold))

# Plot ML vs dvsb estimates of the difference in seroprevalence log odds between
# the two groups 
df_prev_ml_logistic %>%
  filter(param == "risk_factorTRUE") %>%
  mutate(mu_pos = as.numeric(mu_pos),
         method = paste(threshold, "sigma\nFrequentist\nstepwise")) %>%
  bind_rows(df_quantiles %>%
              filter(name == "p_pos_effect_risk_factor") %>%
              rename(mu_pos = mu_pos_true,
                     `2.5 %` = `0.025`,
                     `97.5 %` = `0.975`,
                     Estimate = `0.5`) %>%
              mutate(method = "dvsb")) %>%
  mutate(label = paste(mu_pos + 2.75, "sigma separation of sero +/-")) %>%
  ggplot() +
  geom_errorbar(aes(method, ymin = `2.5 %`, ymax = `97.5 %`)) +
  geom_point(aes(method, Estimate)) +
  geom_hline(yintercept = mastiff::logit(p_pos_high) - mastiff::logit(p_pos_low),
             col = "blue") +
  facet_wrap(~label, nrow = 1) +
  labs(y = "difference in seropositivity log odds between groups")
if (FALSE) ggsave(paste0(outdir, "vary_serostatus_separation_seroprevalence_difference.pdf"), height = 4, width = 12)

