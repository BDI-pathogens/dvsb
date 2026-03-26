library(data.table)
library(tidyverse)
theme_set(theme_classic())

# INPUT ----

# If read_files_from_stan=FALSE, you must run dvsb_do_stan.R first to define
# the required variables in your current R session's memory.
read_files_from_stan <- FALSE
#files_from_stan_basename <- "/Users/cwymant/enable/samples_full_run_2026-03-23-11h59m48" # baseline data downsampled by 2
#files_from_stan_basename <- "/Users/cwymant/enable/samples_full_run_2026-03-23-16h51m56" 
files_from_stan_basename <- "/Users/cwymant/enable/samples_full_run_2026-03-25-18h04m08" 
downsampling_factor_posterior <- 1L
downsampling_factor_prior <- 1L

path_here <- "~/repos/dvsb/"

# READ FILES FROM STAN IF DESIRED ----

if (read_files_from_stan) {
  file_input_code_read_cmdstan <- file.path(path_here, "R", "dvsb_read_cmdstan_out_files.R")
  stopifnot(file.exists(file_input_code_read_cmdstan))
  source(file_input_code_read_cmdstan)
  load(paste0(files_from_stan_basename, ".RData"))
  files_posterior <- Sys.glob(paste0(files_from_stan_basename, "_posterior_chain*.csv"))
  files_prior     <- Sys.glob(paste0(files_from_stan_basename, "_prior_chain*.csv"))
  df_fit_wide_postonly <- read_cmdstan_out_files(
    file_paths = files_posterior,
    params_to_ignore = params_to_ignore,
    downsampling_factor = downsampling_factor_posterior)
  df_ps <- read_cmdstan_out_files(
    file_paths = files_prior,
    params_to_ignore = params_to_ignore,
    downsampling_factor = downsampling_factor_prior)
} 

# WRANGLE STAN OUTPUT ----

data.table::setnames(df_fit_wide_postonly, mastiff::rename_params_cmdstanfile_to_rstan)
data.table::setnames(df_ps, mastiff::rename_params_cmdstanfile_to_rstan)

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

# The distribution of probabilty of being positive
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
  p <- inner_join(df_prob_pos,
             df_sam %>% summarise(.by = id_sam_int, pos = unique(pos)), 
             by = "id_sam_int") %>%
    mutate(pos = if_else(pos, "pos", "neg")) %>%
    ggplot() 
} else {
  p <- df_prob_pos %>%
    ggplot() 
}
p +
  geom_histogram(aes(prob_pos_q_0.5), #y = after_stat(density)),
                 bins = 100) +
  labs(y = "Number of samples",
       x = "Probability sample is positive") +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(limits = c(NA, NA))
ggsave("~/enable/enable_prob_positive_histogram.pdf", height = 5, width = 5)
p + stat_ecdf(aes(prob_pos_q_0.5)) +
  labs(y = "c.d.f. (Proportion of samples\nwhose probability of being\npositive is less than that)",
       x = "Probability sample is positive") +
  scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
ggsave("~/enable/enable_prob_positive_cdf.pdf", height = 5, width = 5)


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
  filter(sample %% 100 == 0) %>%
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
  scale_y_continuous(limits = c(NA, NA))
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
  filter(sample %% 3 == 0) %>%
  select(sample, starts_with("y_sam_sim_conditional[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  mutate(value = log10(value)) %>%
  ggplot() +
  geom_histogram(data = df_sam, 
                 aes(#x = y,
                     x = log10(y),
                     y = after_stat(density)), 
                 bins = 100) +
  geom_density(aes(value, group = sample), alpha = 0.01) +
  coord_cartesian(expand = F) +
  scale_x_continuous(limits = c(-3, 1)) +
  labs(#x = "OD value",
       x = "log10(OD value)",
       y = "probability density") +
  NULL 
#ggsave("~/enable/enable_retrodictive_check_linear.pdf", height = 4, width = 6)
ggsave("~/enable/enable_retrodictive_check_logarithmic.pdf", height = 4, width = 6)

# Plot the posterior distribution of the population level distribution of 
# stochastically redrawn y_sam_sim_unconditional 
df_fit_wide_postonly %>%
  filter(sample %% 30 == 0) %>%
  #filter(sample == 1000 ) %>%
  select(sample, starts_with("y_sam_sim_unconditional[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  tidyr::extract(param, 
                 into = "which_sam_rep", 
                 regex = "y_sam_sim_unconditional\\[([0-9]+)\\]") %>%
  mutate(which_sam_rep = as.integer(which_sam_rep)) %>%
  left_join(df_sam, by = "which_sam_rep") %>%
  #mutate(value = log10(value)) %>%
  ggplot() +
  #geom_density(aes(value, group = sample), alpha = 0.01) +
  #geom_histogram(aes(value, y = after_stat(density)), fill = NA, col = "black", alpha = 0.2, bins = 100) +
  geom_step(aes(x = value,
                y = after_stat(density),
                group = sample),
            stat="bin",
            fill="white", color="black",
            bins = 100,
            alpha = 0.1) +
  geom_step(data = df_sam,
            aes(x = y,
                y = after_stat(density)),
            stat="bin",
            bins = 100,
            fill="white", color="blue") +
  coord_cartesian(expand = F) +
  scale_x_continuous(limits = c(-0.2, 3)) +
  labs(x = "log10(OD value)",
       y = "probability density") +
  NULL 
ggsave("~/enable/enable_retrodictive_check_linear.pdf", height = 4, width = 6)

# Plot plate-level random effects on the stochastic noise in OD values
if (data_was_simulated) {
  df_fit_wide_postonly %>%
    select(sample, starts_with("y_obs_sd_min_multiplier_per_plate")) %>%
    pivot_longer(-sample, names_to = "param") %>%
    tidyr::extract(param, 
                   into = "plate_int", 
                   regex = "y_obs_sd_min_multiplier_per_plate\\[([0-9]+)\\]") %>%
    mutate(plate_int = as.integer(plate_int)) %>%
    ggplot() +
    geom_violin(aes(as.factor(plate_int), log10(value))) +
    geom_point(data = df_plate,
               aes(plate_int, log10(y_obs_sd_min_multiplier_per_plate))) +
    labs(x = "plate index",
         y = "OD stochastic noise multiplier for lower asymptote")
  df_fit_wide_postonly %>%
    select(sample, starts_with("y_obs_sd_jump_multiplier_per_plate")) %>%
    pivot_longer(-sample, names_to = "param") %>%
    tidyr::extract(param, 
                   into = "plate_int", 
                   regex = "y_obs_sd_jump_multiplier_per_plate\\[([0-9]+)\\]") %>%
    mutate(plate_int = as.integer(plate_int)) %>%
    ggplot() +
    geom_violin(aes(as.factor(plate_int), log10(value))) +
    geom_point(data = df_plate,
               aes(plate_int, log10(y_obs_sd_jump_multiplier_per_plate))) +
    labs(x = "plate index",
         y = "OD stochastic noise multiplier for difference in asymptotes")
}
