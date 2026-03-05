# PRELIMINARIES ----

# This script simulates infection events for each person and each time step,
# randomly based only on a dynamic infection risk at each time, then calculates
# antibody trajectories based only on the amount of time since the last
# infection event (forgetting any before that) plus some non-dynamic
# individual-level variation.

# Abbreviations:
# num = number
# df = dataframe
# ab = antibody
# base = baseline ab
# jump = increase in ab upon infection
# wane = rate of ab waning after infection
# traj = trajectory
# exp = expected
# mu = mean
# sd = standard deviation
# pred = predictor
# infec = infection
# tsi = time since infection
# beta = a vector of regression coefficients
# x = design matrix, i.e. x * beta shifts the predicted mean for something 
# poss = possibility or possibilities, referring to infection histories and
#        derived quantities like antibody trajectories

#rm(list = ls())

library(tidyverse)
library(mvtnorm)
options(mc.cores = parallel::detectCores())
theme_set(theme_classic())

# Use rstan or cmdstanr? cmdstanr seems more verbose about problematic infinities,
# and is maybe faster to compile and run.
# WARNING: I haven't added support for parsing cmdstanr output for now. 
# It's messy, and slow due to not having the option to ignore some parameters in
# output. Best stick with rstan.
use_cmdstanr <- FALSE
if (use_cmdstanr) {
  library(cmdstanr)
  check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
} else {
  library(rstan)
  rstan_options(auto_write = TRUE)
}


# INPUT ----

# Stan things
file_input_stan <- "~/Infectious Disease Dropbox/Vaccine Work/Lassa/code_serology_model/lassa_serology_model_crosssectional_v2.stan"
num_mc_chains <- 4
num_mc_iterations <- 200

# Unmodelled aspects of the data-generating process (things we condition on)
num_bat <- 20
num_con_per_bat <- 2 # This many replicates per x value
num_sam_per_bat <- 4
xs <- log(1000 * 2^(-2:3))

y_obs_sd <- 0.1

# The four parameters of the logistic regression (f_1, f_2, f_3, f_4)
# which control the log MFI through
# f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (x - f_4))) 
f <- c(4,
       0.5,
       8,
       log(1000))

# The covariance matrix for the batch-level random effects on f, parameterised
# by the square root of the diagonal entries and the dimensionless correlation
# matrix.
sigma_bat_f <- c(0.1,
                 0.05,
                 0.3,
                 0.3)
Rho_bat <- matrix(c(1, 0, 0, -0.1,
                    0, 1, 0.2, 0,
                    0, 0.2, 1, 0,
                    -0.1, 0, 0, 1),
                  4, 4, byrow = TRUE)
stopifnot(isSymmetric(Rho_bat))
stopifnot(all(diag(Rho_bat) == 1))

# SIMULATE ----

# Make a df with one row per batch
df_bat <- tibble(bat = 1:num_bat)

# Assign batch effects
Sigma_bat <- diag(sigma_bat_f) %*% Rho_bat %*% diag(sigma_bat_f)
f_bat_effects <- rmvnorm(num_bat, c(0, 0, 0, 0), Sigma_bat)
df_bat <- df_bat %>%
  mutate(f_1 = f[[1]] + f_bat_effects[, 1],
         f_2 = f[[2]] + f_bat_effects[, 2],
         f_3 = f[[3]] + f_bat_effects[, 3],
         f_4 = f[[4]] + f_bat_effects[, 4])

# Expand to one row per con (one for each x). Calculate y expected.
df_con <- df_bat %>%
  expand_grid(x = xs, con = 1:num_con_per_bat) %>%
  mutate(y_con_mean = f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (x - f_4))))

# Plot y expected by batch
ggplot(df_con %>% 
         filter(con == 1)) +
  geom_line(aes(exp(x), y_con_mean, group = bat, col = as.factor(bat))) +
  scale_x_log10(breaks = exp(xs)) +
  scale_y_log10() +
  labs(x = "Concentration (e.g. pg/mL)",
       y = "Expected MFI",
       col = "batch") +
  coord_cartesian(expand = F)

# Draw observed y
df_con <- df_con %>%
  mutate(y_con = y_con_mean + y_obs_sd * rnorm(nrow(.)))

# Plot observed y
ggplot(df_con) +
  geom_point(aes(jitter(exp(x)), exp(y_con), col = as.factor(bat))) +
  scale_x_log10(breaks = exp(xs)) +
  scale_y_log10() +
  labs(x = "Concentration",
       y = "Observed MFI",
       col = "batch")

# RUN STAN ----

stan_input_posterior <- list(
  num_bat = num_bat,
  num_obs_con = nrow(df_con),
  which_bat_con = df_con$bat,
  y_con = df_con$y_con,
  x_con = df_con$x,
  
  f_lower = c(2, 0, 6,  min(xs)),
  f_upper = c(6, 1, 10, max(xs)),
  sigma_bat_f_lower = 0 * sigma_bat_f,
  sigma_bat_f_upper = 2 * sigma_bat_f,
  y_obs_sd_lower = 0 * y_obs_sd,
  y_obs_sd_upper = 2 * y_obs_sd,
  
  
  sample_posterior_not_prior = 1L
)
stan_input_prior <- stan_input_posterior
stan_input_prior$sample_posterior_not_prior <- 0L

params_to_ignore <- c(
  "f_bat_effects_unscaled",
  "y_con_mean_per_obs"
)

# Compile the Stan code
model_compiled <- stan_model(file_input_stan)

# Run the Stan code
start_time <- Sys.time()
cat("Started running Stan at")
print(start_time)
samples_posterior <- sampling(model_compiled,
                              data = stan_input_posterior,
                              iter = num_mc_iterations,
                              chains = num_mc_chains,
                              pars = params_to_ignore,
                              include = FALSE)
samples_prior <- sampling(model_compiled,
                          data = stan_input_prior,
                          iter = num_mc_iterations,
                          chains = num_mc_chains,
                          pars = params_to_ignore,
                          include = FALSE)
end_time <- Sys.time()
cat("Finished running Stan at")
print(end_time)
print(end_time - start_time)

# WRANGLE AND PLOT STAN OUTPUT ----

df_fit_wide <- bind_rows(samples_posterior %>%
                           as.data.frame() %>%
                           mutate(density_type = "posterior",
                                  sample = row_number()),
                         samples_prior %>%
                           as.data.frame() %>%
                           mutate(density_type = "prior",
                                  sample = row_number())) %>%
  as_tibble()

# Pivot to long format: one col for all params
df_fit <- df_fit_wide %>%
  pivot_longer(-c("sample", "density_type"), names_to = "param")

# Plot pop-level params: prior, posterior and true value
df_true_pop_params <- tribble(
  ~param, ~value,
  "y_obs_sd", y_obs_sd,
  "f[1]", f[1],
  "f[2]", f[2],
  "f[3]", f[3],
  "f[4]", f[4],
  "sigma_bat_f[1]", sigma_bat_f[1],
  "sigma_bat_f[2]", sigma_bat_f[2],
  "sigma_bat_f[3]", sigma_bat_f[3],
  "sigma_bat_f[4]", sigma_bat_f[4],
  "Rho_bat[1,2]", Rho_bat[1,2],
  "Rho_bat[1,3]", Rho_bat[1,3],
  "Rho_bat[1,4]", Rho_bat[1,4],
  "Rho_bat[2,3]", Rho_bat[2,3],
  "Rho_bat[2,4]", Rho_bat[2,4],
  "Rho_bat[3,4]", Rho_bat[3,4]
)
ggplot() +
  geom_histogram(data = df_fit %>%
                   filter(param %in% df_true_pop_params$param),
                 aes(value, fill = density_type, y = after_stat(density)),
                 alpha = 0.6,
                 position = "identity",
                 bins = 50) +
  geom_vline(data = df_true_pop_params, aes(xintercept = value)) +
  facet_wrap(~param, scales = "free", nrow = 3) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(expand = FALSE) +
  labs(fill = "",
       x = "param value",
       y = "probability density")

# The posteriors for the 4PL function by batch
df_fit %>%
  filter(density_type == "posterior") %>%
  filter(startsWith(param, "f_per_bat[")) %>%
  tidyr::extract(param, 
                 into = c("bat", "f_index"), 
                 regex = "f_per_bat\\[([0-9]+),([0-9]+)\\]") %>%
  mutate(bat = as.integer(bat),
         f_index = as.integer(f_index)) %>%
  pivot_wider(names_from = f_index, names_prefix = "f_") %>%
  expand_grid(x = seq(min(xs), max(xs), length.out = 20)) %>%
  mutate(y = f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (x - f_4)))) %>%
  ggplot() +
  geom_line(aes(x, y, group = sample)) +
  geom_line(data = df_bat %>%
              expand_grid(x = seq(min(xs), max(xs), length.out = 20)) %>%
              mutate(y = f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (x - f_4)))),
            aes(x, y), col = "blue") +
  facet_wrap(~bat)

# A posterior retrodictive check of model fit
df_posterior_retrodictive <- df_fit %>%
  filter(density_type == "posterior") %>%
  filter(startsWith(param, "y_con_sim[")) %>%
  tidyr::extract(param, 
                 into = "index", 
                 regex = "y_con_sim\\[([0-9]+)\\]") %>%
  mutate(index = as.integer(index)) %>%
  left_join(df_con %>%
              select(bat, x) %>%
              mutate(index = row_number()),
            by = "index")
ggplot(df_posterior_retrodictive) +
  geom_violin(aes(x = as.factor(exp(x)), y = exp(value))) +
  geom_point(data = df_con, 
             aes(as.factor(exp(x)), exp(y_con)), 
             col = "blue") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) +
  scale_y_log10() +
  facet_wrap(~bat) +
  labs(x = "Concentration",
       y = "MFI")


