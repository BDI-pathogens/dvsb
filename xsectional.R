# PRELIMINARIES ----

# See abbreviations below.
# This script models the expected y as a 4PL function of x, with batch-level
# random effects on the 4PL parameters, and observation noise. We simulate data
# from this model and then use the associated Stan file to infer the parameters.

# Abbreviations:
# num = number
# y = optical density = OD
# x = log(antibody concentration)
# mu = mean
# sd = standard deviation
# df = dataframe or degrees of freedom, context-dependent
# con = control (has known x)
# sam = sample (has unknown x)
# pos = (sero-)positive 
# neg = (sero-)negative 
# rep = replicate
# param = parameter
# p = prob = probability

#rm(list = ls())

library(rethinking) # Non-CRAN dependencies: install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos"))); install.packages("remotes"); remotes::install_github("rmcelreath/rethinking") 
library(tidyverse)
library(mvtnorm)
library(rstan)
library(ggforce)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
theme_set(theme_classic())

# INPUT ----

set.seed(123)

# Stan things
file_input_stan <- "~/PathogenDynamics Dropbox/Vaccine Work/Lassa/code_serology_model/lassa_serology_model_crosssectional_v5.stan"
sample_prior_manually <- TRUE
num_mc_chains <- 4
num_mc_iterations_posterior <- 700
num_mc_iterations_prior <- 5000

# Switch between normal and student t distributions
use_student_for_obs <- TRUE
student_df_x   <- 4
student_df_obs <- 5

# Unmodelled aspects of the data-generating process (things we condition on)
num_bat <- 15
num_rep_per_con <- 2
num_rep_per_sam <- 2
num_sam_per_bat <- 20
xs <- c(-0.7055697, 0.3930426, 1.4916549, 2.5902672, 3.6888795) # concentrations of cons

y_obs_sd_min <- 0.001
y_obs_sd_jump <- 0.016
x_exp_sam_neg_alpha <- 1
x_exp_sam_pos_alpha_jump <- 1
x_exp_sam_pos_beta <- 2/3 # beta a.k.a. gamma a.k.a rate
x_exp_sam_neg_beta_jump <- 30
p_sam_pos <- 0.5

# The four parameters of the logistic regression (f_1, f_2, f_3, f_4)
# which control the OD, y, through
# f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (x - f_4))) 
f <- c(0.8,
       0,
       3.5,
       2.5)

# The covariance matrix for the batch-level random effects on f, parameterised
# by the square root of the diagonal entries and the dimensionless correlation
# matrix.
sigma_bat_f <- c(0.01,
                 0.005,
                 0.01,
                 0.03)
Rho_bat <- matrix(c(1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1),
                  4, 4, byrow = TRUE)
stopifnot(isSymmetric(Rho_bat))
stopifnot(all(diag(Rho_bat) == 1))

# Upper and lower bounds for priors
df_priors_scalars <- tribble(
  ~param, ~lower, ~upper,
  "x_exp_sam_neg_alpha", 0, x_exp_sam_neg_alpha + 1,
  "x_exp_sam_pos_alpha_jump", 0, x_exp_sam_pos_alpha_jump + 1,
  "x_exp_sam_pos_beta", 0, x_exp_sam_pos_beta * 2,
  "x_exp_sam_neg_beta_jump", 0, x_exp_sam_neg_beta_jump * 2,
  "p_sam_pos", 0, 1,
  "y_obs_sd_min",  0, 2 * y_obs_sd_min,
  "y_obs_sd_jump", 0, 2 * y_obs_sd_jump
)
df_priors_vectors <- tribble(
  ~param, ~lower, ~upper,
  "sigma_bat_f", sigma_bat_f * 0.5, sigma_bat_f * 2,
  "f", f-1.5, f+1.5
)

# The eta parameter of the LKJ prior for Rho_bat
Rho_bat_prior_eta <- 2

# SIMULATE ----

draw_student_or_norm <- function(n, student, student_df) {
  stopifnot(is.logical(student))
  if (student) return(rt(n, df = student_df))
  return(rnorm(n))
}

# Derived params
x_exp_sam_neg_beta <- x_exp_sam_pos_beta + x_exp_sam_neg_beta_jump
x_exp_sam_pos_alpha <- x_exp_sam_neg_alpha + x_exp_sam_pos_alpha_jump
y_obs_sd_max <- y_obs_sd_min + y_obs_sd_jump
x_exp_sam_pos_mu <- x_exp_sam_pos_alpha / x_exp_sam_pos_beta
x_exp_sam_neg_mu <- x_exp_sam_neg_alpha / x_exp_sam_neg_beta


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

PL4 <- function(x, f_1, f_2, f_3, f_4) {
  f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (x - f_4)))
}

# Expand to one row per con (one for each x). Calculate y expected.
df_con <- df_bat %>%
  expand_grid(x = xs, con = 1:num_rep_per_con) %>%
  mutate(which_con = row_number(),
         y_mean = PL4(x, f_1, f_2, f_3, f_4))

# Plot y expected by batch
ggplot(df_con %>% 
         filter(con == 1)) +
  #geom_line(aes(exp(x), y_mean, group = bat, col = as.factor(bat))) +
  geom_line(aes(x, y_mean, group = bat, col = as.factor(bat))) +
  #scale_x_log10(breaks = exp(xs)) +
  labs(x = "x = Ab concentration",
       y = "y = ELISA OD",
       col = "batch") +
  coord_cartesian(expand = F) +
  ylim(0, NA)

# Draw observed y
df_con <- df_con %>%
  mutate(y_obs_sd = PL4(x, f_1, y_obs_sd_min, y_obs_sd_min + y_obs_sd_jump, f_4),
         y = y_mean + y_obs_sd *
           draw_student_or_norm(nrow(.), use_student_for_obs, student_df_obs))

# Plot observed y
ggplot(df_con) +
  geom_point(aes(jitter(exp(x)), y, col = as.factor(bat))) +
  scale_x_log10(breaks = exp(xs)) +
  labs(x = "Concentration",
       y = "OD",
       col = "batch")
ggplot(df_con %>% filter(bat <= 25)) +
  geom_point(aes(exp(x), y)) +
  scale_x_log10(breaks = exp(xs)) +
  labs(x = "Concentration",
       y = "OD",
       col = "batch") +
  facet_wrap(~bat)

# Calculate indepent ML 4PL curves per batch and add to the plot
if (FALSE) {
  library(dr4pl)
  df_ml <- df_con %>%
    split(.$bat) %>%
    map(function(df_) {
      fit <- dr4pl(data = df_, dose = x, response = y)$parameters
      tibble(param = names(fit), value = as.numeric(fit))
    }) %>%
    bind_rows(.id = "bat") %>%
    mutate(bat = as.integer(bat)) %>%
    pivot_wider(names_from = param, values_from = value) %>%
    inner_join(df_bat, by = "bat") %>%
    expand_grid(x = seq(min(xs), max(xs), length.out = 100)) %>%
    mutate(`independent\nmax-likelihood` = theta_1 + (theta_4 - theta_1) / (1 + (x / theta_2)^theta_3),
           truth = PL4(x, f_1, f_2, f_3, f_4)) %>%
    pivot_longer(c("truth", "independent\nmax-likelihood"), names_to = "y", values_to = "value") %>%
    mutate(y = factor(y, levels = c("truth", "independent\nmax-likelihood")))
  ggplot() +
    scale_x_log10(breaks = exp(xs)) +
    labs(x = "x = Ab concentration",
         y = "y = OD",
         col = "batch",
         linetype = "") +
    geom_line(data = df_ml %>% filter(y == "truth"),
              aes(exp(x), value, col = as.factor(bat))) 
  ggsave("~/foo_1.pdf", height = 3.5, width = 4)
  ggplot(df_con) +
    geom_point(aes(exp(x), y, col = as.factor(bat))) +
    scale_x_log10(breaks = exp(xs)) +
    labs(x = "x = Ab concentration",
         y = "y = OD",
         col = "batch",
         linetype = "") +
    facet_wrap(~bat, ncol = 1) +
    geom_line(data = df_ml %>% filter(y == "truth"),
              aes(exp(x), value, col = as.factor(bat))) 
  ggsave("~/foo_2.pdf", height = 8, width = 3.2)
  ggplot(df_con) +
    geom_point(aes(exp(x), y, col = as.factor(bat))) +
    scale_x_log10(breaks = exp(xs)) +
    labs(x = "x = Ab concentration",
         y = "y = OD",
         col = "batch",
         linetype = "") +
    facet_wrap(~bat, ncol = 1) +
    geom_line(data = df_ml, aes(exp(x), value, col = as.factor(bat), linetype = y)) 
  ggsave("~/foo_3.pdf", height = 8, width = 4)
}

# Allocate each unique sample to a bat, draw its x, and calculate its mean y
# using that bat's f parameters...
num_sam_id <- num_bat * num_sam_per_bat
df_sam <- df_bat %>%
  slice(rep(row_number(), num_sam_per_bat)) %>%
  arrange(bat) %>%
  mutate(id_sam = row_number(),
         pos = runif(num_sam_id) < p_sam_pos,
         x = if_else(pos,
                     rgamma(num_sam_id, rate = x_exp_sam_pos_beta, 
                            shape = x_exp_sam_pos_alpha),
                     rgamma(num_sam_id, rate = x_exp_sam_neg_beta, 
                            shape = x_exp_sam_neg_alpha)),
         x = log(x),
         y_mean = PL4(x, f_1, f_2, f_3, f_4))

if (FALSE) {
  ggplot(df_sam) +
    geom_histogram(aes(exp(x)), fill = "grey") +
    scale_x_log10() +
    coord_cartesian(expand = F) +
    labs(x = "x = Ab concentration",
         y = "Number of samples") 
  ggsave("~/foo_4.pdf", height = 4, width = 5)
}

# ... then create the desired number of reps of each sample, and draw their ys
num_sam_tot <- num_sam_id * num_rep_per_sam
df_sam <- df_sam %>%
  slice(rep(row_number(), num_rep_per_sam)) %>%
  arrange(id_sam) %>%
  mutate(which_sam = row_number(),
         y_obs_sd = PL4(x, f_1, y_obs_sd_min, y_obs_sd_min + y_obs_sd_jump, f_4),
         y = y_mean + y_obs_sd *
           draw_student_or_norm(num_sam_tot, use_student_for_obs, student_df_obs)) 

ggplot(df_sam) +
  geom_histogram(aes(y)) +
  scale_x_log10(limits = c(1e-3, 1e1))




# Plot controls and samples by batch
bind_rows(df_con %>% mutate(label = "control") ,
          df_sam %>% mutate(label = if_else(pos, "+ sample", "- sample"))) %>%
  ggplot() +
  geom_point(aes(jitter(exp(x)), y, col = label)) +
  scale_x_log10(breaks = exp(xs)) +
  facet_wrap(~bat) +
  labs(x = "Concentration",
       y = "Observed optical density",
       col = "") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
ggsave("~/lassa_serology_cross-sectional_data.pdf", height = 9, width = 12)  


# RUN STAN ----

stan_input_posterior <- list(
  use_student_for_obs = as.integer(use_student_for_obs),
  student_df_x = student_df_x,
  student_df_obs = student_df_obs,
  num_bat = num_bat,
  num_con_tot = nrow(df_con),
  num_sam_id = num_sam_id,
  num_sam_tot = num_sam_tot,
  which_bat_con = df_con$bat,
  which_bat_sam = df_sam$bat,
  which_id_sam = df_sam$id_sam,
  y_con = df_con$y,
  y_sam = df_sam$y,
  x_con = df_con$x,
  Rho_bat_prior_eta = Rho_bat_prior_eta,
  sample_posterior_not_prior = 1L
)
for (row in 1:nrow(df_priors_scalars)) {
  stan_input_posterior[[paste0(df_priors_scalars$param[[row]], "_lower")]] <-
    df_priors_scalars$lower[[row]]
  stan_input_posterior[[paste0(df_priors_scalars$param[[row]], "_upper")]] <-
    df_priors_scalars$upper[[row]]
}
for (row in 1:nrow(df_priors_vectors)) {
  stan_input_posterior[[paste0(df_priors_vectors$param[[row]], "_lower")]] <-
    df_priors_vectors$lower[[row]]
  stan_input_posterior[[paste0(df_priors_vectors$param[[row]], "_upper")]] <-
    df_priors_vectors$upper[[row]]
}

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
max_treedepth <- 12
samples_posterior <- sampling(model_compiled,
                              data = stan_input_posterior,
                              iter = num_mc_iterations_posterior,
                              chains = num_mc_chains,
                              control = list(max_treedepth = max_treedepth),
                              pars = params_to_ignore,
                              include = FALSE)
if (! sample_prior_manually) {
  samples_prior <- sampling(model_compiled,
                            data = stan_input_prior,
                            iter = num_mc_iterations_prior,
                            chains = num_mc_chains,
                            control = list(max_treedepth = max_treedepth),
                            pars = params_to_ignore,
                            include = FALSE)
}
end_time <- Sys.time()
cat("Finished running Stan at")
print(end_time)
print(end_time - start_time)

# WRANGLE STAN OUTPUT ----

# Ugly code to sample from the prior manually. Sorry programming.
if (sample_prior_manually) {
  
  # Combine df_priors_scalars and df_priors_vectors
  df_priors <- df_priors_scalars
  names(df_priors_vectors$lower[[1]]) <- paste0("sigma_bat_f[", 1:4, "]")
  names(df_priors_vectors$upper[[1]]) <- paste0("sigma_bat_f[", 1:4, "]")
  names(df_priors_vectors$lower[[2]]) <- paste0("f[", 1:4, "]")
  names(df_priors_vectors$upper[[2]]) <- paste0("f[", 1:4, "]")
  
  df_priors <- df_priors %>%
    bind_rows(full_join(df_priors_vectors %>% 
                          select(param, lower) %>%
                          pivot_wider(names_from = param, values_from = lower) %>%
                          unnest_wider(c("sigma_bat_f", "f")) %>%
                          pivot_longer(everything(), names_to = "param", values_to = "lower"),
                        df_priors_vectors %>% 
                          select(param, upper) %>%
                          pivot_wider(names_from = param, values_from = upper) %>%
                          unnest_wider(c("sigma_bat_f", "f")) %>%
                          pivot_longer(everything(), names_to = "param", values_to = "upper"),
                        by = "param"))
  
  # Sample params with simple uniform distributions
  N <- num_mc_iterations_prior
  df_prior_samples <- tibble(sample = 1:N)
  for (row in 1:nrow(df_priors)) {
    df_prior_samples[[df_priors$param[[row]]]] <-
      runif(N, df_priors$lower[[row]], df_priors$upper[[row]])
  }
  
  # Sample other params
  Rho_bat_samples <- rlkjcorr(N, 4, eta = Rho_bat_prior_eta)
  for (i in 1:3) {
    for (j in seq(i+1, length.out = 4-i)) {
      df_prior_samples[[paste0("Rho_bat[", i, ",", j, "]")]] <- 
        purrr::map_dbl(1:N, ~ Rho_bat_samples[.x, i, j])
    }
  }
  df_prior_samples <- df_prior_samples %>%
    mutate(x_exp_sam_pos_alpha = x_exp_sam_neg_alpha + x_exp_sam_pos_alpha_jump,
           x_exp_sam_neg_beta = x_exp_sam_pos_beta + x_exp_sam_neg_beta_jump,
           y_obs_sd_max = y_obs_sd_min + y_obs_sd_jump,
           x_exp_sam_pos_mu = x_exp_sam_pos_alpha / x_exp_sam_pos_beta,
           x_exp_sam_neg_mu = x_exp_sam_neg_alpha / x_exp_sam_neg_beta)
  
} else {
  df_prior_samples <- samples_prior %>%
    as.data.frame() %>%
    mutate(sample = row_number())
}
df_prior_samples$density_type <- "prior"

df_fit_wide <- bind_rows(samples_posterior %>%
                           as.data.frame() %>%
                           mutate(density_type = "posterior",
                                  sample = row_number()),
                         df_prior_samples) %>%
  as_tibble()

# Pivot to long format: one col for all params
df_fit <- df_fit_wide %>%
  pivot_longer(-c("sample", "density_type"), names_to = "param")

# Plot pop-level params: prior, posterior and true value
df_true_pop_params <- tribble(
  ~param, ~value,
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
  "Rho_bat[3,4]", Rho_bat[3,4],
  "x_exp_sam_neg_alpha", x_exp_sam_neg_alpha,
  "x_exp_sam_pos_alpha", x_exp_sam_pos_alpha,
  "x_exp_sam_neg_beta", x_exp_sam_neg_beta,
  "x_exp_sam_pos_beta", x_exp_sam_pos_beta,
  "x_exp_sam_neg_mu", x_exp_sam_neg_mu,
  "x_exp_sam_pos_mu", x_exp_sam_pos_mu,
  "p_sam_pos", p_sam_pos,
  "y_obs_sd_min", y_obs_sd_min,
  "y_obs_sd_max", y_obs_sd_max
)


# PLOT STAN OUTPUT ----

#vec_true_pop_params <- df_true_pop_params$value
#names(vec_true_pop_params) <- df_true_pop_params$param
#mastiff::plot_posterior(df_fit_wide %>% filter(density_type == "posterior"), 
#                        prior_samples = df_fit_wide %>% filter(density_type == "prior"),
#                        true_param_values = vec_true_pop_params,
#                        params_desired = names(vec_true_pop_params),
#                        skip_stanfit_to_dt = TRUE)

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

#ggplot(data = df_fit_wide %>% select(density_type, x_sam_neg_mu, x_sam_pos_mu)) +
#  geom_bin2d(aes(x_sam_neg_mu, x_sam_pos_mu)) +
#  facet_wrap(~density_type)

quantiles <- c(0.025, 0.5, 0.975)
df_sam_x <- df_fit %>%
  filter(density_type == "posterior") %>%
  select(-density_type) %>%
  filter(startsWith(param, "x_exp_sam[")) %>%
  tidyr::extract(param, 
                 into = "id_sam", 
                 regex = "x_exp_sam\\[([0-9]+)\\]") %>%
  mutate(id_sam = as.integer(id_sam)) %>%
  group_by(id_sam) %>%
  reframe(value = quantile(value, probs = quantiles),
          quantile = quantiles) %>%
  pivot_wider(names_from = quantile, names_prefix = "x_exp_q_")
df_sam_x %>%
  left_join(df_sam %>%
              select(id_sam, x) %>%
              distinct(),
            by = "id_sam") %>%
  filter(id_sam %% 15 == 0) %>%
  ggplot() +
  geom_errorbar(aes(exp(x), ymin = x_exp_q_0.025, ymax = x_exp_q_0.975)) +
  geom_point(aes(exp(x), x_exp_q_0.5)) +
  geom_abline() +
  scale_x_log10(breaks = exp(xs)) +
  scale_y_log10(breaks = exp(xs)) +
  labs(x = "True Ab",
       y = "Estimated Ab")
ggsave("~/foo_8.pdf", height = 3.3, width = 3.3)


# TODO: classification plot
df_prob_pos <- df_fit %>%
  filter(density_type == "posterior") %>%
  select(-density_type) %>%
  filter(startsWith(param, "p_sam_is_pos[")) %>%
  tidyr::extract(param, 
                 into = "id_sam", 
                 regex = "p_sam_is_pos\\[([0-9]+)\\]") %>%
  mutate(id_sam = as.integer(id_sam)) %>%
  group_by(id_sam) %>%
  reframe(value = quantile(value, probs = quantiles),
          quantile = quantiles) %>%
  pivot_wider(names_from = quantile, names_prefix = "prob_pos_q_")

inner_join(df_prob_pos,
           df_sam %>% summarise(.by = id_sam, pos = unique(pos)), 
           by = "id_sam") %>%
  mutate(pos = if_else(pos, "pos", "neg")) %>%
  ggplot() +
  #geom_violin(aes(pos, prob_pos_q_0.5)) +
  #geom_sina(aes(pos, prob_pos_q_0.5)) +
  #geom_point(aes(jitter(as.numeric(pos)), prob_pos_q_0.5)) +
  geom_histogram(aes(prob_pos_q_0.5, fill = pos), #y = after_stat(density))
                 position = "identity",
                 alpha = 0.6,
                 bins = 30) +
  labs(y = "Number of samples",
       #x = "Estimated probability of being positive (posterior median)",
       x = "Probability sample is positive",
       fill = "Truth:") +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(limits = c(NA, NA))
ggsave("~/foo_5.pdf", height = 2.7, width = 3.5)

inner_join(df_prob_pos,
           df_sam %>% summarise(.by = id_sam, y = mean(y)), 
           by = "id_sam") %>%
  ggplot() +
  geom_point(aes(y, prob_pos_q_0.5)) +
  labs(x = "Mean observed OD",
       y = "Estimated probability of being positive")  

inner_join(df_sam_x, df_prob_pos, by = "id_sam") %>%
  ggplot() +
  #geom_errorbar(aes(x_q_0.5, ymin = prob_pos_q_0.025, ymax = prob_pos_q_0.975)) +
  #geom_errorbarh(aes(y = prob_pos_q_0.5, xmin = x_q_0.025, xmax = x_q_0.975)) +
  geom_point(aes(x_exp_q_0.5, prob_pos_q_0.5)) +
  labs(x = "Estimated concentration",
       y = "Estimated probability of being positive")

# The posteriors for the 4PL function by batch
df_4pl <- df_fit %>%
  filter(density_type == "posterior") %>%
  filter(startsWith(param, "f_per_bat[")) %>%
  tidyr::extract(param, 
                 into = c("bat", "f_index"), 
                 regex = "f_per_bat\\[([0-9]+),([0-9]+)\\]") %>%
  mutate(bat = as.integer(bat),
         f_index = as.integer(f_index)) %>%
  pivot_wider(names_from = f_index, names_prefix = "f_") %>%
  expand_grid(x = seq(min(xs), max(xs), length.out = 20)) %>%
  mutate(y = PL4(x, f_1, f_2, f_3, f_4))
bats_to_plot <- 1:num_bat
ggplot(df_4pl %>%
         filter(bat %in% bats_to_plot) %>%
         filter(sample %% 10 == 0)) +
  geom_line(aes(exp(x), y, group = sample), alpha = 0.1) +
  geom_line(data = df_bat %>%
              filter(bat %in% bats_to_plot) %>%
              expand_grid(x = seq(min(xs), max(xs), length.out = 20)) %>%
              mutate(y = PL4(x, f_1, f_2, f_3, f_4)),
            aes(exp(x), y), col = "blue", linewidth = 1) +
  geom_point(data = df_con %>%
               filter(bat %in% bats_to_plot),
             aes(exp(x), y)) +
  facet_wrap(~bat) +
  scale_x_log10(breaks = exp(xs)) +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
ggsave("~/foo_6.pdf", height = 3.3, width = 3.3)
ggplot() +
  geom_point(data = df_con %>%
               filter(bat %in% bats_to_plot),
             aes(exp(x), y)) +
  facet_wrap(~bat) +
  scale_x_log10(breaks = exp(xs)) +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
ggsave("~/foo_7.pdf", height = 3.3, width = 3.3)


xs_plot <- seq(from = 0, to = 10,
               length.out = 500)
df_x_distributions_truth <-
  tibble(x_exp = xs_plot,
         `P(x_exp | pos)` = dgamma(x_exp, shape = x_exp_sam_pos_alpha, rate = x_exp_sam_pos_beta),
         `P(x_exp | neg)` = dgamma(x_exp, shape = x_exp_sam_neg_alpha, rate = x_exp_sam_neg_beta),
         `P(x_exp)` = p_sam_pos * `P(x_exp | pos)` + (1 - p_sam_pos) * `P(x_exp | neg)`,
         `P(pos | x_exp)` = p_sam_pos * `P(x_exp | pos)` / `P(x_exp)`) %>%
  pivot_longer(-x_exp)
df_x_distributions <- df_fit_wide %>%
  filter(density_type == "posterior") %>%
  filter(sample %% 10 == 0) %>%
  select("sample", "x_exp_sam_pos_alpha", "x_exp_sam_pos_beta",
         "x_exp_sam_neg_alpha", "x_exp_sam_neg_beta", "p_sam_pos") %>%
  full_join(tibble(x_exp = xs_plot),
            by = character()) %>%
  mutate(`P(x_exp | pos)` = dgamma(x_exp, shape = x_exp_sam_pos_alpha, rate = x_exp_sam_pos_beta),
         `P(x_exp | neg)` = dgamma(x_exp, shape = x_exp_sam_neg_alpha, rate = x_exp_sam_neg_beta),
         `P(x_exp)` = p_sam_pos * `P(x_exp | pos)` + (1 - p_sam_pos) * `P(x_exp | neg)`,
         `P(pos | x_exp)` = p_sam_pos * `P(x_exp | pos)` / `P(x_exp)`) %>%
  select(sample, x_exp, `P(x_exp | pos)`, 
         `P(x_exp | neg)`, `P(x_exp)`, `P(pos | x_exp)`) %>%
  pivot_longer(c("P(x_exp | pos)", "P(x_exp | neg)", "P(x_exp)", "P(pos | x_exp)"))
p <- ggplot() +
  geom_line(data = df_x_distributions,
            aes(x = x_exp, y = value, group = sample), alpha = 0.15) +
  geom_line(data = df_x_distributions_truth,
            aes(x = x_exp, y = value), colour = "blue") +
  facet_wrap(vars(name), scales = "free_y", ncol = 1) +
  labs(x = "x_exp",
       y = "y") +
  scale_x_log10(expand = c(0, 0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, NA))
p
p <- ggplot() +
  geom_line(data = df_x_distributions %>%
              filter(name != "P(pos | x)"),
            aes(x = x, y = value, group = sample), alpha = 0.15) +
  geom_line(data = df_x_distributions_truth %>%
              filter(name != "P(pos | x)"),
            aes(x = x, y = value), colour = "blue", linewidth = 1) +
  facet_wrap(vars(name), scales = "free_y", ncol = 1) +
  labs(x = "x = log(Ab concentration)",
       y = "") +
  scale_x_continuous(expand = c(0, 0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) 
p
ggsave("~/foo_9.pdf", height = 6, width = 3)


# A posterior retrodictive check of model fit
df_posterior_retrodictive <- df_fit %>%
  filter(density_type == "posterior") %>%
  filter(startsWith(param, "y_con_sim[")) %>%
  tidyr::extract(param, 
                 into = "which_con", 
                 regex = "y_con_sim\\[([0-9]+)\\]") %>%
  mutate(which_con = as.integer(which_con)) %>%
  left_join(df_con %>%
              select(bat, x, which_con),
            by = "which_con")
bats_to_plot <- 1:num_bat
ggplot(df_posterior_retrodictive %>% filter(bat %in% bats_to_plot)) +
  geom_violin(aes(x = as.factor(exp(x)), y = value)) +
  geom_point(data = df_con %>% filter(bat %in% bats_to_plot), 
             aes(as.factor(exp(x)), y), 
             col = "blue") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) +
  facet_wrap(~bat, scales = "free") +
  labs(x = "Concentration",
       y = "OD")


