data_was_simulated <- TRUE

# SYNCHRONISE SIMULATED AND REAL DATA PREVIOUS STEPS ----

if (data_was_simulated) {
  df_plate$plate_int <- df_plate$plate
  df_cal$plate_int <- df_cal$plate
} else {
  PL4 <- function(xlog, f_1, f_2, f_3, f_4) {
    f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (xlog - f_4)))
  }
}

# RUN STAN ----

params_to_ignore <- c(
  "f_plate_effects_unscaled",
  "y_cal_mean_per_obs",
  "p_sam_pos_log",
  "p_sam_neg_log",
  "y_cal_mean_per_obs",
  "y_sam_mean_per_obs",
  "y_obs_sd_cal",
  "y_obs_sd_sam"
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
names(df_priors_vectors$lower[[1]]) <- paste0("sigma_f_plate[", 1:4, "]")
names(df_priors_vectors$upper[[1]]) <- paste0("sigma_f_plate[", 1:4, "]")
names(df_priors_vectors$lower[[2]]) <- paste0("f[", 1:4, "]")
names(df_priors_vectors$upper[[2]]) <- paste0("f[", 1:4, "]")

df_priors <- df_priors %>%
  bind_rows(full_join(df_priors_vectors %>% 
            select(param, lower) %>%
            pivot_wider(names_from = param, values_from = lower) %>%
            unnest_wider(c("sigma_f_plate", "f")) %>%
            pivot_longer(everything(), names_to = "param", values_to = "lower"),
          df_priors_vectors %>% 
            select(param, upper) %>%
            pivot_wider(names_from = param, values_from = upper) %>%
            unnest_wider(c("sigma_f_plate", "f")) %>%
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
Rho_plate_samples <- rethinking::rlkjcorr(N, 4, eta = Rho_plate_prior_eta)
for (i in 1:3) {
  for (j in seq(i+1, length.out = 4-i)) {
    df_prior_samples[[paste0("Rho_plate[", i, ",", j, "]")]] <- 
      purrr::map_dbl(1:N, ~ Rho_plate_samples[.x, i, j])
  }
}
df_prior_samples <- df_prior_samples %>%
  mutate(x_sam_pos_alpha = x_sam_neg_alpha + x_sam_pos_alpha_jump,
         x_sam_neg_beta = x_sam_pos_beta + x_sam_neg_beta_jump,
         y_obs_sd_max = y_obs_sd_min + y_obs_sd_jump,
         x_sam_pos_mu = x_sam_pos_alpha / x_sam_pos_beta,
         x_sam_neg_mu = x_sam_neg_alpha / x_sam_neg_beta)

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
if (data_was_simulated) {
df_true_pop_params <- tribble(
  ~param, ~value,
  "f[1]", f[1],
  "f[2]", f[2],
  "f[3]", f[3],
  "f[4]", f[4],
  "sigma_f_plate[1]", sigma_f_plate[1],
  "sigma_f_plate[2]", sigma_f_plate[2],
  "sigma_f_plate[3]", sigma_f_plate[3],
  "sigma_f_plate[4]", sigma_f_plate[4],
  "Rho_plate[1,2]", Rho_plate[1,2],
  "Rho_plate[1,3]", Rho_plate[1,3],
  "Rho_plate[1,4]", Rho_plate[1,4],
  "Rho_plate[2,3]", Rho_plate[2,3],
  "Rho_plate[2,4]", Rho_plate[2,4],
  "Rho_plate[3,4]", Rho_plate[3,4],
  "x_sam_neg_alpha", x_sam_neg_alpha,
  "x_sam_pos_alpha", x_sam_pos_alpha,
  "x_sam_pos_alpha_jump", x_sam_pos_alpha_jump,
  "x_sam_neg_beta_jump", x_sam_neg_beta_jump,
  "x_sam_neg_beta", x_sam_neg_beta,
  "x_sam_pos_beta", x_sam_pos_beta,
  "x_sam_neg_mu", x_sam_neg_mu,
  "x_sam_pos_mu", x_sam_pos_mu,
  "p_sam_pos", p_sam_pos,
  "y_obs_sd_min", y_obs_sd_min,
  "y_obs_sd_jump", y_obs_sd_jump,
  "y_obs_sd_max", y_obs_sd_max
)
}

# PLOT STAN OUTPUT ----

#vec_true_pop_params <- df_true_pop_params$value
#names(vec_true_pop_params) <- df_true_pop_params$param
#mastiff::plot_posterior(df_fit_wide %>% filter(density_type == "posterior"), 
#                        prior_samples = df_fit_wide %>% filter(density_type == "prior"),
#                        true_param_values = vec_true_pop_params,
#                        params_desired = names(vec_true_pop_params),
#                        skip_stanfit_to_dt = TRUE)

p <- ggplot() +
  geom_histogram(data = df_fit %>%
                   filter(param %in% names(df_prior_samples)) %>%
                   mutate(value = if_else(param %in% c("x_sam_pos_mu", "x_sam_neg_mu"),
                                          log(value),
                                          value)),
                 aes(value, fill = density_type, y = after_stat(density)),
                 alpha = 0.6,
                 position = "identity",
                 bins = 50) +
  facet_wrap(~param, scales = "free", nrow = 5) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(expand = FALSE) +
  labs(fill = "",
       x = "param value",
       y = "probability density")
if (data_was_simulated) {
  p <- p + geom_vline(data = df_true_pop_params %>%
                        mutate(value = if_else(param %in% c("x_sam_pos_mu", "x_sam_neg_mu"),
                                               log(value),
                                               value)),
                      aes(xintercept = value))
}
p

# Compare true and estimated sample x 
quantiles <- c(0.025, 0.5, 0.975)
df_sam_x <- df_fit %>%
  filter(density_type == "posterior") %>%
  select(-density_type) %>%
  filter(startsWith(param, "x_sam[")) %>%
  tidyr::extract(param, 
                 into = "id_sam", 
                 regex = "x_sam\\[([0-9]+)\\]") %>%
  mutate(id_sam = as.integer(id_sam)) %>%
  group_by(id_sam) %>%
  reframe(value = quantile(value, probs = quantiles),
          quantile = quantiles) %>%
  pivot_wider(names_from = quantile, names_prefix = "x_q_")
if (data_was_simulated) {
df_sam_x %>%
  left_join(df_sam %>%
              select(id_sam, x) %>%
              distinct(),
            by = "id_sam") %>%
  filter(id_sam %% 15 == 0) %>%
  ggplot() +
  geom_errorbar(aes(x, ymin = x_q_0.025, ymax = x_q_0.975)) +
  geom_point(aes(x, x_q_0.5)) +
  geom_abline() +
  scale_x_log10(breaks = xs) +
  scale_y_log10(breaks = xs) +
  labs(x = "True Ab",
       y = "Estimated Ab")
ggsave("~/foo_8.pdf", height = 3.3, width = 3.3)
}

# Classification plot
quantiles <- c(0.025, 0.5, 0.975)
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
if (data_was_simulated) {
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
}

# Compare prob positivity vs OD
inner_join(df_prob_pos,
           df_sam %>% summarise(.by = id_sam, y = mean(y)), 
           by = "id_sam") %>%
  ggplot() +
  geom_point(aes(y, prob_pos_q_0.5)) +
  labs(x = "Mean observed OD",
       y = "Estimated probability of being positive")  

# Compare prob positivity vs x
inner_join(df_sam_x, df_prob_pos, by = "id_sam") %>%
  ggplot() +
  geom_point(aes(x_q_0.5, prob_pos_q_0.5)) +
  labs(x = "Estimated concentration",
       y = "Estimated probability of being positive")

# Plot cal data by plate  
plate_ints_to_plot <- 1:4
ggplot() +
  geom_point(data = df_cal %>%
               filter(plate_int %in% plate_ints_to_plot),
             aes(x, y)) +
  facet_wrap(~plate) +
  scale_x_log10() +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
ggsave("~/foo_7.pdf", height = 3.3, width = 3.3)

# The posteriors for the 4PL function by plate
xlog_range <- seq(log(0.4), log(40), length.out = 20)
df_4pl <- df_fit %>%
  filter(density_type == "posterior") %>%
  filter(startsWith(param, "f_per_plate[")) %>%
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
plate_ints_to_plot <- 1:5
p <- ggplot(df_4pl %>%
         filter(plate_int %in% plate_ints_to_plot) %>%
         filter(sample %% 10 == 0)) +
  geom_line(aes(x, y, group = sample), alpha = 0.1) +
  geom_point(data = df_cal %>%
               filter(plate_int %in% plate_ints_to_plot),
             aes(x, y)) +
  facet_wrap(~plate) +
  scale_x_log10() +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
if (data_was_simulated) {
  p <- p + 
    geom_line(data = df_plate %>%
                filter(plate_int %in% plate_ints_to_plot) %>%
                expand_grid(xlog = xlog_range) %>%
                mutate(x = exp(xlog)) %>%
                mutate(y = PL4(xlog, f_1, f_2, f_3, f_4)),
              aes(x, y), col = "blue", linewidth = 1)
}
p
ggsave("~/foo_6.pdf", height = 3.3, width = 3.3)


# Plot P(x | pos), P(x | neg), P(x), P(pos | x)
xs_plot <- seq(from = 0, to = 10,
               length.out = 500)
df_x_distributions <- df_fit_wide %>%
  filter(density_type == "posterior") %>%
  filter(sample %% 10 == 0) %>%
  select("sample", "x_sam_pos_alpha", "x_sam_pos_beta",
         "x_sam_neg_alpha", "x_sam_neg_beta", "p_sam_pos") %>%
  full_join(tibble(x = xs_plot),
            by = character()) %>%
  mutate(`P(x | pos)` = dgamma(x, shape = x_sam_pos_alpha, rate = x_sam_pos_beta),
         `P(x | neg)` = dgamma(x, shape = x_sam_neg_alpha, rate = x_sam_neg_beta),
         `P(x)` = p_sam_pos * `P(x | pos)` + (1 - p_sam_pos) * `P(x | neg)`,
         `P(pos | x)` = p_sam_pos * `P(x | pos)` / `P(x)`) %>%
  select(sample, x, `P(x | pos)`, 
         `P(x | neg)`, `P(x)`, `P(pos | x)`) %>%
  pivot_longer(c("P(x | pos)", "P(x | neg)", "P(x)", "P(pos | x)"))
p <- ggplot() +
  geom_line(data = df_x_distributions,
            aes(x = x, y = value, group = sample), alpha = 0.15) +
  facet_wrap(vars(name), scales = "free_y", ncol = 1) +
  labs(x = "x",
       y = "y") +
  #scale_x_continuous(expand = c(0, 0), limits = c(NA, 2.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, NA))
if (data_was_simulated) {
  df_x_distributions_truth <-
    tibble(x = xs_plot,
           xlog = log(x),
           `P(x | pos)` = dgamma(x, shape = x_sam_pos_alpha, rate = x_sam_pos_beta),
           `P(x | neg)` = dgamma(x, shape = x_sam_neg_alpha, rate = x_sam_neg_beta),
           `P(x)` = p_sam_pos * `P(x | pos)` + (1 - p_sam_pos) * `P(x | neg)`,
           `P(pos | x)` = p_sam_pos * `P(x | pos)` / `P(x)`) %>%
    pivot_longer(-c("x", "xlog"))
  p <- p +
    geom_line(data = df_x_distributions_truth,
              aes(x = x, y = value), colour = "blue") 
}
p

# Plot P(x | pos), P(x | neg), P(x), P(pos | x) again but now with logx
xlogs_plot <- log(10) * -30:20 / 10
df_xlog_distributions <- df_fit_wide %>%
  filter(density_type == "posterior") %>%
  filter(sample %% 10 == 0) %>%
  select("sample", "x_sam_pos_alpha", "x_sam_pos_beta",
         "x_sam_neg_alpha", "x_sam_neg_beta", "p_sam_pos") %>%
  full_join(tibble(xlog = xlogs_plot,
                   x = exp(xlog)),
            by = character()) %>%
  mutate(`P(xlog | pos)` = x * dgamma(x, shape = x_sam_pos_alpha, rate = x_sam_pos_beta),
         `P(xlog | neg)` = x * dgamma(x, shape = x_sam_neg_alpha, rate = x_sam_neg_beta),
         `P(xlog)` = p_sam_pos * `P(xlog | pos)` + (1 - p_sam_pos) * `P(xlog | neg)`,
         `P(pos | xlog)` = p_sam_pos * `P(xlog | pos)` / `P(xlog)`) %>%
  select(sample, xlog, `P(xlog | pos)`, 
         `P(xlog | neg)`, `P(xlog)`, `P(pos | xlog)`) %>%
  pivot_longer(c("P(xlog | pos)", "P(xlog | neg)", "P(xlog)", "P(pos | xlog)"))
p <- ggplot() +
  geom_line(data = df_xlog_distributions,
            aes(x = xlog, y = value, group = sample), alpha = 0.15) +
  facet_wrap(vars(name), scales = "free_y", ncol = 1) +
  labs(x = "x",
       y = "y") +
  #scale_x_log10(expand = c(0, 0), limits = c(NA, NA)) +
  scale_x_continuous(expand = c(0, 0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, NA))
if (data_was_simulated) {
  df_x_distributions_truth <-
    tibble(xlog = xlogs_plot,
           x = exp(xlog),
           `P(xlog | pos)` = x * dgamma(x, shape = x_sam_pos_alpha, rate = x_sam_pos_beta),
           `P(xlog | neg)` = x * dgamma(x, shape = x_sam_neg_alpha, rate = x_sam_neg_beta),
           `P(xlog)` = p_sam_pos * `P(xlog | pos)` + (1 - p_sam_pos) * `P(xlog | neg)`,
           `P(pos | xlog)` = p_sam_pos * `P(xlog | pos)` / `P(xlog)`) %>%
    pivot_longer(-c("x", "xlog"))
  p <- p +
    geom_line(data = df_x_distributions_truth,
              aes(x = xlog, y = value), colour = "blue") 
}
p

# Posterior retrodictive check
group_size <- 4
df_plot_group <- df_cal %>%
  select(plate, plate_int) %>%
  distinct() %>%
  mutate(plot_group = 1 + (row_number() - 1) %/% group_size) 
df_plot <- df_fit_wide %>% # TODO: this needed for full dataset, instead of df_fit, due to memory exhaustion
  filter(density_type == "posterior") %>%
  filter(sample %% 10 == 0) %>%
  select(sample, starts_with("y_cal_sim[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  rename(y_sim = value) %>%
  tidyr::extract(param, 
                 into = c("cal_rep"), 
                 regex = "y_cal_sim\\[([0-9]+)\\]") %>%
  mutate(cal_rep = as.integer(cal_rep)) %>%
  left_join(df_cal %>% select(plate, x) %>% mutate(cal_rep = row_number()),
            by = "cal_rep") %>%
  left_join(df_plot_group, by = "plate")
pdf("~/enable_y_dependent_noise.pdf",
    width = 18.5, height = 10.5)
for (group in unique(df_plot_group$plot_group)) {
  cat("Now doing page ", group, " of ", max(df_plot_group$plot_group), "\n")
  df_plot %>%
    filter(plot_group == group) %>%
    mutate(x = round(x, digits = 1)) %>%
    {print(ggplot(.) +
             geom_violin(aes(factor(x), y_sim)) +
             geom_point(data = df_cal %>% 
                          left_join(df_plot_group, by = "plate") %>%
                          filter(plot_group == group) %>%
                          mutate(x = round(x, digits = 1)),
                        aes(factor(x), y), col = "blue") +
             facet_wrap(~plate, nrow = 3, scales = "free") +
             scale_x_discrete(drop = TRUE) +
             labs(x = "log concentration",
                  y = "absorbance")
    )}
}
dev.off()

# Plot the posterior distribution of the population level distribution of point
# estimates of x_sam (not the posterior distribution of the parametric
# pop-level distribution of x_sam)
df_fit_wide %>%
  filter(sample %% 3 == 0) %>%
  filter(density_type == "posterior") %>%
  select(sample, starts_with("x_sam[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  mutate(value = log(value)) %>%
  ggplot() +
  geom_density(aes(value, group = sample), alpha = 0.01) +
  coord_cartesian(expand = F) +
  scale_x_continuous(limits = c(-6, 5)) +
  labs(x = "Ab concentration (log_e)",
       y = "population distribution of point estimates (inverting the 4PL)") +
  NULL 
ggsave("~/enable_PopDistributionOfX_LogScale.pdf", height = 6, width = 6)









  

