# PRELIMINARIES ----

# See abbreviations below.
# This script models the expected y as a 4PL function of x, with plate-level
# random effects on the 4PL parameters, and observation noise. We simulate data
# from this model and then use the associated Stan file to infer the parameters.

# Abbreviations:
# num = number
# y = optical density = OD
# x = log(antibody concentration)
# mu = mean
# sd = standard deviation
# df = dataframe or degrees of freedom, context-dependent
# cal = calibrator (has known x)
# sam = sample (has unknown x)
# pos = (sero-)positive 
# neg = (sero-)negative 
# rep = replicate
# param = parameter
# p = prob = probability

#rm(list = ls())

library(tidyverse)
library(mvtnorm)
library(ggforce)
theme_set(theme_classic())

# INPUT ----

set.seed(123)

# Switch between normal and student t distributions
use_student_for_obs <- FALSE
student_df_obs <- 5

# Unmodelled aspects of the data-generating process (things we condition on)
num_plate <- 50
num_rep_per_cal <- 2
num_rep_per_sam <- 2
num_sam_per_plate <- 20
xlogs <- log(c(0.5, 1.5, 4.5, 13, 40)) #c(-0.7055697, 0.3930426, 1.4916549, 2.5902672, 3.6888795) # concentrations of cals

y_obs_sd_min <- 0.006
y_obs_sd_jump <- 0.015
x_sam_neg_mu <- -2.9
x_sam_neg_sd <- 1
x_sam_pos_mu <- 0.7
x_sam_pos_sd <- 1.3
p_sam_pos <- 0.5

# The four parameters of the logistic regression (f_1, f_2, f_3, f_4)
# which calibrator the OD, y, through
# f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (xlog - f_4))) 
f <- c(0.95,
       0.01,
       3.35,
       2.45)

# The covariance matrix for the plate-level random effects on f, parameterised
# by the square root of the diagonal entries and the dimensionless correlation
# matrix.
sigma_f_plate <- c(0.0065,
                   0.0009,
                   0.09,
                   0.035)
rho <- matrix(c(1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1),
              4, 4, byrow = TRUE)
stopifnot(isSymmetric(rho))
stopifnot(all(diag(rho) == 1))

# f_pred_vars should either be an empty list, or a named list in which 
# each element is a character vector of length at least 2 with no duplicates.
# Each character vector consists of categories that systematically differ in
# their f values. 
# We randomly assign each plate exactly one element from each character vector.
# sigma_f_pred_vars should be a list with the same names as f_pred_vars.
# Each element in sigma_f_pred_vars is a 4-vector of standard deviations of the
# elements of f associated with that predictor variable.
f_pred_vars <- list(
  #OP = c("chris", "anton"),
  #LAB = c("ben", "gui", "lib")
)
sigma_f_pred_vars <- list(
  #OP = 2 * sigma_f_plate,
  #LAB = 3 * sigma_f_plate
)

# SIMULATE ----

draw_student_or_norm <- function(n, student, student_df) {
  stopifnot(is.logical(student))
  if (student) return(rt(n, df = student_df))
  return(rnorm(n))
}

# Check f_pred_vars 
f_pred_vars_names <- names(f_pred_vars)
stopifnot(! any(f_pred_vars_names %in% # avoid name clashes with variables
                  c("plate", "f", "f_1", "f_2", "f_3", "f_4", "label")))
stopifnot(identical(sort(f_pred_vars_names),
                    sort(names(sigma_f_pred_vars))))
num_f_pred_vars <- length(f_pred_vars)
predict_f <- num_f_pred_vars > 0L
if (predict_f) {
  for (name_ in f_pred_vars_names) {
    stopifnot(is.character(f_pred_vars[[name_]]))
    stopifnot(length(f_pred_vars[[name_]]) >= 2L)
    stopifnot(!anyDuplicated(f_pred_vars[[name_]]))
  }
}
num_cat_per_f_pred_var <- map_int(f_pred_vars, length)
num_f_pred_var_cats <- sum(num_cat_per_f_pred_var)

# Derived params
y_obs_sd_max <- y_obs_sd_min + y_obs_sd_jump

xs <- exp(xlogs)

# Make a df with one row per plate
df_plate <- tibble(plate = 1:num_plate)
for (f_pred_var in f_pred_vars_names) {
  df_plate[[f_pred_var]] <- sample(f_pred_vars[[f_pred_var]],
                                   size = num_plate,
                                   replace = TRUE)
}

# Draw plate-level variation in f
Sigma_plate <- diag(sigma_f_plate) %*% rho %*% diag(sigma_f_plate)
f_plate_effects <- rmvnorm(num_plate, c(0, 0, 0, 0), Sigma_plate)
df_plate$f_effect_plate <- map(1:num_plate, ~ f_plate_effects[.x, ])

# Draw variation in f due to f_pred_vars
f_effects_by_pred_var <- list()
for (f_pred_var in f_pred_vars_names) {
  sigma_f_ <- sigma_f_pred_vars[[f_pred_var]]
  Sigma_f_ <- diag(sigma_f_) %*% rho %*% diag(sigma_f_)
  num_cats <- length(f_pred_vars[[f_pred_var]])
  f_effects_ <- rmvnorm(num_cats, c(0, 0, 0, 0), Sigma_f_)
  rownames(f_effects_) <- f_pred_vars[[f_pred_var]]
  f_effects_by_pred_var[[f_pred_var]] <- f_effects_
  df_plate[[paste0("f_effect_", f_pred_var)]] <- map(
    df_plate[[f_pred_var]], ~ f_effects_[.x, ])
}

# Assign f by plate
df_plate$f <- map(1:num_plate, ~ f)
for (plate in 1:num_plate) {
  df_plate$f[[plate]] <- f + df_plate$f_effect_plate[[plate]]
}
for (f_pred_var in f_pred_vars_names) {
  for (plate in 1:num_plate) {
    df_plate$f[[plate]] <- df_plate$f[[plate]] +
      df_plate[[paste0("f_effect_", f_pred_var)]][[plate]]
  }
}
df_plate <- df_plate %>%
  unnest_wider(f, names_sep = "_")

# Label plates for plotting
df_plate$label <- paste("plate", df_plate$plate)
for (f_pred_var in f_pred_vars_names) {
  df_plate$label <- paste0(df_plate$label, ", ", f_pred_var, " ",
                           df_plate[[f_pred_var]])
}
df_plate <- df_plate %>%
  mutate(label = fct_reorder(label, plate))

PL4 <- function(xlog, f_1, f_2, f_3, f_4) {
  f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (xlog - f_4)))
}

# Expand to one row per cal (one for each x). Calculate y expected.
df_cal <- df_plate %>%
  expand_grid(xlog = xlogs, cal = 1:num_rep_per_cal) %>%
  mutate(x = exp(xlog),
         which_cal = row_number(),
         y_mean = PL4(xlog, f_1, f_2, f_3, f_4))

# Plot y expected by plate
ggplot(df_cal %>% 
         filter(cal == 1)) +
  geom_line(aes(xlog, y_mean, group = label, col = label)) +
  labs(x = "x = Ab concentration (log)",
       y = "y = ELISA OD",
       col = "") +
  coord_cartesian(expand = F) +
  ylim(0, NA)

# Draw observed y
df_cal <- df_cal %>%
  mutate(y_obs_sd = PL4(xlog, f_1, y_obs_sd_min, y_obs_sd_max, f_4),
         y = y_mean + y_obs_sd *
           draw_student_or_norm(nrow(.), use_student_for_obs, student_df_obs))

# Plot observed y
ggplot(df_cal) +
  geom_point(aes(jitter(x), y, col = label)) +
  scale_x_log10(breaks = xs) +
  labs(x = "Concentration",
       y = "OD",
       col = "")
ggplot(df_cal %>% filter(plate <= 25)) +
  geom_point(aes(x, y)) +
  scale_x_log10(breaks = xs) +
  labs(x = "Concentration",
       y = "OD",
       col = "plate") +
  facet_wrap(~label)

# Calculate indepent ML 4PL curves per plate and add to the plot
if (FALSE) {
  library(dr4pl)
  df_ml <- df_cal %>%
    split(.$plate) %>%
    map(function(df_) {
      fit <- dr4pl(data = df_, dose = x, response = y)$parameters
      tibble(param = names(fit), value = as.numeric(fit))
    }) %>%
    bind_rows(.id = "plate") %>%
    mutate(plate = as.integer(plate)) %>%
    pivot_wider(names_from = param, values_from = value) %>%
    inner_join(df_plate, by = "plate") %>%
    expand_grid(xlog = seq(min(xlogs), max(xlogs), length.out = 100)) %>%
    mutate(x = exp(xlog)) %>%
    mutate(`independent\nmax-likelihood` = theta_1 + (theta_4 - theta_1) / (1 + (x / theta_2)^theta_3),
           truth = PL4(xlog, f_1, f_2, f_3, f_4)) %>%
    pivot_longer(c("truth", "independent\nmax-likelihood"), names_to = "y", values_to = "value") %>%
    mutate(y = factor(y, levels = c("truth", "independent\nmax-likelihood")))
  ggplot() +
    scale_x_log10(breaks = xs) +
    labs(x = "x = Ab concentration",
         y = "y = OD",
         col = "plate",
         linetype = "") +
    geom_line(data = df_ml %>% filter(y == "truth"),
              aes(x, value, col = as.factor(plate))) 
  ggsave("~/foo_1.pdf", height = 3.5, width = 4)
  ggplot(df_cal) +
    geom_point(aes(x, y, col = as.factor(plate))) +
    scale_x_log10(breaks = xs) +
    labs(x = "x = Ab concentration",
         y = "y = OD",
         col = "plate",
         linetype = "") +
    facet_wrap(~plate) +
    geom_line(data = df_ml %>% filter(y == "truth"),
              aes(x, value, col = as.factor(plate))) 
  ggsave("~/foo_2.pdf", height = 8, width = 3.2)
  ggplot(df_cal) +
    geom_point(aes(x, y, col = as.factor(plate))) +
    scale_x_log10(breaks = xs) +
    labs(x = "x = Ab concentration",
         y = "y = OD",
         col = "plate",
         linetype = "") +
    facet_wrap(~plate) +
    geom_line(data = df_ml, aes(x, value, col = as.factor(plate), linetype = y)) 
  ggsave("~/foo_3.pdf", height = 8, width = 4)
}

# Allocate each unique sample to a plate, draw its x, and calculate its mean y
# using that plate's f parameters...
num_sam_id <- num_plate * num_sam_per_plate
df_sam <- df_plate %>%
  slice(rep(row_number(), num_sam_per_plate)) %>%
  arrange(plate) %>%
  mutate(id_sam = row_number(),
         pos = runif(num_sam_id) < p_sam_pos,
         xlog = if_else(pos,
                        rnorm(num_sam_id, mean = x_sam_pos_mu, 
                              sd = x_sam_pos_sd),
                        rnorm(num_sam_id, mean = x_sam_neg_mu, 
                              sd = x_sam_neg_sd)),
         x = exp(xlog),
         y_mean = PL4(xlog, f_1, f_2, f_3, f_4))

if (FALSE) {
  ggplot(df_sam) +
    geom_histogram(aes(x), fill = "grey") +
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
         y_obs_sd = PL4(xlog, f_1, y_obs_sd_min, y_obs_sd_max, f_4),
         y = y_mean + y_obs_sd *
           draw_student_or_norm(num_sam_tot, use_student_for_obs, student_df_obs)) 

ggplot(df_sam) +
  geom_histogram(aes(y)) +
  scale_x_log10(limits = c(1e-2, 3), expand = c(0, 0))


# Plot calibrators and samples by plate
if (FALSE) {
  bind_rows(df_cal %>% mutate(type = "calibrator") ,
            df_sam %>% mutate(type = if_else(pos, "+ sample", "- sample"))) %>%
    ggplot() +
    geom_point(aes(jitter(x), y, col = type)) +
    scale_x_log10(breaks = xs) +
    facet_wrap(~label) +
    labs(x = "Concentration",
         y = "Observed optical density",
         col = "") +
    theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
  ggsave("~/lassa_serology_cross-sectional_data.pdf", height = 9, width = 12)  
}

# PREPARE DATA FOR STAN ----

stan_input_posterior <- list(
  use_student_for_obs = as.integer(use_student_for_obs),
  student_df_obs = student_df_obs,
  num_plate = num_plate,
  num_cal_tot = nrow(df_cal),
  num_sam_id = num_sam_id,
  num_sam_tot = num_sam_tot,
  which_plate_cal = df_cal$plate,
  which_plate_sam = df_sam$plate,
  which_id_sam = df_sam$id_sam,
  y_cal = df_cal$y,
  y_sam = df_sam$y,
  x_cal = df_cal$x,
  num_f_pred_vars = num_f_pred_vars,
  num_cat_per_f_pred_var = num_cat_per_f_pred_var,
  sample_posterior_not_prior = 1L
)




