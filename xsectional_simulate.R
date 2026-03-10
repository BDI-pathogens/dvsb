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

set.seed(1234567)

# Unmodelled aspects of the data-generating process (things we condition on)
num_plate <- 1
num_sam_per_plate <- 2
num_rep_per_sam <- 2
num_rep_per_cal <- 2
xlogs <- log(c(0, 0.5, 1.5, 4.5, 13, 40)) #c(-0.7055697, 0.3930426, 1.4916549, 2.5902672, 3.6888795) # concentrations of cals

y_obs_sd_cal_min <- 0.01
y_obs_sd_cal_jump <- 0.4
y_obs_sd_sam_min <- 0.002
y_obs_sd_sam_jump <- 0.4
mu_neg <- -3.4
sd_neg <- 1
mu_pos <- 0.8
sd_pos <- 1.1
p_pos <- 0.5
p_blank <- 0.2

# The four parameters of the logistic regression (f_1, f_2, f_3, f_4)
# which control the OD, y, through
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
  #op_ = letters[1:10]  #c("chris", "anton"),
  #lab_ = c("ben", "gui", "lib")
)
sigma_f_pred_vars <- list(
  #op_ = 2 * sigma_f_plate
  #lab_ = 3 * sigma_f_plate
)

# Initialise empty predictor vars for the x mix distribution (unnecessary
# because we'll overwrite them next, but it shows the structure when empty).
x_mix_pred_vars <- list()
x_mix_pred_vars_sds <- list()
x_mix_params <- c("p_pos", "mu_neg", "mu_pos", "sd_neg", "sd_pos")
for (x_mix_pred_var_ in x_mix_params) {
  x_mix_pred_vars[[x_mix_pred_var_]]  <- list()
  x_mix_pred_vars_sds[[x_mix_pred_var_]] <- numeric()
}

# p_pos_pred_vars should either be an empty list, or a named list in which 
# each element is a character vector of length at least 2 with no duplicates.
# Each character vector consists of categories that systematically differ in
# their p_pos values. 
# We randomly assign each sample exactly one element from each character vector.
# sigma_p_pos_pred_vars should be a list with the same names as f_pred_vars.
# Each element in sigma_p_pos_pred_vars is a standard deviation of the
# variability in p_pos (on a logit scale) associated with that predictor variable.
x_mix_pred_vars$p_pos <- list(letter = letters[1:7],
                              int = as.character(1:7),
                              foo = c("bar", "spam"),
                              age = c("0-9", "10-19", "20+"))
x_mix_pred_vars_sds$p_pos <- c(letter = 1,
                               int = 1,
                               foo = 1,
                               age = 1)
x_mix_pred_vars$mu_neg <- list(letter = letters[1:7])
x_mix_pred_vars_sds$mu_neg <- c(letter = 1)
x_mix_pred_vars$mu_pos <- list(letter = letters[1:7])
x_mix_pred_vars_sds$mu_pos <- c(letter = 1)
x_mix_pred_vars$sd_neg <- list(letter = letters[1:7])
x_mix_pred_vars_sds$sd_neg <- c(letter = 1)
x_mix_pred_vars$sd_pos <- list(letter = letters[1:7])
x_mix_pred_vars_sds$sd_pos <- c(letter = 1)

p_pos_binary_pred_vars <- c("boolA", "boolB", "boolC")
p_pos_binary_pred_vars_sd <- 2

# INPUT CHECKS ----

# Check f_pred_vars 
f_pred_vars_names <- names(f_pred_vars)
stopifnot(! any(f_pred_vars_names %in% # avoid name clashes with variables
                  c("plate", "f", "f_1", "f_2", "f_3", "f_4", "label")))
stopifnot(identical(sort(f_pred_vars_names),
                    sort(names(sigma_f_pred_vars))))
if (is.null(f_pred_vars_names)) f_pred_vars_names <- character() # more intuitive
num_f_pred_vars <- length(f_pred_vars)
predict_f <- num_f_pred_vars > 0L
if (predict_f) {
  if (num_plate < 2) stop("2+ plates are needed if plate covariates are used")
  for (name_ in f_pred_vars_names) {
    stopifnot(is.character(f_pred_vars[[name_]]))
    stopifnot(length(f_pred_vars[[name_]]) >= 2L)
    stopifnot(!anyDuplicated(f_pred_vars[[name_]]))
  }
}

# Check x mix pred vars 
num_sam_id <- num_plate * num_sam_per_plate
stopifnot(identical(names(x_mix_pred_vars),     x_mix_params))
stopifnot(identical(names(x_mix_pred_vars_sds), x_mix_params))
for (param_ in x_mix_params) {
  if (! identical(names(x_mix_pred_vars[[param_]]),
                  names(x_mix_pred_vars_sds[[param_]]))) {
    stop(paste0("Different predictor variables were named in x_mix_pred_vars$",
                param_, " and in x_mix_pred_vars$", param_, ":\n",
                paste(names(x_mix_pred_vars[[param_]]), collapse = " "), "\nand\n",
                paste(names(x_mix_pred_vars_sds[[param_]]), collapse = " "), 
                "\nrespectively. These must be identical.\n"))
  }
}
x_mix_pred_vars_names <- map(x_mix_pred_vars, names)
x_mix_pred_vars_nums <- map_int(x_mix_pred_vars_names, length)
if (any(x_mix_pred_vars_nums) && num_sam_id == 0) {
  stop("You need some samples if the x mix parameters are to be predicted")
}
for (param in x_mix_params) {
  if (x_mix_pred_vars_nums[[param]]) {
    for (pred_var in x_mix_pred_vars_names[[param]]) {
      if (!is.character(x_mix_pred_vars[[param]][[pred_var]])) {
        stop(paste0("x_mix_pred_vars$", param, "$", pred_var,
                    " must be a character vector"))
      }
      if (length(x_mix_pred_vars[[param]][[pred_var]]) < 2L) {
        stop(paste0("x_mix_pred_vars$", param, "$", pred_var,
                    " must contain at least two elements; we found it equal to ",
                    x_mix_pred_vars[[param]][[pred_var]]))
      }
      if (anyDuplicated(x_mix_pred_vars[[param]][[pred_var]])) {
        stop(paste0("x_mix_pred_vars$", param, "$", pred_var,
                    " must not contain duplicates; we ", 
                    "found it equal to ", 
                    paste(x_mix_pred_vars[[param]][[pred_var]], collapse = " ")))
      }
    }
  }
}

for (i in seq(1, 4)) {
  param_1 <- x_mix_params[[i]]
  if (is.null(x_mix_pred_vars_names[[param_1]])) next
  for (j in (i+1):5) {
    param_2 <- x_mix_params[[j]]
    if (is.null(x_mix_pred_vars_names[[param_2]])) next
    pred_vars_shared <- x_mix_pred_vars_names[[param_1]][
      x_mix_pred_vars_names[[param_1]] %in% x_mix_pred_vars_names[[param_2]]]
    for (pred_var in pred_vars_shared) {
      if (! identical(sort(x_mix_pred_vars[[param_1]][[pred_var]]),
                      sort(x_mix_pred_vars[[param_2]][[pred_var]]))) {
        stop(paste0(pred_var, " was specified as a predictor variable for both ",
                    param_1, " and ", param_2, ", but different categories were specified: ",
                    paste(x_mix_pred_vars[[param_1]][[pred_var]], collapse = ", "),
                    " for ", param_1, ", and ",
                    paste(x_mix_pred_vars[[param_2]][[pred_var]], collapse = ", "),
                    " for ", param_2,
                    ". The categories must be the same for a given predictor variable."))
      }
    }
  }
}

# Check p_pos_binary_pred_vars
stopifnot(is.character(p_pos_binary_pred_vars))
stopifnot(is.numeric(p_pos_binary_pred_vars_sd))
if (length(p_pos_binary_pred_vars_sd) == 1) {
  stopifnot(p_pos_binary_pred_vars_sd >= 0)
  if (length(p_pos_binary_pred_vars) == 0) {
    stop(paste("If p_pos_binary_pred_vars_sd has length 1, then",
               "p_pos_binary_pred_vars must have length at least 1"))
  }
  for (x_mix_param in x_mix_params) {
    if (is.null(x_mix_pred_vars_names[[x_mix_param]])) next 
    pred_vars_shared <- x_mix_pred_vars_names[[x_mix_param]][
      x_mix_pred_vars_names[[x_mix_param]] %in% p_pos_binary_pred_vars]
    if (length(pred_vars_shared)) {
      stop(paste0("p_pos_binary_pred_vars must not contain any predictor ",
                  "variables that are also used as (non-binary) predictors for ", 
                  x_mix_param, "; found these variables used for both: ", 
                  paste(pred_vars_shared, collapse = " ")))
    }
  }
} else if (length(p_pos_binary_pred_vars_sd) == 0) {
  if (length(p_pos_binary_pred_vars) != 0) {
    stop(paste("If p_pos_binary_pred_vars_sd has length 0, then",
               "p_pos_binary_pred_vars must have length 0; found length",
               length(p_pos_binary_pred_vars)))
  }
  } else {
  stop(paste("p_pos_binary_pred_vars_sd should have length 0 or 1; found length",
             length(p_pos_binary_pred_vars_sd)))
}

# SIMULATE PLATE VARIABILITY AND CALS ----

# Derived params
y_obs_sd_cal_max <- y_obs_sd_cal_min + y_obs_sd_cal_jump
y_obs_sd_sam_max <- y_obs_sd_sam_min + y_obs_sd_sam_jump

xs <- exp(xlogs)

# Make a df with one row per plate.
# Sample each plate's f predictor variables.
# Ensure that we don't randomly sample the same category for every plate.
# Delete any unsampled categories.
df_plate <- tibble(plate = 1:num_plate)
for (f_pred_var in f_pred_vars_names) {
  sampled_pred_vars <- character()
  while(length(sampled_pred_vars) < 2) {
    sampled_pred_vars <- sample(f_pred_vars[[f_pred_var]],
                                size = num_plate,
                                replace = TRUE)
  } 
  if (length(sampled_pred_vars) < length(f_pred_vars[[f_pred_var]])) {
    f_pred_vars[[f_pred_var]] <- sort(unique(sampled_pred_vars))
  } 
  df_plate[[f_pred_var]] <- sampled_pred_vars
}
num_cat_per_f_pred_var <- map_int(f_pred_vars, length)
num_f_pred_var_cats <- sum(num_cat_per_f_pred_var)

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
  f_effects_col_means <- colMeans(f_effects_)
  for (cat_num in 1:num_cats) {
    f_effects_[cat_num, ] <- f_effects_[cat_num, ] - f_effects_col_means
  }
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
  mutate(y_obs_sd = PL4(xlog, f_1, y_obs_sd_cal_min, y_obs_sd_cal_max, f_4),
         y = rnorm(nrow(.), mean = y_mean, sd = y_obs_sd))

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
  ggsave("~/foo_1.pdf", height = 5.5, width = 6)
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
  ggplot(df_cal %>% filter(plate == 18)) +
    geom_point(aes(x, y)) +
    scale_x_log10(breaks = xs) +
    labs(x = "x = Ab concentration",
         y = "y = OD") +
    geom_line(data = df_ml %>% filter(y == "truth") %>% filter(plate == 18),
              aes(x, value)) 
  ggsave("~/foo_2b.pdf", height = 4, width = 4)
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

# SIMULATE SAMS ----

x_mix_baseline <- c(p_pos  = p_pos,
                    mu_pos = mu_pos,
                    sd_pos = sd_pos,
                    mu_neg = mu_neg,
                    sd_neg = sd_neg)

# Allocate each unique sample to a plate and draw its x mix predictors.
# Ensure that we don't randomly sample the same category for every sample.
# Delete any unsampled categories.
# For those pred vars shared by multiple params, sample once only.
df_sam <- df_plate %>%
  slice(rep(row_number(), num_sam_per_plate)) %>%
  arrange(plate) %>%
  mutate(id_sam = row_number())
x_mix_pred_vars_num_cats <- list()
x_mix_pred_vars_num_cats_tots <- integer()
for (param in x_mix_params) {
  for (pred_var in x_mix_pred_vars_names[[param]]) {
    if (pred_var %in% names(df_sam)) {
      # We've already sampled this pred_var for a previous param.
      # Ensure we remove any unsampled cats the same as previously, then skip.
      x_mix_pred_vars[[param]][[pred_var]] <- sort(unique(df_sam[[pred_var]]))
      next
    } 
    sampled_pred_vars <- character()
    while(n_distinct(sampled_pred_vars) < 2) {
      sampled_pred_vars <- sample(x_mix_pred_vars[[param]][[pred_var]],
                                  size = num_sam_id,
                                  replace = TRUE)
    } 
    x_mix_pred_vars[[param]][[pred_var]] <- sort(unique(sampled_pred_vars))
    df_sam[[pred_var]] <- sampled_pred_vars
  }
  x_mix_pred_vars_num_cats[[param]] <-
    map_int(x_mix_pred_vars[[param]], length)
  x_mix_pred_vars_num_cats_tots[[param]] <- 
    sum(x_mix_pred_vars_num_cats[[param]])
}
for (pred_var in p_pos_binary_pred_vars) {
  sampled_pred_vars <- character()
  while(n_distinct(sampled_pred_vars) < 2) {
    sampled_pred_vars <- sample(c(TRUE, FALSE),
                                size = num_sam_id,
                                replace = TRUE)
  } 
  df_sam[[pred_var]] <- sampled_pred_vars
}
all_names_x_mix_pred_vars <- x_mix_pred_vars %>% 
  map(names) %>%
  unlist() %>%
  unique()
if (is.null(all_names_x_mix_pred_vars)) {
  df_sam <- df_sam %>%
    mutate(x_mix_group = NA_character_)
} else {
  df_sam <- df_sam %>%
    unite("x_mix_group", all_of(all_names_x_mix_pred_vars), sep = "_", remove = FALSE)
}
for (pred_var in p_pos_binary_pred_vars) {
 df_sam$x_mix_group <- paste0(df_sam$x_mix_group, "_", pred_var, df_sam[[pred_var]])
}

# Draw effects on the x mix params from each pred var
x_mix_effects <- list()
x_mix_overall <- list()
for (param in x_mix_params) {
  x_mix_effects[[param]] <- list()
  x_mix_overall[[param]] <- list()
  for (pred_var in x_mix_pred_vars_names[[param]]) {
    sigma_ <- x_mix_pred_vars_sds[[param]][[pred_var]]
    num_cats <- x_mix_pred_vars_num_cats[[param]][[pred_var]]
    effects_ <- rnorm(num_cats, 0, sigma_)
    effects_ <- effects_ - mean(effects_)
    names(effects_) <- x_mix_pred_vars[[param]][[pred_var]]
    x_mix_effects[[param]][[pred_var]] <- effects_
    df_sam[[paste0(param, "_effect_", pred_var)]] <- map_dbl(
      df_sam[[pred_var]], ~ effects_[[.x]])
    if (param == "p_pos") {
      x_mix_overall[[param]][[pred_var]] <- 
        mastiff::logistic(mastiff::logit(x_mix_baseline[[param]]) + effects_)
    } else if (param %in% c("sd_pos", "sd_neg")) {
      x_mix_overall[[param]][[pred_var]] <- 
        exp(log(x_mix_baseline[[param]]) + effects_)
    } else {
      x_mix_overall[[param]][[pred_var]] <- 
        x_mix_baseline[[param]] + effects_
    }
  }
}

if (length(p_pos_binary_pred_vars)) {
  p_pos_binary_effects <- rnorm(n = length(p_pos_binary_pred_vars),
                                mean = 0, 
                                sd = p_pos_binary_pred_vars_sd)  
  names(p_pos_binary_effects) <- p_pos_binary_pred_vars
  for (pred_var in p_pos_binary_pred_vars) {
    df_sam[[paste0("p_pos_effect_", pred_var)]] <- 
      p_pos_binary_effects[[pred_var]] * df_sam[[pred_var]]
  }
}


# Calculate each sam's x mix params given its predictors
for (param in x_mix_params) {
  if (param == "p_pos") {
    df_sam[[param]] <- mastiff::logit(x_mix_baseline[[param]]) 
  } else if (param %in% c("sd_pos", "sd_neg")) {
    df_sam[[param]] <- log(x_mix_baseline[[param]])
  } else {
    df_sam[[param]] <- x_mix_baseline[[param]]
  }
  for (pred_var in x_mix_pred_vars_names[[param]]) {
    df_sam[[param]] <- df_sam[[param]] +
      df_sam[[paste0(param, "_effect_", pred_var)]]
  }
  if (param == "p_pos") {
    for (pred_var in p_pos_binary_pred_vars) {
      df_sam[[param]] <- df_sam[[param]] +
        df_sam[[paste0("p_pos_effect_", pred_var)]]
    }
    df_sam[[param]] <- mastiff::logistic(df_sam[[param]]) 
  } else if (param %in% c("sd_pos", "sd_neg")) {
    df_sam[[param]] <- exp(df_sam[[param]])
  } 
}

# For each sam: draw x using x mix params, then calculate its mean y using its plate's
# f parameters...
df_sam <- df_sam %>%
  mutate(pos = runif(num_sam_id) < p_pos,
         xlog = if_else(pos,
                        rnorm(num_sam_id, mean = mu_pos, 
                              sd = sd_pos),
                        rnorm(num_sam_id, mean = mu_neg, 
                              sd = sd_neg)),
         x = exp(xlog),
         y_mean = PL4(xlog, f_1, f_2, f_3, f_4))

ggplot(df_sam %>% mutate(x_mix_group = paste("x mix group =", x_mix_group))) +
  geom_histogram(aes(x), fill = "grey", bins = 30) +
  scale_x_log10(limits = c(NA, NA)) +
  coord_cartesian(expand = F) +
  labs(x = "x = Ab concentration",
       y = "Number of samples") +
  facet_wrap(~x_mix_group, ncol = 8) +
  geom_vline(xintercept = exp(mu_pos)) +
  geom_vline(xintercept = exp(mu_neg))
if (FALSE) {
  ggsave("~/foo_4.pdf", height = 4, width = 5)
}

# ... then create the desired number of reps of each sample, and draw their ys
num_sam_rep <- num_sam_id * num_rep_per_sam
df_sam <- df_sam %>%
  slice(rep(row_number(), num_rep_per_sam)) %>%
  arrange(id_sam) %>%
  mutate(which_sam_rep = row_number(),
         y_obs_sd = PL4(xlog, f_1, y_obs_sd_sam_min, y_obs_sd_sam_max, f_4),
         is_blank = runif(nrow(.)) < p_blank,
         y = if_else(is_blank,
                     rnorm(nrow(.), mean = f_2,    sd = y_obs_sd_sam_min),
                     rnorm(nrow(.), mean = y_mean, sd = y_obs_sd)))

ggplot(df_sam %>% mutate(x_mix_group = paste("x mix group =", x_mix_group))) +
  geom_histogram(aes(y), fill = "grey", bins = 30) +
  scale_x_log10(limits = c(NA, NA)) +
  coord_cartesian(expand = F) +
  labs(x = "x = Ab concentration",
       y = "Number of samples") +
  facet_wrap(~x_mix_group, ncol = 8) 


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
  num_plate = num_plate,
  num_cal_tot = nrow(df_cal),
  num_sam_id = num_sam_id,
  num_sam_rep = num_sam_rep,
  which_plate_cal = df_cal$plate,
  which_plate_sam = df_sam$plate,
  which_id_sam = df_sam$id_sam,
  y_cal = df_cal$y,
  y_sam = df_sam$y,
  x_cal = df_cal$x,
  num_f_pred_vars = num_f_pred_vars,
  num_cat_per_f_pred_var = num_cat_per_f_pred_var %>% as.array(),
  sample_posterior_not_prior = 1L
)




