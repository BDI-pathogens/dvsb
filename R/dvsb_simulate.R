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

PL4 <- function(xlog, f_1, f_2, f_3, f_4) {
  f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (xlog - f_4)))
}

simulate_data <- function(
    seed = 1234567,
    num_plate = 4,
    num_sam_per_plate = 20,
    num_rep_per_sam = 2,
    num_rep_per_cal = 2 
){

# INPUT ----

set.seed(seed)

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
p_blank <- 0

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
x_mix_pred_vars$p_pos <- list(letter = letters[1:4])
x_mix_pred_vars_sds$p_pos <- c(letter = 1)
x_mix_pred_vars$mu_neg <- list(letter = letters[1:4])
x_mix_pred_vars_sds$mu_neg <- c(letter = 1)
x_mix_pred_vars$mu_pos <- list(letter = letters[1:4])
x_mix_pred_vars_sds$mu_pos <- c(letter = 1)
x_mix_pred_vars$sd_neg <- list(letter = letters[1:4])
x_mix_pred_vars_sds$sd_neg <- c(letter = 1)
x_mix_pred_vars$sd_pos <- list(letter = letters[1:4])
x_mix_pred_vars_sds$sd_pos <- c(letter = 1)

p_pos_binary_effects <- c("boolA" = -2,
                          "boolB" = 0,
                          "boolC" = 2)

y_obs_sd_min_log_shift_sd <- 0
y_obs_sd_jump_log_shift_sd <- 0

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
stopifnot(is.numeric(p_pos_binary_effects))
p_pos_binary_pred_vars <- names(p_pos_binary_effects)
if (length(p_pos_binary_effects)) {
  stopifnot(!is.null(p_pos_binary_pred_vars))
  stopifnot(!anyNA(p_pos_binary_pred_vars))
  stopifnot(!anyNA(p_pos_binary_effects))
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

y_obs_sd_min_multiplier_per_plate_unscaled <- rnorm(num_plate)
y_obs_sd_min_multiplier_per_plate_unscaled <-
  y_obs_sd_min_multiplier_per_plate_unscaled -
  mean(y_obs_sd_min_multiplier_per_plate_unscaled)
y_obs_sd_jump_multiplier_per_plate_unscaled <- rnorm(num_plate)
y_obs_sd_jump_multiplier_per_plate_unscaled <-
  y_obs_sd_jump_multiplier_per_plate_unscaled -
  mean(y_obs_sd_jump_multiplier_per_plate_unscaled)
df_plate <- df_plate %>%
  mutate(y_obs_sd_min_multiplier_per_plate =
           exp(y_obs_sd_min_multiplier_per_plate_unscaled * 
                 y_obs_sd_min_log_shift_sd - y_obs_sd_min_log_shift_sd^2 / 2),
         y_obs_sd_jump_multiplier_per_plate =
           exp(y_obs_sd_jump_multiplier_per_plate_unscaled * 
                 y_obs_sd_jump_log_shift_sd - y_obs_sd_jump_log_shift_sd^2 / 2),
         y_obs_sd_cal_min  = y_obs_sd_cal_min  * y_obs_sd_min_multiplier_per_plate,
         y_obs_sd_cal_jump = y_obs_sd_cal_jump * y_obs_sd_jump_multiplier_per_plate,
         y_obs_sd_sam_min  = y_obs_sd_sam_min  * y_obs_sd_min_multiplier_per_plate,
         y_obs_sd_sam_jump = y_obs_sd_sam_jump * y_obs_sd_jump_multiplier_per_plate)

# Label plates for plotting
df_plate$label <- paste("plate", df_plate$plate)
for (f_pred_var in f_pred_vars_names) {
  df_plate$label <- paste0(df_plate$label, ", ", f_pred_var, " ",
                           df_plate[[f_pred_var]])
}
df_plate <- df_plate %>%
  mutate(label = fct_reorder(label, plate))

# Expand to one row per cal (one for each x). Calculate y expected.
df_cal <- df_plate %>%
  expand_grid(xlog = xlogs, cal = 1:num_rep_per_cal) %>%
  mutate(x = exp(xlog),
         which_cal = row_number(),
         y_mean = PL4(xlog, f_1, f_2, f_3, f_4))

# Draw observed y
df_cal <- df_cal %>%
  mutate(y_obs_sd = PL4(xlog, f_1, y_obs_sd_cal_min, y_obs_sd_cal_min + y_obs_sd_cal_jump, f_4),
         y = rnorm(nrow(.), mean = y_mean, sd = y_obs_sd))

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
  mutate(id_sam = as.character(row_number()))
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
for (param in x_mix_params) {
  x_mix_effects[[param]] <- list()
  for (pred_var in x_mix_pred_vars_names[[param]]) {
    sigma_ <- x_mix_pred_vars_sds[[param]][[pred_var]]
    num_cats <- x_mix_pred_vars_num_cats[[param]][[pred_var]]
    effects_ <- rnorm(num_cats, 0, sigma_)
    effects_ <- effects_ - mean(effects_)
    names(effects_) <- x_mix_pred_vars[[param]][[pred_var]]
    x_mix_effects[[param]][[pred_var]] <- effects_
    df_sam[[paste0(param, "_effect_", pred_var)]] <- map_dbl(
      df_sam[[pred_var]], ~ effects_[[.x]])
  }
}

if (length(p_pos_binary_effects)) {
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

# ... then create the desired number of reps of each sample, and draw their ys
num_sam_rep <- num_sam_id * num_rep_per_sam
df_sam <- df_sam %>%
  slice(rep(row_number(), num_rep_per_sam)) %>%
  arrange(id_sam) %>%
  mutate(y_obs_sd = PL4(xlog, f_1, y_obs_sd_sam_min, y_obs_sd_sam_max, f_4),
         is_blank = runif(nrow(.)) < p_blank,
         y = if_else(is_blank,
                     rnorm(nrow(.), mean = f_2,    sd = y_obs_sd_sam_min),
                     rnorm(nrow(.), mean = y_mean, sd = y_obs_sd)))

return(list(
  df_sam = df_sam,
  df_plate = df_plate,
  df_cal = df_cal,
  f_pred_vars_names = f_pred_vars_names,
  x_mix_pred_vars_names = x_mix_pred_vars_names,
  p_pos_binary_pred_vars = p_pos_binary_pred_vars,
  params = list(
  y_obs_sd_cal_min = y_obs_sd_cal_min,
  y_obs_sd_cal_max = y_obs_sd_cal_max,
  y_obs_sd_cal_jump = y_obs_sd_cal_jump,
  y_obs_sd_sam_min = y_obs_sd_sam_min,
  y_obs_sd_sam_max = y_obs_sd_sam_max,
  y_obs_sd_sam_jump = y_obs_sd_sam_jump,
  mu_neg = mu_neg,
  sd_neg = sd_neg,
  mu_pos = mu_pos,
  sd_pos = sd_pos,
  p_pos = p_pos,
  p_blank = p_blank,
  f = f,
  sigma_f_plate = sigma_f_plate,
  rho = rho,
  x_mix_pred_vars_sds = x_mix_pred_vars_sds,
  x_mix_effects = x_mix_effects,
  f_effects_by_pred_var = f_effects_by_pred_var,
  sigma_f_pred_vars = sigma_f_pred_vars,
  p_pos_binary_effects = p_pos_binary_effects,
  y_obs_sd_min_log_shift_sd = y_obs_sd_min_log_shift_sd,
  y_obs_sd_jump_log_shift_sd = y_obs_sd_jump_log_shift_sd
  )))

}
