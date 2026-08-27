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

#' Simulate x and y values for samples and calibrators across different plates
#'
#' @param seed seed used for random number generation (if no seed is provided,
#'   none is used)
#' @param num_plate number of plates
#' @param num_sam_per_plate number of samples per plate
#' @param num_rep_per_sam number of replicates per sample
#' @param num_rep_per_cal number of replicates per calibrator
#' @param x_cals the x values of the set of calibrators on each plate
#' @param y_obs_sd_cal_min the limit, as x tends to zero, of the scale of
#'   stochastic observational noise in y for calibrators
#' @param y_obs_sd_cal_jump the difference between the lower and upper limits,
#'   as x tends to zero and infinity, of the scale of stochastic observational
#'   noise in y for calibrators
#' @param y_obs_sd_sam_min the limit, as x tends to zero, of the scale of
#'   stochastic observational noise in y for samples
#' @param y_obs_sd_sam_jump the difference between the lower and upper limits,
#'   as x tends to zero and infinity, of the scale of stochastic observational
#'   noise in y for samples
#' @param mu_neg the population mean log_e(x) for seronegatives
#' @param sd_neg the standard deviation of the population distribution of
#'   log_e(x) for seronegatives
#' @param mu_pos the population mean log_e(x) for seropositives
#' @param sd_pos the standard deviation of the population distribution of
#'   log_e(x) for seropositives
#' @param p_pos the probability of a sample being seropositive (i.e.
#'   seroprevalence) before addition of subpopulation-specific deviations
#' @param p_blank the probability that any given sample replicate is an
#'   accidental blank. Beware: in the main dvsb inference model this is assumed
#'   to be zero, so values greater than zero introduce model misspecification
#'   for inference. The dvsb_accidental_blanks inference model does not assume
#'   `p_blank` is zero.
#' @param f a vector with the four parameters of the four-parameter logistic
#'   (4PL) function that link log_e(x) to the expected value of y through y =
#'   f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (log_e(x) - f_4)))
#' @param sigma_f_plate a vector with the four scales of normal variability
#'   between plates of the four elements of the f vector
#' @param rho 4x4 correlation matrix for the variability between plates of the
#'   four elements of the f vector
#' @param f_pred_vars a list (whose names are the names of categorical
#'   variables) of character vectors (whose values are the different categories
#'   of a given categorical variable). Each plate has one category randomly
#'   allocated for each categorical variable. Categories differ systematically
#'   in their f vector (in addition to the random variability between plates).
#' @param sigma_f_pred_vars a list (whose names are the names of categorical
#'   variables, matching those of `f_pred_vars`) of length-4 numeric vectors.
#'   Each of these vectors specifies the scales of normal variability in the f
#'   vector between different categories of the corresponding categorical
#'   variable.
#' @param x_mix_pred_vars a list (whose names are the parameters p_pos, mu_neg,
#'   mu_pos, sd_neg, sd_pos) of lists (whose names are the names of categorical
#'   variables) of character vectors (whose values are the different categories
#'   of a given categorical variable). Each parameter has a regression model
#'   specified by additively combining its categorical variables. For each
#'   categorical variable, each sample has one of the categories randomly
#'   allocated.
#' @param x_mix_pred_vars_sds a list (whose names are the parameters p_pos,
#'   mu_neg, mu_pos, sd_neg, sd_pos) of named numeric vectors (whose names must
#'   match the categorical variables named in the inner lists of
#'   `x_mix_pred_vars`). Each \{name, numeric value\} pair within one of these
#'   vectors specifies the scale of variability between the regression
#'   coefficients for the different categories of named categorical variable.
#' @param x_mix_effects Normally you will want to leave this at its default
#'   value of NA, in which case we will randomly draw regression coefficients
#'   for the parameters of the regression model for the x distribution according
#'   to the values specified in `x_mix_pred_vars_sds`. Alternatively, this
#'   argument can be used to specify values for all these coefficients; this may
#'   be useful to keep them fixed over several different simulations of a
#'   dataset, e.g. with different stochastic seeds, to compare inference with
#'   these random effects fixed. If used, this argument should be a list (whose
#'   names are the x mixture parameters p_pos, mu_neg, mu_pos, sd_neg, sd_pos)
#'   of lists (whose names match those of the same parameter in x_mix_pred_vars)
#'   of named numeric vectors (one per categorical variable used in a regression
#'   model for this parameter; the names, one per category of this variable,
#'   must match those for the variable as specified in x_mix_pred_vars).
#' @param p_pos_binary_effects a named numeric vector. Each \{name, numeric
#'   value\} pair within this vector specifies the name of a logical variable and
#'   the additive shift in seroprevalence (on a logit scale) between when this
#'   variable is true and when it is false. Each sample will be randomly
#'   allocated a value of true or false for each such variable.
#' @param y_obs_sd_min_log_shift_sd the standard deviation (on a log scale) of
#'   the multiplicative variability in both y_obs_sd_cal_min and
#'   y_obs_sd_sam_min between plates. Beware: such variability is assumed to be
#'   zero in the inference model, so values greater than the default of zero
#'   will introduce model misspecification for inference (which may be of
#'   interest for testing purposes).
#' @param y_obs_sd_jump_log_shift_sd the standard deviation (on a log scale) of
#'   the multiplicative variability in both y_obs_sd_cal_jump and
#'   y_obs_sd_sam_jump between plates. Beware: such variability is assumed to be
#'   zero in the inference model, so values greater than the default of zero
#'   will introduce model misspecification for inference (which may be of
#'   interest for testing purposes).
#'
#' @returns a named list whose elements are: df_sam (a dataframe with one row
#'   per simulated sample replicate), df_plate (a dataframe with one row per
#'   simulated plate), df_cal (a dataframe with one row per simulated calibrator
#'   replicate), params (a list of the values of parameters used for
#'   simulation), f_pred_vars_names (a character vector of the names of any
#'   variables used for a regression model for the four vector f),
#'   x_mix_pred_vars_names (a named list whose names are the five parameters of
#'   the normal mixture model for x; each one of the elements is a character
#'   vector naming the variables used in a regression model for that parameter),
#'   p_pos_binary_pred_vars (a character vector of the names of any variables
#'   used for a regression model for p_pos, using only binary fixed effects).
#' @importFrom magrittr %>%
#' @export
#'
simulate_data <- function(
    seed = NA,
    num_plate = 10,
    num_sam_per_plate = 40,
    num_rep_per_sam = 2,
    num_rep_per_cal = 2,
    x_cals = c(0, 0.37, 1.1, 3.3, 10, 30),
    y_obs_sd_cal_min = 0.003,
    y_obs_sd_cal_jump = 0.242,
    y_obs_sd_sam_min = 0.004,
    y_obs_sd_sam_jump = 0.718,
    mu_neg = -2.7,
    sd_neg = 1,
    mu_pos = 0.7,
    sd_pos = 1,
    p_pos = 0.5,
    p_blank = 0,
    f = c(0.95,
          0.01,
          3.33,
          2.44),
    sigma_f_plate = c(0.073,
                      0.004,
                      0.704,
                      0.404),
    rho = matrix(c(1, 0, 0, 0,
                   0, 1, 0, 0,
                   0, 0, 1, 0,
                   0, 0, 0, 1),
                 4, 4, byrow = TRUE),
    f_pred_vars = list(),
    sigma_f_pred_vars = list(),
    x_mix_pred_vars = list(
      p_pos  = list(letter = letters[1:4]),
      mu_neg = list(letter = letters[1:4]),
      mu_pos = list(letter = letters[1:4]),
      sd_neg = list(letter = letters[1:4]),
      sd_pos = list(letter = letters[1:4])
    ),
    x_mix_pred_vars_sds = list(
      p_pos = c(letter = 1),
      mu_neg = c(letter = 1),
      mu_pos = c(letter = 1),
      sd_neg = c(letter = 0.2),
      sd_pos = c(letter = 0.2)
    ),
    x_mix_effects = NA,
    p_pos_binary_effects = c("boolA" = -2,
                             "boolB" = 0,
                             "boolC" = 2),
    y_obs_sd_min_log_shift_sd = 0,
    y_obs_sd_jump_log_shift_sd = 0
){
  
  # INPUT CHECKS ----
  
  x_mix_params <- c("mu_neg", "mu_pos", "p_pos", "sd_neg", "sd_pos")
  
  # Check numeric scalars  
  if (! identical(seed, NA)) mastiff::check_numeric(seed)
  mastiff::check_numeric(num_plate, lower = 0)
  mastiff::check_numeric(num_sam_per_plate, lower = 0)
  mastiff::check_numeric(num_rep_per_sam, lower = 0)
  mastiff::check_numeric(num_rep_per_cal, lower = 0)
  mastiff::check_numeric(y_obs_sd_cal_min, lower = 0)
  mastiff::check_numeric(y_obs_sd_cal_jump, lower = 0)
  mastiff::check_numeric(y_obs_sd_sam_min, lower = 0)
  mastiff::check_numeric(y_obs_sd_sam_jump, lower = 0)
  mastiff::check_numeric(mu_neg)
  mastiff::check_numeric(mu_pos, lower = mu_neg, lower_inclusive = FALSE)
  mastiff::check_numeric(sd_neg, lower = 0)
  mastiff::check_numeric(sd_pos, lower = 0)
  mastiff::check_numeric(p_pos, lower = 0, upper = 1)
  mastiff::check_numeric(p_blank, lower = 0, upper = 1)
  mastiff::check_numeric(y_obs_sd_min_log_shift_sd, lower = 0)
  mastiff::check_numeric(y_obs_sd_jump_log_shift_sd, lower = 0)
  
  # Check vectors and matrices
  stopifnot(is.numeric(x_cals))
  stopifnot(!anyDuplicated(x_cals))
  stopifnot(all(x_cals >= 0))
  stopifnot(is.numeric(f))
  stopifnot(length(f) == 4)
  stopifnot(is.numeric(sigma_f_plate))
  stopifnot(length(sigma_f_plate) == 4)
  stopifnot(all(sigma_f_plate >= 0))
  stopifnot(is.matrix(rho))
  stopifnot(is.numeric(rho))
  stopifnot(isSymmetric(rho))
  stopifnot(all(diag(rho) == 1))
  stopifnot(all(rho >= -1))
  stopifnot(all(rho <= 1))
  
  # Then check more complicated objects...
  
  # Check f_pred_vars and sigma_f_pred_vars
  f_pred_vars_names <- names(f_pred_vars)
  stopifnot(! any(f_pred_vars_names %in% # avoid name clashes with variables
                    c("plate", "f", "f_1", "f_2", "f_3", "f_4")))
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
      stopifnot(is.numeric(sigma_f_pred_vars[[name_]]))
      stopifnot(length(sigma_f_pred_vars[[name_]]) == 4)
      stopifnot(all(sigma_f_pred_vars[[name_]] >= 0))
    }
  }
  
  # Check x mix pred vars 
  num_sam_id <- num_plate * num_sam_per_plate
  stopifnot(identical(sort(names(x_mix_pred_vars)),     x_mix_params))
  stopifnot(identical(sort(names(x_mix_pred_vars_sds)), x_mix_params))
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
  x_mix_pred_vars_names <- purrr::map(x_mix_pred_vars, names)
  x_mix_pred_vars_nums <- purrr::map_int(x_mix_pred_vars_names, length)
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
        if (! (is.numeric(x_mix_pred_vars_sds[[param_]][[pred_var]]) &&
               length(x_mix_pred_vars_sds[[param_]][[pred_var]]) == 1 &&
               is.finite(x_mix_pred_vars_sds[[param_]][[pred_var]]) &&
               x_mix_pred_vars_sds[[param_]][[pred_var]] >= 0 )) {
          stop(paste('x_mix_pred_vars_sds[["', param_, '"]][["', pred_var,
                     '"]] must be a non-negative number'))
        }
      }
    }
  }
  
  # Check a manually specified x_mix_effects
  if (! identical(x_mix_effects, NA)) {
    stopifnot(is.list(x_mix_effects))
    stopifnot(identical(sort(names(x_mix_effects)),
                        sort(x_mix_params)))
    for (param in x_mix_params) {
      stopifnot(is.list(x_mix_effects[[param]]))
      stopifnot(identical(sort(names(x_mix_effects[[param]])),
                          sort(names(x_mix_pred_vars[[param]]))))
      for (param_inner in names(x_mix_effects[[param]])) {
        vector_ <- x_mix_effects[[param]][[param_inner]]
        stopifnot(is.numeric(vector_))
        stopifnot(identical(sort(names(vector_)),
                            sort(x_mix_pred_vars[[param]][[param_inner]])))
        stopifnot(all(is.finite(vector_)))
        stopifnot(abs(sum(vector_)) < 1e-5)
      }
    }
  }  
  
  # Check that if the same predictor variable is specified for different params,
  # it has the same categories, to avoid confusing output.
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
  
  if (! identical(seed, NA)) set.seed(seed)
  
  # Derived params
  y_obs_sd_cal_max <- y_obs_sd_cal_min + y_obs_sd_cal_jump
  y_obs_sd_sam_max <- y_obs_sd_sam_min + y_obs_sd_sam_jump
  
  x_cals_log <- log(x_cals)
  
  
  # Make a df with one row per plate.
  # Sample each plate's f predictor variables.
  # Ensure that we don't randomly sample the same category for every plate.
  # Delete any unsampled categories.
  df_plate <- tibble::tibble(plate = 1:num_plate)
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
  num_cat_per_f_pred_var <- purrr::map_int(f_pred_vars, length)
  num_f_pred_var_cats <- sum(num_cat_per_f_pred_var)
  
  # Draw plate-level variation in f
  Sigma_plate <- diag(sigma_f_plate) %*% rho %*% diag(sigma_f_plate)
  f_plate_effects <- mvtnorm::rmvnorm(num_plate, c(0, 0, 0, 0), Sigma_plate)
  df_plate$f_effect_plate <- purrr::map(1:num_plate, ~ f_plate_effects[.x, ])
  
  # Draw variation in f due to f_pred_vars
  f_effects_by_pred_var <- list()
  for (f_pred_var in f_pred_vars_names) {
    sigma_f_ <- sigma_f_pred_vars[[f_pred_var]]
    Sigma_f_ <- diag(sigma_f_) %*% rho %*% diag(sigma_f_)
    num_cats <- length(f_pred_vars[[f_pred_var]])
    f_effects_ <- mvtnorm::rmvnorm(num_cats, c(0, 0, 0, 0), Sigma_f_)
    f_effects_col_means <- colMeans(f_effects_)
    for (cat_num in 1:num_cats) {
      f_effects_[cat_num, ] <- f_effects_[cat_num, ] - f_effects_col_means
    }
    rownames(f_effects_) <- f_pred_vars[[f_pred_var]]
    f_effects_by_pred_var[[f_pred_var]] <- f_effects_
    df_plate[[paste0("f_effect_", f_pred_var)]] <- purrr::map(
      df_plate[[f_pred_var]], ~ f_effects_[.x, ])
  }
  
  # Assign f by plate
  df_plate$f <- purrr::map(1:num_plate, ~ f)
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
    tidyr::unnest_wider(f, names_sep = "_")
  
  y_obs_sd_min_multiplier_per_plate_unscaled <- stats::rnorm(num_plate)
  y_obs_sd_min_multiplier_per_plate_unscaled <-
    y_obs_sd_min_multiplier_per_plate_unscaled -
    mean(y_obs_sd_min_multiplier_per_plate_unscaled)
  y_obs_sd_jump_multiplier_per_plate_unscaled <- stats::rnorm(num_plate)
  y_obs_sd_jump_multiplier_per_plate_unscaled <-
    y_obs_sd_jump_multiplier_per_plate_unscaled -
    mean(y_obs_sd_jump_multiplier_per_plate_unscaled)
  df_plate <- df_plate %>%
    dplyr::mutate(y_obs_sd_min_multiplier_per_plate =
             exp(y_obs_sd_min_multiplier_per_plate_unscaled * 
                   y_obs_sd_min_log_shift_sd - y_obs_sd_min_log_shift_sd^2 / 2),
           y_obs_sd_jump_multiplier_per_plate =
             exp(y_obs_sd_jump_multiplier_per_plate_unscaled * 
                   y_obs_sd_jump_log_shift_sd - y_obs_sd_jump_log_shift_sd^2 / 2),
           y_obs_sd_cal_min  = y_obs_sd_cal_min  * y_obs_sd_min_multiplier_per_plate,
           y_obs_sd_cal_jump = y_obs_sd_cal_jump * y_obs_sd_jump_multiplier_per_plate,
           y_obs_sd_sam_min  = y_obs_sd_sam_min  * y_obs_sd_min_multiplier_per_plate,
           y_obs_sd_sam_jump = y_obs_sd_sam_jump * y_obs_sd_jump_multiplier_per_plate)
  
  # Expand to one row per cal (one for each x). Calculate y expected.
  df_cal <- df_plate %>%
    tidyr::expand_grid(xlog = x_cals_log, cal = 1:num_rep_per_cal) %>%
    dplyr::mutate(x = exp(xlog),
           which_cal = dplyr::row_number(),
           y_mean = PL4(xlog, f_1, f_2, f_3, f_4))
  
  # Draw observed y
  df_cal <- df_cal %>%
    dplyr::mutate(y_obs_sd = PL4(xlog, f_1, y_obs_sd_cal_min, y_obs_sd_cal_min + y_obs_sd_cal_jump, f_4),
           y = stats::rnorm(nrow(.), mean = y_mean, sd = y_obs_sd))
  
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
    dplyr::slice(rep(dplyr::row_number(), num_sam_per_plate)) %>%
    dplyr::arrange(plate) %>%
    dplyr::mutate(id_sam = as.character(dplyr::row_number()))
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
      while(dplyr::n_distinct(sampled_pred_vars) < 2) {
        sampled_pred_vars <- sample(x_mix_pred_vars[[param]][[pred_var]],
                                    size = num_sam_id,
                                    replace = TRUE)
      } 
      x_mix_pred_vars[[param]][[pred_var]] <- sort(unique(sampled_pred_vars))
      df_sam[[pred_var]] <- sampled_pred_vars
    }
    x_mix_pred_vars_num_cats[[param]] <-
      purrr::map_int(x_mix_pred_vars[[param]], length)
    x_mix_pred_vars_num_cats_tots[[param]] <- 
      sum(x_mix_pred_vars_num_cats[[param]])
  }
  for (pred_var in p_pos_binary_pred_vars) {
    sampled_pred_vars <- character()
    while(dplyr::n_distinct(sampled_pred_vars) < 2) {
      sampled_pred_vars <- sample(c(TRUE, FALSE),
                                  size = num_sam_id,
                                  replace = TRUE)
    } 
    df_sam[[pred_var]] <- sampled_pred_vars
  }
  all_names_x_mix_pred_vars <- x_mix_pred_vars %>% 
    purrr::map(names) %>%
    unlist() %>%
    unique()
  if (is.null(all_names_x_mix_pred_vars)) {
    df_sam <- df_sam %>%
      dplyr::mutate(x_mix_group = NA_character_)
  } else {
    df_sam <- df_sam %>%
      tidyr::unite("x_mix_group", tidyselect::all_of(all_names_x_mix_pred_vars), sep = "_", remove = FALSE)
  }
  for (pred_var in p_pos_binary_pred_vars) {
    df_sam$x_mix_group <- paste0(df_sam$x_mix_group, "_", pred_var, df_sam[[pred_var]])
  }
  
  # Draw effects on the x mix params from each pred var
  if (identical(x_mix_effects, NA)) {
    x_mix_effects <- list()
    for (param in x_mix_params) {
      x_mix_effects[[param]] <- list()
      for (pred_var in x_mix_pred_vars_names[[param]]) {
        sigma_ <- x_mix_pred_vars_sds[[param]][[pred_var]]
        num_cats <- x_mix_pred_vars_num_cats[[param]][[pred_var]]
        effects_ <- stats::rnorm(num_cats, 0, sigma_)
        effects_ <- effects_ - mean(effects_)
        names(effects_) <- x_mix_pred_vars[[param]][[pred_var]]
        x_mix_effects[[param]][[pred_var]] <- effects_
      }
    }
  } 
  for (param in x_mix_params) {
    for (pred_var in x_mix_pred_vars_names[[param]]) {
      effects_ <- x_mix_effects[[param]][[pred_var]]
      df_sam[[paste0(param, "_effect_", pred_var)]] <- purrr::map_dbl(
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
    dplyr::mutate(pos = stats::runif(num_sam_id) < p_pos,
           xlog = dplyr::if_else(pos,
                          stats::rnorm(num_sam_id, mean = mu_pos, 
                                sd = sd_pos),
                          stats::rnorm(num_sam_id, mean = mu_neg, 
                                sd = sd_neg)),
           x = exp(xlog),
           y_mean = PL4(xlog, f_1, f_2, f_3, f_4))
  
  # ... then create the desired number of reps of each sample, and draw their ys
  num_sam_rep <- num_sam_id * num_rep_per_sam
  df_sam <- df_sam %>%
    dplyr::slice(rep(dplyr::row_number(), num_rep_per_sam)) %>%
    dplyr::arrange(id_sam) %>%
    dplyr::mutate(y_obs_sd = PL4(xlog, f_1, y_obs_sd_sam_min, y_obs_sd_sam_max, f_4),
           is_blank = stats::runif(nrow(.)) < p_blank,
           y = dplyr::if_else(is_blank,
                       stats::rnorm(nrow(.), mean = f_2,    sd = y_obs_sd_sam_min),
                       stats::rnorm(nrow(.), mean = y_mean, sd = y_obs_sd)))
  
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
