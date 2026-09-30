#' Wrangles inputs into the format required for [run_stan_interfaces()]
#'
#' @param df_sam a dataframe of the format output by [simulate_data()] (inside
#'   its list of outputs): with one row per sample replicate, including columns
#'   id_sam, y, and plate. y should be numeric.
#' @param df_cal a dataframe of the format output by [simulate_data()] (inside
#'   its list of outputs): with one row per calibrator replicate, including
#'   columns x, y, and plate. x and y should be numeric.
#' @param df_priors_scalars,df_priors_vectors dataframes of the format output by
#'   [read_priors()] (inside its list of outputs): each having one row per
#'   parameter, including columns `param`, `lower`, and `upper`. `param` should
#'   be character, `lower` and `upper` should be numeric.
#' @param rho_prior_eta the single positive number 'eta' that is a
#'   hyperparameter of the LKJ distribution used as the prior for the
#'   correlation matrix for variability in the parameters of the x<->y
#'   relationship between plates. See the Stan docs on the LKJ distribution for
#'   the definition of eta.
#' @param x_mix_pred_vars_names a list (whose names are the parameters p_pos,
#'   mu_neg, mu_pos, sd_neg, sd_pos) of character vectors (each one naming the
#'   categorical predictor variables used in the regression model for that
#'   parameter). The default value of a list with all empty character vectors is
#'   appropriate for using no regression model for any of these five parameters.
#' @param p_pos_binary_pred_vars a character vector containing the names of
#'   binary variables used as fixed effects in the regression mode for
#'   seroprevalence. The default value of an empty character vector is
#'   appropriate for using no such regression model.
#' @param f_pred_vars_names a character vector containing the names of
#'   variables used for a regression model for the four vector f. The default
#'   value of an empty character vector is appropriate for using no regression
#'   model.
#'
#' @returns a list whose elements are
#'    * `df_plate`, a dataframe with one row per plate;
#'    * `data_descriptors`, a list of things required as input for [rename_params_from_stan()] and [wrangle_true_params()];
#'    * `stan_input_posterior`, a list of things required as input for [run_stan_interfaces()] if you want to sample from the posterior;
#'    * `stan_input_prior`, a list of things required as input for [run_stan_interfaces()] if you want to sample from the prior (this reduces the number of samples to zero but keeps the same 'shape' of the dataset in terms of covariates used, to more efficiently sample population-level parameters while avoiding sampling individual-level parameters such as x);
#'    * `df_sam` and `df_cal` (as provided as input but with extra columns added for Stan indexing).
#' @export
#'
prepare_data_for_stan <- function(
    df_sam,
    df_cal,
    df_priors_scalars,
    df_priors_vectors,
    rho_prior_eta = 1,
    x_mix_pred_vars_names = list(
      p_pos  = character(),
      mu_pos = character(),
      sd_pos = character(),
      mu_neg = character(),
      sd_neg = character()
    ),
    p_pos_binary_pred_vars = character(),
    f_pred_vars_names = character()
    ) {
  
  x_mix_params <- c("p_pos", "mu_neg", "mu_pos", "sd_neg", "sd_pos")
  
  # CHECK ARGS ----
  
  stopifnot(is.data.frame(df_sam))
  stopifnot(is.data.frame(df_cal))
  stopifnot(is.character(f_pred_vars_names))
  stopifnot("plate" %in% names(df_sam))
  stopifnot("plate" %in% names(df_cal))
  stopifnot(all(f_pred_vars_names %in% names(df_sam)))
  stopifnot(all(f_pred_vars_names %in% names(df_cal)))
  stopifnot(!anyDuplicated(f_pred_vars_names))
  stopifnot(is.list(x_mix_pred_vars_names))
  stopifnot(all(names(x_mix_pred_vars_names) %in% x_mix_params))
  for (param in x_mix_params) {
    stopifnot(is.null(x_mix_pred_vars_names[[param]]) ||
                is.character(x_mix_pred_vars_names[[param]]))
  }
  stopifnot(is.character(p_pos_binary_pred_vars))
  stopifnot(all(p_pos_binary_pred_vars %in% names(df_sam)))
  mastiff::check_numeric(rho_prior_eta, lower = 0)
  stopifnot(is.data.frame(df_priors_scalars))
  stopifnot(is.data.frame(df_priors_vectors))
  cols_expected_sam <- c("y", "plate", "id_sam")
  for (col in cols_expected_sam) {
    if (! col %in% names(df_sam)) stop(paste("Column", col,
                                             "missing from df_sam"))
  }
  cols_expected_cal <- c("x", "y", "plate")
  for (col in cols_expected_cal) {
    if (! col %in% names(df_cal)) stop(paste("Column", col,
                                             "missing from df_cal"))
  }
  cols_unexpected_sam <- c("plate_int", "id_sam_int", "which_sam_rep")
  for (col in cols_unexpected_sam) {
    if (col %in% names(df_sam)) stop(paste(
      "Column", col, "already present in df_sam; please remove it to prevent",
    "overwriting"))
  }
  cols_unexpected_cal <- c("plate_int")
  for (col in cols_unexpected_cal) {
    if (col %in% names(df_cal)) stop(paste(
      "Column", col, "already present in df_cal; please remove it to prevent",
      "overwriting"))
  }
  stopifnot(is.numeric(df_sam$y))
  stopifnot(is.numeric(df_cal$y))
  stopifnot(is.numeric(df_cal$x))
  stopifnot(!anyNA(df_sam$y))
  stopifnot(!anyNA(df_cal$y))
  stopifnot(!anyNA(df_cal$x))
  stopifnot(!anyNA(df_sam$plate))
  stopifnot(!anyNA(df_cal$plate))
  cols_expected_df_priors <- c("param", "lower", "upper")
  for (col in cols_expected_df_priors) {
    if (! col %in% names(df_priors_scalars)) stop(paste(
      "Column", col, "missing from df_priors_scalars"))
    if (! col %in% names(df_priors_vectors)) stop(paste(
      "Column", col, "missing from df_priors_vectors"))
  }
  # TODO: check df_priors_scalars and df_priors_vectors have all the rows needed
  
  # GENERAL WRANGLING ----
  
  df <- dplyr::bind_rows(df_sam, df_cal)
  
  df_plate <- data.frame(plate = unique(df$plate))
  df_plate$plate_int <- seq(from = 1, length.out = nrow(df_plate))
  num_plate <- nrow(df_plate)
  df_sam <- dplyr::left_join(df_sam, df_plate, by = "plate")
  df_cal <- dplyr::left_join(df_cal, df_plate, by = "plate")
  
  num_sam_rep <- nrow(df_sam)
  
  # TODO: set int id cols 
  df_sam <- df_sam %>%
    dplyr::mutate(.by = id_sam, id_sam_int = dplyr::cur_group_id()) %>% 
    dplyr::arrange(id_sam_int) %>%
    dplyr::mutate(which_sam_rep = dplyr::row_number())
  
  num_sam_id <- ifelse(num_sam_rep, max(df_sam$id_sam_int), 0)
  
  # WRANGLE x_mix_pred_vars ----
  
  x_mix_pred_vars <- list()
  x_mix_pred_vars_nums <- integer()
  x_mix_pred_vars_num_cats <- list()
  x_mix_pred_vars_num_cats_tots <- integer()
  for (param in x_mix_params) {
    x_mix_pred_vars[[param]] <- list()
    x_mix_pred_vars_num_cats[[param]] <- numeric()
    x_mix_pred_vars_num_cats_tots[[param]] <- 0L
    if (is.null(x_mix_pred_vars_names[[param]])) {
      x_mix_pred_vars_nums[[param]] <- 0L
      next
    }
    x_mix_pred_vars_nums[[param]] <- length(x_mix_pred_vars_names[[param]])
    for (pred_var in x_mix_pred_vars_names[[param]]) {
      stopifnot(pred_var %in% names(df_sam))
      if (anyNA(df_sam[[pred_var]])) {
        stop(paste("Some samples have missing", pred_var))
      }
      if (dplyr::n_distinct(df_sam[[pred_var]]) < 2L) stop(paste(
        "Fewer than 2 unique values of", pred_var, "found in df_sam;",
        "need at least 2 to use it as an x_mix_pred_var"
      ))
      x_mix_pred_vars[[param]][[pred_var]] <- unique(df_sam[[pred_var]])
    }
    x_mix_pred_vars_num_cats[[param]] <-
      purrr::map_int(x_mix_pred_vars[[param]], length)
    x_mix_pred_vars_num_cats_tots[[param]] <- 
      sum(x_mix_pred_vars_num_cats[[param]])
  }
  names(x_mix_pred_vars_nums) <- x_mix_params
  predict_p_pos  <- x_mix_pred_vars_nums[["p_pos"]]  > 0L
  predict_mu_pos <- x_mix_pred_vars_nums[["mu_pos"]] > 0L
  predict_sd_pos <- x_mix_pred_vars_nums[["sd_pos"]] > 0L
  predict_mu_neg <- x_mix_pred_vars_nums[["mu_neg"]] > 0L
  predict_sd_neg <- x_mix_pred_vars_nums[["sd_neg"]] > 0L
  
  # For each of the 5 x mix params, 0-or-1 encode each cat of each pred var,
  # for every sample
  x_mix_design_matrices <- list()
  for (param in x_mix_params) {
    
    # No pred vars
    if (x_mix_pred_vars_nums[[param]] == 0L) {
      x_mix_design_matrices[[param]] <- matrix(NA_real_, nrow = num_sam_id, ncol = 0)
      
      # Some pred vars
    } else {
      mat <- df_sam %>%
        dplyr::select(id_sam_int, tidyselect::all_of(x_mix_pred_vars_names[[param]])) %>%
        dplyr::distinct() # remove duplicates because of replicates...
      stopifnot(identical(mat$id_sam_int, # ... and check one row per sam_id, in order
                          1:nrow(mat)))
      mat <- mat %>%
        dplyr::select(-id_sam_int) %>%
        dplyr::mutate(dplyr::across(tidyselect::everything(), as.character)) %>%
        dplyr::mutate(dplyr::across(tidyselect::everything(), as.factor)) %>%
        {stats::model.matrix(~ . - 1,
                      data = .,
                      contrasts.arg = lapply(.[, , drop = FALSE],
                                             stats::contrasts, contrasts = FALSE))}
      colnames_expected <-
        purrr::map(x_mix_pred_vars_names[[param]],
            ~ paste0(.x, x_mix_pred_vars[[param]][[.x]])) %>%
        unlist
      stopifnot(identical(sort(colnames(mat)),
                          sort(colnames_expected)))
      mat <- mat[, colnames_expected]
      stopifnot(identical(colnames(mat),
                          colnames_expected))
      x_mix_design_matrices[[param]] <- mat
    }
  }
  
  # Look-ups between int and string encodings of x_mix_pred_vars 
  lookup_pred_var_int <- list()
  lookup_pred_var_cat_int <- list()
  for (param in x_mix_params) {
    lookup_pred_var_int[[param]] <- tibble::tibble(
      pred_var = x_mix_pred_vars_names[[param]],
      pred_var_int = seq(from = 1, length.out = x_mix_pred_vars_nums[[param]]))
    lookup_pred_var_cat_int[[param]] <- tibble::tibble(
      pred_var_cat = colnames(x_mix_design_matrices[[param]]),
      pred_var_cat_int = seq(from = 1, length.out = x_mix_pred_vars_num_cats_tots[[param]]))
    if (nrow(lookup_pred_var_int[[param]]) == 0) {
      lookup_pred_var_int[[param]]$pred_var <- character()
    }
    if (nrow(lookup_pred_var_cat_int[[param]]) == 0) {
      lookup_pred_var_cat_int[[param]]$pred_var_cat <- character()
    }
  }

  # WRANGLE p_pos_binary_pred_vars ----
  
  # Check that we have both TRUE and FALSE values for each p_pos_binary_pred_var
  for (pred_var in p_pos_binary_pred_vars) {
    stopifnot(is.logical(df_sam[[pred_var]]))
    if (dplyr::n_distinct(df_sam[[pred_var]]) < 2L) stop(paste(
      "Fewer than 2 unique values of", pred_var, "found in df_sam;",
      "need at least 2 to use it as a p_pos_binary_pred_var"
    ))
  }
  predict_p_pos_binary <- length(p_pos_binary_pred_vars) > 0L
  
  # Encode the design matrix for p_pos_binary 
  if (! predict_p_pos_binary) {
    p_pos_binary_design_matrix <- matrix(NA_real_, nrow = num_sam_id, ncol = 0)
  } else {
    mat <- df_sam %>%
      dplyr::select(id_sam_int, tidyselect::all_of(p_pos_binary_pred_vars)) %>%
      dplyr::distinct() 
    stopifnot(identical(mat$id_sam_int, 
                        1:nrow(mat)))
    mat <- mat %>%
      dplyr::select(-id_sam_int) %>%
      dplyr::mutate(dplyr::across(tidyselect::everything(), as.integer))
    stopifnot(identical(colnames(mat),
                        p_pos_binary_pred_vars))
    p_pos_binary_design_matrix <- mat
  }
  
  # WRANGLE f_pred_vars ----

  num_f_pred_vars <- length(f_pred_vars_names)
  predict_f <- num_f_pred_vars > 0L
  f_pred_vars <- list()
  for (f_pred_var in f_pred_vars_names) {
    if (anyNA(df[[f_pred_var]])) {
      stop(paste("Some samples or calibrators have missing", f_pred_var))
    }
    if (dplyr::n_distinct(df[[f_pred_var]]) < 2L) {
      stop(paste("Found fewer than 2 unique values for f_pred_var",
                 f_pred_var, "across df_sam and df_cal. At least 2 unique",
                 "values are required."))
    }
    f_pred_vars[[f_pred_var]] <- unique(df[[f_pred_var]])
  }
  num_cat_per_f_pred_var <- purrr::map_int(f_pred_vars, length)
  num_f_pred_var_cats <- sum(num_cat_per_f_pred_var)
  
  # 0-or-1 encode each cat of f_pred_vars for every plate
  if (predict_f) {
    design_matrix_f <- df_plate %>%
      dplyr::select(tidyselect::all_of(f_pred_vars_names)) %>%
      dplyr::mutate(dplyr::across(tidyselect::everything(), as.factor)) %>%
      {stats::model.matrix(~ . - 1,
                    data = .,
                    contrasts.arg = lapply(.[, , drop = FALSE],
                                           stats::contrasts, contrasts = FALSE))}
  } else {
    design_matrix_f <- matrix(NA_real_, nrow = num_plate, ncol = 0)
  }
  
  # Count cats per f pred var. Ensure the col names of design_matrix_f are as 
  # expected.
  design_matrix_f_colnames_expected <-
    purrr::map(f_pred_vars_names, ~ paste0(.x, f_pred_vars[[.x]])) %>%
    unlist
  stopifnot(identical(sort(colnames(design_matrix_f)),
                      sort(design_matrix_f_colnames_expected)))
  if (!is.null(design_matrix_f_colnames_expected)) {
    design_matrix_f <- design_matrix_f[, design_matrix_f_colnames_expected]
  }
  stopifnot(identical(colnames(design_matrix_f),
                      design_matrix_f_colnames_expected))

  # Look-ups between int and string encodings
  df_f_pred_vars <- tibble::tibble(f_pred_var = f_pred_vars_names,
                           f_pred_var_int = seq(1, length.out = num_f_pred_vars))
  df_f_pred_vars_cats <- tibble::tibble(f_pred_var_cat = design_matrix_f_colnames_expected,
                                f_pred_var_cat_int = seq(1, length.out = num_f_pred_var_cats))

  
  # GET POSTERIOR INPUT INTO STAN FORMAT ----
  
  stan_input_posterior <- list(
    num_plate = num_plate,
    num_cal_tot = nrow(df_cal),
    num_sam_id = num_sam_id,
    num_sam_rep = num_sam_rep,
    which_plate_cal = df_cal$plate_int,
    which_plate_sam = df_sam$plate_int,
    which_id_sam = df_sam$id_sam_int,
    y_cal = df_cal$y,
    y_sam = df_sam$y,
    x_cal = df_cal$x,
    num_f_pred_vars = num_f_pred_vars,
    num_cat_per_f_pred_var = num_cat_per_f_pred_var %>% as.array(),
    num_p_pos_binary_pred_vars = length(p_pos_binary_pred_vars),
    design_matrix_f = design_matrix_f,
    design_matrix_p_pos_binary = p_pos_binary_design_matrix,
    sample_posterior_not_prior = 1L
  )
  
  for (param in x_mix_params) {
    stan_input_posterior[[paste0("design_matrix_", param)]] <-
      x_mix_design_matrices[[param]]
    stan_input_posterior[[paste0("num_", param, "_pred_vars")]] <- 
      x_mix_pred_vars_nums[[param]]
    stan_input_posterior[[paste0("num_cat_per_", param, "_pred_var")]] <- 
      x_mix_pred_vars_num_cats[[param]] %>% as.array()
  }
  
  # Hyperparameters of the priors
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
  stan_input_posterior$rho_prior_eta <- rho_prior_eta
  
  # GET PRIOR INPUT INTO STAN FORMAT ----
  
  # Define minimal data to sample the desired number of params from the prior 
  stan_input_prior <- stan_input_posterior
  stan_input_prior$sample_posterior_not_prior <- 0L
  stan_input_prior$num_plate <- 0L
  stan_input_prior$num_cal_tot <- 0L
  stan_input_prior$num_sam_id <- 0L
  stan_input_prior$num_sam_rep <- 0L
  stan_input_prior$which_plate_cal <- integer() %>% as.array()
  stan_input_prior$which_plate_sam <- integer() %>% as.array()
  stan_input_prior$which_id_sam <- integer() %>% as.array()
  stan_input_prior$y_cal <- numeric() %>% as.array()
  stan_input_prior$y_sam <- numeric() %>% as.array()
  stan_input_prior$x_cal <- numeric() %>% as.array()
  stan_input_prior$design_matrix_f <- stan_input_posterior$design_matrix_f[0, ]
  for (param in x_mix_params) {
    stan_input_prior[[paste0("design_matrix_", param)]] <-
      stan_input_posterior[[paste0("design_matrix_", param)]][0, ]
  }
  stan_input_prior$design_matrix_p_pos_binary <- 
    stan_input_posterior$design_matrix_p_pos_binary[0, ]
  
  # RETURN OUTPUT ----
  
  output <- list(
    df_sam = df_sam,
    df_cal = df_cal,
    df_plate = df_plate,
    stan_input_posterior = stan_input_posterior,
    stan_input_prior = stan_input_prior
  ) 
  
  output$data_descriptors <- list(
    num_plate = num_plate,
    num_sam_rep = num_sam_rep,
    num_sam_id = num_sam_id,
    # x mix vars
    x_mix_pred_vars = x_mix_pred_vars,
    x_mix_pred_vars_names = x_mix_pred_vars_names,
    x_mix_pred_vars_nums = x_mix_pred_vars_nums,
    x_mix_pred_vars_num_cats = x_mix_pred_vars_num_cats,
    x_mix_pred_vars_num_cats_tots = x_mix_pred_vars_num_cats_tots,
    predict_p_pos = predict_p_pos,
    predict_mu_pos = predict_mu_pos,
    predict_sd_pos = predict_sd_pos,
    predict_mu_neg = predict_mu_neg,
    predict_sd_neg = predict_sd_neg,
    x_mix_design_matrices = x_mix_design_matrices,
    lookup_pred_var_int = lookup_pred_var_int,
    lookup_pred_var_cat_int = lookup_pred_var_cat_int,
    # p_pos_binary vars
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    predict_p_pos_binary = predict_p_pos_binary,
    p_pos_binary_design_matrix = p_pos_binary_design_matrix,
    # f pred vars
    f_pred_vars_names = f_pred_vars_names,
    num_f_pred_vars = num_f_pred_vars,
    predict_f = predict_f,
    f_pred_vars = f_pred_vars,
    num_cat_per_f_pred_var = num_cat_per_f_pred_var,
    num_f_pred_var_cats = num_f_pred_var_cats,
    design_matrix_f = design_matrix_f,
    df_f_pred_vars = df_f_pred_vars,
    df_f_pred_vars_cats = df_f_pred_vars_cats
  )
  
  return(output)
  
}