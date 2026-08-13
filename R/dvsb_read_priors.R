library(tidyverse)

#' Reads uniform prior distributions (specified with a min and max) from a csv
#'
#' @param path_to_csv path to a csv file with columns parameter, min and max,
#'   and one row for each of the (top-level, fixed-effect) parameters in our
#'   statistical model. Such a file is included in the github repository for
#'   this code, with an additional column explaining the parameter.
#'
#' @returns a list with elements rho_prior_eta (a single number controlling the
#'   prior for the eta parameter), df_priors_scalars and df_priors_vectors
#'   (dataframes containing the priors for scalar and vector quantities
#'   respectively).
#' @export
#'
read_priors <- function(path_to_csv) {
  
  stopifnot(file.exists(path_to_csv))
  
  param_names <- c(
    "mu_neg",
    "sd_neg",
    "sd_pos",
    "mu_pos",
    "p_pos",
    "p_blank",
    "y_obs_sd_cal_min",
    "y_obs_sd_cal_jump",
    "y_obs_sd_sam_min",
    "y_obs_sd_sam_jump",
    "sigma_p_pos_pred_vars",
    "sigma_mu_pos_pred_vars",
    "sigma_sd_pos_pred_vars",
    "sigma_mu_neg_pred_vars",
    "sigma_sd_neg_pred_vars",
    "p_pos_binary_effects",
    "rho_prior_eta",
    "f[1]",
    "f[2]",
    "f[3]",
    "f[4]",
    "sigma_f_plate[1]",
    "sigma_f_plate[2]",
    "sigma_f_plate[3]",
    "sigma_f_plate[4]",
    "sigma_f_pred_vars[1]",
    "sigma_f_pred_vars[2]",
    "sigma_f_pred_vars[3]",
    "sigma_f_pred_vars[4]"
  )
  non_neg_params <- c(
    "sd_neg", 
    "sd_pos", 
    "p_pos", 
    "p_blank",
    "y_obs_sd_cal_min",
    "y_obs_sd_cal_jump", 
    "y_obs_sd_sam_min",
    "y_obs_sd_sam_jump", 
    "sigma_p_pos_pred_vars",
    "sigma_mu_pos_pred_vars",
    "sigma_sd_pos_pred_vars",
    "sigma_mu_neg_pred_vars", 
    "sigma_sd_neg_pred_vars",
    "rho_prior_eta",
    "sigma_f_plate[1]",
    "sigma_f_plate[2]",
    "sigma_f_plate[3]",
    "sigma_f_plate[4]",
    "sigma_f_pred_vars[1]",
    "sigma_f_pred_vars[2]",
    "sigma_f_pred_vars[3]",
    "sigma_f_pred_vars[4]"
  )
  probability_params <- c(
    "p_pos", 
    "p_blank"
  )
  vector_params <- c(
    "f[1]",
    "f[2]",
    "f[3]",
    "f[4]",
    "sigma_f_plate[1]",
    "sigma_f_plate[2]",
    "sigma_f_plate[3]",
    "sigma_f_plate[4]",
    "sigma_f_pred_vars[1]",
    "sigma_f_pred_vars[2]",
    "sigma_f_pred_vars[3]",
    "sigma_f_pred_vars[4]"  
  )
  
  df_prior_limits <- readr::read_csv(path_to_csv, col_types = readr::cols_only(
    parameter = readr::col_character(),
    min = readr::col_double(),
    max = readr::col_double()
  )) %>%
    dplyr::rename(param = parameter)
  
  # Check for duplicated params
  problem_params <- df_prior_limits$param[duplicated(df_prior_limits$param)]
  if (length(problem_params)) {
    stop(paste(
      "The following parameters appeared in more than one row:",
      paste(problem_params, collapse = ", ")))
  }

  # Check for missing params
  if (! all(param_names %in% df_prior_limits$param)) {
      stop(paste(
        "The following expected parameters were missing:",
        paste(param_names[! param_names %in% df_prior_limits$param], 
              collapse = ", ")))
  }
  
  # Check for unexpected params
  if (! all(df_prior_limits$param %in% param_names)) {
    stop(paste(
      "The following parameters were unexpected:",
      paste(df_prior_limits$param[! df_prior_limits$param %in% param_names],
            collapse = ", ")))
  }
  
  # Check for missing mins
  problem_params <- df_prior_limits %>% 
    dplyr::filter(is.na(min)) %>%
    dplyr::pull(param)
  if (length(problem_params)) {
    stop(paste("The following parameters had missing (NA) min:",
               paste(problem_params, collapse = ", ")))
  }
  
  # Check for missing maxs
  problem_params <- df_prior_limits %>% 
    dplyr::filter(is.na(max)) %>%
    dplyr::pull(param)
  if (length(problem_params)) {
    stop(paste("The following parameters had missing (NA) max:",
               paste(problem_params, collapse = ", ")))
  }
  
  # Check max > min
  problem_params <- df_prior_limits %>% 
    dplyr::filter(max <= min,
           param != "rho_prior_eta") %>%
    dplyr::pull(param)
  if (length(problem_params)) {
    stop(paste("The following parameters had a max <= min:",
               paste(problem_params, collapse = ", ")))
  }
  
  # Check rho_prior_eta has min = max
  problem_params <- df_prior_limits %>% 
    dplyr::filter(param == "rho_prior_eta",
           max != min) %>%
    dplyr::pull(param)
  if (length(problem_params)) {
    stop(paste("For parameter rho_prior_eta only, min should equal max (see",
               "parameter explanation in the csv)"))
  }
  
  # Check mins and maxs are in allowable ranges
  problem_params <- df_prior_limits %>% 
    dplyr::filter(param %in% non_neg_params) %>%
    dplyr::filter(min < 0) %>%
    dplyr::pull(param)
  if (length(problem_params)) {
    stop(paste("The following parameters, which can never be negative, had a",
               "negative min value:", paste(problem_params, collapse = ", ")))
  }
  problem_params <- df_prior_limits %>% 
    dplyr::filter(param %in% probability_params) %>%
    dplyr::filter(max > 1) %>%
    dplyr::pull(param)
  if (length(problem_params)) {
    stop(paste("The following parameters, which must always be less than or",
               "equal to 1, had a max value greater than 1:", 
               paste(problem_params, collapse = ", ")))
  }
  mu_neg_min <- df_prior_limits$min[df_prior_limits$param == "mu_neg"]
  mu_pos_min <- df_prior_limits$min[df_prior_limits$param == "mu_pos"]
  mu_neg_max <- df_prior_limits$max[df_prior_limits$param == "mu_neg"]
  mu_pos_max <- df_prior_limits$max[df_prior_limits$param == "mu_pos"]
  if (mu_pos_min < mu_neg_min) stop("mu_pos min must be greater than mu_neg min")
  if (mu_pos_max < mu_neg_max) stop("mu_pos max must be greater than mu_neg max")
  
  df_prior_limits <- df_prior_limits %>%
    dplyr::rename(lower = min,
           upper = max)
  
  df_priors_scalars <- df_prior_limits %>%
    dplyr::filter(! param %in% c(vector_params, "rho_prior_eta"))
  df_priors_vectors <- df_prior_limits %>%
    dplyr::filter(param %in% vector_params)
  rho_prior_eta <- df_prior_limits$lower[df_prior_limits$param == "rho_prior_eta"]
  
  # Coerce separate rows for "f[1]", "f[2]" etc. into a single vector-valued
  # row for "f". It's very ugly.
  df_priors_vectors <- df_priors_vectors %>%
    tidyr::extract(param, 
                   into = c("param", "index"), 
                   regex = "(.*)\\[([0-9]+)\\]$") %>%
    tidyr::pivot_wider(names_from = index, values_from = c("upper", "lower")) 
  df_priors_vectors$lower <- list(numeric(4), numeric(4), numeric(4))
  df_priors_vectors$upper <- list(numeric(4), numeric(4), numeric(4))
  for (index in 1:4) {
    for (row in (1:nrow(df_priors_vectors))) {
      df_priors_vectors$lower[[row]][[index]] <- 
        df_priors_vectors[[paste0("lower_", index)]][[row]]
      df_priors_vectors$upper[[row]][[index]] <- 
        df_priors_vectors[[paste0("upper_", index)]][[row]]
    }
  }
  df_priors_vectors <- df_priors_vectors %>%
    dplyr::select(param, lower, upper)
  
  list(df_priors_scalars = df_priors_scalars,
       df_priors_vectors = df_priors_vectors,
       rho_prior_eta = rho_prior_eta)
}
