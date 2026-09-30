#' Change parameter names from Stan format to be more human readable
#'
#' @param original_names a character vector of parameter names as they are
#'   output by Stan. Tensor parameters are expected to be named as they are in
#'   rstan output, e.g. `my_matrix[1,2]`, not as they are named in cmdstan
#'   output, e.g. `my_matrix.1.2`. The latter can be converted to the former
#'   using [mastiff::rename_params_cmdstanfile_to_rstan()].
#' @param data_descriptors a list of things describing the dataset, of the
#'   format output by [prepare_data_for_stan()] (inside its list of outputs).
#'
#' @returns a character vector containing the renamed parameters, of the same
#'   length and in the same order as `original_names`.
#' @export
#'
rename_params_from_stan <- function(original_names, data_descriptors) {
  
  stopifnot(is.character(original_names))
  stopifnot(is.list(data_descriptors))
  descriptors_expected <- c(
    "predict_f", "df_f_pred_vars", "df_f_pred_vars_cats", 
    "x_mix_pred_vars_nums", "lookup_pred_var_int", "lookup_pred_var_cat_int",
    "predict_p_pos_binary", "p_pos_binary_pred_vars")
  for (descriptor in descriptors_expected) {
    if (! descriptor %in% names(data_descriptors)) stop(paste(
      "data_descriptors should have an element named", descriptor
    ))
  }
  
  d <- data_descriptors

  x_mix_params <- c("p_pos", "mu_neg", "mu_pos", "sd_neg", "sd_pos")
  
  df_param_names <- tibble::tibble(orig = original_names,
                           new = stringr::str_replace(orig,
                                             "sigma_f_plate\\[([0-9]+)\\]",
                                             "sigma_f[\\1]_plate"))
  if (d$predict_f) {
    df_param_names <- df_param_names %>%
      tidyr::extract(orig, 
                     into = c("f_pred_var_int", "which_f_foo"), 
                     regex = "sigma_f_pred_vars\\[([0-9]+),([0-9]+)\\]",
                     remove = FALSE) %>%
      tidyr::extract(orig, 
                     into = c("f_pred_var_cat_int", "which_f_spam"), 
                     regex = "f_effects_by_pred_var_cat\\[([0-9]+),([0-9]+)\\]",
                     remove = FALSE) %>%
      dplyr::mutate(f_pred_var_int = as.integer(f_pred_var_int),
             f_pred_var_cat_int = as.integer(f_pred_var_cat_int)) %>%
      dplyr::left_join(d$df_f_pred_vars, by = "f_pred_var_int") %>%
      dplyr::left_join(d$df_f_pred_vars_cats, by = "f_pred_var_cat_int") %>% 
      dplyr::mutate(new = dplyr::case_when(
        !is.na(f_pred_var_int) ~ paste0("sigma_f[", which_f_foo, "]_", f_pred_var),
        !is.na(f_pred_var_cat_int) ~ paste0("f_effect[", which_f_spam, "]_", f_pred_var_cat),
        TRUE ~ new
      ))
  }
  
  for (param in x_mix_params) {
    if (d$x_mix_pred_vars_nums[[param]] == 0) next
    df_param_names <- df_param_names %>%
      tidyr::extract(orig, 
                     into = "pred_var_int", 
                     regex = paste0("sigma_", param, "_pred_vars\\[([0-9]+)\\]"),
                     remove = FALSE) %>%
      tidyr::extract(orig, 
                     into = "pred_var_cat_int", 
                     regex = paste0(param, "_effects_by_pred_var_cat\\[([0-9]+)\\]"),
                     remove = FALSE) %>%
      dplyr::mutate(pred_var_int = as.integer(pred_var_int),
             pred_var_cat_int = as.integer(pred_var_cat_int)) %>%
      dplyr::left_join(d$lookup_pred_var_int[[param]], by = "pred_var_int") %>%
      dplyr::left_join(d$lookup_pred_var_cat_int[[param]], by = "pred_var_cat_int") %>% 
      dplyr::mutate(new = dplyr::case_when(
        !is.na(pred_var_int) ~ paste0("sigma_", param, "_pred_vars_", pred_var),
        !is.na(pred_var_cat_int) ~ paste0(param, "_effect_", pred_var_cat),
        TRUE ~ new
      )) %>%
      dplyr::select(orig, new)
  }
  
  if (d$predict_p_pos_binary) {
    df_param_names <- df_param_names %>%
      tidyr::extract(orig, 
                     into = "pred_var_int", 
                     regex = paste0("p_pos_binary_effects\\[([0-9]+)\\]"),
                     remove = FALSE) %>%
      dplyr::mutate(pred_var_int = as.integer(pred_var_int)) %>%
      dplyr::left_join(tibble::tibble(pred_var_int = seq(from = 1, length.out = length(d$p_pos_binary_pred_vars)),
                       pred_var = d$p_pos_binary_pred_vars),
                by = "pred_var_int") %>% 
      dplyr::mutate(new = dplyr::case_when(
        !is.na(pred_var_int) ~ paste0("p_pos_effect_", pred_var),
        TRUE ~ new
      )) %>%
      dplyr::select(orig, new)
  }
  
  df_param_names$new
}