#' Gets the true parameter values from simulated data into dataframe format
#'
#' @param param_true_values_list a list with the true parameter values, in the
#'   format output by [simulate_data()] (inside its list of outputs).
#' @param data_descriptors a list of things describing the dataset, of the
#'   format output by [prepare_data_for_stan()] (inside its list of outputs).
#'
#' @returns A dataframe with columns `param` and `value`, with one row per
#'   parameter.
#' @export
#'
#' @examples
wrangle_true_params <- function(param_true_values_list,
                                data_descriptors) {
  
  x_mix_params <- c("p_pos", "mu_neg", "mu_pos", "sd_neg", "sd_pos")
  
  p <- param_true_values_list
  d <- data_descriptors

df_true_pop_params <- tribble(
  ~param, ~value,
  "f[1]", p$f[1],
  "f[2]", p$f[2],
  "f[3]", p$f[3],
  "f[4]", p$f[4],
  "sigma_f_plate[1]", p$sigma_f_plate[1],
  "sigma_f_plate[2]", p$sigma_f_plate[2],
  "sigma_f_plate[3]", p$sigma_f_plate[3],
  "sigma_f_plate[4]", p$sigma_f_plate[4],
  "rho[1,2]", p$rho[1,2],
  "rho[1,3]", p$rho[1,3],
  "rho[1,4]", p$rho[1,4],
  "rho[2,3]", p$rho[2,3],
  "rho[2,4]", p$rho[2,4],
  "rho[3,4]", p$rho[3,4],
  "sd_neg", p$sd_neg,
  "sd_pos", p$sd_pos,
  "mu_neg", p$mu_neg,
  "mu_pos", p$mu_pos,
  "p_pos", p$p_pos,
  "p_blank", p$p_blank,
  "y_obs_sd_cal_min", p$y_obs_sd_cal_min,
  "y_obs_sd_cal_jump", p$y_obs_sd_cal_jump,
  "y_obs_sd_sam_min", p$y_obs_sd_sam_min,
  "y_obs_sd_sam_jump", p$y_obs_sd_sam_jump,
  "y_obs_sd_min_log_shift_sd", p$y_obs_sd_min_log_shift_sd,
  "y_obs_sd_jump_log_shift_sd", p$y_obs_sd_jump_log_shift_sd
) 

if (d$predict_f) {
  df_true_f_effects_by_pred_var <- p$f_effects_by_pred_var %>% 
    map(function(mat) {mat %>%
        as_tibble(.name_repair = "universal_quiet") %>%
        mutate(cat = rownames(mat))}) %>% 
    bind_rows(.id = "f_pred_var") %>%
    mutate(f_pred_var_cat = paste0(f_pred_var, cat)) %>% 
    inner_join(d$df_f_pred_vars_cats, by = "f_pred_var_cat")
  stopifnot(identical(sort(df_true_f_effects_by_pred_var$f_pred_var_cat),
                      sort(d$design_matrix_f_colnames_expected)))
  df_true_f_effects_by_pred_var <- df_true_f_effects_by_pred_var %>%
    pivot_longer(paste0("...", 1:4),
                 names_prefix = "...",
                 names_to = "which_f") %>%
    mutate(param = paste0("f_effects_by_pred_var_cat[", f_pred_var_cat_int,
                          ",", which_f, "]")) %>%
    select(param, value)
  df_true_pop_params <- df_true_pop_params %>%
    bind_rows(df_true_f_effects_by_pred_var,
              tibble(f_pred_var = names(p$sigma_f_pred_vars),
                     value = p$sigma_f_pred_vars) %>% 
                unnest_longer(value, indices_to = "which_f") %>%
                left_join(p$df_f_pred_vars, by = "f_pred_var") %>%
                mutate(param = paste0("sigma_f_pred_vars[", f_pred_var_int,
                                      ",", which_f, "]")) %>%
                select(param, value))
}

df_true_pop_params <- df_true_pop_params %>% bind_rows(
  map(x_mix_params, function(param) {
    if (d$x_mix_pred_vars_nums[[param]] == 0) {
      return(tibble(param = character(), value = numeric()))
    }
    df_true_effects_by_pred_var <- p$x_mix_effects[[param]] %>% 
      map(function(mat) {mat %>%
          as_tibble(.name_repair = "universal_quiet") %>%
          mutate(cat = names(mat))}) %>% 
      bind_rows(.id = "pred_var") %>%
      mutate(pred_var_cat = paste0(pred_var, cat)) %>% 
      inner_join(d$lookup_pred_var_cat_int[[param]], by = "pred_var_cat")
    stopifnot(identical(sort(df_true_effects_by_pred_var$pred_var_cat),
                        sort(colnames(d$x_mix_design_matrices[[param]]))))
    df_true_effects_by_pred_var <- df_true_effects_by_pred_var %>%
      mutate(param = paste0(param, "_effects_by_pred_var_cat[", 
                            pred_var_cat_int, "]")) %>%
      select(param, value)
    df_true_sigma_pred_vars <-
      tibble(pred_var = names(p$x_mix_pred_vars_sds[[param]]),
             value = p$x_mix_pred_vars_sds[[param]]) %>%
      inner_join(d$lookup_pred_var_int[[param]], by = "pred_var") %>%
      mutate(param = paste0("sigma_", param, "_pred_vars[", pred_var_int, "]")) %>%
      select(param, value)
    bind_rows(df_true_effects_by_pred_var, df_true_sigma_pred_vars)
  }) %>%
    bind_rows())

if (d$predict_p_pos_binary) {
  df_true_pop_params <- df_true_pop_params %>%
    bind_rows(tibble(
      param = paste0("p_pos_binary_effects[",
                     seq(from = 1,
                         length.out = length(d$p_pos_binary_pred_vars)), "]"),
      value = p$p_pos_binary_effects)) 
}

df_true_pop_params
}
