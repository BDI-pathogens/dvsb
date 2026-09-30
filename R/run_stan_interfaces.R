#' Run one of rstan, cmdstanr or cmdstan on a file of Stan code
#'
#' @param input_to_stan a list containing all the input the Stan code expects.
#' @param path_to_stan_code the path to the file containing the Stan code.
#' @param interface one of `"rstan"`, `"cmdstanr"` or `"cmdstan"`.
#' @param iter_warmup a positive integer: the number of warmup iterations per
#'   chain (during which the sampling algorithm adapts; these are excluded from
#'   the output).
#' @param iter_sampling a positive integer: the number of sampling iterations
#'   per chain (which are included in the output).
#' @param chains a positive integer: the number of chains used for sampling.
#' @param cores a positive integer: the number cores used in parallel for
#'   computation.
#' @param params_to_ignore a character vector naming parameters to be excluded
#'   from output (if possible; interface dependent).
#' @param downsampling_factor a positive integer: the factor by which to
#'   downsample the posterior. e.g. if a value of 2 is specified, we keep 1 in
#'   every 2 samples. The default of 1 means we keep all samples. TODO:
#'   currently only implemented for cmdstan.
#' @param cmdstan_path_to_installation the path to where cmdstan is installed on
#'   your system; you need to specify this if `interface="cmdstan"`, but not
#'   otherwise. Inside this directory there should be an executable file named
#'   `make` (which we use to compile Stan code).
#' @param cmdstan_path_to_output the path to where we will write output files
#'   from cmdstan; you need to specify this if `interface="cmdstan"`, but not
#'   otherwise. Several files will be created with things appended to this path:
#'   _chain1.csv, _chain2.csv etc.
#' @param cmdstan_path_to_json the path to where we will write a temporary json
#'   file to hold the input for cmdstan; you need to specify this if
#'   `interface="cmdstan"`, but not otherwise.
#' @param cmdstan_overwrite_json a single logical value: should we overwrite a
#'   file at `cmdstan_path_to_json` if it exists already?
#' @param cmdstan_path_to_compiled_model the path where we will create the
#'   compiled version of the Stan code. Some value (such as the default) is
#'   needed if `interface="cmdstan"`, but not otherwise.
#' @param cmdstan_read_output_into_df a single logical value: should we read
#'   cmdstan output files into a dataframe that is returned by this function? If
#'   a value `FALSE` is specified, this function returns a value `NULL`.
#' @param ... additional arguments will be passed to [rstan::sampling()] (if
#'   `interface="rstan"`) or to the `$sample()` method of the
#'   [cmdstanr::CmdStanModel()] object (if `interface="cmdstanr"`).
#'
#' @returns a dataframe with one row per sample from the posterior and one
#'   column per parameter (unless `interface` is set to `cmdstan` and
#'   `cmdstan_read_output_into_df` is set to `FALSE`). Parameters you'll find in
#'   here include all of the parameters in the csv file that you previously read
#'   into a dataframe using [read_priors()] (see that csv or dataframe for
#'   descriptions of these parameters), plus:
#'    * `rho[..., ...]` is the 4x4 dimensionless correlation matrix between the values of the f vector for different plates, with the square brackets containing two indices for the matrix element.
#'    * `xlog_sam[...]` is the log_e antibody level for samples, with the square brackets containing an integer index for which sample.
#'    * `p_pos_effect_bool{name}`, with `{name}` being one of the names inside `p_pos_binary_pred_vars` given as input to [prepare_data_for_stan()], is the effect on the overall `p_pos` parameter due to the boolean variable `{name}` taking the value `TRUE` instead of `FALSE`.
#'    * `f_per_plate[..., ...]` is the f vector that relates x and y for each plate. The first integer indexes the plate and the second indexes one of the four elements of the vector.
#'    * `sigma_p_pos_pred_vars_{name1}`, with `{name1}` being one of the names inside the `p_pos` element of the `x_mix_pred_vars_names` list given as input to [prepare_data_for_stan()] (i.e. the name of one categorical variable used in the regression model for `p_pos`) is the scale of variability in `p_pos` between different categories of this variable.
#'    * `sigma_mu_pos_pred_vars_{name1}`, `sigma_sd_pos_pred_vars_{name1}`, `sigma_mu_neg_pred_vars_{name1}`, `sigma_sd_neg_pred_vars_{name1}` are all defined analogously to `sigma_p_pos_pred_vars_{name1}` but for the other four parameters of the x mixture distribution: `mu_pos`, `sd_pos`, `mu_neg` and `sd_neg` respectively.
#'    * `p_pos_effect_{name1}{name2}`, with `{name1}` matching `{name1}` in `sigma_p_pos_pred_vars_{name1}` and `{name2}` being one of the categories of `{name1}`, is the effect of this category on `p_pos`: adding this parameter to the overall `p_pos` parameter gives the value of `p_pos` for this category (a logit link function is used, and if more than one categorical variable was used for the regression model, all of the associated `{name1}` parameters must be summed over to get to a single subpopulation, e.g. this could look like `logistic(logit(p_pos) + p_pos_effect_Age18-29 + p_pos_effect_JobFarmer)`);
#'    * `mu_pos_effect_{name1}{name2}`, `sd_pos_effect_{name1}{name2}`, `mu_neg_effect_{name1}{name2}`, `sd_neg_effect_{name1}{name2}`, these are all defined analogously to `p_pos_effect_{name1}{name2}` but for the other four parameters of the x mixture distribution: `mu_pos`, `sd_pos`, `mu_neg` and `sd_neg` respectively (with a log link function for the `sd` parameters and no link function for the `mu` parameters).
#'    * `y_cal_sim[...]` is simulated new y values for calibrators (with the square brackets containing an integer index for which calibrator), providing a posterior retrodictive check for the calibrator y values actually observed.
#'    * `y_sam_sim_conditional[...]` is simulated new y values for samples, conditioning on the y values actually observed for that specific sample, i.e. imagining obtaining new measurements from the same individuals.
#'    * `y_sam_sim_unconditional[...]` is simulated new y values for samples, not conditioning on the y values actually observed for that specific sample, i.e. imagining resampling a new individual from the same subpopulation with all the same covariates. This provides a more stringent posterior retrodictive check for subpopulation distributions of y values.
#'    * `xlog_sam_sim_unconditional[...]` is simulated new x values, one per individual (indexed by the integer in the square brackets), not conditioning on the y values actually observed for that specific individual, i.e. imagining resampling a new individual from the same subpopulation with all the same covariates. This allows checking of the consistency between the distribution of estimated x values and the estimated distribution of x values.
#'    * `pos_sam_sim_unconditional[...]` is simulated new serostatuses, one per individual (indexed by the integer in the square brackets), not conditioning on the y values actually observed for that specific sample, i.e. imagining resampling a new individual from the same subpopulation with all the same covariates. Grouping these values together based on any covariate(s) desired, counting the fraction of them that are positive, and then taking the distribution over posterior samples provides a convenient posterior predictive distribution for stratified seroprevalence in a new dataset of individuals with the same distribution of covariates as the one actually sampled. It takes into account the joint effect of all modelled covariates and their correlations, and the sampling uncertainty associated with the number of individuals possessing each combination of covariates.
#'    * `p_sam_is_pos[...]` is the probability that a given sample (indexed by the integer in square brackets) is seropositive, conditional on its observed y values. The best way to understand this is as the fraction of a large hypothetical population of individuals with identical y values that would be positive; it takes continuous values between 0 and 1, and has a posterior distribution capturing its uncertainty.
#' @export
#'
run_stan_interfaces <- function(input_to_stan,
                                path_to_stan_code,
                                interface = c("rstan", "cmdstanr", "cmdstan"),
                                iter_warmup = 250,
                                iter_sampling = 250,
                                chains = 4,
                                cores = parallel::detectCores(),
                                params_to_ignore = c(
                                  "lp__",
                                  "rho[1,1]",
                                  "rho[2,2]",
                                  "rho[3,3]",
                                  "rho[4,4]",
                                  "rho[2,1]",
                                  "rho[3,1]",
                                  "rho[4,1]",
                                  "rho[3,2]",
                                  "rho[4,2]",
                                  "rho[4,3]",
                                  "rho.1.1",
                                  "rho.2.2",
                                  "rho.3.3",
                                  "rho.4.4",
                                  "rho.2.1",
                                  "rho.3.1",
                                  "rho.4.1",
                                  "rho.3.2",
                                  "rho.4.2",
                                  "rho.4.3",
                                  "f_plate_effects_unscaled",
                                  "f_effects_by_pred_var_cat_unscaled",
                                  "p_pos_effects_by_pred_var_cat_unscaled",
                                  "mu_pos_effects_by_pred_var_cat_unscaled",
                                  "sd_pos_effects_by_pred_var_cat_unscaled",
                                  "mu_neg_effects_by_pred_var_cat_unscaled",
                                  "sd_neg_effects_by_pred_var_cat_unscaled",
                                  "y_cal_mean_per_obs",
                                  "y_cal_mean_per_obs",
                                  "y_sam_mean_per_obs",
                                  "y_obs_sd_cal",
                                  "y_obs_sd_sam",
                                  "exp_f_1_mult_f_4_per_plate",
                                  "p_pos_log_per_sam_id",
                                  "p_neg_log_per_sam_id",
                                  "p_pos_per_sam_id",
                                  "mu_pos_per_sam_id",
                                  "mu_neg_per_sam_id",
                                  "sd_pos_per_sam_id",
                                  "sd_neg_per_sam_id",
                                  "f_3_min_f_2_per_plate",
                                  "p_blank_log",
                                  "p_blank_log1m",
                                  "y_obs_sd_min_log_shift_unscaled",
                                  "y_obs_sd_jump_log_shift_unscaled",
                                  "y_obs_sd_sam_min_per_plate",
                                  "y_obs_sd_cal_min_per_plate",
                                  "y_obs_sd_sam_jump_per_plate",
                                  "y_obs_sd_cal_jump_per_plate"
                                ),
                                downsampling_factor = 1L,
                                cmdstan_path_to_installation = NA,
                                cmdstan_path_to_json = NA, 
                                cmdstan_overwrite_json = FALSE, 
                                cmdstan_path_to_output = NA,
                                cmdstan_path_to_compiled_model =
                                  stringr::str_remove(path_to_stan_code, ".stan$"),
                                cmdstan_read_output_into_df = TRUE,
                                ...){
  
  # Check args
  stopifnot(is.list(input_to_stan))
  stopifnot(is.character(path_to_stan_code))
  stopifnot(length(path_to_stan_code) == 1)
  stopifnot(file.exists(path_to_stan_code))
  stopifnot(endsWith(path_to_stan_code, ".stan"))
  mastiff::check_numeric(iter_warmup, lower = 1)
  mastiff::check_numeric(iter_sampling, lower = 1)
  mastiff::check_numeric(chains, lower = 1)
  mastiff::check_numeric(cores, lower = 1)
  stopifnot(is.character(params_to_ignore))
  stopifnot(is.character(interface))
  interface <- match.arg(interface)
  if (interface == "cmdstan") {
    if (identical(cmdstan_path_to_json, NA)) stop(paste(
      "If interface is set to cmdstan, the cmdstan_path_to_json option",
      "must be used"
    ))
    if (identical(cmdstan_path_to_output, NA)) stop(paste(
      "If interface is set to cmdstan, the cmdstan_path_to_output option",
      "must be used"
    ))
    if (identical(cmdstan_path_to_installation, NA)) stop(paste(
      "If interface is set to cmdstan, the cmdstan_path_to_installation option",
      "must be used"
    ))
    stopifnot(is.character(cmdstan_path_to_json))
    stopifnot(length(cmdstan_path_to_json) == 1)
    stopifnot(is.character(cmdstan_path_to_output))
    stopifnot(length(cmdstan_path_to_output) == 1)
    stopifnot(is.character(cmdstan_path_to_installation))
    stopifnot(length(cmdstan_path_to_installation) == 1)
    stopifnot(is.character(cmdstan_path_to_compiled_model))
    stopifnot(length(cmdstan_path_to_compiled_model) == 1)
    stopifnot(dir.exists(cmdstan_path_to_installation))
    cmdstan_path_to_make <- file.path(cmdstan_path_to_installation, "make")
    if (! file.exists(cmdstan_path_to_make)) stop(paste(
      "Could not find a make file inside", cmdstan_path_to_installation))
    if (! cmdstan_overwrite_json && file.exists(cmdstan_path_to_json)) stop(paste(
      cmdstan_path_to_json, 
      "exists already; please move/rename/delete to prevent overwriting,",
      "or run again with cmdstan_overwrite_json set to TRUE"
    ))
  }  
  mastiff::check_logical(cmdstan_read_output_into_df)
  mastiff::check_numeric(downsampling_factor, lower = 1)
  
  # Compile
  if (interface == "rstan") {
    model_compiled <- rstan::stan_model(path_to_stan_code, auto_write = TRUE)
  } else if (interface == "cmdstanr") {
    model_compiled <- cmdstanr::cmdstan_model(path_to_stan_code)
  } else {
    system(paste("cd", cmdstan_path_to_installation, "&& make STAN_THREADS=true", cmdstan_path_to_compiled_model))
  }
  
  start_time <- Sys.time()
  cat("Started running Stan at")
  print(start_time)
  
  if (interface == "rstan") {
    df_samples <- rstan::sampling(model_compiled,
                                  data = input_to_stan,
                                  iter = iter_warmup + iter_sampling,
                                  warmup = iter_warmup,
                                  chains = chains,
                                  cores = cores,
                                  pars = params_to_ignore,
                                  include = FALSE,
                                  ...) %>%
      as.data.frame()
    data.table::setDT(df_samples)
    
  } else if (interface == "cmdstanr") {
    samples <- model_compiled$sample(
      data = input_to_stan,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      parallel_chains = cores,
      ...
    )
    df_samples <- samples$draws(format = "draws_df")
    data.table::setDT(df_samples)

  } else {
    cmdstanr::write_stan_json(input_to_stan, file = cmdstan_path_to_json)
    files_out_stan <- paste0(cmdstan_path_to_output, "_chain", 1:chains, ".csv")
    files_out_stan_profile <- paste0(cmdstan_path_to_output, "_profiles.csv")
    command <- paste0(cmdstan_path_to_compiled_model,
                      " method=sample",
                      " num_chains=", chains,
                      " num_warmup=", iter_warmup,
                      " num_samples=", iter_sampling,
                      " num_threads=", cores,
                      " data file=", cmdstan_path_to_json, 
                      " output file=", paste(files_out_stan, collapse = ","),
                      " profile_file=", files_out_stan_profile)
    print("About to run this command:")
    print(command)
    system(command)
    if (! all(file.exists(files_out_stan))) stop(paste(
      "Internal error: we expected to create all of the following files, but at",
      "least one does not exist:", paste(files_out_stan, collapse = " ")))
    cat(paste("Created the following output files:", 
              paste(files_out_stan, collapse = " "), "\n"))
    if (cmdstan_read_output_into_df) {
      df_samples <- read_cmdstan_out_files(
        file_paths = files_out_stan, 
        params_to_ignore = params_to_ignore,
        downsampling_factor = downsampling_factor)
    } 
  } 
  
  end_time <- Sys.time()
  cat("Finished running Stan at")
  print(end_time)
  print(end_time - start_time)
  
  # There is no dataframe to return in this case:
  if (interface == "cmdstan" && ! cmdstan_read_output_into_df) return(NULL)
  
  # Exclude cols if desired
  keep_col <- rep(TRUE, ncol(df_samples))
  for (param in params_to_ignore) {
    keep_based_on_this_param <- 
      colnames(df_samples) != param &
      ! startsWith(colnames(df_samples), paste0(param, ".")) &
      ! startsWith(colnames(df_samples), paste0(param, "["))
    keep_col <- keep_col & keep_based_on_this_param
  }
  df_samples <- df_samples[, ..keep_col]
  
  # TODO:
  #mastiff::rename_params_cmdstanfile_to_rstan() if cmdstan(r)
  #data.table::setnames(df_samples, function(names) {
  #rename_params_from_stan(names,
  #                        data_descriptors = input_to_stan$data_descriptors)})
  # Add a note to for the user to do that themself if ! cmdstan_read_output_into_df
  
  df_samples
  
}