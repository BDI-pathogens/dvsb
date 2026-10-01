#' Run the dvsb statistical model
#' 
#' Run the Stan code implementing the dvsb statistical model, for Bayesian
#' estimation of the parameters that generated a given dataset.
#' 
#' The first time you call `run_stan()` (after each installation or update of
#' dvsb and for a given combination of `model` and `interface` arguments) the
#' Stan code will get compiled before the model is run. The second and
#' subsequent times, compilation will normally be skipped due to the previously
#' compiled code still being available.
#' 
#' If `distribution` is set to `prior`, we reduce the number of samples to zero
#' but keep the same 'shape' of the dataset (in terms of the definition of
#' regression models with covariates). This lets us sample the prior for
#' population-level parameters, which are the main thing of interest, forgetting
#' about the individual-level parameters (x values), gaining a lot of speed for
#' big datasets. Note that the sampling algorithm used by Stan (NUTS, a type of
#' Hamiltonian Monte Carlo) is far from the most efficient way of sampling from
#' priors: more efficient would be to sequentially randomly draw a new set of
#' parameters independent of the previous draw, not to use fancy Hamiltonian
#' dynamical equations to make successive draws dependent. Nevertheless we
#' re-use the same Stan code to sample the prior so that we and you can ensure
#' that the prior specified in Stan really is what was intended, reducing the
#' potential for silent bugs when sampling from the posterior. The downside is
#' that Stan can have difficulty exploring the geometry of the prior
#' distribution, just like for the posterior distribution, so be careful to run
#' enough iterations to get convergence here too.
#'
#' @param data_wrangled a list containing everything about the data we need, as
#'   output by [wrangle_data()]
#' @param distribution one of "posterior" or "prior" - which distribution should
#'   we sample from? (See the details section of this help.)
#' @param interface one of `"rstan"`, `"cmdstanr"` or `"cmdstan"`.
#' @param model one of `"dvsb"`, `"dvsb_accidental_blanks"` or
#'   `"dvsb_no_mu_x_predictors"`. `dvsb` is the main dvsb model;
#'   `dvsb_accidental_blanks` is the extended model that allows for the
#'   possibility that each sample well accidentally contains nothing at all
#'   (i.e. is blank); `dvsb_no_mu_x_predictors` is the restricted model in which
#'   regression models are not specified for the means of the two antibody
#'   distributions (one for seropositives and one for seronegatives).
#' @param params_to_ignore a character vector naming parameters to be excluded
#'   from output. The default value is those parameters, including derived
#'   quantities that are calculated as part of the statistical model, that we
#'   think will typically not be of interest.
#'
#' @returns a dataframe with one row per sample from the posterior and one
#'   column per parameter (unless `interface` is set to `cmdstan` and
#'   `cmdstan_read_output_into_df` is set to `FALSE`). Parameters you'll find in
#'   here include all of the parameters in the csv file that you previously read
#'   into a dataframe using [get_priors()] (see that csv or dataframe for
#'   descriptions of these parameters), plus:
#'    * `rho[..., ...]` is the 4x4 dimensionless correlation matrix between the values of the f vector for different plates, with the square brackets containing two indices for the matrix element.
#'    * `xlog_sam[...]` is the log_e antibody level for samples, with the square brackets containing an integer index for which sample.
#'    * `p_pos_effect_bool{name}`, with `{name}` being one of the names inside `p_pos_binary_pred_vars` given as input to [wrangle_data()], is the effect on the overall `p_pos` parameter due to the boolean variable `{name}` taking the value `TRUE` instead of `FALSE`.
#'    * `f_per_plate[..., ...]` is the f vector that relates x and y for each plate. The first integer indexes the plate and the second indexes one of the four elements of the vector.
#'    * `sigma_p_pos_pred_vars_{name1}`, with `{name1}` being one of the names inside the `p_pos` element of the `x_mix_pred_vars_names` list given as input to [wrangle_data()] (i.e. the name of one categorical variable used in the regression model for `p_pos`) is the scale of variability in `p_pos` between different categories of this variable.
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
run_stan <- function(data_wrangled,
                     distribution = c("posterior", "prior"),
                     model = c("dvsb", "dvsb_accidental_blanks", "dvsb_no_mu_x_predictors"),
                     interface = c("rstan", "cmdstanr", "cmdstan"),
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
                     ...){
  
  # Check args
  stopifnot(is.list(data_wrangled))
  stopifnot("stan_input_posterior" %in% names(data_wrangled))
  stopifnot("stan_input_prior" %in% names(data_wrangled))
  stopifnot("data_descriptors" %in% names(data_wrangled))
  stopifnot(is.list(data_wrangled[["stan_input_posterior"]]))
  stopifnot(is.list(data_wrangled[["stan_input_prior"]]))
  stopifnot(is.list(data_wrangled[["data_descriptors"]]))
  stopifnot(is.character(params_to_ignore))
  stopifnot(is.character(interface))
  interface <- match.arg(interface)
  stopifnot(is.character(model))
  model <- match.arg(model)
  stopifnot(is.character(distribution))
  distribution <- match.arg(distribution)
  
  path_to_stan_code <- get_stan_file_path(model)
  input_to_stan <- data_wrangled[[paste0("stan_input_", distribution)]]
  result <- mastiff::run_stan_interfaces(path_to_stan_code = path_to_stan_code,
                                         input_to_stan = input_to_stan, 
                                         interface = interface, 
                                         params_to_ignore = params_to_ignore, 
                                         ...)
  
  if (! is.null(result)) {
    data.table::setnames(result, function(names) {
    rename_params_from_stan(
      names, data_descriptors = data_wrangled[["data_descriptors"]])})
  }
  
  result
  
}

get_stan_file_path <- function(model) {
  stan_subdir <- system.file("stan", package = "dvsb")
  stan_path <- file.path(stan_subdir, paste0(model, ".stan"))
  if (! file.exists(stan_path)) {
    stop(paste("No file exists at the expected path of", stan_path))
  }
  stan_path
}
