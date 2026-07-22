#' Run one of rstan, cmdstanr or cmdstan on a file of Stan code
#'
#' @param input_to_stan a list containing all the input the Stan code expects.
#' @param path_to_stan_code the path to the file containing the Stan code.
#' @param interface one of "rstan", "cmdstanr" or "cmdstan".
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
#'   compiled version of the stan code. Some value (such as the default) is
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
#'   `cmdstan_read_output_into_df` is set to `FALSE`).
#' @export
#'
#' @examples
run_stan_interfaces <- function(input_to_stan,
                                path_to_stan_code,
                                interface = c("rstan", "cmdstanr", "cmdstan"),
                                iter_warmup = 250,
                                iter_sampling = 250,
                                chains = 4,
                                cores = parallel::detectCores(),
                                params_to_ignore = character(),
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
  mastiff::check_numeric(iterations, lower = 1)
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
      "exists already; please move/rename/delete to prevent overwriting"
    ))
  }  
  mastiff::check_logical(cmdstan_read_output_into_df)
  mastiff::check_numeric(downsampling_factor, lower = 1)
  
  # Compile
  if (interface == "rstan") {
    #rstan::rstan_options(auto_write = TRUE) # TODO?
    model_compiled <- rstan::stan_model(path_to_stan_code)
  } else if (interface == "cmdstanr") {
    model_compiled <- cmdstanr::cmdstan_model(path_to_stan_code)
  } else {
    system(paste("cd", dir_stan, "&& make STAN_THREADS=true", cmdstan_path_to_compiled_model))
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
    keep_col <- rep(TRUE, ncol(df_samples))
    for (param in params_to_ignore) {
      keep_based_on_this_param <- 
        colnames(df_samples) != param &
        ! startsWith(colnames(df_samples), paste0(param, ".")) &
        ! startsWith(colnames(df_samples), paste0(param, "["))
      keep_col <- keep_col & keep_based_on_this_param
    }
    df_samples <- df_samples[, ..keep_col]

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
  
  if (interface == "cmdstan" && ! cmdstan_read_output_into_df) return(NULL)
  
  df_samples
  
}