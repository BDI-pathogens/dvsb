library(data.table)
library(tidyverse)
theme_set(theme_classic())

data_was_simulated <- TRUE
read_posterior_from_file <- FALSE
#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run-202510141258-*-97b346.csv") # first run with Anton's code debugged
#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2025-10-27-18h07m40_chain*.csv") # v19 on data 2025-10-27
#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2025-10-29-16h26m52_chain*.csv") # v19 on data 2025-10-27 with 15 of my dodgy plates excluded and longer chains
#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2025-11-03-22h29m50_chain*.csv") # v19 on all data
#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2025-11-13-21h38m27_chain*.csv") # v19 on all data with 1500 iter, with sex as p_pos predictor, slow but OK convergence 
files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2025-11-25-18h44m40_chain*.csv") # v20 with age as a predictor for all 5 params, with cluster and site and job for p_pos

#files_out_stan <- Sys.glob("/Users/cwymant/enable/samples_full_run_2025-12-03-10h28m49_chain*.csv") # first restriction to baseline only in some time(!), v20 with age as a predictor for all 5 params, with site and job for p_pos

# INPUT ABOUT STAN ----

file_input_stan <- "~/repos/dvsb/xsectional.stan"
dir_stan <- "~/.cmdstan/cmdstan-2.37.0/"
num_mc_chains <- 5
num_mc_iterations_posterior <- 500 # per chain, half of them warmup
num_mc_iterations_prior <- 10000
# one of: "rstan", "cmdstanr", "cmdstan". cmdstan uses cmdstanr for prior sampling.
stan_method <- "cmdstan" 
file_stan_temp <- "/Users/cwymant/foo.json" # for writing the data for cmdstan
file_out_stan_basename <- "/Users/cwymant/enable/samples_full_run_"

# Upper and lower bounds for priors
if (data_was_simulated) {
  
  df_priors_scalars <- tribble(
    ~param, ~lower, ~upper,
    "mu_neg", -5, -1,
    "mu_pos", -1.5, 3,
    "sd_neg", 0, sd_neg * 2,
    "sd_pos", 0, sd_pos * 2,
    "p_pos", 0, 1,
    "p_blank", 0, 2 * p_blank,
    "y_obs_sd_cal_min",  0, 2 * y_obs_sd_cal_min,
    "y_obs_sd_cal_jump", 0, 2 * y_obs_sd_cal_jump,
    "y_obs_sd_sam_min",  0, 2 * y_obs_sd_sam_min,
    "y_obs_sd_sam_jump", 0, 2 * y_obs_sd_sam_jump,
    "sigma_p_pos_pred_vars",  0, 10,
    "sigma_mu_pos_pred_vars", 0, 10,
    "sigma_sd_pos_pred_vars", 0, 10,
    "sigma_mu_neg_pred_vars", 0, 10,
    "sigma_sd_neg_pred_vars", 0, 10,
    "sigma_p_pos_binary_pred_vars", 0, 2
  )
  df_priors_vectors <- tribble(
    ~param, ~lower, ~upper,
    "sigma_f_plate", sigma_f_plate * 0, sigma_f_plate * 4,
    "f", f-0.5, f+0.5,
    "sigma_f_pred_vars", rep(0, 4), rep(0, 4)) # not needed if ! predict_f
  if (predict_f) {
    df_priors_vectors$upper[[3]] <- 2 * do.call(pmax, sigma_f_pred_vars)
  }
  
} else {
  df_priors_scalars <- tribble(
    ~param, ~lower, ~upper,
    "mu_neg", -5, -1,
    "sd_neg", 0.5, 2,
    "sd_pos", 0, 2.5, 
    "mu_pos", -1.5, 3,
    "p_pos", 0, 1,
    "p_blank", 0, 0.05,
    "y_obs_sd_cal_min", 0, 0.015,
    "y_obs_sd_cal_jump", 0.1, 1,
    "y_obs_sd_sam_min", 0, 0.03,
    "y_obs_sd_sam_jump", 0.1, 1,
    "sigma_p_pos_pred_vars", 0, 3,
    "sigma_mu_pos_pred_vars", 0, 2,
    "sigma_sd_pos_pred_vars", 0, 1,
    "sigma_mu_neg_pred_vars", 0, 2,
    "sigma_sd_neg_pred_vars", 0, 1,
    "sigma_p_pos_binary_pred_vars", 0, 2
  )
  df_priors_vectors <- tribble(
    ~param, ~lower, ~upper,
    "sigma_f_plate", c(0, 0, 0, 0), c(0.25, 0.025, 2.5, 2),
    "f", c(0.8, -0.05, 3, 2.2), c(1.1, 0.05, 5.5, 3.8),
    "sigma_f_pred_vars", c(0, 0, 0, 0), c(0.6, 0.1, 4, 3)
  )
  
}

# The eta parameter of the LKJ prior for rho
rho_prior_eta <- 1

# SYNCHRONISE SIMULATED AND REAL DATA PREVIOUS STEPS ----

if (data_was_simulated) {
  df_plate$plate_int <- df_plate$plate
  df_cal$plate_int <- df_cal$plate
} else {
  library(mvtnorm)
  PL4 <- function(xlog, f_1, f_2, f_3, f_4) {
    f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (xlog - f_4)))
  }
}

# FINISH PREPARING FOR STAN ----

predict_p_pos  <- x_mix_pred_vars_nums[["p_pos"]]  > 0L
predict_mu_pos <- x_mix_pred_vars_nums[["mu_pos"]] > 0L
predict_sd_pos <- x_mix_pred_vars_nums[["sd_pos"]] > 0L
predict_mu_neg <- x_mix_pred_vars_nums[["mu_neg"]] > 0L
predict_sd_neg <- x_mix_pred_vars_nums[["sd_neg"]] > 0L
predict_p_pos_binary <- length(p_pos_binary_pred_vars) > 0L

# 0-or-1 encode each cat of f_pred_vars for every plate
if (predict_f) {
  design_matrix_f <- df_plate %>%
    select(all_of(f_pred_vars_names)) %>%
    mutate(across(everything(), as.factor)) %>%
    {model.matrix(~ . - 1,
                  data = .,
                  contrasts.arg = lapply(.[, , drop = FALSE],
                                         contrasts, contrasts = FALSE))}
} else {
  design_matrix_f <- matrix(NA_real_, nrow = num_plate, ncol = 0)
}

# Count cats per f pred var. Ensure the col names of design_matrix_f are as 
# expected.
design_matrix_f_colnames_expected <-
  map(f_pred_vars_names, ~ paste0(.x, f_pred_vars[[.x]])) %>%
  unlist
stopifnot(identical(sort(colnames(design_matrix_f)),
                    sort(design_matrix_f_colnames_expected)))
if (!is.null(design_matrix_f_colnames_expected)) {
  design_matrix_f <- design_matrix_f[, design_matrix_f_colnames_expected]
}
stopifnot(identical(colnames(design_matrix_f),
                    design_matrix_f_colnames_expected))
stan_input_posterior$design_matrix_f <- design_matrix_f

# Look-ups between int and string encodings
if (predict_f) {
  df_f_pred_vars <- tibble(f_pred_var = f_pred_vars_names,
                           f_pred_var_int = 1:num_f_pred_vars)
  df_f_pred_vars_cats <- tibble(f_pred_var_cat = design_matrix_f_colnames_expected,
                                f_pred_var_cat_int = 1:num_f_pred_var_cats)
}

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
      select(id_sam, all_of(x_mix_pred_vars_names[[param]])) %>%
      distinct() # remove duplicates because of replicates...
    stopifnot(identical(mat$id_sam, # ... and check one row per sam_id, in order
                        1:nrow(mat)))
    mat <- mat %>%
      select(-id_sam) %>%
      mutate(across(everything(), as.character)) %>%
      mutate(across(everything(), as.factor)) %>%
      {model.matrix(~ . - 1,
                    data = .,
                    contrasts.arg = lapply(.[, , drop = FALSE],
                                           contrasts, contrasts = FALSE))}
    colnames_expected <-
      map(x_mix_pred_vars_names[[param]],
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

# Encode the design matrix for p_pos_binary 
if (! predict_p_pos_binary) {
  p_pos_binary_design_matrix <- matrix(NA_real_, nrow = num_sam_id, ncol = 0)
} else {
  mat <- df_sam %>%
    select(id_sam, all_of(p_pos_binary_pred_vars)) %>%
    distinct() 
  stopifnot(identical(mat$id_sam, 
                      1:nrow(mat)))
  mat <- mat %>%
    select(-id_sam) %>%
    mutate(across(everything(), as.integer))
  stopifnot(identical(colnames(mat),
                      p_pos_binary_pred_vars))
  p_pos_binary_design_matrix <- mat
}

for (param in x_mix_params) {
  stan_input_posterior[[paste0("design_matrix_", param)]] <-
    x_mix_design_matrices[[param]]
  stan_input_posterior[[paste0("num_", param, "_pred_vars")]] <- 
    x_mix_pred_vars_nums[[param]]
  stan_input_posterior[[paste0("num_cat_per_", param, "_pred_var")]] <- 
    x_mix_pred_vars_num_cats[[param]] %>% as.array()
}

stan_input_posterior$num_p_pos_binary_pred_vars <- length(p_pos_binary_pred_vars)
stan_input_posterior$design_matrix_p_pos_binary <- p_pos_binary_design_matrix

# Look-ups between int and string encodings of pred vars & their cats
lookup_pred_var_int <- list()
lookup_pred_var_cat_int <- list()
for (param in x_mix_params) {
  lookup_pred_var_int[[param]] <- tibble(
    pred_var = x_mix_pred_vars_names[[param]],
    pred_var_int = seq(from = 1, length.out = x_mix_pred_vars_nums[[param]]))
  lookup_pred_var_cat_int[[param]] <- tibble(
    pred_var_cat = colnames(x_mix_design_matrices[[param]]),
    pred_var_cat_int = seq(from = 1, length.out = x_mix_pred_vars_num_cats_tots[[param]]))
  if (nrow(lookup_pred_var_int[[param]]) == 0) {
    lookup_pred_var_int[[param]]$pred_var <- character()
  }
  if (nrow(lookup_pred_var_cat_int[[param]]) == 0) {
    lookup_pred_var_cat_int[[param]]$pred_var_cat <- character()
  }
}

# Priors
stan_input_posterior$rho_prior_eta <- rho_prior_eta
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

params_to_ignore <- c(
  "f_plate_effects_unscaled",
  "y_cal_mean_per_obs",
  "y_cal_mean_per_obs",
  "y_sam_mean_per_obs",
  "y_obs_sd_cal",
  "y_obs_sd_sam",
  "f_effects_by_pred_var_cat_unscaled",
  "p_pos_effects_by_pred_var_cat_unscaled",
  "mu_pos_effects_by_pred_var_cat_unscaled",
  "sd_pos_effects_by_pred_var_cat_unscaled",
  "mu_neg_effects_by_pred_var_cat_unscaled",
  "sd_neg_effects_by_pred_var_cat_unscaled",
  "p_pos_binary_effects_by_pred_var_unscaled",
  "exp_f_1_mult_f_4_per_plate",
  "p_pos_log_per_sam_id",
  "p_neg_log_per_sam_id",
  "p_pos_per_sam_id",
  "mu_pos_per_sam_id",
  "mu_neg_per_sam_id",
  "sd_pos_per_sam_id",
  "sd_neg_per_sam_id",
  "f_3_min_f_2_per_plate",
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
  "lp__",
  "p_blank_log",
  "p_blank_log1m",
  ".chain",
  ".iteration",
  ".draw",
  # actually interesting, but memory is limited:
  #"xlog_sam",
  #"y_sam_sim_conditional", 
  #"p_sam_is_pos", 
  #"y_sam_loglik_per_obs",
  "xlog_sam_sim_unconditional"
)

# Set up Stan
max_treedepth <- 14
if (stan_method == "cmdstanr") {
  model_compiled <- cmdstanr::cmdstan_model(file_input_stan)
} else if (stan_method == "rstan") {
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  model_compiled <- rstan::stan_model(file_input_stan)
} else if (stan_method != "cmdstan") {
  stop(paste("Unknown value", stan_method, "specified for stan_method"))
}

# GET THE POSTERIOR FROM STAN OR FROM FILE ----

if (! read_posterior_from_file) {
  
  start_time <- Sys.time()
  cat("Started running Stan at")
  print(start_time)
  
  if (stan_method == "cmdstanr") {
    samples_posterior <- model_compiled$sample(
      data = stan_input_posterior,
      iter_warmup = num_mc_iterations_posterior / 2,
      iter_sampling = num_mc_iterations_posterior / 2,
      chains = num_mc_chains,
      max_treedepth = max_treedepth,
      parallel_chains = num_mc_chains
    )
    df_fit_wide_postonly <- samples_posterior$draws(format = "draws_df")
    setDT(df_fit_wide_postonly)
    keep_col <- rep(TRUE, ncol(df_fit_wide_postonly))
    for (param in params_to_ignore) {
      keep_based_on_this_param <- 
        colnames(df_fit_wide_postonly) != param &
        ! startsWith(colnames(df_fit_wide_postonly), paste0(param, ".")) 
      keep_col <- keep_col & keep_based_on_this_param
    }
    df_fit_wide_postonly <- df_fit_wide_postonly[, ..keep_col]
    
  } else if (stan_method == "rstan") {
    samples_posterior <- sampling(model_compiled,
                                  data = stan_input_posterior,
                                  iter = num_mc_iterations_posterior,
                                  chains = num_mc_chains,
                                  control = list(max_treedepth = max_treedepth),
                                  pars = params_to_ignore,
                                  include = FALSE)
    df_fit_wide_postonly <- samples_posterior %>%
      as.data.frame()
    
  } else if (stan_method == "cmdstan") {
    cmdstanr::write_stan_json(stan_input_posterior, file = file_stan_temp)
    stopifnot(endsWith(file_input_stan, ".stan"))
    file_stan_exe <- str_remove(file_input_stan, ".stan$")
    system(paste("cd", dir_stan, " && make STAN_THREADS=true", file_stan_exe))
    time <- format(Sys.time(), "%Y-%m-%d-%Hh%Mm%S")
    files_out_stan <- paste0(file_out_stan_basename, time, "_chain", 1:num_mc_chains, ".csv")
    files_out_stan_summary <- paste0(file_out_stan_basename, time, "_summary.csv")
    files_out_stan_profile <- paste0(file_out_stan_basename, time, "_profiles.csv")
    files_out_r_image <- paste0(file_out_stan_basename, time, ".RData")
    command <- paste0(file_stan_exe,
                      " method=sample",
                      " num_chains=", num_mc_chains,
                      " num_warmup=", round(num_mc_iterations_posterior / 2),
                      " num_samples=", round(num_mc_iterations_posterior / 2),
                      " num_threads=", num_mc_chains,
                      " data file=", file_stan_temp, 
                      " output file=", paste(files_out_stan, collapse = ","),
                      " profile_file=", files_out_stan_profile)
    print("Running this command:")
    print(command)
    system(command)
    system(paste0(dir_stan, "/bin/stansummary ", files_out_stan,
                  " --percentiles 50",
                  " --csv_filename ", files_out_stan_summary))
    save.image(files_out_r_image)
    
  } else {
    stop(paste("Unknown value", stan_method, "specified for stan_method"))
  }
  
  end_time <- Sys.time()
  cat("Finished running Stan at")
  print(end_time)
  print(end_time - start_time)
  
  #samples_posterior$save_output_files("/Users/cwymant/enable/", basename = "samples_full_run")
  # rm(samples_posterior)
  # save.image("/Users/cwymant/enable/samples_full_run_TODO_DATE.RData")
  
  # samples_posterior$profiles() 
  
}

if (read_posterior_from_file || stan_method == "cmdstan") {
  stopifnot(length(files_out_stan) > 0)
  stopifnot(all(file.exists(files_out_stan)))
  df_fit_wide_postonly <- map(files_out_stan, function(file_){
    print(Sys.time())
    cat("Now reading file", file_, "\n")
    df_ <- data.table::fread(cmd = paste("grep -v '^#'", file_))
    keep_col <- rep(TRUE, ncol(df_))
    for (param in params_to_ignore) {
      keep_based_on_this_param <- 
        colnames(df_) != param &
        ! startsWith(colnames(df_), paste0(param, ".")) 
      keep_col <- keep_col & keep_based_on_this_param
    }
    df_ <- df_[, ..keep_col]
    df_ <- df_[seq(1, .N, by = 2)] # Keep every other row for memory management
    df_
  }) %>% data.table::rbindlist()
  
  data.table::setnames(df_fit_wide_postonly, mastiff::rename_params_cmdstanfile_to_rstan)
}

# GET THE PRIOR FROM STAN ----

if (stan_method %in% c("cmdstanr", "cmdstan")) {
  model_compiled <- cmdstanr::cmdstan_model(file_input_stan)
  df_ps <- model_compiled$sample(
    data = stan_input_prior,
    iter_warmup = num_mc_iterations_prior / 2,
    iter_sampling = num_mc_iterations_prior / 2,
    chains = num_mc_chains,
    max_treedepth = max_treedepth,
    parallel_chains = num_mc_chains
  )
  df_ps <- df_ps$draws(format = "draws_df")
  setDT(df_ps)
  keep_col <- rep(TRUE, ncol(df_ps))
  for (param in params_to_ignore) {
    keep_based_on_this_param <- 
      colnames(df_ps) != param &
      ! startsWith(colnames(df_ps), paste0(param, ".")) &
      ! startsWith(colnames(df_ps), paste0(param, "[")) 
    keep_col <- keep_col & keep_based_on_this_param
  }
  df_ps <- df_ps[, ..keep_col]
} else if (stan_method == "rstan") {
  df_ps <- sampling(model_compiled,
                    data = stan_input_prior,
                    iter = num_mc_iterations_prior,
                    chains = num_mc_chains,
                    control = list(max_treedepth = max_treedepth),
                    pars = params_to_ignore,
                    include = FALSE)
  df_ps <- as.data.frame(df_ps)
  setDT(df_ps)
} else {
  stop(paste("Unknown value", stan_method, "specified for stan_method"))
}
df_ps[, sample := 1:nrow(df_ps)]
df_ps[, density_type := "prior"]


# WRANGLE STAN OUTPUT ----

df_fit_wide_postonly[, density_type := "posterior"]
df_fit_wide_postonly[, sample := 1:nrow(df_fit_wide_postonly)]


# Merge prior and posterior samples
desired_cols <- names(df_ps)
df_fit_wide_postandprior <- rbind(df_fit_wide_postonly[,..desired_cols],
                                  df_ps) 

# Record true values of params
if (data_was_simulated) {
  
  df_true_pop_params <- tribble(
    ~param, ~value,
    "f[1]", f[1],
    "f[2]", f[2],
    "f[3]", f[3],
    "f[4]", f[4],
    "sigma_f_plate[1]", sigma_f_plate[1],
    "sigma_f_plate[2]", sigma_f_plate[2],
    "sigma_f_plate[3]", sigma_f_plate[3],
    "sigma_f_plate[4]", sigma_f_plate[4],
    "rho[1,2]", rho[1,2],
    "rho[1,3]", rho[1,3],
    "rho[1,4]", rho[1,4],
    "rho[2,3]", rho[2,3],
    "rho[2,4]", rho[2,4],
    "rho[3,4]", rho[3,4],
    "sd_neg", sd_neg,
    "sd_pos", sd_pos,
    "mu_neg", mu_neg,
    "mu_pos", mu_pos,
    "p_pos", p_pos,
    "p_blank", p_blank,
    "y_obs_sd_cal_min",  y_obs_sd_cal_min,
    "y_obs_sd_cal_jump", y_obs_sd_cal_jump,
    "y_obs_sd_sam_min",  y_obs_sd_sam_min,
    "y_obs_sd_sam_jump", y_obs_sd_sam_jump,
    "y_obs_sd_cal_max", y_obs_sd_cal_max,
    "y_obs_sd_sam_max", y_obs_sd_sam_max
  ) 
  
  if (predict_f) {
    df_true_f_effects_by_pred_var <- f_effects_by_pred_var %>% 
      map(function(mat) {mat %>%
          as_tibble(.name_repair = "universal_quiet") %>%
          mutate(cat = rownames(mat))}) %>% 
      bind_rows(.id = "f_pred_var") %>%
      mutate(f_pred_var_cat = paste0(f_pred_var, cat)) %>% 
      inner_join(df_f_pred_vars_cats, by = "f_pred_var_cat")
    stopifnot(identical(sort(df_true_f_effects_by_pred_var$f_pred_var_cat),
                        sort(design_matrix_f_colnames_expected)))
    df_true_f_effects_by_pred_var <- df_true_f_effects_by_pred_var %>%
      pivot_longer(paste0("...", 1:4),
                   names_prefix = "...",
                   names_to = "which_f") %>%
      mutate(param = paste0("f_effects_by_pred_var_cat[", f_pred_var_cat_int,
                            ",", which_f, "]")) %>%
      select(param, value)
    df_true_pop_params <- df_true_pop_params %>%
      bind_rows(df_true_f_effects_by_pred_var,
                tibble(f_pred_var = names(sigma_f_pred_vars),
                       value = sigma_f_pred_vars) %>% 
                  unnest_longer(value, indices_to = "which_f") %>%
                  left_join(df_f_pred_vars, by = "f_pred_var") %>%
                  mutate(param = paste0("sigma_f_pred_vars[", f_pred_var_int,
                                        ",", which_f, "]")) %>%
                  select(param, value))
  }
  
  df_true_pop_params <- df_true_pop_params %>% bind_rows(
    map(x_mix_params, function(param) {
      if (x_mix_pred_vars_nums[[param]] == 0) {
        return(tibble(param = character(), value = numeric()))
      }
      df_true_effects_by_pred_var <- x_mix_effects[[param]] %>% 
        map(function(mat) {mat %>%
            as_tibble(.name_repair = "universal_quiet") %>%
            mutate(cat = names(mat))}) %>% 
        bind_rows(.id = "pred_var") %>%
        mutate(pred_var_cat = paste0(pred_var, cat)) %>% 
        inner_join(lookup_pred_var_cat_int[[param]], by = "pred_var_cat")
      stopifnot(identical(sort(df_true_effects_by_pred_var$pred_var_cat),
                          sort(colnames(x_mix_design_matrices[[param]]))))
      df_true_effects_by_pred_var <- df_true_effects_by_pred_var %>%
        mutate(param = paste0(param, "_effects_by_pred_var_cat[", 
                              pred_var_cat_int, "]")) %>%
        select(param, value)
      df_true_sigma_pred_vars <-
        tibble(pred_var = names(x_mix_pred_vars_sds[[param]]),
               value = x_mix_pred_vars_sds[[param]]) %>%
        inner_join(lookup_pred_var_int[[param]], by = "pred_var") %>%
        mutate(param = paste0("sigma_", param, "_pred_vars[", pred_var_int, "]")) %>%
        select(param, value)
      bind_rows(df_true_effects_by_pred_var, df_true_sigma_pred_vars)
    }) %>%
      bind_rows())
  
  if (predict_p_pos_binary) {
    df_true_pop_params <- df_true_pop_params %>%
      bind_rows(tibble(
        param = paste0("p_pos_binary_effects_by_pred_var[",
                       1:length(p_pos_binary_pred_vars), "]"),
        value = p_pos_binary_effects)) %>%
      bind_rows(tibble(
        param = "sigma_p_pos_binary_pred_vars[1]",
        value = p_pos_binary_pred_vars_sd
      ))
  }

}

# Rename params for interpretability
rename_params <- function(original_names) {
  df_param_names <- tibble(orig = original_names,
                           new = str_replace(orig,
                                             "sigma_f_plate\\[([0-9]+)\\]",
                                             "sigma_f[\\1]_plate"))
  if (predict_f) {
    df_param_names <- df_param_names %>%
      tidyr::extract(orig, 
                     into = c("f_pred_var_int", "which_f_foo"), 
                     regex = "sigma_f_pred_vars\\[([0-9]+),([0-9]+)\\]",
                     remove = FALSE) %>%
      tidyr::extract(orig, 
                     into = c("f_pred_var_cat_int", "which_f_spam"), 
                     regex = "f_effects_by_pred_var_cat\\[([0-9]+),([0-9]+)\\]",
                     remove = FALSE) %>%
      mutate(f_pred_var_int = as.integer(f_pred_var_int),
             f_pred_var_cat_int = as.integer(f_pred_var_cat_int)) %>%
      left_join(df_f_pred_vars, by = "f_pred_var_int") %>%
      left_join(df_f_pred_vars_cats, by = "f_pred_var_cat_int") %>% 
      mutate(new = case_when(
        !is.na(f_pred_var_int) ~ paste0("sigma_f[", which_f_foo, "]_", f_pred_var),
        !is.na(f_pred_var_cat_int) ~ paste0("f_effect[", which_f_spam, "]_", f_pred_var_cat),
        TRUE ~ new
      ))
  }
  
  for (param in x_mix_params) {
    if (x_mix_pred_vars_nums[[param]] == 0) next
    df_param_names <- df_param_names %>%
      tidyr::extract(orig, 
                     into = "pred_var_int", 
                     regex = paste0("sigma_", param, "_pred_vars\\[([0-9]+)\\]"),
                     remove = FALSE) %>%
      tidyr::extract(orig, 
                     into = "pred_var_cat_int", 
                     regex = paste0(param, "_effects_by_pred_var_cat\\[([0-9]+)\\]"),
                     remove = FALSE) %>%
      mutate(pred_var_int = as.integer(pred_var_int),
             pred_var_cat_int = as.integer(pred_var_cat_int)) %>%
      left_join(lookup_pred_var_int[[param]], by = "pred_var_int") %>%
      left_join(lookup_pred_var_cat_int[[param]], by = "pred_var_cat_int") %>% 
      mutate(new = case_when(
        !is.na(pred_var_int) ~ paste0("sigma_", param, "_pred_vars_", pred_var),
        !is.na(pred_var_cat_int) ~ paste0(param, "_effect_", pred_var_cat),
        TRUE ~ new
      )) %>%
      select(orig, new)
  }
  if (predict_p_pos_binary) {
    df_param_names <- df_param_names %>%
      tidyr::extract(orig, 
                     into = "pred_var_int", 
                     regex = paste0("p_pos_binary_effects_by_pred_var\\[([0-9]+)\\]"),
                     remove = FALSE) %>%
      mutate(pred_var_int = as.integer(pred_var_int)) %>%
      left_join(tibble(pred_var_int = 1:length(p_pos_binary_pred_vars),
                       pred_var = p_pos_binary_pred_vars),
                by = "pred_var_int") %>% 
      mutate(new = case_when(
        !is.na(pred_var_int) ~ paste0("p_pos_effect_", pred_var),
        TRUE ~ new
      )) %>%
      select(orig, new) %>%
      mutate(new = if_else(orig == "sigma_p_pos_binary_pred_vars[1]",
                           "sigma_p_pos_binary_pred_vars",
                           new))
  }
  df_param_names$new
}
setnames(df_fit_wide_postonly, rename_params)
setnames(df_fit_wide_postandprior, rename_params)
setnames(df_ps, rename_params)
if (data_was_simulated) df_true_pop_params$param <- rename_params(df_true_pop_params$param)

# Define x mix params by group, from overall params + effects
for (param in x_mix_params) {
  if (x_mix_pred_vars_nums[[param]] == 0) next
  for (pred_var_cat in lookup_pred_var_cat_int[[param]]$pred_var_cat) {
    if (param == "p_pos") {
      df_fit_wide_postandprior[[paste0(param, "_for_", pred_var_cat)]] <- mastiff::logistic(
        mastiff::logit(df_fit_wide_postandprior[[param]]) +
          df_fit_wide_postandprior[[paste0(param, "_effect_", pred_var_cat)]])
      df_ps[[paste0(param, "_for_", pred_var_cat)]] <- mastiff::logistic(
        mastiff::logit(df_ps[[param]]) + df_ps[[paste0(param, "_effect_", pred_var_cat)]])
    } else if (param %in% c("sd_pos", "sd_neg")) {
      df_fit_wide_postandprior[[paste0(param, "_for_", pred_var_cat)]] <- exp(
        log(df_fit_wide_postandprior[[param]]) + 
          df_fit_wide_postandprior[[paste0(param, "_effect_", pred_var_cat)]])
      df_ps[[paste0(param, "_for_", pred_var_cat)]] <- exp(
        log(df_ps[[param]]) + df_ps[[paste0(param, "_effect_", pred_var_cat)]])
    }
    else {
      df_fit_wide_postandprior[[paste0(param, "_for_", pred_var_cat)]] <-
        df_fit_wide_postandprior[[param]] + 
        df_fit_wide_postandprior[[paste0(param, "_effect_", pred_var_cat)]]
      df_ps[[paste0(param, "_for_", pred_var_cat)]] <-
        df_ps[[param]] + df_ps[[paste0(param, "_effect_", pred_var_cat)]]
    }
  }
}
if (data_was_simulated) {
  df_true_pop_params <- df_true_pop_params %>% bind_rows(
    x_mix_params %>% map(function(param) {
      if (x_mix_pred_vars_nums[[param]] == 0) {
        return(tibble(param = character(), value = numeric()))
      }
      x_mix_overall[[param]] %>%
        map(function(mat) {mat %>%
            as_tibble(.name_repair = "universal_quiet") %>%
            mutate(cat = names(mat))}) %>%
        bind_rows(.id = "pred_var") %>%
        mutate(param = paste0(param, "_for_", pred_var, cat)) %>%
        select(param, value)
    }) 
  )
}

if ("p_pos_for_sex_Male"   %in% names(df_fit_wide_postandprior) &&
    "p_pos_for_sex_Female" %in% names(df_fit_wide_postandprior)) {
  df_fit_wide_postandprior$p_pos_for_sex_Male_min_Female <-
    df_fit_wide_postandprior$p_pos_for_sex_Male -
    df_fit_wide_postandprior$p_pos_for_sex_Female
}

# PLOT STAN OUTPUT ----

param_names_pop <- colnames(df_fit_wide_postandprior)
param_names_pop <- param_names_pop[param_names_pop != "sample"]

regex_for_params_to_plot <- ""
params_desired <- param_names_pop[grepl(regex_for_params_to_plot, param_names_pop)]
params_for_log_transform <- c() # params_desired[grepl("sd_", params_desired)]
list_for_log_transform <- list()
for (param in params_for_log_transform) list_for_log_transform[[param]] <- log
if (data_was_simulated) {
  true_params_to_plot <- df_true_pop_params[grepl(regex_for_params_to_plot,
                                                  df_true_pop_params$param),]$value
  names(true_params_to_plot) <- df_true_pop_params[grepl(regex_for_params_to_plot,
                                                         df_true_pop_params$param),]$param
  p <- mastiff::plot_posterior(
    posterior_samples = df_fit_wide_postandprior %>% select(-sample) %>% filter(density_type == "posterior"),
    prior_samples     = df_fit_wide_postandprior %>% select(-sample) %>% filter(density_type == "prior"),
    params_desired = params_desired,
    true_param_values = true_params_to_plot,
    skip_stanfit_to_dt = TRUE)
} else {
  p <- mastiff::plot_posterior(
    posterior_samples = df_fit_wide_postandprior %>% select(-sample) %>% filter(density_type == "posterior"),
    prior_samples     = df_fit_wide_postandprior %>% select(-sample) %>% filter(density_type == "prior"),
    params_desired = params_desired,
    skip_stanfit_to_dt = TRUE)
}
p

 
# Compare true and estimated sample x 
quantiles <- c(0.025, 0.5, 0.975)
df_sam_x <- df_fit_wide_postonly %>%
  select(sample, starts_with("xlog_sam[")) %>%
  pivot_longer(-sample, names_to = "param") %>%
  tidyr::extract(param, 
                 into = "id_sam", 
                 regex = "xlog_sam\\[([0-9]+)\\]") %>%
  mutate(id_sam = as.integer(id_sam)) %>%
  group_by(id_sam) %>%
  reframe(value = quantile(value, probs = quantiles),
          quantile = quantiles) %>%
  pivot_wider(names_from = quantile, names_prefix = "x_q_")
if (data_was_simulated) {
  df_sam_x %>%
    left_join(df_sam %>%
                select(id_sam, xlog) %>%
                distinct(),
              by = "id_sam") %>%
    ggplot() +
    geom_errorbar(aes(xlog, ymin = x_q_0.025, ymax = x_q_0.975)) +
    geom_point(aes(xlog, x_q_0.5)) +
    geom_abline() +
    labs(x = "True Ab",
         y = "Estimated Ab")
}

# Classification plot
quantiles <- c(0.025, 0.5, 0.975)
df_prob_pos <- df_fit_wide_postonly %>%
  select(sample, starts_with("p_sam_is_pos[")) %>%
  pivot_longer(-sample, names_to = "param") %>%
  tidyr::extract(param, 
                 into = "id_sam", 
                 regex = "p_sam_is_pos\\[([0-9]+)\\]") %>%
  mutate(id_sam = as.integer(id_sam)) %>%
  group_by(id_sam) %>%
  reframe(value = quantile(value, probs = quantiles),
          quantile = quantiles) %>%
  pivot_wider(names_from = quantile, names_prefix = "prob_pos_q_")
if (data_was_simulated) {
  inner_join(df_prob_pos,
             df_sam %>% summarise(.by = id_sam, pos = unique(pos)), 
             by = "id_sam") %>%
    mutate(pos = if_else(pos, "pos", "neg")) %>%
    ggplot() +
    geom_histogram(aes(prob_pos_q_0.5, fill = pos), #y = after_stat(density))
                   position = "identity",
                   alpha = 0.6,
                   bins = 30) +
    labs(y = "Number of samples",
         x = "Probability sample is positive",
         fill = "Truth:") +
    coord_cartesian(expand = FALSE) +
    scale_x_continuous(limits = c(NA, NA))
}

# Compare prob positivity vs OD
inner_join(df_prob_pos,
           df_sam %>% summarise(.by = id_sam, y = mean(y)), 
           by = "id_sam") %>%
  ggplot() +
  geom_point(aes(y, prob_pos_q_0.5)) +
  labs(x = "Mean observed OD",
       y = "Estimated probability of being positive")  

# Compare prob positivity vs x
inner_join(df_sam_x, df_prob_pos, by = "id_sam") %>%
  ggplot() +
  geom_point(aes(x_q_0.5, prob_pos_q_0.5)) +
  labs(x = "Estimated concentration",
       y = "Estimated probability of being positive") 

# Plot cal data by plate  
ggplot() +
  geom_point(data = df_cal,
             aes(x, y)) +
  facet_wrap(~plate) +
  scale_x_log10() +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 

# The posteriors for the 4PL function by plate
xlog_range <- seq(log(0.002), log(40), length.out = 50)
df_4pl <- df_fit_wide_postonly %>%
  select(sample, starts_with("f_per_plate[")) %>%
  filter(sample %% 10 == 0) %>%
  pivot_longer(-sample, names_to = "param") %>%
  tidyr::extract(param, 
                 into = c("plate_int", "f_index"), 
                 regex = "f_per_plate\\[([0-9]+),([0-9]+)\\]") %>%
  mutate(plate_int = as.integer(plate_int),
         f_index = as.integer(f_index)) %>%
  left_join(df_plate %>% select(plate, plate_int), by = "plate_int") %>%
  pivot_wider(names_from = f_index, names_prefix = "f_") %>%
  expand_grid(xlog = xlog_range) %>%
  mutate(x = exp(xlog)) %>%
  mutate(y = PL4(xlog, f_1, f_2, f_3, f_4))
p <- ggplot(df_4pl %>%
              filter(sample %% 10 == 0)) +
  geom_line(aes(x, y, group = sample), alpha = 0.1) +
  geom_point(data = df_cal,
             aes(x, y), col = "blue") +
  facet_wrap(~plate, nrow = 3) +
  scale_x_log10() +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
if (data_was_simulated) {
  p <- p + 
    geom_line(data = df_plate %>%
                expand_grid(xlog = xlog_range) %>%
                mutate(x = exp(xlog)) %>%
                mutate(y = PL4(xlog, f_1, f_2, f_3, f_4)),
              aes(x, y), col = "blue", linewidth = 1)
}
p

# Plot P(x | pos), P(x | neg), P(x), P(pos | x) again but now with logx
xlogs_plot <- log(10) * -90:60 / 30
df_xlog_distributions <- df_fit_wide_postonly %>%
  filter(sample %% 10 == 0) %>%
  select("sample", "mu_pos", "sd_pos",
         "mu_neg", "sd_neg", "p_pos") %>%
  full_join(tibble(xlog = xlogs_plot,
                   x = exp(xlog)),
            by = character()) %>%
  mutate(`P(xlog | pos)` = dnorm(xlog, mean = mu_pos, sd = sd_pos),
         `P(xlog | neg)` = dnorm(xlog, mean = mu_neg, sd = sd_neg),
         `P(xlog)` = p_pos * `P(xlog | pos)` + (1 - p_pos) * `P(xlog | neg)`,
         `P(pos | xlog)` = p_pos * `P(xlog | pos)` / `P(xlog)`) %>%
  select(sample, xlog, `P(xlog | pos)`, 
         `P(xlog | neg)`, `P(xlog)`, `P(pos | xlog)`) %>%
  pivot_longer(c("P(xlog | pos)", "P(xlog | neg)", "P(xlog)", "P(pos | xlog)"))
p <- ggplot() +
  geom_line(data = df_xlog_distributions,
            aes(x = xlog, y = value, group = sample), alpha = 0.15) +
  facet_wrap(vars(name), scales = "free_y", ncol = 1) +
  labs(x = "log_e(Ab concentration)",
       y = "") +
  scale_x_continuous(expand = c(0, 0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, NA))
if (data_was_simulated) {
  df_x_distributions_truth <-
    tibble(xlog = xlogs_plot,
           x = exp(xlog),
           `P(xlog | pos)` = dnorm(xlog, mean = mu_pos, sd = sd_pos),
           `P(xlog | neg)` = dnorm(xlog, mean = mu_neg, sd = sd_neg),
           `P(xlog)` = p_pos * `P(xlog | pos)` + (1 - p_pos) * `P(xlog | neg)`,
           `P(pos | xlog)` = p_pos * `P(xlog | pos)` / `P(xlog)`) %>%
    pivot_longer(-c("x", "xlog"))
  p <- p +
    geom_line(data = df_x_distributions_truth,
              aes(x = xlog, y = value), colour = "blue") 
}
p

# Plot the posterior distribution of the population level distribution of point
# estimates of x_sam and the posterior distribution of the parametric
# pop-level distribution of x_sam
df_fit_wide_postonly %>%
  filter(sample %% 20 == 0) %>%
  select(sample, starts_with("xlog_sam[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  mutate(value = exp(value)) %>%
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = log(1.8)) +
  geom_line(data = df_xlog_distributions %>% filter(name == "P(xlog)"),
            aes(x = xlog, y = value, group = sample), alpha = 0.15, col = "blue") +
  geom_density(aes(log(value) , group = sample), alpha = 0.01) +
  coord_cartesian(expand = F) +
  labs(x = "log_e(Ab concentration)",
       y = "probability density") +
  NULL 

# Plot the posterior distribution of the population level distribution of 
# stochastically redrawn y_sam_sim_conditional 
df_fit_wide_postonly %>%
  filter(sample %% 30 == 0) %>%
  select(sample, starts_with("y_sam_sim_conditional[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  mutate(value = log10(value)) %>%
  ggplot() +
  geom_histogram(data = df_sam, 
                 aes(x = log10(y), y = after_stat(density)), 
                 bins = 60) +
  geom_density(aes(value, group = sample), alpha = 0.01) +
  coord_cartesian(expand = F) +
  scale_x_continuous(limits = c(-3, 1)) +
  labs(x = "log10(OD value)",
       y = "probability density") +
  NULL 

# Plot the posterior distribution of the population level distribution of 
# stochastically redrawn y_sam_sim_unconditional 
df_fit_wide_postonly %>%
  filter(sample %% 100 == 0) %>%
  select(sample, starts_with("y_sam_sim_unconditional[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  tidyr::extract(param, 
                 into = "which_sam_rep", 
                 regex = "y_sam_sim_unconditional\\[([0-9]+)\\]") %>%
  mutate(which_sam_rep = as.integer(which_sam_rep)) %>%
  left_join(df_sam, by = "which_sam_rep") %>%
  mutate(value = log10(value)) %>%
  ggplot() +
  geom_density(aes(value, group = sample), alpha = 0.01) +
  geom_histogram(data = df_sam, 
                 aes(x = log10(y), y = after_stat(density)),
                 col = "blue", fill = NA,
                 bins = 60) +
  coord_cartesian(expand = F) +
  scale_x_continuous(limits = c(-3, 1)) +
  labs(x = "log10(OD value)",
       y = "probability density") +
  NULL 



