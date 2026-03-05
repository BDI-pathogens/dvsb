library(data.table)

data_was_simulated <- FALSE
read_samples_from_file <- TRUE
#files_samples <- Sys.glob("~/enable/samples_full_run-202509160846-*-76908a.csv") # non-baseline, no observation mixture model
#files_samples <- Sys.glob("~/enable/samples_full_run-202509180821-*-54ceae.csv") # baseline only, mixture observation with p_sam_rep_is_blank and correct y_sim
#files_samples <- Sys.glob("~/enable/samples_full_run-202509192120-*-212f13.csv") # non-baseline, mixture observation without p_sam_rep_is_blank and with y_sim excluding blank possibility
#files_samples <- Sys.glob("~/enable/samples_full_run-202509231402-*-34f4eb.csv") # all calibration only, including all blanks to identify which are outliers
files_samples <- Sys.glob("/Users/cwymant/enable/samples_full_run-202509240958-*-8dbbf0.csv") # non-baseline, mixture observation without p_sam_rep_is_blank and with y_sim excluding blank possibility, after removing some outliers

# INPUT ABOUT STAN ----

file_input_stan <- "~/PathogenDynamics Dropbox/Vaccine Work/Lassa/code_serology_model/Xsectional_v16.stan"
sample_prior_manually <- TRUE
num_mc_chains <- 4
num_mc_iterations_posterior <- 500
num_mc_iterations_prior <- 4000
use_cmdstanr <- TRUE

# Upper and lower bounds for priors
if (data_was_simulated) {
  
  df_priors_scalars <- tribble(
    ~param, ~lower, ~upper,
    "x_sam_neg_mu", x_sam_neg_mu - 1, x_sam_neg_mu + 1,
    "x_sam_pos_mu", x_sam_pos_mu - 1, x_sam_pos_mu + 1,
    "x_sam_neg_sd", 0, x_sam_neg_sd * 2,
    "x_sam_pos_sd", 0, x_sam_pos_sd * 2,
    "p_pos", 0, 1,
    "p_blank", 0, 2 * p_blank,
    "y_obs_sd_cal_min",  0, 2 * y_obs_sd_cal_min,
    "y_obs_sd_cal_jump", 0, 2 * y_obs_sd_cal_jump,
    "y_obs_sd_sam_min",  0, 2 * y_obs_sd_sam_min,
    "y_obs_sd_sam_jump", 0, 2 * y_obs_sd_sam_jump,
    "sigma_p_pos_pred_vars", 0, 2 * max(sigma_p_pos_pred_vars)
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
    "x_sam_neg_mu", -5, -1,
    "x_sam_neg_sd", 0, 2,
    "x_sam_pos_sd", 0, 2.5, 
    "x_sam_pos_mu", -1, 3,
    "p_pos", 0, 1,
    "p_blank", 0, 0.05,
    "y_obs_sd_cal_min", 0, 0.015,
    "y_obs_sd_cal_jump", 0.1, 1,
    "y_obs_sd_sam_min", 0, 0.03,
    "y_obs_sd_sam_jump", 0.1, 1,
    "sigma_p_pos_pred_vars", 0, 2
  )
  df_priors_vectors <- tribble(
    ~param, ~lower, ~upper,
    "sigma_f_plate", c(0, 0, 0, 0), c(0.25, 0.025, 2.5, 2),
    "f", c(0.8, -0.05, 3, 2.2), c(1.1, 0.05, 5.5, 3.8),
    "sigma_f_pred_vars", c(0, 0, 0, 0), c(0.6, 0.1, 4, 3)
  )
  
}

# The eta parameter of the LKJ prior for rho
rho_prior_eta <- 2

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

# FINISH PREPARING INPUT TO STAN ----

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
  design_matrix_f <- matrix(nrow = num_plate, ncol = 0)
}

# Count cats per f pred var. Ensure the col names of design_matrix_f are as 
# expected.
design_matrix_f_colnames_expected <-
  map(f_pred_vars_names, ~ paste0(.x, f_pred_vars[[.x]])) %>%
  unlist
stopifnot(identical(sort(colnames(design_matrix_f)),
                    sort(design_matrix_f_colnames_expected)))
design_matrix_f <- design_matrix_f[, design_matrix_f_colnames_expected]
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

# 0-or-1 encode each cat of p_pos_pred_vars for every sample
if (predict_p_pos) {
  design_matrix_p_pos <- df_sam %>%
    select(id_sam, all_of(p_pos_pred_vars_names)) %>%
    distinct()
  stopifnot(identical(design_matrix_p_pos$id_sam,
                      1:nrow(design_matrix_p_pos)))
  design_matrix_p_pos <- design_matrix_p_pos %>%
    select(-id_sam) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), as.factor)) %>%
    {model.matrix(~ . - 1,
                  data = .,
                  contrasts.arg = lapply(.[, , drop = FALSE],
                                         contrasts, contrasts = FALSE))}
} else {
  design_matrix_p_pos <- matrix(nrow = num_sam_id, ncol = 0)
}

# Count cats per p_pos pred var. Ensure the col names of design_matrix_p_pos are as 
# expected.
design_matrix_p_pos_colnames_expected <-
  map(p_pos_pred_vars_names, ~ paste0(.x, p_pos_pred_vars[[.x]])) %>%
  unlist
stopifnot(identical(sort(colnames(design_matrix_p_pos)),
                    sort(design_matrix_p_pos_colnames_expected)))
design_matrix_p_pos <- design_matrix_p_pos[, design_matrix_p_pos_colnames_expected]
stopifnot(identical(colnames(design_matrix_p_pos),
                    design_matrix_p_pos_colnames_expected))
stan_input_posterior$design_matrix_p_pos <- design_matrix_p_pos

# Look-ups between int and string encoding
if (predict_p_pos) {
  df_p_pos_pred_vars <- tibble(p_pos_pred_var = p_pos_pred_vars_names,
                               p_pos_pred_var_int = 1:num_p_pos_pred_vars)
  df_p_pos_pred_vars_cats <- tibble(p_pos_pred_var_cat = design_matrix_p_pos_colnames_expected,
                                    p_pos_pred_var_cat_int = 1:num_p_pos_pred_var_cats)
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
stan_input_prior <- stan_input_posterior
stan_input_prior$sample_posterior_not_prior <- 0L

params_to_ignore <- c(
  "f_plate_effects_unscaled",
  "y_cal_mean_per_obs",
  "y_cal_mean_per_obs",
  "y_sam_mean_per_obs",
  "y_obs_sd_cal",
  "y_obs_sd_sam",
  "f_effects_by_pred_var_cat_unscaled",
  "p_pos_effects_by_pred_var_cat_unscaled",
  "exp_f_1_mult_f_4_per_plate",
  "p_pos_log_per_sam_id",
  "p_neg_log_per_sam_id",
  "f_3_min_f_2_per_plate"
)

# RUN STAN ----

if (read_samples_from_file) {
  
  df_fit_wide_postonly <- map(files_samples, function(file_){
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
    df_
  }) %>% data.table::rbindlist()
  
  data.table::setnames(df_fit_wide_postonly, mastiff::rename_params_cmdstanfile_to_rstan)
  #colnames(df_fit_wide_postonly) <- 
  #  mastiff::rename_params_cmdstanfile_to_rstan(colnames(df_fit_wide_postonly))
  
} else {
  
  if (use_cmdstanr) {
    library(cmdstanr)
  } else {
    library(rstan)
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
  }
  
  # Compile the Stan code
  if (use_cmdstanr) {
    model_compiled <- cmdstan_model(file_input_stan)
  } else {
    model_compiled <- stan_model(file_input_stan)
  }
  
  # Run the Stan code
  start_time <- Sys.time()
  cat("Started running Stan at")
  print(start_time)
  max_treedepth <- 12
  if (use_cmdstanr) {
    
    samples_posterior <- model_compiled$sample(
      data = stan_input_posterior,
      iter_warmup = num_mc_iterations_posterior / 2,
      iter_sampling = num_mc_iterations_posterior / 2,
      chains = num_mc_chains,
      max_treedepth = max_treedepth,
      parallel_chains = num_mc_chains
    )
    df_fit_wide_postonly <- samples_posterior$draws(format = "draws_df")
    if (! sample_prior_manually) {
      samples_prior <- model_compiled$sample(
        data = stan_input_prior,
        iter_warmup = num_mc_iterations_posterior / 2,
        iter_sampling = num_mc_iterations_posterior / 2,
        chains = num_mc_chains,
        max_treedepth = max_treedepth,
        parallel_chains = num_mc_chains
      )
    }
    
  } else {
    
    samples_posterior <- sampling(model_compiled,
                                  data = stan_input_posterior,
                                  iter = num_mc_iterations_posterior,
                                  chains = num_mc_chains,
                                  control = list(max_treedepth = max_treedepth),
                                  pars = params_to_ignore,
                                  include = FALSE)
    df_fit_wide_postonly <- samples_posterior %>%
      as.data.frame()
    
    if (! sample_prior_manually) {
      samples_prior <- sampling(model_compiled,
                                data = stan_input_prior,
                                iter = num_mc_iterations_prior,
                                chains = num_mc_chains,
                                control = list(max_treedepth = max_treedepth),
                                pars = params_to_ignore,
                                include = FALSE)
    }
  }
  end_time <- Sys.time()
  cat("Finished running Stan at")
  print(end_time)
  print(end_time - start_time)
  
  #samples_posterior$save_output_files("~/enable/", basename = "samples_full_run")
  # rm(samples_posterior)
  # save.image("~/enable/samples_full_run_TODO_DATE.RData")
  
}

# WRANGLE STAN OUTPUT ----

setDT(df_fit_wide_postonly)
df_fit_wide_postonly[, density_type := "posterior"]
df_fit_wide_postonly[, sample := 1:nrow(df_fit_wide_postonly)]

# Ugly code to sample from the prior manually. Sorry programming.
if (sample_prior_manually) {
  
  # Combine df_priors_scalars and df_priors_vectors
  df_priors <- df_priors_scalars
  df_pri_vec_sampling <- df_priors_vectors %>% 
    mutate(replicates_needed = case_when(
      param == "sigma_f_pred_vars" ~ num_f_pred_vars,
      TRUE ~ 1)) %>%
    uncount(replicates_needed) 
  names(df_pri_vec_sampling$lower[[1]]) <- paste0("sigma_f_plate[", 1:4, "]")
  names(df_pri_vec_sampling$lower[[2]]) <- paste0("f[", 1:4, "]")
  if (predict_f) {
    for (which_f_pred_var in 1:num_f_pred_vars) {
      names(df_pri_vec_sampling$lower[[2 + which_f_pred_var]]) <-
        paste0("sigma_f_pred_vars[", which_f_pred_var, ",", 1:4, "]")
    }
  }
  for (row in 1:nrow(df_pri_vec_sampling)) {
    names(df_pri_vec_sampling$upper[[row]]) <- names(df_pri_vec_sampling$lower[[row]])
  }
  df_pri_vec_sampling$param <- 1:nrow(df_pri_vec_sampling) # anything unique
  df_priors <- df_priors %>%
    bind_rows(full_join(df_pri_vec_sampling %>% 
                          select(param, lower) %>%
                          pivot_wider(names_from = param, values_from = lower) %>%
                          unnest_wider(everything()) %>%
                          pivot_longer(everything(), names_to = "param", values_to = "lower"),
                        df_pri_vec_sampling %>% 
                          select(param, upper) %>%
                          pivot_wider(names_from = param, values_from = upper) %>%
                          unnest_wider(everything()) %>%
                          pivot_longer(everything(), names_to = "param", values_to = "upper"),
                        by = "param"))
  df_priors <- df_priors %>%
    mutate(replicates_needed = case_when(
      param == "sigma_p_pos_pred_vars" ~ num_p_pos_pred_vars,
      TRUE ~ 1)) %>%
    uncount(replicates_needed) 
  if (predict_p_pos) {
    df_priors <- bind_rows(
      df_priors %>% 
        filter(param != "sigma_p_pos_pred_vars"),
      df_priors %>% 
        filter(param == "sigma_p_pos_pred_vars") %>%
        mutate(param = paste0(param, "[", row_number(), "]")))
  }
  
  # Sample params with simple uniform distributions.
  # df_ps = a df with prior samples. Short name due to heavy usage.
  N <- num_mc_iterations_prior
  df_ps <- tibble(sample = 1:N)
  for (row in 1:nrow(df_priors)) {
    df_ps[[df_priors$param[[row]]]] <-
      runif(N, df_priors$lower[[row]], df_priors$upper[[row]])
  }
  
  # Sample other params
  rho_samples <- rethinking::rlkjcorr(N, 4, eta = rho_prior_eta)
  for (i in 1:3) {
    for (j in seq(i+1, length.out = 4-i)) {
      df_ps[[paste0("rho[", i, ",", j, "]")]] <- 
        purrr::map_dbl(1:N, ~ rho_samples[.x, i, j])
    }
  }
  df_ps <- df_ps %>%
    mutate(y_obs_sd_cal_max = y_obs_sd_cal_min + y_obs_sd_cal_jump,
           y_obs_sd_sam_max = y_obs_sd_sam_min + y_obs_sd_sam_jump)
  x_sam_pos_mu_lower <- df_priors %>% filter(param == "x_sam_pos_mu") %>% pull(lower)
  x_sam_pos_mu_upper <- df_priors %>% filter(param == "x_sam_pos_mu") %>% pull(upper)
  df_ps <- df_ps %>%
    mutate(x_sam_pos_mu = runif(N,
                                pmax(x_sam_pos_mu_lower, x_sam_neg_mu),
                                x_sam_pos_mu_upper),
           density_type = "prior")
  
  if (predict_f) {
    # Independently draw num_f_pred_var_cats 4-vectors with mean zero, then scale,
    # mirroring these parts of the Stan code:
    # f_effects_by_pred_var_unscaled ~ multi_normal(zeros_for_f_pred_vars, rho);
    #...
    #int cat_current = 1;
    #for (f_pred_var in 1:num_f_pred_vars) {
    #  int num_cat_this_f_pred_var = num_cat_per_f_pred_var[f_pred_var];
    #  for (cat in cat_current:(cat_current + num_cat_this_f_pred_var - 1)) {
    #    f_effects_by_pred_var[cat, ] = f_effects_by_pred_var_unscaled[cat] .* 
    #      sigma_f_pred_vars[f_pred_var];
    #  }
    #  cat_current += num_cat_this_f_pred_var;
    #}
    f_effects_by_pred_var_samples <- 
      array(dim = c(num_mc_iterations_prior, num_f_pred_var_cats, 4))
    for (sample in 1:num_mc_iterations_prior) {
      rho_ <- matrix(
        c(1, df_ps$`rho[1,2]`[[sample]], df_ps$`rho[1,3]`[[sample]], df_ps$`rho[1,4]`[[sample]],
          df_ps$`rho[1,2]`[[sample]], 1, df_ps$`rho[2,3]`[[sample]], df_ps$`rho[2,4]`[[sample]],
          df_ps$`rho[1,3]`[[sample]], df_ps$`rho[2,3]`[[sample]], 1, df_ps$`rho[3,4]`[[sample]],
          df_ps$`rho[1,4]`[[sample]], df_ps$`rho[2,4]`[[sample]], df_ps$`rho[3,4]`[[sample]], 1),
        4, 4, byrow = TRUE)
      f_effects_by_pred_var_ <- rmvnorm(num_f_pred_var_cats, c(0, 0, 0, 0), rho_)
      cat_current <- 1
      for (f_pred_var_int in 1:num_f_pred_vars) {
        num_cat_this_f_pred_var = num_cat_per_f_pred_var[f_pred_var_int]
        sigma_f_pred_vars_ <- df_ps[sample, paste0(
          "sigma_f_pred_vars[", f_pred_var_int, ",", 1:4, "]")] %>% unlist
        for (cat in cat_current:(cat_current + num_cat_this_f_pred_var - 1)) {
          f_effects_by_pred_var_[cat, ] = f_effects_by_pred_var_[cat, ] * 
            sigma_f_pred_vars_
        }
        cat_current <- cat_current + num_cat_this_f_pred_var
      }
      f_effects_by_pred_var_samples[sample, , ] <- f_effects_by_pred_var_
    }
    for (f_pred_var_cat_int in 1:num_f_pred_var_cats) {
      for (which_f in 1:4) {
        df_ps[[paste0("f_effects_by_pred_var_cat[", f_pred_var_cat_int, ",",
                      which_f, "]" )]] <-
          f_effects_by_pred_var_samples[ , f_pred_var_cat_int, which_f]
      }
    }
  }
  
  if (predict_p_pos) {
    # Independently draw num_f_pred_var_cats 4-vectors with mean zero, then scale,
    # mirroring these parts of the Stan code:
    # p_pos_effects_by_pred_var_cat_unscaled ~ std_normal();
    #...
    #vector[tot_cat_per_p_pos_pred_var] p_pos_effects_by_pred_var_cat;
    #{
    #  int cat_current = 1;
    #  for (p_pos_pred_var in 1:num_p_pos_pred_vars) {
    #    int num_cat_this_p_pos_pred_var = num_cat_per_p_pos_pred_var[p_pos_pred_var];
    #    for (cat in cat_current:(cat_current + num_cat_this_p_pos_pred_var - 1)) {
    #      p_pos_effects_by_pred_var_cat[cat] = p_pos_effects_by_pred_var_cat_unscaled[cat] * 
    #        sigma_p_pos_pred_vars[p_pos_pred_var]; 
    #    }
    #    cat_current += num_cat_this_p_pos_pred_var;
    #  }
    #}
    p_pos_effects_by_pred_var_cat_samples <- # define unscaled, then scale
      matrix(rnorm(n = num_mc_iterations_prior * num_p_pos_pred_var_cats),
             nrow = num_mc_iterations_prior,
             ncol = num_p_pos_pred_var_cats)
    for (sample in 1:num_mc_iterations_prior) {
      cat_current <- 1
      for (p_pos_pred_var_int in 1:num_p_pos_pred_vars) {
        num_cat_this_p_pos_pred_var = num_cat_per_p_pos_pred_var[p_pos_pred_var_int]
        sigma_p_pos_pred_vars_ <- df_ps[[paste0(
          "sigma_p_pos_pred_vars[", p_pos_pred_var_int, "]")]][[sample]]
        cats <- cat_current:(cat_current + num_cat_this_p_pos_pred_var - 1)
        p_pos_effects_by_pred_var_cat_samples[sample, cats] <-
          p_pos_effects_by_pred_var_cat_samples[sample, cats] *
          sigma_p_pos_pred_vars_
        cat_current <- cat_current + num_cat_this_p_pos_pred_var
      }
    }
    for (p_pos_pred_var_cat_int in 1:num_p_pos_pred_var_cats) {
      df_ps[[paste0("p_pos_effects_by_pred_var_cat[", p_pos_pred_var_cat_int, "]" )]] <-
        p_pos_effects_by_pred_var_cat_samples[, p_pos_pred_var_cat_int]
    }
  }
  
} else {
  df_ps <- samples_prior %>%
    as.data.frame() %>%
    mutate(sample = row_number())
}
df_ps$density_type <- "prior"
setDT(df_ps)

# Merge prior and posterior samples
desired_cols <- names(df_ps)
df_fit_wide_postandprior <- rbind(df_fit_wide_postonly[,..desired_cols],
                                  df_ps) 

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
    "x_sam_neg_sd", x_sam_neg_sd,
    "x_sam_pos_sd", x_sam_pos_sd,
    "x_sam_neg_mu", x_sam_neg_mu,
    "x_sam_pos_mu", x_sam_pos_mu,
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
  
  if (predict_p_pos) {
    df_true_p_pos_effects_by_pred_var <- p_pos_effects_by_pred_var %>% 
      map(function(mat) {mat %>%
          as_tibble(.name_repair = "universal_quiet") %>%
          mutate(cat = names(mat))}) %>% 
      bind_rows(.id = "p_pos_pred_var") %>%
      mutate(p_pos_pred_var_cat = paste0(p_pos_pred_var, cat)) %>% 
      inner_join(df_p_pos_pred_vars_cats, by = "p_pos_pred_var_cat")
    stopifnot(identical(sort(df_true_p_pos_effects_by_pred_var$p_pos_pred_var_cat),
                        sort(design_matrix_p_pos_colnames_expected)))
    df_true_p_pos_effects_by_pred_var <- df_true_p_pos_effects_by_pred_var %>%
      mutate(param = paste0("p_pos_effects_by_pred_var_cat[", 
                            p_pos_pred_var_cat_int, "]")) %>%
      select(param, value)
    df_true_sigma_p_pos_pred_vars <-
      tibble(p_pos_pred_var = names(sigma_p_pos_pred_vars),
             value = sigma_p_pos_pred_vars) %>%
      inner_join(df_p_pos_pred_vars, by = "p_pos_pred_var") %>%
      mutate(param = paste0("sigma_p_pos_pred_vars[", p_pos_pred_var_int, "]")) %>%
      select(param, value)
    df_true_pop_params <- df_true_pop_params %>%
      bind_rows(df_true_p_pos_effects_by_pred_var,
                df_true_sigma_p_pos_pred_vars)
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
  if (predict_p_pos) {
    df_param_names <- df_param_names %>%
      tidyr::extract(orig, 
                     into = "p_pos_pred_var_int", 
                     regex = "sigma_p_pos_pred_vars\\[([0-9]+)\\]",
                     remove = FALSE) %>%
      tidyr::extract(orig, 
                     into = "p_pos_pred_var_cat_int", 
                     regex = "p_pos_effects_by_pred_var_cat\\[([0-9]+)\\]",
                     remove = FALSE) %>%
      mutate(p_pos_pred_var_int = as.integer(p_pos_pred_var_int),
             p_pos_pred_var_cat_int = as.integer(p_pos_pred_var_cat_int)) %>%
      left_join(df_p_pos_pred_vars, by = "p_pos_pred_var_int") %>%
      left_join(df_p_pos_pred_vars_cats, by = "p_pos_pred_var_cat_int") %>% 
      mutate(new = case_when(
        !is.na(p_pos_pred_var_int) ~ paste0("sigma_p_pos_pred_vars_", p_pos_pred_var),
        !is.na(p_pos_pred_var_cat_int) ~ paste0("p_pos_effect_", p_pos_pred_var_cat),
        TRUE ~ new
      ))
  }
  df_param_names$new
}
setnames(df_fit_wide_postonly, rename_params)
setnames(df_fit_wide_postandprior, rename_params)
setnames(df_ps, rename_params)
if (data_was_simulated) df_true_pop_params$param <- rename_params(df_true_pop_params$param)

# Define p_pos by group, from overall p_pos and p_pos effects
if (predict_p_pos) {
  for (p_pos_pred_var_cat in df_p_pos_pred_vars_cats$p_pos_pred_var_cat) {
    df_fit_wide_postandprior[[paste0("p_pos_for_", p_pos_pred_var_cat)]] <- mastiff::logistic(
      mastiff::logit(df_fit_wide_postandprior$p_pos) + df_fit_wide_postandprior[[paste0("p_pos_effect_", p_pos_pred_var_cat)]])
    df_ps[[paste0("p_pos_for_", p_pos_pred_var_cat)]] <- mastiff::logistic(
      mastiff::logit(df_ps$p_pos) + df_ps[[paste0("p_pos_effect_", p_pos_pred_var_cat)]])
  }
  if (data_was_simulated) {
    df_true_pop_params <- df_true_pop_params %>%
      bind_rows(p_pos_overall_by_pred_var %>% 
                  map(function(mat) {mat %>%
                      as_tibble(.name_repair = "universal_quiet") %>%
                      mutate(cat = names(mat))}) %>%
                  bind_rows(.id = "p_pos_pred_var") %>%
                  mutate(param = paste0("p_pos_for_", p_pos_pred_var, cat)) %>%
                  select(param, value))
  }
}

# PLOT STAN OUTPUT ----

# Plot prior vs posterior for all main params, except f_effects
p <- ggplot() +
  geom_histogram(data = df_fit_wide_postandprior %>%
                   select(!matches("f_effect")) %>%
                   select(!matches("p_pos_")) %>%
                   select(!matches("_jump")) %>%
                   pivot_longer(-c("sample", "density_type"), names_to = "param"),
                 aes(value, fill = density_type, y = after_stat(density)),
                 alpha = 0.6,
                 position = "identity",
                 bins = 50) +
  facet_wrap(~param, scales = "free", nrow = 5) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(expand = FALSE) +
  labs(fill = "",
       x = "param value",
       y = "probability density")
if (data_was_simulated) {
  p <- p + geom_vline(data = df_true_pop_params %>%
                        filter(!str_detect(param, "f_effect"),
                               !str_detect(param, "p_pos_")),
                      aes(xintercept = value))
}
p
ggsave("~/enable/enable_posteriors.pdf", height = 8, width = 9.4)

p <- ggplot() +
  geom_histogram(data = df_fit_wide_postandprior %>%
                   select(density_type, sample, matches("p_pos_for_")) %>%
                   pivot_longer(-c("sample", "density_type"), names_to = "param"),
                 aes(value, fill = density_type, y = after_stat(density)),
                 alpha = 0.6,
                 position = "identity",
                 bins = 50) +
  facet_wrap(~param, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(expand = FALSE) +
  labs(fill = "",
       x = "param value",
       y = "probability density")
if (data_was_simulated) {
  p <- p + geom_vline(data = df_true_pop_params %>%
                        filter(str_detect(param, "p_pos_for_")),
                      aes(xintercept = value))
}
p
ggsave("~/enable/enable_posteriors_p_pos_by_cluster.pdf", height = 20, width = 32)

p <- ggplot() +
  geom_histogram(data = df_fit_wide_postandprior %>%
                   select(density_type, sample, matches("p_pos_effect_")) %>%
                   pivot_longer(-c("sample", "density_type"), names_to = "param"),
                 aes(value, fill = density_type, y = after_stat(density)),
                 alpha = 0.6,
                 position = "identity",
                 bins = 50) +
  facet_wrap(~param, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(expand = FALSE) +
  labs(fill = "",
       x = "param value",
       y = "probability density")
if (data_was_simulated) {
  p <- p + geom_vline(data = df_true_pop_params %>%
                        filter(str_detect(param, "p_pos_effect_")),
                      aes(xintercept = value))
}
p

p <- ggplot() +
  geom_violin(data = df_fit_wide_postonly %>%
                select(matches("p_pos_effect_age_group_")) %>%
                pivot_longer(everything(), names_to = "param") %>%
                mutate(param = str_remove_all(param, "^p_pos_effect_age_group_"),
                       param = factor(param, levels = c(
                         "2-5y", "6-10y", "11-17y", "18-50y", ">50y"))),
              aes(param, value)) +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(expand = FALSE) +
  labs(fill = "",
       x = "Age group",
       y = "Increase in log odds of being positive")
p
ggsave("~/enable/enable_positivity_predictors_age.pdf", height = 4, width = 6)  

p <- ggplot() +
  geom_violin(data = df_fit_wide_postonly %>%
                select(matches("p_pos_effect_job_")) %>%
                pivot_longer(everything(), names_to = "param") %>%
                mutate(param = str_remove_all(param, "^p_pos_effect_job_")),
              aes(param, value)) +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(expand = FALSE) +
  labs(fill = "",
       x = "",
       y = "Increase in log odds of being positive")
p
ggsave("~/enable/enable_positivity_predictors_job.pdf", height = 4, width = 6)  


p <- ggplot() +
  geom_violin(data = df_fit_wide_postonly %>%
                select(matches("p_pos_effect_country_")) %>%
                pivot_longer(everything(), names_to = "param") %>%
                mutate(param = str_remove_all(param, "^p_pos_effect_country_")),
              aes(param, value)) +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(expand = FALSE) +
  labs(fill = "",
       x = "",
       y = "Increase in log odds of being positive")
p
ggsave("~/enable/enable_positivity_predictors_country.pdf", height = 4, width = 6)  

# Plot the posterior median & CI estimates for f_effects.
# NB almost all CIs from country and lot_id overlap zero - suggests don't need
# to model this variation. 
quantiles <- c(0.025, 0.5, 0.975)
df_fit_wide_postonly %>%
  select(sample, , contains("f_effect[")) %>%
  pivot_longer(-sample, names_to = "param") %>%
  group_by(param) %>%
  reframe(value = quantile(value, probs = quantiles),
          quantile = quantiles) %>%
  pivot_wider(names_from = quantile, names_prefix = "q_") %>%
  tidyr::extract(param, remove = FALSE,
                 into = "which_f", 
                 regex = "f_effect\\[([0-9]+)\\]") %>%
  tidyr::extract(param, 
                 into = "param", 
                 regex = "f_effect\\[[0-9]+\\]_(.+)") %>%
  ggplot() +
  geom_point(aes(param, q_0.5)) +
  geom_errorbar(aes(param, ymin = q_0.025, ymax = q_0.975)) +
  geom_hline(yintercept = 0) +
  facet_wrap(~which_f, ncol = 1, scales = "free_y") +
  labs(x = "Plate category",
       y = "Effect on 4 parameters of Ab<->OD logistic function") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0))
ggsave("~/enable/enable_f_effects_of_lab_and_lot.pdf", height = 8, width = 8)  


# Compare true and estimated sample x 
quantiles <- c(0.025, 0.5, 0.975)
df_sam_x <- df_fit_wide_postonly %>%
  select(sample, starts_with("x_sam[")) %>%
  pivot_longer(-sample, names_to = "param") %>%
  tidyr::extract(param, 
                 into = "id_sam", 
                 regex = "x_sam\\[([0-9]+)\\]") %>%
  mutate(id_sam = as.integer(id_sam)) %>%
  group_by(id_sam) %>%
  reframe(value = quantile(value, probs = quantiles),
          quantile = quantiles) %>%
  pivot_wider(names_from = quantile, names_prefix = "x_q_")
if (data_was_simulated) {
  df_sam_x %>%
    left_join(df_sam %>%
                select(id_sam, x) %>%
                distinct(),
              by = "id_sam") %>%
    filter(id_sam %% 15 == 0) %>%
    ggplot() +
    geom_errorbar(aes(x, ymin = x_q_0.025, ymax = x_q_0.975)) +
    geom_point(aes(x, x_q_0.5)) +
    geom_abline() +
    scale_x_log10(breaks = xs) +
    scale_y_log10(breaks = xs) +
    labs(x = "True Ab",
         y = "Estimated Ab")
  ggsave("~/foo_8.pdf", height = 3.3, width = 3.3)
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
    #geom_violin(aes(pos, prob_pos_q_0.5)) +
    #geom_sina(aes(pos, prob_pos_q_0.5)) +
    #geom_point(aes(jitter(as.numeric(pos)), prob_pos_q_0.5)) +
    geom_histogram(aes(prob_pos_q_0.5, fill = pos), #y = after_stat(density))
                   position = "identity",
                   alpha = 0.6,
                   bins = 30) +
    labs(y = "Number of samples",
         #x = "Estimated probability of being positive (posterior median)",
         x = "Probability sample is positive",
         fill = "Truth:") +
    coord_cartesian(expand = FALSE) +
    scale_x_continuous(limits = c(NA, NA))
  ggsave("~/foo_5.pdf", height = 2.7, width = 3.5)
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
       y = "Estimated probability of being positive") +
  scale_x_log10()

# Plot cal data by plate  
plate_ints_to_plot <- 1:4
ggplot() +
  geom_point(data = df_cal %>%
               filter(plate_int %in% plate_ints_to_plot),
             aes(x, y)) +
  facet_wrap(~plate) +
  scale_x_log10() +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
ggsave("~/foo_7.pdf", height = 3.3, width = 3.3)

# The posteriors for the 4PL function by plate
xlog_range <- seq(log(0.2), log(40), length.out = 50)
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
plate_ints_to_plot <- 1:12
p <- ggplot(df_4pl %>%
              filter(plate_int %in% plate_ints_to_plot) %>%
              filter(sample %% 10 == 0)) +
  geom_line(aes(x, y, group = sample), alpha = 0.1) +
  geom_point(data = df_cal %>%
               filter(plate_int %in% plate_ints_to_plot),
             aes(x, y), col = "blue") +
  facet_wrap(~plate, nrow = 3) +
  scale_x_log10() +
  labs(x = "x = Ab concentration",
       y = "y = OD") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
if (data_was_simulated) {
  p <- p + 
    geom_line(data = df_plate %>%
                filter(plate_int %in% plate_ints_to_plot) %>%
                expand_grid(xlog = xlog_range) %>%
                mutate(x = exp(xlog)) %>%
                mutate(y = PL4(xlog, f_1, f_2, f_3, f_4)),
              aes(x, y), col = "blue", linewidth = 1)
}
p
ggsave("~/enable/enable_calibrator_posterior_curves.pdf", height = 4.5, width = 6)


# Plot P(x | pos), P(x | neg), P(x), P(pos | x)
xs_plot <- seq(from = 0, to = 4,
               length.out = 500)
df_x_distributions <- df_fit_wide_postonly %>%
  filter(sample %% 10 == 0) %>%
  select("sample", "x_sam_pos_alpha", "x_sam_pos_beta",
         "x_sam_neg_alpha", "x_sam_neg_beta", "p_pos") %>%
  full_join(tibble(x = xs_plot),
            by = character()) %>%
  mutate(`P(x | pos)` = dgamma(x, shape = x_sam_pos_alpha, rate = x_sam_pos_beta),
         `P(x | neg)` = dgamma(x, shape = x_sam_neg_alpha, rate = x_sam_neg_beta),
         `P(x)` = p_pos * `P(x | pos)` + (1 - p_pos) * `P(x | neg)`,
         `P(pos | x)` = p_pos * `P(x | pos)` / `P(x)`) %>%
  select(sample, x, `P(x | pos)`, 
         `P(x | neg)`, `P(x)`, `P(pos | x)`) %>%
  pivot_longer(c("P(x | pos)", "P(x | neg)", "P(x)", "P(pos | x)"))
p <- ggplot() +
  geom_line(data = df_x_distributions,
            aes(x = x, y = value, group = sample), alpha = 0.15) +
  facet_wrap(vars(name), scales = "free_y", ncol = 1) +
  labs(x = "x = Ab concentration",
       y = "") +
  #scale_x_continuous(expand = c(0, 0), limits = c(NA, 2.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, NA))
if (data_was_simulated) {
  df_x_distributions_truth <-
    tibble(x = xs_plot,
           xlog = log(x),
           `P(x | pos)` = dgamma(x, shape = x_sam_pos_alpha, rate = x_sam_pos_beta),
           `P(x | neg)` = dgamma(x, shape = x_sam_neg_alpha, rate = x_sam_neg_beta),
           `P(x)` = p_pos * `P(x | pos)` + (1 - p_pos) * `P(x | neg)`,
           `P(pos | x)` = p_pos * `P(x | pos)` / `P(x)`) %>%
    pivot_longer(-c("x", "xlog"))
  p <- p +
    geom_line(data = df_x_distributions_truth,
              aes(x = x, y = value), colour = "blue") 
}
p

# Plot P(x | pos), P(x | neg), P(x), P(pos | x) again but now with logx
xlogs_plot <- log(10) * -90:60 / 30
df_xlog_distributions <- df_fit_wide_postonly %>%
  filter(sample %% 10 == 0) %>%
  select("sample", "x_sam_pos_mu", "x_sam_pos_sd",
         "x_sam_neg_mu", "x_sam_neg_sd", "p_pos") %>%
  full_join(tibble(xlog = xlogs_plot,
                   x = exp(xlog)),
            by = character()) %>%
  mutate(`P(xlog | pos)` = dnorm(xlog, mean = x_sam_pos_mu, sd = x_sam_pos_sd),
         `P(xlog | neg)` = dnorm(xlog, mean = x_sam_neg_mu, sd = x_sam_neg_sd),
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
  #scale_x_log10(expand = c(0, 0), limits = c(NA, NA)) +
  scale_x_continuous(expand = c(0, 0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, NA))
if (data_was_simulated) {
  df_x_distributions_truth <-
    tibble(xlog = xlogs_plot,
           x = exp(xlog),
           `P(xlog | pos)` = x * dgamma(x, shape = x_sam_pos_alpha, rate = x_sam_pos_beta),
           `P(xlog | neg)` = x * dgamma(x, shape = x_sam_neg_alpha, rate = x_sam_neg_beta),
           `P(xlog)` = p_pos * `P(xlog | pos)` + (1 - p_pos) * `P(xlog | neg)`,
           `P(pos | xlog)` = p_pos * `P(xlog | pos)` / `P(xlog)`) %>%
    pivot_longer(-c("x", "xlog"))
  p <- p +
    geom_line(data = df_x_distributions_truth,
              aes(x = xlog, y = value), colour = "blue") 
}
p
ggsave("~/enable/enable_parametric_prob_pos.pdf", height = 5, width = 6)



# Posterior retrodictive check
group_size <- 15
df_plot_group <- df_cal %>%
  select(plate, plate_int) %>%
  distinct() %>%
  mutate(plot_group = 1 + (row_number() - 1) %/% group_size) 
df_plot <- df_fit_wide_postonly %>% 
  filter(sample %% 10 == 0) %>%
  select(sample, starts_with("y_cal_sim[")) %>%
  pivot_longer(-c("sample"), names_to = "param") %>%
  rename(y_sim = value) %>%
  tidyr::extract(param, 
                 into = c("cal_rep"), 
                 regex = "y_cal_sim\\[([0-9]+)\\]") %>%
  mutate(cal_rep = as.integer(cal_rep)) %>%
  left_join(df_cal %>% select(plate, x) %>% mutate(cal_rep = row_number()),
            by = "cal_rep") %>%
  left_join(df_plot_group, by = "plate")
pdf("~/enable/enable_y_dependent_noise.pdf",
    width = 18.5, height = 10.5)
for (group in unique(df_plot_group$plot_group)) {
  cat("Now doing page ", group, " of ", max(df_plot_group$plot_group), "\n")
  df_plot %>%
    filter(plot_group == group) %>%
    mutate(x = round(x, digits = 1)) %>%
    {print(ggplot(.) +
             geom_violin(aes(x, y_sim, group = x)) +
             geom_point(data = df_cal %>% 
                          left_join(df_plot_group, by = "plate") %>%
                          filter(plot_group == group) %>%
                          mutate(x = round(x, digits = 1)),
                        aes(x, y), col = "blue") +
             facet_wrap(~plate, nrow = 3, scales = "free") +
             scale_x_discrete(drop = TRUE) +
             labs(x = "concentration",
                  y = "absorbance")
    )}
}
dev.off()

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
  # 
  #scale_x_continuous(limits = c(0, 4)) +
  #scale_x_log10(limits = c(1e-3, 100)) + 
  labs(x = "log_e(Ab concentration)",
       #y = "population distribution of point estimates (inverting the 4PL)") +
       y = "probability density") +
  #geom_line(data = df_gamma_, aes(x, p), col = "blue") +
  #geom_density(data = df_x_sam_point, aes(x), col = "blue") +
  NULL 
ggsave("~/enable/enable_PopDistributionOfX_LogScale.pdf", height = 4, width = 6.5)

df_x_sam_point <- df_fit_wide_postonly %>%
  select(starts_with("x_sam[")) %>%
  map(median) %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "param", values_to = "x") %>%
  mutate(x = exp(x)) %>%
  tidyr::extract(param, 
                 into = c("which_sam_rep"), 
                 regex = "x_sam\\[([0-9]+)\\]") %>%
  mutate(which_sam_rep = as.integer(which_sam_rep)) %>%
  left_join(df_sam %>% 
              select(sample_id, id_sam, which_sam_rep),
            by = "which_sam_rep")
df_x_sam_point %>%
  summarise(.by = c("id_sam", "sample_id"),
            x = if_else(all(is.infinite(x)),
                        Inf,
                        mean(x[is.finite(x)]))) %>%
  write_csv("~/enable/enable/enable_SampleX_empirical_v9.csv")

# Plot the posterior distribution of the population level distribution of 
# stochastically redrawn y_sam 
df_fit_wide_postonly %>%
  filter(sample %% 30 == 0) %>%
  select(sample, starts_with("y_sam_sim[")) %>%
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
ggsave("~/enable/enable_PopDistributionOfODs.pdf", height = 4, width = 6.5)








