test_that("simulated data has expected size and indices", {
  num_simulations <- 100
  for (i in 1:num_simulations) {
    num_plate <- sample(1:3, 1)
    num_sam_per_plate <- sample(0:2, 1)
    num_cal_per_plate <- sample(0:2, 1)
    x_cals <- runif(n = num_cal_per_plate)
    num_rep_per_sam <- sample(0:2, 1)
    num_rep_per_cal <- sample(0:2, 1)
    data <- simulate_data(num_plate = num_plate, 
                          num_sam_per_plate = num_sam_per_plate, 
                          x_cals = x_cals,
                          num_rep_per_sam = num_rep_per_sam,
                          num_rep_per_cal = num_rep_per_cal,
                          x_mix_pred_vars = list(),
                          x_mix_pred_vars_sds = list(),
                          p_pos_binary_effects = numeric()
                          )
    expect_equal(data$df_sam$plate,
                 rep(1:num_plate, each = num_sam_per_plate * num_rep_per_sam))
    expect_equal(data$df_sam$id_sam,
                 as.character(rep(seq_len(num_sam_per_plate * num_plate), each = num_rep_per_sam)))
    expect_equal(data$df_cal$plate,
                 rep(1:num_plate, each = num_cal_per_plate * num_rep_per_cal))
    expect_equal(data$df_cal$x,
                 rep(rep(x_cals, each = num_rep_per_cal), times = num_plate))
  }
})

test_that("All combinations of x_mix parameters can have regression models", {
  x_mix_params <- c("mu_neg", "mu_pos", "p_pos", "sd_neg", "sd_pos")
  bools <- list()
  for (bool_mu_ in c(FALSE, TRUE)) {
      for (bool_sd_neg in c(FALSE, TRUE)) {
        for (bool_sd_pos in c(FALSE, TRUE)) {
          for (bool_p_pos in c(FALSE, TRUE)) {
            
            bools <- list(mu_neg = bool_mu_,
                          mu_pos = bool_mu_,
                          sd_neg = bool_sd_neg,
                          sd_pos = bool_sd_pos,
                          p_pos  = bool_p_pos)
            x_mix_pred_vars = list()
            x_mix_pred_vars_sds = list()
            for (param in x_mix_params) {
              if (bools[[param]]) {
                x_mix_pred_vars[[param]] <- list(a = letters[1:3],
                                                 b = LETTERS[1:4])
                x_mix_pred_vars_sds[[param]] <- c(a = 1, b = 1)
              }
            }
            
            data <- simulate_data(x_mix_pred_vars = x_mix_pred_vars,
                                  x_mix_pred_vars_sds = x_mix_pred_vars_sds)
            data_wrangled <- wrangle_data(df_sam = data$df_sam, 
                                          df_cal = data$df_cal, 
                                          priors_list = get_priors(),
                                          x_mix_pred_vars_names = 
                                            data$x_mix_pred_vars_names)
            expect_no_error(wrangle_true_params(data$params, 
                                                data_wrangled$data_descriptors))
            
            for (param in x_mix_params) {
              if (! bools[[param]]) {
                expect_identical(data$params$x_mix_effects[[param]], list())
              } else {
                expect_identical(purrr::map(data$params$x_mix_effects[[param]], 
                                            names),
                                 x_mix_pred_vars[[param]])
              }
            }
          }
        }
      }
    
  }
})

test_that("A regression model for only one of mu_neg or mu_pos gives warning", {
  expect_warning(simulate_data(
    x_mix_pred_vars = list(mu_pos = list(a = letters[1:3])),
    x_mix_pred_vars_sds = list(mu_pos = c(a = 1))),
    regexp = paste("You have chosen to simulate using a regression for only",
                   "one of mu_pos and mu_neg"))
  expect_warning(simulate_data(
    x_mix_pred_vars = list(mu_neg = list(a = letters[1:3])),
    x_mix_pred_vars_sds = list(mu_neg = c(a = 1))),
    regexp = paste("You have chosen to simulate using a regression for only",
                   "one of mu_pos and mu_neg"))
})

test_that("simulate_data works with f_pred_vars", {
  data <- simulate_data(num_plate = 100,
                        f_pred_vars = list(age = c("1", "2"),
                                           place = c("A", "B", "C")),
                        sigma_f_pred_vars = list(age = rep(1, 4), 
                                                 place = rep(1, 4)))
  expect_all_equal(sort(unique(data$df_plate$age)) == c("1", "2"), TRUE)
  expect_all_equal(sort(unique(data$df_plate$place)) == c("A", "B", "C"), TRUE)
})



test_that("All combinations of x_mix parameters work in x_mix_effects", {
  x_mix_params <- c("mu_neg", "mu_pos", "p_pos", "sd_neg", "sd_pos")
  for (bool_mu_ in c(FALSE, TRUE)) {
    for (bool_sd_neg in c(FALSE, TRUE)) {
      for (bool_sd_pos in c(FALSE, TRUE)) {
        for (bool_p_pos in c(FALSE, TRUE)) {
          
          bools <- c(mu_neg = bool_mu_,
                     mu_pos = bool_mu_,
                     sd_neg = bool_sd_neg,
                     sd_pos = bool_sd_pos,
                     p_pos  = bool_p_pos)
          if (! any(bools)) next
          x_mix_effects = list()
          x_mix_pred_vars <- list()
          for (param in x_mix_params) {
            if (bools[[param]]) {
              x_mix_effects[[param]] <- list(a = c(a = 1, b = 2, d = -3),
                                             b = c(e = 1, f = 2, g = -3))
              x_mix_pred_vars[[param]] <- list(a = c("a", "b", "d"),
                                               b = c("e", "f", "g"))
            }
          }
          
          data <- simulate_data(x_mix_effects = x_mix_effects, num_plate = 100)
          data_wrangled <- wrangle_data(df_sam = data$df_sam, 
                                        df_cal = data$df_cal, 
                                        priors_list = get_priors(),
                                        x_mix_pred_vars_names = 
                                          data$x_mix_pred_vars_names)
          expect_no_error(wrangle_true_params(data$params, 
                                              data_wrangled$data_descriptors))
          for (param in x_mix_params) {
            if (! bools[[param]]) {
              expect_identical(data$params$x_mix_effects[[param]], list())
            } else {
              expect_identical(purrr::map(data$params$x_mix_effects[[param]], 
                                          names),
                               x_mix_pred_vars[[param]])
            }
          }
        }
      }
    }
    
  }
})

test_that("Expected error messages in simulate_data are triggered", {
  
  expect_error(simulate_data(num_plate = 0), regexp = 
    "num_plate must be at least 1")
  
  expect_error(simulate_data(seed = "foo"), regexp = 
                 "is\\.numeric\\(x\\) is not TRUE")
  
  expect_error(simulate_data(p_pos_exact = TRUE,
                             x_mix_pred_vars = list(p_pos = list(a = letters[1:3])),
                             x_mix_pred_vars_sds = list(p_pos = c(a = 1))),
                             regexp = paste(
    "If you set p_pos_exact to TRUE, x_mix_pred_vars\\$p_pos must be an",
    "empty list"
  ))
  
  expect_error(simulate_data(
    x_mix_pred_vars = list(p_pos = list(age = 1:4)),
    x_mix_pred_vars_sds = list(p_pos = c(age = 1))), regexp = 
                 "x_mix_pred_vars\\$p_pos\\$age must be a character vector")
  
  expect_error(simulate_data(
    x_mix_pred_vars = list(p_pos = list(age = "foo")),
    x_mix_pred_vars_sds = list(p_pos = c(age = 1))), regexp = 
      "x_mix_pred_vars\\$p_pos\\$age must contain at least two elements")
  
  expect_error(simulate_data(num_plate = 1, num_sam_per_plate = 1,
    x_mix_pred_vars = list(p_pos = list(age = letters[1:4])),
    x_mix_pred_vars_sds = list(p_pos = c(age = 1))), regexp = paste(
      "x_mix_pred_vars\\$p_pos\\$age must not contain more elements than the",
      "total number of samples"))
  
  expect_error(simulate_data(
    x_mix_pred_vars = list(p_pos = list(age = c("A", "A"))),
    x_mix_pred_vars_sds = list(p_pos = c(age = 1))), regexp =
      "x_mix_pred_vars\\$p_pos\\$age must not contain duplicates")
  
  expect_error(simulate_data(
    x_mix_pred_vars = list(p_pos = list(age = letters[1:4])),
    x_mix_pred_vars_sds = list(p_pos = c(age = -1))), regexp =
      "x_mix_pred_vars_sds\\$p_pos\\$age must be a non-negative number")
  
  expect_error(simulate_data(
    x_mix_pred_vars = list(p_pos  = list(age = letters[1:4]),
                           mu_pos = list(age = letters[1:5])),
    x_mix_pred_vars_sds = list(p_pos = c(age = 1), mu_pos = c(age = 1))),
    regexp = paste("age was specified as a predictor variable for both mu_pos",
                   "and p_pos, but different categories were specified"))
  
  expect_error(simulate_data(
    x_mix_pred_vars = list(mu_pos = list(age = letters[1:4])),
    x_mix_pred_vars_sds = list(mu_pos = c(age = 1)),
    p_pos_binary_effects = c(age = 1)), regexp = paste(
      "p_pos_binary_pred_vars must not contain any predictor variables that",
      "are also used as \\(non-binary\\) predictors for "))
  
  expect_error(simulate_data(num_plate = 1, num_sam_per_plate = 1,
    x_mix_pred_vars = list(),
    x_mix_pred_vars_sds = list(),
    p_pos_binary_effects = c(age = 1)), regexp = 
      "At least 2 samples are needed to use p_pos_binary_effects")
  
  expect_error(simulate_data(num_plate = 1,
                             f_pred_vars = list(age = c("1", "2"),
                                   place = c("A", "B", "C")),
                sigma_f_pred_vars = list(age = rep(1, 4), 
                                         place = rep(1, 4))), regexp = 
                 "2\\+ plates are needed if plate covariates are used")
  
  expect_error(simulate_data(
    x_mix_pred_vars = list(foo = list(a = letters[1:4])),
    x_mix_pred_vars_sds = list(foo = 1)),
    regexp = paste("Unexpected variable name foo in the elements of the",
                   "x_mix_pred_vars list")) 
  
  expect_error(simulate_data(
    x_mix_pred_vars = list(
      p_pos = list(age = c("young", "old"), height = c("tall", "short"))),
    x_mix_pred_vars_sds = list(p_pos = c(age = 1))),
    regexp = paste("Different predictor variables were named in",
                              "x_mix_pred_vars\\$p_pos and in",
                              "x_mix_pred_vars_sds\\$p_pos:"))
  
  expect_error(simulate_data(p_pos_exact = TRUE,
                             p_pos_binary_effects = c(a = 1)),
               regexp = paste(
                 "If you set p_pos_exact to TRUE,",
                 "p_pos_binary_effects must be an empty numeric vector"
               ))
  
  expect_error(simulate_data(p_pos_exact = TRUE,
                             x_mix_pred_vars = list(p_pos = character(),
                                                    mu_pos = character(),
                                                    sd_pos = character(),
                                                    mu_neg = character(),
                                                    sd_neg = character()),
                             x_mix_pred_vars_sds = list(p_pos = numeric(),
                                                        mu_pos = numeric(),
                                                        sd_pos = numeric(),
                                                        mu_neg = numeric(),
                                                        sd_neg = numeric()),
                             p_pos_binary_effects = numeric(),
                             num_plate = 1,
                             num_sam_per_plate = 1,
                             p_pos = 0.5),
               regexp = paste(
                 "If you set p_pos_exact to TRUE, the number of samples",
                 "\\(i\\.e\\. num_plate \\* num_sam_per_plate\\) times p_pos",
                 "must be an integer"
               ))
  
})

test_that("p_pos_exact works", {
  data <- simulate_data(p_pos_exact = TRUE,
                x_mix_pred_vars = list(p_pos = character(),
                                       mu_pos = character(),
                                       sd_pos = character(),
                                       mu_neg = character(),
                                       sd_neg = character()),
                x_mix_pred_vars_sds = list(p_pos = numeric(),
                                           mu_pos = numeric(),
                                           sd_pos = numeric(),
                                           mu_neg = numeric(),
                                           sd_neg = numeric()),
                p_pos_binary_effects = numeric(),
                num_plate = 2,
                num_sam_per_plate = 1,
                p_pos = 0.5)
  expect_equal(sum( data$df_sam$pos), 2)
  expect_equal(sum(!data$df_sam$pos), 2)
})

test_that("wrangling true params works on simulated data", {
  
  priors_list <- get_priors()
  
  # Simulate lots of sams to ensure no categories are stochastically absent
  data <- simulate_data(num_plate = 100,
                        num_sam_per_plate = 100,
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
                        p_pos_binary_effects = c("boolA" = -2,
                                                 "boolB" = 0,
                                                 "boolC" = 2))
  
  data_wrangled <- wrangle_data(
  df_sam = data$df_sam, 
  df_cal = data$df_cal, 
  priors_list = priors_list,
  x_mix_pred_vars_names = data$x_mix_pred_vars_names, 
  p_pos_binary_pred_vars = data$p_pos_binary_pred_vars,
  f_pred_vars_names = data$f_pred_vars_names
  )  
  
  wrangled_names <- names(
  wrangle_true_params(param_true_values_list = data$params,
                      data_descriptors = data_wrangled$data_descriptors))
  expect_identical(sort(wrangled_names),
                         sort(c(
  "f[1]", "f[2]", "f[3]", "f[4]", "sigma_f[1]_plate", "sigma_f[2]_plate", "sigma_f[3]_plate", "sigma_f[4]_plate", "rho[1,2]", "rho[1,3]", "rho[1,4]", "rho[2,3]", "rho[2,4]", "rho[3,4]", "sd_neg", "sd_pos", "mu_neg", "mu_pos", "p_pos",  "p_blank",  "y_obs_sd_cal_min", "y_obs_sd_cal_jump", "y_obs_sd_sam_min", "y_obs_sd_sam_jump", "y_obs_sd_min_log_shift_sd", "y_obs_sd_jump_log_shift_sd", "p_pos_effect_lettera", "p_pos_effect_letterb", "p_pos_effect_letterc", "p_pos_effect_letterd", "sigma_p_pos_pred_vars_letter", "mu_neg_effect_lettera", "mu_neg_effect_letterb", "mu_neg_effect_letterc", "mu_neg_effect_letterd", "sigma_mu_neg_pred_vars_letter", "mu_pos_effect_lettera", "mu_pos_effect_letterb", "mu_pos_effect_letterc", "mu_pos_effect_letterd", "sigma_mu_pos_pred_vars_letter", "sd_neg_effect_lettera", "sd_neg_effect_letterb", "sd_neg_effect_letterc", "sd_neg_effect_letterd", "sigma_sd_neg_pred_vars_letter", "sd_pos_effect_lettera", "sd_pos_effect_letterb", "sd_pos_effect_letterc", "sd_pos_effect_letterd", "sigma_sd_pos_pred_vars_letter", "p_pos_effect_boolA", "p_pos_effect_boolB", "p_pos_effect_boolC"
  )))
})

test_that("Expected error messages in wrangle_true_params are triggered", {
  
  priors_list <- get_priors()
  data <- simulate_data()
  data_wrangled <- wrangle_data(
    df_sam = data$df_sam, 
    df_cal = data$df_cal, 
    priors_list = priors_list,
    x_mix_pred_vars_names = data$x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = data$p_pos_binary_pred_vars,
    f_pred_vars_names = data$f_pred_vars_names
  )  
  
  params_problem <- data$params
  params_problem$p_pos_binary_effects <- NULL
  expect_error(
    wrangle_true_params(param_true_values_list = params_problem,
                        data_descriptors = data_wrangled$data_descriptors),
    regexp = paste("param_true_values_list does not have an element named",
                   "p_pos_binary_effects"))
  
  data_descriptors_problem <- data_wrangled$data_descriptors
  data_descriptors_problem$lookup_pred_var_cat_int <- NULL
  expect_error(
    wrangle_true_params(param_true_values_list = data$params,
                        data_descriptors = data_descriptors_problem),
    regexp = paste("data_descriptors does not have an element named",
                   "lookup_pred_var_cat_int"))
  
  # TODO: all the others!

})
