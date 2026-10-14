test_that("Expected error messages in simulate_data are triggered", {
  
  priors_list <- get_priors()
  data <- simulate_data(x_mix_pred_vars = list(p_pos = list(age = letters[1:4])),
                        x_mix_pred_vars_sds = list(p_pos = c(age = 1)))
  df_sam <- data$df_sam 
  df_cal <- data$df_cal 
  x_mix_pred_vars_names <- data$x_mix_pred_vars_names 
  p_pos_binary_pred_vars <- data$p_pos_binary_pred_vars
  f_pred_vars_names <- data$f_pred_vars_names
  
  df_sam_problem <- df_sam
  df_sam_problem$y <- NULL
  expect_error(wrangle_data(
    df_sam = df_sam_problem, 
    df_cal = df_cal, 
    priors_list = priors_list,
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = "Column y missing from df_sam")
  
  df_cal_problem <- df_cal
  df_cal_problem$y <- NULL
  expect_error(wrangle_data(
    df_sam = df_sam, 
    df_cal = df_cal_problem, 
    priors_list = priors_list,
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = "Column y missing from df_cal")
  
  df_sam_problem <- df_sam
  df_sam_problem$plate_int <- NA
  expect_error(wrangle_data(
    df_sam = df_sam_problem, 
    df_cal = df_cal, 
    priors_list = priors_list, 
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = paste("Column plate_int already present in df_sam; please remove",
  "it to prevent overwriting"))
  
  df_cal_problem <- df_cal
  df_cal_problem$plate_int <- NA
  expect_error(wrangle_data(
    df_sam = df_sam, 
    df_cal = df_cal_problem, 
    priors_list = priors_list, 
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = paste("Column plate_int already present in df_cal; please remove",
                    "it to prevent overwriting"))
  
  priors_list_problem <- priors_list
  priors_list_problem$df_priors_scalars$param <- NULL
  expect_error(wrangle_data(
    df_sam = df_sam, 
    df_cal = df_cal, 
    priors_list = priors_list_problem,
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = "Column param missing from priors_list\\$df_priors_scalars")
  
  priors_list_problem <- priors_list
  priors_list_problem$df_priors_vectors$param <- NULL
  expect_error(wrangle_data(
    df_sam = df_sam, 
    df_cal = df_cal, 
    priors_list = priors_list_problem,
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = "Column param missing from priors_list\\$df_priors_vectors")
  
  priors_list_problem <- priors_list
  priors_list_problem$df_priors_scalars <-
    priors_list_problem$df_priors_scalars %>%
    dplyr::filter(param != "p_pos")
  expect_error(wrangle_data(
    df_sam = df_sam, 
    df_cal = df_cal, 
    priors_list = priors_list_problem,
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = "Parameter p_pos missing from priors_list\\$df_priors_scalars")
  
  priors_list_problem <- priors_list
  priors_list_problem$df_priors_vectors <-
    priors_list_problem$df_priors_vectors %>%
    dplyr::filter(param != "sigma_f_plate")
  expect_error(wrangle_data(
    df_sam = df_sam, 
    df_cal = df_cal, 
    priors_list = priors_list_problem,
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = "Parameter sigma_f_plate missing from priors_list\\$df_priors_vectors")

  df_sam_problem <- df_sam
  df_sam_problem$age <- "a"
  expect_error(wrangle_data(
    df_sam = df_sam_problem, 
    df_cal = df_cal, 
    priors_list = priors_list,
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = paste("Fewer than 2 unique values of age found in df_sam; need",
                    "at least 2 to use it as an x_mix_pred_var"))
  
  df_sam_problem <- df_sam
  df_sam_problem$age[[1]] <- NA
  expect_error(wrangle_data(
    df_sam = df_sam_problem, 
    df_cal = df_cal, 
    priors_list = priors_list,
    x_mix_pred_vars_names = x_mix_pred_vars_names, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = paste("Some samples have missing age"))
  
  x_mix_pred_vars_names_problem <- x_mix_pred_vars_names
  x_mix_pred_vars_names_problem$mu_pos <- "age"
  expect_error(wrangle_data(
    df_sam = df_sam, 
    df_cal = df_cal, 
    priors_list = priors_list,
    x_mix_pred_vars_names = x_mix_pred_vars_names_problem, 
    p_pos_binary_pred_vars = p_pos_binary_pred_vars,
    f_pred_vars_names = f_pred_vars_names
  ), regexp = paste("You cannot specify a regression for mu_pos but not",
                    "mu_neg, or vice versa: it's both or neither"))  
    
})
