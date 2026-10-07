test_that("can get priors", {
  priors_list <- get_priors()
  expect_true(is.list(priors_list))
  expect_all_true(c("df_priors_scalars", "df_priors_vectors", "rho_prior_eta")
                  %in% names(priors_list))
  expect_true(is.data.frame(priors_list$df_priors_scalars))
  expect_true(is.data.frame(priors_list$df_priors_vectors))
  expect_true(is.numeric(priors_list$rho_prior_eta))
  expect_all_true(c("param", "lower", "upper")
                  %in% names(priors_list$df_priors_scalars))
  expect_all_true(c("param", "lower", "upper")
                  %in% names(priors_list$df_priors_vectors))
})

test_that("Expected error messages are triggered", {
  
  local({nonexistent_path <<- withr::local_tempfile(lines = c("x,y", "1,2"))})
  expect_error(get_priors(nonexistent_path), regexp =
                 paste("No file exists at the path", nonexistent_path))
  
  priors_df <- readr::read_csv(
    system.file("input_priors/priors.csv", package = "dvsb"))
  withr::with_tempfile("priors_file", fileext = ".csv", {
    
    priors_df_problem <- rbind(priors_df, list(
      parameter = "f[1]", min = 0, max = 1, explanation = ""))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp =
                   "The following parameters appeared in more than one row: f\\[1\\]")
    
    priors_df_problem <- priors_df %>%
      dplyr::filter(parameter != "f[1]")
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp =
                   "The following expected parameters were missing: f\\[1\\]")
    
    priors_df_problem <- rbind(priors_df, list(
      parameter = "foo", min = 0, max = 1, explanation = ""))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp =
                   "The following parameters were unexpected: foo")
    
    priors_df_problem <- priors_df %>%
      dplyr::mutate(min = dplyr::if_else(parameter == "f[1]", NA_real_, min))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp =
                   "The following parameters had missing \\(NA\\) min: f\\[1\\]")
    
    priors_df_problem <- priors_df %>%
      dplyr::mutate(max = dplyr::if_else(parameter == "f[1]", NA_real_, max))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp =
                   "The following parameters had missing \\(NA\\) max: f\\[1\\]")
    
    priors_df_problem <- priors_df %>%
      dplyr::mutate(min = dplyr::if_else(parameter == "f[1]", Inf, min))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp =
                   "The following parameters had a max <= min: f\\[1\\]")
    
    priors_df_problem <- priors_df %>%
      dplyr::mutate(max = dplyr::if_else(parameter == "rho_prior_eta", 2e5, max))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp =
                   "For parameter rho_prior_eta only, min should equal max")
    
    priors_df_problem <- priors_df %>%
      dplyr::mutate(min = dplyr::if_else(parameter == "sd_neg", -Inf, min))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp = paste(
      "The following parameters, which can never be negative, had",
      "a negative min value: sd_neg"))
    
    priors_df_problem <- priors_df %>%
      dplyr::mutate(max = dplyr::if_else(parameter == "p_pos", 2, max))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp = paste(
                   "The following parameters, which must always be less than",
                   "or equal to 1, had a max value greater than 1: p_pos"))
    
    priors_df_problem <- priors_df %>%
      dplyr::mutate(max = dplyr::if_else(parameter == "mu_neg", 2e5, max))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp =
                   "mu_pos max must be greater than mu_neg max")
    
    priors_df_problem <- priors_df %>%
      dplyr::mutate(min = dplyr::if_else(parameter == "mu_neg", 1e5, min),
                    max = dplyr::if_else(parameter == "mu_neg", 2e5, max))
    readr::write_csv(priors_df_problem, priors_file)
    expect_error(get_priors(priors_file), regexp =
                   "mu_pos min must be greater than mu_neg min")
    
    
    
  })
})

