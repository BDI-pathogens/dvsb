read_cmdstan_out_files <- function(file_paths,
                                   params_to_ignore = character(),
                                   downsampling_factor = 1L,
                                   verbose = TRUE) {
  
  stopifnot(is.character(file_paths))
  stopifnot(length(file_paths) > 0)
  stopifnot(all(file.exists(file_paths)))
  stopifnot(is.character(params_to_ignore))
  
  df_fit_wide_postonly <- map(file_paths, function(file_){
    if (verbose) {
      print(Sys.time())
      cat("Now reading file", file_, "\n")
    }
    df_ <- data.table::fread(cmd = paste("grep -v '^#'", file_))
    if (length(params_to_ignore)) {
      keep_col <- rep(TRUE, ncol(df_))
      for (param in params_to_ignore) {
        keep_based_on_this_param <- 
          colnames(df_) != param &
          ! startsWith(colnames(df_), paste0(param, ".")) 
        keep_col <- keep_col & keep_based_on_this_param
      }
      df_ <- df_[, ..keep_col]
    }
    if (downsampling_factor > 1L) {
      df_ <- df_[seq(1, .N, by = downsampling_factor)] 
    }
    df_
  }) %>% data.table::rbindlist()
  
}