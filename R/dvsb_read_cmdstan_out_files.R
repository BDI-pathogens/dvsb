read_cmdstan_out_files <- function(file_paths,
                                   params_to_ignore = character(),
                                   downsampling_factor = 1L,
                                   verbose = TRUE) {
  
  stopifnot(is.character(file_paths))
  stopifnot(length(file_paths) > 0)
  stopifnot(all(file.exists(file_paths)))
  stopifnot(is.character(params_to_ignore))
  mastiff::check_numeric(downsampling_factor, lower = 1, upper_inclusive = FALSE)

  map(file_paths, function(file_){
    if (verbose) {
      print(Sys.time())
      cat("Now reading file", file_, "\n")
    }
    if (endsWith(file_, "gz")) {
      cmd <- paste("gunzip -c", file_, "| grep -v '^#'")
    } else {
      cmd <- paste("grep -v '^#'", file_)
    }
    df_ <- data.table::fread(cmd = cmd)
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
  }) %>% data.table::rbindlist(use.names = TRUE)

}