#' Reads posterior output files from cmdstan into a dataframe
#'
#' @param file_paths the paths of the cmdstan output files. If the file names
#'   end "gz" then we gunzip them (into a pipe, leaving the files themselves
#'   zipped).
#' @param params_to_ignore a character vector naming parameters to be excluded
#'   from output. If you name a tensor parameter (e.g. "my_vector"), we exclude
#'   all elements of it (e.g. "my_vector.1", "my_vector.2").
#' @param downsampling_factor a positive integer: the factor by which to
#'   downsample the posterior. e.g. if a value of 2 is specified, we keep 1 in
#'   every 2 samples. The default of 1 means we keep all samples.
#' @param comment_lines a single logical value: do the cmdstan output files
#'   contain comment lines (beginning with #)? Normally they do; set this
#'   argument to FALSE only if you have already removed such lines (for slightly
#'   faster parsing of the files).
#' @param verbose a single logical value: should we update on progress reading
#'   in the files?
#'
#' @returns a dataframe binding the contents of all the output files: one row
#'   per sample from the posterior, one column per parameter (and an extra
#'   column indicating which chain that sample came from, if that can be
#'   detected from the file names).
#' @export
#'
read_cmdstan_out_files <- function(file_paths,
                                   params_to_ignore = character(),
                                   downsampling_factor = 1L,
                                   comment_lines = TRUE,
                                   verbose = TRUE) {
  
  stopifnot(is.character(file_paths))
  stopifnot(length(file_paths) > 0)
  stopifnot(all(file.exists(file_paths)))
  stopifnot(is.character(params_to_ignore))
  mastiff::check_numeric(downsampling_factor, lower = 1, upper_inclusive = FALSE)

  purrr::map(file_paths, function(file_){
    if (verbose) {
      print(Sys.time())
      cat("Now reading file", file_, "\n")
    }
    if (comment_lines) {
      if (endsWith(file_, "gz")) {
        cmd <- paste("gunzip -c", file_, "| grep -v '^#'")
      } else {
        cmd <- paste("grep -v '^#'", file_)
      }
      df_ <- data.table::fread(cmd = cmd)
    } else {
      df_ <- data.table::fread(file_)
    }
    if (length(params_to_ignore)) {
      keep_col <- rep(TRUE, ncol(df_))
      for (param in params_to_ignore) {
        keep_based_on_this_param <- 
          colnames(df_) != param &
          ! startsWith(colnames(df_), paste0(param, ".")) 
        keep_col <- keep_col & keep_based_on_this_param
      }
      df_ <- df_[, ..keep_col]
      if (verbose) {
        cat("Reduced from", length(keep_col), "to", ncol(df_),
            "columns due to params_to_ignore\n")
      }
    }
    if (downsampling_factor > 1L) {
      df_ <- df_[seq(1, .N, by = downsampling_factor)] 
    }
    chain <- stringr::str_match(file_, "chain_([0-9+])")[,2]
    if (!is.na(chain)) {
      df_[, chain := as.integer(chain)]
    }
    df_
  }) %>% data.table::rbindlist(use.names = TRUE)

}