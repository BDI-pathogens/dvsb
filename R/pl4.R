#' Four-parameter logistic function
#'
#' @param xlog a number
#' @param f_1 a number
#' @param f_2 a number
#' @param f_3 a number
#' @param f_4 a number
#'
#' @returns a number
#' 
#' @export
#' @keywords internal
PL4 <- function(xlog, f_1, f_2, f_3, f_4) {
  f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (xlog - f_4)))
}