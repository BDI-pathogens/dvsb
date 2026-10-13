#' Four-parameter logistic function
#'
#' @param xlog a number
#' @param f_1 a number
#' @param f_2 a number
#' @param f_3 a number
#' @param f_4 a number
#'
#' @returns y (a number)
#' 
#' @export
#' @keywords internal
PL4 <- function(xlog, f_1, f_2, f_3, f_4) {
  f_2 + (f_3 - f_2) / (1 + exp(-f_1 * (xlog - f_4)))
}

#' Inverse of the four-parameter logistic function
#'
#' @param y a number
#' @param f_1 a number
#' @param f_2 a number
#' @param f_3 a number
#' @param f_4 a number
#'
#' @returns xlog (a number)
#' 
#' @export
#' @keywords internal
PL4_inv <- function(y, f_1, f_2, f_3, f_4) {
  stopifnot(y >= f_2)
  stopifnot(y <= f_3)
  f_4 - log((f_3 - f_2) / (y - f_2) - 1) / f_1
}