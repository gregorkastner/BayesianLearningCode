#' Rounding with Trailing Zeros
#'
#' Returns vector of strings with the rounded input, keeping trailing zeros.
#'
#' @param x vector of numbers to be rounded.
#' @param digits number of digits after the decimal separator.
#'
#' @returns \code{roun} returns the character vector of rounded numbers.
#'
#' @name roun
NULL

#' @rdname roun
#' @export
roun <- function(x, digits = 0L) {
  sprintf(paste0("%.", digits, "f"), x)
}
