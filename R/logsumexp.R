#' The Log-sum-exp function and the Log-mean-exp function
#'
#' Numerically stable variant to compute the log of the sum (mean) of exponentials
#'
#' @param x sequence of real numbers
#'
#' @returns a single real number
#'
#' @name Logsumexp
NULL

#' @rdname Logsumexp
#' @export
logsumexp <- function(x) {
  maxx <- max(x)
  maxx + log(sum(exp(x - maxx)))
}

#' @rdname Logsumexp
#' @export
logmeanexp <- function(x) {
  logsumexp(x) - log(length(x))
}
