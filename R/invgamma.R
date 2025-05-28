#' The Inverse Gamma Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the inverse gamma distribution with parameters a and b.
#'
#' For details, please see \code{\link[stats]{GammaDist}}.
#'
#' @param x vector of quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param a positive parameter.
#' @param b positive parameter.
#' @param ndraws number of observations.
#' @param log logical; if \code{TRUE}, densities are returned as logarithms.
#'
#' @returns \code{dinvgamma} gives the density, \code{pinvgamma} gives the
#' distribution function, \code{qinvgamma} gives the quantile function, and
#' \code{rinvgamma} generates random deviates.
#'
#' @name InvGammaDist
NULL

#' @rdname InvGammaDist
#' @export
#' @importFrom stats dgamma
dinvgamma <- function(x, a, b, log = FALSE) {
  logdens <- dgamma(1/x, a, b, log = TRUE) - 2 * log(x)
  if (log) logdens else exp(logdens)
}

#' @rdname InvGammaDist
#' @export
#' @importFrom stats pgamma
pinvgamma <- function(q, a, b) {
  1 - pgamma(1 / q, a, b)
}

#' @rdname InvGammaDist
#' @export
#' @importFrom stats qgamma
qinvgamma <- function(p, a, b) {
  1 / qgamma(1 - p, a, b)
}

#' @rdname InvGammaDist
#' @export
#' @importFrom stats rgamma
rinvgamma <- function(ndraws, a, b) {
  1 / rgamma(ndraws, a, b)
}
