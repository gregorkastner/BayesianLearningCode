#' The Non-standardized Student t Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the location-scale-transformed Student t distribution.
#'
#' The non-standardized Student t distribution with location \eqn{\mu}, scale
#' \eqn{\tau}, and df \eqn{\nu > 2} has mean \eqn{\mu} and variance
#' \eqn{\tau^2 \frac{\nu}{\nu-2}}.
#'
#' @param x vector of quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param location parameter.
#' @param scale positive parameter.
#' @param df positive parameter.
#' @param ndraws number of observations.
#' @param log logical; if \code{TRUE}, densities are returned as logarithms.
#'
#' @returns \code{dstudt} gives the density, \code{pstudt} gives the
#' distribution function, \code{qstudt} gives the quantile function, and
#' \code{rstudt} generates random deviates.
#'
#' @name StudTDist
NULL

#' @rdname StudTDist
#' @export
dstudt <- function(x, location = 0, scale = 1, df, log = FALSE) {
  logdens <- dt((x - location) / scale, df = df, log = TRUE) - log(scale)
  if (log) logdens else exp(logdens)
}

#' @rdname StudTDist
#' @export
pstudt <- function(x, location = 0, scale = 1, df) {
  pt((x - location) / scale, df = df)
}

#' @rdname StudTDist
#' @export
qstudt <- function(p, location = 0, scale = 1, df) {
  location + scale * qt(p, df = df)
}

#' @rdname StudTDist
#' @export
rstudt <- function(n, location = 0, scale = 1, df) {
  location + scale * rt(n, df = df)
}
