#' The Dirichlet Distribution
#'
#' Random generation for the Dirichlet distribution with parameter
#' vector alpha.
#'
#' @param alpha positive parameter vector.
#' @param ndraws number of observations.
#'
#' @returns \code{rdirichlet} generates random deviates.
#'
#' @name Dirichlet
NULL

#' @rdname Dirichlet
#' @export
#' @importFrom stats rgamma
rdirichlet <- function(ndraws, alpha) {
    x <- matrix(0, ndraws, length(alpha))
    for (j in seq_along(alpha)) {
        x[,j] <- rgamma(ndraws, shape = alpha[j])
    }
    x <- x / rowSums(x)
    x
}

