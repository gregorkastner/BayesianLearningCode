#' Accidents Data
#'
#' Seriously injured or killed pedestrians per month in Linz, Austria, from
#' January 1987 to December 2002.
#'
#' @format A monthly time series object with 192 rows and 4 columns
#' \describe{
#'   \item{children_accidents}{number of children aged 6 to 10 years seriously injured or
#'   killed}
#'   \item{children_exposure}{estimated number of children exposed}
#'   \item{seniors_accidents}{senior persons above 65 seriously injured or killed}
#'   \item{seniors_exposure}{estimated number of seniors exposed}
#' }
#' @source TBD
"accidents"

#' Eye Tracking Data
#'
#' Count data on eye tracking anomalies in 101 schizophrenic patients,
#' sorted by size
#'
#' @format A data frame with 101 observations of 1 variable
#' \describe{
#'   \item{anomalies}{number of eye tracking anomalies}
#' }
#'
#' @source Escobar & West (1998)
"eyetracking"

#' Movie Open Box Office Data
#'
#' TBD
#'
#' @format A data frame with 94 observations of 29 variables
#' \describe{
#'    \item{\code{OpenBoxOffice}}{a numeric vector}
#'    \item{\code{Action}}{a binary vector}
#'    \item{\code{Adventure}}{a binary vector}
#'    \item{\code{Animation}}{a binary vector}
#'    \item{\code{Comedy}}{a binary vector}
#'    \item{\code{Crime}}{a binary vector}
#'    \item{\code{Drama}}{a binary vector}
#'    \item{\code{Family}}{a binary vector}
#'    \item{\code{Fantasy}}{a binary vector}
#'    \item{\code{Mystery}}{a binary vector}
#'    \item{\code{Romance}}{a binary vector}
#'    \item{\samp{Sci-Fi}}{a binary vector}
#'    \item{\code{Thriller}}{a binary vector}
#'    \item{\code{PG}}{a binary vector}
#'    \item{\code{PG13}}{a binary vector}
#'    \item{\code{R}}{a binary vector}
#'    \item{\code{Budget}}{a numeric vector}
#'    \item{\code{Weeks}}{a numeric vector}
#'    \item{\code{Screens}}{a numeric vector}
#'    \item{\samp{S-21-27}}{a numeric vector}
#'    \item{\samp{S-14-20}}{a numeric vector}
#'    \item{\samp{S-7-13}}{a numeric vector}
#'    \item{\samp{S-4-6}}{a numeric vector}
#'    \item{\samp{S-1-3}}{a numeric vector}
#'    \item{\samp{Vol-21-27}}{a numeric vector}
#'    \item{\samp{Vol-14-20}}{a numeric vector}
#'    \item{\samp{Vol-7-13}}{a numeric vector}
#'    \item{\samp{Vol-4-6}}{a numeric vector}
#'    \item{\samp{Vol-1-3}}{a numeric vector}
#' }
#' @source TBD: \doi{10.7910/DVN/JS3AVM}
"movies"

#' Words Data
#'
#' List of the 1000 most common English headwords (i.e., the dictionary or
#' encyclopedia entries under which a set of related words appear) according to
#' the British National Corpus and the Corpus of Contemporary American English
#' (BNC/COCA).
#'
#' @format A data frame with 1000 observations of 2 variables
#' \describe{
#'   \item{word}{an English (head)word}
#'   \item{frequency}{number of occurrences}
#' }
#' @source Sourced from <https://www.eapfoundation.com/vocab/general/bnccoca/>
#' on March 21, 2025
"words"
