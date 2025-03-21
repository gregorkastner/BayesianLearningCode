#' Accidents Data
#'
#' Seriously injured or killed pedestrians per month in Linz, Austria, from
#' January 1987 to December 2002.
#'
#' @format A monthly time series object with 192 rows and 4 columns
#' \describe{
#'   \item{children}{number of children aged 6 to 10 years seriously injured or
#'   killed}
#'   \item{seniors}{senior persons above 65 seriously injured or killed}
#'   \item{children_exposure}{estimated number of children exposed}
#'   \item{seniors_exposure}{estimated number of seniors exposed}
#' }
#' @source TBD
"accidents"

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
