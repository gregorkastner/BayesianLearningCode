#' Accidents Data
#'
#' Seriously injured or killed pedestrians per month in Linz, Austria, from
#' January 1987 to December 2002.
#'
#' @format A monthly time series object with 192 rows and 4 columns
#' \describe{
#'   \item{children}{number of children aged 6 to 10 years seriously injured or killed}
#'   \item{seniors}{senior persons above 65 seriously injured or killed}
#'   \item{children_exposure}{estimated number of children exposed}
#'   \item{seniors_exposure}{estimated number of seniors exposed}
#' }
#' @source TBD
"accidents"

#' Words Data
#'
#' List of the 100 most common English lemmas (i.e., words whose different
#' forms have been grouped together), according to the Corpus of Contemporary
#' American English (COCA).
#'
#' @format A data frame with 100 observations of 2 variables
#' \describe{
#'   \item{word}{an English word (lemma)}
#'   \item{frequency}{number of occurrences}
#' }
#' @source Sourced from <https://www.wordfrequency.info> on March 20, 2025
"words"
