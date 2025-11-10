#' Inflation Data
#'
#' Harmonized Index of Consumer Prices (HICP) in the Euro area (changing
#' composition). Monthly data from January 1997 to June 2025.
#'
#' @format A sequence of length 342.
#' @source European Central Bank Data Portal, retrieved from
#' \url{https://data.ecb.europa.eu/data/datasets/ICP/ICP.M.U2.N.000000.4.ANR},
#' July 25, 2025.
"inflation"

#' GDP data
#'
#' Gross domestic product of the USA. Quarterly data from 1947Q1 to 2025Q1.
#'
#' Gross domestic product (GDP), the featured measure of U.S. output,
#' is the market value of the goods and services produced by labor and property
#' located in the United States. The units are Billions of Dollars,
#' Seasonally Adjusted Annual Rate.
#'
#' @format A sequence of length 313
#' @source U.S. Bureau of Economic Analysis, Gross Domestic Product (GDP),
#' retrieved from FRED, Federal Reserve Bank of St. Louis;
#' \url{https://fred.stlouisfed.org/series/GDP}, July 24, 2025.
"gdp"

#' Labor Market Data
#'
#' The Austrian social security authority collects detailed data for
#' all workers. This data set is a random sample of 4376 persons from
#' the birth cohort 1921-1980 with their work status recorded for the
#' years 1986-1998.
#'
#' @format A data frame with 4376 observations (persons) on 41 variables:
#' \describe{
#'   \item{female}{logical, TRUE if the person is female}
#'   \item{birthyear}{integer, year of birth}
#'   \item{income_1986, income_1987, income_1988, income_1989,
#'   income_1990, income_1991, income_1992, income_1993, income_1994,
#'   income_1995, income_1996, income_1997, income_1998}{ordered
#'   factor, income status based on the gross monthly wage on May 31st
#'   of the specific year with levels `zero` (no income) and `<= 1st
#'   quintile`, `<= 2nd quintile`, `<= 3rd quintile`, `<= 4th
#'   quintile` and `<= 5th quintile` according to the yearly wage
#'   distribution}
#'   \item{employerchanges_1986, employerchanges_1987,
#'   employerchanges_1988, employerchanges_1989, employerchanges_1990,
#'   employerchanges_1991, employerchanges_1992, employerchanges_1993,
#'   employerchanges_1994, employerchanges_1995, employerchanges_1996,
#'   employerchanges_1997, employerchanges_1998}{integer, number of
#'   times the person changed employer in the respective year}
#'   \item{wcollar_1986, wcollar_1987,
#'   wcollar_1988, wcollar_1989, wcollar_1990,
#'   wcollar_1991, wcollar_1992, wcollar_1993,
#'   wcollar_1994, wcollar_1995, wcollar_1996,
#'   wcollar_1997, wcollar_1998}{logical, TRUE if the person is
#'   classified as white collar employee and not a blue collar
#'   employee in the respective year}
#' }
#'
#' @source Zweimüller, J., R. Winter-Ebmer, R. Lalive, A. Kuhn,
#' J.-P. Wuellrich, O. Ruf, and S. Büchi (2009): The Austrian
#' Social Security Database (ASSD). Working Paper 0903, NRN: The
#' Austrian Center for Labor Economics and the Analysis of the
#' Welfare State, Linz, Austria.
"labor"

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
