
# BayesianLearningCode

<!-- badges: start -->
<!-- badges: end -->

BayesianLearningCode provides R code to reproduce (most of the) examples
in a future book.

## Typical use case

Click on **Articles** on top of the page to view the R code.

## Installation

If for some reason (e.g., because you would like to have access to the data
sets) you'd like to have the entire package on your computer,
you can do so by installing the development version of BayesianLearningCode
from [GitHub](https://github.com/gregorkastner/BayesianLearningCode/) with:

``` r
install.packages("pak")
pak::pak("gregorkastner/BayesianLearningCode")
```

Another variant (for Sylvia's course)
``` r
install.packages(c("corrplot", "robustbase", "stochvol", "mvtnorm", "pgdraw", "coda", "numDeriv", "RColorBrewer", "gsl", "knitr", "rmarkdown"))
install.packages("https://github.com/gregorkastner/BayesianLearningCode/releases/download/v0.0.1/BayesianLearningCode_0.0.1.tar.gz", repos = NULL, type = "source")
```
