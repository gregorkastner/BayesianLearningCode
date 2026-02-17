# Chapter 10: Bayesian Model Selection

## Section 9.1: The Foundations of Bayesian Model Selections

### Table 3.1: (Log) Bayes factor and model probabilities

``` r
logBF <- -7:7
BF <- exp(logBF)
PrM1 <- BF / (1 + BF)
PrM2 <- 1 - PrM1
#xtable::xtable(cbind(BF, logBF, PrM1, PrM2), digits = c(0, 3, 1, 3, 3))
knitr::kable(cbind(BF, logBF, PrM1, PrM2))
```

|           BF | logBF |      PrM1 |      PrM2 |
|-------------:|------:|----------:|----------:|
|    0.0009119 |    -7 | 0.0009111 | 0.9990889 |
|    0.0024788 |    -6 | 0.0024726 | 0.9975274 |
|    0.0067379 |    -5 | 0.0066929 | 0.9933071 |
|    0.0183156 |    -4 | 0.0179862 | 0.9820138 |
|    0.0497871 |    -3 | 0.0474259 | 0.9525741 |
|    0.1353353 |    -2 | 0.1192029 | 0.8807971 |
|    0.3678794 |    -1 | 0.2689414 | 0.7310586 |
|    1.0000000 |     0 | 0.5000000 | 0.5000000 |
|    2.7182818 |     1 | 0.7310586 | 0.2689414 |
|    7.3890561 |     2 | 0.8807971 | 0.1192029 |
|   20.0855369 |     3 | 0.9525741 | 0.0474259 |
|   54.5981500 |     4 | 0.9820138 | 0.0179862 |
|  148.4131591 |     5 | 0.9933071 | 0.0066929 |
|  403.4287935 |     6 | 0.9975274 | 0.0024726 |
| 1096.6331584 |     7 | 0.9990889 | 0.0009111 |

### Example 3.5: Testing for zero mean in the CHF-USD log returns

We load the data, compute the required parameters, and evaluate the two
marginal likelihoods.

``` r
data(exrates, package = "stochvol")
y <- 100 * diff(log(exrates$USD / exrates$CHF))

c0 <- 1
C0 <- 1
N0 <- 10^(-5:3)

N <- length(y)
cN <- c0 + N  / 2
CN_M1 <- C0 + sum(y^2) / 2
CN_M2 <- C0 + 0.5 * N * (mean((y - mean(y))^2) + N0 / (N0 + N) * mean(y)^2)

(logmarglikM1 <- lgamma(cN) + c0 * log(C0) -
   lgamma(c0) - cN * log(CN_M1) - 0.5 * N * log(2 * pi))
#> [1] -3456.478
(logmarglikM2 <- lgamma(cN) + c0 * log(C0) + 0.5 * log(N0) -
   lgamma(c0) - cN * log(CN_M2) - 0.5 * N * log(2 * pi) - 0.5 * log(N0 + N))
#> [1] -3465.344 -3464.192 -3463.041 -3461.890 -3460.739 -3459.588 -3458.440
#> [8] -3457.329 -3456.493
```

The “direct” formula for the log Bayes factor is even simpler, and we
can check whether we get the same answer.

``` r
(logBF <- 0.5 * (log(N0 + N) - log(N0)) + cN * (log(CN_M2) - log(CN_M1)))
#> [1] 8.86553606 7.71424355 6.56295141 5.41166293 4.26041101 3.10952463 1.96228321
#> [8] 0.85048010 0.01501194
all.equal(logBF, logmarglikM1 - logmarglikM2)
#> [1] TRUE
```
