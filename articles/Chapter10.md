# Chapter 10: Bayesian Model Selection

## Section 10.1: The Foundations of Bayesian Model Selections

### Table 10.1: (Log) Bayes factor and model probabilities

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

## Section 10.2: Bayesian Testing of Hypotheses

### Example 10.5: Testing for zero mean in the CHF-USD log returns

Before producing Figure 10.1, we begin with a concrete example. We load
the data, compute the required parameters, and evaluate the two marginal
likelihoods.

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

Now we are ready to investigate the sensitivity of the log Bayes factor
with respect to prior hyperparameter choices.

``` r
C0 <- c(0, 0, 0.001, 1, 100)
N0 <- 10^seq(1, 10, by = 0.1)

logBF <- matrix(NA_real_, nrow = length(N0), ncol = length(c0))
for (i in seq_along(c0)) {
  cN <- c0[i] + N  / 2
  CN_M1 <- C0[i] + sum(y^2) / 2
  CN_M2 <- C0[i] + 0.5 * N * (mean((y - mean(y))^2) + N0 / (N0 + N) * mean(y)^2)
  logBF[, i] <- 0.5 * (log(N0 + N) - log(N0)) + cN * (log(CN_M2) - log(CN_M1))
}
matplot(N0, logBF, log = "x", type = "l",
        xlab = expression(paste(N[0], " [log scale]")), ylab = "log BF")
abline(h = 0, lty = 3)
legend("topright", paste0("IG(", c0, ",", C0, ")"), col = seq_along(c0),
       lty = seq_along(c0))
title("Log Bayes factor in favor of the zero-mean model")
```

![](Chapter10_files/figure-html/unnamed-chunk-5-1.png)

### Example 10.6: Testing for heterogeneity of the no-income risk

First, we re-load the data from Chapter 3.

``` r
data("labor", package = "BayesianLearningCode")
labor <- subset(labor,
                income_1997 != "zero" & female,
                c(income_1998, wcollar_1986))
labor <- with(labor,
              data.frame(unemployed = income_1998 == "zero",
                         wcollar = wcollar_1986))
```

Next, we compute the marginal likelihoods for the homogeneity model.

``` r
N <- length(labor$unemployed)
SN <- sum(labor$unemployed)
hN <- SN / N
N0 <- c(2, 10)
a0 <- c(1, 0.5, hN * N0) 
b0 <- c(1, 0.5, N0 * (1 - hN))

aN <- a0 + sum(labor$unemployed)
bN <- b0 + N - sum(labor$unemployed)

logmarglikM1 <- lbeta(aN, bN) - lbeta(a0, b0)
```

And for the heterogeneity model.

``` r
N1 <- with(labor, sum(wcollar))
SN1 <- with(labor, sum(wcollar & unemployed))
N2 <- with(labor, sum(!wcollar))
SN2 <- with(labor, sum(!wcollar & unemployed))
hN <- (SN1 + SN2) / (N1 + N2) # redundant, same as above
a01 <- a02 <- c(1, 0.5, hN * N0) # redundant, same as a0
b01 <- b02 <- c(1, 0.5, N0 * (1 - hN)) # redundant, same as b0

aN1 <- a01 + SN1
bN1 <- b01 + N1 - SN1
aN2 <- a02 + SN2
bN2 <- b02 + N2 - SN2

logmarglikM2 <- lbeta(aN1, bN1) + lbeta(aN2, bN2) -
   lbeta(a01, b01) - lbeta(a02, b02)
```

We can now compute the model probabilities and the log BFs.

``` r
logBF <- logmarglikM1 - logmarglikM2
PrM2 <- 0.5 / (0.5 + 0.5 * exp(logBF))
PrM1 <- 1 - PrM2
#print(xtable::xtable(res, digits = c(0, 2, 2, 5, 5, 5), row.names = FALSE))
knitr::kable(res <- cbind(logmarglikM1, logmarglikM2, logBF, PrM1, PrM2))
```

| logmarglikM1 | logmarglikM2 |     logBF |      PrM1 |      PrM2 |
|-------------:|-------------:|----------:|----------:|----------:|
|    -387.5820 |    -385.9187 | -1.663244 | 0.1593270 | 0.8406730 |
|    -387.5560 |    -385.8867 | -1.669302 | 0.1585173 | 0.8414827 |
|    -387.2933 |    -385.3980 | -1.895305 | 0.1306408 | 0.8693592 |
|    -386.2587 |    -383.4054 | -2.853317 | 0.0545101 | 0.9454899 |

### Example 10.7: Testing for heterogeneity of the mortality rate

We proceed exactly as above, just with different data.

``` r
N1  <- 1668 # number at risk in City A
SN1 <- 2    # cancer deaths in City A
N2  <- 583  # number at risk in City B
SN2 <- 1    # cancer deaths in City B
```

We again begin with the homogeneity model.

``` r
N <- N1 + N2
SN <- SN1 + SN2
hN <- SN / N
N0 <- c(2, 10)
a0 <- c(1, 0.5, hN * N0) 
b0 <- c(1, 0.5, N0 * (1 - hN))

aN <- a0 + SN
bN <- b0 + N - SN

logmarglikM1 <- lbeta(aN, bN) - lbeta(a0, b0)
```

Followed by the heterogeneity model.

``` r
a01 <- a02 <- a0
b01 <- b02 <- b0

aN1 <- a01 + SN1
bN1 <- b01 + N1 - SN1
aN2 <- a02 + SN2
bN2 <- b02 + N2 - SN2

logmarglikM2 <- lbeta(aN1, bN1) + lbeta(aN2, bN2) -
   lbeta(a01, b01) - lbeta(a02, b02)
```

We can now compute the model probabilities and the log BFs.

``` r
logBF <- logmarglikM1 - logmarglikM2
PrM2 <- 0.5 / (0.5 + 0.5 * exp(logBF))
PrM1 <- 1 - PrM2
#print(xtable::xtable(res, digits = c(0, 2, 2, 3, 3, 3), row.names = FALSE))
knitr::kable(res <- cbind(logmarglikM1, logmarglikM2, logBF, PrM1, PrM2))
```

| logmarglikM1 | logmarglikM2 |    logBF |      PrM1 |      PrM2 |
|-------------:|-------------:|---------:|----------:|----------:|
|    -29.08387 |    -34.30308 | 5.219212 | 0.9946175 | 0.0053825 |
|    -26.95877 |    -30.22452 | 3.265757 | 0.9632352 | 0.0367648 |
|    -28.40707 |    -33.09584 | 4.688776 | 0.9908859 | 0.0091141 |
|    -26.84585 |    -29.97906 | 3.133212 | 0.9582421 | 0.0417579 |

### Example 10.8: Testing for a structural break in the road safety data

We load the data.

``` r
data("accidents", package = "BayesianLearningCode")
```

We define the priors hyperparameters and compute the log marginal
likelihood for the exchangable model.

``` r
y <- accidents[, "children_accidents"]
a0_tmp <- c(0.01, 0.1, 0.5, 1, 2)
m0_tmp <- c(1, mean(y), 5, 7, 10)
grid <- expand.grid(a0 = a0_tmp, m0 = m0_tmp)
a0 <- grid$a0
b0 <- grid$a0 / grid$m0
aN <- sum(y) + a0
bN <- length(y) + b0
logmarglikM1 <- matrix(a0 * log(b0) + lgamma(aN) -
                       aN * log(bN) - lgamma(a0) -
                       sum(lgamma(y + 1)),
                       nrow = length(a0_tmp), ncol = length(m0_tmp),
                       dimnames = list(a0 = a0_tmp, m0 = round(m0_tmp, 2)))
knitr::kable(round(logmarglikM1, 2))
```

|      |       1 |    1.84 |       5 |       7 |      10 |
|:-----|--------:|--------:|--------:|--------:|--------:|
| 0.01 | -320.55 | -320.55 | -320.55 | -320.55 | -320.56 |
| 0.1  | -318.50 | -318.47 | -318.51 | -318.53 | -318.56 |
| 0.5  | -317.43 | -317.31 | -317.49 | -317.61 | -317.75 |
| 1    | -317.12 | -316.89 | -317.26 | -317.49 | -317.77 |
| 2    | -316.96 | -316.51 | -317.24 | -317.70 | -318.26 |

And now the same for the model with structural break.

``` r
accidents1 <- window(accidents, end = c(1994, 9))
accidents2 <- window(accidents, start = c(1994, 10))

a01 <- a02 <- a0
b01 <- b02 <- b0

aN1 <- a01 + sum(accidents1[, "children_accidents"])
aN2 <- a02 + sum(accidents2[, "children_accidents"])
bN1 <- b01 + length(accidents1[, "children_accidents"])
bN2 <- b02 + length(accidents2[, "children_accidents"])

logmarglikM2 <- matrix(a01 * log(b01) + lgamma(aN1) +
                       a02 * log(b02) + lgamma(aN2) -
                       aN1 * log(bN1) - lgamma(a01) -
                       aN2 * log(bN2) - lgamma(a02) -
                       sum(lgamma(y + 1)),
                       nrow = length(a0_tmp), ncol = length(m0_tmp),
                       dimnames = list(a0 = a0_tmp, m0 = round(m0_tmp, 2)))
knitr::kable(round(logmarglikM2, 2))
```

|      |       1 |    1.84 |       5 |       7 |      10 |
|:-----|--------:|--------:|--------:|--------:|--------:|
| 0.01 | -321.40 | -321.40 | -321.40 | -321.41 | -321.41 |
| 0.1  | -317.30 | -317.25 | -317.33 | -317.37 | -317.43 |
| 0.5  | -315.17 | -314.94 | -315.31 | -315.54 | -315.81 |
| 1    | -314.58 | -314.12 | -314.85 | -315.31 | -315.86 |
| 2    | -314.30 | -313.38 | -314.83 | -315.75 | -316.86 |

We can now compute log Bayes factors and corresponding model
probabilities (under uniform prior probabilities).

``` r
logBF <- logmarglikM1[, 2] - logmarglikM2[, 2]
PrM2 <- 0.5 / (0.5 + 0.5 * exp(logBF))
PrM1 <- 1 - PrM2
knitr::kable(round(rbind(PrM1, PrM2), 3))
```

|      |  0.01 |   0.1 |   0.5 |     1 |     2 |
|:-----|------:|------:|------:|------:|------:|
| PrM1 | 0.701 | 0.228 | 0.086 | 0.059 | 0.042 |
| PrM2 | 0.299 | 0.772 | 0.914 | 0.941 | 0.958 |

### Example 10.9: Testing normal versus t

After loading the data, we define prior hyperparameters and compute log
marginal likelihoods. This is straightforward for the normal model.

``` r
nu <- 7
y <- 100 * diff(log(exrates$USD / exrates$CHF))
N <- length(y)

c10 <- c30 <- 1
C10 <- 2
C30 <- (nu - 2) / nu * C10

c1N <- c10 + N / 2
C1N <- C10 + sum(y^2) / 2

logmarglikM1 <- lgamma(c1N) + c10 * log(C10) -
   lgamma(c10) - c1N * log(C1N) - 0.5 * N * log(2 * pi)
```

For the Student $t$ model, we need, e.g., numerical integration. Note
that we apply the “log-sum-exp” trick here to normalize the integrand so
that its maximum is 1.

``` r
integrand_nonvec <- function(sigma2, y, c0, C0, nu, const = 0, log = FALSE) {
  N <- length(y)
  logint <- -(N/2 + c0 + 1) * log(sigma2) -
            0.5 * (nu + 1) * sum(log(1 + y^2 / (nu * sigma2))) -
            C0 / sigma2
  if (log) logint + const else exp(logint + const)
}
integrand <- Vectorize(integrand_nonvec, "sigma2")

resolution <- 1000
grid <- seq(0.1, 0.6, length.out = resolution + 1)

tmp <- integrand(grid, y, c30, C30, nu, log = TRUE)
const <- -max(tmp)
tmp <- integrand(grid, y, c30, C30, nu, const = const)
plot(grid, tmp, type = 'l')
```

![](Chapter10_files/figure-html/unnamed-chunk-19-1.png)

``` r
logarea <- log(sum(diff(grid) * .5 * (head(tmp, -1) + tail(tmp, -1)))) - const

logmarglikM3 <- c30 * log(C30) +
                N * lgamma((nu + 1) / 2) -
                lgamma(c30) -
                N * lgamma(nu / 2) -
                N / 2 * log(nu * pi) +
                logarea
```

We can now compute log Bayes factors and corresponding model
probabilities (under equal prior model probabilities).

``` r
(logBF <- logmarglikM1 - logmarglikM3)
#> [1] -150.6869
PrM2 <- 0.5 / (0.5 + 0.5 * exp(logBF))
PrM1 <- 1 - PrM2
knitr::kable(round(rbind(PrM1, PrM2), 3))
```

|      |     |
|:-----|----:|
| PrM1 |   0 |
| PrM2 |   1 |
