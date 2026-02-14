# Chapter 9: Bayesian Predictive Analysis

## Section 9.1

### Example 9.1: Road Safety Data

We load the data and extract the observations for the senior people in
Linz. We then plot the pdf and cdf for the predictive distribution which
corresponds under the flat prior to
$$y_{f}|\mathbf{y} \sim \mathcal{N}B\left( N\bar{y} + 1,N \right).$$

``` r
data("accidents", package = "BayesianLearningCode")
y <- accidents[, "seniors_accidents"]
(aN <- sum(y) + 1)
#> [1] 1009
(bN <- length(y))
#> [1] 192
mu <- aN / bN
yf <- 0:20
plot(yf, dnbinom(yf, size = aN, mu = mu),
     type = "h", xlab = bquote(y[f]), ylab = "")
plot(yf, pnbinom(yf, size = aN, mu = mu),
     type = "h", xlab = bquote(y[f]), ylab = "")
probs <- c(0.025, 0.975)
abline(h = probs, lty = 3)
mtext(probs, side = 2, at = probs, adj = c(0, 1), cex = .8, col = "dimgrey")
```

![](Chapter09_files/figure-html/unnamed-chunk-3-1.png)

### Example 9.2: Exchange Rate Data

We load the data and then plot the pdf and cdf for the predictive
distribution which corresponds under the improper prior to
$$y_{f}|\mathbf{y} \sim {\mathcal{t}}_{2c_{N}}\left( \bar{y},\frac{Ns_{y}^{2}}{N - 1}\left( 1 + \frac{1}{N} \right) \right).$$

``` r
library("BayesianLearningCode")
data("exrates", package = "stochvol")
y <- 100 * diff(log(exrates$USD / exrates$CHF))
ybar <- mean(y)
N <- length(y)
b0 <- 0
N0 <- 0
c0 <- -1/2
C0 <- 0
BN <- 1 / (N + N0)
bN <- BN * (N * ybar + N0 * b0)
cN <- c0 + N/2
CN <- C0 + .5 * sum((y - ybar)^2) + N0 * N / (2 * (N0 + N)) * (b0 - ybar)^2
x <- seq(-3, 3, length.out = 200)
scale <- sqrt(CN / cN * (BN + 1))
plot(x, dstudt(x, location = bN, scale = scale, df = 2 * cN),
     type = "l", xlab = bquote(y[f]), ylab = "")
plot(x, pstudt(x, location = bN, scale = scale, df = 2 * cN),
     type = "l", xlab = bquote(y[f]), ylab = "")
probs <- c(0.025, 0.975)
abline(h = probs, lty = 3)
mtext(probs, side = 2, at = probs, adj = c(0, 1), cex = .8, col = "dimgrey")
```

![](Chapter09_files/figure-html/unnamed-chunk-4-1.png)

We inspect the parameters of the Student-t distribution.

``` r
round(c(mu = bN, sigma2 = scale^2, df = 2 * cN), digits = 3)
#>       mu   sigma2       df 
#>    0.018    0.528 3138.000
```

### Example 9.3: Exchange Rate Data cont’d

To compare with a method that ignores parameter uncertainty, we now plot
the posterior predictive alongside the “classical” forecasting
distribution for varying $N$.

``` r
Ns <- c(5, 10, 30, N)

for (i in seq_along(Ns)) {
  N <- Ns[i]
  yshort <- y[1:N]
  ybar <- mean(yshort)
  BN <- 1 / (N + N0)
  bN <- BN * (N * ybar + N0 * b0)
  cN <- c0 + N/2
  CN <- C0 + .5 * sum((yshort - ybar)^2) + N0 * N / (2 * (N0 + N)) * (b0 - ybar)^2
  scale <- sqrt(CN / cN * (BN + 1))
  plot(x, dnorm(x, ybar, sd(yshort)), lty = 2, type = "l",
       xlab = bquote(y[f]), ylab = "", main = paste("N =", N))
  lines(x, dstudt(x, location = bN, scale = scale, df = 2 * cN))
}
```

![](Chapter09_files/figure-html/unnamed-chunk-6-1.png)

### Example 9.4: Road Safety Data cont’d

We verify that a 95% predictive interval is given by \[1, 9\] using the
cdf and compute the effective coverage.

``` r
pnbinom(9, size = aN, mu = mu) - pnbinom(0, size = aN, mu = mu)
#> [1] 0.9522248
```

### Example 9.5: Exchange Rate Data cont’d

We determine the 95% predictive interval using the quantile function.

``` r
round(qstudt(probs, location = bN, scale = scale, df = 2 * cN), digits = 3)
#> [1] -1.408  1.443
```

We now work out the effective coverage of the naive interval which
ignores parameter uncertainty.

``` r
coverage <- rep(NA_real_, length(Ns))
names(coverage) <- Ns
for (i in seq_along(Ns)) {
  N <- Ns[i]
  yshort <- y[1:N]
  quants <- qnorm(probs, mean(yshort), sd(yshort))
  coverage[i] <- diff(pstudt(quants, mean(yshort), sd(yshort) * sqrt(1 + 1/N),
                             2 * (c0 + N/2)))
}
knitr::kable(t(round(coverage, 4)))
```

|      5 |     10 |     30 |   3139 |
|-------:|-------:|-------:|-------:|
| 0.8519 | 0.9055 | 0.9363 | 0.9499 |
