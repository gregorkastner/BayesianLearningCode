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
aN <- sum(y) + 1
bN <- length(y)
mu <- aN / bN
yf <- 0:20
plot(yf, dnbinom(yf, size = aN, mu = mu),
     type = "h", xlab = "", ylab = "")
plot(yf, pnbinom(yf, size = aN, mu = mu),
     type = "h", xlab = "", ylab = "")
probs <- c(0.025, 0.975)
abline(h = probs)
mtext(probs, side = 2, at = probs, adj = c(0, 1), cex = .8, col = "dimgrey")
```

![](Chapter09_files/figure-html/unnamed-chunk-3-1.png)

We inspect the parameters of the negative binomial distributation and
verify that a 95% predictive interval is given by \[1, 9\] using the
cdf.

``` r
c(aN = aN, bN = bN)
#>   aN   bN 
#> 1009  192
pnbinom(9, size = aN, mu = mu) - pnbinom(0, size = aN, mu = mu)
#> [1] 0.9522248
```

### Example 9.2: Exchange Rate Data

We load the data and then plot the pdf and cdf for the predictive
distribution which corresponds under the improper prior to
$$y_{f}|\mathbf{y} \sim {\mathcal{t}}_{2c_{N}}\left( \bar{y},\frac{Ns_{y}^{2}}{N - 1}\left( 1 + \frac{1}{N} \right) \right).$$

``` r
par(mfrow = c(1, 2))
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
CN <- C0 + 1/2 * sum((y - ybar)^2) + N0 * N / (2 * (N0 + N)) * (b0 - ybar)^2
x <- seq(-5, 5, length.out = 200)
scale <- CN / cN * (BN + 1)
plot(x, dstudt(x, location = bN, scale = scale, df = 2 * cN),
     type = "l", xlab = "", ylab = "")
plot(x, pstudt(x, location = bN, scale = scale, df = 2 * cN),
     type = "l", xlab = "", ylab = "")
probs <- c(0.025, 0.975)
abline(h = probs)
mtext(probs, side = 2, at = probs, adj = c(0, 1), cex = .8, col = "dimgrey")
```

![](Chapter09_files/figure-html/unnamed-chunk-5-1.png)

We inspect the parameters of the negative binomial distributation and
determine the 95% predictive interval using the quantile function.

``` r
round(c(location = bN, scale = scale, df = 2 * cN), digits = 3)
#> location    scale       df 
#>    0.018    0.528 3138.000
round(qstudt(probs, location = bN, scale = scale, df = 2 * cN), digits = 3)
#> [1] -1.018  1.053
```
