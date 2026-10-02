# Chapter 12: Model Assessment and Predictive Model Checking

## Section 12.1: Beyond classical Bayesian model comparison

### Section 12.1.1: Schwarz criterion and BIC

#### Example 12.1: Labor market data: Bayesian variable selection using BIC

We return to the probit analysis of Example 11.11 and run variable
selection using the Schwarz criterion instead of the marginal
likelihood. We first implement this approach.

``` r

modsel_probit_BIC <- function(y, X, burnin = 1000L, M = 5000L){
  n <- dim(X)[1]
  p <- dim(X)[2] # change to X later

  gamma_post <- matrix(ncol = p, nrow = M)
  acc <- numeric(length = M)

  gamma <- matrix(rep(1, p), nrow = 1)
  mod0 <- glm(y ~ X, family = binomial(link = "probit"))
  SC_old <- logLik(mod0) - log(n) * (p + 1) / 2

  for (m in seq_len(burnin + M)) {
    gamma_proposed <- gamma_old <- gamma

    j <- sample(1:p, size = 1)
    gamma_proposed[j] <- 1 - gamma_old[j]

    X_gamma <- X[, gamma_proposed == 1]
    model <- glm(y ~ X_gamma, family = binomial(link = "probit"))
    SC_proposed <- logLik(model) - log(n) * (sum(gamma_proposed) + 1) / 2
    
    # compute acceptance probability and decide on acceptance
    log_acc <- SC_proposed - SC_old

    if (log(runif(1)) < log_acc) {
      gamma <- gamma_proposed
      SC_old <- SC_proposed
      accept <- 1
    } else {
      gamma <- gamma_old
      accept <- 0
    }

    if (m > burnin) {
      gamma_post[m - burnin, ] <- gamma
      acc[m - burnin] <- accept
    }
  }
  return(list(gamma_post = gamma_post, acc = acc))
}
```

We next load the data and as in exercise 11.11. add 5 irrelevant,
standard normal regressors. Then we run the algorithm.

``` r

library("BayesianLearningCode")
data("labor", package = "BayesianLearningCode")
y <- labor$income_1998 == "zero"
N <- length(y)
X_unemp <- with(labor, cbind(female = female,
                             age18 = 1998 - birthyear - 18,
                             wcollar = wcollar_1997,
                             unemp97 = income_1997 == "zero")) # regressor matrix
p_irrel <- 5
set.seed(seed)
X <- cbind(X_unemp, matrix(rnorm(p_irrel * N), ncol = p_irrel))
colnames(X)[4 + (1:p_irrel)] <-
  c("irrel1", "irrel2", "irrel3", "irrel4", "irrel5")

p <- dim(X)[2]
M <- 20000 / mcmcspeedup
res <- modsel_probit_BIC(y, X, M = M)
```

We determine the number of models visited during MCMC, the average model
complexity, i.e. the average number of covariates included in the model
and show the posterior inclusion probabilities.

``` r

if (pdfplots) {
  pdf("12-1_1.pdf", width = 5, height = 6)
}
par(mfrow = c(1, 1), mar = c(2.5, 2.5, 1.5, .5), mgp = c(1.5, .5, 0))

gammas_unique <- unique(res$gamma_post)
num_mod <- dim(gammas_unique)[1]
print(num_mod)
#> [1] 21

k_gamma <- rowSums(res$gamma_post)
print(mean(k_gamma))
#> [1] 3.7335

# PIPs
print(colMeans(res$gamma_post))
#> [1] 0.8995 1.0000 0.7020 1.0000 0.0320 0.0085 0.0065 0.0710 0.0140

barplot(colMeans(res$gamma_post), col = "blue", names.arg = 1:p,
        xlab = "Covariate", ylab = "PIP")
```

![](Chapter12_files/figure-html/unnamed-chunk-4-1.png)

To determine the model visited most often we re-use a function which we
defined in Chapter 11.

``` r

number_draws <- function(gamma_post, models) {
  nmod <- dim(models)[1]
  freq <- rep(NA, nmod)
  
  for (j in (1:nmod)) {
    freq[j] <- sum(apply(gamma_post, 1, function(x) identical(x, models[j,])))
  }
  freq
}

freq_gammas <- number_draws(res$gamma_post, gammas_unique)
io <- order(freq_gammas, decreasing = TRUE)

knitr::kable(cbind(gammas_unique, freq_gammas / M)[io[1:5], ],
             digits = cbind(rep(0, p), 4))
```

|     |     |     |     |     |     |     |     |     |        |
|----:|----:|----:|----:|----:|----:|----:|----:|----:|-------:|
|   1 |   1 |   1 |   1 |   0 |   0 |   0 |   0 |   0 | 0.5915 |
|   1 |   1 |   0 |   1 |   0 |   0 |   0 |   0 |   0 | 0.2060 |
|   0 |   1 |   0 |   1 |   0 |   0 |   0 |   0 |   0 | 0.0580 |
|   1 |   1 |   1 |   1 |   0 |   0 |   0 |   1 |   0 | 0.0380 |
|   0 |   1 |   1 |   1 |   0 |   0 |   0 |   0 |   0 | 0.0255 |

### Section 12.1.2: Perspectives on Bayesian model comparison

\[no code\]

## Section 12.2: Comparative model assessment through predictive information criteria

### Section 12.2.1: Measuring predictive loss

\[no code\]

### Section 12.2.2: Akaike information criterion (AIC)

#### Example 12.2: U.S. GDP data: Choosing the model order via AIC and BIC

We need to re-load the function which creates the design matrix.

``` r

ARdesignmatrix <- function(dat, p = 1, conditioninglength = p) {
  d <- p + 1
  N <- length(dat) - p

  Xy <- matrix(NA_real_, N, d)
  Xy[, 1] <- 1
  for (i in seq_len(p)) {
    Xy[, i + 1] <- dat[(p + 1 - i) : (length(dat) - i)]
  }
  Xy[(1 + conditioninglength - p):N, , drop = FALSE]
}
```

Let’s load the data and compute log returns.

``` r

data(gdp, package = "BayesianLearningCode")
logret <- diff(log(gdp))
```

Now we compute five OLS estimates and the corresponding AICs and BICs.

``` r

AIC <- BIC <- rep(NA_real_, 5)
for (i in seq_along(AIC)) {
  y <- tail(logret, -4)
  p <- i - 1     # lag length
  d <- i + 1     # total number of parameters
  n <- length(y) # number of data points used for estimation
  Xy <- ARdesignmatrix(logret, p, conditioninglength = 4)
  fit <- lm(y ~ Xy - 1)
  betahat <- coef(fit)
  sigma2hat <- sum(resid(fit)^2) / length(y)
  minus2timesloglikmax <- n * (1 + log(2 * pi)) + n * log(sigma2hat)
  AIC[i] <- minus2timesloglikmax + 2 * d      # equals bult-in AIC(fit)
  BIC[i] <- minus2timesloglikmax + log(n) * d # equals built-in BIC(fit)
}
knitr::kable(rbind(AIC, BIC), digits = 2)
```

|     |          |          |          |          |          |
|:----|---------:|---------:|---------:|---------:|---------:|
| AIC | -1808.61 | -1828.87 | -1839.65 | -1837.66 | -1836.24 |
| BIC | -1801.15 | -1817.68 | -1824.72 | -1819.01 | -1813.86 |
