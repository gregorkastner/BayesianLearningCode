# The Inverse Gamma Distribution

Density, distribution function, quantile function and random generation
for the inverse gamma distribution with parameters a and b.

## Usage

``` r
dinvgamma(x, a, b, log = FALSE)

pinvgamma(q, a, b)

qinvgamma(p, a, b)

rinvgamma(ndraws, a, b)
```

## Arguments

- x:

  vector of quantiles.

- a:

  positive parameter.

- b:

  positive parameter.

- log:

  logical; if `TRUE`, densities are returned as logarithms.

- q:

  vector of quantiles.

- p:

  vector of probabilities.

- ndraws:

  number of observations.

## Value

`dinvgamma` gives the density, `pinvgamma` gives the distribution
function, `qinvgamma` gives the quantile function, and `rinvgamma`
generates random deviates.

## Details

For details, please see
[`GammaDist`](https://rdrr.io/r/stats/GammaDist.html).
