# The Non-standardized Student t Distribution

Density, distribution function, quantile function and random generation
for the location-scale-transformed Student t distribution.

## Usage

``` r
dstudt(x, location = 0, scale = 1, df, log = FALSE)

pstudt(q, location = 0, scale = 1, df)

qstudt(p, location = 0, scale = 1, df)

rstudt(n, location = 0, scale = 1, df)
```

## Arguments

- x:

  vector of quantiles.

- location:

  parameter.

- scale:

  positive parameter.

- df:

  positive parameter.

- log:

  logical; if `TRUE`, densities are returned as logarithms.

- q:

  vector of quantiles.

- p:

  vector of probabilities.

- n:

  number of observations.

## Value

`dstudt` gives the density, `pstudt` gives the distribution function,
`qstudt` gives the quantile function, and `rstudt` generates random
deviates.

## Details

The non-standardized Student t distribution with location \\\mu\\, scale
\\\tau\\, and df \\\nu \> 2\\ has mean \\\mu\\ and variance \\\tau^2
\frac{\nu}{\nu-2}\\.
