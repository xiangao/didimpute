# Summarise a did_impute result

Returns a tidy `data.frame` with one row per estimand (effects,
pre-trends, and control coefficients combined).

## Usage

``` r
# S3 method for class 'did_impute'
summary(object, ...)
```

## Arguments

- object:

  a `did_impute` object.

- ...:

  ignored.

## Value

A `data.frame` with columns `term`, `estimate`, and `std.error`.

## Examples

``` r
set.seed(1)
panel <- expand.grid(i = 1:4, t = 1:6)
panel$Ei <- ifelse(panel$i <= 2, 4L, NA_integer_)
panel$y  <- 0.3 * (!is.na(panel$Ei) & panel$t >= panel$Ei) + rnorm(nrow(panel), sd = 0.1)
res <- did_impute(panel, y = "y", i = "i", t = "t", Ei = "Ei",
                  horizons = 0:1, pretrends = 2, minn = 0)
#> The number of treated entities for 'wtr0' is too small for some cohorts. Standard Errors may be wrong; consider using avgeffectsby option, averaging the the effect by treated X post variable.
#> The number of treated entities for 'wtr1' is too small for some cohorts. Standard Errors may be wrong; consider using avgeffectsby option, averaging the the effect by treated X post variable.
summary(res)
#>   term    estimate  std.error
#> 1 tau0  0.18004283 0.07299948
#> 2 tau1  0.35146589 0.05508206
#> 3 pre1 -0.02143844 0.19970161
#> 4 pre2 -0.02571258 0.14878493
```
