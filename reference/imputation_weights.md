# Compute BJS imputation weights by iterative orthogonalization

Faithful port of Python `imputation_weights` (did_imputation.py lines
100-143). Adds `copy<w>` columns to `dt` (in-place style) and returns
`dt`.

## Usage

``` r
imputation_weights(
  dt,
  wei = "wei",
  fe,
  wtr = character(),
  controls = character(),
  tol = 1e-06,
  maxit = 1000L
)
```

## Arguments

- dt:

  data.table (modified in-place)

- wei:

  name of observation-weight column

- fe:

  character vector of fixed-effect column names

- wtr:

  character vector of treatment-weight columns

- controls:

  character vector of control variable names

- tol:

  convergence tolerance

- maxit:

  maximum iterations
