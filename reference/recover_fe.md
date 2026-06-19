# Recover individual fixed-effect levels by alternating weighted means

Faithful port of Python `recover_fixed_effects_iterative`. Decomposes a
combined fixed-effects column into per-dimension level effects so they
can be extrapolated to treated observations.

## Usage

``` r
recover_fe(df, fe_cols, combined_col, weight_col, max_iter = 100L, tol = 1e-06)
```

## Arguments

- df:

  data.frame containing `fe_cols`, `combined_col`, `weight_col`.

- fe_cols:

  character vector of fixed-effect column names.

- combined_col:

  name of the combined-effect column.

- weight_col:

  name of the weight column.

- max_iter, tol:

  convergence controls.

## Value

named list; `out[[fe]]` is a named numeric vector keyed by level.
