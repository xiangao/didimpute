# Compute cluster-robust standard errors via BJS smartweights

Faithful port of Python `compute_effect_se` (did_imputation.py lines
615-711, no-controls path). Operates on a COPY of dt (caller must pass
one).

## Usage

``` r
compute_effect_se(
  dt,
  wtr,
  y,
  cluster,
  avgeffectsby,
  gr_var,
  fe,
  controls = character()
)
```

## Arguments

- dt:

  data.table (will be modified; pass a copy)

- wtr:

  character vector of original treatment-weight column names

- y:

  name of outcome column

- cluster:

  name of cluster column

- avgeffectsby:

  character vector of columns defining the estimand group

- gr_var:

  character vector = union of cluster + avgeffectsby

- fe:

  character vector of fixed-effect columns (passed to
  imputation_weights)

- controls:

  character vector of control variable names (default empty)

## Value

list with `se` (named numeric), `group_sums` matrix, `V` matrix

## Details

For each weight column in `wtr`:

- Fill NaN copy with original w (non-iterated weights)

- Compute smartweights on treated rows

- Residual: y - y_hat (all rows), overridden to tau - avg_taus (treated
  rows)

- product = resid \* copy \* wei

- group_sums = sum(product) per cluster Then V = G'G where G =
  cbind(group_sums); ses = sqrt(diag(V)).
