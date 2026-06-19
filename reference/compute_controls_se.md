# Compute cluster-robust SEs for control-variable coefficients (influence function)

Faithful port of Python `compute_controls_se` (did_imputation.py lines
351-404, 691-695). For each control c: demean c and the other controls
on FE (untreated subsample, after dropping singleton FE levels); regress
demeaned c on demeaned others; form per-observation ctrlweight from
residuals; cluster- sum ctrlweight \* imput_resid; sandwich V = M'M.

## Usage

``` r
compute_controls_se(dt, controls, fe, cluster, aw, df_a, beta)
```

## Arguments

- dt:

  data.table containing untreated and treated rows (needs imput_resid on
  untreated rows)

- controls:

  character vector of control column names

- fe:

  character vector of FE column names

- cluster:

  name of cluster column

- aw:

  name of analytic weight column, or NULL

- df_a:

  absorbed degrees of freedom (from .compute_df_a)

- beta:

  named numeric vector of WLS coefficients for each control

## Value

named numeric vector of cluster-robust SEs, one per control
