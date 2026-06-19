# Compute cluster-robust standard errors for pretrend coefficients

Faithful port of Python pretrends block (did_imputation.py lines
477-566, 700-704). Operates on the untreated subsample after dropping
singleton FE levels. For each h=1..pretrends: residualise pretrendvar_h
on FE (and other pretrend/control vars); WLS to get preresid; form
per-observation preweight; normalise; cluster-sum preweight \* preresid
\* sqrt(dof_adj). Sandwich V = M'M.

## Usage

``` r
compute_pretrends(dt, y, fe, cluster, controls, pretrends, aw)
```

## Arguments

- dt:

  data.table (all rows, including treated — untreated filter applied
  inside)

- y:

  name of outcome column

- fe:

  character vector of FE column names

- cluster:

  name of cluster column

- controls:

  character vector of control column names

- pretrends:

  integer number of pretrend horizons (h = 1..pretrends)

- aw:

  name of analytic weight column, or NULL

## Value

list with `coefs` (named numeric) and `ses` (named numeric), plus
`list_pre_weps` (matrix of per-cluster influence vectors)
