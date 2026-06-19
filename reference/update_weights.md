# Orthogonalize weight columns against one variable/FE stratum

Faithful port of Python `update_weights` (did_imputation.py lines
62-95). Modifies `dt` in-place and returns it.

## Usage

``` r
update_weights(
  dt,
  varlist = NULL,
  w = character(),
  wei = "wei",
  d = "untreated",
  denom = "",
  by = character()
)
```

## Arguments

- dt:

  data.table

- varlist:

  column name to project against (character scalar); NULL = intercept

- w:

  character vector of weight columns to update

- wei:

  name of observation-weight column

- d:

  name of untreated indicator column

- denom:

  name of denominator column

- by:

  character vector of grouping columns (empty = global sum)
