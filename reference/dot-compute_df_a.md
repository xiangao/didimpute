# Compute absorbed degrees of freedom matching pyhdfe pairwise method

Mirrors the pyhdfe `algorithm.degrees` logic (used with cluster_ids): an
FE is dropped if it is nested within (i.e., its levels are a refinement
of) any cluster group. The remaining FEs contribute their unique level
counts, minus the number of connected components across pairs (pairwise
method).

## Usage

``` r
.compute_df_a(dt, fe, cluster)
```

## Arguments

- dt:

  data.table restricted to the untreated subsample

- fe:

  character vector of FE column names

- cluster:

  name of cluster column

## Value

integer absorbed degrees of freedom

## Details

For the standard case (fe = c(i, t), cluster = i): the i-FE is nested
within the i-cluster, so only the t-FE contributes =\> df_a =
uniqueN(t).
