# Difference-in-differences imputation estimator (BJS 2024)

Implements the imputation estimator of Borusyak, Jaravel and Spiess
(2024, AER) for staggered difference-in-differences designs. The
estimator imputes each treated unit's counterfactual outcome from a
two-way fixed effects model fit on untreated observations, then
aggregates the resulting unit-level treatment effects into
user-specified estimands. Standard errors are cluster-robust (clustered
by unit `i` by default) using the smartweights influence-function
approach of the original Python package.

## Usage

``` r
did_impute(
  df,
  y,
  i,
  t,
  Ei,
  controls = character(),
  fe = NULL,
  timecontrols = character(),
  aw = NULL,
  unitcontrols = character(),
  wtr = character(),
  sum = FALSE,
  horizons = NULL,
  allhorizons = FALSE,
  hbalance = FALSE,
  hetby = character(),
  project = character(),
  minn = 30,
  saveweights = FALSE,
  shift = 0,
  pretrends = 0,
  cluster = "",
  avgeffectsby = character(),
  leaveoneout = FALSE,
  nose = FALSE,
  delta = NULL,
  seed = 1L
)
```

## Arguments

- df:

  A data frame (or data.table) in long format.

- y:

  Character. Name of the outcome column.

- i:

  Character. Name of the unit identifier column.

- t:

  Character. Name of the time column (must be integer-valued or
  coercible to integer).

- Ei:

  Character. Name of the treatment-timing column. `NA` indicates a
  never-treated unit; a finite integer is the first treated period.

- controls:

  Character vector of time-varying control variable names (default:
  none).

- fe:

  Character vector of fixed-effect column names. Defaults to `c(i, t)`
  (the standard two-way FE).

- timecontrols:

  Not yet implemented; passing a non-empty vector raises an error.

- aw:

  Character. Name of an analytic-weight column, or `NULL` (default).

- unitcontrols:

  Not yet implemented; passing a non-empty vector raises an error.

- wtr:

  Character vector of user-supplied treatment-weight column names.
  Cannot be combined with `horizons`, `allhorizons`, or `hbalance`.

- sum:

  Logical. If `TRUE`, sum effects rather than averaging (weights are not
  normalised to sum to 1 over treated rows). Default `FALSE`.

- horizons:

  Numeric vector of relative-time horizons at which to estimate
  event-study effects (e.g. `0:3`).

- allhorizons:

  Logical. If `TRUE`, estimate effects at all observed relative-time
  horizons (inferred from the data). Cannot be combined with `horizons`.

- hbalance:

  Logical. If `TRUE` and `horizons` is specified, restrict to cohorts
  that are observed at every horizon in `horizons` (horizon-balanced
  sample).

- hetby:

  Not yet implemented; passing a non-empty vector raises an error.

- project:

  Not yet implemented; passing a non-empty vector raises an error.

- minn:

  Minimum effective sample size (default 30). An estimand is suppressed
  (weight set to zero) when its effective *N* falls below this
  threshold. Set to 0 to disable.

- saveweights:

  Logical. If `TRUE`, include the per-observation treatment weights in
  the returned object (not yet used by downstream methods; included for
  API parity with the Python package).

- shift:

  Integer. Time shift applied when constructing relative time and the
  untreated indicator (default 0).

- pretrends:

  Non-negative integer. Number of pre-treatment horizons at which to
  estimate placebo ("pre-trend") coefficients. Default 0 (none).

- cluster:

  Character. Name of the cluster column for cluster-robust standard
  errors. Defaults to `i` (cluster by unit).

- avgeffectsby:

  Character vector. Columns that define the estimand grouping used for
  the smartweights influence function. Defaults to `c(Ei, t)`.

- leaveoneout:

  Not yet implemented; setting to `TRUE` raises an error.

- nose:

  Logical. If `TRUE`, skip standard error computation. Estimates are
  still returned. Default `FALSE`.

- delta:

  Integer or `NULL`. Time-step size. If `NULL` (default) the step is
  detected automatically as the modal first difference of the observed
  time values.

- seed:

  Integer seed passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html) before the randomised
  rank / collinearity check on the control variables. Setting a fixed
  seed makes this check reproducible across calls. This is the one
  deliberate deviation from the upstream Python package, which uses a
  non-seeded random draw. Default `1L`.

## Value

An S3 object of class `"did_impute"` with the following components:

- `estimates`:

  Named list of point estimates. The baseline case (no horizons, no user
  `wtr`) produces a single element `tau_ate`. Horizon estimates are
  named `tau0`, `tau1`, ...; user-weight estimates use the weight column
  names.

- `std_errors`:

  Named list of cluster-robust standard errors, matching `estimates`.
  `NULL` when `nose = TRUE`.

- `pretrends_estimates`:

  Named list of pre-trend (placebo) coefficients (`pre1`, `pre2`, ...).
  `NULL` when `pretrends = 0`.

- `pretrends_std_errors`:

  Named list of cluster-robust SEs for the pre-trend coefficients.
  `NULL` when `pretrends = 0` or `nose = TRUE`.

- `controls_estimates`:

  Named list of WLS coefficients on `controls`. `NULL` when `controls`
  is empty.

- `controls_std_errors`:

  Named list of cluster-robust SEs for the control coefficients. `NULL`
  when `controls` is empty or `nose = TRUE`.

- `n_obs`:

  Total observation count entering the estimand computation (untreated
  rows plus treated rows with non-zero weight).

- `V`:

  A square covariance matrix of dimension \\K \times K\\, where \\K\\ is
  the total number of reported estimands (effects, then pre-trends, then
  controls, in that order). Row and column names match the estimand
  names. `sqrt(diag(V))` reproduces the reported per-estimand standard
  errors to numerical precision. This is a deliberate improvement over
  the upstream Python package, whose `V` was a scalar (sum of squared
  SEs). `NULL` when `nose = TRUE` or when no standard-error components
  are present.

- `weights`:

  Reserved for future use (`NULL` unless `saveweights = TRUE` is
  requested; the field is included for API parity with the Python
  package).

## Details

**Absorbed degrees of freedom.** The cluster-robust SE scaling factor
uses absorbed degrees of freedom (\\df_a\\) computed by a port of the
*pyhdfe* pairwise method: an FE that is nested within the cluster
variable is dropped; the remaining FEs contribute their unique-level
counts minus the number of bipartite connected components shared across
pairs. This is exact for the default two-way FE / cluster-by-unit case
and for any configuration with at most two non-nested fixed effects. For
more than two non-nested fixed effects the value is an approximation and
a warning is emitted.

**Unimplemented options.** The arguments `timecontrols`, `unitcontrols`,
`leaveoneout`, `hetby`, and `project` match the Python package API but
are not yet implemented; they raise an error if supplied.

## References

Borusyak, K., Jaravel, X., and Spiess, J. (2024). Revisiting Event-Study
Designs: Robust and Efficient Estimation. *Review of Economic Studies*,
91(6), 3253–3285.

## See also

[`event_plot`](https://xiangao.github.io/didimpute/reference/event_plot.md),
[`summary.did_impute`](https://xiangao.github.io/didimpute/reference/summary.did_impute.md),
[`print.did_impute`](https://xiangao.github.io/didimpute/reference/print.did_impute.md)

## Examples

``` r
# Minimal synthetic staggered panel (4 units, 6 periods)
set.seed(42)
n_units <- 4; n_t <- 6
panel <- expand.grid(i = seq_len(n_units), t = seq_len(n_t))
panel$Ei <- ifelse(panel$i <= 2, 4L, NA_integer_)  # units 1-2 treated at t=4
panel$y  <- 0.5 * (!is.na(panel$Ei) & panel$t >= panel$Ei) +
              rnorm(nrow(panel), sd = 0.2)

# Baseline ATT estimate
res <- did_impute(panel, y = "y", i = "i", t = "t", Ei = "Ei")
#> WARNING: suppressing wtr, consider lower minn or minn=0.
#> The number of treated entities is too small for some cohorts. Standard Errors may be wrong, consider using avgeffectsby option, averaging the the effect by treated X post variable.
print(res)
#> <didimpute did_impute result>
#>   Effects: tau_ate=0 
#>   SE     : tau_ate=0 
#>   N obs: 18 

# Event-study with horizons 0 and 1 and one pre-period
res2 <- did_impute(panel, y = "y", i = "i", t = "t", Ei = "Ei",
                   horizons = 0:1, pretrends = 1, minn = 0)
#> The number of treated entities for 'wtr0' is too small for some cohorts. Standard Errors may be wrong; consider using avgeffectsby option, averaging the the effect by treated X post variable.
#> The number of treated entities for 'wtr1' is too small for some cohorts. Standard Errors may be wrong; consider using avgeffectsby option, averaging the the effect by treated X post variable.
summary(res2)
#>   term    estimate std.error
#> 1 tau0  0.38111364 0.1955443
#> 2 tau1  0.41610601 0.2809816
#> 3 pre1 -0.09815805 0.1807401
```
