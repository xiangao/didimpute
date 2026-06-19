# didimpute — project instructions

## What this is
`didimpute` is an R package: a **faithful port of the Python `did_imputation` package**
(Georgii Marinichev's implementation, https://github.com/Gmarinichev/did_imputation) of the
Borusyak, Jaravel & Spiess (2024) difference-in-differences **imputation estimator**, plus an
`event_plot`. Faithfulness is to *that Python implementation's numerical output* — not to the
original Stata package, and independent of Kyle Butts's separate CRAN `didimputation`.

## Layout
- Package (shippable): `~/projects/software/didimpute` (this repo).
- Port workspace (NOT part of the package): `~/projects/claude/didimpute_port` — Python `.venv`,
  the golden-grid generator (`code/gen_golden.py`, `code/make_panels.py`), saved golden CSVs, the
  design spec, the implementation plan, and the SDD progress ledger. Mirrors the `drlate_port` precedent.

## Engine
- `fixest::demean()` is the analog of Python's `pyhdfe.residualize()` (HDFE absorption on the
  untreated subsample). `stats::lm.wfit` recovers control coefficients. The BJS-specific pieces
  (FE-level recovery, imputation-weight iteration, smartweight/influence-function SEs) are
  hand-translated — no R library implements them.
- `data.table` throughout. NOTE the recurring traps that are already handled: `get(col)` mis-resolves
  inside `dt[...]` when a column shares the name (capture column-name strings / use `dt[[col]]`);
  never rely on merge row order (re-extract vectors after each merge); a never-treated row's
  `Rel_time` is `NA` in R (Python `NaN`) — every Rel_time-derived indicator MUST guard `NA -> 0`.

## Files (R/)
- `did_impute.R` — orchestrator: prep, delta detection, untreated mask, imputation, aggregation,
  `df_a`, SE wiring, combined `$V` assembly. (Largest file.)
- `recover_fe.R` — iterative FE-level recovery (alternating weighted means).
- `imputation_weights.R` — `update_weights`, `imputation_weights` (BJS weight iteration).
- `standard_errors.R` — `compute_effect_se` (smartweight), `compute_controls_se`,
  `compute_pretrends`, `.compute_df_a`.
- `output.R` — S3 `print`/`summary`.
- `event_plot.R` — ggplot2 event-study plot.

## Validation (do not weaken)
- Correctness is enforced by a **Python golden grid**: 14 configurations matched against the Python
  package — estimates AND standard errors — to **~1e-6** (most to 1e-8). Tests live in
  `tests/testthat/`; the golden CSVs are in `tests/testthat/fixtures/`.
- To regenerate the golden grid: in `~/projects/claude/didimpute_port`, `source .venv/bin/activate`
  then `python code/gen_golden.py`, and copy `output/golden/*_result.csv` (+ any `*_data.csv`) into
  `tests/testthat/fixtures/`. Each `did_imputation(...)` call MUST be immediately preceded by
  `np.random.seed(1)` (matches the R `seed=1` so the one randomized rank-check path agrees).
- The shared panel is **150 units** so `minn=30` does not suppress horizon estimands. Do not shrink it.
- Run tests: `Rscript -e 'devtools::load_all(); testthat::test_local()'`. `R CMD check` is clean
  (the only warning is the vignette build needing Pandoc on PATH, which is environmental).

## Deliberate decisions / known limitations
- `$V` is the **full cluster-robust covariance matrix** of all reported estimands (effects,
  pretrends, controls), `sqrt(diag(V))` reproduces the SEs. This is a *deliberate improvement* over
  the upstream Python package, whose `V` was a scalar (sum of squared SEs).
- The one deliberate deviation from Python: the randomized rank/collinearity check is seeded via the
  `seed` argument (default `1L`) for reproducibility.
- Absorbed degrees of freedom (`df_a`) for the cluster-robust SEs are **exact** for the default
  two-way FE / cluster-by-unit case and for ≤2 fixed effects; for >2 non-nested fixed effects the
  value is **approximate** and a warning is emitted.
- Unimplemented options (faithful to Python — they error): `timecontrols`, `unitcontrols`,
  `leaveoneout`, `hetby`, `project`. `fe=NULL` maps to the two-way FE `c(i,t)`; the Python
  intercept-only (`_const_`) path is not implemented (no golden case needs it).

## Workflow
- This package was built via subagent-driven TDD against the golden grid. Any change that touches the
  estimator MUST keep the golden grid matching to 1e-6 — if it diverges, debug against Python
  (temporary dumps in `gen_golden.py`), never loosen the tolerance to pass.
- Docs deploy via pkgdown; README links to the pkgdown site. Re-run `devtools::document()` after
  roxygen edits.
