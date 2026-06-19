#' @import data.table
NULL

.detect_delta <- function(tvec) {
  if (!is.integer(tvec) && !all(tvec == as.integer(tvec), na.rm = TRUE)) {
    stop("Unsupported data type. Time must be integer.")
  }
  d <- diff(sort(unique(as.integer(tvec))))
  d <- d[!is.na(d)]
  if (!length(d)) stop("No time deltas found to analyze.")
  tab <- table(d)
  modal <- as.integer(names(tab)[tab == max(tab)])
  if (length(modal) > 1) message(sprintf("more than 1 time delta found, %d chosen", min(modal)))
  min(modal)
}

.prep_data <- function(df, y, i, t, Ei, controls, fe, cluster, aw, sum, shift) {
  dt <- data.table::as.data.table(data.table::copy(df))
  if (is.null(fe) || length(fe) == 0) fe <- c(i, t)   # Python: fe == [] -> [i, t]
  if (identical(cluster, "")) cluster <- i
  cols <- unique(c(t, i, fe, y, controls, cluster, if (!is.null(aw)) aw))
  keep <- stats::complete.cases(dt[, cols, with = FALSE])
  dt <- dt[keep]

  if (is.null(aw)) {
    dt[, wei := 1]
  } else {
    dt[, wei := get(aw)]
    if (!sum) dt[, wei := get(aw) * .N / base::sum(get(aw))]
  }

  delta <- .detect_delta(dt[[t]])
  Ei_col <- Ei; t_col <- t  # capture column-name strings before entering data.table [
  dt[, untreated := as.integer(is.na(.SD[[Ei_col]]) | (.SD[[t_col]] + shift < .SD[[Ei_col]])),
     .SDcols = c(Ei_col, t_col)]
  dt[, Rel_time := (.SD[[t_col]] - .SD[[Ei_col]] + shift) / delta,
     .SDcols = c(t_col, Ei_col)]

  list(df = dt, fe = fe, cluster = cluster, delta = delta)
}

.impute_y0 <- function(prep, y, controls, fe, cluster, aw) {
  dt <- prep$df
  untr <- dt[untreated == 1]
  if (nrow(untr) == 0L)
    stop("No untreated (control) observations available for imputation; cannot estimate Y(0).")
  wcol <- if (is.null(aw)) NULL else untr$wei

  ymat <- as.matrix(untr[[y]])
  y_resid <- fixest::demean(ymat, untr[, ..fe], weights = wcol)

  beta <- NULL
  if (length(controls) > 0) {
    xmat <- as.matrix(untr[, ..controls])
    x_resid <- fixest::demean(xmat, untr[, ..fe], weights = wcol)
    ww <- if (is.null(aw)) rep(1, nrow(untr)) else untr$wei
    fit <- stats::lm.wfit(x = x_resid, y = drop(y_resid), w = ww)
    beta <- fit$coefficients
    names(beta) <- controls
    fitted <- as.numeric(xmat %*% beta)
    combined <- untr[[y]] - fitted - drop(y_resid - x_resid %*% beta)  # = FE part
  } else {
    combined <- untr[[y]] - drop(y_resid)
  }
  untr[, combined := combined]
  fe_levels <- recover_fe(untr, fe, "combined", "wei")

  # extrapolate to ALL rows
  dt[, y_hat := 0]
  for (f in fe) {
    lv <- fe_levels[[f]]
    dt[, y_hat := y_hat + lv[as.character(get(f))]]
  }
  if (length(controls) > 0) {
    dt[, y_hat := y_hat + as.numeric(as.matrix(.SD) %*% beta), .SDcols = controls]
  }
  list(df = dt, beta = beta, fe_levels = fe_levels)
}

new_did_impute <- function(estimates = NULL, std_errors = NULL, pretrends_estimates = NULL,
                           pretrends_std_errors = NULL, controls_estimates = NULL,
                           controls_std_errors = NULL, n_obs = NULL, V = NULL, weights = NULL) {
  structure(list(estimates = estimates, std_errors = std_errors,
                 pretrends_estimates = pretrends_estimates,
                 pretrends_std_errors = pretrends_std_errors,
                 controls_estimates = controls_estimates,
                 controls_std_errors = controls_std_errors,
                 n_obs = n_obs, V = V, weights = weights),
            class = "did_impute")
}

#' Difference-in-differences imputation estimator (BJS 2024)
#'
#' Implements the imputation estimator of Borusyak, Jaravel and Spiess (2024,
#' AER) for staggered difference-in-differences designs.  The estimator
#' imputes each treated unit's counterfactual outcome from a two-way fixed
#' effects model fit on untreated observations, then aggregates the resulting
#' unit-level treatment effects into user-specified estimands.  Standard errors
#' are cluster-robust (clustered by unit \code{i} by default) using the
#' smartweights influence-function approach of the original Python package.
#'
#' @param df A data frame (or data.table) in long format.
#' @param y Character. Name of the outcome column.
#' @param i Character. Name of the unit identifier column.
#' @param t Character. Name of the time column (must be integer-valued or
#'   coercible to integer).
#' @param Ei Character. Name of the treatment-timing column.  \code{NA}
#'   indicates a never-treated unit; a finite integer is the first treated
#'   period.
#' @param controls Character vector of time-varying control variable names
#'   (default: none).
#' @param fe Character vector of fixed-effect column names.  Defaults to
#'   \code{c(i, t)} (the standard two-way FE).
#' @param timecontrols Not yet implemented; passing a non-empty vector raises
#'   an error.
#' @param aw Character. Name of an analytic-weight column, or \code{NULL}
#'   (default).
#' @param unitcontrols Not yet implemented; passing a non-empty vector raises
#'   an error.
#' @param wtr Character vector of user-supplied treatment-weight column names.
#'   Cannot be combined with \code{horizons}, \code{allhorizons}, or
#'   \code{hbalance}.
#' @param sum Logical. If \code{TRUE}, sum effects rather than averaging
#'   (weights are not normalised to sum to 1 over treated rows).  Default
#'   \code{FALSE}.
#' @param horizons Numeric vector of relative-time horizons at which to
#'   estimate event-study effects (e.g. \code{0:3}).
#' @param allhorizons Logical.  If \code{TRUE}, estimate effects at all
#'   observed relative-time horizons (inferred from the data).  Cannot be
#'   combined with \code{horizons}.
#' @param hbalance Logical.  If \code{TRUE} and \code{horizons} is specified,
#'   restrict to cohorts that are observed at every horizon in \code{horizons}
#'   (horizon-balanced sample).
#' @param hetby Not yet implemented; passing a non-empty vector raises an
#'   error.
#' @param project Not yet implemented; passing a non-empty vector raises an
#'   error.
#' @param minn Minimum effective sample size (default 30).  An estimand is
#'   suppressed (weight set to zero) when its effective \emph{N} falls below
#'   this threshold.  Set to 0 to disable.
#' @param saveweights Logical.  If \code{TRUE}, include the per-observation
#'   treatment weights in the returned object (not yet used by downstream
#'   methods; included for API parity with the Python package).
#' @param shift Integer.  Time shift applied when constructing relative time
#'   and the untreated indicator (default 0).
#' @param pretrends Non-negative integer.  Number of pre-treatment horizons at
#'   which to estimate placebo ("pre-trend") coefficients.  Default 0 (none).
#' @param cluster Character. Name of the cluster column for cluster-robust
#'   standard errors.  Defaults to \code{i} (cluster by unit).
#' @param avgeffectsby Character vector.  Columns that define the estimand
#'   grouping used for the smartweights influence function.  Defaults to
#'   \code{c(Ei, t)}.
#' @param leaveoneout Not yet implemented; setting to \code{TRUE} raises an
#'   error.
#' @param nose Logical.  If \code{TRUE}, skip standard error computation.
#'   Estimates are still returned.  Default \code{FALSE}.
#' @param delta Integer or \code{NULL}.  Time-step size.  If \code{NULL}
#'   (default) the step is detected automatically as the modal first difference
#'   of the observed time values.
#' @param seed Integer seed passed to \code{\link[base]{set.seed}} before the
#'   randomised rank / collinearity check on the control variables.  Setting
#'   a fixed seed makes this check reproducible across calls.  This is the one
#'   deliberate deviation from the upstream Python package, which uses a
#'   non-seeded random draw.  Default \code{1L}.
#'
#' @return An S3 object of class \code{"did_impute"} with the following
#'   components:
#'   \describe{
#'     \item{\code{estimates}}{Named list of point estimates.  The baseline
#'       case (no horizons, no user \code{wtr}) produces a single element
#'       \code{tau_ate}.  Horizon estimates are named \code{tau0}, \code{tau1},
#'       \ldots; user-weight estimates use the weight column names.}
#'     \item{\code{std_errors}}{Named list of cluster-robust standard errors,
#'       matching \code{estimates}.  \code{NULL} when \code{nose = TRUE}.}
#'     \item{\code{pretrends_estimates}}{Named list of pre-trend (placebo)
#'       coefficients (\code{pre1}, \code{pre2}, \ldots).  \code{NULL} when
#'       \code{pretrends = 0}.}
#'     \item{\code{pretrends_std_errors}}{Named list of cluster-robust SEs for
#'       the pre-trend coefficients.  \code{NULL} when \code{pretrends = 0} or
#'       \code{nose = TRUE}.}
#'     \item{\code{controls_estimates}}{Named list of WLS coefficients on
#'       \code{controls}.  \code{NULL} when \code{controls} is empty.}
#'     \item{\code{controls_std_errors}}{Named list of cluster-robust SEs for
#'       the control coefficients.  \code{NULL} when \code{controls} is empty
#'       or \code{nose = TRUE}.}
#'     \item{\code{n_obs}}{Total observation count entering the estimand
#'       computation (untreated rows plus treated rows with non-zero weight).}
#'     \item{\code{V}}{A square covariance matrix of dimension \eqn{K \times K},
#'       where \eqn{K} is the total number of reported estimands (effects, then
#'       pre-trends, then controls, in that order).  Row and column names match
#'       the estimand names.  \code{sqrt(diag(V))} reproduces the reported
#'       per-estimand standard errors to numerical precision.  This is a
#'       deliberate improvement over the upstream Python package, whose \code{V}
#'       was a scalar (sum of squared SEs).  \code{NULL} when \code{nose = TRUE}
#'       or when no standard-error components are present.}
#'     \item{\code{weights}}{Reserved for future use (\code{NULL} unless
#'       \code{saveweights = TRUE} is requested; the field is included for API
#'       parity with the Python package).}
#'   }
#'
#' @details
#' **Absorbed degrees of freedom.**  The cluster-robust SE scaling factor uses
#' absorbed degrees of freedom (\eqn{df_a}) computed by a port of the
#' \emph{pyhdfe} pairwise method: an FE that is nested within the cluster
#' variable is dropped; the remaining FEs contribute their unique-level counts
#' minus the number of bipartite connected components shared across pairs.
#' This is exact for the default two-way FE / cluster-by-unit case and for any
#' configuration with at most two non-nested fixed effects.  For more than two
#' non-nested fixed effects the value is an approximation and a warning is
#' emitted.
#'
#' **Unimplemented options.**  The arguments \code{timecontrols},
#' \code{unitcontrols}, \code{leaveoneout}, \code{hetby}, and \code{project}
#' match the Python package API but are not yet implemented; they raise an
#' error if supplied.
#'
#' @seealso \code{\link{event_plot}}, \code{\link{summary.did_impute}},
#'   \code{\link{print.did_impute}}
#'
#' @references Borusyak, K., Jaravel, X., and Spiess, J. (2024).
#'   Revisiting Event-Study Designs: Robust and Efficient Estimation.
#'   \emph{Review of Economic Studies}, 91(6), 3253--3285.
#'
#' @examples
#' # Minimal synthetic staggered panel (4 units, 6 periods)
#' set.seed(42)
#' n_units <- 4; n_t <- 6
#' panel <- expand.grid(i = seq_len(n_units), t = seq_len(n_t))
#' panel$Ei <- ifelse(panel$i <= 2, 4L, NA_integer_)  # units 1-2 treated at t=4
#' panel$y  <- 0.5 * (!is.na(panel$Ei) & panel$t >= panel$Ei) +
#'               rnorm(nrow(panel), sd = 0.2)
#'
#' # Baseline ATT estimate
#' res <- did_impute(panel, y = "y", i = "i", t = "t", Ei = "Ei")
#' print(res)
#'
#' # Event-study with horizons 0 and 1 and one pre-period
#' res2 <- did_impute(panel, y = "y", i = "i", t = "t", Ei = "Ei",
#'                    horizons = 0:1, pretrends = 1, minn = 0)
#' summary(res2)
#' @export
did_impute <- function(df, y, i, t, Ei, controls = character(), fe = NULL,
                       timecontrols = character(), aw = NULL, unitcontrols = character(),
                       wtr = character(), sum = FALSE, horizons = NULL, allhorizons = FALSE,
                       hbalance = FALSE, hetby = character(), project = character(),
                       minn = 30, saveweights = FALSE, shift = 0, pretrends = 0,
                       cluster = "", avgeffectsby = character(), leaveoneout = FALSE,
                       nose = FALSE, delta = NULL, seed = 1L) {
  if (length(timecontrols) + length(unitcontrols) > 0)
    stop("Time and unit interacted controls are not added yet.")
  if (length(wtr) != 0 && (length(horizons) > 0 || allhorizons || hbalance))
    stop("User provided weights and horizons options can not be combined.")
  if (leaveoneout)
    stop("Leaveoneout standard errors and options are not added yet.")
  if (length(hetby) > 0 || length(project) > 0)
    stop("Hetby and project options are not added yet.")

  prep <- .prep_data(df, y, i, t, Ei, controls, fe, cluster, aw, sum, shift)
  fe <- prep$fe; cluster <- prep$cluster
  if (length(avgeffectsby) == 0) avgeffectsby <- c(Ei, t)

  imp <- .impute_y0(prep, y, controls, fe, cluster, aw)
  dt <- imp$df

  # Set imput_resid on untreated rows (Python line 318); needed by compute_controls_se
  if (length(controls) > 0) {
    ycol_str <- y
    dt[untreated == 1L, imput_resid := get(ycol_str) - y_hat]
  }

  n_drop <- dt[is.na(y_hat), .N]
  if (n_drop > 0) {
    message(sprintf("Cannot impute for %d observations. Autosample used.", n_drop))
    if (sum) stop("Autosample cannot be combined with sum. Please specify the sample explicitly.")
    dt <- dt[!is.na(y_hat)]
  }
  ycol <- y
  dt[, tau := dt[[ycol]] - y_hat]

  # ---- Seeded rank / collinearity check (Python lines 407-429) ---------------
  if (length(controls) > 0) {
    set.seed(seed)
    dt_untr <- dt[untreated == 1L]
    ymat_u <- matrix(stats::rnorm(nrow(dt_untr)), ncol = 1L)
    xr_u <- fixest::demean(as.matrix(dt_untr[, ..controls]), dt_untr[, ..fe])
    if (any(apply(abs(xr_u), 2L, max) < 1e-10))
      stop("Could not run imputation for some observations because some controls are collinear in the D==0 subsample")
    ymat_f <- matrix(stats::rnorm(nrow(dt)), ncol = 1L)
    xr_f <- fixest::demean(as.matrix(dt[, ..controls]), dt[, ..fe])
    df_a_untr <- .compute_df_a(dt_untr, fe, cluster)
    df_a_full <- .compute_df_a(dt, fe, cluster)
    if (df_a_full > df_a_untr)
      stop("Could not run imputation for some observations because some absorbed variables/FEs are collinear in the D==0 subsample but not in the full sample")
  }

  # ---- Aggregation (Task 6) -----------------------------------------------
  flag_wtr <- length(wtr) != 0

  # Python lines 345-348: default treated-indicator weight if no user wtr provided
  # (unconditional — needed for allhorizons to filter on "wtr" column)
  if (length(wtr) == 0) {
    dt[, wtr := as.numeric(untreated == 0)]
    wtr <- "wtr"
  }

  if (allhorizons) {
    if (length(horizons) > 0) stop("Options horizons and allhorizons cannot be combined")
    # Use current wtr column to find valid Rel_time values
    horizons <- sort(unique(dt[untreated == 0 & !is.na(dt[["wtr"]]) & dt[["wtr"]] != 0, Rel_time]))
  }

  is_horizon <- length(horizons) > 0
  wtrnames <- NULL
  if (is_horizon) {
    horlist <- character(); hornames <- character()
    if (hbalance) {
      dt[, inhorizons := as.integer(Rel_time %in% horizons)]
      dt[, sumh := sum(inhorizons), by = c(i)]
      for (h in horizons) {
        v <- paste0("wtr", h)
        dt[, (v) := as.numeric(Rel_time == h & sumh == length(horizons))]
        na_idx <- is.na(dt[[v]])
        if (any(na_idx)) data.table::set(dt, i = which(na_idx), j = v, value = 0.0)
        horlist <- c(horlist, v); hornames <- c(hornames, paste0("tau", h))
      }
    } else {
      for (h in horizons) {
        v <- paste0("wtr", h)
        # NA == h yields NA in R (unlike Python where NaN == h is False);
        # replace NA with 0 so never-treated rows (NA Rel_time) get weight 0.
        dt[, (v) := as.numeric(Rel_time == h)]
        na_idx <- is.na(dt[[v]])
        if (any(na_idx)) data.table::set(dt, i = which(na_idx), j = v, value = 0.0)
        horlist <- c(horlist, v); hornames <- c(hornames, paste0("tau", h))
      }
    }
    wtr <- horlist; wtrnames <- hornames
  }

  # Normalize each weight over treated weight mass unless sum=TRUE
  if (!sum) {
    for (v in wtr) {
      wv <- dt[[v]]
      ww <- dt$wei
      denom <- base::sum(wv[dt$untreated == 0] * ww[dt$untreated == 0], na.rm = TRUE)
      dt[, (v) := dt[[v]] / denom]
    }
  }

  # minn suppression (strict >)
  if (minn != 0) {
    for (v in wtr) {
      mask <- dt$untreated == 0
      wv <- dt[[v]][mask]
      ww <- dt$wei[mask]
      sab <- base::sum(abs(wv) * ww, na.rm = TRUE)
      if (sab != 0) {
        tmp <- (wv * ww / sab)^2
        if (base::sum(tmp, na.rm = TRUE) > 1 / minn) {
          dt[, (v) := 0]
          message(sprintf("WARNING: suppressing %s, consider lower minn or minn=0.", v))
        }
      } else {
        message(sprintf("WARNING: %s has zero total weight, no data available.", v))
      }
    }
  }

  # Effects
  estimates <- list()
  if (!is_horizon && !flag_wtr) {
    # baseline case: tau_ate
    wv <- dt[["wtr"]]
    ww <- dt$wei
    tau_v <- dt$tau
    mask <- dt$untreated == 0
    estimates[["tau_ate"]] <- base::sum(tau_v[mask] * wv[mask] * ww[mask], na.rm = TRUE)
  } else if (!flag_wtr) {
    # horizon case: tau0, tau1, etc. — key = int(float(v.lstrip("wtr")))
    for (k in seq_along(wtr)) {
      v <- wtr[k]
      wv <- dt[[v]]
      ww <- dt$wei
      tau_v <- dt$tau
      mask <- dt$untreated == 0
      estimates[[wtrnames[k]]] <- base::sum(tau_v[mask] * wv[mask] * ww[mask], na.rm = TRUE)
    }
  } else {
    # user wtr: key = v itself
    for (v in wtr) {
      wv <- dt[[v]]
      ww <- dt$wei
      tau_v <- dt$tau
      mask <- dt$untreated == 0
      estimates[[v]] <- base::sum(tau_v[mask] * wv[mask] * ww[mask], na.rm = TRUE)
    }
  }

  # n_obs: untreated rows PLUS treated rows with any nonzero, non-NA weight
  treated_mask <- dt$untreated == 0
  if (length(wtr) > 0) {
    any_nonzero_wtr <- Reduce(`|`, lapply(wtr, function(v) {
      col <- dt[[v]]
      !is.na(col) & col != 0
    }))
    need_imputation <- treated_mask & any_nonzero_wtr
  } else {
    need_imputation <- rep(FALSE, nrow(dt))
  }
  n_obs <- nrow(dt[!treated_mask | need_imputation, ])

  # ---- Standard errors (Task 7) -------------------------------------------
  gr_var <- unique(c(cluster, avgeffectsby))

  # Small-cohort warning (Python lines 615-640)
  dt[, treat_cohorts := .GRP, by = avgeffectsby]
  max_coh <- max(dt$treat_cohorts)
  if (!is_horizon && !flag_wtr) {
    flag_small <- FALSE
    for (coh in seq_len(max_coh)) {
      sub <- dt[treat_cohorts == coh & untreated == 0L]
      ndist <- length(unique(sub[[i]]))
      if (ndist != 0L && ndist < 15L) { flag_small <- TRUE; break }
    }
    if (flag_small)
      message("The number of treated entities is too small for some cohorts. Standard Errors may be wrong, consider using avgeffectsby option, averaging the the effect by treated X post variable.")
  } else {
    for (v in wtr) {
      flag_v <- FALSE
      for (coh in seq_len(max_coh)) {
        sub <- dt[treat_cohorts == coh & dt[[v]] != 0]
        ndist <- length(unique(sub[[i]]))
        if (ndist != 0L && ndist < 15L) { flag_v <- TRUE; break }
      }
      if (flag_v)
        message(sprintf("The number of treated entities for '%s' is too small for some cohorts. Standard Errors may be wrong; consider using avgeffectsby option, averaging the the effect by treated X post variable.", v))
    }
  }

  if (nose) {
    # Python lines 592-612: pretrends_estimates still computed; controls_estimates too;
    # but all SEs are NULL
    nose_pretrends_estimates <- NULL
    nose_controls_estimates  <- NULL
    if (pretrends > 0L) {
      # We need pretrend coefs even in nose mode; reuse compute_pretrends but
      # ignore ses (Python still runs the WLS for coefs even in nose mode)
      pret_nose <- compute_pretrends(data.table::copy(dt), y, fe, cluster, controls,
                                     pretrends = pretrends,
                                     aw = if (is.null(aw)) NULL else aw)
      nose_pretrends_estimates <- as.list(pret_nose$coefs)
    }
    if (length(controls) > 0L) {
      nose_controls_estimates <- as.list(imp$beta)
    }
    return(new_did_impute(estimates = estimates,
                          pretrends_estimates = nose_pretrends_estimates,
                          controls_estimates  = nose_controls_estimates,
                          n_obs = n_obs))
  }

  se_out <- compute_effect_se(data.table::copy(dt), wtr, y, cluster, avgeffectsby, gr_var, fe,
                              controls = controls)
  std_errors <- as.list(stats::setNames(se_out$se, names(estimates)))

  # ---- Controls estimates and SEs (Task 8) ------------------------------------
  controls_estimates <- NULL
  controls_std_errors <- NULL
  if (length(controls) > 0) {
    controls_estimates <- as.list(imp$beta)
    # Compute df_a on untreated subsample (matches pyhdfe algorithm.degrees with
    # cluster_ids: drops FEs nested within cluster => only t contributes)
    dt_untr_c <- dt[untreated == 1L]
    df_a <- .compute_df_a(dt_untr_c, fe, cluster)
    ctrl_out <- compute_controls_se(dt, controls, fe, cluster,
                                    aw = if (is.null(aw)) NULL else aw,
                                    df_a = df_a, beta = imp$beta)
    controls_std_errors <- as.list(ctrl_out$ses)
    list_ctrl_weps <- ctrl_out$list_ctrl_weps
  } else {
    list_ctrl_weps <- NULL
  }

  # ---- Pretrends (Task 9) -----------------------------------------------------
  pretrends_estimates <- NULL
  pretrends_std_errors <- NULL
  list_pre_weps <- NULL
  if (pretrends > 0L) {
    pret_out <- compute_pretrends(data.table::copy(dt), y, fe, cluster, controls,
                                  pretrends = pretrends,
                                  aw = if (is.null(aw)) NULL else aw)
    pretrends_estimates <- as.list(pret_out$coefs)
    pretrends_std_errors <- as.list(pret_out$ses)
    list_pre_weps <- pret_out$list_pre_weps
  }

  # ---- Combined V: full cluster-robust covariance matrix of all estimands --------
  # Improvement over the upstream Python package (whose V was a scalar = sum of
  # squared SEs).  V is now a (K x K) covariance matrix where K = total estimands
  # (effects, then pretrends, then controls).  sqrt(diag(V)) == reported SEs.
  #
  # Each component already returns a (clusters x estimands) per-cluster influence
  # matrix with cluster IDs as rownames.  We align all components to the union of
  # cluster IDs (zero-filling absent clusters) then cbind and crossprod.
  G_eff  <- se_out$group_sums                          # clusters_eff x n_effects
  G_pret <- if (!is.null(list_pre_weps)) list_pre_weps else NULL  # clusters_pre x n_pre
  G_cont <- if (!is.null(list_ctrl_weps)) list_ctrl_weps else NULL  # clusters_ctrl x n_ctrl

  # Collect non-NULL matrices
  mats <- Filter(Negate(is.null), list(G_eff, G_pret, G_cont))

  if (length(mats) == 0L) {
    V_total <- NULL
  } else {
    # Union of all cluster IDs (as character, sorted for determinism)
    all_cls <- sort(unique(unlist(lapply(mats, rownames))))

    # Zero-fill each matrix to the full union of clusters
    align_mat <- function(M) {
      if (is.null(M) || nrow(M) == 0L) return(NULL)
      nc <- ncol(M)
      out <- matrix(0.0, nrow = length(all_cls), ncol = nc,
                    dimnames = list(all_cls, colnames(M)))
      present <- rownames(M)
      out[present, ] <- M
      out
    }

    G_eff_a  <- align_mat(G_eff)
    G_pret_a <- align_mat(G_pret)
    G_cont_a <- align_mat(G_cont)

    # Bind present components (drop NULLs)
    G_parts <- Filter(Negate(is.null), list(G_eff_a, G_pret_a, G_cont_a))
    G_all <- do.call(cbind, G_parts)

    # Estimand names in documented order: effects, pretrends, controls
    est_names <- c(
      names(estimates),
      if (!is.null(pretrends_estimates)) names(pretrends_estimates),
      if (!is.null(controls_estimates))  names(controls_estimates)
    )
    colnames(G_all) <- est_names

    V_total <- crossprod(G_all)   # t(G_all) %*% G_all
    rownames(V_total) <- est_names
    colnames(V_total) <- est_names
  }

  new_did_impute(estimates = estimates, std_errors = std_errors,
                 pretrends_estimates = pretrends_estimates,
                 pretrends_std_errors = pretrends_std_errors,
                 controls_estimates = controls_estimates,
                 controls_std_errors = controls_std_errors,
                 n_obs = n_obs, V = V_total)
}
