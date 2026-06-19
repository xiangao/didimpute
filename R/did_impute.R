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

  n_drop <- dt[is.na(y_hat), .N]
  if (n_drop > 0) {
    message(sprintf("Cannot impute for %d observations. Autosample used.", n_drop))
    if (sum) stop("Autosample cannot be combined with sum. Please specify the sample explicitly.")
    dt <- dt[!is.na(y_hat)]
  }
  ycol <- y
  dt[, tau := dt[[ycol]] - y_hat]

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
    return(new_did_impute(estimates = estimates, n_obs = n_obs))
  }

  se_out <- compute_effect_se(data.table::copy(dt), wtr, y, cluster, avgeffectsby, gr_var, fe)
  std_errors <- as.list(stats::setNames(se_out$se, names(estimates)))
  new_did_impute(estimates = estimates, std_errors = std_errors, n_obs = n_obs,
                 V = se_out$V)
}
