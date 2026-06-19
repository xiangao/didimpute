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

  # default wtr = treated indicator, normalized over treated wei
  flag_wtr <- length(wtr) != 0
  if (!flag_wtr && length(horizons) == 0 && !allhorizons) {
    dt[, wtr := as.numeric(untreated == 0)]
    if (!sum) dt[, wtr := wtr / base::sum(dt[untreated == 0, wtr * wei])]
    ate <- dt[untreated == 0, base::sum(tau * wtr * wei)]
    estimates <- list(tau_ate = ate)
  } else {
    stop("horizons / wtr path arrives in Task 6")
  }
  new_did_impute(estimates = estimates, n_obs = dt[, .N])
}
