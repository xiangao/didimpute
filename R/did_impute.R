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
  stop("not implemented yet")
}
