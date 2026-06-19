#' Recover individual fixed-effect levels by alternating weighted means
#'
#' Faithful port of Python `recover_fixed_effects_iterative`. Decomposes a
#' combined fixed-effects column into per-dimension level effects so they can be
#' extrapolated to treated observations.
#'
#' @param df data.frame containing `fe_cols`, `combined_col`, `weight_col`.
#' @param fe_cols character vector of fixed-effect column names.
#' @param combined_col name of the combined-effect column.
#' @param weight_col name of the weight column.
#' @param max_iter,tol convergence controls.
#' @return named list; `out[[fe]]` is a named numeric vector keyed by level.
#' @keywords internal
recover_fe <- function(df, fe_cols, combined_col, weight_col,
                       max_iter = 100L, tol = 1e-6) {
  keys <- lapply(fe_cols, function(col) as.character(df[[col]]))
  names(keys) <- fe_cols
  combined <- df[[combined_col]]
  w <- df[[weight_col]]

  effects <- lapply(fe_cols, function(col) {
    lv <- unique(keys[[col]])
    stats::setNames(numeric(length(lv)), lv)
  })
  names(effects) <- fe_cols

  for (iter in seq_len(max_iter)) {
    prev <- effects
    for (cur in fe_cols) {
      others <- setdiff(fe_cols, cur)
      if (length(others)) {
        other_sum <- Reduce(`+`, lapply(others, function(o) effects[[o]][keys[[o]]]))
      } else {
        other_sum <- 0
      }
      diff <- combined - other_sum
      num <- tapply(diff * w, keys[[cur]], sum)
      den <- tapply(w, keys[[cur]], sum)
      effects[[cur]] <- c(num / den)          # c() drops tapply's dim; keeps names
    }
    max_change <- max(vapply(fe_cols, function(col) {
      max(abs(effects[[col]] - prev[[col]][names(effects[[col]])]))
    }, numeric(1)))
    if (is.finite(max_change) && max_change < tol) break
  }
  effects
}
