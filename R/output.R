#' Print a did_impute result
#'
#' @param x a \code{did_impute} object.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @examples
#' set.seed(1)
#' panel <- expand.grid(i = 1:4, t = 1:6)
#' panel$Ei <- ifelse(panel$i <= 2, 4L, NA_integer_)
#' panel$y  <- 0.3 * (!is.na(panel$Ei) & panel$t >= panel$Ei) + rnorm(nrow(panel), sd = 0.1)
#' res <- did_impute(panel, y = "y", i = "i", t = "t", Ei = "Ei")
#' print(res)
#' @export
print.did_impute <- function(x, ...) {
  cat("<didimpute did_impute result>\n")
  if (!is.null(x$estimates))
    cat("  Effects:", paste(sprintf("%s=%.4g", names(x$estimates), unlist(x$estimates)), collapse = "  "), "\n")
  if (!is.null(x$std_errors))
    cat("  SE     :", paste(sprintf("%s=%.4g", names(x$std_errors), unlist(x$std_errors)), collapse = "  "), "\n")
  if (!is.null(x$pretrends_estimates))
    cat("  Pretrends:", paste(sprintf("%s=%.4g", names(x$pretrends_estimates), unlist(x$pretrends_estimates)), collapse = "  "), "\n")
  if (!is.null(x$controls_estimates))
    cat("  Controls:", paste(sprintf("%s=%.4g", names(x$controls_estimates), unlist(x$controls_estimates)), collapse = "  "), "\n")
  cat("  N obs:", x$n_obs, "\n")
  invisible(x)
}

#' Summarise a did_impute result
#'
#' Returns a tidy \code{data.frame} with one row per estimand (effects,
#' pre-trends, and control coefficients combined).
#'
#' @param object a \code{did_impute} object.
#' @param ... ignored.
#' @return A \code{data.frame} with columns \code{term}, \code{estimate}, and
#'   \code{std.error}.
#' @examples
#' set.seed(1)
#' panel <- expand.grid(i = 1:4, t = 1:6)
#' panel$Ei <- ifelse(panel$i <= 2, 4L, NA_integer_)
#' panel$y  <- 0.3 * (!is.na(panel$Ei) & panel$t >= panel$Ei) + rnorm(nrow(panel), sd = 0.1)
#' res <- did_impute(panel, y = "y", i = "i", t = "t", Ei = "Ei",
#'                   horizons = 0:1, pretrends = 2, minn = 0)
#' summary(res)
#' @export
summary.did_impute <- function(object, ...) {
  terms    <- names(object$estimates)
  ests     <- unlist(object$estimates)
  se_vals  <- if (is.null(object$std_errors)) rep(NA_real_, length(terms)) else unlist(object$std_errors)
  df_out   <- data.frame(term = terms, estimate = ests, std.error = se_vals,
                         row.names = NULL, stringsAsFactors = FALSE)

  # Append pretrends rows if present
  if (!is.null(object$pretrends_estimates)) {
    pre_terms <- names(object$pretrends_estimates)
    pre_ests  <- unlist(object$pretrends_estimates)
    pre_se    <- if (is.null(object$pretrends_std_errors)) rep(NA_real_, length(pre_terms)) else unlist(object$pretrends_std_errors)
    df_out <- rbind(df_out,
                    data.frame(term = pre_terms, estimate = pre_ests, std.error = pre_se,
                               row.names = NULL, stringsAsFactors = FALSE))
  }

  # Append controls rows if present
  if (!is.null(object$controls_estimates)) {
    ctrl_terms <- names(object$controls_estimates)
    ctrl_ests  <- unlist(object$controls_estimates)
    ctrl_se    <- if (is.null(object$controls_std_errors)) rep(NA_real_, length(ctrl_terms)) else unlist(object$controls_std_errors)
    df_out <- rbind(df_out,
                    data.frame(term = ctrl_terms, estimate = ctrl_ests, std.error = ctrl_se,
                               row.names = NULL, stringsAsFactors = FALSE))
  }

  df_out
}
