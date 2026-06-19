#' Event-study plot for a did_impute result
#'
#' Produces a ggplot2 event-study figure from a \code{did_impute} result object
#' or from explicitly supplied coefficient and standard-error lists.
#' Pre-trend estimates keyed \code{pre<k>} are placed at \eqn{x = -k}; horizon
#' estimates keyed \code{tau<k>} at \eqn{x = +k}.  Confidence intervals use the
#' critical value \eqn{z_{1-\alpha/2}}.
#'
#' @param results_obj A \code{did_impute} result object, or \code{NULL}.
#' @param pretrends Named list of pre-trend estimates (keys \code{pre<k>}).
#'   Used only when \code{results_obj} is \code{NULL}.
#' @param pretrends_std Named list of pre-trend standard errors.
#' @param effects Named list of effect estimates (keys \code{tau<k>}).
#' @param effects_std Named list of effect standard errors.
#' @param plot_type \code{"rcap"} (error bars, default) or \code{"rarea"}
#'   (shaded ribbon).
#' @param significance_level Significance level for CIs (default 0.05 → 95\%).
#' @param together Logical; if \code{TRUE} pre-trends and effects are combined
#'   into a single series (default \code{FALSE}).
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param title Optional plot title.
#' @param pretrends_color Colour for pre-trend series (default \code{"blue"}).
#' @param effects_color Colour for effects series (default \code{"red"}).
#' @param ... Currently unused.
#'
#' @return A \code{ggplot} object.  The underlying data frame (accessible as
#'   \code{p$data}) contains columns \code{x}, \code{coef}, \code{err}, and
#'   \code{group}.
#'
#' @examples
#' \dontrun{
#' res <- did_impute(df, y = "y", i = "i", t = "t", Ei = "Ei",
#'                   horizons = 0:2, pretrends = 2)
#' event_plot(res)
#' event_plot(res, plot_type = "rarea", together = TRUE)
#' }
#' @export
event_plot <- function(results_obj = NULL,
                       pretrends = NULL, pretrends_std = NULL,
                       effects = NULL, effects_std = NULL,
                       plot_type = c("rcap", "rarea"),
                       significance_level = 0.05,
                       together = FALSE,
                       xlab = "Time Relative to Treatment",
                       ylab = "Coefficient",
                       title = NULL,
                       pretrends_color = "blue",
                       effects_color = "red",
                       ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for event_plot(). ",
         "Install it with install.packages('ggplot2').")
  }

  plot_type <- match.arg(plot_type)

  # Extract from results object or use supplied lists
  if (!is.null(results_obj)) {
    pre <- results_obj$pretrends_estimates
    pse <- results_obj$pretrends_std_errors
    eff <- results_obj$estimates
    ese <- results_obj$std_errors
  } else {
    pre <- pretrends
    pse <- pretrends_std
    eff <- effects
    ese <- effects_std
  }

  # Critical value: z_{1 - alpha/2}
  cv <- stats::qnorm(1 - significance_level / 2)

  # Build a tidy data.frame for one series.
  # sign = -1 for pre-trends (pre<k> -> x = -k), +1 for effects (tau<k> -> x = +k)
  .mk <- function(est, se, sign, grp) {
    if (is.null(est) || length(est) == 0L) return(NULL)
    k_str  <- names(est)
    # Keep only keys whose suffix after the prefix is a non-negative integer
    # (e.g. "pre2", "tau0"). This excludes aggregate keys like "tau_ate".
    suffix <- sub("^(pre|tau)", "", k_str)
    valid  <- grepl("^[0-9]+$", suffix)
    k_str  <- k_str[valid]
    if (length(k_str) == 0L) return(NULL)
    x_vals <- sign * as.integer(sub("^(pre|tau)", "", k_str))
    coef_v <- as.numeric(unlist(est[k_str]))
    err_v  <- if (!is.null(se) && length(se) > 0L)
                cv * as.numeric(unlist(se[k_str]))
              else
                rep(NA_real_, length(x_vals))
    data.frame(x = x_vals, coef = coef_v, err = err_v,
               group = grp, stringsAsFactors = FALSE)
  }

  df_pre <- .mk(pre, pse, -1L, "Pre-trends")
  df_eff <- .mk(eff, ese, +1L, "Effects")

  if (together) {
    # Merge into one series; use effects colour/label
    df_all <- rbind(df_pre, df_eff)
    if (!is.null(df_all)) {
      # Zero-fill missing-SE side (Python parity)
      # Identify which side has errors
      pre_has_err  <- !is.null(df_pre) && any(!is.na(df_pre$err))
      eff_has_err  <- !is.null(df_eff) && any(!is.na(df_eff$err))

      # If only one side has errors, zero-fill the other side
      if (pre_has_err && !eff_has_err) {
        # Effects missing SEs: set effects err to 0
        df_all$err[is.na(df_all$err)] <- 0.0
      } else if (!pre_has_err && eff_has_err) {
        # Pretrends missing SEs: set pretrends err to 0
        df_all$err[is.na(df_all$err)] <- 0.0
      }
      # If both have errors or neither has errors, leave as-is

      df_all$group <- "Effects"
    }
    df <- df_all
    color_map <- c("Effects" = effects_color)
  } else {
    df <- rbind(df_pre, df_eff)
    color_map <- c("Pre-trends" = pretrends_color,
                   "Effects"    = effects_color)
  }

  if (is.null(df) || nrow(df) == 0L) {
    stop("No data to plot: supply pretrends, effects, or a did_impute result.")
  }

  df <- df[order(df$x), ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = coef, color = group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "black", alpha = 0.5) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df$x))) +
    ggplot2::labs(x = xlab, y = ylab, title = title, color = NULL) +
    ggplot2::theme_minimal()

  if (plot_type == "rarea") {
    # Ribbon: fill mapped to group so colours match points/lines
    fill_map <- if (together)
      c("Effects" = effects_color)
    else
      c("Pre-trends" = pretrends_color, "Effects" = effects_color)

    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = coef - err, ymax = coef + err, fill = group),
        alpha = 0.3, color = NA, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = fill_map, drop = FALSE)
  } else {
    # rcap: error bars
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = coef - err, ymax = coef + err),
        width = 0.1)
  }

  p
}
