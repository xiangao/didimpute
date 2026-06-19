#' @import data.table
NULL

#' Compute cluster-robust standard errors via BJS smartweights
#'
#' Faithful port of Python `compute_effect_se` (did_imputation.py lines 615-711,
#' no-controls path). Operates on a COPY of dt (caller must pass one).
#'
#' For each weight column in `wtr`:
#'   - Fill NaN copy<w> with original w (non-iterated weights)
#'   - Compute smartweights on treated rows
#'   - Residual: y - y_hat (all rows), overridden to tau - avg_taus (treated rows)
#'   - product = resid * copy<w> * wei
#'   - group_sums = sum(product) per cluster
#' Then V = G'G where G = cbind(group_sums); ses = sqrt(diag(V)).
#'
#' @param dt data.table (will be modified; pass a copy)
#' @param wtr character vector of original treatment-weight column names
#' @param y name of outcome column
#' @param cluster name of cluster column
#' @param avgeffectsby character vector of columns defining the estimand group
#' @param gr_var character vector = union of cluster + avgeffectsby
#' @param fe character vector of fixed-effect columns (passed to imputation_weights)
#' @return list with `se` (named numeric), `group_sums` matrix, `V` matrix
#' @keywords internal
compute_effect_se <- function(dt, wtr, y, cluster, avgeffectsby, gr_var, fe) {
  # Run imputation weights iteration (adds copy<w> columns)
  dt <- imputation_weights(dt, wei = "wei", fe = fe, wtr = wtr)

  weights <- paste0("copy", wtr)
  group_sums_list <- list()

  # Pre-extract constant vectors to avoid data.table j-env column-name shadowing
  untreated_vec <- dt[["untreated"]]
  m <- untreated_vec == 0L           # treated mask (logical)
  y_vec    <- dt[[y]]
  y_hat    <- dt[["y_hat"]]
  tau_vec  <- dt[["tau"]]
  wei_vec  <- dt[["wei"]]
  cl_vec   <- dt[[cluster]]

  for (k in seq_along(weights)) {
    w_col  <- weights[k]    # "copy<wtr>"
    bw_col <- wtr[k]        # original "<wtr>"

    # Extract copy<w> vector; fill NA with original <wtr>
    w_vec <- dt[[w_col]]
    bw_vec <- dt[[bw_col]]
    na_mask <- is.na(w_vec)
    if (any(na_mask)) w_vec[na_mask] <- bw_vec[na_mask]
    # Update the column in dt (needed for merge-based groupby)
    data.table::set(dt, j = w_col, value = w_vec)

    # ---- clusterweight: sum(copy<w> * wei) per gr_var group (treated rows) ----
    # Python: df_sub = df[mask]; df_sub["w_times_wei"] = w * wei
    #         clwei_df = df_sub.groupby(gr_var)["w_times_wei"].sum()
    w_times_wei <- w_vec * wei_vec
    # Build gr_var key for grouping — use data.table
    if ("clusterweight" %in% names(dt)) data.table::set(dt, j = "clusterweight", value = NULL)
    # Set w_times_wei temporarily
    data.table::set(dt, j = ".wtw_", value = w_times_wei)
    cl_agg <- dt[m, .(clusterweight = sum(.wtw_)), by = gr_var]
    data.table::set(dt, j = ".wtw_", value = NULL)
    dt <- merge(dt, cl_agg, by = gr_var, all.x = TRUE, sort = FALSE)

    # ---- smartdenom: sum(clusterweight * copy<w> * wei) per avgeffectsby (treated rows) ----
    clw_vec <- dt[["clusterweight"]]
    temp_weighted <- clw_vec * w_vec * wei_vec
    # Note: w_vec, wei_vec indices match dt order; clusterweight was merged so it also matches
    if ("smartdenom" %in% names(dt)) data.table::set(dt, j = "smartdenom", value = NULL)
    data.table::set(dt, j = ".tw_", value = temp_weighted)
    sm_agg <- dt[m, .(smartdenom = sum(.tw_)), by = avgeffectsby]
    data.table::set(dt, j = ".tw_", value = NULL)
    dt <- merge(dt, sm_agg, by = avgeffectsby, all.x = TRUE, sort = FALSE)

    # Re-extract after merge (row order may differ)
    w_vec    <- dt[[w_col]]
    wei_vec  <- dt[["wei"]]
    tau_vec  <- dt[["tau"]]
    y_vec    <- dt[[y]]
    y_hat    <- dt[["y_hat"]]
    untreated_vec <- dt[["untreated"]]
    m        <- untreated_vec == 0L
    cl_vec   <- dt[[cluster]]
    clw_vec  <- dt[["clusterweight"]]
    smd_vec  <- dt[["smartdenom"]]

    # ---- smartweight (treated rows only) ----
    smartweight <- rep(0.0, nrow(dt))
    denom_safe  <- smd_vec
    smartweight[m] <- (clw_vec[m] * w_vec[m] * wei_vec[m]) / denom_safe[m]
    smartweight[is.na(smartweight)] <- 0.0

    # ---- avg_taus: smartweight-weighted mean of tau per avgeffectsby group ----
    # Python: df['temp_product'] = tau * smartweight
    #         df['avg_taus'] = df.groupby(avgeffectsby)['temp_product'].transform('sum')
    # smartweight=0 for untreated => transform broadcasts sum (incl. untreated) per group
    temp_product <- tau_vec * smartweight
    data.table::set(dt, j = ".tp_", value = temp_product)
    dt[, avg_taus := sum(.tp_), by = avgeffectsby]
    data.table::set(dt, j = ".tp_", value = NULL)
    avg_taus <- dt[["avg_taus"]]

    # ---- residuals ----
    # All rows: y - y_hat; then treated rows overridden: tau - avg_taus
    resid_vec <- y_vec - y_hat
    resid_vec[m] <- tau_vec[m] - avg_taus[m]

    # ---- product and cluster sum ----
    product <- resid_vec * w_vec * wei_vec
    # Sum product per cluster (all rows)
    gs_dt <- data.table::data.table(.cl_ = cl_vec, .prod_ = product)
    gs_agg <- gs_dt[, .(s = sum(.prod_)), by = .cl_]
    gs_agg <- gs_agg[order(.cl_)]
    group_sums_list[[bw_col]] <- gs_agg$s
  }

  G <- do.call(cbind, group_sums_list)
  V <- t(G) %*% G
  ses <- sqrt(diag(V))
  list(se = ses, group_sums = G, V = V)
}
