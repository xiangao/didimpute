#' @import data.table
NULL

#' Compute absorbed degrees of freedom matching pyhdfe pairwise method
#'
#' Mirrors the pyhdfe `algorithm.degrees` logic (used with cluster_ids): an FE
#' is dropped if it is nested within (i.e., its levels are a refinement of) any
#' cluster group.  The remaining FEs contribute their unique level counts, minus
#' the number of connected components across pairs (pairwise method).
#'
#' For the standard case (fe = c(i, t), cluster = i): the i-FE is nested within
#' the i-cluster, so only the t-FE contributes => df_a = uniqueN(t).
#'
#' @param dt data.table restricted to the untreated subsample
#' @param fe character vector of FE column names
#' @param cluster name of cluster column
#' @return integer absorbed degrees of freedom
#' @keywords internal
.compute_df_a <- function(dt, fe, cluster) {
  # For each FE f, check if f is nested within cluster (each cluster contains
  # only one value of f, i.e., within-cluster variance of f is 0).
  is_nested <- vapply(fe, function(f) {
    # A FE is nested within cluster if every cluster has exactly 1 unique level
    tab <- dt[, data.table::uniqueN(get(f)), by = cluster]
    all(tab$V1 == 1L)
  }, logical(1L))

  active_fe <- fe[!is_nested]
  if (length(active_fe) == 0L) return(0L)

  # Each active FE contributes its unique level count
  counts <- vapply(active_fe, function(f) data.table::uniqueN(dt[[f]]), integer(1L))
  degrees <- sum(counts)

  # Pairwise subtraction: for each active FE beyond the first, subtract the
  # number of connected components it shares with prior FEs.
  # For the typical case (single active FE = t), no subtraction needed.
  if (length(active_fe) > 1L) {
    # Warn if more than 2 active FEs: the pairwise max(comp_counts) approach
    # is an approximation of pyhdfe's global bipartite connected-components count,
    # exact only for <=2 fixed effects.
    if (length(active_fe) > 2L) {
      warning("compute_df_a: absorbed degrees of freedom for more than two non-nested fixed effects is approximate; cluster-robust standard errors may be mis-scaled in this configuration.")
    }
    for (j in seq(2L, length(active_fe))) {
      f_j <- active_fe[j]
      prior <- active_fe[seq_len(j - 1L)]
      # connected components: max over prior FEs of the number of components
      # f_j shares with each prior FE.
      # Approximate: use the number of distinct (f_prior, f_j) pairs relative
      # to unique(f_j) — i.e. subtract (unique(f_j) - connected_components).
      # For a balanced panel this equals min(unique(f_j), unique(prior)) - 1 per pair.
      # Use the same formula pyhdfe uses: subtract components shared with each prior.
      comp_counts <- vapply(prior, function(fp) {
        # number of connected components between fp and f_j (graph coloring)
        # For simplicity: n_unique(f_j) - (n_unique(f_j) - n_unique(fp)) when fp is smaller
        # The correct formula uses union-find but for our cases this is simple:
        # components = unique rows of (fp, fj) - unique(fj) + unique(fp) is NOT right.
        # Use rank-based: treat as bipartite graph; components = n_i + n_j - rank(biadj)
        n_i <- data.table::uniqueN(dt[[fp]])
        n_j <- data.table::uniqueN(dt[[f_j]])
        # rank of incidence matrix
        cross <- unique(dt[, c(fp, f_j), with = FALSE])
        # number of edges
        n_edges <- nrow(cross)
        # components = n_i + n_j - rank(incidence)
        # rank(incidence) = n_i + n_j - components  =>  components = n_i + n_j - rank
        # For a connected bipartite graph: rank = n_i + n_j - 1
        # Use spanning-tree approach: components = n_i + n_j - rank via SVD is expensive;
        # use union-find instead
        parents <- seq_len(n_i + n_j)
        find_p <- function(x) { while (parents[x] != x) x <- parents[x]; x }
        i_ids <- as.integer(factor(cross[[fp]]))
        j_ids <- as.integer(factor(cross[[f_j]])) + n_i
        for (e in seq_len(nrow(cross))) {
          ri <- find_p(i_ids[e]); rj <- find_p(j_ids[e])
          if (ri != rj) parents[ri] <- rj
        }
        roots <- vapply(seq_len(n_i + n_j), find_p, integer(1L))
        length(unique(roots))
      }, integer(1L))
      degrees <- degrees - max(comp_counts)
    }
  }
  as.integer(degrees)
}


#' Compute cluster-robust standard errors for pretrend coefficients
#'
#' Faithful port of Python pretrends block (did_imputation.py lines 477-566,
#' 700-704).  Operates on the untreated subsample after dropping singleton FE
#' levels.  For each h=1..pretrends: residualise pretrendvar_h on FE (and other
#' pretrend/control vars); WLS to get preresid; form per-observation preweight;
#' normalise; cluster-sum preweight * preresid * sqrt(dof_adj).  Sandwich V = M'M.
#'
#' @param dt data.table (all rows, including treated — untreated filter applied
#'   inside)
#' @param y name of outcome column
#' @param fe character vector of FE column names
#' @param cluster name of cluster column
#' @param controls character vector of control column names
#' @param pretrends integer number of pretrend horizons (h = 1..pretrends)
#' @param aw name of analytic weight column, or NULL
#' @return list with \code{coefs} (named numeric) and \code{ses} (named numeric),
#'   plus \code{list_pre_weps} (matrix of per-cluster influence vectors)
#' @keywords internal
compute_pretrends <- function(dt, y, fe, cluster, controls, pretrends, aw) {
  # Build pretrendvar_h = 1{Rel_time == -h} on the full dataset (matching Python).
  # NA Rel_time (never-treated units) => 0, not NA (Python NaN == -h is False).
  for (h in seq_len(pretrends)) {
    vname <- paste0("pretrendvar_", h)
    rel_h <- dt[["Rel_time"]]
    dt[, (vname) := as.numeric(!is.na(rel_h) & rel_h == -h)]
  }
  pretrendvars <- paste0("pretrendvar_", seq_len(pretrends))

  # Restrict to untreated and drop singleton FE levels (Python lines 486-491)
  contr <- data.table::copy(dt[untreated == 1L])
  for (f in fe) {
    cnt <- contr[, .N, by = f]
    keep <- cnt[N > 1L][[f]]
    contr <- contr[get(f) %in% keep]
  }

  wcol <- if (is.null(aw)) NULL else contr$wei

  # Residualise [pretrendvars + controls] jointly on FE, with weights (Python 493-503)
  all_xvars <- c(pretrendvars, controls)
  ymat     <- as.matrix(contr[[y]])
  y_resid  <- drop(fixest::demean(ymat, contr[, ..fe], weights = wcol))

  if (length(all_xvars) > 0L) {
    xmat    <- as.matrix(contr[, ..all_xvars])
    X_resid <- fixest::demean(xmat, contr[, ..fe], weights = wcol)
  } else {
    X_resid <- matrix(1.0, nrow(contr), 1L)
  }

  ww <- if (is.null(aw)) rep(1.0, nrow(contr)) else contr$wei

  # Compute df_a using cluster_ids (same as effect SE: FE nested in cluster dropped)
  df_a <- .compute_df_a(contr, fe, cluster)

  # WLS for pretrend coefs (Python line 508-514)
  fit_main <- stats::lm.wfit(x = X_resid, y = y_resid, w = ww)
  coef_all <- fit_main$coefficients  # length = pretrends + length(controls)
  preresid <- y_resid - drop(X_resid %*% coef_all)

  coefs <- stats::setNames(coef_all[seq_len(pretrends)], paste0("pre", seq_len(pretrends)))

  # Standard (non-cluster) SEs from WLS fit — used for skip condition below.
  # Python results.bse is the standard statsmodels WLS SE (same quantity).
  sigma2_main <- if (fit_main$df.residual > 0L)
    sum(fit_main$weights * fit_main$residuals^2) / fit_main$df.residual else 0.0
  XtWX_inv <- tryCatch(solve(crossprod(X_resid * sqrt(ww))), error = function(e) NULL)
  main_bse <- if (!is.null(XtWX_inv))
    sqrt(pmax(0.0, sigma2_main * diag(XtWX_inv))) else rep(0.0, length(coef_all))

  ncl <- contr[, data.table::uniqueN(get(cluster))]
  n_contr <- nrow(contr)

  list_pre_weps <- list()
  for (h in seq_len(pretrends)) {
    pv <- paste0("pretrendvar_", h)
    # Python line 522: skip (set weps=0) only when coef==0 AND se==0 (standard WLS SE).
    if (isTRUE(coefs[[paste0("pre", h)]] == 0) && isTRUE(main_bse[[h]] == 0)) {
      list_pre_weps[[h]] <- rep(0.0, ncl)
      next
    }

    # Demean pretrendvar_h on FE (no cluster_ids — Python uses create(ids) not create(ids, cluster_ids))
    Y_h      <- as.matrix(contr[[pv]])
    Y_h_resid <- drop(fixest::demean(Y_h, contr[, ..fe], weights = wcol))

    # Other pretrend vars + controls as RHS
    # Guard: paste0("pretrendvar_", integer(0)) yields "pretrendvar_" (length 1) not
    # character(0), so we must check the index before pasting.
    others <- setdiff(seq_len(pretrends), h)
    other_pretrendvars <- if (length(others) > 0L) paste0("pretrendvar_", others) else character(0)
    rhs_vars <- c(controls, other_pretrendvars)
    if (length(rhs_vars) > 0L) {
      X_h      <- as.matrix(contr[, ..rhs_vars])
      X_h_resid <- fixest::demean(X_h, contr[, ..fe], weights = wcol)
    } else {
      X_h_resid <- matrix(1.0, nrow(contr), 1L)
    }

    # WLS -> preweight = residual (Python line 549)
    fit_h    <- stats::lm.wfit(x = X_h_resid, y = Y_h_resid, w = ww)
    preweight <- Y_h_resid - drop(X_h_resid %*% fit_h$coefficients)

    # Normalise by sum over rows where pretrendvar_h == 1 (Python lines 551-556)
    wei_vec <- if (is.null(aw)) rep(1.0, nrow(contr)) else contr$wei
    preweight_weighted <- preweight * wei_vec
    mask_h <- contr[[pv]] == 1
    sumpreweight <- sum(preweight_weighted[mask_h])
    if (sumpreweight != 0) {
      preweight_normalized <- preweight_weighted / sumpreweight
    } else {
      preweight_normalized <- rep(0.0, length(preweight_weighted))
    }

    # product = normalized * preresid (Python line 558)
    contr[, product := preweight_normalized * preresid]

    # dof_adj: subtract pretrends in denominator (Python line 561)
    dof_adj <- (n_contr - 1L) / (n_contr - length(controls) - pretrends - df_a + 1L) *
               ncl / (ncl - 1L)

    # per-cluster sum * sqrt(dof_adj) (Python line 564)
    sums_dt <- contr[, .(v = sum(product)), by = cluster]
    sums_dt <- sums_dt[order(get(cluster))]
    list_pre_weps[[h]] <- sums_dt$v * sqrt(dof_adj)
  }

  # Build per-cluster influence matrix (each column = one pretrend horizon)
  M    <- do.call(cbind, list_pre_weps)
  V    <- t(M) %*% M
  ses  <- stats::setNames(sqrt(diag(V)), paste0("pre", seq_len(pretrends)))

  list(coefs = coefs, ses = ses, list_pre_weps = M)
}


#' Compute cluster-robust SEs for control-variable coefficients (influence function)
#'
#' Faithful port of Python \code{compute_controls_se} (did_imputation.py lines
#' 351-404, 691-695).  For each control c: demean c and the other controls on
#' FE (untreated subsample, after dropping singleton FE levels); regress demeaned
#' c on demeaned others; form per-observation ctrlweight from residuals; cluster-
#' sum ctrlweight * imput_resid; sandwich V = M'M.
#'
#' @param dt data.table containing untreated and treated rows (needs imput_resid
#'   on untreated rows)
#' @param controls character vector of control column names
#' @param fe character vector of FE column names
#' @param cluster name of cluster column
#' @param aw name of analytic weight column, or NULL
#' @param df_a absorbed degrees of freedom (from .compute_df_a)
#' @param beta named numeric vector of WLS coefficients for each control
#' @return named numeric vector of cluster-robust SEs, one per control
#' @keywords internal
compute_controls_se <- function(dt, controls, fe, cluster, aw, df_a, beta) {
  contr <- data.table::copy(dt[untreated == 1L])

  # Drop singleton FE levels (value_counts > 1 in Python means keep levels
  # that appear more than once)
  for (f in fe) {
    cnt <- contr[, .N, by = f]
    keep <- cnt[N > 1L][[f]]
    contr <- contr[get(f) %in% keep]
  }

  ncl <- contr[, data.table::uniqueN(get(cluster))]
  wcol <- if (is.null(aw)) NULL else contr$wei

  list_ctrl_weps <- list()

  for (c in controls) {
    # Python: if cont_coefs[c] == 0 and cont_se[c] == 0 -> zeros
    # Brief says: if beta[c] == 0 -> weps = 0
    if (isTRUE(beta[[c]] == 0)) {
      list_ctrl_weps[[c]] <- rep(0.0, ncl)
      next
    }

    other_controls <- setdiff(controls, c)

    # Demean c on FE
    Y_h <- as.matrix(contr[[c]])
    Y_h_resid <- drop(fixest::demean(Y_h, contr[, ..fe], weights = wcol))

    if (length(other_controls) > 0L) {
      X_h <- as.matrix(contr[, ..other_controls])
      X_h_resid <- fixest::demean(X_h, contr[, ..fe], weights = wcol)
    } else {
      X_h_resid <- matrix(1.0, nrow(contr), 1L)
    }

    ww <- if (is.null(aw)) rep(1.0, nrow(contr)) else contr$wei
    fit_h <- stats::lm.wfit(x = X_h_resid, y = Y_h_resid, w = ww)
    resid_h <- Y_h_resid - drop(X_h_resid %*% fit_h$coefficients)

    # dof adjustment (Python line 389)
    n_contr <- nrow(contr)
    dof_adj <- (n_contr - 1L) / (n_contr - length(controls) - df_a + 1L) *
               ncl / (ncl - 1L)

    # ctrlweight
    wei_h <- if (is.null(aw)) rep(1.0, nrow(contr)) else contr$wei
    cw <- resid_h * wei_h
    sp <- sum(cw * contr[[c]])
    if (sp != 0) {
      cw <- cw / sp
    } else {
      cw <- rep(0.0, length(cw))
    }
    contr[, ctrlweight := cw]

    # cluster sums: sum(ctrlweight * imput_resid) per cluster, scaled by sqrt(dof_adj)
    sums_dt <- contr[, .(v = sum(ctrlweight * imput_resid)), by = cluster]
    # order clusters consistently
    sums_dt <- sums_dt[order(get(cluster))]
    list_ctrl_weps[[c]] <- sums_dt$v * sqrt(dof_adj)
  }

  M <- do.call(cbind, list_ctrl_weps)
  V_cont <- t(M) %*% M
  se_cont <- sqrt(diag(V_cont))
  list(ses = stats::setNames(se_cont, controls), list_ctrl_weps = M)
}


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
#' @param controls character vector of control variable names (default empty)
#' @return list with `se` (named numeric), `group_sums` matrix, `V` matrix
#' @keywords internal
compute_effect_se <- function(dt, wtr, y, cluster, avgeffectsby, gr_var, fe,
                              controls = character()) {
  # Run imputation weights iteration (adds copy<w> columns)
  # Pass controls so that weights are also orthogonalized w.r.t. them (Python line 643)
  dt <- imputation_weights(dt, wei = "wei", fe = fe, wtr = wtr, controls = controls)

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

    # Re-extract after first merge (row order may differ)
    w_vec   <- dt[[w_col]]
    wei_vec <- dt[["wei"]]

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
