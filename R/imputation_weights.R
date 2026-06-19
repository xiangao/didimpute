#' @import data.table
NULL

#' Orthogonalize weight columns against one variable/FE stratum
#'
#' Faithful port of Python `update_weights` (did_imputation.py lines 62-95).
#' Modifies `dt` in-place and returns it.
#'
#' @param dt data.table
#' @param varlist column name to project against (character scalar); NULL = intercept
#' @param w character vector of weight columns to update
#' @param wei name of observation-weight column
#' @param d name of untreated indicator column
#' @param denom name of denominator column
#' @param by character vector of grouping columns (empty = global sum)
#' @keywords internal
update_weights <- function(dt, varlist = NULL, w = character(), wei = "wei",
                           d = "untreated", denom = "", by = character()) {
  drop_const <- FALSE
  if (is.null(varlist)) {
    data.table::set(dt, j = "_const", value = 1.0)
    varlist <- "_const"
    drop_const <- TRUE
  }

  # Extract vectors outside data.table j to avoid column-name shadowing
  # (data.table resolves bare names against column names first, so a parameter
  #  named "wei" would find the "wei" column as a vector, not the string "wei")
  for (wj in w) {
    wei_vec  <- dt[[wei]]
    wj_vec   <- dt[[wj]]
    var_vec  <- dt[[varlist]]
    den_vec  <- dt[[denom]]
    d_vec    <- dt[[d]]

    if (length(by)) {
      # sum within each by-group level
      grp <- dt[[by]]   # single FE column (by is length-1 for our use)
      # compute per-group sums then broadcast
      grp_keys <- unique(grp)
      sumw_vec <- numeric(nrow(dt))
      for (gk in grp_keys) {
        idx <- grp == gk
        sumw_vec[idx] <- sum(wei_vec[idx] * wj_vec[idx] * var_vec[idx])
      }
      data.table::set(dt, j = "sumw", value = sumw_vec)
    } else {
      sumw_val <- sum(wei_vec * wj_vec * var_vec)
      data.table::set(dt, j = "sumw", value = sumw_val)
    }

    # update: only untreated rows where denom != 0 and not NA
    cond <- d_vec == 1L & !is.na(den_vec) & den_vec != 0
    # new wj value: wj - sumw * varlist / denom
    sumw_v <- dt[["sumw"]]
    new_wj <- wj_vec
    new_wj[cond] <- wj_vec[cond] - sumw_v[cond] * var_vec[cond] / den_vec[cond]
    data.table::set(dt, j = wj, value = new_wj)
  }
  data.table::set(dt, j = "sumw", value = NULL)
  if (drop_const) data.table::set(dt, j = "_const", value = NULL)
  dt
}


#' Compute BJS imputation weights by iterative orthogonalization
#'
#' Faithful port of Python `imputation_weights` (did_imputation.py lines 100-143).
#' Adds `copy<w>` columns to `dt` (in-place style) and returns `dt`.
#'
#' @param dt data.table (modified in-place)
#' @param wei name of observation-weight column
#' @param fe character vector of fixed-effect column names
#' @param wtr character vector of treatment-weight columns
#' @param controls character vector of control variable names
#' @param tol convergence tolerance
#' @param maxit maximum iterations
#' @keywords internal
imputation_weights <- function(dt, wei = "wei", fe, wtr = character(),
                               controls = character(), tol = 1e-6, maxit = 1000L) {
  # Step 1: copy wtr columns to copy<w>
  weights <- character()
  for (w in wtr) {
    v <- paste0("copy", w)
    data.table::set(dt, j = v, value = dt[[w]])
    weights <- c(weights, v)
  }

  # Step 2: demean controls on untreated rows
  wei_vec <- dt[[wei]]
  untr_mask <- dt[["untreated"]] == 1L
  for (v in controls) {
    v_vec <- dt[[v]]
    wm <- stats::weighted.mean(v_vec[untr_mask], wei_vec[untr_mask])
    dm_name <- paste0("dm_", v)
    data.table::set(dt, j = dm_name, value = v_vec - wm)
    dm_vec <- dt[[dm_name]]
    sumval <- sum(wei_vec[untr_mask] * dm_vec[untr_mask]^2)
    data.table::set(dt, j = paste0("denom_", v), value = sumval)
  }

  # Step 3: FE denominators — sum of wei within FE level ON UNTREATED ROWS ONLY
  # Python: df[df["untreated"]==1].groupby(f)["wei"].transform('sum')
  # Untreated rows get the per-level sum; treated rows get NA (excluded by denom != 0 check).
  for (f in fe) {
    dn <- paste0("weight", f)
    f_vec   <- dt[[f]]
    # initialize all NA
    dn_vec <- rep(NA_real_, nrow(dt))
    # compute per-level sums from untreated rows
    f_untr  <- f_vec[untr_mask]
    w_untr  <- wei_vec[untr_mask]
    grp_keys <- unique(f_untr)
    for (gk in grp_keys) {
      idx_untr <- f_untr == gk
      s        <- sum(w_untr[idx_untr])
      # assign to all untreated rows with this FE level
      idx_all  <- untr_mask & (f_vec == gk)
      dn_vec[idx_all] <- s
    }
    data.table::set(dt, j = dn, value = dn_vec)
  }

  # Step 4: iterate until convergence
  it <- 0L
  while (it < maxit) {
    # capture current copy<w> values on untreated rows (as a named list)
    w_copy <- lapply(weights, function(wc) dt[[wc]][untr_mask])
    names(w_copy) <- weights

    for (cont in controls) {
      dt <- update_weights(dt, varlist = paste0("dm_", cont), w = weights,
                           wei = wei, d = "untreated",
                           denom = paste0("denom_", cont), by = character())
    }

    for (f in fe) {
      dt <- update_weights(dt, varlist = NULL, w = weights,
                           wei = wei, d = "untreated",
                           denom = paste0("weight", f), by = f)
    }

    # check convergence: which weights still changing (L1 norm > tol)?
    newkeep <- character()
    for (wc in weights) {
      cur    <- dt[[wc]][untr_mask]
      sumdif <- sum(abs(w_copy[[wc]] - cur), na.rm = TRUE)
      if (sumdif > tol) newkeep <- c(newkeep, wc)
    }
    weights <- newkeep
    it <- it + 1L
    if (!length(weights)) break
  }

  dt
}
