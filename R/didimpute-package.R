#' @keywords internal
"_PACKAGE"

#' @importFrom stats setNames
NULL

## Suppress R CMD check "no visible binding for global variable" notes that
## arise from data.table NSE column references inside [.data.table calls,
## and for the data.table special symbols `.` and `N`.
utils::globalVariables(c(
  # data.table specials
  ".", "N",
  # column names used in NSE inside [.data.table
  "untreated", "wei", "Rel_time", "y_hat", "tau", "wtr",
  "inhorizons", "sumh", "treat_cohorts", "combined", "imput_resid",
  "ctrlweight", "product", "clusterweight", "smartdenom",
  "avg_taus", "sumw",
  # ..col references (data.table column-by-reference)
  "..fe", "..controls", "..other_controls", "..all_xvars", "..rhs_vars",
  # temporary columns used inside compute_effect_se
  ".wtw_", ".tw_", ".tp_", ".cl_", ".prod_",
  # ggplot2 column names referenced in aes()
  "x", "coef", "err", "group"
))
