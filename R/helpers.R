# ============================================================
# helpers.R
#
# Internal helpers used by both the basis path and the GAM path
# of dcame_aipe(). Not exported.
# ============================================================

#' Place interior knots for a spline basis
#' @keywords internal
#' @noRd
place_knots <- function(d_values, n_knots, strategy = "quantile") {
  if (n_knots == 0) return(numeric(0))
  if (strategy == "quantile") {
    probs <- seq(0, 1, length.out = n_knots + 2)[-c(1, n_knots + 2)]
    as.numeric(stats::quantile(d_values, probs = probs, na.rm = TRUE))
  } else {
    rng <- range(d_values, na.rm = TRUE)
    seq(rng[1], rng[2], length.out = n_knots + 2)[-c(1, n_knots + 2)]
  }
}

#' Compute boundary knots
#' @keywords internal
#' @noRd
compute_boundary_knots <- function(d_values) {
  range(d_values, na.rm = TRUE)
}

#' Create a basis matrix (polynomial or B-spline)
#' @keywords internal
#' @noRd
create_basis <- function(d_values, spec, knots_ref = NULL, bknots_ref = NULL) {
  if (spec$type == "poly") {
    B <- stats::poly(d_values, degree = spec$degree, raw = TRUE)
    colnames(B) <- paste0("poly_", 1:spec$degree)
    return(list(basis = B, knots_interior = NULL, bknots = NULL, spec = spec))
  }
  if (is.null(knots_ref)) {
    strategy <- if (!is.null(spec$knot_strategy)) spec$knot_strategy else "quantile"
    knots_interior <- place_knots(d_values, spec$knots, strategy)
    bknots <- compute_boundary_knots(d_values)
  } else {
    knots_interior <- knots_ref
    bknots <- bknots_ref
  }
  if (spec$type %in% c("spline", "bspline")) {
    B <- splines::bs(d_values, knots = knots_interior, degree = spec$degree,
                     Boundary.knots = bknots, intercept = FALSE)
    colnames(B) <- paste0("bs_", 1:ncol(B))
    return(list(basis = B, knots_interior = knots_interior,
                bknots = bknots, spec = spec))
  }
  stop(sprintf("Unknown spec type: '%s'", spec$type))
}

#' Complexity of a specification (used for tie-breaking in CV/BIC)
#' @keywords internal
#' @noRd
spec_complexity <- function(s) {
  if (s$type == "poly") return(s$degree)
  s$knots + s$degree
}

#' Human-readable label for a specification
#' @keywords internal
#' @noRd
spec_label <- function(s) {
  strat <- if (!is.null(s$knot_strategy)) paste0(",", substr(s$knot_strategy, 1, 1)) else ""
  switch(s$type,
         "poly"    = paste0("Poly(", s$degree, ")"),
         "bspline" = ,
         "spline"  = paste0("BS(k=", s$knots, ",d=", s$degree, strat, ")"),
         paste0(s$type, "(k=", s$knots, ")")
  )
}

#' Parse user-supplied spline string "knots,degree"
#' @keywords internal
#' @noRd
parse_spline_spec <- function(s, type = "bspline") {
  p <- as.integer(strsplit(s, ",")[[1]])
  list(type = type, knots = p[1], degree = p[2], knot_strategy = "quantile")
}

# ============================================================
# KERNEL DENSITY HELPERS (for AIPE-targeted CV weights)
# ============================================================

#' Weighted kernel density estimator
#' @keywords internal
#' @noRd
weighted_kde <- function(d_values, wgts = NULL, bw_method = "SJ") {
  if (is.null(wgts)) wgts <- rep(1, length(d_values))
  wgts <- wgts / sum(wgts) * length(d_values)
  bw <- tryCatch(
    stats::bw.SJ(d_values),
    error = function(e) stats::bw.nrd0(d_values)
  )
  function(x_eval) {
    n <- length(d_values)
    sapply(x_eval, function(x0) {
      sum(wgts * stats::dnorm((x0 - d_values) / bw)) / (n * bw)
    })
  }
}

#' Compute AIPE-targeted CV weights
#'
#' For each observation i in group g, returns w_i = f_pool(D_i) / f_own_g(D_i),
#' where f_pool and f_own_g are KDE-estimated densities of D on the full sample
#' and within group g respectively. Weights are capped at quantile(0.98) * max_ratio
#' to bound the influence of low-density tails. Generalizes the binary x in {0,1}
#' case in the paper appendix to K groups (g in 0..K-1).
#' @keywords internal
#' @noRd
compute_targeted_cv_weights <- function(d, g, survey_wgts = NULL,
                                        max_ratio = 20) {
  n <- length(d)
  groups <- sort(unique(g[!is.na(g)]))
  if (length(groups) < 2L) {
    # Degenerate: only one group. Return uniform weights.
    return(rep(1, n))
  }
  kde_pool <- weighted_kde(d, survey_wgts)
  f_pool_all <- kde_pool(d)

  w_target <- numeric(n)
  eps <- 1e-10
  for (gg in groups) {
    idx_g <- which(g == gg)
    if (length(idx_g) < 2L) {
      w_target[idx_g] <- 1
      next
    }
    wgts_g <- if (!is.null(survey_wgts)) survey_wgts[idx_g] else NULL
    kde_g <- weighted_kde(d[idx_g], wgts_g)
    f_own_g <- pmax(kde_g(d[idx_g]), eps)
    w_target[idx_g] <- f_pool_all[idx_g] / f_own_g
  }

  cap <- stats::quantile(w_target[w_target > 0], probs = 0.98, na.rm = TRUE)
  cap <- max(cap, 1)
  w_target <- pmin(w_target, cap * max_ratio)
  w_target[w_target <= 0] <- 1
  w_target <- w_target / sum(w_target) * n

  w_target
}

# ============================================================
# FIXED EFFECTS (DUMMY) HELPER
# ============================================================

#' Build fixed-effect dummy matrix from a set of variables.
#'
#' Each variable is coerced to a factor, expanded via model.matrix, and the
#' first level dropped to avoid collinearity with the intercept. Columns are
#' named FE_{var}_{level}. Columns with zero variance in the analysis sample
#' are dropped (can happen after tercile/common-support trimming).
#' @keywords internal
#' @noRd
build_fe_dummies <- function(data, fe_vars) {
  if (is.null(fe_vars) || length(fe_vars) == 0) return(NULL)
  miss <- setdiff(fe_vars, names(data))
  if (length(miss) > 0)
    stop(sprintf("FE variable(s) not found in data: %s",
                 paste(miss, collapse = ", ")))
  pieces <- list()
  for (v in fe_vars) {
    vals <- data[[v]]
    if (any(is.na(vals)))
      stop(sprintf("FE variable '%s' contains NA values.", v))
    f <- factor(vals)
    lv <- levels(f)
    if (length(lv) < 2) {
      warning(sprintf("FE variable '%s' has only one level; skipping.", v))
      next
    }
    mm <- stats::model.matrix(~ f)[, -1, drop = FALSE]  # drops first level
    colnames(mm) <- paste0("FE_", v, "_", lv[-1])
    pieces[[v]] <- mm
  }
  if (length(pieces) == 0) return(NULL)
  do.call(cbind, pieces)
}

