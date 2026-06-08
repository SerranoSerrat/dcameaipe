# ============================================================
# dcame_aipe.R
#
# Public-facing function dcame_aipe() (basis path + GAM dispatch),
# the shared results renderer .render_results(), S3 print/summary
# methods, and the export-table helper dcame_aipe_table().
#
# Companion to: "The Pitfalls of Assuming Linear Treatment Effects
#                in Interaction Tests" (Serrano-Serrat).
#
# Model: Y = mu + alpha'f_0(D)*(1-X) + alpha'f_1(D)*X + eta*X
#            + delta'Z + e
#
# This code has been enhanced and optimized using Claude Opus 4.7 and 4.8.
# ============================================================

# Silence R CMD check NOTEs for ggplot2 NSE column references.
utils::globalVariables(c(
  "D", "Yhat", "lo", "up", "Group", "se", "IPE",
  "Specification", "Estimate", "y_start", "y_end",
  "y_label", "label", "D_X"
))

#' Estimate DCAME and AIPE for interactive treatment effects
#'
#' Implements the empirical approach from Serrano-Serrat,
#' *The Pitfalls of Assuming Linear Treatment Effects in Interaction Tests*.
#' Fits the K-group joint model
#' \deqn{Y_i = \sum_{k=0}^{K-1} \big\{ \mu_k + f_k(D_i) + \alpha_k X_i \big\} G_{k,i}
#'              + \delta' Z_i + \mathrm{FE} + \varepsilon_i,}
#' where \eqn{G_{k,i} = \mathbf{1}\{g(X_i) = k\}} is an indicator for the
#' \eqn{k}-th moderator group and \eqn{f_k(\cdot)} is a group-specific
#' polynomial or B-spline (basis path) / penalised smooth (GAM path) of \eqn{D}.
#' For binary \eqn{X}, \eqn{K = 2} and the \eqn{\alpha_k X G_k} term is
#' collinear with \eqn{G_k}, so it is omitted. For continuous \eqn{X},
#' \eqn{K = n\_parts} groups defined by quantile binning, and middle groups
#' are kept (they contribute to fitting nuisance parameters but the headline
#' D-CAME/AIPE contrasts use only \eqn{G_0} vs. \eqn{G_{K-1}}).
#'
#' Returns the Average Interactive Partial Effect (AIPE) and/or the
#' Difference in Conditional Average Marginal Effects (D-CAME), with
#' standard errors, confidence intervals, and diagnostic plots.
#'
#' @param data A `data.frame` containing all variables referenced by name.
#' @param Y Character. Name of the outcome variable.
#' @param D Character. Name of the (continuous) treatment variable.
#' @param X Character. Name of the moderator. Can be binary (0/1), numeric
#'   continuous (auto-binned via `n_parts`), or a two-level
#'   character/factor (auto-coded to 0/1).
#' @param Z Character vector of control variable names. Default `NULL`.
#' @param model Either `"basis"` (default, polynomial/B-spline) or `"GAM"`
#'   (penalised splines via `mgcv`).
#' @param estimand `"AIPE"` or `"DCAME"`. Must be specified.
#' @param vce `"robust"` (default, HC1) or `"cluster"` (requires `cluster`).
#'   Note: for `model = "GAM"` with `inference = "regular"`, the analytic
#'   mgcv posterior is used regardless of `vce`; cluster SEs require
#'   `inference = "bootstrap"`.
#' @param cluster Character. Name of the clustering variable (used when
#'   `vce = "cluster"` for the basis path, and for cluster bootstrap in both
#'   paths).
#' @param wgt Character. Optional name of a survey-weights variable. Honoured
#'   throughout the basis path and the GAM path (fit, prediction Z-bar,
#'   FE modal-level selection, AIPE/CAME averaging, bootstrap resampling).
#' @param Z_interact Logical. If `TRUE`, controls in `Z` are interacted with
#'   the K-group factor (i.e., \eqn{Z \cdot G_k}) in both the basis and GAM
#'   paths. For binary X this is equivalent to \eqn{Z \cdot X}.
#' @param Z_FE Character vector of fixed-effect variable names (additive dummies).
#' @param Z_FE_absorb Logical. If `TRUE`, FE in `Z_FE` are absorbed via
#'   `fixest::feols` instead of expanded into dummies. Use this when one or
#'   more FE variables has high cardinality (school, person, fine grid) and
#'   would blow up memory. Requires the `fixest` package. **Only supported in
#'   the basis path**; ignored (with a warning) under `model = "GAM"`.
#'   Under analytic inference (`"regular"` / `"crossfit"`), predicted-value
#'   CIs reflect within-FE-group uncertainty only; use `inference = "bootstrap"`
#'   for PV CIs comparable to `Z_FE_absorb = FALSE`. Default `FALSE`.
#' @param user_spec_0,user_spec_1,user_spec_2,user_spec_3,user_spec_4
#'   Optional user-specified functional forms, one per group. Numeric =>
#'   polynomial degree; character (e.g. `"3,3"`) => B-spline with knots,degree.
#'   For binary X, only `user_spec_0`/`user_spec_1` are used. For continuous X
#'   with `n_parts = K`, supply the first K (`user_spec_0..user_spec_{K-1}`).
#'   Either supply ALL K or leave ALL NULL (no partial pre-specification).
#' @param compare `FALSE` (default), `"LINEAR"`, `"DCAME"`, `"AIPE"`, or `"ALL"`.
#'   Adds the requested side-by-side comparison(s). The `"LINEAR"` bar is the
#'   joint linear-per-group D-CAME contrast (slope at group K-1 minus slope at
#'   group 0), fitted as `Y ~ D * factor(g) + X * factor(g) + Z + FE` on the
#'   same trimmed sample. It is exactly the basis estimate when
#'   `user_spec_0 = ... = user_spec_{K-1} = 1`.
#' @param selection Model selection criterion for the basis path:
#'   `"CV"` (default, 10-fold per-group CV), `"BIC"` (per-group), or
#'   `"CV_targeted"` (10-fold CV with AIPE-targeted density-ratio weights
#'   `f_pool(D)/f_own_k(D)`; only meaningful when `estimand = "AIPE"` —
#'   falls back to plain CV with a warning otherwise).
#' @param inference `"crossfit"` (default), `"bootstrap"`, or `"regular"`.
#'   GAM path uses analytic (mgcv) or bootstrap inference.
#' @param B_boot Number of bootstrap replicates (default 500).
#' @param n_grid Number of evaluation points on the D-grid (default 1000).
#' @param level_pv Confidence level for predicted-value bands (default 0.84).
#' @param level_est Confidence level for estimands (default 0.95).
#' @param n_parts Number of bins for continuous X (2, 3, 4, or 5; default 3).
#' @param smooth GAM smoothness selection: `"REML"` (default) or `"GCV"`.
#' @param k_s,k_t GAM basis dimensions for `s()` and tensor `te()/ti()` smooths.
#' @param inter_gam GAM continuous-X interaction form: `"ti"` (default) or `"te"`.
#' @param max_poly Highest polynomial degree to search (basis path).
#' @param max_knots Maximum number of interior knots to search (basis path).
#' @param joint_select Logical. Joint vs. independent CV selection across X groups.
#' @param hist_pv Rug/percentile annotation on the predicted-values plot:
#'   `"rug"` (default), `"perc"`, or `"rug_1000"`.
#' @param plot_came_all Logical. If `TRUE` and X is continuous, also produce
#'   `plot_came_all` showing CAME at every level (`n_parts` groups), with CIs.
#'   The K per-group CAMEs come directly from the joint K-group model
#'   (basis path) or from evaluating the fitted GAM at K quantile midpoints
#'   (GAM path). Default `FALSE`. Ignored when X is binary.
#' @param verbose Logical. Print progress messages.
#'
#' @return An object of class `dcame_aipe`. Key fields:
#'   \describe{
#'     \item{`aipe`, `se_aipe`, `ci_aipe`}{Average Interactive Partial Effect.}
#'     \item{`dcame`, `se_dcame`, `ci_dcame`}{Difference in CAMEs.}
#'     \item{`came0`, `se_came0`, `came1`, `se_came1`}{CAME at low / high
#'       moderator group endpoints.}
#'     \item{`linear_dx`, `se_linear_dx`}{Linear D-CAME contrast from the joint
#'       linear-per-group fit (when `compare` includes `"LINEAR"`/`"ALL"`).}
#'     \item{`d_grid`, `yhat_X0`, `yhat_X1`, `se_yhat_X0`, `se_yhat_X1`}{
#'       Predicted-value curves at low/high moderator levels with bands.}
#'     \item{`slope_grid_X0`, `slope_grid_X1`, `ipe_grid`, `se_ipe_grid`}{
#'       Pointwise slope and IPE curves along D.}
#'     \item{`came_all`}{Data frame of per-group CAMEs when `plot_came_all = TRUE`.}
#'     \item{`common_support`, `n_trimmed`, `n_used`}{Common-support trim info.}
#'     \item{`K`, `tercile_cuts`, `Z_FE_absorb`, `compare`}{Reproducibility metadata.}
#'     \item{`plot`, `plot_predicted`, `plot_ipe`, `plot_coef`, `plot_came_all`}{
#'       Ready-to-print `ggplot` / `gtable` objects.}
#'     \item{`fit`, `vcov`}{Underlying fit object (when stored) and vcov.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage with a simulated data frame
#' res <- dcame_aipe(data = mydata, Y = "Y", D = "D", X = "X",
#'                   estimand = "AIPE", compare = "ALL")
#' summary(res)
#' res$plot
#' }
#'
#' @import stats
#' @export
dcame_aipe <- function(data, Y, D, X, Z = NULL,
                       model = "basis",
                       estimand,
                       vce = "robust", cluster = NULL,
                       wgt = NULL,
                       Z_interact = FALSE,
                       Z_FE = NULL,
                       Z_FE_absorb = FALSE,
                       user_spec_0 = NULL, user_spec_1 = NULL,
                       user_spec_2 = NULL, user_spec_3 = NULL,
                       user_spec_4 = NULL,
                       compare = FALSE,
                       selection = "CV",
                       inference = "crossfit",
                       B_boot = 500,
                       n_grid = 1000,
                       level_pv = 0.84,
                       level_est = 0.95,
                       n_parts = 3,
                       smooth = "REML",
                       k_s = NULL,
                       k_t = NULL,
                       inter_gam = "ti",
                       max_poly = 5,
                       max_knots = 4,
                       joint_select = FALSE,
                       hist_pv = "rug",
                       plot_came_all = FALSE,
                       verbose = TRUE) {

  # ══════════════════════════════════════════════════════════════════════════
  # ARGUMENT VALIDATION
  # ══════════════════════════════════════════════════════════════════════════

  if (missing(estimand))
    stop("`estimand` must be specified explicitly: \"AIPE\" or \"DCAME\".")

  model       <- match.arg(model, c("basis", "GAM"))
  smooth      <- match.arg(smooth, c("REML", "GCV"))
  selection   <- match.arg(selection,  c("CV", "BIC", "CV_targeted"))
  inference   <- match.arg(inference,  c("crossfit", "bootstrap", "regular"))
  hist_pv     <- match.arg(hist_pv, c("rug", "perc", "rug_1000"))
  estimand    <- match.arg(estimand, c("AIPE", "DCAME"))
  inter_gam   <- match.arg(inter_gam, c("ti", "te"))

  if (!is.null(k_s)) {
    if (!is.numeric(k_s) || length(k_s) != 1 || k_s < 3 || k_s != as.integer(k_s))
      stop("`k_s` must be NULL or a single integer >= 3.")
    k_s <- as.integer(k_s)
  }
  if (!is.null(k_t)) {
    if (!is.numeric(k_t) || !(length(k_t) %in% c(1, 2)) ||
        any(k_t < 3) || any(k_t != as.integer(k_t)))
      stop("`k_t` must be NULL or an integer vector of length 1 or 2, each >= 3.")
    k_t <- as.integer(k_t)
    if (length(k_t) == 1) k_t <- rep(k_t, 2)
  }

  # ── Hard-coded internals ──────────────────────────────────────────────
  spline_type <- "bspline"
  FE          <- NULL
  FE_vbles    <- NULL

  # Collect user_spec_0..user_spec_4 into a list (length 5).
  user_spec_all <- list(user_spec_0, user_spec_1, user_spec_2, user_spec_3, user_spec_4)

  # Helper: convert a user_spec value to a spec list, or NULL.
  .parse_user_spec_value <- function(v) {
    if (is.null(v)) return(NULL)
    if (is.numeric(v)) return(list(type = "poly", degree = as.integer(v), knots = 0))
    if (is.character(v)) return(parse_spline_spec(v, type = spline_type))
    stop("user_spec_* must be NULL, a single integer (polynomial degree), or a 'knots,degree' string.")
  }

  if (Z_FE_absorb && !requireNamespace("fixest", quietly = TRUE)) {
    stop("Z_FE_absorb = TRUE requires the 'fixest' package. ",
         "Install it with install.packages(\"fixest\").")
  }
  if (Z_FE_absorb && is.null(Z_FE)) {
    warning("Z_FE_absorb = TRUE but Z_FE is NULL: there are no fixed effects to absorb. ",
            "Setting Z_FE_absorb = FALSE.", call. = FALSE)
    Z_FE_absorb <- FALSE
  }
  stopifnot(is.numeric(level_pv), length(level_pv) == 1,
            level_pv > 0, level_pv < 1)
  stopifnot(is.numeric(level_est), length(level_est) == 1,
            level_est > 0, level_est < 1)

  stopifnot(is.numeric(n_parts), length(n_parts) == 1,
            n_parts %in% c(2, 3, 4, 5))

  if (!is.null(Z_FE)) {
    if (!is.character(Z_FE))
      stop("Z_FE must be a character vector of variable names.")
    miss_fe <- setdiff(Z_FE, names(data))
    if (length(miss_fe) > 0)
      stop(sprintf("Z_FE variable(s) not found in data: %s",
                   paste(miss_fe, collapse = ", ")))
    if (length(Z_FE) == 0) Z_FE <- NULL
  }

  if (model == "GAM") {
    if (inference == "crossfit") {
      warning("model = 'GAM' does not support cross-fitting. ",
              "Inference has been overridden to 'regular' (analytic mgcv posterior). ",
              "Set inference = 'bootstrap' for a case/cluster bootstrap.",
              call. = FALSE)
      inference <- "regular"
    }
    if (isTRUE(Z_FE_absorb)) {
      warning("model = 'GAM' does not support Z_FE_absorb (no FE absorption in mgcv). ",
              "Z_FE variables will be included as additive factor terms in the GAM ",
              "(every level becomes a coefficient). For high-cardinality FE, ",
              "consider model = 'basis' with Z_FE_absorb = TRUE.",
              call. = FALSE)
      Z_FE_absorb <- FALSE
    }
  }

  # Warn about predicted-value CI semantics when absorbing FE under analytic inference.
  # Slopes / CAME / AIPE / D-CAME CIs are unaffected (FWL), but predicted-value
  # bands omit the variance of the average FE contribution. Bootstrap captures it.
  # Placed AFTER the GAM check so it doesn't fire spuriously for GAM (which has
  # already coerced Z_FE_absorb to FALSE above).
  if (isTRUE(Z_FE_absorb) && inference != "bootstrap") {
    warning(
      "Z_FE_absorb = TRUE: predicted-value CIs from analytic inference (",
      inference, ") reflect within-FE-group uncertainty only and will be ",
      "narrower than under Z_FE_absorb = FALSE. CAME / AIPE / D-CAME CIs ",
      "are unaffected (FWL). For PV CIs comparable to Z_FE_absorb = FALSE, ",
      "use inference = \"bootstrap\".",
      call. = FALSE
    )
  }

  # Warn if B_boot was explicitly set but won't be used.
  if (!missing(B_boot) && inference != "bootstrap") {
    warning(sprintf(
      "B_boot = %d was supplied but inference = '%s' does not use bootstrap; ",
      B_boot, inference),
      "B_boot will be ignored. Set inference = 'bootstrap' to use it.",
      call. = FALSE)
  }

  if (isTRUE(compare)) {
    compare <- if (estimand == "AIPE") "DCAME" else "AIPE"
  }
  if (!isFALSE(compare)) {
    if (is.character(compare)) {
      valid_compare <- c("LINEAR", "DCAME", "AIPE", "ALL")
      compare <- toupper(compare)
      compare <- match.arg(compare, valid_compare)
    } else {
      stop("compare must be FALSE, 'LINEAR', 'DCAME', 'AIPE', or 'ALL'")
    }
  }

  # CV_targeted only makes sense for AIPE (the density-ratio weights target the
  # pooled-D distribution that AIPE integrates over). For DCAME, fall back to
  # plain CV and inform the user.
  if (selection == "CV_targeted" && estimand != "AIPE") {
    warning("selection = 'CV_targeted' is only meaningful when estimand = 'AIPE' ",
            "(the targeted weights align CV with the AIPE estimand). ",
            "Falling back to plain 10-fold CV for estimand = '", estimand, "'.",
            call. = FALSE)
    selection <- "CV"
  }
  targeted_cv <- selection == "CV_targeted"

  y <- data[[Y]]; d <- data[[D]]; x_raw <- data[[X]]

  .strip_labels <- function(v) {
    if (inherits(v, c("haven_labelled", "labelled"))) as.numeric(v) else v
  }
  y     <- .strip_labels(y)
  d     <- .strip_labels(d)
  x_raw <- .strip_labels(x_raw)

  x_original_labels <- NULL
  if (is.character(x_raw) || is.factor(x_raw)) {
    x_lvls <- if (is.factor(x_raw)) levels(droplevels(as.factor(x_raw)))
    else sort(unique(x_raw[!is.na(x_raw)]))
    if (length(x_lvls) == 0)
      stop(sprintf("Moderator '%s' contains only NA values.", X))
    if (length(x_lvls) == 1)
      stop(sprintf("Moderator '%s' has only one unique level ('%s'). A moderator requires at least two levels.",
                   X, x_lvls[1]))
    if (length(x_lvls) != 2)
      stop(sprintf(paste0("Character/factor moderator '%s' must have exactly 2 levels, ",
                          "found %d: %s. Consider recoding into a binary variable or ",
                          "using a numeric moderator with n_parts."),
                   X, length(x_lvls), paste(x_lvls, collapse = ", ")))
    x_original_labels <- x_lvls
    x_raw <- as.integer(x_raw == x_lvls[2])
    if (verbose) cat(sprintf("Moderator X auto-coded: '%s' = 0, '%s' = 1\n",
                             x_original_labels[1], x_original_labels[2]))
  }

  if (!is.numeric(y))
    stop(sprintf("Outcome '%s' must be numeric. Found class: %s. %s",
                 Y, paste(class(y), collapse = "/"),
                 if (is.character(y)) "If this is a categorical outcome, recode to numeric first."
                 else ""))
  if (!is.numeric(d))
    stop(sprintf("Treatment '%s' must be numeric. Found class: %s. %s",
                 D, paste(class(d), collapse = "/"),
                 if (is.character(d)) "Character treatments are not supported; recode to numeric."
                 else ""))
  if (!is.numeric(x_raw))
    stop(sprintf("Moderator '%s' must be numeric or a two-level character/factor. Found class: %s.",
                 X, paste(class(x_raw), collapse = "/")))
  if (vce == "cluster" && is.null(cluster))
    stop("cluster variable name required when vce = 'cluster'.")
  if (!is.null(cluster) && vce != "cluster" && inference != "bootstrap") {
    warning("`cluster` was supplied but vce = '", vce,
            "' and inference = '", inference,
            "'. The cluster variable will be ignored for standard errors. ",
            "Set vce = 'cluster' (analytic) or inference = 'bootstrap' (cluster bootstrap) to use it.",
            call. = FALSE)
  }
  if (!is.null(cluster) && vce == "cluster") {
    n_clusters_check <- length(unique(data[[cluster]][!is.na(data[[cluster]])]))
    if (n_clusters_check < 30L) {
      warning(sprintf(
        "Only %d unique cluster(s) detected (variable '%s'). ",
        n_clusters_check, cluster),
        "Liang-Zeger / CR1 cluster-robust SEs can be badly biased downward when ",
        "the number of clusters is small (rule of thumb: G >= 30). ",
        "Consider inference = 'bootstrap' (cluster bootstrap), which has better ",
        "finite-sample behavior.",
        call. = FALSE)
    }
  }

  n_total <- length(y)
  n_na_y <- sum(is.na(y)); n_na_d <- sum(is.na(d)); n_na_x <- sum(is.na(x_raw))
  est_complete <- !is.na(y) & !is.na(d) & !is.na(x_raw)
  if (!is.null(Z)) {
    z_check <- data[, Z, drop = FALSE]
    z_na <- !complete.cases(z_check)
    est_complete <- est_complete & !z_na
  }
  if (!is.null(wgt)) {
    w_tmp <- .strip_labels(data[[wgt]])
    est_complete <- est_complete & !is.na(w_tmp)
  }
  if (!is.null(cluster)) est_complete <- est_complete & !is.na(data[[cluster]])
  if (!is.null(Z_FE)) est_complete <- est_complete & complete.cases(data[, Z_FE, drop = FALSE])
  if (!is.null(FE) && !is.null(FE_vbles))
    est_complete <- est_complete & complete.cases(data[, FE_vbles, drop = FALSE])

  n_na_total <- sum(!est_complete)
  if (n_na_total > 0 && verbose) {
    parts <- character(0)
    if (n_na_y > 0) parts <- c(parts, sprintf("%d in Y", n_na_y))
    if (n_na_d > 0) parts <- c(parts, sprintf("%d in D", n_na_d))
    if (n_na_x > 0) parts <- c(parts, sprintf("%d in X", n_na_x))
    n_other <- n_na_total - sum(!(!is.na(y) & !is.na(d) & !is.na(x_raw)))
    if (n_other > 0) parts <- c(parts, sprintf("%d in Z/weights/cluster/FE", n_other))
    cat(sprintf("Note: %d/%d observations excluded from estimation due to missing values (%s).\n",
                n_na_total, n_total, paste(parts, collapse = ", ")))
  }
  if (sum(est_complete) < 10)
    stop(sprintf("Only %d complete observations available for estimation. Need at least 10.",
                 sum(est_complete)))

  d_full_for_aipe <- d[!is.na(d)]
  survey_wgts_full_for_aipe <- if (!is.null(wgt)) {
    w_tmp <- .strip_labels(data[[wgt]])
    w_tmp[!is.na(d) & !is.na(w_tmp)]
  } else NULL
  if (!is.null(survey_wgts_full_for_aipe) &&
      length(survey_wgts_full_for_aipe) != length(d_full_for_aipe)) {
    survey_wgts_full_for_aipe <- NULL
  }

  y     <- y[est_complete]
  d     <- d[est_complete]
  x_raw <- x_raw[est_complete]
  data  <- data[est_complete, , drop = FALSE]

  survey_wgts <- NULL
  if (!is.null(wgt)) {
    survey_wgts <- data[[wgt]]
    survey_wgts <- .strip_labels(survey_wgts)
    stopifnot(is.numeric(survey_wgts), all(survey_wgts > 0))
    if (verbose) cat(sprintf("Survey weights from '%s' (range: [%.2f, %.2f])\n",
                             wgt, min(survey_wgts), max(survey_wgts)))
  }

  # ══════════════════════════════════════════════════════════════════════════
  # DISPATCH: GAM PATH or BASIS PATH
  # ══════════════════════════════════════════════════════════════════════════

  if (model == "GAM") {
    if (!requireNamespace("mgcv", quietly = TRUE))
      stop("Package 'mgcv' is required for model = 'GAM'. Install it with install.packages('mgcv').")

    return(.dcame_aipe_gam(y = y, d = d, x_raw = x_raw, data = data,
                           X_name = X, Z = Z,
                           cluster = cluster, vce = vce,
                           survey_wgts = survey_wgts,
                           wgt = wgt,
                           estimand = estimand, compare = compare,
                           B_boot = B_boot,
                           n_grid = n_grid,
                           level_pv = level_pv, level_est = level_est,
                           n_parts = n_parts,
                           hist_pv = hist_pv, verbose = verbose,
                           Z_FE = Z_FE,
                           FE = FE, FE_vbles = FE_vbles,
                           inference = inference,
                           smooth = smooth,
                           k_s = k_s, k_t = k_t,
                           inter_gam = inter_gam,
                           Z_interact = Z_interact,
                           plot_came_all = plot_came_all,
                           x_original_labels = x_original_labels))
  }
  # ══════════════════════════════════════════════════════════════════════════
  # BASIS PATH (K-GROUP JOINT MODEL)
  # ══════════════════════════════════════════════════════════════════════════
  #
  # Model: Y = sum_{k=0..K-1} alpha_k' f_k(D)*I(g=k)
  #            + sum_{k=1..K-1} eta_k * I(g=k)
  #            + delta'Z + FE + epsilon
  #
  # Group 0 = lowest, group K-1 = highest. Headline estimands compare those
  # two endpoints; middle groups contribute to fitting nuisance parameters
  # (Z, FE) but not directly to AIPE/DCAME. Each group has its own basis
  # f_k() — different functional forms across groups are the rule, not the
  # exception. For binary X (auto-detected), K = 2 and g = x_raw.

  x_is_binary <- all(x_raw %in% c(0, 1))
  tercile_cuts <- NULL
  x_continuous_raw <- NULL

  if (x_is_binary) {
    K <- 2L
    g <- as.integer(x_raw)
    if (verbose) cat("Moderator X detected as binary.\n")
    if (!missing(n_parts) && n_parts != 2L) {
      warning(sprintf(
        "n_parts = %d was supplied but moderator '%s' is binary (0/1). ",
        n_parts, X),
        "n_parts is ignored; X is used as the group indicator directly (K = 2).",
        call. = FALSE)
    }
  } else {
    K <- as.integer(n_parts)
    part_label <- switch(as.character(K),
                         "2" = "halves (median split)",
                         "3" = "terciles",
                         "4" = "quartiles",
                         "5" = "quintiles")
    if (verbose) cat(sprintf("Moderator X detected as continuous. Binning into %s (n_parts=%d).\n",
                             part_label, K))

    cut_probs <- seq(0, 1, length.out = K + 1L)[-c(1L, K + 1L)]
    tercile_cuts <- as.numeric(stats::quantile(x_raw, probs = cut_probs, na.rm = TRUE))

    g <- rep(NA_integer_, length(x_raw))
    g[x_raw <= tercile_cuts[1]] <- 0L
    if (K > 2L) {
      for (kk in 2:(K - 1L)) {
        g[x_raw > tercile_cuts[kk - 1L] & x_raw <= tercile_cuts[kk]] <- as.integer(kk - 1L)
      }
    }
    g[x_raw > tercile_cuts[length(tercile_cuts)]] <- as.integer(K - 1L)

    x_continuous_raw <- x_raw
    if (verbose) {
      cat(sprintf("  Cut-points: %s\n",
                  paste(sprintf("%.3f", tercile_cuts), collapse = ", ")))
      sizes <- vapply(0:(K - 1L), function(kk) sum(g == kk, na.rm = TRUE), integer(1))
      cat(sprintf("  Group sizes (low..high, no obs dropped): %s\n",
                  paste(sizes, collapse = ", ")))
    }
  }

  # ── User-spec validation: all-or-nothing across K groups ─────────────
  if (length(user_spec_all) > K) {
    for (extra_i in (K + 1L):length(user_spec_all)) {
      if (!is.null(user_spec_all[[extra_i]]))
        stop(sprintf("user_spec_%d is non-NULL but only %d groups are needed (n_parts=%d or binary X). Set it to NULL.",
                     extra_i - 1L, K, n_parts))
    }
  }
  active_idx <- which(!vapply(user_spec_all[1:K], is.null, logical(1)))
  if (length(active_idx) > 0L && length(active_idx) < K) {
    missing_idx <- setdiff(1:K, active_idx) - 1L
    stop(sprintf("Partial pre-specification not allowed. Either set ALL user_spec_0..user_spec_%d, or leave all NULL for CV. Missing: %s.",
                 K - 1L,
                 paste(sprintf("user_spec_%d", missing_idx), collapse = ", ")))
  }

  user_specs <- lapply(user_spec_all[1:K], .parse_user_spec_value)
  all_user_specified <- length(active_idx) == K

  if (all_user_specified && inference == "crossfit") {
    warning("All ", K, " user_spec_* arguments are supplied, so no model selection is needed; ",
            "cross-fitting (which exists to separate selection from inference) is unnecessary. ",
            "Inference has been overridden to 'regular'.",
            call. = FALSE)
    inference <- "regular"
  }

  if (joint_select && K > 3L)
    stop(sprintf("joint_select = TRUE is supported only for K <= 3 (you have K=%d). Use joint_select = FALSE for per-group independent CV.", K))

  # ── Covariate matrix Z (and continuous-X linear control) ─────────────
  d_all_for_aipe <- d_full_for_aipe
  survey_wgts_all_for_aipe <- survey_wgts_full_for_aipe

  z_mat <- if (!is.null(Z)) {
    m <- as.matrix(data[, Z, drop = FALSE])
    storage.mode(m) <- "double"
    m[, apply(m, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  } else NULL

  if (!x_is_binary) {
    x_ctrl <- matrix(x_continuous_raw, ncol = 1)
    colnames(x_ctrl) <- "X_continuous"
    z_mat <- if (!is.null(z_mat)) cbind(x_ctrl, z_mat) else x_ctrl
  }

  # ── Fixed effects: either expand to dummies, or hold raw for absorption
  fe_mat <- NULL
  fe_absorb_factors <- NULL
  if (!is.null(Z_FE)) {
    if (Z_FE_absorb) {
      fe_absorb_factors <- lapply(Z_FE, function(v) factor(data[[v]]))
      names(fe_absorb_factors) <- Z_FE
      if (verbose) cat(sprintf("Z_FE: will absorb %s via fixest::feols.\n",
                               paste(Z_FE, collapse = ", ")))
    } else {
      fe_mat <- build_fe_dummies(data, Z_FE)
      if (!is.null(fe_mat)) {
        sds_fe <- apply(fe_mat, 2, sd, na.rm = TRUE)
        fe_mat <- fe_mat[, sds_fe > 0, drop = FALSE]
        if (ncol(fe_mat) == 0) fe_mat <- NULL
        if (!is.null(fe_mat) && verbose)
          cat(sprintf("Z_FE: added %d fixed-effect dummies from %s\n",
                      ncol(fe_mat), paste(Z_FE, collapse = ", ")))
      }
    }
  }

  cluster_var <- if (!is.null(cluster)) data[[cluster]] else NULL
  n <- length(y)
  zc_pv  <- qnorm(1 - (1 - level_pv) / 2)
  zc_est <- qnorm(1 - (1 - level_est) / 2)

  # ── Common support: intersection of group 0 and group K-1 D-ranges ───
  d_g0 <- d[g == 0L]
  d_gK <- d[g == (K - 1L)]
  cs_lo <- max(min(d_g0, na.rm = TRUE), min(d_gK, na.rm = TRUE))
  cs_hi <- min(max(d_g0, na.rm = TRUE), max(d_gK, na.rm = TRUE))
  if (is.na(cs_lo) || is.na(cs_hi) || cs_lo >= cs_hi)
    stop(sprintf("No common support between group 0 (D: [%.3g, %.3g]) and group K-1 (D: [%.3g, %.3g]).",
                 min(d_g0, na.rm = TRUE), max(d_g0, na.rm = TRUE),
                 min(d_gK, na.rm = TRUE), max(d_gK, na.rm = TRUE)))

  keep <- (d >= cs_lo) & (d <= cs_hi) & !is.na(g)
  nt <- sum(!keep)
  if (nt > 0) {
    n_pre <- length(y)
    share_total <- 100 * nt / n_pre
    # Per-group share of the trim
    trim_msg_parts <- character(0)
    for (kk in 0:(K - 1L)) {
      n_g <- sum(g == kk, na.rm = TRUE)
      n_trim_g <- sum(!keep & !is.na(g) & g == kk)
      lbl <- if (kk == 0L) "low-X" else if (kk == K - 1L) "high-X" else sprintf("mid-X(%d)", kk)
      trim_msg_parts <- c(trim_msg_parts,
                          sprintf("%s: %.1f%%", lbl,
                                  if (n_g > 0) 100 * n_trim_g / n_g else 0))
    }
    warning(sprintf(
      "Common-support trimming: %d / %d observations trimmed (%.1f%% of sample). ",
      nt, n_pre, share_total),
      "Breakdown by moderator group (share of within-group obs trimmed): ",
      paste(trim_msg_parts, collapse = "; "), ". ",
      "Trimming keeps only D in [", sprintf("%.3g", cs_lo), ", ", sprintf("%.3g", cs_hi),
      "], the overlap between group 0 and group K-1. ",
      "Large trims (>5%) can shift estimands materially. ",
      "Note: the Linear comparison reported below uses the same trimmed sample for consistency.",
      call. = FALSE)
    y <- y[keep]; d <- d[keep]; g <- g[keep]
    if (!is.null(z_mat)) z_mat <- z_mat[keep, , drop = FALSE]
    if (!is.null(fe_mat)) fe_mat <- fe_mat[keep, , drop = FALSE]
    if (!is.null(fe_absorb_factors))
      fe_absorb_factors <- lapply(fe_absorb_factors, function(f) factor(f[keep]))
    if (!is.null(cluster_var)) cluster_var <- cluster_var[keep]
    if (!is.null(survey_wgts)) survey_wgts <- survey_wgts[keep]
    if (!x_is_binary) x_continuous_raw <- x_continuous_raw[keep]
    n <- length(y)
  }

  d_all_for_aipe <- pmin(pmax(d_all_for_aipe, cs_lo), cs_hi)

  if (!is.null(fe_mat)) {
    sds_fe <- apply(fe_mat, 2, sd, na.rm = TRUE)
    if (any(sds_fe == 0)) {
      fe_mat <- fe_mat[, sds_fe > 0, drop = FALSE]
      if (ncol(fe_mat) == 0) fe_mat <- NULL
    }
  }

  # When Z_FE_absorb = TRUE, FE are absorbed via fixest::feols in the
  # fitting steps below. We only need to track df lost (for reporting).
  absorbed_df_lost <- if (Z_FE_absorb && !is.null(fe_absorb_factors))
    sum(vapply(fe_absorb_factors, function(f) nlevels(f) - 1L, integer(1)))
  else 0L

  d_obs_by_group <- lapply(0:(K - 1L), function(kk) d[g == kk])
  d_obs_X0 <- d_obs_by_group[[1L]]
  d_obs_X1 <- d_obs_by_group[[K]]

  mod_label <- if (x_is_binary) "binary" else sprintf("continuous (n_parts=%d)", K)
  cat(sprintf("=== DCAME/AIPE Estimation (model = basis, estimand = %s, inference = %s, moderator = %s) ===\n",
              estimand, inference, mod_label))
  cat(sprintf("N = %d | Common support [g=0 vs g=K-1]: [%.3f, %.3f]", n, cs_lo, cs_hi))
  if (nt > 0) cat(sprintf(" | %d trimmed", nt))
  cat(sprintf(" | Z_interact: %s", Z_interact))
  if (Z_FE_absorb) cat(sprintf(" | FE absorbed: %s", paste(Z_FE, collapse = ", ")))
  else if (!is.null(fe_mat)) cat(sprintf(" | Z_FE: %d dummies", ncol(fe_mat)))
  if (!is.null(wgt)) cat(" | weighted")
  cat("\n")

  # ══════════════════════════════════════════════════════════════════════════
  # INTERNAL HELPERS
  # ══════════════════════════════════════════════════════════════════════════

  .get_knots_one <- function(d_sub, spec) {
    if (spec$type == "poly") return(list(k = NULL, bk = NULL, ok = TRUE))
    deg <- if (!is.null(spec$degree)) spec$degree else 3
    if (length(d_sub) <= spec$knots + deg + 1)
      return(list(k = NULL, bk = NULL, ok = FALSE))
    strategy <- if (!is.null(spec$knot_strategy)) spec$knot_strategy else "quantile"
    list(k = place_knots(d_sub, spec$knots, strategy),
         bk = compute_boundary_knots(d_sub),
         ok = TRUE)
  }

  # Build K-group joint design matrix.
  # Implements: Y = sum_k { mu_k + f_k(D) + alpha_k * X } * G_k + delta'Z + FE + e
  # The alpha_k * X * G_k block is included whenever x_raw_loc is supplied
  # (continuous X). For binary X, X equals G and the block would be collinear,
  # so x_raw_loc = NULL skips it.
  .build_design_K <- function(d_vec, g_vec, z_local, specs_list, kns_list,
                              fe_local = NULL, x_raw_loc = NULL) {
    nn <- length(d_vec)
    basis_blocks <- vector("list", K)
    for (kk_i in 0:(K - 1L)) {
      bk <- create_basis(d_vec, specs_list[[kk_i + 1L]],
                         kns_list[[kk_i + 1L]]$k, kns_list[[kk_i + 1L]]$bk)
      block <- bk$basis * as.numeric(g_vec == kk_i)
      colnames(block) <- paste0("B_g", kk_i, "_", colnames(bk$basis))
      basis_blocks[[kk_i + 1L]] <- block
    }
    design <- do.call(cbind, basis_blocks)

    if (K > 1L) {
      gd <- matrix(0, nn, K - 1L)
      for (kk_i in 1:(K - 1L)) gd[, kk_i] <- as.numeric(g_vec == kk_i)
      colnames(gd) <- paste0("G_", 1:(K - 1L))
      design <- cbind(design, gd)
    }

    # Continuous X: add alpha_k * X * I(g=k) blocks (group-specific linear X).
    # Drop X_continuous from z_local to avoid double-counting.
    z_eff <- z_local
    if (!is.null(x_raw_loc)) {
      if (!is.null(z_eff)) {
        drop_idx <- which(colnames(z_eff) == "X_continuous")
        if (length(drop_idx) > 0L)
          z_eff <- z_eff[, -drop_idx, drop = FALSE]
        if (ncol(z_eff) == 0L) z_eff <- NULL
      }
      x_blocks <- matrix(0, nn, K)
      for (kk_i in 0:(K - 1L))
        x_blocks[, kk_i + 1L] <- as.numeric(x_raw_loc) * as.numeric(g_vec == kk_i)
      colnames(x_blocks) <- paste0("X_g", 0:(K - 1L))
      design <- cbind(design, x_blocks)
    }

    if (!is.null(z_eff) && ncol(z_eff) > 0) {
      if (Z_interact) {
        for (kk_i in 0:(K - 1L)) {
          zk <- z_eff * as.numeric(g_vec == kk_i)
          colnames(zk) <- paste0("Z_g", kk_i, "_", colnames(z_eff))
          design <- cbind(design, zk)
        }
      } else {
        zt <- z_eff
        colnames(zt) <- paste0("Z_", colnames(z_eff))
        design <- cbind(design, zt)
      }
    }
    if (!is.null(fe_local) && ncol(fe_local) > 0) {
      design <- cbind(design, fe_local)
    }
    design
  }

  .make_candidate_specs <- function() {
    specs <- list()
    for (dd in 1:max_poly)
      specs[[length(specs) + 1L]] <- list(type = "poly", degree = dd, knots = 0)
    for (kk in 1:max_knots)
      specs[[length(specs) + 1L]] <- list(type = "bspline", knots = kk, degree = 3,
                                          knot_strategy = "quantile")
    specs
  }

  # Single fit helper: lm() + sandwich for the default path; fixest::feols
  # with absorbed FE when fe_factors is supplied (Z_FE_absorb = TRUE).
  .fit_lm <- function(y_loc, design_mat, survey_w = NULL, cl_var = NULL,
                      vce_type = "robust", fe_factors = NULL) {
    if (!is.null(fe_factors)) {
      dfx <- data.frame(.Y = y_loc, design_mat, check.names = FALSE)
      for (v in names(fe_factors)) dfx[[v]] <- fe_factors[[v]]
      design_rhs <- paste0("`", colnames(design_mat), "`", collapse = " + ")
      fe_rhs     <- paste(paste0("`", names(fe_factors), "`"), collapse = " + ")
      fml <- stats::as.formula(sprintf(".Y ~ %s | %s", design_rhs, fe_rhs))
      if (vce_type == "cluster" && !is.null(cl_var)) {
        dfx$.CL_VAR <- cl_var
        vcov_arg <- stats::as.formula("~ .CL_VAR")
      } else {
        vcov_arg <- "hetero"
      }
      fit <- fixest::feols(fml, data = dfx,
                           weights = survey_w,
                           vcov    = vcov_arg,
                           warn    = FALSE, notes = FALSE)
      V  <- as.matrix(stats::vcov(fit))
      bh <- stats::coef(fit)
      return(list(fit = fit, V = V, bh = bh))
    }
    dfx <- data.frame(.Y = y_loc, design_mat, check.names = FALSE)
    fit <- if (!is.null(survey_w))
      stats::lm(.Y ~ ., data = dfx, weights = survey_w)
    else
      stats::lm(.Y ~ ., data = dfx)
    V <- if (vce_type == "cluster" && !is.null(cl_var))
      sandwich::vcovCL(fit, cluster = cl_var)
    else sandwich::vcovHC(fit, type = "HC1")
    bh <- stats::coef(fit)
    list(fit = fit, V = V, bh = bh)
  }

  # ══════════════════════════════════════════════════════════════════════════
  # MODEL SELECTION (per-group independent CV by default; joint up to K=3)
  # ══════════════════════════════════════════════════════════════════════════

  .select_model_K <- function(y_loc, d_loc, g_loc, z_loc, user_specs_loc,
                              sel_method, survey_w = NULL, verbose_sel = FALSE,
                              fe_loc = NULL, fe_factors_loc = NULL,
                              x_raw_loc = NULL) {
    if (all(!vapply(user_specs_loc, is.null, logical(1)))) return(user_specs_loc)

    asp <- .make_candidate_specs()
    sp_lists <- lapply(user_specs_loc, function(u) if (is.null(u)) asp else list(u))

    is_bic <- sel_method == "BIC"

    nfolds <- 10L
    .old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(123)
    fid <- sample(rep(1:nfolds, length.out = length(y_loc)))
    if (!is.null(.old_seed))
      assign(".Random.seed", .old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv))
      rm(".Random.seed", envir = .GlobalEnv)

    # Fit + predict helper for CV: uses feols with absorbed FE if fe_train
    # is provided, else falls back to lm. df_tr/df_te already include the
    # basis (+ z) columns. Returns numeric predictions or NULL on failure.
    .fit_predict_sel <- function(df_tr, df_te, w_tr, fe_train, fe_test) {
      if (!is.null(fe_train)) {
        for (v in names(fe_train)) {
          df_tr[[v]] <- fe_train[[v]]
          df_te[[v]] <- fe_test[[v]]
        }
        design_cols <- setdiff(colnames(df_tr), c("Y_out", names(fe_train)))
        design_rhs  <- paste0("`", design_cols, "`", collapse = " + ")
        fe_rhs      <- paste(paste0("`", names(fe_train), "`"), collapse = " + ")
        fml <- stats::as.formula(sprintf("Y_out ~ %s | %s", design_rhs, fe_rhs))
        ft <- tryCatch(fixest::feols(fml, data = df_tr, weights = w_tr,
                                     warn = FALSE, notes = FALSE),
                       error = function(e) NULL)
        if (is.null(ft)) return(NULL)
        return(tryCatch(stats::predict(ft, newdata = df_te),
                        error = function(e) NULL))
      }
      ft <- tryCatch({
        if (!is.null(w_tr))
          stats::lm(Y_out ~ ., data = df_tr, weights = w_tr)
        else stats::lm(Y_out ~ ., data = df_tr)
      }, error = function(e) NULL)
      if (is.null(ft)) return(NULL)
      tryCatch(stats::predict(ft, newdata = df_te), error = function(e) NULL)
    }

    # BIC helper for the BIC selection branch.
    .bic_fit_sel <- function(df_g, w_g, fe_g) {
      if (!is.null(fe_g)) {
        for (v in names(fe_g)) df_g[[v]] <- fe_g[[v]]
        design_cols <- setdiff(colnames(df_g), c("Y_out", names(fe_g)))
        design_rhs  <- paste0("`", design_cols, "`", collapse = " + ")
        fe_rhs      <- paste(paste0("`", names(fe_g), "`"), collapse = " + ")
        fml <- stats::as.formula(sprintf("Y_out ~ %s | %s", design_rhs, fe_rhs))
        ft <- tryCatch(fixest::feols(fml, data = df_g, weights = w_g,
                                     warn = FALSE, notes = FALSE),
                       error = function(e) NULL)
        if (is.null(ft)) return(Inf)
        return(stats::BIC(ft))
      }
      ft <- tryCatch({
        if (!is.null(w_g))
          stats::lm(Y_out ~ ., data = df_g, weights = w_g)
        else stats::lm(Y_out ~ ., data = df_g)
      }, error = function(e) NULL)
      if (is.null(ft)) return(Inf)
      stats::BIC(ft)
    }

    # Helper to subset a list of FE factors to a row-index set.
    .fe_factors_subset <- function(fe_full, idx) {
      if (is.null(fe_full)) return(NULL)
      lapply(fe_full, function(f) factor(as.character(f[idx])))
    }

    # ── AIPE-targeted CV weights (Appendix CV).
    # When sel_method == "CV_targeted", w_target[i] = f_pool(D_i)/f_own_g(D_i)
    # re-weights the per-fold MSE so it integrates over the pooled-D distribution
    # rather than the within-group one. Combined with survey weights when both
    # are present. Computed once on the full sample passed to .select_model_K().
    w_target_full <- NULL
    if (sel_method == "CV_targeted") {
      w_target_full <- tryCatch(
        compute_targeted_cv_weights(d_loc, g_loc, survey_w),
        error = function(e) {
          if (verbose_sel)
            cat("    Warning: targeted CV weights failed; falling back to uniform.\n")
          NULL
        })
      if (!is.null(w_target_full) &&
          (any(!is.finite(w_target_full)) || any(w_target_full < 0))) {
        if (verbose_sel)
          cat("    Warning: targeted CV weights non-finite; falling back to uniform.\n")
        w_target_full <- NULL
      }
      if (!is.null(w_target_full) && verbose_sel)
        cat("    AIPE-targeted CV: density-ratio weights computed.\n")
    }
    # Helper: build the effective per-row CV weight in a fold.
    .cv_fold_weight <- function(rows) {
      sw <- if (!is.null(survey_w)) survey_w[rows] else NULL
      wt <- if (!is.null(w_target_full)) w_target_full[rows] else NULL
      if (!is.null(sw) && !is.null(wt)) return(sw * wt)
      if (!is.null(sw)) return(sw)
      if (!is.null(wt)) return(wt)
      NULL
    }

    # ── Per-group CV: each group's basis selected independently ────────
    cv_for_group <- function(kk, sp_list_k) {
      idx_g <- which(g_loc == kk)
      if (length(idx_g) < 10L)
        return(list(idx = 1L, cvm = rep(NA_real_, length(sp_list_k))))
      fid_g <- fid[idx_g]
      fmm_g <- matrix(Inf, length(sp_list_k), nfolds)
      for (m in seq_along(sp_list_k)) {
        s <- sp_list_k[[m]]
        for (kf in 1:nfolds) {
          tr_g <- idx_g[fid_g != kf]; te_g <- idx_g[fid_g == kf]
          if (length(tr_g) < 5 || length(te_g) < 1) next
          kn <- .get_knots_one(d_loc[tr_g], s)
          if (!kn$ok) next
          B_tr <- tryCatch(create_basis(d_loc[tr_g], s)$basis, error = function(e) NULL)
          if (is.null(B_tr)) next
          B_te <- tryCatch(create_basis(d_loc[te_g], s,
                                        knots_ref = kn$k, bknots_ref = kn$bk)$basis,
                           error = function(e) NULL)
          if (is.null(B_te)) next
          df_tr <- data.frame(Y_out = y_loc[tr_g], B_tr)
          df_te <- data.frame(Y_out = y_loc[te_g], B_te)
          if (!is.null(z_loc) && ncol(z_loc) > 0) {
            df_tr <- cbind(df_tr, z_loc[tr_g, , drop = FALSE])
            df_te <- cbind(df_te, z_loc[te_g, , drop = FALSE])
          }
          if (!is.null(fe_loc) && ncol(fe_loc) > 0) {
            df_tr <- cbind(df_tr, fe_loc[tr_g, , drop = FALSE])
            df_te <- cbind(df_te, fe_loc[te_g, , drop = FALSE])
          }
          fe_tr <- .fe_factors_subset(fe_factors_loc, tr_g)
          fe_te <- .fe_factors_subset(fe_factors_loc, te_g)
          # Training uses survey weights only (not targeted): targeted weights
          # are an AIPE evaluation criterion, not a re-weighted fit.
          yp <- .fit_predict_sel(df_tr, df_te,
                                 if (!is.null(survey_w)) survey_w[tr_g] else NULL,
                                 fe_tr, fe_te)
          if (is.null(yp)) next
          resid2 <- (y_loc[te_g] - yp)^2
          fold_w <- .cv_fold_weight(te_g)
          fmm_g[m, kf] <- if (!is.null(fold_w))
            stats::weighted.mean(resid2, w = fold_w, na.rm = TRUE)
          else mean(resid2, na.rm = TRUE)
        }
      }
      cvm_g <- rowMeans(fmm_g, na.rm = TRUE)
      valid <- which(is.finite(cvm_g))
      bi_g <- if (length(valid) == 0L) 1L else valid[which.min(cvm_g[valid])]
      list(idx = bi_g, cvm = cvm_g)
    }

    if (is_bic) {
      out <- vector("list", K)
      for (kk in 0:(K - 1L)) {
        if (!is.null(user_specs_loc[[kk + 1L]])) {
          out[[kk + 1L]] <- user_specs_loc[[kk + 1L]]; next
        }
        idx_g <- which(g_loc == kk)
        bv <- rep(Inf, length(sp_lists[[kk + 1L]]))
        for (m in seq_along(sp_lists[[kk + 1L]])) {
          s <- sp_lists[[kk + 1L]][[m]]
          kn <- .get_knots_one(d_loc[idx_g], s)
          if (!kn$ok) next
          B <- tryCatch(create_basis(d_loc[idx_g], s)$basis, error = function(e) NULL)
          if (is.null(B)) next
          df_g <- data.frame(Y_out = y_loc[idx_g], B)
          if (!is.null(z_loc) && ncol(z_loc) > 0) df_g <- cbind(df_g, z_loc[idx_g, , drop = FALSE])
          if (!is.null(fe_loc) && ncol(fe_loc) > 0)
            df_g <- cbind(df_g, fe_loc[idx_g, , drop = FALSE])
          fe_g <- .fe_factors_subset(fe_factors_loc, idx_g)
          bv[m] <- .bic_fit_sel(df_g,
                                if (!is.null(survey_w)) survey_w[idx_g] else NULL,
                                fe_g)
        }
        bi_g <- which.min(bv)
        if (length(bi_g) == 0L || !is.finite(bv[bi_g])) bi_g <- 1L
        out[[kk + 1L]] <- sp_lists[[kk + 1L]][[bi_g]]
      }
      if (verbose_sel)
        cat(sprintf("    Selected (BIC, per-group): %s\n",
                    paste(vapply(0:(K - 1L),
                                 function(kk) sprintf("g%d:%s", kk, spec_label(out[[kk + 1L]])),
                                 character(1)), collapse = " | ")))
      return(out)
    }

    # ── Joint CV (K <= 3): grid search over all groups simultaneously ──
    if (joint_select) {
      groups_to_search <- which(vapply(user_specs_loc, is.null, logical(1)))
      grid_list <- lapply(groups_to_search, function(gi) seq_along(sp_lists[[gi]]))
      names(grid_list) <- paste0("g", groups_to_search - 1L)
      grids <- do.call(expand.grid, c(grid_list, list(KEEP.OUT.ATTRS = FALSE)))
      nm_grid <- nrow(grids)
      fmm <- matrix(Inf, nm_grid, nfolds)
      for (mi in 1:nm_grid) {
        sp_try <- user_specs_loc
        for (gi in groups_to_search)
          sp_try[[gi]] <- sp_lists[[gi]][[grids[mi, paste0("g", gi - 1L)]]]
        for (kf in 1:nfolds) {
          tr <- fid != kf; te <- fid == kf
          kns_try <- lapply(0:(K - 1L), function(kk)
            .get_knots_one(d_loc[tr & g_loc == kk], sp_try[[kk + 1L]]))
          if (any(!vapply(kns_try, function(z) z$ok, logical(1)))) next
          dm_tr <- tryCatch(.build_design_K(d_loc[tr], g_loc[tr],
                                            if (!is.null(z_loc)) z_loc[tr, , drop = FALSE] else NULL,
                                            sp_try, kns_try,
                                            if (!is.null(fe_loc)) fe_loc[tr, , drop = FALSE] else NULL,
                                            x_raw_loc = if (!is.null(x_raw_loc)) x_raw_loc[tr] else NULL),
                            error = function(e) NULL)
          if (is.null(dm_tr)) next
          dm_te <- tryCatch(.build_design_K(d_loc[te], g_loc[te],
                                            if (!is.null(z_loc)) z_loc[te, , drop = FALSE] else NULL,
                                            sp_try, kns_try,
                                            if (!is.null(fe_loc)) fe_loc[te, , drop = FALSE] else NULL,
                                            x_raw_loc = if (!is.null(x_raw_loc)) x_raw_loc[te] else NULL),
                            error = function(e) NULL)
          if (is.null(dm_te)) next
          df_tr_j <- data.frame(Y_out = y_loc[tr], dm_tr)
          df_te_j <- data.frame(Y_out = y_loc[te], dm_te)
          fe_tr_j <- .fe_factors_subset(fe_factors_loc, which(tr))
          fe_te_j <- .fe_factors_subset(fe_factors_loc, which(te))
          yp <- .fit_predict_sel(df_tr_j, df_te_j,
                                 if (!is.null(survey_w)) survey_w[tr] else NULL,
                                 fe_tr_j, fe_te_j)
          if (is.null(yp)) next
          resid2 <- (y_loc[te] - yp)^2
          fold_w <- .cv_fold_weight(which(te))
          fmm[mi, kf] <- if (!is.null(fold_w))
            stats::weighted.mean(resid2, w = fold_w, na.rm = TRUE)
          else mean(resid2, na.rm = TRUE)
        }
      }
      cvm <- rowMeans(fmm, na.rm = TRUE)
      valid <- which(is.finite(cvm))
      bi <- if (length(valid) == 0L) 1L else valid[which.min(cvm[valid])]
      out <- user_specs_loc
      for (gi in groups_to_search)
        out[[gi]] <- sp_lists[[gi]][[grids[bi, paste0("g", gi - 1L)]]]
      if (verbose_sel)
        cat(sprintf("    Selected (joint CV): %s\n",
                    paste(vapply(0:(K - 1L),
                                 function(kk) sprintf("g%d:%s", kk, spec_label(out[[kk + 1L]])),
                                 character(1)), collapse = " | ")))
      return(out)
    }

    # ── Default: per-group independent CV ──────────────────────────────
    out <- vector("list", K)
    for (kk in 0:(K - 1L)) {
      if (!is.null(user_specs_loc[[kk + 1L]])) {
        out[[kk + 1L]] <- user_specs_loc[[kk + 1L]]
      } else {
        r <- cv_for_group(kk, sp_lists[[kk + 1L]])
        out[[kk + 1L]] <- sp_lists[[kk + 1L]][[r$idx]]
      }
    }
    if (verbose_sel)
      cat(sprintf("    Selected (per-group CV): %s\n",
                  paste(vapply(0:(K - 1L),
                               function(kk) sprintf("g%d:%s", kk, spec_label(out[[kk + 1L]])),
                               character(1)), collapse = " | ")))
    out
  }

  # ══════════════════════════════════════════════════════════════════════════
  # ESTIMATE GIVEN K SPECS
  # ══════════════════════════════════════════════════════════════════════════

  .estimate_K <- function(y_loc, d_loc, g_loc, z_loc, specs_list,
                          vce_type, cl_var, survey_w = NULL,
                          fe_loc = NULL,
                          fe_factors = NULL,
                          x_raw_loc = NULL,
                          d_all_for_aipe_local = NULL,
                          survey_w_all_for_aipe_local = NULL) {
    kns_list <- lapply(0:(K - 1L),
                       function(kk) .get_knots_one(d_loc[g_loc == kk], specs_list[[kk + 1L]]))

    design_mat <- .build_design_K(d_loc, g_loc, z_loc, specs_list, kns_list,
                                  fe_loc, x_raw_loc = x_raw_loc)

    # Group-specific mean of X used to evaluate the alpha_k * X * G_k terms
    # at "average composition" for predictions (only when continuous X).
    x_grp_means <- if (!is.null(x_raw_loc))
      vapply(0:(K - 1L),
             function(kk) mean(x_raw_loc[g_loc == kk], na.rm = TRUE),
             numeric(1))
    else NULL
    fit_out <- .fit_lm(y_loc, design_mat, survey_w = survey_w,
                       cl_var = cl_var, vce_type = vce_type,
                       fe_factors = fe_factors)
    fit <- fit_out$fit; V <- fit_out$V; bh_full <- fit_out$bh
    keep_coef <- !is.na(bh_full)
    bh <- bh_full[keep_coef]
    if (ncol(V) != length(bh)) {
      common <- intersect(names(bh), colnames(V))
      bh <- bh[common]
      V  <- V[common, common, drop = FALSE]
    }
    p_dim <- length(bh)
    dg <- seq(cs_lo, cs_hi, length.out = n_grid)
    # If we added X*G blocks, drop X_continuous from the Z-mean vector to match
    # the design built in .build_design_K (avoids double counting X).
    z_loc_eff <- z_loc
    if (!is.null(x_raw_loc) && !is.null(z_loc_eff)) {
      drop_idx <- which(colnames(z_loc_eff) == "X_continuous")
      if (length(drop_idx) > 0L)
        z_loc_eff <- z_loc_eff[, -drop_idx, drop = FALSE]
      if (ncol(z_loc_eff) == 0L) z_loc_eff <- NULL
    }
    zm  <- if (!is.null(z_loc_eff) && ncol(z_loc_eff) > 0) colMeans(z_loc_eff, na.rm = TRUE) else NULL
    fem <- if (!is.null(fe_loc) && ncol(fe_loc) > 0) colMeans(fe_loc, na.rm = TRUE) else NULL

    # Prediction-row matrix at (d_vec, group = kk).
    # For the alpha_k * X * I(g=k) blocks, predict at the group-kk mean of X
    # (i.e., "average composition within the group"). Only the X_g{kk} column
    # is nonzero on this prediction surface; others stay 0 because we are
    # plotting the group-kk curve.
    .bpm_at_group <- function(d_vec, kk) {
      nd <- length(d_vec)
      basis_blocks <- vector("list", K)
      for (kj in 0:(K - 1L)) {
        bk <- create_basis(d_vec, specs_list[[kj + 1L]],
                           kns_list[[kj + 1L]]$k, kns_list[[kj + 1L]]$bk)
        block <- bk$basis * as.numeric(kj == kk)
        colnames(block) <- paste0("B_g", kj, "_", colnames(bk$basis))
        basis_blocks[[kj + 1L]] <- block
      }
      mt <- do.call(cbind, basis_blocks)
      if (K > 1L) {
        gd <- matrix(0, nd, K - 1L)
        for (kj in 1:(K - 1L)) gd[, kj] <- as.numeric(kj == kk)
        colnames(gd) <- paste0("G_", 1:(K - 1L))
        mt <- cbind(mt, gd)
      }
      if (!is.null(x_raw_loc)) {
        x_blk <- matrix(0, nd, K)
        x_blk[, kk + 1L] <- x_grp_means[kk + 1L]
        colnames(x_blk) <- paste0("X_g", 0:(K - 1L))
        mt <- cbind(mt, x_blk)
      }
      if (!is.null(zm)) {
        if (Z_interact) {
          for (kj in 0:(K - 1L)) {
            zk <- matrix(rep(zm * as.numeric(kj == kk), each = nd), nrow = nd)
            colnames(zk) <- paste0("Z_g", kj, "_", colnames(z_loc_eff))
            mt <- cbind(mt, zk)
          }
        } else {
          zz <- matrix(rep(zm, each = nd), nrow = nd)
          colnames(zz) <- paste0("Z_", colnames(z_loc_eff))
          mt <- cbind(mt, zz)
        }
      }
      if (!is.null(fem)) {
        fe_blk <- matrix(rep(fem, each = nd), nrow = nd)
        colnames(fe_blk) <- names(fem)
        mt <- cbind(mt, fe_blk)
      }
      mt
    }
    .align_to_bh <- function(mt, bn) {
      al <- matrix(0, nrow(mt), length(bn)); colnames(al) <- bn
      sh <- intersect(colnames(mt), bn); al[, sh] <- mt[, sh, drop = FALSE]; al
    }
    .add_intercept <- function(mt, bn) {
      if ("(Intercept)" %in% bn) cbind(`(Intercept)` = 1, mt) else mt
    }

    # Predicted Y and SE at the K group "endpoints" on the D grid.
    yhat_by_group <- vector("list", K)
    se_yhat_by_group <- vector("list", K)
    for (kk in 0:(K - 1L)) {
      mt_full <- .align_to_bh(.add_intercept(.bpm_at_group(dg, kk), names(bh)), names(bh))
      yhat_by_group[[kk + 1L]] <- as.numeric(mt_full %*% bh)
      se_yhat_by_group[[kk + 1L]] <- sqrt(pmax(0, rowSums((mt_full %*% V) * mt_full)))
    }

    # When feols absorbed FE, bh has no intercept and the curve sits at
    # "FE = 0". Shift it to the average FE composition so the plotted Y
    # is at the sample mean. The shift is a constant => slopes, CAME,
    # AIPE, DCAME, and the SE bands are unaffected.
    if (!is.null(fe_factors)) {
      design_in <- .align_to_bh(.add_intercept(design_mat, names(bh)), names(bh))
      pred_in   <- as.numeric(design_in %*% bh)
      if (!is.null(survey_w)) {
        sw_n  <- survey_w / sum(survey_w)
        shift <- sum(sw_n * y_loc) - sum(sw_n * pred_in)
      } else {
        shift <- mean(y_loc) - mean(pred_in)
      }
      for (kk in 0:(K - 1L))
        yhat_by_group[[kk + 1L]] <- yhat_by_group[[kk + 1L]] + shift
    }

    # Slope (finite difference) at each grid point for each group.
    h <- diff(range(dg)) / 1e5
    Sg_by_group <- vector("list", K)
    slope_by_group <- vector("list", K)
    se_slope_by_group <- vector("list", K)
    for (kk in 0:(K - 1L)) {
      m1 <- .align_to_bh(.add_intercept(.bpm_at_group(dg + h, kk), names(bh)), names(bh))
      m0 <- .align_to_bh(.add_intercept(.bpm_at_group(dg - h, kk), names(bh)), names(bh))
      Sg <- (m1 - m0) / (2 * h)
      Sg_by_group[[kk + 1L]] <- Sg
      slope_by_group[[kk + 1L]] <- as.numeric(Sg %*% bh)
      se_slope_by_group[[kk + 1L]] <- sqrt(pmax(0, rowSums((Sg %*% V) * Sg)))
    }

    # IPE: high-minus-low slope contrast.
    Cg <- Sg_by_group[[K]] - Sg_by_group[[1L]]
    ipeg <- as.numeric(Cg %*% bh)
    seig <- sqrt(pmax(0, rowSums((Cg %*% V) * Cg)))

    # Interpolation onto arbitrary D values along an n_grid-row matrix.
    .interp_to_d <- function(Gmat, d_vals) {
      n_d <- length(d_vals); p_g <- ncol(Gmat)
      Co <- matrix(0, n_d, p_g)
      for (i in 1:n_d) {
        pos <- (d_vals[i] - cs_lo) / (cs_hi - cs_lo) * (n_grid - 1) + 1
        pos <- max(1, min(n_grid, pos))
        lo_idx <- floor(pos); hi_idx <- ceiling(pos)
        if (lo_idx == hi_idx || lo_idx >= n_grid) {
          Co[i, ] <- Gmat[min(lo_idx, n_grid), ]
        } else {
          w <- pos - lo_idx
          Co[i, ] <- (1 - w) * Gmat[lo_idx, ] + w * Gmat[hi_idx, ]
        }
      }
      Co
    }

    # AIPE: integrate IPE over full-sample D distribution (clamped to CS).
    d_for_aipe <- if (!is.null(d_all_for_aipe_local)) d_all_for_aipe_local else d_loc
    sw_for_aipe <- if (!is.null(d_all_for_aipe_local) && !is.null(survey_w_all_for_aipe_local))
      survey_w_all_for_aipe_local else survey_w
    d_for_aipe_c <- pmin(pmax(d_for_aipe, cs_lo), cs_hi)
    Co_ipe_all <- .interp_to_d(Cg, d_for_aipe_c)
    if (!is.null(sw_for_aipe)) {
      sw_n <- sw_for_aipe / sum(sw_for_aipe)
      Cb_aipe <- as.numeric(t(Co_ipe_all) %*% sw_n)
    } else {
      Cb_aipe <- colMeans(Co_ipe_all)
    }
    aipe <- as.numeric(Cb_aipe %*% bh)
    seA <- sqrt(max(0, as.numeric(t(Cb_aipe) %*% V %*% Cb_aipe)))

    # CAME_k for each group: average slope over group k's D distribution.
    came_by_group <- numeric(K)
    se_came_by_group <- numeric(K)
    Cb_came_by_group <- vector("list", K)
    for (kk in 0:(K - 1L)) {
      idx_kk <- which(g_loc == kk)
      if (length(idx_kk) == 0L) {
        came_by_group[kk + 1L] <- NA_real_
        se_came_by_group[kk + 1L] <- NA_real_
        Cb_came_by_group[[kk + 1L]] <- rep(0, p_dim)
        next
      }
      Co_kk <- .interp_to_d(Sg_by_group[[kk + 1L]], d_loc[idx_kk])
      if (!is.null(survey_w)) {
        sw_kk <- survey_w[idx_kk]; sw_kk <- sw_kk / sum(sw_kk)
        Cb_kk <- as.numeric(t(Co_kk) %*% sw_kk)
      } else {
        Cb_kk <- colMeans(Co_kk)
      }
      came_by_group[kk + 1L] <- as.numeric(Cb_kk %*% bh)
      se_came_by_group[kk + 1L] <- sqrt(max(0, as.numeric(t(Cb_kk) %*% V %*% Cb_kk)))
      Cb_came_by_group[[kk + 1L]] <- Cb_kk
    }

    came0 <- came_by_group[1L]; came1 <- came_by_group[K]
    se_came0 <- se_came_by_group[1L]; se_came1 <- se_came_by_group[K]
    Cb_dcame <- Cb_came_by_group[[K]] - Cb_came_by_group[[1L]]
    dcame <- came1 - came0
    se_dcame <- sqrt(max(0, as.numeric(t(Cb_dcame) %*% V %*% Cb_dcame)))

    list(fit = fit, vcov = V, coef = bh,
         d_grid = dg,
         yhat_X0 = yhat_by_group[[1L]], yhat_X1 = yhat_by_group[[K]],
         se_yhat_X0 = se_yhat_by_group[[1L]], se_yhat_X1 = se_yhat_by_group[[K]],
         slope_grid_X0 = slope_by_group[[1L]], slope_grid_X1 = slope_by_group[[K]],
         se_slope_X0 = se_slope_by_group[[1L]], se_slope_X1 = se_slope_by_group[[K]],
         ipe_grid = ipeg, se_ipe_grid = seig,
         aipe = aipe, se_aipe = seA,
         came0 = came0, came1 = came1,
         se_came0 = se_came0, se_came1 = se_came1,
         dcame = dcame, se_dcame = se_dcame,
         came_by_group = came_by_group,
         se_came_by_group = se_came_by_group,
         specs_list = specs_list)
  }

  # ══════════════════════════════════════════════════════════════════════════
  # BOOTSTRAP and CROSSFIT WRAPPERS
  # ══════════════════════════════════════════════════════════════════════════

  .bootstrap_K <- function(y_loc, d_loc, g_loc, z_loc, user_specs_loc,
                           sel_method, vce_type, cl_var, survey_w,
                           fe_loc, fe_factors_loc = NULL,
                           x_raw_loc_full = NULL, verbose_boot = FALSE) {
    n_bt <- length(y_loc)
    if (verbose_boot) cat(sprintf("--- Bootstrap inference (%d replicates) ---\n", B_boot))

    specs_full <- .select_model_K(y_loc, d_loc, g_loc, z_loc, user_specs_loc,
                                  sel_method, survey_w = survey_w,
                                  verbose_sel = FALSE, fe_loc = fe_loc,
                                  fe_factors_loc = fe_factors_loc,
                                  x_raw_loc = x_raw_loc_full)
    est_full <- .estimate_K(y_loc, d_loc, g_loc, z_loc, specs_full,
                            vce_type, cl_var, survey_w, fe_loc = fe_loc,
                            fe_factors = fe_factors_loc,
                            x_raw_loc = x_raw_loc_full,
                            d_all_for_aipe_local = d_all_for_aipe,
                            survey_w_all_for_aipe_local = survey_wgts_all_for_aipe)

    boot_aipe   <- numeric(B_boot)
    boot_dcame  <- numeric(B_boot)
    boot_came0  <- numeric(B_boot)
    boot_came1  <- numeric(B_boot)
    boot_came_all <- matrix(NA_real_, B_boot, K)
    boot_ipe    <- matrix(NA, B_boot, n_grid)
    boot_yh0    <- matrix(NA, B_boot, n_grid)
    boot_yh1    <- matrix(NA, B_boot, n_grid)
    boot_slope0 <- matrix(NA, B_boot, n_grid)
    boot_slope1 <- matrix(NA, B_boot, n_grid)

    for (b in 1:B_boot) {
      if (verbose_boot && b %% 50 == 0)
        cat(sprintf("  Bootstrap replicate %d/%d\n", b, B_boot))
      if (!is.null(cl_var)) {
        ucl <- unique(cl_var)
        idx <- unlist(lapply(sample(ucl, replace = TRUE),
                             function(cl) which(cl_var == cl)))
      } else {
        idx <- sample(n_bt, replace = TRUE)
      }
      yb <- y_loc[idx]; db <- d_loc[idx]; gb <- g_loc[idx]
      zb <- if (!is.null(z_loc)) z_loc[idx, , drop = FALSE] else NULL
      cb <- if (!is.null(cl_var)) cl_var[idx] else NULL
      wb <- if (!is.null(survey_w)) survey_w[idx] else NULL
      feb <- if (!is.null(fe_loc)) fe_loc[idx, , drop = FALSE] else NULL
      if (!is.null(feb) && ncol(feb) > 0) {
        sds_feb <- apply(feb, 2, sd, na.rm = TRUE)
        feb <- feb[, sds_feb > 0, drop = FALSE]
        if (ncol(feb) == 0) feb <- NULL
      }
      fefb <- if (!is.null(fe_factors_loc))
        lapply(fe_factors_loc, function(f) factor(as.character(f[idx])))
      else NULL
      xrb <- if (!is.null(x_raw_loc_full)) x_raw_loc_full[idx] else NULL

      specs_b <- tryCatch(
        .select_model_K(yb, db, gb, zb, user_specs_loc, sel_method,
                        survey_w = wb, verbose_sel = FALSE, fe_loc = feb,
                        fe_factors_loc = fefb, x_raw_loc = xrb),
        error = function(e) NULL)
      if (is.null(specs_b)) next

      est_b <- tryCatch(
        .estimate_K(yb, db, gb, zb, specs_b, vce_type, cb, wb, fe_loc = feb,
                    fe_factors = fefb, x_raw_loc = xrb),
        error = function(e) NULL)
      if (is.null(est_b)) next

      boot_aipe[b]    <- est_b$aipe
      boot_dcame[b]   <- est_b$dcame
      boot_came0[b]   <- est_b$came0
      boot_came1[b]   <- est_b$came1
      boot_came_all[b, ] <- est_b$came_by_group
      boot_ipe[b, ]   <- est_b$ipe_grid
      boot_yh0[b, ]   <- est_b$yhat_X0
      boot_yh1[b, ]   <- est_b$yhat_X1
      boot_slope0[b, ] <- est_b$slope_grid_X0
      boot_slope1[b, ] <- est_b$slope_grid_X1
    }

    vld <- !is.na(boot_aipe)
    est_full$se_aipe   <- sd(boot_aipe[vld])
    est_full$se_dcame  <- sd(boot_dcame[vld])
    est_full$se_came0  <- sd(boot_came0[vld])
    est_full$se_came1  <- sd(boot_came1[vld])
    est_full$se_came_by_group <- apply(boot_came_all[vld, , drop = FALSE], 2, sd, na.rm = TRUE)
    est_full$se_ipe_grid <- apply(boot_ipe[vld, , drop = FALSE], 2, sd, na.rm = TRUE)
    est_full$se_yhat_X0 <- apply(boot_yh0[vld, , drop = FALSE], 2, sd, na.rm = TRUE)
    est_full$se_yhat_X1 <- apply(boot_yh1[vld, , drop = FALSE], 2, sd, na.rm = TRUE)
    est_full$se_slope_X0 <- apply(boot_slope0[vld, , drop = FALSE], 2, sd, na.rm = TRUE)
    est_full$se_slope_X1 <- apply(boot_slope1[vld, , drop = FALSE], 2, sd, na.rm = TRUE)
    est_full
  }

  .crossfit_K <- function(y_loc, d_loc, g_loc, z_loc, user_specs_loc,
                          sel_method, vce_type, cl_var, survey_w,
                          fe_loc, fe_factors_loc = NULL,
                          x_raw_loc_full = NULL, verbose_cf = FALSE) {
    n_cf <- length(y_loc)
    if (verbose_cf) cat("--- Cross-fitting inference (all quantities) ---\n")

    .old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(42)
    half1 <- integer(0)
    for (kk in 0:(K - 1L)) {
      idx_kk <- which(g_loc == kk)
      s_kk <- sample(idx_kk)
      half1 <- c(half1, s_kk[seq_len(floor(length(s_kk) / 2))])
    }
    half1 <- sort(half1)
    half2 <- setdiff(1:n_cf, half1)
    if (!is.null(.old_seed))
      assign(".Random.seed", .old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv))
      rm(".Random.seed", envir = .GlobalEnv)

    results <- list()
    for (fold in 1:2) {
      sel_idx <- if (fold == 1) half1 else half2
      est_idx <- if (fold == 1) half2 else half1

      if (verbose_cf) cat(sprintf("  Fold %d: selection N=%d, estimation N=%d\n",
                                  fold, length(sel_idx), length(est_idx)))

      fe_sel <- if (!is.null(fe_loc)) fe_loc[sel_idx, , drop = FALSE] else NULL
      fe_est <- if (!is.null(fe_loc)) fe_loc[est_idx, , drop = FALSE] else NULL
      if (!is.null(fe_sel) && ncol(fe_sel) > 0) {
        sds <- apply(fe_sel, 2, sd, na.rm = TRUE)
        fe_sel <- fe_sel[, sds > 0, drop = FALSE]
        if (ncol(fe_sel) == 0) fe_sel <- NULL
      }
      if (!is.null(fe_est) && ncol(fe_est) > 0) {
        sds <- apply(fe_est, 2, sd, na.rm = TRUE)
        fe_est <- fe_est[, sds > 0, drop = FALSE]
        if (ncol(fe_est) == 0) fe_est <- NULL
      }
      fefac_sel <- if (!is.null(fe_factors_loc))
        lapply(fe_factors_loc, function(f) factor(as.character(f[sel_idx])))
      else NULL
      fefac_est <- if (!is.null(fe_factors_loc))
        lapply(fe_factors_loc, function(f) factor(as.character(f[est_idx])))
      else NULL
      xr_sel <- if (!is.null(x_raw_loc_full)) x_raw_loc_full[sel_idx] else NULL
      xr_est <- if (!is.null(x_raw_loc_full)) x_raw_loc_full[est_idx] else NULL

      specs_fold <- .select_model_K(
        y_loc[sel_idx], d_loc[sel_idx], g_loc[sel_idx],
        if (!is.null(z_loc)) z_loc[sel_idx, , drop = FALSE] else NULL,
        user_specs_loc, sel_method,
        survey_w = if (!is.null(survey_w)) survey_w[sel_idx] else NULL,
        verbose_sel = verbose_cf, fe_loc = fe_sel,
        fe_factors_loc = fefac_sel, x_raw_loc = xr_sel)

      est <- .estimate_K(
        y_loc[est_idx], d_loc[est_idx], g_loc[est_idx],
        if (!is.null(z_loc)) z_loc[est_idx, , drop = FALSE] else NULL,
        specs_fold, vce_type,
        if (!is.null(cl_var)) cl_var[est_idx] else NULL,
        if (!is.null(survey_w)) survey_w[est_idx] else NULL,
        fe_loc = fe_est,
        fe_factors = fefac_est,
        x_raw_loc = xr_est,
        d_all_for_aipe_local = d_all_for_aipe,
        survey_w_all_for_aipe_local = survey_wgts_all_for_aipe)

      results[[fold]] <- list(est = est, specs = specs_fold)
    }

    w1 <- length(half2) / n_cf; w2 <- length(half1) / n_cf
    r1 <- results[[1]]$est; r2 <- results[[2]]$est

    list(
      aipe = w1 * r1$aipe + w2 * r2$aipe,
      se_aipe = sqrt(w1^2 * r1$se_aipe^2 + w2^2 * r2$se_aipe^2),
      dcame = w1 * r1$dcame + w2 * r2$dcame,
      se_dcame = sqrt(w1^2 * r1$se_dcame^2 + w2^2 * r2$se_dcame^2),
      came0 = w1 * r1$came0 + w2 * r2$came0,
      came1 = w1 * r1$came1 + w2 * r2$came1,
      se_came0 = sqrt(w1^2 * r1$se_came0^2 + w2^2 * r2$se_came0^2),
      se_came1 = sqrt(w1^2 * r1$se_came1^2 + w2^2 * r2$se_came1^2),
      came_by_group = w1 * r1$came_by_group + w2 * r2$came_by_group,
      se_came_by_group = sqrt(w1^2 * r1$se_came_by_group^2 + w2^2 * r2$se_came_by_group^2),
      d_grid = r1$d_grid,
      yhat_X0 = w1 * r1$yhat_X0 + w2 * r2$yhat_X0,
      yhat_X1 = w1 * r1$yhat_X1 + w2 * r2$yhat_X1,
      se_yhat_X0 = sqrt(w1^2 * r1$se_yhat_X0^2 + w2^2 * r2$se_yhat_X0^2),
      se_yhat_X1 = sqrt(w1^2 * r1$se_yhat_X1^2 + w2^2 * r2$se_yhat_X1^2),
      ipe_grid = w1 * r1$ipe_grid + w2 * r2$ipe_grid,
      se_ipe_grid = sqrt(w1^2 * r1$se_ipe_grid^2 + w2^2 * r2$se_ipe_grid^2),
      slope_grid_X0 = w1 * r1$slope_grid_X0 + w2 * r2$slope_grid_X0,
      slope_grid_X1 = w1 * r1$slope_grid_X1 + w2 * r2$slope_grid_X1,
      se_slope_X0 = sqrt(w1^2 * r1$se_slope_X0^2 + w2^2 * r2$se_slope_X0^2),
      se_slope_X1 = sqrt(w1^2 * r1$se_slope_X1^2 + w2^2 * r2$se_slope_X1^2),
      specs_fold1 = results[[1]]$specs, specs_fold2 = results[[2]]$specs
    )
  }

  # ══════════════════════════════════════════════════════════════════════════
  # DISPATCH
  # ══════════════════════════════════════════════════════════════════════════

  fit_obj <- NULL; V_obj <- NULL

  if (inference == "crossfit") {
    cf <- .crossfit_K(y, d, g, z_mat, user_specs, selection, vce, cluster_var,
                      survey_w = survey_wgts, fe_loc = fe_mat,
                      fe_factors_loc = fe_absorb_factors,
                      x_raw_loc_full = x_continuous_raw,
                      verbose_cf = verbose)
    dg <- cf$d_grid
    aipe <- cf$aipe; seA <- cf$se_aipe
    dcame <- cf$dcame; se_dcame <- cf$se_dcame
    came0 <- cf$came0; came1 <- cf$came1
    se_came0 <- cf$se_came0; se_came1 <- cf$se_came1
    came_all <- cf$came_by_group; se_came_all <- cf$se_came_by_group
    ipeg <- cf$ipe_grid; seig <- cf$se_ipe_grid
    yh0 <- cf$yhat_X0; yh1 <- cf$yhat_X1
    se0 <- cf$se_yhat_X0; se1 <- cf$se_yhat_X1
    slope0g <- cf$slope_grid_X0; slope1g <- cf$slope_grid_X1
    se_slope0g <- cf$se_slope_X0; se_slope1g <- cf$se_slope_X1
    bX0 <- cf$specs_fold1[[1L]]; bX1 <- cf$specs_fold1[[K]]
    bX0_label <- sprintf("%s / %s",
                         spec_label(cf$specs_fold1[[1L]]),
                         spec_label(cf$specs_fold2[[1L]]))
    bX1_label <- sprintf("%s / %s",
                         spec_label(cf$specs_fold1[[K]]),
                         spec_label(cf$specs_fold2[[K]]))

  } else if (inference == "bootstrap") {
    br <- .bootstrap_K(y, d, g, z_mat, user_specs, selection, vce, cluster_var,
                       survey_w = survey_wgts, fe_loc = fe_mat,
                       fe_factors_loc = fe_absorb_factors,
                       x_raw_loc_full = x_continuous_raw,
                       verbose_boot = verbose)
    dg <- br$d_grid
    aipe <- br$aipe; seA <- br$se_aipe
    dcame <- br$dcame; se_dcame <- br$se_dcame
    came0 <- br$came0; came1 <- br$came1
    se_came0 <- br$se_came0; se_came1 <- br$se_came1
    came_all <- br$came_by_group; se_came_all <- br$se_came_by_group
    ipeg <- br$ipe_grid; seig <- br$se_ipe_grid
    yh0 <- br$yhat_X0; yh1 <- br$yhat_X1
    se0 <- br$se_yhat_X0; se1 <- br$se_yhat_X1
    slope0g <- br$slope_grid_X0; slope1g <- br$slope_grid_X1
    se_slope0g <- br$se_slope_X0; se_slope1g <- br$se_slope_X1
    bX0 <- br$specs_list[[1L]]; bX1 <- br$specs_list[[K]]
    bX0_label <- spec_label(bX0); bX1_label <- spec_label(bX1)
    fit_obj <- br$fit; V_obj <- br$vcov

  } else {
    if (verbose) cat("--- Model selection (regular) ---\n")
    specs_final <- .select_model_K(y, d, g, z_mat, user_specs, selection,
                                   survey_w = survey_wgts, verbose_sel = verbose,
                                   fe_loc = fe_mat,
                                   fe_factors_loc = fe_absorb_factors,
                                   x_raw_loc = x_continuous_raw)
    if (verbose) cat("--- Final estimation ---\n")
    est <- .estimate_K(y, d, g, z_mat, specs_final, vce, cluster_var, survey_wgts,
                       fe_loc = fe_mat,
                       fe_factors = fe_absorb_factors,
                       x_raw_loc = x_continuous_raw,
                       d_all_for_aipe_local = d_all_for_aipe,
                       survey_w_all_for_aipe_local = survey_wgts_all_for_aipe)
    dg <- est$d_grid
    aipe <- est$aipe; seA <- est$se_aipe
    dcame <- est$dcame; se_dcame <- est$se_dcame
    came0 <- est$came0; came1 <- est$came1
    se_came0 <- est$se_came0; se_came1 <- est$se_came1
    came_all <- est$came_by_group; se_came_all <- est$se_came_by_group
    ipeg <- est$ipe_grid; seig <- est$se_ipe_grid
    yh0 <- est$yhat_X0; yh1 <- est$yhat_X1
    se0 <- est$se_yhat_X0; se1 <- est$se_yhat_X1
    slope0g <- est$slope_grid_X0; slope1g <- est$slope_grid_X1
    se_slope0g <- est$se_slope_X0; se_slope1g <- est$se_slope_X1
    bX0 <- est$specs_list[[1L]]; bX1 <- est$specs_list[[K]]
    bX0_label <- spec_label(bX0); bX1_label <- spec_label(bX1)
    fit_obj <- est$fit; V_obj <- est$vcov
  }

  # ── Per-group CAME table (used by plot_came_all and reporting) ───────
  came_all_df <- NULL
  if (!x_is_binary && isTRUE(plot_came_all)) {
    bin_labels <- vapply(0:(K - 1L), function(kk) {
      if (kk == 0L)
        sprintf("X <= %.3g", tercile_cuts[1])
      else if (kk == K - 1L)
        sprintf("X > %.3g", tercile_cuts[length(tercile_cuts)])
      else
        sprintf("%.3g < X <= %.3g", tercile_cuts[kk], tercile_cuts[kk + 1L])
    }, character(1))
    x_mid_per_group <- vapply(0:(K - 1L), function(kk)
      stats::median(x_continuous_raw[g == kk], na.rm = TRUE), numeric(1))
    came_all_df <- data.frame(
      group = 0:(K - 1L),
      group_label = bin_labels,
      X_mid = x_mid_per_group,
      n = vapply(0:(K - 1L), function(kk) sum(g == kk), integer(1)),
      CAME = came_all,
      SE = se_came_all,
      CI_lo = came_all - zc_est * se_came_all,
      CI_hi = came_all + zc_est * se_came_all,
      p_value = 2 * stats::pnorm(-abs(came_all / se_came_all)),
      stringsAsFactors = FALSE
    )
  } else if (x_is_binary && isTRUE(plot_came_all) && verbose) {
    cat("Note: plot_came_all is ignored when the moderator is binary.\n")
  }

  # ══════════════════════════════════════════════════════════════════════════
  # RENDER
  # ══════════════════════════════════════════════════════════════════════════

  n_fe_total <- 0L
  if (!is.null(fe_mat)) n_fe_total <- ncol(fe_mat)
  if (Z_FE_absorb && !is.null(Z_FE)) n_fe_total <- absorbed_df_lost

  .render_results(
    model = "basis",
    estimand = estimand, compare = compare,
    level_pv = level_pv, level_est = level_est,
    zc_pv = zc_pv, zc_est = zc_est,
    aipe = aipe, seA = seA,
    dcame = dcame, se_dcame = se_dcame,
    came0 = came0, came1 = came1,
    se_came0 = se_came0, se_came1 = se_came1,
    dg = dg, yh0 = yh0, yh1 = yh1, se0 = se0, se1 = se1,
    slope0g = slope0g, slope1g = slope1g,
    se_slope0g = se_slope0g, se_slope1g = se_slope1g,
    ipeg = ipeg, seig = seig,
    d_obs_X0 = d_obs_X0, d_obs_X1 = d_obs_X1,
    d_all = d,
    y = y, x = g, z_mat = z_mat,
    x_raw_full = x_continuous_raw,
    cluster_var = cluster_var, survey_wgts = survey_wgts, vce = vce,
    x_is_binary = x_is_binary, tercile_cuts = tercile_cuts,
    n_dropped_t2 = 0L,
    bX0_label = bX0_label, bX1_label = bX1_label,
    bX0 = bX0, bX1 = bX1,
    fit_obj = fit_obj, V_obj = V_obj,
    cs_lo = cs_lo, cs_hi = cs_hi,
    nt = nt, n = n,
    n_grid = n_grid, hist_pv = hist_pv,
    inference = inference, selection = selection,
    Z_interact = Z_interact,
    FE = FE, FE_vbles = FE_vbles,
    wgt = wgt,
    Z_FE = Z_FE, Z_FE_absorb = Z_FE_absorb,
    fe_absorb_factors = fe_absorb_factors,
    n_fe_dummies = n_fe_total,
    fe_mat = fe_mat,
    x_original_labels = x_original_labels,
    came_all_df = came_all_df,
    plot_came_all = isTRUE(plot_came_all),
    K = K
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SHARED RENDERER (summary, plots, return object)
# ══════════════════════════════════════════════════════════════════════════════

#' Shared renderer producing the final S3 object returned by dcame_aipe().
#' @keywords internal
#' @noRd
.render_results <- function(model, estimand, compare,
                            level_pv, level_est,
                            zc_pv, zc_est,
                            aipe, seA, dcame, se_dcame, came0, came1,
                            se_came0, se_came1,
                            dg, yh0, yh1, se0, se1,
                            slope0g, slope1g, se_slope0g, se_slope1g,
                            ipeg, seig,
                            d_obs_X0, d_obs_X1, d_all,
                            y, x, z_mat, cluster_var, survey_wgts, vce,
                            x_raw_full = NULL,
                            x_is_binary, tercile_cuts, n_dropped_t2,
                            bX0_label, bX1_label, bX0, bX1,
                            fit_obj, V_obj, cs_lo, cs_hi, nt, n,
                            n_grid, hist_pv,
                            inference, selection,
                            Z_interact, FE, FE_vbles = NULL,
                            wgt,
                            Z_FE = NULL, Z_FE_absorb = FALSE,
                            fe_absorb_factors = NULL,
                            n_fe_dummies = 0L, fe_mat = NULL,
                            x_original_labels = NULL,
                            came_all_df = NULL,
                            plot_came_all = FALSE,
                            K = 2L) {

  do_linear <- !isFALSE(compare) && compare %in% c("LINEAR", "ALL")
  do_dcame_compare <- !isFALSE(compare) && compare %in% c("DCAME", "ALL")
  do_aipe_compare  <- !isFALSE(compare) && compare %in% c("AIPE", "ALL")

  if (estimand == "DCAME") {
    cat(sprintf("\n  D-CAME = %.4f (SE = %.4f) %d%% CI: [%.4f, %.4f]\n",
                dcame, se_dcame, round(level_est * 100),
                dcame - zc_est * se_dcame, dcame + zc_est * se_dcame))
    cat(sprintf("  CAME(X=0) = %.4f (SE = %.4f) | CAME(X=1) = %.4f (SE = %.4f)\n",
                came0, se_came0, came1, se_came1))
  } else {
    cat(sprintf("\n  AIPE = %.4f (SE = %.4f) %d%% CI: [%.4f, %.4f]\n",
                aipe, seA, round(level_est * 100),
                aipe - zc_est * seA, aipe + zc_est * seA))
  }

  aiL <- NULL; seL <- NULL

  if (do_aipe_compare && estimand == "DCAME") {
    cat(sprintf("  AIPE = %.4f (SE = %.4f) %d%% CI: [%.4f, %.4f]\n",
                aipe, seA, round(level_est * 100),
                aipe - zc_est * seA, aipe + zc_est * seA))
  }
  if (do_dcame_compare && estimand == "AIPE") {
    cat(sprintf("  D-CAME = %.4f (SE = %.4f) %d%% CI: [%.4f, %.4f]\n",
                dcame, se_dcame, round(level_est * 100),
                dcame - zc_est * se_dcame, dcame + zc_est * se_dcame))
    cat(sprintf("  CAME(X=0) = %.4f (SE = %.4f) | CAME(X=1) = %.4f (SE = %.4f)\n",
                came0, se_came0, came1, se_came1))
  }

  if (do_linear) {
    cat("--- Linear comparison (joint linear-per-group; D-CAME = slope_high - slope_low) ---\n")
    # Use ALL observations. `x` here is the group index g in {0,...,K-1}
    # (no NA: middle groups are kept; same binning the basis model uses).
    # Model fitted (matches Y = sum_j {mu_j + beta_j*D + alpha_j*X} * G_j + delta'Z + FE + e):
    #   Y = mu_0 + beta_0 * D + sum_{k>=1} (beta_k - beta_0) * D * I(g=k)
    #         + sum_{k>=1} eta_k * I(g=k)
    #         + [continuous X only:] alpha_0 * X + sum_{k>=1}(alpha_k - alpha_0) * X * I(g=k)
    #         + delta'Z + FE + e
    # Contrast reported = beta_{K-1} - beta_0 = coef on D:g_fac{K-1}.
    # This makes the Linear bar IDENTICAL to the basis-model D-CAME when
    # user_spec_0 = ... = user_spec_{K-1} = 1 (linear in every group).
    # For binary X, X equals the group index, so X*g_fac is collinear and is omitted.
    # Build the K-group index LOCALLY so the Linear comparison is identical
    # for model = "basis" and model = "GAM" (the GAM path passes x = {0,1,NA},
    # which would otherwise drop middle observations and break the K-1 contrast).
    if (x_is_binary) {
      g_lin <- as.integer(x)            # binary X: x already in {0,1}
    } else if (!is.null(x_raw_full)) {
      cut_probs <- seq(0, 1, length.out = K + 1L)[-c(1L, K + 1L)]
      cuts_lin  <- as.numeric(stats::quantile(x_raw_full, probs = cut_probs,
                                              na.rm = TRUE))
      g_lin <- rep(NA_integer_, length(x_raw_full))
      g_lin[x_raw_full <= cuts_lin[1]] <- 0L
      if (K > 2L) {
        for (kk in 2:(K - 1L)) {
          g_lin[x_raw_full > cuts_lin[kk - 1L] &
                  x_raw_full <= cuts_lin[kk]] <- as.integer(kk - 1L)
        }
      }
      g_lin[x_raw_full > cuts_lin[length(cuts_lin)]] <- as.integer(K - 1L)
    } else {
      # Continuous X with no raw vector available: fall back to whatever was passed.
      g_lin <- as.integer(x)
    }
    # The inputs y, d_all, x_raw_full, z_mat, etc. arriving here are ALREADY
    # the common-support-trimmed samples (trimmed in the dispatcher before
    # .render_results() was called). So the Linear comparison uses exactly
    # the same observations as the basis/GAM fits. Middle moderator groups
    # are kept for all K-group bins.
    keep_lin <- !is.na(g_lin)
    dl_y <- y[keep_lin]
    dl_d <- d_all[keep_lin]
    dl_g <- factor(g_lin[keep_lin], levels = 0:(K - 1L))
    # Strip the X_continuous control from z_mat (we add it back, interacted with g_fac).
    z_for_lin <- if (!is.null(z_mat)) {
      drop_col <- which(colnames(z_mat) == "X_continuous")
      if (length(drop_col) > 0L) z_mat[, -drop_col, drop = FALSE] else z_mat
    } else NULL
    dl_z <- if (!is.null(z_for_lin)) z_for_lin[keep_lin, , drop = FALSE] else NULL
    dl_w <- if (!is.null(survey_wgts)) survey_wgts[keep_lin] else NULL
    dl_c <- if (!is.null(cluster_var)) cluster_var[keep_lin] else NULL
    dl_fe <- if (!is.null(fe_mat)) fe_mat[keep_lin, , drop = FALSE] else NULL
    dl_fefac <- if (isTRUE(Z_FE_absorb) && !is.null(fe_absorb_factors))
      lapply(fe_absorb_factors, function(f) factor(as.character(f[keep_lin])))
    else NULL

    dl <- data.frame(Y_out = dl_y, D_treat = dl_d, g_fac = dl_g)
    # Continuous X only: add raw X (will be interacted with g_fac in the formula).
    has_xraw <- !is.null(x_raw_full) && !x_is_binary
    if (has_xraw) dl$X_raw <- x_raw_full[keep_lin]
    if (!is.null(dl_z) && ncol(dl_z) > 0) dl <- cbind(dl, dl_z)
    if (!is.null(dl_fe) && ncol(dl_fe) > 0) dl <- cbind(dl, dl_fe)

    contrast_name <- paste0("D_treat:g_fac", K - 1L)

    # RHS pieces: D*g_fac (group-specific slopes & intercepts) + optional X*g_fac
    # (group-specific linear X effects) + other Z controls + FE.
    rhs_core <- if (has_xraw) "D_treat * g_fac + X_raw * g_fac" else "D_treat * g_fac"

    if (!is.null(dl_fefac)) {
      for (v in names(dl_fefac)) dl[[v]] <- dl_fefac[[v]]
      reserved <- c("Y_out", "D_treat", "g_fac",
                    if (has_xraw) "X_raw" else character(0),
                    names(dl_fefac))
      other_cols <- setdiff(colnames(dl), reserved)
      base_rhs <- rhs_core
      if (length(other_cols) > 0L)
        base_rhs <- paste(base_rhs, "+",
                          paste0("`", other_cols, "`", collapse = " + "))
      fe_rhs <- paste(paste0("`", names(dl_fefac), "`"), collapse = " + ")
      fml <- stats::as.formula(sprintf("Y_out ~ %s | %s", base_rhs, fe_rhs))
      if (vce == "cluster" && !is.null(dl_c)) {
        dl$.CL_VAR <- dl_c
        vcov_arg <- stats::as.formula("~ .CL_VAR")
      } else {
        vcov_arg <- "hetero"
      }
      fl <- fixest::feols(fml, data = dl, weights = dl_w,
                          vcov = vcov_arg, warn = FALSE, notes = FALSE)
      Vl <- as.matrix(stats::vcov(fl))
      coefs_fl <- stats::coef(fl)
    } else {
      reserved <- c("Y_out", "D_treat", "g_fac",
                    if (has_xraw) "X_raw" else character(0))
      other_cols <- setdiff(colnames(dl), reserved)
      base_rhs <- rhs_core
      if (length(other_cols) > 0L)
        base_rhs <- paste(base_rhs, "+",
                          paste0("`", other_cols, "`", collapse = " + "))
      fml <- stats::as.formula(paste("Y_out ~", base_rhs))
      fl <- if (!is.null(dl_w))
        lm(fml, data = dl, weights = dl_w)
      else
        lm(fml, data = dl)
      Vl <- if (vce == "cluster" && !is.null(dl_c))
        sandwich::vcovCL(fl, cluster = dl_c)
      else sandwich::vcovHC(fl, type = "HC1")
      coefs_fl <- coef(fl)
    }

    if (!(contrast_name %in% names(coefs_fl)) ||
        !(contrast_name %in% rownames(Vl))) {
      warning(sprintf(
        "Linear comparison: the contrast coefficient '%s' was dropped from the fit, ",
        contrast_name),
        "usually because of collinearity in the design (rank-deficient X). ",
        "Common causes: (1) the highest moderator group is empty or near-empty ",
        "after common-support trimming, (2) one of the FE variables has a level ",
        "with all observations in a single moderator group, ",
        "(3) a column in Z is collinear with the group dummies. ",
        "Try inspecting the group sizes (n per K-group), removing the offending ",
        "FE / Z column, or reducing n_parts. Reporting Linear = NA.",
        call. = FALSE)
      aiL <- NA_real_; seL <- NA_real_
    } else {
      aiL <- as.numeric(coefs_fl[contrast_name])
      seL <- sqrt(as.numeric(Vl[contrast_name, contrast_name]))
    }
    cat(sprintf("  Linear D-CAME (slope_high - slope_low) = %.4f (SE = %.4f)\n",
                aiL, seL))
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PLOTS
  # ══════════════════════════════════════════════════════════════════════════
  cat("--- Building plot objects ---\n")

  if (x_is_binary) {
    if (!is.null(x_original_labels)) {
      grp0_label <- x_original_labels[1]
      grp1_label <- x_original_labels[2]
    } else {
      grp0_label <- "X = 0"
      grp1_label <- "X = 1"
    }
  } else if (model == "GAM") {
    # GAM passes the representative low/high moderator values (x_lo, x_hi).
    grp0_label <- sprintf("X = %.3g (low)",  tercile_cuts[1])
    grp1_label <- sprintf("X = %.3g (high)", tercile_cuts[length(tercile_cuts)])
  } else {
    # Basis path passes the K-1 internal cut-points. The low group is the
    # bottom bin (X <= first cut); the high group is the top bin (X > last cut).
    grp0_label <- sprintf("X <= %.3g (low)", tercile_cuts[1])
    grp1_label <- sprintf("X > %.3g (high)", tercile_cuts[length(tercile_cuts)])
  }

  dfp <- data.frame(
    D = rep(dg, 2),
    Yhat = c(yh0, yh1),
    se = c(se0, se1),
    Group = rep(c(grp0_label, grp1_label), each = n_grid))
  dfp$lo <- dfp$Yhat - zc_pv * dfp$se
  dfp$up <- dfp$Yhat + zc_pv * dfp$se

  p1 <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = dfp,
                         ggplot2::aes(x = D, ymin = lo, ymax = up, fill = Group), alpha = .15) +
    ggplot2::geom_line(data = dfp,
                       ggplot2::aes(x = D, y = Yhat, color = Group, linetype = Group), linewidth = 1)

  if (hist_pv == "perc") {
    pctl_probs <- c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)
    pctl_labels <- c("p1", "p10", "p25", "p50", "p75", "p90", "p99")

    pctl_X0 <- quantile(d_obs_X0, probs = pctl_probs, na.rm = TRUE)
    pctl_X1 <- quantile(d_obs_X1, probs = pctl_probs, na.rm = TRUE)

    pctl_df <- data.frame(
      D = c(as.numeric(pctl_X0), as.numeric(pctl_X1)),
      Group = c(rep(grp0_label, length(pctl_probs)),
                rep(grp1_label, length(pctl_probs))),
      label = rep(pctl_labels, 2),
      stringsAsFactors = FALSE
    )

    y_range <- range(c(dfp$lo, dfp$up), na.rm = TRUE)
    y_span <- diff(y_range)
    tick_len <- y_span * 0.03

    pctl_ticks <- pctl_df
    pctl_ticks$y_start <- ifelse(pctl_ticks$Group == grp1_label,
                                 y_range[2], y_range[1])
    pctl_ticks$y_end <- ifelse(pctl_ticks$Group == grp1_label,
                               y_range[2] - tick_len, y_range[1] + tick_len)
    pctl_ticks$y_label <- ifelse(pctl_ticks$Group == grp1_label,
                                 y_range[2] + tick_len * 0.3,
                                 y_range[1] - tick_len * 0.3)

    p1 <- p1 +
      ggplot2::geom_segment(data = pctl_ticks,
                            ggplot2::aes(x = D, xend = D, y = y_start, yend = y_end),
                            color = ifelse(pctl_ticks$Group == grp1_label, "black", "gray50"),
                            linewidth = 0.4, alpha = 0.7) +
      ggplot2::geom_text(data = pctl_ticks,
                         ggplot2::aes(x = D, y = y_label, label = label),
                         size = 1.8,
                         color = ifelse(pctl_ticks$Group == grp1_label, "black", "gray50"),
                         vjust = ifelse(pctl_ticks$Group == grp1_label, 0, 1))

    p1_caption <- paste0("Tick marks: top = ", grp1_label,
                         " | bottom = ", grp0_label,
                         " (p1, p10, p25, p50, p75, p90, p99)")

  } else {
    if (hist_pv == "rug_1000") {
      .old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
        get(".Random.seed", envir = .GlobalEnv) else NULL
      set.seed(999)
      max_rug <- 1000
      d_rug_X0 <- if (length(d_obs_X0) > max_rug) sample(d_obs_X0, max_rug) else d_obs_X0
      d_rug_X1 <- if (length(d_obs_X1) > max_rug) sample(d_obs_X1, max_rug) else d_obs_X1
      if (!is.null(.old_seed))
        assign(".Random.seed", .old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      d_rug_X0 <- d_obs_X0
      d_rug_X1 <- d_obs_X1
    }

    .rug_alpha <- function(n_grp) {
      if (n_grp <= 750) return(0.3)
      max(0.02, 0.3 * (750 / n_grp)^0.5)
    }
    alpha_X0 <- .rug_alpha(length(d_rug_X0))
    alpha_X1 <- .rug_alpha(length(d_rug_X1))

    rug_df_X1 <- data.frame(D = d_rug_X1, Group = grp1_label)
    rug_df_X0 <- data.frame(D = d_rug_X0, Group = grp0_label)

    p1 <- p1 +
      ggplot2::geom_rug(data = rug_df_X1, ggplot2::aes(x = D), sides = "t",
                        color = "black", alpha = alpha_X1,
                        length = ggplot2::unit(0.03, "npc")) +
      ggplot2::geom_rug(data = rug_df_X0, ggplot2::aes(x = D), sides = "b",
                        color = "gray50", alpha = alpha_X0,
                        length = ggplot2::unit(0.03, "npc"))

    rug_suffix <- if (hist_pv == "rug_1000") " (max 1000 per group)" else ""
    p1_caption <- paste0("Rug: top = ", grp1_label, " | bottom = ", grp0_label, rug_suffix)
  }

  p1 <- p1 +
    ggplot2::scale_fill_manual(values = stats::setNames(c("gray50", "black"),
                                                 c(grp0_label, grp1_label))) +
    ggplot2::scale_color_manual(values = stats::setNames(c("gray50", "black"),
                                                  c(grp0_label, grp1_label))) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "dashed"),
                                                     c(grp0_label, grp1_label))) +
    ggplot2::labs(x = "D", y = "Fitted values",
                  title = "Predicted Values",
                  subtitle = sprintf("(%d%% CI)", round(level_pv * 100)),
                  caption = p1_caption) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   legend.title = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = .5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = .5, size = 9, color = "gray40"),
                   plot.caption = ggplot2::element_text(size = 8, color = "gray50"),
                   plot.margin = ggplot2::margin(t = 15, r = 5, b = 15, l = 5))

  p2 <- NULL
  dip <- NULL
  if (estimand == "AIPE") {
    pcts <- c(.01, .10, .25, .50, .75, .90, .99)
    plb  <- c("p1", "p10", "p25", "p50", "p75", "p90", "p99")
    dpv  <- quantile(d_all, probs = pcts, na.rm = TRUE)
    ipv  <- approx(dg, ipeg, xout = dpv, rule = 2)$y
    spv  <- approx(dg, seig, xout = dpv, rule = 2)$y

    dip <- data.frame(D = as.numeric(dpv), label = plb,
                      IPE = ipv, lo = ipv - zc_est * spv, up = ipv + zc_est * spv)

    dfipe <- data.frame(D = dg, IPE = ipeg,
                        lo = ipeg - zc_est * seig,
                        up = ipeg + zc_est * seig)

    p2 <- ggplot2::ggplot(dfipe, ggplot2::aes(x = D)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = up), alpha = 0.15, fill = "black") +
      ggplot2::geom_line(ggplot2::aes(y = IPE), linewidth = 1.3, color = "black") +
      ggplot2::geom_point(data = dip, ggplot2::aes(x = D, y = IPE),
                          shape = 16, size = 2.5, color = "gray30") +
      ggplot2::scale_x_continuous(
        sec.axis = ggplot2::dup_axis(name = NULL,
                                     breaks = as.numeric(dpv), labels = plb)) +
      ggplot2::labs(x = "D", y = "Interactive Effects",
                    title = "Interactive Partial Effects",
                    subtitle = sprintf("(%d%% CI)", round(level_est * 100))) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = .5, size = 9, color = "gray40"),
                     axis.text.x.top = ggplot2::element_text(size = 7, color = "gray40"),
                     axis.ticks.x.top = ggplot2::element_line(color = "gray40"))
  }

  has_any_compare <- !isFALSE(compare)

  if (estimand == "DCAME" && !has_any_compare) {
    da <- data.frame(
      Specification = factor(c("CAME(X=0)", "CAME(X=1)", "D-CAME"),
                             levels = c("CAME(X=0)", "CAME(X=1)", "D-CAME")),
      Estimate = c(came0, came1, dcame),
      lo = c(came0 - zc_est * se_came0, came1 - zc_est * se_came1, dcame - zc_est * se_dcame),
      up = c(came0 + zc_est * se_came0, came1 + zc_est * se_came1, dcame + zc_est * se_dcame)
    )
  } else if (estimand == "DCAME" && has_any_compare) {
    da <- data.frame(Specification = "D-CAME", Estimate = dcame,
                     lo = dcame - zc_est * se_dcame, up = dcame + zc_est * se_dcame)
    if (do_aipe_compare) {
      da <- rbind(da,
                  data.frame(Specification = "AIPE", Estimate = aipe,
                             lo = aipe - zc_est * seA, up = aipe + zc_est * seA))
    }
    if (do_linear && !is.null(aiL)) {
      da <- rbind(data.frame(Specification = "Linear", Estimate = aiL,
                             lo = aiL - zc_est * seL, up = aiL + zc_est * seL), da)
    }
    spec_levels <- c("Linear", "D-CAME", "AIPE")
    da$Specification <- factor(da$Specification,
                               levels = intersect(spec_levels, da$Specification))
  } else if (estimand == "AIPE") {
    da <- data.frame(Specification = "AIPE", Estimate = aipe,
                     lo = aipe - zc_est * seA, up = aipe + zc_est * seA)
    if (do_dcame_compare) {
      da <- rbind(data.frame(Specification = "D-CAME", Estimate = dcame,
                             lo = dcame - zc_est * se_dcame, up = dcame + zc_est * se_dcame),
                  da)
    }
    if (do_linear && !is.null(aiL)) {
      da <- rbind(data.frame(Specification = "Linear", Estimate = aiL,
                             lo = aiL - zc_est * seL, up = aiL + zc_est * seL), da)
    }
    spec_levels <- c("Linear", "D-CAME", "AIPE")
    da$Specification <- factor(da$Specification,
                               levels = intersect(spec_levels, da$Specification))
  }

  p3 <- ggplot2::ggplot(da, ggplot2::aes(x = Specification, y = Estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lo, ymax = up),
                             size = 1.2, linewidth = 1.5) +
    ggplot2::labs(x = "Estimand", y = "Marginal Effects",
                  title = "Interactive Effects",
                  subtitle = sprintf("(%d%% CI)", round(level_est * 100)),
                  caption = sprintf("model: %s", model)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = .5, size = 9, color = "gray40"),
                   plot.caption = ggplot2::element_text(hjust = .5, size = 8, color = "gray50"))

  if (estimand == "DCAME") {
    cp <- gridExtra::arrangeGrob(p1, p3, nrow = 1, widths = c(1.2, .8))
  } else {
    cp <- gridExtra::arrangeGrob(p1, p2, p3, nrow = 1, widths = c(1.2, 1, .8))
  }
  class(cp) <- c("dcame_aipe_grob", class(cp))

  # ── Optional plot: CAME at every level of X (plot_came_all) ──────────
  p_came_all <- NULL
  if (isTRUE(plot_came_all) && !is.null(came_all_df) && nrow(came_all_df) > 0L) {
    da_all <- came_all_df
    da_all$group_factor <- factor(da_all$group_label, levels = da_all$group_label)
    p_came_all <- ggplot2::ggplot(da_all,
                                  ggplot2::aes(x = group_factor, y = CAME)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = CI_lo, ymax = CI_hi),
                               size = 1.0, linewidth = 1.3) +
      ggplot2::labs(x = "Moderator level", y = "CAME",
                    title = sprintf("CAME at each of %d moderator levels", K),
                    subtitle = sprintf("(%d%% CI)", round(level_est * 100)),
                    caption = "Each CAME is the average dY/dD within that group, from the joint K-group model.") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title    = ggplot2::element_text(hjust = .5, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = .5, size = 9, color = "gray40"),
                     plot.caption  = ggplot2::element_text(hjust = .5, size = 8, color = "gray50"),
                     axis.text.x   = ggplot2::element_text(angle = 25, hjust = 1))
  }

  res <- list(
    model            = model,
    estimand         = estimand,
    inference_method = inference,
    selection_method = selection,
    Z_interact       = Z_interact,
    FE_method        = FE,
    FE_vbles         = FE_vbles,
    Z_FE             = Z_FE,
    n_fe_dummies     = n_fe_dummies,
    weighted         = !is.null(wgt),
    moderator_type   = if (x_is_binary) "binary" else "continuous",
    tercile_cuts     = tercile_cuts,
    n_dropped_middle = n_dropped_t2,
    spec_X0 = bX0, spec_X1 = bX1,
    spec_X0_label = bX0_label, spec_X1_label = bX1_label,
    fit = fit_obj, vcov = V_obj,

    level_pv = level_pv, level_est = level_est,

    d_grid = dg,
    yhat_X0 = yh0, yhat_X1 = yh1,
    se_yhat_X0 = se0, se_yhat_X1 = se1,
    slope_grid_X0 = slope0g, slope_grid_X1 = slope1g,
    se_slope_X0 = se_slope0g, se_slope_X1 = se_slope1g,
    ipe_grid = ipeg, se_ipe_grid = seig,
    ipe_at_percentiles = dip,

    aipe = aipe, se_aipe = seA,
    ci_aipe = c(aipe - zc_est * seA, aipe + zc_est * seA),

    dcame = dcame, se_dcame = se_dcame,
    ci_dcame = c(dcame - zc_est * se_dcame, dcame + zc_est * se_dcame),
    came0 = came0, came1 = came1,
    se_came0 = se_came0, se_came1 = se_came1,

    linear_dx = aiL, se_linear_dx = seL,

    common_support = c(cs_lo, cs_hi),
    n_trimmed = nt, n_used = n,

    K = K,
    came_all = came_all_df,
    Z_FE_absorb = isTRUE(Z_FE_absorb),
    compare = compare,

    plot = cp, plot_predicted = p1, plot_ipe = p2, plot_coef = p3,
    plot_came_all = p_came_all
  )
  class(res) <- "dcame_aipe"
  invisible(res)
}


# ══════════════════════════════════════════════════════════════════════════════
# S3 print methods
# ══════════════════════════════════════════════════════════════════════════════

#' Print the combined plot grob, redrawing the figure
#' @param x a `dcame_aipe_grob` object
#' @param ... unused
#' @export
print.dcame_aipe_grob <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}

#' Print method for `dcame_aipe` results
#' @param x a `dcame_aipe` object returned by [dcame_aipe()].
#' @param ... unused
#' @export
print.dcame_aipe <- function(x, ...) {
  cat("\n=== DCAME/AIPE Results ===\n")
  model_str <- if (!is.null(x$model)) x$model else "basis"
  cat(sprintf("Model: %s | Estimand: %s | Inference: %s | Selection: %s\n",
              model_str,
              x$estimand,
              x$inference_method,
              x$selection_method))
  cat(sprintf("Moderator: %s\n", x$moderator_type))
  if (x$moderator_type == "continuous") {
    if (model_str == "GAM") {
      cat(sprintf("  Contrast points: %.3f (low), %.3f (high)\n",
                  x$tercile_cuts[1], x$tercile_cuts[2]))
    } else {
      cat(sprintf("  Cut-points: %s | n_parts = %d (no obs dropped)\n",
                  paste(sprintf("%.3f", x$tercile_cuts), collapse = ", "),
                  if (!is.null(x$K)) x$K else (length(x$tercile_cuts) + 1L)))
    }
  }
  cat(sprintf("Z_interact: %s", x$Z_interact))
  if (!is.null(x$FE_method)) {
    cat(sprintf(" | FE: %s (%s)", x$FE_method,
                paste(x$FE_vbles, collapse = ", ")))
  }
  if (!is.null(x$Z_FE) && length(x$Z_FE) > 0)
    cat(sprintf(" | Z_FE: %s (%d dummies)", paste(x$Z_FE, collapse = ","),
                x$n_fe_dummies))
  if (x$weighted) cat(" | weighted")
  cat("\n")
  cat(sprintf("X=0: %s | X=1: %s\n", x$spec_X0_label, x$spec_X1_label))
  cat(sprintf("N = %d", x$n_used))
  if (x$n_trimmed > 0) cat(sprintf(" (%d trimmed)", x$n_trimmed))
  cat(sprintf("\nSupport: [%.3f, %.3f]\n", x$common_support[1], x$common_support[2]))
  cat(sprintf("CI levels: predicted values = %d%% | estimands = %d%%\n",
              round(x$level_pv * 100), round(x$level_est * 100)))

  # Respect the original `compare` setting so the same estimand is not shown
  # twice (once during call, once at print time). With compare = FALSE,
  # only the headline estimand is printed.
  cmp <- x$compare
  show_dcame  <- x$estimand == "DCAME" ||
    (!isFALSE(cmp) && is.character(cmp) && cmp %in% c("DCAME", "ALL"))
  show_aipe   <- x$estimand == "AIPE" ||
    (!isFALSE(cmp) && is.character(cmp) && cmp %in% c("AIPE", "ALL"))
  show_linear <- !isFALSE(cmp) && is.character(cmp) && cmp %in% c("LINEAR", "ALL") &&
    !is.null(x$linear_dx)

  if (show_dcame) {
    cat(sprintf("D-CAME = %.4f (SE = %.4f) CI: [%.4f, %.4f]\n",
                x$dcame, x$se_dcame, x$ci_dcame[1], x$ci_dcame[2]))
    cat(sprintf("  CAME(X=0) = %.4f (SE = %.4f) | CAME(X=1) = %.4f (SE = %.4f)\n",
                x$came0, x$se_came0, x$came1, x$se_came1))
  }
  if (show_aipe) {
    cat(sprintf("AIPE   = %.4f (SE = %.4f) CI: [%.4f, %.4f]\n",
                x$aipe, x$se_aipe, x$ci_aipe[1], x$ci_aipe[2]))
  }
  if (show_linear) {
    cat(sprintf("Linear = %.4f (SE = %.4f)\n", x$linear_dx, x$se_linear_dx))
  }

  if (!is.null(x$ipe_at_percentiles) && "label" %in% names(x$ipe_at_percentiles)) {
    cat("\nIPE at percentiles:\n")
    print(x$ipe_at_percentiles[, c("label", "D", "IPE", "lo", "up")],
          row.names = FALSE, digits = 4)
  }
  if (!is.null(x$came_all) && nrow(x$came_all) > 0L) {
    cat("\nCAME at each moderator level (plot_came_all):\n")
    print(x$came_all[, c("group", "group_label", "n", "CAME", "SE",
                         "CI_lo", "CI_hi", "p_value")],
          row.names = FALSE, digits = 4)
  }
  cat("\n")
  invisible(x)
}


# ══════════════════════════════════════════════════════════════════════════════
# Summary method
# ══════════════════════════════════════════════════════════════════════════════

#' Summary table for a `dcame_aipe` result
#'
#' Returns a clean data.frame with AIPE, D-CAME, CAME(X=0), CAME(X=1)
#' (and Linear if requested via `compare = "ALL"`/`"LINEAR"`), each with
#' SE, confidence interval, and two-sided Wald p-value (H0: effect = 0).
#'
#' @param object a `dcame_aipe` object returned by [dcame_aipe()].
#' @param ... unused
#' @return Invisibly returns the summary `data.frame`.
#' @export
summary.dcame_aipe <- function(object, ...) {
  x <- object
  zc_est <- qnorm(1 - (1 - x$level_est) / 2)
  zc_pv  <- qnorm(1 - (1 - x$level_pv) / 2)

  rows <- list()

  p_aipe <- 2 * pnorm(-abs(x$aipe / x$se_aipe))
  rows[["AIPE"]] <- data.frame(
    Estimand  = "AIPE",
    Estimate  = x$aipe,
    SE        = x$se_aipe,
    CI_lo     = x$aipe - zc_est * x$se_aipe,
    CI_hi     = x$aipe + zc_est * x$se_aipe,
    p_value   = p_aipe,
    stringsAsFactors = FALSE
  )

  p_dcame <- 2 * pnorm(-abs(x$dcame / x$se_dcame))
  rows[["D-CAME"]] <- data.frame(
    Estimand  = "D-CAME",
    Estimate  = x$dcame,
    SE        = x$se_dcame,
    CI_lo     = x$dcame - zc_est * x$se_dcame,
    CI_hi     = x$dcame + zc_est * x$se_dcame,
    p_value   = p_dcame,
    stringsAsFactors = FALSE
  )

  p_came0 <- 2 * pnorm(-abs(x$came0 / x$se_came0))
  rows[["CAME(X=0)"]] <- data.frame(
    Estimand  = "CAME(X=0)",
    Estimate  = x$came0,
    SE        = x$se_came0,
    CI_lo     = x$came0 - zc_est * x$se_came0,
    CI_hi     = x$came0 + zc_est * x$se_came0,
    p_value   = p_came0,
    stringsAsFactors = FALSE
  )

  p_came1 <- 2 * pnorm(-abs(x$came1 / x$se_came1))
  rows[["CAME(X=1)"]] <- data.frame(
    Estimand  = "CAME(X=1)",
    Estimate  = x$came1,
    SE        = x$se_came1,
    CI_lo     = x$came1 - zc_est * x$se_came1,
    CI_hi     = x$came1 + zc_est * x$se_came1,
    p_value   = p_came1,
    stringsAsFactors = FALSE
  )

  if (!is.null(x$linear_dx)) {
    p_lin <- 2 * pnorm(-abs(x$linear_dx / x$se_linear_dx))
    rows[["Linear"]] <- data.frame(
      Estimand  = "Linear",
      Estimate  = x$linear_dx,
      SE        = x$se_linear_dx,
      CI_lo     = x$linear_dx - zc_est * x$se_linear_dx,
      CI_hi     = x$linear_dx + zc_est * x$se_linear_dx,
      p_value   = p_lin,
      stringsAsFactors = FALSE
    )
  }

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  ci_label <- sprintf("%d%%", round(x$level_est * 100))
  names(tab)[4:5] <- paste0("CI_", ci_label, c("_lo", "_hi"))

  cat("\n=== Summary: DCAME/AIPE Estimates ===\n")
  cat(sprintf("Model: %s | Inference: %s | N = %d\n",
              if (!is.null(x$model)) x$model else "basis",
              x$inference_method, x$n_used))
  cat(sprintf("CI level (estimands): %s | CI level (predicted values): %d%%\n",
              ci_label, round(x$level_pv * 100)))
  cat(sprintf("Support: [%.3f, %.3f]\n\n", x$common_support[1], x$common_support[2]))

  print(tab, digits = 4, row.names = FALSE)

  if (!is.null(x$came_all) && nrow(x$came_all) > 0L) {
    cat("\nCAME at each moderator level:\n")
    print(x$came_all[, c("group", "group_label", "n", "CAME", "SE",
                         "CI_lo", "CI_hi", "p_value")],
          row.names = FALSE, digits = 4)
  }
  cat("\n")

  invisible(tab)
}


# ══════════════════════════════════════════════════════════════════════════════
# Export-ready table of estimates
# ══════════════════════════════════════════════════════════════════════════════

#' Export estimates table to data.frame, .xlsx, .tex, or .csv
#'
#' Produces a clean `data.frame` with all estimands (AIPE, D-CAME,
#' CAME(X=0), CAME(X=1), and Linear if requested via `compare = "ALL"`/
#' `"LINEAR"`) that can be written to spreadsheet, LaTeX, or CSV.
#'
#' @param object a `dcame_aipe` object returned by [dcame_aipe()].
#' @param file Optional path. Extension determines format
#'   (`.xlsx`/`.xls` requires `openxlsx`; `.tex`; `.csv`).
#' @param digits Number of decimals in LaTeX/xlsx output (default 4).
#' @param caption Optional LaTeX caption.
#' @param label Optional LaTeX label.
#' @return The estimates `data.frame` (invisibly if `file` is set).
#' @examples
#' \dontrun{
#' tab <- dcame_aipe_table(res)                       # data.frame
#' dcame_aipe_table(res, file = "out.xlsx")           # Excel
#' dcame_aipe_table(res, file = "out.tex")            # LaTeX
#' }
#' @export
dcame_aipe_table <- function(object, file = NULL, digits = 4,
                             caption = NULL, label = NULL) {
  if (!inherits(object, "dcame_aipe"))
    stop("`object` must be a result returned by dcame_aipe().")

  zc_est <- qnorm(1 - (1 - object$level_est) / 2)

  .row <- function(name, est, se) {
    if (is.null(est) || is.null(se) || is.na(est) || is.na(se))
      return(NULL)
    data.frame(
      Estimand = name,
      Estimate = as.numeric(est),
      SE       = as.numeric(se),
      CI_lo    = as.numeric(est - zc_est * se),
      CI_hi    = as.numeric(est + zc_est * se),
      p_value  = as.numeric(2 * pnorm(-abs(est / se))),
      stringsAsFactors = FALSE
    )
  }

  rows <- list(
    .row("AIPE",      object$aipe,      object$se_aipe),
    .row("D-CAME",    object$dcame,     object$se_dcame),
    .row("CAME(X=0)", object$came0,     object$se_came0),
    .row("CAME(X=1)", object$came1,     object$se_came1),
    .row("Linear",    object$linear_dx, object$se_linear_dx)
  )
  if (!is.null(object$came_all) && nrow(object$came_all) > 0L) {
    ca <- object$came_all
    for (i in seq_len(nrow(ca))) {
      rows[[length(rows) + 1L]] <- .row(
        sprintf("CAME[%s]", ca$group_label[i]),
        ca$CAME[i], ca$SE[i])
    }
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL

  ci_lab <- sprintf("CI_%d%%", round(object$level_est * 100))
  names(tab)[4:5] <- paste0(ci_lab, c("_lo", "_hi"))

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("openxlsx", quietly = TRUE))
        stop("Package 'openxlsx' is required to write .xlsx files. ",
             "Install with install.packages('openxlsx').")
      openxlsx::write.xlsx(tab, file = file, rowNames = FALSE)
    } else if (ext == "tex") {
      .fmt <- function(v) formatC(v, format = "f", digits = digits)
      body <- apply(tab, 1, function(r) {
        paste0(r[[1]], " & ",
               .fmt(as.numeric(r[[2]])), " & ",
               .fmt(as.numeric(r[[3]])), " & [",
               .fmt(as.numeric(r[[4]])), ", ",
               .fmt(as.numeric(r[[5]])), "] & ",
               .fmt(as.numeric(r[[6]])), " \\\\")
      })
      ci_hdr <- sprintf("%d\\%% CI", round(object$level_est * 100))
      tex <- c(
        "\\begin{table}[ht]",
        "\\centering",
        if (!is.null(caption)) sprintf("\\caption{%s}", caption),
        if (!is.null(label))   sprintf("\\label{%s}",   label),
        "\\begin{tabular}{lrrcrr}",
        "\\hline",
        paste("Estimand & Estimate & SE &", ci_hdr,
              "& $p$-value \\\\"),
        "\\hline",
        body,
        "\\hline",
        "\\end{tabular}",
        "\\end{table}"
      )
      writeLines(tex, con = file)
    } else if (ext == "csv") {
      utils::write.csv(tab, file = file, row.names = FALSE)
    } else {
      stop(sprintf("Unsupported file extension '%s'. Use .xlsx, .tex, or .csv.", ext))
    }
    invisible(tab)
  } else {
    tab
  }
}
