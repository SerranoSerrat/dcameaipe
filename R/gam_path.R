# ============================================================
# gam_path.R
#
# GAM-based estimation path for dcame_aipe().
#
# Contains the GAM fit/predict engine, two inference methods
# (analytic posterior + case/cluster bootstrap), and the GAM
# dispatcher .dcame_aipe_gam() which is called from the main
# dcame_aipe() function when model = "GAM".
#
# All functions are internal.
# ============================================================

#' Fit GAM and return all estimands of interest on a grid
#' @keywords internal
#' @noRd
.gam_fit_and_predict <- function(y, d, x_raw, z_mat, d_grid, x_lo, x_hi,
                                 survey_w = NULL, fe_data = NULL,
                                 fe_vars = NULL,
                                 d_all_for_aipe = NULL,
                                 survey_w_all_for_aipe = NULL,
                                 smooth = "REML",
                                 k_s = NULL, k_t = NULL,
                                 inter_gam = "ti",
                                 compute_se = TRUE,
                                 z_interact = FALSE,
                                 g_loc = NULL,    # K-group index (continuous X)
                                 n_parts = 3L) {  # number of bins (continuous X)
  gam_method <- if (smooth == "GCV") "GCV.Cp" else "REML"
  # Assemble data frame
  df_fit <- data.frame(Y_out = y, D_treat = d, X_mod = x_raw)
  if (!is.null(z_mat) && ncol(z_mat) > 0) {
    df_fit <- cbind(df_fit, z_mat)
    z_names <- colnames(z_mat)
  } else {
    z_names <- character(0)
  }

  # Attach FE variables (as-is, will be wrapped with factor() in formula)
  fe_names <- character(0)
  if (!is.null(fe_vars) && length(fe_vars) > 0 && !is.null(fe_data)) {
    for (v in fe_vars) {
      df_fit[[v]] <- factor(fe_data[[v]])
    }
    fe_names <- fe_vars
  }

  x_is_bin <- all(x_raw %in% c(0, 1))
  K_g <- if (x_is_bin) 2L else as.integer(n_parts)

  # Continuous X: attach the K-group factor for use in Z*g_fac interactions
  # (4G: controls interact with moderator GROUPS, not raw X).
  if (!x_is_bin && !is.null(g_loc)) {
    df_fit$g_fac <- factor(g_loc, levels = 0:(K_g - 1L))
  }

  # Build k = ... fragments. When NULL we omit `k` entirely so mgcv uses its
  # internal defaults (10 for s(), 5 per marginal for te()/ti()).
  k_arg_s <- if (!is.null(k_s)) sprintf(", k = %d", k_s) else ""
  k_arg_t <- if (!is.null(k_t)) sprintf(", k = c(%d, %d)", k_t[1], k_t[2]) else ""

  if (x_is_bin) {
    df_fit$x_factor <- factor(x_raw, levels = c(0, 1))
    rhs <- sprintf("x_factor + s(D_treat, by = x_factor%s)", k_arg_s)
  } else if (inter_gam == "te") {
    rhs <- sprintf("te(D_treat, X_mod%s)", k_arg_t)
  } else {
    rhs <- sprintf(
      "s(D_treat%s) + s(X_mod%s) + ti(D_treat, X_mod%s)",
      k_arg_s, k_arg_s, k_arg_t)
  }

  if (length(z_names) > 0) {
    if (z_interact && x_is_bin) {
      z_rhs <- paste(sprintf("`%s` * x_factor", z_names), collapse = " + ")
    } else if (z_interact && !x_is_bin) {
      # 4G: For continuous X, interact controls with the K-group factor g_fac
      # (matches the basis path's Z_g{k} blocks), not with raw X_mod.
      if (is.null(df_fit$g_fac)) {
        warning("Z_interact = TRUE with continuous X but no group factor available; ",
                "falling back to additive Z without interaction.", call. = FALSE)
        z_rhs <- paste(paste0("`", z_names, "`"), collapse = " + ")
      } else {
        z_rhs <- paste(sprintf("`%s` * g_fac", z_names), collapse = " + ")
      }
    } else {
      z_rhs <- paste(paste0("`", z_names, "`"), collapse = " + ")
    }
    rhs <- paste0(rhs, " + ", z_rhs)
  }
  if (length(fe_names) > 0) {
    fe_rhs <- paste(paste0("`", fe_names, "`"), collapse = " + ")
    rhs <- paste0(rhs, " + ", fe_rhs)
  }
  form <- stats::as.formula(paste0("Y_out ~ ", rhs))

  fit <- if (!is.null(survey_w)) {
    mgcv::gam(form, data = df_fit, method = gam_method, weights = survey_w)
  } else {
    mgcv::gam(form, data = df_fit, method = gam_method)
  }

  # 2E: weighted means of Z for prediction-frame Z-bar.
  .wmean_col <- function(v, w) {
    if (is.null(w)) mean(v, na.rm = TRUE)
    else stats::weighted.mean(v, w = w, na.rm = TRUE)
  }
  z_bar <- if (length(z_names) > 0)
    vapply(seq_len(ncol(z_mat)),
           function(j) .wmean_col(z_mat[, j], survey_w),
           numeric(1)) else NULL
  if (!is.null(z_bar)) names(z_bar) <- z_names

  # 2E: weighted mode for FE prediction-frame levels.
  fe_modes <- list()
  if (length(fe_names) > 0) {
    for (v in fe_names) {
      if (!is.null(survey_w)) {
        agg <- tapply(survey_w, df_fit[[v]], sum, default = 0)
        fe_modes[[v]] <- names(agg)[which.max(agg)]
      } else {
        tb <- table(df_fit[[v]])
        fe_modes[[v]] <- names(tb)[which.max(tb)]
      }
    }
  }

  # ── Build prediction frame ─────────────────────────────────────────────
  # gv: optional group index (0..K-1) used to set g_fac in the prediction frame.
  # If NULL, default to bin matching xv via quantile cut (back-compat).
  .make_pred_df <- function(dv, xv, gv = NULL) {
    out <- data.frame(D_treat = dv, X_mod = rep(xv, length(dv)))
    if (x_is_bin) {
      out$x_factor <- factor(rep(xv, length(dv)), levels = c(0, 1))
    } else if (!is.null(df_fit$g_fac)) {
      gv_eff <- if (!is.null(gv)) gv else {
        # Best-effort default: assign by xv against group medians.
        if (xv <= stats::median(x_raw[g_loc == 0], na.rm = TRUE)) 0L
        else if (xv >= stats::median(x_raw[g_loc == K_g - 1L], na.rm = TRUE)) K_g - 1L
        else floor(K_g / 2)
      }
      out$g_fac <- factor(rep(gv_eff, length(dv)), levels = 0:(K_g - 1L))
    }
    if (length(z_names) > 0) {
      for (nm in z_names) out[[nm]] <- z_bar[nm]
    }
    if (length(fe_names) > 0) {
      for (nm in fe_names) {
        out[[nm]] <- factor(rep(fe_modes[[nm]], length(dv)),
                            levels = levels(df_fit[[nm]]))
      }
    }
    out
  }

  # 3A/3D: identify low- and high-bin observations using the K-group factor
  # (n_parts-based binning), not hardcoded terciles.
  if (x_is_bin) {
    idx_low <- which(x_raw == 0)
    idx_high <- which(x_raw == 1)
    g_low <- 0L; g_high <- 1L
  } else {
    if (is.null(g_loc))
      stop("Internal error: g_loc missing for continuous X in .gam_fit_and_predict")
    idx_low <- which(g_loc == 0L)
    idx_high <- which(g_loc == (K_g - 1L))
    g_low <- 0L; g_high <- K_g - 1L
  }

  # ── Predicted values on grid ────────────────────────────────────────────
  if (compute_se) {
    pred_X0 <- stats::predict(fit, newdata = .make_pred_df(d_grid, x_lo, g_low),  se.fit = TRUE)
    pred_X1 <- stats::predict(fit, newdata = .make_pred_df(d_grid, x_hi, g_high), se.fit = TRUE)
    yhat_X0 <- as.numeric(pred_X0$fit); se_yhat_X0 <- as.numeric(pred_X0$se.fit)
    yhat_X1 <- as.numeric(pred_X1$fit); se_yhat_X1 <- as.numeric(pred_X1$se.fit)
  } else {
    yhat_X0 <- as.numeric(stats::predict(fit, newdata = .make_pred_df(d_grid, x_lo,  g_low)))
    yhat_X1 <- as.numeric(stats::predict(fit, newdata = .make_pred_df(d_grid, x_hi, g_high)))
    se_yhat_X0 <- rep(NA_real_, length(yhat_X0))
    se_yhat_X1 <- rep(NA_real_, length(yhat_X1))
  }

  # ── Slopes on grid ─────────────────────────────────────────────────────
  h_num <- diff(range(d_grid)) / 1e5

  .slope_at <- function(dv, xv, gv) {
    p_up <- stats::predict(fit, newdata = .make_pred_df(dv + h_num, xv, gv))
    p_dn <- stats::predict(fit, newdata = .make_pred_df(dv - h_num, xv, gv))
    (as.numeric(p_up) - as.numeric(p_dn)) / (2 * h_num)
  }

  slope_grid_X0 <- .slope_at(d_grid, x_lo,  g_low)
  slope_grid_X1 <- .slope_at(d_grid, x_hi, g_high)
  ipe_grid <- slope_grid_X1 - slope_grid_X0

  # ── AIPE, CAME, D-CAME ─────────────────────────────────────────────────
  d_aipe <- if (!is.null(d_all_for_aipe)) d_all_for_aipe else d
  sw_aipe <- if (!is.null(d_all_for_aipe) && !is.null(survey_w_all_for_aipe)) {
    survey_w_all_for_aipe
  } else {
    survey_w
  }

  # AIPE: average IPE over the pooled D distribution. Use the contrast between
  # x_hi (group K-1) and x_lo (group 0) evaluated at every D in d_aipe.
  slope_aipe_X0 <- .slope_at(d_aipe, x_lo,  g_low)
  slope_aipe_X1 <- .slope_at(d_aipe, x_hi, g_high)
  ipe_aipe <- slope_aipe_X1 - slope_aipe_X0

  # CAME slopes evaluated at the observation's own D, with X held at the
  # group's representative level (x_lo or x_hi) and the group factor set
  # appropriately (so Z*g_fac interactions evaluate correctly).
  slope_obs_X0 <- .slope_at(d, x_lo,  g_low)
  slope_obs_X1 <- .slope_at(d, x_hi, g_high)

  if (!is.null(sw_aipe)) {
    sw_norm <- sw_aipe / sum(sw_aipe)
    aipe <- sum(sw_norm * ipe_aipe)
    sw_lo <- survey_w[idx_low]  / sum(survey_w[idx_low])
    sw_hi <- survey_w[idx_high] / sum(survey_w[idx_high])
    came0 <- sum(sw_lo * slope_obs_X0[idx_low])
    came1 <- sum(sw_hi * slope_obs_X1[idx_high])
  } else {
    aipe <- mean(ipe_aipe)
    {
      idx_lo <- idx_low
      idx_hi <- idx_high
      came0 <- mean(slope_obs_X0[idx_lo])
      came1 <- mean(slope_obs_X1[idx_hi])
    }
  }
  dcame <- came1 - came0

  list(
    fit = fit,
    aipe = aipe, dcame = dcame, came0 = came0, came1 = came1,
    yhat_X0 = yhat_X0, yhat_X1 = yhat_X1,
    se_yhat_X0 = se_yhat_X0, se_yhat_X1 = se_yhat_X1,
    slope_grid_X0 = slope_grid_X0, slope_grid_X1 = slope_grid_X1,
    ipe_grid = ipe_grid
  )
}


#' GAM: analytic inference using mgcv posterior covariance.
#' @keywords internal
#' @noRd
.gam_analytic <- function(y, d, x_raw, z_mat, d_grid, x_lo, x_hi,
                          survey_w = NULL, fe_data = NULL, fe_vars = NULL,
                          d_all_for_aipe = NULL,
                          survey_w_all_for_aipe = NULL,
                          smooth = "REML",
                          k_s = NULL, k_t = NULL,
                          inter_gam = "ti",
                          z_interact = FALSE,
                          g_loc = NULL, n_parts = 3L,
                          verbose = TRUE) {

  if (verbose) cat(sprintf("--- GAM fit with analytic (mgcv) inference [smooth=%s] ---\n", smooth))
  est <- .gam_fit_and_predict(y, d, x_raw, z_mat, d_grid, x_lo, x_hi,
                              survey_w, fe_data = fe_data,
                              fe_vars = fe_vars,
                              d_all_for_aipe = d_all_for_aipe,
                              survey_w_all_for_aipe = survey_w_all_for_aipe,
                              smooth = smooth,
                              k_s = k_s, k_t = k_t,
                              inter_gam = inter_gam,
                              z_interact = z_interact,
                              g_loc = g_loc, n_parts = n_parts)

  fit <- est$fit
  n_obs <- length(y)
  n_grid <- length(d_grid)

  # Bayesian posterior covariance
  Vp <- fit$Vp

  x_is_bin <- all(x_raw %in% c(0, 1))
  K_g <- if (x_is_bin) 2L else as.integer(n_parts)
  z_names <- if (!is.null(z_mat) && ncol(z_mat) > 0) colnames(z_mat) else character(0)
  fe_names <- if (!is.null(fe_vars) && length(fe_vars) > 0) fe_vars else character(0)

  # 2E: weighted Z-bar and FE modes
  .wmean_col2 <- function(v, w) {
    if (is.null(w)) mean(v, na.rm = TRUE)
    else stats::weighted.mean(v, w = w, na.rm = TRUE)
  }
  z_bar <- if (length(z_names) > 0)
    vapply(seq_len(ncol(z_mat)),
           function(j) .wmean_col2(z_mat[, j], survey_w),
           numeric(1)) else NULL
  if (!is.null(z_bar)) names(z_bar) <- z_names

  df_fit_ref <- data.frame(Y_out = y, D_treat = d, X_mod = x_raw)
  if (length(z_names) > 0) df_fit_ref <- cbind(df_fit_ref, z_mat)
  if (length(fe_names) > 0) {
    for (v in fe_names) df_fit_ref[[v]] <- factor(fe_data[[v]])
  }
  if (!x_is_bin && !is.null(g_loc))
    df_fit_ref$g_fac <- factor(g_loc, levels = 0:(K_g - 1L))

  fe_modes <- list()
  if (length(fe_names) > 0) {
    for (v in fe_names) {
      if (!is.null(survey_w)) {
        agg <- tapply(survey_w, df_fit_ref[[v]], sum, default = 0)
        fe_modes[[v]] <- names(agg)[which.max(agg)]
      } else {
        tb <- table(df_fit_ref[[v]])
        fe_modes[[v]] <- names(tb)[which.max(tb)]
      }
    }
  }

  .make_pred_df_local <- function(dv, xv, gv = NULL) {
    out <- data.frame(D_treat = dv, X_mod = rep(xv, length(dv)))
    if (x_is_bin) out$x_factor <- factor(rep(xv, length(dv)), levels = c(0, 1))
    if (!x_is_bin && !is.null(df_fit_ref$g_fac) && !is.null(gv))
      out$g_fac <- factor(rep(gv, length(dv)), levels = 0:(K_g - 1L))
    if (length(z_names) > 0) for (nm in z_names) out[[nm]] <- z_bar[nm]
    if (length(fe_names) > 0) {
      for (nm in fe_names)
        out[[nm]] <- factor(rep(fe_modes[[nm]], length(dv)), levels = levels(df_fit_ref[[nm]]))
    }
    out
  }

  h_num <- diff(range(d_grid)) / 1e5

  .get_slope_Lp <- function(dv, xv, gv) {
    Lp_up <- stats::predict(fit, newdata = .make_pred_df_local(dv + h_num, xv, gv), type = "lpmatrix")
    Lp_dn <- stats::predict(fit, newdata = .make_pred_df_local(dv - h_num, xv, gv), type = "lpmatrix")
    (Lp_up - Lp_dn) / (2 * h_num)
  }

  g_low  <- 0L
  g_high <- K_g - 1L
  Sg0 <- .get_slope_Lp(d_grid, x_lo,  g_low)
  Sg1 <- .get_slope_Lp(d_grid, x_hi, g_high)
  Cg <- Sg1 - Sg0

  se_slope_X0 <- sqrt(pmax(0, rowSums((Sg0 %*% Vp) * Sg0)))
  se_slope_X1 <- sqrt(pmax(0, rowSums((Sg1 %*% Vp) * Sg1)))
  se_ipe_grid <- sqrt(pmax(0, rowSums((Cg %*% Vp) * Cg)))

  d_aipe <- if (!is.null(d_all_for_aipe)) d_all_for_aipe else d
  sw_aipe <- if (!is.null(d_all_for_aipe) && !is.null(survey_w_all_for_aipe)) {
    survey_w_all_for_aipe
  } else {
    survey_w
  }

  Sg0_aipe <- .get_slope_Lp(d_aipe, x_lo,  g_low)
  Sg1_aipe <- .get_slope_Lp(d_aipe, x_hi, g_high)
  Cg_aipe <- Sg1_aipe - Sg0_aipe

  if (!is.null(sw_aipe)) {
    sw_norm <- sw_aipe / sum(sw_aipe)
    Cb_aipe <- as.numeric(t(Cg_aipe) %*% sw_norm)
  } else {
    Cb_aipe <- colMeans(Cg_aipe)
  }
  se_aipe <- sqrt(max(0, as.numeric(t(Cb_aipe) %*% Vp %*% Cb_aipe)))

  Sg0_obs <- .get_slope_Lp(d, x_lo,  g_low)
  Sg1_obs <- .get_slope_Lp(d, x_hi, g_high)

  # 3A/3D: bin observations by g_loc (n_parts-based), not by hardcoded terciles.
  if (x_is_bin) {
    idx_low  <- which(x_raw == 0)
    idx_high <- which(x_raw == 1)
  } else {
    idx_low  <- which(g_loc == 0L)
    idx_high <- which(g_loc == (K_g - 1L))
  }
  if (!is.null(survey_w)) {
    sw_lo_a <- survey_w[idx_low]  / sum(survey_w[idx_low])
    sw_hi_a <- survey_w[idx_high] / sum(survey_w[idx_high])
    Cb_came0 <- as.numeric(t(Sg0_obs[idx_low,  , drop = FALSE]) %*% sw_lo_a)
    Cb_came1 <- as.numeric(t(Sg1_obs[idx_high, , drop = FALSE]) %*% sw_hi_a)
  } else {
    Cb_came0 <- colMeans(Sg0_obs[idx_low,  , drop = FALSE])
    Cb_came1 <- colMeans(Sg1_obs[idx_high, , drop = FALSE])
  }
  Cb_dcame <- Cb_came1 - Cb_came0

  se_came0 <- sqrt(max(0, as.numeric(t(Cb_came0) %*% Vp %*% Cb_came0)))
  se_came1 <- sqrt(max(0, as.numeric(t(Cb_came1) %*% Vp %*% Cb_came1)))
  se_dcame <- sqrt(max(0, as.numeric(t(Cb_dcame) %*% Vp %*% Cb_dcame)))

  if (verbose) cat(sprintf("  Analytic: AIPE SE=%.4f, D-CAME SE=%.4f\n",
                           se_aipe, se_dcame))

  est$se_aipe      <- se_aipe
  est$se_dcame     <- se_dcame
  est$se_came0     <- se_came0
  est$se_came1     <- se_came1
  est$se_yhat_X0   <- est$se_yhat_X0
  est$se_yhat_X1   <- est$se_yhat_X1
  est$se_slope_X0  <- se_slope_X0
  est$se_slope_X1  <- se_slope_X1
  est$se_ipe_grid  <- se_ipe_grid
  est$n_boot_valid <- NA_integer_
  est$d_grid       <- d_grid

  est
}


#' GAM: case/cluster bootstrap for SEs.
#' @keywords internal
#' @noRd
.gam_bootstrap <- function(y, d, x_raw, z_mat, d_grid, x_lo, x_hi,
                           cl_var = NULL, survey_w = NULL,
                           B_boot = 500,
                           verbose = TRUE, fe_data = NULL, fe_vars = NULL,
                           d_all_for_aipe = NULL,
                           survey_w_all_for_aipe = NULL,
                           smooth = "REML",
                           k_s = NULL, k_t = NULL,
                           inter_gam = "ti",
                           z_interact = FALSE,
                           g_loc = NULL, n_parts = 3L) {

  if (verbose) cat(sprintf("--- GAM fit on full sample [smooth=%s] ---\n", smooth))
  est_full <- .gam_fit_and_predict(y, d, x_raw, z_mat, d_grid, x_lo, x_hi,
                                   survey_w, fe_data = fe_data,
                                   fe_vars = fe_vars,
                                   d_all_for_aipe = d_all_for_aipe,
                                   survey_w_all_for_aipe = survey_w_all_for_aipe,
                                   smooth = smooth,
                                   k_s = k_s, k_t = k_t,
                                   inter_gam = inter_gam,
                                   compute_se = FALSE,
                                   z_interact = z_interact,
                                   g_loc = g_loc, n_parts = n_parts)

  n_obs <- length(y)
  n_grid <- length(d_grid)

  if (verbose) cat(sprintf("--- Case%s bootstrap (%d replicates) ---\n",
                           if (!is.null(cl_var)) "/cluster" else "", B_boot))

  boot_aipe   <- numeric(B_boot)
  boot_dcame  <- numeric(B_boot)
  boot_came0  <- numeric(B_boot)
  boot_came1  <- numeric(B_boot)
  boot_yh0    <- matrix(NA_real_, B_boot, n_grid)
  boot_yh1    <- matrix(NA_real_, B_boot, n_grid)
  boot_slp0   <- matrix(NA_real_, B_boot, n_grid)
  boot_slp1   <- matrix(NA_real_, B_boot, n_grid)
  boot_ipe    <- matrix(NA_real_, B_boot, n_grid)
  boot_aipe[] <- NA_real_; boot_dcame[] <- NA_real_
  boot_came0[] <- NA_real_; boot_came1[] <- NA_real_

  for (b in 1:B_boot) {
    if (verbose && b %% 50 == 0)
      cat(sprintf("  Case bootstrap replicate %d/%d\n", b, B_boot))

    if (!is.null(cl_var)) {
      ucl <- unique(cl_var)
      idx <- unlist(lapply(sample(ucl, replace = TRUE),
                           function(cl) which(cl_var == cl)))
    } else {
      idx <- sample(n_obs, replace = TRUE)
    }
    yb <- y[idx]; db <- d[idx]; xb <- x_raw[idx]
    zb <- if (!is.null(z_mat)) z_mat[idx, , drop = FALSE] else NULL
    wb <- if (!is.null(survey_w)) survey_w[idx] else NULL
    feb <- if (!is.null(fe_data)) fe_data[idx, , drop = FALSE] else NULL
    gb <- if (!is.null(g_loc)) g_loc[idx] else NULL

    est_b <- tryCatch(
      .gam_fit_and_predict(yb, db, xb, zb, d_grid, x_lo, x_hi, wb,
                           fe_data = feb, fe_vars = fe_vars,
                           d_all_for_aipe = d_all_for_aipe,
                           survey_w_all_for_aipe = survey_w_all_for_aipe,
                           smooth = smooth,
                           k_s = k_s, k_t = k_t,
                           inter_gam = inter_gam,
                           compute_se = FALSE,
                           z_interact = z_interact,
                           g_loc = gb, n_parts = n_parts),
      error = function(e) NULL)
    if (is.null(est_b)) next

    boot_aipe[b]   <- est_b$aipe
    boot_dcame[b]  <- est_b$dcame
    boot_came0[b]  <- est_b$came0
    boot_came1[b]  <- est_b$came1
    boot_yh0[b, ]  <- est_b$yhat_X0
    boot_yh1[b, ]  <- est_b$yhat_X1
    boot_slp0[b, ] <- est_b$slope_grid_X0
    boot_slp1[b, ] <- est_b$slope_grid_X1
    boot_ipe[b, ]  <- est_b$ipe_grid
  }

  vld <- !is.na(boot_aipe)
  n_valid <- sum(vld)

  se_aipe  <- stats::sd(boot_aipe[vld])
  se_dcame <- stats::sd(boot_dcame[vld])
  se_came0 <- stats::sd(boot_came0[vld])
  se_came1 <- stats::sd(boot_came1[vld])
  se_yhat_X0 <- apply(boot_yh0[vld, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
  se_yhat_X1 <- apply(boot_yh1[vld, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
  se_slope_X0 <- apply(boot_slp0[vld, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
  se_slope_X1 <- apply(boot_slp1[vld, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
  se_ipe_grid <- apply(boot_ipe[vld, , drop = FALSE], 2, stats::sd, na.rm = TRUE)

  if (verbose) cat(sprintf("  Case bootstrap: AIPE SE=%.4f, D-CAME SE=%.4f (%d/%d valid)\n",
                           se_aipe, se_dcame, n_valid, B_boot))

  est_full$se_aipe      <- se_aipe
  est_full$se_dcame     <- se_dcame
  est_full$se_came0     <- se_came0
  est_full$se_came1     <- se_came1
  est_full$se_yhat_X0   <- se_yhat_X0
  est_full$se_yhat_X1   <- se_yhat_X1
  est_full$se_slope_X0  <- se_slope_X0
  est_full$se_slope_X1  <- se_slope_X1
  est_full$se_ipe_grid  <- se_ipe_grid
  est_full$n_boot_valid <- n_valid
  est_full$d_grid <- d_grid

  est_full
}


# ============================================================
# GAM DISPATCHER — called from dcame_aipe() when model = "GAM"
# ============================================================

#' GAM dispatcher used by dcame_aipe(model = "GAM").
#' @keywords internal
#' @noRd
.dcame_aipe_gam <- function(y, d, x_raw, data, X_name, Z,
                            cluster, vce, survey_wgts,
                            wgt = NULL,
                            estimand, compare,
                            B_boot,
                            n_grid,
                            level_pv, level_est,
                            n_parts,
                            hist_pv, verbose,
                            Z_FE = NULL,
                            FE = NULL, FE_vbles = NULL,
                            inference = "bootstrap",
                            smooth = "REML",
                            k_s = NULL, k_t = NULL,
                            inter_gam = "ti",
                            Z_interact = FALSE,
                            plot_came_all = FALSE,
                            x_original_labels = NULL) {

  zc_pv  <- stats::qnorm(1 - (1 - level_pv) / 2)
  zc_est <- stats::qnorm(1 - (1 - level_est) / 2)

  x_is_binary <- all(x_raw %in% c(0, 1))
  if (x_is_binary) {
    x_lo <- 0; x_hi <- 1
    g_full <- as.integer(x_raw)  # 0/1
  } else {
    p_lo <- 1 / (2 * n_parts)
    p_hi <- (2 * n_parts - 1) / (2 * n_parts)
    qs <- stats::quantile(x_raw, probs = c(p_lo, p_hi), na.rm = TRUE)
    x_lo <- as.numeric(qs[1]); x_hi <- as.numeric(qs[2])
    # 3A/3D: K-group factor by n_parts quantiles (matches the basis path).
    cut_probs_gam <- seq(0, 1, length.out = as.integer(n_parts) + 1L)[
      -c(1L, as.integer(n_parts) + 1L)]
    cuts_gam <- as.numeric(stats::quantile(x_raw, probs = cut_probs_gam, na.rm = TRUE))
    g_full <- rep(NA_integer_, length(x_raw))
    g_full[x_raw <= cuts_gam[1]] <- 0L
    if (as.integer(n_parts) > 2L)
      for (kk in 2:(as.integer(n_parts) - 1L))
        g_full[x_raw > cuts_gam[kk - 1L] & x_raw <= cuts_gam[kk]] <- as.integer(kk - 1L)
    g_full[x_raw > cuts_gam[length(cuts_gam)]] <- as.integer(as.integer(n_parts) - 1L)
  }

  d_all_for_aipe <- d
  survey_wgts_all_for_aipe <- survey_wgts

  # Z matrix
  z_mat <- if (!is.null(Z)) {
    m <- as.matrix(data[, Z, drop = FALSE])
    storage.mode(m) <- "double"
    m[, apply(m, 2, stats::sd, na.rm = TRUE) > 0, drop = FALSE]
  } else NULL

  fe_data <- NULL
  fe_vars <- NULL

  all_factor_vars <- character(0)
  all_factor_data <- NULL

  if (!is.null(Z_FE) && length(Z_FE) > 0) {
    all_factor_vars <- Z_FE
    all_factor_data <- data[, Z_FE, drop = FALSE]
  }

  if (!is.null(FE) && FE == "Conventional") {
    new_vars <- setdiff(FE_vbles, all_factor_vars)
    all_factor_vars <- c(all_factor_vars, new_vars)
    if (is.null(all_factor_data)) {
      all_factor_data <- data[, new_vars, drop = FALSE]
    } else {
      all_factor_data <- cbind(all_factor_data,
                               data[, new_vars, drop = FALSE])
    }
    if (verbose) cat(sprintf("Conventional FE (GAM): adding %s as factor terms\n",
                             paste(FE_vbles, collapse = ", ")))
  }

  if (length(all_factor_vars) > 0) {
    fe_vars <- all_factor_vars
    fe_data <- all_factor_data
  }

  cluster_var <- if (!is.null(cluster)) data[[cluster]] else NULL

  cs_lo <- min(d, na.rm = TRUE); cs_hi <- max(d, na.rm = TRUE)

  if (x_is_binary) {
    cs_lo <- max(min(d[x_raw == 0], na.rm = TRUE), min(d[x_raw == 1], na.rm = TRUE))
    cs_hi <- min(max(d[x_raw == 0], na.rm = TRUE), max(d[x_raw == 1], na.rm = TRUE))
    if (is.na(cs_lo) || is.na(cs_hi) || cs_lo >= cs_hi)
      stop("No common support between X groups in D. Check data.")
    if (cs_lo >= cs_hi) stop("No common support")
    keep <- (d >= cs_lo) & (d <= cs_hi)
    nt <- sum(!keep)
    if (nt > 0) {
      n_pre <- length(y)
      n0 <- sum(x_raw == 0); n1 <- sum(x_raw == 1)
      nt0 <- sum(!keep & x_raw == 0); nt1 <- sum(!keep & x_raw == 1)
      warning(sprintf(
        "Common-support trimming (GAM, binary X): %d / %d obs (%.1f%%). By X group: X=0: %.1f%%, X=1: %.1f%%.",
        nt, n_pre, 100 * nt / n_pre,
        if (n0 > 0) 100 * nt0 / n0 else 0,
        if (n1 > 0) 100 * nt1 / n1 else 0),
        call. = FALSE)
      y <- y[keep]; d <- d[keep]; x_raw <- x_raw[keep]
      if (!is.null(z_mat)) z_mat <- z_mat[keep, , drop = FALSE]
      if (!is.null(cluster_var)) cluster_var <- cluster_var[keep]
      if (!is.null(survey_wgts)) survey_wgts <- survey_wgts[keep]
      if (!is.null(fe_data)) fe_data <- fe_data[keep, , drop = FALSE]
      g_full <- g_full[keep]
    }
  } else {
    nt <- 0L
  }

  n <- length(y)
  plo_pct <- round(100 / (2 * n_parts))
  phi_pct <- round(100 * (2 * n_parts - 1) / (2 * n_parts))
  mod_label <- if (x_is_binary) "binary" else sprintf("continuous (p%d=%.3f, p%d=%.3f)",
                                                      plo_pct, x_lo, phi_pct, x_hi)

  cat(sprintf("=== DCAME/AIPE Estimation (model = GAM, estimand = %s, inference = %s, moderator = %s) ===\n",
              estimand, inference, mod_label))
  cat(sprintf("N = %d | D range: [%.3f, %.3f]", n, cs_lo, cs_hi))
  if (nt > 0) cat(sprintf(" | %d trimmed", nt))
  if (!is.null(cluster)) cat(sprintf(" | cluster: %s", cluster))
  if (!is.null(survey_wgts)) cat(" | weighted")
  if (!is.null(FE)) cat(sprintf(" | FE: %s (%s)", FE, paste(FE_vbles, collapse = ", ")))
  if (!is.null(fe_vars)) cat(sprintf(" | factor terms: %s", paste(fe_vars, collapse = ", ")))
  cat("\n")

  d_grid <- seq(cs_lo, cs_hi, length.out = n_grid)

  if (inference == "bootstrap") {
    result <- .gam_bootstrap(y = y, d = d, x_raw = x_raw, z_mat = z_mat,
                             d_grid = d_grid, x_lo = x_lo, x_hi = x_hi,
                             cl_var = cluster_var, survey_w = survey_wgts,
                             B_boot = B_boot,
                             verbose = verbose,
                             fe_data = fe_data, fe_vars = fe_vars,
                             d_all_for_aipe = d_all_for_aipe,
                             survey_w_all_for_aipe = survey_wgts_all_for_aipe,
                             smooth = smooth,
                             k_s = k_s, k_t = k_t,
                             inter_gam = inter_gam,
                             z_interact = Z_interact,
                             g_loc = g_full, n_parts = n_parts)
  } else {
    # 2C: GAM analytic SEs use mgcv's Bayesian/frequentist posterior, which
    # ignores within-cluster dependence. Surface this as a real warning so
    # users running with verbose = FALSE still see it.
    if (!is.null(cluster_var)) {
      warning("Analytic GAM standard errors (mgcv posterior) do not account for ",
              "clustering. The `cluster` variable you supplied is being ignored ",
              "for SEs. Use inference = 'bootstrap' for cluster-bootstrap SEs.",
              call. = FALSE)
    }
    result <- .gam_analytic(y = y, d = d, x_raw = x_raw, z_mat = z_mat,
                            d_grid = d_grid, x_lo = x_lo, x_hi = x_hi,
                            survey_w = survey_wgts,
                            fe_data = fe_data, fe_vars = fe_vars,
                            d_all_for_aipe = d_all_for_aipe,
                            survey_w_all_for_aipe = survey_wgts_all_for_aipe,
                            smooth = smooth,
                            k_s = k_s, k_t = k_t,
                            inter_gam = inter_gam,
                            z_interact = Z_interact,
                            g_loc = g_full, n_parts = n_parts,
                            verbose = verbose)
  }

  if (x_is_binary) {
    d_obs_X0 <- d[x_raw == 0]
    d_obs_X1 <- d[x_raw == 1]
    x_plot <- x_raw
  } else {
    # 3A/3D: use the n_parts K-group factor (g_full) consistently.
    d_obs_X0 <- d[g_full == 0L]
    d_obs_X1 <- d[g_full == (as.integer(n_parts) - 1L)]
    # x_plot for the renderer's Linear comparison uses g_full directly (the
    # full K-group index), so the K-group factor in the Linear regression
    # matches what the GAM used.
    x_plot <- g_full
  }

  if (x_is_binary) {
    bX0_label <- "GAM s(D, by=X) @ X=0"
    bX1_label <- "GAM s(D, by=X) @ X=1"
  } else {
    formula_lbl <- if (inter_gam == "te") "te(D,X)" else "s(D)+s(X)+ti(D,X)"
    bX0_label <- sprintf("GAM %s @ X=%.3g (p%d)", formula_lbl, x_lo, plo_pct)
    bX1_label <- sprintf("GAM %s @ X=%.3g (p%d)", formula_lbl, x_hi, phi_pct)
  }
  bX0 <- list(type = "gam", x_value = x_lo)
  bX1 <- list(type = "gam", x_value = x_hi)

  tercile_cuts <- if (!x_is_binary) c(x_lo, x_hi) else NULL
  tercile_cuts_raw <- if (!x_is_binary) stats::quantile(x_raw, probs = c(1/3, 2/3), na.rm = TRUE) else NULL
  n_dropped_t2 <- if (!x_is_binary) sum(is.na(x_plot)) else 0L

  n_fe_display <- 0L
  if (!is.null(fe_vars)) {
    n_fe_display <- sum(sapply(fe_vars, function(v) {
      length(unique(fe_data[[v]])) - 1L
    }))
  }

  # ── Per-X-level CAME for plot_came_all (continuous X only) ───────────
  came_all_df <- NULL
  if (!x_is_binary && isTRUE(plot_came_all)) {
    K_gam <- as.integer(n_parts)
    fit <- result$fit
    # Reproduce the prediction-frame construction used inside .gam_fit_and_predict.
    z_names_g <- if (!is.null(z_mat) && ncol(z_mat) > 0) colnames(z_mat) else character(0)
    fe_names_g <- if (!is.null(fe_vars)) fe_vars else character(0)
    z_bar_g <- if (length(z_names_g) > 0) colMeans(z_mat, na.rm = TRUE) else NULL
    df_fit_ref_g <- data.frame(Y_out = y, D_treat = d, X_mod = x_raw)
    if (length(z_names_g) > 0) df_fit_ref_g <- cbind(df_fit_ref_g, z_mat)
    if (length(fe_names_g) > 0)
      for (v in fe_names_g) df_fit_ref_g[[v]] <- factor(fe_data[[v]])
    fe_modes_g <- list()
    if (length(fe_names_g) > 0)
      for (v in fe_names_g) {
        tb <- table(df_fit_ref_g[[v]])
        fe_modes_g[[v]] <- names(tb)[which.max(tb)]
      }
    .pred_df_g <- function(dv, xv) {
      out <- data.frame(D_treat = dv, X_mod = rep(xv, length(dv)))
      if (length(z_names_g) > 0)
        for (nm in z_names_g) out[[nm]] <- z_bar_g[nm]
      if (length(fe_names_g) > 0)
        for (nm in fe_names_g)
          out[[nm]] <- factor(rep(fe_modes_g[[nm]], length(dv)),
                              levels = levels(df_fit_ref_g[[nm]]))
      out
    }
    h_g <- diff(range(d)) / 1e5
    .slope_at_g <- function(dv, xv) {
      p_up <- stats::predict(fit, newdata = .pred_df_g(dv + h_g, xv))
      p_dn <- stats::predict(fit, newdata = .pred_df_g(dv - h_g, xv))
      (as.numeric(p_up) - as.numeric(p_dn)) / (2 * h_g)
    }

    # Theoretical quantile midpoints at (2k+1)/(2K) for k = 0..K-1.
    qprobs_mid <- (2 * (0:(K_gam - 1L)) + 1) / (2 * K_gam)
    x_mids <- as.numeric(stats::quantile(x_raw, probs = qprobs_mid, na.rm = TRUE))

    # Cut points for grouping observations by bin.
    cut_probs_K <- seq(0, 1, length.out = K_gam + 1L)[-c(1L, K_gam + 1L)]
    cuts_K <- as.numeric(stats::quantile(x_raw, probs = cut_probs_K, na.rm = TRUE))
    g_assign <- rep(NA_integer_, length(x_raw))
    g_assign[x_raw <= cuts_K[1]] <- 0L
    if (K_gam > 2L)
      for (kk in 2:(K_gam - 1L))
        g_assign[x_raw > cuts_K[kk - 1L] & x_raw <= cuts_K[kk]] <- as.integer(kk - 1L)
    g_assign[x_raw > cuts_K[length(cuts_K)]] <- as.integer(K_gam - 1L)

    came_vec <- numeric(K_gam); se_vec <- numeric(K_gam)
    bin_labels <- character(K_gam)
    for (kk in 0:(K_gam - 1L)) {
      idx_kk <- which(g_assign == kk)
      bin_labels[kk + 1L] <- if (kk == 0L) sprintf("X <= %.3g", cuts_K[1])
        else if (kk == K_gam - 1L) sprintf("X > %.3g", cuts_K[length(cuts_K)])
        else sprintf("%.3g < X <= %.3g", cuts_K[kk], cuts_K[kk + 1L])
      if (length(idx_kk) == 0L) {
        came_vec[kk + 1L] <- NA_real_
        se_vec[kk + 1L] <- NA_real_
        next
      }
      # Point estimate: average slope at x = x_mids[kk+1], over group k's D obs.
      slope_obs_kk <- .slope_at_g(d[idx_kk], x_mids[kk + 1L])
      if (!is.null(survey_wgts)) {
        sw_kk <- survey_wgts[idx_kk] / sum(survey_wgts[idx_kk])
        came_vec[kk + 1L] <- sum(sw_kk * slope_obs_kk)
      } else {
        came_vec[kk + 1L] <- mean(slope_obs_kk)
      }
      # SE via posterior Lp matrix (analytic) — if not available, leave NA.
      se_vec[kk + 1L] <- NA_real_
      if (!is.null(fit$Vp)) {
        Lp_up <- stats::predict(fit, newdata = .pred_df_g(d[idx_kk] + h_g, x_mids[kk + 1L]),
                                type = "lpmatrix")
        Lp_dn <- stats::predict(fit, newdata = .pred_df_g(d[idx_kk] - h_g, x_mids[kk + 1L]),
                                type = "lpmatrix")
        Sg_kk <- (Lp_up - Lp_dn) / (2 * h_g)
        Cb_kk <- if (!is.null(survey_wgts)) {
          sw_kk <- survey_wgts[idx_kk] / sum(survey_wgts[idx_kk])
          as.numeric(t(Sg_kk) %*% sw_kk)
        } else colMeans(Sg_kk)
        se_vec[kk + 1L] <- sqrt(max(0, as.numeric(t(Cb_kk) %*% fit$Vp %*% Cb_kk)))
      }
    }
    came_all_df <- data.frame(
      group = 0:(K_gam - 1L),
      group_label = bin_labels,
      X_mid = x_mids,
      n = vapply(0:(K_gam - 1L), function(kk) sum(g_assign == kk, na.rm = TRUE), integer(1)),
      CAME = came_vec,
      SE = se_vec,
      CI_lo = came_vec - zc_est * se_vec,
      CI_hi = came_vec + zc_est * se_vec,
      p_value = ifelse(is.na(se_vec) | se_vec == 0, NA_real_,
                       2 * stats::pnorm(-abs(came_vec / se_vec))),
      stringsAsFactors = FALSE
    )
  } else if (x_is_binary && isTRUE(plot_came_all) && verbose) {
    cat("Note: plot_came_all is ignored when the moderator is binary.\n")
  }

  .render_results(
    model = "GAM",
    estimand = estimand, compare = compare,
    level_pv = level_pv, level_est = level_est,
    zc_pv = zc_pv, zc_est = zc_est,
    aipe = result$aipe, seA = result$se_aipe,
    dcame = result$dcame, se_dcame = result$se_dcame,
    came0 = result$came0, came1 = result$came1,
    se_came0 = result$se_came0, se_came1 = result$se_came1,
    dg = d_grid,
    yh0 = result$yhat_X0, yh1 = result$yhat_X1,
    se0 = result$se_yhat_X0, se1 = result$se_yhat_X1,
    slope0g = result$slope_grid_X0, slope1g = result$slope_grid_X1,
    se_slope0g = result$se_slope_X0, se_slope1g = result$se_slope_X1,
    ipeg = result$ipe_grid, seig = result$se_ipe_grid,
    d_obs_X0 = d_obs_X0, d_obs_X1 = d_obs_X1,
    d_all = d,
    y = y, x = x_plot, z_mat = z_mat,
    cluster_var = cluster_var, survey_wgts = survey_wgts, vce = vce,
    x_is_binary = x_is_binary, tercile_cuts = tercile_cuts,
    n_dropped_t2 = n_dropped_t2,
    bX0_label = bX0_label, bX1_label = bX1_label,
    bX0 = bX0, bX1 = bX1,
    fit_obj = result$fit, V_obj = NULL,
    cs_lo = cs_lo, cs_hi = cs_hi,
    nt = nt, n = n,
    n_grid = n_grid, hist_pv = hist_pv,
    inference = inference,
    selection = paste0("GAM:", smooth),   # 3C: reflect actual smooth criterion
    Z_interact = Z_interact,
    FE = FE, FE_vbles = FE_vbles,
    wgt = wgt,
    Z_FE = Z_FE, n_fe_dummies = n_fe_display,
    fe_mat = NULL,
    x_original_labels = x_original_labels,
    came_all_df = came_all_df,
    plot_came_all = isTRUE(plot_came_all),
    K = if (x_is_binary) 2L else as.integer(n_parts),
    x_raw_full = if (!x_is_binary) x_raw else NULL
  )
}
