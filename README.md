# dcameaipe

[![R-CMD-check](https://github.com/SerranoSerrat/dcameaipe/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SerranoSerrat/dcameaipe/actions/workflows/R-CMD-check.yaml)

Estimation of the **Average Interactive Partial Effect (AIPE)** and the
**Difference in Conditional Average Marginal Effects (D-CAME)** for interactive
treatment effects, using flexible (polynomial / B-spline / GAM) functional
forms for the treatment-by-moderator interaction.

Companion package to Serrano-Serrat,
*The Pitfalls of Assuming Linear Treatment Effects in Interaction Tests*,
*Sociological Science*, 2026.

> This package was enhanced and optimized using Claude Opus 4.7 and 4.8.

## Installation

```r
# install.packages("remotes")
remotes::install_github("SerranoSerrat/dcameaipe")
```

Optional dependency `fixest` (only needed for `Z_FE_absorb = TRUE`):

```r
install.packages("fixest")
```

## What the function fits

`dcame_aipe()` estimates the K-group joint model

```
Y_i = sum_{k=0..K-1} { mu_k + f_k(D_i) + alpha_k * X_i } * G_{k,i}
      + delta' Z_i + FE + epsilon_i
```

where `G_{k,i}` is an indicator for the k-th moderator group and `f_k(D)` is a
group-specific polynomial / B-spline (`model = "basis"`) or penalised smooth
(`model = "GAM"`). For continuous `X`, `K = n_parts` quantile bins (default
`n_parts = 3`). Middle groups are kept (they help fit `Z` / FE) but the
headline D-CAME and AIPE contrast group 0 against group K-1.

## Quick start

```r
library(dcameaipe)

# Basic usage: binary X
res <- dcame_aipe(
  data     = mydata,
  Y        = "wage",
  D        = "training_hours",   # continuous treatment
  X        = "union",            # binary moderator (0/1)
  Z        = c("age", "female"), # controls
  estimand = "AIPE",
  compare  = "ALL"               # also show D-CAME and Linear bar
)

summary(res)        # data.frame of all estimands
res$plot            # combined PV / IPE / coef plot
```

```r
# Continuous moderator with fixed effects absorbed via fixest
res2 <- dcame_aipe(
  data         = mydata,
  Y            = "y",
  D            = "treatment",
  X            = "education",         # continuous; binned into terciles
  Z            = c("age", "female"),
  Z_FE         = c("country", "year"),
  Z_FE_absorb  = TRUE,                # high-cardinality FE -> feols
  estimand     = "AIPE",
  compare      = "ALL",
  inference    = "bootstrap",         # cluster-bootstrap if cluster supplied
  cluster      = "country"
)
```

## Key arguments

| Argument | Default | Notes |
|---|---|---|
| `model` | `"basis"` | `"basis"` (poly/B-spline) or `"GAM"` (mgcv penalised). |
| `estimand` | — required — | `"AIPE"` or `"DCAME"`. |
| `inference` | `"crossfit"` | `"crossfit"`, `"bootstrap"`, or `"regular"`. GAM does not support cross-fit. |
| `selection` | `"CV"` | `"CV"`, `"BIC"`, `"CV_targeted"` (AIPE-only). |
| `n_parts` | `3` | Bins for continuous moderator (2–5). Ignored for binary X. |
| `vce` | `"robust"` | `"robust"` (HC1) or `"cluster"`. GAM analytic ignores cluster (use bootstrap). |
| `wgt` | `NULL` | Optional survey weights. |
| `Z_FE` | `NULL` | Additive fixed-effect variables. |
| `Z_FE_absorb` | `FALSE` | Use `fixest::feols` instead of dummies. Basis path only. |
| `Z_interact` | `FALSE` | If `TRUE`, controls interact with the K-group factor. |
| `compare` | `FALSE` | `"LINEAR"`, `"DCAME"`, `"AIPE"`, or `"ALL"`. |

See `?dcame_aipe` for the full argument list.

## Three inference paths

* `inference = "crossfit"` *(basis only; default)* — splits the sample, does
  model selection on one half and estimation on the other; combines folds.
  Best when CV-based spec selection is in play.
* `inference = "bootstrap"` — case (or cluster, if `cluster` is set) bootstrap
  with `B_boot` replicates. Recommended for cluster SEs and for predicted-value
  CIs under `Z_FE_absorb = TRUE`.
* `inference = "regular"` — analytic delta-method (basis) or mgcv posterior
  (GAM). Fast; ignores cluster dependence under GAM.

## Common pitfalls the package warns about

The function emits explicit warnings (not just verbose notes) for:

* `Z_FE_absorb = TRUE` with `model = "GAM"` (silently ignored otherwise).
* `Z_FE_absorb = TRUE` with non-bootstrap inference (PV CIs differ).
* `cluster` set with `vce = "robust"` and no bootstrap (cluster ignored).
* Number of unique clusters < 30 (CR1 SEs badly biased).
* `B_boot` supplied with non-bootstrap inference (ignored).
* `n_parts` supplied with binary X (ignored).
* `selection = "CV_targeted"` with `estimand = "DCAME"` (falls back to plain CV).
* Common-support trimming, with per-group breakdown.
* `model = "GAM"` with `inference = "crossfit"` (overridden to `"regular"`).
* All user_spec_\* supplied with `inference = "crossfit"` (overridden).
* Linear comparison contrast dropped due to collinearity (with diagnostic hints).

## Citation

If you use this package, please cite:

> Serrano-Serrat, J. (2026). The Pitfalls of Assuming Linear Treatment Effects in
> Interaction Tests. *Sociological Science*.

## License

MIT. See `LICENSE`.
