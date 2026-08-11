###############################################################################
##  Does the kernel estimator recover the CAME when the treatment is not
##  conditionally normal?
##
##  Companion code for "Further Discussion" in the dcameaipe documentation.
##  https://serranoserrat.github.io/dcameaipe/discussion.html
##
##  The outcome equation is the one used by Hainmueller, Liu, Liu, Mummolo and
##  Xu (2025, Figure 2), following Simonsohn (2024):
##
##        Y = D^2 - 0.5 D + e
##
##  Only the distribution of D changes across the three designs below. In every
##  design corr(D, X) = 0.5, E[D] = 0 and Var(D) = 1, so the designs differ in
##  shape alone -- not in location, scale, or strength of dependence.
###############################################################################

library(MASS)
library(interflex)

set.seed(20260811)

n   <- 20000
rho <- 0.5


## ---------------------------------------------------------------------------
## 1. Three data-generating processes
## ---------------------------------------------------------------------------
##
##  (a) "normal"       D ~ N(0, 1).  The benchmark. D | X is normal, hence
##                     symmetric at every x.
##
##  (b) "copula"       D ~ Uniform(-sqrt3, sqrt3) marginally -- same mean and
##                     variance as (a), and not skewed unconditionally. But
##                     because D is a monotone transform of a normal that is
##                     correlated with X, the *conditional* distribution D | X
##                     is squeezed against whichever boundary X pushes it
##                     toward: right-skewed for low x, left-skewed for high x,
##                     symmetric only at x = 0.
##
##  (c) "sym_uniform"  D = rho X + uniform noise. Non-normal, but D | X is
##                     symmetric at every x. This is the placebo: it isolates
##                     conditional *symmetry* from normality.

draw <- function(n, design = c("normal", "copula", "sym_uniform"), rho = 0.5) {

  design <- match.arg(design)

  S  <- matrix(c(1, rho, rho, 1), 2, 2)
  WX <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = S)
  W  <- WX[, 1]
  X  <- WX[, 2]

  D <- switch(design,
    normal      = W,
    copula      = sqrt(3) * (2 * pnorm(W) - 1),
    sym_uniform = {
      half <- sqrt(3 * (1 - rho^2))          # keeps Var(D) = 1
      rho * X + runif(n, -half, half)
    })

  data.frame(Y = D^2 - 0.5 * D + rnorm(n), D = D, X = X)
}


## ---------------------------------------------------------------------------
## 2. The true CAME
## ---------------------------------------------------------------------------
##
##  dY/dD = 2D - 0.5, so
##
##        CAME(x) = E[dY/dD | X = x] = 2 E[D | X = x] - 0.5
##
##  (a) and (c):  E[D | X = x] = rho * x
##  (b)           W | X = x ~ N(rho x, 1 - rho^2) and, for W ~ N(mu, s2),
##                E[Phi(W)] = Phi(mu / sqrt(1 + s2)). Hence
##                E[D | X = x] = sqrt(3) * (2 Phi(rho x / sqrt(2 - rho^2)) - 1),
##                which is S-shaped in x rather than linear.

came_true <- function(x, design = c("normal", "copula", "sym_uniform"),
                      rho = 0.5) {

  design <- match.arg(design)

  ED <- switch(design,
    normal      = rho * x,
    copula      = sqrt(3) * (2 * pnorm(rho * x / sqrt(2 - rho^2)) - 1),
    sym_uniform = rho * x)

  2 * ED - 0.5
}


## ---------------------------------------------------------------------------
## 3. What the kernel estimator converges to
## ---------------------------------------------------------------------------
##
##  As the bandwidth shrinks, local linear regression of Y on D within a
##  neighbourhood of x returns Cov(D, Y | X = x) / Var(D | X = x). With
##  Y = D^2 - 0.5 D and E[D^3 | X] = mu^3 + 3 mu sigma^2 + gamma1 sigma^3,
##
##        plim  =  2 mu(x) - 0.5  +  gamma1(x) sigma(x)
##                 \_____________/    \_______________/
##                    true CAME             bias
##
##  where gamma1(x) is the skewness of D | X = x. The bias vanishes whenever
##  D | X is conditionally symmetric -- under normality, but not only there.
##  This function evaluates it by simulating the conditional distribution.

kernel_plim <- function(x, design, rho = 0.5, m = 2e6) {

  d <- draw(m, design, rho)                       # unconditional draw
  keep <- abs(d$X - x) < 0.02                     # thin slice around x
  Dx <- d$D[keep]

  mu <- mean(Dx); s <- sd(Dx)
  g1 <- mean((Dx - mu)^3) / s^3
  c(true = 2 * mu - 0.5, plim = 2 * mu - 0.5 + g1 * s, skew = g1)
}


## ---------------------------------------------------------------------------
## 4. Estimate with interflex and compare
## ---------------------------------------------------------------------------

run_design <- function(design) {

  d <- draw(n, design)

  out <- interflex(estimator = "kernel",
                   data      = d,
                   Y         = "Y",
                   D         = "D",
                   X         = "X",
                   Ylabel    = "Y",
                   Dlabel    = "D",
                   Xlabel    = "X",
                   vartype   = "bootstrap",
                   main      = paste("Kernel model --", design))

  print(out)

  ## interflex returns the estimated CAME on a grid of X. The slot is
  ## est.kernel in current versions -- run names(out) if this errors.
  est  <- as.data.frame(out$est.kernel)
  grid <- est[[1]]
  me   <- est[["ME"]]

  ## Representative moderator values: the median of the bottom and of the top
  ## tercile of X, as recommended for continuous moderators.
  x_lo <- unname(quantile(d$X, 1/6))
  x_hi <- unname(quantile(d$X, 5/6))

  hat  <- approx(grid, me, xout = c(x_lo, x_hi))$y
  tru  <- came_true(c(x_lo, x_hi), design)

  data.frame(
    design    = design,
    tercile   = c("bottom", "top"),
    x         = round(c(x_lo, x_hi), 3),
    came_true = round(tru, 3),
    came_hat  = round(hat, 3),
    row.names = NULL)
}

results <- do.call(rbind, lapply(c("normal", "copula", "sym_uniform"),
                                 run_design))
print(results)


## D-CAME: the contrast between the two terciles, true versus estimated.
dcame <- do.call(rbind, lapply(split(results, results$design), function(r) {
  data.frame(design      = r$design[1],
             dcame_true  = diff(r$came_true),
             dcame_hat   = diff(r$came_hat),
             recovered   = round(diff(r$came_hat) / diff(r$came_true), 3))
}))
print(dcame)


## ---------------------------------------------------------------------------
## 5. The same conclusion analytically, with no estimation noise
## ---------------------------------------------------------------------------

for (dg in c("normal", "copula", "sym_uniform")) {
  cat("\n---", dg, "---\n")
  for (x in c(-0.967, 0, 0.967)) {
    v <- kernel_plim(x, dg)
    cat(sprintf("  x = %6.3f   skew(D|x) = %6.3f   true = %7.4f   plim = %7.4f\n",
                x, v["skew"], v["true"], v["plim"]))
  }
}
