# bass_jeffreys_oracle.R — defect fixture for R BASS's Jeffreys mode (port brief, fitting segment, D03 "R's Jeffreys
# mode has no single target"; F-O6 (b)).
#
# Usage (cwd = r/):  Rscript oracles/bass_jeffreys_oracle.R <design.csv> <out.json> [sigma2] [g1] [g2]
#   design: headerless CSV, y in the first column, then three columns B with the intercept first (k_max = 3, since
#     the reference integrates over a numerically). Births add columns in order: B[, 1:k] to
#     B[, 1:(k+1)]. The change move replaces column 2 by column 3.
#   sigma2: the noise variance the moves condition on (default 1); g1, g2: R's IG(g1, g2) on sigma^2 (default 0, 0).
#
# The reference is quadrature. For each k and sigma^2, I(k, sigma^2) = int exp(-(||y - B_k a||^2 - ||r_k||^2)/(2 sigma^2)) da
# is taken by a tensor Gauss-Legendre rule (8 panels of 16 nodes per coordinate) on a box of +-10 marginal standard
# deviations around the least-squares fit; the box and scales only place the grid. r_k is the least-squares
# residual, computed directly. The rule agrees with the Gaussian integral's closed form, (2 pi sigma^2)^{k/2}
# |G|^{-1/2}, to rounding on the stage-1 fixture (log_I reports both). Every coefficient prior
# compared here is constant in a, so each prior's marginal is its constant times I:
#   flat:      p(a) = 1;                         per-term factor c (2 pi sigma^2)^{1/2} (|G'|/|G|)^{-1/2}, c = 1
#   jeffreys:  p(a | sigma^2) = |G|^{1/2} sigma^-k; per-term factor (2 pi)^{1/2}
#   r_moves:   sigma^{-2k} |G|^{3/4}, the prior D03 says R's moves imply
#   r_sigma:   sigma^{-(k-1)}, the prior D03 says R's sigma^2 draw implies
# A per-term factor is the marginal ratio of a birth with exp(Delta Q / (2 sigma^2)) divided out, Delta Q from the
# direct residuals; in terms of I it is the prior constants' ratio times I(k+1)/I(k).
#
# R's side is evaluated from the pinned source as written: the `alpha` expressions of birth_des, death_des and
# change_des (rjmcmc_des.R) with beta.prior = 'jeffreys' (beta.jprior.ind = 1, tau = beta.prec = 0), unit inverse
# temperature, and the structure terms neutralized (lambda = nc, lbmcmp = 0, equal move probabilities), with R's
# getQf for the factors; and updateMCMC's `s2.a` (mcmc_update.R) for the sigma^2 shape.
#
# Emits { n, k_max, sigma2, g1, g2, log_I: [ { k, sigma2_grid, quadrature, closed_form } ], birth: [ { k, r, r_death,
#         flat, jeffreys, r_moves_prior, closed_flat, closed_jeffreys, closed_r } ], change: { r, flat, jeffreys,
#         r_moves_prior, closed_flat, closed_r },
#         sigma_shape: { k, r, flat, jeffreys, r_sigma_prior, closed_flat, closed_jeffreys, rate_fit, rate_closed } }.
# All factors are on the log scale.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_jeffreys_oracle.R <design.csv> <out.json> [sigma2] [g1] [g2]")
D <- read_matrix(args[1])
sigma2 <- if (length(args) >= 3) as.numeric(args[3]) else 1
g1 <- if (length(args) >= 4) as.numeric(args[4]) else 0
g2 <- if (length(args) >= 5) as.numeric(args[5]) else 0
y <- D[, 1]; B <- D[, -1, drop = FALSE]
n <- length(y); k_max <- ncol(B)
if (k_max != 3) stop("need exactly three columns after y: intercept, then two more")
bass <- bass_source(c("mcmc_funcs.R", "rjmcmc_des.R", "mcmc_update.R"))

resid_sq <- function(Bk) sum(qr.resid(qr(Bk, tol = 0), y)^2)

# I(k, s2) by a tensor Gauss-Legendre rule, centred and scaled per coordinate.
log_I <- function(Bk, s2, panels = 8, nodes = 16, half_width = 10) {
  k <- ncol(Bk)
  ahat <- qr.coef(qr(Bk, tol = 0), y)
  rk <- resid_sq(Bk)
  scale <- sqrt(s2 * diag(solve(crossprod(Bk))))
  gl <- gauss_legendre(nodes)
  edges <- seq(-half_width, half_width, length.out = panels + 1)
  mid <- (edges[-1] + edges[-(panels + 1)]) / 2; half <- diff(edges) / 2
  z <- c(outer(gl$x, half) + rep(mid, each = nodes)); w <- c(outer(gl$w, half))
  points <- as.matrix(expand.grid(rep(list(seq_along(z)), k)))
  total <- 0
  for (chunk in split(seq_len(nrow(points)), ceiling(seq_len(nrow(points)) / 20000))) {
    idx <- points[chunk, , drop = FALSE]
    A <- ahat + t(matrix(z[idx], ncol = k)) * scale                  # k x points
    q <- (colSums((y - Bk %*% A)^2) - rk) / (2 * s2)
    total <- total + sum(exp(-q) * apply(matrix(w[idx], ncol = k), 1, prod))
  }
  log(total) + sum(log(scale))
}

log_det_G <- function(Bk) as.numeric(determinant(crossprod(Bk))$modulus)
log_prior_const <- function(prior, Bk, s2) {
  k <- ncol(Bk)
  switch(prior,
         flat = 0,
         jeffreys = 0.5 * log_det_G(Bk) - k / 2 * log(s2),
         r_moves = -k * log(s2) + 0.75 * log_det_G(Bk),
         r_sigma = -(k - 1) / 2 * log(s2))
}

grid <- sigma2 * c(0.5, 1, 2)
logI <- lapply(1:k_max, function(k) {
  Bk <- B[, 1:k, drop = FALSE]
  q <- vapply(grid, function(s2) log_I(Bk, s2), 0)
  closed <- k / 2 * log(2 * pi * grid) - 0.5 * log_det_G(Bk)
  list(k = k, sigma2_grid = grid, quadrature = q, closed_form = closed)
})
I_at <- function(k, s2) logI[[k]]$quadrature[match(s2, grid)]

# R's alpha expressions, evaluated as written.
alpha_birth <- source_assignment(bass$birth_des, "alpha")
alpha_death <- source_assignment(bass$death_des, "alpha")
alpha_change <- source_assignment(bass$change_des, "alpha")
r_env <- function(Bcurr, Bnew, nbasis, s2) {
  XtX <- crossprod(Bcurr); Xty <- crossprod(Bcurr, y)
  cur <- bass$getQf(XtX, Xty)
  new <- bass$getQf(crossprod(Bnew), crossprod(Bnew, y))
  if (is.null(cur) || is.null(new)) stop("R's getQf refused a fixture design")
  list2env(list(
    curr = list(s2 = s2, qf = cur$qf, beta.prec = 0, lam = nbasis + 1, nc = nbasis + 1, nbasis = nbasis, R = cur$R,
                temp.ind = 1),
    data = list(itemp.ladder = 1, death.prob.next = 1 / 3, birth.prob = 1 / 3, birth.prob.last = 1 / 3,
                death.prob = 1 / 3),
    prior = list(beta.jprior.ind = 1),
    cand.des = list(lbmcmp = 0), lpbmcmp = 0, qf.cand.list = new, delta_qf = c(new$qf - cur$qf)))
}

birth <- lapply(1:(k_max - 1), function(k) {
  Bk <- B[, 1:k, drop = FALSE]; Bk1 <- B[, 1:(k + 1), drop = FALSE]
  e <- r_env(Bk, Bk1, k - 1, sigma2)
  r <- eval(alpha_birth, e) - e$delta_qf / (2 * sigma2)
  # The death from k+1 columns back to k: curr is the larger model, lam = nbasis so that the count terms vanish.
  ed <- r_env(Bk1, Bk, k, sigma2); ed$curr$lam <- k
  r_death <- eval(alpha_death, ed) - ed$delta_qf / (2 * sigma2)
  factor <- function(prior) log_prior_const(prior, Bk1, sigma2) - log_prior_const(prior, Bk, sigma2) +
    I_at(k + 1, sigma2) - I_at(k, sigma2)
  dG <- log_det_G(Bk1) - log_det_G(Bk)
  list(k = k, r = c(r), r_death = c(r_death), flat = factor("flat"), jeffreys = factor("jeffreys"),
       r_moves_prior = factor("r_moves"),
       closed_flat = 0.5 * log(2 * pi * sigma2) - 0.5 * dG, closed_jeffreys = 0.5 * log(2 * pi),
       closed_r = 0.5 * log(2 * pi) - 0.5 * log(sigma2) + 0.25 * dG)
})

# The change move: column 2 replaced by column 3, so k stays 2.
B_old <- B[, 1:2]; B_new <- B[, c(1, 3)]
ec <- r_env(B_old, B_new, 1, sigma2)
r_change <- eval(alpha_change, ec) - ec$delta_qf / (2 * sigma2)
q_new <- log_I(B_new, sigma2)
change_factor <- function(prior) log_prior_const(prior, B_new, sigma2) - log_prior_const(prior, B_old, sigma2) +
  q_new - I_at(2, sigma2)
change <- list(r = c(r_change), flat = change_factor("flat"), jeffreys = change_factor("jeffreys"),
               r_moves_prior = change_factor("r_moves"),
               closed_flat = -0.5 * (log_det_G(B_new) - log_det_G(B_old)), closed_r = 0.25 * (log_det_G(B_new) - log_det_G(B_old)))

# sigma^2's shape for the k = 2 model: fit log p(sigma^2 | y) = c - (shape + 1) log s - rate / s on the three-point
# grid, where log p = log IG kernel(g1, g2) + log marginal (prior constant, -n/2 log(2 pi s), -||r||^2/(2 s), log I).
k <- 2; Bk <- B[, 1:k]; rk <- resid_sq(Bk)
fit_shape <- function(prior) {
  lp <- vapply(grid, function(s2) -(g1 + 1) * log(s2) - g2 / s2 + log_prior_const(prior, Bk, s2) -
                 n / 2 * log(2 * pi * s2) - rk / (2 * s2) + I_at(k, s2), 0)
  coef <- solve(cbind(1, -log(grid), -1 / grid), lp)
  c(shape = coef[2] - 1, rate = coef[3])
}
fits <- lapply(c(flat = "flat", jeffreys = "jeffreys", r_sigma = "r_sigma"), fit_shape)
s2_shape_r <- eval(source_assignment(bass$updateMCMC, "s2.a"),
                   list2env(list(prior = list(g1 = g1, beta.jprior.ind = 1), data = list(n = n))))
sigma_shape <- list(k = k, r = s2_shape_r, flat = unname(fits$flat["shape"]), jeffreys = unname(fits$jeffreys["shape"]),
                    r_sigma_prior = unname(fits$r_sigma["shape"]), closed_flat = g1 + (n - k) / 2,
                    closed_jeffreys = g1 + n / 2, rate_fit = unname(fits$flat["rate"]), rate_closed = g2 + rk / 2)

emit_json(list(n = n, k_max = k_max, sigma2 = sigma2, g1 = g1, g2 = g2, log_I = logI, birth = birth,
               change = change, sigma_shape = sigma_shape), args[2])
