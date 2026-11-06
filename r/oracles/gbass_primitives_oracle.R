# gbass_primitives_oracle.R — GBASS's distributional primitives against their laws (E1's build brief, slice 0).
#
# Usage (cwd = r/):  Rscript oracles/gbass_primitives_oracle.R <cases.csv> <out.json> [draws] [seed]
#   cases: headerless CSV, one row per case: kind, then three parameters.
#     kind 1, GIG(p, a, b), density ∝ x^(p-1) exp(-(a x + b/x)/2): GIGrvg::rgig(lambda = p, psi = a, chi = b) and GBASS's rgig2,
#       which takes 1/rgamma(-p, rate b/2) at a = 0, beside GBASS's mu_gig and var_gig as written.
#     kind 2, the modified half-normal MHN(alpha, beta, gamma), density ∝ x^(alpha-1) exp(-beta x^2 + gamma x) on x > 0 (Sun, Kong
#       and Pal 2021): GBASS's rmhn, with the branch it takes read from its own selectors.
#     kind 3, the modified Bessel function K_nu(x), as (x, nu, 0): base besselK(x, nu, expon.scaled = TRUE).
#   draws: per sampler and case (default 100000); seed: set.seed before each case's draws (default 20261105).
#
# The reference is each law's own integral, never a sampler: the normalized density by adaptive quadrature in logs, split at the
# mode, for the mass, the mean and variance where they exist, and the CDF at the draws' deciles; K_nu by its integral
# representation, ln K_nu(x) = -x + ln ∫_0^∞ exp(-x (cosh t - 1)) cosh(nu t) dt, also in logs. A sampler is right where its deciles
# carry the law's probabilities to within the binomial error.
#
# Emits { gig: [...], mhn: [...], bessel: [...] }, one record per case in input order within its kind:
#   gig: p, a, b, moments (mean_exists, var_exists, mean, var by quadrature; mean_bessel, var_bessel from K where a, b > 0),
#     mu_gig, var_gig (GBASS as written, or the error it raised), and for rgig (where a > 0; null otherwise) and rgig2: deciles
#     (the draws' sample deciles), cdf (the law's CDF at them), n, seconds.
#   mhn: alpha, beta, gamma, branch (negative, sqrt_gamma, normal, series), log_mass, mean, var by quadrature, and the draws' record.
#   bessel: x, nu, log_k_r (ln besselK(x, nu, TRUE) - x), log_k_integral, warned.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: gbass_primitives_oracle.R <cases.csv> <out.json> [draws] [seed]")
cases <- read_matrix(args[1])
draws <- if (length(args) >= 3) as.integer(args[3]) else 100000L
seed <- if (length(args) >= 4) as.integer(args[4]) else 20261105L
ns <- asNamespace("GBASS")

# ∫ exp(lk(x) - lk(mode)) g(x) dx over (0, ∞) by adaptive quadrature in u = ln x, split at the mode, so far ranges keep their digits.
# lk is the log kernel; returns the log of the integral of the kernel times g, for g >= 0.
log_integral <- function(lk, mode, g = function(x) 1) {
  lm <- lk(mode)
  f <- function(u) { x <- exp(u); v <- exp(lk(x) - lm) * g(x) * x; v[!is.finite(v)] <- 0; v }
  u0 <- log(mode)
  # Reach out until the integrand is negligible on each side.
  lo <- u0; while (f(lo - 1) > 1e-300 && lo > u0 - 800) lo <- lo - 1
  hi <- u0; while (f(hi + 1) > 1e-300 && hi < u0 + 800) hi <- hi + 1
  left <- integrate(f, lo - 1, u0, rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value
  right <- integrate(f, u0, hi + 1, rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value
  lm + log(left + right)
}

# The law's moments and CDF from its log kernel and mode.
law <- function(lk, mode, mean_exists, var_exists) {
  lz <- log_integral(lk, mode)
  m1 <- if (mean_exists) exp(log_integral(lk, mode, function(x) x) - lz) else NA_real_
  m2 <- if (var_exists) exp(log_integral(lk, mode, function(x) x^2) - lz) else NA_real_
  cdf <- function(q) vapply(q, function(qq) {
    lm <- lk(mode)
    f <- function(u) { x <- exp(u); v <- exp(lk(x) - lm) * x; v[!is.finite(v)] <- 0; v }
    lo <- log(mode); while (f(lo - 1) > 1e-300 && lo > log(mode) - 800) lo <- lo - 1
    below <- integrate(f, lo - 1, log(qq), rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value
    min(1, exp(lm + log(below) - lz))
  }, 0)
  list(log_mass = lz, mean = m1, var = if (var_exists) m2 - m1^2 else NA_real_, cdf = cdf)
}

sampled <- function(draw, cdf) {
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  x <- draw(draws)
  secs <- proc.time()[["elapsed"]] - t0
  d <- unname(quantile(x, probs = seq(0.1, 0.9, by = 0.1), type = 1))
  list(n = length(x), deciles = d, cdf = cdf(d), seconds = secs)
}

or_error <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))

gig_case <- function(p, a, b) {
  lk <- function(x) (p - 1) * log(x) - (a * x + b / x) / 2
  # The mode of x^(p-1) exp(-(a x + b/x)/2): a x^2 - 2(p - 1) x - b = 0.
  # At b = 0 and p <= 1 the density falls from 0, so the gamma's mean anchors the quadrature instead.
  mode <- if (a > 0) ((p - 1) + sqrt((p - 1)^2 + a * b)) / a else b / (2 * (1 - p))
  if (!(mode > 0)) mode <- 2 * p / a
  mean_exists <- a > 0 || p < -1
  var_exists <- a > 0 || p < -2
  L <- law(lk, mode, mean_exists, var_exists)
  bes <- if (a > 0 && b > 0) {
    s <- sqrt(a * b); k <- function(nu) besselK(s, nu, expon.scaled = TRUE)
    list(mean = sqrt(b / a) * k(p + 1) / k(p), var = (b / a) * (k(p + 2) / k(p) - (k(p + 1) / k(p))^2))
  } else list(mean = NA_real_, var = NA_real_)
  list(p = p, a = a, b = b,
       moments = list(mean_exists = mean_exists, var_exists = var_exists, mean = L$mean, var = L$var, log_mass = L$log_mass,
                      mean_bessel = bes$mean, var_bessel = bes$var),
       mu_gig = or_error(ns$mu_gig(p, a, b)), var_gig = or_error(ns$var_gig(p, a, b)),
       rgig = if (a > 0) sampled(function(n) GIGrvg::rgig(n, lambda = p, chi = b, psi = a), L$cdf) else NULL,
       rgig2 = sampled(function(n) vapply(seq_len(n), function(i) ns$rgig2(p, a, b), 0), L$cdf))
}

mhn_branch <- function(alpha, beta, gamma) {
  if (gamma <= 0) return("negative")
  if (alpha > 1) return(if (ns$.mhn_log_k2(alpha, beta, gamma) <= ns$.mhn_log_k1(alpha, beta, gamma)) "sqrt_gamma" else "normal")
  "series"
}

mhn_case <- function(alpha, beta, gamma) {
  lk <- function(x) (alpha - 1) * log(x) - beta * x^2 + gamma * x
  # The mode solves 2 beta x^2 - gamma x - (alpha - 1) = 0; below alpha = 1 the density falls from 0, so take a scale point.
  mode <- if (alpha > 1) (gamma + sqrt(gamma^2 + 8 * beta * (alpha - 1))) / (4 * beta) else max(gamma / (2 * beta), 1 / sqrt(2 * beta))
  L <- law(lk, mode, TRUE, TRUE)
  list(alpha = alpha, beta = beta, gamma = gamma, branch = mhn_branch(alpha, beta, gamma), log_mass = L$log_mass, mean = L$mean,
       var = L$var, rmhn = sampled(function(n) ns$rmhn(n, alpha, beta, gamma), L$cdf))
}

bessel_case <- function(x, nu) {
  warned <- FALSE
  r <- withCallingHandlers(besselK(x, nu, expon.scaled = TRUE), warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
  # h(t) = -x (cosh t - 1) + ln cosh(nu t), with ln cosh written to hold its digits for large nu t.
  h <- function(t) -x * (cosh(t) - 1) + ifelse(nu * t > 20, nu * t - log(2), log(cosh(nu * t)))
  # Its maximum: x sinh t = nu tanh(nu t), at 0 when nu^2 <= x.
  tstar <- if (nu^2 <= x) 0 else uniroot(function(t) x * sinh(t) - nu * tanh(nu * t), c(1e-12, asinh(nu / x) + 1), tol = 1e-14)$root
  hs <- h(tstar)
  f <- function(t) exp(h(t) - hs)
  hi <- tstar + 1; while (f(hi) > 1e-300) hi <- hi + 1
  parts <- if (tstar > 0) integrate(f, 0, tstar, rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L)$value else 0
  parts <- parts + integrate(f, tstar, hi, rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L)$value
  list(x = x, nu = nu, log_k_r = log(r) - x, log_k_integral = -x + hs + log(parts), warned = warned)
}

gig <- mhn <- bessel <- list()
for (i in seq_len(nrow(cases))) {
  k <- cases[i, 1]; u <- cases[i, 2]; v <- cases[i, 3]; w <- cases[i, 4]
  if (k == 1) gig[[length(gig) + 1]] <- gig_case(u, v, w)
  else if (k == 2) mhn[[length(mhn) + 1]] <- mhn_case(u, v, w)
  else if (k == 3) bessel[[length(bessel) + 1]] <- bessel_case(u, v)
  else stop("unknown kind ", k)
}
emit_json(list(gig = gig, mhn = mhn, bessel = bessel), args[2], null = "null", na = "null")
