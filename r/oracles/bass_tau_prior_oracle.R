# bass_tau_prior_oracle.R — reference for D03's tau hyperprior (port brief, fitting segment): the unit-information
# Gamma family, tau ~ Gamma(a, rate a*N), whose mean is 1/N for every shape a. Zellner-Siow is a = 1/2. Francom et al.
# (2018) state the same intent, "centered over the unit information prior", with Gamma(1/N, 1), the member a = 1/N:
# a shape that shrinks with N spreads ln tau over about N nats, so the prior is centred on 1/N in its mean only.
# The correction keeps the family and fixes the shape independently of N.
#
# Usage (cwd = r/):  Rscript oracles/bass_tau_prior_oracle.R <fixture.csv> <out.json>
#   fixture: headerless CSV, one case per row: N, k, R2. Each case is a k-column structure with coefficient of
#            determination R2 on N observations, under the Liang convention (flat intercept, centred columns).
#
# Emits, per prior ("zs", a = 1/2; "francom", a = 1/N) and per case:
#   log_mean, log_sd       E[ln tau] and sd(ln tau) by quadrature of the density on u = ln tau;
#   log_mean_closed, ...   the same in closed form, psi(shape) - ln(rate) and sqrt(psi'(shape));
#   log_median             ln of the median, from qgamma on the log scale;
#   log_density_unit       the density of ln tau at tau = 1/N;
#   ln_bf                  ln of the Bayes factor of the structure against the intercept-only one,
#                          the integral over u = ln g (g = 1/tau) of p(g) (1+g)^((N-1-k)/2) (1 + g(1-R2))^(-(N-1)/2);
#   ln_bf_wide             the same integral over twice the range, the convergence check.
# Under the Liang convention the intercept-only structure carries no tau, so ln_bf is exact, not up to a constant.
# The reference for the moments is quadrature, checked against the closed forms; the Bayes factor is quadrature
# only, so the checks read its signs and trends, not its digits.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_tau_prior_oracle.R <fixture.csv> <out.json>")
cases <- read_matrix(args[1])

# Each prior is a shape rule for the family Gamma(a, rate a*N).
shapes <- list(
  zs      = function(N) 0.5,
  francom = function(N) 1 / N
)

# Density of u = ln tau for tau ~ Gamma(shape, rate): exp(shape*(u + ln rate) - rate*e^u) / Gamma(shape).
log_density_u <- function(u, a, r) a * (u + log(r)) - r * exp(u) - lgamma(a)

# Moments of u by quadrature. For small shapes the density is exponential in u with rate `shape` below its mode,
# so the range runs 60/shape below it, in pieces one 1/shape wide, which integrate() handles.
log_moments <- function(a, r) {
  lo <- -log(r) - 60 / a; hi <- -log(r) + 5
  edges <- unique(c(seq(lo, hi - 10, by = max(1 / a, 1)), hi - 10, hi))
  piece <- function(f) sum(vapply(seq_len(length(edges) - 1), function(i)
    integrate(f, edges[i], edges[i + 1], rel.tol = 1e-12, abs.tol = 0, subdivisions = 1000L)$value, 0))
  p <- function(u) exp(log_density_u(u, a, r))
  m0 <- piece(p)
  m1 <- piece(function(u) u * p(u)) / m0
  m2 <- piece(function(u) (u - m1)^2 * p(u)) / m0
  c(mass = m0, mean = m1, sd = sqrt(m2))
}

ln_bf <- function(N, k, R2, a, r, half_width) {
  # Over u = ln g, g = 1/tau: the prior density of u is the density of ln tau at -u.
  f <- function(u) exp(log_density_u(-u, a, r) + ((N - 1 - k) / 2) * log1p(exp(u)) -
                         ((N - 1) / 2) * log1p(exp(u) * (1 - R2)))
  # The log integrand's peak sets a scale so the integral stays representable.
  grid <- seq(-half_width, half_width, length.out = 4001)
  lf <- log_density_u(-grid, a, r) + ((N - 1 - k) / 2) * log1p(exp(grid)) - ((N - 1) / 2) * log1p(exp(grid) * (1 - R2))
  peak <- max(lf)
  g <- function(u) exp(log_density_u(-u, a, r) + ((N - 1 - k) / 2) * log1p(exp(u)) -
                         ((N - 1) / 2) * log1p(exp(u) * (1 - R2)) - peak)
  edges <- seq(-half_width, half_width, length.out = 41)
  total <- sum(vapply(1:40, function(i) integrate(g, edges[i], edges[i + 1], rel.tol = 1e-12, abs.tol = 0,
                                                  subdivisions = 1000L)$value, 0))
  peak + log(total)
}

out <- list()
for (name in names(shapes)) {
  rows <- lapply(seq_len(nrow(cases)), function(i) {
    N <- unname(cases[i, 1]); k <- unname(cases[i, 2]); R2 <- unname(cases[i, 3])
    a <- shapes[[name]](N); r <- a * N
    mom <- log_moments(a, r)
    list(N = N, k = k, R2 = R2, shape = a, rate = r,
         mass = mom[["mass"]], log_mean = mom[["mean"]], log_sd = mom[["sd"]],
         log_mean_closed = digamma(a) - log(r), log_sd_closed = sqrt(trigamma(a)),
         log_median = log(qgamma(0.5, a, rate = r)),
         log_density_unit = log_density_u(log(1 / N), a, r),
         ln_bf = ln_bf(N, k, R2, a, r, 100), ln_bf_wide = ln_bf(N, k, R2, a, r, 200))
  })
  out[[name]] <- rows
}
emit_json(out, args[2])
