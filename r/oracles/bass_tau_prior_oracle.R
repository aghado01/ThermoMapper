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
#
# The shared family's other members (decided 2026-10-29), each with a = 3: hyper-g (Liang, Paulo, Molina, Clyde and
# Berger 2008) on g = 1/tau, p(g) = ((a - 2)/2)(1 + g)^(-a/2), and hyper-g/n, the same law on g/N. For each, per case:
#   mass, log_mean, log_sd   the density of u = ln tau by quadrature (Jacobian included);
#   beta_mean                E[g'/(1 + g')] by quadrature, against its closed form 2/a (g'/(1 + g') ~ Beta(1, a/2 - 1));
#   ln_bf, ln_bf_wide        the Bayes factor against the intercept-only structure, by quadrature as above;
#   ln_bf_closed             hyper-g's closed form, Liang et al.'s eq. (17): (a - 2)/(k + a - 2) 2F1((N - 1)/2, 1;
#                            (k + a)/2; R2), its series summed in logs; null for hyper-g/n, which has none this simple.

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

ln_bf <- function(N, k, R2, a, r, half_width, density = function(v) log_density_u(v, a, r)) {
  # Over u = ln g, g = 1/tau: the prior density of u is the density of ln tau at -u.
  # The log integrand's peak sets a scale so the integral stays representable.
  grid <- seq(-half_width, half_width, length.out = 4001)
  lf <- density(-grid) + ((N - 1 - k) / 2) * log1p(exp(grid)) - ((N - 1) / 2) * log1p(exp(grid) * (1 - R2))
  peak <- max(lf)
  g <- function(u) exp(density(-u) + ((N - 1 - k) / 2) * log1p(exp(u)) -
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
# hyper-g's density of u = ln tau, with g' = s/tau: s = 1 (hyper-g) or 1/N (hyper-g/n).
log_density_u_hg <- function(u, a, s) log((a - 2) / 2) - (a / 2) * log1p(s * exp(-u)) + log(s) - u

# ln 2F1(alpha, 1; gamma; z) for 0 <= z < 1, its series summed in logs: term k is (alpha)_k / (gamma)_k z^k.
log_2f1_unit <- function(alpha, gamma, z) {
  k <- 0:200000
  lt <- lgamma(alpha + k) - lgamma(alpha) - lgamma(gamma + k) + lgamma(gamma) + k * log(z)
  m <- max(lt)
  if (lt[length(lt)] > m - 50) stop("2F1's series has not converged")
  m + log(sum(exp(lt - m)))
}

hg_moments <- function(a, s) {
  centre <- log(s); lo <- centre - 60 / (a / 2 - 1); hi <- centre + 60
  edges <- seq(lo, hi, length.out = ceiling(hi - lo) + 1)
  piece <- function(f) sum(vapply(seq_len(length(edges) - 1), function(i)
    integrate(f, edges[i], edges[i + 1], rel.tol = 1e-12, abs.tol = 0, subdivisions = 1000L)$value, 0))
  p <- function(u) exp(log_density_u_hg(u, a, s))
  m0 <- piece(p)
  m1 <- piece(function(u) u * p(u)) / m0
  m2 <- piece(function(u) (u - m1)^2 * p(u)) / m0
  b1 <- piece(function(u) s / (exp(u) + s) * p(u)) / m0
  c(mass = m0, mean = m1, sd = sqrt(m2), beta_mean = b1)
}

hyper <- list(hyperg = function(N) 1, hypergn = function(N) 1 / N)
a_hg <- 3
for (name in names(hyper)) {
  out[[name]] <- lapply(seq_len(nrow(cases)), function(i) {
    N <- unname(cases[i, 1]); k <- unname(cases[i, 2]); R2 <- unname(cases[i, 3])
    s <- hyper[[name]](N)
    mom <- hg_moments(a_hg, s)
    density <- function(v) log_density_u_hg(v, a_hg, s)
    list(N = N, k = k, R2 = R2, a = a_hg, scale = s, mass = mom[["mass"]], log_mean = mom[["mean"]], log_sd = mom[["sd"]],
         beta_mean = mom[["beta_mean"]], beta_mean_closed = 2 / a_hg,
         ln_bf = ln_bf(N, k, R2, a_hg, NA, 100, density), ln_bf_wide = ln_bf(N, k, R2, a_hg, NA, 200, density),
         ln_bf_closed = if (name == "hyperg") log((a_hg - 2) / (k + a_hg - 2)) + log_2f1_unit((N - 1) / 2, (k + a_hg) / 2, R2) else NA)
  })
}
emit_json(out, args[2])
