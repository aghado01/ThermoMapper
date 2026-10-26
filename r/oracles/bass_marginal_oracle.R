# bass_marginal_oracle.R — reference for BASS's collapsed g-prior marginal over tau (port brief, D03 propriety).
#
# Usage (cwd = r/):  Rscript oracles/bass_marginal_oracle.R <fixture.csv> <out.json> [convention] [g1] [g2]
#   fixture: headerless CSV; column 1 is y, the rest are the design B. Under convention "r" B includes the
#            intercept column; under "liang" it holds only the other columns, which the oracle centres.
#   convention: "r" (default): the intercept inside the g-prior, k = M + 1, as R BASS does.
#               "liang": a flat intercept and the g-prior on the centred columns, k = M (Liang et al. 2008).
#   g1, g2: the IG(g1, g2) prior on sigma^2; the default 0 0 is the improper 1/sigma^2.
#
# Emits { convention, N, N_eff, k, Q, rss, exact_fit, lower, log_integral, slope, predicted_slope }.
# log_integral[i] is the log of the integral over log tau, from lower[i] to 1e4, of
#   p(tau) (tau/(1+tau))^(k/2) (g2 + S(tau)/2)^-(g1 + N_eff/2),   S(tau) = rss + tau/(1+tau) Q,
# with tau ~ Gamma(1/2, rate N/2), the Zellner-Siow prior: the structure's marginal after a and sigma^2 are
# integrated out, up to constants every structure shares. Under "liang" N_eff = N - 1 and y and the columns are
# centred; in g = 1/tau this is Liang et al.'s (1+g)^((N-1-M)/2) (1 + g(1-R^2))^(-(N-1)/2).
# A sequence that settles is a finite marginal. An exact fit (rss = 0) under the improper prior grows at
# predicted_slope = (N_eff - k - 1)/2 per unit of log(1/lower), the divergence that makes the posterior improper.
# rss is computed from the residual itself, never as y'y - Q, which cancels; exact_fit means
# |r| <= 1e3 * eps * |y|, the working-precision zero.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_marginal_oracle.R <fixture.csv> <out.json> [convention] [g1] [g2]")

F <- read_matrix(args[1])
convention <- if (length(args) >= 3) args[3] else "r"
g1 <- if (length(args) >= 4) as.numeric(args[4]) else 0
g2 <- if (length(args) >= 5) as.numeric(args[5]) else 0
if (!convention %in% c("r", "liang")) stop("convention must be 'r' or 'liang'")

y <- F[, 1]
B <- F[, -1, drop = FALSE]
N <- length(y)
if (convention == "liang") {
  y <- y - mean(y)
  B <- scale(B, center = TRUE, scale = FALSE)
  N_eff <- N - 1
} else {
  N_eff <- N
}
k <- ncol(B)

fit <- lm.fit(B, y)
if (fit$rank < k) stop("design is rank deficient")
rss <- sum(fit$residuals^2)
Q <- sum(fit$fitted.values^2)
exact_fit <- sqrt(rss) <= 1e3 * .Machine$double.eps * sqrt(sum(y^2))
if (exact_fit) rss <- 0

log_integrand <- function(tau) {
  dgamma(tau, shape = 0.5, rate = N / 2, log = TRUE) + k / 2 * log(tau / (1 + tau)) -
    (g1 + N_eff / 2) * log(g2 + (rss + tau / (1 + tau) * Q) / 2)
}

log_integral <- function(lower, upper = 1e4) {
  u <- seq(log(lower), log(upper), length.out = 4001)
  shift <- max(log_integrand(exp(u)) + u)                 # the grid maximum, so exp() never overflows
  f <- function(u) exp(log_integrand(exp(u)) + u - shift)
  log(integrate(f, log(lower), log(upper), subdivisions = 4000L, rel.tol = 1e-8)$value) + shift
}

lower <- 10^-c(4, 8, 16, 32, 64)
vals <- vapply(lower, log_integral, numeric(1))
n <- length(lower)
slope <- (vals[n] - vals[n - 1]) / (log(lower[n - 1]) - log(lower[n]))

emit_json(list(
  convention      = convention,
  N               = N,
  N_eff           = N_eff,
  k               = k,
  Q               = Q,
  rss             = rss,
  exact_fit       = exact_fit,
  lower           = lower,
  log_integral    = vals,
  slope           = slope,
  predicted_slope = (N_eff - k - 1) / 2
), args[2])
