# bass_intercept_oracle.R — reference for D03's intercept convention (port brief, fitting segment): the Liang et al.
# convention (flat intercept, centred columns) against level centring (R's convention on y - mean(y): the intercept
# a column inside the g-prior, the prior centred on the level), with R's own zero centre beside them. Fixture F-O9.
# All under Zellner-Siow, g = 1/tau ~ IG(1/2, N/2), and the Jeffreys limit of sigma^2's prior, under which Bayes
# factors between structures on the same data are defined.
#
# Usage (cwd = r/):  Rscript oracles/bass_intercept_oracle.R <fixture.csv> <out.json> <mode>
#   grid:   headerless CSV of N, k, R2 (k centred columns with coefficient of determination R2 on N observations).
#           Per case, ln Bayes factors against the intercept-only structure, over u = ln g:
#             liang  int p(g) (1+g)^((N-1-k)/2) (1 + g(1-R2))^(-(N-1)/2)
#             level  int p(g) (1+g)^((N-1-k)/2) (1 + g(1-R2))^(-N/2)  -  ln int p(g) (1+g)^(-1/2)
#           and each over twice the range, the convergence check.
#   design: headerless CSV of x, y. Structures: the intercept alone, then one hinge (x - 0.3)+, then that and
#           (0.6 - x)+. For each offset added to y (0, 10, -250), ln Bayes factors against the intercept-only
#           structure, from the design itself:
#             level      R's tau-form marginal on y - mean(y), with the uncentred Gram over [1, hinges]:
#                        int p(tau) (tau/(1+tau))^((k+1)/2) (y'y - Q/(1+tau))^(-N/2), Q the fit energy;
#             zero       the same on y as given, R's zero centre;
#             liang      the g-form above, with R2 from the centred hinge columns;
#             level_g    the g-form level formula with that R2, which the explicit level marginal must equal.
# The reference is quadrature; the identity level = level_g is the derivation the D03 alternative rests on.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: bass_intercept_oracle.R <fixture.csv> <out.json> <grid|design>")
fixture <- read_matrix(args[1])
mode <- args[3]

log_ig_u <- function(u, N) 0.5 * log(N / 2) - lgamma(0.5) - 0.5 * u - (N / 2) * exp(-u)   # density of ln g
log_gamma_v <- function(v, N) 0.5 * (v + log(N / 2)) - (N / 2) * exp(v) - lgamma(0.5)      # density of ln tau

# ln of the integral of exp(f) over [-w, w], in 40 pieces about the peak's scale.
log_integral <- function(f, w) {
  grid <- seq(-w, w, length.out = 4001)
  peak <- max(f(grid))
  edges <- seq(-w, w, length.out = 41)
  total <- sum(vapply(1:40, function(i) integrate(function(u) exp(f(u) - peak), edges[i], edges[i + 1],
                                                  rel.tol = 1e-12, abs.tol = 0, subdivisions = 1000L)$value, 0))
  peak + log(total)
}

g_form <- function(N, k, R2, power, w) {
  log_integral(function(u) log_ig_u(u, N) + (N - 1 - k) / 2 * log1p(exp(u)) - power * log1p(exp(u) * (1 - R2)), w)
}
level_null <- function(N, w) log_integral(function(u) log_ig_u(u, N) - 0.5 * log1p(exp(u)), w)

liang_bf <- function(N, k, R2, w) g_form(N, k, R2, (N - 1) / 2, w)
level_bf <- function(N, k, R2, w) g_form(N, k, R2, N / 2, w) - level_null(N, w)

if (mode == "grid") {
  out <- lapply(seq_len(nrow(fixture)), function(i) {
    N <- unname(fixture[i, 1]); k <- unname(fixture[i, 2]); R2 <- unname(fixture[i, 3])
    list(N = N, k = k, R2 = R2,
         liang = liang_bf(N, k, R2, 60), liang_wide = liang_bf(N, k, R2, 120),
         level = level_bf(N, k, R2, 60), level_wide = level_bf(N, k, R2, 120))
  })
  emit_json(out, args[2])
} else if (mode == "design") {
  x <- fixture[, 1]; y0 <- fixture[, 2]; N <- length(x)
  hinges <- list(matrix(numeric(0), N, 0), cbind(pmax(x - 0.3, 0)), cbind(pmax(x - 0.3, 0), pmax(0.6 - x, 0)))

  # R's tau-form marginal of a response under the uncentred Gram over [1, H], as a log, up to the constants every
  # structure shares on the same response. S(tau) = y'y - Q/(1+tau) is computed as |r|^2 + Q*tau/(1+tau), with the
  # residual taken directly: far from zero, y'y and Q agree to many digits (D03's cancellation).
  tau_form <- function(y, H) {
    Z <- cbind(1, H)
    fit <- Z %*% solve(crossprod(Z), crossprod(Z, y))
    Q <- sum(fit * fit); rr <- sum((y - fit)^2); k1 <- ncol(Z)
    log_integral(function(v) log_gamma_v(v, N) + (k1 / 2) * (v - log1p(exp(v))) - (N / 2) * log(rr + Q * plogis(v)), 60)
  }
  centred_r2 <- function(y, H) {
    if (ncol(H) == 0) return(0)
    Hc <- scale(H, scale = FALSE); yc <- y - mean(y)
    sum(yc * (Hc %*% solve(crossprod(Hc), crossprod(Hc, yc)))) / sum(yc * yc)
  }

  out <- lapply(c(0, 10, -250), function(offset) {
    y <- y0 + offset
    yt <- y - mean(y)
    level <- vapply(hinges, function(H) tau_form(yt, H), 0)
    zero <- vapply(hinges, function(H) tau_form(y, H), 0)
    r2 <- vapply(hinges, function(H) centred_r2(y, H), 0)
    ks <- vapply(hinges, ncol, 0)
    list(offset = offset, k = ks, r2 = r2,
         level = level[-1] - level[1], zero = zero[-1] - zero[1],
         liang = vapply(2:3, function(j) liang_bf(N, ks[j], r2[j], 60), 0),
         level_g = vapply(2:3, function(j) level_bf(N, ks[j], r2[j], 60), 0))
  })
  emit_json(out, args[2])
} else stop("mode must be grid or design")
