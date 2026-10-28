# bass_moments_oracle.R — reference moments for BASS's sensitivity analysis (port brief, outputs segment, D08; O-O3).
#
# Usage (cwd = r/):  Rscript oracles/bass_moments_oracle.R <cases.csv> <out.json> [continuous|categorical]
#   continuous (default): headerless CSV, one case per row: t1, s1, t2, s2, q. With s2 = 0 the case is the first
#     moment E[h1]; otherwise the cross moment E[h1 h2], where h_i(x) = (s_i (x - t_i))_+^q is the unscaled hinge
#     (R divides each hinge by const(); that scaling is a convention outside these values).
#   categorical: headerless CSV, first row the level probabilities (L columns, zeros allowed), then one row per
#     subset as 0/1 indicators over the levels.
#
# Continuous references, under each measure below (all truncated to an interval and renormalized):
#   quadrature: adaptive quadrature (integrate, rel.tol 1e-13) over the measure's interval split at the knots, the
#     normalizer taken the same way;
#   gauss_legendre: a composite Gauss-Legendre rule, 64 subpanels of 20 nodes between consecutive breakpoints, as a
#     second, independent quadrature; max_rel_quadrature_vs_gl is their largest relative difference.
#   r: R's own moments, C1 and C2 (sobol_normTmix.R), times const() to undo R's scaling, where R supports the
#     measure; r_error records why not (R refuses normal degrees above 1, and its Student moments need hypergeo,
#     which is not installed). max_rel_r_vs_quadrature over the supported cases.
# Measures: uniform_unit U(0, 1); uniform_inner U(0.25, 0.75); normal_mix, weights (0.4, 0.6), means (0.3, 0.7), sds
# (0.1, 0.2) on [0, 1], weights before truncation and the whole mixture renormalized (R's convention, D08's default);
# normal_mix_conditional, the same with each component renormalized on [0, 1] (D08's other named form);
# normal_narrow, N(0.5, 0.3^2) on [0.4, 0.6]; student_nu1, student_nu2, student_nu3, location 0.5 and scale 0.15
# (R's location-scale `sd`) on [0, 1]; student_mix, nu (1, 3), locations (0.3, 0.7), scales (0.1, 0.2), weights
# (0.5, 0.5) on [0, 1].
#
# Categorical references, by enumeration: P(D_i) and P(D_i and D_j) as sums of the supplied probabilities, beside
# R's subset size (or intersection size) over the level count, which ignores the probabilities.
#
# Emits continuous: { measures, cases, quadrature, gauss_legendre, r, r_error, max_rel_quadrature_vs_gl,
#         max_rel_r_vs_quadrature } (matrices are cases x measures; null where R does not support the measure);
#       categorical: { probabilities, subsets, P, P_joint, r_P, r_P_joint }.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_moments_oracle.R <cases.csv> <out.json> [continuous|categorical]")
C <- read_matrix(args[1])
mode <- if (length(args) >= 3) args[3] else "continuous"

if (mode == "categorical") {
  w <- C[1, ]; Dm <- C[-1, , drop = FALSE] == 1
  L <- length(w); nD <- nrow(Dm)
  P <- vapply(seq_len(nD), function(i) sum(w[Dm[i, ]]), 0)
  joint <- outer(seq_len(nD), seq_len(nD), Vectorize(function(i, j) sum(w[Dm[i, ] & Dm[j, ]])))
  r_joint <- outer(seq_len(nD), seq_len(nD), Vectorize(function(i, j) sum(Dm[i, ] & Dm[j, ]) / L))
  emit_json(list(probabilities = w, subsets = Dm * 1, P = P, P_joint = joint, r_P = rowSums(Dm) / L, r_P_joint = r_joint),
            args[2])
  quit(save = "no")
}

bass <- bass_source(c("basis_funcs.R", "sobol_normTmix.R"))
mixture <- function(name, dist, lo, hi, weights = NULL, mean = NULL, sd = NULL, df = NULL, conditional = FALSE)
  list(name = name, dist = dist, trunc = c(lo, hi), weights = weights, mean = mean, sd = sd, df = df, conditional = conditional)
measures <- list(
  mixture("uniform_unit", "uniform", 0, 1),
  mixture("uniform_inner", "uniform", 0.25, 0.75),
  mixture("normal_mix", "normal", 0, 1, c(0.4, 0.6), c(0.3, 0.7), c(0.1, 0.2)),
  mixture("normal_mix_conditional", "normal", 0, 1, c(0.4, 0.6), c(0.3, 0.7), c(0.1, 0.2), conditional = TRUE),
  mixture("normal_narrow", "normal", 0.4, 0.6, 1, 0.5, 0.3),
  mixture("student_nu1", "student", 0, 1, 1, 0.5, 0.15, 1),
  mixture("student_nu2", "student", 0, 1, 1, 0.5, 0.15, 2),
  mixture("student_nu3", "student", 0, 1, 1, 0.5, 0.15, 3),
  mixture("student_mix", "student", 0, 1, c(0.5, 0.5), c(0.3, 0.7), c(0.1, 0.2), c(1, 3)))

component_density <- function(mu, k) switch(mu$dist,
  normal = function(x) dnorm(x, mu$mean[k], mu$sd[k]),
  student = function(x) dt((x - mu$mean[k]) / mu$sd[k], mu$df[k]) / mu$sd[k])
# The unnormalized density on the interval; components renormalized one by one when conditional.
raw_density <- function(mu) {
  if (mu$dist == "uniform") return(function(x) rep(1, length(x)))
  comps <- lapply(seq_along(mu$weights), function(k) component_density(mu, k))
  scale <- if (mu$conditional) vapply(comps, function(f) integrate(f, mu$trunc[1], mu$trunc[2], rel.tol = 1e-13)$value, 0)
           else rep(1, length(comps))
  function(x) rowSums(vapply(seq_along(comps), function(k) mu$weights[k] * comps[[k]](x) / scale[k], numeric(length(x))))
}
hinge <- function(x, t, s, q) pmax(s * (x - t), 0)^q
integrand <- function(cs) function(x) {
  v <- hinge(x, cs[1], cs[2], cs[5])
  if (cs[4] != 0) v <- v * hinge(x, cs[3], cs[4], cs[5])
  v
}
breakpoints <- function(mu, cs) {
  b <- sort(unique(c(mu$trunc, cs[1], if (cs[4] != 0) cs[3])))
  b[b >= mu$trunc[1] & b <= mu$trunc[2]]
}
by_quadrature <- function(f, b) sum(vapply(seq_len(length(b) - 1), function(i)
  integrate(f, b[i], b[i + 1], rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L)$value, 0))
gl <- gauss_legendre(20)
by_gauss_legendre <- function(f, b) sum(vapply(seq_len(length(b) - 1), function(i) {
  edges <- seq(b[i], b[i + 1], length.out = 65)
  x <- c(outer(gl$x, diff(edges) / 2) + rep((edges[-1] + edges[-65]) / 2, each = 20))
  w <- c(outer(gl$w, diff(edges) / 2))
  sum(w * f(x))
}, 0))

# R's moments: the measure as sobol() hands it to C1/C2 (normal weights divided by the truncated mass), and the
# result times const() to undo R's hinge scaling.
r_prior <- function(mu) {
  if (mu$dist == "uniform") return(structure(list(trunc = mu$trunc), class = "uniform"))
  z <- pnorm(mu$trunc[2], mu$mean, mu$sd) - pnorm(mu$trunc[1], mu$mean, mu$sd)
  structure(list(trunc = mu$trunc, mean = mu$mean, sd = mu$sd, weights = mu$weights / sum(mu$weights * z)), class = "normal")
}
r_moment <- function(mu, cs) {
  tl <- list(q = cs[5], t = matrix(c(cs[1], cs[3]), 2, 1), s = matrix(c(cs[2], cs[4]), 2, 1), prior = list(r_prior(mu)))
  if (cs[4] == 0) return(bass$C1(1, 1, tl) * bass$const(cs[2], cs[1], cs[5]))
  bass$C2(1, 2, 1, tl) * bass$const(c(cs[2], cs[4]), c(cs[1], cs[3]), cs[5])
}
r_supported <- function(mu) mu$dist == "uniform" || (mu$dist == "normal" && !mu$conditional)

cases <- C
quad <- glm_ <- r <- matrix(NA_real_, nrow(cases), length(measures))
r_error <- rep(NA_character_, length(measures))
for (j in seq_along(measures)) {
  mu <- measures[[j]]
  dens <- raw_density(mu)
  Zq <- by_quadrature(dens, mu$trunc); Zg <- by_gauss_legendre(dens, mu$trunc)
  for (i in seq_len(nrow(cases))) {
    cs <- cases[i, ]
    f <- integrand(cs); g <- function(x) f(x) * dens(x)
    b <- breakpoints(mu, cs)
    quad[i, j] <- by_quadrature(g, b) / Zq
    glm_[i, j] <- by_gauss_legendre(g, b) / Zg
    if (r_supported(mu)) {
      val <- tryCatch(r_moment(mu, cs), error = function(e) conditionMessage(e))
      if (is.character(val)) { if (is.na(r_error[j])) r_error[j] <- val } else r[i, j] <- val
    }
  }
  if (!r_supported(mu)) r_error[j] <- if (mu$dist == "student") "Student moments need hypergeo (not installed)" else
    "R has no form with weights conditional on truncation"
}
rel <- function(a, b) { d <- abs(a - b); s <- pmax(abs(a), abs(b)); ifelse(s == 0, 0, d / s) }
emit_json(list(measures = lapply(measures, function(mu) mu[!vapply(mu, is.null, TRUE)]), cases = cases,
               quadrature = quad, gauss_legendre = glm_, r = r, r_error = r_error,
               max_rel_quadrature_vs_gl = max(rel(quad, glm_)),
               max_rel_r_vs_quadrature = max(rel(r, quad), na.rm = TRUE)), args[2], na = "null")
