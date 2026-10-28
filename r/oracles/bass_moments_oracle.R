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
#   r: R's own moments, C1 and C2 (sobol_normTmix.R), times const() to undo R's scaling, for every measure R has
#     (all but the conditional mixture); r_error records the first error R raises for the measure. C1 and C2 return
#     0 without evaluating anything when a hinge misses the measure's interval, so R "covers" such cases trivially.
#     r_covered and r_max_rel count R's finite values and give its largest relative error per measure, for first and
#     cross moments at degrees 1-3 (null where R has none); max_rel_r_vs_quadrature is the largest over the uniform
#     and normal measures.
#   r_repaired: R's Student forms (intx1Student, intx2Student) with only the 2F1 evaluator replaced. As written,
#     robust2f1 calls hypergeo::f15.3.8, the Abramowitz-Stegun 15.3.8 transformation: for |x| < 1 it returns a
#     complex value, which fails C1's and C2's sign checks, and for |x| >= 1 R passes it the Pfaff-transformed
#     argument in (0, 1), where it returns NA (the commented-out code used gsl::hyperg_2F1 there; hypergeo's general
#     hypergeo() does not converge at these arguments either). The repair keeps R's branches, Pfaff transform
#     included, and sums the Gauss series, which converges for |x| < 1, to relative 1e-17. It separates R's
#     evaluator from its formula: repaired_covered and repaired_max_rel as above.
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
# Emits continuous: { measures, cases, quadrature, gauss_legendre, r, r_error, groups, r_covered, r_max_rel,
#         r_repaired, repaired_covered, repaired_max_rel, max_rel_quadrature_vs_gl, max_rel_r_vs_quadrature }
#         (matrices are cases x measures; the summaries are measures x groups, a group being first or cross moments
#         at one degree; null where R has no value);
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
repaired <- bass_source(c("basis_funcs.R", "sobol_normTmix.R"))
gauss_series <- function(a, b, c, x) {                         # 2F1(a, b; c; x) for |x| < 1
  term <- 1; total <- 1; n <- 0
  repeat {
    term <- term * (a + n) * (b + n) / ((c + n) * (n + 1)) * x
    total <- total + term; n <- n + 1
    if (abs(term) <= 1e-17 * abs(total) || n >= 1e6) break
  }
  total
}
assign("robust2f1", function(a, b, c, x)
  if (abs(x) < 1) gauss_series(a, b, c, x) else gauss_series(a, c - b, c, 1 - 1 / (1 - x)) / (1 - x)^a, envir = repaired)
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

# R's moments: the measure as sobol() hands it to C1/C2 (weights divided by the truncated mass), and the result times
# const() to undo R's hinge scaling.
r_prior <- function(mu) {
  if (mu$dist == "uniform") return(structure(list(trunc = mu$trunc), class = "uniform"))
  z <- if (mu$dist == "normal") pnorm(mu$trunc[2], mu$mean, mu$sd) - pnorm(mu$trunc[1], mu$mean, mu$sd) else
    pt((mu$trunc[2] - mu$mean) / mu$sd, mu$df) - pt((mu$trunc[1] - mu$mean) / mu$sd, mu$df)
  structure(list(trunc = mu$trunc, mean = mu$mean, sd = mu$sd, df = mu$df, weights = mu$weights / sum(mu$weights * z)),
            class = mu$dist)
}
r_moment <- function(mu, cs, env = bass) {
  tl <- list(q = cs[5], t = matrix(c(cs[1], cs[3]), 2, 1), s = matrix(c(cs[2], cs[4]), 2, 1), prior = list(r_prior(mu)))
  if (cs[4] == 0) return(env$C1(1, 1, tl) * env$const(cs[2], cs[1], cs[5]))
  env$C2(1, 2, 1, tl) * env$const(c(cs[2], cs[4]), c(cs[1], cs[3]), cs[5])
}
r_supported <- function(mu) !mu$conditional

cases <- C
quad <- glm_ <- r <- rep_ <- matrix(NA_real_, nrow(cases), length(measures))
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
      val <- tryCatch(suppressWarnings(r_moment(mu, cs)), error = function(e) conditionMessage(e))
      if (is.complex(val)) val <- "a complex value"
      if (is.character(val)) { if (is.na(r_error[j])) r_error[j] <- val } else r[i, j] <- val
    }
    if (mu$dist == "student") {
      rep_[i, j] <- tryCatch(r_moment(mu, cs, repaired), error = function(e) NA_real_)
    }
  }
  if (!r_supported(mu)) r_error[j] <- "R has no form with weights conditional on truncation"
}
rel <- function(a, b) { d <- abs(a - b); s <- pmax(abs(a), abs(b)); ifelse(s == 0, 0, d / s) }
groups <- paste0(ifelse(cases[, 4] == 0, "first", "cross"), "_q", cases[, 5])
group_names <- sort(unique(groups))
measure_names <- vapply(measures, `[[`, "", "name")
summarize <- function(values) {
  finite <- ifelse(is.finite(values), values, NA)
  err <- rel(finite, quad)
  cell <- function(f) {
    out <- t(vapply(seq_along(measures), function(j) vapply(group_names, function(g) f(err[groups == g, j]), 0),
                    numeric(length(group_names))))
    dimnames(out) <- list(measure_names, group_names)
    out
  }
  list(values = finite, covered = cell(function(v) sum(!is.na(v))),
       max_rel = cell(function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)), err = err)
}
rs <- summarize(r); rr <- summarize(rep_)
uniform_normal <- vapply(measures, function(mu) mu$dist %in% c("uniform", "normal"), TRUE)
emit_json(list(measures = lapply(measures, function(mu) mu[!vapply(mu, is.null, TRUE)]), cases = cases,
               quadrature = quad, gauss_legendre = glm_, r = rs$values, r_error = r_error,
               groups = group_names, r_covered = rs$covered, r_max_rel = rs$max_rel, r_repaired = rr$values,
               repaired_covered = rr$covered, repaired_max_rel = rr$max_rel,
               max_rel_quadrature_vs_gl = max(rel(quad, glm_)),
               max_rel_r_vs_quadrature = max(rs$err[, uniform_normal], na.rm = TRUE)), args[2], na = "null")
