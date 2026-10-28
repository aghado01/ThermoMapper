# bass_sensitivity_defects_oracle.R — defect fixtures for R BASS's sensitivity analysis (port brief, outputs segment,
# D08; O-O2).
#
# Usage (cwd = r/):  Rscript oracles/bass_sensitivity_defects_oracle.R <basis.csv> <out.json>
#   basis: headerless CSV, one row per output location: y_m, y_s, then the basis columns K_1, K_2 (m x 4).
#
# Each section runs the pinned R code (sobol_normTmix.R, bassPCA.R, basis_funcs.R) as written and records R's result
# beside the reference:
#   category_weights: R's categorical moments are subset sizes over level counts. add_tlCat, C2Cat (direct path) and
#     C1Basis, C2Basis (basis path) on one variable with levels (a, b, c) and subsets D1 = {a}, D2 = {a, b}, against
#     enumeration under supplied level probabilities, (0, 1/2, 1/2) (a zero-probability level) and (0.9, 0.05, 0.05);
#     and the inventory's two-level case, P(a) = 0.9, where E[1(z in {a})^2] is 0.9 and R gives 0.5. static: whether
#     each of those functions reads prior.cat, where sobol stores the supplied weights.
#   truncation_bounds: C22Basis (basis path) integrates hinge products over (t2, 1), (t1, t2) or (0, t1) whatever the
#     measure's interval; C2 (direct path) clips to it. Both against adaptive quadrature over the clipped interval,
#     under Uniform(0.25, 0.75), for sign pairs (+, +), (+, -), (-, -) and knot pairs (0.2, 0.4), (0.3, 0.6). All on
#     R's scale, the hinge products divided by R's const.
#   decomposition_units: sobolBasis on a mock two-component fit whose dat$y_s is not 1. The components' surfaces, on
#     p = 2 inputs with R's scaled degree-1 hinges and the default uniform measure on [0, 1]^2:
#       w_1 = 0.5 + 1.2 h(x1; +, 0.3) - 0.8 h(x2; -, 0.6)
#       w_2 = -0.2 + 0.9 h(x1; +, 0.3) h(x2; +, 0.4) + 0.5 h(x2; +, 0.2).
#     The reference is a tensor Gauss-Legendre rule split at the knots (exact for these piecewise polynomials) for
#     the variance of y_m + y_s * (K w) at each location, its main effects and the interaction: R's Var.tot and S.var
#     against the decomposition-unit and the physical-unit references, and R's indices S against the unit-free ones.
#   effects: sobol_des's effect functions and returned mean use C1.temp, the first moment of a scaled hinge under
#     Uniform(0, 1), whatever the declared measure; its effect curves call makeBasis with degree 1. C1.temp
#     evaluated as written against R's own measure-aware C1 and quadrature, under Uniform(0.25, 0.75) and a
#     truncated normal mixture, for (s, t) = (+, 0.3) and (-, 0.6); and the degree argument of every makeBasis call
#     in sobol_des, with a degree-2 hinge's effect curve evaluated as written against the degree-2 hinge.
#
# Emits { category_weights, truncation_bounds, decomposition_units, effects }.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_sensitivity_defects_oracle.R <basis.csv> <out.json>")
Bm <- read_matrix(args[1])
y_m <- Bm[, 1]; y_s <- Bm[, 2]; K <- Bm[, 3:4]
m <- nrow(Bm)
bass <- bass_source(c("misc.R", "basis_funcs.R", "sobol_normTmix.R", "bassPCA.R"))
uniform <- function(lo, hi) structure(list(dist = "uniform", trunc = c(lo, hi)), class = "uniform")

# ---- category weights -----------------------------------------------------------------------------------------------
levels3 <- c("a", "b", "c")
subsets <- list(D1 = "a", D2 = c("a", "b"))
cat_mod <- list(pcat = 1, pdes = 0, nlevels = 3, vars.cat = array(1, c(1, 2, 1)), sub.size = array(c(1, 2), c(1, 2, 1)),
                sub.list = list(list(list("a"), list(c("a", "b")))))
tl <- bass$add_tlCat(list(M = 2, pfunc = 0), cat_mod, 1, 1)
enumerate <- function(w, D) sum(w[levels3 %in% D])
category_case <- function(w) list(
  weights = w,
  enumerated = list(P_D1 = enumerate(w, "a"), P_D2 = enumerate(w, c("a", "b")), P_D1_and_D2 = enumerate(w, "a")))
pc_cat <- list(c(cat_mod, list(vars.cat = array(1, c(1, 2, 1)))))
direct <- list(sub_cnt = c(tl$sub.cnt), C2Cat_D1_D2 = bass$C2Cat(1, 1, 2, tl), C2Cat_D1_D1 = bass$C2Cat(1, 1, 1, tl))
basis_path <- list(C1Basis_D1 = bass$C1Basis(list(), pc_cat, 1, 1, 1, 1), C1Basis_D2 = bass$C1Basis(list(), pc_cat, 1, 2, 1, 1),
                   C2Basis_D1_D2 = bass$C2Basis(list(), pc_cat, 1, 1, 2, 1, 1, 1, 1))
tl2 <- bass$add_tlCat(list(M = 1, pfunc = 0),
                      list(pcat = 1, pdes = 0, nlevels = 2, vars.cat = array(1, c(1, 1, 1)), sub.size = array(1, c(1, 1, 1)),
                           sub.list = list(list(list("a")))), 1, 1)
reads_prior_cat <- vapply(c("add_tlCat", "C2Cat", "C1Basis", "C2Basis"),
                          function(f) "prior.cat" %in% all.names(body(get(f, envir = bass))), TRUE)
category_weights <- list(
  r_direct = direct, r_basis = basis_path,
  cases = list(category_case(c(0, 0.5, 0.5)), category_case(c(0.9, 0.05, 0.05))),
  two_level = list(weights = c(0.9, 0.1), r_second_moment = bass$C2Cat(1, 1, 1, tl2), enumerated = 0.9),
  static_reads_prior_cat = as.list(reads_prior_cat))

# ---- truncation bounds ----------------------------------------------------------------------------------------------
pos <- function(v) pmax(v, 0)
measure <- uniform(0.25, 0.75)
bound_case <- function(t1, t2, s1, s2) {
  cc <- bass$const(signs = c(s1, s2), knots = c(t1, t2), degree = 1)
  f <- function(x) pos(s1 * (x - t1)) * pos(s2 * (x - t2)) / (measure$trunc[2] - measure$trunc[1])
  brk <- sort(unique(c(measure$trunc, t1, t2)))
  brk <- brk[brk >= measure$trunc[1] & brk <= measure$trunc[2]]
  ref <- sum(vapply(seq_len(length(brk) - 1), function(i)
    integrate(f, brk[i], brk[i + 1], rel.tol = 1e-13, abs.tol = 0)$value, 0)) / cc
  tlb <- list(q = 1, t = matrix(c(t1, t2), 2, 1), s = matrix(c(s1, s2), 2, 1), prior = list(measure))
  list(t = c(t1, t2), s = c(s1, s2), r_basis = bass$C22Basis(measure, t1, t2, s1, s2, 1, 1, 2, 1, 1),
       r_direct = bass$C2(1, 1, 2, tlb), quadrature = ref)
}
truncation_bounds <- list(measure = measure$trunc, cases = lapply(
  list(c(0.2, 0.4, 1, 1), c(0.2, 0.4, 1, -1), c(0.2, 0.4, -1, -1), c(0.3, 0.6, 1, 1), c(0.3, 0.6, 1, -1),
       c(0.3, 0.6, -1, -1)), function(a) bound_case(a[1], a[2], a[3], a[4])))

# ---- decomposition units --------------------------------------------------------------------------------------------
# Terms per component: rows of (variable, sign, knot) factors; the intercept is beta[1].
components <- list(
  list(beta = c(0.5, 1.2, -0.8), terms = list(list(v = 1, s = 1, t = 0.3), list(v = 2, s = -1, t = 0.6))),
  list(beta = c(-0.2, 0.9, 0.5), terms = list(list(v = c(1, 2), s = c(1, 1), t = c(0.3, 0.4)), list(v = 2, s = 1, t = 0.2))))
knot_values <- sort(unique(unlist(lapply(components, function(cmp) unlist(lapply(cmp$terms, `[[`, "t"))))))
xx_des <- cbind(knot_values, knot_values)            # knotInd indexes rows; the same table serves both inputs
mock_component <- function(cmp) {
  M <- length(cmp$terms)
  vars <- signs <- kind <- array(NA, c(1, M, 2))
  for (j in seq_len(M)) {
    tm <- cmp$terms[[j]]; k <- length(tm$v)
    vars[1, j, 1:k] <- tm$v; signs[1, j, 1:k] <- tm$s; kind[1, j, 1:k] <- match(tm$t, knot_values)
  }
  list(s2 = 1, p = 2, pdes = 2, pcat = 0, des = TRUE, cat = FALSE, func = FALSE, nbasis = M, model.lookup = 1,
       beta = matrix(cmp$beta, 1), vars.des = vars, signs = signs, knotInd.des = kind, xx.des = xx_des, degree = 1)
}
mod <- structure(list(mod.list = lapply(components, mock_component), dat = list(basis = K, y.m = y_m, y.s = y_s)),
                 class = "bassBasis")
sob <- bass$sobolBasis(mod, int.order = 2, mcmc.use = 1, verbose = FALSE)

# Reference: a tensor Gauss-Legendre rule on [0, 1]^2, panels split at every knot.
panels <- sort(unique(c(0, knot_values, 1)))
gl <- gauss_legendre(6)
nodes <- unlist(lapply(seq_len(length(panels) - 1), function(i) (panels[i] + panels[i + 1]) / 2 + diff(panels[i + 0:1]) / 2 * gl$x))
weights <- unlist(lapply(seq_len(length(panels) - 1), function(i) diff(panels[i + 0:1]) / 2 * gl$w))
surface <- function(cmp, x1, x2) {
  X <- rbind(x1, x2)
  out <- rep(cmp$beta[1], length(x1))
  for (j in seq_along(cmp$terms)) {
    tm <- cmp$terms[[j]]
    out <- out + cmp$beta[j + 1] * bass$makeBasis(tm$s, tm$v, tm$t, X, 1)
  }
  out
}
g <- expand.grid(i = seq_along(nodes), j = seq_along(nodes))
W <- matrix(sapply(components, function(cmp) surface(cmp, nodes[g$i], nodes[g$j])), ncol = 2)
wt <- weights[g$i] * weights[g$j]
F_dec <- W %*% t(K)                                   # grid points x locations, decomposition units
mean_dec <- colSums(wt * F_dec)
var_dec <- colSums(wt * sweep(F_dec, 2, mean_dec)^2)
cond_mean <- function(by) {                           # E[F | x_by] at each node of x_by, then its variance
  idx <- if (by == 1) g$i else g$j
  M <- rowsum(wt * F_dec, idx) / rowsum(wt, idx)[, 1]
  colSums(weights * sweep(M, 2, mean_dec)^2)
}
main_dec <- cbind(cond_mean(1), cond_mean(2))
inter_dec <- var_dec - rowSums(main_dec)
ref_dec <- cbind(main_dec, inter_dec)
r_var <- t(matrix(sob$S.var, ncol = m))               # locations x effects ("1", "2", "1x2" and "other")
r_S <- t(matrix(sob$S, ncol = m))
decomposition_units <- list(
  names = sob$names.ind, y_s = y_s,
  r_var_total = sob$Var.tot, reference_var_total_decomposition = var_dec, reference_var_total_physical = y_s^2 * var_dec,
  max_rel_r_vs_decomposition = max(abs(sob$Var.tot / var_dec - 1)),
  max_rel_r_vs_physical = max(abs(sob$Var.tot / (y_s^2 * var_dec) - 1)),
  max_rel_effects_r_vs_decomposition = max(abs(r_var[, 1:3] - ref_dec)) / max(var_dec),
  max_rel_effects_r_vs_physical = max(abs(r_var[, 1:3] - y_s^2 * ref_dec)) / max(y_s^2 * var_dec),
  max_abs_indices_r_vs_reference = max(abs(r_S[, 1:3] - ref_dec / var_dec)))

# ---- effects --------------------------------------------------------------------------------------------------------
sobol_des <- bass$sobol_des
c1_temp <- source_assignment(sobol_des, "C1.temp")
normal_mix <- structure(list(dist = "normal", trunc = c(0, 1), mean = c(0.3, 0.7), sd = c(0.1, 0.2),
                             weights = c(0.4, 0.6)), class = "normal")
normal_mix$weights <- normal_mix$weights /
  sum(normal_mix$weights * (pnorm(1, normal_mix$mean, normal_mix$sd) - pnorm(0, normal_mix$mean, normal_mix$sd)))
density_of <- function(mu) if (inherits(mu, "uniform")) function(x) dunif(x, mu$trunc[1], mu$trunc[2]) else
  function(x) rowSums(sapply(seq_along(mu$weights), function(k) mu$weights[k] * dnorm(x, mu$mean[k], mu$sd[k])))
effect_case <- function(mu, s, t) {
  tle <- list(q = 1, t = matrix(t), s = matrix(s), prior = list(mu))
  cc <- bass$const(s, t, 1)
  f <- function(x) pos(s * (x - t)) / cc * density_of(mu)(x)
  ref <- sum(vapply(list(c(mu$trunc[1], t), c(t, mu$trunc[2])), function(r)
    if (r[2] > r[1]) integrate(f, r[1], r[2], rel.tol = 1e-13, abs.tol = 0)$value else 0, 0))
  list(measure = mu$dist, trunc = mu$trunc, s = s, t = t, r_c1_temp = c(eval(c1_temp, list2env(list(tl = tle)))),
       r_c1 = bass$C1(1, 1, tle), quadrature = ref)
}
calls <- source_calls(sobol_des, function(e) identical(e[[1]], as.name("makeBasis")))
degree_args <- vapply(calls, function(e) paste(deparse(e[[6]]), collapse = ""), "")
xg <- seq(0, 1, length.out = 101)
effect_env <- list2env(list(tl = list(s = matrix(1), t = matrix(0.3), q = 2), m = 1, pp.use = 1, xxt = t(xg),
                            makeBasis = bass$makeBasis))
r_curve <- eval(calls[[1]], effect_env)
effects <- list(
  first_moment = list(effect_case(uniform(0.25, 0.75), 1, 0.3), effect_case(uniform(0.25, 0.75), -1, 0.6),
                      effect_case(normal_mix, 1, 0.3), effect_case(normal_mix, -1, 0.6),
                      effect_case(uniform(0, 1), 1, 0.3)),
  makebasis_degree_args = degree_args,
  degree2_curve = list(r_max_abs_vs_degree2 = max(abs(r_curve - bass$makeBasis(1, 1, 0.3, t(xg), 2))),
                       r_max_abs_vs_degree1 = max(abs(r_curve - bass$makeBasis(1, 1, 0.3, t(xg), 1)))))

emit_json(list(category_weights = category_weights, truncation_bounds = truncation_bounds,
               decomposition_units = decomposition_units, effects = effects), args[2])
