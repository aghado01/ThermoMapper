# gbass_defects_oracle.R — E1's defects of GBASS 2.0.1 run against the installed package, each confirmed or struck (E1's build brief,
# slice 0). The read list is E1's (issues/regression/bass/extensions-brief.md); each section evaluates GBASS's own code as written
# beside a reference that shares nothing with it.
#
# Usage (cwd = r/):  Rscript oracles/gbass_defects_oracle.R <data.csv> <out.json> [nmcmc] [nburn] [seed]
#   data: headerless CSV, inputs in [0, 1] then y in the last column, for the GBP run (section gbp_step).
#
# Sections:
#   nw_kurtosis: nw_est_mom's fourth moment, its m4 expression as written, against the Normal-Wald's (NIG's) kurtosis from the
#     mixture itself, X = mu + beta V + sqrt(V) Z with V ~ IG(mean delta/gamma, shape delta^2): E[(beta (V - m) + sqrt(V) Z)^k] by
#     quadrature over V's density, at (alpha, beta, delta) = (2, 1, 1) and two more points. corrected is 3 + 3 (1 + 4 beta^2/alpha^2) /
#     (delta gamma).
#   nw_gamma_prior: its two probabilities as written at its starting values (100, 100) and default (q1, q2) = (0.1, 0.9), against
#     the definitions P(xi < q1) = P(|gamma| > q1^-2 - 1) and P(xi > q2) = P(|gamma| < q2^-2 - 1) for gamma ~ N(m, s), xi = (1 +
#     |gamma|)^(-1/2), integrated by quadrature; then the prior it calibrates as written and with the one stray -1 removed (the body
#     deparsed, the expression replaced, parsed again), each with the probabilities it reaches.
#   w_step: the w block of gbass's sweep as written, at a frozen state, 4000 times per arm with the floor w_prior$lower_bound set at
#     the conditional's median: the share of draws below the floor at beta = 0 and beta != 0; and the rMHN arguments at w_prior$a = 0
#     and 5, which the stated conditional separates (exp(-a w / 2) in w).
#   npart: the birth's and the mutation's admission conditions as written, on a column with exactly npart positive rows.
#   gbp_step: hbass on the data, its v step's adapted prop_sigma against the posterior sd of each ln v_i over the kept draws (a
#     random walk's step is about 2.4 times the target's sd in one dimension), the sd over draws of mean(ln v), which the
#     adaptation reads, and the acceptance rate.
#   gamma_sign: nwbass's gamma draw as written at a small state: the normal's mass below 0, where its kernel is not gamma's
#     conditional (the Wald's normalizer is e^|gamma|), which is the distance in total variation to the conditional under a prior on
#     gamma > 0.
#   text_slips: the paper's two slips (E1): r = yhat - y, and Table 1's t row GIG(-nu/2, 0, nu/2). The residual's assignment as
#     written, and tbass's v_prior as written, with the marginal of sqrt(v) zeta (c = 1) by quadrature over v against dt.
#   gig_moments: mu_gig and var_gig as written at a = 0 where the GIG's mean or variance is infinite, and gbass's s2 and bias
#     expressions as written at tbass's v_prior for df in (1, 1.5, 2, 3), with beta = 0 as tbass holds it.

source("oracles/_common.R")
suppressPackageStartupMessages(library(GBASS))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: gbass_defects_oracle.R <data.csv> <out.json> [nmcmc] [nburn] [seed]")
D <- read_matrix(args[1])
nmcmc <- if (length(args) >= 3) as.integer(args[3]) else 3000L
nburn <- if (length(args) >= 4) as.integer(args[4]) else 2000L
seed <- if (length(args) >= 5) as.integer(args[5]) else 20261105L
ns <- asNamespace("GBASS")
gb <- ns$gbass

# ---- nw_kurtosis ----------------------------------------------------------------------------------------------------------------
m4_code <- source_assignment(ns$nw_est_mom, "m4", 1)
nig_kurtosis <- function(alpha, beta, delta) {
  gamma <- sqrt(alpha^2 - beta^2)
  m <- delta / gamma; lambda <- delta^2
  f <- function(v) sqrt(lambda / (2 * pi * v^3)) * exp(-lambda * (v - m)^2 / (2 * m^2 * v))
  E <- function(g) integrate(function(v) f(v) * g(v), 0, Inf, rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value
  c2 <- E(function(v) beta^2 * (v - m)^2 + v)
  c4 <- E(function(v) beta^4 * (v - m)^4 + 6 * beta^2 * (v - m)^2 * v + 3 * v^2)
  list(alpha = alpha, beta = beta, delta = delta, mixture = c4 / c2^2,
       code = eval(m4_code, list(alpha = alpha, beta = beta, delta = delta, gamma = gamma)),
       corrected = 3 + 3 * (1 + 4 * beta^2 / alpha^2) / (delta * gamma))
}
nw_kurtosis <- list(nig_kurtosis(2, 1, 1), nig_kurtosis(3, -2, 0.5), nig_kurtosis(1.5, 0.3, 4))

# ---- nw_gamma_prior ------------------------------------------------------------------------------------------------------------
lhs1_code <- source_assignment(ns$nw_gamma_prior, "lhs1", 1)
lhs2_code <- source_assignment(ns$nw_gamma_prior, "lhs2", 1)
probabilities <- function(m, s, q1, q2) {
  dens <- function(g) dnorm(g, m, s)
  inside <- function(t) integrate(dens, -t, t, rel.tol = 1e-13, abs.tol = 0)$value
  list(below_q1 = 1 - inside(q1^-2 - 1), above_q2 = inside(q2^-2 - 1))
}
at_start <- list(xx = c(100, 100), q1 = 0.1, q2 = 0.9)
fixed_src <- sub("(1 - q2^-2 - 1 - xx[1])", "(1 - q2^-2 - xx[1])", paste(deparse(ns$nw_gamma_prior, width.cutoff = 500L), collapse = "\n"), fixed = TRUE)
nw_gamma_prior_fixed <- eval(parse(text = fixed_src))
environment(nw_gamma_prior_fixed) <- ns
calibrated <- function(fn) {
  par <- fn()
  c(list(m_gamma = unname(par[1]), s_gamma = unname(par[2])), probabilities(par[1], par[2], 0.1, 0.9))
}
nw_gamma_prior_rec <- list(
  code = list(lhs1 = eval(lhs1_code, at_start), lhs2 = eval(lhs2_code, at_start)),
  definition = probabilities(100, 100, 0.1, 0.9),
  corrected_changed = !identical(fixed_src, paste(deparse(ns$nw_gamma_prior, width.cutoff = 500L), collapse = "\n")),
  targets = list(p1 = 0.5, p2 = 0.05),
  as_written = calibrated(ns$nw_gamma_prior), corrected = calibrated(nw_gamma_prior_fixed))

# ---- w_step ---------------------------------------------------------------------------------------------------------------------
w_block <- source_calls(gb, function(e) identical(e[[1]], as.name("if")) && identical(deparse(e[[2]]), "w_prior$type == \"GIG\""))[[1]]
set.seed(seed)
Nw <- 30L; Mw <- 2L
wenv <- function(bet, a_w, lb) {
  e <- new.env(parent = ns)
  assign("w_prior", list(type = "GIG", p = -0.1, a = a_w, b = 0.1, lower_bound = lb), e)
  assign("N", Nw, e); assign("M", Mw, e); assign("scale", 1, e); assign("bet", bet, e)
  assign("shape_shift", (Nw + Mw + 1L) / 2, e); assign("quad_term", 0.1 + 12 + 0.5, e); assign("sum_r", 3, e)
  assign("w", 1e6, e)   # an old w far from the conditional, so a kept old value is plain to see
  e
}
# The conditional's median at beta = 0 by quadrature, GIG(p - (N + M + 1)/2, 0, quad_term), an inverse gamma.
shape0 <- (Nw + Mw + 1L) / 2 + 0.1; scale0 <- (0.1 + 12 + 0.5) / 2
floor_w <- scale0 / qgamma(0.5, shape0)
arm <- function(bet, a_w) {
  out <- vapply(seq_len(4000), function(i) { e <- wenv(bet, a_w, floor_w); eval(w_block, e); get("w", e) }, 0)
  list(beta = bet, a_w = a_w, kept_old = mean(out == 1e6), below_floor = mean(out < floor_w))
}
mhn_args <- function(a_w) {
  e <- wenv(0.4, a_w, floor_w)
  holder <- function() NULL; body(holder) <- w_block
  rc <- source_calls(holder, function(x) identical(x[[1]], as.name("rMHN")))[[1]]
  list(a_w = a_w, alpha = eval(rc$alpha, e), beta = eval(rc$beta, e), gamma = eval(rc$gamma, e))
}
w_step <- list(floor = floor_w, arms = list(arm(0, 0), arm(0.4, 0)), mhn = list(mhn_args(0), mhn_args(5)))

# ---- npart ----------------------------------------------------------------------------------------------------------------------
conds <- source_calls(gb, function(e) identical(e[[1]], as.name("if")) && any(grepl("B_new > 0", deparse(e[[2]]))))
npart_rec <- local({
  npart <- 6; B_new <- c(rep(0.5, 6), rep(0, 14))
  list(npart = npart, positive_rows = sum(B_new > 0), birth_condition = deparse(conds[[1]][[2]]), birth_admits = eval(conds[[1]][[2]]),
       mutation_condition = deparse(conds[[2]][[2]]), mutation_admits = eval(conds[[2]][[2]]))
})

# ---- gbp_step -------------------------------------------------------------------------------------------------------------------
X <- D[, -ncol(D), drop = FALSE]; y <- D[, ncol(D)]
set.seed(seed)
t0 <- proc.time()[["elapsed"]]
fit <- hbass(X, y, nmcmc = nmcmc, nburn = nburn, verbose = FALSE)
secs <- proc.time()[["elapsed"]] - t0
lv <- log(fit$v)
gbp_step <- list(N = length(y), nmcmc = nmcmc, nburn = nburn, seconds = secs, prop_sigma = fit$v_prior$prop_sigma,
                 sd_log_v = unname(apply(lv, 2, sd)), sd_mean_log_v = sd(rowMeans(lv)), acc_rate_v = fit$acc_rate_v)

# ---- gamma_sign -----------------------------------------------------------------------------------------------------------------
gc <- source_assignment(ns$nwbass, "gam", 2)
gamma_sign <- local({
  e <- list(N = 5L, v = c(2, 6, 3, 5, 4), s_gamma = 3, m_gamma = 0.1)
  mu <- eval(gc[[3]], e); sd <- eval(gc[[4]], e)
  list(N = e$N, sum_v = sum(e$v), m_gamma = e$m_gamma, s_gamma = e$s_gamma, mean = mu, sd = sd, mass_below_zero = pnorm(0, mu, sd))
})

# ---- text_slips -----------------------------------------------------------------------------------------------------------------
t_marginal <- function(x, nu, b) {
  # sqrt(v) zeta with v ~ GIG(-nu/2, 0, b) = IG(nu/2, b/2), by quadrature over ln v.
  integrate(function(u) { v <- exp(u); dnorm(x, 0, sqrt(v)) * exp(log(dgamma(1 / v, nu / 2, rate = b / 2)) - 2 * u) * v },
            -40, 40, rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value
}
vp_code <- source_assignment(ns$tbass, "v_prior", 1)
text_slips <- local({
  nu <- 5; vp <- eval(vp_code, list(df = nu)); xs <- c(0, 0.7, 2.5)
  list(residual = deparse(source_assignment(gb, "r", 1)), tbass_v_prior = vp[c("p", "a", "b")],
       x = xs, dt = dt(xs, nu), code_marginal = vapply(xs, t_marginal, 0, nu = nu, b = vp$b),
       table1_marginal = vapply(xs, t_marginal, 0, nu = nu, b = nu / 2))
})

# ---- gig_moments ----------------------------------------------------------------------------------------------------------------
s2_code <- source_assignment(gb, "s2", 1); bias_code <- source_assignment(gb, "bias", 1)
gig_moments <- lapply(c(1, 1.5, 2, 3), function(df) {
  vp <- eval(vp_code, list(df = df)); sc <- ifelse(df > 2, (df - 2) / df, 1)
  mu_v <- tryCatch(ns$mu_gig(vp$p, vp$a, vp$b), error = function(e) NA_real_)
  s2_v <- tryCatch(ns$var_gig(vp$p, vp$a, vp$b), error = function(e) NA_real_)
  e <- list(scale = sc, w = 1, bet = 0, mu_v = mu_v, s2_v = s2_v)
  list(df = df, p = vp$p, mean_exists = -vp$p > 1, var_exists = -vp$p > 2, mu_gig = mu_v, var_gig = s2_v,
       s2 = eval(s2_code, e), bias = eval(bias_code, e))
})

emit_json(list(nw_kurtosis = nw_kurtosis, nw_gamma_prior = nw_gamma_prior_rec, w_step = w_step, npart = npart_rec, gbp_step = gbp_step,
               gamma_sign = gamma_sign, text_slips = text_slips, gig_moments = gig_moments),
          args[2], null = "null", na = "string")
