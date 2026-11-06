# gbass_conditionals_oracle.R — GBASS's full conditionals at a frozen state, as written, against its stated joint (E1's build brief,
# slice 0).
#
# Usage (cwd = r/):  Rscript oracles/gbass_conditionals_oracle.R <state.csv> <out.json> key=value ...
#   state: headerless CSV, one row per observation: y_i, v_i, then the basis row B_i0 .. B_iM (B_i0 = 1, the intercept, inside the
#     prior as GBASS holds it).
#   keys (each a number; a is ';'-separated): a (the M + 1 coefficients), w, beta, tau, gamma, scale (c), w_p, w_a, w_b (w's GIG
#     prior), v_p, v_a, v_b (v's GIG prior), a_tau, b_tau, m_beta, s_beta, m_gamma, s_gamma, lower_bound.
#
# The model GBASS states (Rumsey, Francom and Shen, arXiv 2306.01911, with the code's own choices where it has them):
#   y_i = (B a)_i + sqrt(w) (beta v_i + sqrt(c v_i) zeta_i), zeta_i ~ N(0, 1);  a ~ N(0, w tau I_{M+1});  v_i ~ GIG(v_p, v_a, v_b);
#   w ~ GIG(w_p, w_a, w_b);  tau ~ IG(a_tau, b_tau);  beta ~ N(m_beta, s_beta^2);
#   nwbass: v_i ~ GIG(-1/2, gamma^2, 1), the Wald with normalizer e^|gamma| / sqrt(2 pi), and gamma ~ N(m_gamma, s_gamma^2).
# GIG(p, a, b) ∝ x^(p-1) exp(-(a x + b/x)/2).
#
# Each conditional's parameters are GBASS's own expressions, found in the bodies of gbass and nwbass and evaluated as written
# (source_assignment, source_calls). The reference is the log joint above, evaluated at eight points of the one coordinate with the
# rest held: the conditional is right exactly when the joint less its log kernel is constant over them, so each record carries the
# spread of that difference (in nats; rounding is 1e-12 or so). For w at beta != 0 the coordinate is u = w^(-1/2), the joint in u
# carrying the Jacobian 2 u^(-3). For gamma the points lie on each side of 0. For a, the code's mean and covariance against the
# Gaussian conditional by dense algebra. For the structure, the likelihood part of a birth's log acceptance as written, against the
# difference of log N(z | 0, w c V + w tau B Bᵀ) between the design with the last column and without it.
#
# Emits { r, pen, sum_r, rss_over_v, quad_term, shape_shift, v: {p, a, b, spread}, beta: {mean, sd, spread}, tau: {shape, rate,
#   spread}, w_gig: {p, a, b, spread}, w_mhn: {alpha, beta, gamma, spread}, gamma: {mean, sd, spread_positive, spread_negative,
#   mass_below_zero}, a: {mean_code, mean_dense, cov_max_diff}, birth: {code, dense} }.

source("oracles/_common.R")
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: gbass_conditionals_oracle.R <state.csv> <out.json> key=value ...")
S <- read_matrix(args[1])
kv <- strsplit(args[-(1:2)], "=", fixed = TRUE)
par <- setNames(lapply(kv, function(p) as.numeric(strsplit(p[2], ";", fixed = TRUE)[[1]])), vapply(kv, `[`, "", 1))
need <- c("a", "w", "beta", "tau", "gamma", "scale", "w_p", "w_a", "w_b", "v_p", "v_a", "v_b", "a_tau", "b_tau", "m_beta", "s_beta",
          "m_gamma", "s_gamma", "lower_bound")
if (length(setdiff(need, names(par)))) stop("missing keys: ", paste(setdiff(need, names(par)), collapse = ", "))

ns <- asNamespace("GBASS")
y <- S[, 1]; v <- S[, 2]; B <- S[, -(1:2), drop = FALSE]
N <- nrow(B); M <- ncol(B) - 1L
a <- par$a; w <- par$w; bet <- par$beta; tau <- par$tau; gam <- par$gamma; scale <- par$scale
if (length(a) != M + 1L) stop("a holds ", length(a), " values for ", M + 1L, " columns")

# ---- the stated joint ------------------------------------------------------------------------------------------------------------
lgig <- function(x, p, A, Bp) (p - 1) * log(x) - (A * x + Bp / x) / 2   # unnormalized; normalizers enter where a coordinate moves them
log_joint <- function(a = par$a, v = S[, 2], w = par$w, bet = par$beta, tau = par$tau, gam = par$gamma, nw = FALSE) {
  mu <- as.vector(B %*% a) + sqrt(w) * bet * v
  ll <- sum(dnorm(y, mu, sqrt(w * scale * v), log = TRUE))
  la <- sum(dnorm(a, 0, sqrt(w * tau), log = TRUE))
  lv <- if (nw) sum(abs(gam) - 0.5 * log(2 * pi) + lgig(v, -1/2, gam^2, 1)) else sum(lgig(v, par$v_p, par$v_a, par$v_b))
  lw <- lgig(w, par$w_p, par$w_a, par$w_b)
  lt <- lgig(tau, -par$a_tau, 0, 2 * par$b_tau)
  lb <- if (par$s_beta > 0) dnorm(bet, par$m_beta, par$s_beta, log = TRUE) else 0
  lg <- if (nw && par$s_gamma > 0) dnorm(gam, par$m_gamma, par$s_gamma, log = TRUE) else 0
  ll + la + lv + lw + lt + lb + lg
}
spread <- function(points, joint_at, kernel_at) {
  d <- vapply(points, function(x) joint_at(x) - kernel_at(x), 0)
  max(d) - min(d)
}

# ---- GBASS's expressions, as written ------------------------------------------------------------------------------------------
gb <- ns$gbass
env <- new.env(parent = ns)
assign("v_prior", list(type = "GIG", p = par$v_p, a = par$v_a, b = par$v_b), env)
assign("w_prior", list(type = "GIG", p = par$w_p, a = par$w_a, b = par$w_b, lower_bound = par$lower_bound), env)
for (k in c("N", "M", "B", "y", "v", "a", "w", "bet", "tau", "scale")) assign(k, get(k), env)
for (k in c("a_tau", "b_tau", "m_beta", "s_beta")) assign(k, par[[k]], env)
ev <- function(e) eval(e, env)
assign_as_written <- function(name, n, fn = gb) assign(name, ev(source_assignment(fn, name, n)), env)
# The sweep's own order after a's draw: yhat, r, pen, then rss_over_v, sum_r, quad_term and shape_shift.
assign_as_written("yhat", 1); assign_as_written("r", 1); assign_as_written("pen", 1)
for (k in c("rss_over_v", "sum_r", "quad_term", "shape_shift")) assign_as_written(k, 1)
call_args <- function(fn, head) {
  found <- source_calls(fn, function(e) identical(e[[1]], as.name(head)))
  if (!length(found)) stop("no call to ", head)
  found[[1]]
}
r <- get("r", env)

# v: the rgig2.vec call's arguments.
vc <- call_args(gb, "rgig2.vec")
v_par <- list(p = ev(vc$p), a = ev(vc$a), b = ev(vc$b))
i1 <- which.max(abs(r))   # the observation with the largest residual, where a wrong b would show most
v_pts <- v[i1] * exp(seq(-1.5, 1.5, length.out = 8))
v_rec <- c(v_par, list(
  spread = spread(v_pts, function(x) { vv <- v; vv[i1] <- x; log_joint(v = vv) },
                  function(x) lgig(x, v_par$p, v_par$a, v_par$b[i1]))))

# beta: the second assignment, rnorm(1, mean, sd).
bc <- source_assignment(gb, "bet", 2)
b_mean <- ev(bc[[3]]); b_sd <- ev(bc[[4]])
b_pts <- bet + seq(-2, 2, length.out = 8)
b_rec <- list(mean = b_mean, sd = b_sd,
              spread = if (par$s_beta > 0) spread(b_pts, function(x) log_joint(bet = x), function(x) dnorm(x, b_mean, b_sd, log = TRUE)) else NA_real_)

# tau: 1/rgamma(1, shape, rate).
tc <- source_assignment(gb, "tau", 2)[[3]]
t_shape <- ev(tc[[3]]); t_rate <- ev(tc[[4]])
t_pts <- tau * exp(seq(-1, 1, length.out = 8))
t_rec <- list(shape = t_shape, rate = t_rate,
              spread = spread(t_pts, function(x) log_joint(tau = x), function(x) lgig(x, -t_shape, 0, 2 * t_rate)))

# w at beta = 0: rgig2(p, a, b), its first call in the w block; the joint at beta = 0.
wc <- source_calls(gb, function(e) identical(e[[1]], as.name("rgig2")))[[1]]
wg <- list(p = ev(wc$p), a = ev(wc$a), b = ev(wc$b))
w_pts <- w * exp(seq(-1, 1, length.out = 8))
wg_rec <- c(wg, list(spread = spread(w_pts, function(x) log_joint(w = x, bet = 0), function(x) lgig(x, wg$p, wg$a, wg$b))))

# w at beta != 0: rMHN(n, alpha, beta, gamma) on u = w^(-1/2).
mc <- call_args(gb, "rMHN")
wm <- list(alpha = ev(mc$alpha), beta = ev(mc$beta), gamma = ev(mc$gamma))
u0 <- 1 / sqrt(w); u_pts <- u0 * exp(seq(-0.5, 0.5, length.out = 8))
wm_rec <- c(wm, list(spread = spread(u_pts, function(u) log_joint(w = u^-2) + log(2) - 3 * log(u),
                                     function(u) (wm$alpha - 1) * log(u) - wm$beta * u^2 + wm$gamma * u)))

# gamma (nwbass): the second assignment, rnorm(1, mean, sd), with N and sum(v) as nwbass holds them.
gc <- source_assignment(ns$nwbass, "gam", 2)
genv <- new.env(parent = ns)
for (k in c("N", "v")) assign(k, get(k), genv)
assign("s_gamma", par$s_gamma, genv); assign("m_gamma", par$m_gamma, genv)
g_mean <- eval(gc[[3]], genv); g_sd <- eval(gc[[4]], genv)
g_kernel <- function(x) dnorm(x, g_mean, g_sd, log = TRUE)
g_joint <- function(x) log_joint(gam = x, nw = TRUE)
g_rec <- list(mean = g_mean, sd = g_sd,
              spread_positive = spread(abs(g_mean) * exp(seq(-0.4, 0.4, length.out = 8)) + 1e-3, g_joint, g_kernel),
              spread_negative = spread(-abs(g_mean) * exp(seq(-0.4, 0.4, length.out = 8)) - 1e-3, g_joint, g_kernel),
              mass_below_zero = pnorm(0, g_mean, g_sd))

# a: U and U2 as the sweep forms them after the latent draws, then TCP(U, U2) and scale * w * U Uᵀ.
assign("z", ev(source_assignment(gb, "z", 2)), env)
assign("Vinv", ev(source_assignment(gb, "Vinv", 2)), env)
assign("U", ev(source_assignment(gb, "U", 5)), env)
assign("U2", ev(source_assignment(gb, "U2", 5)), env)
a_code <- as.vector(ev(call("TCP", as.name("U"), as.name("U2"))))
a_cov_code <- as.matrix(scale * w * tcrossprod(get("U", env)))
z <- y - bet * sqrt(w) * v
P <- crossprod(B, B / v) / (w * scale) + diag(1 / (w * tau), M + 1L)
a_dense <- as.vector(solve(P, crossprod(B, z / v) / (w * scale)))
a_rec <- list(mean_code = a_code, mean_dense = a_dense, cov_max_diff = max(abs(a_cov_code - solve(P))), cov_scale = max(abs(solve(P))))

# A birth's likelihood part as written, from the design without the last column (U, U2) to the design with it (U_cand, U2_cand).
birth <- if (M >= 1) {
  benv <- new.env(parent = ns)
  Bp <- B[, -(M + 1L), drop = FALSE]
  for (k in c("scale", "tau", "w", "z")) assign(k, get(k, env), benv)
  assign("v", v, benv); assign("Vinv", Matrix::Diagonal(x = 1 / v), benv)
  assign("B", Bp, benv); assign("M", M - 1L, benv); assign("B_cand", B, benv)
  assign("U", eval(source_assignment(gb, "U", 1), benv), benv)
  assign("U2", eval(source_assignment(gb, "U2", 1), benv), benv)
  assign("U_cand", eval(source_assignment(gb, "U_cand", 1)[[2]], benv), benv)   # the tryCatch's first argument
  assign("U2_cand", eval(source_assignment(gb, "U2_cand", 1), benv), benv)
  lap <- source_assignment(gb, "log_accept_prob_B", 2)   # the first is the refusal's -Inf
  code <- as.numeric(eval(lap, benv)) + log(scale) / 2 - log(tau) / 2   # the second line's likelihood terms
  dens <- function(X) {
    C <- w * scale * diag(v) + w * tau * tcrossprod(X)
    as.numeric(-0.5 * determinant(C)$modulus - 0.5 * crossprod(z, solve(C, z)))
  }
  list(code = code, dense = dens(B) - dens(Bp))
} else NULL

emit_json(list(r = as.vector(r), pen = get("pen", env), sum_r = get("sum_r", env), rss_over_v = get("rss_over_v", env),
               quad_term = get("quad_term", env), shape_shift = get("shape_shift", env),
               v = v_rec, beta = b_rec, tau = t_rec, w_gig = wg_rec, w_mhn = wm_rec, gamma = g_rec, a = a_rec, birth = birth),
          args[2], null = "null", na = "null")
