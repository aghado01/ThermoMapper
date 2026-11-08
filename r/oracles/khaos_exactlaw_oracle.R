# khaos_exactlaw_oracle.R — a khaos sampler's structure law against the exact posterior on a term space small enough to list (E19's
# build brief, slice 0): the check that the corrected samplers (khaos_corrected, the readouts' corrected reference) target the
# posterior, with the samplers as written beside them.
#
# Usage (cwd = r/):  Rscript oracles/khaos_exactlaw_oracle.R <data.csv> <out.json> <prior> <variant> <nmcmc> <seed>
#   data: headerless CSV, four inputs in [0, 1] then y. prior: ridge or gprior. variant: corrected or written.
#
# The space: p = 4, degree 3, order 2, so A holds 30 multi-indices (generate_A), and structures of at most max_basis = 2 distinct terms,
# 466 in all. The degree stays above the order because khaos as written cannot draw a degree at q = d_max: sample(q:d_max, 1, prob)
# reads a single q as 1:q. Each structure's posterior weight, with lambda integrated against its Gamma(a, b) prior and the terms
# uniform on A:
#   ln p(s) = lgamma(a + M) - M ln(b + 1) - M ln|A|  (M! orderings of M distinct terms, each ordered list |A|^-M),
#   ln p(y | s): the ridge's collapsed marginal at tau2 = 10, sigma^2's prior improper and common to every structure (khaos's g1 =
#   g2 = 0): -(1/2) ln|B'B + I/tau2| - (M + 1) ln tau - (n/2) ln(y'y - bhat'(B'B + I/tau2) bhat); the g-prior's with g0^2 integrated
#   against its IG(a_g, b_g) prior by Simpson's rule in ln g0^2 over [10^-6, 10^16], split at the CDF's points, the integrand through
#   the pencil (pencil_conditional, checked against dense_marginal: pencil_check), sigma^2 as khaos defaults it (a_sigma = b_sigma = 0).
# The basis is khaos's own make_basis. The chain: the variant at these settings, nmcmc iterations from set.seed(seed), the first tenth
# discarded, the g-prior's g0^2 by its "mh" branch (as written khaos's independence step against dgsq_full; corrected, one slice step on
# ln g0^2 against the derived conditional).
#
# Emits { prior, variant, terms (A's rows), exact: {M (P(M = 0..2)), structures, inclusion (each term's), g0sq_cuts and g0sq_cdf
#   (P(g0^2 < 300, 1000, 10^4), g-prior), pencil_check (g-prior)}, chain: {M, inclusion, g0sq_cdf, each with its Monte Carlo standard error by batch means over
#   30 batches (M_se, inclusion_se, g0sq_cdf_se), outside (the share of draws holding a term outside A), draws, seconds, accept,
#   propose, ess_M, ess_log_g0sq and g0sq_moves (the share of iterations g0^2 moves)} }.

source("oracles/khaos_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("usage: khaos_exactlaw_oracle.R <data.csv> <out.json> <prior> <variant> <nmcmc> <seed>")
D <- read_matrix(args[1]); prior <- args[3]; variant <- args[4]
nmcmc <- as.integer(args[5]); seed <- as.integer(args[6])
if (!prior %in% c("ridge", "gprior") || !variant %in% c("corrected", "written")) stop("prior is ridge or gprior; variant corrected or written")
X <- D[, -ncol(D), drop = FALSE]; y <- D[, ncol(D)]
n <- length(y); p <- ncol(X)
if (p != 4) stop("the space is set for four inputs")
degree <- 3; order <- 2; max_basis <- 2; tau2 <- 10
A <- kns$generate_A(p, degree, order)
nA <- nrow(A)
key_of <- function(alpha) paste(alpha, collapse = "")
akeys <- apply(A, 1, key_of)
column <- function(alpha) { v <- which(alpha > 0); kns$make_basis(v, alpha[v], X) }
cols <- lapply(seq_len(nA), function(r) column(A[r, ]))
gm_of <- function(alpha) { q <- sum(alpha > 0); d <- sum(alpha); complexity_g(q, d, 1) }

fm <- formals(if (prior == "ridge") kns$adaptive_khaos_ridge else kns$adaptive_khaos_gprior)
a_lam <- eval(if (prior == "ridge") fm$h1 else fm$a_M)
b_lam <- eval(if (prior == "ridge") fm$h2 else fm$b_M, list(y = y))
a_g <- eval(fm$a_g); b_g <- eval(fm$b_g)

# g0^2's log conditional given a structure, with beta and sigma^2 integrated, through the pencil: with D = D(g)^-1 B'B D(g)^-1 = L L' and
# L^-1 B'B L^-T = V diag(lambda) V', the precision B'B + D/g0^2 is L (diag(lambda) + I/g0^2) L' in V's basis, so any g0^2 costs O(p).
# Written apart from the chain's own (khaos_corrected's kc_log_g0sq), and checked against the dense multivariate t (pencil_check), so
# the chain's g0^2 step is checked against an implementation it does not share; vectorized over g0^2.
pencil_conditional <- function(gm, BtB, v, ssy, a_sigma, b_sigma) {
  D <- BtB / tcrossprod(gm)
  L <- t(chol(D)); Li <- forwardsolve(L, diag(nrow(D)))
  e <- eigen(Li %*% BtB %*% t(Li), symmetric = TRUE)
  lam <- e$values; w2 <- as.vector(crossprod(e$vectors, Li %*% v))^2
  ldD <- 2 * sum(log(diag(L))); k <- length(gm)
  function(theta) {
    shifted <- outer(1 / theta, lam, "+")
    quad <- ssy - as.vector(shifted^-1 %*% w2)
    -(a_g + 1) * log(theta) - b_g / theta - 0.5 * (ldD + rowSums(log(shifted))) + 0.5 * (ldD - k * log(theta)) -
      (a_sigma + n / 2) * log(b_sigma + 0.5 * quad)
  }
}
# Simpson's rule for the integral of exp(h(u) - top) over [a, b], 2000 panels.
simpson <- function(h, a, b, top) {
  u <- seq(a, b, length.out = 2001); f <- exp(h(u) - top); s <- (b - a) / 2000
  s / 3 * (f[1] + f[2001] + 4 * sum(f[seq(2, 2000, 2)]) + 2 * sum(f[seq(3, 1999, 2)]))
}
env <- corrected_env()
log_marginal <- function(set) {
  B <- cbind(rep(1, n), do.call(cbind, cols[set]))
  BtB <- crossprod(B); v <- crossprod(B, y); ssy <- sum(y^2)
  if (prior == "ridge") {
    Vinv <- BtB + diag(ncol(B)) / tau2
    bhat <- solve(Vinv, v)
    return(list(lm = -0.5 * as.numeric(determinant(Vinv, logarithm = TRUE)$modulus) - ncol(B) * 0.5 * log(tau2) -
                  (n / 2) * log(ssy - sum(bhat * v)), lg = NA_real_))
  }
  gm <- c(1, vapply(set, function(s) gm_of(A[s, ]), 0))
  cond <- pencil_conditional(gm, BtB, v, ssy, 0, 0)
  f <- function(u) cond(exp(u)) + u
  bounds <- c(log(1e-6), log(g0sq_cuts), log(1e16))
  top <- max(f(seq(bounds[1], bounds[length(bounds)], length.out = 2001)))
  parts <- vapply(seq_len(length(bounds) - 1), function(k) simpson(f, bounds[k], bounds[k + 1], top), 0)
  list(lm = top + log(sum(parts)), lg = cumsum(parts)[seq_along(g0sq_cuts)] / sum(parts))
}
# g0^2's posterior CDF at three points, bounded readings whose Monte Carlo error batch means estimates well; its mean in ln g0^2 has a
# long right tail when no term holds it (the intercept alone, under the vague IG(0.001, 1000)).
g0sq_cuts <- c(300, 1000, 10000)
sets <- c(list(integer(0)), unlist(lapply(1:max_basis, function(m) combn(nA, m, simplify = FALSE)), recursive = FALSE))
weights <- lapply(sets, function(s) {
  M <- length(s); mlk <- log_marginal(s)
  list(lw = lgamma(a_lam + M) - M * log(b_lam + 1) - M * log(nA) + mlk$lm, lg = mlk$lg)
})
lw <- vapply(weights, `[[`, 0, "lw"); post <- exp(lw - max(lw)); post <- post / sum(post)
# The pencil route against the dense multivariate t (dense_marginal) at five structures and three g0^2, under a proper sigma^2 prior
# IG(2, 1): the largest gap after the constants, which both routes must share.
pencil_check <- if (prior == "gprior") max(unlist(lapply(list(integer(0), 1L, 7L, c(2L, 11L), c(5L, 20L)), function(s) {
  B <- cbind(rep(1, n), do.call(cbind, cols[s])); gm <- c(1, vapply(s, function(t) gm_of(A[t, ]), 0))
  cond <- pencil_conditional(gm, crossprod(B), crossprod(B, y), sum(y^2), 2, 1)
  vapply(c(50, 500, 5000), function(th) abs(cond(th) + (a_g + 1) * log(th) + b_g / th + lgamma(2 + n / 2) - lgamma(2) - (n / 2) * log(2 * pi) -
                                              dense_marginal(y, B, gprior_s0(B, gm, th), 2, 1)), 0)
}))) else NULL
Ms <- lengths(sets)
exact <- list(M = vapply(0:max_basis, function(m) sum(post[Ms == m]), 0), structures = length(sets),
              inclusion = vapply(seq_len(nA), function(t) sum(post[vapply(sets, function(s) t %in% s, TRUE)]), 0),
              g0sq_cuts = if (prior == "gprior") g0sq_cuts else NULL,
              g0sq_cdf = if (prior == "gprior") colSums(post * do.call(rbind, lapply(weights, `[[`, "lg"))) else NULL,
              pencil_check = pencil_check)

sampler <- if (variant == "corrected") khaos_corrected(prior) else if (prior == "ridge") kns$adaptive_khaos_ridge else kns$adaptive_khaos_gprior
set.seed(seed)
t0 <- proc.time()[["elapsed"]]
extra <- if (prior == "ridge") list(tau2 = tau2) else list()
fit <- do.call(sampler, c(list(X = X, y = y, degree = degree, order = order, max_basis = max_basis, nmcmc = nmcmc, nburn = nmcmc %/% 10,
                                verbose = FALSE), extra))
secs <- proc.time()[["elapsed"]] - t0
draws <- length(fit$nbasis)
held_at <- matrix(FALSE, draws, nA); outside <- logical(draws)
for (i in seq_len(draws)) {
  for (m in seq_len(fit$nbasis[i])) {
    k <- fit$nint[i, m]; vv <- fit$vars[i, m, seq_len(k)]; dd <- fit$degs[i, m, seq_len(k)]
    if (anyDuplicated(vv)) { outside[i] <- TRUE; next }
    alpha <- integer(p); alpha[vv] <- dd
    hit <- match(key_of(alpha), akeys)
    if (is.na(hit)) outside[i] <- TRUE else held_at[i, hit] <- TRUE
  }
}
# A trace's Monte Carlo standard error by batch means over 30 batches, and its effective sample size.
batch_se <- function(x) { b <- 30; m <- length(x) %/% b; means <- vapply(seq_len(b), function(k) mean(x[(k - 1) * m + seq_len(m)]), 0)
  sqrt(var(means) / b) }
ess <- function(x) { v <- var(x); if (v == 0) NA_real_ else v / batch_se(x)^2 }
indicators <- lapply(0:max_basis, function(m) as.numeric(fit$nbasis == m))
lg <- if (prior == "gprior") log(fit$g2) else NULL
below <- if (prior == "gprior") lapply(g0sq_cuts, function(cut) as.numeric(fit$g2 < cut)) else NULL
chain <- list(M = vapply(indicators, mean, 0), M_se = vapply(indicators, batch_se, 0),
              inclusion = colMeans(held_at), inclusion_se = apply(held_at, 2, function(x) batch_se(as.numeric(x))),
              g0sq_cdf = if (!is.null(below)) vapply(below, mean, 0) else NULL,
              g0sq_cdf_se = if (!is.null(below)) vapply(below, batch_se, 0) else NULL,
              outside = mean(outside), draws = draws, seconds = secs,
              accept = unname(fit$count_accept), propose = unname(fit$count_propose),
              ess_M = ess(fit$nbasis), g0sq_moves = if (!is.null(lg)) mean(diff(fit$g2) != 0) else NULL,
              ess_log_g0sq = if (!is.null(lg)) ess(lg) else NULL)

emit_json(list(prior = prior, variant = variant, n = n, terms = lapply(seq_len(nA), function(r) A[r, ]), exact = exact, chain = chain),
          args[2], null = "null", na = "null")
