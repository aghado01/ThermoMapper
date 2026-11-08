# khaos_defects_oracle.R — the brief's defects of khaos 2.1.1 that are not frozen-state evaluations, each run against the pinned package
# beside a reference that shares nothing with its code (E19's build brief, slice 0). The frozen-state rows (the g-prior's structure
# ratio at a general design, g0^2's conditional, the birth's forward probability, the variable swap's constructed state) are
# khaos_conditionals_oracle.R's.
#
# Usage (cwd = r/):  Rscript oracles/khaos_defects_oracle.R <data.csv> <out.json> [nmcmc] [seed]
#   data: headerless CSV, inputs in [0, 1] then y in the last column, for the seeded chain (section swap_chain).
#
# Sections:
#   normalizer: on an exactly orthogonal design, Sylvester's 64 x 64 Hadamard matrix's first five columns (the intercept and four
#     terms) with y orthogonal to a sixth, adaptive_khaos_gprior's birth of the sixth as written, the code's per-term penalty
#     -llik.alpha, at g0^2 in (1, 10, 100, 1000, 10^4) and (q, d) in ((1, 1), (2, 5), (3, 12)) at zeta = 1; beside it the reference
#     -Delta dense_marginal (sigma^2 ~ IG(2, 1/2)), (1/2) ln n and (1/2) ln(1 + g0^2 g_m^2).
#   evidence: log_evidence_fast_core and log_evidence_full_core as written at the nested structures 0..4 of the same design and a
#     response that loads on each column (nu0 = 4, s0_sq = 1/4, so nu0 s0_sq / 2 = 1/2), beside dense_marginal and
#     (1/2) ln|S0^-1|; and the namespace's functions whose bodies call either core (its callers).
#   partition: random_partition's law (partition_law) at (d, q) in ((4, 2), (5, 3), (6, 3), (7, 4), (8, 3), (6, 1)), beside the
#     uniform 1 / C(d - 1, q - 1) and the ratios' 1 / C(d, q).
#   swap_chain: adaptive_khaos_gprior on the data (nmcmc, nmcmc / 2 burn-in, seed), the kept draws' terms that repeat a variable; the
#     run is khaos_gprior_probed, the sampler as written with a probe that draws nothing, so it also reads, before every move, how far
#     the current state's log marginal as the ratios read it (its Sigma_n and Q from the g0^2 of its last acceptance) is from its value
#     at the g0^2 the iteration holds (stale: the mean and greatest |gap|, and the shares of iterations over 0.1 and 1 nat).
#   normal: the death's P(q > order | q0) per q0 as written (term_q_j, its normal with the continuity correction) beside the
#     Poisson-binomial by enumeration, and P(q = 0 | q0), at four (eta, order) cases.
#   typos: the paper's slips in the code's own reading: legendre_poly(x, 1); the birth's coin draw (deparsed); the support of the
#     birth's degree draw as written (5000 draws at q = 2, degree 6, p = 3); the degree weights as written at s_d = 1 beside the
#     paper's d^(-s_d); the count in the degree mutation's ratio (deparsed).
#   degree_draw: the birth's degree draw as written at q = d_max = 2 (the error it raises) and the corrected draw's value.
#   directions: each factor's exponent in the code's lprop.alpha for the birth, the death and the degree mutation, as written: the
#     change of lprop.alpha when one factor is scaled, over the change of that factor's log.
#   floor: s2_lower's default in each sampler and the clamp its sigma^2 draw applies (deparsed).

source("oracles/khaos_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: khaos_defects_oracle.R <data.csv> <out.json> [nmcmc] [seed]")
D <- read_matrix(args[1])
nmcmc <- if (length(args) >= 3) as.integer(args[3]) else 1500L
seed <- if (length(args) >= 4) as.integer(args[4]) else 20261107L
gp <- kns$adaptive_khaos_gprior
rd <- kns$adaptive_khaos_ridge
ldet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)
cp <- khaos_coin_pars
weights_at <- function(eta, q0) kns$make_weights(eta, q0, cp[[3]], cp[[2]], cp[[4]])

# ---- normalizer ----------------------------------------------------------------------------------------------------------------
H <- matrix(1, 1, 1); for (r in 1:6) H <- rbind(cbind(H, H), cbind(H, -H))
n <- nrow(H); Bh <- H[, 1:6]
set.seed(seed)
yh <- as.numeric(1 + Bh[, 2:5] %*% c(0.7, -0.4, 0.25, 0.1) + rnorm(n, 0, 0.3))
yo <- yh - sum(yh * Bh[, 6]) / n * Bh[, 6]
a0 <- 2; b0 <- 0.5
penalty <- function(g0sq, qn, dn) {
  e <- code_env(gp, y = yo, n = n, ssy = sum(yo^2), a_sigma = a0, b_sigma = b0, zeta = 1, g2 = c(g0sq, g0sq), i = 2L, B.curr = Bh[, 1, drop = FALSE])
  for (nm in c("BtB.curr", "g_vec.curr", "G.curr", "Sigma.curr", "v.curr", "quad.curr", "Q.curr", "ldet.curr")) as_written(gp, e, nm)
  for (k in 1:5) {
    assign("qq", if (k < 5) 1 else qn, e); assign("dd", if (k < 5) 1 else dn, e); assign("B.cand", Bh[, 1:(k + 1), drop = FALSE], e)
    for (nm in c("g_vec.cand", "G.cand", "BtB.cand", "v.cand", "Sigma.cand", "quad.cand", "Q.cand", "ldet.cand", "llik.alpha")) as_written(gp, e, nm)
    if (k == 5) break
    for (nm in c("B", "BtB", "v", "G", "g_vec", "Sigma", "Q", "ldet")) assign(paste0(nm, ".curr"), get(paste0(nm, ".cand"), e), e)
  }
  gm <- complexity_g(qn, dn, 1)
  g5 <- c(1, rep(1, 4)); g6 <- c(g5, gm)
  ref <- dense_marginal(yo, Bh[, 1:5], gprior_s0(Bh[, 1:5], g5, g0sq), a0, b0) - dense_marginal(yo, Bh, gprior_s0(Bh, g6, g0sq), a0, b0)
  list(g0sq = g0sq, q = qn, d = dn, g_m = gm, code = -as.numeric(e$llik.alpha), reference = ref, half_log_n = 0.5 * log(n),
       stated = 0.5 * log(1 + g0sq * gm^2))
}
normalizer <- list()
for (g0sq in c(1, 10, 100, 1000, 1e4)) for (qd in list(c(1, 1), c(2, 5), c(3, 12))) normalizer[[length(normalizer) + 1]] <- penalty(g0sq, qd[1], qd[2])

# ---- evidence ------------------------------------------------------------------------------------------------------------------
gh <- c(1, complexity_g(c(1, 2, 1, 3), c(2, 4, 3, 5), 1))
evidence_at <- function(k, g0sq) {
  Bk <- Bh[, 1:(k + 1), drop = FALSE]; gk <- gh[1:(k + 1)]
  BtB <- crossprod(Bk); Bty <- crossprod(Bk, yh)
  S0 <- gprior_s0(Bk, gk, g0sq)
  list(k = k, fast = kns$log_evidence_fast_core(sum(yh^2), as.numeric(Bty), gk, g0sq, 2 * a0, b0 / a0, n),
       full = kns$log_evidence_full_core(BtB, solve(BtB), Bty, gk, g0sq, 2 * a0, b0 / a0, sum(yh^2), n),
       dense = dense_marginal(yh, Bk, S0, a0, b0), half_ldet_prior_precision = -0.5 * ldet(S0))
}
callers <- Filter(function(f) { fn <- get(f, envir = kns); is.function(fn) &&
  any(c("log_evidence", "log_evidence_fast_core", "log_evidence_full_core") %in% all.names(body(fn))) &&
  !f %in% c("log_evidence") }, ls(kns, all.names = TRUE))
evidence <- list(g0sq = 50, structures = lapply(0:4, evidence_at, g0sq = 50), callers = callers)

# ---- partition -----------------------------------------------------------------------------------------------------------------
partition <- lapply(list(c(4, 2), c(5, 3), c(6, 3), c(7, 4), c(8, 3), c(6, 1)), function(dq) {
  law <- partition_law(dq[1], dq[2])
  list(d = dq[1], q = dq[2], compositions = lapply(seq_len(nrow(law$compositions)), function(i) law$compositions[i, ]), law = law$prob,
       uniform = 1 / choose(dq[1] - 1, dq[2] - 1), code_count = 1 / choose(dq[1], dq[2]))
})

# ---- swap_chain ----------------------------------------------------------------------------------------------------------------
X <- D[, -ncol(D), drop = FALSE]; y <- D[, ncol(D)]
probed <- khaos_gprior_probed()
set.seed(seed)
t0 <- proc.time()[["elapsed"]]
fit <- probed$fit(X, y, nmcmc = nmcmc, nburn = nmcmc %/% 2, verbose = FALSE)
secs <- proc.time()[["elapsed"]] - t0
stale <- probed$stale()
terms <- repeated <- 0; draws_with <- 0
for (i in seq_along(fit$nbasis)) {
  any_rep <- FALSE
  for (m in seq_len(fit$nbasis[i])) {
    vv <- fit$vars[i, m, seq_len(fit$nint[i, m])]
    terms <- terms + 1
    if (anyDuplicated(vv)) { repeated <- repeated + 1; any_rep <- TRUE }
  }
  if (any_rep) draws_with <- draws_with + 1
}
swap_chain <- list(nmcmc = nmcmc, seed = seed, seconds = secs, terms = terms, repeated = repeated, draws = length(fit$nbasis),
                   draws_with_repeat = draws_with, swaps_proposed = fit$count_propose[[4]], swaps_accepted = fit$count_accept[[4]],
                   stale = list(iterations = length(stale), mean_abs = mean(abs(stale)), greatest = max(abs(stale)),
                                over_tenth = mean(abs(stale) > 0.1), over_one = mean(abs(stale) > 1)))

# ---- normal --------------------------------------------------------------------------------------------------------------------
term_q <- source_assignment(rd, "term_q_j", 2)
normal <- lapply(list(list(c(3, 0, 1, 0), 2), list(c(3, 0, 1, 0), 3), list(c(5, 0, 0, 2, 1, 0), 3), list(c(4, 4, 1, 0, 0, 0, 0, 0, 0, 0), 2)), function(cs) {
  eta <- cs[[1]]; order <- cs[[2]]; p <- length(eta)
  J <- cp[[1]](seq_len(order)); Hs <- sum(J)
  per <- lapply(seq_len(order), function(q0) {
    w <- weights_at(eta, q0)
    e <- list2env(list(wts.j = w, jj = q0, J_probs = J, Hjs = Hs, order = order), parent = kns)
    as_written(rd, e, "mu_chi"); as_written(rd, e, "sig_chi")
    pb <- pb_enumerate(w)
    list(q0 = q0, weights = w, normal_over = exp(eval(term_q, e) - log(J[q0]) + log(Hs)), exact_over = sum(pb[-(1:(order + 1))]), zero = pb[1])
  })
  list(eta = eta, order = order, J = J / Hs, per_q0 = per)
})

# ---- typos ---------------------------------------------------------------------------------------------------------------------
xs <- c(-0.5, 0.25, 0.9)
te <- code_env(rd, nint.cand = 2, degree = 6, degree_penalty = 0, p = 3)
as_written(rd, te, "d_probs", 1); as_written(rd, te, "d_probs", 2)
degree_draw <- source_assignment(rd, "dtot.cand", 1)
support <- sort(unique(vapply(seq_len(5000), function(t) eval(degree_draw, te), 0)))
we <- code_env(rd, nint.cand = 2, degree = 6, degree_penalty = 1)
as_written(rd, we, "d_probs", 1); as_written(rd, we, "d_probs", 2)
typos <- list(p1 = list(x = xs, code = kns$legendre_poly(xs, 1)),
              coin_draw = paste(deparse(source_assignment(rd, "chi.cand", 2)), collapse = ""),
              degree_support = list(q = 2, degree = 6, p = 3, support = support),
              degree_weights = list(q = 2, degree = 6, s_d = 1, code = we$d_probs, paper = (2:6)^-1 / sum((2:6)^-1)),
              mutate_count = paste(deparse(source_assignment(rd, "lprop.alpha", 3)), collapse = ""))

# ---- degree_draw ---------------------------------------------------------------------------------------------------------------
# The degree draw as written where the term's order is the largest degree: sample(q:d_max, 1, prob = d_probs) with q = d_max hands
# sample a single number, which it reads as 1:q, so the one probability no longer fits. The corrected draw indexes the support.
dd_env <- code_env(rd, nint.cand = 2, degree = 2, degree_penalty = 0)
as_written(rd, dd_env, "d_probs", 1); as_written(rd, dd_env, "d_probs", 2)
degree_draw_section <- list(q = 2, degree = 2,
                            written = tryCatch({ eval(degree_draw, dd_env); "" }, error = function(e) conditionMessage(e)),
                            corrected = eval(quote((nint.cand:degree)[sample.int(length(d_probs), 1, prob = d_probs)]), dd_env))

# ---- directions ----------------------------------------------------------------------------------------------------------------
# lprop.alpha as written at a base point, and the exponent of one factor: (lprop with the factor moved - lprop at the base) / (the
# change of the factor's log). The base: an order-2 term of degree 4 (degree 6 for the mutation's candidate) under degree 8, so
# d_probs holds 7 equal weights; M = 3; delayed_reject_term ln 0.2; proposal_nint 0.2; proposal_reject 0.1. A weight is moved by
# doubling d_probs at the degree's index, 1 + d - q; a count by the degree moving up one, which leaves the equal weights alone.
lp <- function(k, moved = list()) {
  e <- code_env(rd, move_probs = rep(1/3, 3), i = 2L, nbasis = c(3, NA), d_probs = rep(1/7, 7), dtot.cand = if (k == 3) 6 else 4,
                nint.cand = 2, dtot.curr = 4, nint.curr = 2, delayed_reject_term = log(0.2), proposal_nint = 0.2, proposal_reject = 0.1)
  for (nm in names(moved)) assign(nm, moved[[nm]], e)
  as.numeric(eval(source_assignment(rd, "lprop.alpha", k), e))
}
doubled <- function(i) { v <- rep(1/7, 7); v[i] <- 2 / 7; v }
exponent <- function(k, moved, dlog) (lp(k, moved) - lp(k)) / dlog
up <- function(d) lchoose(d + 1, 2) - lchoose(d, 2)
directions <- list(
  birth = list(d_weight = exponent(1, list(d_probs = doubled(3)), log(2)), count = exponent(1, list(dtot.cand = 5), up(4)),
               terms = exponent(1, list(nbasis = c(4, NA)), log(5) - log(4)), coins = exponent(1, list(delayed_reject_term = log(0.4)), log(2))),
  death = list(d_weight = exponent(2, list(d_probs = doubled(3)), log(2)), count = exponent(2, list(dtot.cand = 5), up(4)),
               terms = exponent(2, list(nbasis = c(4, NA)), log(4) - log(3)), coins = exponent(2, list(proposal_nint = 0.4), log(2))),
  mutate_degree = list(d_weight_curr = exponent(3, list(d_probs = doubled(3)), log(2)),
                       d_weight_cand = exponent(3, list(d_probs = doubled(5)), log(2)),
                       count_curr = exponent(3, list(dtot.curr = 5), up(4)),
                       count_cand = exponent(3, list(dtot.cand = 7), up(6))))

# ---- floor ---------------------------------------------------------------------------------------------------------------------
clamp <- function(fn) paste(deparse(source_calls(fn, function(e) identical(e[[1]], as.name("max")) && "s2_lower" %in% all.names(e))[[1]]), collapse = "")
floor_section <- list(gprior_default = eval(formals(gp)$s2_lower), ridge_default = eval(formals(rd)$s2_lower),
                      gprior_clamp = clamp(gp), ridge_clamp = clamp(rd))

emit_json(list(normalizer = normalizer, evidence = evidence, partition = partition, swap_chain = swap_chain, normal = normal,
               typos = typos, degree_draw = degree_draw_section, directions = directions, floor = floor_section), args[2], null = "null", na = "null")
