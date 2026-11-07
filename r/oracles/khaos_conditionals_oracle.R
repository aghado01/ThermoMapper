# khaos_conditionals_oracle.R — KHAOS's collapsed marginals, conditionals and proposal probabilities at frozen states, the pinned
# khaos 2.1.1's expressions as written beside references that share nothing with its code (E19's build brief, slice 0).
#
# Usage (cwd = r/):  Rscript oracles/khaos_conditionals_oracle.R <state.csv> <out.json> <proposals.csv> key=value ...
#   state: headerless CSV, one row per observation: y_i, then the basis row B_i0 .. B_iM, B_i0 = 1 the intercept, as khaos holds it
#     (inside the prior). The structures are nested: structure k holds columns 0..k, so a birth from k - 1 to k adds column k.
#   proposals: headerless CSV, one row per case, kind first, padded with NaN:
#     kind 1, (1, order, degree, eta_1 .. eta_p): a birth from a state whose term counts are eta, and the death back to it.
#     kind 2, (2, degree, s_d, q, d_curr, d_cand): a degree mutation of an order-q term, s_d the degree_penalty.
#     kind 3, (3, p, q, vars_1 .. vars_q, eta_1 .. eta_p): a variable swap in a term over vars, from a state with term counts eta.
#     kind 4, (4, K, lambda, A): the birth's and death's prior ratios for K = 1..K, and the move probabilities at M = 0.
#   keys (numbers; vectors ';'-separated): q, d (each column's order and degree, columns 1..M), zeta, a_sigma, b_sigma (sigma^2 ~
#     IG(a, b)), a_g, b_g (g0^2 ~ IG(a_g, b_g)), g0sq, tau2 (the ridge's prior variance), beta (M + 1 frozen coefficients), s2 (a frozen
#     sigma^2), theta (a grid of g0^2), s2grid (a grid of sigma^2), draws, seed.
#
# The model (Rumsey, Francom et al., arXiv 2510.25036, the KHAOS model and its coefficient priors): y ~ N(B beta, sigma^2 I);
# beta | sigma^2 ~ N(0, sigma^2 S0), S0 = g0^2 D(g) (B'B)^(-1) D(g) under the modified g-prior (g_0 = 1 for the intercept, g_m =
# (1 + q_m (d_m + q_m - 2))^(-zeta/2)), S0 = tau2 I under the ridge; sigma^2 ~ IG(a_sigma, b_sigma); g0^2 ~ IG(a_g, b_g). The ridge's
# code reads its sigma^2 prior as g1 = a_sigma and g2 = 2 b_sigma in its marginal, d = g2 + y'y - ..., and as g2 in its sigma^2 draw.
#
# Sections:
#   marginal: per structure k, ln p(y | k) under each prior by the closed form through Sigma_n = (B'B + S0^-1)^-1 and by
#     dense_marginal (y a multivariate t over n x n); half_ldet_prior_precision = (1/2) ln|S0^-1|. Per birth k - 1 -> k, the code's
#     llik.alpha as written (adaptive_khaos_gprior's and adaptive_khaos_ridge's, from their initial state through each accepted birth);
#     and the death's llik.alpha as written from the full structure, its last column removed.
#   g0sq: at the full structure, over theta: joint_integrated (dense ln p(y | theta) + ln IG(theta)), joint_given (ln N(y; B beta,
#     s2 I) + ln N(beta; 0, s2 S0(theta)) + ln IG(theta)), code (ln dgsq_full as written), data_term (-(a_sigma + n/2) ln Q(theta)),
#     closed_integrated (the closed form + ln IG); given_shape, given_rate (the conjugate IG given beta and sigma^2); laplace (the
#     IG proposal rg0sq_laplace_full fits, as written).
#   sigma2: at the full structure and frozen beta, over s2grid, each prior's joint ln N(y; B beta, s2 I) + ln N(beta; 0, s2 S0) +
#     ln IG(s2; a_sigma, b_sigma); the code's draw's IG(shape, rate) as written; the stated conditional's IG shape and rate.
#   birth: per kind-1 case, J (q0's weights), weights (make_weights' coins per q0), and every valid inclusion vector chi
#     (1 <= q <= order): exact (the forward probability of chi under the code's
#     own draw, sum_q0 J(q0)/H P(chi | w(q0)) / P(valid | w(q0)), validity by enumeration over all 2^p vectors), freq (the code's own
#     draw, run as written `draws` times), the code's birth term delayed_reject_term over the draws ending at chi (its least, mean
#     and greatest, and its value at chi's q0-marginal, no_redraw), and the death's reverse as written (ln proposal_nint -
#     ln(1 - proposal_reject)), from the state the birth leads to. Also the coins' validity per q0, exact and as the death writes it.
#   mutate_degree: per kind-2 case, the code's lprop.alpha as written; the log ratio under the code's own proposal (its degree weights
#     and random_partition's law) over every pair of compositions, least and greatest; and under uniform compositions.
#   mutate_vars: per kind-3 case, per position and new variable, the code's lprop.alpha as written, the log ratio of its own swap
#     kernel's reverse to forward, the forward probability, and whether the new variable is already in the term.
#   prior: per kind-4 case, per K, the code's birth log prior ratio at M = K - 1 and its death's at M = K, as written, and their sum;
#     the move-type draw as written at M = 0 and M = 1 (the shares of birth, death and change over `draws`), and the move parts of
#     the birth from M = 0 and the death to M = 0 as written (lprop.alpha with every other factor neutral).

source("oracles/khaos_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: khaos_conditionals_oracle.R <state.csv> <out.json> <proposals.csv> key=value ...")
S <- read_matrix(args[1])
P <- read_matrix(args[3])
kv <- strsplit(args[-(1:3)], "=", fixed = TRUE)
par <- setNames(lapply(kv, function(p) as.numeric(strsplit(p[2], ";", fixed = TRUE)[[1]])), vapply(kv, `[`, "", 1))
need <- c("q", "d", "zeta", "a_sigma", "b_sigma", "a_g", "b_g", "g0sq", "tau2", "beta", "s2", "theta", "s2grid", "draws", "seed")
if (length(setdiff(need, names(par)))) stop("missing keys: ", paste(setdiff(need, names(par)), collapse = ", "))

y <- S[, 1]; B <- S[, -1, drop = FALSE]
n <- length(y); M <- ncol(B) - 1L
a <- par$a_sigma; b <- par$b_sigma; g0sq <- par$g0sq; tau2 <- par$tau2
g <- c(1, complexity_g(par$q, par$d, par$zeta))
if (length(par$q) != M || length(par$d) != M || length(par$beta) != M + 1) stop("q, d and beta must match the ", M, " columns")
gp <- kns$adaptive_khaos_gprior
rd <- kns$adaptive_khaos_ridge
ldet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)
log_ig <- function(x, shape, rate) shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - rate / x

closed_marginal <- function(Bk, S0) {
  Sn <- solve(crossprod(Bk) + solve(S0)); v <- crossprod(Bk, y)
  quad <- as.numeric(sum(y^2) - t(v) %*% Sn %*% v)
  lgamma(a + n / 2) - lgamma(a) + a * log(b) - (n / 2) * log(2 * pi) + 0.5 * ldet(Sn) - 0.5 * ldet(S0) - (a + n / 2) * log(b + quad / 2)
}

# ---- marginal ------------------------------------------------------------------------------------------------------------------
structures <- lapply(0:M, function(k) {
  Bk <- B[, 1:(k + 1), drop = FALSE]
  Sg <- gprior_s0(Bk, g[1:(k + 1)], g0sq); Sr <- tau2 * diag(k + 1)
  list(k = k, gprior = list(closed = closed_marginal(Bk, Sg), dense = dense_marginal(y, Bk, Sg, a, b), half_ldet_prior_precision = -0.5 * ldet(Sg)),
       ridge = list(closed = closed_marginal(Bk, Sr), dense = dense_marginal(y, Bk, Sr, a, b), half_ldet_prior_precision = -0.5 * ldet(Sr)))
})

# The g-prior sampler's state from its initial assignments, then each birth k - 1 -> k as written, accepted.
eg <- code_env(gp, y = y, n = n, ssy = sum(y^2), a_sigma = a, b_sigma = b, zeta = par$zeta, g2 = c(g0sq, g0sq), i = 2L, B.curr = B[, 1, drop = FALSE])
for (nm in c("BtB.curr", "g_vec.curr", "G.curr", "Sigma.curr", "v.curr", "quad.curr", "Q.curr", "ldet.curr")) as_written(gp, eg, nm)
er <- code_env(rd, y = y, n = n, ssy = sum(y^2), g1 = a, g2 = 2 * b, tau2 = tau2, i = 2L, B.curr = B[, 1, drop = FALSE])
for (nm in c("BtB.curr", "BtBi.curr", "Vinv.curr", "bhat", "d.curr")) as_written(rd, er, nm)
births <- lapply(seq_len(M), function(k) {
  assign("qq", par$q[k], eg); assign("dd", par$d[k], eg); assign("B.cand", B[, 1:(k + 1), drop = FALSE], eg)
  for (nm in c("g_vec.cand", "G.cand", "BtB.cand", "v.cand", "Sigma.cand", "quad.cand", "Q.cand", "ldet.cand", "llik.alpha")) as_written(gp, eg, nm)
  lg <- as.numeric(eg$llik.alpha)
  for (nm in c("B", "BtB", "v", "G", "g_vec", "Sigma", "Q", "ldet")) assign(paste0(nm, ".curr"), get(paste0(nm, ".cand"), eg), eg)
  assign("nbasis", c(k - 1, NA), er); assign("B.cand", B[, 1:(k + 1), drop = FALSE], er)
  for (nm in c("BtB.cand", "BtBi.cand", "Vinv.cand", "bhat.cand", "d.cand", "llik.alpha")) as_written(rd, er, nm)
  lr <- as.numeric(er$llik.alpha)
  assign("B.curr", er$B.cand, er); assign("BtB.curr", er$BtB.cand, er); assign("BtBi.curr", er$BtBi.cand, er)
  assign("Vinv.curr", er$Vinv.cand, er); assign("bhat", er$bhat.cand, er); assign("d.curr", er$d.cand, er)
  list(k = k, gprior_code = lg, ridge_code = lr)
})
# The death of the last column, from the full structure, as written.
assign("tokill", M, eg)
for (nm in c("B.cand", "BtB.cand", "g_vec.cand", "G.cand", "v.cand", "Sigma.cand", "quad.cand", "Q.cand", "ldet.cand")) as_written(gp, eg, nm, 2)
as_written(gp, eg, "llik.alpha", 2)
assign("tokill", M, er); assign("nbasis", c(M, NA), er)
for (nm in c("B.cand", "BtB.cand", "BtBi.cand", "Vinv.cand", "bhat.cand", "d.cand", "llik.alpha")) as_written(rd, er, nm, 2)
marginal <- list(structures = structures, births = births,
                 death = list(gprior_code = as.numeric(eg$llik.alpha), ridge_code = as.numeric(er$llik.alpha)))

# ---- g0sq ----------------------------------------------------------------------------------------------------------------------
BtB <- crossprod(B); v <- crossprod(B, y); beta <- par$beta; s2 <- par$s2
resid <- as.numeric(y - B %*% beta)
Dinv <- diag(1 / g, M + 1)
prior_energy_g <- as.numeric(t(beta) %*% Dinv %*% BtB %*% Dinv %*% beta)   # beta' S0^-1 beta at g0^2 = 1
ll_y <- function(s2v) sum(dnorm(y, as.numeric(B %*% beta), sqrt(s2v), log = TRUE))
g0 <- lapply(par$theta, function(th) {
  S0 <- gprior_s0(B, g, th)
  Sn <- solve(BtB + solve(S0))
  Q <- b + 0.5 * as.numeric(sum(y^2) - t(v) %*% Sn %*% v)
  list(theta = th,
       joint_integrated = dense_marginal(y, B, S0, a, b) + log_ig(th, par$a_g, par$b_g),
       joint_given = ll_y(s2) + mvtnorm::dmvnorm(beta, rep(0, M + 1), s2 * S0, log = TRUE) + log_ig(th, par$a_g, par$b_g),
       code = log(kns$dgsq_full(th, par$a_g, par$b_g, g, BtB)),
       data_term = -(a + n / 2) * log(Q),
       closed_integrated = closed_marginal(B, S0) + log_ig(th, par$a_g, par$b_g))
})
lap <- kns$rg0sq_laplace_full(-1, par$a_g, par$b_g, g, BtB)
g0sq_section <- list(grid = g0, given_shape = par$a_g + (M + 1) / 2, given_rate = par$b_g + prior_energy_g / (2 * s2),
                     laplace = list(a = lap$a, b = lap$b), columns = M + 1)

# ---- sigma2 --------------------------------------------------------------------------------------------------------------------
rss <- sum(resid^2)
draw_args <- function(fn, marker, env) {
  calls <- source_calls(fn, function(e) identical(e[[1]], quote(stats::rgamma)) && grepl(marker, paste(deparse(e), collapse = "")))
  if (!length(calls)) stop("no sigma^2 draw marked ", marker)
  mc <- match.call(stats::rgamma, calls[[1]])
  list(shape = eval(mc$shape, env), rate = eval(mc$rate, env))
}
code_g <- draw_args(gp, "sum_sq", list2env(list(n = n, a_sigma = a, b_sigma = b, sum_sq = c(NA, rss), i = 2L)))
code_r <- draw_args(rd, "resid", list2env(list(n = n, g1 = a, g2 = 2 * b, resid = resid)))
Sg <- gprior_s0(B, g, g0sq); Sr <- tau2 * diag(M + 1)
sigma2_joint <- function(S0) vapply(par$s2grid, function(s) ll_y(s) + mvtnorm::dmvnorm(beta, rep(0, M + 1), s * S0, log = TRUE) + log_ig(s, a, b), 0)
sigma2 <- list(grid = par$s2grid, rss = rss,
               gprior = list(joint = sigma2_joint(Sg), code_shape = code_g$shape, code_rate = code_g$rate,
                             stated_shape = a + n / 2 + (M + 1) / 2, stated_rate = b + rss / 2 + prior_energy_g / (2 * g0sq)),
               ridge = list(joint = sigma2_joint(Sr), code_shape = code_r$shape, code_rate = code_r$rate,
                            stated_shape = a + n / 2 + (M + 1) / 2, stated_rate = b + rss / 2 + sum(beta^2) / (2 * tau2)))

# ---- proposals -----------------------------------------------------------------------------------------------------------------
set.seed(par$seed)
draws <- as.integer(par$draws)
cp <- khaos_coin_pars
weights_at <- function(eta, q0) kns$make_weights(eta, q0, cp[[3]], cp[[2]], cp[[4]])
birth_case <- function(order, degree, eta) {
  p <- length(eta)
  env <- code_env(rd, p = p, order = order, degree = degree, eta_vec = eta, coin_pars = cp)
  as_written(rd, env, "J_probs"); as_written(rd, env, "Hjs")
  J <- env$J_probs / env$Hjs
  X <- all_chi(p)
  valid <- rowSums(X) >= 1 & rowSums(X) <= order
  W <- lapply(seq_len(order), function(q0) weights_at(eta, q0))
  pv <- vapply(W, function(w) sum(apply(X[valid, , drop = FALSE], 1, coin_prob, w = w)), 0)
  marginal_chi <- function(chi) sum(vapply(seq_len(order), function(q0) J[q0] * coin_prob(chi, W[[q0]]), 0))
  exact <- vapply(which(valid), function(r) sum(vapply(seq_len(order), function(q0) J[q0] * coin_prob(X[r, ], W[[q0]]) / pv[q0], 0)), 0)
  # The code's own draw, as written: q0, the coins, and the redraw loop with its delayed_reject_term.
  rhs <- lapply(c("q0", "wts", "chi.cand", "delayed_reject_term"), function(nm) source_assignment(rd, nm))
  loop <- source_calls(rd, function(e) identical(e[[1]], as.name("while")))[[1]]
  key <- apply(X, 1, paste, collapse = "")
  hit <- character(draws); drt <- numeric(draws)
  for (t in seq_len(draws)) {
    for (j in seq_along(rhs)) assign(c("q0", "wts", "chi.cand", "delayed_reject_term")[j], eval(rhs[[j]], env), env)
    eval(loop, env)
    chi <- integer(p); chi[env$vars.cand] <- 1L
    hit[t] <- paste(chi, collapse = ""); drt[t] <- env$delayed_reject_term
  }
  # The death's reverse, as written, from the state the birth leads to: its coins read eta.cand, the counts with the term removed.
  death_loop <- source_for(rd, "jj")[[2]]
  death <- function(chi) {
    de <- code_env(rd, p = p, order = order, eta.cand = eta, vars.cand = which(chi == 1), J_probs = env$J_probs, Hjs = env$Hjs, coin_pars = cp)
    for (nm in c("proposal_nint", "proposal_nint_0", "proposal_nint_q")) assign(nm, 0, de)
    eval(death_loop, de)
    as_written(rd, de, "proposal_reject")
    list(log_reverse = log(de$proposal_nint) - log(1 - de$proposal_reject), nint_q = de$proposal_nint_q, nint_0 = de$proposal_nint_0)
  }
  per_chi <- lapply(seq_along(which(valid)), function(s) {
    r <- which(valid)[s]; chi <- X[r, ]
    mine <- drt[hit == key[r]]
    dth <- death(chi)
    list(chi = chi, exact = exact[s], freq = length(mine) / draws, n = length(mine),
         code_birth = if (length(mine)) list(least = min(mine), mean = mean(mine), greatest = max(mine)) else NULL,
         no_redraw = log(marginal_chi(chi)), death_reverse = dth$log_reverse)
  })
  # The coins' validity per q0: exact by enumeration, and the parts the death writes (q = 0 exactly, q > order by its normal).
  dth0 <- death(X[which(valid)[1], ])
  list(order = order, degree = degree, eta = eta, J = J, weights = W, valid_exact = pv, q0_marginal_invalid_exact = sum(J * (1 - pv)),
       q0_marginal_zero = dth0$nint_0, q0_marginal_over_normal = dth0$nint_q,
       q0_marginal_over_exact = sum(vapply(seq_len(order), function(q0) J[q0] * sum(pb_enumerate(W[[q0]])[-(1:(order + 1))]), 0)),
       chis = per_chi, draws = draws)
}

mutate_degree_case <- function(degree, s_d, q, d_curr, d_cand) {
  env <- code_env(rd, degree = degree, degree_penalty = s_d, nint.curr = q, dtot.curr = d_curr, dtot.cand = d_cand)
  as_written(rd, env, "d_probs", 5); as_written(rd, env, "d_probs", 6)
  as_written(rd, env, "lprop.alpha", 3)
  pd <- env$d_probs; at <- function(dd) pd[1 + dd - q]
  lc <- partition_law(d_curr, q); ln <- partition_law(d_cand, q)
  pairs <- outer(log(at(d_curr) * lc$prob), log(at(d_cand) * ln$prob), "-")
  list(degree = degree, s_d = s_d, q = q, d_curr = d_curr, d_cand = d_cand, code = as.numeric(env$lprop.alpha),
       actual_least = min(pairs), actual_greatest = max(pairs),
       uniform = log(at(d_curr)) - lchoose(d_curr - 1, q - 1) - log(at(d_cand)) + lchoose(d_cand - 1, q - 1), d_probs = pd)
}

mutate_vars_case <- function(p, vars, eta) {
  q <- length(vars)
  stmts <- source_block(rd, "ind_tochange")
  lhs <- function(s) if (is.call(s) && as.character(s[[1]])[1] %in% c("<-", "=")) deparse(s[[2]]) else ""
  from <- which(vapply(stmts, lhs, "") == "ind_tochange"); to <- which(vapply(stmts, lhs, "") == "newvar.cand")
  out <- list()
  for (ind in seq_len(q)) for (nv in seq_len(p)) {
    if (nv == vars[ind]) next
    env <- code_env(rd, p = p, eta_vec = eta, coin_pars = cp, vars.curr = vars, nint.curr = q)
    for (s in stmts[from:to]) {
      nm <- lhs(s)
      if (nm == "ind_tochange") assign(nm, ind, env) else if (nm == "newvar.curr") assign(nm, nv, env) else eval(s, env)
    }
    as_written(rd, env, "lprop.alpha", 4)
    fwd <- env$newvar_probs.curr[nv]; rev <- env$newvar_probs.cand[vars[ind]]
    out[[length(out) + 1]] <- list(position = ind, new = nv, code = as.numeric(env$lprop.alpha), exact = log(rev) - log(fwd),
                                   forward = fwd, duplicate = nv %in% vars[-ind], candidate = env$vars.cand)
  }
  list(p = p, vars = vars, eta = eta, swaps = out)
}

prior_case <- function(K, lam, A) {
  env <- code_env(rd, lam = c(lam, NA), A_total = A, i = 2L)
  ratios <- lapply(seq_len(K), function(k) {
    assign("nbasis", c(k - 1, NA), env); as_written(rd, env, "lprior.alpha", 1); lb <- env$lprior.alpha
    assign("nbasis", c(k, NA), env); as_written(rd, env, "lprior.alpha", 2); ld <- env$lprior.alpha
    list(K = k, birth = lb, death = ld, sum = lb + ld)
  })
  # The move-type draw, the loop's first three statements as written.
  head3 <- source_block(rd, "move.type")[1:3]
  share <- function(m) {
    me <- code_env(rd, nbasis = c(m, NA), i = 2L, move_probs = rep(1/3, 3), max_basis = 1000)
    tab <- table(factor(vapply(seq_len(draws), function(t) { for (s in head3) eval(s, me); me$move.type }, ""), levels = c("birth", "death", "change")))
    as.numeric(tab) / draws
  }
  # The move parts of the birth from M = 0 and the death from M = 1 as written, every other factor neutral (one degree, one
  # composition, the coins' terms 0).
  me <- code_env(rd, move_probs = rep(1/3, 3), i = 2L, d_probs = 1, dtot.cand = 1, nint.cand = 1, delayed_reject_term = 0,
                 proposal_nint = 1, proposal_reject = 0)
  assign("nbasis", c(0, NA), me); as_written(rd, me, "lprop.alpha", 1); lb <- me$lprop.alpha
  assign("nbasis", c(1, NA), me); as_written(rd, me, "lprop.alpha", 2); ld <- me$lprop.alpha
  list(ratios = ratios, moves_at_0 = share(0), moves_at_1 = share(1), birth_from_0_code = lb, death_to_0_code = ld)
}

birth <- mutate_degree <- mutate_vars <- prior <- list()
for (r in seq_len(nrow(P))) {
  row <- P[r, ]; row <- row[!is.nan(row)]
  k <- row[1]
  if (k == 1) birth[[length(birth) + 1]] <- birth_case(row[2], row[3], row[-(1:3)])
  else if (k == 2) mutate_degree[[length(mutate_degree) + 1]] <- mutate_degree_case(row[2], row[3], row[4], row[5], row[6])
  else if (k == 3) { p <- row[2]; q <- row[3]; mutate_vars[[length(mutate_vars) + 1]] <- mutate_vars_case(p, row[4:(3 + q)], row[(4 + q):(3 + q + p)]) }
  else if (k == 4) prior[[length(prior) + 1]] <- prior_case(row[2], row[3], row[4])
  else stop("unknown kind ", k)
}

emit_json(list(n = n, M = M, g = g, marginal = marginal, g0sq = g0sq_section, sigma2 = sigma2,
               birth = birth, mutate_degree = mutate_degree, mutate_vars = mutate_vars, prior = prior),
          args[2], null = "null", na = "null")
