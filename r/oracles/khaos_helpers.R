# khaos_helpers.R — what the khaos_*_oracle.R modes share (E19's build brief, slice 0): the pinned khaos 2.1.1 namespace, evaluation of
# its expressions as written, and the references that share nothing with its code: compositions and the Poisson-binomial by
# enumeration, the collapsed marginal as a dense multivariate t, and random_partition's law by enumerating its random draws.

source("oracles/_common.R")
suppressPackageStartupMessages(library(khaos))
kns <- asNamespace("khaos")
if (as.character(utils::packageVersion("khaos")) != "2.1.1") stop("khaos is not at the pinned 2.1.1")

# ---- khaos's own expressions, as written ---------------------------------------------------------------------------------------
# An environment whose parent is the namespace, holding the sampler's local safe_solve, so its expressions evaluate as they do there.
code_env <- function(fn, ...) {
  env <- new.env(parent = kns)
  assign("safe_solve", eval(source_assignment(fn, "safe_solve", 1), env), env)
  vals <- list(...)
  for (k in names(vals)) assign(k, vals[[k]], env)
  env
}
# The n-th assignment to name in fn's body, evaluated in env and bound there.
as_written <- function(fn, env, name, n = 1) assign(name, eval(source_assignment(fn, name, n), env), env)
# The block ({ ... }) of fn's body holding an assignment to `name`, as a list of statements in source order.
source_block <- function(fn, name) {
  blocks <- source_calls(fn, function(e) identical(e[[1]], as.name("{")) &&
    any(vapply(as.list(e)[-1], function(s) is.call(s) && as.character(s[[1]])[1] %in% c("<-", "=") && identical(s[[2]], as.name(name)), TRUE)))
  if (!length(blocks)) stop("no block assigns ", name)
  as.list(blocks[[1]])[-1]
}
# The calls `for (var in ...)` in fn's body, in source order.
source_for <- function(fn, var) source_calls(fn, function(e) identical(e[[1]], as.name("for")) && identical(e[[2]], as.name(var)))

# ---- references ------------------------------------------------------------------------------------------------------------------
# Every composition of d into q positive parts, one per row, in lexicographic order.
compositions <- function(d, q) {
  if (q == 1) return(matrix(d, 1, 1))
  do.call(rbind, lapply(seq_len(d - q + 1), function(a) cbind(a, compositions(d - a, q - 1))))
}

# random_partition's own law at (d, q), by enumerating its draw: sample(0:(d - q), q - 1, replace = TRUE) is uniform over the
# (d - q + 1)^(q - 1) ordered tuples, so the function is run as written once per tuple with sample answering that tuple.
partition_law <- function(d, q) {
  comps <- compositions(d, q)
  key <- function(v) paste(v, collapse = ",")
  prob <- setNames(numeric(nrow(comps)), apply(comps, 1, key))
  r <- d - q
  tuples <- if (q == 1) matrix(integer(0), 1, 0) else as.matrix(expand.grid(rep(list(0:r), q - 1)))
  f <- kns$random_partition
  env <- new.env(parent = kns)
  environment(f) <- env
  for (t in seq_len(nrow(tuples))) {
    tup <- unname(tuples[t, ])
    assign("sample", function(x, size, replace = FALSE, prob = NULL) tup, env)
    k <- key(f(d, q))
    if (!k %in% names(prob)) stop("random_partition(", d, ", ", q, ") gave ", k, ", not a composition")
    prob[k] <- prob[k] + 1 / nrow(tuples)
  }
  list(compositions = unname(comps), prob = unname(prob))
}

# P(chi | w) for an inclusion vector chi under independent coins w.
coin_prob <- function(chi, w) prod(ifelse(chi == 1, w, 1 - w))
# All inclusion vectors over p coins, one per row.
all_chi <- function(p) unname(as.matrix(expand.grid(rep(list(0:1), p))))
# The Poisson-binomial pmf of sum(chi) by enumeration over the 2^p vectors.
pb_enumerate <- function(w) {
  X <- all_chi(length(w))
  pr <- apply(X, 1, coin_prob, w = w)
  vapply(0:length(w), function(k) sum(pr[rowSums(X) == k]), 0)
}

# ln p(y | B, S0) with beta ~ N(0, sigma^2 S0), sigma^2 ~ IG(a, b): y is multivariate t with 2a degrees of freedom and scale
# (b/a)(I + B S0 B'), evaluated densely through the n x n Cholesky, with no Sigma_n, no g-form and no ridge algebra.
dense_marginal <- function(y, B, S0, a, b) {
  n <- length(y)
  V <- diag(n) + B %*% S0 %*% t(B)
  L <- chol(V)
  z <- backsolve(L, y, transpose = TRUE)
  lgamma(a + n / 2) - lgamma(a) + a * log(b) - (n / 2) * log(2 * pi) - sum(log(diag(L))) - (a + n / 2) * log(b + sum(z^2) / 2)
}

# The modified g-prior's prior covariance S0 = g0^2 D(g) (B'B)^(-1) D(g) (Rumsey, Francom et al., arXiv 2510.25036, the KHAOS model).
gprior_s0 <- function(B, g, g0sq) g0sq * diag(g, length(g)) %*% solve(crossprod(B)) %*% diag(g, length(g))
# g_m = (1 + q (d + q - 2))^(-zeta/2); the intercept's is 1.
complexity_g <- function(q, d, zeta) (1 + q * (d + q - 2))^(-zeta / 2)

# The pinned defaults the samplers take.
khaos_coin_pars <- eval(formals(kns$adaptive_khaos_ridge)$coin_pars)

# ---- khaos's samplers with slice 0's corrections ---------------------------------------------------------------------------------
# Each patch replaces one expression of the pinned body, matched by identical(), and states how many times it must match, so a
# change in the pinned source stops the oracle rather than leaving a patch unapplied: `patch`, `patch_body` and `replace_lang` in
# _common.R, sourced above.

# The helpers the corrected bodies call, in an environment whose parent is the namespace.
corrected_env <- function() {
  env <- new.env(parent = kns)
  # The move-type probability the sampler's own draw gives at M terms: a birth is forced at 0, and at max_basis the death and the
  # change share the draw.
  env$kc_move_prob <- function(type, M, move_probs, max_basis) {
    if (M == 0) return(c(1, 0, 0)[type])
    w <- if (M == max_basis) c(0, move_probs[2:3]) else move_probs
    (w / sum(w))[type]
  }
  # The Poisson-binomial pmf of the coins' sum, by dynamic programming over the p coins.
  env$kc_pb <- function(w) { pmf <- 1; for (x in w) pmf <- c(pmf * (1 - x), 0) + c(0, pmf * x); pmf }
  # ln of the birth's exact forward probability of the inclusion vector over vars, from term counts eta:
  # sum_q0 J(q0) P(chi | w(q0)) / P(1 <= q <= order | w(q0)), q0 drawn once and the coins redrawn until valid.
  env$kc_log_forward <- function(vars, eta, order, p, J_probs, coin_pars) {
    J <- J_probs / sum(J_probs); chi <- seq_len(p) %in% vars
    total <- 0
    for (q0 in seq_along(J)) {
      w <- kns$make_weights(eta, q0, coin_pars[[3]], coin_pars[[2]], coin_pars[[4]])
      pmf <- env$kc_pb(w)
      total <- total + J[q0] * exp(sum(log(ifelse(chi, w, 1 - w)))) / sum(pmf[2:(min(order, p) + 1)])
    }
    log(total)
  }
  # A composition of d into q positive parts, uniform over the C(d - 1, q - 1): q - 1 distinct cuts from {1, ..., d - 1}.
  env$kc_partition <- function(d, q) if (q == 1) d else diff(c(0, sort(sample.int(d - 1, q - 1)), d))
  # ln|S0^-1| under the modified g-prior, S0^-1 = D(g)^-1 B'B D(g)^-1 / g0^2.
  env$kc_ldet_prior <- function(BtB, g, g0sq) as.numeric(determinant(BtB, logarithm = TRUE)$modulus) - 2 * sum(log(g)) - length(g) * log(g0sq)
  # ln of g0^2's conditional with beta and sigma^2 integrated, up to a constant: the IG(a, b) prior, |Sigma_n|^(1/2) |S0|^(-1/2) and
  # the data term.
  # One slice step (Neal 2003: stepping out by w, at most m steps, then shrinkage) from x0 under the log density logf.
  env$kc_slice <- function(x0, logf, w = 1, m = 50) {
    level <- logf(x0) - stats::rexp(1)
    left <- x0 - w * stats::runif(1); right <- left + w
    j <- floor(m * stats::runif(1)); k <- m - 1 - j
    while (j > 0 && logf(left) > level) { left <- left - w; j <- j - 1 }
    while (k > 0 && logf(right) > level) { right <- right + w; k <- k - 1 }
    repeat {
      x <- left + (right - left) * stats::runif(1)
      if (logf(x) > level) return(x)
      if (x < x0) left <- x else right <- x
    }
  }
  # The density of ln g0^2 up to a constant, at one structure: the decomposition is made once, after which each g0^2 costs O(p).
  # With D(g)^-1 B'B D(g)^-1 = R'R (Cholesky) and R^-T B'B R^-1 = V diag(lambda) V', the posterior precision B'B + R'R/g0^2 has
  # ln-determinant ln|R'R| + sum ln(lambda + 1/g0^2), and v' (precision)^-1 v = sum w^2 / (lambda + 1/g0^2), w = V' R^-T v.
  env$kc_log_g0sq <- function(gm, BtB, v, ssy, a, b, a_sigma, b_sigma, n) {
    R <- chol(BtB / tcrossprod(gm))
    e <- eigen(backsolve(R, t(backsolve(R, BtB, transpose = TRUE)), transpose = TRUE), symmetric = TRUE)
    lam <- e$values; w2 <- as.vector(crossprod(e$vectors, backsolve(R, v, transpose = TRUE)))^2
    k <- length(gm)
    function(u) {
      inv <- exp(-u)
      quad <- ssy - sum(w2 / (lam + inv))
      -a * u - b * inv - 0.5 * sum(log(lam + inv)) - 0.5 * k * u - (a_sigma + n / 2) * log(b_sigma + 0.5 * quad)
    }
  }
  env
}

# The patches both samplers take: the birth's and death's proposal terms with the exact forward probability, uniform compositions
# counted by C(d - 1, q - 1) and the move probabilities the draw gives at the boundaries; the degree mutation's count; the degree draw
# where q = d_max (sample(x, 1, prob) reads a single x as 1:x); the death's prior ratio; the variable swap onto inactive variables
# only; the mutation type held at 1/2.
shared_patches <- list(
  patch(quote(lprop.alpha <- ((log(move_probs[2]) + log(1/(nbasis[i - 1] + 1))) - (log(move_probs[1]) + log(d_probs[1 + dtot.cand - nint.cand]) - lchoose(dtot.cand, nint.cand) + delayed_reject_term))),
        quote(lprop.alpha <- ((log(kc_move_prob(2, nbasis[i - 1] + 1, move_probs, max_basis)) + log(1/(nbasis[i - 1] + 1))) - (log(kc_move_prob(1, nbasis[i - 1], move_probs, max_basis)) + log(d_probs[1 + dtot.cand - nint.cand]) - lchoose(dtot.cand - 1, nint.cand - 1) + kc_log_forward(vars.cand, eta_vec, order, p, J_probs, coin_pars))))),
  patch(quote(lprop.alpha <- ((log(move_probs[1]) + log(proposal_nint) - log(1 - proposal_reject) + log(d_probs[1 + dtot.cand - nint.cand]) - lchoose(dtot.cand, nint.cand)) - (log(move_probs[2]) + log(1/nbasis[i - 1])))),
        quote(lprop.alpha <- ((log(kc_move_prob(1, nbasis[i - 1] - 1, move_probs, max_basis)) + kc_log_forward(vars.cand, eta.cand, order, p, J_probs, coin_pars) + log(d_probs[1 + dtot.cand - nint.cand]) - lchoose(dtot.cand - 1, nint.cand - 1)) - (log(kc_move_prob(2, nbasis[i - 1], move_probs, max_basis)) + log(1/nbasis[i - 1]))))),
  patch(quote(lprop.alpha <- log(d_probs[1 + dtot.curr - nint.curr]) - log(d_probs[1 + dtot.cand - nint.curr]) + (-lchoose(dtot.curr, nint.curr) + lchoose(dtot.cand, nint.curr))),
        quote(lprop.alpha <- log(d_probs[1 + dtot.curr - nint.curr]) - log(d_probs[1 + dtot.cand - nint.curr]) + (-lchoose(dtot.curr - 1, nint.curr - 1) + lchoose(dtot.cand - 1, nint.curr - 1)))),
  patch(quote(sample(nint.cand:degree, 1, prob = d_probs)), quote((nint.cand:degree)[sample.int(length(d_probs), 1, prob = d_probs)])),
  patch(quote(sample(nint.curr:degree, 1, prob = d_probs)), quote((nint.curr:degree)[sample.int(length(d_probs), 1, prob = d_probs)])),
  patch(quote(random_partition(dtot.cand, nint.cand)), quote(kc_partition(dtot.cand, nint.cand))),
  patch(quote(random_partition(dtot.cand, nint.curr)), quote(kc_partition(dtot.cand, nint.curr))),
  patch(quote(lprior.alpha <- (-log(lam[i - 1]) + log(nbasis[i - 1]) + log(A_total) - log(nbasis[i - 1] + 1))),
        quote(lprior.alpha <- (-log(lam[i - 1]) + log(nbasis[i - 1]) + log(A_total) - log(nbasis[i - 1])))),
  patch(quote(newvar_probs.curr[var_tochange] <- 0), quote(newvar_probs.curr[vars.curr] <- 0)),
  patch(quote(newvar_probs.cand[newvar.curr] <- 0), quote(newvar_probs.cand[vars.cand] <- 0)),
  patch(quote(mutate_prob <- mutate_eps + (1 - 2 * mutate_eps) * stats::pnorm(-diff(pmin(5, pmax(-5, stats::qnorm(success_rates)))))),
        quote(mutate_prob <- 0.5)))

# The g-prior's own: the structure ratios with (1/2) Delta ln|S0^-1| (the swap changes B'B, the degree mutation g and B'B); the current
# state's Sigma_n, Q and ln|Sigma_n| refreshed at the g0^2 the iteration holds; g0^2 against its derived conditional, in its default
# "mh" branch by one slice step on ln g0^2, the kernel E19 decided for the engine, in place of khaos's independence step, whose
# Laplace proposal is fitted to dgsq_full and so has lighter tails than the derived conditional (valid, but slow in the tail: the
# exact-law check read P(g0^2 < 10^4) 0.869 against 0.826 with it); sigma^2 drawn given the structure, then beta given sigma^2. The
# other g2_sample branches draw from Laplace approximations or target dgsq_full or dgsq_orth, so the corrected sampler takes only "mh"
# (the default) and "f" (g0^2 held).
gprior_mh_block <- as.call(c(as.name("{"), source_block(kns$adaptive_khaos_gprior, "lap_pars")))
gprior_patches <- list(
  patch(gprior_mh_block,
        quote({
          g2[i] <- exp(kc_slice(log(g2[i - 1]), kc_log_g0sq(g_vec.curr, BtB.curr, v.curr, ssy, a_g, b_g, a_sigma, b_sigma, n)))
        })),
  patch(quote(llik.alpha <- as.numeric(0.5 * (ldet.cand - ldet.curr) + (a_sigma + n/2) * (log(Q.curr) - log(Q.cand)))),
        quote(llik.alpha <- as.numeric(0.5 * (ldet.cand - ldet.curr) + (a_sigma + n/2) * (log(Q.curr) - log(Q.cand))) +
                0.5 * (kc_ldet_prior(BtB.cand, g_vec.cand, g2[i - 1]) - kc_ldet_prior(BtB.curr, g_vec.curr, g2[i - 1]))), n = 4),
  patch(quote(nbasis[i] <- nbasis[i - 1]),
        quote({
          nbasis[i] <- nbasis[i - 1]
          G.curr <- 1 + tcrossprod(1/g_vec.curr)/g2[i - 1]
          Sigma.curr <- solve(BtB.curr * G.curr)
          quad.curr <- ssy - t(v.curr) %*% Sigma.curr %*% v.curr
          Q.curr <- b_sigma + 0.5 * quad.curr
          ldet.curr <- determinant(Sigma.curr, logarithm = TRUE)$modulus
        })),
  patch(quote(mu_n <- Sigma.curr %*% v.curr),
        quote({
          s2[i] <- 1/stats::rgamma(1, n/2 + a_sigma, rate = as.numeric(Q.curr))
          mu_n <- Sigma.curr %*% v.curr
        })),
  patch(quote(beta_curr <- as.vector(mvtnorm::rmvnorm(1, mu_n, s2[i - 1] * Sigma.curr))),
        quote(beta_curr <- as.vector(mvtnorm::rmvnorm(1, mu_n, s2[i] * Sigma.curr)))),
  patch(quote(s2[i] <- max(s2_lower, 1/stats::rgamma(1, n/2 + a_sigma, rate = b_sigma + 0.5 * sum_sq[i]))), NULL))

# The ridge's own: sigma^2 drawn given the structure from the marginal its ratios read, IG(g1 + n/2, d/2), then beta given sigma^2.
ridge_patches <- list(
  patch(quote(Sigma_n <- s2[i - 1] * Lambda_i_n),
        quote({
          s2[i] <- 1/stats::rgamma(1, n/2 + g1, rate = as.numeric(d.curr)/2)
          Sigma_n <- s2[i] * Lambda_i_n
        })),
  patch(quote(s2[i] <- max(s2_lower, 1/stats::rgamma(1, n/2 + g1, rate = g2 + 0.5 * sum(resid^2)))), NULL))

# adaptive_khaos_ridge or adaptive_khaos_gprior with slice 0's corrections. The swap reaches only inactive variables, so where it can
# run (p > 3) the order must stay below p.
khaos_corrected <- function(prior) {
  fn <- if (prior == "ridge") kns$adaptive_khaos_ridge else kns$adaptive_khaos_gprior
  own <- if (prior == "ridge") ridge_patches else gprior_patches
  corrected <- patch_body(fn, c(shared_patches, own), corrected_env())
  function(X, y, ...) {
    order <- list(...)$order; if (is.null(order)) order <- eval(formals(fn)$order)
    if (ncol(X) > 3 && order >= ncol(X)) stop("the corrected swap needs order < p where p > 3")
    g2_sample <- list(...)$g2_sample
    if (!is.null(g2_sample) && !g2_sample %in% c("mh", "f")) stop("the corrected g-prior takes g2_sample \"mh\" or \"f\"")
    corrected(X, y, ...)
  }
}

# The g-prior as written with one statement added: before each move, the stale part of the current state's log marginal, as the
# structure ratio reads it (1/2 ln|Sigma_n| - (a + n/2) ln Q at the g0^2 of the state's last acceptance), against its value at the
# g0^2 the iteration holds. The probe draws nothing, so the chain is the one written.
khaos_gprior_probed <- function() {
  env <- corrected_env()
  env$kc_stale <- numeric(0)
  env$kc_record <- function(x) env$kc_stale <- c(env$kc_stale, x)
  fn <- patch_body(kns$adaptive_khaos_gprior, list(patch(quote(nbasis[i] <- nbasis[i - 1]), quote({
    nbasis[i] <- nbasis[i - 1]
    kc_record(local({
      Gf <- 1 + tcrossprod(1/g_vec.curr)/g2[i - 1]
      Sf <- solve(BtB.curr * Gf)
      qf <- ssy - t(v.curr) %*% Sf %*% v.curr
      lf <- determinant(Sf, logarithm = TRUE)$modulus
      as.numeric(0.5 * (ldet.curr - lf) - (a_sigma + n/2) * (log(Q.curr) - log(b_sigma + 0.5 * qf)))
    }))
  }))), env)
  list(fit = fn, stale = function() env$kc_stale)
}
