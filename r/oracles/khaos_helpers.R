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
