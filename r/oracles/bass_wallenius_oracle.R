# bass_wallenius_oracle.R — reference for BASS's variable-set proposal probabilities (port brief, fitting segment,
# D01's birth proposal; F-O5).
#
# Usage (cwd = r/):  Rscript oracles/bass_wallenius_oracle.R <z.csv> <out.json> [orders] [max_subsets] [draws]
#   z: headerless CSV, one row: the p variable weights z_u (R's z.vec, w2 + appearances; normalization is irrelevant).
#   orders: comma-separated interaction orders K (default 1,2,3,4,5; each at most p).
#   max_subsets: subsets evaluated per order (default 2000). Below C(p, K) a seeded sample of distinct subsets is
#     taken, and the sum over subsets is not reported.
#   draws: if given, that many draws of R's own sampler, sample(1:p, K, prob = z, replace = FALSE) (genCandBasis),
#     for every order whose subsets are all enumerated.
#
# The reference is brute force: a subset S's probability is the sum, over the K! orders of its variables, of the
# product of sequential selections without replacement, prod_k z_{s_k} / sum_{u not yet drawn} z_u. Every term and
# every denominator is a sum of positive numbers, so the reference keeps its relative precision. Against it, per
# subset:
#   r_dmwnch: R's dmwnchBass (mcmc_funcs.R), inclusion-exclusion over the 2^K subsets of S of Wallenius' integral
#     int_0^1 prod_{i in S} (1 - t^alpha_i) dt, alpha_i = z_i / sum_{k not in S} z_k. R calls it only for K >= 2
#     (logProbChangeMod's K = 1 branch scores the uniform 1/p it samples); at K = 1 its loop runs over 1:0 and
#     returns twice the probability, which is recorded here as unreachable.
#   quadrature: the same integral as a quadrature of its positive integrand, D01's proposal above K = 3, on
#     t = exp(-u): int_0^Inf prod_i (-expm1(-alpha_i u)) exp(-u) du.
# With equal weights every subset of size K has probability 1 / C(p, K), an exact check of the brute force itself.
#
# Emits { p, z, orders: [ { K, n_subsets, enumerated, subsets, brute, r_dmwnch, quadrature, max_rel_r, worst_r,
#         max_rel_quadrature, sum_brute, equal_weights, max_rel_equal, draws, max_z_score } ] }; sum_brute, max_rel_equal,
#         and the draw fields are null where they do not apply.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_wallenius_oracle.R <z.csv> <out.json> [orders] [max_subsets] [draws]")
z <- as.numeric(read_matrix(args[1])[1, ])
orders <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else 1:5
max_subsets <- if (length(args) >= 4) as.integer(args[4]) else 2000L
n_draws <- if (length(args) >= 5) as.integer(args[5]) else 0L
p <- length(z)
if (any(z <= 0)) stop("weights must be positive")
if (any(orders > p)) stop("an order exceeds p")
bass <- bass_source("mcmc_funcs.R")

permutations <- function(v) {
  if (length(v) <= 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) for (rest in permutations(v[-i])) out[[length(out) + 1]] <- c(v[i], rest)
  out
}

# Each denominator sums the weights still in the urn, rather than subtracting the drawn ones from the total, which
# cancels when a large weight is drawn before small ones.
brute_force <- function(S) sum(vapply(permutations(S), function(o) {
  left <- rep(TRUE, p); pr <- 1
  for (u in o) { pr <- pr * z[u] / sum(z[left]); left[u] <- FALSE }
  pr
}, 0))

quadrature <- function(S) {
  alpha <- z[S] / sum(z[-S])
  integrate(function(u) vapply(u, function(ui) prod(-expm1(-alpha * ui)), 0) * exp(-u), 0, Inf,
            rel.tol = 1e-13, abs.tol = 0, subdivisions = 1000L)$value
}

subsets_of <- function(K) {
  n <- choose(p, K)
  if (n <= max_subsets) return(list(sets = t(combn(p, K)), enumerated = TRUE))
  set.seed(K)
  sets <- unique(t(replicate(3L * max_subsets, sort(sample.int(p, K)))))
  list(sets = sets[seq_len(min(max_subsets, nrow(sets))), , drop = FALSE], enumerated = FALSE)
}

equal <- all(z == z[1])
results <- lapply(orders, function(K) {
  sub <- subsets_of(K)
  sets <- sub$sets
  brute <- apply(sets, 1, brute_force)
  r <- apply(sets, 1, function(S) bass$dmwnchBass(z, S))
  quad <- apply(sets, 1, quadrature)
  rel_r <- abs(r / brute - 1)
  res <- list(K = K, n_subsets = nrow(sets), enumerated = sub$enumerated, subsets = sets, brute = brute,
              r_dmwnch = r, quadrature = quad, max_rel_r = max(rel_r), worst_r = sets[which.max(rel_r), ],
              max_rel_quadrature = max(abs(quad / brute - 1)),
              sum_brute = if (sub$enumerated) sum(brute) else NA,
              equal_weights = equal,
              max_rel_equal = if (equal) max(abs(brute * choose(p, K) - 1)) else NA,
              draws = NA, max_z_score = NA)
  if (n_draws > 0 && sub$enumerated) {
    set.seed(100 + K)
    key <- function(S) paste(sort(S), collapse = ",")
    counts <- table(factor(replicate(n_draws, key(sample(1:p, size = K, prob = z, replace = FALSE))),
                           levels = apply(sets, 1, key)))
    freq <- as.numeric(counts) / n_draws
    res$draws <- n_draws
    res$max_z_score <- max(abs(freq - brute) / sqrt(brute * (1 - brute) / n_draws))
  }
  res
})

emit_json(list(p = p, z = z, orders = results), args[2], na = "null")
