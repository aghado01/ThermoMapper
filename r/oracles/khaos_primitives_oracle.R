# khaos_primitives_oracle.R — KHAOS's primitives in the pinned khaos 2.1.1, each beside a reference that shares nothing with its code
# (E19's build brief, slice 0).
#
# Usage (cwd = r/):  Rscript oracles/khaos_primitives_oracle.R <cases.csv> <out.json> [draws] [seed]
#   cases: headerless CSV, one row per case, kind first, padded with NaN to a common width:
#     kind 1, (1, jmax): ss_legendre_poly, the standardized shifted Legendre psi_j(x) = sqrt(2j + 1) P_j(2x - 1), for j = 0..jmax, at
#       the (jmax + 1)-point Gauss-Legendre nodes on [0, 1] and at 0, 1/3, 1/2, 0.9, 1. The references: orthonormality, the Gram
#       sum_k w_k psi_i(x_k) psi_j(x_k) under that rule (exact to degree 2 jmax + 1), and psi_jmax's zeros at the jmax-point nodes
#       (Golub-Welsch, an eigenproblem). legendre_poly recurses on both j - 1 and j - 2 without memory, so one call costs about
#       phi^j calls: seconds is psi_jmax's own time.
#     kind 2, (2, p, d, q): |A_{p,d,q}|, the multi-indices of order 1..q and degree 1..d: A_size as written; nrow(generate_A), the
#       package's enumerator, where A_size <= 5000; and a brute count over {0..d}^p where (d + 1)^p <= 2e6.
#     kind 3, (3, q0, eta_1 .. eta_p): make_weights(eta, q0) under the samplers' default coin_pars (epsilon 2, alpha 1, three passes),
#       the coin probabilities a birth flips, meant to sum to q0; beside them the exact solution of sum_j v_j^beta = q0 over the same
#       v_j = (eta_j + 2) / sum(eta + 2) * q0 / p, by uniroot.
#     kind 4, (4, d, q): random_partition(d, q)'s law, by enumerating its draw (partition_law), and `draws` calls of the function
#       itself from set.seed(seed); the compositions of d into q positive parts, each with its probability and its frequency.
#   draws (default 100000), seed (default 20261107).
#
# Emits { legendre: [...], a_size: [...], weights: [...], partition: [...] }, one record per case in input order within its kind:
#   legendre: jmax, x, values (rows j = 0..jmax over x), nodes (the Gram's rule on [0, 1]), gram_diag_error, gram_offdiag_max,
#     zeros_max (max |psi_jmax| at the jmax-point nodes, over sqrt(2 jmax + 1)), seconds.
#   a_size: p, d, q, code, generate (or null), brute (or null).
#   weights: q0, eta, code (the weights), code_sum, exact (the weights at the exact beta), beta_exact.
#   partition: d, q, compositions (rows), law, freq, n.

source("oracles/khaos_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: khaos_primitives_oracle.R <cases.csv> <out.json> [draws] [seed]")
cases <- read_matrix(args[1])
draws <- if (length(args) >= 3) as.integer(args[3]) else 100000L
seed <- if (length(args) >= 4) as.integer(args[4]) else 20261107L

legendre_case <- function(jmax) {
  rule <- gauss_legendre(jmax + 1, 0, 1)
  x <- c(rule$x, 0, 1/3, 1/2, 0.9, 1)
  values <- t(vapply(0:jmax, function(j) kns$ss_legendre_poly(x, j), numeric(length(x))))
  at_nodes <- values[, seq_along(rule$x), drop = FALSE]
  gram <- at_nodes %*% diag(rule$w, length(rule$w)) %*% t(at_nodes)
  zeros <- gauss_legendre(jmax, 0, 1)$x
  t0 <- proc.time()[["elapsed"]]
  last <- kns$ss_legendre_poly(zeros, jmax)
  secs <- proc.time()[["elapsed"]] - t0
  list(jmax = jmax, x = x, values = lapply(seq_len(nrow(values)), function(i) values[i, ]), nodes = rule$x,
       gram_diag_error = max(abs(diag(gram) - 1)), gram_offdiag_max = max(abs(gram[upper.tri(gram)])),
       zeros_max = max(abs(last)) / sqrt(2 * jmax + 1), seconds = secs)
}

a_size_case <- function(p, d, q) {
  code <- kns$A_size(p, d, q)
  generate <- if (code <= 5000) nrow(kns$generate_A(p, d, q)) else NULL
  brute <- if ((d + 1)^p <= 2e6) {
    G <- as.matrix(expand.grid(rep(list(0:d), p)))
    s <- rowSums(G); o <- rowSums(G > 0)
    sum(s >= 1 & s <= d & o <= q)
  } else NULL
  list(p = p, d = d, q = q, code = code, generate = generate, brute = brute)
}

weights_case <- function(q0, eta) {
  cp <- khaos_coin_pars
  w <- kns$make_weights(eta, q0, cp[[3]], cp[[2]], cp[[4]])
  p <- length(eta)
  v <- (eta^cp[[2]] + cp[[3]]) / sum(eta^cp[[2]] + cp[[3]]) * q0 / p
  beta <- if (q0 < p) uniroot(function(b) sum(v^b) - q0, c(1e-9, 1e3), tol = 1e-15)$root else 0
  list(q0 = q0, eta = eta, code = w, code_sum = sum(w), exact = v^beta, beta_exact = beta)
}

partition_case <- function(d, q) {
  law <- partition_law(d, q)
  set.seed(seed)
  keys <- apply(law$compositions, 1, paste, collapse = ",")
  hits <- table(factor(vapply(seq_len(draws), function(i) paste(kns$random_partition(d, q), collapse = ","), ""), levels = keys))
  list(d = d, q = q, compositions = lapply(seq_len(nrow(law$compositions)), function(i) law$compositions[i, ]), law = law$prob,
       freq = as.numeric(hits) / draws, n = draws)
}

legendre <- a_size <- weights <- partition <- list()
for (i in seq_len(nrow(cases))) {
  row <- cases[i, ]; row <- row[!is.nan(row)]
  k <- row[1]
  if (k == 1) legendre[[length(legendre) + 1]] <- legendre_case(row[2])
  else if (k == 2) a_size[[length(a_size) + 1]] <- a_size_case(row[2], row[3], row[4])
  else if (k == 3) weights[[length(weights) + 1]] <- weights_case(row[2], row[-(1:2)])
  else if (k == 4) partition[[length(partition) + 1]] <- partition_case(row[2], row[3])
  else stop("unknown kind ", k)
}
emit_json(list(legendre = legendre, a_size = a_size, weights = weights, partition = partition), args[2], null = "null", na = "null")
