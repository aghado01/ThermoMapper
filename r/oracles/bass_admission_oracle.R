# bass_admission_oracle.R — reference for BASS's numerical admission gate (port brief, D06; F-O2, F-O7).
#
# Usage (cwd = r/):  Rscript oracles/bass_admission_oracle.R <design.csv> <out.json> [mode] [H.csv]
#   design: headerless CSV; the columns of B in insertion order, the intercept first as in R BASS.
#   mode: uncentred (default; R's intercept convention), centred (D03's Liang default) or hadamard (the Liang
#     default for direct functional fits). In hadamard mode the design is A, the case-side factors (n x (M+1)), and
#     H.csv the output-side factors (m x (M+1)), each with the intercept's column of ones first; column j of the
#     design is then the Kronecker product of A_j and H_j, rows ordered (case, location) with the case outer.
#
# uncentred emits { mode, positive_definite, loo, min_loo, in_order, min_in_order, getqf_ratio, getqf_refuses }:
#   loo[j] = 1 / (G_jj (G^-1)_jj), column j's squared sine to the span of all the other columns: the reciprocal
#     of its uncentred variance inflation factor, independent of column order (D06's gate).
#   in_order[j] = d_j / G_jj, the Cholesky pivots in the given order over the diagonal (FactorStatus's
#     MinRelativePivot); each is at least loo[j], and the last equals it.
#   getqf_ratio = max_{j >= 2} r_jj / min_j r_jj for R = chol(G), and getqf_refuses = getqf_ratio > 1e3: R BASS's
#     admission rule in getQf (mcmc_funcs.R), which depends on column norms; a defect fixture, not ported (D06).
#
# centred emits { mode, positive_definite, rho, loo, min_loo, in_order, min_in_order, loo_uncentred,
#                 loo_explicit, explicit_rel_error, loo_explicit_recursive, explicit_recursive_rel_error,
#                 loo_downdate, downdate_rel_error } over the columns after the intercept, which the Liang convention
#   keeps out of G:
#   loo[j] = column j's squared sine to the span of the others and 1, over its centred squared norm: the reciprocal
#     of its usual (centred) variance inflation factor. The reference, by QR of the explicitly centred design.
#   rho[j] = |b_j - mean(b_j)|^2 / |b_j|^2, and loo_uncentred[j] the same column's uncentred angle with the
#     intercept a column, also by QR, so that loo_uncentred = loo * rho (D06).
#   in_order: the pivots of the explicitly centred Gram in the given order, over its diagonal.
#   loo_explicit, loo_downdate: s from an equilibrated Cholesky of the centred Gram formed two ways, from explicitly
#     centred columns and by D03's rank-one downdate B'B - N bbar bbar', each cross-product summed pairwise;
#     *_rel_error = the largest relative difference from loo. The downdate loses about eps / rho_j on near-constant
#     columns, whatever the summation.
#   loo_explicit_recursive: explicitly centred columns again, summed left to right in double. On a column with many
#     equal centred values after one large one, the rounding accumulates systematically, so this depends on the
#     row order.
#
# hadamard emits { mode, n, m, positive_definite, rho, loo, min_loo, loo_three_term, three_term_rel_error,
#                  loo_downdate, downdate_rel_error, gram_three_term_error, gram_downdate_error }:
#   loo, rho: as in centred mode, from the materialized design (fixtures only; the port never forms it).
#   loo_three_term: s from the centred Hadamard form, with A_c, H_c the centred factors and abar, hbar their means
#     (cross-products summed pairwise throughout),
#     G_c = (A_c'A_c) o (H_c'H_c) + n (abar abar') o (H_c'H_c) + m (hbar hbar') o (A_c'A_c);
#   loo_downdate: s from the downdate (A'A) o (H'H) - n m bbar bbar', with bbar_j = abar_j hbar_j.
#   gram_*_error: the largest entry difference of each Gram from the explicitly centred one, over its largest entry.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_admission_oracle.R <design.csv> <out.json> [uncentred|centred|hadamard] [H.csv]")
mode <- if (length(args) >= 3) args[3] else "uncentred"
if (!mode %in% c("uncentred", "centred", "hadamard")) stop("unknown mode: ", mode)

# Two-pass centring: subtract the mean, then the mean of what is left.
centre <- function(X) {
  Xc <- sweep(X, 2, colMeans(X))
  sweep(Xc, 2, colMeans(Xc))
}
# Order-free angles by QR: each column's residual against the others, over its squared norm. tol = 0 keeps R's
# LINPACK QR from dropping nearly dependent columns.
loo_qr <- function(X) {
  vapply(seq_len(ncol(X)), function(j) {
    sum(qr.resid(qr(X[, -j, drop = FALSE], tol = 0), X[, j])^2) / sum(X[, j]^2)
  }, numeric(1))
}
# Order-free angles from a Gram, by an equilibrated Cholesky (the port's route); NA when it is not positive definite.
loo_gram <- function(G) {
  d <- 1 / sqrt(diag(G))
  R <- tryCatch(chol(G * tcrossprod(d)), error = function(e) NULL)
  if (is.null(R)) rep(NA_real_, ncol(G)) else 1 / diag(chol2inv(R))
}
rel_error <- function(s, ref) max(abs(s / ref - 1))
# Cross-products X'Y in plain double arithmetic, so that the results do not depend on the platform's long double:
# pairwise (a balanced tree over the rows) or recursive (left to right).
pairwise_sum <- function(v) {
  while (length(v) > 1) {
    if (length(v) %% 2 == 1) v <- c(v, 0)
    v <- v[c(TRUE, FALSE)] + v[c(FALSE, TRUE)]
  }
  v
}
recursive_sum <- function(v) Reduce(`+`, v, accumulate = FALSE)
cross <- function(X, Y = X, summation = pairwise_sum) {
  G <- matrix(0, ncol(X), ncol(Y))
  for (i in seq_len(ncol(X))) for (j in seq_len(ncol(Y))) G[i, j] <- summation(X[, i] * Y[, j])
  G
}
strip_intercept <- function(X, what) {
  if (ncol(X) < 2 || any(X[, 1] != 1)) stop(what, ": the first column must be the intercept's ones")
  X[, -1, drop = FALSE]
}

B <- read_matrix(args[1])

if (mode == "uncentred") {
  G <- crossprod(B)
  R <- tryCatch(chol(G), error = function(e) NULL)
  if (is.null(R)) {
    emit_json(list(mode = mode, positive_definite = FALSE), args[2])
  } else {
    r <- diag(R)
    loo <- 1 / (diag(G) * diag(chol2inv(R)))
    in_order <- r^2 / diag(G)
    ratio <- if (length(r) > 1) max(r[-1]) / min(r) else 1
    emit_json(list(
      mode              = mode,
      positive_definite = TRUE,
      loo               = loo,
      min_loo           = min(loo),
      in_order          = in_order,
      min_in_order      = min(in_order),
      getqf_ratio       = ratio,
      getqf_refuses     = ratio > 1e3
    ), args[2])
  }
} else if (mode == "centred") {
  X <- strip_intercept(B, "design")
  N <- nrow(X)
  Xc <- centre(X)
  G_explicit <- cross(Xc)
  R <- tryCatch(chol(G_explicit), error = function(e) NULL)
  if (is.null(R)) {
    emit_json(list(mode = mode, positive_definite = FALSE), args[2])
  } else {
    loo <- loo_qr(Xc)
    loo_explicit <- loo_gram(G_explicit)
    loo_recursive <- loo_gram(cross(Xc, summation = recursive_sum))
    loo_downdate <- loo_gram(cross(X) - N * tcrossprod(colMeans(X)))
    in_order <- diag(R)^2 / diag(G_explicit)
    emit_json(list(
      mode               = mode,
      positive_definite  = TRUE,
      rho                = colSums(Xc^2) / colSums(X^2),
      loo                = loo,
      min_loo            = min(loo),
      in_order           = in_order,
      min_in_order       = min(in_order),
      loo_uncentred      = loo_qr(cbind(1, X))[-1],
      loo_explicit       = loo_explicit,
      explicit_rel_error = rel_error(loo_explicit, loo),
      loo_explicit_recursive       = loo_recursive,
      explicit_recursive_rel_error = rel_error(loo_recursive, loo),
      loo_downdate       = loo_downdate,
      downdate_rel_error = rel_error(loo_downdate, loo)
    ), args[2], na = "null")
  }
} else {
  if (length(args) < 4) stop("hadamard mode needs H.csv")
  A <- strip_intercept(B, "A")
  H <- strip_intercept(read_matrix(args[4]), "H")
  if (ncol(A) != ncol(H)) stop("A and H need the same number of columns")
  n <- nrow(A); m <- nrow(H)
  X <- vapply(seq_len(ncol(A)), function(j) kronecker(A[, j], H[, j]), numeric(n * m))
  X <- matrix(X, n * m)
  Xc <- centre(X)
  G_explicit <- cross(Xc)
  R <- tryCatch(chol(G_explicit), error = function(e) NULL)
  if (is.null(R)) {
    emit_json(list(mode = mode, n = n, m = m, positive_definite = FALSE), args[2])
  } else {
    Ac <- centre(A); Hc <- centre(H)
    abar <- colMeans(A); hbar <- colMeans(H)
    HcHc <- cross(Hc); AcAc <- cross(Ac)
    G_three <- AcAc * HcHc + n * tcrossprod(abar) * HcHc + m * tcrossprod(hbar) * AcAc
    G_down <- cross(A) * cross(H) - n * m * tcrossprod(abar * hbar)
    loo <- loo_qr(Xc)
    loo_three <- loo_gram(G_three)
    loo_down <- loo_gram(G_down)
    scale <- max(abs(G_explicit))
    emit_json(list(
      mode                  = mode,
      n                     = n,
      m                     = m,
      positive_definite     = TRUE,
      rho                   = colSums(Xc^2) / colSums(X^2),
      loo                   = loo,
      min_loo               = min(loo),
      loo_three_term        = loo_three,
      three_term_rel_error  = rel_error(loo_three, loo),
      loo_downdate          = loo_down,
      downdate_rel_error    = rel_error(loo_down, loo),
      gram_three_term_error = max(abs(G_three - G_explicit)) / scale,
      gram_downdate_error   = max(abs(G_down - G_explicit)) / scale
    ), args[2], na = "null")
  }
}
