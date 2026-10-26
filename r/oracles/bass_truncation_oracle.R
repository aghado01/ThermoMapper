# bass_truncation_oracle.R — reference for BASS's truncation residuals (port brief, outputs segment, D07; O-O4).
#
# Usage (cwd = r/):  Rscript oracles/bass_truncation_oracle.R <Y.csv> <out.json> [perc_var] [n_pc]
#   Y: headerless CSV, one case per row and one output location per column (n x m).
#   perc_var: the rank target in percent (default 99). r is the smallest rank whose cumulative variance share is at
#     least the target, capped at the numerical rank of the centred data (D07; R's strict '>' fails at 100).
#   n_pc: an explicit rank, overriding perc_var as in R's bassPCAsetup, capped the same way.
#
# Centres by the column mean with y_s = 1 (R's scale = FALSE); scaling is outside this oracle. The fixture's
# singular values must be distinct at the cut (d_r > d_{r+1}), or the rank-r subspace is not defined.
#
# Emits { n, m, centred_rank, r, full_retention, zero_tolerance, in_sample, in_sample_max_abs, in_sample_rms,
#         loo, loo_rank, loo_rms, ratio, coords_vs_output_in_sample, coords_vs_output_loo }:
#   in_sample[i, ] = y_ci - U_r U_r' y_ci, observed minus reconstructed (D07's orientation; R stores the negative).
#     At full retention (r = centred_rank) it is zero to within zero_tolerance = max(n, m) * eps * d_1, the
#     tolerance that judges every numerical rank here, d_1 being the largest singular value of all the data.
#   loo[i, ] = the same for case i against the basis refitted without it: centred by the other cases' mean, of rank
#     loo_rank[i] = min(r, the other cases' centred rank). Computed by n explicit refits in output space: the
#     reference. With n = 2 the other case alone has rank 0, so loo[i, ] is the difference of the two cases.
#   coords_vs_output_*: the largest absolute difference between those residuals and the same ones computed in case
#     coordinates from the one thin SVD, Y_c' = U D V' with C = D V', at O(n^4) cost independent of m. Every
#     leave-one-out quantity lies in span(U): case i centred by the others' mean is n/(n-1) c_i, and case j is
#     c_j + c_i/(n-1).
#   ratio = loo_rms / in_sample_rms, D07's basis-stability diagnostic; null when the in-sample residual is zero.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_truncation_oracle.R <Y.csv> <out.json> [perc_var] [n_pc]")

Y <- read_matrix(args[1])
perc_var <- if (length(args) >= 3) as.numeric(args[3]) else 99
n_pc <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
n <- nrow(Y); m <- ncol(Y)
if (n < 2) stop("need at least two cases")

choose_rank <- function(d, rho) {
  if (!is.na(n_pc)) return(min(n_pc, rho))
  if (rho == 0) return(0L)
  share <- cumsum(d^2) / sum(d^2)
  r <- which(share[seq_len(rho)] >= perc_var / 100)[1]
  if (is.na(r)) rho else r
}
# Residual of the columns of V against the orthonormal columns of U (an empty U leaves V unchanged).
residual <- function(V, U) if (ncol(U) == 0) V else V - U %*% crossprod(U, V)

# The fit on all cases, in output space.
Yc <- sweep(Y, 2, colMeans(Y))
s <- svd(t(Yc))                                            # m x n, as R's bassPCAsetup
# One tolerance, at the scale of all the data, judges every rank below, so that a leave-one-out matrix that is zero
# in exact arithmetic (n = 2) has rank 0 in both computations.
zero_tolerance <- max(n, m) * .Machine$double.eps * s$d[1]
numerical_rank <- function(d) sum(d > zero_tolerance)
rho <- numerical_rank(s$d)
r <- choose_rank(s$d, rho)
in_sample <- t(residual(t(Yc), s$u[, seq_len(r), drop = FALSE]))

# Leave-one-case-out by explicit refits in output space.
loo <- matrix(0, n, m)
loo_rank <- integer(n)
for (i in seq_len(n)) {
  Yo <- Y[-i, , drop = FALSE]
  mo <- colMeans(Yo)
  so <- svd(t(sweep(Yo, 2, mo)))
  loo_rank[i] <- min(r, numerical_rank(so$d))
  loo[i, ] <- residual(Y[i, ] - mo, so$u[, seq_len(loo_rank[i]), drop = FALSE])
}

# The same residuals in case coordinates, from the one SVD.
in_sample_c <- matrix(0, n, m)
loo_c <- matrix(0, n, m)
if (rho > 0) {
  Uq <- s$u[, seq_len(rho), drop = FALSE]
  C <- diag(s$d[seq_len(rho)], rho) %*% t(s$v[, seq_len(rho), drop = FALSE])        # Yc' = Uq C
  tail_k <- setdiff(seq_len(rho), seq_len(r))
  if (length(tail_k) > 0) in_sample_c <- t(Uq[, tail_k, drop = FALSE] %*% C[tail_k, , drop = FALSE])
  for (i in seq_len(n)) {
    Z <- C[, -i, drop = FALSE] + C[, i] / (n - 1)                                   # the others, centred
    sz <- svd(Z)
    ri <- min(r, numerical_rank(sz$d))
    loo_c[i, ] <- Uq %*% residual(C[, i] * n / (n - 1), sz$u[, seq_len(ri), drop = FALSE])
  }
}

in_sample_max_abs <- max(abs(in_sample))
in_sample_rms <- sqrt(mean(in_sample^2))
loo_rms <- sqrt(mean(loo^2))
emit_json(list(
  n                          = n,
  m                          = m,
  centred_rank               = rho,
  r                          = r,
  full_retention             = r == rho,
  zero_tolerance             = zero_tolerance,
  in_sample                  = in_sample,
  in_sample_max_abs          = in_sample_max_abs,
  in_sample_rms              = in_sample_rms,
  loo                        = loo,
  loo_rank                   = loo_rank,
  loo_rms                    = loo_rms,
  ratio                      = if (in_sample_max_abs > zero_tolerance) loo_rms / in_sample_rms else NA_real_,
  coords_vs_output_in_sample = max(abs(in_sample_c - in_sample)),
  coords_vs_output_loo       = max(abs(loo_c - loo))
), args[2], na = "null")
