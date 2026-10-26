# bass_admission_oracle.R — reference for BASS's numerical admission gate (port brief, D06).
#
# Usage (cwd = r/):  Rscript oracles/bass_admission_oracle.R <design.csv> <out.json>
#   design: headerless CSV; the columns of B in insertion order, the intercept first as in R BASS.
#
# Emits { positive_definite, loo, min_loo, in_order, min_in_order, getqf_ratio, getqf_refuses }:
#   loo[j] = 1 / (G_jj (G^-1)_jj), column j's squared sine to the span of all the other columns: the reciprocal
#     of its uncentred variance inflation factor, independent of column order (D06's gate).
#   in_order[j] = d_j / G_jj, the Cholesky pivots in the given order over the diagonal (FactorStatus's
#     MinRelativePivot); each is at least loo[j], and the last equals it.
#   getqf_ratio = max_{j >= 2} r_jj / min_j r_jj for R = chol(G), and getqf_refuses = getqf_ratio > 1e3: R BASS's
#     admission rule in getQf (mcmc_funcs.R), which depends on column norms; D06's BassLiteral option.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_admission_oracle.R <design.csv> <out.json>")

B <- read_matrix(args[1])
G <- crossprod(B)
R <- tryCatch(chol(G), error = function(e) NULL)

if (is.null(R)) {
  emit_json(list(positive_definite = FALSE), args[2])
} else {
  r <- diag(R)
  loo <- 1 / (diag(G) * diag(chol2inv(R)))
  in_order <- r^2 / diag(G)
  ratio <- if (length(r) > 1) max(r[-1]) / min(r) else 1
  emit_json(list(
    positive_definite = TRUE,
    loo               = loo,
    min_loo           = min(loo),
    in_order          = in_order,
    min_in_order      = min(in_order),
    getqf_ratio       = ratio,
    getqf_refuses     = ratio > 1e3
  ), args[2])
}
