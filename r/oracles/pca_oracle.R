# pca_oracle.R — reference PCA for cross-checking Maths.LinAlg.Pca.
#
# Usage (cwd = r/):  Rscript oracles/pca_oracle.R <fixture.csv> <out.json> [k]
# Emits { mean, eigenvalues, components } where components is k x d (rows = PCs).
#
# NOTE: prcomp eigenvalues use the (n-1) denominator; the C# Pca uses the MLE (n).
# Compare eigenVECTORS sign-agnostically (subspace distance) and rescale eigenvalues.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: pca_oracle.R <fixture.csv> <out.json> [k]")

X <- read_matrix(args[1])
k <- if (length(args) >= 3) as.integer(args[3]) else ncol(X)
idx <- seq_len(k)

pc <- prcomp(X, center = TRUE, scale. = FALSE)

emit_json(list(
  mean        = pc$center,
  eigenvalues = (pc$sdev^2)[idx],
  components  = t(pc$rotation[, idx, drop = FALSE])
), args[2])
