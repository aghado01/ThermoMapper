# mom_oracle.R — Grassmann geometric median via Riemann::riem.median.
#
# Usage (cwd = r/): Rscript oracles/mom_oracle.R <frames.csv> <out.json> <ambient_dim> <rank> [maxiter] [eps]
# Each CSV row is one Grassmann frame flattened column-major as a d x k basis.
# Emits { median, variation } where median is a d x k matrix.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: mom_oracle.R <frames.csv> <out.json> <ambient_dim> <rank> [maxiter] [eps]")
}

X <- read_matrix(args[1])
d <- as.integer(args[3])
k <- as.integer(args[4])
maxiter <- if (length(args) >= 5) as.integer(args[5]) else 200L
eps <- if (length(args) >= 6) as.double(args[6]) else 1e-5

if (ncol(X) != d * k) {
  stop(sprintf("frame width %d does not match ambient_dim * rank = %d", ncol(X), d * k))
}

frames <- vector("list", nrow(X))
for (i in seq_len(nrow(X))) {
  frames[[i]] <- matrix(as.numeric(X[i, ]), nrow = d, ncol = k)
}

riem <- Riemann::wrap.grassmann(frames)
med <- Riemann::riem.median(riem, geometry = "intrinsic", maxiter = maxiter, eps = eps)

emit_json(list(
  median = med$median,
  variation = med$variation
), args[2])
