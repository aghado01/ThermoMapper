# T4transport / lpSolve parity oracle for DiagramMetrics (sliced + Sinkhorn + exact Hungarian).
#
# Fixture: headerless CSV, one row per finite bar: side (0 = A, 1 = B), birth, death.
# Args:    <fixture.csv> <out.json> <p> <num_proj> <seed> <epsilon>
#
# The diagonal-augmented balanced geometry is reconstructed here, independently of the C#
# implementation, from the raw bars (side X = own bars + diagonal projections of the other side):
#   - swdist   : T4transport sliced Wasserstein on the augmented 2-D clouds. Slicing only ever
#                sees 1-D projected distances, so the L-inf vs Euclidean ground-metric question
#                does not arise — this is a full-semantics check of the screening metric.
#                Two T4transport quirks are reconciled here rather than C#-side: its scalar
#                `distance` is mean(W_p) over projections (not the documented (mean W_p^p)^(1/p)),
#                so we recombine from `projdist` in the documented convention; and its per-slice
#                1-D transport uses linearly interpolated ecdf quantiles on a 1000-point grid, a
#                smoothing that biases values at small cloud sizes — the fixture is sized (s ~ 100)
#                so the residual is a few percent, absorbed by the C# tolerance.
#   - lp_cost  : lpSolve::lp.assign optimum on the L-inf^p augmented cost matrix (unit masses) —
#                an external exact oracle for the Hungarian path.
#   - sinkhorn : T4transport::sinkhornD on the same L-inf augmented distance matrix, with
#                lambda = epsilon * max-finite-cost (the C# side normalizes costs by that max, so
#                its dimensionless epsilon maps to lambda exactly; plans then agree up to the
#                probability-vs-unit mass scale s).
#
# T4transport masses are probability (1/s per point); the C# metrics use unit mass per point, so
# swdist and sinkhorn distances relate to the C# values by the factor s^(1/p).

source("oracles/_common.R")
suppressMessages({
  library(T4transport)
  library(lpSolve)
})

args    <- commandArgs(trailingOnly = TRUE)
fixture <- args[[1]]
outp    <- args[[2]]
p       <- as.numeric(args[[3]])
nproj   <- as.integer(args[[4]])
seed    <- as.integer(args[[5]])
epsilon <- as.numeric(args[[6]])

bars <- read_matrix(fixture)
A <- bars[bars[, 1] == 0, 2:3, drop = FALSE]
B <- bars[bars[, 1] == 1, 2:3, drop = FALSE]
n <- nrow(A)
m <- nrow(B)
s <- n + m

diagproj <- function(P) {
  mid <- (P[, 1] + P[, 2]) / 2
  cbind(mid, mid)
}

Aaug <- rbind(A, diagproj(B))
Baug <- rbind(B, diagproj(A))

set.seed(seed)
sw <- swdist(Aaug, Baug, p = p, num_proj = nproj)
sw_rd <- mean(sw$projdist^p)^(1 / p)   # documented SW_p combination, from the per-slice values

# L-inf distances with the balanced structure: real-vs-real Dinf; each real bar may escape only
# to its own diagonal phantom (others forbidden); phantom-vs-phantom is free.
dinf  <- function(u, v) max(abs(u - v))
diagA <- (A[, 2] - A[, 1]) / 2
diagB <- (B[, 2] - B[, 1]) / 2

D <- matrix(0, s, s)
for (i in seq_len(n))
  for (j in seq_len(m))
    D[i, j] <- dinf(A[i, ], B[j, ])
# Prohibitive but tempered: big^p = (maxd+1)^p * (s+1) exceeds any legitimate assignment cost
# (<= s * maxd^p) while keeping the LP's dynamic range small enough for lpSolve's pivot
# tolerances — mirroring the C# BuildCost sizing. A 1e6-style sentinel loses ~2e-4 relative
# accuracy at p = 2 from the ~1e16 cost range.
maxd <- max(D, diagA, diagB)
big  <- (maxd + 1) * (s + 1)^(1 / p)
for (i in seq_len(n))
  for (k in seq_len(n))
    D[i, m + k] <- if (k == i) diagA[i] else big
for (j in seq_len(m))
  for (jj in seq_len(m))
    D[n + j, jj] <- if (jj == j) diagB[j] else big

lp <- lp.assign(D^p)

cmax   <- max(D[D < big])^p
lambda <- epsilon * cmax
sink   <- sinkhornD(D, p = p, lambda = lambda, maxiter = 20000, abstol = 1e-12)

emit_json(list(
  swdist   = sw_rd,
  lp_cost  = lp$objval,
  sinkhorn = sink$distance,
  s        = s
), outp)
