# tda_oracle.R — reference persistence diagram via Ripser (TDAstats::calculate_homology), the exact
# engine TDAkit::diagRips wraps. Full Vietoris-Rips on a point cloud, for cross-checking the C#
# RipsFiltration (complete graph) -> PersistentHomology pipeline against the gold-standard Ripser.
#
# Usage (cwd = r/): Rscript oracles/tda_oracle.R <cloud.csv> <out.json> [maxdim] [threshold]
# Emits { dimension[], birth[], death[] }; death = -1 marks an essential (infinite) bar.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: tda_oracle.R <cloud.csv> <out.json> [maxdim] [threshold]")

X <- read_matrix(args[1])
maxdim <- if (length(args) >= 3) as.integer(args[3]) else 1L
thr    <- if (length(args) >= 4) as.double(args[4]) else 1e9   # large => full Rips, no cutoff

h <- TDAstats::calculate_homology(X, dim = maxdim, threshold = thr, format = "cloud", return_df = TRUE)
h <- as.data.frame(h)
colnames(h) <- c("dimension", "birth", "death")

death <- as.double(h$death)
death[!is.finite(death)] <- -1   # JSON-safe sentinel for essential / infinite bars

emit_json(list(
  dimension = as.integer(h$dimension),
  birth     = as.double(h$birth),
  death     = death
), args[2])
