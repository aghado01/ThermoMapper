# khaos_sobol_oracle.R — sobol_khaos on given draws (E19's build brief, slice 4): the engine's draws' multi-indices, coefficients and
# s2, packed as khaos 2.1.1's adaptive_khaos object holds them (beta, vars, degs, nint, X, s2), so sobol_khaos reads exactly the
# engine's coefficients and the comparison is draw for draw. sobol_khaos divides each index by the explained variance plus s2.
#
# Usage (cwd = r/):  Rscript oracles/khaos_sobol_oracle.R <draws.csv> <out.json> <p>
#   draws: headerless numeric CSV rows kind,draw,term,a,b, draws and terms counted from 1, term 0 the intercept:
#     1,draw,term,coefficient,0     one per coefficient, the intercept's included;
#     2,draw,term,variable,degree   one per factor, variables counted from 1;
#     3,draw,0,s2,0                 one per draw.
# Emits { p, draws, labels, S (draws × labels, by row), T (draws × p, by row), leftover } from sobol_khaos(obj, plot_it = FALSE).

source("oracles/khaos_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: khaos_sobol_oracle.R <draws.csv> <out.json> <p>")
D <- read_matrix(args[1]); p <- as.integer(args[3])

coef <- D[D[, 1] == 1, , drop = FALSE]; fac <- D[D[, 1] == 2, , drop = FALSE]; s2r <- D[D[, 1] == 3, , drop = FALSE]
n_iter <- as.integer(max(D[, 2]))
M <- max(1L, as.integer(max(coef[, 3])))
q <- max(1L, if (nrow(fac)) max(table(paste(fac[, 2], fac[, 3]))) else 1L)

beta <- matrix(NA_real_, n_iter, M + 1)
nint <- matrix(NA_integer_, n_iter, M)
vars <- array(NA_integer_, c(n_iter, M, q))
degs <- array(NA_integer_, c(n_iter, M, q))
for (r in seq_len(nrow(coef))) beta[coef[r, 2], coef[r, 3] + 1] <- coef[r, 4]
for (r in seq_len(nrow(fac))) {
  i <- fac[r, 2]; m <- fac[r, 3]
  k <- if (is.na(nint[i, m])) 1L else nint[i, m] + 1L
  nint[i, m] <- k; vars[i, m, k] <- fac[r, 4]; degs[i, m, k] <- fac[r, 5]
}
s2 <- rep(NA_real_, n_iter); s2[s2r[, 2]] <- s2r[, 4]

obj <- structure(list(beta = beta, vars = vars, degs = degs, nint = nint, X = matrix(0, 1, p), s2 = s2), class = "adaptive_khaos")
sob <- sobol_khaos(obj, plot_it = FALSE)

emit_json(list(p = p, draws = n_iter, labels = names(sob$S), S = unname(as.matrix(sob$S)), T = unname(as.matrix(sob$T)),
               leftover = sob$leftover),
          args[2], null = "null", na = "null")
