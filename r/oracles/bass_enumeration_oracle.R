# bass_enumeration_oracle.R — structure-space enumeration for BASS's birth, death and change moves (port brief, fitting
# segment): the reference for R's tempering and swap density (D05; F-O6 (a)) and for the count law each birth ratio
# implies (D01; F-O6 (c)).
#
# Usage (cwd = r/):  Rscript oracles/bass_enumeration_oracle.R <data.csv> <out.json> <mode> [key=value ...]
#   data: headerless CSV, y in the first column, then the inputs x_1..x_p in [0, 1] (R's scaled inputs).
#   mode: swap | count | laws | posterior.
#   keys (defaults): M_max (2), K_max (2), npart (3), w1 (5), w2 (5), lambda (2), sigma2 (0.1), tau (1/N),
#     rungs (1,0.5,0.2), g1 (0), g2 (0), h1 (10), h2 (10), s_min (1e-10), c (0.4).
#   posterior mode's own keys, refusing any it does not know: coef (sampled | fixed), centre (mean | flat | zero),
#     sigma2_prior (scaled | jeffreys), sigma2_shape (1e-3), tau_shape (0.5), likelihood (collapsed | constant),
#     exact_fit_c (1e3), rungs (0,0.05,0.2,0.5,1), categorical (0), L_max (3), subsets (size | uniform),
#     form (g | ridge); tau is the held tau (1/N, unit information, under the g-form; declared under the ridge).
#
# The structure space. A term is an order K <= K_max, K distinct variables, a sign and a knot per variable, the knots
# ranging over each variable's distinct values (D01's knot measure; R's unique.ind). Its prior is R's unconstrained
# measure, pi(t) = (1/K_max) C(p, K)^-1 (1/2)^K prod 1/n_v. A state is a set of M <= M_max distinct terms; under the
# count prior Pois(M; lambda) times the M! orderings, its prior weight is lambda^M prod pi(t), restricted to F (D02's
# global truncation). Birth proposals are R's: the order from I(K) = (w1 + order-K count) normalized, one variable
# uniformly for K = 1 and Wallenius with z_u = w2 + appearances for K >= 2 (brute force over draw orders), signs and
# knots uniformly. A change redraws one factor's sign and knot. Move-type probabilities are state dependent: birth
# only at M = 0, death and change 1/2 each at M_max, 1/3 each in between, as R selects them (mcmc_update.R).
# Every kernel's stationary law is found from its transition matrix on the enumerated space, independently of any
# target formula.
#
# mode = swap (F-O6 (a)). F = terms with at least npart nonzero values (R's support rule), designs [1, B_S] that R's
# getQf admits; getqf_mismatches counts states where getQf's verdict differs from the strict rank in either column
# order (0 makes F order-free). The continuous parameters are fixed at (sigma2, tau, lambda), the slice on which R's
# structure moves act. Kernels:
#   r_isolated: R's acceptance, the `alpha` expressions of birth_des, death_des and change_des (rjmcmc_des.R)
#     evaluated as written with R's getQf and logProbChangeMod, at inverse temperature b, but with the state's actual
#     move probabilities in data$birth.prob etc., so that R's fixed thirds (F-O4's defect) are kept apart;
#   r_literal: the same with R's fixed thirds, reported at b = 1 only;
#   corrected: D05's likelihood path, log p_b(new) - log p_b(old) plus the untempered prior and proposal ratios,
#     log p_b = -(Nb/2) log(2 pi sigma2) + (k/2) log(tau/(b + tau)) - b S_b/(2 sigma2), S_b = ||r||^2 + Q tau/(b + tau),
#     with the residual r computed directly.
# Candidate tempered targets at rung b (log scale, up to a constant): like_int, the prior times p_b (D05's default);
# like_lin, the prior times p_1^b; pow_coll, (prior times p_1)^b; pow_joint, the joint over (structure, a) powered
# in a's coordinates and then integrated. R's swap accepts with (b1 - b2)(lp(x2) - lp(x1)), lp evaluated by R's own
# `lp` (mcmc_funcs.R). It satisfies detailed balance for a pair of rung laws (s1, s2) exactly when
# D(S) = log s1(S) - log s2(S) - (b1 - b2) lp(S) is constant over S; swap_spread reports max D - min D for R's own
# stationary laws and for each candidate family (0 to rounding is balance). D05's swap uses log p_b in place of
# lp, and its spread is reported against the corrected kernels' stationary laws. reversibility_r and
# reversibility_corrected measure each kernel's detailed balance with its own stationary law; with p >= 3 the
# Wallenius factor of R's proposal enters R's tempered acceptance, raised to the power b - 1 (rungs must include 1).
#
# mode = count (F-O6 (c)). A constant likelihood and F = distinct terms with M <= M_max (no support or rank rule).
# Kernels: d01, a birth ratio of lambda pi(t) / (M + 1) times d_{M+1} / (b_M q(t)) (R's and D01's); khaos, the same
# with a further 1/(M + 1), the KHAOS manuscript's double count. The count law of each kernel's stationary law
# against the stated prior's, by enumeration of lambda^M prod pi(t) over F, and against that law times 1/M!.
#
# mode = laws (F-O3, F-O4), D01 and D02 as decided (2026-10-27). A constant likelihood; F = terms with at least npart
# nonzero values, and sets whose centred columns pass D06's order-free angle gate, s_j = 1/(G_jj (G^-1)_jj) >= s_min
# on the centred Gram (the intercept orthogonal to them). The count prior has lambda integrated out: a state's weight
# is Gamma(h1 + M) / (h2 + 1)^M times prod pi(t), restricted to F (global truncation), so a birth from M carries
# (M + h1)/(h2 + 1) for the count. Kernels, each from its own transition matrix:
#   d01: D01's move probabilities, equal over the moves available, in selection and ratio;
#   cancel: the prior-cancelling schedule, b_M = c min(1, r_M), d_{M+1} = c min(1, 1/r_M), the rest change, with
#     r_M = (M + h1)/((M + 1)(h2 + 1)), in selection and ratio (the option D01 keeps for F-M1);
#   thirds: D01's selection with R's fixed thirds in the ratios (F-O4's defect).
# Emits the count law and the joint law of term orders (n_1, ..., n_Kmax) of each against the target's enumeration;
# the boundary ratios pi(1)/pi(0) and pi(M_max)/pi(M_max - 1), and pi(2)/pi(1) when M_max >= 3, each kernel's over
# the target's; the admitted prior mass of each order, sum of pi(t | K) over single-term states in F, against the
# kernels' order law at M = 1 (the tilt D01 reports); and lambda's marginal at retention, the mixture of
# Gamma(h1 + M, rate h2 + 1) over the count law, against the prior Gamma(h1, rate h2) tilted by F in closed form.
#
# mode = posterior (stage 4's S6 and S8). The exact posterior over the laws mode's F given y. A state's prior is the
# laws target, Gamma(h1 + M)/(h2 + 1)^M prod pi(t), normalized over F. Its likelihood is D03's collapsed marginal,
# a and sigma^2 integrated out and tau integrated against its prior by composite Gauss-Legendre on u = log tau
# (16 nodes per unit panel, from log(1/N) - 37/a to log(1/N) + log(80/a)). Every intercept convention reads the same
# column space, so the marginal is computed on its own design: centre = mean (default, GPriorCentre.Mean) reads
# y~ = y - mean(y) on [1, Z_c], k = M + 1; zero (R's, GPriorCentre.Zero) reads y on [1, Z], k = M + 1; flat
# (Liang, GPriorCentre.Flat) reads y~ on Z_c, k = M, its flat intercept integrated with density 1. With r the
# residual (computed directly, never as y'y - Q) and Q the fit energy of that response on that design, and at
# inverse temperature b (b = 1 is the posterior), D05's integrated form is
#   log p_b(y | S, tau) = C - (n_b/2) log(2 pi) - [flat] log(bN)/2 + lgamma(d) + (k/2) log(tau/(tau + b))
#                         - d log(g2 + b S_b/2),   S_b = |r|^2 + Q tau/(tau + b),   d = g1 + n_b/2,
# with n_b = Nb (Nb - 1 under flat) and sigma^2 ~ IG(g1, g2), C = g1 log g2 - lgamma(g1): scaled (default),
# g1 = sigma2_shape, g2 = g1 s^2 with s^2 the mean square about the level (mean(y) under mean and flat, 0 under
# zero); jeffreys, p(sigma^2) = 1/sigma^2 as written (g1 = g2 = 0, C = 0), whose F excludes exact fits,
# |r| <= exact_fit_c eps |y~| (D03; refusals counted). tau ~ Gamma(tau_shape, rate tau_shape N), Zellner-Siow at
# 1/2; coef = fixed holds tau at the tau key (1/N, unit information). Only the defaults' proper priors make the
# marginal absolute; under flat or jeffreys log_Z is null and log_marginal is relative. likelihood = constant puts
# log p_b = 0, the laws mode's target. The scale of a column (R's basis normalization) does not move the marginal:
# a g-prior reads the column space only.
# form = ridge (the normalized ridge, decided 2026-10-29): beta | sigma^2, tau ~ N(0, sigma^2 (tau D)^-1) on the centred
# columns, D = diag(Z_c'Z_c)/N, so sigma^2 I/tau on the columns scaled to unit root mean square; the intercept, a unit-RMS
# column already, takes sigma^2/tau about its level inside it (mean, zero) or stays flat. Each state's design ([1, Z_c],
# or Z_c under flat) is scaled to unit RMS and its Gram decomposed by eigen(); with lambda_j its eigenvalues and q_j the
# fit energy of the convention's response along eigenvector j, (k/2) log(tau/(tau + b)) becomes
# sum_j (1/2) log(tau/(tau + b lambda_j)) and Q tau/(tau + b) becomes sum_j q_j tau/(tau + b lambda_j). tau's prior is
# Gamma(tau_shape, rate tau_shape), its unit information 1; a held ridge's tau must be declared (D03). The explicit
# check reads the prior precision tau diag(1, d) on [1, Z_c] (tau diag(d) on Z_c beside the flat intercept).
# Each closed form is checked in the oracle at the posterior's most probable state of each size, at b = 1 and 1/2
# and three tau (the fixed one under fixed-g): against quadrature over log sigma^2 of the Gaussian integral over a
# on the explicit uncentred design [1, Z] with y as given (prior mean mean(y) e_0 under mean; the flat intercept
# and the centred columns' g-prior under flat), and the composite rule against integrate() over log tau.
# Rungs (D05, for S8): the integrated form p_b above, and the linear form, the collapsed marginal p_1(y | S, tau)
# powered by b with tau's prior untempered; each rung's law over F is prior(S) int p(tau) exp(l_b) dtau. A rung is
# undefined (defined: false) where p_b's sigma^2 integral diverges: b = 0 under jeffreys or flat, b <= (1 - 2 g1)/N
# under flat. mean_slope and var_slope are E_b and Var_b of dl_b/db over (S, tau): the closed-form derivative of
# log p_b on the integrated form, log p_1 on the linear one, so that d log Z_b / db = mean_slope.
# Categorical predictors (D01's categorical part): categorical = n reads the last n input columns as factor codes,
# their levels the distinct codes the cases take, ascending (a one-level factor is refused). A term is then a
# continuous part of order J in {0..K_max} times a categorical part of order L in {0..L_max}, the product
# prod_l 1(z_l in D_l) over L distinct factors with each D_l a proper nonempty subset of its levels; J = L = 0 is the
# intercept, outside F. Its prior: J and L each uniform over their range, D01's continuous measure given J, the
# factors uniform over the C(p_cat, L) sets, and each subset by the size-uniform measure (default: the size uniform
# on {1..n_l - 1}, then the subset uniform given the size; R's) or uniform over the 2^n_l - 2 proper subsets (the
# paper's; the two agree up to three levels). D01 leaves the mixed order prior open and R's has none (its zero-order
# part scores lbmcmp = 0, its prior then the proposal's weight): the uniform order on {0..cap} for each part is this
# oracle's reading of "each with its own cap", flagged for the decision. Support counts the nonzero values of the
# whole product (R's birth_des_cat); D06's gate as before, which refuses a factor's main effect on D beside one on
# its complement (they span 1 with the intercept). Wallenius' variable-set probability is the proposal's at L >= 2
# and does not enter the posterior.
# Emits, beside the shared fields, in the C# conventions (0-based throughout): terms as { cont: [[var, sign, knot]],
# cat: [[var, [levels]]], func: [] }, var the position among the continuous or the categorical predictors, knot an
# index into the variable's ascending distinct values (BassInputs' knots), levels indices into its levels, and states
# as sorted indices into terms; supported (per term); data { y, x }; response { centre, level, s2, sigma2_shape,
# sigma2_rate, n_eff, absolute }; categorical { count, levels, L_max, subsets } (null without one); tau_prior
# { shape, rate } or { fixed }; n_gate_refused, n_exact_fit_refused;
# per state: M, k, rss, Q, log_prior, log_marginal, posterior, tau_mean (E[tau | S, y]); count { M, prior,
# posterior }; order { types (each term type's [J, L], L = 0 without categorical predictors), levels (the
# per-type term counts; n_1..n_Kmax without them), prior, posterior }; inclusion (per term);
# tau { mean, log_mean, log_sd, levels, quantiles } (null under fixed-g); log_Z; quadrature { states, sigma2_gap,
# tau_gap }; tempered { rungs, integrated, linear }, each rung { b, defined, log_Z, law, count, order, tau_mean,
# mean_slope, var_slope }.
#
# Emits { mode, n, p, n_terms, n_states, ... } with the mode's fields; states are listed as sorted term indices into
# `terms` (K, vars, signs, knots).

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: bass_enumeration_oracle.R <data.csv> <out.json> <swap|count|laws|posterior> [key=value ...]")
D <- read_matrix(args[1])
mode <- args[3]
y <- D[, 1]; X <- D[, -1, drop = FALSE]
N <- length(y); p <- ncol(X)
opt <- list(M_max = 2, K_max = 2, npart = 3, w1 = 5, w2 = 5, lambda = 2, sigma2 = 0.1, tau = 1 / N, rungs = "1,0.5,0.2",
            g1 = 0, g2 = 0, h1 = 10, h2 = 10, s_min = 1e-10, c = 0.4)
if (mode == "posterior")
  opt[c("coef", "centre", "sigma2_prior", "sigma2_shape", "tau_shape", "likelihood", "exact_fit_c", "rungs",
        "categorical", "L_max", "subsets")] <-
    list("sampled", "mean", "scaled", 1e-3, 0.5, "collapsed", 1e3, "0,0.05,0.2,0.5,1", 0, 3, "size")
if (mode == "posterior") opt$form <- "g"
string_keys <- c("rungs", "coef", "centre", "sigma2_prior", "likelihood", "subsets", "form")
tau_given <- any(startsWith(args[-(1:3)], "tau="))
for (kv in args[-(1:3)]) {
  parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
  if (mode == "posterior" && !parts[1] %in% names(opt)) stop("unknown key ", parts[1])
  opt[[parts[1]]] <- if (parts[1] %in% string_keys) parts[2] else as.numeric(parts[2])
}
n_cat <- if (mode == "posterior") opt$categorical else 0   # the last n_cat input columns are factor codes
if (n_cat > 0) {
  if (n_cat >= p) stop("the categorical predictors need a continuous one beside them")
  Z_cat <- X[, p - n_cat + seq_len(n_cat), drop = FALSE]
  X <- X[, seq_len(p - n_cat), drop = FALSE]; p <- ncol(X)
}
M_max <- opt$M_max; K_max <- min(opt$K_max, p)
bass <- bass_source(c("basis_funcs.R", "mcmc_funcs.R", "rjmcmc_des.R"))

# ---- terms -------------------------------------------------------------------------------------------------------
knot_values <- lapply(seq_len(p), function(v) if (mode == "posterior") sort(unique(X[, v])) else unique(X[, v]))
n_v <- lengths(knot_values)
terms <- list()
for (K in seq_len(K_max)) for (vars in split(combn(p, K), col(combn(p, K)))) {
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), K)))
  knots <- as.matrix(expand.grid(lapply(vars, function(v) seq_len(n_v[v]))))
  for (i in seq_len(nrow(signs))) for (j in seq_len(nrow(knots))) {
    kn <- vapply(seq_len(K), function(f) knot_values[[vars[f]]][knots[j, f]], 0)
    terms[[length(terms) + 1]] <- list(K = K, vars = vars, signs = unname(signs[i, ]), knot_index = unname(knots[j, ]),
                                       knots = kn)
  }
}
n_terms <- length(terms)
term_key <- vapply(terms, function(t) paste(t$K, paste(t$vars, collapse = ","), paste(t$signs, collapse = ","),
                                            paste(t$knot_index, collapse = ","), sep = "|"), "")
term_index <- new.env()
for (i in seq_len(n_terms)) assign(term_key[i], i, envir = term_index)
basis <- sapply(terms, function(t) bass$makeBasis(t$signs, t$vars, t$knots, t(X), 1))
supported <- if (mode %in% c("swap", "laws", "posterior")) colSums(basis != 0) >= opt$npart else rep(TRUE, n_terms)
log_prior_term <- vapply(terms, function(t) -log(K_max) - lchoose(p, t$K) - t$K * log(2) - sum(log(n_v[t$vars])), 0)
if (n_cat > 0) {   # every continuous part (J = 0 included) times every categorical part, less the intercept
  cat_levels <- lapply(seq_len(n_cat), function(v) sort(unique(Z_cat[, v])))
  n_lev <- lengths(cat_levels)
  if (any(n_lev < 2)) stop("a factor with one level is refused (D01)")
  if (!opt$subsets %in% c("size", "uniform")) stop("subsets must be size or uniform")
  L_max <- min(opt$L_max, n_cat)
  columns_of <- function(m) split(m, col(m))
  proper_subsets <- function(n) unlist(lapply(seq_len(n - 1), function(s) columns_of(combn(n, s))), recursive = FALSE)
  log_subset <- function(n, size) if (opt$subsets == "uniform") -log(2^n - 2) else -log(n - 1) - lchoose(n, size)
  cat_parts <- list(list(vars = integer(0), subsets = list()))
  for (L in seq_len(L_max)) for (vars in columns_of(combn(n_cat, L))) {
    choices <- lapply(vars, function(v) proper_subsets(n_lev[v]))
    pick <- as.matrix(expand.grid(lapply(choices, seq_along)))
    for (r in seq_len(nrow(pick)))
      cat_parts[[length(cat_parts) + 1]] <- list(vars = vars, subsets = lapply(seq_along(vars), function(l) choices[[l]][[pick[r, l]]]))
  }
  cat_basis <- sapply(cat_parts, function(cp) {
    col <- rep(1, N)
    for (l in seq_along(cp$vars)) col <- col * (Z_cat[, cp$vars[l]] %in% cat_levels[[cp$vars[l]]][cp$subsets[[l]]])
    col
  })
  cat_lp <- vapply(cat_parts, function(cp) -log(L_max + 1) - lchoose(n_cat, length(cp$vars)) +
                     sum(vapply(seq_along(cp$vars), function(l) log_subset(n_lev[cp$vars[l]], length(cp$subsets[[l]])), 0)), 0)
  cont_parts <- c(list(list(K = 0, vars = integer(0), signs = numeric(0), knot_index = integer(0), knots = numeric(0))), terms)
  cont_lp <- c(-log(K_max + 1), log_prior_term + log(K_max) - log(K_max + 1))   # J uniform on {0..K_max}
  pairs <- expand.grid(c = seq_along(cont_parts), k = seq_along(cat_parts))[-1, ]
  terms <- lapply(seq_len(nrow(pairs)), function(r) c(cont_parts[[pairs$c[r]]],
                                                      list(cat_vars = cat_parts[[pairs$k[r]]]$vars, cat_subsets = cat_parts[[pairs$k[r]]]$subsets)))
  basis <- cbind(1, basis)[, pairs$c, drop = FALSE] * cat_basis[, pairs$k, drop = FALSE]
  log_prior_term <- cont_lp[pairs$c] + cat_lp[pairs$k]
  n_terms <- length(terms)
  supported <- colSums(basis != 0) >= opt$npart   # the whole product's nonzero values (R's birth_des_cat)
}

# ---- proposals ----------------------------------------------------------------------------------------------------
permutations <- function(v) {
  if (length(v) <= 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) for (rest in permutations(v[-i])) out[[length(out) + 1]] <- c(v[i], rest)
  out
}
wallenius <- function(z, S) sum(vapply(permutations(S), function(o) {
  left <- rep(TRUE, length(z)); pr <- 1
  for (u in o) { pr <- pr * z[u] / sum(z[left]); left[u] <- FALSE }   # denominators sum the weights left
  pr
}, 0))
counts_of <- function(S) {
  I_star <- rep(opt$w1, K_max); z_star <- rep(opt$w2, p)
  for (i in S) { t <- terms[[i]]; I_star[t$K] <- I_star[t$K] + 1; z_star[t$vars] <- z_star[t$vars] + 1 }
  list(I = I_star / sum(I_star), z = z_star / sum(z_star))
}
log_q_term <- function(i, cnt) {   # birth proposal of term i given the counts
  t <- terms[[i]]
  v <- if (t$K == 1) 1 / p else wallenius(cnt$z, t$vars)
  log(cnt$I[t$K]) + log(v) - t$K * log(2) - sum(log(n_v[t$vars]))
}
schedule <- "equal"   # laws mode also builds the prior-cancelling schedule
count_ratio <- function(M) (M + opt$h1) / ((M + 1) * (opt$h2 + 1))   # p(M+1)/p(M) with lambda integrated out
move_probs <- function(M) {
  if (schedule == "cancel") {
    b <- if (M < M_max) opt$c * min(1, count_ratio(M)) else 0
    d <- if (M > 0) opt$c * min(1, 1 / count_ratio(M - 1)) else 0
    return(c(birth = b, death = d, change = 1 - b - d))
  }
  if (M == 0) return(c(birth = 1, death = 0, change = 0))
  if (M == M_max) return(c(birth = 0, death = 0.5, change = 0.5))
  c(birth = 1, death = 1, change = 1) / 3
}

# ---- states -------------------------------------------------------------------------------------------------------
candidates <- which(supported)
state_sets <- list(integer(0))
for (M in seq_len(M_max)) if (length(candidates) >= M)
  state_sets <- c(state_sets, split(combn(candidates, M), col(combn(candidates, M))))
state_sets <- lapply(state_sets, function(S) sort(as.integer(S)))

design <- function(S) cbind(1, basis[, S, drop = FALSE])
ssy <- sum(y^2)
describe <- function(S) {
  B <- design(S)
  if (mode == "count") return(list(valid = TRUE))
  if (mode %in% c("laws", "posterior")) {   # D06's numerical gate on the centred columns
    if (length(S) == 0) return(list(valid = TRUE))
    Bc <- scale(basis[, S, drop = FALSE], scale = FALSE)
    if (qr(Bc, tol = 1e-12)$rank < length(S)) return(list(valid = FALSE))
    G <- crossprod(Bc)
    s_j <- 1 / (diag(G) * diag(solve(G)))
    return(list(valid = all(s_j >= opt$s_min)))
  }
  qf_list <- bass$getQf(crossprod(B), crossprod(B, y))
  Brev <- B[, c(1, rev(seq_len(ncol(B))[-1])), drop = FALSE]
  rev_ok <- !is.null(bass$getQf(crossprod(Brev), crossprod(Brev, y)))
  strict <- qr(B, tol = 1e-9)$rank == ncol(B)
  mismatch <- (!is.null(qf_list) != strict) || (rev_ok != strict)
  if (is.null(qf_list) || !strict) return(list(valid = FALSE, mismatch = mismatch))
  list(valid = TRUE, mismatch = mismatch, qf_list = qf_list, qf = c(qf_list$qf),
       rss = sum(qr.resid(qr(B, tol = 0), y)^2), ldet_half = sum(log(diag(qf_list$R))))
}
info <- lapply(state_sets, describe)
getqf_mismatches <- sum(vapply(info, function(s) isTRUE(s$mismatch), TRUE))
keep <- vapply(info, function(s) s$valid, TRUE)
state_sets <- state_sets[keep]; info <- info[keep]
n_states <- length(state_sets)
state_key <- vapply(state_sets, function(S) paste0("s", paste(S, collapse = "-")), "")
state_index <- new.env()
for (i in seq_len(n_states)) assign(state_key[i], i, envir = state_index)
lookup <- function(S) { k <- paste0("s", paste(sort(S), collapse = "-")); if (exists(k, envir = state_index, inherits = FALSE)) get(k, envir = state_index) else NA_integer_ }
M_of <- lengths(state_sets)
log_prior <- vapply(state_sets, function(S) length(S) * log(opt$lambda) + sum(log_prior_term[S]), 0)

# ---- transitions ----------------------------------------------------------------------------------------------------
# Each proposal is a row (from, to, probability, type, term, log q of the birth it reverses or makes); rejections
# stay on the diagonal. type: 1 birth, 2 death, 3 change, 4 a change that redraws the same term.
type_names <- c("birth", "death", "change", "self")
rows_of <- function(i) {
  S <- state_sets[[i]]; M <- length(S); mp <- move_probs(M); cnt <- counts_of(S)
  out <- matrix(numeric(0), 0, 6)
  add <- function(j, prob, type, t, lq) out <<- rbind(out, c(i, j, prob, type, t, lq))
  if (mp["birth"] > 0) for (t in seq_len(n_terms)) {
    if (!supported[t] || t %in% S) next
    j <- lookup(c(S, t)); if (is.na(j)) next
    lq <- log_q_term(t, cnt)
    add(j, mp[["birth"]] * exp(lq), 1, t, lq)
  }
  if (mp["death"] > 0) for (t in S) add(lookup(setdiff(S, t)), mp[["death"]] / M, 2, t, log_q_term(t, counts_of(setdiff(S, t))))
  if (mp["change"] > 0) for (t in S) {
    tt <- terms[[t]]
    for (f in seq_len(tt$K)) for (sg in c(-1, 1)) for (kn in seq_len(n_v[tt$vars[f]])) {
      new_signs <- tt$signs; new_signs[f] <- sg; new_knots <- tt$knot_index; new_knots[f] <- kn
      t2 <- get(paste(tt$K, paste(tt$vars, collapse = ","), paste(new_signs, collapse = ","),
                      paste(new_knots, collapse = ","), sep = "|"), envir = term_index)
      prob <- mp[["change"]] / (M * tt$K * 2 * n_v[tt$vars[f]])
      if (t2 == t) { add(i, prob, 4, t, 0); next }
      if (!supported[t2] || t2 %in% S) next
      j <- lookup(c(setdiff(S, t), t2)); if (is.na(j)) next
      add(j, prob, 3, t2, 0)
    }
  }
  out
}
if (mode != "posterior") {   # the posterior is enumerated directly, without a kernel
  tr <- do.call(rbind, lapply(seq_len(n_states), rows_of))
  from <- as.integer(tr[, 1]); to <- as.integer(tr[, 2]); prob <- tr[, 3]
  type <- type_names[tr[, 4]]; term <- as.integer(tr[, 5]); lq_fwd <- tr[, 6]
}

kernel <- function(log_accept) {
  a <- ifelse(type == "self", 1, pmin(1, exp(log_accept)))
  agg <- rowsum(prob * a, (to - 1) * n_states + from)
  P <- matrix(0, n_states, n_states)
  P[as.numeric(rownames(agg))] <- agg
  diag(P) <- diag(P) + 1 - rowSums(P)
  P
}
# Swap mode compares log ratios of stationary laws, so it uses GTH elimination (Grassmann, Taksar and Heyman 1985),
# which keeps every component's relative accuracy; count mode's laws are absolute quantities, and a dense solve
# serves its larger space.
stationary <- function(P) {
  if (mode %in% c("count", "laws")) {
    A <- t(P) - diag(n_states); A[n_states, ] <- 1
    s <- solve(A, c(rep(0, n_states - 1), 1))
    return(s / sum(s))
  }
  stationary_gth(P)
}
normalize_log <- function(l) { w <- exp(l - max(l)); w / sum(w) }
tv <- function(a, b) 0.5 * sum(abs(a - b))

# The prior-and-proposal part of each move's log ratio, untempered, with the state's actual move probabilities:
# birth rate(M) pi(t) d_{M+1} / ((M + 1) b_M q(t)); death its reciprocal from the smaller state; change 0. The count
# factor is lambda by default, (M + h1)/(h2 + 1) with lambda integrated out; ratio_probs supplies the move
# probabilities the ratio reads, which a defect may hold fixed while selection follows the state.
prior_proposal <- function(extra_birth = function(M) 0, log_rate = function(M) log(opt$lambda), ratio_probs = move_probs) {
  vapply(seq_along(from), function(r) {
    if (type[r] %in% c("self", "change")) return(0)
    Mf <- M_of[from[r]]
    if (type[r] == "birth") {
      mp_to <- ratio_probs(Mf + 1)
      return(log_rate(Mf) + log_prior_term[term[r]] + log(mp_to[["death"]] / (Mf + 1)) -
               log(ratio_probs(Mf)[["birth"]]) - lq_fwd[r] + extra_birth(Mf))
    }
    mp_to <- ratio_probs(Mf - 1)
    -(log_rate(Mf - 1) + log_prior_term[term[r]] + log(ratio_probs(Mf)[["death"]] / Mf) - log(mp_to[["birth"]]) -
        lq_fwd[r] + extra_birth(Mf - 1))
  }, 0)
}

result <- list(mode = mode, n = N, p = p, options = opt, n_terms = n_terms, n_supported = sum(supported),
               n_states = n_states, terms = lapply(terms, function(t) t[c("K", "vars", "signs", "knots")]),
               states = state_sets)

if (mode == "count") {
  pp <- prior_proposal()
  laws <- list()
  count_law <- function(s) vapply(0:M_max, function(M) sum(s[M_of == M]), 0)
  stated <- normalize_log(log_prior)
  s_d01 <- stationary(kernel(pp))
  s_khaos <- stationary(kernel(prior_proposal(function(M) -log(M + 1))))
  stated_over_factorial <- normalize_log(log_prior - lfactorial(M_of))
  result$count <- list(M = 0:M_max, stated = count_law(stated), d01 = count_law(s_d01), khaos = count_law(s_khaos),
                       stated_over_factorial = count_law(stated_over_factorial),
                       tv_d01_stated = tv(s_d01, stated), tv_khaos_stated = tv(s_khaos, stated),
                       tv_khaos_over_factorial = tv(s_khaos, stated_over_factorial))
}

if (mode == "laws") {
  h1 <- opt$h1; h2 <- opt$h2
  nb_rate <- function(M) log((M + h1) / (h2 + 1))
  thirds <- function(M) c(birth = 1, death = 1, change = 1) / 3
  log_target <- vapply(state_sets, function(S) lgamma(h1 + length(S)) - length(S) * log(h2 + 1) + sum(log_prior_term[S]), 0)
  target <- normalize_log(log_target)
  count_law <- function(s) vapply(0:M_max, function(M) sum(s[M_of == M]), 0)
  orders <- t(vapply(state_sets, function(S) tabulate(vapply(S, function(i) terms[[i]]$K, 0), K_max), numeric(K_max)))
  order_key <- apply(orders, 1, paste, collapse = ",")
  order_levels <- unique(order_key[order(rowSums(orders), order_key)])
  order_law <- function(s) vapply(order_levels, function(k) sum(s[order_key == k]), 0)

  s_d01 <- stationary(kernel(prior_proposal(log_rate = nb_rate)))
  s_thirds <- stationary(kernel(prior_proposal(log_rate = nb_rate, ratio_probs = thirds)))
  schedule <- "cancel"   # rebuild the transitions with the prior-cancelling selection probabilities
  tr <- do.call(rbind, lapply(seq_len(n_states), rows_of))
  from <- as.integer(tr[, 1]); to <- as.integer(tr[, 2]); prob <- tr[, 3]
  type <- type_names[tr[, 4]]; term <- as.integer(tr[, 5]); lq_fwd <- tr[, 6]
  s_cancel <- stationary(kernel(prior_proposal(log_rate = nb_rate)))
  schedule <- "equal"
  # The prior-cancelling selection probabilities are themselves proper: every one in [0, 1], summing to 1.
  cancel_probs <- t(vapply(0:M_max, function(M) { schedule <<- "cancel"; v <- move_probs(M); schedule <<- "equal"; v }, numeric(3)))

  ratio_at <- function(s, M) { cl <- count_law(s); cl[M + 1] / cl[M] }
  tc <- count_law(target)
  boundary <- function(s) c(low = ratio_at(s, 1) / ratio_at(target, 1), high = ratio_at(s, M_max) / ratio_at(target, M_max),
                            middle = if (M_max >= 3) ratio_at(s, 2) / ratio_at(target, 2) else NA)

  # The tilt: the admitted prior mass of each order among single-term states, against the order law at M = 1.
  single <- which(M_of == 1)
  single_K <- vapply(state_sets[single], function(S) terms[[S]]$K, 0)
  admitted <- vapply(seq_len(K_max), function(K) sum(exp(log_prior_term[unlist(state_sets[single][single_K == K])] + log(K_max))), 0)
  order_at_one <- function(s) vapply(seq_len(K_max), function(K) sum(s[single][single_K == K]), 0) / sum(s[single])

  # lambda at retention: the mixture over the count law against the prior tilted by F, in closed form.
  c_M <- vapply(0:M_max, function(M) sum(exp(vapply(state_sets[M_of == M], function(S) sum(log_prior_term[S]), 0))), 0)
  log_Z <- log(sum(c_M * exp(lgamma(h1 + 0:M_max) - (h1 + 0:M_max) * log(h2 + 1)))) + h1 * log(h2) - lgamma(h1)
  lambda_grid <- c(0.25, 0.5, 1, 2, 5, 10, 20)
  tilted <- vapply(lambda_grid, function(l) exp(h1 * log(h2) - lgamma(h1) + (h1 - 1) * log(l) - (h2 + 1) * l +
                                                  log(sum(c_M * l^(0:M_max))) - log_Z), 0)
  mixture <- vapply(lambda_grid, function(l) sum(count_law(s_d01) * dgamma(l, h1 + 0:M_max, rate = h2 + 1)), 0)

  result$laws <- list(
    M = 0:M_max, n_gate_refused = sum(!keep),
    count = list(target = tc, d01 = count_law(s_d01), cancel = count_law(s_cancel), thirds = count_law(s_thirds)),
    tv = list(d01 = tv(s_d01, target), cancel = tv(s_cancel, target), thirds = tv(s_thirds, target)),
    order = list(levels = order_levels, target = order_law(target), d01 = order_law(s_d01), cancel = order_law(s_cancel)),
    tv_order = list(d01 = tv(order_law(s_d01), order_law(target)), cancel = tv(order_law(s_cancel), order_law(target))),
    boundary = list(d01 = boundary(s_d01), cancel = boundary(s_cancel), thirds = boundary(s_thirds)),
    cancel_probs = cancel_probs,
    admitted_mass = admitted, order_at_one = list(target = order_at_one(target), d01 = order_at_one(s_d01)),
    lambda = list(grid = lambda_grid, mixture = mixture, tilted = tilted))
}

if (mode == "swap") {
  rungs <- as.numeric(strsplit(opt$rungs, ",")[[1]])
  s2 <- opt$sigma2; tau <- opt$tau; lam <- opt$lambda
  qf <- vapply(info, `[[`, 0, "qf"); rss <- vapply(info, `[[`, 0, "rss"); ldh <- vapply(info, `[[`, 0, "ldet_half")
  k <- M_of + 1
  log_pb <- function(b) -(N * b / 2) * log(2 * pi * s2) + (k / 2) * log(tau / (b + tau)) -
    b * (rss + qf * tau / (b + tau)) / (2 * s2)
  # R's lp, as written, at each state on the slice.
  prior_r <- list(g1 = opt$g1, g2 = opt$g2, h1 = opt$h1, h2 = opt$h2, beta.gprior.ind = 1, a.beta.prec = 0.5,
                  b.beta.prec = 2 / N)
  lp_r <- vapply(seq_len(n_states), function(i) {
    ql <- info[[i]]$qf_list
    curr <- list(s2.rate = c(0.5 * (ssy - crossprod(ql$bhat, crossprod(design(state_sets[[i]]), y)) / (1 + tau))),
                 s2 = s2, nbasis = M_of[i], R = ql$R, lam = lam, beta.prec = tau)
    c(bass$lp(curr, prior_r, list(n = N)))
  }, 0)
  # R's alpha, evaluated as written at b = 1; the move probabilities are the state's (isolated) or R's thirds.
  alpha_expr <- list(birth = source_assignment(bass$birth_des, "alpha"), death = source_assignment(bass$death_des, "alpha"),
                     change = source_assignment(bass$change_des, "alpha"))
  r_alpha <- function(r, b, literal) {
    i <- from[r]; j <- to[r]; S <- state_sets[[i]]; M <- M_of[i]
    probs <- if (literal) list(birth.prob = 1 / 3, death.prob.next = 1 / 3, death.prob = 1 / 3, birth.prob.last = 1 / 3)
             else list(birth.prob = move_probs(M)[["birth"]], death.prob.next = move_probs(M + 1)[["death"]],
                       death.prob = move_probs(M)[["death"]], birth.prob.last = move_probs(max(M - 1, 0))[["birth"]])
    lbmcmp <- 0
    if (type[r] %in% c("birth", "death")) {
      t <- terms[[term[r]]]
      cnt <- counts_of(if (type[r] == "birth") S else setdiff(S, term[r]))
      lbmcmp <- bass$logProbChangeMod(t$K, t$vars, cnt$I, cnt$z, p, n_v, K_max, 0)
    }
    e <- list2env(list(
      curr = list(s2 = s2, qf = info[[i]]$qf, beta.prec = tau, lam = lam, nc = M + 1, nbasis = M, R = info[[i]]$qf_list$R,
                  temp.ind = 1),
      data = c(list(itemp.ladder = b), probs), prior = list(beta.jprior.ind = 0),
      cand.des = list(lbmcmp = lbmcmp), lpbmcmp = lbmcmp, qf.cand.list = info[[j]]$qf_list))
    c(eval(alpha_expr[[type[r]]], e))
  }
  moving <- which(type != "self")
  alpha1_iso <- vapply(seq_along(from), function(r) if (type[r] == "self") 0 else r_alpha(r, 1, FALSE), 0)
  alpha1_lit <- vapply(seq_along(from), function(r) if (type[r] == "self") 0 else r_alpha(r, 1, TRUE), 0)
  # R's alpha is b times its b = 1 value: checked on every 25th moving transition at b = 0.5.
  probe <- moving[seq(1, length(moving), by = 25)]
  linear_gap <- max(abs(vapply(probe, function(r) r_alpha(r, 0.5, FALSE), 0) - 0.5 * alpha1_iso[probe]))

  pp <- prior_proposal()
  dl <- function(f) f[to] - f[from]
  posterior <- normalize_log(log_prior + log_pb(1))
  targets <- function(b) list(
    like_int = log_prior + log_pb(b),
    like_lin = log_prior + b * log_pb(1),
    pow_coll = b * (log_prior + log_pb(1)),
    pow_joint = b * (log_prior - N / 2 * log(2 * pi * s2) + k / 2 * log(tau / (2 * pi * s2)) + ldh) -
      b * (rss + qf * tau / (1 + tau)) / (2 * s2) + k / 2 * log(2 * pi * s2 / (b * (1 + tau))) - ldh)
  # Detailed balance of a kernel with its own stationary law: the largest |s_x P_xy - s_y P_yx| relative to the
  # largest flow s_x P_xy between distinct states. Zero to rounding means the kernel is reversible.
  reversibility_gap <- function(P, st) {
    Fl <- st * P; diag(Fl) <- 0
    max(abs(Fl - t(Fl))) / max(Fl)
  }
  rung <- lapply(rungs, function(b) {
    P_r <- kernel(b * alpha1_iso); P_c <- kernel(dl(log_pb(b)) + pp)
    s_r <- stationary(P_r); s_c <- stationary(P_c)
    tg <- lapply(targets(b), normalize_log)
    list(b = b, r_isolated = s_r, corrected = s_c, tv_corrected_target = tv(s_c, tg$like_int),
         tv_r_to = lapply(tg, function(g) tv(s_r, g)),
         reversibility_r = reversibility_gap(P_r, s_r), reversibility_corrected = reversibility_gap(P_c, s_c))
  })
  b1_rung <- rung[[match(1, rungs)]]
  spread <- function(d) max(d) - min(d)
  pairs <- lapply(seq_len(length(rungs) - 1), function(a) {
    b1 <- rungs[a]; b2 <- rungs[a + 1]
    t1 <- targets(b1); t2 <- targets(b2)
    list(b1 = b1, b2 = b2,
         r_swap_vs_r_kernels = spread(log(rung[[a]]$r_isolated) - log(rung[[a + 1]]$r_isolated) - (b1 - b2) * lp_r),
         r_swap_vs_family = lapply(setNames(names(t1), names(t1)), function(f) spread(t1[[f]] - t2[[f]] - (b1 - b2) * lp_r)),
         corrected_swap_vs_corrected_kernels = spread(log(rung[[a]]$corrected) - log(rung[[a + 1]]$corrected) -
                                                        (log_pb(b1) - log_pb(b2))))
  })
  result$swap <- list(
    getqf_mismatches = getqf_mismatches, rungs = rungs, n_transitions = length(moving), alpha_linear_in_b_gap = linear_gap,
    tv_r_isolated_b1_posterior = tv(b1_rung$r_isolated, posterior),
    tv_r_literal_b1_posterior = tv(stationary(kernel(alpha1_lit)), posterior),
    rung = lapply(rung, function(x) x[c("b", "tv_corrected_target", "tv_r_to", "reversibility_r", "reversibility_corrected")]),
    corrected_target = lapply(rungs, function(b) normalize_log(log_prior + log_pb(b))),
    pairs = pairs, lp_r = lp_r)
}

if (mode == "posterior") {
  centre <- opt$centre; flat <- centre == "flat"; jeffreys <- opt$sigma2_prior == "jeffreys"
  fixed_g <- opt$coef == "fixed"; constant <- opt$likelihood == "constant"
  if (!centre %in% c("mean", "flat", "zero")) stop("centre must be mean, flat or zero")
  if (!opt$sigma2_prior %in% c("scaled", "jeffreys")) stop("sigma2_prior must be scaled or jeffreys")
  if (!opt$coef %in% c("sampled", "fixed")) stop("coef must be sampled or fixed")
  if (!opt$likelihood %in% c("collapsed", "constant")) stop("likelihood must be collapsed or constant")
  if (!opt$form %in% c("g", "ridge")) stop("form must be g or ridge")
  ridge <- opt$form == "ridge"
  if (ridge && fixed_g && !tau_given) stop("a held ridge's tau is part of the model and must be declared (D03): pass tau=")
  n0 <- if (ridge) 1 else N   # the form's unit-information precision is 1/n0
  rungs <- as.numeric(strsplit(opt$rungs, ",")[[1]])
  level <- if (centre == "zero") 0 else mean(y)
  if (!constant && (if (centre == "zero") all(y == 0) else all(y == y[1])))
    stop("a zero response (zero centre) or a constant one (mean, flat) is refused (D03)")
  yt <- y - level
  s2 <- mean(yt^2)
  g1 <- if (jeffreys) 0 else opt$sigma2_shape
  g2 <- g1 * s2
  ig_const <- if (jeffreys) 0 else g1 * log(g2) - lgamma(g1)
  lse <- function(l) { m <- max(l); m + log(sum(exp(l - m))) }
  row_lse <- function(A) { m <- A[cbind(seq_len(nrow(A)), max.col(A, "first"))]; m + log(rowSums(exp(A - m))) }

  # Each state's fit of the response its convention reads, on that convention's design.
  design_of <- function(S) {
    Z <- basis[, S, drop = FALSE]
    if (centre == "zero") return(cbind(1, Z))
    Zc <- Z - rep(colMeans(Z), each = N)
    if (flat) Zc else cbind(1, Zc)
  }
  fit_stats <- t(vapply(state_sets, function(S) {
    Bs <- design_of(S)
    if (ncol(Bs) == 0) return(c(rss = sum(yt^2), Q = 0))
    fitted <- qr.fitted(qr(Bs, tol = 0), yt)
    c(rss = sum((yt - fitted)^2), Q = sum(fitted^2))
  }, c(rss = 0, Q = 0)))
  exact_fit <- if (jeffreys && !constant) sqrt(fit_stats[, "rss"]) <= opt$exact_fit_c * .Machine$double.eps * sqrt(sum(yt^2))
               else rep(FALSE, n_states)
  state_sets <- state_sets[!exact_fit]; fit_stats <- fit_stats[!exact_fit, , drop = FALSE]
  n_states <- length(state_sets); M_of <- lengths(state_sets)
  rss <- fit_stats[, "rss"]; Q <- fit_stats[, "Q"]; k <- M_of + !flat

  # The ridge's spectrum per state: the convention's design scaled to unit RMS, its Gram's eigenvalues, and the fit
  # energy of the convention's response along each eigenvector.
  ridge_stats <- if (ridge) lapply(state_sets, function(S) {
    Z <- basis[, S, drop = FALSE]
    Zc <- Z - rep(colMeans(Z), each = N)
    Xr <- if (flat) Zc else cbind(1, Zc)
    if (ncol(Xr) == 0) return(list(lambda = numeric(0), q = numeric(0)))
    Xn <- sweep(Xr, 2, sqrt(colSums(Xr^2) / N), "/")
    e <- eigen(crossprod(Xn), symmetric = TRUE)
    w <- c(crossprod(e$vectors, crossprod(Xn, yt)))
    list(lambda = e$values, q = w^2 / e$values)
  }) else NULL

  log_target <- vapply(state_sets, function(S) lgamma(opt$h1 + length(S)) - length(S) * log(opt$h2 + 1) +
                         sum(log_prior_term[S]), 0)
  log_prior_F <- log_target - lse(log_target)

  # tau's quadrature on u = log tau, with the prior's density on u folded into the weights; fixed-g is a point.
  a_tau <- opt$tau_shape; r_tau <- a_tau * n0
  log_density_u <- function(u) a_tau * (u + log(r_tau)) - r_tau * exp(u) - lgamma(a_tau)
  if (fixed_g) {
    u <- log(opt$tau); logw <- 0; panel <- 1
  } else {
    lo <- log(1 / n0) - 37 / a_tau; hi <- log(1 / n0) + log(80 / a_tau)
    if (hi - lo > 4000) stop("tau's prior is too diffuse for the quadrature")
    edges <- seq(lo, hi, length.out = ceiling(hi - lo) + 1)
    gl <- gauss_legendre(16)
    width <- diff(edges)
    u <- as.vector(outer((gl$x + 1) / 2, width) + rep(edges[-length(edges)], each = 16))
    logw <- log(as.vector(outer(gl$w / 2, width))) + log_density_u(u)
    panel <- rep(seq_along(width), each = 16)
  }
  tau <- exp(u)

  # log p_b(y | S, tau) and its b-derivative, at the states idx and the points tt (states by points).
  defined <- function(b) constant || (g1 + (N * b - flat) / 2 > 0 && (b > 0 || !(jeffreys || flat)))
  lp <- function(b, tt = tau, idx = seq_len(n_states)) {
    if (constant) return(matrix(0, length(idx), length(tt)))
    d <- g1 + (N * b - flat) / 2
    if (ridge) {
      base <- ig_const - (N * b - flat) / 2 * log(2 * pi) - (if (flat) 0.5 * log(b * N) else 0) + lgamma(d)
      out <- matrix(0, length(idx), length(tt))
      for (r in seq_along(idx)) {
        st <- ridge_stats[[idx[r]]]
        h <- outer(st$lambda * b, tt, function(bl, t0) t0 / (t0 + bl))   # columns by tau
        Tb <- b * (rss[idx[r]] + colSums(st$q * h))
        out[r, ] <- base + 0.5 * colSums(log(h)) - d * log(g2 + Tb / 2)
      }
      return(out)
    }
    shrink <- tt / (tt + b)
    Tb <- b * (rss[idx] + outer(Q[idx], shrink))
    ig_const - (N * b - flat) / 2 * log(2 * pi) - (if (flat) 0.5 * log(b * N) else 0) + lgamma(d) +
      outer(k[idx] / 2, log(shrink)) - d * log(g2 + Tb / 2)
  }
  dlp <- function(b, tt = tau) {
    if (constant) return(matrix(0, n_states, length(tt)))
    d <- g1 + (N * b - flat) / 2
    if (ridge) {
      out <- matrix(0, n_states, length(tt))
      for (i in seq_len(n_states)) {
        st <- ridge_stats[[i]]
        h <- outer(st$lambda * b, tt, function(bl, t0) t0 / (t0 + bl))
        Tb <- b * (rss[i] + colSums(st$q * h)); dTb <- rss[i] + colSums(st$q * h^2)
        dld <- -0.5 * colSums(st$lambda * h / rep(tt, each = length(st$lambda)))   # -(1/2) sum lambda/(tau + b lambda)
        out[i, ] <- -(N / 2) * log(2 * pi) - (if (flat) 1 / (2 * b) else 0) + (N / 2) * digamma(d) + dld -
          (N / 2) * log(g2 + Tb / 2) - d * (dTb / 2) / (g2 + Tb / 2)
      }
      return(out)
    }
    shrink <- tt / (tt + b)
    Tb <- b * (rss + outer(Q, shrink)); dTb <- rss + outer(Q, shrink^2)
    -(N / 2) * log(2 * pi) - (if (flat) 1 / (2 * b) else 0) + (N / 2) * digamma(d) - outer(k / 2, 1 / (tt + b)) -
      (N / 2) * log(g2 + Tb / 2) - d * (dTb / 2) / (g2 + Tb / 2)
  }

  count_law <- function(s) vapply(0:M_max, function(M) sum(s[M_of == M]), 0)
  JL <- t(vapply(terms, function(t) c(t$K, length(t$cat_vars)), c(0, 0)))   # each term's (J, L)
  types <- unique(JL); types <- types[order(rowSums(types), types[, 1]), , drop = FALSE]
  type_of <- match(paste(JL[, 1], JL[, 2]), paste(types[, 1], types[, 2]))
  orders <- matrix(vapply(state_sets, function(S) tabulate(type_of[S], nrow(types)), numeric(nrow(types))),
                   ncol = nrow(types), byrow = TRUE)
  order_key <- apply(orders, 1, paste, collapse = ",")
  order_levels <- unique(order_key[order(rowSums(orders), order_key)])
  order_law <- function(s) vapply(order_levels, function(kk) sum(s[order_key == kk]), 0, USE.NAMES = FALSE)
  absolute <- constant || !(jeffreys || flat)

  # One rung: the per-state tempered marginal, the law over F, and the slope's moments over (S, tau).
  rung_at <- function(b, form) {
    if (form == "integrated" && !defined(b)) return(list(b = b, defined = FALSE, log_Z = NA))
    l1 <- if (form == "linear") lp(1) else NULL
    A <- sweep(if (form == "linear") b * l1 else lp(b), 2, logw, "+")
    logm <- row_lse(A)
    W <- exp(A - logm)                       # tau | S, y at rung b, as weights over the nodes
    law <- exp(log_prior_F + logm - lse(log_prior_F + logm))
    slope <- if (form == "linear") l1 else dlp(b)
    s1 <- rowSums(W * slope); s2_ <- rowSums(W * slope^2)
    mean_slope <- sum(law * s1); var_slope <- sum(law * s2_) - mean_slope^2
    # At b = 0 the integrated slope carries -(k/2)/tau, whose prior moments E[1/tau] and E[1/tau^2] are infinite for
    # tau shapes up to 1 and 2: null there, not the truncated quadrature's number.
    singular <- form == "integrated" && b == 0 && !fixed_g && !constant
    if (singular && a_tau <= 1) mean_slope <- NA
    if (singular && a_tau <= 2) var_slope <- NA
    list(b = b, defined = TRUE, log_Z = if (absolute) lse(log_prior_F + logm) else NA, law = law,
         count = count_law(law), order = order_law(law), tau_mean = sum(law * rowSums(W * rep(tau, each = n_states))),
         mean_slope = mean_slope, var_slope = var_slope, log_marginal = logm, W = W)
  }

  post <- rung_at(1, "integrated")
  pst <- post$law; logm <- post$log_marginal; W <- post$W
  tau_given_S <- rowSums(W * rep(tau, each = n_states))
  incl <- numeric(n_terms)
  incl_sum <- rowsum(rep(pst, M_of), unlist(state_sets))
  incl[as.integer(rownames(incl_sum))] <- incl_sum[, 1]

  tau_summary <- NA   # null in the JSON
  if (!fixed_g) {
    mix <- colSums(pst * W)                  # the posterior of u on the nodes, mixed over states
    log_mean <- sum(mix * u)
    mix_density <- function(uu) colSums(pst * exp(sweep(lp(1, exp(uu)), 2, log_density_u(uu), "+") - logm))
    cum <- c(0, cumsum(tapply(mix, panel, sum)))
    cdf <- function(uq) {
      j <- findInterval(uq, edges, rightmost.closed = TRUE)
      nodes <- edges[j] + (gl$x + 1) / 2 * (uq - edges[j])
      cum[j] + sum(gl$w / 2 * (uq - edges[j]) * mix_density(nodes))
    }
    levels <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    quantiles <- vapply(levels, function(q) {
      j <- findInterval(q, cum, rightmost.closed = TRUE)
      exp(uniroot(function(uq) cdf(uq) - q, c(edges[j], edges[j + 1]), tol = 1e-13)$root)
    }, 0)
    tau_summary <- list(mean = sum(pst * tau_given_S), log_mean = log_mean, log_sd = sqrt(sum(mix * (u - log_mean)^2)),
                        levels = levels, quantiles = quantiles)
  }

  # The closed forms against quadrature over sigma^2 of the explicit Gaussian integral over a, on [1, Z] with y as
  # given, and the composite rule over tau against integrate(), at the posterior's most probable state of each size.
  explicit_lpb <- function(S, b, t0) {
    Z <- basis[, S, drop = FALSE]; Xd <- cbind(1, Z); kx <- ncol(Xd); G <- crossprod(Xd); Xty <- crossprod(Xd, y)
    if (ridge) {
      # The ridge's prior precision on the centred design [1, Z_c]: t0 diag(1, d), d = colSums(Z_c^2)/N, the intercept's
      # flat under the flat intercept.
      Zc <- Z - rep(colMeans(Z), each = N)
      Xd <- cbind(1, Zc); G <- crossprod(Xd); Xty <- crossprod(Xd, y)
      P0 <- diag(t0 * c(1, colSums(Zc^2) / N), kx)
      if (flat) P0[1, 1] <- 0
      mu0 <- c(if (flat) 0 else level, rep(0, kx - 1))
      prior_ld <- function(s2v) {
        keep <- if (flat) -1 else seq_len(kx)
        if (length(diag(P0)[keep]) == 0) return(0)
        0.5 * sum(log(diag(P0)[keep] / s2v)) - length(diag(P0)[keep]) / 2 * log(2 * pi) -
          0.5 * c(t(mu0) %*% P0 %*% mu0) / s2v
      }
    } else if (flat) {
      P0 <- matrix(0, kx, kx)
      if (length(S)) P0[-1, -1] <- t0 * crossprod(Z - rep(colMeans(Z), each = N))
      mu0 <- rep(0, kx)
      prior_ld <- function(s2v) if (length(S)) 0.5 * c(determinant(P0[-1, -1, drop = FALSE] / s2v)$modulus) -
        (kx - 1) / 2 * log(2 * pi) else 0
    } else {
      P0 <- t0 * G
      mu0 <- c(level, rep(0, kx - 1))
      prior_ld <- function(s2v) 0.5 * c(determinant(P0 / s2v)$modulus) - kx / 2 * log(2 * pi) -
        0.5 * c(t(mu0) %*% P0 %*% mu0) / s2v
    }
    f <- function(v) vapply(v, function(vv) {
      s2v <- exp(vv)
      A <- (b * G + P0) / s2v; h <- (b * Xty + P0 %*% mu0) / s2v
      ll <- -(N * b / 2) * log(2 * pi * s2v) - b * sum(y^2) / (2 * s2v) + kx / 2 * log(2 * pi) -
        0.5 * c(determinant(A)$modulus) + 0.5 * c(crossprod(h, solve(A, h))) + prior_ld(s2v)
      ll + (if (jeffreys) -vv else g1 * log(g2) - lgamma(g1) - (g1 + 1) * vv - g2 / s2v) + vv
    }, 0)
    grid <- seq(log(s2) - 60, log(s2) + 40, length.out = 401)
    peak_at <- grid[which.max(f(grid))]; peak <- f(peak_at)
    pieces <- seq(peak_at - 25, peak_at + 25, length.out = 11)
    total <- sum(vapply(1:10, function(i) integrate(function(v) exp(f(v) - peak), pieces[i], pieces[i + 1],
                                                    rel.tol = 1e-12, abs.tol = 0, subdivisions = 1000L)$value, 0))
    peak + log(total)
  }
  quadrature <- NA
  if (!constant) {
    reps <- vapply(0:M_max, function(M) { i <- which(M_of == M); if (length(i)) i[which.max(pst[i])] else NA_integer_ }, 0L)
    reps <- reps[!is.na(reps)]
    tau_points <- if (fixed_g) opt$tau else c(0.1, 1, 10) / n0
    bs <- Filter(defined, c(1, 0.5))
    sigma2_gap <- max(unlist(lapply(reps, function(i) lapply(bs, function(b) vapply(tau_points, function(t0)
      abs(lp(b, t0, i)[1, 1] - explicit_lpb(state_sets[[i]], b, t0)), 0)))))
    tau_gap <- 0
    if (!fixed_g) tau_gap <- max(unlist(lapply(reps, function(i) lapply(bs, function(b) {
      g <- function(uu) c(lp(b, exp(uu), i)) + log_density_u(uu)
      peak <- max(g(u))
      pieces <- seq(lo, hi, length.out = 41)
      total <- sum(vapply(1:40, function(j) integrate(function(uu) exp(g(uu) - peak), pieces[j], pieces[j + 1],
                                                      rel.tol = 1e-13, abs.tol = 0, subdivisions = 1000L)$value, 0))
      abs(row_lse(sweep(lp(b, tau, i), 2, logw, "+")) - (peak + log(total)))
    }))))
    quadrature <- list(states = I(as.integer(reps - 1L)), sigma2_gap = sigma2_gap, tau_gap = tau_gap)
  }

  strip <- function(r) r[setdiff(names(r), c("log_marginal", "W"))]
  cs_term <- function(t) list(cont = lapply(seq_len(t$K), function(f) as.integer(c(t$vars[f] - 1, t$signs[f], t$knot_index[f] - 1))),
                              cat = lapply(seq_along(t$cat_vars), function(l) list(as.integer(t$cat_vars[l] - 1),
                                                                                   I(as.integer(t$cat_subsets[[l]] - 1)))),
                              func = list())
  result$terms <- lapply(terms, cs_term)
  result$states <- lapply(state_sets, function(S) I(as.integer(S - 1L)))
  result$n_states <- n_states
  result$posterior <- list(
    supported = supported, data = if (n_cat > 0) list(y = y, x = X, z = Z_cat) else list(y = y, x = X),
    response = list(centre = centre, level = level, s2 = s2, sigma2_shape = g1, sigma2_rate = g2,
                    n_eff = N - flat, absolute = absolute),
    form = opt$form,
    categorical = if (n_cat > 0) list(count = n_cat, levels = lapply(cat_levels, I), L_max = L_max, subsets = opt$subsets) else NA,
    tau_prior = if (fixed_g) list(fixed = opt$tau) else list(shape = a_tau, rate = r_tau),
    n_gate_refused = sum(!keep), n_exact_fit_refused = sum(exact_fit),
    M = M_of, k = k, rss = rss, Q = Q, log_prior = log_prior_F, log_marginal = logm, posterior = pst,
    tau_mean = tau_given_S,
    count = list(M = 0:M_max, prior = count_law(exp(log_prior_F)), posterior = count_law(pst)),
    order = list(types = types, levels = lapply(strsplit(order_levels, ","), function(v) I(as.integer(v))),
                 prior = order_law(exp(log_prior_F)), posterior = order_law(pst)),
    inclusion = incl, tau = tau_summary, log_Z = post$log_Z, quadrature = quadrature,
    tempered = list(rungs = rungs,
                    integrated = lapply(rungs, function(b) strip(rung_at(b, "integrated"))),
                    linear = lapply(rungs, function(b) strip(rung_at(b, "linear")))))
}

emit_json(result, args[2])
