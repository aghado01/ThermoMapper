# bass_enumeration_oracle.R — structure-space enumeration for BASS's birth, death and change moves (port brief, fitting
# segment): the reference for R's tempering and swap density (D05; F-O6 (a)) and for the count law each birth ratio
# implies (D01; F-O6 (c)).
#
# Usage (cwd = r/):  Rscript oracles/bass_enumeration_oracle.R <data.csv> <out.json> <mode> [key=value ...]
#   data: headerless CSV, y in the first column, then the inputs x_1..x_p in [0, 1] (R's scaled inputs).
#   mode: swap | count.
#   keys (defaults): M_max (2), K_max (2), npart (3), w1 (5), w2 (5), lambda (2), sigma2 (0.1), tau (1/N),
#     rungs (1,0.5,0.2), g1 (0), g2 (0), h1 (10), h2 (10).
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
# Emits { mode, n, p, n_terms, n_states, ... } with the mode's fields; states are listed as sorted term indices into
# `terms` (K, vars, signs, knots).

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: bass_enumeration_oracle.R <data.csv> <out.json> <swap|count> [key=value ...]")
D <- read_matrix(args[1])
mode <- args[3]
y <- D[, 1]; X <- D[, -1, drop = FALSE]
N <- length(y); p <- ncol(X)
opt <- list(M_max = 2, K_max = 2, npart = 3, w1 = 5, w2 = 5, lambda = 2, sigma2 = 0.1, tau = 1 / N, rungs = "1,0.5,0.2",
            g1 = 0, g2 = 0, h1 = 10, h2 = 10)
for (kv in args[-(1:3)]) {
  parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
  opt[[parts[1]]] <- if (parts[1] == "rungs") parts[2] else as.numeric(parts[2])
}
M_max <- opt$M_max; K_max <- min(opt$K_max, p)
bass <- bass_source(c("basis_funcs.R", "mcmc_funcs.R", "rjmcmc_des.R"))

# ---- terms -------------------------------------------------------------------------------------------------------
knot_values <- lapply(seq_len(p), function(v) unique(X[, v]))
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
supported <- if (mode == "swap") colSums(basis != 0) >= opt$npart else rep(TRUE, n_terms)
log_prior_term <- vapply(terms, function(t) -log(K_max) - lchoose(p, t$K) - t$K * log(2) - sum(log(n_v[t$vars])), 0)

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
move_probs <- function(M) {
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
tr <- do.call(rbind, lapply(seq_len(n_states), rows_of))
from <- as.integer(tr[, 1]); to <- as.integer(tr[, 2]); prob <- tr[, 3]
type <- type_names[tr[, 4]]; term <- as.integer(tr[, 5]); lq_fwd <- tr[, 6]

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
  if (mode == "count") {
    A <- t(P) - diag(n_states); A[n_states, ] <- 1
    s <- solve(A, c(rep(0, n_states - 1), 1))
    return(s / sum(s))
  }
  stationary_gth(P)
}
normalize_log <- function(l) { w <- exp(l - max(l)); w / sum(w) }
tv <- function(a, b) 0.5 * sum(abs(a - b))

# The prior-and-proposal part of each move's log ratio, untempered, with the state's actual move probabilities:
# birth lambda pi(t) d_{M+1} / ((M + 1) b_M q(t)); death its reciprocal from the smaller state; change 0.
prior_proposal <- function(extra_birth = function(M) 0) {
  vapply(seq_along(from), function(r) {
    if (type[r] %in% c("self", "change")) return(0)
    Mf <- M_of[from[r]]
    if (type[r] == "birth") {
      mp_to <- move_probs(Mf + 1)
      return(log(opt$lambda) + log_prior_term[term[r]] + log(mp_to[["death"]] / (Mf + 1)) -
               log(move_probs(Mf)[["birth"]]) - lq_fwd[r] + extra_birth(Mf))
    }
    mp_to <- move_probs(Mf - 1)
    -(log(opt$lambda) + log_prior_term[term[r]] + log(move_probs(Mf)[["death"]] / Mf) - log(mp_to[["birth"]]) -
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

emit_json(result, args[2])
