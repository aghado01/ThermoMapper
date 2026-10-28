# diagnostics_oracle.R — reference convergence and PSIS-LOO diagnostics (OQ9, N13).
#
# Cross-checks Maths.Samplers.Ensemble.ChainDiagnostics (rank-normalized split R-hat, bulk and tail ESS, ESS of the
# mean) against the posterior package, and Maths.Regression.PredictiveChecks.PsisLoo (k-hat, r_eff, the elpd's Monte
# Carlo error) against the loo package.
#
# Usage (cwd = r/):
#   Rscript oracles/diagnostics_oracle.R <draws.csv> <out.json> chains
#     draws.csv: iterations x chains. Emits rhat, rhat_bulk, rhat_tail, ess_bulk, ess_tail, ess_q05, ess_q95, ess_mean,
#     plus rhat_tail_split and ess_tail_split: the tail diagnostics with the median and the 5% and 95% quantiles taken
#     over the split chains (an odd middle draw dropped), the convention of the C# port. posterior takes them over all
#     draws; the two agree whenever the chains have an even length.
#   Rscript oracles/diagnostics_oracle.R <loglik.csv> <out.json> loo <chains>
#     loglik.csv: draws x units, each chain's draws contiguous; chains = 0 for independent draws. Emits elpd, se, p_loo,
#     the pointwise elpd_i, mcse_i and pareto_k, r_eff, threshold, and mcse (null when a k-hat is above the threshold).

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: diagnostics_oracle.R <fixture.csv> <out.json> chains|loo [chains]")
X <- read_matrix(args[1])
mode <- args[3]

if (mode == "chains") {
  # rhat.default's two halves, as it computes them.
  split <- posterior:::.split_chains(X)
  rhat_bulk <- posterior:::.rhat(posterior:::z_scale(split))
  rhat_tail <- posterior:::.rhat(posterior:::z_scale(posterior:::.split_chains(posterior:::fold_draws(X))))

  # The C# convention: the median and the tail quantiles over the split chains.
  folded_split <- abs(split - median(split))
  rhat_tail_split <- posterior:::.rhat(posterior:::z_scale(folded_split))
  ess_q_split <- sapply(c(0.05, 0.95), function(p) posterior:::.ess(split <= quantile(split, p)))

  q <- posterior::ess_quantile(X, c(0.05, 0.95))
  emit_json(list(
    rhat            = posterior::rhat(X),
    rhat_bulk       = rhat_bulk,
    rhat_tail       = rhat_tail,
    ess_bulk        = posterior::ess_bulk(X),
    ess_tail        = posterior::ess_tail(X),
    ess_q05         = unname(q[1]),
    ess_q95         = unname(q[2]),
    ess_mean        = posterior::ess_mean(X),
    rhat_tail_split = rhat_tail_split,
    ess_tail_split  = min(ess_q_split)
  ), args[2], na = "null")
} else if (mode == "loo") {
  chains <- if (length(args) >= 4) as.integer(args[4]) else 0L
  S <- nrow(X)
  r_eff <- if (chains > 0) {
    loo::relative_eff(exp(X), chain_id = rep(seq_len(chains), each = S / chains))
  } else {
    rep(1, ncol(X))
  }
  fit <- loo::loo(X, r_eff = r_eff)
  pw <- fit$pointwise
  emit_json(list(
    elpd      = fit$estimates["elpd_loo", "Estimate"],
    se        = fit$estimates["elpd_loo", "SE"],
    p_loo     = fit$estimates["p_loo", "Estimate"],
    elpd_i    = unname(pw[, "elpd_loo"]),
    mcse_i    = unname(pw[, "mcse_elpd_loo"]),
    pareto_k  = unname(pw[, "influence_pareto_k"]),
    r_eff     = unname(r_eff),
    threshold = loo:::ps_khat_threshold(S),
    mcse      = loo:::mcse_loo(fit)
  ), args[2], na = "null")
} else {
  stop("mode must be 'chains' or 'loo'")
}
