# BASS port: shared infrastructure not yet built (chip brief)

Date: 2026-10-27. Issued from the BASS planning session, which continues the D01–D06 design discussion while this runs. A second chip, [reference-fixtures-brief.md](reference-fixtures-brief.md), runs at the same time in the same checkout and works in R only; the file ownership below keeps the two apart. This chip is the only one that runs `dotnet`.

## Goal

Three shared pieces that BASS's decision drafts lean on and nothing in the repo provides yet. Each belongs below BASS (Samplers, the shared g-prior block, the R oracle bridge), serves existing engines too, and does not wait on the design pass. Done means each has landed with its tests, its numerics-checklist item (N11–N13 in [port-numerics-checklist.md](port-numerics-checklist.md)) is ticked with a report, the docs describe it, and this brief's Report section says what landed.

Suggested order: ask for S3's download approval first (the user can answer while the rest runs), then S2 (small), then S1.

## S1 (N11): evidence over a tempering ladder

D03's draft reports evidence from proper-prior runs by "the SMC normalizer or path sampling over a tempering ladder", and D05's says "Evidence: path sampling over the same rungs or stages, with ∂_b ℓ_b in closed form". The SMC core already estimates its normalizer (`SmcRun`'s `LogNormalizerRatio`). Nothing estimates evidence from a ladder.

- **What to build, in `Maths.Samplers`:** an estimator of log Z₁ − log Z₀ from per-rung draws on a ladder 0 = b₀ < … < b_K = 1:
  - *path sampling* (thermodynamic integration; Gelman and Meng 1998; Friel and Pettitt 2008), from each rung's draws of ∂_b log q_b;
  - *stepping-stone* (Xie, Lewis, Fan, Kuo and Chen 2011), from each rung's draws of log q_{b_{k+1}} − log q_{b_k};
  - each with a Monte Carlo error that uses the effective sample size of each rung's trace (reuse the repo's ESS code; do not write another), and a discretization diagnostic for path sampling (for example Friel, Hurn and Wyse 2014's corrected trapezoid, which uses each rung's variance).

  It knows nothing about models: it takes the ladder and the per-rung scalar traces. Make its input fit what `ParallelTempering` and the SMC core produce.
- **Per-rung readouts.** The estimator needs a scalar trace from every rung, not just the cold one. If `ParallelTempering` does not record them, add a generic, opt-in per-replica readout recorded after burn-in (the target supplies the scalar). It stays model-neutral; the boundary rules are in the [hub](../bass-full-port-brief.md#placement-and-dependency-boundaries).
- **Which estimator is the default.** Make the call with a stated rationale from the literature and from your tests, keep the other selectable, and record the alternatives and trade-offs in `docs/Maths/Samplers.md`. No default from a partial study: if the tests do not separate the two, say so and pick on the literature's grounds.
- **Tests, in the Samplers test host, on targets independent of Regression:**
  - a conjugate Normal prior and likelihood, whose tempered normalizer is closed-form at every b, sampled exactly at each rung, so the estimator's bias and its reported error can be read separately from sampler error;
  - an enumerable discrete target, sampled by a chain, so the effective-sample-size correction is exercised;
  - the discretization error falls as rungs are added, and the reported Monte Carlo error covers the truth at about its nominal rate across seeds.
- **Out of scope:** wiring an engine to report evidence. Normal BARS could, under its proper default σ² prior and its exact reference rung, but whether and how its result carries evidence is a consumer decision. Name it in the report as a follow-up.

## S2 (N12): the g-prior's log marginal at (τ, b) as one scalar

`CollapsedGPrior.At(τ, b)` builds a `GPriorConditional`, which allocates the p-vector of conditional means on every call. The log marginal needs only S, the determinant term and the row count. BASS's τ slice sampler (D04) and every tempered swap (D05) evaluate the marginal many times per sweep and never read the mean.

- **What to build:** a path that returns the log marginal at (τ, b) under a σ² prior without allocating the mean. Unweighted only; the weighted block keeps its single (τ, b = 1).
- **Callers:** move the existing callers that use `At(…)` only for the log marginal onto it (the BARS tempered marginal and its relatives in `src/maths/regression/spline/bars/`).
- **Check:** the new path returns exactly the same double as `At(τ, b).LogMarginal(prior)` over Jeffreys, `ScaledToResponse()` and a fixed IG prior, b ∈ {0, 0.3, 1} (b = 0 under proper priors only), and τ from 10⁻⁶ to 10³. The Regression suite stays green (426 tests at `663ae7f`, plus yours).

## S3 (N13): an external reference for the rank-normalized and PSIS diagnostics (OQ9)

The rank-normalized R̂ and the bulk and tail ESS (step 3.4) are checked against theory only, and so are the PSIS k̂, r_eff and Monte Carlo error built for OQ17 (`PredictiveChecks`), whose algorithms were matched to R `loo`'s source by reading it. D10 makes automatic stopping and performance claims wait on an external reference. The row is OQ9 in the [preparation brief](../bass-preparation-brief.md).

- **Approval first.** R `posterior` and `loo` are not in the project's renv. Installing them is a CRAN download. Before anything else in this chip, ask the user in your own session, naming the packages and their dependencies, the source (CRAN through renv) and the approximate download size, and wait for a clear yes. If the answer is no, record OQ9 as still open and skip S3.
- **What to build, on approval:** add the packages to `r/DESCRIPTION` and snapshot `r/renv.lock` (see `r/README.md`); a new oracle `r/oracles/diagnostics_oracle.R` that emits R̂, bulk ESS and tail ESS (`posterior`) and k̂, r_eff and the elpd Monte Carlo error (`loo`) for supplied draws; and parity tests in `tests/oracle/` through `ROracle.cs`, in the style of the existing parity tests. Fixtures: AR(1) chains with known autocorrelation, chains with mismatched scales, heavy tails (the existing theory checks' cases), and a pointwise log-density matrix with one influential unit.
- **Where they disagree:** find out why before changing either side. A difference in convention (split chains, folding, the tail quantiles, the tail length) is recorded with its source, not silently aligned.

## How to work

- Read [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md) first. Generic machinery stays in Samplers or LinAlg and never inspects model meaning. Docstrings stay lean.
- `docs/Maths/*.md` use CRLF line endings; briefs, checklists and C# use LF. For multi-line splices into the docs, the repo's habit is a Python script with raw strings, checked with `git diff --stat`.
- Run the Samplers and Regression suites after S1 and S2. Run the oracle tests after S3 (cwd and R activation as in `r/README.md`).

## File ownership

This chip writes only:
- `src/maths/samplers/**`, `src/maths/regression/gprior/CollapsedGPrior.cs`, and the BARS files whose calls move onto S2's path;
- tests under the Samplers, Regression and Oracle test roots;
- for S3, `r/DESCRIPTION`, `r/renv.lock`, and new files `r/oracles/diagnostics_*.R` and `r/scripts/diagnostics_*.R`;
- `docs/Maths/Samplers.md`, and `docs/Maths/BayesianRegression.md` only where S2 needs a line;
- [port-numerics-checklist.md](port-numerics-checklist.md) (N11–N13 and their reports, prepended under "## Reports"), OQ9's row in the preparation brief, and this brief's Report section.

It does not edit `r/oracles/_common.R` or `r/scripts/bass_stage1_checks.R` (the other chip's), the hub or the segment briefs (under discussion in the planning session; put any finding that bears on a draft in your report instead), or the fitting, outputs and calibration checklists.

## Commits

Work on `main` with targeted commits of your own paths only: `git add <your new files>`, then `git commit -m "..." -- <your paths>`, so the other chip's staged work is never swept in. Never `git add -A` or `git add .`. If `.git/index.lock` exists, the other session is committing: wait and retry. Leave the user's untracked notes alone (`issues/regression/clumped-problems.md`, `issues/regression/error-capture-proposal.md`, `issues/spc/opus-SMC-notes.md`, `opus-gvdv-ideas.md`, `issues/spred/pca-degeneracy.md`). End each commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Report

Append the report here when the work lands, newest first: date, commits, test evidence, the default chosen in S1 and why, S3's outcome (or why it was skipped), limits, and follow-ups.
