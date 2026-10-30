# BASS port: shared numerics and sampler contracts: checklist

Companion to [port-numerics-brief.md](port-numerics-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions

- [x] **N-D10** Adaptation policy for adaptive kernels (note): adopted 2026-10-27 in the user's design pass; D10's stopping semantics was drafted and decided on 2026-10-29 (the brief's D10 section).

### Build (stages 2–3)

- [x] **N1** Bordered insertion on the strict dense factor, with candidate work private until acceptance: `BorderedCholesky`.
- [x] **N2** Givens deletion, and a changed term as a staged delete and insert: `BorderedCholesky`.
- [x] **N3** Strict rebuild: periodic and on failure, with drift checks against fresh factorizations: `BorderedCholesky.Rebuild`, `FactorDrift`.
- [x] **N4** diag(G⁻¹) maintained through births and deaths, for D06's order-free gate: `BorderedCholesky.LeaveOneOutPivot`.
- [x] **N5** Stable thin SVD with numerical-rank handling in LinAlg: `ThinSvd`.
- [x] **N6** Diagonal-plus-low-rank solves and log-determinants, without m × m inverses: `DiagonalPlusLowRank`.
- [x] **N7** Face-splitting design contraction with the Cartesian (Hadamard) shortcut: `KhatriRao`, both intercept conventions.
- [x] **N8** Generic SMC core in `samplers/smc` (optional; specified in the calibration segment): `SmcRun`, `SmcSampler`.
- [x] **N9** Qualified baseline rerun after the shared changes, before the first BASS integration tests: requalified at `3a04601`, and again at `cf37584` before stage 4 (2026-10-28).
- [x] **N10** Cross-products summed pairwise or with compensation: `PairwiseSum`. F-O7's near-constant design is the acceptance case: its left-to-right error is 10⁻³ relative in s.
- [x] **N11** Evidence over a tempering ladder: path sampling and stepping-stone from per-rung traces, with Monte Carlo error (D03, D05). Chip: [shared-infrastructure-brief.md](shared-infrastructure-brief.md). `LadderEvidence` with stepping-stone as the default; `ParallelTempering.TraceRungs` (`84c173f`).
- [x] **N12** The g-prior's log marginal at (τ, b) as one scalar, without the mean vector (D04's τ slice, D05's swaps). Same chip. `CollapsedGPrior.MarginalAt` → `GPriorMarginal`, bit-identical to `At`'s (`0afe33d`).
- [x] **N13** External reference for the rank-normalized and PSIS diagnostics against R `posterior` and `loo` (OQ9; D10's gate). Same chip; needs the user's approval for the download. Approved; `diagnostics_oracle.R` and `DiagnosticsParityTests` agree to 10⁻⁹ (`5eecbc1`, `ee8b97f`).
- [x] **N14** D10's convergence verdict: read at the end of every run in the shared ensemble, with the engine's agreement groups and a declared stopping scope; BASS's three groups, the log posterior among its monitored quantities (decided 2026-10-29).
- [ ] **N15** `DiagonalPlusLowRank` gains the capacitance matrix's Cholesky factor, O(m·k² + k³), beside its per-call Jacobi thin SVD. Found by stage 7 on 2026-10-29: 1.4 ms per σ² proposal at m = 20 with 39 columns. Built when real sizes make it binding.
- [ ] **N16** The SMC core's conveniences, both stage 8 chips' notes (2026-10-29): `TemperToPosterior` returning at its stage budget rather than throwing, `SmcRun` exposing its `MaxStages`, and a first-step bisection on cached increments in `SmcWeights` beside `ConditionalEss`.

### Room for extensions

- [ ] **N-X1** The face-splitting contraction covers tensor designs (ABMRS).
- [ ] **N-X2** `ISpdFactor` admits a wide-band factor beside `DenseCholesky`.
- [ ] **N-X3** `ConjugateUpdates` expressible as limits of a GIG update.

## Reports

### N14 verified: the convergence verdict — 2026-10-29

The N14 chip built it ([its report](n14-convergence-verdict-brief.md#reports)), and it was read and verified here.
- *Commits:* `e8547d3` (the shared verdict, committed first for stage 7), `e88b33e` (BASS), `1c26dcd` (BARS, BGP, BAPS, `HierarchicalPSpline`), `6fe653a` (the report).
- *Read here:* one rule, `ConvergenceVerdict.Of`, called by the stop between rounds and by the verdict at the end, with `MayStop` deleted. Under the default single group it is the old rule.
- *Suites rerun here:* Samplers 204/204, Regression 1,085/1,085, the BASS oracle classes 63/63 (with stage 7's work in progress in the tree). The eight bit-identity hashes hold.
- *As built against the brief:*
  - frozen and unavailable readings fail, as D10's rules say; the brief's API bullet that called them undetermined was wrong;
  - the log posterior equals the chain's own density bit for bit;
  - BARS declares two groups, predictive and structural, and its stop is unchanged;
  - BGP, BAPS and `HierarchicalPSpline` carry the verdict on one group.
- *Changed behaviour:*
  - every BASS fingerprint moved, since the config's JSON now carries `StoppingScope`, so a checkpoint from before `e88b33e` will not resume (pre-release, no shim);
  - a fit with an R̂ target may stop later, since the log posterior now gates it. `BassOutputComponentsTests` now sets the predictive scope to keep exercising unequal stops.
- *On S10's fits* every group is not converged, with log-posterior gaps of 21.8–42.8 nats untempered and 27.0 tempered on Friedman.
- *Docs:* `docs/Maths/Samplers.md` (the verdict, the API and source rows) and `docs/Maths/BayesianRegression.md` (BASS's diagnostics and configuration, BARS's verdict, the API row).
- *Not built:* a stopping scope on BARS, and any stop on the mode-separation readings.

### D10's stopping semantics decided — 2026-10-29

The user took the draft as written, with its four points:
- the reference target is 1.01;
- the log posterior joins the monitored quantities;
- every group gates the stop by default;
- R's budget stays the provisional default.

N14 builds it. The verdict belongs to the shared ensemble, so BARS gets it with BASS.

### Baseline requalified before stage 4 — 2026-10-28

**Qualified revision:** `cf37584`, N9's baseline (`3a04601`) plus 14 engine commits, all in `regression` and `samplers`; `linalg`, `distributions` and `rng` are unchanged since. The toolchain is N9's: .NET SDK 10.0.201, `net10.0`, Release, on the same machine.

**Suites** (`scripts/fact-harness.ps1 -Project`, each case its own `dotnet test` process):

| Suite | Cases | Passed | Run |
|---|---|---|---|
| `Maths.Regression.Tests` | 438 | 438 | `20261028_013448602` |
| `Maths.Samplers.Tests` | 137 | 137 | `20261028_014037203` |
| `Maths.LinAlg.Tests` | 118 | 118 | `20261028_014345840` |
| `Maths.Distributions.Tests` | 93 | 93 | `20261028_014552996` |

**Downstream consumers** of the changed code: `Synthetic.Tests` 38 of 38 and `Maths.Oracle.Tests` 16 of 16 under plain `dotnet test`; `Clustering.Graphical.SPC` builds. `Maths.Geometry.Tests` reads `linalg` alone, unchanged, and was not rerun.

**The drift record.** Each of the 14 commits has a recorded account of what it changes:
- *One declared path change:* `c7f8626` makes the scaled σ² prior Normal BARS's default ([BayesianRegression](../../../docs/Maths/BayesianRegression.md), the σ² prior). Every default Normal BARS fit's marginal moves by a near-constant amount across configurations, so its accept decisions can differ.
- *Options and readouts that leave the default paths:*
  - `fe70b32` adds the σ² prior with Jeffreys still the default;
  - `92540de` adds exact reference rungs, and ladders without one keep their streams (preparation brief, OQ2);
  - `00eb9da` and `1fdaacf` report the evidence only on a likelihood-path ladder to a reference rung;
  - `128fd1f` reconciles the PSIS and tail-diagnostic conventions; the bulk diagnostics are untouched;
  - `a1971e1` adds the genealogy error as a readout;
  - `cf37584` (F-C2) keeps the σ² prior's factories and its fingerprint spelling.
- *New code beside the engines:* `84c173f` (N11), `0afe33d` (N12, bit-identical to `At`'s), `ae82374` (SMC traces and `FromSmc`).
- *Paths outside the sweep's claims:*
  - `70ea7b5` deletes Poisson's BIC-only code (the score, its ridge, its MLE start), leaving the g-prior path as it was, and no claim uses Poisson;
  - `d798961` makes SMC resampling multinomial by default, which changes SMC runs' draws, and no claim uses SMC.
- `66f6171` changes a doc comment and the bimodal fixture.

**Seed sweep** (`SEED_SWEEP=50`):
- **Hold counts.** Regression's `20261028_015051357` matches N9's `20261027_023326165` claim by claim. Samplers' `20261028_015056249` matches `20261027_023344333` row for row.
- **Per-seed values.** The regression rows changed only in the claims that run Normal BARS:
  - the adaptive kernel: 101 of 101 rows;
  - Student-t weights: 101 of 101;
  - the model-preference claim: 51 of 51;
  - contaminated-normal weights: 65 of 101;
  - EBIC: 15 of 101, its mean knot counts being coarse;
  - one Gibbs-loss row, through its least-squares comparator, which is Normal BARS (MSE 0.0825453 → 0.0825456).
- **Split at `c7f8626`.** To attribute those changes, the sweep was rerun on the source of `c7f8626^` (`3690061`) and of `c7f8626`, extracted with `git archive` (scratch). The comparisons:
  - `3a04601` → `3690061`: no row changed;
  - `3690061` → `c7f8626`: every change above;
  - `c7f8626` → `cf37584`: no row changed.
  So the default σ² prior accounts for all the drift, and the other 13 commits leave every claim's per-seed values exactly as they were.
- **Sizes.**
  - *The adaptive kernel:* the chain axis's worst ratio went from 0.67 to 0.76, and its worst margin from 14.4 to 9.3 in ESS, still 50 of 50.
  - *The robust claims' margins* moved in the fourth digit.
  - *The model-preference claim:* the margin held at 3.05 nats on the data axis, while the median ratio went from 0.934 to 0.969. The complete marginal adds the same constant to both sides, which moves the ratio but not the difference.

**Gate status:** closed at `cf37584`. Stage 4's first integration tests measure against this baseline.

### The SMC normalizer's genealogy error (Lee and Whiteley 2018) — 2026-10-27

At the user's request, after D05's fourth point (the Kish error under-covers at one move per stage).
- **Built.** `SmcWeights.GenealogyRelativeVariance`: 1 − (N/(N−1))^{R+1}·Σ_{E_i ≠ E_j} W_i W_j over the particles' initial ancestors, which `SmcRun` already tracked as families. Every `SmcStageReport` carries it as `NormalizerRelativeVariance`, read after the reweighting and before any resampling; `SmcResult.LogNormalizerError` gives √ln(1 + v); `LadderEvidence.FromSmc`'s stepping-stone takes it in place of the stages' weights.
- **Multinomial only.** The probe across schemes (200 seeds, 1000 particles) found it calibrated under multinomial resampling (coverage 0.950 at one move per stage, 0.965 at five) and too small under systematic, stratified and residual (0.90 at one move, 0.795 at five under systematic). So it is NaN from the stage after the first non-multinomial resampling. The full table is in `Samplers.md`.
- **The kernel.** A population-tuned kernel (`RandomWalkKernel`) adds a dependence the genealogy does not see: 5% low at N = 16, invisible at N = 1000. The unbiasedness check uses a fixed-step kernel.
- **Checks** (`SmcTests`, `LadderEvidenceTests`): a hand case and the one-family case; importance sampling's (N·ΣW² − 1)/(N − 1) without resampling; unbiasedness over 20,000 replicates at N = 16 (0.020 ± 0.013, with the exponent R rejected at 0.081 ± 0.012); coverage at one and five moves under multinomial resampling; NaN after a systematic resampling, kept through a resume taken just after it (the checkpoint carries `GenealogyHolds`, since that stage's own report is still finite). The SMC evidence test now resamples multinomially, so the path-sampling numbers `Samplers.md` quotes moved (−2.273, diagnostic +0.149). Samplers 137/137, Regression 432/432.
- **Not built.** An error for path sampling over SMC stages, which a sum of stage means needs; the fixed-lag variant for long runs whose lineages collapse (Olsson and Douc 2019).
- **The default is multinomial** (the user's call, 2026-10-27, revising the same day's "keep systematic", `173b8b2`). The reason is that a default run should report its normalizer's error, which only multinomial resampling makes estimable, and that the SMC limit theory is multinomial's; systematic lacks the general convergence theory the others have (Gerber, Chopin and Whiteley 2019). On this target the lower-variance schemes did not lower the spread of ln Ẑ, and multinomial costs the same O(N). Systematic, stratified and residual stay selectable for runs that need only posterior expectations. Posterior moments under the two schemes were not measured; a consumer that finds systematic pays in its regime opts in.

### The linear form's path sampling, and the evidence's priors — 2026-10-27

The user's call after the seam was read against D05 (`port-fitting-brief.md`): rehearse the linear-form route on BARS before BASS relies on it, and carry the priors as D05 asks. Commit `1fdaacf`.

- **Path sampling on the linear form.**
  - `ITemperedMarginal.LogMarginalLikelihoodSlope` returns ℓ₁ at every b on the linear form. Its default refuses, and the integrated form keeps that default.
  - Runs on `TemperedForm.Linear` therefore report `BarsEvidence.PathSampling` beside stepping-stone, over the same draws.
  - ℓ_b is linear in b there, so the Friel–Hurn–Wyse diagnostic reads the integrand's true derivative. Path sampling minus stepping-stone reads the trapezoid's error against an estimator that has none.
- **The integrated form stays stepping-stone only.** Its ∂_b ℓ_b has a closed form, but its diagnostic would also need E_b[∂²ℓ_b].
- **Priors carried** (D05). `BarsEvidencePriors` holds:
  - the count prior and the support gate;
  - the centre and its resolved level;
  - the σ² prior as specified and as resolved, through `WeightedNormalModel.ResolvedPriors`.

  Both resolved pieces read the response, so the default evidence is empirical-Bayes. A comparison is licensed only on the same data, where they agree.
- **The gate names the Normal model.** It is now `WeightedNormalModel { HasCompleteMarginal: true }`, since only the Normal model's priors are described. `IObservationModel.HasCompleteMarginal` is gone.
- **Checks** (`BarsEvidenceTests`). One knot, linear form, 10 levels: stepping-stone −1.6313 ± 0.0047 and path sampling −1.6322 ± 0.0049, against a quadrature of −1.6231, with a diagnostic of 0.0010. The test also covers the priors carried, and the integrated form's missing path sampling.
- **Tests.** Regression 432/432.
- **For D05:** folded into the D05 draft's evidence bullet ([port-fitting-brief.md](port-fitting-brief.md#d05-draft-tempering-and-smc-target-paths)) by the planning session, 2026-10-27. Four points: the estimator follows the form; the linear form makes path sampling free; BASS's priors record would add τ's prior, the term measure, M's count prior and the interaction cap; and the SMC error under-covers at one sweep per stage.

### Follow-ups: conventions reconciled, SMC stage traces, BARS evidence — 2026-10-27

The user's follow-up to N11–N13: shore up the gaps, surface and reconcile differences from R instead of leaving them silent, give SMC per-stage readouts, and have Normal BARS report evidence by default behind a knob. Commits: `128fd1f` (conventions), `ae82374` (SMC), `00eb9da` (BARS).

- **posterior's Pareto guard.** `FitGeneralizedPareto` returns (NaN, NaN) when the tail's first quartile is its minimum, as `gpdfit` returns NA, and uses log1p as `gpdfit` does. PSIS reads that as `NoFit`, k̂ = +∞, as loo does. It was +∞ before too, but by NaN arithmetic. The parity `nofit` fixture agrees with loo: +∞ on both sides, and both withhold the total's error.
- **Tail treatment surfaced.** `PredictiveCheck.ParetoTails` gives each unit's `ParetoTail`: `Fitted`, `Tied`, `TooShort` or `NoFit`.
  - *The tied tail keeps k̂ = −∞, after weighing loo's +∞.* A tie means at least M draws hold the largest ratio: the sample's tail is bounded and no draw dominates, which is what k̂ detects. Ties are the norm for the collapsed content. On the three-structure case the unsmoothed estimate is the exact harmonic mean, which loo would refuse.
  - *The divergence is asserted.* The parity `tied` fixture checks loo's +∞ and its NA total against the C#'s −∞ and its finite total.
- **Odd lengths aligned.** The tail R̂'s median and the tail ESS's quantiles are taken over all draws, as posterior and ArviZ take them. At 4 × 999 the C# now matches posterior to rounding (it differed by 2.8·10⁻⁵ before). The oracle's split-chain variant is gone.
- **SMC per-stage readouts.**
  - `SmcOptions.TraceStages` records an `SmcStageTrace` per stage: weights, slopes, and the log ratios to the next stage. They ride in the result, the run and the checkpoint, and resume is exact.
  - `LadderEvidence.FromSmc`: stepping-stone reproduces the run's `LogNormalizerRatio` to 10⁻¹⁰, and path sampling works over the same stages.
  - `RungDraws.FromPopulation` drops zero-weight draws.
  - *Settled by a probe* (200 seeds, 1000 particles). The Kish-based error covers 0.935 (z² 1.14) with 5 moves per stage, and 0.795 (z² 2.15) with 1 move: shared ancestry the weights miss. A genealogy-based estimate (Lee and Whiteley 2018) is the fix, not built.
- **Normal BARS evidence.** `BarsResult.Evidence` (`BarsEvidence`) is on by default (`evidence = true`, `ReportsEvidence`) wherever a run can give one: a likelihood-path ladder to an exact reference rung, and `IObservationModel.HasCompleteMarginal`, the Normal model under a proper σ² prior.
  - *How.* Each chain's ladder traces its rungs after burn-in. The estimate is the mean of the chains' stepping-stone estimates, since each chain adapts its own rungs, with the per-chain estimates kept.
  - *Checks.* With no knots it is the empty design's marginal exactly. With at most one, −1.6165 ± 0.0059 against a quadrature of −1.6231.
  - *The default BARS run is untempered,* and the reference rung is opt-in (OQ2), so the default reports none. Making evidence default for every run would mean tempering every run, which is the user's call.
- **Tests.**
  - Samplers 131/131 and Regression 431/431.
  - Oracle `DiagnosticsParityTests` 11/11 (odd length direct, `tied`, `nofit`).
- **Open.**
  - Path sampling for BARS needs the Normal model's ∂_b ℓ_b, which has a closed form.
  - D05's diagnostic needs E_b[∂²ℓ_b]; there is no seam for it yet.
  - The SMC genealogy variance.
  - The `r/README.md` oracle row.

### External reference for the diagnostics (N13, OQ9) — 2026-10-27

From the shared-infrastructure chip. The user approved the CRAN download. Commits: `5eecbc1` (environment), `ee8b97f` (oracle, parity tests, docs, OQ9 row).

**Environment.** R `posterior` 1.7.0 and `loo` 2.10.1 are in `r/DESCRIPTION` and `r/renv.lock`, with eight dependencies not pinned before: abind, distributional, matrixStats, numDeriv, pillar, tensorA, tibble and utf8. They are Windows binaries, 6.1 MB, and the lockfile diff is those ten packages alone.

**Oracle.** `r/oracles/diagnostics_oracle.R` has two modes.
- *`chains`* emits posterior's R̂, bulk and tail R̂, bulk ESS, tail ESS (with q05 and q95) and `ess_mean`. Bulk and tail R̂ are taken from `rhat.default`'s two halves through posterior's internals. It also emits a split-chain variant of the tail diagnostics, for the convention below.
- *`loo`* emits elpd, SE, p_loo, per-unit elpd, k̂ and Monte Carlo error, `relative_eff` (from `exp(ll)` and the chain ids), the threshold and `mcse_loo`.

**Parity** (`tests/oracle/DiagnosticsParityTests.cs`, 9 cases; Oracle suite 14 of 14). The tolerance is 10⁻⁹ relative, with k̂ also allowed 10⁻⁹ absolute; observed differences are of order 10⁻¹⁴ relative.
- *Convergence:* AR(1) at φ = 0.9 and −0.3, one chain at three times the scale, one shifted by half a standard deviation, Cauchy draws, and AR(1) at an odd length (4 × 999).
- *PSIS-LOO:* 30 Normal units with one influential unit, 4000 draws of (μ, σ).
  - Independent draws: largest k̂ 0.58.
  - Four AR(1) chains: r_eff below 1, k̂ 0.59, and the total's Monte Carlo error compared.
  - An outlier at 9: k̂ 1.45; both withhold the total's error.

**Where they differ, and why.**
- **Tail-diagnostic convention.** `posterior` takes the median (`fold_draws`) and the 5% and 95% quantiles (`.ess_quantile`) over all draws. The C# takes them over the split chains, which drop an odd middle draw.
  - At even lengths the two coincide.
  - At 4 × 999 the tail R̂ is 1.0121491 here and 1.0121767 in `posterior`. The C# matches the oracle's split-chain variant to 3·10⁻¹⁶.
  - It is kept and recorded in `docs/Maths/Samplers.md`, not aligned.
- **Not exercised by these fixtures, read from the sources.**
  - Constant draws are NA in `posterior` and `Constant` here.
  - A PSIS tail that is flat to rounding: `loo` leaves k̂ = +∞ and warns, while the C# sets −∞ and leaves the tail unsmoothed, as documented in `PredictiveChecks`.
  - `posterior::gpdfit` returns no fit when the tail's first quartile equals its minimum; the C# has no such guard.

**Limits.**
- R stays opt-in: without the toolchain the parity tests return early.
- `r/README.md`'s oracle map has no row for this oracle yet. That file is outside this chip's ownership, so the row is a follow-up.

### Ladder evidence (N11) and the g-prior's mean-free marginal (N12) — 2026-10-27

From the shared-infrastructure chip ([shared-infrastructure-brief.md](shared-infrastructure-brief.md)). Commits: `0afe33d` (N12), `84c173f` (N11).

**N12: `CollapsedGPrior.MarginalAt(τ, b)`.** It returns a `GPriorMarginal` value holding S, the determinant term and m, with `LogMarginal(prior)` and `Sigma2Posterior(prior)`, and it allocates nothing.
- `At(τ, b)` builds its conditional from the same value, and `GPriorConditional`'s `S`, `LogDeterminantTerm`, `LogMarginal` and `Sigma2Posterior` read through it, so the two paths cannot drift apart. The domain is `At`'s: unweighted at any (τ, b), weighted at its own τ and b = 1.
- Normal BARS's tempered marginal (`WeightedNormalModel.NormalMarginal`, at b = 1 and at every rung) now reads `MarginalAt`. `NormalConditional` keeps `At`, since it needs the mean.
- **Evidence.**
  - A before/after dump of 1,233 values was bit-identical: block and BARS marginals over 3 centres, 3 knot sets, 4 σ² priors, τ ∈ [10⁻⁶, 10³] and b ∈ {0, 10⁻⁶, 0.3, 0.7, 1}, the weighted b = 1 and the exception types included.
  - The new test checks `MarginalAt` against `At`'s log marginal bit for bit under the Jeffreys, `ScaledToResponse()` and IG(2, 0.05) priors, at b ∈ {0, 0.3, 1} and 11 values of τ from 10⁻⁶ to 10³, and on the weighted block. It also checks zero allocation (`GC.GetAllocatedBytesForCurrentThread`).
  - Regression: 427 of 427.

**N11: `Maths.Samplers.Evidence.LadderEvidence`.** It estimates ln Z(b_max) − ln Z(b_min) from per-rung draws (`RungDraws`): ordered chains, whose ESS is `ChainDiagnostics.MeanEss`'s, or a weighted SMC population, whose ESS is Kish's.
- *Stepping-stone* (the default) and *path sampling* (selectable), with a Monte Carlo error. Path sampling also carries Friel, Hurn and Wyse's correction as its discretization diagnostic.
- *The error* is the larger of the rungs-as-independent reading and the joint reading, which sums each scan across rungs (rationale and calibration table in `docs/Maths/Samplers.md`, "Evidence over a ladder").
- *The target contract.* `ITargetEvaluation.LogLikelihoodSlope(b)` is a default interface member that refuses. `SeparableEvaluation` returns ℓ₁, and `Temperature.LogDensitySlope` gives ∂_b ln q_b on either path.
- *Per-rung readouts.* `ParallelTempering.TraceRungs()` is opt-in, after burn-in. Each scan it records every rung's slope and its log ratio to the next colder rung, from the evaluations the rungs already hold. `RungTraces` and `TracesSlopes` expose them, the snapshot carries them, and `AdaptLadder` refuses while tracing. `FromLadders(ladders, method)` reads one or more traced ladders.
- **Evidence** (`LadderEvidenceTests`, 7 tests; Samplers 129 of 129, Regression 427 of 427):
  - *Exact draws, conjugate Normal* (16 intervals of (k/K)⁵, 1000 draws per rung, 400 replicates):
    - stepping-stone: coverage 0.963, mean z² 0.98, bias −0.0023 against an error of 0.030;
    - path sampling: coverage 0.953, mean z² 0.96 against the exact trapezoid, whose own bias is −0.084.
  - *Discretization.* The trapezoid's error falls about 4× per doubling from 4 to 64 intervals. The diagnostic leaves 0.18× the error at 8 intervals and 0.003× at 64, and overshoots at 4.
  - *Traced ladders, enumerable target* with a likelihood that is not linear in b (`ReferenceLadder(8, 0.01)`, two ladders × 3000 scans, slope ESS down to 1–3% of the draws). At 200 replicates:
    - stepping-stone: coverage 0.955, mean z² 0.90;
    - path sampling: coverage 0.975, mean z² 1.01.
  - *Each reading alone,* at 400 replicates: the independent one under-reads path sampling (z² 1.23) and the joint one stepping-stone (z² 1.22); the larger reads 0.95 and 1.06.
- **Default: stepping-stone.** It reads ℓ_b only, which every target supplies, and it has no discretization bias: −0.004 against the trapezoid's +0.669 on the nonlinear target. It is also the SMC normalizer's estimator. The trade-offs are documented.
- **Limits.**
  - The diagnostic reads V_k as the mean slope's derivative, so on a collapsed likelihood it misses E_b[∂²ℓ_b]. On the test's nonlinear target it read 0.094 against a bias of 0.669.
  - A population's Kish ESS ignores shared ancestry.
  - `SmcRun` exposes no per-stage readouts.
- **Follow-ups.**
  - Whether Normal BARS reports an evidence under its proper σ² prior and exact reference rung is a consumer decision. Its evaluations supply no slope, so it would use stepping-stone, which `MarginalAt` makes O(1) per rung.
  - For D05's path sampling, BASS should supply ∂_b ℓ_b and, for a valid diagnostic, E_b[∂²ℓ_b] (the latter has no seam yet).

### Baseline requalified (N9) — 2026-10-27

**Qualified revision:** `3a04601`, the step 3.5 baseline (`2efa134`) plus 84 engine commits. The toolchain is step 3.5's: .NET SDK 10.0.201, `net10.0`, Release, on the same machine.

**Suites** (`scripts/fact-harness.ps1 -Project`, each case its own `dotnet test` process):

| Suite | Cases | Passed | Run |
|---|---|---|---|
| `Maths.Regression.Tests` | 406 | 406 | `20261027_020935664` |
| `Maths.Samplers.Tests` | 119 | 119 | `20261027_021721150` |
| `Maths.LinAlg.Tests` | 118 | 118 | `20261027_021924460` |
| `Maths.Distributions.Tests` | 93 | 93 | `20261027_022127437` |

**Downstream consumers:**
- `Maths.Geometry.Tests` passes 51 of 51 and `Synthetic.Tests` 38 of 38, under plain `dotnet test`.
- `Clustering.Graphical.SPC` builds.
- Earlier the same day the Pca migration also ran Mapper 28/28, Oracle 5/5, VizCore 240 (2 skipped by design) and DimReduction 83/83.

**The drift record.** Every one of the 84 engine commits since `2efa134` (in `regression`, `samplers`, `linalg`, `distributions` and `rng`) has a recorded account of what it changes:
- 76 are named in the preparation brief's or this segment's reports, each with its probe or its declared path change;
- 7 are this segment's N-series, recorded by commit range above, which add new code beside the engines;
- `87397fb` changed only a doc comment in `Pointwise.cs`; its report is in the preparation brief under its title.

**Seed sweep** (`SEED_SWEEP=50`):
- **Hold counts.** Regression's `20261027_023326165` matches the last sweep (`20261022_182007796`) claim by claim, except the adaptive-kernel comparison on the data axis, which went from 49 to 50 of 50. Samplers' `20261027_023344333` is identical to `20261021_195902630`.
- **Per-seed values.** To trace those, the sweep was rerun on the source of `6bc5134^` (`bc51934`, extracted with `git archive`, scratch) and split there:
  - *Oct 22 → `bc51934`:* every `HierarchicalPSpline` row changed, and nothing else. That is `d2db820`'s declared path change, the ensemble shell with overdispersed starts, whose agreement in distribution with the single chain is in its report.
  - *`bc51934` → HEAD:* 20–56% of the rows of the five BARS claims changed, and nothing else: adaptive kernel, EBIC, Gibbs loss, Student-t and contaminated-normal weights. That is `6bc5134`, which sheds knots from a start the posterior refuses. The later commits in that span change a doc comment, a PSIS readout and code outside the BARS sampler.
  - *Unchanged throughout:* BAPS, the anisotropic and tensor splines, the penalty and REML claims and the changepoint claim.
- **The 49 → 50.** It is seed 46 (data seed 372476911). On Oct 22 the too-tight fixed kernel reported a peak-location bulk ESS of 4,586, the signature of a chain held at a refused start. It now reports 9.7, against the adaptive kernel's 359.

**Unlike step 3.5,** this rerun did not rebuild a field-by-field output probe across revisions. Its drift evidence is each step's own recorded probe, plus the seed sweep's per-seed table split at `6bc5134^`.

**Gate status:** closed at `3a04601`. The first BASS integration tests measure against this baseline.

### The generic SMC core (N8) — 2026-10-27

`3a04601`, in `src/maths/samplers/smc/`. No existing code changed. The calibration segment's list, item by item:
- **Contracts.**
  - *`SmcStage<TState>`* is a target at a temperature, so a tempered run walks one target from b = 0 to 1 and a cut run walks different conditional targets.
  - *`ISmcParticles<TState>`* draws from the first stage, which must be normalized, and forks a particle for resampling's offspring.
  - *`ISmcKernel<TState>`* is a stage-invariant move. It gets a `Prepare` call with the population before the stage's moves, so it can tune from the past while staying fixed during the stage, as D10's policy asks.
  - Densities come through `ISamplingTarget`: moving to another temperature of the same target queries each particle's evaluation instead of repeating its work.
- **Weights** (`SmcWeights`): log-sum-exp normalization, ESS, the conditional ESS (Zhou, Johansen and Aston 2016), and the log-normalizer increment.
- **Resampling:** multinomial, stratified, systematic and residual. A zero-weight particle is never chosen.
- **Driver** (`SmcRun`, stage by stage): reweight, resample when ESS/N falls below a fraction, then mutate, each particle on its own `SeedTree` stream. `SmcSampler` runs a fixed sequence, or tempers adaptively by choosing each next b through a bisection on the conditional ESS.
- **Diagnostics, kept apart from chain R̂ and autocorrelation ESS:** per stage, the ESS, conditional ESS, whether the population resampled, distinct ancestors, surviving families (distinct initial lineages), acceptance, and the log-normalizer increment.
- **Checkpoints** resume the uninterrupted run bit for bit, because every stream is derived by position.
- **`RandomWalkKernel`:** Gaussian random-walk Metropolis for real vectors, scaled per coordinate from the population's spread.

**Tests** (`tests/samplers/SmcTests.cs`, 10):
- the weights at extreme scales;
- each resampling scheme's offspring averaging N·W_i over 20,000 repetitions, with systematic and residual keeping their deterministic parts;
- a conjugate Gaussian: the posterior moments and the closed-form evidence. Over 20 seeds (scratch) the log-evidence error had mean −0.0003 and sd 0.039;
- a bimodal likelihood, whose modes' masses (0.3 and 0.7) are recovered;
- without resampling or moves, the evidence telescoping to plain importance sampling exactly;
- results independent of parallelism;
- checkpoint resume.

**Not done:** BASS's fitting and cut-calibration adapters (stage 8), the cut SMC's bridges, OQ2's exact reference rung for parallel tempering (and the OQ1 rerun it enables), and SPC population annealing.

### N6, N7, and the SVD callers moved — 2026-10-27

Commits `09572c6`..`28fec55`.

- **`DiagonalPlusLowRank`** (N6): Σ = D + U·diag(s)·Uᵀ.
  - The factor of I + WᵀW (W = D^(−½)·U·diag(√s)) is taken as the thin SVD of W. That makes bᵀΣ⁻¹b and ln|Σ| sums of nonnegative terms, and nothing m × m is formed.
  - Solves and draws from N(0, Σ) and N(0, Σ⁻¹) each cost O(m·r).
  - Along a heavy direction (weight 10¹⁶) it keeps the quadratic form to 10⁻¹⁰ relative. The Woodbury subtraction loses all of it; the test checks both.
- **`KhatriRao`** (N7): the direct functional design X = A ⊙ H, through its factors.
  - It provides the Hadamard Gram, birth borders and response cross-products, for R's convention and centred for the Liang convention (the three-term form), so N7 does not wait on the intercept call.
  - It also accumulates over scattered (case, location) samples through a pairwise cascade.
  - *Terminology:* BASS's design is the column-wise Khatri–Rao product. Its rows are Hadamard products of factor rows, and its transpose is the face-splitting product. ABMRS's rows are face-splitting products themselves. The brief's "Room for extensions" paragraph now says so.
- **`ThinSvd` fixed.** The Grassmann exp map gave it two parallel columns of norms 1.8·10⁻¹⁶¹ and 0.115. ζ² overflowed, the rotation became the identity, and the Jacobi never converged. Three changes fix it:
  - rotation angles are now computed through the ratio of the norms, and norms and cosines are rescaled outside the normal range;
  - columns below ε·‖R‖_F are left alone, and their left vectors are completed;
  - non-finite input is refused.

  Regression tests cover that case and extreme column scales. `PairwiseSum`'s accumulators are now vector lanes, bit-identical to the scalar form.
- **Callers moved onto `ThinSvd`:**
  - *`Pca`* now reads the SVD of the centred data, at 8.0 s against 22.3 s on a random 6,238 × 617 matrix, with eigenvalues agreeing to 2·10⁻¹³. Wide inputs are slower: 0.53 s against 0.32 s at 300 × 2,000. Components are now signed so that their largest entry is positive.
  - *`StiefelManifold` and `GrassmannManifold`* use it in place of their private Gram-route SVDs. They keep their rule of zero left vectors where σ < 10⁻¹²: at the cut locus that direction's sign is arbitrary, and the corrupted-block tests of MoMPCA and SPRED depend on the two choices cancelling.
- **Suites:**
  - LinAlg 118/118, Geometry 49/49, Mapper 28/28, Oracle 5/5, VizCore 240 passed with its 2 by-design skips, DimReduction 83/83.
  - Before the fix, the manifold change broke 8 tests: 2 MoMPCA, 4 DistributedSpred and 2 SpredCylinder.
  - SPRED's ISOLET pilot reruns its warm start to 10⁻¹⁴. Its annealed objectives move at 10⁻⁵ relative, because the annealer's trajectory diverges from bit-level differences. The pilot JSONs that the suite rewrites in `issues/spred/pilot/` were restored to their committed versions.
- **Noticed:** Maths.Oracle failed 2 of 5 once at baseline and then passed on a rerun with no change in between.

### Decision-free stage 2 numerics (N1–N5, N10) — 2026-10-26

Three new LinAlg types, in `src/maths/linalg/`. No existing code changed, and no caller is migrated yet.

- **`PairwiseSum`** (N10): `Of` and `Dot` sum pairwise. The tree runs over blocks of 128 terms, and each block is summed in four interleaved accumulators. It uses plain double arithmetic in a fixed order, so the bits are the same on every machine. `MatrixOps.Dot` and `ColumnGramMatrix` keep their SIMD left-to-right sums for their existing callers.
- **`BorderedCholesky`** (N1–N4), an `ISpdFactor` that keeps A and diag(A⁻¹) beside L.
  - *Births* border the factor with a new last row, in `DenseCholesky`'s arithmetic. A factor grown one row at a time is bit-identical to the strict one.
  - *Deaths* drop a row of L and restore the triangle with Givens rotations.
  - *Changes* are a death and then a birth, staged together on a copy. **The new column goes last**, so the BASS design layer reorders its terms to match.
  - *Proposals* (`ProposeBirth`, `ProposeDeath`, `ProposeChange`) return a `FactorProposal` without touching the accepted factor. It carries positive definiteness, the relative pivot, the log-determinant change and the smallest leave-one-out pivot. `ProposedInverseQuadraticForm` gives bᵀA′⁻¹b for the candidate, and `Accept` commits it. Deaths and changes are staged on a second copy and swapped in, which costs an O(k²) copy per proposal.
  - *diag(A⁻¹)* gains w·wᵀ/d on a birth. On a death it loses v·vᵀ/v_m, and any entry that keeps less than 10⁻³ of itself is recomputed from the new factor, because that subtraction cancels when a near-duplicate dies. `LeaveOneOutPivot(j)` = 1/(A_jj·(A⁻¹)_jj) is D06's s_j.
  - *Rebuilds* (`Rebuild`, and automatically every `rebuildInterval` accepted updates) refactor the kept A strictly and recompute diag(A⁻¹). They report the drift removed as a `FactorDrift`, and keep the incremental factor if the strict one fails. The default interval of 100 is a placeholder, not a measured choice.
- **`ThinSvd`** (N5): Householder QR, then one-sided Jacobi rotations on R, without forming AᵀA. U is completed where σ is exactly zero. The rank tolerance is max(m, n)·ε·σ₁, the same as the O-O4 oracle's.

**Tests** (`tests/linalg/`, 33 new): `Maths.LinAlg.Tests` 98 of 98 under the fact harness (`20261026_222717257`).
- *`PairwiseSumTests`*: errors are measured against exact BigInteger sums and stay within (32 + ⌈log₂(N/128)⌉ + 3)·u·Σ|terms|, including on F-O7's near-constant column.
- *`BorderedCholeskyTests`*:
  - births are bit-identical to `DenseCholesky`, condition estimate included;
  - over 400 random births, deaths and changes with no rebuild, the factor tracks a fresh factorization to 10⁻¹²·√A_ii, diag(A⁻¹) to 10⁻¹¹ relative, and the log-determinant, the proposals' predictions and bᵀA⁻¹b likewise;
  - the final rebuild finds drift ≤ 10⁻¹²;
  - rejected proposals leave every bit of the accepted factor unchanged;
  - leave-one-out pivots match the squared sine against the others, computed through `ThinSvd`, to 10⁻⁸;
  - D06's order counterexample: the newest pivot moves between 10⁻⁶ and 10⁻¹² with insertion order, while the leave-one-out vector does not;
  - after a near-duplicate's death, diag(A⁻¹) holds to 10⁻¹³. With the cancellation guard disabled, this test fails.
- *`ThinSvdTests`*:
  - reconstruction and orthonormality, for tall, wide and square matrices;
  - a constructed spectrum from 1 to 10⁻¹¹, recovered to 10·m·ε absolute;
  - rank 7 from a known product;
  - centred curves lose one rank;
  - zero matrices, and a zero column completed in the middle.

**Not done here:**
- N6 (low-rank solves), N7 (face-splitting contraction, which waits on the intercept call), N8 and N9;
- moving `Pca` and the geometry manifolds' private Gram-route SVDs (`StiefelManifold`, `GrassmannManifold`) onto `ThinSvd`;
- benchmarks.

### Cross-product summation — 2026-10-26

New item N10, from D06's fixture F-O7. On a design with a near-constant column, left-to-right summation of the Gram gave 7·10⁻⁴ to 2·10⁻³ relative error in s ≈ 10⁻⁹. Pairwise summation gave 2·10⁻⁷ to 3·10⁻⁷, which is ε·κ. Nothing in `src/maths/linalg` sums with compensation yet.

### Adaptation policy drafted — 2026-10-25

D10's note written as a proposal (`53a4876`, carried into this segment at the split).

