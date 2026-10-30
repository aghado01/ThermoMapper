# BASS port: calibration and SMC backends: checklist

Companion to [port-calibration-brief.md](port-calibration-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions (stage 1)

- [x] **C-D09** Calibration target and covariance formulations: drafted 2026-10-25 as a proposal; decided by the user on 2026-10-29 (the brief's D09, amendments marked).

### Reference fixtures

- [x] **C-O1** R defect fixtures: the stale likelihood after a discrepancy draw, the branch chosen by shape, normalized returned coordinates, the unit mismatch, the all-groups `s2.df` switch: `r/oracles/bass_calibration_defects_oracle.R`.
- [x] **C-O2** Dense covariance against the low-rank path; scaled and unscaled y give the same physical posterior: `r/oracles/bass_covariance_oracle.R`.

### Build (stages 7–8)

- [x] **C1** Typed covariance formulations: functional and multivariate (stage 7 chip, `24014cd`).
- [x] **C2** Variance groups, discrepancy, and categorical calibration parameters (`50565b9`; the conditional sampler `d3af6a3`).
- [x] **C3** Nested conditional reference sampler, with inner-run error measured (`ab14608`).
- [x] **C4** Fitting SMC adapter on D05's sequence (stage 8 chip C4: `d0a0fd8`, `9e8c349`, `b0978d5`, `73352f0`).
- [x] **C5** Cut-calibration SMC adapter, with ordering and induced-density bridges (stage 8 chip C5: `86c5969`, `2e77db3`, `eb05485`, `cb0ecef`).

### Room for extensions

- [x] **C-X1** The measurement-error model is a declared component. Built in stage 7 itself (D09, 2026-10-29; `3e935b5`).
- [x] **C-X2** Emulator draws behind an interface. Built in stage 7 itself (D09, 2026-10-29; `066a054`).

## Reports

### Stage 8 verified: fitting SMC (C4) and cut SMC (C5) — 2026-10-29

The two stage 8 chips built them in parallel ([C4's report](stage8-fitting-smc-brief.md#reports), [C5's](stage8-cut-smc-brief.md#reports)), and they were read and verified here.
- *Commits:* C4 `d0a0fd8`, `9e8c349`, `b0978d5`, `73352f0` and report `9036a38`; C5 `86c5969`, `2e77db3`, `eb05485`, `cb0ecef` and report `e3ddc02`. `src/maths/samplers/smc/` is untouched.
- *Read here:*
  - `DrawSample` is Akinshin's window over the pooled step function, with n* counted over the individual draws;
  - the incremental weights are the core's `Temperature.LogDensity` differences;
  - each cut-SMC visit pools at exactly 1/S.
- *Suites rerun here:* Samplers 206/206, Regression 1,169/1,169 and the BASS oracle classes 71/71 at the chips' last commits; after the presentation fix here (`ebf0697`), its classes 40/40 and the full Regression suite 1,170/1,170.
- *A fix here,* the gap C4 named in stage 9's files: an SMC fit's `Traces()` now refuses (its draws are a population with no order), and its `BassRunLengths` carry its `SmcSpec` with the sweep counts null.
- *What the readings show:*
  - fitting SMC finds the best MCMC representation on the mixed and functional fits, but not on Friedman;
  - lineages collapse to one or two families per replicate, so the genealogy error saturates and the replicates' spread is the error;
  - cut SMC's tour halves the bridge stages against a random order, while C3 stays cheaper at S = 6;
  - every C3 inner run converged under the default gate, which closes stage 7's open item.
- *Docs:*
  - `docs/Maths/BayesianRegression.md` gains "BASS: SMC" and cut SMC in the calibration section, the new weighted quantile, the weighted pointwise forms, and the API, configuration, status and source rows;
  - `docs/Maths/Samplers.md` gains `DrawSample`'s rule.
- *Open:*
  - Friedman's worse modes under SMC: the lineages, not the particle count, bind;
  - cut SMC at larger S;
  - C3's short-run bias, which its gate does not see;
  - N15, the covariance factor's cost;
  - N16, the core's conveniences.

### Stage 7 verified: calibration — 2026-10-29

The stage 7 chip built it ([its report](stage7-calibration-brief.md#reports)), and it was read and verified here.
- *Commits:* `066a054` (C-X2), `3e935b5` (C-X1), `24014cd` (C1), `50565b9` (C2), `d3af6a3` (the conditional sampler), `ab14608` (C3) and `c92255e` (the report). Outside `calibration/`, it added `Maths.Distributions.HalfCauchy` and extended the two R oracles. No accessor was added to BASS's fits.
- *Read here:*
  - the variance move's Jacobian stays untempered (the log target is the current density plus u);
  - weighted outer draws are resampled by weight, then each run is weighted equally;
  - the half-Cauchy's density, CDF and quantile are right;
  - the cut holds structurally: each inner run fixes its emulator draw.
- *Suites rerun here:* Regression 1,101/1,101, Distributions 148/148, the BASS and calibration oracle classes 63/63, and R's stage-1 checks, 97 verdicts passing (run from r/).
- *As built, where D09's text was loose* (recorded in D09 as built):
  - the multivariate correlation is over the outputs, m × m, as R's `meas.error.cor` is;
  - the multivariate truncation term defaults to the residual set's second moment about zero, like the functional one;
  - the diagonal option is R's centred n − 1 variance;
  - above the cap, single-site steps are Metropolis;
  - `Pooled` also leaves out undetermined runs.
- *Open, for later:*
  - how often inner runs converge under the default gate (R̂ 1.01, ESS 100 per chain): the recovery runs declared 1.05 and 50, and 3–6 of 6 converged per seed;
  - `DiagonalPlusLowRank`'s per-proposal thin SVD (N15);
  - a per-site Gibbs draw above the cap.
- *Docs:* `docs/Maths/BayesianRegression.md` gains "BASS: calibration" and its API, configuration, status and source rows; `docs/Maths/ProbabilityAndInformation.md` gains `HalfCauchy`.
- *Next:* stage 8, C4 and C5.

### D09 decided — 2026-10-29

The user took D09 with argon's resolutions and four amendments from the main session, each marked *(2026-10-29)* in [the brief](port-calibration-brief.md#d09-decided-2026-10-29-calibration-target-and-covariance-formulations).
- *Adopted as drafted:*
  - the cut posterior;
  - the declared covariance formulations, functional on `DiagonalPlusLowRank` and multivariate on the dense factor;
  - per-group variance priors with R's inverse gamma as the default;
  - categorical calibration parameters;
  - the nested conditional reference sampler.
- *Argon's resolutions:*
  - C-X2's emulator interface and C-X1's declared measurement model are built in stage 7 itself;
  - shapes and positive definiteness are validated at construction;
  - categorical enumeration is capped at 1,024 combinations, with single-site updates above the cap.
- *The four amendments:*
  - an emulator without a basis adds σ²_e·I on the diagonal;
  - the functional formulation carries D07's truncation residuals as a low-rank term, with R's diagonal as an option;
  - the heavy-tailed variance prior is a half-Cauchy on σ_g with scale ŝ_g, and R's improper, unit-dependent kernel is a defect fixture;
  - the discrepancy's scale s_v is a required declaration, s_v = 1 reproducing R.
- *Stage 7 is open.* C1–C3 build against this, with C-X1 and C-X2.

### R's calibration tempering (C-O1) — 2026-10-27

From the reference-fixtures chip's follow-up, `f7acf2a`: a `tempering` section in `r/oracles/bass_calibration_defects_oracle.R`. One variance group, IG(5, 0.05), with θ at θ* and the discrepancy at 0, on the ladder 1, 2, 5. All checks pass.
- **R's moves are exact for a powered target.** The variance move's `alpha`, evaluated as written, is exactly b times (log likelihood + log prior kernel + log Jacobian). On a 61-point grid of log σ², R's kernel reaches exp(b(ℓ + log prior + log σ²)), the target powered in log σ² coordinates, to 5·10⁻¹⁶ in total variation. The discrepancy draw's covariance, `discrep.S / itl[t]`, is exactly that target's conditional, S/b.
- **R's swap leaves that target.** It omits the Jacobian and the discrepancy prior. Its largest gap from detailed balance against its own rung laws is (b₁ − b₂) times the grid's log σ² range, 4.61 and 2.76. Between states that differ only in the discrepancy, the omitted prior costs up to 1.06.
- **Against D09's likelihood path** R is far off: total variation 0.12 and 0.27 at b = 0.5 and 0.2, swap gaps of 234 and 140, and discrepancy covariances 1.4% and 5.2% off (I/τ + b DᵀΣ⁻¹D)⁻¹. D09's kernels reach the likelihood path to 1·10⁻¹⁵, and its swap, (b₁ − b₂)(ℓ(x₂) − ℓ(x₁)), balances on them to 1·10⁻¹³.
- **For D09:** tempering the Jacobian is not by itself the defect. R's variance move and discrepancy draw are consistent with a target powered in log σ², which is what D05 calls the powered path in log coordinates; the swap alone breaks it. D09's correction to the likelihood path stands on D05's reasons, not on that inconsistency.
- **Not covered:**
  - The emulator draw shared across temperatures.
  - Several variance groups on the grid.
  - The stale likelihood's interaction with the swap.

### Calibration defects and covariance references (C-O1, C-O2) — 2026-10-27

From the reference-fixtures chip ([brief](reference-fixtures-brief.md)), `a6960f0`. Two oracles, run by `r/scripts/bass_stage1_checks.R`, all checks passing.
- **C-O1,** `r/oracles/bass_calibration_defects_oracle.R`. R's `calibrate.bassBasis` (`2671014`) runs one iteration as written on a mock emulator served through R's own `predict.bassBasis`. The setup: m = 20, two components, y_s from 0.5 to 2.5, θ on [10, 20] × [−1, 1], eight rungs.
  - *Stale likelihood.* On the one rung whose θ move was rejected, the returned state's cached log-likelihood differs by 0.674 from that of its own prediction, discrepancy and covariance. That difference entered the θ acceptance test. Every other rung matches to rounding.
  - *Branch by shape.* Type 1 with a square discrepancy matrix runs the multivariate branch, and type 2 with a two-column one runs the functional branch.
  - *Coordinates.* R returns θ₁ in [0.15, 0.96] while the emulator was called on [11.52, 19.56]. R's starting θ is 18.6 away from the inputs of the emulator's first call, which predicts at it. Mapped to physical coordinates it matches them exactly.
  - *Units.* R's log-determinant equals the dense one in R's units to rounding, and is 0.45 below the physical one on every rung.
  - *The `s2.df` switch.* With `s2.df` = (0, 10), the second group's log prior kernel is −log(1 + σ²): −0.001 to −0.693 for σ² from 0.001 to 1. Its IG(5, 0.05) kernel gives −8.55, 22.63, 13.32 and −0.05. With (10, 10), R uses the IG.
- **C-O2,** `r/oracles/bass_covariance_oracle.R`.
  - *Functional formulation.* The reference is the dense Cholesky solve and log-determinant, checked by a symmetric eigendecomposition. R's `swm` and `swm.ldet`, and D09's low-rank route (the factor of I + E^{1/2}UᵀD⁻¹UE^{1/2}), match it to 1·10⁻¹⁴ at κ = 1·10³ and 4·10⁻¹⁰ at κ = 9·10⁷.
  - *At κ = 9·10¹³* the two dense algorithms disagree by 1·10⁻³, so no dense reference certifies beyond that. Every route, dense or low-rank, is backward stable there (structured residual ≤ 1.4·10⁻¹⁵), and their forward gaps of about 5·10⁻⁴ are the conditioning.
  - *Multivariate formulation:* an AR(1) group correlation, a squared-exponential discrepancy and the full truncation covariance, κ = 4·10³. Cholesky and eigendecomposition agree to 1·10⁻¹³.
  - *Units.* The likelihood computed in standardized coordinates is shifted by exactly Σ log y_s (spread 2·10⁻¹⁴), so the θ posteriors agree to 9·10⁻¹⁷. With R's units, for contrast, the posterior moves (total variation 0.96, mean 19.9 against 16.0).
- **Finding for D09:** at these sizes R's Woodbury helpers are as accurate as the low-rank route. D09's case against m × m inverses rests on cost and storage, not accuracy.
- **Not covered:**
  - The tempering defects (variance Jacobians tempered, swap terms omitted).
  - The discrepancy's fixed prior precision.
  - More than one iteration of R's chain.
  - The cut target and the nested sampler (C3).

### D09 drafted — 2026-10-25

Proposal, not a decision: `0844bdc`; categorical calibration parameters added in `a32c701`.

