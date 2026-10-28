# BASS port: calibration and SMC backends: checklist

Companion to [port-calibration-brief.md](port-calibration-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions (stage 1)

- [ ] **C-D09** Calibration target and covariance formulations: drafted 2026-10-25 as a proposal; settles at the user's design pass.

### Reference fixtures

- [x] **C-O1** R defect fixtures: the stale likelihood after a discrepancy draw, the branch chosen by shape, normalized returned coordinates, the unit mismatch, the all-groups `s2.df` switch: `r/oracles/bass_calibration_defects_oracle.R`.
- [x] **C-O2** Dense covariance against the low-rank path; scaled and unscaled y give the same physical posterior: `r/oracles/bass_covariance_oracle.R`.

### Build (stages 7–8)

- [ ] **C1** Typed covariance formulations: functional and multivariate.
- [ ] **C2** Variance groups, discrepancy, and categorical calibration parameters.
- [ ] **C3** Nested conditional reference sampler, with inner-run error measured.
- [ ] **C4** Fitting SMC adapter on D05's sequence.
- [ ] **C5** Cut-calibration SMC adapter, with ordering and induced-density bridges.

### Room for extensions

- [ ] **C-X1** The measurement-error model is a declared component.
- [ ] **C-X2** Emulator draws behind an interface.

## Reports

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

