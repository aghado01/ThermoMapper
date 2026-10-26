# BASS port: calibration and SMC backends: checklist

Companion to [port-calibration-brief.md](port-calibration-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions (stage 1)

- [ ] **C-D09** Calibration target and covariance formulations: drafted 2026-10-25 as a proposal; settles at the user's design pass.

### Reference fixtures

- [ ] **C-O1** R defect fixtures: the stale likelihood after a discrepancy draw, the branch chosen by shape, normalized returned coordinates, the unit mismatch, the all-groups `s2.df` switch.
- [ ] **C-O2** Dense covariance against the low-rank path; scaled and unscaled y give the same physical posterior.

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

### D09 drafted — 2026-10-25

Proposal, not a decision: `0844bdc`; categorical calibration parameters added in `a32c701`.

