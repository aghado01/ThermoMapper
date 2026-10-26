# BASS port: response bases, prediction and sensitivity: checklist

Companion to [port-outputs-brief.md](port-outputs-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions (stage 1)

- [ ] **O-D07** Response units and uncertainty: drafted 2026-10-25 as a proposal; settles at the user's design pass.
- [ ] **O-D08** Sensitivity measures and moments: drafted 2026-10-25 as a proposal; settles at the user's design pass.

### Reference fixtures

- [ ] **O-O1** R defect fixtures: the truncation sampler's indexing, recycled noise, the residual's sign and units, 100% retention, constant columns, draw order and dropped axes, interval quantiles.
- [ ] **O-O2** Sensitivity defect fixtures: ignored category weights, unclipped truncation bounds, decomposition units, hard-coded effects.
- [ ] **O-O3** Moments against quadrature and category enumeration, including higher degrees and truncated mixtures.
- [ ] **O-O4** Truncation residuals: the in-sample residual is zero at full retention; leave-one-case-out equals explicit refits, and the computation in case coordinates equals the one in output space.

### Build (stages 5, 6, 9)

- [ ] **O1** PCA and supplied-basis decomposition on the thin SVD.
- [ ] **O2** Component execution and joint draws across components.
- [ ] **O3** Prediction: latent and predictive draws, component noise, truncation residuals, physical units.
- [ ] **O4** Pointwise density for supplied-basis fits.
- [ ] **O5** One moment service for direct and basis sensitivity.
- [ ] **O6** Direct sensitivity, including both functional forms.
- [ ] **O7** Basis sensitivity with cross-component terms and a subset budget.
- [ ] **O8** Effect functions for the declared measure and degree.
- [ ] **O9** Presentation data above maths (stage 9).

### Room for extensions

- [ ] **O-X1** The decomposed surface is a declared input.
- [ ] **O-X2** The moment service is written over piecewise polynomials.
- [ ] **O-X3** The basis contract admits robust bases, an uncertainty slot and a declared truncation model.
- [ ] **O-X4** Term derivatives are available to sensitivity.
- [ ] **O-X5** The response representation carries a declared reconstruction map (linear now; composition for elastic emulation).

## Reports

### Truncation residual held out of the basis — 2026-10-26

From reading Andros et al. (arXiv 2406.18751, 2602.09351). D07's draft now proposes leave-one-case-out truncation residuals by default, with R's in-sample residual as an option and their ratio as a basis-stability diagnostic. The scratch check behind it is exploration (one toy family), not a decision; O-O4 is the fixture to build. O-X5's reconstruction map now also declares where it can be evaluated.

### Prep for elastic emulation — 2026-10-26

"Room for extensions" gains the representation with a declared reconstruction map (O-X5), from Francom, Tucker et al.'s elastic calibration (arXiv 2305.08834).

### D07 and D08 drafted — 2026-10-25

Proposals, not decisions: `0844bdc`; categorical moments' source added in `a32c701`.

