# BASS port: shared numerics and sampler contracts: checklist

Companion to [port-numerics-brief.md](port-numerics-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions

- [ ] **N-D10** Adaptation policy for adaptive kernels (note): drafted 2026-10-25 as a proposal; settles at the user's design pass.

### Build (stages 2–3)

- [ ] **N1** Bordered insertion on the strict dense factor, with candidate work private until acceptance.
- [ ] **N2** Givens deletion, and a changed term as a staged delete and insert.
- [ ] **N3** Strict rebuild: periodic and on failure, with drift checks against fresh factorizations.
- [ ] **N4** diag(G⁻¹) maintained through births and deaths, for D06's order-free gate.
- [ ] **N5** Stable thin SVD with numerical-rank handling in LinAlg.
- [ ] **N6** Diagonal-plus-low-rank solves and log-determinants, without m × m inverses.
- [ ] **N7** Face-splitting design contraction with the Cartesian (Hadamard) shortcut.
- [ ] **N8** Generic SMC core in `samplers/smc` (optional; specified in the calibration segment).
- [ ] **N9** Qualified baseline rerun after the shared changes, before the first BASS integration tests.

### Room for extensions

- [ ] **N-X1** The face-splitting contraction covers tensor designs (ABMRS).
- [ ] **N-X2** `ISpdFactor` admits a wide-band factor beside `DenseCholesky`.
- [ ] **N-X3** `ConjugateUpdates` expressible as limits of a GIG update.

## Reports

### Adaptation policy drafted — 2026-10-25

D10's note written as a proposal (`53a4876`, carried into this segment at the split).

