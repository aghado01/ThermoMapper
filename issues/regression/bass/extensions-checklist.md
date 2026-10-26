# BASS extensions: checklist

Companion to [extensions-brief.md](extensions-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Entries

- [ ] **E1** Generalized-hyperbolic noise (GBMARS), with OQ20's BARS side.
- [ ] **E2** Robust centring and bases (OQ16).
- [ ] **E3** Free-knot tensor-product surfaces (ABMRS).
- [ ] **E4** Additional proposal kernels.
- [ ] **E5** KHAOS's modified g-prior.
- [ ] **E6** Basis-estimation uncertainty.
- [ ] **E7** Moments for piecewise polynomials.
- [ ] **E8** Active-subspace analysis.
- [ ] **E9** Cross-engine evidence.
- [ ] **E10** Direct-functional BASS for families of SPC sweeps.
- [ ] **E11** Elastic emulation of misaligned functional outputs.
- [ ] **E12** Residuals correlated along the output index.

## Reports

### GBMARS's global scale derived — 2026-10-26

From the user's derivation. With u = √w and zᵢ = u·vᵢ, u has a GIG full conditional for any GIG mixing. The prior must go on u (Gamma, inverse-gamma or 1/u) and the coefficient prior's variance must scale as u, or not at all. Variance ∝ w leaves the family, and variance ∝ u puts y's units into Σ₀ or τ (OQ19). The scale-equivariant alternative is log-concave in log u and is left to a slice step. A numerical check is in E1's brief entry. E1 stays deferred.

### Andros et al. read — 2026-10-26

From arXiv 2406.18751 and 2602.09351: new entry E12 (correlated residuals along the output index); E9 gains held-out proper scores (the energy and interval scores); E10 notes the functional GP as a comparison point. Sketching and the functional GP are not adopted, with reasons in the hub's register.

### Elastic emulation added — 2026-10-26

New entry from Francom, Tucker, Huerta, Shuler and Ries (arXiv 2305.08834), with its prep in the outputs segment (O-X5).
