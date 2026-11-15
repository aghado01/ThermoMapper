# GBASS noise, slices 0 and 1: landed (chip, 2026-11-05)

Commits on `main`, nothing pushed: `ee2cce6` (slice 0, the oracle modes and parity tests), `e978092` (slice 0's report), `c6e0baf` (slice 1, Student t on the core, code and tests), `c45dd3b` (slice 1's docs, ledger and report). Reports: `issues/regression/bass/gbass-noise-brief.md`, under "Reports".

For whoever builds slice 2 next:
- `gbass_conditionals_oracle.R` and `gbass_primitives_oracle.R` already hold the skewed conditionals (w's modified half-normal, β's normal) and the GIG/MHN laws; the engine's GIG sampler and MHN meet them there.
- The carrier: `ObservationLatents` (v, β) rides `WorkingResponse.Latents`; σ² and a ride `BassStructureMoves.SweepFit` within a sweep (not the record; reason in the report).
- `StudentTMixing.SampleConditional` and `VarianceMeanMixture.At` refuse a drift with a pointer to slices 2 and 3; the latent block reads `working.Values`, which under skew is y − σβv, so it must add the offset back for the total residual.
- Retention: 2N doubles a draw under t; the "v alone on the draw" lever is proposed for slice 2, where it quarters.
- Geweke harness: `BassStudentTQualificationTests`; Zellner–Siow τ makes the successive chain sticky, so the cases use τ ~ Gamma(5, rate 5N) or held.
- Running gates: `cmd start /belownormal` is denied in this sandbox; a PowerShell that lowers its own priority (dotnet test and testhost inherit) works.

`artifacts\chip-gbass-noise\` was deleted.
