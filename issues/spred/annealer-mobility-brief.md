# SubspaceAnnealer mobility — proposal redesign brief

**Status:** ready to execute
**Date:** 2026-07-16
**Home:** engine layer — `src/maths/geometry/dim-reduction/SubspaceAnnealer.cs` + `tests/geometry/dim-reduction/`
**Seed evidence:** [chip-grassman-median.md](chip-grassman-median.md) iteration-budget probe (`2c84cac`)

## The finding this executes against

At I=1000 on ISOLET S0 (Gr(617, 30), intrinsic dimension 30×587 = 17,610), seven of eight blocks end
bit-identical to the PCA warm start — one accepted move in 8,000 block-iterations, and that one is the
very first step at full temperature. Compute is a non-constraint (corrected marginal 0.092 s/iter,
~34 min full S0); the block is proposal design. Mechanism: an isotropic random horizontal tangent's
directional derivative shrinks like 1/√dim while the finite-step curvature penalty stays O(step²), so
the improving fraction of fixed-length isotropic proposals vanishes at high codimension, and the
`0.99^iter` cooling extinguishes uphill acceptance by ~iteration 500. More budget provably does not
help. The cylinder-scale proposal design does not transfer.

## Design (settled)

1. **Two-plane Givens proposals — the primary move.** Pick a retained column i ∈ [k] and a unit
   direction v orthogonal to span(Y) (draw a random d-vector, project off span(Y) — the existing
   `HorizontalTangent` projection specialized to one column — and normalize); rotate
   `y_i ← y_i·cosθ + v·sinθ`. This is exactly the Grassmann geodesic for the rank-1 horizontal
   tangent `θ·v·e_iᵀ` (its SVD has a single σ = θ), so it stays inside the design.md
   "geodesic SA on Grassmann" story — `ExpMap` may be reused as-is or specialized closed-form.
   Rationale: restores O(1)-dimensional moves whose improving fraction does not vanish with
   codimension; the classic subspace-search move.
2. **Acceptance-adaptive step scale.** Adapt the θ scale toward ~25% acceptance (multiplicative or
   windowed update) with a floor and ceiling, replacing the fixed `temp·0.1` magnitude. Cooling
   continues to govern the Metropolis temperature only. Adaptive-alone is known-insufficient: it
   buys acceptance by shrinking toward diffusion-scale steps; it is the controller inside the
   two-plane move, not the fix itself.
3. **Isotropic fraction.** A small isotropic-proposal mixture (default candidate 10%) is optional
   ergodicity insurance — decide during implementation; pure two-plane is acceptable if the cylinder
   validates.
4. **Config surface.** A declarative `SubspaceAnnealerOptions` record (proposal kind/mixture, target
   acceptance, step floor/ceiling, cooling), defaults reproducing sensible current behavior;
   `Compute` takes the options. One Xoshiro stream; adaptation is state-dependent but deterministic
   given the seed — same-seed bit-identical runs remain the contract.

## Validation

- **New engine-level mobility fact (the key test):** a synthetic quadratic-on-Grassmann objective at
  high codimension (e.g. d=200, k=5, minimize a distance/trace objective to a target subspace) —
  the current isotropic proposals stall from a deliberately offset start, the new proposals descend.
  Fast, deterministic, no PH dependency.
- Cylinder fixture still beats every axis-aligned view (`SpredCylinderTests`) — no small-d regression.
- Same-seed determinism facts updated for the options surface.
- The S0 acceptance/objective re-probe belongs to the pilot lane (`SpredIsoletPilotTests` — owned by
  the concurrent pilot session; coordinate, do not collide).

## Tracked, not executed

- **Gradient-informed Riemannian search** (PH subgradients via matched pairs) — the bigger lift;
  revisit only if two-plane + adaptive stall on S0.
- Restart ensembles / early-stop on plateau (dev-sequence P3 remainder) — orthogonal to mobility;
  cheap to add later on top of the options surface.

## Report-back

Append implementation and probe results here, per chip conventions.
