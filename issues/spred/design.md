# SPRED — Design

**Updated:** 2026-07-03

SPRED is a *topology-preserving linear* dimensionality reduction: it optimizes a linear projection
over the Grassmann manifold against a **persistent-homology objective** — "PCA keeps variance, UMAP
keeps neighborhoods, SPRED keeps the persistence." (Kisung You, arXiv:2106.02096; the project's
reference — see the `project_spred` memory note.)

This doc holds the architecture and the rationale behind each layer. For the ordered build plan and
current status, see [dev-sequence.md](dev-sequence.md). It **supersedes** the point-in-time
`project_spred` / `project_kisungyou_dr_track` memory notes where they disagree (those predate the
engine rename).

## The three-layer split

The track is deliberately split so the optimizer stays a reusable geometry primitive and the
topology lives with the barcodes:

| Layer | What | Home | Depends on |
|---|---|---|---|
| 1. Engine | `SubspaceAnnealer` — simulated annealing over Gr(k, d) for any subspace cost | `Maths.Geometry.DimReduction` | LinAlg, Rng, Geometry |
| 2. PH cost adapter | `projection ↦ W_p(PH(project(X)), PH(X))` | `TDA.Ph` | PersistenceBarcode, DiagramMetrics, RipsFiltration |
| 3. SPRED driver | wires engine + PH cost into "reduce while preserving persistence" | new consumer above geometry+ph | Maths.Geometry **and** TDA.Ph |

**Why the split.** The engine, as written, has *no homology dependency* — the objective enters
through an injected `SubspaceCostFunction` delegate. So it belongs with the manifold machinery, not
under a PH module: anyone wanting the Stiefel/Grassmann annealer for a non-topological cost must not
be forced to depend on `TDA.Ph`. The *cost* is what legitimately reaches into topology, and every
ingredient it needs (barcodes, diagram metrics, Rips filtration) already lives in `TDA.Ph` — so the
adapter lives there. The *driver* composes an engine from `Maths.Geometry` with a cost from
`TDA.Ph`, so it can sit in neither (geometry must never depend on ph); it belongs in a consumer
layer above both. Per engine-first, the driver + adapter are deferrable consumers — the engine is
the primitive to make sound first.

## Design notes

**Grassmann over Stiefel.** The objective is subspace-invariant — rotating within the target
subspace is an isometry of the projected cloud, so its persistence is unchanged — which makes the
Grassmann quotient the correct substrate. `GrassmannManifold.ExpMap` is a true Edelman geodesic;
`StiefelManifold.ExpMap` is only an add-tangent-then-orthonormalize retraction (its own comments
flag it "simplified"), so it is the wrong engine substrate. Reserve Stiefel/Grassmann *estimators*
for the distributed-aggregation target space (dev-sequence, distributed track), not the anneal step.

**Proposal + step semantics.** A proposal is a random isotropic ambient Gaussian projected onto the
horizontal tangent space at the current subspace (`Δ ← Δ − Y(YᵀΔ)`), scaled to a geodesic step, then
retracted via `GrassmannManifold.ExpMap`. The step is a genuine geodesic distance (principal-angle
radians): `step = temp · 0.1`, `temp = 0.99^iter`. This is a behavioral change from the old
per-entry perturbation magnitude; convergence values are not pinned by tests.

**Naming seam.** `SubspaceAnnealer` (engine) reserves the name *SPRED/Spred* for the driver.
`SubspaceCostFunction` is the seam between the two: the engine optimizes any scalar cost, and the
topological identity lives entirely in the injected delegate. `SubspaceAnnealer` +
`SubspaceCostFunction` are meant to read as a coherent pair.

## References
- Kisung You, arXiv:2106.02096 — the project's SPRED reference.
- Memory: `project_spred`, `project_kisungyou_dr_track`, `project_intrinsic_geometry_northstar`.
- Not to be confused with the hyperbolic / HoroPCA-style DR (arXiv 2604.24895) — SPRED proper is
  Euclidean linear DR with a PH objective, a separate track from the negative-curvature branch.
