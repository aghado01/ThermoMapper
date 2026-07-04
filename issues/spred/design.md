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
| 1. Engine | `SubspaceAnnealer` — simulated annealing over Gr(k, d) for any subspace objective | `Maths.Geometry.DimReduction` | LinAlg, Rng, Geometry |
| 2. PH objective | `PersistenceObjective` — `projection ↦ W_p(PH(project X), PH X)`, built via the injected GraphCompiler recipe | `TDA.DimReduction` | Graphs, TDA.Ph, LinAlg |
| 3. SPRED driver | `Spred` — wires objective + engine into "reduce while preserving persistence" | `TDA.DimReduction` | Maths.Geometry + row 2 |

**Why the split.** The engine, as written, has *no homology dependency* — the objective enters
through an injected `SubspaceObjectiveFunction` delegate. So it belongs with the manifold machinery,
not under a PH module: anyone wanting the Stiefel/Grassmann annealer for a non-topological objective
must not be forced to depend on `TDA.Ph`. The *objective* is what legitimately reaches into topology;
its barcode ingredients live in `TDA.Ph`, but assembling it also needs the `GraphCompiler` (to build
the projected graph), which `TDA.Ph` does not reference — so `PersistenceObjective` lives one layer
up, in the `TDA.DimReduction` consumer. The *driver* composes an engine from `Maths.Geometry` with
that objective, so it can sit in neither engine-level assembly (geometry must never depend on ph); it
lives beside the objective in `TDA.DimReduction`. Per engine-first, the objective + driver are
deferrable consumers — the engine is the primitive to make sound first.

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
`SubspaceObjectiveFunction` is the seam between the two: the engine optimizes any scalar objective,
and the topological identity lives entirely in the injected delegate (`PersistenceObjective.Evaluate`
is the concrete SPRED objective). `SubspaceAnnealer` + `SubspaceObjectiveFunction` are meant to read
as a coherent pair. Vocabulary convention: *objective* (never *cost*) throughout — the value the SA
minimizes, the delegate, the config, and the concrete type all say "objective".

## Source grounding & deviations from Yu & You (2106.02096)

Reviewed v3, 2026-07-03. Objective (Eq. 1): `min_{P∈St(n,k)} W_p(D_j(d_X), D_j(d_Y))`, `Y = XP` — the
Wasserstein distance between the persistent diagrams of the **ambient** and **projected** clouds, each
built from its *own* Rips (or Čech) filtration. §6 gives the multi-order form
`f_λ = λ·W_p(D_0) + (1−λ)·W_p(D_1)`. This confirms the recompile-per-proposal model, maps `Dimensions`
onto the paper's `(λ, 1−λ)`, and settles that `B_ref` uses the **same construction** as the projected
side (one recipe both sides; a different `ReferenceGraph` recipe is an opt-in experiment).

Two deliberate, credited deviations:

1. **Grassmann-geodesic SA, not the paper's Stiefel-QR perturbation.** The paper's basic PERTURB
   (§3.1) adds i.i.d. Gaussian noise to every entry of `P` and re-orthonormalizes by QR — a numerical
   fallback it adopts *because the Stiefel manifold has no closed-form exponential map*. §3.2 then
   switches to the **Grassmann** manifold for the median, justified by exactly the invariance the
   engine was rebuilt on: *"the persistent homology induced by a projection matrix is invariant under
   the right action of the orthogonal group, which justifies Grassmannian geometry."* Our
   `GrassmannManifold` has the closed-form Edelman/Bendokat exp map, so `SubspaceAnnealer` uses
   Grassmann geodesics for the SA too — unifying the walk and the median on the geometrically-correct
   manifold the paper could only afford for the median. Horizontal-only proposals also skip the wasted
   in-subspace O(k) rotations the paper's all-entry perturbation spends (unchanged objective by the same
   invariance). Trade-off: exact oracle parity with Rdimtools/the paper wants a Stiefel-QR PERTURB
   *mode* — deferred to validation (dev-sequence 3b).

2. **Graph-restricted, recipe-configurable Rips, not full Rips/Čech.** The paper builds the full
   Rips/Čech complex on the point cloud; we build a neighborhood-graph-restricted Rips via the injected
   `GraphCompilerConfig`. Full Rips is recoverable as a dense recipe (ε-ball at diam scale, or k=n−1);
   kNN recipes are the scalable generalization the paper flags as future work (§6). Čech is deferred
   (the `nerves/` path exists but is unwired here); ground metric is fixed L∞ by `DiagramMetrics` (the
   paper allows a general q-norm).

The §4 topological-equivalence measures (μ_quasi-iso, μ_equiv) are *post-hoc* quality measures via the
filtration homomorphism, **not** the SA objective — a separate validation deliverable (μ_quasi-iso is
barcode-computable per Prop 4.3; μ_equiv needs π₁ of a quotient, hard). See dev-sequence 3b.

## References
- Kisung You, arXiv:2106.02096 — the project's SPRED reference.
- Memory: `project_spred`, `project_kisungyou_dr_track`, `project_intrinsic_geometry_northstar`.
- Not to be confused with the hyperbolic / HoroPCA-style DR (arXiv 2604.24895) — SPRED proper is
  Euclidean linear DR with a PH objective, a separate track from the negative-curvature branch.
