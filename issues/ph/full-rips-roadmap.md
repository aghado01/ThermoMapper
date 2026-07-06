# Full Vietoris–Rips → `src` — roadmap

**Status:** P1 landed 2026-07-06. `TDA.Ph.FullRips.Build` is a first-class,
threshold-bounded full-Rips density API over the existing 2-skeleton materializer. **Additive** — the
existing graph-restricted path is untouched; production (SPRED et al.) stays on it. Seed: the SPRED
oracle parity test (`tests/oracle/TdaParityTests.cs`) needed a full Rips and initially built one ad-hoc.

## Motivation

A first-class, threshold-bounded **full Vietoris–Rips** complex — the flag complex of the *complete*
proximity graph (all pairs ≤ threshold), as opposed to the flag complex of a *provided* skeleton. Uses:
oracle validation (Ripser parity), small-cloud exact topology, and general engine capability. It is a
sibling of the graph-restricted builder, not a replacement.

## Current state in `src`

The Rips machinery is capped at the 2-skeleton (H0/H1) — triangles are the top enumerated simplex.
There are now two density entry points over that same materializer:

- `FullRips.Build(points, maxDimension, threshold, label)` — complete Euclidean proximity graph
  (all pairs ≤ threshold) → `GraphRips`; for validation and small-cloud exact topology.
- `RipsFiltration.GraphRips(g, weights, maxDim)` — materializes a `SimplicialFiltration`; vertices +
  edges + triangles (via `FlagComplex.Triangles`) when `maxDim ≥ 2`. Graph-restricted and still the
  production SPRED path. **No tetrahedra.**
- `LazyRipsFiltration` — the lazy/cohomology-oriented `IFiltration`; also discovers **only triangles**
  (`CommonNeighbors`). Note its boundary/cofacet logic (`RemoveAt`, `PackKey`) is already
  *dimension-general* — only the enumeration is capped.
- `FlagComplex` — enumerates **only 2-simplices**; no k-clique enumeration.

Consumers, all graph-restricted / H0–H1: `PersistenceObjective` (SPRED), `ConditionedFiltration`,
`H1CycleEdges`. The graph-restriction is **load-bearing** for SPRED (recompile-a-sparse-kNN-graph per
proposal is the scalability story).

## What the parity test's "full Rips" is

`TdaParityTests` now uses `FullRips.Build(points, maxDimension: 2)` and hands the resulting filtration
to `PersistentHomology.Compute`. So it is "full" in **density** (all pairs vs a kNN skeleton) and still
the same **2-skeleton** (H0/H1), reusing all existing machinery.

## Two orthogonal axes

- **Density** — complete graph (all pairs, threshold-bounded) vs kNN-restricted skeleton. *P1 landed*
  for complete Euclidean density via `FullRips.Build`.
- **Dimension** — 2-skeleton (H0/H1) vs higher (H2 voids, …). *The genuinely missing capability;* capped
  in the **enumeration**, not the reducer. Deferred.

## Naming

- **`FullRips`** — the dense/complete builder. "Full Rips" is the standard TDA term; unambiguous.
- **`GraphRips`** (or `SkeletonRips`) — the current graph-restricted builder.
- **Avoid `SparseRips`** — it collides with Sheehy's *sparse-Rips approximation* (a specific linear-size
  construction with interleaving guarantees). Ours is *exact on the provided skeleton*, not an
  approximation. `DenseRips` isn't a standard term either.
- **`RipsFromGraph` → `GraphRips` landed 2026-07-06.** This keeps the graph-restricted API symmetric
  with `FullRips` without adding a pre-release compatibility alias.

## API / design scoping

Landed surface (P1), in `TDA.Ph`:

```
// FullRips.cs (new). No new project deps — CsrGraph/Edge already referenced by TDA.Ph.
public static SimplicialFiltration FullRips.Build(
    double[][] points,           // n×d cloud
    int maxDimension = 2,        // 2 = H0+H1 today; higher gated on the dimension axis
    double threshold = double.PositiveInfinity,  // include edges with dist ≤ threshold
    string label = "FullRips")
```

Internals: enumerate pairs `i<j`, keep those with `dist(i,j) ≤ threshold` as `Edge(i,j,dist)` →
`CsrGraph.FromEdges` → delegate to `GraphRips`. Returns the same `SimplicialFiltration` currency, so
it drops straight into `PersistentHomology.Compute` and every existing consumer/test.

- **Threshold semantics.** Bounds the whole complex (triangle births = max edge birth ≤ threshold).
  Default ∞ = true full Rips (all pairs). P2: offer Ripser's **enclosing-radius** default.
- **Distance metric.** Euclidean inline for P1 (keeps TDA.Ph dep-free). P2: optional `IDistanceMetric`
  (adds a `Graphs.Distance` edge) for non-Euclidean clouds.
- **Home.** Settled: `FullRips` is a dedicated static class; the graph-restricted builder remains a
  method on `RipsFiltration` named `GraphRips`.

## Performance

- **Size.** Unbounded full Rips: O(n²) edges, **O(n³) triangles**; the reduction then dominates. With
  threshold `r`: O(n² · avg-degree-within-`r`) — dramatically smaller. **Threshold is the primary lever.**
- **Reducer.** The standard `PersistentHomology.Compute` (SortedSet symmetric-difference) is ~O(m³)
  worst-case and slow on dense complexes. Route dense full-Rips reductions through the
  **cohomology/clearing path** already in `src` (`PersistentCohomology` + `PersistenceClearing`;
  `PersistentCohomologyRipsParityTests` already proves they agree with the standard path). Secondary lever.
- The n=30 dense-complex cost is **real but separate** from the parked R-env hang in the parity test.
- P2: benchmark full Rips at increasing n and record the tractable-n envelope (threshold on/off,
  standard vs cohomology reducer).

## Dimension axis (deferred)

To reach H_k, build the (k+1)-skeleton = all cliques up to size k+2:
- Generalize `FlagComplex` (triangles) → **k-clique enumeration** (bounded Bron–Kerbosch or degeneracy
  ordering), and extend the filtration builders to emit higher simplices (birth = max face birth).
- **`PersistentHomology.Compute` needs no change** (dimension-agnostic; reads boundary indices), and
  `LazyRipsFiltration`'s boundary/cofacet logic is already general — the whole gap is enumeration.
- Clique count blows up (O(n^{k+1}) dense), so higher dims are tractable only with a threshold and/or
  sparse graphs. **Build only when voids (H2) are actually needed** by the project's data.

## Phased plan

- **P0 — prereq ✔ 2026-07-05/06.** Rebuild the `r/` renv at the new path, unskip
  `TdaParityTests`, get it **green**. This anchors correctness (our PH + the full-Rips-via-complete-graph
  vs Ripser) before we build the API.
- **P1 — ball-rolling ✔ 2026-07-06.** Add `FullRips.Build` (threshold-bounded, reusing `GraphRips`); refactor
  `TdaParityTests.FullRipsBarcode` onto it; add a small self-contained unit test (known cloud → known
  diagram). Minimal, non-disruptive.
- **P2 — full integration.** Add the `IDistanceMetric` option and the enclosing-radius default; wire
  the cohomology/clearing reducer path + benchmark the envelope; scope (not necessarily build) the
  dimension axis.

## Non-disruption

SPRED's per-proposal recompile stays on `GraphRips` (sparse kNN skeleton); `ConditionedFiltration` and
`H1CycleEdges` stay graph-restricted. `FullRips` is purely additive — validation, exact small-cloud
topology, capability.

## Open decisions

- Threshold default: ∞ vs Ripser enclosing radius?
- Euclidean-only vs `IDistanceMetric` abstraction, and when?
- Does any project dataset actually need H2 (voids) — i.e., is the dimension axis ever on the critical path?

## References

- `src/tda/ph/{RipsFiltration,LazyRipsFiltration,FlagComplex,PersistentHomology,PersistentCohomology,PersistenceClearing}.cs`
- `tests/oracle/TdaParityTests.cs` (the `FullRips.Build` + Ripser oracle)
- SPRED tracker item: `issues/spred/dev-sequence.md` → "Full Vietoris–Rips → `src`" (pointer; this file is canonical)
