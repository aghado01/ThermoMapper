# Dev sequencing — conditioned-persistence engine + graph-engine renovation

**Status:** live sequencing record (updated 2026-07-03). **Tactical / transient** — lives in local `issues/`
while this arc is active, and **archives to the MarkBrain vault when the P0–P2 arc closes**. The design briefs it
links are the stable canonical reference, already in MarkBrain (relative `../../../MarkBrain/…`, both repos under
`D:\aghado01\`). Ties the two open fronts — the **PH engine** (`tda-purification/persistent-homology/`) and the
**graph-engine renovation** (`graph-engine-expansion/`) — into one ordered plan. Each brief keeps its own
detail; this page owns *order*, not design.

## The two-front map

**Locked / settled** (design done, execute against it — don't re-litigate):
- **Conditioned Persistence P0–P4** — [opus-brief-conditioned-persistence-synthesis.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/opus-brief-conditioned-persistence-synthesis.md).
  One prior, one directed skeleton, read twice (Φ homology ⊥ Ψ sheaf). Reuse-maximal: P0–P2 ride built primitives.
- **Placement** — [tda-placement.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/tda-placement.md), authoritative.

**Open / not locked** (design still live):
- **Graph sculptor + compiler modularity** — [opus-brief-graph-sculptor-modularity.md](../../../MarkBrain/ThermoMapper/issues/graph-engine-expansion/opus-brief-graph-sculptor-modularity.md).
- **Representation axis** (expansion ⟂ reduction) — [opus-brief-representation-axis.md](../../../MarkBrain/ThermoMapper/issues/graph-engine-expansion/opus-brief-representation-axis.md).

## DONE — the forced first move (commit `31ca346`, 2026-07-03)

Run as **one pass** because the fold and the injection fix touch the same seam (`H1CycleEdges` +
`GraphCompiler`), so the files get touched once, not twice:

- **`primitives → ph` fold** executed — `TDA.Primitives` retired, all 10 PH files now `TDA.Ph`/`TDA.Ph.Nerves`.
- **`graphs → tda` inversion cleared** — `GraphCompiler.Build` takes a `ProtectedEdgeSource` delegate; the H1
  protect-set is computed above the seam and injected by the UserRepl callers. `PreserveH1Cycles` and
  `Graphs.Proximity`'s TDA reference are gone.
- **`AGENTS.md` invariant written** — the construction-layering wall, with the cleared `GraphCompiler` violation
  as the worked example.
- **Gate:** solution builds; TDA.Ph.Tests 293/293, TDA.Mapper.Tests 28/28, RepoAudit 10/10. VizCore's 15
  failures pre-date the pass.

This unblocks **P0** — born in `tda/ph` against a clean tree, which was the whole reason to sequence it first.

## Main line — conditioned persistence (carries the momentum)

Engine-first, P0 → P1 → P2. This is the track with the first genuinely new *science* (P2's M2).

- **P0 · `ConditionedFiltration`** — ✅ **landed & green** (`src/tda/ph/ConditionedFiltration.cs`;
  TDA.Ph.Tests **298/298**, +5 P0 tests). Union convention (`BuildGraph`) + convenience `ComputeBarcode`;
  SIFTS `τ≡0` degenerate proven, barcode through the existing reducer unchanged. *Commit pending* (selective,
  atop the Spred session's tree). Brief: [p0-conditioned-filtration-brief.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/p0-conditioned-filtration-brief.md).
- **P1 · prior generalizes + multiparameter frame** — **← NEXT (scoped).**
  - **P1a** (residual band, monotone δ) — ✅ **landed** (`c91143d`; `src/tda/ph/ResidualPrior.cs`, TDA.Ph.Tests
    303/303). `ResidualEdges` feeds P0's `BuildGraph`; `τ≡0` subsumes P0's similarity, a nonzero prior shifts a
    return's birth. Brief: [p1a-residual-prior-brief.md](p1a-residual-prior-brief.md).
  - **P1b** (Δ reach axis) — **← NEXT (scoped).** brief: [p1b-delta-reach-brief.md](p1b-delta-reach-brief.md).
    A reach bound `|τ| ≤ Δ` making `(δ,Δ)` a monotone bifiltration read by PH slices. Tight, reuse-maximal.
  - **P1c** (deferred) — the non-monotone **zigzag** reader (persistent-Mapper-over-`T`); grounds in the built
    zigzag stack, verify against git before building — multi-track area
    ([zigzag-frontier.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/zigzag-engine/zigzag-frontier.md)).
  - **full multiparameter module** → **Z6**.
- **P2 · directedness v1.5** — directed flag complex (combinatorial) + the built magnetic Laplacian (spectral);
  experiments M1–M3. **M2 = two returns equal in persistence, separated by flux — the gate that earns the phase.**

**Load-bearing convergence:** the already-built `MagneticLaplacianOperator` serves P2's flux grader, the `λ_q`
order parameter, **and** the sculptor's criterion-C-as-magnetic-resistance (§4.5b). One object in
`graphs/spectral`, three threads — main-line work de-risks the sculptor's research bet for free.

## Second front — gated behind the main line

Compiler renovation (payload-aware two-tier `AdjacencyView` + vocabulary/type neutralization) → **sculptor
Phase 0** (static-rebuild sculptor, prove it improves SPC before any PMA).

- The **minimal decouple already landed** (the injection seam). The full renovation earns its slot when you
  actually build sculptor Phase 0 — not before. Open-design + big refactor is the classic project-swallower;
  both briefs are self-aware about this (Phase 0 before PMA; seam now, layout later).

## Cheap floaters (palate-cleansers between P-phases)

- **Toy frustrated-graph validation of C ≡ A-on-phases** (sculptor §4.5b) — highest information-per-hour in the
  portfolio; decides whether criterion C exists as a *separate* scorer at all. Best run **right after P2's
  M1–M3**, while the magnetic machinery is warm.
- **RFF as pure `maths` code** (RNG + matmul + cos, no PH/Stiefel) — buildable anytime; *wiring* it as the
  Stage-0 `Representation` bookend waits on the decoupled compiler, so no urgency.

## Parked (each behind an explicit gate in its own doc)

P3/P4 sheaf cohomology + persistent-sheaf-Laplacian stack (net-new; the barcode never depends on it) · the
PMA / BP-CSR dynamic layout (earned by measured iteration cost, not assumed) · criterion C's full design ·
the full multiparameter module (Z6 quasi-zigzag) · GLMY directed path homology (v2) · the representation-axis
§6 feed-direct-vs-SPRED-first branch.

## Open decisions — where they landed (2026-07-03)

1. **Interleave appetite** → **serial through P0** at minimum, *then* the renovation may start. P0 is short,
   and the fold-pass just destabilized the exact files the renovation would touch.
2. **Renovation scope, when it comes** → **bundle the vocabulary cleanup** (`Edge.J`, `EdgeWeightKind.Coupling`)
   with the `AdjacencyView` seam. The whole failure mode is assistants re-fusing *via the vocabulary*, so the
   naming neutralization is the anti-refusion medicine, not optional polish.
3. **Toy-C experiment timing** → **after P2's M1–M3**, on the warm magnetic machinery.

> One forced first move, done. Then: ship Φ up the main line (P0 → P1 → P2), earn Ψ; let the compiler
> renovation and sculptor wait behind P0; keep the toy-C check and RFF as cheap detours. The built magnetic
> operator is the pivot both fronts turn on.
