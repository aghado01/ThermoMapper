# Brief — P1c: the non-monotone reader (H0 zigzag over a conditioned sweep)

**Status:** build ticket — **local `issues/ph` (transient; archive to MarkBrain when P1 closes).** Third increment
of synthesis §6 **P1** — the non-monotone slice. **Scope in one line:** read a *non-monotone* sweep of the
conditioned filtration (edges enter **and** leave — a trajectory / sliding window / `T`-sweep) through the built
**graph-zigzag** oracle, giving **H0 components-over-`T`** (birth / death / merge / split) — the
persistent-Mapper-over-`T` telos, seated on the zigzag engine proper. **The filled-H1 case is the harder half,
deferred; the full module stays Z6.**

**Background (MarkBrain — canonical):** synthesis §6 P1 (non-monotone slice → zigzag; "the honest carrier for
persistent-Mapper-over-`T`") — [opus-brief-conditioned-persistence-synthesis.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/opus-brief-conditioned-persistence-synthesis.md).
Zigzag-engine status — [zigzag-frontier.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/zigzag-engine/zigzag-frontier.md).
Builds on P1a/P1b ([p1a](p1a-residual-prior-brief.md), [p1b](p1b-delta-reach-brief.md)). **Sequencing:** [dev-sequencing.md](dev-sequencing.md).

> **⚠ Zigzag-engine discipline (from the frontier).** That engine is multi-track and its per-rung briefs *lag
> committed code*. **Before building: `git log`, glob the target types, and trust the code over any prose.** This
> brief's reuse claims were grounded at HEAD `66d6f24` (`ZigzagFiltration`, `GraphZigzag.Compute`,
> `ZigzagBarcodeNaive`/`FastZigzag` present; TDA.Ph.Tests 306/306) — re-verify `GraphZigzag.Compute`'s exact
> input signature at build time.

---

## 1. Why non-monotone needs a zigzag (the boundary with P1b)

P1a/P1b are **monotone** — growing `δ` or `Δ` only *adds* edges, so every slice is ordinary PH. A non-monotone
driver — a **trajectory** where the field `{t_i}` evolves, a **sliding window**, a **`T`-sweep** — makes edges
**enter and leave** between frames. A component can *split* as well as merge. That is exactly zigzag persistence,
and it is where the built zigzag engine becomes the reader (not a re-derivation).

## 2. The reuse surface (grounded at HEAD)

- **`GraphZigzag.Compute`** — the graph-zigzag oracle: **H0 always** (parity-complete, cross-validated vs the
  Z1/Z2 oracles), H1 folded when `maxDimension ≥ 1`. Dim-0 (vertices) + dim-1 (edges) — graphs are 1-complexes.
  **Its H1 is graph cycle space — NO triangle filling** (unlike P0/P1's Rips-H1). This is the crux of §3's split.
- **`ZigzagFiltration`** — the general cell-level currency (`Add(cellId, boundary)` / `Delete(cellId)`,
  starts/ends empty), read by `ZigzagBarcodeNaive` / `FastZigzag` for the full *filled* simplicial case.
- **`PersistentMapper.BuildFiltration` + `NerveDiff`** — the **started** persistent-Mapper-over-`T`, but it tracks
  nerve topology **frame-to-frame, not through a zigzag oracle**. P1c seats a conditioned sweep on the oracle
  proper — the open work the frontier names. **Complement `NerveDiff`; do not rebuild it.**

## 3. The split — the H0 / H1 faithfulness boundary

- **P1c = H0 components over a conditioned sweep**, via `GraphZigzag` H0. Unambiguous (no filling question),
  telos-core (clusters birth / merge / split over `T`), reuse-maximal. **This brief.**
- **Filled-H1 over the zigzag (deferred — the harder half).** `GraphZigzag`'s folded H1 is *unfilled* cycle
  space; P0/P1's H1 fills loops with triangles. Reconciling that over a zigzag needs the full simplicial
  `ZigzagFiltration` (vertices + edges + **triangles**, cell-level, stable global ids) read by
  `ZigzagBarcodeNaive` / `FastZigzag`. **Its key decision:** does a return-over-`T` die by triangle-fill (Rips)
  or by edge-departure (graph)? — deferred, not dodged. Formalize as P1d if/when taken up.
- **Full multiparameter module** → **Z6**.

## 4. The one new thing — the sweep → churn driver

Everything else is `GraphZigzag`. The new code is a **driver** turning a *sequence of conditioned graphs* into the
graph-zigzag insert/delete sequence:
- **input:** a sequence of conditioned edge-sets (each a P1a/P1b residual-edge set at a sweep step — an evolving
  field `{t_i}`, or a sliding window over a trajectory);
- **diff** consecutive frames → edges added / removed, with **stable vertex ids** (`0..n-1`) and per-edge cell
  ids; an edge that leaves and returns is a re-entry (`GraphZigzag` is multigraph-safe);
- **output:** the churn sequence `GraphZigzag.Compute` consumes → the H0 zigzag barcode.

**Placement:** `src/tda/ph/ConditionedZigzag.cs`, `TDA.Ph`. Functional identifiers; the sweep parameter is `T` in
prose.

## 5. Tests / exit (`TDA.Ph.Tests`)

1. **Monotone sweep = ordinary PH (the bridge).** A sweep that only *adds* edges → the H0 zigzag barcode equals
   the monotone-PH H0 of the final graph. Ties P1c back to P1a/P1b.
2. **A split the monotone reader can't see (the payoff).** A component that *merges then splits* across the sweep
   (an edge enters, then leaves) → the H0 zigzag records the split as a distinct bar; assert the merge/split
   structure. Monotone PH cannot represent this.
3. **Re-entry.** An edge that leaves and returns keeps the component bookkeeping correct (`GraphZigzag`
   multigraph-safe).
4. **Green**; the barcode comes from the built `GraphZigzag` oracle unchanged — the driver adds no reduction code.

## 6. Out of scope

- **Filled-H1 over the zigzag** (triangle filling; the Rips-vs-graph death question) → deferred (P1d).
- **Full multiparameter module** → **Z6**.
- **Directedness / magnetic flux** → **P2**.
- **`NerveDiff` frame-to-frame tracking** — already built; P1c is the oracle-seated reading, not a rewrite.

> H0 first, on the parity-complete graph-zigzag. The new code is the sweep → churn driver; the oracle is built.
> If P1c reaches for triangles, that's the filled-H1 half (P1d) — a separate seam, and its own faithfulness call.
