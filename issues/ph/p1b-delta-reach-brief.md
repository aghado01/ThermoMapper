# Brief — P1b: the Δ reach axis (the prior gains a tunable horizon)

**Status:** build ticket — **local `issues/ph` (transient; archive to MarkBrain when P1 closes).** Second
increment of synthesis §6 **P1**. **Scope in one line:** add the **Δ reach bound** (`|τ| ≤ Δ`) as a second,
*monotone* filtration axis over P1a's residual band, making `K_{δ,Δ}` a genuine bifiltration read by **monotone
Δ-slices** (ordinary PH). **The non-monotone zigzag reader is split out to P1c; the full module stays Z6.**

**Background (MarkBrain — canonical):** synthesis §2 (the `τ ≤ Δ` cap), §6 P1, and the axes section
(`Δ` = prior reach / horizon) — [opus-brief-conditioned-persistence-synthesis.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/opus-brief-conditioned-persistence-synthesis.md).
Builds on **P1a** — [p1a-residual-prior-brief.md](p1a-residual-prior-brief.md) (committed `c91143d`).
**Sequencing:** [dev-sequencing.md](dev-sequencing.md).

---

## 1. The split (decided in scoping — the synthesis bundled three sizes)

Synthesis §6 P1 reads as "Δ axis + non-monotone zigzag reader + full module." Those are three very different
builds, so P1 decomposes:

- **Δ reach axis** — a monotone second parameter: a filter on the prior + PH slices. Tight, reuse-maximal.
  **This brief (P1b).**
- **Non-monotone zigzag reader** — building a `ZigzagFiltration` over a splitting/merging sweep (trajectory,
  sliding window, `T`-sweep) and reading it through the built zigzag oracle. Bigger, telos-connected
  (persistent-Mapper-over-`T`), and it touches the **multi-track zigzag engine** → **its own brief, P1c.**
- **Full multiparameter module** (`K_{δ,Δ}` rank invariant / fibered barcode) → **Z6** (unchanged).

**P1b = the Δ axis only.** It completes P1a's single-`δ` filtration into the `(δ,Δ)` bifiltration, read by
monotone slices — no zigzag, no module.

## 2. What Δ is — pin it, "reach" is ambiguous

Per synthesis §2, Δ bounds the **prior magnitude**: an edge `(i,j,τ)` is admitted only when `|τ| ≤ Δ`. Δ is the
prior's **horizon** — how far ahead a prediction is allowed to reach. It is **not** the backbone span `|i−j|`,
and **not** the residual (that is `δ`). Growing Δ **adds** edges (longer-reach predictions), so Δ is
monotone-increasing, exactly like `δ`.

- small Δ → only short-reach predictions admitted (local priors);
- large Δ → long-reach predictions too.

`(δ, Δ)` is a bifiltration, **both axes monotone**. A return that relies on a long-reach prior edge appears
**only above that edge's Δ threshold** — the reach axis turns "how far the prior had to reach to explain this
return" into a filtration coordinate.

## 3. Reuse-maximal — a reach filter on P1a's producer

The one new thing is a **reach filter** applied before P1a's residual computation:
- `ResidualPrior.ResidualEdges` gains a `reachBound` (drop `(i,j,τ)` with `|τ| > reachBound`); a **Δ-slice at
  fixed δ is P1a with the filtered prior** → the existing barcode path, unchanged;
- the bifiltration is exposed as a **slice family**: a monotone Δ-grid → a nested sequence of residual-edge sets,
  each read by the existing reducer.

**Placement:** extend `ResidualPrior` (a `reachBound` parameter) and/or a thin `ReachAxis` companion in
`src/tda/ph/`. Functional identifiers; `Δ`/`δ` live in prose.

## 4. API (shape, not gospel)

```csharp
// A — a bound on the existing producer (default +inf recovers P1a exactly):
ResidualPrior.ResidualEdges(observations, prior, symmetry,
                            double reachBound = double.PositiveInfinity);   // admit only |tau| <= reachBound

// B — the slice family for the (delta, Delta) bifiltration:
ReachAxis.Slices(observations, prior, IReadOnlyList<double> reachGrid)
    -> IReadOnlyList<(double reach, IReadOnlyList<(int,int,double)> edges)>;
```

## 5. Tests / exit (`TDA.Ph.Tests`)

1. **`reachBound = +∞` subsumes P1a.** The unbounded producer returns P1a's residual edges exactly — same
   barcode. Δ generalizes, doesn't perturb (the same subsumption discipline as τ≡0 → P0).
2. **A reach-gated return.** Two chords: one short-reach (`small |τ|`), one long-reach (`large |τ|`). At small Δ
   only the short-reach return appears; raising Δ past the long-reach edge's `|τ|` admits its return — assert
   `β₁` steps up at the Δ threshold.
3. **Monotone in Δ.** Growing Δ only *adds* edges → the admitted sets are nested → the slice family is monotone
   (no edge ever leaves as Δ grows). This is the property that keeps monotone-slice = ordinary PH, and it is the
   exact boundary with P1c: the non-monotone case comes from a *time-like* parameter, never from Δ.
4. **Green**; slices route through the existing reducer unchanged.

## 6. Out of scope (P1c, Z6, P2)

- **Non-monotone / zigzag reader** (trajectories, sliding windows, `T`-sweeps that split *and* merge) → **P1c**.
  Grounds in the built zigzag stack — currency `ZigzagFiltration`; oracles `ZigzagBarcodeNaive` / `FastZigzag` /
  `ZigzagMapBarcode` (arbitrary cell maps). **Verify against the zigzag FRONTIER doc + `git log`/glob before
  building** — that engine is multi-track and its per-rung briefs lag committed code
  ([zigzag-frontier.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/zigzag-engine/zigzag-frontier.md)).
- **Full multiparameter module** (`K_{δ,Δ}` rank invariant) → **Z6**.
- **Directedness / magnetic flux** → **P2**.

> A reach bound is a filter; a Δ-slice is P1a. If P1b needs more than the reach filter + slice family + a test,
> the zigzag reader (P1c) has leaked in — that is the non-monotone case, and it is a separate engine seam.
