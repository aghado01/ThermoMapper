# Ascent well-posedness — the "conditioning pays twice" stack

> Provenance: the 2026-06-12 resolution-layer discourse (modal ascent's two applications; this
> stack covers **application #2** — the graph mode-seek periphery policy, the alternative to
> Domany's peripheral capture). The selector-side walk math (per-leaf landscape integrals on the
> dendrogram) is a separate future note once R1 lands. Siblings: `fable-BARS-lemma.md`,
> `fable-bifiltration-lemmas.md`, `PKWangAudit.md` (the critical register), and
> `CumulativeField/design.md` (the recovered construction).
> Map context: `../hierarchy/dendrogram-integration-map.md` § R1 refinements.

A different genre than the A/B register. Those are "surprising exact reduction, apparatus rendered
vestigial." This stack is **preconditions made theorems**: the graph engine's knobs (mutual filter,
MST repair, deterministic tie-breaks, gauge discipline) were each justified for the *dynamics*; the
lemmas below show the same knobs deliver well-posedness for the *resolution stage* — termination,
coverage-with-honest-abstains, orientation-freedom, gauge robustness. The composite claim:
**conditioning pays twice**, and that's provable, not asserted. Everything here is finite
combinatorics + order theory — zero Mathlib infrastructure gaps, which makes it the cheapest fish
in the queue (cheaper than the BARS pair, which at least needed Fermat).

## Setup — the ascent operator, defined so the lemmas are forced

The landscape codomain needs only a **linear order** — that's not a convenience, it's L4's thesis
baked into the signature (ascent is an *ordinal* consumer; it cannot even see cardinal structure).

```lean
import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Order.Lexicographic

variable {n : ℕ} (G : SimpleGraph (Fin n)) [DecidableRel G.Adj]
variable {α : Type*} [LinearOrder α] (L : Fin n → α)

/-- Tie-broken height key: landscape first, lower index wins a plateau.
    Injective by the second component — the plateau guard IS the order. -/
def key (i : Fin n) : α ×ₗ (Fin n)ᵒᵈ := (L i, OrderDual.toDual i)

/-- One ascent step: the key-maximum of the closed neighborhood.
    Deterministic, seed-free — the resolution stage adds zero randomness. -/
def step (i : Fin n) : Fin n :=
  ((insert i (G.neighborFinset i)).image id).max' (by simp) -- argmax under key; sorry-shaped sugar
  -- real def: Finset.max' over the closed nbhd ordered by key; uniqueness from key injectivity
```

(Plateaus: two nodes with equal `L` are ordered by index, so `key` is a strict total order on any
subset and `max'` is unique. The same lowest-index discipline as the SW giant tie-break — one
convention, used twice, and no RNG anywhere in resolution.)

## L1 — termination, forest, basins (the foundation)

```lean
/-- Along a non-trivial step the key strictly increases. -/
theorem key_strictMono_step (i : Fin n) (h : step G L i ≠ i) :
    key L i < key L (step G L i) := by
  sorry -- step is the closed-nbhd max; if it differs from i, its key beats i's

/-- Ascent terminates: at most n steps reach a fixed point. Basins are well-defined. -/
theorem ascent_terminates (i : Fin n) : ∃ k ≤ n, step G L (step G L)^[k] i = (step G L)^[k] i := by
  sorry -- strictly increasing key sequence in a finite linear order; pigeonhole

def basin (i : Fin n) : Fin n := sorry -- σ^[n] i, or via WellFounded.fix on key
```

**Obligations:** (a) `key` injective ⇒ unique argmax; (b) `key_strictMono_step`; (c) pigeonhole
termination; (d) the trajectory relation `i → step i` is a forest rooted at the fixed points
(modes); (e) `basin` is total ⇒ **basins partition V** — coverage needs no connectivity.
**Knob consumed:** the deterministic tie-break. **Why load-bearing:** quick-shift-family papers
hand-wave this ("follow density uphill"); on plateaus and ties the naive operator cycles or forks.
The tie-break-as-order makes the operator a function and the proof mechanical.

## L2 — abstains are landscape phenomena (connectivity's real payoff)

Two statements. First, the free one:

```lean
/-- Trajectories never descend: ascent cannot cross a valley. -/
theorem ascent_never_descends (i : Fin n) (k : ℕ) :
    key L i ≤ key L ((step G L)^[k] i) := by
  sorry -- monotone composition of key_strictMono_step / refl
```

Second, what `EnsureConnected` (MST repair) actually buys — NOT coverage (L1 gives that), but
**informativeness of failure**: in a connected graph, if a point's basin terminates at an
*unlabeled* mode (abstain), a route to every labeled core exists — so the abstain certifies that
every such route crosses a watershed (a basin-boundary edge), i.e. a genuine valley separates the
point from every core. Disconnection can never masquerade as a landscape feature.

```lean
theorem abstain_certifies_watershed (hconn : G.Connected)
    (i c : Fin n) (hbasin : basin G L i ≠ basin G L c) :
    ∀ w : G.Walk i c, ∃ e ∈ w.edges, -- some edge of every walk crosses basins
      basin G L e.fst ≠ basin G L e.snd := by
  sorry -- basins partition V (L1.e); a walk between different classes crosses at an edge;
        -- walk existence is hconn. Finite induction along the walk.
```

**Knob consumed:** `RepairKind.MstMin` / `EnsureConnected`. **Why load-bearing:** it's the formal
version of "ascent fails closed" — an abstain is a statement about the landscape (frustration
valley on every route), never an artifact of a torn graph. Capture has no analog of either theorem.

## L3 — symmetry: the precondition L1 quietly assumed

`SimpleGraph` is symmetric by type; the engine-side content is that **both** filter knobs produce
it from the asymmetric KNN relation:

```lean
/-- Mutual (∩) and OR-rule (∪) filters both symmetrize a directed relation. -/
theorem filter_symmetrizes (R : Fin n → Fin n → Prop) :
    Symmetric (fun i j => R i j ∧ R j i) ∧ Symmetric (fun i j => R i j ∨ R j i) := by
  sorry -- one line each
```

**Knob consumed:** `FilterKind.MutualKnn` / OR-rule. **Why it's here at all:** the ascent operator
is only *constructible* on a symmetric adjacency (orientation-free flow); raw KNN is not one. The
lemma is trivial; naming it keeps the stack honest about where `SimpleGraph` came from.

## L4 — gauge invariance (the sleeper)

```lean
/-- Basins depend on the landscape only through its order:
    any strictly monotone regauging changes nothing. -/
theorem basin_gauge_invariant {β : Type*} [LinearOrder β]
    (φ : α → β) (hφ : StrictMono φ) :
    basin G (φ ∘ L) = basin G L := by
  sorry -- StrictMono preserves < and (via injectivity) =; key comparisons unchanged;
        -- step unchanged pointwise; basin unchanged by funext + induction
```

**Knob consumed:** the sink/gauge choice (the metric–measure factorization's measure side).
**Why load-bearing:** it splits the resolution family into **ordinal consumers** (ascent — provably
invariant under per-sink normalizations, α-family rescalings, any order-preserving gauge) and
**cardinal consumers** (threshold cuts, EOM mass integrals — NOT invariant). LocalField-as-α0
"re-ranks, doesn't rescale" lands exactly here: re-ranking is the *only* thing ascent can see. Any
A/B between policies must therefore declare gauges for the cardinal side only — the ordinal side is
gauge-free by theorem, one less knob to litigate.

## L5 — capture ≁ ascent: the one-witness separation

The formal teeth behind the P5 A/B (same genre as the BARS note's `argmax_expectation_noncommute` —
an explicit finite witness, near-`decide`-able):

```lean
/-- There is a configuration where edge-greedy capture crosses a watershed
    that height-greedy ascent provably cannot. The policies are inequivalent. -/
theorem capture_crosses_where_ascent_cannot :
    ∃ (L : Fin 4 → ℚ) (g : Sym2 (Fin 4) → ℚ) (i : Fin 4),
      -- path graph 0–1–2–3; L = [3, 1, ½, 4]; g(0,1)=0.2, g(1,2)=0.9, g(2,3)=0.3
      capture g i ∉ basinSet L i ∧ step _ L i ∈ basinSet L i := by
  sorry -- node 1: max-g neighbor is 2 (0.9) → capture unions across the valley floor;
        -- ascent from 1 goes to 0 (key-max of {0,2}); basins {0,1} | {2,3}; watershed (1,2)
```

**Why load-bearing:** it certifies the periphery A/B is theoretically contentful, with the failure
mode named — a strong saddle edge is sufficient for capture to fuse across a frustration valley,
and `ascent_never_descends` (L2) is the proof ascent can't. Prediction it pins for the experiment:
capture/ascent disagreements concentrate on watershed nodes and correlate with the fuzzy-entropy
frustration meter.

## Knob → property map

| graph-engine / config knob | lemma | property delivered to resolution |
|---|---|---|
| deterministic tie-break (index order; no RNG) | L1 | termination, forest, unique basins, plateau-safe |
| `FilterKind.MutualKnn` / OR-rule | L3 | symmetric adjacency — the operator is constructible |
| `EnsureConnected` (`RepairKind.MstMin`) | L2 | abstains certify valleys, never disconnection |
| sink + gauge choice (measure side) | L4 | basins gauge-free; ordinal/cardinal consumer split |
| — (no knob; intrinsic) | L5 | capture/ascent inequivalence, witnessed |

## Non-goals (the circularity guard)

"Well-conditioned graph ⇒ no spurious modes of L" is **not in the stack** and must not sneak in:
it's a theorem about the generative setting (data → graph faithfulness), not a retro-formalizable
invariant of the code. Proving it from the engine's own assumptions would be the circular-validation
sin — the stack certifies the resolution is *well-posed given* the graph; whether the graph is
*faithful* answers to external oracles (the parity battery) and cross-sink triangulation, never to
these lemmas.

## Fishability verdict

All five are finite combinatorics over a linear order — no analysis, no measure theory, no Mathlib
gaps. Closeability ranking: **L3, L4** (one-liners once `key` is defined) → **L5** (explicit
witness over `ℚ`, decidable arithmetic) → **L1** (mechanical pigeonhole/well-founded, the only real
work is defining `step` cleanly via `Finset.max'`) → **L2** (walk induction over the L1 partition).
Lead with L1 (everything else imports its definitions), but L4 is the statement worth having
*first* if the order is forced — it retires the gauge-litigation for half the resolution family in
four lines.
