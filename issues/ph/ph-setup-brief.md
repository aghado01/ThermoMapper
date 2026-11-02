# Brief: PH setup, the relative backbone reading and the frame builder

**Date:** 2026-11-01. **Status:** built; scope and checks were fixed before the first run (`14ff3d0`), and the
[report](#report-2026-11-01) is below.
**Asked by the user:** "some setup work on ph stuff like the backbone persistence update", in this session.
Engine work only: no consumer is wired, and nothing here assumes one application.

## Scope

1. **The backbone update: the relative reading.** `ConditionedFiltration` documents a "return relative to the
   prior" but computes absolute persistence of the union, which matches only when K₀ has no cycles
   ([backbone-persistence.md](../backbone-persistence.md)). The design has named both flavours since 2026-06:
   backbone-at-base H(K_δ) and backbone-quotiented H(K_δ, K₀). This builds the second.
2. **P1c: the frame builder.** A sequence of graphs over stable vertex ids becomes a `ZigzagFiltration`, and the
   barcode reads back on the frames. Nothing outside `src/tda` calls the zigzag engines today.
3. **Status.** Refresh the stale rows in [dev-sequencing.md](dev-sequencing.md) and the docs in `docs/TDA`.

Out of scope: filled H₁ over the zigzag (P1d), any SPC or other consumer, monotonizing transforms, the
ultrametricity test, relative zigzag.

## Design calls

- **K₀ is the flag complex of the backbone edges, on their endpoints.** That is what the union filtration holds at
  ε₀, so the two readings describe one object. A content edge at residual 0 stays outside K₀ even though it enters
  at the same value, so the prior's predicted returns stay visible. An empty backbone gives K₀ = ∅.
- **The relative complex needs no new reduction code.** Its cells are K's cells outside K₀, in K's order, with faces
  in K₀ dropped from each boundary. `PersistentHomology`, the explicit column reducer, reads any `IFiltration`. The
  clearing and involuted paths assume absolute structure (H₀ by union-find on edges), so they are not used.
- **ε₀ is inert in the relative reading.** K₀ is present throughout the pair, so the entry point takes no `baseWeight`.
- **The convenience entry reports only exact degrees** (below the skeleton's top), so it does not carry the H₂
  truncation artifact `ComputeBarcode` documents. `RelativeFiltration` itself reports every degree.
- **Correction to the design notes: this is not VSMJ2011.** The notes say the relative flavour "needs the VSMJ2011
  machinery". Their relative modules (§2.4) are H(Xₙ, Xᵢ), the whole complex against the moving stage, with surjective
  chain maps. H(K_δ, K₀) against a fixed K₀ has injective ones and is ordinary persistence of the quotient, so it costs
  far less than those notes assumed. (Written from memory first; checked against the paper afterward, see below.)
- **Frames carry vertex sets**, because a sliding window has points entering and leaving. The convenience builders
  keep every vertex in every frame.
- **The between-frames order is an option, `FrameTransition { Union, Intersection }`, default `Union`.** Correction to
  what this session told the user on 2026-10-30, that the two orders "give different barcodes": by the Diamond
  Principle (Carlsson & de Silva 2010, §5) they should agree on the frames and differ only in transients. Check D1
  tests this.
- **The per-frame threshold is a list.** A constant θ, a warped path and a non-monotone ridge are all the caller's
  input, so the builder takes no position on γ. That frame stays the working option it is
  ([project memory, 2026-10-25]).
- **Names.** The P1c brief sketched `ConditionedZigzag`. The driver is not specific to conditioning, so it is
  `GraphFrameZigzag`, and the conditioned sweep is a short `ConditionedFiltration.BuildSweep` on top of it.
- **P1c's H₀-first scoping needs a caveat.** With a connected spanning backbone, H₀ over any sweep is one bar. A
  conditioned sweep says something in H₁ (already folded by `GraphZigzag`) or with a non-spanning backbone. H₀ is the
  payoff for frame sequences with no backbone.

## Checks, fixed before the first run

All are exact equalities of multisets, so no tolerance is needed. "Significant" means `Death > Birth` or infinite;
ties among equal births can reorder zero-length pairs only.

**Relative reading** (`tests/TDA.Ph.Tests/RelativeFiltrationTests.cs`):

| # | Check | Expected | Independent of the code under test because |
|---|---|---|---|
| R1 | Empty backbone | Relative equals absolute in every degree | K₀ = ∅ |
| R2 | Path backbone, chords | H₀ empty; H₁ equals absolute H₁ | H₁(K₀) = 0 |
| R3 | 6-cycle backbone + chord (0,3) at d | Absolute H₁ {[0,∞), [d,∞)}; relative H₁ {[d,∞)}, generator the chord | Hand |
| R4 | 4-cycle backbone + chord (0,2) at d | Absolute H₁ {[0,d)}; relative H₁ none significant, H₂ {[d,∞)} | Hand; the long exact sequence's connecting map |
| R5 | Backbone on {0,1,2}, edge (2,3) at a, separate pair (4,5) at b | Relative H₀ {[0,a), [0,b), [0,∞)} | Hand |
| R6 | Content at residual 0 closing a loop | Relative H₁ {[0,∞)}: a predicted return stays | Hand |
| R7 | Random graphs, connected spanning K₀ with cycles, content d > 0 (200 trials) | Relative H₁ equals absolute H₁ bars born after 0 | The user's identity H₁(K,K₀) ≅ H₁(K)/im H₁(K₀), read by the involuted engine on the union |
| R8 | Random graphs, arbitrary non-empty K₀ (200 trials) | Relative degrees 0–2 equal reduced homology of K ∪ cone(K₀) | H̃ₙ(K ∪ CK₀) ≅ Hₙ(K, K₀), an absolute computation on a different complex |
| R9 | A non-subcomplex predicate | Throws | Guard |

**Frame builder** (`tests/TDA.Ph.Tests/GraphFrameZigzagTests.cs`):

| # | Check | Expected | Independent because |
|---|---|---|---|
| Z1 | n=2, edges ∅, {e}, ∅ | H₀ frame bars {[0,2], [0,0], [2,2]} | Hand decomposition |
| Z2 | n=2, edges {e}, ∅, {e} (re-entry) | H₀ frame bars {[0,2], [1,1]} | Hand decomposition |
| Z3 | Random nested frames (100 trials) | Frame bars equal `PersistentHomology` on the edge-birth filtration; H₀ deaths equal union-find merge frames | Ordinary persistence; Kruskal |
| Z4 | Random churn with re-entry (30 trials, n=5, 5 frames) | `GraphZigzag` Reference and Fast, `FastZigzag`, `ZigzagBarcodeNaive` agree, step-level bars, degrees 0–1 | Four engines |
| Z5 | Random churn with vertices entering and leaving (30 trials) | As Z4 | Four engines |
| Z6 | Single frame | Each component one bar [0,0], each independent cycle one H₁ bar [0,0] | Hand |
| Z7 | Conditioned sweep: path backbone, chord present in frames 1, 2, 4 of 5 | H₀ {[0,4]}, H₁ {[1,2], [4,4]} | Hand |
| Z8 | Edge field on one substrate with varying θ | Admits by `≥ θ` on j>i slots, ignores mirror slots | Hand |
| Z9 | `FromGraphs` and `FromEdgeSets` on the same frames | Identical steps | Same input |
| Z10 | Bad frames (edge endpoint absent, id out of range) | Throw | Guard |

**Diamond Principle** (same file):

| # | Check | Expected |
|---|---|---|
| D1 | Random churn (200 trials): `Union` and `Intersection` | Identical frame bars |
| D2 | Same trials: transients | H₁ transients of `Union` at t match H₀ transients of `Intersection` at t; `Union` has no H₀ transient, `Intersection` no H₁ transient |
| D3 | Hand: {(0,1),(1,2)} then {(0,2),(1,2)} | One H₁ transient under `Union`, one H₀ transient under `Intersection`, frame bars H₀ {[0,1]} under both |

D1's bijection is the theorem as sourced. D2's degree shift is from memory of the theorem: what I could read confirms
the bijection and the section, not the shift. It is a prediction, and a failure gets recorded as one.

## Notes

The brief's hand expectations were worked before the code exists. If one is wrong, the report says so and why,
rather than the expectation being edited quietly.

## Report (2026-11-01)

**What landed.** `RelativeFiltration` and `ConditionedFiltration.BuildRelative` / `ComputeRelativeBarcode` (the
backbone update); `GraphFrameZigzag`, `GraphFrameFiltration` and `ConditionedFiltration.BuildSweep` (P1c). Docs:
`docs/TDA/PersistentHomology.md`, `docs/TDA/ZigzagPersistence.md`. No engine changed, and nothing is wired to a
consumer. `ComputeBarcode` keeps its behaviour.

**Results.** Every check above passed on its first run, with no hand expectation changed, including D2's degree shift,
which was only a prediction. `TDA.Ph.Tests` is 375/375 and `TDA.Mapper.Tests` 28/28.

**The sweeps reach what they claim.** Counts over the trials, after the first run raised the question:

| Check | Reached |
|---|---|
| R7 | 130 backbone bars quotiented; 201 relative H₁ bars, 107 of them finite |
| R8 | Bars in H₀/H₁/H₂: 169/155/167; 114 non-spanning backbones; 90 with K₀ triangles |
| Z3 | 210 H₀ merges, 270 H₁ bars |
| Z11 | 200 frame bars, 18 transients |
| D1, D2 | 791 transients |

**The tests can fail.** Two deliberate breaks, run once and reverted:

- K₀ losing its triangles failed R8 alone. That is the expected reach: R7 reads H₁, which K₀'s triangles don't touch,
  so the cone is the only guard on this choice.
- The read-back's level window off by one failed Z1, Z2, Z3, Z11 and D1.

**What the review changed, after the first green run:**

- **A bug the checks missed.** `GraphRips` stops at triangles for any `maxDimension ≥ 2`, so `ComputeRelativeBarcode`
  at `maxDimension: 3` would have reported the truncated H₂ it promises to leave out. Fixed; R11 pins it.
- **The read-back had no independent reference beyond hand cases and nested frames.** Z11 now decomposes the module on
  frames and between-frames complexes directly with `ZigzagMapBarcode`, with no single-cell refinement and no read-back.
  It agrees for both transitions, frame bars and transients.
- **R10** pins a finite bar's representative, in the ambient's indices. By hand: path 0–4 with chords (0,4), (0,2),
  (2,4) gives H₁ [1,3) whose relative cycle is the chord (0,4) alone.
- **Z4 and Z5 run 10 trials, not 30.** They ran green at 30, at 110 s against a suite of about 20 s, almost all of it
  the naive oracle. D1 and Z11 cover the same ground far more cheaply.

Checks added after the first run, so not pre-registered:

| # | Check | Expected |
|---|---|---|
| R10 | Finite relative bar, hand | [1,3), cycle {(0,4)} |
| R11 | `ComputeRelativeBarcode` at `maxDimension` 3 and 1 | No H₂; H₀ only |
| Z11 | Random frames with vertex churn (40 trials, both transitions) | Read-back equals direct decomposition of the coarse module |

## Checked against the sources afterward

The user pointed me at their library (`graveyard/bibliotecha/bibliotecha/compendia/ph`), which settled both claims this
brief had left unread.

- **VSMJ2011's relative modules** are H(Xₙ, Xᵢ) with i moving (§2.4), and the paper's remark in that section says their
  chain maps are surjective. Confirmed.
- **The Diamond Principle's degree shift.** Carlsson, de Silva and Morozov 2009 (`CDSM2009.md`) pair a transient of degree
  p at the intersection with one of degree p+1 at the union, and every other bar keeps its degree. Confirmed. The interval
  rules are in Theorem 7 of Dey & Hou 2022 (`2204.11080v2.md`).

## What VSMJ's Prop 2.4 suggested: R12

Their Prop 2.4 says their absolute and relative barcodes carry the same information, with a degree shift for finite bars.
A fixed K₀ has an analogue, which I derived after reading it and added as a check after the fact:

| Absolute bar, run with K₀ strictly first | Relative bar |
|---|---|
| Born outside K₀ | The same bar |
| Born in K₀, killed outside it at δ | [δ, ∞), one degree up |
| Born in K₀, essential or killed within K₀ | None |

**The argument.** With K₀'s rows lowest in the boundary matrix, deleting them leaves every pivot outside K₀ as it was. A
column whose pivot was in K₀ reduces to zero, so it becomes a birth. That cell was negative in the absolute run, so no
column has it as a pivot, and the new class is essential.

**The check.** R12 passes on 200 random trials: 427 bars kept, 64 moved up, 1184 dropped. The complexes stop at triangles,
so no higher cell could have killed a moved-up class, and R12 alone can't test that one.

**The proof.** The user asked for a formalization, since this was derived and not cited:
`lean/Lemmas/RelativePersistence.lean`, no apologies, standard axioms only, for any field. It turned out simpler than the
argument above: any reduction D = R·U restricts to a reduction of the relative matrix, which gives the table directly.

- `relative_barcode` is the statement, for any reduction of either matrix.
- `isLow_iff_of_reductions` shows the pairing doesn't depend on the reduction.
- `column_zero_of_isLow` is the only place D·D = 0 enters, and it covers the "essential" claim R12 couldn't.
- `low_degree` gives the degree.
- A worked example over `ZMod 2` shows the hypotheses can be met.

It covers the mathematics, not the C#: `RelativeFiltration` is still checked by tests. I haven't looked for whether this
is stated in the literature, so it makes no claim to being new.

**What it means for the design notes.**

- The notes call absolute and relative "two questions", the relative one giving features "orthogonal to" the backbone.
  With K₀ first they are one barcode read two ways. The relative reading adds a reading, not information. That is not an
  objection to it: it says which bars the prior already explains.
- `ComputeBarcode` does not run K₀ first. K₀ ties at 0 with the vertices outside it and with content at residual 0, so
  its output alone doesn't determine the relative barcode. A run with K₀ strictly first does.
- A second reduction is not strictly needed, so a faster path exists in principle. It is not built: the involuted path's
  H₀ union-find inverts the elder rule when vertex births differ (`docs/TDA/PersistentHomology.md`), which K₀-first needs.

## Status of the PH-for-SPC ideas

Nothing below was decided except where it says so. Application framing stays out of the engine.

| Idea | Status | Decision |
|---|---|---|
| Backbone persistence, relative reading | **Built** | User, 2026-11-01: do it, here |
| Frames → zigzag (P1c) | **Built**, no consumer | Same request, as setup work |
| SPC over temperature as zigzag | Not wired. The builder takes `CoMembership.G` as a field on the coupling graph as it stands. At fixed θ with J ≥ 0, the exact stack only gains edges, so zigzag reads the raw estimate, non-nesting cuts, a θ*(T) ridge or signed couplings | Open |
| (β, θ) bifiltration and transforms | Not built; `lean/Enthymemata/Bifiltration.lean` is empty. "Not a bifiltration" holds for the cumulative field, where (T, θ) reduces to τ = T·κ(θ), and not for the sampler, whose edge curves cross | Open |
| Ultrametricity test | Not built; placement undecided | Open |
| Per-sweep partitions for PH | Not built. A zigzag over thinned draws is order-dependent; resampling them per temperature is the better route to barcodes with uncertainty | Open |

**Next on the engine side:** filled H₁ over frames landed the same day as P1d
([p1d-filled-h1-brief.md](p1d-filled-h1-brief.md)). A relative reading over frames is the open seam.

**To forward to the Bayesian clustering session:** θ read as a Binder cost ratio applies to each edge's decision, while
the θ-graph's components are the closure, which is single linkage and not the Binder optimum. The robust window and H₀
bar length both measure an interval of stability, one per partition and one per cluster.

**Older notes this corrects, left as they are:** the MarkBrain synthesis brief and `backbone-conditioned-persistence.md`
still say relative persistence needs the VSMJ2011 machinery, and still show it unbuilt.
