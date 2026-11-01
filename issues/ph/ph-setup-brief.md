# Brief: PH setup, the relative backbone reading and the frame builder

**Date:** 2026-11-01. **Status:** scope and checks fixed before the first run; the report is appended below.
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
  machinery". Their relative objects are H(K, K_i), against the moving stage of the filtration. H(K_δ, K₀) against a
  fixed K₀ is ordinary persistence of the quotient, so it costs far less than those notes assumed. (Their definition
  is from memory; the abstract I could read is consistent with it and does not state it.)
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
