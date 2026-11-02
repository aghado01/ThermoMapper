# Brief: P1d, filled H₁ over frames

**Date:** 2026-11-01. **Status:** built; scope and checks were fixed before the first run (`72763fb`), and the
[report](#report-2026-11-01) is below.
**Asked by the user:** "Yes, build P1d", after [ph-setup-brief.md](ph-setup-brief.md) named it as the next engine step.
Engine work only: no consumer is wired.

## Scope

`GraphFrameZigzag` gains `FrameComplex { Graph, Flag }`. `Flag` makes each frame its flag complex through triangles,
read by `FastZigzag`, so a cycle can die because triangles fill it, and a hole can open because they leave. `Graph`
stays the default and its steps do not change.

Out of scope: dimensions above triangles, relative reading over frames, any consumer, new engine code.

## The question P1c left open is not a choice

[P1c §3](p1c-nonmonotone-zigzag-brief.md) deferred "does a return die by triangle-fill or by edge-departure?" as P1d's key
decision. It needs no decision. Both happen, and the arrow that ends a bar says which: `FrameBar` reports it per bar.

| Event | As a birth | As a death |
|---|---|---|
| `Arrival` (cells enter) | a vertex enters; an edge closes a cycle | a merge; triangles fill a cycle |
| `Departure` (cells leave) | an edge leaves and splits a component; triangles leave and open a hole | an edge leaves and breaks a cycle; a vertex leaves |
| `Unobserved` | the bar is at the first frame | the bar is at the last frame |

A bar at the sequence's edge has no observed start or end there, so it says so instead of reporting build-up or teardown.

**A prediction this rests on.** The event should not depend on `FrameTransition`. Under the Diamond Principle's interval
map, a bar ending at frame t under `Union` by a forward arrow (t → union) ends at the intersection under `Intersection`,
again by a forward arrow, and likewise for backward arrows. I derived that here and have not read it in the theorem, so
D1 tests it.

## Design calls

- **`Union` is Flag(G_t) ∪ Flag(G_{t+1}), the union of the two complexes.** The alternative, Flag(G_t ∪ G_{t+1}), adds
  triangles neither frame holds, so a class could be filled by a triangle that exists at no frame. Union of complexes
  also keeps Mayer–Vietoris, since Flag(G_t) ∩ Flag(G_{t+1}) = Flag(G_t ∩ G_{t+1}), so the Diamond Principle still
  applies. The alternative is a different zigzag and is not built.
- **Through triangles only.** H₀ and H₁ are exact. H₂ has nothing above to fill it, so it is the same truncation artifact
  `docs/TDA/PersistentHomology.md` documents: read `maxDimension: 1`.
- **The graph engines do not read `Flag`.** They assume a 1-complex. Z12 pins that they refuse it.
- **The event comes from the filtration's own step,** not from the engine's `IntervalEnd`, so it is the same for every
  engine.
- **Representatives are not addressed.** `FastZigzag` gives none for bars born by a deletion, and its docs mark them
  unverified above H₀. So a filled reading has no complete "which cycle" yet.

## Checks, fixed before the first run

All are exact equalities of multisets. Tests: `tests/TDA.Ph.Tests/GraphFrameZigzagTests.cs`. C4 is the 4-cycle
(0,1), (1,2), (2,3), (0,3); the chord is (0,2), which makes two triangles. Hand cases run under both transitions.

| # | Frames | Expected (degree, first, last, start, end; U, A, D for `Unobserved`, `Arrival`, `Departure`) |
|---|---|---|
| F1 | C4, C4 + chord; `Flag` | H₀ (0,0,1,U,U); H₁ (1,0,0,U,A): filled |
| F2 | C4, C4 − (0,3); `Flag` | H₀ (0,0,1,U,U); H₁ (1,0,0,U,D): broken |
| F3 | C4 + chord, C4; `Flag` | H₀ (0,0,1,U,U); H₁ (1,1,1,D,U): a hole opens |
| F4 | C4, C4 + chord, C4; `Flag` | H₀ (0,0,2,U,U); H₁ (1,0,0,U,A), (1,2,2,D,U) |
| F4′ | Same; `Graph` | H₀ (0,0,2,U,U); H₁ (1,0,2,U,U), (1,1,1,A,D): the same sequence, read unfilled |
| F5 | `BuildSweep`, path 0–3, chord (0,3) in frames 0–2, chord (0,2) in frame 1; `Flag` | H₀ (0,0,3,U,U); H₁ (1,0,0,U,A), (1,2,2,D,D): one return dies both ways |

| # | Check | Expected | Independent because |
|---|---|---|---|
| Z12 | `GraphZigzag`, Reference and Fast, on `Flag` frames | Throw | Guard |
| N1 | Random nested frames, `Flag` (100 trials) | Frame bars in H₀, H₁ equal `PersistentHomology` on `GraphRips` with each edge born at its first frame | Ordinary persistence |
| B1 | Random churn, both complexes (100 trials) | At every frame, bars covering it equal the Betti numbers of that frame's own complex | Each frame reduced alone, no zigzag |
| E1 | Random churn with re-entry, `Flag` (10 trials) | `FastZigzag`, `ZigzagBarcodeNaive`, `ZigzagMapBarcode.ComputeFromZigzag` agree at step level, H₀–H₂ | Three engines |
| C1 | Random frames with vertex churn, both complexes, both transitions (40 trials) | Frame bars, transients and events equal a direct decomposition of the module on frames and between-frames complexes | No single-cell refinement, no read-back, triangles enumerated by the test |
| D1 | Random churn, `Flag` (200 trials) | Frame bars with events identical across transitions | Prediction above |
| D2 | Same trials | A transient of degree d under `Union` matches one of degree d−1 under `Intersection` at the same position; `Union` has none in H₀, `Intersection` none in H₂ | Diamond Principle, H₁ ↔ H₀ and H₂ ↔ H₁ |
| H1 | Random churn (100 trials) | H₀ frame bars with events, `Flag` through `FastZigzag`, equal `Graph` through `GraphZigzag` | Triangles do not change connectivity |

The Z11 check from [ph-setup-brief.md](ph-setup-brief.md) becomes C1, and now also covers events.

If a hand expectation turns out wrong, the report says so and why.

## Report (2026-11-01)

**What landed.** `FrameComplex { Graph, Flag }` on every `GraphFrameZigzag` builder and on `BuildSweep`, and
`FrameBar.Start` / `End` as `FrameEvent { Unobserved, Arrival, Departure }`. No engine changed, and `Graph` frames build
the same steps as before. Docs: `docs/TDA/ZigzagPersistence.md`, `docs/TDA/PersistentHomology.md`.

**Results.** `TDA.Ph.Tests` is 394/394 and `TDA.Mapper.Tests` 28/28. Every hand expectation held on its first run. The
prediction held too: a bar's start and end events are the same under both transitions.

**The first run had one failure, and it changed the tests.** C1's reach guard failed under `Flag`: the read-back agreed
everywhere, but over 40 trials no cycle was filled and no hole opened. Independent random frames share too few edges for
a cycle to survive into the next frame, so it can't be filled there. D1's prediction was passing without being tested on
the events it is about.

The fix was `EvolvingFrames`, where each frame is the last with a few vertices and pairs flipped. C1 now runs 80 trials
and D1 keeps 200, each half independent and half evolving. No expectation changed, only what the trials reach:

| Check | Reached after the fix |
|---|---|
| D1 | H₁ bars ended by arrival (filled): 37; by departure (broken): 179; started by departure (hole opened): 38 |
| D2 | Union transients: H₁ 287, H₂ 32 |
| C1, `Flag` | 224 frame bars, 106 transients, 5 fills, 6 holes opened |
| C1, `Graph` | 577 frame bars, 106 transients, 0 fills, 0 holes opened |
| N1 | 38 H₁ bars ended before the last frame |

The zeros under `Graph` are asserted, not just observed: a graph cycle can't be filled, and no cell leaving can open one.

**The tests can fail.** Two deliberate breaks, run once and reverted:

- The events read from the wrong step failed F2, F4, C1 (both complexes) and D1.
- Tearing the whole frame down under `Intersection`, so the between-frames complex is empty, failed F1–F5 under
  `Intersection`, C1 and D1.

**What the review changed.** `FrameTransition`'s docstrings said a transient is a cycle under `Union` and a component
under `Intersection`. With triangles it can also be a void and a cycle. Fixed.

**Cost.** The suite took 51 s on one run and 20 s on the next, at `xUnit.MaxParallelThreads=4`, against 19 s before these
checks. The spread is machine load, so the checks add little.

**Checked against the source afterward.** Theorem 7 of Dey & Hou 2022 (`2204.11080v2.md`, in the user's library) gives the
interval rules, for a diamond at position j. Under them a bar ending at j−1 through the union ends at j through the
intersection, and a bar ending at j ends at j−1. The first pair both die on a forward arrow and the second on a backward
one, and starts map the same way. So the prediction is a direct reading of the theorem, not only my derivation.

## Not verified

- **Representatives.** `FastZigzag` gives none for bars born by a deletion, and its docs mark representatives above H₀
  unverified. So a `Flag` bar reports that it was filled, but cannot yet reliably name the cycle.
