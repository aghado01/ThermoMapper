# Brief: closing the PH setup's loose ends

**Date:** 2026-11-02. **Status:** chipped, not started.
**Asked by the user:** "lets chip off to close the loose ends", the four entries under "PH setup loose ends" in
[the Shelf](../SHELF.md). All four are decided: build them. Engine work only, with no consumer wired, and nothing here assumes
one application.

**Read first:** [ph-setup-brief.md](ph-setup-brief.md) (the relative reading, its R-checks and the proved table),
[p1d-filled-h1-brief.md](p1d-filled-h1-brief.md), `docs/TDA/ZigzagPersistence.md` (the `FastZigzag` section and Status) and
`lean/Lemmas/RelativePersistence.lean`. Briefs lag code, so verify against git and grep before trusting one.

## Order

One concern per commit, in this order, since each step leans on the one before it or is cheap enough to clear first.

### 1. The older notes

`MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/`: `backbone-conditioned-persistence.md` (lines 34-36,
80-82, 90) and `opus-brief-conditioned-persistence-synthesis.md`. They say relative needs VSMJ2011 and show it unbuilt.

- A dated note goes beside each stale claim, with the original text kept. Fixed K₀ is ordinary persistence of the quotient, built in
  `fd311e8`, and VSMJ2011's H(Xₙ, Xᵢ) is the moving case.
- "Absolute vs relative = two questions" is the user's interpretive framing and stays. The note adds that with K₀ first they
  are one reduction read two ways, so the two questions differ in what they ask, not in what they compute.
- MarkBrain's commits are the user's. Leave the edits uncommitted there and name them in the report, along with the "Frame
  sequences" row already uncommitted in `zigzag-engine/zigzag-frontier.md`.

### 2. The H₀ elder rule

`PersistenceClearing.ComputeH0` (`src/tda/ph/PersistenceClearing.cs:74`) keeps the component with the larger birth, an engine
bug documented at `docs/TDA/PersistentHomology.md:157`. It blocks step 3, and it is wrong for any filtration with distinct vertex
births. `PersistentInvolutedHomology` and `PersistentCohomology` both read it.

Checks:
- Both engines equal `PersistentHomology` as multisets of significant bars, every degree, on random simplicial filtrations with
  distinct vertex births. A reach guard counts the trials where the inverted rule would have differed.
- Filtrations with all vertices born at 0 keep their bars unchanged. Any test that pinned the inversion gets reported.
- The differential may turn up further non-Rips deviations in these engines. Fix those that share this cause, report the others.
- Remove the deviation from `docs/TDA/PersistentHomology.md` (lines 37 and 157).

### 3. The relative barcode from one absolute run

The table is proved: `relative_barcode`, `low_degree`. Build the reading as a code path, with K₀ strictly first.
`RelativeFiltration.ComputeBarcode` stays as the reference.

Three things to settle while building:
- **Strictly first.** K₀ ties at 0 with the other vertices and with content at residual 0. R12
  (`RelativeBarcode_IsARereadingOfTheAbsoluteRunWithK0First`, `tests/TDA.Ph.Tests/RelativeFiltrationTests.cs:205`) already
  builds this order, so start there.
- **Row dispatch.** The table needs to know each bar's birth and death cell, or to recover it from the order. `ComputeH0`'s bars
  carry values only.
- **Representatives.** For a moved-up class the relative cycle is the chain V_j, with ∂V_j = R_j lying in K₀. `Cycle` holds R_j, so
  it is not that chain. A bars-only path with representatives left to the reference is a fair answer if V isn't at hand. State
  which it is.

Checks:
- Exact multiset equality with `RelativeFiltration.ComputeBarcode` across the R-suite's families, every degree that isn't truncated.
- Each table row reached, with counts.
- One timing against the reference at a few sizes, with the load conditions stated.

Its gain is unproven, since nothing has measured the relative reduction as slow. Its other value is both readings from one run.
Report the measured ratio plainly. If it isn't faster, say so and let the user rule on whether it stays on the public surface.

### 4. `FastZigzag` H₁ representatives

`FastZigzag.Compute` (`src/tda/ph/FastZigzag.cs:154`) restricts the coned cycle to original cells for Add-born bars, gives null for
Delete-born ones, and its docs mark p ≥ 1 unverified: a cycle might include an incarnation of a re-entering cell deleted before the birth
step (`docs/TDA/ZigzagPersistence.md:141`). This is the real mathematics in this chip.

- **Contract first.** Say what `Bar.Cycle` means for a step bar, including closed and open birth ends, and put it in the docs.
- **Delete-born bars.** They carry a cone dimension shift. Derive what the coned cycle gives, or read it: the library has
  `graveyard/bibliotecha/bibliotecha/compendia/ph/references/2412.02591v2.md`, which mentions zigzag representatives (not yet
  read). Dey, Hou and Morozov's *A Fast Algorithm for Computing Zigzag Representatives* (SODA 2025) is the likely source. That
  citation is from memory, and the paper is not in the library as far as I could find.
- **The oracle stays outside the engine**, exact ℤ₂ linear algebra over the complex at the birth step:
  - The cycle lies in that complex and has zero boundary.
  - It is nontrivial there.
  - Add-born: its class is outside the image of the previous step's homology.
  - Delete-born: its class is in the kernel of the map into the previous step's homology.
  - These are necessary conditions. Say whether a stronger one is available, such as a representative sequence along the bar.
- **Re-entry** is the case to aim at, since it is what the unverified claim is about. Reach guards count such bars.
- **Through the frames.** A `Flag` frame bar that was filled names its cycle via `StepBar`, in frame terms, on a hand-built case and
  on the random churn.
- Update `docs/TDA/ZigzagPersistence.md` (lines 138-141, 349, 392).

If a representative turns out to be unavailable for some class of bar, that is a finding: document which, and don't paper over it.

### 5. The relative reading over frames

H(K_t, K₀) over a frame sequence, with K₀ in every frame (refuse otherwise). Union and intersection frames both contain K₀, so
relative Mayer-Vietoris holds and the Diamond Principle should carry over. That is my reasoning, so test it.

Route is the chip's call. Dropping K₀'s cells leaves cells with missing faces, which the zigzag engines may assume away. One option:
H_p(K, K₀) ≅ H̃_p(K ∪ CK₀), so coning K₀ to a new apex lets the unmodified engines read it, with one essential H₀ bar dropped.
Unverified, an idea only. Compare the extra cells against the quotient.

Checks:
- Monotone frames equal `RelativeFiltration.ComputeBarcode`, the bridge between the two readings.
- `FastZigzag`, `ZigzagBarcodeNaive` and `ZigzagMapBarcode.ComputeFromZigzag` agree at step level.
- `Union` and `Intersection` give the same frame bars and events.
- K₀ with its own cycle: those cycles vanish. K₀ a spanning tree: relative H₁ frame bars equal the absolute ones.

The API goes with `GraphFrameFiltration` and `ConditionedFiltration.BuildSweep`. Take the obvious shape, and bring it to the user only
if it isn't obvious.

### 6. Close

Remove the Shelf entries as they close, update the `zigzag-frontier.md` row, append a dated report to this brief (what was built,
what was measured, what stayed open and why), and delete the build folder.

## How to work

- **Each step ends with a fresh-eyes review**, correctness and coverage gaps, fixed before reporting. Fix each step's checks before its
  first run, and break the code on purpose to confirm they catch it.
- **New writing takes the module's symbols.** Log splits in `issues/notation-ledger.md`, rename nothing.
- **Tests and the one timing in step 3 are in scope.** Other experiments wait for the user to ask.
- **Decisions beyond these four items are the user's**, and the user may talk to you directly. In open discussion, talk it through
  without forced-choice cards. Pressure-test ideas, including this brief's.
- **Keep to `src/tda/ph`, `tests/TDA.Ph.Tests`, `docs/TDA`, `issues/ph` and `issues/SHELF.md`.** `src/maths`, `issues/regression` and
  `issues/spc-BARS/stopping-theory-primer.md` belong to other sessions. Check `git status` and `git log` before touching
  `Clustering.Dendrograms` or shared SPC files.
- **No sub-agents** without the user's permission.

## Repo practice

- Work directly on `main`. Commit your own files only, by explicit path (`git commit -- paths`), and leave nothing staged: the index
  is shared.
- Build and test with `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\ph-loose-ends\` set inline on each command, and delete
  that folder when done. VizCore tests run in Release.
- The machine is shared: `-m:2`, `-- xUnit.MaxParallelThreads=4`, heavy runs at below-normal priority. Run anything over about two
  minutes in the background.
- `docs/TDA` is CRLF, so splice it with Python scripts kept in your scratchpad. C#, briefs and Lean are LF. Check `git diff --stat` after.
- The Grep tool is denied by a hook here, so use `rg` through Bash. `python` is a `.bat` shim, so put scripts in files.

## Report
