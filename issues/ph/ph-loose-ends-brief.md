# Brief: closing the PH setup's loose ends

**Date:** 2026-11-02. **Status:** built; the [report](#report) is below. The user ruled `ComputeBothReadings` internal. Three items left over from it landed the same day, at the report's end.
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

### 2026-11-02

**What landed.** All four Shelf entries are closed. No consumer is wired. On `main`:

| Commit | What |
|---|---|
| `039f649` | The H₀ elder rule fixed, with a differential against `PersistentHomology` |
| `98ae2e3` | `RelativeFiltration.ComputeBothReadings`, the proved table as a code path |
| `37bc6a6` | `FastZigzag` names every bar's class where it begins; the contract in the docs; `ZigzagModuleOracle`; `GraphFrameFiltration.VerticesOf` |
| `8119c7e` | `GraphFrameZigzag.BuildRelative`, `ConditionedFiltration.BuildRelativeSweep`, `RelativeCycle` |
| `464f549` | Two other non-Rips deviations reported in the docs, not fixed |
| `a4a0d44` | The timing as measured, in the docstring and docs; the differential extended below zero |

**Left uncommitted in MarkBrain, as asked:** `persistent-homology/backbone-conditioned-persistence.md` and
`persistent-homology/opus-brief-conditioned-persistence-synthesis.md` (dated notes beside each stale claim, original text kept),
and `zigzag-engine/zigzag-frontier.md` (the "Frame sequences" row already uncommitted there, now updated, and the
"Representatives" row). All under `MarkBrain/ThermoMapper/issues/tda-purification/`.

**Suites.** `TDA.Ph.Tests` 433/433, `TDA.Mapper.Tests` 28/28.

**One thing that went wrong in the process.** Before committing step 4 I ran the full suite with `--no-build`, which read the
binary from the last deliberate break. It showed five failures, all that break's, and my command committed anyway because the
commit wasn't gated on the result. A rebuild straight after was 418/418 on the committed code. The step 5 commit was gated.

#### Step 1: the older notes

A dated note beside each stale claim, with the original text kept. "Absolute vs relative = two questions" stays as the user's framing, and
the note adds that with K₀ first they are one reduction read two ways. Whether Li's loop is VSMJ2011, which the synthesis brief also says,
I did not recheck, and its note says so.

#### Step 2: the elder rule

`ComputeH0` kept the larger birth. Flipped both comparisons, so ties go to the earlier cell.

| Check | Result |
|---|---|
| Random non-Rips filtrations through tetrahedra (400 trials), a quarter below zero; `maxDimension` all, 0, 1, 2 | Both engines equal `PersistentHomology`, every degree |
| Reach | The inverted rule differs on 400/400 trials. Finite H₀/H₁/H₂ 1759/884/21; essential 485/761/17 |
| All vertices at 0 (200 trials) | Bars equal the reference, and the old rule and the new give the same H₀, so no Rips result moved |
| Tests that pinned the inversion | None |
| Involuted cycles, non-negative values (200 trials) | 455 finite and 366 essential close and lie in the birth complex |

Watched failing on the old code first: the two hand cases and the differential, with the zero-birth check passing.

The differential turned up no further deviation in bars, in any degree. Reading the code and probing turned up two others with
different causes, reported and not fixed, and now in `docs/TDA/PersistentHomology.md`:

- **`PersistentInvolutedHomology` with `representatives` hangs on an essential H₁ bar when an edge value is negative.**
  `CycleReconstruction` runs a shortest path weighted by filtration value, and a negative undirected edge loops it. Three small
  cases hung past a 4 s limit. Rips values are distances, so nothing reaches this today, but a K₀-first run through that path
  would have, which is one reason step 3 goes through cohomology.
- **`PersistentCohomology` throws on a `CellFiltration`,** though it takes any `IFiltration`: its H₀ path reads `GetVertices`,
  which that type leaves empty.

#### Step 3: the relative barcode from one run

`RelativeFiltration.ComputeBothReadings(maxDimension)` returns `(Absolute, Relative)` from one run of the clearing engine.

- **Strictly first.** An `IFiltration` view of the ambient with K₀'s cells first, reading just below every ambient value. R12's
  order, without rebuilding a `SimplicialFiltration`.
- **Row dispatch.** Cell identity, not values. `PersistentCohomology.ComputePairs` returns each pair's birth and death cell, and
  a cell is in K₀ or not. For that to be the filtration's own pairing, `ComputeH0` now walks edges in filtration order instead of
  an unstable sort by value. That is a behaviour change at exact ties: which tied edge a bar names can differ on lists over 16
  edges, so an LMP protect-set could move at tied distances. The bars do not, and no test did.
- **Representatives.** Bars and generators only, so this is the bars-only answer. `Generator` matches the reference, including a
  moved-up bar's killer cell. `Cycle` is null: a moved-up class's relative cycle is V_j, and no engine keeps V. The finite kept
  bars could have one from an involuted run (R_j without K₀'s rows), but the reference has none for essential bars either, so
  `ComputeBarcode` stays the place to ask.

Checks are stronger than the brief's: equal pairs with generators, zero-length bars included, at `maxDimension` all, 0, 1, 2.
The one miss on the first run was a reach guard, with every equality holding: the any-filtration family never moved a class from
H₂ to H₃, since K₀ rarely held a closed 2-sphere. A denser third of trials fixed it, with no expectation changed.

| Family | Kept finite | Kept essential | Moved up H₀/H₁/H₂ | Dropped finite | Dropped essential |
|---|---|---|---|---|---|
| R7 and R8 backbones (200 trials) | 829 | 340 | 7/65/0 | 1253 | 248 |
| Any filtration, any subcomplex, through tetrahedra (300 trials) | 2183 | 1385 | 261/23/541 | 2134 | 341 |

In the second family, 207 moved-up classes had cells above them that could have killed them, the claim R12 couldn't reach, and
K₀ sat above the base in 289 trials. `Absolute` equals the ambient's barcode with K₀ at the base, built by hand. Plus the hand
cases R1–R6 and R10, and K₀ a 2-sphere that fills, relative H₃ [d, ∞).

**The timing, corrected.** My first write-up said to read CPU time because it holds up under load, and that "the one run is about 3 times
slower". The user pointed out that a busy machine probably confounds this, and it limits what these figures can say. CPU time drops the
waiting for a core, not contention for cache, memory bandwidth and hyperthreads. The machine stayed at 80–100% across 18 logical processors
through every run, with other sessions' test runs and long jobs going, so there is still **no quiet-machine timing**. I added allocated
bytes, which barely depend on load and in these reducers roughly track the toggles of column entries. They measure work and not time.

All degrees, in Release, at below-normal priority. K₀ is a spanning tree with a few extra edges, and content is a uniform draw. n=800 was
in a first version of the harness and ran 30 minutes under this load without finishing, so I dropped it.

| n | Edges | Triangles | MB: relative, `ComputeBarcode` | MB: one run | MB: absolute, column reducer | MB: absolute, clearing | CPU s, best of 3: relative | CPU s, best of 3: one run |
|---|---|---|---|---|---|---|---|---|
| 100 | 1567 | 5155 | 12.6 | 5.5 | 25.0 | 3.5 | 0.016 | 0.016 |
| 200 | 4300 | 13314 | 76.6 | 39.5 | 304 | 36.9 | 0.453 | 0.531 |
| 400 | 12626 | 42302 | 1203 | 1646 | 8392 | 1632 | 9.55 | 28.5 |

What survives the load:

- **At n=400 the relative barcode alone costs less through `ComputeBarcode`.** The one run allocates 1.37 times as much, and took 2.8–3.0 times the
  CPU on all three passes. The sign holds; **the size of the gap isn't established.** The two figures differ by a factor of 2. The one run costs
  about 2.3 times as much CPU per byte allocated, at both n=200 and n=400. That could be larger columns and deeper trees, and it could be cache
  misses, which are what contention inflates. Only a quiet run separates them.
- **The ratio grows with size.** In allocations it goes 0.44, 0.52, 1.37.
- **K₀ first is not what costs.** The one run and the clearing engine on the ambient allocate the same, 1646 and 1632 MB at n=400.
- **The relative complex is a smaller problem.** Quotienting K₀ cuts the column reducer's allocations 7 times at n=400 (time: 11 times).
- **For both readings the one run is the cheapest route measured.** At n=400 it allocates 1.7 times less than `ComputeBarcode` plus the clearing
  engine and 5.8 times less than `ComputeBarcode` plus the column reducer (CPU: 1.3 and 4.2 times).

What doesn't survive:

- **"Faster past small complexes."** At n≤200 the sign is open. Allocations favour the one run by about 2 times, and single-pass CPU ratios at n=200
  ran from 0.38 to 1.29, which is noise on this machine. n=100's CPU figures are at the clock's 16 ms tick besides.

This is one K₀ shape, small against the content. A K₀ holding more of the complex could move all of these.

Breaks, run once and reverted: K₀ not first failed 2 tests, the moved-up degree off by one 4, the elder tie-break flipped 3.

#### Step 4: `FastZigzag` representatives

**The contract,** now in `docs/TDA/ZigzagPersistence.md`: `Cycle` is a ℤ₂ cycle of the complex at the bar's first level whose
class there is the bar's own. It names the class where the bar begins, not a chain that lasts.

**Bars born by a deletion turned out to be the easy case.** The kernel of the backward map is spanned by [∂τ], so the class is
forced and the cone's dimension shift never has to be read. `Cycle` is ∂τ. I derived this, then found it stated: Definition 5 of
Dey and Hou's graph paper, `2103.07353v1.md` in the library ("the non-zero element in ker").

**The unverified claim is settled, in the direction the docs hoped, and it's a proof.** A closed birth's reduced column holds only
incarnations alive at the birth step, in every degree (argument in the docs): any column summed in has a pivot after σ, so it was
added later or its cell is deleted later. Re-entry raises no question, since incarnations are separate cells in Ê. The engine's
closed-birth path needed no change beyond the docs.

**The oracle** is `tests/TDA.Ph.Tests/ZigzagModuleOracle.cs`, exact ℤ₂ linear algebra with no engine code: the brief's necessary
conditions, and a stronger one. It answers to the barcode first, reading the same bars as `ZigzagBarcodeNaive`.

**A stronger condition is available, and it is built.** The right filtration at a level (Carlsson and de Silva, mirrored) says where
each class ends. A class labelled d has a compatible sequence to level d that ends there as Definition 5 asks, so this is the "set
of representatives" condition without building the sequence. It has teeth: of 1389 perturbed cycles passing every birth condition,
it rejected 447. The labels are my reading of their construction, mirrored, and the oracle's agreement with the naive engine is what checks it.

| Random churn with re-entry (120 trials) | H₀ | H₁ | H₂ |
|---|---|---|---|
| Add-born, all passing | 1360 | 1162 | 446 |
| Delete-born, all passing | 762 | 667 | 353 |
| Add-born through a re-entered cell | 798 | 920 | 378 |
| Delete-born through a re-entered cell | 526 | 543 | 307 |

Breaks: the wrong deleted cell failed 7 tests, dropping the pivot 8, keeping ω 4, and a break in the oracle's own walk 5.

**Through the frames.** `VerticesOf` reads a cycle in frame terms. By hand: C4 filled names C4's four edges, under both transitions.
On churn (300 sequences): 55 filled, 129 broken and 53 opened H₁ frame bars, all representatives. Under `Graph`, a class is its one
chain, and the 1407 H₁ bars' cycles are in every frame they cover.

**A finding, as the brief asked, not papered over.** Under `Flag`, the chain of a hole that opens is never in its first frame: 53 of
53, and it has to be, since a triangle leaves a flag complex only when an edge of it does. The chain is right where the bar begins,
and homologous there to the cycle at the frame, but a consumer asking for the hole at that frame doesn't get it from `Cycle`. Pinned
by a test. Bars born by an arrival stayed in their first frame in all 184 cases, which is a count and not a theorem: under `Union`
a chain could hold an edge a leaving triangle covers.

**Not built here: compatible cycles along a bar** (built later the same day, at the end of this report). That is what closes the finding above. Dey, Hou and Morozov's arXiv:2410.20565 exists
as cited; I read the abstract only, and it isn't in the library (`2412.02591v2.md` only cites it, as [26]). An idea, derived and
**not verified**: the coned run's V may give the sequence directly, with z_ℓ = R + Σ ∂ρ over the Add columns of V added before ℓ for
a closed birth, and z_ℓ = Σ ∂τ over the cone columns of cells deleted before ℓ for an open one. V isn't kept, and testing this is an
experiment, so it waits for the user.

#### Step 5: the relative reading over frames

`GraphFrameZigzag.BuildRelative(n, frames, subcomplex, transition, complex)` and `ConditionedFiltration.BuildRelativeSweep`, which
refuse a frame without K₀. `IsRelative` and `RelativeCycle(bar)` read the result.

**The route is the cone,** the brief's option. A quotient isn't available as things stand: a chord between two K₀ vertices has an empty relative
boundary and every engine infers dimension from the boundary, so it would read as a vertex. The apex enters first and leaves last, so
its component is the one bar born at step 0, and the views drop it.

**The Diamond reasoning holds, more simply than the brief's.** The union and intersection of two coned frames are the coned union and
intersection, so it is the absolute theorem applied to the coned complexes, and relative Mayer-Vietoris isn't needed. Tested below.

**Cost against the quotient.** |K₀|+1 extra cells where a quotient would save |K₀|−1, so 4|K₀| columns more in `FastZigzag`, once across
the sequence. On a 5-cycle backbone and two chords: 46 steps against 24 absolute, and a quotient would be 6. Not timed.

| Check | Result |
|---|---|
| Monotone frames equal `RelativeFiltration` (100 trials, both complexes, both transitions) | Equal, degrees through the top: 532 bars, 100 ending before the last frame, 89 H₂ under `Flag` |
| Against the relative chain complex's own decomposition, with no cone, through the oracle (300 trials) | Frame bars, events and transients equal, H₀–H₂: 1384 bars, H₁ 963 with 19 filled and 32 opened, H₂ 383, 15 transients, 104 K₀ triangles coned |
| Each frame bar's `RelativeCycle` | A representative in the relative complex, 1384 of 1384 |
| `FastZigzag`, `ZigzagBarcodeNaive`, `ZigzagMapBarcode.ComputeFromZigzag`, step level (10 trials) | Agree on 157 bars |
| `Union` and `Intersection` (200 trials) | Same frame bars and events on 670 bars; 30 transients one degree apart |
| K₀ with its own cycle, by hand | The cycle vanishes, a chord stays; filled, it moves up to H₂ |
| K₀ a spanning tree (100 trials) | Relative H₁ equals absolute on 425 bars, with events; no H₀ |
| K₀ empty (100 trials) | Equals the absolute reading |
| Frames without K₀; the graph engines on a cone over an edge | Refused; throw |

Everything passed on the first run. Reach was thin in churn at first (7 filled, 5 opened), so I raised the trials before the figures
above. Breaks: the apex's bar kept failed 8 tests, no cone over edges 10, no cone over K₀'s triangles 2, `RelativeCycle` keeping K₀ 2.

The graph engines can't read this, since the cone over an edge is a triangle, so a relative graph sweep runs through `FastZigzag`.

#### What stayed open, and why

- **`ComputeBothReadings` is internal, by the user's ruling (2026-11-02).** I recommended keeping it public and the user chose internal, at a point
  when my figure was still "about 3 times slower". Asked again on the softened figures above, the user kept it internal. `ComputeBarcode` is the one
  public path. The reading and its tests stay, with nothing outside the assembly calling it.
- **A timing on a quiet machine.** Not available today. The harness is ready; the command is on the Shelf.
- **Compatible cycles along a bar.** Above. The user's call; built later the same day, below.
- **The two reported deviations.** On the Shelf, per the brief; both fixed later the same day, below.
- **The tie change in `ComputeH0`.** Worth knowing if a tied-distance LMP run moves. I didn't measure one, since that's an experiment.

### 2026-11-02, later: three of the items left over

The user asked for three of the four Shelf items. The quiet-machine timing stays there. `TDA.Ph.Tests` 459/459, `TDA.Mapper.Tests` 28/28.

**`PersistentCohomology` reads any cell complex (`5af778d`).** The Shelf entry was the smaller half of this. The crash came from reading `GetVertices`,
but two more limits sat behind it, and they would have lost bars with no error: H₀ assumed every 1-cell has two faces, and the walk reached a degree only
through cofaces, so a cell with no faces was never read (the minimal torus's 2-cell, a relative cell inside K₀) and an empty degree ended the walk.

- A 1-cell now reads through its faces. Two join; one joins its vertex to a root older than every cell, since its far end lies in a quotiented
  subcomplex; none closes a cycle. A general chain complex's 1-cell with more than two faces has no union-find reading, so that case reduces the columns
  of ∂₁. Each degree's cells are found by dimension.
- Checked against the column reducer at every cutoff degree: 400 random filtered chain complexes (2,096 cells above degree 0 with no faces; 409 1-cells with
  one face, 347 with over two; finite H₀–H₂ 962/165/507, essential H₀–H₃ 1685/220/664/896), and 300 relative quotients (335 empty boundaries, 6,020 bars).
- **This reverses something above.** I wrote that the clearing engine reaches the relative barcode only through `ComputeBothReadings`. It now reads the
  quotient directly, which is a third route for the relative barcode alone. I haven't measured it; the timing's Shelf entry says to add it.

**Involuted representatives (`0cc8224`).** Negative values read by edge count, since a value is not then a length and an undirected negative edge has no
shortest path. That ends the hang.

- **The sweep found an older bug the Shelf didn't know about.** The path could run through an edge tied with the birth edge and after it, and with a tied
  triangle that makes the "cycle" a boundary, a class that is zero. All four value ranges failed on it when I first asked whether the cycle's class was
  the bar's own; a test that only checked the cycle closes, like the one above, can't see this. The path now runs only through edges before the birth edge.
- **Behaviour change at exact ties.** A Rips filtration with tied distances can now name a different cycle, the right one, so an LMP protect set
  (`H1CycleEdges`) at tied distances could shift. No bars move. I didn't measure a run, since that's an experiment.
- Lengths whose sum overflows reported "no path"; they now find it. Swept at scale 5e307, at shift −1e308, −20 and −1.5.

**Compatible cycles (`70fe9be`).** The idea at step 4 holds, and passed on its first run. `FastZigzag.ComputeCycles` returns `ZigzagCycles`, with
`CycleAt(bar, level)`; `GraphFrameFiltration` adds `CycleAtFrame`, `RelativeCycleAtFrame`, `LevelOfFrame` and `LevelBetweenFrames`.

- **The formula.** With the bar's death column V = ω·S + W: born by an addition, z_ℓ = S + ∂(W's cells added at step ≥ ℓ); born by a deletion, z_ℓ = ∂(S's cells
  deleted at step < ℓ). A bar keeps one chain, and a cycle is a view of it. `docs/TDA/ZigzagPersistence.md` has the argument: each z_ℓ lies in K_ℓ, consecutive
  ones are homologous in the larger complex, the ends are as the bar says, and induction on the level gives a basis at each level. It takes from Dey and
  Hou only that their pairs are the right bars. That argument is in prose; there is no Lean proof.
- **Checked at every level of every bar,** against the oracle, which gained `CheckSequence` and `CheckBasis`. On churn with re-entry (120 sequences, 23,360 levels):
  closed-closed/closed-open/open-closed/open-open in H₀ 655/705/421/341, H₁ 495/667/394/273, H₂ 197/249/160/193; bars whose chain moves 212/171/38; 3,560 through
  a re-entered cell. On frames (200 sequences, both complexes, both transitions, 51,173 levels), and on relative frames (423 bars, 1,157 frame cycles).
- **The check can refuse.** Of 204 bars that move, 112 were refused with their first cycle held still and 154 with the right cycles one level off; the rest stay
  valid up to homology. All 753 sequences with a level swapped for its sum with another bar alive at that level and the next were refused. A bar alive at that
  level alone would leave a valid sequence, since representatives aren't unique, so the test leaves those out.
- **The step-4 finding closes.** A hole that opens under `Flag`: 34 of 34 birth chains miss the first frame, as before, and every `CycleAtFrame` lies in its frame.
- **Cost, in allocated bytes** (Release, random flag frames, n = 20, 40, 80, at below-normal priority on a busy machine):

  | Steps | Bars | MB: bars | MB: `representatives` | MB: `ComputeCycles` | Ratio to `representatives` | Reading |
  |---|---|---|---|---|---|---|
  | 856 | 253 | 1.6 | 1.8 | 2.7 | 1.55 | 2,220 cycles, 0.5 MB |
  | 5,120 | 961 | 23.4 | 24.4 | 43.3 | 1.78 | 8,656 cycles, 1.8 MB |
  | 37,528 | 3,883 | 612.5 | 619.4 | 1,194.7 | 1.93 | 34,465 cycles, 7.0 MB |

  n = 160 did not finish inside the run's ten-minute limit, so it has no figure. The ratio rises with size, so 1.93 is not a ceiling. Mean cycle size there is 3, since a flag complex's cycles are short; other complexes will read larger. Set
  `ZIGZAG_CYCLES_COST` to a file and run `ComputeCycles_Cost` in Release to redo it.
- **Against Dey, Hou and Morozov.** I read their introduction and definitions, not their algorithm. They compute the same object with wires and bundles in O(m²n). This costs
  what the coned reduction costs, O(m³ log m) worst case with sorted-set columns, and I have not compared the two on any input. Their paper has no reading from
  a coned run in the part I read; I don't know whether one appears elsewhere, so the docs call the formula derived here and claim nothing beyond that.

**Stayed open.**

- **Essential bars' cycles.** `PersistentHomology` names none, `RelativeFiltration.ComputeBarcode` inherits that, and the involuted engine gives H₂ and up the
  birth simplex alone, which isn't a cycle. An unpaired column's V is its cycle, and the reducer now keeps V for paired ones. On the Shelf.
- **A Lean proof of the formula.** Not asked for, and not done. It would suit the harness.
- **MarkBrain's commits are the user's.** One more uncommitted edit there: `zigzag-engine/zigzag-frontier.md`, the Representatives row.
