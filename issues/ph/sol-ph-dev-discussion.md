# So im trying to sift through old markdown files and try to…

> So im trying to sift through old markdown files and try to isolate unfinished development in topology/persistent homology by reconciling discussions against current source code and doc strings. these are two directories that are close to recent dev work `D:\aghado01\MarkBrain\ThermoMapper\issues\tda-purification\persistent-homology` `D:\aghado01\MarkBrain\ThermoMapper\issues\tda-purification\zigzag-engine`

<details><summary>9 previous messages</summary>

> I'll start by mapping out what's in these two discussion directories, then reconcile against current source. Let me explore the structure first.
>
> I see the scope — the user pointed me at the two most active subdirs (`persistent-homology` with 18 files, `zigzag-engine` with 10). Per project memory, `zigzag-frontier.md` is the canonical cross-track status and the per-rung briefs _lag_ committed code, so I should anchor on the status/plan docs first, then verify against source.
>
> Let me gauge file sizes and modification times before reading, to spend tokens wisely.
>
> All modified in the same 2026-06-23/24 window (likely a vault sync — mod time won't discriminate). Note `perplexity-sifts-scaffolded 1.md` is a byte-identical duplicate of `perplexity-sifts-scaffolded.md`.
>
> Let me read the canonical status/plan trackers first — these enumerate planned vs done, and per memory the per-rung briefs lag actual code.
>
> This gives me a clear planned-vs-claimed picture. Now the critical step per the frontier doc's own rule — _"if a row and the code disagree, the code wins"_ — is to verify these against the **current source**, since `phased-plan.md` is dated 2026-06-19 and may be stale. Let me locate the source tree and confirm which "unfinished" items are genuinely still missing vs. already landed.
>
> Already a key reconciliation signal: `CombinatorialIndex.cs`, `IFiltration.cs`, and `LazyRipsFiltration.cs` **exist** — so Phase C's C1/C2/C3 (which `phased-plan.md` lists as not-done in "Track β") have in fact landed since the doc's 2026-06-19 date. Notably, there's **no `tda/pipelines/sifts`** and no Mapper-over-T consumer in the tree. Let me run a single content sweep to pin down each candidate precisely.
>
> Good signal. Several "Track β"/Phase C items have clearly landed (C1/C2/C3/C4), while genuine gaps emerge: **SIFTS pipeline**, **Z/p coefficients**, **alpha/cubical/merge-tree/circular-coord features**, **Z6 quasi-zigzag**, and **zigzag moving-PCD consumer** are all ABSENT. A few hits need disambiguation (some are prose false-positives). Let me verify the ambiguous ones directly: the 🟡 zigzag in-progress items (ReflectionZigzag backward arrows, RuVineyard), `PersistentMapper`, `DiagramMetrics.Bottleneck`, apparent-pairs, and what `StrongCollapse` actually is.
>
> Sharp findings. Three confirmations: **`DiagramMetrics.Bottleneck` is a `throw new NotImplementedException(...)` stub** (D2 genuinely unfinished); **apparent-pairs (C5) is absent** from the reducer (the earlier hit was prose); and **ReflectionZigzag** has advanced past the 2026-06-23 frontier snapshot — stage 2b-i (reverse teardown) landed, but stage 2b-ii (arbitrary removal via vineyard transpose-to-last) still throws. Also **`PersistentMapper` landed** but builds a _monotone_ `NerveFiltration` over T, not the _zigzag_ filtration the frontier's consumer item called for.
>
> Let me close the last gaps: RuVineyard's completeness, what `StrongCollapse` actually is (E2 edge-collapse vs strong-collapse), confirm no SIFTS pipeline dir, and whether SPRED's SA path consumes the PH engine.
>
> Clear now. **RuVineyard** is a complete oracle-grade primitive (`VineSwap` validated in isolation, no stubs) but is _not yet wired_ into ReflectionZigzag for arbitrary removal — so stage 2b-ii remains the open seam. **StrongCollapse** is genuinely Boissonnat–Pritam _strong_ collapse (arXiv:1809.10945, vertex domination) — note its docstring flags a "fast nerve-based algorithm" and "§4 persistence version" as follow-ons. **SIFTS pipeline is confirmed absent** (no `tda/pipelines`), and **SPRED does not reference the PH engine at all** (D1 wiring unmet). Let me make one last check on whether StrongCollapse's fast/persistence variants and the distinct flag-complex _edge_ collapse exist.
>
> Strong resolution. **StrongCollapse**: only the slow ground-truth `Core` is implemented; the §3 Algorithm-1 _fast_ nerve-based collapse is still a docstring-flagged follow-on. But **StrongCollapseZigzag** does implement the §4 core-assembly over a `ZigzagFiltration` — so the frontier's "✅ verify" scope question resolves _yes, it covers the zigzag path_. The flag-complex **edge** collapse (the different Boissonnat–Pritam paper, 1810.11602) that Phase E2 wanted is not present.
>
> I've anchored on the two authoritative status docs, which index the construction briefs. Before I synthesize, let me make sure no _live_ unbuilt thread hides in the idea-docs I haven't opened — especially `void-inference-three-engines.md` (the freshest file, 2026-06-23, largest). I'll extract just headers + status lines to stay economical.

</details>

I now have a source-verified reconciliation. The idea-docs surfaced two more live threads the phased-plan spine doesn't track (the "three engines" Engine 3, and a construction-layer reorg). Here's the consolidated picture.

## Method

I anchored on the two authoritative status docs — `zigzag-frontier.md` (which by its own rule says _"if a row and the code disagree, the code wins"_) and `phased-plan.md` (dated 2026-06-19, so suspect) — then verified every "unfinished" candidate against the actual tree in [src/tda/ph](src/tda/ph), [src/tda/primitives](src/tda/primitives), and [src/maths/topology](src/maths/topology). The 2026-06-19 plan is materially stale, so several "next" items are already done and a few "done/landed" rows hide real stubs.

## Genuinely unfinished — discussed, source-confirmed open

| #   | Item                                                                                   | Source state (verified)                                                                                                                                                                                                                             | Discussion doc                                                                            |
| --- | -------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| 1   | **SIFTS pipeline** (text → Time-Skeleton → PH → Barcode)                               | **No `tda/pipelines` dir at all.** Fully unbuilt despite ~5 docs.                                                                                                                                                                                   | `opus-brief-sifts-integration.md`, `perplexity-sifts-scaffolded.md`, `Integrating SIFTS…` |
| 2   | **`DiagramMetrics.Bottleneck`** (D2)                                                   | **Stub:** `throw new NotImplementedException("Bottleneck … is P1")` at [DiagramMetrics.cs:75](src/tda/ph/DiagramMetrics.cs#L75). Gate still runs on Wasserstein.                                                                                    | `phased-plan.md` D2                                                                       |
| 3   | **Apparent pairs (C5)**                                                                | Absent from [FiltrationCohomologyReducer.cs](src/tda/ph/FiltrationCohomologyReducer.cs). (C4 emergent pairs _did_ land.)                                                                                                                            | `phased-plan.md` Phase C                                                                  |
| 4   | **SPRED ↔ pCoh wiring (D1)**                                                           | [Spred.cs](src/maths/geometry/dim-reduction/Spred.cs) has **zero** references to Barcode/Persistent/Topology. The SA cost is not on pCoh; Track β exit unmet.                                                                                       | `phased-plan.md` Track β / D1                                                             |
| 5   | **ReflectionZigzag stage 2b-ii** (3rd general oracle)                                  | Stage 2b-i (reverse teardown) landed, but arbitrary removal **throws** `NotSupportedException` at [ReflectionZigzag.cs:91](src/tda/ph/ReflectionZigzag.cs#L91). `RuVineyard` primitive is complete but **not wired in**.                            | `z3-reflection-zigzag-brief.md`                                                           |
| 6   | **StrongCollapse §3 fast path**                                                        | Only the slow ground-truth `Core` exists; the "nerve-based fast algorithm (§3 Algorithm 1)" is still a docstring-flagged follow-on at [StrongCollapse.cs:15](src/tda/ph/StrongCollapse.cs#L15).                                                     | engine docstring                                                                          |
| 7   | **Zigzag consumers / telos**                                                           | Moving-PCD, embedding/transformer-layer trajectories, and a true **zigzag** nerve-over-T builder all absent. `PersistentMapper` exists but is **monotone** (`NerveFiltration`, non-decreasing T) — not the non-monotone zigzag the frontier wanted. | `zigzag-frontier.md` "consumers", `plan.md` §consumers                                    |
| 8   | **Z6 quasi-zigzag** (multiparameter × zigzag)                                          | Not started — confirmed absent.                                                                                                                                                                                                                     | `zigzag-frontier.md`                                                                      |
| 9   | **"Three Engines" Engine 3** (GA / population-MCMC for void / minimal-cycle inference) | Concept synthesis only (§4.4 "buildable spec", §7 "open decisions before building"). GA papers staged as raw PDFs, _not even ingested_. Engines 1 (phaseless-skeleton/Hodge) and 2 (RJMCMC) also conceptual.                                        | `void-inference-three-engines.md` (freshest doc, 2026-06-22)                              |
| 10  | **Construction-layer reorg**                                                           | Proposed `tda/complexes` / `tda/nerve` re-seat, triangle-enum sink to `Graphs.Primitives`, "deferred ripples from the tda/ph move" still ⬜. Tree is still flat `src/tda/primitives`.                                                               | `rips-flag-nerve-brief.md`                                                                |
| —   | **Phase E design-for-later**                                                           | E1 Z/p, E2 flag-complex _edge_ collapse (1810.11602 — distinct from the strong collapse that _did_ land), E3 alpha, E4 cubical, E5 merge trees, E6 packed coeffs, E7 circular coords — **all absent.**                                              | `phased-plan.md` Phase E                                                                  |

## Docs say "next/unfinished" but the code already landed (don't rebuild)

The 2026-06-19 plan lists these under open "Track β" / 🟡, but they're in the tree now:

- **C1 combinatorial index** → [CombinatorialIndex.cs](src/maths/topology/CombinatorialIndex.cs)
- **C2 `IFiltration`** → [IFiltration.cs](src/tda/ph/IFiltration.cs) · **C3 lazy Rips** → [LazyRipsFiltration.cs](src/tda/primitives/LazyRipsFiltration.cs) · **C4 emergent pairs** → in the reducer
- **ReflectionZigzag** is past "stage 1 only / all backward arrows throw" — reverse-teardown works now.
- **Strong-collapse zigzag scope** (frontier's "✅ verify"): confirmed — [StrongCollapseZigzag.cs](src/tda/ph/StrongCollapseZigzag.cs) implements the §4 core-assembly over a `ZigzagFiltration`.

## Doc hygiene / discrepancies

- **Duplicate**: `perplexity-sifts-scaffolded 1.md` is byte-identical to `perplexity-sifts-scaffolded.md` (both 38.1 KB) — delete the ` 1` copy.
- **`phased-plan.md` is stale** — it lacks the table-thin "code wins" discipline that makes `zigzag-frontier.md` reliable. The PH side has no equivalent live dashboard.
- **Placement drift (real, worth a decision):** the plan's discriminator says _pure → `Maths.Topology`_, but pure types `DiagramMetrics`, `Barcode`, the reducers, and the zigzag engine all live in `TDA.Ph` under `src/tda/ph`. Only `Simplex`/`CombinatorialIndex`/`Z2LinearAlgebra` actually floored to `Maths.Topology`. Either the floor move stalled or the placement decision changed — the docs and tree disagree, and `rips-flag-nerve-brief.md` flags the same "needs a scrub."

The genuinely open, collision-free, high-leverage targets are **#1 SIFTS** (most-discussed, zero code), **#2 Bottleneck** (small, unblocks the gate's metric), and **#5 ReflectionZigzag 2b-ii** (the primitive is already sitting there unwired).

Want me to (a) write a thin `ph-frontier.md` status dashboard mirroring the zigzag one and flip the stale `phased-plan.md` rows, (b) delete the duplicate, and/or (c) spin off any of these (SIFTS, Bottleneck) as a scoped follow-up task?

> So im trying to sift through old markdown files and try to isolate unfinished development in topology/persistent homology by reconciling discussions against current source code and doc strings. these are two directories that are close to recent dev work `D:\aghado01\MarkBrain\ThermoMapper\issues\tda-purification\persistent-homology` `D:\aghado01\MarkBrain\ThermoMapper\issues\tda-purification\zigzag-engine`

<details><summary>4 previous messages</summary>

> The user re-sent the same task, now with Ultracode on. Per your CLAUDE.md note, I'm reading that as a request for **maximum primary-agent depth and completeness** — not a license to fan out into Workflow/subagents (your repo policy forbids that, and you flagged this exact conflict). So I'll do an exhaustive _sequential_ pass: digest **every** doc in both directories for unfinished-development signal, then verify each named component against source — no anchoring-and-skipping this time.
>
> Let me start by digesting all the docs I haven't yet read in full, extracting headers + every "unfinished/deferred/stub/open" line.
>
> The digest sharpens everything and surfaces precise specs for the open items. Key refinements: `z3-reflection-zigzag-brief.md` pins the exact remaining ReflectionZigzag work (stage 2b-ii arbitrary removal + stage 3 churn-sweep); `z5a` was a parity-fix task on `GraphZigzag` (was WIP 13/18); `ph-directive-addenda.md` confirms WRD2025/DBK2023 deferred consumers; `implementation_plan.md` + `opus-brief-sifts-integration.md` define the SIFTS structural reorg; and `void-inference-three-engines.md` is a three-engine concept doc with Engine 3 explicitly "the one to build."
>
> Let me read the high-signal unbuilt-thread docs in full and run a second verification sweep for the newly-named components.
>
> This is decisive. The SIFTS brief and `implementation_plan.md` reveal that the **`tda/primitives` → `tda/ph` reorg is the spine of the SIFTS workstream** — and the tree shows it's _half-executed_: `src/tda/ph/` is populated (namespace `TDA.Ph`), yet `src/tda/primitives/` still exists with `FlagComplex.cs`, a second `Barcode.cs`, `NerveDiff.cs`, `Filtration.cs`. Let me verify the exact partial-reorg state plus the remaining newly-named components (nerve relocation, ph→viz adapter, void-inference engines, complex-analytic Laplacian, GraphZigzag WIP status, SIFTS project).
>
> Now the picture is fully resolved, and the biggest finding is structural. Let me make one final precise check — the two `Barcode.cs` files (one in `TDA.Ph`, one in `TDA.Primitives`) — to characterize whether that's a genuine duplicate or the `plan.md` pure/graph-H0 split, before I synthesize.

</details>

Fully resolved — and the two `Barcode.cs` files are **not** a duplicate: `TDA.Ph` has the pure `Barcode` currency, `TDA.Primitives` has the graph-H0 `PersistenceBarcode.ComputeH0`. That's exactly the `plan.md` split, executed correctly. The real structural finding is bigger and sits underneath everything.

---

# TDA / PH unfinished-development reconciliation

_Method: anchored on the two authoritative status docs, then digested **all 28** discussion files for unfinished-signal lines and verified every named component against the current tree (`src/tda/ph`, `src/tda/primitives`, `src/maths/topology`, `src/graphs/spectral`). "Source verdict" = what the code actually shows, which overrides the docs where they disagree._

## The meta-finding: the "purification floor" is half-executed (and the authoritative doc says it's done)

`plan.md`, `phased-plan.md`, and `engine-port-plan.md` all decree: pure engine → `Maths.Topology`. `phased-plan.md` even states it is _"already executed."_ **The code contradicts this:**

- `Maths.Topology` (`src/maths/topology`) contains only `Simplex`, [CombinatorialIndex](src/maths/topology/CombinatorialIndex.cs), `Z2LinearAlgebra`.
- The actual pure engine — `PersistentHomology`, `PersistentCohomology`, `PersistentInvolutedHomology`, [Barcode](src/tda/ph/Barcode.cs), [DiagramMetrics](src/tda/ph/DiagramMetrics.cs), `SimplicialFiltration`, `IFiltration`, `FiltrationCohomologyReducer` — lives in **`src/tda/ph`, namespace `TDA.Ph`**. These are pure (no `CsrGraph`) yet were never floored.
- `src/tda/primitives` (project `TDA.Primitives`) **still exists** with 10 files, all still `namespace TDA.Primitives`: `FlagComplex`, `PersistenceBarcode`, `Filtration`(nerve)+`NerveDiff`, `RipsFiltration`, `LazyRipsFiltration`, the codim-one duals, `EmbeddedComplexZigzag`, `H1CycleEdges`.

So the code actually followed the **SIFTS brief's** `tda/ph`/`TDA.Ph` destination — which `phased-plan.md` explicitly calls _"stale, superseded by Maths.Topology."_ The doc and the tree disagree on the single most load-bearing architectural decision in this whole issue folder, and three `projects/` (`TDA.Ph`, `TDA.Primitives`, `TDA.Mapper`) coexist where the plan wanted a clean `TDA.PH` + `TDA.Pipelines.Sifts`. **This is the root unfinished item; most of the SIFTS workstream below is blocked on resolving it.**

## A · PH engine internals — genuinely open

| Item                                                                               | Discussed in                                                       | Source verdict                                                                                                    |
| ---------------------------------------------------------------------------------- | ------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------- |
| **C5 apparent pairs**                                                              | `phased-plan` C5, `ripserer-notes`, `engine-port-plan`             | **Absent** from [FiltrationCohomologyReducer](src/tda/ph/FiltrationCohomologyReducer.cs) (C4 emergent _did_ land) |
| **D2 `DiagramMetrics.Bottleneck`**                                                 | `phased-plan` D2, `ph-directives` Spec 1                           | **Stub** — `throw new NotImplementedException` at [DiagramMetrics.cs:75](src/tda/ph/DiagramMetrics.cs#L75)        |
| **D1 SPRED SA → pCoh wiring**                                                      | `phased-plan` D1, `ripserer-notes` §perf                           | [Spred.cs](src/maths/geometry/dim-reduction/Spred.cs) has **zero** topology refs — Track β exit unmet             |
| **R4 bitvector/XOR column rep**                                                    | `opus-brief-sifts` R4, `hpc-considerations`                        | Reducer still `SortedSet<int>`, not bit-packed — deferred, unbuilt                                                |
| **StrongCollapse §3 fast nerve algo**                                              | [StrongCollapse.cs:15](src/tda/ph/StrongCollapse.cs#L15) docstring | Only slow ground-truth `Core`; fast Algorithm-1 unbuilt                                                           |
| **E1 Z/p coefficients**                                                            | `engine-port-plan`, `phased-plan` E1                               | Absent (Z/2 only)                                                                                                 |
| **E2 flag-complex _edge_ collapse** (Boissonnat–Pritam **1810.11602**)             | `engine-port-plan`, `phased-plan` E2                               | Absent — only _strong_ collapse (1809.10945) landed; distinct paper                                               |
| **E3 alpha · E4 cubical · E5 merge trees · E6 packed coeffs · E7 circular coords** | `engine-port-plan` §94, `phased-plan` E                            | All absent (design-for-later)                                                                                     |

## B · Zigzag engine — genuinely open

| Item                                                                       | Discussed in       | Source verdict                                                                                                                                                                                          |
| -------------------------------------------------------------------------- | ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **ReflectionZigzag stage 2b-ii** (arbitrary removal via transpose-to-last) | `z3-brief` §4      | Throws `NotSupportedException` at [ReflectionZigzag.cs:91](src/tda/ph/ReflectionZigzag.cs#L91); [RuVineyard](src/tda/ph/RuVineyard.cs) primitive is **ready but unwired**. Live frontier (`gudhi port`) |
| **ReflectionZigzag stage 3** (fold into churn-sweep, "5th corner")         | `z3-brief` §5      | Not integrated                                                                                                                                                                                          |
| **Z6 quasi-zigzag** (multiparameter × zigzag)                              | `frontier`, `plan` | Absent — not started                                                                                                                                                                                    |
| `strictRepresentatives:true`                                               | `z5d` §1.4         | Throws `NotSupportedException` in [GraphZigzag](src/tda/ph/GraphZigzag.cs)/`GraphZigzagFast` — known omission                                                                                           |
| Worst-case-bound faithfulness gaps                                         | `z5d` §1.6         | Correctness intact; perf bounds not the paper's — documented, open                                                                                                                                      |

## C · Consumers / the telos — almost entirely unbuilt

| Item                                                                | Discussed in                                  | Source verdict                                                                                                                            |
| ------------------------------------------------------------------- | --------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| **Zigzag** nerve-over-T builder (non-monotone)                      | `frontier` consumers, `plan`                  | [PersistentMapper](src/tda/mapper/PersistentMapper.cs) is **monotone only** (requires non-decreasing T) — the _zigzag_ version is unbuilt |
| Neural moving-PCD / embedding-layer trajectories                    | `frontier` consumers                          | Absent                                                                                                                                    |
| **`ph→viz` adapter** (engine geometry → viz layers)                 | `engine-port-plan` §107                       | Absent                                                                                                                                    |
| **WRD2025 localized gate (P2)** — which vertices/loops drive change | `ph-directive-addenda`, `ph-engine-placement` | Absent — deferred consumer; R_eff gate (DBK2023) _did_ ship                                                                               |

## D · SIFTS workstream — the largest unbuilt application

| Item                                                                             | Discussed in                                  | Source verdict                                                                                         |
| -------------------------------------------------------------------------------- | --------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| Dissolve `tda/primitives` → `tda/ph` + `mapper/nerves`                           | `opus-brief-sifts` §49, `implementation_plan` | **Half done** — `tda/ph` stood up; `tda/primitives` (10 files, `TDA.Primitives`) never dissolved       |
| **R2 nerve relocation** → `tda/mapper/nerves` (`TDA.Mapper.Nerves`)              | `opus-brief-sifts` R2                         | **Absent** — `Filtration.cs`+`NerveDiff.cs` still in primitives                                        |
| **`TDA.Pipelines.Sifts` project + `SiftsPipeline`**                              | `opus-brief-sifts` §54, `implementation_plan` | **Absent** — no `tda/pipelines` dir, no project                                                        |
| **R1 Time Skeleton** (the "T" in SIFTS)                                          | `opus-brief-sifts` R1                         | Absent                                                                                                 |
| **SIFTS experimental design** (nursery-rhyme β₁ tests + Time-Skeleton unit test) | `opus-brief-sifts` R3, `implementation_plan`  | Absent (no Sifts test project)                                                                         |
| R3 canonical topological zoo (S¹/2-circles/T²/S²/fig-8, _known_ Betti)           | `opus-brief-sifts` R3                         | Partial — `phased-plan` A1 has H2/torus/2-sphere; the SIFTS-named zoo + nursery rhymes are not present |

## E · Research / exploratory tracks (concept-only)

| Track                                                                                                         | Discussed in                                           | Source verdict                                                                                                                                                                                              |
| ------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Void-inference "three engines"** (1 phaseless-sync/Hodge · 2 RJMCMC region-count · 3 GA/pop-MCMC min-cycle) | `void-inference-three-engines` (freshest, 2026-06-22)  | **Concept synthesis only.** The RJMCMC _substrate_ exists (`src/maths/samplers/rjmcmc/*`, BARS) but the void engines are unbuilt; Engine 3 = "the one to build", GA papers **not even ingested**            |
| **Backbone-conditioned persistence** (SIFTS reframe; neural manifolds on task-event backbones)                | `backbone-conditioned-persistence`                     | Design realization, "future research" — unbuilt                                                                                                                                                             |
| **Complex-analytic PH** (Hermitian/persistent Laplacian, cellular sheaves, conformal Mapper lenses)           | `Gemini-Complex-Analysis`, `gemini-magnetic-laplacian` | Research ideation. [MagneticLaplacianOperator](src/graphs/spectral/MagneticLaplacianOperator.cs) **landed** in graphs/spectral, but the _persistent/Hermitian-Laplacian_ + sheaf TDA constructs are unbuilt |

## Stale "done/next" — actually landed, don't rebuild

C1/C2/C3/C4 (`CombinatorialIndex`, `IFiltration`, `LazyRipsFiltration`, emergent pairs) · ReflectionZigzag stage 1 + **2a** (`RuVineyard` + vineyard transposition) + **2b-i** (reverse teardown) · z5a `GraphZigzag` now **"parity-complete (H0)"** (was WIP 13/18) · all of z5d codim-one (A1, A2 p=2/p≥3 `VoidBoundary`, item C/F, multigraph-safe H0, E-builder hardening) · `StrongCollapseZigzag` §4 core-assembly (resolves the frontier "✅ verify") · DBK2023 `EffectiveResistance` gate · monotone `PersistentMapper` · the `Barcode`/`PersistenceBarcode` split.

## Doc hygiene & contradictions

1. **`phased-plan.md` (2026-06-19) is materially stale** — wrongly claims the `Maths.Topology` floor is "executed", and lists landed C1–C4 as open. No PH-side equivalent of the well-maintained `zigzag-frontier.md` dashboard, which is _why_ it drifts.
2. **Unresolved placement contradiction** (the meta-finding) — docs say `Maths.Topology`; code says `TDA.Ph`. This needs _your_ decision: finish the floor, or ratify `TDA.Ph` and correct the docs.
3. `perplexity-sifts-scaffolded 1.md` is a byte-identical duplicate of `perplexity-sifts-scaffolded.md` — delete.
4. `ph-directives.md` is marked "retired (A2 ✅)" but its Spec 1 is the **live, unimplemented** Bottleneck spec — not fully retired.

**Highest-leverage open targets:** (1) resolve the `Maths.Topology` vs `TDA.Ph` placement and finish/ratify the floor (unblocks SIFTS) → (2) `Bottleneck` (small, declared, gate-relevant) → (3) ReflectionZigzag 2b-ii (primitive already sitting unwired) → (4) the SIFTS `SiftsPipeline` itself.

Want me to (a) write a thin, `zigzag-frontier`-style **`ph-frontier.md`** dashboard encoding this reconciliation and flip the stale `phased-plan.md` rows, (b) delete the duplicate file, and/or (c) spin off `Bottleneck` or the placement-floor reconciliation as scoped follow-up task chips? I won't touch docs or commit anything without your go-ahead.

> maths.topology and tda.ph are two sides of the same coin. this entire project is organized around things being reusable by posing things in general terms and building application layers bottom up. so topology primitive concepts live in maths, ph application specifics live in ph etc

Got it — that collapses my meta-finding, and I'm retracting it. `Maths.Topology` ↔ `TDA.Ph` isn't a half-executed floor; it's the bottom-up split working as designed. General topology/algebra primitives (`Simplex`, combinatorial index, Z/2 algebra) belong in `maths/topology`; PH-application specifics (the reducers, `Barcode`/`DiagramMetrics`, the abstract filtrations, the zigzag oracles) belong in `tda/ph`. My error was treating the plan's `CsrGraph` discriminator as the _whole_ rule — "pure / no `CsrGraph`" is necessary but not sufficient for flooring to maths; domain-generality is the real criterion, and persistence/barcodes/diagram-metrics are PH-domain, so `TDA.Ph` is exactly right.

By the same logic `tda/primitives` is consistent too — it's the graph-coupled construction/implementor layer (`FlagComplex`, `Rips`/`LazyRips`, nerve + diff, the dual-graph builders, `EmbeddedComplexZigzag`). So: general primitives in `maths`, the PH engine in `ph`, graph-coupled builders in `primitives`. The SIFTS brief's proposed regroup (FlagComplex→ph, nerves→mapper/nerves) is organizational polish, not missing function — I'm pulling it out of the "unfinished" list.

What this does **not** change — the genuine backlog is all _functionality_, and it stands as verified:

- **PH engine:** apparent pairs (C5), the `Bottleneck` stub (D2), SPRED↔pCoh wiring (D1), bit-packed columns (R4), StrongCollapse's fast nerve path, and the E1–E7 feature set.
- **Zigzag:** ReflectionZigzag stage 2b-ii + churn-sweep stage 3, Z6 quasi-zigzag, `strictRepresentatives`.
- **Consumers:** zigzag nerve-over-T (vs. today's monotone `PersistentMapper`), moving-PCD / embedding trajectories, the `ph→viz` adapter, WRD2025 localized gate.
- **SIFTS proper:** `SiftsPipeline` + Time Skeleton + experimental-design tests — the file moves were never the substance.
- **Research tracks:** void-inference three engines, backbone-conditioned persistence, Hermitian/persistent Laplacian (the magnetic-Laplacian operator already landed in `graphs/spectral`).

Doc hygiene shrinks accordingly: the byte-identical `perplexity-sifts-scaffolded 1.md` duplicate, and `phased-plan.md` being stale on Phase C (C1–C4 have landed). Its "floored to `Maths.Topology`" line is loose phrasing for the right design — a wording tighten, not a decision to make.

Want me to draft the `ph-frontier.md` dashboard (functionality-only now) and/or delete the duplicate?

> SIFTS turns out to be a rather simplistic case of something more general that i'm reaching for. so its not going to be 'sifts pipeline'. you still have some more reading to do before we attack that

<details><summary>2 previous messages</summary>

> You're right — I header-digested those but never actually read them. The generalization is exactly what I parked in the "research" bucket without reading: `backbone-conditioned-persistence.md` (literally subtitled "reframing the SIFTS workstream"), the scaffold-concept ideation in `perplexity-sifts-scaffolded.md`, and `void-inference-three-engines.md`. Let me read those three in full now.
>
> The picture is sharpening fast. `backbone-conditioned-persistence.md` is explicit: **SIFTS-proper is a "basic seed"; the real object is the backbone as a general primitive — _prior-conditioned persistence_**. The mechanism (build a filtration with a base subcomplex K₀) floors to the engine; the _choice_ of K₀ (reading-order, sweep-axis, spine, MST, task-events) is the application. And the load-bearing unification is **embeddings ≡ neurons** — token-position-on-embeddings and task-events-on-neurons as the _same_ backbone-conditioned-trajectory problem, with your Bichot 2019 PFC→V4 data as the flagship.
>
> The third doc is your own voice on the scaffold concept. Let me read it, plus the canonical backbone brief in the sibling `mapper-integrations/` folder.

</details>

Read them — and the canonical brief in `mapper-integrations/` (`opus-brief-backbone-conditioned-persistence.md`) is the keystone. I see the shape now, and why "SIFTS pipeline" was the wrong target. Let me reflect it back so we're aligned before attacking.

## What SIFTS is a special case of

The genus is **backbone-conditioned (prior-conditioned) persistence**: one primitive, `BackboneConditionedFiltration(B, X, metric)` — a _mandatory_ backbone `B` laid in at birth ε₀≈0, similarity edges born at distance, → `SimplicialFiltration` → the existing PH/zigzag engine. The backbone is the whole trick: **it converts a _similarity_ structure into a _return_ structure** — Rips alone says "what's alike"; add the backbone and the same chords become "where the sequence loops back." H₁ = a return. SIFTS is the degenerate instance: backbone = reading order, metric = cosine on TF-IDF (which your own thread correctly identifies as the two most-aged assumptions). It floors to `tda/ph` as the shared core; the _choice_ of backbone (reading-order = SIFTS, delay-graph = Li, task-events = COSYNE) stays in thin `pipelines/*` clients. That refines — doesn't contradict — the SIFTS brief's R1.

## Your proximal target (the thing we actually attack)

Not neural manifolds — that's the telos. The near object is the **two time-like axes** for NLP/embeddings and the **warped path between them**:

- **Axis 1 — positional backbone:** reading / token position (zero-cost scaffold edges).
- **Axis 2 — thermal:** the SPC temperature sweep, where the coupling kernel `J_ij = f(d(p_i,p_j))` _is_ the filtration weight — so the thermal axis and the similarity filtration are **coupled by construction**, not bolted together. The filtration threshold _is_ the temperature.
- **The object:** a warped path `s ↦ (λ(s), τ(s))`, topology computed _along the path_ — explicitly **not** a 2-parameter persistence module. ("Scaffolded latent persistence" / a thermal barcode over a sequentially-scaffolded graph.)
- **Metric:** the principled geodesics you've already built — Fisher–Rao (simplex/PMF _and_ hyperbolic) and Poincaré — with the contrastive question (Euclidean vs simplex-FR vs hyperbolic) being a genuine research probe of whether curvature reveals structure Euclidean mangles, not a convenience.
- **DR front-end:** SPRED (Betti-preserving), precisely _because_ it's not UMAP's "loops might be layout artifacts" — define the object upstream on the graph, project only to look.
- **Temporal carrier = zigzag.** Slide a window; H₁ classes born/die as context drifts; the zigzag bar-lengths _are_ cycle-persistence over time. This is exactly where the zigzag engine I inventoried cashes out — it's not a separate track, it's the time dimension of this same object.

And the conceptual spine that makes the thermo⋈topology weld non-arbitrary: the **measure-zero / Lebesgue argument** — a pun, a tie-back, a re-entrant cycle live on a measure-zero set, invisible to thermodynamic averaging but visible to homology. Thermo weighs the bulk; topology reads the measure-zero skeleton; each is the other's blind spot.

## Correcting my earlier reconciliation

Strike "**SIFTS pipeline — fully unbuilt**" as a deliverable. The genuinely-unbuilt thing is the **`BackboneConditionedFiltration` primitive** (canonical brief's E1: factor it into `tda/ph`, reproduce the nursery-rhyme β₁ as a _special case_ of backbone=path) plus the **two-axis thermal warped-path construction** on top. SIFTS-proper (reading-order + cosine) becomes the degenerate baseline and the engine's first real-data validation vehicle — not the goal. Per the canonical brief's own src survey, the undirected case is _nearly free today_ (`CsrGraph = skeleton ∪ kNN-similarity` → `RipsFromGraph`), so most of what's "missing" is the union convention + the warped-path-over-T abstraction + ingestion, not new engine math.

## The design tension I'd want to settle before building

The sharp fork in "warped path between two axes": **is it one parameter or two?**

- If position is the _birth-0 backbone_ and T is the _single filtration parameter_ (via the coupling kernel), the object collapses to a clean **1-parameter backbone-conditioned filtration over T** — ordinary PH, or **zigzag** if the T-sweep is non-monotone. This is essentially persistent-Mapper-over-T (ThermoMapper layer C) with a backbone prior.
- If you genuinely traverse _both_ non-trivially — position advancing _and_ temperature changing along `s` — it's a path through 2-D `(T, s)` space, which is a real modeling choice (what is the warp?) and pulls toward zigzag along that curve.

You've ruled out the full bifiltration; what's underdetermined is whether the two axes are one-knob-with-a-prior or a chosen coupled trajectory. That choice decides whether this is "PH along a path," "zigzag along a path," and where the warp `s↦(λ,τ)` is specified. I won't pick it for you — but it's the first thing I'd pin.

## What I think I still owe before we attack

I've read the persistent-homology + zigzag spine and the backbone brief. Still unread and likely load-bearing: `mapper-integrations/thermomapper-concept-map.md` (the C/D layers this warped-path feeds, and the "reverse arrow"), and the Zhu 2013 / Li 2025 compendia the brief cites. Want me to read the concept-map next so I have the ThermoMapper layer semantics straight — or would you rather drive the warped-path-vs-zigzag question first and have me pull references on demand?

> Does it help that I already have fast union find in under src/graphs? i wrote it for SPC a while ago

<details><summary>2 previous messages</summary>

> Good question — and yes, but let me be precise about _where_ it helps and where it doesn't, because the answer forks exactly on the monotone-vs-zigzag distinction we just raised. Let me check the actual state of it in the tree first rather than guess.
>
> `src/graphs/primitives/UnionFind.cs` is clearly the canonical one and it's already load-bearing across SPC (`SwendsenWang`, `UnionFindLabeler`, `LineagePersistence`), MST (`Kruskal`/`Boruvka`), dendrograms, HDBSCAN, and the TDA side (`NerveDiff`, `CodimensionOneDualGraph`, `EmbeddedComplexZigzag`). Notably it never floored to `maths/structures` as `plan.md` proposed — everyone just reaches into `graphs` for it. Let me read its surface to see exactly how much of the H0 elder-rule it already gives us.

</details>

Yes — and materially. That file is exactly the kind of asset the project's "reuse, don't rebuild, floor to a shared primitive" discipline is built around, and it's _already_ load-bearing on both sides of this work: SW labeling ([SwendsenWang](src/clustering/graphical/spc/runtime/core/sampler/SwendsenWang.cs)), MST/single-linkage, HDBSCAN, and the TDA side ([NerveDiff](src/tda/primitives/NerveDiff.cs), `MapperConnectedComponents`, the codim-one duals). So reusing it for backbone-conditioned persistence isn't opportunistic — it's the same primitive the rest of the stack already trusts.

**Where it directly serves the backbone-conditioned object:**

- **The H0 pass + clearing partition.** Backbone-conditioned persistence anchors H0 (the backbone connects everything at birth 0) — so H0 itself is cheap, but the _clearing partition_ (`to_skip`/`to_reduce`) that hands columns to the H1 reducer comes from exactly a union-find elder-rule pass. The interesting object (returns = H1) rides on top of this. One small extension: it tracks `parent`+`size` but not **birth** — the elder rule needs a birth stamp per root to assign deaths. `plan.md` already anticipated this ("extend `UnionFind` → birth-tracking DSU"); it's a parallel `birth[]` array, not a rewrite. And [Reroot](src/graphs/primitives/UnionFind.cs:72) already gives you the merge-tree internal-node shape (it's the HDBSCAN/Kruskal mechanism), which is the same structure H0-elder-rule and merge-trees want.
- **The thermal axis.** [Reset](src/graphs/primitives/UnionFind.cs:84) — zero-realloc reuse across SW sweeps — is a real win for the T-sweep: each temperature re-labels the lattice cheaply, which is precisely persistent-Mapper-over-T's inner loop. You built that reuse _for_ SPC; it pays off again as the thermal carrier.
- **Mapper/nerve node construction.** [GetLabels](src/graphs/primitives/UnionFind.cs:97)/[WriteRootSizesTo](src/graphs/primitives/UnionFind.cs:119) are zero-alloc component extraction — Mapper nodes _are_ connected components of the pullback cover.

**The honest limit — and it maps onto last turn's fork.** This UF is **monotone**: `Union` only ever merges, there's no delete. That's perfect for the _absolute_ half — backbone-at-base, fixed-T clustering, a one-directional T-sweep. But the **zigzag / temporal half** (slide a window; cells _leave_; the non-monotone warp) needs deletions, which a plain union-find structurally can't do. That's exactly why `src/tda/ph/dynamic/` exists — `DynamicConnectivity`, `EulerTourTree`, `LinkCutTree`, `DynamicMsf` — the dynamic-connectivity family that `GraphZigzag` rides instead of this UF.

So the union-find actually _answers_ the design fork I flagged: if the warped path is **monotone in T** (position = birth-0 backbone, T = single increasing filtration parameter), this UF + birth-tracking is nearly the whole H0/clearing engine, and the sweep reuse is free. If the warp is **non-monotone** (genuine zigzag along the `(T,s)` curve), you're in dynamic-CC territory. Either way you already own both substrates — the static one you wrote for SPC, the dynamic family already in the zigzag stack. The choice of warp decides which one carries the load; it doesn't require building a new one.

(Footnote: `plan.md` proposed floating `UnionFind` up to `maths/structures` under the second-pure-consumer rule. That never fired — and it doesn't need to here. `BackboneConditionedFiltration` lands in `tda/ph`, `tda` already depends on `graphs`, so it consumes `Graphs.Primitives.UnionFind` directly. No float, no duplication.)

> does the dynamic stuff work on csrgraphs?

<details><summary>1 previous message</summary>

> Let me check what those dynamic structures actually operate on rather than infer it — the CSR question is architecturally loaded (CSR is a _static_ layout, and dynamic connectivity needs edge insert/delete).

</details>

No — and deliberately. They don't touch CSR at all. Every one of them is allocated with just a vertex count and driven by edge deltas over `int` ids, with zero `Graphs.*` / `Csr` dependency (the only `using`s are `System` / `Collections.Generic`):

- [`DynamicConnectivity(int n)`](src/tda/ph/dynamic/DynamicConnectivity.cs:30) → `Insert(int u, int v)` / `Delete(int u, int v)` / `Connected(u,v)`
- [`DynamicMsf(int n)`](src/tda/ph/dynamic/DynamicMsf.cs:28) → `Insert(int u, int v, int w)` / `Delete(int u, int v)` — note the **weight** `w`, i.e. the weighted MSF you need for the H0 elder-rule _with deletions_
- [`EulerTourTree(int n)`](src/tda/ph/dynamic/EulerTourTree.cs:65) → `Link` / `Cut` / `Connected`
- [`LinkCutTree(int n)`](src/tda/ph/dynamic/LinkCutTree.cs:28) → path-max + `Connected`

They maintain their **own** mutable forests (splay trees / Euler tours / the HDT spanning-forest hierarchy). That's not an oversight — it's forced: `CsrGraph` is a _static_, compressed adjacency layout; you can't insert or delete an edge in it without rebuilding the row offsets. Dynamic connectivity structurally requires mutable structures, so a CSR-backed dynamic oracle is a contradiction in terms.

**How they compose with CSR, then:** CSR is the static snapshot; you bridge by _streaming_ its edges (or the zigzag filtration's add/remove events) into `Insert`/`Delete`. The flow is `CsrGraph → edge-event stream → dynamic oracle`, never "oracle reads CSR." For the warped-path/zigzag-over-T object, that's exactly the wiring: build the backbone+similarity graph as a CSR snapshot per T (where `UnionFind` rides for the absolute side), and feed the _deltas_ between T-steps into `DynamicMsf`/`DynamicConnectivity` for the non-monotone side.

**This is the precise gap the void-inference doc §8 earmarks BP-CSR for.** Its backend table is the same split we keep hitting: SW union-find (fixed lattice) → static **CSR**; filtration/**zigzag** persistence → **BP-CSR** (because the complex grows _and_ shrinks in batches). BP-CSR (Wheatman–Burns–Xu 2024) is precisely "a CSR-shaped store that _is_ mutable — batch insert+delete at CSR-traversal speed." So the three-way picture:

- **Static `CsrGraph`** — backs the absolute/fixed-T side (`UnionFind`, SW, MST).
- **The `TDA.Ph.Dynamic` family** — backs the zigzag/temporal side, fed by edge deltas, _not_ CSR-backed. Already built, already CSR-independent, so it composes with whatever emits the edge stream.
- **BP-CSR** — the _unbuilt_ unifying mutable-CSR backend (C++ ref; the doc locks it as a `graphs`-layer primitive, port cost lands there) if/when you want one storage seam that serves both static and dynamic consumers.

One placement note in passing: because these dynamic structures import no `Graphs.*`, they're "pure" by the discriminator — they're general dynamic-graph primitives currently homed in `TDA.Ph.Dynamic` only because zigzag is the sole consumer (the defer-abstraction rule). The day a second consumer appears — void-merge tracking is the obvious one — they float cleanly to `maths/structures` or `graphs`, same `maths`↔application logic you described for `topology`/`ph`.

> and what about the Xin Li stuff that i connected to SIFTS in my generalized temporal backbone?

<details><summary>1 previous message</summary>

> Let me ground this rather than riff on the brief's one-line gloss — I want to see whether the Li primary source is captured in your compendia and whether the engine has any of the directed/path-homology machinery the Li-faithful version would need.

</details>

Grounded answer first: **directed path-homology is absent from `src`** (no `PathHomology`/`GLMY`/`DeltaHomology`/directed-flag anything), so the "Li-faithful" version is genuinely new engine math, not a reuse. And there's no primary-source digest of Li in the tree that I can find — the term cluster (`polychronous`, `delta-homology`, `path homology`, `PNG`) only shows up _inside the briefs_, and tellingly it co-occurs with the **complex-analytic / magnetic-Laplacian** docs, not just the backbone one. So I have the briefs' gloss of Li, not the paper. With that caveat:

**Where Li sits in your generalized backbone.** It's the third domain — text/SIFTS, spikes/Li, decision/COSYNE — and it's not just "another instance." It's the one that **introduces direction**. SIFTS (reading order) and COSYNE (belief/pupil proximity) use the backbone _undirectedly_ — similarity is symmetric, the return is a plain H₁ loop. Li's backbone is **axonal-delay temporal consistency**, and a PNG is a _delay-locked, directed_ firing chain — so its "return" is a **re-entrant directed cycle**, and throwing away the arrow throws away the science. Li is the stress test that breaks the undirected-Rips assumption the other two get away with.

**That forces the fork the brief names (now source-confirmed):**

- **v1 — undirected H₁ + orientation post-check on `CycleReconstruction`.** Reuses the existing engine entirely; "covers the science" for detecting the loop, just not its causal orientation natively. Nearly free.
- **v2 — directed path homology (Grigor'yan–Lin–Muranov–Yau).** A _different_ homology theory (path complex over directed graphs, not the simplicial PH the engine computes). **Not in the tree** — this is a real new build, not a port-and-reuse.

**The cross-link I'd flag as the interesting part.** The directed side isn't a dead-end spikes-only detour, and it isn't only reachable via GLMY. There are _two_ routes to "directed," and you already have substrate for one:

- **Combinatorial:** GLMY path homology (unbuilt).
- **Spectral:** the **magnetic / Hermitian Laplacian** — which encodes edge direction as phase, and whose operator _is already built_ ([`MagneticLaplacianOperator`](src/graphs/spectral/MagneticLaplacianOperator.cs)). That's exactly why the Li terms co-occur with the complex-analytic-thread and magnetic-Laplacian briefs — directed cycles ↔ phase-encoded Laplacian ↔ directed homology is one thread, and you've half-built its spectral leg.

And the payoff isn't confined to neurons: your **proximal NLP target inherits this too.** A SIFTS tie-back (conclusion returns to intro) _is_ a directed return — undirected H₁ sees the loop but discards whether the discourse went forward-and-back or the reverse. So directed/path homology would give the _text_ barcode a causal/discourse-direction reading, not just the spike barcode. The Li-faithful machinery upgrades all three domains, which is more reason to treat it as v2 rather than a one-off.

The honest next step, since I'm working off the brief and not the paper: I can pull `opus-brief-complex-analytic-thread.md` + `roadmap-complex-analytic-ph-papers.md` (where the Li/path-homology/magnetic material actually lives) and, if you have it, the 2508.11646 source — before we decide whether the temporal backbone commits to v1 (ship undirected now) or earmarks v2 (directed, and which route). Want me to go read that thread, or is the directed-vs-undirected call something you'd rather make first?

> yea i want to drop the anchoring in neurons. thats a downstream applicaiton. i want to keep this framed as persistent homology ground work. but here's a raw ingest of the xin-li paper `D:\aghado01\codex-scientiae\ingestion\compendia\ph\xin-li\2508.11646v1`

<details><summary>2 previous messages</summary>

> Good — reframing to PH groundwork, neurons as downstream. Let me read the actual paper. First let me see what's in the ingest directory.
>
> This is dense and substantive. I've read through Section III (the core topological machinery); Sections IV–VI are the CCUP/cognitive model you want to drop. Let me map the remaining formal content surgically rather than read the cognitive framing in full.

</details>

Read it (Sections II–III in full — the topological machinery; IV–VI are the CCUP/cognitive model you're dropping, and I mapped their formal content without dwelling on the cognition). Stripped of neurons, here's what the paper actually is as PH groundwork — and there's one correction that changes the build cost materially.

## What the construction is, in neutral terms

Strip "spikes/delays/memory" and Li's object is: a **directed flag complex built from a vertex scalar field + an edge-level generative prior, filtered by prediction residual.**

- Vertices carry a scalar `t_i`; edges carry an offset/prior `τ_ij`. A directed edge `i→j` is admitted iff `|t_j − (t_i + τ_ij)| < δ` (Def 1/2). So `t_i + τ_ij` is a **prediction** of where `j` should sit, and the filtration parameter δ is the **residual slack** you allow. Higher simplices enter when all ordered faces are consistent.
- Increasing δ → more edges → a monotone filtration → persistent H₁. A cycle = a **return relative to the prior**; the "delta-homology" framing is just: the meaningful features are the sparse, irreducible H₁ generators (`ker ∂₁ \ im ∂₂`), Z₂ throughout.
- There's a **second parameter** Δ (a cap on `τ_ij`), so `K_δ^Δ` is genuinely a **bifiltration** in (tolerance, prior-magnitude).

## The correction that matters for the engine

The earlier brief said "Li-faithful = directed **path homology** (GLMY)." **That's not what the paper does.** Its reference (line 932) is **Lütgehetmann–Govc–Smith–Levi 2020, "Computing persistent homology of directed flag complexes"** — i.e., ordinary simplicial PH over a complex whose simplices are _ordered_. The "directed" content is in the **construction** (respect edge orientation, build ordered simplices), and the homology is computed by the **same Z₂ column reduction the engine already has.** So the directed capability is a **`FlagComplex` build variant**, not a new homology theory. That de-risks the whole "v2 needs new algebra" worry — it's a builder on top of the existing reducer, not GLMY.

(Current [`FlagComplex`](src/tda/primitives/FlagComplex.cs) is undirected; the gap is a directed/ordered variant, which the cell-generic engine — [`AbstractCell`](src/tda/ph/AbstractCell.cs)/`CellFiltration` — is already shaped to host.)

## How it upgrades backbone-conditioned persistence

This is the real payoff, and it generalizes your backbone primitive one level up:

- **SIFTS backbone** = "zero-distance path edges at birth 0" — a fixed `K₀` subcomplex.
- **Li backbone** = "an edge-level **generative prior** `t_i + τ_ij`, with the filtration measuring **deviation from the prediction**." The backbone stops being a static subcomplex and becomes a _model that predicts where the next point should be_; persistence reads how much slack closes a cycle.

So `BackboneConditionedFiltration` generalizes from "include `K₀` at ε₀" to "**admit edges by residual `|observed − predicted| < δ` against a prior**" — strictly more expressive, and SIFTS falls out as the degenerate prior (`τ ≡ 0`, undirected). And note the shape: Li's `(δ, Δ)` two-parameter complex is _exactly_ your "two time-like axes / warped path, not a full bifiltration" question, recurring verbatim. The fork I flagged two turns ago — warped path vs genuine bifiltration — is the same fork Li sidesteps by mostly working 1-parameter in δ.

## The second half the paper adds: a "context" layer dual to the homology

Beyond the construction, Li bolts a **cellular-sheaf / cohomology layer** onto the homology (Def 4, Thm 2, and the appendix): stalks per cell, restriction maps along covering relations, **coherence = local sections glue into a global section**, and failure = a cohomological obstruction `H¹(F) ≠ 0` (line 1096). Framed as Φ (homology = "content") ⊥ Ψ (cohomology/sheaf/filtration = "context") with a pairing `⟨α, γ⟩` (line 655, 1162).

Engine-wise that splits cleanly:

- The **Φ / homology half** and the absolute-vs-relative pairing is your existing pCoh + the **VSMJ2011 relative-(co)homology** feature already scoped in the port — not new.
- The **Ψ / cellular-sheaf-cohomology half is genuinely new capability** — and it's the _same_ thing the complex-analytic thread calls "complex-valued cellular sheaves / persistent sheaf Laplacian," whose spectral cousin ([`MagneticLaplacianOperator`](src/graphs/spectral/MagneticLaplacianOperator.cs)) you've already built. So Li is the construction that makes the magnetic-Laplacian/sheaf thread _concrete_ rather than speculative.

## Net for the groundwork inventory

As PH groundwork (neurons dropped), Li contributes four primitives, three of them cheap reuse:

| Primitive                                                      | Engine status                                                                   |
| -------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| Directed/ordered flag complex (Lütgehetmann)                   | **build variant** on `FlagComplex` → existing Z₂ reducer; _not_ GLMY            |
| Prediction-residual ("band") filtration = generalized backbone | new `IFiltration` implementor; generalizes `BackboneConditionedFiltration`      |
| Cell-poset homology (combinatorial compression)                | largely the existing **cell-generic** engine (`AbstractCell`/`CellFiltration`)  |
| Cellular sheaf cohomology + global-section/`H¹` obstruction    | **genuinely new**; the bridge to the complex-analytic/magnetic-Laplacian thread |

So Li doesn't push you toward exotic new algebra — it says the directed, prior-conditioned construction rides your existing reducer, and the only net-new engine capability it demands is the cellular-sheaf "context" layer, which you already have a foothold on. The cleanest framing: **`BackboneConditionedFiltration` should take a _prior over edges_ (residual-band admission), with reading-order/zero-prior as the SIFTS degenerate, directedness as a flag-complex build flag, and the sheaf layer as the optional Ψ companion** — that single generalization absorbs SIFTS and Li as instances and leaves the warped-path-vs-bifiltration call as the one open modeling fork.

Want me to pull the complex-analytic-thread + magnetic-Laplacian briefs next so the sheaf/Ψ half is fully specced against what's already built, or hold here?

> yea lets review those sources and write a new brief synthesis of the two as you describing with the new material gleaned from xi and integrated in to a plan

<details><summary>3 previous messages</summary>

> Let me pull the complex-analytic thread sources and check the current state of the magnetic-Laplacian code, so the synthesis is specced against what's actually built.
>
> The complex-analytic thread is well-developed and — crucially — it's the **spectral** half of exactly the object Li gives the **combinatorial** half of. Seam A's "a U(1) magnetic Laplacian _is_ a cellular sheaf (ℂ stalks, restriction maps = e^{iθ})" is the _same sheaf_ Li builds over his cell poset (Def 4, coherence = global section, obstruction `H¹(F)≠0`). The synthesis writes itself. Let me verify the three code anchors the plan will lean on, then write the brief.
>
> Code-grounded, and the picture is clean:
>
> - **`MagneticLaplacianOperator`** — `L^(q) = D − A∘exp(i·2πq·Θ)` — is **built**, matrix-free over `CsrGraph`, with a real-`2n×2n` embedding so the existing **`LOBPCG`** solves the complex-Hermitian eigenproblem unchanged. So Seam A's v1.5 spectral operator already exists.
> - **`LOBPCG`** + dense/fast eigensolvers exist; **FEAST/contour (Seam B's exact harmonic projector) does not** — λ_q is computable today via LOBPCG, contour-integral is the later optimization.
> - **Sheaf / persistent / connection Laplacian: absent** — confirmed. The Ψ layer (Li's sheaf ⟷ Wei–Wei Persistent Sheaf Laplacians umbrella) is the genuine new build.
>
> That's everything I need. Writing the synthesis now.

</details>

Written — [opus-brief-conditioned-persistence-synthesis.md](../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/opus-brief-conditioned-persistence-synthesis.md) in `persistent-homology/` (I put it there, not `mapper-integrations/`, because you reframed this as PH groundwork — neurons/ThermoMapper-layer-D are downstream; the cross-links into the complex-analytic cluster are relative). Not staged, not committed.

**The synthesis claim** (the spine of the brief): the two threads aren't siblings that feed each other — they're the **homological (Φ) and spectral (Ψ) projections of one primitive**, _Conditioned Persistence_: a filtration admitted by **residual against a prior** (the backbone), built on a **directed** skeleton, **read twice** — homology (which returns are real) and a connection/sheaf (whether they cohere, how far they reach, their continuous order parameter). Li is what fuses them, because his paper builds _both halves over one object_ and nothing else.

**The load-bearing identity** that makes it one object, not two bolted together: _a U(1) magnetic Laplacian **is** a cellular sheaf_ (ℂ stalks, restriction maps `e^{iθ}`) — so the `MagneticLaplacianOperator` already in your tree is the first concrete instance of **Li's contextual sheaf Ψ**, and Wei–Wei's Persistent Sheaf Laplacians is the umbrella that makes it persistent. Li gives the combinatorial sheaf; the complex-analytic thread gives its spectral engine; they're the same Ψ.

**Two corrections that _lower_ the build cost** (both from actually reading Li):

- Li is **directed-flag PH (Lütgehetmann), not GLMY** — directedness is a `FlagComplex` build flag riding the existing Z₂ reducer, not a new homology theory. GLMY stays the _true_ v2, now optional.
- **SIFTS is the `τ≡0`, undirected degenerate** — the general primitive subsumes the seed rather than sitting beside it.

**The plan** is five reuse-maximal phases: **P0** conditioned core (Φ, undirected — reproduce SIFTS β₁ as the degenerate) → **P1** prior generalizes to the residual band (SIFTS/Li/thermal as prior _choices_) → **P2** directedness, _both_ readings (directed flag + the already-built magnetic flux/`λ_q`, v1.5 before GLMY v2) → **P3** the sheaf layer (Li coherence ⟷ sheaf Laplacian, magnetic = U(1) instance) → **P4** `λ_q` as the continuous observable + the one **genuinely open research term**: the magnetic/sheaf _persistent_ Laplacian `dim ker(α,q)` pairing, obtained as a _specialization_ of the persistent-sheaf-Laplacian umbrella. P0–P2 ride entirely on built primitives; only P3–P4 are net-new, and they're additive — the barcode never depends on them. **Ship Φ, earn Ψ.**

Three things I'd flag for your pressure-test:

- The **warped-path-vs-bifiltration fork** is recorded, not resolved (`(δ,Δ)` / `(T,position)`) — still your call, and P1 is where it gets made.
- I mapped **Li's homology–cohomology dual loop → VSMJ2011 relative-(co)homology (content) + sheaf (context)**. That's my interpretive bridge; it's clean but it's the seam most worth your eye.
- I asserted **Botnan–Dong's near-linear PL for _non-branching_ complexes is a direct hit** because your backbones are path-like — a real engine win if it holds, worth confirming against the actual complex shapes P1 produces.

Next, your pick: ingest the two staged anchors (Wei–Wei _Persistent Sheaf Laplacians_ 2112.10906, Botnan–Dong 2408.16741) so P3–P4 are specced against the real math; start **P0** against the built engine; or push on any of those three flags first.

> ask and ye shall receive `D:\aghado01\codex-scientiae\compendia\ph\2408.16741v2.md`, `D:\aghado01\codex-scientiae\compendia\ph\2508.11646v1.md` , `D:\aghado01\codex-scientiae\compendia\ph\2112.10906v4.md`

<details><summary>3 previous messages</summary>

> The two anchors I asked for, plus the canonical Li copy. Let me read the two new ones — Botnan–Dong (the near-linear PL engine claim) and Wei–Wei (the sheaf-Laplacian umbrella) — since those are what P3–P4 must be specced against.
>
> Botnan–Dong is clear enough to correct my brief already — and it's exactly the flag you told me to check. Let me note the key finding before moving on, then read Wei–Wei.
>
> **The Botnan–Dong applicability is narrower than my brief claimed, in an important way:** "non-branching" means each _q_-simplex bounds **≤2** (_q+1_)-simplices — i.e. **pseudomanifolds and cubical/image complexes**, _not_ an arbitrary path-like graph. Concretely:
>
> - A **pure backbone** (a path/sequence) _is_ 0-non-branching, so the H0 up-persistent-Laplacian is near-linear there — via **weak column reduction = union-find on the incidence graph** (literally the DSU you already own).
> - But the **conditioned complex is backbone ∪ similarity chords**, and the moment chords raise a vertex's degree past 2 (H0) or put an edge in >2 triangles (H1), it is **no longer non-branching** — so the near-linear guarantee does _not_ transfer to the full flag complex. It holds cleanly for **cubical/image (E4)** data and the pure skeleton, not the chorded returns.
> - The durable gifts regardless: **nullity(Δ_q^{K,L}) = persistent Betti number** (the formal basis for "λ_q is the continuous shadow of β_q"), and eigenvalues via **SVD/Golub–Kahan–Lanczos** on the sparse boundary matrices.
>
> Now Wei–Wei — the sheaf umbrella that P3 actually specializes.
>
> Both anchors confirm the synthesis — and Wei–Wei even more strongly than I expected. Two findings that sharpen P3–P4:
>
> 1. **Wei–Wei Example 2.3 _is_ the bridge, verbatim.** They define a cellular sheaf on a **path complex** (integers as vertices, `[n,n+1]` as edges — a backbone) with **`e^{±iθ}` unitary restriction maps**. That's a U(1)/magnetic sheaf. So "the magnetic Laplacian is a cellular sheaf, and its persistent version is a PSL specialization" isn't my extrapolation — it's the paper's own example. And `ker Δ_q^{X,Y} ≅ persistent sheaf Betti`, so the magnetic-persistent `(α,q)` operator has a _defined_ construction (instantiate ℂ stalks + `e^{iθ}` restriction maps in their §3.2); only the pairing/stability rule is genuinely open.
> 2. **Botnan–Dong is narrower than I wrote** (the flag you called): "non-branching" = pseudomanifold/cubical, and the near-linear reduction _is_ union-find — a real win for the **pure backbone** (a path is 0-non-branching) and **cubical/image** data, but **not** the chorded conditioned complex. The durable gift is `nullity(Δ_q) = persistent Betti`, which grounds `λ_q` as the literal shadow of `β_q`.
>
> Let me fold these into the brief so P3–P4 are specced against the real math.

</details>

Done — the brief now specs P3–P4 against the actual math, with the two findings folded in and the overstatement flagged. [opus-brief-conditioned-persistence-synthesis.md](../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/opus-brief-conditioned-persistence-synthesis.md).

**What the anchors changed:**

- **Wei–Wei _shrank_ the open research term.** Their Example 2.3 — a path complex with `e^{±iθ}` restriction maps — _is_ a backbone magnetic sheaf, so the identification isn't mine to defend, it's the paper's. And §3.2 gives the persistent sheaf Laplacian with `ker ≅ persistent sheaf Betti`. So the magnetic-persistent `(α,q)` operator is now a **defined construction** (ℂ stalks + `e^{iθ}` restriction maps in their PSL), not an "umbrella to aim at." What's genuinely open collapsed from _"define the operator"_ to _"define the `(α,q)` pairing/stability rule"_ — a real narrowing. Their **labeled-sheaf** mechanism (per-vertex quantity riding the restriction maps) also turns out to be _exactly_ Li's context-modulation Ψ — the prior enters through the restriction maps, not the coefficients.

- **Botnan–Dong got corrected** (your flag was right). "Non-branching" = pseudomanifold/cubical, and the near-linear step is **union-find** (the DSU you own) — a genuine win for the **pure backbone** (a path is 0-non-branching) and **cubical/image** data, but it **fails** the moment similarity chords branch the conditioned complex. I pulled the "direct hit for backbones" claim and replaced it with a scoped caveat. The durable gift survives: `nullity(Δ_q) = persistent Betti`, which grounds `λ_q` as the _literal_ shadow of `β_q` rather than an analogy.

**Net effect on the plan:** P3–P4 are no longer "new research stack" — they're "specialize a defined theory + own the eigensolve," with only the two-parameter `(α,q)` pairing rule as open contribution (and the roadmap's Wolf–Fan–Monod 2509.20220 already warns which Laplacian components keep monotonicity/stability — the lead on that rule). P0–P2 remain buildable on what's in the tree today.

Two things now sit in front of you:

1. **The P1 warped-path-vs-bifiltration fork** is still your call, and it's the one that gates everything above P0 — Li's `(δ,Δ)` and your `(T,position)` are the same two-parameter object, and whether you commit to a _chosen path_ or a genuine bifiltration decides whether P2's temporal reads are PH-along-a-path or zigzag.
2. **P0 is buildable now** — `ConditionedFiltration = backbone(prior, birth ε₀) ∪ similarity(metric)` → `RipsFromGraph` → existing reducer, with SIFTS β₁ as the degenerate test. No new engine math.

Want me to resolve the P1 fork with you first (I'd argue it, not just list options), or scope P0 into a concrete build ticket against the existing engine?

> regarding 1. I don't see why i have to choose. I think I want both since they probably have different applications. A whole other part of the PH unfinished work is my zigzag capabilities and fast implementations. Also I want to step back on my two parameters (T, position) feels kind of bespoke. My conditioned backbone isn't specifically about time, or temperature, or reading order. The frameing there feels of

> [Request interrupted by user]

> regarding 1. I don't see why i have to choose. I think I want both since they probably have different applications. A whole other part of the PH unfinished work is my zigzag capabilities and fast implementations. Also I want to step back on my two parameters (T, position) feels kind of bespoke. My conditioned backbone isn't specifically about time, or temperature, or reading order. The frameing there feels off

You're right on all three, and the third one is the key that fixes the other two. Let me not use Workflow here (your CLAUDE.md standing rule overrides the Ultracode nudge) — this is a framing question to think through, not fan out.

## Why `(T, position)` felt off: `position` isn't an axis

In "`(T, position)`," **`position` is the backbone prior `K₀`, not a persistence parameter.** The backbone sits at birth-0 — it's the _thing you condition on_, the same slot as reading-order / delay-prediction / task-events. Promoting it to a second "time-like axis" is a category error, and that's the bespoke smell. The backbone is a _choice of prior_, exactly parallel to "which metric." It never was an axis.

So the genuine axes are not `(T, position)`. They are:

- **one content axis** — the similarity/residual **filtration scale** (the homology direction; this is always there), and
- **one-or-more conditioning parameters** — whatever modulates the prior/context.

This isn't my coinage — it's your own **Li Φ/Ψ** (content vs context) and **Wei–Wei's "geometric vs non-geometric information"** said twice. The conditioning is the sheaf/label side; the filtration is the geometry side. `(T, position)` collapsed a prior and a parameter into a fake symmetric pair.

## Why you don't have to choose: the conditioning parameter comes in two kinds

Here's the distinction I think dissolves the fork — and it's the one worth pressure-testing:

- **Filtration-type conditioning** (grows the complex, monotone): Li's `Δ` (max prior magnitude), a backbone-inclusion scale, an SPC `T` where clusters merge. These are _genuine multiparameter-persistence axes_. `K_{ε, c}` is then a real **bifiltration**, and —
  - a **monotone slice/path** through it → ordinary PH,
  - a **non-monotone slice/path** → **zigzag** (`T`-sweeps are non-monotone — clusters split as well as merge — which is _why_ persistent-Mapper-over-T is honestly a zigzag, not a clean bifiltration line),
  - the **whole thing** → the full multiparameter module.
- **Gauge-type conditioning** (fixes the complex, changes the Ψ operator): the magnetic **charge `q`**, a connection phase. A `q`-sweep is a _family of sheaf-Laplacians on one complex_ → a **spectral surface `λ_q`**, not a bifiltration at all. Orthogonal to the above.

So "path vs bifiltration" was never a choice — **the multiparameter module is the object; the path (monotone → PH, non-monotone → zigzag) and the full module are just how you read it, and the gauge sweep is a third, orthogonal reading.** You get all of them because they're one module and its slices, plus a spectral family over it. Different applications, exactly as you said: a canonical traversal → a slice; no canonical traversal → the module; directedness/order-parameter → the gauge surface.

## This unifies your "other part" — zigzag isn't separate, it's the slice-reader

The zigzag engine and its fast implementations aren't a parallel workstream to conditioned persistence — **zigzag _is_ the engine that reads the non-monotone slices** of the conditioned multiparameter filtration. The frontier work I inventoried (ReflectionZigzag stage 2b-ii, the fast graph-zigzag paths, the dynamic-connectivity substrate) is precisely the machinery this object needs the moment a slice isn't monotone. And the thing that reads **both at once** — a non-monotone traversal _through_ a multiparameter space — is exactly what `zigzag-frontier.md` already names as the far, unstarted rung: **Z6 "quasi-zigzag = multiparameter × zigzag."** Your "I want both" _is_ Z6. The two "parts" of the unfinished PH work were one program the whole time; Z6 is their intersection.

So the corrected picture, de-bespoked:

| reading                       | when                                                    | engine                                        |
| ----------------------------- | ------------------------------------------------------- | --------------------------------------------- |
| monotone slice                | a canonical increasing traversal                        | ordinary PH (built)                           |
| non-monotone slice            | trajectories, sliding windows, `T`-sweeps               | **zigzag** (built + frontier)                 |
| full module                   | no canonical path; joint content×conditioning invariant | multiparameter persistence (mostly unbuilt)   |
| gauge sweep                   | directedness / order parameter                          | magnetic/sheaf `λ_q` surface (operator built) |
| non-monotone × multiparameter | both at once                                            | **Z6 quasi-zigzag** (unstarted far rung)      |

I'd rewrite **P1** to drop the "pick a warp" fork entirely and instead state: _the conditioned filtration is multiparameter (content scale × conditioning), read via monotone slices (PH), non-monotone slices (zigzag), or the full module — with gauge parameters as an orthogonal spectral family; the backbone is the prior `K₀`, never an axis._ And fold the zigzag frontier into the plan as the non-monotone-slice reader rather than a separate track.

Two things I'd want your read on before I edit: (1) does the **filtration-axis vs gauge-axis** split land for you, or do you want them treated uniformly (some multiparameter-persistence people would fold the gauge into a "parameter" too)? (2) do you want the brief to name the general object in **Li's Φ/Ψ** vocabulary, **Wei–Wei's geometric/non-geometric**, or the plain **multiparameter (content × conditioning)** — you draw vocab by fitness, so I'd rather you pick the register than have me impose one.

> I think i'd keep guage is its own axis. I'm inclined towards Li's notation and general conceptualization , so potentially inclined to adopt `Φ/Ψ` and `δ,Δ`, thoughts? alot of li's more lofty stuff i dont necessarily like but some of his operational abstractions kind of appeal

Yes to Li's operational layer, with two scrubs and one sharpening — and gauge as its own axis fits cleanly into it. Here's my read on adopting `Φ/Ψ` and `δ,Δ`:

## The split I'd make in Li's vocabulary

| Li's abstraction                                                       | verdict                                  | why                                                                                                                                                                                                                                                |
| ---------------------------------------------------------------------- | ---------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`Φ/Ψ`** (content ⊥ context)                                          | **adopt** — as the dual reading          | Clean, compact, and it _is_ Wei–Wei's geometric/non-geometric said twice. `Φ` = the homology (which returns are real), `Ψ` = the sheaf/connection (the conditioning). Your instinct is right.                                                      |
| **`δ, Δ`** (the two filtration params)                                 | **adopt the symbols, scrub the meaning** | `δ` = residual/slack (the primary content scale), `Δ` = prior reach/horizon (the second filtration axis). Keep the letters; **drop "timing tolerance / axonal delay."** They're not about time.                                                    |
| **cell poset** (combinatorial compression)                             | keep the idea, no new name               | Subsumed by your existing `CellFiltration`/`AbstractCell`.                                                                                                                                                                                         |
| sheaf coherence / global section / `H¹(F)`                             | **adopt**                                | It's Wei–Wei's sheaf cohomology anyway — real, not lofty.                                                                                                                                                                                          |
| **"delta-homology"**                                                   | **drop**                                 | It reifies nothing: it's just persistent H₁ over the `δ`-filtration. Naming it a "homology" invites the mistake that it's a new theory. The insight (sparse, irreducible, path-dependent returns) is real — credit it in prose, don't coin a term. |
| CCUP / content-context _uncertainty principle_, low/high-entropy story | **drop**                                 | The lofty overlay. `Φ/Ψ` survive it intact — you keep Li's _assignment_ (Φ=homology, Ψ=sheaf) and lose his _justification_.                                                                                                                        |
| memory-amortized inference / perception–action loop                    | **drop**                                 | The cognitive narrative; the operational residue is just the Φ↔Ψ pairing, which you already have.                                                                                                                                                  |

## The one sharpening — pin `Ψ`, don't inherit Li's grab-bag

Li lets `Ψ` mean _"a filtration, **or** a cohomology class, **or** a sheaf"_ — three different things. That muddiness is worth fixing on adoption: **the filtration `(δ,Δ)` is `Φ`-side** (it produces the content homology); **`Ψ` is specifically the sheaf/connection.** Don't let `Ψ` mean "filtration." That's a real precision gain over Li, and it's exactly what makes gauge fall out cleanly:

- **`Φ`** — content homology, read off the `(δ, Δ)` filtration axes.
- **`Ψ`** — the context sheaf/connection, carrying the **gauge `q`** as _its_ axis.

So the object is indexed `(δ, Δ ; q)`: two `Φ`-filtration axes and one `Ψ`-gauge axis — your "gauge is its own axis," placed precisely. And the axes appear _as the structure demands_, not as a fixed bespoke tuple: **`δ` always; add `Δ` when the prior has a tunable reach; add `q` when there's a connection.** SIFTS is `δ` only (τ≡0, no reach, no gauge); Li lights up `Δ`; the directed/magnetic case lights up `q`. That's what de-bespokes it — no forced two-axis pair, just a content scale plus whatever conditioning the prior actually carries.

## One discipline note (so the vocab doesn't over-reify)

Keep `Φ/Ψ/δ/Δ` as **notation in the prose and docstrings** — the locating narrative — but let the **code identifiers stay functional** (`ContentFiltration`/`ContextSheaf`/`residual`/`priorReach`, not a `PhiPsiModule` class). Same discipline you already apply to math vocabulary: it's reasoning scaffolding, first-class in the narrative, not reified into the type system.

Net: you get Li's genuinely good operational abstractions (`Φ/Ψ` duality, the `δ,Δ` residual-band construction, the context-as-sheaf) credited to him, with the CCUP/memory/"delta-homology" scaffolding scrubbed and `Ψ` tightened from three meanings to one.

Want me to fold this into the brief — rename the axes to `(δ, Δ ; q)` with `Φ/Ψ`, add a short "adopted from Li, with these scrubs" provenance note (crediting the operational seed, naming the dropped heuristics honestly), and delete the residual `(T, position)` framing? That'd also let me correct P1 to the multiparameter-module-plus-slices picture in the same pass.

> yes to all of the above

<details><summary>2 previous messages</summary>

> Applying all of it. Let me read the current brief to edit precisely against its actual state.
>
> All edits confirmed against the current text. Applying them.

</details>

Done — all four folded in:

- **Axes renamed to `(δ, Δ ; q)`** with `Φ/Ψ`, and a new **"Notation & axes"** subsection in §2 pins them: `δ` = residual/slack (Φ content scale), `Δ` = prior reach (Φ, appears only when the prior scales), `q` = gauge (Ψ-side, its own orthogonal axis). Plus the explicit rule that the **backbone `K₀` is not an axis** and axes appear _as the structure demands_ (SIFTS = `δ` only; Li lights `Δ`; magnetic lights `q`).
- **Provenance + scrubs** recorded honestly: credits Li's operational layer (`Φ/Ψ` duality, `δ`-residual-band, context-as-sheaf); drops the timing/delay reading of `δ,Δ`, "delta-homology," CCUP/entropy, and memory-amortized inference — noting `Φ/Ψ` survive because we keep his _assignment_ and drop his _justification_. Includes the **Ψ sharpening** (pinned to the sheaf, not Li's filtration-or-cohomology-or-sheaf grab-bag) and the **code-stays-functional** discipline.
- **`(T, position)` deleted** everywhere (§2, P1, §7) — it was the bespoke framing; `T` now only appears as one _instance_ of a non-monotone slice.
- **P1 rewritten** to the multiparameter-module-plus-slices picture: `K_{δ,Δ}` read as monotone slice (PH) / non-monotone slice (**zigzag** — with the zigzag engine + frontier named as _the reader_, not a separate track) / full module (deferred to **Z6 quasi-zigzag** = "both at once"), gauge orthogonal. §7 defers updated to match.

The brief now says, in one line: _one prior `K₀`, a content filtration `(δ,Δ)` read by PH or zigzag or the full module, and a context sheaf `Ψ` swept by gauge `q` — Li's operational abstractions, his loftiness scrubbed._ And it structurally unifies your two "unfinished PH" fronts: conditioned persistence and the zigzag/fast-implementation work are the same program, meeting at Z6.

Natural next step is to scope **P0** into a concrete build ticket against the existing engine (`ConditionedFiltration` = `backbone ∪ similarity` → `RipsFromGraph` → reducer, SIFTS β₁ as the degenerate test — no new math), since it's the one phase that's fully buildable today and gates nothing. Want me to write that ticket, or keep refining the design (e.g., how the filtration-axis-vs-gauge split shows up concretely in the `IFiltration`/`ContextSheaf` seam)?

> lets write that up i guess

> [Request interrupted by user]

> lets write that up i guess

Written — [p0-conditioned-filtration-brief.md](../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/p0-conditioned-filtration-brief.md). It's deliberately small: the only new code is the **union convention** (`backbone @ ε₀ ∪ similarity @ distance` → hand to the existing `RipsFromGraph`); everything downstream is the engine you already have, named to exact files. Two tests — a synthetic path+chord (`β₁`: 1 with the tie-back chord, born at its distance; 0 without) and the SIFTS nursery-rhyme degenerate — and a hard "if you're writing more than a builder + a test, stop, you're rebuilding" guard.

One small correction I baked in: I placed it in `src/tda/primitives` beside `RipsFiltration` (it builds a `CsrGraph`, so it's graph-coupled _construction_, not the pure `ph` engine) — the parent brief's "lands in tda/ph" was loose. Flagging so the two docs don't drift.

On the "i guess" — no pressure to build it now. P0 is genuinely the smallest step and it gates nothing, so it's a safe checkpoint to have specced whether or not you pick it up next; the design thread is equally live if you'd rather keep going there.

When you do want it built, the clean move per your repo rules is a **spawn_task chip** — it lands the ticket in its own hook-instrumented session/worktree rather than an in-process subagent. Want me to spin that up (I'd point it at the ticket + the named reuse), or leave it as a spec for the Copilot lane and get back to design?

> yea this ongoing confusion about what belongs in tda/ph is getting really annoying. tda/primitives is itself supposed to be reorganized into tda/ph because literally everything in primitives in PH related

<details><summary>1 previous message</summary>

> You're right, and I'm the one who added the churn — two turns ago I "corrected" the synthesis brief's `tda/ph` _to_ `primitives`, inventing a "graph-coupled builder tier" that was never supposed to exist. Killing that now. The clean architecture is **two tiers, not three**:
>
> - **`maths/topology`** — general topology primitives, reusable beyond PH (`Simplex`, `CombinatorialIndex`, `Z2LinearAlgebra`).
> - **`tda/ph`** — the PH engine, _everything_ PH-specific, including the graph-coupled builders (`RipsFiltration`, `FlagComplex`, the dual graphs). It's allowed to consume `CsrGraph` because `tda → graphs` is already a downward dep.
>
> **`tda/primitives` is a legacy bucket, not a layer** — everything in it is PH, and it's meant to dissolve into `tda/ph`. So new PH code goes straight to `ph`, and the synthesis brief was right the first time. Let me fix the ticket.

</details>

Fixed. And the reason this keeps biting you is worth naming plainly: **the confusion _is_ the unfinished `primitives → ph` fold.** It's the same structural item my first reconciliation flagged (the SIFTS-brief reorg: "dissolve `tda/primitives` → `tda/ph`") — as long as those ~10 PH files sit in a `TDA.Primitives` bucket, every new file re-litigates "ph or primitives?". Doing the fold is the permanent fix; until then the answer is just _"always ph."_

That fold is a mechanical refactor (move ~8 files `primitives → ph`, rename `TDA.Primitives → TDA.Ph`, update the `.csproj` globs and references, run RepoAudit) — exactly the kind of broad dotnet change your repo rules say to hand to the executor with a clear spec rather than have me grind inline. I can write it up the same way I did P0.

One thing the fold-spec has to decide, so flagging it rather than silently routing: the two **nerve** files (`Filtration.cs`/`NerveFiltration`, `NerveDiff.cs`). The SIFTS brief R2 wanted those in `tda/mapper/nerves` (Mapper-owned, since `PersistentMapper` is their only consumer), not `ph`. You just said "literally everything in primitives is PH" — so either they go to `ph` with the rest, or they're the one exception to `mapper/nerves`. Your call; the other 8 are unambiguously `ph`.

Want me to write the `primitives → ph` fold ticket (and you decide the nerve-files line in it), or leave the reorg for now and keep the rule "new PH code → `tda/ph`, primitives is legacy" in force?

> i think my thinking would be to move mapper/nerves at least in part into something like ph/nerves because conceptually nerves are a ph thing not a mapper thing right? then whatever interfaces for ph stuff under mapper though i dont know what those would. this project is getting so confusing

Your instinct is right, and it _simplifies_ — it removes a distinction rather than adding one.

**Nerves are PH.** A nerve is the Nerve-Theorem construction: turn a cover into a simplicial complex that captures the space's homotopy type. That's algebraic topology; Mapper is just one _consumer_. So `NerveFiltration`/`NerveDiff` belong in `ph` — a `ph/nerves` subfolder is a fine cohesive home (or flat in `ph`; the repo leans flat until a concern grows, and this is ~2–3 files, so it's a judgment call, not a principle).

**What stays in mapper — this answers your "what interfaces" hesitation: nothing PH does.** The seam is a _cover_:

- **Mapper produces a cover** — lens → pullback → cluster → a set of vertex subsets. That's Mapper's whole job.
- **`ph/nerves` turns a cover into a nerve filtration** — cover → complex → `IFiltration`. PH's job.
- **the PH engine computes persistence.**

So Mapper depends _downward_ on `ph/nerves` and hands it a cover — there's no PH interface living "under mapper." The one concrete decoupling: `NerveFiltration` should take an **abstract cover** (`int[][]` member-indices), not a `MapperResult`; the `PersistentMapper` _orchestrator_ stays in mapper (because "run Mapper per parameter" is a Mapper app) and feeds the cover down. And you're basically already there — `PersistentMapper.BuildFiltration` already reduces each nerve to `Nodes[i].MemberIndices`, so the cover-as-`int[][]` seam just needs to become the _input type_ instead of something dug out of `MapperResult`.

That collapses everything to **one rule, one direction:**

> `maths/topology` (pure primitives) → `tda/ph` (PH engine: reducers, currencies, builders, **nerves**) → `tda/mapper` (cover production + orchestration) → `tda/pipelines` (apps)

No "graph-coupled builder tier," no "nerves are Mapper." (The _abstract_ nerve could float to `maths/topology` later if it grows a representation-agnostic core — same deferred call as `FlagComplex` — but it's CSR-coupled today, so `ph`.)

On "this project is getting so confusing" — the architecture isn't the problem; the **scatter** is. The placement rules live in `plan.md`'s split table, `ph-engine-placement.md`, `rips-flag-nerve-brief.md`, _and_ my own drift across this thread, none of them agreeing, so every new file re-derives them. The fix is one **canonical placement doc** — the single source of truth, the way `zigzag-frontier.md` is for zigzag status — stating the one rule above, with the others getting a "superseded → see X" banner.

Want me to write that canonical placement doc and fold the nerve decision into the `primitives → ph` fold ticket? That trades the current five-doc scatter for **one rule + one reorg spec** — which is the actual exit from the confusion, rather than us re-adjudicating placement every few files.

> sure write it up

<details><summary>2 previous messages</summary>

> Two artifacts, then. Before I write the fold ticket, one cheap check — whether the nerve files already take an abstract cover or still reference Mapper types (that decides whether the nerve-move is purely mechanical or needs the decoupling first).
>
> Confirmed — the nerve files already take an abstract cover (`CsrGraph Nerve` + `int[][] NodeMemberIndices`) with **zero** `MapperResult`/`TDA.Mapper` type dependency; "Mapper" appears only in docstrings. So the nerve move is purely mechanical, no decoupling. Writing both artifacts.

</details>

Both written:

- **[tda-placement.md](../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/tda-placement.md)** — the one rule. A four-tier ladder (`maths/topology → tda/ph → tda/mapper → tda/pipelines`, one-way), a four-step discriminator for placing any new file, and the three clarifications that kept biting written down flat (primitives is legacy, nerves are PH, "pure" ≠ maths/topology). It's `zigzag-frontier`-style: table-thin, "if a file and this disagree, one is the bug."
- **[primitives-to-ph-fold-brief.md](../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/primitives-to-ph-fold-brief.md)** — the reorg spec. Per-file move table (8 → `ph`, 2 → `ph/nerves`), the one real gotcha (the `Barcode.cs` filename collision → rename to `PersistenceBarcode.cs`), the verified-mechanical nerve move, no back-compat shims, and a RepoAudit + `TDA.Ph.Tests`-green gate. Scope guard: a pure relocation must not move a single test result.

One sequencing note I baked into the ticket: **do the fold before P0**, so `ConditionedFiltration` is born in `tda/ph` against a clean tree instead of chasing a moving target — otherwise P0 lands in the very bucket we're dissolving.

Two finishing moves, your call on each:

1. **Add the "superseded → see tda-placement.md" banners** to `plan.md`, `ph-engine-placement.md`, and `rips-flag-nerve-brief.md` — that's what actually collapses the scatter (the placement doc names what it supersedes, but those docs don't yet point back). Low-risk, ~3 one-line edits; I left it for you to okay rather than editing three existing docs unprompted.
2. **Spawn the fold as a task chip** for the executor lane (the broad dotnet move) whenever you want it built — I'd point it at the ticket + `tda-placement.md`.

That's the exit from the confusion: one rule, one reorg spec, and after the fold there's no `primitives` left to re-litigate.

> but dont those other docs like plan.md, ph-engine-placement etc cover things other than "where do we put files" question

<details><summary>1 previous message</summary>

> Good catch — you're right, and it changes the approach. A top-of-doc "SUPERSEDED" banner would be wrong: `plan.md` is mostly the purification _reasoning_ (the SPC→persistence cycle-avoidance, the `UnionFind`/`BoundedMinHeap` floor decisions, the `IFiltration` cut), `ph-engine-placement.md` carries the gate/SPRED flow-across-the-seam and the cocycle→conditioning notes, and `rips-flag-nerve-brief.md` has the naming survey + validation. Only the _placement tables_ in each are what `tda-placement.md` now owns. Blanket-retiring the whole docs would torch live content — exactly the "don't conflate scope" trap.
>
> So the pointer has to be **scoped to the placement section**, not the top of the doc: a one-line "for placement, defer to `tda-placement.md`; the rest stands" _at_ each doc's split-table / placement-table / homes-by-concern heading. The placement doc's `Superseded` wording implied whole-doc supersession — my sloppiness, not the design. Fixing it now:

</details>

Fixed — the placement doc's section is now "**Placement authority — these docs defer here on placement only**," and it spells out, per doc, exactly what stays authoritative (plan.md's floor/`IFiltration` reasoning, ph-engine-placement's gate/SPRED seam, rips-flag-nerve's naming + validation). It also carries the exact scoped pointer to drop _at the placement section_ of each — not the top.

So the corrected finishing move is three **surgical** edits, each a single blockquote line placed right at the relevant heading:

- `plan.md` → at `## Split table`
- `ph-engine-placement.md` → at its placement table
- `rips-flag-nerve-brief.md` → at `## The discipline — homes by concern`

Each says only "placement moved to `tda-placement.md`; this section kept for its reasoning" — nothing above or below those sections is touched, so all the live content (the why, the floor decisions, the seams, the naming) stays put. That's the difference between _consolidating the placement scatter_ and _retiring three docs_, which is what you were right to stop.

Want me to drop those three scoped pointers in?

> sure

<details><summary>3 previous messages</summary>

> Getting the exact heading text in each doc so the pointers land precisely.
>
> All three headings confirmed. Dropping the scoped pointers in.
>
> Those two files need a Read before editing (I'd only grepped them). Reading the heading regions, then re-applying.

</details>

Done — all three scoped pointers are in, each a single blockquote at the placement section only; nothing above or below is touched, so plan.md's floor reasoning, ph-engine-placement's gate/SPRED seam, and rips-flag-nerve's naming/validation all stay intact.

One substantive thing I caught while placing the rips-flag-nerve pointer: that doc proposes `tda/complexes` and `tda/nerve` as homes — which your "everything PH → `tda/ph`" ruling overrides. So I made its pointer explicit: _"the fold lands the units in `tda/ph`, not a separate `tda/complexes`/`tda/nerve`."_ That's the one place the old doc's specifics actively disagreed with the new authority, now reconciled rather than left to mislead.

Net state of the `persistent-homology/` folder now:

- **`tda-placement.md`** — the placement authority (the one rule).
- **`primitives-to-ph-fold-brief.md`** — the reorg that makes it true (dissolve `primitives`, nerves → `ph/nerves`).
- **`opus-brief-conditioned-persistence-synthesis.md`** + **`p0-conditioned-filtration-brief.md`** — the design + first build step.
- The old placement docs now defer cleanly instead of drifting.

The placement churn is closed. Whenever you want to move from docs to code, the order is fold → P0 (both are executor-lane tickets ready to spawn); until then the rule "new PH code → `tda/ph`" holds. Nothing else needs a decision from you right now.

> well ive got a draft ingestion of this guy `D:\aghado01\codex-scientiae\compendia\ph\2509.20220v2.md` which i think factors into the planning

<details><summary>2 previous messages</summary>

> Yes — that's Wolf–Fan–Monod, _Generalized Persistent Laplacians_, and it's the exact lead I flagged for P4's open term (which Laplacian keeps monotonicity/stability). Let me read it. (Staying primary-agent per your CLAUDE.md, not fanning out.)
>
> That's the load-bearing content — I have what the plan needs (the hook's right that I don't need the other ~150KB of operator-theory proofs; the abstract + contributions + the persistent-Laplacian definition carry the result). Here's how it factors in, and I'll fold it into the brief.
>
> **What Wolf–Fan–Monod settles:** the _classical full_ persistent Laplacian fails monotonicity + stability, but `Δ_k^{K,L} = up + down`, and the **up- and down-components each satisfy both**, and — the key line — **their spectra fully determine the full's**, so the components are "not only preferable but sufficient." Plus a sufficient condition for when the _full_ one is monotone. Generalized to any chain complex of Hilbert spaces (finite **and** infinite-dim, incl. de Rham).
>
> **How it factors into P4 — it de-risks the open term by more than half:**
>
> 1. **It answers "which Laplacian to compute": the up- (and down-) components, separately — not the full.** And `λ_q` should be the **up-component's** gap (the stable, monotone one). This also aligns the whole PL stack on the up-Laplacian: Botnan–Dong's near-linear (non-branching) _and_ Wei–Wei's homotopy-continuation eigenvalues _and_ Wolf–Fan–Monod's stability all target the same object. Consistent directive, three papers.
> 2. **The magnetic/sheaf persistent Laplacian's stability along the filtration `α` is now _inherited_, not open.** A cellular-sheaf cochain complex (ℂ stalks, `e^{iθ}` restriction maps) is a finite-dim **Hilbert complex** with the isometric-inclusion filtration Wei–Wei's PSL already sets up — i.e. _exactly_ Wolf–Fan–Monod's framework. So the up-persistent sheaf Laplacian's `α`-monotonicity + stability come free from their theorem.
> 3. **So the genuinely-open contribution narrows sharply** — from "define the `(α,q)` operator _and_ its stability" to just **the gauge `q` axis**: `q` is not a filtration axis, so filtration-monotonicity says nothing about the `q`-sweep or the joint `(α,q)`. That's the residual novel bit.
>
> Folding it in:

</details>

Folded in — [synthesis brief](../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/opus-brief-conditioned-persistence-synthesis.md) updated in three places: P4 now leads with the "compute up/down components, not the full" directive and the narrowed open term; the status-table row records `α`-stability as _inherited_ with only gauge `q` open; and Wolf–Fan–Monod is in the Inputs.

**Net effect on the plan:** P4's research contribution shrank from _"define the `(α,q)` operator and prove its stability"_ to _"characterize the gauge `q` axis of an already-defined, already-`α`-stable operator."_ Three of your four Laplacian anchors now point the same way — up-component-centric — which is a clean, non-accidental convergence: Botnan–Dong computes it fast, Wei–Wei gives its eigenvalues and the sheaf umbrella, Wolf–Fan–Monod proves it's the stable/monotone one that determines the full.

Two honest flags on what I did and didn't do:

- **The inheritance claim is _my_ bridge, not a cited result.** Neither paper discusses the magnetic/sheaf case in the other's framework — I'm asserting that a finite-dim cellular-sheaf cochain complex with Wei–Wei's isometric-inclusion filtration lands inside Wolf–Fan–Monod's "chain complex of Hilbert spaces" setting, so their up-component stability transfers. It's sound (finite-dim Hilbert complex + isometric inclusions = their base case), but it's the one link to sanity-check before you lean on it — the brief says "exactly Wolf–Fan–Monod's framework," which is a synthesis I'm confident in but that no referee has signed off.
- **I read the abstract/intro/contributions + the persistent-Laplacian definition, not the detailed §Spectral-Properties theorems or the "sufficient condition for full monotonicity."** I have the load-bearing claims with high confidence, but if you ever want the _full_ Laplacian (not just components), that sufficient condition is the thing to have — I can pull just those statements via `ctx_execute_file` without dragging the other 150KB into context.

Want me to extract the precise monotonicity/stability theorem + the full-monotonicity condition (surgical, via ctx), or is the components-are-sufficient takeaway enough to leave it there and keep the gauge-`q` axis as the marked-open term?

> i dont think it would be wise to simply skim material like this and declare victory confidently. lets take a closer look

<details><summary>4 previous messages</summary>

> You're right — that was a skim-and-declare, exactly the thing to not do with a stability result whose whole content is in its hypotheses. Let me actually read the monotonicity/stability section and check my "inheritance" claim against the real theorem statements and their assumptions, not the hand-wave. First, locate the relevant sections.
>
> Good — the structure is clear, and there's a directly relevant tell: they treat **Persistent Cosheaf Laplacians** as an example (line 1212), so sheaf/cosheaf-type Laplacians are explicitly in their framework by their own examples. But I need the actual hypotheses of the monotonicity/stability theorems, not just that a cosheaf example exists. Reading the load-bearing ranges: the persistent-chain-Laplacian definition, the cosheaf example, the monotonicity definitions + theorems, and the stability theorem + full-monotonicity condition.
>
> This changes the picture in an important way — the closer look shows my claim is _directionally_ right but I was sloppy about the load-bearing hypothesis. The containment (Def, line 1013) **requires the inclusions to be isometries** — condition (b), `⟨v,w⟩^P = ⟨ι_k v, ι_k w⟩^Q` — and the filtration (line 1580) is _by_ isometric inclusions. And crucially it's a **cohomology-side** instance (their Example at line 1027 handles cohomology via `ι = π*` with _compatible_ inner products), not the naive chain-inclusion I was hand-waving. So the real question isn't "is it a Hilbert complex" but "does the sheaf filtration have compatible inner products across α, at fixed q." Let me read the actual stability theorem, the full-monotonicity condition, and the cosheaf example to nail it.
>
> Now I can give you a real answer, and it matters: **the closer look confirms the conclusion but shows my _mechanism_ was wrong.** Not "it's a Hilbert complex, therefore covered" — that was the skim talking. The actual situation:
>
> **What the theorems actually say (and where my hand-wave was load-bearing):**
>
> 1. **Containment requires isometry.** WFM's `P ⊂ Q` (line 1013) demands condition (b): `⟨v,w⟩^P = ⟨ι_k v, ι_k w⟩^Q` — the inclusions must be **isometries** — and a filtration (line 1580) is _by_ such isometric inclusions. That hypothesis is the whole game; I glossed it.
> 2. **But the sheaf case is covered — more explicitly than I claimed.** WFM devote a section to **Persistent Cosheaf Laplacians** (line 1212) and prove (line 1257) that **any Hermitian PSD operator is realizable as one**. `L^(q)` is Hermitian PSD, so it's genuinely in their framework — the _right_ justification, which I didn't have.
> 3. **Up/down stability (Thm 25) is real** — `d_I(λ_±, λ_±) ≤ d_I(P,Q)` — but strictly w.r.t. the **filtration** interleaving distance.
> 4. **The full Laplacian is definitively out**: they give a counterexample where it fails monotonicity _and_ stability **even for `q=1` simplicial complexes** (line 2188–2279). So "compute up/down components, not the full" is **necessary, not tidy**.
>
> **The verdict on my claim:** it holds, but **contingent on a hypothesis I have to state, not assume** — the `α`-filtration of magnetic-sheaf cochain complexes must be by isometric inclusions, which for the magnetic case means **fixed `q` with a filtration-consistent connection** (the phase fixed on the ambient complex, restricted to sub-complexes). The magnetic operator's phase-only/unit-magnitude/standard-inner-product setup satisfies this cleanly; a general _weighted_ sheaf would additionally need weight-consistency. So `α`-stability is inherited **under that condition**, and the gauge `q` remains genuinely open. Correcting the brief to say that rather than hand-wave:

</details>

You were right to stop me. Here's what the closer look actually changed, plainly:

**My skim gave the wrong _mechanism_.** I said "a cellular-sheaf cochain complex is a Hilbert complex, so it's covered" — hand-waving past the hypothesis that carries the whole result. The theorems don't care that it's a Hilbert complex; they care that the **filtration is by isometric inclusions** (containment condition (b)). Had that propagated into the plan, it would have been a confident non-sequitur.

**The conclusion survives — but only under a hypothesis I'd have missed, and via a better justification I didn't have:**

- The _right_ reason the magnetic case is in-scope: WFM explicitly treat **persistent (co)sheaf Laplacians** and prove **any Hermitian-PSD operator is one**. `L^(q)` is Hermitian PSD — so it's genuinely covered, cleanly.
- But `α`-stability is inherited **contingent on the `α`-filtration being isometric**, which for the magnetic operator means **fixed `q` + a filtration-consistent connection** (phase fixed on the ambient complex, restricted down). The built operator's phase-only/unit-magnitude/standard-inner-product setup _does_ satisfy this — but a general _weighted_ sheaf would need weight-consistency too. That's a real condition, not automatic.

**Two things came out stronger than the skim, both usefully:**

- "Compute up/down components, not the full" went from _preferable_ to _necessary_ — WFM give a counterexample where the full Laplacian fails monotonicity **and** stability _even for `q=1` ordinary simplicial complexes_.
- The gauge `q` being open is now precise, not vibes: Thm 25's stability is w.r.t. the _filtration_ interleaving distance, and `q` isn't a filtration parameter, so the theorem provably says nothing about the `q`-sweep.

**Net for the plan:** the open work now cleanly splits into two tiers — (a) a _checkable condition_: confirm the built `MagneticLaplacianOperator`'s `α`-filtration gives isometric inclusions at fixed `q` (a bounded math/engineering check, not research); and (b) the _genuinely open_ bit: the gauge-`q` axis and joint `(α,q)`. The brief now states the hypothesis and both tiers instead of asserting inheritance for free.

That's the difference the closer read bought: same headline, but with the load-bearing condition surfaced rather than buried — and it's the kind of thing that, left as a skim, quietly rots a plan from the inside.
