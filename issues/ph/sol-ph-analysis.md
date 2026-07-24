# So I want to continue trying to shore up outstanding work a…

> So I want to continue trying to shore up outstanding work and current state of my persistent homology engine from before i did my migration to ThermoMapper. There's definitely the conditioned backbone concept that's still taking shape `D:\aghado01\MarkBrain\ThermoMapper\issues\tda-purification\persistent-homology\opus-brief-conditioned-persistence-synthesis.md`, there's zigzag frontier `D:\aghado01\MarkBrain\ThermoMapper\issues\tda-purification\zigzag-engine\zigzag-frontier.md` which I believe includes the void boundary work outstanding `D:\aghado01\MarkBrain\ThermoMapper\issues\tda-purification\persistent-homology\void-inference-three-engines.md` . There were some loose ends with the Gudhi assimilation effort though I can't find the latest report on that atm. then there's the emerging persistent laplacian angle and complex analysis inspired concepts that very much interest me, which have some ground broken on `D:\aghado01\MarkBrain\ThermoMapper\issues\tda-purification\persistent-homology\Gemini-Complex-Analysis-PH-ideas.md` "magnetic laplacian" `D:\aghado01\MarkBrain\ThermoMapper\issues\tda-purification\persistent-homology\gemini-magnetic-laplacian.md`

<details><summary>3 previous messages</summary>

> Good — this is a "re-orient on the PH engine corpus and map outstanding work" pass. Per your memory, zigzag-frontier.md is the canonical cross-track status doc and per-rung briefs lag committed code, so I'll read the corpus first, hunt down the missing Gudhi assimilation report, then verify against the actual repo state. Reading the named docs now (sequentially in-session, no sub-agents).
>
> The frontier doc loaded fine but Glob came back empty — likely a path-resolution quirk. Let me read the named docs directly and list the directory another way.
>
> You've hit your session limit · resets 12:30am (America/Los_Angeles)

</details>

No response requested.

> picking up where we left off

<details><summary>2 previous messages</summary>

> Resuming exactly where we stopped: HEAD confirmed as `513c4bd`, all verification done. Two steps left — one honesty-check grep (does the strong-collapse validation actually ride `ZigzagMapBarcode`?), then truing up the frontier dashboard and delivering the state map.
>
> Table rows flipped. Now the contention map and the frontier list:

</details>

The frontier dashboard is now current. Here's the full picture, verified against code rather than briefs.

## State of the PH engine, verified (HEAD `513c4bd`, tests green: TDA.Ph 293/293, TDA.Mapper 28/28)

**The Gudhi report you couldn't find doesn't exist as a separate file** — I grepped all of MarkBrain and the repo. The latest Gudhi-assimilation status lives inside [z3-reflection-zigzag-brief.md](../MarkBrain/ThermoMapper/issues/tda-purification/zigzag-engine/z3-reflection-zigzag-brief.md) as dated status sections (2026-06-22), and the code confirms them: **stages 1, 2a (`RuVineyard`), and 2b-i (reverse-teardown backward arrows) all landed and green.** What remains of the port is **2b-ii** (re-seat the backward arrow on `RuVineyard` streaming append/transpose/truncate so non-last deletes stop throwing) and **stage 3** (the churn 5th corner). The frontier table had it worse than reality ("stage 1 only") while the contention map had it better ("landing stages 2–3") — both now corrected in [zigzag-frontier.md](../MarkBrain/ThermoMapper/issues/tda-purification/zigzag-engine/zigzag-frontier.md), along with a migration note: the squash means every pre-migration hash and track prefix cited in your briefs no longer resolves in git.

**The migration commit's "new work on PH" turns out to be the consumers rung starting** — this was undocumented anywhere: [PersistentMapper.cs](src/tda/mapper/PersistentMapper.cs) (Mapper-per-T → `NerveFiltration` → barcode), [NerveDiff.cs](src/tda/primitives/NerveDiff.cs) (birth/death/merge/split events between nerve frames), and [ZigzagMapBarcode.cs](src/tda/ph/ZigzagMapBarcode.cs) (a zigzag oracle over _arbitrary cell maps_, built to validate strong-collapse core-assembly). All three have green tests. One caveat I confirmed: `ZigzagMapBarcode` is only exercised by its own tests — the strong-collapse tests don't consume it yet, so the Z4 scope confirmation is _enabled_, not _done_. Also note `PersistentMapper` tracks nerve topology frame-to-frame rather than going through the zigzag oracles — seating it on `ZigzagFiltration` proper is still open.

**Conditioned persistence:** the synthesis brief is the authoritative head of that whole thread — it subsumes both Gemini docs you linked (those are raw seeds, already distilled through the two opus briefs in `mapper-integrations/`). Its code-status table checks out against the tree: the Φ engine, `MagneticLaplacianOperator`, and `Lobpcg` are all built as claimed; **none of P0–P4 has started** — no `ConditionedFiltration`, no directed-flag build flag, no sheaf stack.

**Void inference:** pure concept-stage — the three-engines fold is settled on paper (population-MCMC is "the one to build"), but the GA corpus is still staged pre-ingestion and the §7 open decisions (sampler/optimizer/both, orientation co-evolution, code home) are unresolved. No code anywhere.

## The outstanding-work menu, roughly ordered by readiness

1. **ReflectionZigzag 2b-ii + stage 3** — closes the third oracle; the mule-ready spec exists and `RuVineyard` is validated in isolation. The most shovel-ready item.
2. **Wire `ZigzagMapBarcode` into strong-collapse validation** — small, closes the Z4 question.
3. **P0 `ConditionedFiltration`** — opens the whole conditioned-persistence arc, and rides entirely on built primitives (the brief's sequencing discipline: ship Φ, earn Ψ).
4. **Seat `PersistentMapper` on the zigzag engine** — connects the new consumer work to the oracles it was built to feed.
5. **Void engines** — furthest from code; needs corpus ingestion and the §7 decisions first.

One housekeeping note: the vault note behind `project_zigzag_engine_z5` still says "ReflectionZigzag 3rd-oracle stage-1" — since the memory index is sync-generated from MarkBrain, that's an edit to the vault note plus a sync re-run (I didn't touch it). Happy to dig into whichever thread you want to advance first.

> ok there's also the broader context of the TDA purification work and reorganizing persistent homology related concepts to be more first class. originally PH was being built more like an appendage of mapper, and then work started on things like migrating topology primitives to maths.topology, tda/ph organization, etc. With the persistent laplacian and sheaf work, it feels like we kind of forking conceptually between persistent homology proper and related graph persistent laplacian/sheaf concepts. even if the latter isn't very mature yet i feel its worth thinking through what the target organization emerging is

<details><summary>1 previous message</summary>

> Good — this is the right moment to think that through, and there's actually a placement corpus from the purification arc I should read before opining (per your own read-the-corpus discipline). Let me ground in what's settled: the placement briefs, plus what actually lives in `maths/topology` and `graphs` today.

</details>

Read the placement corpus — [tda-placement.md](../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/tda-placement.md) is the authoritative ladder, and it turns out the codebase has already started answering your question in a way worth naming. Here's my read.

## Where the corpus already stands

The settled ladder is `maths/topology → tda/ph → tda/mapper → tda/pipelines`, everything allowed to lean downward on `graphs`/`maths`. The floor today holds exactly three files (`Simplex`, `CombinatorialIndex`, `Z2LinearAlgebra`), and the discriminator is domain-generality, not mere purity. Two pieces of pending mechanical debt matter here: the `tda/primitives` → `tda/ph` fold is specced but not executed (the directory and csproj still exist), and — this is the telling one — `AbstractCell`/`CellFiltration` physically sit in `tda/ph` but declare the `Maths.Topology` **namespace**. That drift isn't sloppiness; it's the architecture trying to tell you something.

## The seam, named precisely

I'd argue the fork you're sensing is real, but it's neither "PH vs graphs" nor "homology vs Laplacian." The honest cut is **coefficient structure**: what you decorate the incidence data with.

- **The barcode reading**: chain complexes over **Z₂** — exact, interval-decomposable, combinatorial. Value type: `Barcode`. This is `tda/ph` as it exists — hardened, zero eigensolver dependencies.
- **The spectral reading**: cochain complexes over **ℝ/ℂ with an inner product and a connection** (stalks + restriction maps — the sheaf; magnetic is the U(1) instance). Value type: spectra / `λ_q` surfaces. Kernel dimensions recover the barcode's Betti numbers, but the machinery, the dependencies (`maths/linalg`, complex arithmetic), and the maturity are entirely different.

Both are _readers over the same currency_ — the filtered complex. That's the same architectural move as `project-thermomapper-architecture` (Potts as model, inference strategies as consumers): the filtration is the model; barcode reduction and persistent Laplacians are two observables on it. And it's the structural expression of your Hodge-bridge frame — the spectral wing literally _is_ the bridge layer, so its position between `graphs/spectral` and `tda/ph` is not an accident of the Li/Wei–Wei reading.

## The target picture

```
maths/topology      incidence floor + THE CURRENCY (cells, complexes, filtration order)
maths/linalg        eigensolvers (Lobpcg; later FEAST)
graphs/spectral     operators on a FIXED graph — MagneticLaplacianOperator, R_eff (unchanged)
tda/ph              reader 1: Z₂ chain reading → barcodes (reducers, zigzag, nerves)
tda/spectral (new)  reader 2: inner-product cochain reading → sheaf cochains,
                    sheaf/persistent Laplacians (up/down), λ_q observables
tda/mapper          cover production + orchestration (unchanged)
tda/pipelines       apps (unchanged)
```

Three consequences fall out:

**1. The currency floors.** The corpus's own second-consumer rule, applied to the filtration currency itself: `tda-placement` correctly kept `Barcode` in `tda/ph` because it's PH-domain — but the moment a _second engine_ consumes filtrations, `AbstractCell`/`CellFiltration`/the filtration order stop being PH-domain and become domain-general. The prophetic namespace gets its file move. (`Barcode` and `DiagramMetrics` stay in `ph` — they're barcode-reading-specific.)

**2. `graphs/spectral` stays the fixed-graph operator shop.** The precedent is already in `ph-engine-placement`: R_eff metric → `graphs/spectral`, the Rips adapter consuming it → `tda`. Same seam one level up: the magnetic operator stays put as the fast 1-skeleton instance; `tda/spectral` owns everything filtration-indexed, pair-indexed (Δ^{K,L}), and q≥2. Dependency is legally downward.

**3. The oracle seam lives in tests.** `nullity(Δ_q^{K,L}) = persistent Betti` is a cross-engine validation identity — exactly your multi-oracle discipline extended across the fork. Test-level coupling, no src-level dependency between the readers.

## The genuine forks still open, and my leanings

**Sibling vs wing** (`tda/spectral` vs `tda/ph/spectral`). The wing keeps one engine and matches the discriminator's "any topology-engine machinery → ph" reading. But I lean **sibling**, for an engineering reason rather than a conceptual one: `TDA.Ph` today needs no eigensolver, and the SPRED/gate hot loops ride on that lean Φ core. Dragging `Maths.LinAlg` and complex arithmetic into the same project to accommodate a research-grade P3–P4 stack couples maturity profiles that want to stay decoupled. The pressure-test against me: the doc once killed a fake tier ("there is no graph-coupled builder tier"), and I'm proposing a new box — but a sibling _within_ the ladder's engine rank isn't a new tier, it's a second occupant of an existing one.

**Naming.** Survey: `tda/spectral` (parallel to `graphs/spectral`; "spectral TDA" / Wei's "persistent spectral theory" are the actual field terms), `tda/laplacians` (concrete, but narrower than the sheaf content), `tda/sheaves` (the true umbrella per Wei–Wei, but the near-term deliverables are Laplacian-flavored). I'd take `spectral` — field-idiomatic and the parallelism with `graphs/spectral` encodes the fixed-object/filtered-object split for free. Per the synthesis brief's own discipline, Φ/Ψ stay in prose.

**Timing.** Nothing forces the directory into existence today — P2-spectral (flux, M1–M3) consumes the graph-level magnetic operator directly, and the sheaf stack only materializes at P3. What I'd do now is settle the _decision_: add the spectral row + currency-floors clarification to `tda-placement.md` (it's the placement authority; this is exactly the kind of drift it exists to prevent), and sequence the `primitives → ph` fold **before** the new wing lands so the fold isn't chasing a moving target.

The part I'd most want your push-back on is the sibling-vs-wing call — if you see the oracle fabric as wanting source-level proximity rather than test-level, the wing reading gets stronger.

> ok, i'll buy the wing fork. lets try to capture this in one of the synthesis docs to centralize holistic design thinking. There's still more i want to ask about and some convergent lines of work from outside of the PH docs after that, such as BP-CSR integration into graph compiler internals and the connection to zigzag.

> [Request interrupted by user]
