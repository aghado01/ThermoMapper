# Docs owed: proposed additions (draft)

Compiled 2026-11-12 from two checks: `AGENTS-draft.md` and `CONTRIBUTING-draft.md` against the current `AGENTS.md` and `CONTRIBUTING.md`, and the live memory cards against both drafts. Each bullet gives proposed wording, then its source.

- `[cur]`: present in the current docs, absent from the drafts, with no other home found.
- `[mem: name]`: a memory card.

**Sorting rule.** A line goes to AGENTS or CONTRIBUTING when it can be stated without naming a type, file, dataset or incident and holds in every area. Anything else goes beside the thing it describes.

## AGENTS-draft.md

### §2 Environment

- Python work runs through `scripts/env/portable-python.ps1` and `scripts/env/venv-boostrap.ps1`. `[cur]` (the current doc gives the path as `scripts/`; the files are under `scripts/env/`)

### §3 Toolchain

- Build row: in a shared checkout each session builds into its own root, `-p:ArtifactsRoot=<repo>\artifacts\<session>\`, on every `dotnet build` and `dotnet test`. `[mem: artifacts-root-chains]`
- Lean row in §1: `Theorems` is an optional promotion for consequential verified results, and `Archeion` is the non-building exit for superseded material. The draft reads as one linear chain. `[cur]`, matches the header of `lean/scripts/meta-ci.ps1`

### §4 Execution, repair and commit

- Commit in one step by pathspec, `git commit -m "…" -- <paths>`, and leave nothing staged between commands. Sessions share one git index. `[mem: shared-index-commits]`
- A step ends with a review pass: reread the whole diff as a first-time reviewer, look for options and refusal paths no test runs, fix, and repeat until a pass finds nothing. The report records what each pass found. `[mem: self-review-cycle]`
- Removing an option, rule or output is put to the user as its own question, stating what is lost and what replaces it. A brief that says "check X" covers reporting on X. `[mem: removals-own-question]`
- Heavy or long runs yield the machine: below-normal priority and capped workers. Long harnesses append each result as it ends and skip finished rows on a rerun. `[mem: shared-machine-load]`, `[cur]` for the test-select settings

### §5 DocOps and coordination

- Chips message each other directly when they need to coordinate. `[cur]`
- A hand-over entry states what is uncommitted, what is running, and what comes next. `[cur]`
- Pages under `issues/` reference each other by wikilink or explicit relative link. `[cur]`

## CONTRIBUTING-draft.md

### §1 Engine first

- The engine serves many applications. State a capability by the problem class it covers, and judge methods and order of work on engine grounds (seam generality, cost, correctness, external validation). A named application appears as one example among several. `[mem: application-agnostic-engine]`
- One-line reason: a primitive shaped by its first consumer gets rebuilt when the second one needs it exposed differently. `[cur]`
- Circle back on the corrections applications expose once they exist. `[cur]`

### §4 Naming, shape and documentation

- Notation follows the module: a concept keeps the symbol the module already uses, whatever the source paper writes. A notation fix scopes to the text being fixed, and a split goes to the notation ledger. `[mem: module-notation]`
- A name states what the thing is or does, in the term its field already uses for that object or operation: `H1CycleEdges` holds the edges of H1 cycles, `repo-audit` audits the repo. The test: a reader who knows the field and has not read the brief can say what it is. `[mem: artificial-memes]` (would absorb the draft's "Functional names over project prefixes" bullet)
- The fitness test fails in both directions: a term that is intuitive but narrowing, and a term that is general but inert. `[cur]` (the draft keeps the first half)
- The locating clause in a docstring is a clause, with one example to copy: "accumulates the per-draw bond events and reduces them to the `Affinities` currency." `[cur]`
- A decided statistical or design choice is stated on its docs page with the alternatives considered, the trade-offs, and the implications for users. The brief keeps the provenance. `[mem: document-tradeoffs]`

### §5 Refactoring

- One-line reason for clean breaks: nothing outside the repo depends on the current shape, so a shim is cost with no benefit. Revisit when external consumers exist. `[cur]`
- One-line reason for deferred deletes: the build stays green mid-flight, the old code stays available for comparison, and "build the new shape" and "delete the old shape" land as separate commits. `[cur]`

### §6 Verification

- Validation ground truth originates outside the estimator's own model: the real metric and data construction, or an independently derived reference. A test whose generator shares the estimator's forward model checks inversion only. `[mem: validation-independence]`
- An oracle is a parity target only where it is right. `[cur]`

## Other docs

### docs/Architecture.md

- Dependency layers: the topology stack runs `maths/topology → tda/ph → tda/mapper → tda/pipelines`, and all persistence code lives in `tda/ph`. `[cur]` (the current doc also points at `tda-placement.md` in MarkBrain as the placement authority)
- Cross-module contracts, beside the existing `ProtectedEdgeSource` row: `GraphCompiler.Build(config, n, metric, ProtectedEdgeSource? protectedEdges = null)` illustrates substrate decoupling—the caller (`UserRepl.Commands` / upper layer) computes topology-scored criteria above the seam and passes them into graph construction (`src/graphs/**`, referencing only `graphs` and `maths`) as data/delegates. Later scored passes (persistence or spectral criteria) enter the same way. `[cur]`, moved from `CONTRIBUTING-draft.md`
- Substrate litmus check (`CsrGraph`): `CsrGraph` round-trips to disk via `WriteTo`/`FromBinary` alongside its provenance manifest; construction stays decoupled when non-SPC readers (pure PH, Mapper nerve, spectral embedding) consume the output unchanged without SPC-specific construction hooks. `[cur]`, moved from `AGENTS.md`
- Responsibilities by area, `Clustering.Evaluation` row or `Clustering/ClusterEvaluation.md`: evaluation holds indices of a result's quality, external against a reference and internal from the data or graph. A functional that compares any two partitions is an observable and lives in `Graphs.Observables`. `[mem: evaluation-is-quality]`

### docs/Clustering/SpcEngine.md (and/or docs/Concepts.md)

- Function-algebra vocabulary & docstring locating clause: weave the SPC/Swendsen–Wang vocabulary (*field / accumulate / reduce / observable / currency / sampler / solver / affinity / form-degree*) as one clause locating the type in the recurring shape (e.g., *"accumulates the per-draw bond events and reduces them to the `Affinities` currency"*) with `<see cref="…"/>` links. Reason in general terms (functionals, maps), reify concrete domain types (`IReduction`, not `IFunctional`). `[cur]`, moved from `CONTRIBUTING-draft.md`

### docs/Design.md

- New principle, retained generative state: each draw keeps its structure, parameters, noise scale and weights. Derived curves are recomputed from that state, and streaming accumulators are online views of it. A seam returns what it computed. `[mem: retain-dont-discard]`
- New principle, calibration as a selected value: the API exposes the capability space established from first principles, and a choice that came from matching one dataset or paper is one selectable point in it. The method owns its calibration procedures and consumers call them. `[mem: calibration-not-api]` (the card names `SusceptibilityKind` and `--bandwidth-strategy` as existing examples; not re-verified)
- §1 Declarative config: why chained ergonomics stay at the CLI and REPL (`UserRepl.Commands`). A chain implies sequence, while the engine evaluates the config as a whole (LMP rescaling makes the global bandwidth dead). The shell validates and warns before the engine runs. `[cur]`
- §6 Reproducible randomness: concrete RNG/sweep types moved from `CONTRIBUTING-draft.md`—derive child streams via `SeedTree` / `SplitMix64`, prefer `Xoshiro256PlusPlus` over `System.Random` for fixture stability, and register stochastic claims in `ComparativeSweep` (`SEED_SWEEP=<n>`). A cache whose path differs in the low bits from the uncached path is switched by a static rule (size, configuration), so a restored or handed-over run continues bit for bit. `[cur]` + `[mem: pure-caches]`
- §10 Capability = nullability: when a payload makes a path expensive, the path that does not read it opts out (`maxCycleDimension`) and the output itself stays. `[mem: cost-on-a-payload]`

### docs/Graphs/GraphConstruction.md

- Experiments and tests state each graph as `GraphBuildResult.Resolved` reports it: topology and K, filter, repair and its edge count, what `Auto` became, kernel, bandwidth rule and value. Findings are scoped to the graphs measured, and a test names and pins the graph it asserts on. `[mem: graph-settings-explicit]`

### A development and testing page (no home yet)

Nothing under `docs/` mentions `test-select`, `SEED_SWEEP` or `ArtifactsRoot` today. The current CONTRIBUTING is the only description besides the script headers.

- How selection follows names: a changed file reaches every file in the compile closure that names one of its top-level types, to the fixed point. A csproj or `Directory.Build` change selects every class. `[cur]`
- Cadence: `Edit` scope at each edit, `Commit` scope at each commit, `Slice` scope once at a slice's end, qualification gates as each slice states. `[cur]`
- Run settings: Release, below-normal priority, xUnit capped at 6 threads, build servers off. `[cur]`
- Tiers: a new pin takes the `Pinned` trait. The slow list holds each method that summed to 10 s or more in a whole-assembly profile, and a new profile rewrites it. `[cur]`
- Comparative assertions: an assertion that fails often on the chain axis is replaced. `[cur]` (the draft keeps the rest)
- Numerical primitives get a seeded sweep over parameters and intervals across many orders of magnitude, checked against an independent property and kept as a test (`tests/distributions/DistributionSweepTests.cs`). New far-range cases run under a hang timeout. `[mem: parameter-sweep]`
- Measuring on a busy machine: read allocated bytes beside times, pair blocks across clean `git archive` builds, and treat counts as the part that stands under load. An R-oracle `TimeoutException` is rerun at normal priority before it is read as a failure. `[mem: shared-machine-load]`

### Shell and build quirks (development page, or CLAUDE.md)

`CLAUDE.md` reserves itself for behaviors no other agent shares. These apply to any agent on this machine, so the development page fits better.

- Run `dotnet` with `-p:ArtifactsRoot=` through PowerShell with the property quoted. Git Bash drops the backslashes and the build writes stray `bin/` and `obj/` beside the csproj. `[mem: artifacts-root-chains]`
- Line endings are mixed, and the `docs/` pages are CRLF. `sed -i` and Python universal-newline writes convert the whole file, so edit in place and check `git diff --stat`. `[mem: script-edits-line-endings]`
- A file restored with `Copy-Item` keeps the backup's timestamp, so MSBuild treats the stale DLL as current. Touch the file after restoring. `[mem: restore-keeps-old-mtime]`

## Fixes owed (not additions)

- `docs/Architecture.md`: 23 relative page links are broken. The pages moved into subfolders, `./DesignPrinciples.md` is now `./Design.md`, and `./BayesianRegression.md` is archived. The `Maths.Regression` row predates the `Regression/` folder.
- `docs/Design.md`: at least 17 relative page links are broken for the same reason.
- *(Resolved in drafts)* Both `CONTRIBUTING-draft.md` and `AGENTS-draft.md` now state the same `docs/` ownership rule (`docs/` is currently frozen during cleanup; sessions log updates under `Doc changes owed`, applied in curated batches by rewriting in place).
- *(Resolved in drafts)* `CONTRIBUTING-draft.md` §6 no longer links into `docs/Maths/Regression/` or `issues/regression/`; it states the general source-correction rule instead. (Note: the live `CONTRIBUTING.md` link to `docs/Maths/BayesianRegression.md` remains dead until promoted.)

## Considered and left out

- Probes need a decision, scrub preserves intent, names follow methods: guidance on how an assistant works with you. They stay in memory.
- Coherent upgrades: you declined a CONTRIBUTING entry for it on 2026-11-02.
- Constants derivation home: the card names `src/maths/constants/`, which does not exist.
- Automated discovery (trim safety for reflection-based discovery): narrow, and unassessed against the code.
- Notation ledger timing: already stated in the header of `issues/notation-ledger.md`.
- Rules new in the drafts with no source in the current docs, the repo's practice or memory: the three-cycle repair cap, the closeout template, the sibling sweep, "do not claim enforcement by a checker that does not exist". They are yours to keep or cut.
