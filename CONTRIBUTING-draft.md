# CONTRIBUTING.md — Design Principles & Engineering Standards

Normative engineering rules for ThermoMapper: architecture invariants, API shape, naming, refactoring discipline, docstrings, and verification standards. Repository layout, toolchain invocations, runtime environment rules, and `issues/` coordination live in [AGENTS.md](AGENTS.md).

## 1. Engine-First Strategy

Build and validate engine primitives to mathematical and numerical soundness **before** wiring downstream applications or UI consumers.

- State capabilities by problem class and evaluate methods and work sequencing strictly on engine grounds (seam generality, cost, correctness, external validation) to avoid consumer-coupled rebuilds. Cite specific applications only as examples.
- Validate in-progress engine work against external oracles, brute-force references, and property checks rather than wiring a consumer to test it.
- Sequence engine depth (exact bounds, paper-faithful extensions, reference twins) ahead of integration. Incorporate corrections applications expose once consumers exist.

## 2. Placement & Acyclic Architecture

Repository layout and concrete dependency layers are catalogued in [docs/Architecture.md](docs/Architecture.md). Placement follows strict acyclic layering:

- **Downward-only dependencies.** Modules sit by abstraction tier (foundations and neutral substrates below consumers and pipelines). Dependency edges form a DAG pointing strictly downward toward foundations; avoid lateral or upward imports that introduce circular dependencies.
- **Substrate decoupling (dependency inversion across seams).** Lower-tier substrates never import higher-tier consumers. When a substrate stage requires criteria computed by a consumer or upper-tier engine, compute the criterion above the seam and inject it downward as typed data or a delegate contract.
- **Consumer-neutral artifacts (litmus test).** Substrate outputs and serialized artifacts must remain self-contained so any unrelated downstream consumer can read them unchanged without consumer-specific hooks. If a consumer requires construction changes that leak its concerns into the substrate, re-cut the seam rather than threading consumer types downward.

## 3. Strict Core, Fluent Shell

Recurring implementation patterns are documented in [docs/Design.md](docs/Design.md). For multi-stage pipelines with interdependent configuration (overrides, auto-selections, fallbacks, conditional skips):

- **Declarative core.** Backends expose stateless static entry points (`public static XxxResult Build(XxxConfig config, …)`) taking one deeply-immutable, JSON-serializable `XxxConfig` record. The engine resolves auto-selections and fallbacks holistically before executing.
- **Config shape.** Use `init`-only properties (`required` only for non-defaultable inputs; nullable or defaulted otherwise) and nested per-stage sub-configs pairing a strategy enum (`Auto`, `Explicit`, …) with its parameters. Live objects and delegates ride the method call, never the config DTO.
- **Provenance manifests.** Emit a manifest snapshotting the resolved config, diagnostic reports, and automation decisions (which auto-picks fired and why).
- **Shell-only builders.** Fluent or chained ergonomics live exclusively at the CLI/REPL boundary, validating and normalizing input into clean config DTOs before invoking the engine. Never put fluent builders in backend engines.
- **Complexity trigger.** When a static entry point accumulates ≥6 parameters encoding discrete pipeline stages or strategy selectors, replace the parameter list with a config record.

## 4. Naming, Shape & Documentation

### Naming & Vocabulary
- **Self-describing domain names.** Name types, commands, scripts, namespaces, and assemblies for what they are or do using established terminology from the owning field rather than project prefixes or invented labels. A domain reader unfamiliar with the brief must be able to infer the role from the name alone.
- **Cross-field precision.** Borrow exact terms from the fittest discipline without anchoring to a single field. A borrowed term fails if it is intuitive-but-narrowing (imports unwanted constraints) or general-but-inert (carries no structural weight).
- **Reason abstractly, reify concretely.** Use general mathematical abstractions while designing, then implement concrete domain types rather than generic mathematical interface wrappers.
- **Module-consistent notation.** Follow the module's existing symbol conventions regardless of source-paper notation. Scope notation fixes to the text under edit; log unresolved splits in [issues/notation-ledger.md](issues/notation-ledger.md).
- **One unified shape.** Render one simplest unified structure across sketched patterns and sibling modules. Do not introduce asymmetric layouts for incidental implementation differences; surface genuine misfits as explicit design questions.

### Docstrings & Canonical Documentation
- **Terse operational docstrings.** State what the symbol does, attribution, pertinent behavioral invariants, and one concise clause locating its role in the module's pipeline with `<see cref="…"/>` links to referenced concrete types.
- **No persuasive or archival prose in code.** Docstrings must not contain manifesto language, academic critiques, full derivations, or prior-art reviews. Park derivations and critiques in `issues/` design notes.
- **Documentation ownership & trade-offs.** `docs/` holds current canonical specifications rewritten in place using repo-relative links without personal context, machine-specific paths, or chronological logs, and is updated in curated offline batches (`docs/` is currently frozen during cleanup). Record settled statistical or design choices—including alternatives considered, trade-offs, and user implications—under **Doc changes owed** in `issues/` reports while keeping decision provenance in the governing brief.

## 5. Refactoring & Change Discipline

- **Pre-release clean breaks.** While pre-release, nothing outside the repository depends on current APIs (revisit once external consumers exist). Do not add compatibility shims, `[Obsolete]` markers, re-export stubs, or schema-version migration branches; delete superseded APIs, schemas, and single-caller helpers, and migrate all callers.
- **Staged multi-step refactors.** Keep the build green mid-flight and retain old code for comparison: write replacement files directly in their permanent locations alongside old ones, migrate callers incrementally, and isolate all file deletions into a single dedicated cleanup commit.
- **Sibling sweep.** When fixing a bug, API misuse, or invariant violation in one location, search the codebase for sibling occurrences and fix them in the same change.
- **No hard-coded paths.** Keep user-specific or machine-specific paths out of source files, tests, and scripts. Confine generated artifacts to `artifacts/`.

## 6. Verification & Testing Standards

- **Oracles and source skepticism.** Validate numerical and algorithmic correctness against external reference oracles, brute-force references, and property checks. External papers and reference implementations are checked, not trusted: an oracle is a parity target only where it is right. Correct proven flaws (incoherent definitions, broken intent, silent model changes, or dominated methods) rather than preserving them as options; retain regime-dependent trade-offs with measured defaults, and record each correction with its evidence in the module's source-changes register and governing brief.
- **Independent ground truth.** Validation ground truth must originate outside the estimator's own model (real metrics and data, or independently derived references). Tests whose generator shares the estimator's forward model check inversion only.
- **Tiered verification & gates.** Routine verification targets the change's compile-closure selection (excluding slow tests) on edits, adds stored-expectation pinned tests at commits (acting as the net for reflection and non-lexical reachability; same-build consistency checks stay unmarked), and runs the full assembly plus slow tests at slice completion alongside declared qualification gates. Runner invocations live in [AGENTS.md](AGENTS.md).
- **RNG hygiene.** Derive child seed streams structurally (seed trees / splitting), never via ad-hoc arithmetic (`seed + k`, `seed ^ C`) that aliases streams. Use pinned generators over platform defaults where cross-version fixture stability matters, and treat pinned-fixture expectation regeneration as an explicit, documented act.
- **Sweeps & comparative assertions.** Sweep numerical primitives across parameter ranges and orders of magnitude against independent properties under a hang timeout. Evaluate stochastic superiority claims across both sampler-seed (chain) and data-seed axes: replace assertions that fail frequently on the chain axis, qualify dataset-dependent claims, and rerun sweeps whenever trajectories change.
- **Load-invariant benchmarking & resumable harnesses.** On shared hosts, evaluate performance using allocated bytes and operation counts alongside wall-clock times (pairing blocks across clean snapshot builds). Long-running harnesses append each completed result incrementally and skip finished rows on rerun.

## 7. Antipatterns

- Backwards-compatibility shims, `[Obsolete]` aliases, or schema-version migration branches (pre-release).
- Fluent builders or sequence-dependent state in backend engines.
- Wiring downstream consumers to validate an in-progress engine.
- Upward or circular dependencies between substrates and consumers.
- Deleting superseded files mid-refactor instead of in the final cleanup commit.
- Partial fixes that leave identical sibling instances unfixed.
- Project-prefixed or invented names instead of standard domain terms; asymmetric pattern renderings.
- Persuasive prose, academic critiques, or long derivations in docstrings.
- Arithmetic seed derivation (`seed + k`, `seed ^ C`) or unpinned RNG in stable fixtures.
- Hard-coded local paths or artifact writes outside `artifacts/`.
