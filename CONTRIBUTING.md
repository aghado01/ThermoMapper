# CONTRIBUTING.md — development practices & design principles

Source of truth for *how* development work is done here: strategy, design
principles, conventions, antipatterns, and testing standards. Project
orientation, layering invariants, tooling, and agent-facing behavioral
guidance live in [AGENTS.md](AGENTS.md).

## Strategy: engine first

Build the engine and its primitives from first principles to soundness
**before** building applications or consumers. Circle back on corrections the
applications expose once they exist.

**Why:** locking primitives to a premature consumer is a known failure mode —
it has repeatedly forced painful rebuilding once another consumer needed the
primitives exposed differently or more completely. The underlying theory is
developed in parallel; sound bedrock makes applications proposable without
speculation.

**In practice:**

- Mid-engine, do not wire up a consumer "to validate" — validate with oracles,
  brute-force references, and property checks (the engine's own correctness
  spec; see Testing standards below).
- Sequence further engine depth (representatives, paper-faithful extensions,
  exact bounds) ahead of consumers. Integration goals get surfaced when
  raised, then deferred.

## Architecture: strict core, fluent shell

Multi-stage pipeline backends (graph construction, clustering, mapper, …) with
interdependent logic — interrupts, overrides, auto-picks, conditional skips —
keep the **backend purely declarative**: one deeply-immutable `XxxConfig`
record passed to a `Build(config, …)` entry point. The engine evaluates the
whole config holistically, resolves auto-picks and fallbacks internally, then
executes.

**Never put a fluent builder in the backend.** Fluent/chained ergonomics live
at the CLI/REPL boundary (e.g. `UserRepl.Commands`), where they translate
messy human input into the clean config DTO and can validate/warn/throw before
the engine ever runs.

**Why:**

- A fluent chain implies sequence; the engine's semantics are holistic (e.g.
  LMP rescaling renders global bandwidth dead — a chain misrepresents that).
- "Define the state of the universe and hit run" is the intended mental model.
- Reproducibility is trivial when **the manifest IS the config record**:
  serialize the DTO, deserialize later, get the same manifold.
- The engine stays a deterministic mathematical kernel.

**Shape conventions:**

- Top-level config: one record with `init`-only properties, `required` for the
  few non-defaultable fields, everything else nullable or sensibly defaulted.
- Per-stage sub-configs: nested records carrying a strategy enum (`Auto`,
  `Explicit…`) plus that strategy's parameters.
- Entry point: `public static XxxResult Build(XxxConfig config, …)` — no
  instance state, no method chaining. Delegates and live objects ride the
  call, never the config (configs stay JSON-serializable).
- Emit a manifest snapshotting the config, diagnostic reports, and the
  automation log (which auto-picks fired and why).

**The smell that triggers this shape:** a static method whose parameter list
has grown to encode discrete pipeline stages (≥6 params, several of them
enum/strategy selectors). Don't add parameters or sibling statics — make the
stages first-class in a config record and keep any chained ergonomics at the
shell.

## Refactoring discipline

Two principles — one for the end state, one for the journey. Keep them
distinct.

### End state: clean breaks, no compat bandaids

Pre-release code has nothing in the wild depending on its current shape, so
compat shims are pure cost: they accrete, delay the real cleanup, and keep the
old shape cluttering the new one ("I keep chasing bandaids and compromises").

- Replacing an API? The old surface is deleted and every caller migrated — no
  `[Obsolete]` attributes, no shims routing old calls into the new structure.
- Changing on-disk serialization shapes (manifests, presets)? Rewrite the
  schema without versioning logic — no `SchemaVersion` branch handling.
- Moving files? Move them and let compile errors guide caller updates — no
  re-export stubs at the old path.
- Helpers whose only callers used the old API get deleted too, not kept as
  "potentially useful primitives."
- The changelog describes the clean break without apology — it's intentional.

This holds as long as the project is pre-release; revisit when external
consumers exist.

### Journey: defer deletes to a cleanup sweep

During a multi-step refactor, write new files alongside the old ones and leave
the old ones in place until a dedicated cleanup pass at the end.

**Why:** the build never breaks mid-flight (unmigrated callers still resolve);
the old code is right there for A/B comparison or fallback; "build new shape"
and "delete old shape" land as separate commits, which reads better in
history. Git is the safety net — nothing to lose.

**In practice:** add new files in their permanent locations (no transitional
dirs), migrate callers incrementally, and make the final cleanup task the
*only* place deletions happen — as its own discrete commit.

## Naming & vocabulary

### Functional names over project prefixes

Name CLI commands, scripts, namespaces, and assemblies by what they do
(`repo-audit`, `file-ownership`), not which project owns them (`spcx-audit`,
`ps-core-x`). Functional names are self-describing and transferable.

### Cross-field vocabulary, chosen by fitness

Vocabulary is drawn deliberately across math, physics, topology, CS, and
cybernetics — the fittest, most precise term per spot, with **no single field
as the anchor**. Field terms are first-class and often exact, not lossy
metaphors: *observables* (physics: `field → value`), *Betti* (topology:
`b₁ = E−V+C`), *currency*, *kernel*.

**Fitness test for a borrowed term:** intuitive + native to a serious
community that studies the object + carries its weight without loss of
generality. Failing either way — intuitive-but-narrowing (drags baggage) or
general-but-inert — disqualifies the borrow.

**Reason in general terms, reify concrete ones.** Math vocabulary
(functionals, maps, fields-over-indexes) is prized as reasoning scaffolding
for its generality, but it is not reified: talk in terms of functionals on the
way to writing `IReduction` — there is no `IFunctional`. Build the concrete
thing, named from whatever field fits best.

### One unified shape per pattern

When a structural pattern (namespaces, file layouts, API shapes) is sketched
with examples, render the **simplest unified shape** across them. Do not
derive asymmetric treatments from incidental implementation differences
(stateful vs pure, etc.) — that multiplies decisions and obscures the pattern.
If a detail genuinely doesn't fit, surface it as a single question, not a fork
in the proposal.

## Documentation in code

Docstrings are lean operational references: what it does + attribution + the
one or two behavioral facts that matter (e.g. *"mean-field formulation due to
Wang 2020; computes the closed form instead of Monte Carlo"*).

- **No academic reviews in docstrings** — critiques, full derivations, lemma
  proofs, and prior-art quibbles are parked in `issues/` design notes (see
  DocOps in [AGENTS.md](AGENTS.md)).
- **Locating narrative:** weave the project's function-algebra vocabulary
  (*field / accumulate / reduce / observable / currency / sampler / solver /
  affinity / form-degree*) as one clause locating the thing in the recurring
  shape — e.g. *"accumulates the per-draw bond events and reduces them to the
  `Affinities` currency."* This is load-bearing prose: it lets a reader
  scanning docstrings see the parallel structure across tiers. A clause, not a
  paragraph.
- **cref discipline:** `<see cref="…"/>` the vocabulary nouns that are real
  types (`Affinities`, `Accumulator`, `IReduction`, …) so the parallel
  structure is navigable, not just narrated.

## Testing standards

- **Oracles over consumers.** Engine correctness is established against its
  own spec: R oracles (`r/oracles`), brute-force references, and property
  checks — not by wiring up a downstream application.
- **Test layout.** xUnit test projects live in `tests/` mirroring the
  `projects/` split (`Maths.LinAlg.Tests`, …). The parallel fact runner
  (`src/test-harness`, invoked per [AGENTS.md](AGENTS.md)) drops per-suite
  manifests and a `summary.json` under `artifacts/test-runs/`.
- **RNG hygiene.** Derive seed streams structurally (`SeedTree` /
  `SplitMix64`), never by ad-hoc arithmetic (`seed + k·block`, `seed ^ C`) —
  arithmetic derivation aliases streams. Prefer pinned generators
  (`Xoshiro256PlusPlus`) over `System.Random` where cross-version stream
  stability matters. Regenerating pinned fixture expectations is a deliberate,
  called-out act, not a side effect.

## PowerShell conventions

- Default to `-LiteralPath` over `-Path` for all filesystem cmdlets
  (`Get-Content`, `Test-Path`, `Move-Item`, …). `-Path` treats `[ ] * ?` as
  wildcards, so user-named or crawled files silently fail or mismatch. Use
  `-Path` only when wildcard expansion is explicitly intended.

## Antipatterns — quick reference

- Compat shims, `[Obsolete]` markers, schema-version branches (pre-release).
- Fluent builders in the backend; sequence-dependent backend semantics.
- Wiring a consumer to validate an in-progress engine.
- Deleting superseded files mid-refactor instead of in the final sweep.
- Project-prefixed names for functional tools.
- Asymmetric renderings of a single sketched pattern ("oversplitting").
- Academic critique in docstrings.
- Arithmetic seed derivation; unpinned RNG where streams are fixtures.
- `-Path` on filesystem cmdlets operating on discovered paths.

---

*Distilled 2026-07-21 from the feedback memories staged in `issues/docs/`.*
