# lean — the rigor harness

Retrospective formalization of the project's load-bearing math claims:
work backward from implemented engineering to the implicit mathematical
contracts it rests on, then verify them in Lean 4 + mathlib. Theory-level
unit tests, sibling of `tests/`.

## Taxonomy

Four active stages, promotion-driven, plus a curated archive:

| Stage | Where | Contract |
| --- | --- | --- |
| **protolemmata** | `Protolemmata/*.md` | Informal briefs, conjectures, design arguments, and candidate statements. Expected to change; not compiled. |
| **enthymemata** | `Enthymemata/*.lean` | Statements are real and **compile**; proofs may still apologize (`sorry`). An enthymema suppresses a premise — these suppress proof steps. |
| **lemmas** | `Lemmas/*.lean` | Verified, reusable results with no apologies and no dependency on `Enthymemata`. |
| **theorems** | `Theorems/*.lean` | Verified results whose consequence and stability merit treatment as named project deliverables. The distinction from lemmas is curatorial, not a Lean keyword distinction. |
| **archeion** | `Archeion/` | Superseded or retired material preserved with provenance. A side exit from the active path, not a maturity stage; not compiled. |

An enthymema becomes **eligible for promotion** when it stops apologizing.
After its statement and dependencies pass semantic review, move the file (or
the finished declarations) to `Lemmas/`. Declaration names and namespaces
never change on promotion — stage is location, not identity. Promotion from
`Lemmas/` to `Theorems/` is optional and deliberate: most verified results
remain lemmas. Material may instead leave any active stage for `Archeion/`
when it is superseded or retired.

## meta-CI

`scripts/meta-ci.ps1` enforces the taxonomy:

The script prefers the portable Lake executable under
`$env:PORTABLE_ROOT\elan\bin`, then falls back to `lake` on `PATH`. When it
uses the portable executable, it supplies a process-local `ELAN_HOME` if the
calling process inherited a stale environment snapshot.

1. `lake build` green — all active formal stages must compile (enthymemata owe
   proofs, never statements).
2. Each active aggregate (`Lemmas.lean`, `Theorems.lean`, and
   `Enthymemata.lean`) imports every `.lean` file in its stage, so a draft
   cannot silently evade the build.
3. Active Lean source uses scoped Mathlib modules; the `import Mathlib`
   umbrella is forbidden.
4. The verified stages never apologize: no `sorry` token in their Lean source
   (including comments), and no `import Enthymemata.*` — an import could
   launder a sorried dependency into a verified result without the token
   appearing. `Lemmas` also cannot depend upward on `Theorems`.
5. Enthymema ledger: each file is `unstated` (no declarations yet),
   `apologizing(n)`, or `PROOF-CLOSED` (declarations, zero sorries — a
   candidate for semantic review and promotion, not an automatic endorsement).

Flags: `-Validate` makes ledger notices fail (strict mode); `-NoBuild` skips
the compile gate (CI runs it after lean-action has already built).

## Import discipline

Lean modules import the narrowest stable public Mathlib modules that own the
declarations, notation, and tactics they actually use. Active source must not
use the `Mathlib` umbrella; `scripts/meta-ci.ps1` enforces that boundary.

- Add a scoped import when a new dependency first arises rather than loading
  likely future dependencies in advance. Import tactics explicitly, such as
  `Mathlib.Tactic.NormNum`.
- Prefer a declaration's owning public module (or a coherent public feature
  module) over an unrelated transitive import. The goal is an intentional,
  stable dependency boundary, not the smallest accidental import closure.
- Compile the affected module directly after changing its imports. When code
  stops using a feature, remove the corresponding import and compile again.
- Apply the same discipline to Lean examples in `Protolemmata/`, so
  promotion into an active stage does not reintroduce the umbrella.
- When a new import lacks cached artifacts, fetch that module specifically,
  for example `lake exe cache get Mathlib.Analysis.Calculus.Deriv.Polynomial`.
  The cache tool follows its transitive imports. A bare `lake exe cache get`
  starts from `Mathlib.lean` and restores the entire library, so use it only
  when a complete Mathlib cache is intentional.

## Build notes

- Toolchain: Lean 4.32.2 + mathlib v4.32.2. Fetch missing artifacts by scoped
  module as described above; artifacts unavailable from the cache build
  locally.
- The `Mathlib` umbrella was measured at 5–12 min to elaborate cold and about
  a minute warm in this harness; scoped imports also make dependency changes
  reviewable.
- `SimpleGraph.symm` and `SimpleGraph.loopless` are `Std.Symm` and
  `Std.Irrefl` structures; tactic proofs need `constructor` before `intro`.
- The mathlib standard linter set is on; the copyright-header linter is off
  (`lakefile.toml`).
- The `.github/workflows/` here are inert while `lean/` lives inside the
  monorepo — GitHub only honors workflows at a repo root. They activate if
  this directory is ever split into its own repository; until then,
  `scripts/meta-ci.ps1` is the gate.
