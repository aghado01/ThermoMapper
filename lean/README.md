# lean — the rigor harness

Retrospective formalization of the project's load-bearing math claims:
work backward from implemented engineering to the implicit mathematical
contracts it rests on, then verify them in Lean 4 + mathlib. Theory-level
unit tests, sibling of `tests/`.

## Taxonomy

Three tiers, promotion-driven:

| Tier | Where | Contract |
| --- | --- | --- |
| **proto-lemmas** | `proto-lemmas/*.md` | Informal lemma briefs (chat-exported). Informal, markdown format. Not compiled. |
| **enthymemes** | `Enthymemes/*.lean` | Statements are real and **compile**; proofs still apologize (`sorry`). An enthymeme suppresses a premise — these suppress proof steps. |
| **lemmas** | `Lemmas/*.lean` | No apologies. Complete as far as current concerns go. |

An enthymeme is **promoted** when it stops apologizing: move the file (or the
finished declarations) to `Lemmas/`. Declaration names and namespaces never
change on promotion — tier is location, not identity.

## meta-CI

`scripts/meta-ci.ps1` enforces the taxonomy:

1. `lake build` green — both tiers must compile (enthymemes owe proofs, never
   statements).
2. The lemmas tier never apologizes: no `sorry` token under `Lemmas/` (not
   even in prose), and no `import Enthymemes.*` — an import could launder a
   sorried dependency into a "proved" lemma without the token appearing.
3. Enthymeme ledger: each file is `unstated` (no declarations yet),
   `apologizing(n)`, or `PROMOTION-READY` (declarations, zero sorries —
   move it).

Flags: `-Validate` makes ledger notices fail (strict mode); `-NoBuild` skips
the compile gate (CI runs it after lean-action has already built).

## Build notes

- Toolchain: Lean 4.30.0 + mathlib v4.30.0 (cache fetched — never build
  mathlib from source; `lake exe cache get` if it ever goes missing).
- A file importing `Mathlib` elaborates ~8 min on first touch, ~25 s warm.
- `SimpleGraph.loopless` is `Std.Irrefl` — a structure; tactic proofs need
  `constructor` before `intro`.
- The mathlib standard linter set is on; the copyright-header linter is off
  (`lakefile.toml`).
- The `.github/workflows/` here are inert while `lean/` lives inside the
  monorepo — GitHub only honors workflows at a repo root. They activate if
  this directory is ever split into its own repository; until then,
  `scripts/meta-ci.ps1` is the gate.
