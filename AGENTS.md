# AGENTS.md — project orientation & working guidance

Project-specific guidance for agents: how the repo is laid out, the invariants
that govern placement and dependencies, the tooling, and how to behave while
working here. Development practices, design principles, and testing standards
are the province of [CONTRIBUTING.md](CONTRIBUTING.md) — read it before
writing code.

## Project orientation

| Path        | What lives there |
| ----------- | ---------------- |
| `src/`      | C# source, one folder per functional area (`maths`, `graphs`, `clustering`, `tda`, `viz`, `viz-core`, `synthetic`, `hashish`, `archivory`, `user-repl`, `repo-audit`, `test-harness`) |
| `projects/` | `.csproj` files, one per assembly (`Maths.LinAlg`, `Graphs.Observables`, `Clustering.Graphical.SPC`, …) pointing into `src/` |
| `tests/`    | xUnit test projects mirroring `projects/` (`*.Tests`), plus harness plumbing |
| `scripts/`  | PowerShell entry points (`repo-audit.ps1`, `fact-harness.ps1`, `portable-python.ps1`, `venv-boostrap.ps1`, `parse-lean-docs.ps1`) |
| `r/`        | R oracle package (renv-managed) — reference implementations that validate the C# engine (`r/oracles`) |
| `lean/`     | Lean 4 rigor harness (`Protolemmata` → `Enthymemata` → `Lemmas`, with consequential results in `Theorems` and retired material in `Archeion`) |
| `datasets/` | Benchmark datasets (`iris`, `isolet`, `landsat`, …) + prep and reference material |
| `presets/`  | Named run configurations (e.g. `presets/spc`) |
| `artifacts/`| Build output, test-run manifests, health reports — generated, never hand-edited |
| `issues/`   | Design notes and work coordination — see DocOps below |

## The one-way layer order

```
maths  →  graphs  →  { tda, clustering, viz, user-repl }
maths/topology  →  tda/ph  →  tda/mapper  →  tda/pipelines
```

Dependencies point left (downward). `graphs` is a neutral substrate: it emits a
weighted `CsrGraph` + provenance manifest, and SW/SPC/PH/Mapper are _consumers_
of that artifact. Placement authority for the topology stack lives in
`tda-placement.md` (MarkBrain: `ThermoMapper/issues/tda-purification/persistent-homology/`).
There is no `tda/primitives` — everything PH lives in `tda/ph`.

## The construction invariant

Every graph-construction stage (`src/graphs/**`) references only `graphs` and
`maths` namespaces. Consumer and engine semantics enter construction as **data
through typed contracts** — config values, injected delegates, precomputed edge
scores — never as an import of a consumer or PH type.

Worked example (violation cleared 2026-07-03): `GraphCompiler` used to import
the PH engine and run involuted persistence mid-build
(`H1CycleEdges.FromDistanceGraph`) to compute its LMP protect-set. The fix is
the pattern to repeat: `GraphCompiler.Build(..., ProtectedEdgeSource? protectedEdges)`
— the caller (UserRepl commands) computes the H1 protect-set above the seam and
injects it as a delegate over the distance graph. Future topology-scored passes
(the graph sculptor's persistence/spectral criteria) flow in the same way:
scores in, neutral substrate out.

## The litmus test - under construction

`CsrGraph` round-trips to disk (`WriteTo`/`FromBinary`). The compiler is
decoupled exactly when a non-SPC reader — a pure PH run, a Mapper nerve, a
spectral embedding — consumes the output unchanged. Needing anything
SPC-specific from construction means fusion has crept back in: re-cut the seam
instead of threading the consumer concern into a build stage.

## Keep the boundary artifacts boring

- The compile boundary emits an immutable, serializable artifact + manifest;
  runs snapshot their full parameter set (config-artifact provenance).
- Configs are declarative, JSON-serializable DTOs. Delegates and live objects
  ride the call, not the config.
- Fluent/chained ergonomics live at the CLI/REPL boundary; the backend consumes
  the pure DTO. (Rationale and shape conventions: *strict core, fluent shell*
  in [CONTRIBUTING.md](CONTRIBUTING.md).)
- Pre-release discipline: no back-compat shims, no `[Obsolete]` aliases —
  superseded surfaces are deleted wholesale.

## How to work here — behavioral guidance

- **Engine first.** Don't push toward wiring up consumers/applications to
  validate in-progress engine work — validate with oracles and property
  checks, and sequence engine depth ahead of integration (see
  [CONTRIBUTING.md](CONTRIBUTING.md)).
- **Match the sketched shape.** When a structural pattern is proposed with
  examples, render one unified shape; don't invent asymmetries the examples
  didn't ask for. Genuine misfits get raised as a single question.
- **Refactors:** new files land beside old ones; deletions happen only in the
  final cleanup sweep, as their own commit. End state carries no compat shims.
- **Commits:** work lands directly on `main` as targeted, per-concern commits
  during interactive sessions.
- **Escalate tool inaccuracies.** When repo-audit (or other project tooling)
  produces false positives/negatives, verify root cause in source, then flag
  it to the user as a candidate bug fix or enhancement — don't silently work
  around it.

## Environment

- Toolchains (dotnet, R, Lean, python, llama-cpp, …) come from the portable
  environment (pdenv) and reach shells via the **user-registry PATH**;
  `$env:PORTABLE_ROOT` points at the install root. The old
  bootstrap-profile/`SHARED_ENV` mechanism is dead — never gate on
  `SHARED_ENV_LOADED` or assume its aliases exist.
- A registry PATH change only reaches processes started afterward, so a
  spawned tool shell may hold a **stale PATH snapshot**. A failed
  `Get-Command` probe therefore cannot distinguish "not on PATH yet" from
  "not installed" — don't conclude from the probe; use an absolute path under
  `$env:PORTABLE_ROOT`, or ask.
- Python work goes through the local venv (`scripts/portable-python.ps1`,
  `scripts/venv-boostrap.ps1`). Convenience-alias loader scripts for R /
  dotnet / Lean are planned ([TODO.md](TODO.md)).
- PowerShell filesystem cmdlets: `-LiteralPath` by default (see
  [CONTRIBUTING.md](CONTRIBUTING.md)).

## DocOps — conventions for `issues/`

`issues/` is the design-note and work-coordination space, organized by
workstream subfolder (`ph`, `spred`, `viz`, `docs`, …). Canonical file kinds:

- **Briefs** (`*-brief.md`) — scoped work assignments, JIRA-ticket-like:
  context, scope, constraints, sequencing. Written before delegated or
  deferred work begins.
- **Reports** — writeups of work performed: what was done, open ends, caveats,
  follow-up items. A report on delegated work is appended to (or sits beside)
  the brief that commissioned it, closing the loop.
- **Design notes** (`design.md`, `dev-sequence.md`, roadmaps) — settled
  decisions with rationale. Also the designated parking lot for academic
  critique, derivations, and prior-art discussion kept out of docstrings.
- **Discussion digests / transcripts** — distilled ideation and reviews
  (`*-thread.md`, `*-transcript-*.md`, retrospectives).
- **[The Shelf](issues/SHELF.md)** (`issues/SHELF.md`) — action items jotted
  organically mid-work; entries incubate and may be amended before being
  promoted to a brief.

Docs reference each other by wikilink or explicit relative link. Generated
analysis output goes to `artifacts/`, not `issues/`.
