# AGENTS.md — Runtime Orientation, Toolchain & Agent Protocol

Operational guidance for agents working in ThermoMapper: repository layout, environment resolution, standard toolchain invocations, execution/repair limits, and `issues/` DocOps. Engineering principles, architecture invariants, and testing policy live in [CONTRIBUTING.md](CONTRIBUTING.md) — read it before writing code. Harness-specific notes live in [CLAUDE.md](CLAUDE.md).

## 1. Repository Layout

| Path | Contents |
|---|---|
| `src/` | C# source by functional area (`archivory`, `clustering`, `graphs`, `hashish`, `maths`, `repo-audit`, `storage-audit`, `synthetic`, `tda`, `test-harness`, `user-repl`, `viz`, `viz-core`) |
| `projects/` | `.csproj` files (`Area.Module`, e.g., `Maths.LinAlg`, `Graphs.Observables`, `Clustering.Graphical.SPC`) linking files from `src/` |
| `tests/` | xUnit test assemblies (`*.Tests`) and smoke projects mirroring `projects/` |
| `scripts/` | PowerShell entry points (`repo-audit.ps1`, `storage-audit.ps1`, `fact-harness.ps1`, `test-select.ps1`) |
| `r/` | R validation oracle package (`renv`-managed): `oracles/`, `scripts/oracle-ci.ps1`, `scripts/r-session.ps1`, `setup.ps1` |
| `lean/` | Lean 4 rigor harness (`Protolemmata` → `Enthymemata` → `Lemmas`, optional promotion to `Theorems`, non-building `Archeion` archive), plus `scripts/meta-ci.ps1` and `parse-lean-docs.ps1` |
| `datasets/` | Benchmark datasets (`iris`, `isolet`, `landsat`, …) and prep material |
| `presets/` | Named run configurations (e.g., `presets/spc`) |
| `docs/` | Canonical evergreen documentation ([Architecture](docs/Architecture.md), [Design](docs/Design.md), module specs) — currently frozen during cleanup |
| `artifacts/` | Generated build outputs, test-run manifests, audit reports, and probe outputs — never hand-edited |
| `issues/` | Work coordination, briefs, design notes, and reports (see [DocOps](#5-docops--multi-session-coordination)) |

## 2. Environment & Shell Invariants

- **Portable toolchain resolution.** External toolchains (`dotnet`, `R`, `lean`/`lake`, `python`, `llama-cpp`) come from `$env:PORTABLE_ROOT` and reach shells via the user-registry `PATH`. Do not gate on or assume legacy `SHARED_ENV_LOADED` profile variables/aliases. Because spawned subshells may inherit a stale `PATH` snapshot, never treat a failed `Get-Command` probe as proof a tool is absent — resolve directly under `$env:PORTABLE_ROOT` or ask.
- **Python virtual environments.** Execute Python tasks inside a project `.venv` backed by the portable Python runtime under `$env:PORTABLE_ROOT`, rather than global `site-packages`.
- **PowerShell path safety.** Default to `-LiteralPath` over `-Path` for all filesystem cmdlets (`Get-Content`, `Test-Path`, `Move-Item`, `Remove-Item`, etc.). `-Path` expands `[ ] * ?` as wildcards and silently fails or mismatches literal paths. Use `-Path` only when globbing is explicitly intended.
- **Line-ending preservation.** Files under `docs/` use CRLF and repository line endings are mixed. Edit files in place without wholesale LF conversions (e.g., `sed -i` or universal-newline rewrites) and verify `git diff --stat` before committing.
- **Restore timestamps.** Restoring a file with `Copy-Item` preserves the backup's modification timestamp, causing MSBuild to skip recompilation. Touch restored source files so MSBuild rebuilds dependent assemblies.
- **Lean import discipline.** Active Lean modules import only declaration- or feature-owning Mathlib modules (`import Mathlib` is prohibited and rejected by `lean/scripts/meta-ci.ps1`).
- **Artifact confinement.** Write all build outputs, test logs, and probe data under `artifacts/`. Never write notes, temporary files, or probes outside the repository (including external app side-session folders).

## 3. Standard Toolchain & Utility Invocations

| Domain | Task | Standard Command | Options & Notes |
|---|---|---|---|
| **C# Build / Test** | Build or test with isolated outputs | `dotnet build projects/<Name>/<Name>.csproj -v minimal --nologo "-p:ArtifactsRoot=<repo>\artifacts\<session>\"` | In shared checkouts, pass `-p:ArtifactsRoot=<repo>\artifacts\<session>\` on every `dotnet build` and `dotnet test`. Invoke from PowerShell with the property quoted and trailing backslash intact (Git Bash drops backslashes and writes stray `bin/`/`obj/` beside `.csproj`). Add `-c Release` for perf parity; `-p:ExtraDefineConstants=<FLAG>` selects reference solvers (e.g., `EIGEN_REFERENCE`). |
| **Targeted Tests** | Change-selected run (`Edit` / `Commit` / `Slice`) | `pwsh scripts/test-select.ps1 -Scope Edit -Changed <paths> -ArtifactsRoot artifacts/<run>` | Runs in Release at below-normal priority, xUnit capped at 6 threads, build servers off. Pass `-Changed` (comma-separated in bash) in shared checkouts to exclude other sessions' uncommitted work. `-Scope Commit` adds `[Trait("Tier", "Pinned")]`; `-Scope Slice` runs the full assembly. `-ListOnly -Explain` previews selection; `-IncludeSlow` includes `slow-tests.txt`; `-WriteSlowList <trx> -SlowSeconds 10` updates `slow-tests.txt`; `-SourceRoot <dir>` targets a snapshot. |
| **Parallel Tests** | Per-case Release fact runner | `pwsh scripts/fact-harness.ps1 -Project tests/<Suite>/<Suite>.csproj` | Defaults to `VizCore.Tests` in `Release`. Options: `-Configuration Release\|Debug`, `-MaxWorkers <n>`, `-ListOnly`. Writes suite manifests and `summary.json` under `artifacts/test-runs/`. |
| **Qualification / Sweeps** | Slice gates & stochastic stability | `$env:BASS_QUALIFICATION=1; dotnet test tests/<Suite>/<Suite>.csproj -c Release --filter <Expr>`<br>`$env:SEED_SWEEP=<n>; dotnet test ...` | `SEED_SWEEP` evaluates registered `ComparativeSweep` claims across `n` seeds per axis and writes tables under `artifacts/test-runs/SeedSweep/`. |
| **Repo Audit** | Layering & invariant validation | `pwsh scripts/repo-audit.ps1 -Validate -NoDisplay` | Builds/runs `RepoAudit.exe`. Options: `-Impact <path>`, `-NoGit`, `-Rebuild`. |
| **Storage Audit** | Artifact disk profile & cleanup | `pwsh scripts/storage-audit.ps1 profile -NoDisplay`<br>`pwsh scripts/storage-audit.ps1 plan`<br>`pwsh scripts/storage-audit.ps1 apply -PlanPath <plan> -Execute` | Previews unless `apply` and `-Execute` are supplied; see [src/storage-audit/README.md](src/storage-audit/README.md). |
| **Lean Rigor** | Taxonomy & build gate | `pwsh lean/scripts/meta-ci.ps1` | Enforces `lake build`, barrel import reachability, scoped Mathlib imports, and zero `sorry` in `Lemmas`/`Theorems`. |
| **Lean Build / Docs** | Compile Lean / convert Verso HTML | `lake build` *(in `lean/`)*<br>`pwsh lean/parse-lean-docs.ps1` | Compile individual affected modules with `lake build <Module>` during iteration. |
| **R Oracles** | Oracle bridge & `renv` sync check | `pwsh r/scripts/oracle-ci.ps1` | Resolves R via `r/scripts/r-session.ps1`, verifies `renv.lock` sync, and smoke-tests a base oracle. Pass `-Restore` to run `renv::restore`. An R-oracle `TimeoutException` under heavy load is rerun at normal priority before treating it as a failure. |

**Shared machine load.** Run heavy or long-running benchmarks and harnesses at below-normal priority with capped workers so concurrent interactive work remains responsive.

## 4. Execution, Repair & Commit Protocol

- **Follow the engineering contract.** Apply [CONTRIBUTING.md](CONTRIBUTING.md) rules during execution: engine-first validation, substrate decoupling, unified pattern shapes, sibling sweeps on bug fixes, and staged refactors (new files beside old; deletions in a dedicated final commit).
- **Surface removals; respect read-only scopes.** Removing an option, rule, or output must be put to the user as an explicit question stating what is lost and what replaces it. When a brief or prompt asks to "check X," report findings on X without modifying code.
- **Pre-commit diff review.** End each step with a review pass: read the full diff as a first-time reviewer, check for options and refusal paths that no test exercises, fix gaps, and repeat until a pass finds no defects. Record each pass's findings in the step report.
- **Bounded repair loops.** Run an initial check plus at most **three** targeted repair/rerun cycles for the same failure. If still failing, record the exact command, failure output, and remaining blockers in the task report and stop or continue independent scoped work. Never weaken contracts, numerical tolerances, or test assertions merely to force a pass.
- **Escalate tooling defects.** When `repo-audit`, `test-select`, or another project utility produces false positives or false negatives, verify the root cause in source and report it to the user as a bug/enhancement candidate rather than silently working around it.
- **Single-step pathspec commits & hygiene.** Because sessions share one Git index, commit in a single step by explicit pathspec (`git commit -m "<subject>" -m "<body>" -- <paths>`) and leave nothing staged between commands. In interactive sessions, land targeted per-concern commits directly on `main` unless directed otherwise. Write an imperative subject and a concise body stating what rule/behavior changed and why. Do not include AI co-authorship trailers (`Co-Authored-By: ...`). Do not claim enforcement by a checker or policy file that does not exist.

## 5. DocOps & Multi-Session Coordination

- **Documentation ownership.** `docs/` holds current canonical specifications, rewritten in place without chronological work logs, personal context, or machine-specific paths; `docs/` is currently frozen during cleanup and updated offline in curated batches. During feature and fix work, sessions leave `docs/` untouched and record required updates (including settled design/statistical trade-offs, alternatives considered, and user implications) under **Doc changes owed** in their report. `issues/` holds active work coordination organized by workstream subfolder (`ph`, `regression`, `spred`, `viz`, `docs`, …).
- **Canonical `issues/` artifacts.** Cross-reference pages under `issues/` using wikilinks or explicit relative links.
  - **Briefs (`*-brief.md`)**: scoped assignments stating status, dependencies, scope, acceptance checks, and reversal before delegated or deferred work starts.
  - **Reports**: dated implementation records appended to or placed beside the governing brief, closing the loop.
  - **Design notes (`design.md`, `dev-sequence.md`)**: settled decisions with rationale, and the designated home for academic critiques, derivations, and prior-art reviews kept out of docstrings.
  - **Discussion digests (`*-thread.md`, `*-transcript-*.md`)**: distilled ideation and reviews.
  - **[The Shelf](issues/SHELF.md)** (`issues/SHELF.md`): candidate action items incubating prior to promotion into a brief.
  - **[The notation ledger](issues/notation-ledger.md)** (`issues/notation-ledger.md`): terminology/identifier splits logged newest-first for scheduled unification passes (record only; do not rename inline).
- **Closeout & report format.** Mark a brief/task `done` only when all required scope and acceptance checks pass; use `implemented, verification pending` if required checks are unverified or skipped. Append dated entries in this structure:
  ```markdown
  ### YYYY-MM-DD Implementation or follow-up
  - Outcome: status (`done` | `implemented, verification pending`) and pending gates.
  - Changes: behavior, owned paths, contracts, and deviations.
  - Verification: exact commands, environments, passed/failed/skipped counts, and `artifacts/` paths.
  - Doc changes owed: `docs/` updates required by the change (trade-offs, behavior, specs).
  - Remaining work: destination and prerequisite for each unresolved item.
  - Closeout: scratch cleanup, commit hash (or no-commit reason), and remaining uncommitted changes.
  - Reversal: undo steps preserving durable data and evidence.
  ```
- **Multi-session ("chip") coordination.** Keep handover state (`In progress` entry stating what is uncommitted, what is running, and what comes next) in the governing brief's report, updated in place and folded into the slice report on completion. Chips message each other directly when they need to coordinate; inter-session messages stay terse with only needed context and point to the governing in-repo report. Sessions do not read each other's transcripts.
