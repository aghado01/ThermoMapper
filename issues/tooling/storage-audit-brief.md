# storage-audit: profiling this repository's disk use, and a cleanup the user approves (chip brief)

Date: 2026-10-31. The D: drive is 99% full (2.8 GB free of 239 GB on 2026-10-31), and this repository's work has been eating it steadily: build outputs, chips' artifacts roots, test and benchmark runs, snapshots, and git's objects. The user is formatting a thumb drive to take what is worth keeping.
- *What this chip does:*
  - surveys this repository's disk use and the temp data its work produces;
  - proposes a cleanup, item by item, which the user approves;
  - builds a reusable tool, `storage-audit`, so that profiling and cleanup can run again as routine.
- *Scope:* this repository, `D:\aghado01\ThermoMapper`, and the temp data its sessions write, nothing else on the drive.

Read first: [AGENTS.md](../../AGENTS.md) (its tooling section and `repo-audit`, the pattern to follow) and [CONTRIBUTING.md](../../CONTRIBUTING.md).

## Safety, before anything else

- **Nothing is deleted or moved without the user's explicit approval in chat, batch by batch.** The survey and the plan are read-only. Show each batch with its paths, sizes and category, and wait for a clear yes. Approval of one batch does not extend to the next.
- **Never touch:**
  - tracked files;
  - the user's notes, the untracked files under `issues/` included;
  - source, briefs or docs;
  - another session's working files.
- **Chips are running beside you,** each with its own artifacts root: `artifacts\chip-fh3\` and later `chip-fh4\` (the flat-histogram chip), `artifacts\chip-fr1\` (F-R1) and `artifacts\chip-ext\` (the extension seams). Treat these, and any path written to in the last 24 hours, as active, and leave them alone. They are cleaned by their own chips.
- **Git:** propose `git gc`, `git worktree prune` or a reflog expiry with the space each would free, and run one only on approval. Never rewrite history.
- **Disk space is nearly gone,** so keep your own footprint small: build into `artifacts\chip-storage\` and delete it when you finish, and avoid writing large intermediate files.

## The steps

1. **The survey and the quick wins** (read-only; aim to report within the first half hour).
   - *Sizes,* by top-level directory, then by kind:
     - build outputs (`bin/`, `obj/`, `artifacts\bin`, `artifacts\obj`);
     - chips' and sessions' artifacts roots (`artifacts\chip-*`, `artifacts\main-*`);
     - test, smoke and benchmark runs (`artifacts\test-runs`, `smoke-runs`, `bench`, `benchmarks`, `reports`, and the like);
     - caches (`graph_cache`, `doc-dive`, `gallery`, …);
     - reposnapshot payloads (`src/.snapshot` and any others);
     - git: `.git`'s packs and loose objects (`git count-objects -vH`), worktrees (`git worktree list`, and any `.claude\worktrees` or orphaned worktree folders), and the reflog;
     - large untracked or ignored files (`git status --ignored`).
   - *The temp data* the work writes: the sessions' scratchpads and task outputs under `C:\Users\azrie\AppData\Local\Temp\claude\D--aghado01-ThermoMapper\`, and any temp written to D:. Note which drive each is on: only D:'s frees the D: drive.
   - *Each item gets a category,* with its reason:
     - **regenerable:** a build output or a cache that a build or a run recreates, safe to delete;
     - **archivable:** worth keeping but not here (old run outputs, old snapshots), to move to the thumb drive;
     - **keep;**
     - **ask:** unclear, so ask the user.
   - *The quick-win list:* the regenerable, inactive items, largest first, with the space they free. Present it to the user for approval, and on approval delete those batches and report the space freed.
2. **The plan for the rest.** The archivable and ask items, grouped, with sizes and a proposed action each. Archiving waits until the user gives the thumb drive's path. Moves are verified (sizes and file counts compared) before the source is removed, and only on approval.
3. **The tool, `storage-audit`,** on the pattern of `repo-audit`: C# under `src/storage-audit/`, its project under `projects/`, its tests under `tests/`, and an entry point `scripts/storage-audit.ps1`. It needs no dependency beyond the repo's own.
   - `profile`: the survey above as a report (markdown and JSON), with sizes by category and by directory, the largest items, the worktrees, git's object statistics, and the age of each item's newest file.
   - *A rules file in the repo,* mapping paths and globs to categories with reasons. Its defaults come from the survey, and the user can edit it.
   - `plan`: the proposed actions with their sizes, written to a plan file the user reviews.
   - `apply`: executes a reviewed plan file only.
     - It is a dry run unless `--execute` is given.
     - It moves archivables to `--archive-to <path>`, with a manifest, and verifies before removing the source.
     - It deletes regenerables.
     - It refuses tracked files, paths outside the repo (temp roots excepted, when listed in the rules), and paths written to within a recent window.
     - It logs every action to a ledger.
   - *Tests* on a temporary fixture tree: the rules' classification, the plan, the dry run, the archive's verification, and every refusal.
   - *Docs:* a usage section in the tool's own README, and a one-line entry in AGENTS.md's tooling table.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Switch no branch, push nothing.
- **Your files:** `src/storage-audit/`, its project and tests, `scripts/storage-audit.ps1`, the rules file, the tool's README, the AGENTS.md tooling line, and this brief's Reports section.
- **Commit your own files only, by explicit path,** in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. The user's untracked notes and other sessions' edits stay unstaged; before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.
- **Each step ends with a fresh-eyes review,** for correctness and for coverage gaps (the refusals above all), fixed before the next.

## Report

Under "Reports" below, newest first:
- the survey's table;
- what was cleaned on approval, and the space freed;
- the plan for the rest;
- the tool's commands and rules;
- its tests with their counts;
- anything the user must decide.

## Reports

Newest first.

### 2026-10-31: survey, quick-win cleanup, the tool (chip `storage-audit`)

**The survey** (read-only, 14:10). The repository held about 4.9 GB of D:'s 236 GB used, so cleaning it frees at most about 5 GB. Whatever fills the rest of the drive is outside this chip's scope.

| Path | Size | Newest | Category |
|---|---:|---|---|
| `lean/.lake/packages/*/.lake/build` (Mathlib 2,398 MB, aesop, batteries, Qq, …) | 2,649 MB | 10-24 | regenerable (`lake exe cache get`) |
| `lean/.lake/packages` sources (git clones, Mathlib's `.git` 510 MB) | ~770 MB | 10-24 | keep (re-clone needs the network) |
| `artifacts/bin` (`UserRepl.Publish` 168 MB) | 738 MB | 10-29 | regenerable |
| `r/renv/library` | 205 MB | 10-27 | keep (R oracles in use) |
| `artifacts/obj` | 93 MB | 10-29 | regenerable |
| `.git` (16 MB loose in 3,856 objects, 29 MB in 16 packs; no worktrees) | 60 MB | today | keep; gc proposed |
| `artifacts/chip-ext`, `chip-fh3`, `chip-fr1` | 23–198 MB | today | ask (active chips; left alone) |
| `artifacts/test-runs` (`TestHarness.Runner`: 151 runs, 30.8 MB) | 51 MB | today | regenerable runs; rest keep |
| `lean/.lake/build` (the project's own) | 31 MB | 10-26 | regenerable |
| `projects/UserRepl/notebook` (**tracked** pdb/xml/dll) | 28 MB | Jul | keep (tracked) |
| `artifacts/smoke-runs`, `spc`, `gallery`, `repo-audit`, `bench` | 27, 25, 10.5, 7, 2 MB | May–Sep | regenerable / archivable |
| `.codex` (chat-export 14 MB, tmp 8 MB) | 22 MB | Aug | ask |
| `src/.snapshot` (4 payloads) | 5.7 MB | 08-27 | archivable (newest kept) |
| Sessions' temp, `%TEMP%\claude\D--aghado01-ThermoMapper` | 371 MB | today | on **C:**: frees nothing on D: |

**Cleaned on approval** (batches A, B, C; the user said "lets do A B C"). Each path's recency was re-checked just before deletion, and 165 targets were deleted. D: went from 2.59 GB free to 6.72 GB free, **4.13 GB freed** (3.57 GB of file bytes plus cluster slack).
- A, .NET build output: `artifacts/bin` and `artifacts/obj`, 830 MB.
- B, Lean build caches: every `lean/.lake/packages/*/.lake/build` and `lean/.lake/build`, 2,680 MB. To restore: `lake exe cache get` in `lean/` (Mathlib's compressed cache is on C: at `~/.cache/mathlib`), then `lake build`.
- C, old run logs, 65 MB: 150 of 151 `TestHarness.Runner` runs (the newest kept), `smoke-runs/VizCoreSmoke`, `repo-audit` reports from May–June, and `graph_cache`.

**The plan for the rest** (the tool's own `plan`, now the record of it). Archiving waits for the thumb drive's path. Moves go through `apply -ArchiveTo`, which verifies size, SHA-256 and file count before removing a source.
- *Archivable, 41 MB on D:*:
  - `artifacts/spc` (25.2 MB);
  - `artifacts/gallery` (10.5 MB);
  - `artifacts/bench` (2.0 MB);
  - the 3 older `src/.snapshot` payloads (3.7 MB);
  - `artifacts/benchmarks` (46 KB).
- *Regenerable:* `artifacts/doc-dive` (3 KB). Session temp older than a week is on C:.
- *Ask:*
  - `.codex` (22 MB);
  - the loose files at `artifacts/` top level (89 KB: the June `p2-*` files, the tsvs, `README_old.md`);
  - `artifacts/artifacts/` (4 KB; a run from the wrong working directory);
  - the chip roots, which their chips clean.
- *Larger levers, for the user:*
  - `lean/.lake` (~770 MB of package sources) could move to C: (443 GB free) behind a junction; Lean keeps working.
  - `git gc` would pack 3.9k loose objects (~16 MB) and merge 16 packs. Run it when no chip is committing.
  - The 28 MB of tracked binaries in `projects/UserRepl/notebook` could be untracked. That's a git change, and history keeps them anyway.

**The tool, `storage-audit`.**
- *Files:*
  - the source, in `src/storage-audit/`: `PathGlob`, `StorageRules`, `RepoGit`, `StorageProfiler`, `StoragePlanner`, `PlanApplier` and `ReportWriter`;
  - the projects: `projects/StorageAudit` and `tests/StorageAudit.Tests`;
  - the entry point `scripts/storage-audit.ps1`, the rules `presets/storage-audit/rules.jsonc`, and the usage notes `src/storage-audit/README.md`.
- *Commands:*
  - `profile` writes `profile.md` and `profile.json`: sizes by category, drive and top-level entry, the largest items with their age, files no rule claims, git statistics, worktrees, and orphaned `.claude/worktrees` folders.
  - `plan` adds `plan.md` and `plan.json`: one delete per regenerable item, one archive per archivable item, what is held back and why, and git advice that the tool never runs.
  - `apply <plan.json>` is a dry run unless `--execute` is given, and it logs every line to `artifacts/storage-audit/ledger.jsonl`.
- *Rules:*
  - Each rule is a glob with a category and a reason, plus an optional `kind`, `root`, `keepNewest` and `recentWindowHours`. A per-rule window may only lengthen the global one; session temp gets 168 h, because idle sessions resume.
  - A deeper rule carves its item out of a shallower one, and for the same path the first rule listed wins.
  - `tempRoots` are named folders with `%VAR%` and `{repoSlug}` expanded, and they must not overlap the repo. `protect` lists `issues/` and `docs/`; `.git/` is always protected.
- *Refusals* (shared by `plan` and `apply`; `apply` re-profiles before acting):
  - a plan made for another repository is refused whole;
  - path problems: the path escapes its root, names an unlisted root, no rule classifies it, or it no longer exists;
  - the rules now give a different category;
  - the item was written within the window, or holds a tracked file, a file git does not ignore, a link, a protected folder, `.git`, or the running binary;
  - it grew more than 10% (and 1 MB) since the plan;
  - archive problems: no `--archive-to`, an archive inside a root, an existing destination, or a drive too full.
- *During deletion,* each file's write time is re-read, so a file written since the profile is kept and the action reports `partial`. The script runs the tool from a copy in `%TEMP%`, so it can delete `artifacts/bin` without deleting itself.
- *Tests:* 61 xUnit tests on throwaway git repositories under `%TEMP%`, all passing:
  - `PathGlobTests` 6, `StorageRulesTests` 15, `StorageProfilerTests` 13 and `StoragePlannerTests` 6;
  - `PlanApplierTests` 3 (dry run, delete with a carve-out, skip) and `ArchiveTests` 3 (copy, verify, manifest, verifier, archive refusals);
  - `ApplyRefusalTests` 15 (each refusal, the mid-apply write, temp roots).
  
  The run takes about 2 minutes, because each git call takes about 1 s on this machine.

**For the user to decide.**
- The thumb drive's path, for the 41 MB of archivables. They are small; the drive's real lever is outside this repo.
- ~~Whether to move `lean/.lake` to C: behind a junction.~~ Decided: leave it where it is.
- `git gc`, once the chips finish.
- `.codex`, and the loose `artifacts/` files.
- Whether to add `StorageAudit` and `StorageAudit.Tests` to `ThermoMapper.sln`. The solution file wasn't in this chip's list of files.
