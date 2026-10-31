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
