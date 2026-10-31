# storage-audit

Profiles this repository's disk use and the temp data its sessions write, proposes a cleanup
you review, and applies only what you reviewed. Built for the 2026-10-31 disk-space crunch
([brief](../../issues/tooling/storage-audit-brief.md)); it runs on the pattern of `repo-audit`.

## Usage

```powershell
scripts/storage-audit.ps1 profile                         # the survey: profile.md + profile.json
scripts/storage-audit.ps1 plan                            # the survey plus plan.md + plan.json to review
scripts/storage-audit.ps1 apply <plan.json>               # dry run: what apply would do, and what it refuses
scripts/storage-audit.ps1 apply <plan.json> -Execute      # act; add -ArchiveTo <path> for archivable items
```

Reports land in `artifacts/storage-audit/runs/<stamp>/`. Every `apply` line, dry runs included,
goes to the ledger, `artifacts/storage-audit/ledger.jsonl`. Options: `-Rules <file>` for another
rules file, `-ArtifactsRoot <dir>` to build the tool elsewhere (a chip's own root), `-Rebuild`,
`-NoDisplay`. The script runs the tool from a copy under `%TEMP%`, so `apply` can delete
`artifacts/bin` without deleting itself.

Review the plan, set `"enabled": false` on any action in `plan.json` to skip it, dry-run, then
execute.

## Categories

| Category | Meaning | Plan proposes |
|---|---|---|
| `regenerable` | a build output or cache a build or run recreates | delete |
| `archivable` | worth keeping, not here | archive (copy, verify, then remove the source) |
| `keep` | stays | nothing |
| `ask` | unclear | nothing; listed for you to decide |

Files no rule claims are rolled up per top-level entry (tracked, untracked, ignored-no-rule) and
reported only; nothing acts on a rollup. An ignored bulk there is a candidate for a new rule.

## Rules

[`presets/storage-audit/rules.jsonc`](../../presets/storage-audit/rules.jsonc), JSON with comments:

- `rules`: `glob` (relative to its root, `/` separators; `*` and `?` within a segment, `**`
  across segments, a trailing `/` matches folders, otherwise files), `category`, `reason`,
  optional `kind`, `root`, `keepNewest` (keep the newest n of the rule's matches under one
  parent) and `recentWindowHours` (lengthen the window for this rule's items). A deeper rule
  carves its item out of a shallower one; for the same path, the first listed wins.
- `tempRoots`: named folders outside the repository the audit may reach (`%VAR%` and
  `{repoSlug}` expand). Rules name them with `root`. They walk without git's protections, so
  they must not overlap the repository or each other.
- `protect`: repository folders never acted on, nor anything holding them. `.git/` always is.
- `recentWindowHours`: a path written within it is active and left alone (at least 1).

## What apply refuses

`apply` executes a reviewed plan file and nothing else. A plan made for another repository is
refused whole. It profiles afresh with the repository's current rules and re-checks each action,
refusing when:

- the path escapes its root, or names a root the rules do not list (outside the repository, only
  listed temp roots are reachable);
- no rule classifies the path now, or the rules no longer give the action's category;
- anything the item owns was written within the recent window;
- it holds a tracked file, or a repository file git does not ignore (notes, work in progress);
- it holds a junction or symbolic link (never followed, never deleted through);
- it is or holds a protected folder or `.git/`, or the running tool's own binary;
- it grew by more than 10% (and 1 MB) since the plan;
- for an archive: no `--archive-to` when executing, an archive inside a root, a destination that
  already exists, or too little free space on the archive drive.

While deleting, it reads each file's write time again and keeps any written since the profile
(the action reports `partial`). Archiving copies each owned file while hashing it, verifies every copy's size and SHA-256 and
the file count, writes `storage-audit-manifest-<stamp>.json` into the archive, and only then
removes the source. A failed verification keeps the source. Git maintenance (`gc`, worktree
prune, reflog expiry) is proposed in the plan with its numbers and never run.
