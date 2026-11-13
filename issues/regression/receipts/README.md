# Regression receipts

Tracked copies of the small reproduction inputs behind the verification claims in the
regression briefs and reports. Each file sits at the same relative path it has under the
git-ignored `artifacts/` folder, so `artifacts/<path>` maps to `receipts/<path>`.

What is here: run settings, baseline and commit pins, source manifests, live-state
records, patches and probe sources used in investigations, and small result summaries
(provenance, coverage tables, coverage figure).

What stays local in `artifacts/`: TRX logs, console logs, source tarballs, the slice
implementation patch (git history holds the source) and bulky raw outputs such as
`e34-slice3/results/coverage.csv`. The documents quote the counts from those logs;
rerunning a tier with the tracked run settings on the pinned commit reproduces them.

The `verify-*` folders hold what a coordinating session's verification ran on. The scripts that
produce them are in `scripts/`: `verify-csharp.ps1` (the Commit selection from a clean archive of
HEAD) and `verify-lean.ps1` (the checks, module build, axioms and meta-ci from a clean archive of a
commit's `lean/`). Their headers give the usage; both write their receipts under `artifacts/`.
