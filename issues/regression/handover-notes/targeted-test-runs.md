# Targeted test runs (2026-11-07)

Routine runs of Maths.Regression.Tests no longer need the whole assembly. The rule is in CONTRIBUTING.md, "Targeted runs", and the GBASS, KHAOS and conformal briefs' common-gate and working-rule lines point at it.

- Per edit: `scripts/test-select.ps1 -Changed <your paths> -ArtifactsRoot <your root>`. That is the classes your change can reach, without the slow tier.
- Per commit: add `-Scope Commit`, which also runs the pinned tier, `[Trait("Tier", "Pinned")]` (31 cases, about 1 min).
- At a slice's end: `-Scope Slice`, the whole assembly once.
- Pass your own paths with `-Changed` (commas from bash). The default takes in everyone's uncommitted work.
- Building from a snapshot: add `-SourceRoot <snapshot>` and `-NoBuild` after your build.
- A new test that pins a hash or a fingerprint takes `[Trait("Tier", "Pinned")]`.

Commits: 39c1008 (traits), 42de39d (script), 82e91eb (slow list), 1d8ffa1 (docs). The report is issues/tooling/targeted-test-runs.md.
