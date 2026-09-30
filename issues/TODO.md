# TODO — canonical documentation shore-up

Organizing and shoring up canonical project documentation files (AGENTS.md,
CLAUDE.md, README.md, CONTRIBUTING.md), with additional ones to be added
later.

## Done (2026-07-21)

- [x] Delineated the doc responsibilities and codified the staged feedback
      memories (`issues/docs/feedback_*.md`) into the right homes:
  - **AGENTS.md** — project orientation (root-directory map), layering
    invariants, tooling (repo-audit, test-harness), environment (pdenv
    registry PATH, stale-snapshot caveat), agent behavioral guidance,
    DocOps conventions for `issues/`.
  - **CONTRIBUTING.md** — source of truth for development practices:
    engine-first strategy, strict-core/fluent-shell, refactor discipline
    (clean breaks + deferred deletes), naming & vocabulary principles,
    docstring conventions, testing standards, antipatterns.
  - **CLAUDE.md** — thin pointer to the two above (resolution of the "what
    is Claude-specific?" question: nothing yet; agent-agnostic guidance
    belongs in the shared docs).
- [x] DocOps first pass — canonical `issues/` file kinds (briefs, reports,
      design notes, digests, SHELF) codified in AGENTS.md. Deepen as the
      conventions firm up.

Resolved en route: `feedback_fluent_apis`'s `[Obsolete]`-shim migration-window
advice is superseded by no-compat-bandaids + strict-core/fluent-shell; only
its "param lists encoding stages → make stages first-class" signal was kept.

## Done (2026-09-15)

- [x] README landing page with a code map, EyeTorusToy portrait, and a working
      fixture-render command. Camera settings and five screenshots are in
      [the gallery](docs/gallery.md).

## Remaining

- [ ] README.md — author's conceptual overview and broader vision. The current
      landing page covers the code and fixture gallery.
- [ ] Convenience aliases for R / dotnet / Lean etc. — add PowerShell loader
      scripts (activate local python venv, load aliases), then surface them
      in AGENTS.md § Environment.
- [ ] Language-specific development guidelines — e.g. distillations of C#
      SDK guidelines; likely `issues/docs/` or a `docs/` subtree, referenced
      from CONTRIBUTING.md.
- [ ] Decide the fate of the staged memory files in `issues/docs/` now that
      they're codified (keep as provenance, or sweep once confident).
