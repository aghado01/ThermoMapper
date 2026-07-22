---
name: feedback-no-compat-bandaids
description: "Pre-release codebase — no backwards-compat shims, no [Obsolete] markers, no migration windows; delete old surfaces wholesale when refactoring"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ef8e192e-720a-4222-addb-e184edf1d3ba
---

When refactoring APIs in ps.core.pwshspc (and the user's other in-development
projects), DELETE the old surface wholesale and migrate every caller in the
same pass. Do not add `[Obsolete]` attributes, do not leave shims that route
old calls into the new structure, do not add schema-version migration paths
for serialized formats.

**Why:** The user has not released or seriously used any of this code yet —
nothing in the wild depends on the current shape. The user explicitly named
the failure mode this avoids: "I keep chasing bandaids and compromises."
Compat shims here are pure cost — they accrete, they delay the real cleanup,
and they make the new shape harder to read because the old shape is still
present.

**How to apply:**
- When introducing a new API to replace an old one, delete the old one in the
  same change (or in the immediately-following commit), don't `[Obsolete]` it
- When changing on-disk serialization shapes (RunManifest, presets, etc.),
  rewrite the schema without versioning logic — no SchemaVersion 1.0 → 1.1
  branch handling
- When moving files (e.g. MRPGraph.cs → primitives/mst/Prim.cs), don't leave
  a stub at the old path that re-exports — just move it and let the compile
  errors guide caller updates
- When deleting helpers whose only callers used the old API, delete them too
  rather than leaving them as "potentially useful primitives"
- Update changelog.md to describe the clean break, but don't apologize for
  the breakage — it's intentional, not collateral damage

This applies as long as the project is pre-release. If/when a release ships
and external consumers appear, revisit. See also [[feedback_design_before_build]]
and [[feedback_fluent_apis]].
