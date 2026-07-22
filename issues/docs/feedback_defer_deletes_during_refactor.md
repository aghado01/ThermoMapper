---
name: feedback-defer-deletes-during-refactor
description: "During multi-step refactors, write new files alongside the old ones and defer deletion of superseded code to a dedicated cleanup pass at the end; git history is the safety net"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ef8e192e-720a-4222-addb-e184edf1d3ba
---

When executing a multi-step refactor in ps.core.pwshspc (and other
pre-release projects), write the new files first and **leave the old
ones in place** until a dedicated cleanup sweep at the end of the
refactor. Do not delete superseded files in the same pass that
introduces their replacements.

**Why:** The user said explicitly: "instead of deleting old files,
why don't you just write the new ones and we can clean up old files
later? nothing to lose." The reasoning:

1. **The build never breaks mid-flight** — callers that haven't been
   migrated yet still find what they need. Each in-progress phase
   stays buildable, which keeps the loop tight.
2. **Reversibility within a phase** — if the new code has a bug, the
   old code is right there for A/B comparison or to fall back on
   without consulting git history.
3. **Separation of concerns in commit history** — "build new shape"
   and "delete old shape" land as separate commits, which makes
   review and post-hoc archaeology cleaner.
4. **Git history is the safety net** — deleted files can always be
   restored from history later, so "nothing to lose" by keeping them
   around through the refactor.

**How to apply:**
- Add new files in their permanent locations (no transitional dirs)
- Migrate callers one-by-one to the new surface as needed
- Resist the urge to `rm` or `git rm` the superseded files until the
  entire refactor is feature-complete and verified
- The final "cleanup" task in a refactor plan should be the *only*
  place that does file deletions — and it should be its own discrete
  commit

This pairs with [[feedback_no_compat_bandaids]] — that one says
*end state* has no `[Obsolete]` shims or migration paths; this one
says *path to the end state* keeps the old code alive in tree until
the sweep. End state and journey, separate principles.
