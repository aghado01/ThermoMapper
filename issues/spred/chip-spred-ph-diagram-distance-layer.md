I'll work this solo and sequentially per your no-subagent rule (ultracode → depth in the primary agent). Starting by reading the disposition doc and the target file.

Now let me look at the existing tests and the test project layout.

Now Sinkhorn's body and its doc sentence:

Now build and run the test suite:

Done — all three findings are fixed, tested, and committed as `f01949d`, scoped to exactly [DiagramMetrics.cs](src/tda/ph/DiagramMetrics.cs), the two test files, and the disposition doc. TDA.Ph.Tests passes **343/343** (328 pre-existing + 15 new facts). No SPRED-lane files were touched, and the out-of-lane `PersistenceObjectiveDistanceKindTests` uses only finite fixtures, so nothing ripples there.

**Finding 1 (essential births ignored).** Essentials now match by birth in all three backends — with death = ∞ on both sides, the L∞ ground metric between two essentials reduces to |Δbirth|. I chose the assignment route over sorted-order pairing: sorted order is exact for equal counts (1-D transport with convex ground distance is monotone), but under `FinitePenalty` with unequal counts it can't decide _which_ surplus bars go unmatched. So `MatchEssentialBirths` reuses the existing KMN `MinAssignment` on a birth-only matrix padded with zero-slack rows/columns — one exact path for both cases, and essential counts are Betti-number scale so O(k³) costs nothing. Sliced adds the term once, slice-independent (an essential has no finite death to project; its exact transport is already 1-D); Sinkhorn adds it exactly, never smoothed. The shifted-loop pair (1, ∞) vs (100, ∞) now scores exactly 99 in all three backends (pinned by a cross-backend equality test), surplus is still charged at `perBar^p`, and `InfiniteOnMismatch` still returns +∞ on count mismatch.

**Finding 3 (Sinkhorn forbidden cells).** Sentinel cells are masked as log-kernel −∞ inside both LSE sweeps and skipped in the transport sum, so the plan is confined to the diagonal-augmented support at every ε. The discriminating test uses a forced-support fixture — two bars vs empty, where the only admissible cells force the identity plan — so the entropic value must equal the exact 101 at any ε; before the mask it drifted to ≈138 (ε=1) and ≈323 (ε=10), values I verified via the 2×2 closed form. One honest caveat, recorded in the disposition note: "stays close to exact at large ε" is a theorem only on forced-support fixtures — on general diagrams, large ε smears mass across _admissible_ cells (ordinary entropic bias), so the general large-ε test asserts LP feasibility (finite, never below the exact optimum) instead.

**Finding 4 (validation holes).** `ValidateP` rejects NaN; `EssentialPolicy.FinitePenalty` requires finite `perBar` ≥ 0. I allowed zero deliberately — it's a coherent degenerate policy that disables the surplus charge (tested), whereas a negative penalty would reward essential-count mismatch (the exception message says so). Rejection tests cover NaN p and negative/NaN/+∞ perBar.

The one-paragraph result note is appended to the disposition section of [gpt-review-20261717.md](issues/spred/gpt-review-20261717.md). The R transport oracle was not touched or run, per instructions.
