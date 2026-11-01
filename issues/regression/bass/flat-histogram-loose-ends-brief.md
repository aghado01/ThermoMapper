# The flat-histogram sampler's loose ends (chip brief)

Date: 2026-10-31. Three items Phase 2 left and one the warm start's report did, for the chip that made the warm start robust ([its brief](flat-histogram-warmstart-brief.md)), to run after that brief is reported. Each lands with its checks and ends with a fresh-eyes review, fixed before the next. Its report goes under "Reports" below, newest first. `docs/Maths` and the checklists are the main session's.

Read first: the design's "Phase 1, as built" and "Phase 2, as built" ([flat-histogram-design.md](flat-histogram-design.md)), Phase 2's "For the main session" ([its report](flat-histogram-phase2-brief.md#reports)), and the warm start's ([its report](flat-histogram-warmstart-brief.md#reports)).

## The items

1. **The window's top.** A pilot that never found Friedman's posterior basin reads the window's top below the posterior's ℓ₁ (−778.3 on seed 8017).
   - *What it costs today:* production put 1–13% of each chain's draws above the window, where the slope-1 extension holds the posterior's own law and the b = 1 weights read them exactly. But the learning spent 8–29% of its visits there, which SAMC's update leaves inert.
   - *What the replications add* (2026-10-31): on seeds 8217 and 8417 the pilots' tops sat higher (−773.6 and −771.7), learning spent 1.9–3.5% of its visits above the window, and the b = 1 Kish ESS fell to 22–32% of seed 8017's. So the top trades learning's inert visits against production's b = 1 draws; read both.
   - *Build* an option, not a default: the window's top read again from the learning's own visits, extended at the same node spacing while the gain is still large, within burn-in (D10). A wider margin above the pilot's b = 1 quantile is the simpler alternative; build it too if cheap.
   - *Check:* on the toys, an option left off is bit for bit today's run; on Friedman seed 8017, read each option against today's window: the learning's visits above the window, production's round trips, the b = 1 Kish ESS, ln Ẑ(1), a1's share and the time. One run each, a reading, recommending no default.
2. **The seed-sweep claims.** Fold `tests/samplers/FlatHistogramSweep.cs` and `tests/regression/BassFlatHistogramSweep.cs` into `ComparativeSweep`, as Phase 1's report proposed.
   - *Check:* the registered claims run under `SEED_SWEEP` as before, and nothing else changes.
3. **SMC's lost mass below the lower split** (a reading). On the linear form, fitting SMC from the prior read ln Z(0.10) about 5.5 nats low against the flat histogram's curve, losing mass between b = 0.05 and 0.10 as its lineages fell to 172–220 of 4,000.
   - D05's SMC path is the integrated form, so read the same stages there: ln Z at b = 0.03, 0.05, 0.07 and 0.10 on both forms (4 replicates, 4,000 particles), against the flat histogram's curve on seed 8017 and exact importance sampling from the prior.
   - Report the families at each b, and whether more mutation steps per stage (5 and 20) close the gap.
   - Recommend no default.

4. **The warm start's readout and adapter** (the warm start's "For the main session").
   - *A `TiedNodes` field* on `BassFlatHistogramPilotReadout`: the nodes the warm start tied, empty where the fixed point accepted ĝ, so a fit's readout says whether its start was tied. Correct `BassDiagnostics`' doc of `NodeDraws`, which says only that an empty node is interpolated.
   - *The adapter's `BassFlatHistogramRun.NodeDraws`* duplicates `Wham.NodeDraws`: call the engine's.
   - *The refusal's message names the failing segment:* move the hull reading (`FlatHistogramLimit.Infeasibility`, test code) into the engine as a diagnostic, and have a refusal say at which segment flat occupancy leaves the hull, and by how much, rather than only how far Newton strayed.
   - *Check:* `BassBitIdentityTests`' hashes and every pinned start hold; the field on a tied start (seed 8217's pilot, or a synthetic edge pilot) and on an accepted one; the message on the gap pilots and on a pilot too thin to tie.

## Working rules

As the warm-start brief's, with `artifacts\chip-fh4\` as the artifacts root. Your files are those of the warm-start brief, plus `tests/samplers/ComparativeSweep.cs` and the two sweep files above, and for item 4 `src/maths/regression/spline/bass/BassDiagnostics.cs` (the readout's field and doc) and `BassFlatHistogramRun.cs` (the count and the field's fill).

## Reports

Newest first.
