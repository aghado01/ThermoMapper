# The flat-histogram warm start, robust across seeds, and Friedman replicated (chip brief)

Date: 2026-10-30. N17 Phase 2 found Friedman's posterior on one seed ([its brief and reports](flat-histogram-phase2-brief.md#reports)). On master seed 8017, from the WHAM warm start, every chain tunnelled, the chains agreed within 0.34 nats, and ln Ẑ(1) = −841.4 ± 0.2 reached the −864.6 bound, with the posterior on a1 (four x1·x2 terms).
- *A second seed (8217) failed at its start.* SAMC's fixed point refused the pilot's ĝ over [−1552.56, −773.559], K = 128, b_lo = 0. No node was empty, but a node's occupancy stood 125 from flat, against the 10⁻³ the solve must reach.
- *A flat start there, at the same budget, learned nothing:* no round trip in 100,000 iterations, b = 1 Kish ESS 8 of 500,000.
- *What this chip does:*
  - make the warm start hold across seeds, without changing a start the fixed point already accepts;
  - replicate the Friedman reading on seed 8217 as pre-registered, and on a third seed;
  - read whether a flat start can learn Friedman's window with a larger gain scale.
- *Not in scope:*
  - any default change;
  - the window's top, which a pilot that never found a1 sets below the posterior's ℓ₁ (a later option);
  - SPC (Phase 3);
  - `docs/Maths` and the checklists, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- [the design](flat-histogram-design.md), with its "Phase 1, as built" and "Phase 2, as built";
- [Phase 1's brief](flat-histogram-phase1-brief.md), its reports (the `FixedPoint` fix `9f561d0` above all);
- [Phase 2's brief](flat-histogram-phase2-brief.md), its reports and the main session's verification;
- `docs/Maths/Samplers.md`, "Flat-histogram sampling", and `docs/Maths/BayesianRegression.md`, "BASS: the flat-histogram sampler".

In the code:
- `src/maths/samplers/flathistogram/SamcLearning.cs` (`FixedPoint`: the free nodes, the empty-node rule, the start fitted per node cell, Newton, the stall rule and the 10⁻³ refusal);
- `Wham.cs` (`Estimate`, `WarmStart`, `Window`, `Nodes`) and `TemperingPilot.cs`;
- `src/maths/regression/spline/bass/BassFlatHistogramRun.cs` (the pilot, `Pooled`, the start's shape, the refusal's message);
- `tests/regression/BassFlatHistogramQualificationTests.cs`: arm (i), `Friedman_ArmOne_SecondSeed_TwiceTheSweeps` (seed 8217, `Thin` 400, its criteria written before its run) and `Friedman_ArmOne_SecondSeed_FlatStart`.

## The steps

Each step lands with its checks and ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next. Name the passes and their fixes in the report.

1. **Diagnose seed 8217's refusal.**
   - Rebuild its pilot (the pilot is deterministic from the master seed; it takes about a minute), then WHAM's ĝ and the start's shape.
   - Find the node that departs and why:
     - its hat's mass and where in its cell that mass sits;
     - its neighbours' masses;
     - the pilot's draws per segment;
     - Newton's path: the iterations, the step sizes, where it stalled or ran out.
   - Say which of these it is, with the numbers:
     - an ill-posed limit, a hat seeing mass only at an edge of its support, so θ* is extreme;
     - a well-posed limit Newton fails to reach from its start;
     - something else.
   - Do the same for seeds 7801–7803, where Phase 2's exploration also saw refusals on fixture A.
2. **Make the start robust.**
   - *Build* the remedy the diagnosis supports. The candidates:
     - widen the empty-node rule to thin nodes, a hat whose effective pilot mass falls below a stated threshold, interpolated as empty nodes are;
     - interpolate ln ĝ log-linearly across segments the pilot left empty, before the fixed point;
     - coarsen K until the fixed point flattens;
     - a better Newton start or line search, if the limit is well posed.
   - Choose by measured robustness, and state the rule and its threshold.
   - *The invariant:* a start the fixed point accepts today is unchanged, bit for bit. Check it on seed 8017's Friedman pilot, by the start's weights' hash taken before the change, and on every Phase 1 and Phase 2 test.
     - So the remedy runs where the fixed point refuses, or under a rule that no accepted start triggers.
     - Seed 8017's accepted pilot already had nodes under 10 draws (0–2 and 153–155), so a thin-node threshold must not reach them, or must apply only on refusal.
   - *Check:*
     - on a sweep of Friedman pilots (at least 16 master seeds, declared before the sweep, 8017 and 8217 among them), the refusal rate before and after;
     - on the toys and fixture A, the start's accuracy against SAMC's limit where it was accepted, unchanged, and where it was refused, now within Phase 1's tolerance (0.25 RMS);
     - the refusal kept for a ĝ with no usable mass, with its message.
3. **Replicate Friedman** (`BASS_QUALIFICATION` = 1; S10's Friedman fixture; Release; one fit at a time).
   - *Seed 8217* by `Friedman_ArmOne_SecondSeed_TwiceTheSweeps` as written: the WHAM start, the default gain start, `Thin` 400, its criteria as committed in `56c66f3`, unchanged.
   - *A third seed,* declared before any run and not used by any harness: arm (i) at the same settings, with the same criteria, in a test written and committed before it runs.
   - *What would count,* for each seed, as committed:
     - every production chain tunnels;
     - the chains' ln Ẑ_c(1) SD is at most 1 nat, and each labelled basin at 5% or more passes the χ² at p ≥ 0.01;
     - the pooled ln Ẑ(1) is at or above −864.6;
     - against seed 8017, ln Ẑ(1) lies within 1 nat plus both runs' errors of −841.36 ± 0.17, and a1's pooled share within 5 points of 65.9%.
   - *A reading, not ruled:* the flat start on seed 8217 over its pilot's window with t₀ = 100·K (Liang et al.'s range runs to 100K), the rest as `Friedman_ArmOne_SecondSeed_FlatStart`. It reads whether a flat start can learn the window when the warm start is unavailable.
   - *Read per run,* by the harness's `Read`: the window, K, t₀ and the gain's start; the pilot and the fixed point; the learning's round trips, flatness and trace; production's round trips, flatness and b = 1 Kish ESS; ln Ẑ(1) per chain and pooled; the basins with their χ²; σ's interval; D10's verdict; the time.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Other sessions may be working in this checkout,** an SPC session among them. If a build fails in a file you do not own, wait and retry; never edit it.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-fh3\` to every `dotnet build` and `dotnet test`. Build in Release for the runs. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-fh3\` when you finish.
- **Your files:**
  - `src/maths/samplers/flathistogram/` (the warm start's code) and its tests under `tests/samplers/`;
  - `BassFlatHistogramRun.cs`, only where the start's shape or message must follow the change;
  - `tests/regression/BassFlatHistogramQualificationTests.cs` and `BassFlatHistogramTests.cs`;
  - this brief's Reports section.
- **Commit your own files only, by explicit path,** one commit per step, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked notes under `issues/regression/` and other sessions' edits stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
  - Commit each run's test before running it.
- **Line endings:** briefs and C# are LF.
- **Bit-identity:** `BassBitIdentityTests`' hashes hold, and so does every start the fixed point accepts today.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - the remedy leaves an accepted start bit for bit as it was;
  - a thin or interpolated node's weight is a start only, which SAMC's learning corrects;
  - the refusal kept where nothing usable remains;
  - each run's criteria committed before it ran.
- *Coverage:*
  - a pilot with an empty segment beside a drawn node;
  - a thin node at the window's edge;
  - a pilot whose ĝ spans more than 745 nats;
  - a seed whose replication fails a criterion, reported as it falls.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the diagnosis, with its numbers;
- the remedy and its rule, with the refusal rates before and after;
- each replication: what would count as committed, then its tables, then what they show, kept apart;
- the flat-start reading;
- the suites with their counts (Samplers, Regression, the BASS oracle classes);
- the review passes and their fixes;
- anything the main session must know.

Recommend no default.

## Reports

Newest first.

### The warm start made robust, and Friedman replicated (chip, 2026-10-31)

**What landed**, on `main`, nothing pushed:
- `694dfa9`, step 1: the diagnosis's checks. `FlatHistogramLimit.Infeasibility` (test code) reads whether flat hat occupancy lies inside the hull of a discrete measure's hat vectors; the gap pilots on fixture A are refused exactly where it does not, and seed 8217's Friedman diagnosis, written before its run, rebuilds the second arm's two pilots.
- `4f0dd6f`, step 2: the hashes of every start the fixed point accepted, pinned on the unchanged engine before the remedy.
- `46332d7`, step 2: the remedy and its checks, and the Friedman sweep's harness, written before it ran, with the third seed declared.
- `ca091f8`, step 2: the sweep read one seed a case (its first run lost its test host about 45 minutes in, before any row was written; no crash was logged).
- `7fc684e`, step 3: the third seed's test and the flat-start reading, written before they ran.
- Source files changed: `SamcLearning.cs` and `Wham.cs` only. `MulticanonicalRun` and `BassFlatHistogramRun.cs` are untouched; the adapter's refusal message already names the window, K, b_lo and the empty nodes, and carries the warm start's reason.

**The diagnosis** (`GapPilots_AreRefused_ExactlyWhereNoLimitFlattensTheirHats`, `Friedman_SecondSeed_WhyTheFixedPointRefusedItsPilot`; Newton's path from a traced copy of `FixedPoint`, run on the dumped ĝ and since deleted).
- *What fails is the limit, not Newton.* SAMC's limit on a discrete measure exists only where flat hat occupancy lies inside the convex hull of the atoms' hat vectors φ(ℓ_j), the means its tilts reach. Read along the window, the hull's cumulative mass x_m must satisfy x_m + (1 − u_max)P_m ≤ C_m ≤ x_m + (1 − u_min)P_m for each segment's atoms, C_m the cumulative target. Where it fails, G is unbounded below and no finite θ* exists.
- *Seed 8217*, window [−1552.56, −773.559], K 128, spacing 6.134, b_lo 0, 78 pooled rungs, 168,000 draws (168,000 atoms), WHAM residual 1.6e−9:
  - node 0 holds 0.1 pilot draws under its hat: one atom, at u = 0.869 of segment 0 (ln mass under the hat −11.04). Node 1 holds 3.0 (segment 1: 8 atoms, u 0.254–0.976), node 2 28.8, node 3 179.2.
  - Node 0's flat share, 1/254, would need that atom to carry 1/(254 × 0.131) = 3.0% of the tilted mass, which hands node 1 2.6%, past its whole share of 0.79%: flat occupancy lies outside the hull at segment 1, by 0.0183.
  - Newton from the start (worst departure 10.1 at node 127): steps t = 0.125, 1, 1, 1 with |direction|∞ 511, 2,070, 3.9e5, 8.8e11; G −21.0, −21.2, −24.0, −958, −2.1e9; then no descent in 40 halvings at |direction|∞ 2.6e17. Every θ but node 0's had run to about 7.7e11; the node "125 from flat" (node 3) is where the mass ran to, not the cause.
  - So: an ill-posed limit, a hat seeing ĝ only at the far edge of its support. A better Newton start or line search cannot help.
- *Seed 8017*, accepted: node 0 holds 0.4 draws, segment 0 two atoms at u 0.634 and 0.999, segment 1 none; flat occupancy inside the hull. It was accepted by where its few atoms fell, not by a margin.
- *Fixture A's gap pilots*, K = 20 over [−37.03, −15.80], τ held (about 130 atoms):
  - 7801 refused (17.1 from flat): node 13 holds 0.1 draws, one atom at u = 0.873 between empty segments 12 and 14; outside the hull at segment 13, by 0.0736. Newton's first step had |direction|∞ 6.2e13.
  - 7802 refused (13.9): nodes 13 and 14 empty, node 15 holds 0.2 draws, one atom at u = 0.761; outside at segment 16, by 0.1273. Steps 387, 7,381, then 1.9e25.
  - 7803 accepted, inside the hull (node 12 sees 34 draws at u 0.666–0.987 beside the gap).
- *The same at the window's edges, generally.* Fixture A, τ sampled, 48 pilots over the shapes they resolve themselves (K 39–44): 7 refused, each outside the hull at segment 1 or 2, the hot rung's tail the window's margin reaches. The double well's sparse pilots over their own shapes likewise: 28 of 32 refused at 5,000 scans, all 32 at 2,000.

**The remedy, its rule, and the refusal rates.**
- *The rule.* Where SAMC's fixed point refuses ĝ as it stands, the nodes holding fewer than 10 of the pilot's draws under their hat (`Wham.NodeDraws`, the readout's own count) are tied: a tied node's ln w is held on the line through its nearest held nodes, and past the outermost held node continues at the extensions' slopes (b_lo below, 1 above). The limit on that coarser interpolant (knot removal: the merged hats brought to their own flat occupancy, `SamcLearning.FixedPoint(shape, levels, logMasses, tied)`) is the start; SAMC's learning then moves every node from it. Where even that is refused, or fewer than two nodes hold ten draws, the start is refused, saying which. `Wham.WarmStart(pilot, shape, thinDraws, out tied)` exposes the threshold and the nodes tied; `Wham.ThinDraws` = 10.
- *The invariant.* The untied solve runs first, its arithmetic unchanged (moved into `Solve` and `LogLinearStart`, nothing reordered), and its accepted starts are returned as they were. Pinned before the change (`4f0dd6f`) and holding after: the double well's 16 starts `e3aa60d4d15c7b42`, the Potts lattice's 8 `1d49177c45e38b17`, fixture A's 8 pilots `95ad874b5dce0d4e`, its 41 accepted of 48 `e09cb5b3b470f47e`, the gap pilot 7803 `33f442de917867dc`, Friedman seed 8017 `86f45b016b3f3f47`. The sweep asserts the warm start equals the fixed point bit for bit on its 14 accepted seeds; `BassBitIdentityTests`' hashes hold.
- *Chosen by measured robustness*, on the refused pilots read before the Friedman sweep (fixture A's 7, its gap pilots at K = 20, 30 and 40 (8 refused), the double well at 2,000 and 5,000 scans (60 refused), Friedman 8217), still refused after each candidate:

  | candidate | double well, 2k / 5k scans | gap pilots, K 20 / 30 / 40 | fixture A's 7 | start's RMS to θ* (fixture A's 7) |
  | --- | --- | --- | --- | --- |
  | empty the thin nodes (atoms set aside, the empty-node rule), 1 / 2 draws | 25 / 19 of 32 at 2k; 5k not read | K 20: 0; K 40: 3 of 3 (at 1 and 3 draws); K 30 not read | 1 / 0 (4 at 0.5 draws) | 0.20–0.78 (0.19–0.38 on the drawn nodes) |
  | coarsen K by halving | not read | 0 at K 20 and 40 (down to K 6–21); K 30 not read | 0 (K 12–23) | 0.33–1.18 |
  | tie, 1 draw | 12 / 11 | 0 / 3 / 3 | 0 | 0.20–0.43 |
  | tie, 2 draws | 0 / 2 | 0 / 3 / 3 | 0 | 0.20–0.42 |
  | tie, 5 draws | 0 / 0 | 0 / 3 / 1 | 0 | 0.20–0.31 |
  | tie, 10 draws | 0 / 0 | 0 / 0 / 0 | 0 | 0.22–0.34 |

  Emptying makes the thin nodes' neighbours one-sided, with an interior node's target, the configuration that refused in the first place. Coarsening always ends accepted but loses resolution everywhere for a local fault. Tying merges the thin hats into their neighbours, so the atoms stay and the coarser hats see them; 10 draws is the least threshold measured at which every refused pilot was accepted. The tie measurements above were made with the tied edges held level; on the extensions' slopes, as adopted, feasibility is the same (only fixed offsets change) and fixture A's numbers are identical (its ties are at the bottom edge, b_lo = 0), while a smooth Gaussian ln g read 13.8 RMS level and 1.1 on the slopes, and the double well 0.26 and 0.32 (its top is ℓ's maximum, where g is singular).
- *The Friedman sweep* (`Friedman_WarmStart_RefusalsBeforeAndAfter_OverSixteenPilots`; 16 master seeds declared before it: 8017, 8217 and 8301–8314; arm (i)'s pilots over the shapes their runs resolve, K 104–204):
  - *before*, the fixed point refused 2 of 16: 8217 (125 from flat; outside the hull at segment 1, by 0.0183) and 8313 (167; at segment 2, by 0.0011);
  - *after*, the warm start refused none: both tie nodes 0 and 1. The 14 accepted keep their starts bit for bit.
  - Readings on the two refused: ties at 1 draw refused 8313; at 2 and 5, both accepted; emptying at 2 draws refused 8313; coarsening accepted at K 65 and 88.
  - Every seed's window starts thin: node 0 under 10 draws on all 16, nodes 0 and 1 on 10, and three seeds had a thin interior node (8017's 153–155, 8301's 103, 8305's 106).
- *Fixture A and the toys*:
  - fixture A, 48 pilots: 7 refused before, none after; the 41 accepted unchanged;
  - the gap pilots: 7801 and 7802 now run (7802's fit completes with its evidence, tied [0, 7, 10, 11, 13, 14, 15, 16]); 7803 unchanged;
  - the double well's 32 sparse pilots: 28 refused before, none after; the 4 accepted unchanged.
- *Accuracy where it was refused*, the centred RMS of the start to SAMC's limit on the exact density of states, pooled: fixture A's 7 tied 0.290, against 0.838 for its 41 accepted as they stood (median 0.35; the accepted starts' window-edge nodes sit up to 17 nats off); the double well's 28 tied 0.323, against 0.709 for its 4 accepted. *Phase 1's 0.25 is not met*, by the tied starts or by the accepted ones on the same shapes: it was read on denser pilots over a narrower window. The bounds the tests assert (0.35) were set after these readings.
- *The refusal kept*: a pilot too thin for two nodes to hold ten draws is refused, its message saying so (fixture A's gap pilot traced 3 scans, 18 draws over 20 nodes, through `Bass.Fit`, naming the window, K, b_lo and the empty nodes; a synthetic four-draw pilot). A Gaussian pilot of 10 draws a rung (210 over 48–56 nodes) is refused the same way: at this threshold the rule needs a pilot of a few thousand draws.

**Seed 8217, replicated** (`Friedman_ArmOne_SecondSeed_TwiceTheSweeps` as committed in `56c66f3`; Release, alone on the fixture).

*What would count, as committed:* (1) every production chain tunnels; the chains' ln Ẑ_c(1) SD at most 1 nat, a1 (or any labelled basin at 5% or more) passing the χ² at p ≥ 0.01; the pooled ln Ẑ(1) at or above −864.6; (2) against seed 8017, ln Ẑ(1) within 1 nat plus both runs' errors of −841.36 ± 0.17, and a1's pooled share within 5 points of 65.9%.

| | seed 8217 | seed 8417 | seed 8017 (Phase 2, Thin 200) |
| --- | --- | --- | --- |
| window (width), K, b_lo, t₀ | [−1552.6, −773.6] (779), 128, 0, 1,280 | [−1550.9, −771.7] (779), 124, 0, 1,240 | [−1551.0, −778.3] (773), 204, 0, 2,040 |
| the warm start | the fixed point refused ĝ; nodes 0 and 1 tied (the sweep's row for this pilot; the run's start level over nodes 0–2, as tied at b_lo = 0) | accepted as it stood (nodes 0 and 1 under 10 draws) | accepted as it stood |
| the gain counted from | 42,000 | 42,500 | 42,000 |
| pilot: rungs per ladder; pooled; draws; round trips per traced scan; barrier; worst pair rejection; g | 21 × 4; 78; 168,000; 0; 8.52; 1.000; 1.1–1,244 | 21, 21, 21, 22; 79; 170,000; 0; 8.17; 0.993; 1.0–1,513 | 21 × 4; 78; 168,000; 0; 8.30; 0.975; 1.1–1,345 |
| learning: round trips per chain; flatness; visits above the window | 11, 10, 12, 13; 2.21; 14,184 of 400,000 (3.5%) | 11, 7, 11, 12; 2.27; 7,746 of 400,000 (1.9%) | 5, 9, 5, 4; 2.80; 29% |
| the weights' centred RMS change, start–10%, 10–50%, 50–90%, 90–100% | 3.68, 2.02, 0.13, 0.069 | 3.97, 1.76, 0.11, 0.063 | 0.88, 3.42, 0.15, 0.093 |
| production: round trips per chain (all) | 17, 17, 14, 14 (62) | 6, 16, 21, 14 (57) | 17, 9, 12, 9 (47) |
| production flatness, pooled (per chain) | 0.58 (1.05, 0.48, 0.35, 0.46) | 0.34 (0.60, 0.85, 0.74, 0.34) | 0.39 |
| production above the window, per chain of 125,000 | 1,970; 1,724; 1,464; 1,650 | 311; 1,173; 986; 757 | 8,500; 4,889; 8,970; 9,726 |
| b = 1 Kish ESS of 500,000 (per chain) | 10,937 (3,300; 2,703; 2,345; 2,590) | 7,684 (702; 2,650; 2,466; 1,867) | 34,256 |
| per chain ln Ẑ_c(1) | −841.37, −841.57, −841.74, −841.54 | −842.16, −840.78, −840.86, −841.21 | −841.30, −841.90, −841.20, −841.19 |
| their SD (range) | 0.15 (0.37) | 0.63 (1.38) | 0.34 (0.71) |
| pooled ln Ẑ(1) | −841.55 ± 0.07 | −841.14 ± 0.32 | −841.36 ± 0.17 |
| ln Ẑ(b) at 0.1, 0.2, 0.25, 0.3, 0.4, 0.6, 0.8 | −126.0, −210.3, −251.4, −292.0, −372.4, −530.8, −686.9 | −125.9, −210.2, −251.3, −291.9, −372.3, −530.6, −686.5 | −125.9, −210.3, −251.4, −292.0, −372.3, −530.8, −686.7 |
| σ's 95% interval | [0.9765, 1.1125] | [0.9774, 1.1143] | [0.9757, 1.1139] |
| D10's verdict on π_w's functionals | converged, every group | converged, every group | converged, every group |
| time: pilot, learning, production, all; sweeps a second a chain | 53, 1,013, 1,292, 2,358 s; 39,000 | 39, 934, 1,213, 2,186 s; 41,900 | 67, 337, 375, 779 s; 63,000 |

| basin | seed 8217: pooled (log mass); chains; χ² p | seed 8417 |
| --- | --- | --- |
| a1, 4 × x1·x2 + 2 × x3 + x4 + x5 | 66.6% (−842.0); 71.1, 58.5, 68.3, 67.8%; p 0.613 | 68.4% (−841.5); 62.7, 68.4, 67.5, 71.9%; p 0.848 |
| a1 + x6·x8 (unlabelled) | 9.8% (−843.9); 9.2, 9.9, 8.8, 11.3%; p 0.687 | 10.1% (−843.4); 13.8, 11.4, 10.3, 6.7%; p 0.000 (unlabelled, outside the criteria) |
| B | 0.9% (−846.2) | 0.5% (−846.4) |
| a2 | 0.4% (−847.2) | 0.2% (−847.4) |
| A | under 0.05% (−851.2), 6,954 draws visiting | under 0.05% (−851.4), 6,327 draws visiting |
| T, a5; G, S | under 0.05% (T 33 draws, a5 482); none | under 0.05% (T 18 draws, a5 453); G none, S 1 draw |
| x1·x2 terms 3, 4, 5 | 3.3, 96.6, 0.1% (each chain 95.6–97.5% at 4) | 2.3, 97.3, 0.4% (each chain 96.2–97.8% at 4) |

The rest of 8217's weight sits on a1 with one term added (x2·x9·x10 2.3%, x6·x8·x10 1.3%, x8 1.0%, a third x3 0.9%) or one x1·x2 dropped (0.8%). Every chain visits a1 and a2.

*Against what would count (seed 8217):* (1) every chain tunnels, fewest 14 round trips; the chains' SD is 0.15 nats; a1, the one labelled basin at 5% or more, passes at p 0.613; ln Ẑ(1) −841.55 is above −864.6. (2) It lies 0.19 nats from −841.36, within 1 + 0.07 + 0.17 = 1.24; a1's pooled share, 66.6%, is 0.7 points from 65.9%. Every criterion is met.

**Seed 8417, a third seed** (`Friedman_ArmOne_ThirdSeed_TwiceTheSweeps`; declared in `46332d7`, its test committed in `7fc684e` before it ran; its columns are in the tables above).
- *What would count:* the second seed's criteria, as above, unchanged.
- *Its warm start*, read before the run from the pilot rebuilt from the seed: the fixed point accepted ĝ as it stood (nodes 0 and 1 under ten draws, none tied), and the run's start equals it bit for bit (asserted).
- The rest of its weight sits on a1 with one term added (x2·x9·x10 1.6%, x6·x8·x10 1.2%, x8 1.2%, a third x3 1.0%, a second x5 0.8%). Every chain visits a1 and a2.

*Against what would count (seed 8417):* (1) every chain tunnels, fewest 6 round trips (chain 0, its b = 1 ESS 702); the chains' SD is 0.63 nats; a1 passes at p 0.848; ln Ẑ(1) −841.14 is above −864.6. (2) It lies 0.22 nats from −841.36, within 1 + 0.32 + 0.17 = 1.49; a1's share, 68.4%, is 2.5 points from 65.9%. Every criterion is met. Beside them, the unlabelled a1 + x6·x8 holds 10.1% and its chains disagree (13.8 to 6.7%, χ² p < 0.001 against their batch errors); the criteria read labelled basins only, so it is reported here, not ruled.

*What the replications show* (kept apart from the tables):
- *Friedman's reading now stands on three master seeds.* ln Ẑ(1) is −841.36 ± 0.17, −841.55 ± 0.07 and −841.14 ± 0.32 on 8017, 8217 and 8417, within 0.41 nats of one another and 23 nats above the 99.9% bound. a1 holds 65.9, 66.6 and 68.4%, a1 + x6·x8 10.6, 9.8 and 10.1%, four x1·x2 terms 96.0, 96.6 and 97.3%; A, B, G, S and T together 1.3% or less on each. The ln Ẑ(b) curves agree within 0.4 nats at every b read.
- *The remedy's start reached the same posterior as an accepted start*, on the one seed where it fired: 8217's start tied nodes 0 and 1, and its run agrees with 8017's and 8417's, whose starts the fixed point accepted. One seed is all the remedy's start has here.
- *Twice the sweeps bought more crossings, not 20 a chain.* 62 and 57 production round trips against 8017's 47 at half the sweeps, but per chain 14–17 and 6–21: one chain of the eight reached the 20 the lengths were sized for.
- *The b = 1 Kish ESS fell*, 10,937 and 7,684 of 500,000 against 8017's 34,256. A reading: these windows' tops (−773.6 and −771.7) sit above 8017's (−778.3), so fewer production draws lay above the window, where the extension holds the posterior's own law (0.2–1.6% a chain against 4–8%), and more of the posterior sits inside the window at a flat share.
- *8417's chains spread more* (SD 0.63 against 0.15 and 0.34): its chain 0 made 6 round trips, at b = 1 ESS 702, and sits a nat below the others. The same chain holds the most a1 + x6·x8 (13.8%).
- σ's interval covers 1, and D10's verdict on π_w's functionals is converged in every group, on both runs.
- *Strength:* three master seeds, four chains each; each run's evidence error rests on its chains' spread, 3 degrees of freedom.

**The flat start at t₀ = 100·K** (a reading, not ruled; `Friedman_ArmOne_SecondSeed_FlatStart_GainScaleOneHundredK`):
- *As run:* seed 8217's pilot window [−1552.6, −773.6] (779 nats), K 128, b_lo 0, a flat start with the gain from 0 and t₀ = 100·K = 12,800, the second seed's lengths and Thin 400 (90M sweeps a chain); no pilot, no fixed point.
- *Learning:* round trips per chain 13, 12, 12, 16; flatness 0.12; visits below / above the window 2,328 / 2,871 of 400,000; the weights' centred RMS change 16.39 (start to 10%), 1.34 (10–50%), 0.41 (50–90%), 0.408 (90–100%).
- *Production:* round trips 14, 17, 20, 18 (69); flatness 1.36 pooled (2.15, 1.36, 1.89, 0.38); above the window 3,274, 2,324, 2,958, 974 of 125,000; b = 1 Kish ESS 15,473 of 500,000 (5,187; 3,806; 4,765; 1,718).
- *Evidence:* ln Ẑ_c(1) −841.05, −841.45, −841.14, −842.24, SD 0.54 (range 1.18); pooled −841.38 ± 0.27. ln Ẑ(b) at 0.1, 0.2, 0.25, 0.3, 0.4, 0.6, 0.8: −125.9, −210.2, −251.3, −291.9, −372.3, −530.7, −686.7.
- *Basins:* a1 67.1% (−841.8), chains 67.3, 68.0, 68.1, 61.6%, χ² p 0.910; a1 + x6·x8 9.6% (p 0.552); B 0.8%, a2 0.2%, A under 0.05% (6,527 draws visiting); x1·x2 terms 3, 4, 5: 2.9, 96.8, 0.2%. σ's interval [0.9779, 1.1139]; D10's verdict converged in every group.
- *Time:* learning 1,254 s, production 793 s, all 2,052 s, at 44,000 sweeps a second a chain.
- *Against the first flat start* (t₀ = 10·K, `699bd8b`, the main session's): no round trip in learning or production, flatness 0.99 and 0.98, the weights still moving 1.5 nats RMS between 10% and 50%, b = 1 ESS 8, ln Ẑ(1) −843.64 ± 0.46, a1 74.7% failing the χ².

*What it shows* (a reading, not ruled): at t₀ = 100·K a flat start learns Friedman's window in the warm start's budget. Every chain crosses in learning and in production, and its evidence and posterior agree with the three warm-started runs (−841.38 against −841.14 to −841.55; a1 67.1%). The price is a gain still at t₀/t = 0.128 when learning stops: the weights moved 0.41 nats RMS over its last tenth, and production's flatness, 1.36 (one chain 2.15), is two to four times the warm starts' 0.34–0.58. One run, on one seed; it shows a flat start can stand in where the warm start is unavailable at this budget, not which start is the better.

**The suites** (Release, built into `artifacts\chip-fh3\`):
- *Samplers:* 277 passed, 0 failed, in 53 s. This chip's: the tied solve's three tests and the warm start's four, and the double-well and Potts start pins in two existing tests; the rest of the growth since Phase 2's 262 is other sessions' (SMC).
- *Regression:* 1,403 passed, 0 failed, in 7 min 4 s, at `94723f0` with another session's untracked cost-qualification files in the tree. This chip's: three new fixture-A tests (the gap pilots' diagnosis, the 48 pilots, the pilot too thin to tie), the gap test reworked (its two cases), three start pins in existing tests; and in the qualification harness 20 cases that run only under `BASS_QUALIFICATION` (the diagnosis, seed 8017's pin, the sweep's 16, the third seed, the flat start).
- *The BASS oracle classes* (`Maths.Oracle.Tests.Bass*`): 90 passed, 0 failed, in 6 min 39 s.
- `BassBitIdentityTests`' hashes hold, checked at step 2 and again on the shared tree before the fits.
- *The qualification runs* (`BASS_QUALIFICATION` = 1, Release, one at a time): the diagnosis 1 min 55 s; seed 8017's pin 1 min 3 s; the sweep 35 s – 4 min 17 s a seed; seed 8217 39 min 18 s; seed 8417 36 min 26 s and its warm start's read 101 s; the flat start 34 min 12 s.

**The review passes and their fixes.**
- *Step 1.*
  - The infeasibility helper mirrors the fixed point's targets (an empty node's share spread over the rest); it agrees with the engine's verdict on the five pilots its tests read, on fixture A's 7 refused (asserted in the 48-pilot test) and on the sweep's two.
  - The node report's "1 draws, 1 atoms" fixed in step 2's commit.
- *Step 2.*
  - The first remedy built in exploration, emptying the thin nodes, failed the double well (19–25 of 32 still refused): replaced by the tie before anything was committed.
  - Level edge ties put a smooth ln g's tied top nodes 47 and 94 nats off (13.8 RMS): the edges now continue at the extensions' slopes.
  - `ThinDraws`' docstring claimed every refused pilot read was accepted, which the small Gaussian pilots belie: it now names the populations measured and the small pilot's refusal. `WarmStart`'s said its refusal names the tied nodes, which one of its two refusals does not: corrected.
  - The tied solve's comment omitted the offsets; one signature reflowed; a splice in the test file had clipped the next test's header, restored before any build ran.
  - The sweep's first run lost its host with nothing written: one case a seed (`ca091f8`).
  - *Coverage added:* an empty segment beside a drawn node (a synthetic gap pilot, and 7802), a thin node at the window's edge (synthetic, 8217, 8313 and fixture A's 7), ĝ spanning more than 745 nats (a Gaussian ln g over 1,099 nats), the refusal kept (two pilots, one through `Bass.Fit`), the tied solve's own algebra (ln g recovered to 2.6e−9 where it is linear across tied nodes; edges; empty nodes; refusals).
- *Step 3.* The harness's pilot line said the fixed point accepted ĝ, which a tied start would belie: it now names the warm start. Its refusal line likewise. The third seed's warm start is read before its run and the run's start checked against it.
- *Final, against the brief's list:* (correctness) accepted starts bit for bit, by six pins and the sweep's 14; a tied node's weight is a start only, SAMC's learning moving every node (7802's fit, the replications); the refusal kept where too few draws remain; each run's criteria committed before it ran (`46332d7`, `ca091f8`, `7fc684e`; 8217's in `56c66f3`). (Coverage) the four items above, and the replication that fails a criterion, reported as it falls: none failed (8417's unlabelled a1 + x6·x8 disagreeing across chains is reported beside its criteria).

**For the main session.**
- *Friedman's a1 posterior now rests on three seeds*, each meeting what would count as committed. `docs/Maths/BayesianRegression.md` ("BASS: the flat-histogram sampler") and the design's "Phase 2, as built" still say one seed and "not yet robust"; `docs/Maths/Samplers.md`'s fixed-point text says a refused ĝ is refused. Neither is touched here.
- *The rule's threshold is absolute*, 10 pilot draws. A pilot too small for two nodes to hold ten draws is refused (a few thousand draws are needed); a threshold relative to the pilot's size is an option, not measured.
- *Accepted starts can sit far from SAMC's limit at the window's edges.* On fixture A's pilots over their own shapes the accepted starts' lowest nodes sit up to 17 nats off θ* (pooled RMS 0.84, against 0.29 for the tied starts). Tying thin nodes always, not only on refusal, would likely be more accurate, but would change accepted starts, which this chip's invariant forbade; it is a later option, with the window's lower margin (which reaches into the hot rung's sparse tail on every Friedman seed read) another.
- *Phase 1's 0.25 is not met* by the refused-then-tied starts (0.29 and 0.32 pooled) nor by the accepted starts on the same shapes (0.84 and 0.71); the brief's check is reported as it fell.
- *The adapter duplicates `Wham.NodeDraws`*: `BassFlatHistogramRun.NodeDraws` computes the same sum (identical arithmetic, so the readout's `NodeDraws` is the rule's count), left in place as the file was outside this chip's remit beyond its message. `BassDiagnostics`' doc of `NodeDraws` says a node holding none is interpolated; a thin node is now tied where the fixed point refuses ĝ.
- *The harness's readout cannot say whether a start was tied*: the BASS readout carries no field for it (`BassDiagnostics.cs` was outside this chip). The step-3 harness reads it before the run from the rebuilt pilot; a `TiedNodes` on `BassFlatHistogramPilotReadout` would carry it.
- *`FlatHistogramLimit.Infeasibility`* (test code) decides exactly whether SAMC's limit exists on a discrete measure; as an engine diagnostic it could say in the refusal's message which segment fails, rather than how far Newton strayed.
- *The sweep's first run lost its test host* about 45 minutes in, with no crash logged (other sessions' test hosts were running); the rerun completed, one case a seed. The fits ran alone, at 39,000–42,000 sweeps a second a chain against Phase 2's 52,000–63,000, the machine shared.
- *Other sessions edited BASS's core* (an extension-seams chip) during this chip; `BassBitIdentityTests`' hashes were checked on the shared tree before the fits, and held.
- `docs/Maths` and the checklists are untouched. `artifacts\chip-fh3\` is deleted.

Recommend no default.
