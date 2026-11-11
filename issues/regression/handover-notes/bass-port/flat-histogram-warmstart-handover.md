# Flat-histogram warm start, made robust; Friedman replicated (N17): hand-over, 2026-10-31

Full report: `issues/regression/bass/flat-histogram-warmstart-brief.md`, under Reports ("The warm start made robust, and Friedman replicated").

- Commits on `main`, none pushed: `694dfa9` (step 1, diagnosis checks), `4f0dd6f` (accepted-start pins taken on the unchanged engine), `46332d7` (remedy + checks + sweep harness, third seed 8417 declared), `ca091f8` (sweep one case a seed), `7fc684e` (step 3 tests, before their runs). The report is `3990c63`. On the way it was swept into the cost-qualification session's `fa39a1d` by a shared-index race and taken back out by that session in `50650e6`; my own first report commit had swept in that session's staged files, undone at once by a soft reset of my commit only.
- Diagnosis: SAMC's limit on WHAM's discrete ĝ does not exist where a node's hat sees ĝ only through atoms far from the node (flat occupancy outside the hull of the atoms' hat vectors). 8217: node 0 one atom at u = 0.869, outside at segment 1 by 0.0183; Newton diverges (|dir| to 2.6e17), "125 from flat" at node 3 is where the mass ran. Fixture A 7801/7802 the same beside empty segments. Not Newton's fault.
- Remedy: where the fixed point refuses ĝ, nodes under 10 pilot draws are tied (knot removal; interior on the line through held neighbours, edges continuing at b_lo below / 1 above); `SamcLearning.FixedPoint(…, tied)`, `Wham.WarmStart(pilot, shape, thinDraws, out tied)`, `Wham.NodeDraws`, `Wham.ThinDraws` = 10. Accepted starts bit for bit (six pins + the sweep's 14).
- Friedman sweep, 16 declared seeds: fixed point refused 2 (8217, 8313), warm start 0.
- Replications (criteria of 56c66f3, all met): 8217 ln Ẑ(1) −841.55 ± 0.07, a1 66.6%, SD 0.15; 8417 (start accepted as it stood) −841.14 ± 0.32, a1 68.4%, SD 0.63 (unlabelled a1 + x6·x8 disagrees across chains, reported). With 8017: three seeds within 0.41 nats.
- Flat start at t₀ = 100·K (reading): learns the window, 69 production round trips, ln Ẑ(1) −841.38 ± 0.27, a1 67.1%; frozen weights noisier (flatness 1.36).
- Not met: Phase 1's 0.25 (tied starts 0.29/0.32, accepted starts on the same shapes 0.84/0.71). Tiny pilots (< two nodes with ten draws) refused.
- Suites: Samplers 277/277, Regression 1,403/1,403, BASS oracle classes 90/90; BassBitIdentityTests hold.
- `artifacts\chip-fh3\` deleted. docs/Maths and checklists untouched. No default recommended.
