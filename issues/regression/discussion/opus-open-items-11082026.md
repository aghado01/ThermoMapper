Here's everything open, grouped by what it's waiting on.

**In progress**
- **E8 slice 0:** in the working tree, waiting on its full regression run. I'll verify it next. The chip already found four defects in `coactivity` that it inherited from `concordance`; we correct them on our side.
- **The Grassmann/Stiefel discussion:** waiting on you in [its pane](#local_58551eda-672f-47a3-8439-28c0db7e8f81).

**Waiting on your decision**
1. **Conformal slice 5:** candidate seams for BARS, BAPS, the hierarchical P-spline, BGP and changepoints. E35's C# build follows it; its Lean is done, with six gates recorded.
2. **GBASS's flat histogram under a noise:** have the walker carry the coefficients and σ², move each latent by Metropolis, or leave it refused, since ladders and SMC already handle noisy fits.
3. **Quantile BASS's intervals undercover** (63.8% for nominal 90%). This probably folds into E15, calibrated intervals for loss-based posteriors.
4. **Parallel rungs:**
   - whether the rung-parallelism default changes, which needs a timing on free cores;
   - the rung chip's finding that fits are GC-bound (server GC cut a 4-chain fit from 2.46 s to 0.93 s), which it proposed as its own session.
5. **Lean promotion:** your semantic review of the proof-complete modules (conformal, guarded estimators, four BASS, five BGP, Bifiltration, PottsPK) before they move to `Lemmas/`.
6. **Archiving the finished chip sessions:** KHAOS, GBASS, adaptive SMC, parallel rungs, the conformal Lean chip and lean-leads.

**Build queue, in order**
- **E8 slices 1–4,** then **E34** (its build brief is written once E8 lands), then **E22** on top of E34.
- **Conformal slice 5,** then **E35's C# build.**
- **E25** (dependent inputs) follows E8.
- **Rotate-and-refit and BPPR** wait for E8 slice 4 and the Grassmann discussion.

**Recorded follow-ups, nothing pending** (pick up when useful)
- **GBASS:**
  - the interweaving default for t and other scale mixtures;
  - a ln γ prior for the Normal-Wald;
  - five speed and memory items;
  - the `hbass` half-Cauchy joint-law gate;
  - the skewed horseshoe's hot-rung limit.
- **SMC:**
  - the three trap remedies not built (a step cap, a residual-informed birth, a start nearer the posterior);
  - an adaptive count that sees continuous coordinates;
  - the mutation-steps default, 5 against 2.
- **KHAOS:** the modified g-prior on hinge terms and under a weighted prior Gram. For random intercepts: a sampled spread, nested groupings and random slopes, plus running them under noises and conformal, which are refused today.
- **Conformal:** a per-unit bound that leaves out degenerate orderings.
- **Lean:** the rest of the contraction line, and the old open proofs (BARS ×6, Stability, CumulativeField/Local).
- **Geometry:** the Stiefel distance and log-map rewrite (when something needs frames), and sine-based Grassmann angles.
- **Test tooling:**
  - gated tests reporting a skip rather than a pass (the old card);
  - pinning option-on config fingerprints;
  - splitting `BassFlatHistogramTests`;
  - the 24-minute GBP/NVM distributions tests.
- **Unscheduled extensions:** E15 and E18, among others.
- **Minor:** the path-profile brief doesn't say which σ² prior it ran.