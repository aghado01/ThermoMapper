# Regression coordination: handoff (2026-11-09)

Written by the coordinating session (`local_b6732445`) for a new parent session to take over. It gathers what is in flight, what waits on the user, and where each record lives. The briefs are the records of truth; this file points at them and lists what crosses them.

## How the coordination has run

- **The parent coordinates chips (sessions) and verifies each slice before the next starts.** Verifying means reading the commits' diffs and rerunning the tests from a snapshot of HEAD.
  - *Snapshot:* `git archive HEAD | tar -x -C artifacts/verify-main/snap --exclude=r`. The oracle tests find the real `r/` by walking up from the snapshot.
  - *Selection:* `scripts/test-select.ps1 -SourceRoot <snap> -Scope Commit -Changed <files since the last verified commit> -ArtifactsRoot artifacts/verify-main/bin/`. The totals are in the `.trx`'s `<Counters>`.
  - *Oracles:* each oracle class is run with `dotnet test tests/Maths.Oracle.Tests/... --filter <class>`.
  - *Hygiene:* run at below-normal priority with `MSBUILDDISABLENODEREUSE=1`, and delete `artifacts/verify-main` after.
- **Chips report slice by slice and wait for a verification note.** Decisions beyond a brief come back to the parent, and anything that changes a default or removes an option goes to the user. A removal is always its own question.
- **Repository rules.** Work on `main`, commit by pathspec (`git commit -m … -- <paths>`), leave nothing staged, push nothing.
  - Commit messages carry no attribution lines (CLAUDE.md).
  - Leave out of every commit `issues/regression/bass/basin-hold-validation.csv` and `path-profile-validation.csv`, and the user's own notes (`issues/regression/discussion/opus-open-items-11082026.md`, `issues/geometry/gemini-abliteralus-thermomapper.md`).
- **The chips' hand-over notes** are in `C:\Users\azrie\AppData\Roaming\Claude\side-session-notes\local_b6732445-a9c0-404a-81ad-1ab68d7b1dd2\`, one per chip and slice (for example `e34-multi-output-slice-1.md` and `gbass-noise-slices-6-7.md`). That folder belongs to this session, so a new parent should read from it and tell the chips where to write next.
- **The extension index:** [extensions-checklist.md](bass/extensions-checklist.md), [extensions-brief.md](bass/extensions-brief.md), and the build order in [lineage-review.md](bass/lineage-review.md#build-order).

## Repository state at handoff

- **HEAD before this file's commit:** `40ee799` on `main`. No worktrees, no other branches, no stashes, so there was nothing to merge.
- **Uncommitted work, all from chips stopped at the usage limit.** Leave it for its owner to finish:

| Files | Owner | State | Next |
|---|---|---|---|
| `BassConfig.cs`, `BassLadder.cs`, `BassModel.cs`, `BassObservationModel.cs`, `BassPriorSampler.cs`, `BassRandomIntercepts.cs`, `BassResolved.cs`, `BassSmc.cs`, `BassState.cs`, `BassSweep.cs`, `posterior/BassPrediction.cs`, `tests/regression/RandomEffectsTests.cs` | GBASS chip: ridge blocks, slice 1 (s drawn) | Source written and compiling, tests written, none run, nothing committed | Run `RandomEffectsTests`, the held-s pins and `RandomInterceptTests`; run the Commit scope; commit; then the qualification gates and the docs (its hand-over note lists them) |
| `tests/regression/MultiOutputLawTests.cs` | E34 chip | Two Geweke seed rows, uncommitted until the ln τ defect below is fixed | See "Open defects" |

## The sessions

Links open each session in the app.

| Session | Role | Brief | State | Next |
|---|---|---|---|---|
| [E34: several outputs on one basis](claude://claude.ai/epitaxy/local_26b4e7ee-1be4-482b-a19d-995477d32291) | E34 build | [multi-output-brief.md](bass/multi-output-brief.md) | Slices 0 and 1 verified (to `2b5b60e`). Stopped at the limit during the Geweke recheck | Fix the ln τ defect, then slice 2: samplers, evidence, Huang and Wand as an option |
| [GBASS slices](claude://claude.ai/epitaxy/local_add4be38-7ddb-4981-8b14-6f8ac64f7dab) | E1 built; now E21's remaining work (ridge blocks) | [ridge-blocks-brief.md](bass/ridge-blocks-brief.md), [gbass-noise-brief.md](bass/gbass-noise-brief.md) | Ridge-blocks slice 0 verified (`7cf1761`). Slice 1 uncommitted (above) | Finish slice 1, then slices 2 to 5 |
| [KHAOS slices](claude://claude.ai/epitaxy/local_d665e4ce-cd64-427f-9c2a-8c94d96497dd) | KHAOS (E19) built; now E4's kernels | [proposal-kernels-brief.md](bass/proposal-kernels-brief.md), [khaos-brief.md](bass/khaos-brief.md) | Slice 1 (history birth) verified. Slice 2 (residual birth) landed. Its code to `ff2ec0d` passed the coordinating session's 2,408-test run. `95e5af7`'s probability-cache change and the slice report **are not yet verified**. The session ended on a limit error | Verify `95e5af7`, then put slice 2's decisions to the user (below) |
| [Conformal slices](claude://claude.ai/epitaxy/local_69bae1e5-39ad-45c0-8df2-e8e705b21e42) | E16 conformal layer; E35 built | [conformal-layer-brief.md](conformal-layer-brief.md) | E35 verified. The user's three calls recorded (`40ee799`) and sent. Slice 7, seams that name the forms they serve, landed after this file was written (`e3071cf`, `941c1b1`, `f4e6b83`; the whole assembly 2,736 of 2,736 by the chip) and **is not yet verified**. Weighted BASS candidates turned out to have no reachable fit to serve | Verify slice 7, then slice 8; 9, 10 and 12 wait for their Lean |
| [Conformal validity in Lean](claude://claude.ai/epitaxy/local_21649695-833f-41ac-8d7c-326b4f691dda) | Lean for the conformal layer, plus SG-4(b), SC-5 and CW-1 | [conformal-chip-brief.md](../../lean/Protolemmata/Conformal/conformal-chip-brief.md), [lean-leads-chip-brief.md](../../lean/Protolemmata/BASS/lean-leads-chip-brief.md) | Idle, all verified | Commission it for the Lean of conformal slices 9, 10 and 12 (Lean first, the user's call) |
| [E8: active subspaces](claude://claude.ai/epitaxy/local_66d45fa0-e3f5-4429-b588-1786ae46c0ef) | E8 build | [active-subspace-brief.md](bass/active-subspace-brief.md) | Complete, slices 0 to 5 verified | None. Its follow-ups are listed below |
| [Grassmann and Stiefel geometry](claude://claude.ai/epitaxy/local_58551eda-672f-47a3-8439-28c0db7e8f81) | The user's discussion | [grassmann-stiefel-subspaces-brief.md](discussion/grassmann-stiefel-subspaces-brief.md) | Its items are settled: E8 slice 5, E37, E36's direction prior | Whatever the user takes up there |
| [IRLS solver scale fix](claude://claude.ai/epitaxy/local_f817795b-ddff-4126-9061-e2e154da0004) | Geometry fixes, started from the discussion | [irls-scale-brief.md](../geometry/irls-scale-brief.md) | Idle. It landed the IRLS scaling, MoMPCA's calibration, stable principal angles (`0b4ed4e`, `db0d13f`), asinh hyperbolic models, and test-select's two fixes (`e36330f`). **Its commits were not rerun by the coordinating session** | Its open findings below |
| Older, idle | [Parallel rungs](claude://claude.ai/epitaxy/local_3d04d451-1151-47ba-ba9a-cd0db122380b), [adaptive SMC](claude://claude.ai/epitaxy/local_47a5acd7-ea2c-4afa-acaa-9fd1c4cb06a1), [targeted tests](claude://claude.ai/epitaxy/local_f6922d89-48ce-4c93-b35e-9cb0de18c1d0), [lean leads](claude://claude.ai/epitaxy/local_198ac1bc-ddd9-4159-a180-ad1ed5a8de1b), [path profile](claude://claude.ai/epitaxy/local_8e1db685-72b7-4ee4-ba95-1ef7407166d6), [efficiency pass](claude://claude.ai/epitaxy/local_610246f0-d98b-41d3-a671-4cd2d636a6ec), [cost profile](claude://claude.ai/epitaxy/local_3b899cb2-2bde-499b-9d1a-83341dceecd9), [spectral gaps](claude://claude.ai/epitaxy/local_3153bd3a-eea1-4ee8-80b5-0a58899f9322), [WHAM refusal](claude://claude.ai/epitaxy/local_9c9a9e89-927e-4217-b89e-1e3caf3dc14e) | Their briefs and reports in `issues/regression/bass/` and `issues/tooling/` | Done | Archive when the user wishes |

## Waiting on the user

1. **SMC's trap remedy** (proposal kernels slice 2's report).
   - *What was measured:* the residual birth reaches the posterior in 16 of 16 replicates, but only with a step cap of 1.25. Under the engine's schedule it reaches it in 0 of 16.
   - *The question:* whether the pair of cap and residual birth becomes the remedy for KHAOS's SMC trap, beside the flags.
   - *What the pair does not fix:* the evidence still spreads 21 nats. A reference nearer the posterior is what the evidence would need.
2. **The residual birth as a candidate for KHAOS under the raw ridge.** About 2× ESS per second at the same cost, and a 12× shorter climb from an empty structure. It is one fixture, so no default is proposed. On hinge fits it costs 19 to 88× the CPU with no gain, and it stays opt-in.
3. **E34's default prior on Σ.** The default is the IW that reduces to scalar BASS's at K = 1; it puts a correlation's prior mass near ±1 at small n/K. Slice 3 measures coverage against Huang and Wand at ν = 2 for n/K from 3 to 30, and the default question comes back with those numbers.
4. **GBASS's open decisions** ([gbass-noise-brief.md](bass/gbass-noise-brief.md)):
   - the flat histogram under a noise: the walker carries a and σ², a Metropolis step on each latent, or keep it refused;
   - quantile BASS's intervals undercover (63.8% at nominal 90%), which belongs to E15;
   - interweaving's z′ = σ²·v default for scale mixtures, which needs a measurement;
   - a prior on ln γ for the Normal-Wald.
5. **Parallel rungs** ([rung-parallelism-report.md](bass/rung-parallelism-report.md)):
   - `RungWorkers`' default waits for a wall-clock reading on free cores;
   - the rung chip's GC finding (server GC took a 4-chain fit from 2.46 to 0.93 s) was proposed as an allocation or GC pass.
6. **Lean promotion.** The proof-complete enthymemata (`ConformalValidity`, `GuardedEstimator`, the BASS and BGP modules, `Bifiltration`, `PottsPK`) wait for the user's semantic review before moving to `Lemmas/`.
7. **E8's follow-ups:**
   - rotate-and-refit, and BPPR (E36);
   - pairing a basis fit with another fit, a co-active projection mean, and output locations beyond the training ones;
   - caching the matrix Poincaré check so it can default on.
8. **LinAlg's `DenseEigen`** uses an absolute 10⁻¹² threshold, which is inaccurate on small-scale matrices (E8's slice 4 found it; E8 rescales around it). A relative criterion could move other consumers' bits, so it is the user's call.
9. **The next wave of extensions.** The coordinating session recommended E25 (dependent inputs, which follows E8) and E15 (calibrated intervals, including quantile BASS) once E34 and the conformal slices land. E3 (ABMRS) is queued after the current queue, by the user's word.
10. **Test tooling's older follow-ups** ([targeted-test-runs.md](../tooling/targeted-test-runs.md)):
    - gated tests that report a skip, not a pass;
    - pinning option-on config fingerprints;
    - splitting `BassFlatHistogramTests`;
    - the GBP and NVM distributions tests, which took 24 minutes in one filter.

## Open defects and investigations

- **E34: a ln τ bias under the K × K form, on one Geweke seed** (20261123, z −7.25 at the scalar tests' length; three other seeds pass).
  - *What is ruled out:* the target matches the matrix-t reference at all 600,000 visited states to 2.8·10⁻¹³.
  - *Where the chain falls short:* it under-reaches ln τ's lower tail (0.1%, 1% and 5% quantiles at −10.4, −8.5 and −6.8, against the prior's −15.9, −11.3 and −8.0).
  - *Next* (its hand-over note): evaluate the target over ln τ from −25 to 0 at fixed structures for −∞, NaN or a refusal; test τ's step against quadrature; fix; rerun the four seeds; commit the seed rows.
- **The IRLS chip's open findings** ([irls-scale-brief.md](../geometry/irls-scale-brief.md)):
  - Weiszfeld's slow rate near a data point that barely fails the optimality condition;
  - `EuclideanMedian` stalling when warm-started on a data point;
  - GMM's robust initialization untested;
  - `DiagramMetricsTransportParityTests` timing out in R.
- **KHAOS's SMC evidence** stays unreliable on the trap fixture under every remedy measured. The flags report it.

## LinAlg work queued after the handoff (2026-11-09)

The [eigensolver brief](../linalg/eigensolver-performance-brief.md) goes to the next available LinAlg owner, independently of E21, E34 and conformal slice 8. First fix the capped LOBPCG residual and EigenFast's no-rotation sweeps. Then qualify a production Ritz solver on actual block sizes and measure reuse of `AX` separately. Reserve the shared machine's quiet window for timing comparisons. The E8 tolerance question above remains a separate decision.

## Recorded follow-ups, no decision waiting

Each brief's report lists its own, with the decision each would change:
- [active-subspace-brief.md](bass/active-subspace-brief.md)
- [multi-output-brief.md](bass/multi-output-brief.md)
- [proposal-kernels-brief.md](bass/proposal-kernels-brief.md): Nott and Kohn's best linear predictor, adaptation past burn-in, history weights in the other moves
- [ridge-blocks-brief.md](bass/ridge-blocks-brief.md)
- [conformal-layer-brief.md](conformal-layer-brief.md)
- [khaos-brief.md](bass/khaos-brief.md)
- [gbass-noise-brief.md](bass/gbass-noise-brief.md): five speed and memory items; the half-Cauchy σ² gate; the skewed horseshoe's hot-rung limit
- [smc-adaptive-mutation-report.md](bass/smc-adaptive-mutation-report.md)
- the Lean contraction line in [contraction-walk.md](../../lean/Protolemmata/BASS/contraction-walk.md): the local restatement, σ²'s upper tail, the fast rate at K ≥ 2

## Decisions taken on 2026-11-08 and 2026-11-09

Each is recorded in its brief, with the user's words where given:
- E34 ordered after E8, then commissioned.
- The move along β·v is the skewed horseshoe's default.
- SMC keeps its fixed mutation count.
- The conformal margin is opt-in.
- The KHAOS slice order was 5, 4, 3.
- The SMC trap is flagged (the coordinating session's call, which the user can reverse).
- E8 was commissioned, and slice 5 added from the Grassmann discussion: per-draw k opt-in, the output side, basis adequacy, stable angles.
- E36's direction prior: independent by default, with D06's gate on ridge columns and the determinant as an option with c_M stated.
- E37 recorded.
- E4's kernels went to the KHAOS chip, E21's remaining work to the GBASS chip, E16's completion to the conformal chip, and E3 was queued.
- The three conformal calls (`40ee799`): slice 7's refusal of the hierarchical seam, Lean first, and slices 7 to 12 now.

## What crossed between sessions without passing the user

- **The Grassmann discussion** commissioned E8's slice 5 and wrote E37, from the user's settlement in that session (`c5859ad`, `07c60ec`). It started the IRLS chip, which fixed test-select's two defects (`e36330f`), among them the earlier misreading of a record parameter as a type.
- **E8's chip** answered the user's question about the Poincaré ball and hyperboloid against the Poincaré inequality, in its own session.
- **KHAOS and E34** agreed an order for `BassConfig.cs` and `BassSweep.cs`. GBASS and E34 agreed the whitening's move onto the working response.
- **E34 corrected the coordinating session's brief:** the intercept's count is n, not n − 1 (`8196778`).
- **E34 found that decision 2's prior** puts a correlation's mass at ±1 (now in "Waiting on the user").
- **Defects found in the sources and oracles:**
  - the truncated normal's second moment in 2307.11241 and `coactivity`;
  - three slips in 2311.18146 §3.1.1;
  - `coactivity`'s converters;
  - `gcbeck/bppr`'s determinant prior missing its normalizer in the birth and death ratios.

  Each was corrected on the engine's side, never ported, and is recorded in the briefs and in `docs/Maths/BayesianRegression.md`'s corrections table.

## Resources

- **Oracle pins** in `r/renv.lock`: GBASS 2.0.1, khaos 2.1.1, BASS 1.3.1 (source), and `coactivity` at `knrumsey/coactivity@35e4573` with `lhs` and `zipfR`.
- **Read-only clones** in `D:\aghado01\packages\bayesian-related\`: `coactivity`, `concordance` (`b34e7c6`), `bppr-gcbeck` (`1b0d47a`; GPL-3.0, ideas only).
- **Papers:**
  - `D:\aghado01\codex-scientiae\supellex-private\gauntlet\bayesian\`: 2210.09181 (BPPR), 2307.11241, 2311.18146 and the conformal-martingale set;
  - `…\PDFs\Bayesian\`: `jrssb98.pdf` (Brown, Vannucci and Fearn), `HuangWand2013.pdf`, `NottKohn2005_Biometrika92.pdf`, `Nott2005_StatsComp15.pdf`, `BayesPPR.pdf`.
