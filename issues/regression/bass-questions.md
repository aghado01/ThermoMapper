# N14 verified (records in 2d1ba7d):

One rule, ConvergenceVerdict.Of, now drives both the stop between rounds and the verdict at the end; the old MayStop is gone.
I reran the suites and all pass: Samplers 204/204, Regression 1,085/1,085, BASS oracle 63/63. Those counts include stage 7’s in-progress files, which were in the tree at the time. The eight draw hashes are unchanged.
Where the chip departed from my brief: frozen readings fail, as D10 says. My brief’s API section had wrongly called them “undetermined”.
Two consequences:
Every BASS fingerprint moved, because the config’s JSON now carries StoppingScope. Checkpoints from before the change won’t resume, which is fine before release.
A fit with an R̂ target can now stop later, because the log posterior is part of what it waits on.
On S10’s fits, every group reads not converged. The gap in log posterior between chains is 42.8 nats on Friedman, 21.8 on mixed, 30.9 on functional, and 27.0 on tempered Friedman.
Docs: docs/Maths/Samplers.md and docs/Maths/BayesianRegression.md now describe the verdict as built. N14 is ticked, and D10 is marked built.


# Stage 7, calibration

is verified and documented (2232dc3). With it, everything in the BASS port is built except stage 8, SMC.

What I checked

Suites: Regression 1,101/1,101, Distributions 148/148 (for the new HalfCauchy), and the BASS and calibration oracle tests 63/63.
R’s stage-1 checks: 97 verdicts pass. My first run failed on the first oracle it launched, but that was my mistake, not the code: the script looks up its oracles relative to r/, and I had run it from the repo root. From r/ it passes.
Code against D09:
The ln σ² move’s Jacobian stays out of the tempering.
Weighted emulator draws are resampled by weight, then each inner run counts equally.
The half-Cauchy’s formulas are right.
The cut holds: each inner run keeps its emulator draw fixed.
Recovery: on a synthetic simulator, the 95% intervals cover the true θ in 5/5 seeds for each parameter with a basis emulator. With the direct emulator, one parameter is covered in 4/5 seeds and the other in 5/5; the direct emulator is the less accurate of the two. Coverage is 5/5 again with a declared discrepancy and with two variance groups.
Inner-run length: at 25 kept draws per chain the pooled mean is off by about 0.2 posterior sd and its spread 15% low. From 100 draws on, the error is Monte Carlo noise.

Docs and records

docs/Maths/BayesianRegression.md has a new “BASS: calibration” section, plus API, configuration, status and source rows.
docs/Maths/ProbabilityAndInformation.md lists HalfCauchy.
C1–C3, C-X1 and C-X2 are ticked.
D09 now records three points its text left loose, as built:
the multivariate correlation runs over outputs, as in R;
the multivariate truncation term defaults to the residuals’ second moment, matching the functional one;
the diagonal option uses R’s centred n − 1 variance.

## Open, for later

How often inner runs converge under the default gate is unmeasured. The recovery tests used a looser 1.05 / 50 gate, and 3–6 of 6 inner runs converged per seed.
A new numerics item, N15: DiagonalPlusLowRank runs a full thin SVD on every σ² proposal. A Cholesky factor of the small capacitance matrix is the cheaper route if real sizes make that cost bite.
Above the categorical enumeration cap, single-site moves are Metropolis steps; a per-site Gibbs draw would avoid rejections.


## Stage 8 planning 

I've got the generic SMC core, D05's sequence, and D09's cut-SMC shape in place, with `BassFit` already carrying draw weights. Before finalizing the brief, I need to verify the weighted type-7 quantile reference.

Both references check out, and Mathews et al. is already in the paper library. Drafting the two stage 8 briefs now: C4, fitting SMC, and C5, cut SMC.

I've planned stage 8 as two chips that can run in parallel without touching the same files. Both briefs are committed (`542c442`), and nothing has started yet.

The briefs are [stage8-fitting-smc-brief.md](issues/regression/bass/stage8-fitting-smc-brief.md) (C4) and [stage8-cut-smc-brief.md](issues/regression/bass/stage8-cut-smc-brief.md) (C5).

**What's already fixed.** D05 sets C4's sequence: start from the prior, temper on the likelihood path, use D04's sweep as the mutation, place stages by conditional ESS, and resample multinomially. D09 sets C5's shape after Mathews et al.: carry particles through each emulator draw's conditional, order the draws, and bridge between them through their induced densities. The generic SMC core (N8) is built, and `BassFit` already carries draw weights into its views.

**My recommendations beyond those, for your call:**
1. **Two chips in parallel.** C4 owns BASS's fitting files and `DrawSample`; C5 owns `calibration/`.
2. **C4's settings (`SmcSpec` on `BassSampling`), provisional and not tuned:**
   - 1,000 particles, a conditional-ESS target of 0.99, and 5 sweeps per stage;
   - `Chains` counts independent SMC replicates, 4 by default;
   - `Tempering` and an R̂ target are refused alongside SMC, because SMC tempers itself and ends when b reaches 1.
3. **Weights flow through.** Each replicate's final weighted population becomes the fit's draws. `PointwiseLogDensity` puts the weights into its importance ratios, checked against an equal-weight resample.
4. **`DrawSample` becomes one continuous weighted type-7 quantile.** It matches type 7 exactly at equal weights, varies smoothly with the weights, and returns the min and max at p = 0 and p = 1. It's built from Akinshin's *Weighted quantile estimators* (arXiv 2304.07265). BARS stays bit-identical.
5. **Checks across replicates:**
   - R̂ across replicates as a between-run reading, kept apart from each replicate's ESS and families;
   - the evidence pooled by `LadderEvidence.Pool`'s rule;
   - SMC fits resume from a checkpoint.
6. **C4's qualification answers D10's stage 8 question.** On S10's fixtures at 250, 1,000 and 4,000 particles, it asks whether the population finds and keeps Friedman's better representation. It reads only and changes no default.
7. **C5's design:**
   - the first leg tempers from the prior;
   - the emulator draws are toured nearest-neighbour then 2-opt, on the distance between their predictions at 64 space-filling θ points, with a random order as the comparison;
   - between draws it uses geometric bridges, with t placed by conditional ESS;
   - each block of stage 7's sampler is aimed at the bridge (the discrepancy's draw stays Gaussian);
   - each draw's population is pooled with weight 1/S;
   - it's checked against stage 7's nested reference and timed against it at matched error.

Both papers check out: Mathews et al. is arXiv 2406.07555, with a local copy in the gauntlet library, and Akinshin is arXiv 2304.07265.

Sources:
- [Akinshin, *Weighted quantile estimators*, arXiv 2304.07265](https://arxiv.org/abs/2304.07265v1)
- [Mathews et al., *Sequential Monte Carlo for Cut-Bayesian Posterior Computation*, arXiv 2406.07555v2](https://arxiv.org/abs/2406.07555v2)

Going through them in order. Re-reading the briefs with fresh eyes, I'd change four of my own recommendations; I've marked those as revisions.

**1. Two chips in parallel.** I still recommend it: C4 and C5 share no dependency except the weighted quantile, and C5's checks already avoid that by comparing means and variances.
- *A flaw in both briefs:* each lets its chip add accessors in `src/maths/samplers/smc/`, so two writers could touch the same folder.
- *Revision:* neither chip edits `smc/`. Its public API (`SmcRun.Start`, `SmcSampler.RunTempered`, `NextInverseTemperature`) should be enough, and a chip that needs more names it in its report.

**2. C4's settings: 1,000 particles, conditional ESS 0.99, 5 sweeps per stage.**
- *0.99 makes for many gentle stages.* It's the safer setting for the evidence and for the question that matters here: whether reweighting through Friedman's steep band, b ≈ 0.16–0.32, sorts the representations rather than collapsing onto one. 0.9 would be roughly an order of magnitude cheaper, with more resampling noise.
- *5 sweeps per stage* matches D05's evidence check. The genealogy error was calibrated at 1 and 5 moves, and the Kish error under-covered at 1.
- *Cost at Release speeds:* roughly a minute per fit on Friedman at n = 500, and a few minutes at n = 2,000.
- These are provisional and read by the qualification, never tuned from it. I'd keep them.

**3. `Chains` as the replicate count.** *Revision.* This breaks the names-follow-methods rule: "chains" isn't true of independent SMC runs. Better is `SmcSpec.Replicates`, 4 by default, with `Chains` refused when `Smc` is set, as `Tempering` is.

**4. Weights flow through.** Keeping the weighted population is the retain-don't-discard default, and the views already take `DrawWeights`. If the last stage resampled, the weights come out equal anyway.
- *What the brief misses:* particles that share ancestry aren't independent draws. So `PointwiseLogDensity`'s Monte Carlo errors would read optimistic, as would any ESS computed as though the draws were independent.
- *Revision:* the chip reports the families beside those errors and doesn't present them as calibrated. The replicates are the honest error.

**5. `DrawSample` as one continuous weighted type 7.** I'd keep this in stage 8, because unequal weights arrive here. Without it, SMC's intervals would come from the midpoint rule and MCMC's from type 7: two samplers of one posterior reporting quantiles by different rules.
- *The risk:* I haven't read Akinshin's construction. The brief already lets the chip use any form it shows equivalent, as long as it meets the properties (type 7 at equal weights, continuous, Q(0) = min, Q(1) = max).

**6. R̂ across replicates.** *Revision.* Rank-normalized R̂ splits each chain in half, and a population has no order, so split R̂ means nothing here. Unsplit R̂ would invite comparison with D10's 1.01, which assumes autocorrelated chains.
- *Cleaner:* for each monitored quantity, the spread of the replicates' posterior means in posterior-sd units, beside the evidence's between-replicate error. That says directly whether the runs agree to within a fraction of the posterior's width.

**7. The evidence:** each replicate's normalizer with its genealogy error, pooled by `LadderEvidence.Pool`'s larger-of-two rule, and reported only under D03's proper defaults. I'd keep it as is.

**8. Resume.** `SmcCheckpoint` exists, so wiring it into `Bass.Resume` keeps the fit's contract uniform. But SMC runs are short and the wiring isn't free. It's the most droppable item: keep it if cheap, and defer it if it means reworking the checkpoint type.

**9. C4's qualification** (S10's fixtures, 250, 1,000 and 4,000 particles, 4 replicates). This answers D10's stage 8 question. If populations reliably reach Friedman's best representation where MCMC chains stall 24–42 nats below, that opens a later question: should SMC become the default for structured fits? That would be a pre-registered pass of its own, not a reading from this one. I'd keep it as exploration.

**10. C5's ordering.** *Revision.* I proposed the L₂ distance between predictions at 64 fixed points in the unit θ box. That measures closeness where the posterior has no mass, which isn't where a bridge loses ESS.
- *Better:* the spread of the log-density difference between two emulator draws, ℓ(θ; e_s) − ℓ(θ; e_s′), over the first leg's particles. Conditional ESS reads exactly that quantity, so the tour minimizes what the bridges pay.
- *Cost:* S × N predictions once, then S² × N arithmetic: about 200,000 predictions for 200 draws and 1,000 particles.
- *Validity is unaffected.* Any order gives a valid estimator, since each μ_s is targeted exactly, so letting the ordering depend on the population costs nothing in validity. The space-filling design stays as the fallback.

**11. C5's bridges, mutation and pooling:** the geometric bridge, with t placed by conditional ESS. Each block of stage 7's sampler is aimed at the bridge: the discrepancy stays Gaussian, the categorical block enumerates the bridged conditional, and the θ walk's scale is fixed per stage from the population. Each draw's population is pooled with weight 1/S. I'd keep all of this; it's D09 carried out.

So the four revisions are no edits to `smc/` (1), `Replicates` in place of `Chains` (3), the spread of the replicates' means in place of R̂ (6), and C5's distance on the log-density scale (10). The note about ancestry-optimistic errors (4) is a smaller addition to the C4 brief. Everything else stands.