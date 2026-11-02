# Performance qualification, the statistical half (chip brief)

Date: 2026-10-31. Status: briefed for a chip; the decisions below are the main session's recommendations, taken when the user starts the chip. This brief is the hub's performance qualification ([its section](../bass-full-port-brief.md#validation-and-performance-qualification)) for what a run buys: statistical error per unit time per backend, read with each backend's own error report beside it, and the run lengths D10 left as a budget. What a run costs is the cost half ([its brief and report](performance-cost-brief.md#reports)). This brief's report goes under "Reports" below, newest first. `docs/Maths` and the checklists are the main session's.

## Why

The hub asks for statistical error per unit time with method-appropriate diagnostics, and warns that particle-weight ESS and autocorrelation ESS are not one currency. The defaults' lengths, 9,000 burn-in sweeps and 1,000 kept draws, are R's and provisional (`BassSampling`); D10 made run lengths a budget and left the budget open.

What stands:
- *The cost half* (2026-10-31): per-unit costs on every axis, and the cost model at the cap. It reads cost only.
- *F-R1* (2026-10-31): untempered BASS failed the predictive R̂ gate on every deterministic benchmark at S10's lengths, with about 5 effective draws in 8,000. There the self-report was honest. A chain held in one basin can report a confident ESS: S10's untempered Friedman chains each held a narrow band of M.
- *Friedman:* only the flat histogram meets its pre-registered criteria, now on three master seeds (ln Z(1) −841.4, a1 at 66–68%). SMC's 16,000-particle replicates disagreed by 6.5 nats, and tempering completes no round trip.

So the question, per backend and regime, is three-part:
- how much work buys a stated error on the functionals a user reads;
- whether the backend's own error report can be trusted;
- whether D10's verdict, and its early stop, say so.

## The decisions (recommended; starting the chip takes them)

1. **One currency: replicate effective draws.** A backend's effective draws for a functional ĝ at work W are V / MSE(W).
   - MSE(W) is taken over R independent replicate fits against the functional's reference, and V is g's posterior variance (p(1 − p) for an indicator).
   - Every backend is scored the same way.
   - Each backend's own report is read beside it as a self-report, checked against the replicates, and never used as the currency: the autocorrelation ESS for MCMC and tempering's cold rung; the particle ESS and the populations' spread for SMC; the chains' spread and the Kish ESS for the flat histogram.
2. **Error per unit work, converted to time once.**
   - Work is counted in each backend's own units: sweeps; rung-sweeps and swap rounds; particle mutation steps, reweightings and resamplings; iterations. The flat histogram's pilot counts too.
   - The work is converted to µs by quiet cost blocks at each fixture's own state, by M1's method, with the cost half's readers.
   - So replicate fits may run concurrently, one thread each, and the grid needs no quiet machine. Each fit's thread CPU time is read beside, as a check.
3. **References that stand outside the scored backend.**
   - *On small fixtures,* the reference is exact, by enumeration.
   - *On S10's fixtures,* a functional's reference pools the longest runs of the backends that pass D10's verdict in every group on at least 7 of 8 replicates and agree within three combined errors.
   - A backend is never scored against a reference it contributed to.
   - Where no reference stands for a backend, it is scored on variance alone, flagged, and its bias is reported as unmeasured. A reference whose own error exceeds a third of the target is "reference-limited" and scores variance alone.
4. **Two tiers of fixtures.**
   - *Tier E, enumerable:* fixture A (`BassSweepTests`' F, 2,263 structures, continuous), and `BassKernelEnumeration`'s categorical and functional fixtures.
     - Each is enumerated by `BassFlatHistogramSpace`, τ integrated over its u grid; extend that class test-side where a fixture needs it.
     - f's exact posterior mean comes from each structure's closed-form conditional mean, weighed by the enumeration.
     - Fixtures B and C add ln Z against the oracle's exact −28.589 and −30.903.
   - *Tier S:* S10's three fixtures as `BassMixing.Get` builds them (Friedman, mixed categorical, direct functional), at S10's sizes.
   - F-R1's deterministic simulators are not in. Their re-run is the user's call.
5. **Backends at their defaults:**
   - untempered MCMC (D01–D06);
   - tempering (D05's self-sizing ladder);
   - SMC for fitting (D05: 1,000 particles in each of 4 populations, 5 mutation steps);
   - the flat histogram (the WHAM warm start, its spec as Friedman's qualification sets it, lengths matched to the grid).

   F-M2's move options are not raced again (its pass, 2026-10-29).
6. **Each grid point is an ordinary fit, scaled as a user would scale one,** with 8 replicates a point. MCMC and tempering scale by s: `BurnIn` = 9,000·s, `Samples` = 1,000 and `Thin` = s for s ≥ 1, and 2,250 + 250 at Thin 1 for s = ¼. So kept draws stay at 1,000 a chain while the path lengthens.
7. **Targets, fixed before the runs:** 100 replicate effective draws for every scalar functional (Monte Carlo error at most 0.1 posterior sd), and 0.25 nats RMSE for ln Z. A length reaches a functional's target when its effective draws reach 100, read with the replicate MSE's own uncertainty.
8. **A reading, not a ruling.** No default changes, and the chip recommends none. It reports what the readings below show, per backend and regime. The run-length defaults, the early stop's default and any regime guidance are the user's call afterwards.
9. **A budget:** projected from a 1/10 trial, up to 16 CPU-hours of fits and 1 hour of quiet cost blocks, with a declared trim order.

## The measurements, fixed before any run

Declare every value below in the harness class before the first run. A change after a run is reported as a deviation, with its reason.

**Seeds.** Master seeds 9101–9108, used by no harness before. The same eight run at every length and backend.

**The lengths:**

| Backend | Grid | Notes |
| --- | --- | --- |
| Untempered MCMC | s ∈ {¼, 1, 4, 16, 64, 256} | Friedman adds s = 1,024 and 4,096 (41M sweeps a chain, the flat histogram's Friedman budget) |
| Tempering | s ∈ {¼, 1, 4, 16} | per replica; Friedman adds s = 64 |
| SMC | 250, 1,000 and 4,000 particles a population | Friedman adds 16,000 |
| Flat histogram | learning plus production sweeps a chain matched to MCMC's s = 1 and 16 on tier E, and s = 16 and 256 on tier S, split between learning and production as Friedman's qualification splits them | Friedman: its qualification's spec at `Thin` 50 and 200 |
| The early stop | MCMC and tempering at the default burn-in (9,000), with `RHatTarget` 1.01, `BatchSize` 100 and a budget of 100,000 kept draws at `Thin` 1 | `StoppingScope` All; the stop ends only the sampling, since burn-in is spent whole |

**The functionals,** grouped as D10 groups them, plus what a user reads beyond them:
- *Predictive:* σ²'s posterior mean, and f's posterior mean at D10's monitored points (`BassDiagnostics.DefaultPoints`) and at 20 test points drawn uniformly over the inputs (seed declared). For the direct functional, at the middle output location and one other.
- *Level:* the log posterior's mean.
- *Structural:* E[M], P(M = m) at the modal m and its neighbours, E[τ] where τ is sampled, each predictor's inclusion probability, and on Friedman the labelled basins' shares (a1, a2, A, B, G, S, T, by `Fixture.Basin`).
- *Sensitivity:* each predictor's total Sobol index, posterior mean (`BassSensitivity.Of`).
- *Evidence:* ln Z(1) wherever a backend reports one.
- *V,* each functional's posterior variance: exact on tier E, and from the reference pool's draws on tier S. A functional with V below 10⁻⁴ (an inclusion probability within 10⁻⁴ of 0 or 1) is reported by its absolute error alone.
- *Friedman's ln Z(1)* also stands against R2's 99.9% lower bound, −864.6, which no backend here supplied.

**What each fit records:**
- its work in its own units, and its thread CPU time and wall time;
- each functional's estimate, from the fit's own draw weights where it has them (the flat histogram), and the backend's own error and ESS for it;
- D10's verdict per group, as the fit reports it;
- for the early stop, the draws at which it stopped;
- f's posterior mean at the test points, against the fixture's true f (tier S).

**The cost blocks:**
- *Per (fixture, backend),* µs per work unit by the cost half's readers (`tests/regression/BassCostSweep.cs` and `BassCostQualification.cs`), at the end state of s = 1 and of the backend's largest length. Call them as they are; if one needs a change, copy it into your helpers.
- *Where M differs between the two,* a run's time is its work times the per-unit cost interpolated at the run's mean M.
- *The quiet rule:* no other test host, and other processes' busy cores at most 1.5, at each block's start and end. The cost half saw 1.6–6.9 cores busy with no test host running. A block that fails the rule is flagged and re-timed once quiet; if no quiet window comes, report it with its load.
- *The sensitivity moment work alone* (2026-11-01; the completion audit's G-8, the benchmark list's other open measure): a few points of `BassMoments` timed on their own, at the n, M and interaction orders the cost half read Sobol at, in the same quiet blocks.
- *When* (2026-11-01): the blocks are the one part that needs a quiet machine. Run them apart from the grid, in any quiet hour after it, keeping what each needs (the run's checkpoint, or the seeds that rebuild its end state). A later efficiency pass would leave the grid's error per unit work standing and rerun only these.

## The readings, fixed before any run

1. **Effective draws per µs, and per unit work, along the grid,** per functional, backend and fixture.
2. **The work and µs to reach the target,** per functional, and the length at which every functional in a D10 group reaches it.
3. **A bias floor:** a functional whose MSE stops falling over two grid steps while above its target.
4. **D10's verdict against the error:** for each run, its verdict per group beside whether that group's functionals reach their targets. The rates of "converged" with an error above twice the target, and of "not converged" with every target met.
5. **The self-reports' calibration:** the coverage of the reference by the estimate ± 1.96 × the backend's own error, over the replicates, and the backend's own ESS over the replicate effective draws. *Added 2026-11-01, before any run:* for the flat histogram on tier S, the replicate effective draws against the b = 1 Kish ESS, against round trips, and against crossings between two levels of ℓ₁ fixed per fixture (the reweighted medians at b = 0.10 and b = 1, read from the trial). Round trips count at each window's own ends, and the top differs by pilot, so only the fixed levels compare across runs. On seed 8017 the Kish ESS moved 4.6× across three window tops while the reported errors moved less, which four chains cannot resolve ([the loose ends' report](flat-histogram-loose-ends-brief.md#reports)). Say which of the three tracks the error.
6. **The early stop:** where it stopped, its error there, and its work against the first grid length that reached the targets.
7. **The inferential plateau** (tier S): the RMSE of f's posterior mean to the true f over the test points, per length, and the shortest length within 5% of the longest's. That is how long is long enough for predictions, beside how long is long enough for the posterior.

## The budget

- **Order:** tier E first, then the 1/10 trial on tier S (a tenth of every length, the replicates whole), projected; then tier S.
- **Up to 16 CPU-hours of fits,** concurrent up to 5 at a time, one thread each (`MaxDegreeOfParallelism` 1). Up to 1 hour of quiet cost blocks.
- **Above that, trim in this order and record it:**
  1. Friedman's MCMC at s = 4,096;
  2. tempering at s = 64 on Friedman;
  3. SMC at 16,000 particles on Friedman;
  4. the replicates from 8 to 6 at each backend's largest length;
  5. s = 256 on the mixed and direct functional fixtures.
- **No need to ask:** the user may be away, and the brief's rules decide.

## Checks before the runs

- *The currency on exact draws:* on fixture A, i.i.d. draws from the enumerated posterior give effective draws equal to the draw count, within the replicate MSE's error. That validates the estimator before any sampler is scored.
- *The references:* tier E's exact values against the enumeration's own sums; ln Z on fixtures B and C against the oracle within 10⁻³.
- *The conversion:* on one quiet fit, its thread CPU time against its work × µs per unit, within 15%.
- *The verdict as reported:* D10's verdict is read from each fit's own diagnostics, never recomputed.
- *Draws unchanged:* the eight bit-identity hashes (`BassBitIdentityTests`) pass with the harness in the tree.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Tests only.** The harness is `tests/regression/BassStatisticalQualificationTests.cs`, with helpers in `BassStatisticalQualification.cs`, behind `BASS_QUALIFICATION` = 1. If a reading truly needs a production hook, stop and say so in the report rather than add one.
- **A run survives a restart.** A Windows update's reboot ended the first attempt on 2026-11-01, before any harness existed. Append and flush each fit's row to the results CSV as the fit ends, and on a rerun skip the rows already there, keyed by fixture, backend, length and replicate, so a reboot costs only the fits in flight. A fit's seeds fix its draws, so a fit run again is the fit it would have been.
- **Other chips may work beside you.**
  - Build against the checkout as it is; if a commit changes an API your harness uses, adapt your harness; never edit another chip's files.
  - Another chip's timed readings would suffer your load: run tier E and the trial first, and start tier S's concurrent grid once no other session's test host is running.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-perf2\` to every `dotnet build` and `dotnet test`, on the command itself rather than through a variable set earlier in a chain. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-perf2\` when you finish.
- **Your files:**
  - the harness and its helpers;
  - `performance-statistical-results.csv` beside this brief, one row per fit;
  - `performance-statistical-summary.csv`, one row per fixture, backend, length and functional;
  - this brief's Reports section.
- **Commit your own files only, by explicit path** (`git commit -- paths`), never leaving anything staged, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Commit the harness, with the seeds and the grid declared, before any run. Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - no backend scored against a reference it contributed to;
  - effective draws computed against the reference, not against the replicates' own mean, wherever a reference stands;
  - V taken from the reference, not from the scored backend's draws;
  - every run's time converted at its own fixture's cost, and the conversion check reported;
  - each run's verdict as the fit reported it.
- *Coverage:*
  - a functional reference-limited on every backend;
  - a backend with no length reaching a target, reported with its bias floor;
  - a fixture where the stop never fired;
  - a hub measure this run does not read.

## Report

Under "Reports" below, newest first:
- **per fixture, the curves first:** effective draws per µs along the grid, per backend, for each D10 group's functionals and the others;
- **the work and µs to the targets,** per functional and group, with the bias floors;
- **D10's verdict against the error,** and **the self-reports' calibration,** per backend;
- **the early stop:** where it stopped, its error there, and its work;
- **the inferential plateau** on tier S;
- **the references:** which backends stood for which functionals, and what was scored on variance alone;
- **the costs:** the cost blocks, their quiet readings, and the conversion check;
- **the hub's list:** what this closes, and what is left;
- **the trims, deviations, review passes, and anything the main session must know.**

Recommend no default.

## Not in this brief

- **F-R1's re-run** on tempered ladders or longer chains, the user's call.
- **The efficiency pass,** which acts on the cost half's ranking.
- **F-M2's move options,** decided by its pass.
- **Calibration's samplers** (stage 7), whose targets need references of their own.
- **The completion audit** (deferred by the user).

## Reports

Newest first.
