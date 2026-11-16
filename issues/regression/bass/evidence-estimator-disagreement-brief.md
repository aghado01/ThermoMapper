# BASS evidence estimates from SMC and tempering

Opened 2026-11-15. This investigation concerns the marginal likelihood used to compare models. The BPPR sampler decision can proceed independently.

## What needs explaining

On the mixed generating function with a pure K = 4 ridge basis, SMC reported ln Z = 197.8 ± 0.7 and the tempered ladder's stepping-stone estimate was 194.6 ± 0.5. Neither run raised a flag. Under independent runs, the 3.2-nat gap is about 3.7 combined reported standard errors. On the smooth axis function with directional hinges, the gap was 7.0 nats (167.8 against 174.8), and SMC flagged overdispersion. The source report identifies the data, matched proper priors, budgets, flags and probe. [BPPR pure-basis readings](bppr-brief.md#pure-basis-readings-before-slice-6-2026-11-15)

## Investigation

1. Run both estimators on a small BASS model with independently known ln Z. The existing finite-structure enumeration and dense closed form supply a reference. Across independent seeds, compare signed error, observed run-to-run spread and reported Monte Carlo error, including how often the intervals cover the reference. [Existing reference](../../../tests/regression/EngineEvidenceTests.cs)
2. Reproduce both BPPR cases on the original data and priors. Repeat whole runs with fresh seeds. Increase particle count, independent SMC replicates and mutation effort separately. For stepping-stone, vary independent ladders, rung spacing and run length. Record elapsed cost, rung overlap and round trips, SMC genealogy and replicate spread, and every evidence flag. Compare the estimated log-normalizer increments on common temperature intervals where the paths match. [Evidence diagnostics](../../../src/maths/regression/spline/bass/BassEvidence.cs)
3. For the matched bases, report both each model's ln Z and each estimator's log Bayes factor. Give the uncertainty of the contrast using the actual replication design. Identify model rankings that survive changes in estimator and budget, and any for which the data leave the ranking unresolved.
4. Carry each estimator's vector of log evidences through conformal Bayesian model averaging on identical training scores, candidate grids and model priors. Record the central posterior model weights, their reported ranges, and the accepted grid values for each case. Compare the cases that change sets between estimators with each run's `EvidenceReading.Fragile` flags. In particular, count changes on cases neither run flagged. The current diagnostic tilts log evidences by twice each model's reported Monte Carlo error by default. A common shift cancels from model weights; the test must vary the contrasts. Include a close model comparison as well as one with a large evidence margin. [Model averaging and its evidence tilts](../../../src/maths/regression/conformal/ModelAveraging.cs)

## Completion

Report where the discrepancy first appears and whether it is attributable to the estimate, its error bar, the flags, or remains unresolved at the measured budgets. A proposed correction needs a reference or a reproducible failure that distinguishes the cause. Leave the measurements and seeds in repository artifacts, append the report here, and list doc changes owed.

## Report — matched contrasts, diagnostics and conformal propagation (2026-11-15–16)

The evidence chat consumed the defects worker's live `artifacts/chip-defects/evidence/{step1.txt,step2.txt,step2.txt.curves,runs.txt}` without rerunning its budget batch. Its separate investigation lives in `artifacts/evidence-contrasts/`, with a clean source archive, distinct build folder, below-normal priority and `DOTNET_PROCESSOR_COUNT=2`. It completed six whole-run replicates per estimator for the matched degree-1 axis-hinge companion on each original dataset, a dense closed-form calibration, and fixed posterior density tables for the close axis directional/axis hinge comparison and the large-margin mixed ridge/axis hinge comparison. Equal model priors, test inputs, grids and score rows stayed fixed across evidence vectors. The consumed target receipts were frozen at **2026-11-15 23:51:26 −08:00**, before conformal propagation; subsequent live target results are not silently substituted into this report. Shared checkout changes outside this brief were left alone.

The source snapshot is `a83db0c60c3486453b72956303043247379fd594`. A file-hash comparison against the defects worker's snapshot found only the sampler-default change among regression, samplers, linear algebra and RNG source files: five changed files and the added `BassSamplerDefault.cs`. Every measured fit explicitly declares SMC or tempering, which the new resolver passes through; the evidence, target, proposal and numerical implementations match. The four score fits have finished. Each uses seed 9001, two 40-rung ladders, burn 5,000, 1,000 retained draws per chain, thin 2. Density tables retain 128 evenly spaced rows across those 2,000 draws. Test inputs are eight uniform points generated with seed 9002, with each grid the known generating mean ±0.8 at step 0.025 (65 candidates). These are fixed-score propagation measurements, not predictive-coverage qualification.

The companion runs use whole-run seeds 101–106, the original data and priors, and each estimator's original budget. Matched-model log Bayes factors use paired whole-run differences by seed and their empirical standard error, retaining any covariance from shared master seeds. Individual runs' quadrature errors are only indicative because they omit that covariance. Budget means compared to the baseline companion use their separate empirical mean errors: target seeds 201 onward and companion seeds 101–106 are distinct master-seed sets. The Gamma-τ enumeration and the independent dense held-τ reference are reported separately. `.codex` denied scratch-environment creation; the old portable virtualenvs also refer to moved Python installs. A dependency-free investigation virtualenv lives with the probe artifacts; this did not affect sampling.

**Receipt gap closed for the originals:** the defects worker's existing rows retain genealogy-based errors but not the full `BassSmcReadout` (surviving families and stage reports). The evidence chat's probe writes `<case>-<variant>-<seed>.diagnostics.json` for its later companion fits. As requested by the supervisor, four isolated diagnostic replays recovered both estimators' original runs; all four normalizers and reported errors matched the saved values exactly. No existing target budget batch was changed or restarted by this chat.

### Known-reference calibration

These are different references. The defects worker's first fixture has 26,189 finite structures and integrates Gamma τ numerically through `BassTemperingTests.Exact`; its reference is −16.742427708392842. The independent dense fixture has 37 structures, holds τ = 0.7, uses zero prior centre and σ² ~ IG(3, 2), and evaluates each design's dense multivariate-t density through `EngineEvidenceTests.DenseMarginal`. Its reference is −23.631769838043557; the engine enumeration gives −23.631769838043553. Neither reference is a known answer for the continuous directional models below.

Each row uses 30 whole runs (master seeds 1–30). Bias means the average estimated log evidence minus the reference, spread is the sample SD across whole runs, and reported error is the RMS of the runs' reported Monte Carlo errors. Coverage counts intervals estimate ±1.96 reported errors; it is an observed count, not a qualification of nominal 95% coverage.

| Reference | Estimator | Signed bias | Empirical spread | RMS reported error | Reference covered | Flags |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Gamma-τ enumeration | SMC | +0.002120 | 0.017079 | 0.019514 | 27/30 | 0/30 |
| Gamma-τ enumeration | Ladder | −0.007985 | 0.046748 | 0.063897 | 28/30 | 0/30 |
| Held-τ dense formula | SMC | +0.000257 | 0.003941 | 0.004264 | 29/30 | 0/30 |
| Held-τ dense formula | Ladder | +0.000052 | 0.003215 | 0.004069 | 28/30 | 0/30 |

The empirical t intervals for mean signed bias are respectively [−0.004257, 0.008498], [−0.025441, 0.009471], [−0.001215, 0.001729] and [−0.001149, 0.001252]. None excludes zero. The measured small fixtures do not expose an evidence-constant defect or establish that either estimator's error is calibrated on the much larger directional models. Receipts: `artifacts/evidence-contrasts/{consumed-step1.txt,dense-calibration.txt,calibration-summary.json}`; probes: `DenseCalibrationProbe.cs` here and the worker's `EvidenceCalibrationProbe.cs`.

### Matched-model log Bayes factors

These use the original n = 300 datasets and matched proper priors, including the same Zellner–Siow τ prior. The default centre and σ² scale read the training response: these are the original empirical-Bayes settings, held identical between the models on a dataset, rather than a newly chosen prior. The first model is directional hinges on the axis function, and K = 4 ridge blocks on the mixed function; the second is degree-1 axis hinges. Equal model prior probabilities apply.

All entries below use the six matched whole-run seeds 101–106. Each ± value is an empirical standard error of the mean over whole runs. The log Bayes-factor error comes from the six paired differences themselves, preserving covariance from the shared master seeds. It is not the quadrature of the individual model errors. The ten-run target-ladder summaries remain available separately; the table deliberately uses only the six runs with companion fits.

| Function | Estimator | First model ln Z | Axis hinges ln Z | ln BF (first / axis hinges) | Empirical 95% t interval for ln BF |
| --- | --- | ---: | ---: | ---: | --- |
| Axis | SMC | 171.297 ± 0.681 | 174.729 ± 0.352 | −3.432 ± 0.535 | [−4.807, −2.056] |
| Axis | Ladder | 175.215 ± 1.020 | 177.944 ± 0.936 | −2.729 ± 1.589 | [−6.813, 1.356] |
| Mixed | SMC | 197.254 ± 0.218 | 148.674 ± 0.159 | +48.581 ± 0.313 | [47.777, 49.384] |
| Mixed | Ladder | 196.982 ± 0.190 | 151.629 ± 0.733 | +45.354 ± 0.581 | [43.860, 46.847] |

The mixed-function ranking survives either estimator at these budgets, although its log Bayes factor changes by several nats. The close axis comparison's two means favor axis hinges; the default ladder does not resolve that ranking across whole runs. On the matched six runs the original target's SMC-minus-ladder gap is +0.272 ± 0.292 on mixed ridges (95% t interval [−0.477, 1.021]), and −3.918 ± 1.264 on axis directional hinges ([−7.167, −0.669]). The 3.2-nat original mixed-ridge gap is not a stable offset over these new seeds. The axis discrepancy persists, and the companion's evidence also shifts between estimators. These comparisons cannot determine signed bias on either continuous model because neither has a known ln Z.

Receipts: `artifacts/evidence-contrasts/{companions.txt,estimator-summary.json,estimator-gaps.json,contrasts.json}`. `per_run` in `contrasts.json` also retains individual contrasts and indicative quadrature errors, explicitly distinguished from the empirical paired uncertainty above.

### Budgets and the original diagnostic replays

The worker varied one named control at a time. Baseline SMC is 500 particles × 4 replicates × 12 mutation steps; the alternatives are respectively 2,000 particles, 16 replicates, and 48 mutation steps. Baseline stepping-stone uses two self-sizing ladders, burn 5,000, samples 1,000, thin 2; the alternatives are eight ladders, burn/samples 20,000/4,000, or fixed 40 rungs for axis hinges and 24 for mixed ridges. The fixed-rung alternative retains adaptation of rung locations, so this is a measured denser ladder rather than a claim about a predetermined uniform spacing. Fresh baseline target seeds are 101–106 for SMC and 101–110 for ladders; raised budgets use 201–202 for SMC and 201–203 for ladders.

| Target | Setting | Whole runs | Mean ln Z ± empirical SE | Flagged runs |
| --- | --- | ---: | ---: | ---: |
| Axis directional hinges | SMC baseline | 6 | 171.297 ± 0.681 | 5 |
| | SMC particles | 2 | 174.323 ± 0.276 | 2 |
| | SMC replicates | 2 | 172.081 ± 0.911 | 1 |
| | SMC mutation | 2 | 176.244 ± 1.011 | 0 |
| | Ladder baseline | 10 | 175.286 ± 0.843 | 6 |
| | Ladder chains | 3 | 173.787 ± 0.246 | 3 |
| | Ladder length | 3 | 177.053 ± 1.055 | 1 |
| | Ladder 40 rungs | 3 | 179.121 ± 0.369 | 0 |
| Mixed ridges | SMC baseline | 6 | 197.254 ± 0.218 | 0 |
| | Ladder baseline | 10 | 196.908 ± 0.133 | 3 |
| | Ladder chains | 3 | 197.103 ± 0.167 | 3 |
| | Ladder length | 3 | 197.625 ± 0.398 | 0 |
| | Ladder 24 rungs | 3 | 197.891 ± 0.335 | 0 |

These ± values are sampling errors of the means at the measured settings, not bounds on finite-budget bias. Two raised-budget SMC runs give very weak information about tails. The axis measurements distinguish mutation/particle effort from simply adding replicates, and denser/longer ladders from adding short ladders. They do not establish a common limiting ln Z. The mixed ladder settings move the mean toward the SMC readings; the mixed SMC particle/replicate/mutation batch remained with the original worker and was not complete in this frozen receipt set. No conclusion about those unmeasured settings is implied.

The original replays use axis seed 3763 and mixed seed 3771. Axis SMC ended with **2, 3, 8 and 2 initial families**, versus particle ESS 466.75, 314.96, 471.66 and 342.38. Mixed SMC ended with **16, 18, 36 and 25 families**, versus ESS 279.41, 493.25, 353.12 and 250.25. Every original SMC stage after initialization had conditional ESS at least about 495/500; final steps started above b = 0.92. Thus good local reweighting ESS coexists with substantial cumulative ancestry loss; these are different diagnostics. The axis original did not reach exactly one surviving family. In the broader receipts, the repeated error 0.832554611… has a code-level explanation: `SmcWeights.GenealogyRelativeVariance` becomes 1 when only one family remains, and `LogNormalizerError` maps it to √ln(2). This ceiling is a limitation of that genealogy reading, not an identified posterior-draw defect.

The original axis ladders used 19/18 rungs, max pair rejection 0.63531 and zero round trips per scan; the mixed ladders used 11/12 rungs, max pair rejection 0.80379 and round trips per scan 0.0005. Full per-pair rejection/acceptance and rung locations are retained. Across all 16 original replicate/chain pairings the **only exact shared temperature endpoints are b = 0 and 1**. Summing the saved increments on that common interval reproduces every component's normalizer. There is no exact shared interior interval here on which to locate the first divergence; interpolating cumulative curves would not establish that location or its error. No single bad interior rung is established by these receipts.

Receipts: `diagnostics-replays.txt`, `replay-parity.json`, the four original `*.diagnostics.json` files, `original-diagnostics-summary.json`, and `common-interval-increments.json`, all under `artifacts/evidence-contrasts/`. The engine's fit readout exposes stage-level ancestor/family counts, not each particle's full ancestry tree; no such tree is claimed. Every run's elapsed time is in `runs.txt`, with per-fit time in the scalar rows. Timings include varying shared-host contention and are not an isolated speed comparison.

### A reproducible pooling mechanism, with a bounded correction proposal

`LadderEvidence.Pool` averages component **log** normalizers. That reduces variation of the mean log estimate as component count grows, but it does not remove the difference between E[log Ẑ] and log E[Ẑ]. This matters separately from mutation and overlap.

`PoolingReferenceProbe.cs` isolates the rule on the independent dense BASS fixture. Draw S independently from its exactly enumerated, normalized finite prior and let Ẑ = p(y | S), evaluated by the dense formula. Then E[Ẑ] is exactly the known evidence. The same finite enumeration gives:

- log E[Ẑ] = −23.631769838043557;
- E[log Ẑ] = −23.6711640205547;
- their Jensen gap = **0.0393941825111 nats**.

The probe supplies each component's exactly enumerated log SD as its error and also lets the current pool measure the between-component error. Across 30 whole experiments with nested independent-draw prefixes, raising component count from 32 to 256 to 2,048 reduced mean reported errors to about 0.04989, 0.01726 and 0.006054, but reference coverage fell from **28/30 to 9/30 to 0/30**. At 2,048 components the mean log pool is −23.66998655; the logarithm of the mean normalizer is −23.63035328, whose whole-experiment mean interval [−23.632793, −23.627914] contains the dense reference. This is a reproducible failure to recover the evidence by merely increasing the number of noisy unbiased normalizer inputs, isolated from adaptive SMC and posterior mixing.

The saved actual components show this is potentially material: over the six fresh axis SMC whole runs, mean log normalizer is 171.297 while log mean normalizer over their 24 component estimates is 175.155. Mixed-ridge values are 197.254 and 197.613. These alternate reductions are **sensitivity measurements**, not replacements validated against a known continuous-model evidence. Stepping-stone components have their own stationarity, cross-rung dependence and finite-sample biases; changing the generic pool for every caller is not justified by this witness. Adaptive SMC may have additional finite-particle error too.

The bounded candidate is an **SMC-specific normalizer-space aggregate**, with uncertainty computed on normalizers before taking the log and calibrated across whole runs against known references. It needs separate qualification from the ladder aggregate. No correction was implemented, and no posterior code changed. This mechanism cannot alone explain why the original mixed SMC reading was higher than its ladder reading. The original mixed gap's instability, the axis budget sensitivity, and the pooling mechanism should not be collapsed into a single proven cause.

Receipts: `pooling-reference.json`, `pooling-reference-summary.json`, `pooling-sensitivity.json`, and `results/pooling-reference.trx`. The earlier exploratory zero-within-error receipt is retained as `pooling-reference-zero-within.json`; the final witness above uses the exactly known component log SD.

### Fixed-score conformal propagation

`ContrastProbe.cs` carried **61 frozen evidence vectors** through the production `ModelAveraging.FullStep`, α = 0.1 and default error multiple 2. All use the same eight test cases, 65-point grids, equal model priors, and 128 retained score rows per model described above. Per-run vectors retain each engine's reported error. Mean vectors use empirical standard errors across the available whole runs; a single whole run never supplies a zero empirical error for a mean vector. Mean target ladders use all ten available baseline runs, while the paired log-BF table above deliberately uses six. Reported weight ranges are sensitivity ranges over ±2-error contrast tilts, **not confidence bounds**. A common +100 shift of both log evidences left the sets unchanged in the checked per-run vector on each dataset.

| Mean vector | Central weight | Reported range for that weight | Fragile cases / 8 |
| --- | ---: | --- | ---: |
| Axis, SMC: directional hinges | 0.031325 | [0.004077, 0.203493] | 3 |
| Axis, ladder: directional hinges | 0.065506 | [0.001990, 0.711298] | 5 |
| Mixed, SMC: axis-hinge competitor | 7.974×10⁻²² | [3.749×10⁻²², 1.696×10⁻²¹] | 0 |
| Mixed, ladder: axis-hinge competitor | 2.165×10⁻²⁰ | [3.827×10⁻²¹, 1.225×10⁻¹⁹] | 0 |

The ridge model's complementary mixed-function probability rounds to 1 in double precision; that is not a claim of certainty. Every vector's central weights, ranges, flags and **all accepted grid values** are retained in `propagation.jsonl`, not merely interval endpoints.

Across the six paired baseline seed comparisons, the axis sets changed on **14/48 case-readings**, and the mixed sets on **0/48**. No baseline change occurred where both readings had `Fragile = false` (18 such axis readings, 48 mixed). The 16 mixed baseline readings with neither engine flag also had no set change. Counts reuse eight fixed cases across run comparisons; they are not independent predictive-coverage trials.

Crossing the four available axis SMC budget means with the four ladder budget means gives 36 changes in 128 case-readings, including **one change among the 55 readings where neither `Fragile` flag warned**. The four mixed mean comparisons give zero changes in 32 readings. This is a **synthesized vector-of-means witness, not a pair of individual engine runs**. Each model's supplied error is its sample SD of whole-run log evidences divided by √R, not an average of the engines' reported errors and not an allowance for bias. Here R = 2 for the 16-replicate SMC target, R = 3 for the 40-rung ladder target, and R = 6 for each estimator's axis-hinge companion. The supplied error vectors are [0.910903, 0.352135] and [0.368678, 0.936346]. The default ±2-error tilts on these empirical mean errors are especially weakly estimated at R = 2–3.

The missed change is case index 2: directional-model central weights 0.066139 and 0.764290, with reported ranges [0.005632, 0.469665] and [0.192522, 0.977826]. Both readings mark this case non-fragile. SMC accepts the 13 grid values from 1.32246846120823 through 1.62246846120823 in steps of 0.025; the ladder additionally accepts **1.6474684612082307**. Each central weight lies outside the other reading's tilt range. This demonstrates a missed estimator/budget change in those synthesized inputs, not an error in evaluating the stated tilt vertices or a measured false-negative rate for individual runs. The mean vectors carry the union of engine flags from constituent fits, including overdispersion or saturation; absence of `Fragile` and absence of engine flags are explicitly different counts. The separate individual-run comparison counts remain 14/48 axis changes and 0/48 mixed changes, with zero changes escaping both fragility readings in those paired single-run comparisons.

Against each estimator's own baseline mean, the raised-budget axis vectors change 14/48 case-readings, with zero changes escaping both fragility readings; mixed vectors change 0/24. Receipts: `vectors.json`, `propagation.jsonl`, `propagation-comparisons.json`, `propagation-counts.json`, and the four `*.scores` files. This qualifies propagation of evidence changes through fixed scores; it does not qualify the posterior score tables or predictive coverage.

### Disposition, validation, and handoff

The independently known small-model references do not expose a constant error. The discrepancy first appears in the larger continuous-model readings and remains budget-sensitive there. The original mixed-ridge 3.2-nat difference does not survive fresh-seed averaging as a stable offset. The axis disagreement persists; genealogy loss, mutation/particle effort and rung/length changes distinguish relevant finite-budget limits, but do not supply its true ln Z or a unique defective rung. A separate, reference-based pooling-bias mechanism is established, and its correction proposal is bounded above. Flags and ±2-error tilts do not cover every observed estimator/budget change. No posterior-draw defect was established.

Validation on the archived production source: **8/8 existing `EngineEvidenceTests` and `ModelEvidenceTests`**, 30+30 dense-calibration whole runs, four score fits, 24 companion whole runs, four original diagnostic replays with exact scalar parity, the 61-vector conformal probe, and the finite-prior pooling witness. All invoked test filters passed. `results/*.trx`, `runs.txt`, configurations/resolved priors, seeds and SHA-256/size receipts in `manifest.json` identify the work. The four `*.config.json` exports are base model configurations; `effective-configs.json` expands each measured case/variant/master seed with its full sampling settings, using the engine-emitted original-run sampling templates and the recorded one-control overrides. The archived probe sources specify the dense and score configurations. No production source changed, so this is evidence-specific investigation validation, not a qualification of the other sessions' current mixed-basis/sensitivity edits or a rerun of the entire regression assembly.

The original worker still owns the mixed-ridge raised-SMC budget batch. Its results after the frozen cut are owed a new explicit receipt snapshot and, if they change vectors, a new propagation pass; `analyze.py --refresh` makes that update explicit. The remaining scientific uncertainty is finite-budget evidence bias/precision on the continuous models, including the SMC-specific pooling proposal's behavior under adaptive schedules and the ladder components' own bias. Those are not covered by the finite-prior witness or by the strong mixed-model ranking alone.

**Doc changes owed**

- `docs/Maths/Regression/BASS.md`: distinguish finite-structure/Gamma-τ enumeration from the independent held-τ dense reference; separate whole-run spread and bias from reported Monte Carlo error; state the measured budgets and conditional limits of the model rankings.
- The evidence/sampler documentation: explain genealogy error versus particle ESS and surviving families, and distinguish pooling log normalizers from pooling normalizers; qualify any SMC-specific correction separately from stepping-stone.
- Conformal model-averaging documentation: preserve the fixed-score missed-fragility example, distinguish engine flags from `EvidenceReading.Fragile`, and describe the weight ranges as sensitivity readings rather than bounds. The counterexample above does not authorize enlarging error bars by an arbitrary multiplier.
