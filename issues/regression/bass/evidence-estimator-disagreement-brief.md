# BASS evidence estimates from SMC and tempering

Opened 2026-11-15. This investigation concerns the marginal likelihood used to compare models. The BPPR sampler decision can proceed independently.

## What needs explaining

On the mixed generating function with a pure K = 4 ridge basis, SMC reported ln Z = 197.8 ± 0.7 and the tempered ladder's stepping-stone estimate was 194.6 ± 0.5. Neither run raised a flag. Under independent runs, the 3.2-nat gap is about 3.7 combined reported standard errors. On the smooth axis function with directional hinges, the gap was 7.0 nats (167.8 against 174.8), and SMC flagged overdispersion. The source report identifies the data, matched proper priors, budgets, flags and probe. [BPPR pure-basis readings](bppr-brief.md#pure-basis-readings-before-slice-6-2026-11-15)

## Investigation

1. Run both estimators on a small BASS model with independently known ln Z. The existing finite-structure enumeration and dense closed form supply a reference. Across independent seeds, compare signed error, observed run-to-run spread and reported Monte Carlo error, including how often the intervals cover the reference. [Existing reference](../../../tests/regression/EngineEvidenceTests.cs)
2. Reproduce both BPPR cases on the original data and priors. Repeat whole runs with fresh seeds. Increase particle count, independent SMC replicates and mutation effort separately. For stepping-stone, vary independent ladders, rung spacing and run length. Record elapsed cost, rung overlap and round trips, SMC genealogy and replicate spread, and every evidence flag. Compare the estimated log-normalizer increments on common temperature intervals where the paths match. [Evidence diagnostics](../../../src/maths/regression/spline/bass/BassEvidence.cs)
3. For the matched bases, report both each model's ln Z and each estimator's log Bayes factor. Give the uncertainty of the contrast using the actual replication design. Identify model rankings that survive changes in estimator and budget, and any for which the data leave the ranking unresolved.

## Completion

Report where the discrepancy first appears and whether it is attributable to the estimate, its error bar, the flags, or remains unresolved at the measured budgets. A proposed correction needs a reference or a reproducible failure that distinguishes the cause. Leave the measurements and seeds in repository artifacts, append the report here, and list doc changes owed.
