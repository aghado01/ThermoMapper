# Where from-prior SMC leaves the flat histogram's curve on Friedman (chip brief)

Date: 2026-11-01. A reading, tests only: it changes no code under `src/`. It recommends no default. `docs/Maths` and the checklists are the main session's, which verifies and records the report.

- *Why.* The flat histogram (N17) computes Friedman's structural posterior on three master seeds, ln Z(1) ≈ −841, where from-prior SMC's best replicate reads about −855. The loose ends' item 3 ([its report](flat-histogram-loose-ends-brief.md#reports)) read SMC only up to b = 0.10, and found:
  - on the linear form SMC loses 5.6 nats by then (5 sweeps; 3.9 at 20), its families down to 188–217 of 4,000;
  - on the integrated form, D05's path, about 0.4 nats (−150.59 ± 0.14 and −150.58 ± 0.03 against −150.17 ± 0.02), against one flat-histogram run reweighted.
- So most of SMC's shortfall on its own path accrues above b = 0.10, and nobody has read where. The earlier work puts Friedman's band at b ≈ 0.16–0.32 ([BayesianRegression.md](../../../docs/Maths/BayesianRegression.md), "Friedman's band, measured").
- *The question:* on Friedman, at which b does from-prior SMC fall off the flat histogram's ln Z(b) curve, and is the 0.4 nats it loses by b = 0.10 on the integrated form SMC's or one seed's?

Read first: item 3's test, `Friedman_SmcBelowTheLowerSplit_OnBothForms` (`tests/regression/BassFlatHistogramQualificationTests.cs`), and its report.

## The reading

Every part is a test behind `BASS_QUALIFICATION` = 1 in `BassFlatHistogramQualificationTests.cs`, committed with its marks and seeds before any run.

- *Marks* (`BandMarks`): b = 0.03, 0.05, 0.07, 0.10, 0.13, 0.16, 0.20, 0.25, 0.32, 0.40, 0.50, 0.75, 1, on both forms.
- *SMC* (`Friedman_SmcBand_FromThePrior`): from the prior, a stage forced at each mark, 4 replicates × 4,000 particles at C4's defaults (master seed 8039, as item 3), at 5 and at 20 mutation sweeps a stage. Each replicate's ln Ẑ(b) and families at each mark, and the stages and acceptance between marks. To say what forcing the marks costs, an unforced run to b = 1 on the same seeds, on each form at each sweep count.
- *The reference* (`Friedman_SmcBand_FlatHistogramOnBothForms`): the flat histogram's arm (i) on the three master seeds the qualification stands on, at their settings: 8017 at `Thin` 200, 8217 and 8417 at `Thin` 400. Each must reproduce its ln Ẑ(1), −841.36, −841.55 and −841.14. Each run's draws are reweighted to both forms as item 3 does: Ẑ(b) = Σ m_b(S, τ)/w(ℓ₁) over Σ 1/w(ℓ₁), with m_b = e^{b·ℓ₁} on the linear form and the integrated form's tempered marginal on the other.
- *Importance sampling from the prior* (`Friedman_SmcBand_PriorImportanceSampling`): item 3's 200,000 exact draws at the same marks on both forms, as context, with the Kish ESS of its weights.

*Checks* (an assertion fails the part):
- each flat-histogram rerun reproduces its seed's ln Ẑ(1) to two places (the full value is printed);
- each draw's rebuilt ℓ₁ is its level (the scoring throws otherwise);
- on the linear form the reweighting equals the run's curve at every mark, to 6 places;
- at b = 1 the two forms' reweightings coincide;
- on seed 8017 the first four marks reproduce item 3's flat-histogram column on both forms, to two places;
- each SMC replicate's ln Ẑ at b = 1 is its run's final one.

*A change to item 3's test code, behaviour-preserving:* its reweighting moved into a static `FlatOnForm`, and `DrawsOnTheIntegratedForm` takes the marks (item 3's by default). The seed-8017 check above is what shows it held.

## How the numbers will be read

Written before the runs, so the reading is not fitted to them. The tables carry everything; these are conventions for the sentences.

- *Reference.* On each form, ref(b) is the mean over the three seeds of their reweighted ln Ẑ(b); the seeds' range at b is its spread, and each seed's chains give an error. A replicate's deficit is d(b) = its ln Ẑ(b) − ref(b); the pooled deficit is the mean over the four.
- *Falls off at mark m:* the first mark from which the pooled deficit stays below −0.25 nats at every later mark. A replicate's own is read on its own deficit.
- *Where a deficit accrues:* the increments Δd between consecutive marks (the first from b = 0, where d = 0, to 0.03), which sum to the final deficit D = d(1). For a replicate with D below −0.25: *one drop* if one interval's increment is at least 2/3 of D; *spread* if the largest is under 1/3; *between* otherwise. Intervals inside the band (0.16–0.32) are named as such. Shares are read as they come; a recovery (a positive increment) is a recovery.
- *Do 20 sweeps move it:* at a mark, if the pooled deficits at 5 and 20 sweeps differ by more than 0.25 nats and by more than twice the two standard errors added in quadrature.
- *The 0.4 nats at b = 0.10 on the integrated form:* for each seed s, G_s = its ln Ẑ(0.10) − SMC's pooled value (at each sweep count). *SMC's* if all three G_s exceed 0.25; *one seed's* if one does and the other two stay under 0.15; otherwise mixed, and said so with the three values.
- *What forcing the marks costs:* the forced and unforced pooled ln Ẑ(1) compared against their errors added in quadrature, and replicate by replicate (the same seeds); stages and final families beside. "A little" means within those errors; if not, say so.
- *Where importance sampling still reads:* the marks at which its value lies within three of its errors of ref(b) on both forms, with the Kish ESS beside. Context only.
- *Kept apart:* what the runs show goes in one list, what is made of it in another. No default is recommended.

## Working rules

- Work on `main`. Commit by explicit path (`git commit -- paths`), leaving nothing staged: the index is shared. Push nothing. Commits end with the `Co-Authored-By` line the task gave.
- Build and run into `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-smc-band\`, set on each `dotnet` command itself. Delete that folder when done. No subagents, no worktrees.
- Do not run it while the statistical-half chip times its cost blocks.
- Leave the user's uncommitted notes under `issues/regression/` alone.
- Ends with a fresh-eyes review, fixed before the report.

## Reports

Newest first.
