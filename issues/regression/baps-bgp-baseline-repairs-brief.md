# BAPS and BGP baseline repairs — chip brief

Date: 2026-10-19. Status: chip brief, not started. Start once step 0.1 of the [BASS preparation brief](bass-preparation-brief.md) has landed on `main`. Append the report below when done.

## Context

ThermoMapper is preparing to add BASS to `Maths.Regression`. Before BASS is integrated, the existing engines must be repaired and qualified, so that pre-existing defects cannot be mistaken for BASS defects. The [adaptation prerequisites review](bass-adaptation-prerequisites-review.md) records the findings and groups them into batches B01–B06.

This chip owns **B04 (BAPS)** and **B05 (BGP)**. Neither has a code path into BASS, so this work runs in parallel with the main preparation track, which owns B01, B02, B03 and B06. The chip's closure report is an input to the main track's baseline qualification (its step 0.6).

Read first:

- [AGENTS.md](../../AGENTS.md) and [CONTRIBUTING.md](../../CONTRIBUTING.md).
- From the prerequisites review: R02, R03, R06 and R10 in the findings table; the sections "BAPS tail and Gibbs conditioning" and "BGP scaling qualification"; rows B04 and B05 of the remediation register, with its evidence paragraph.

## Why step 0.1 comes first

Step 0.1 changes how non-finite R̂ is treated in [ChainEnsemble](../../src/maths/samplers/ensemble/ChainEnsemble.cs) and [ChainDiagnostics](../../src/maths/samplers/ensemble/ChainDiagnostics.cs). [BgpSampler](../../src/maths/regression/bgp/BgpSampler.cs#L92) runs on `ChainEnsemble`, and [BapsSampler](../../src/maths/regression/spline/baps/BapsSampler.cs#L160) reports `ChainDiagnostics.RHat`. Starting after that change means the chip's checks run against the repaired diagnostics.

## Scope

### B05: BGP

1. **R06 protection (test only).** Add concurrent-versus-serial checks for [GpRegression.Fit](../../src/maths/regression/bgp/GpRegression.cs#L80) with immutable inputs and a stateless kernel. The local scratch fix is already in place; the check keeps it there.
2. **R10 support.**
   - Reproduce the empty support: two one-dimensional inputs at 0 and 100 with the default `gamma1 = 0.25` give a lower bandwidth bound of 2500 against the fixed upper bound of 1 in [EmpiricalBayesBandwidthPrior](../../src/maths/regression/bgp/EmpiricalBayesBandwidthPrior.cs#L67).
   - Then reject empty or non-finite support before the chain starts, with a clear error.
   - Define and document the input units, and one fitted transform shared by training and prediction.
   - Cover boundary and constant-coordinate inputs.
3. **R10 wording.** In the BGP section of [BayesianRegression.md](../../docs/Maths/BayesianRegression.md) and in [bgp-twcd-lemmas.md](../../lean/Protolemmata/BARS/bgp-twcd-lemmas.md#L67), replace the claim that coordinate scaling breaks the approximation theorem:
   - A fixed nonsingular affine map preserves the rates and changes only the constants.
   - Empirical min-max transforms depend on sample size and need further control.

   Claim no new proof.

### B04: BAPS

4. **R03 hierarchical scan.**
   - Reproduce the stale statistic in [HierarchicalPSpline](../../src/maths/regression/spline/baps/HierarchicalPSpline.cs#L118): the deviation sum is formed with the old population coefficients, those coefficients are redrawn at line 133, and the deviation variance is then updated from the stale sum at line 136.
   - Repair by recomputing the deviation statistic after the population draw, or by another order verified against the full joint density.
   - Check each full conditional against its closed form on a tiny model, and the stationary behavior against an independent reference.
5. **R02 `MarginalEvidence` tail.**
   - Reproduce the improper target: on a small design, the profiled evidence of [PenalizedSpline](../../src/maths/regression/spline/baps/PenalizedSpline.cs#L81) approaches a constant as `log λ` grows, and [BapsSampler](../../src/maths/regression/spline/baps/BapsSampler.cs#L197) samples it under a flat prior on `log λ`.
   - Then write a decision memo and stop for the user's choice. Options:
     - a declared proper prior on `log λ`;
     - a declared bounded support, stated as part of the prior;
     - a separate evidence-optimization mode beside the Bayesian sampler.
   - Implement the chosen option. Because the target is one-dimensional in `log λ`, check it against quadrature.
   - `Gibbs`, `AdaptiveGibbs` and `MarginalEvidence` use different priors, so agreement between them is not a criterion.
6. **B01 policy in BAPS display.** Make BAPS's reported R̂ follow the step 0.1 policy: unavailable values are shown as unavailable and never as converged.

## Order of work

1. Write the failing checks and reproductions for items 2, 4 and 5, plus the R06 protection test, and record each failure at its parent revision.
2. Send the R02 decision memo. Include the R10 transform policy as well if the public behavior changes.
3. Make the repairs for R03 (item 4) and R10 (items 2 and 3), and apply item 6.
4. Repair R02 (item 5) once the decision is made.
5. Run the suite with `scripts/fact-harness.ps1 -Project tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj`, then write the report.

## Boundaries

**May edit:**

- `src/maths/regression/spline/baps/**` and `src/maths/regression/bgp/**`.
- Their tests under `tests/regression/`. New test files take a `Baps`, `Bgp`, `PenalizedSpline` or `HierarchicalPSpline` prefix.
- The BAPS and BGP sections of `docs/Maths/BayesianRegression.md`, including their bullets under "Status and caveats".
- `lean/Protolemmata/BARS/bgp-twcd-lemmas.md`.

**Leave untouched** (owned by the main track):

- `src/maths/samplers/**`, `linalg/**`, `distributions/**`, `rng/**`, `regression/spline/bars/**` and the shared spline files.
- `docs/Maths/Samplers.md`.

**Factor calls.** Keep the existing `BandCholesky` and `CholeskyDecomposition` calls as they are. The main track adds a status API beside them and migrates these call sites in its cleanup sweep. If a repair needs a change outside these boundaries, record it in the report as a request instead.

## Working rules

- **Check first.** Each repair begins with an independent check or counterexample, and the check and the fix land in one commit.
- **Commits.** One concern per commit, on `main`. The main track works in the same checkout, so stage files by explicit path and commit nothing outside these boundaries.
- **Semantics before tests.** Settle statistical semantics with the user before changing any test expectation or public behavior. Known-incorrect outputs become defect fixtures, not gold standards.
- **Stay in scope.** No BASS work, no unrelated subsystem fixes, and no subagent dispatch.

## Done when

- B04 and B05 meet the closure evidence in the remediation register.
- The regression suite passes at a recorded revision.
- The report below is complete.

## Report

For each item, record:

- the reproduction or counterexample;
- the chosen specification;
- affected paths;
- commits;
- independent checks with results;
- remaining limitations;
- requests for the main track.
