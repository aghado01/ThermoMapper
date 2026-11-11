# E34 slice 1: the core (chip, 2026-11-09)

Commits: `8196778` (brief corrections from the slice 0 verification), `ad72e2f` (InverseWishart, OutputMarginal),
`19dda31` (K-output core), `2b907aa` (fits, Σ's prior options), `878fd19` (laws, Geweke, recovery, cost), `b03d926`
(refusal gate names KHAOS's forms), `2b5b60e` (docs/Maths/MultiOutput.md, pointer in BayesianRegression.md,
notation-ledger entry, the brief's slice 1 report). Whole regression assembly on 878fd19: 2,699 of 2,699.
The report with gates, decisions beyond the brief and follow-ups is in `issues/regression/bass/multi-output-brief.md`.

## Where things sit, for slice 2 and other chips

- **Scalar BASS untouched by construction.** `BassDesign` has a K-output side (`DesignOutputs` in `BassOutputs.cs`); every
  K-output path is gated on `_outputs is not null`. Scalar tables, pencils, follow records and spectra are untouched.
  Anyone editing `BassDesign.Decide`/`ProposeDeath`/`Refit`/`Snapshot` should keep the `_outputs` branches.
- **Target:** `BassTarget.Outputs` non-null for several outputs; `Marginal`, `Given` and `LogLikelihoodSlope` refuse
  them (slope in b is slice 2). `Defines(b)` has decision 3's rule.
- **Draws:** `BassDraw.Sigma` / `OutputCoefficients` (rows terms in ids' order, columns outputs, uncentred); `Sigma2` NaN,
  `Coefficients` empty. Readers of one response (`BassPointwise`, predict, sensitivity, conformal) refuse such draws.
- **Config:** `BassConfig.OutputCovariance` (`OutputCovariancePrior`: ScaledToResponse, Standardized, InverseWishart),
  null in JSON when unset, so existing fingerprints hold. Huang and Wand at ν = 2 is slice 2.
- **Refused in slice 1:** ridge/KHAOS priors, noise, random intercepts, ladder/SMC/flat histogram, functional response,
  offsets/weights/latents/reference/priorGram on the working response, missing outputs.
- **Engine pieces others can use:** `Maths.Distributions.InverseWishart` (density, mean, Bartlett draw through Ψ's
  factor), `SpecialFunctions.LogMultivariateGamma`, `BorderedCholesky.ProposedWhiten`.

## Geweke recheck after verification (2026-11-09) — SUPERSEDED

Finished elsewhere: `ac9726b` (fixed-structure τ grid and τ block against quadrature, proposal fits against dense refits,
12,000 independent paths as the committed gate; no sampler defect, the single successive path's slow fitted-level
exploration explains seed 20261123), then slice 2 (`41727ac`) and slice 3. See the brief's "Slice 1 Geweke recheck"
report. The text below is the chip's record up to its stop; its "next" and "uncommitted" lines no longer apply.
`artifacts/chip-e34/` keeps only the probe output and run logs (the checkpoint review links the probe); bin, obj and
the snapshot are deleted.

Worst |z| per seed (bound 4.0). The 4× length (400k draws, 600k steps, 50 batches) is the scalar Geweke gates' own convention:

| seed | 1× (100k/150k) | 4× (400k/600k) |
| --- | --- | --- |
| 20261122 | 1.78 (ln τ −0.72) | 1.13 (ln τ −0.89) |
| 20261123 | 3.38 (ln τ −3.38) | **7.25 (ln τ −7.25)**, ln S11 3.99 |
| 20261124 | **4.54** (x0 in; M 3.98; ln τ 1.69) | 2.38 (x0 in −0.20, ln τ −0.53) |
| 20261125 | 3.54 (ln S22; ln τ −2.27) | 1.92 (ln τ −1.36) |

Reading: seed 24's 1× failure shrank at 4× (signs flipped) → stickiness/undersized batch error at 1×. Seed 23's ln τ
doubled from −3.38 to −7.25, the √4 signature of a bias; ln τ z negative (successive τ too high) on 6 of 8 runs. So the
main session's trigger holds: τ's step under the K × K form needs examining before slice 2.

Ruled out by reading: `BassOutputDraws` (A rows get L_Σ·v_j, Cov = Σ·(G_c⁻¹)_jj; intercept Σ·shrink/N);
`CoefficientConstant` = 1 for multi (shared initializer, reference refused); τ's block is a plain slice step on
`target.Evaluate`.

Running when stopped: `artifacts/chip-e34/snap/tests/regression/GewekeProbe.cs` (snapshot only, partial of
MultiOutputLawTests; env GEWEKE_PROBE=<out path>, GEWEKE_SEED) reproduces seed 20261123's 4× run and writes
`artifacts/chip-e34/probe-20261123.txt`: worst |engine − matrix-t reference| over visited states, ln τ z under 50/20/10
batches, ln τ and M per 10k block, ln τ quantiles. Next: if the reference disagrees at some state → target defect there;
if not → suspect proposal forms (R′, F′ by bordered updates) at near-singular N = 12 designs with interactions, which
the enumeration law (one-way, n = 24) never reached; check a proposal's forms against explicit ones in the probe.

**Probe result (seed 20261123, 4×):** the target equals the matrix-t reference at all 600,000 visited states (worst
relative 2.8e-13), so the target is right wherever the chain goes. The bias is ln τ's LEFT TAIL: successive quantiles
0.1%/1%/5% = −10.4/−8.5/−6.8 against the prior's −15.9/−11.3/−8.0; the upper quantiles agree; z stays −7 to −10 under
10, 20 or 50 batches. The chain under-reaches small τ. Next: at frozen states, evaluate `target.Evaluate(fit, u)` over
u ∈ [−25, 0] and look for −∞/NaN/refusal below some τ (OffSupport, LogTauDensity, OutputMarginal at tiny τ, SliceSampler
stepping-out with Width 1, MaxSteps 32), M = 0 and M > 0 alike, then τ's step alone against quadrature there.

Uncommitted in the repo: `tests/regression/MultiOutputLawTests.cs` adds seeds 20261124 and 20261125 (InlineData). Do
not commit until the gate holds; then set draws/steps to 400_000/600_000 to match the scalar gates. The scratch runner is
`geweke-seed.ps1` in this session's scratchpad. `artifacts/chip-e34/` holds the snapshot and logs; delete at the end.

## Open with the main session

- Ridge/KHAOS forms for several outputs: main session's call was g-form first, the others as a second commit if
  clean, else slice 2. They are not clean (spectrum and pencil must carry K crosses), so slice 2; refused until then,
  the refusal gated in `MultiOutputCoreTests.Refuses_WhatLaterSlicesBring_AndOutputsThatCannotBeFitted`.
- Geweke ln τ z = −3.38 on seed 20261123 (−0.72 on 20261122): a third seed would settle it.
- E4's residual birth (KHAOS chip, ff2ec0d) refuses several outputs (`data.Response is not ScalarResponse`): its score
  reads a scalar residual. A K-output score (e.g. on the K × N residual matrix) is a later item, outside slice 2's row.
