# BASS completion report

Date: 2026-11-01, updated 2026-11-02 when the gaps closed. Written by the audit in [bass/completion-audit-brief.md](bass/completion-audit-brief.md), against the contract in [the hub](bass-full-port-brief.md#scope-and-completion-contract). The row-by-row detail, with each test and reference, is in [bass/completion-audit-matrix.md](bass/completion-audit-matrix.md).

## Is the contract fulfilled?

**Not yet in full, and only by two benchmark measures. Every capability P01–P15 exists and its central behaviour is verified, and every item of R's surface is now built or left out by a recorded decision: no row of the matrix is a gap.** No statistical defect was found. All 215 rows of the matrix have a verdict; the probes this audit added found no wrong number.

What stands between: **the benchmark list's two open measures,** statistical error per unit time and the sensitivity moment work alone. Both are in [the statistical half](bass/performance-statistical-brief.md), not run, and wait on a quiet window.

What the audit first found open, and how it closed on the user's rulings (the gap list has each commit):

1. **Built:** the `other` member of basis sensitivity (P05 names it), conditioning on one of several output variables (P04, `func.var`), selecting draws for sensitivity (`mcmc.use`), a declared ladder as a list (`temp.ladder`), and `predict.bassBasis`'s two switches.
2. **Left out by decision:** the swap start (`start.temper`).

Rendering and text formatting remain above maths; the contract lists them as such, so they are not gaps.

**My reading, apart from what I found.** On R's default paths the port has each capability, and it corrects the inventory's defects (18 rows in the matrix) with fixtures on both sides. The largest honest limit is not in this list: Friedman's posterior, where the default chains do not agree at practical lengths (below). That limit is measured and stated, and the contract does not promise otherwise.

## The audited revision and the suites

- **Audited revision:** `8f4d2da0081c7f17ecb524d5514a27f994191d02`. R: BASS 1.3.1, `26710146c623a258854162367a939041efedd091`.
- **Suites at the audit,** once at the end, one at a time, at below-normal priority with `-m:2` and `xUnit.MaxParallelThreads=4`, `BASS_QUALIFICATION` and `SEED_SWEEP` unset:
  - `Maths.Samplers.Tests`: 284 of 284, 1.3 minutes.
  - `Maths.Regression.Tests`: 1,507, of which 1,506 passed and 1 is skipped on purpose (G-2), 5.1 minutes.
  - The BASS classes of `Maths.Oracle.Tests` (`FullyQualifiedName~Maths.Oracle.Tests.Bass`): 90 of 90, 6.1 minutes.
- **Later commits touching BASS:** this audit's `348ad80` (probes) and `ea6bb46` (matrix), and another session's `cb0400b`, a flat-histogram refusal: a b = 1 rung on one level refuses K instead of overflowing it (`BassFlatHistogramRun.cs`, `Wham.cs`, tests). It touches no row's verdict: the flat histogram is optional and beyond the contract. The other commits since `8f4d2da` are SPC, graphs and persistent homology.
- **The tree they ran on:** the main checkout, not a clean snapshot, so that the R oracle's library was present. It began at `121005e` with `cb0400b`'s five files as another session's uncommitted edits, which it committed during the run; it ended clean in those paths at `8614a95`. So these counts are the audited revision plus the probes and that refusal. Under `src/maths/` they differ from the audited revision by `cb0400b` alone.
- **Suites with the gaps closed, at `b844f53`,** the same way: `Maths.Samplers.Tests` 284 of 284; `Maths.Regression.Tests` 1,544 of 1,544, none skipped; the BASS classes of `Maths.Oracle.Tests` 90 of 90, R present. Since the audited revision, BASS's source has changed only by `cb0400b` above, one comment in `8ac7b05`, and the gaps' own commits.

`oracle/` tests return early without R; R was present for this run (the parity tests' durations show them running).

## The qualified existing-engine revision and B01–B06

**Qualified revision:** `cf37584` (2026-10-28), recorded in the [numerics checklist](bass/port-numerics-checklist.md#baseline-requalified-before-stage-4--2026-10-28): Regression 438, Samplers 137, LinAlg 118, Distributions 93, all passing, each commit since the last baseline accounted for. Those counts are records; this audit read the closure tests below.

| Batch | Findings | Closed at | Test read |
|---|---|---|---|
| B01 Diagnostics | R09 | `4c97de2` | `DiagnosticsPolicyTests.RHat_ChainsFrozenAtDifferentValues_IsFrozen`, `RHat_OneChainOrOneDraw_IsUnavailable`; `ConvergenceVerdictTests.FrozenAndUnavailableReadings_Fail` |
| B02 Mutable state | R05, R07 | `4c97de2` | `KernelOwnershipTests`; `AdaptationAcrossResamplingTests.AdaptiveTau_IsUnaffectedByAResamplerThatKeepsTheWeights` |
| B03 BARS contracts | R01, R04, R08, R12, R13, R14 | `4c97de2` | `CountLawTests.KnotCount_FollowsTheCountPrior`; `BarsDrawTests.Fit_IsTheMeanOfTheDrawsConditionalMeans_UnderRobustWeights`; `DesignAdmissionTests.NearlyCollinearDesign_IsRejected`; `WeightedNormalMarginalTests.Marginal_ZeroResponse_IsRejected_WhileAnInterpolatingFitIsNot`; `SpanCoverageTests.NestedPeaks_CoverageStaysAProbability` |
| B04 BAPS | R02, R03 | `fbe6d0d` | `BapsMarginalEvidenceTests.ProfiledEvidence_Plateaus_SoAFlatPriorIsImproper`; `HierarchicalPSplineConditionalTests.VarianceDraws_MatchTheirFullConditionals` |
| B05 BGP | R06, R10 | `fbe6d0d` | `BgpConcurrencyTests.ConcurrentFits_MatchSerialFits_BitForBit`; `BgpSupportTests.EmptySupport_IsRejected_BeforeAnyChainStarts` |
| B06 Numerical | R08, R11 | `4c97de2` | `DenseCholeskyTests.NonPositiveDefinite_StopsAtTheFailingPivot_WithoutFlooring`; `BorderedCholeskyTests.RejectedProposals_LeaveTheAcceptedFactorUntouched` |

**Baseline repairs and port changes stay apart:** the repairs above are the commits up to `cf37584`. Everything under `spline/bass/` and its shared additions is the port. **Thirteen commits under the qualified engines' folders follow `cf37584`** (N14's verdict and BARS's groups, `DrawSample`'s two quantile fixes, the ladder's resize, the seams, the inverse temperature's rename), each with its own report or bit-identity pin, and no requalification has run since. The hub asks for one before the first integration tests only, so this is a note. `Maths.LinAlg.Tests` and `Maths.Distributions.Tests` were not among this audit's suites.

## Fulfilled capability IDs

| ID | R entry | Verdict | What limits it |
|---|---|---|---|
| P01 | `bass` | Fulfilled, limited | Tempering and log-posterior records aggregated, by the user's ruling. `start.temper` and `s2.lower` left out by decision |
| P02 | `bassPCA` | Fulfilled | |
| P03 | `bassBasis` | Fulfilled | |
| P04 | `sobol` | Fulfilled | |
| P05 | `sobolBasis` | Fulfilled | |
| P06 | `calibrate.bassBasis` | Fulfilled | Cold rung's draws only, by the user's ruling; the adaptation scale is on the block |
| P07 | `predict.bass` | Fulfilled | |
| P08 | `predict.bassBasis` | Fulfilled | |
| P09–P11 | three `plot` methods | Fulfilled, data | Rendering remains |
| P12–P15 | `print`, `summary` | Fulfilled, data | Formatting remains |

Three response representations crossed with three input kinds, and each output, have a test in every cell that applies (the matrix, 1i), the supplied basis apart from PCA. Six of the 12 rows had no test on these outputs before this audit's probes: categorical-only inputs under each of the four representations, mixed inputs under a direct response and under a supplied basis.

## The resolved statistical decisions

| ID | As built | Evidence (test read) |
|---|---|---|
| D01 | Unordered terms, one ordering factor; knots on distinct values; global F; move probabilities per state; Wallenius by a positive recursion; R's subset measure, the paper's selectable | `BassSweepTests.ConstantLikelihood_SamplesThePriorOnF`; `BassProposalTests.BirthProbability_ByHand` |
| D02 | Global truncation, λ integrated out, λ given M at retention | `BassProposalTests.CountPrior_IsTheNegativeBinomial` |
| D03 | τ's unit-information family, a = ½; level centring; σ² scaled; g-form and normalized ridge | `BassTargetTests.TauPrior_IsTheUnitInformationGamma`; `oracle/BassTargetParityTests.Marginals_MatchTheOraclesPosterior`; `BassRidgeTests.Marginal_MatchesTheDenseComputation` |
| D04 | Full collapse, τ by slice, composition at retention | `BassSweepTests.FixedStructure_TauAndComposition_MatchQuadratureAndAJointGibbsReference` |
| D05 | Likelihood path, integrated by default; DEO; reference rung opt-in; stepping-stone; SMC on one sequence | `oracle/BassTemperingParityTests.Swap_BalancesOnTheOraclesRungLaws`; `oracle/BassSmcParityTests.Fit_NormalizerMatchesTheOraclesLogZ` |
| D06 | Order-free angle gate on the centred Gram, s_min 10⁻¹⁰; support count; a ridge is a prior | `BassDesignTests.TheGate_IsOrderFree_WhereANewestColumnRuleIsNot`; `oracle/BassAdmissionParityTests.ANearConstantIndicator_KeepsItsAngle` |
| D07 | Leave-one-case-out residuals by default; noise per case; physical units; fixed axes | `BassOutputDecompositionTests.LeaveOneOut_EqualsExplicitRefits`; `BassPredictBasisTests.Covariance_InPhysicalUnits_MatchesDenseComputations` |
| D08 | Uniform default, empirical option; Student by quadrature; zero variance undefined | `BassSensitivityMomentTests.TermMoments_MatchTensorQuadratureOfTheBasisColumns`; `BassSensitivityTests.ZeroVariance_LeavesTheIndicesUndefined` |
| D09 | Cut posterior; declared formulations; half-Cauchy; declared discrepancy scale; categorical parameters enumerated | `BassCalibrationConditionalTests.Target_IsTheDenseGaussianWithEveryPriorWrittenOut`; `BassCalibrationCutSmcTests.Pooled_IsTheExactEqualMixture_ByQuadrature` |
| D10 | A verdict on every run; three groups; a declared stopping scope; lengths as a budget; adaptation frozen after burn-in | `BassVerdictTests.Fit_CarriesTheVerdict_AndTheModeSeparation`; `ConvergenceVerdictTests.FrozenAndUnavailableReadings_Fail` |

**The defaults, mechanically:** all 28 BASS rows of the Configuration table match the code (the matrix, 3b). `docs/Maths/Samplers.md` has no BASS row, and the rows BASS reads match.

**The register's Corrected entries** (five), and the inventory's defects (18 rows): each keeps R's value in the oracle and asserts the corrected one in the C# (the matrix, 2b and 3c).

## The remaining optional backends and extensions

- **Built and optional, tracked apart from capability completion:** SMC for fitting (C4), SMC for cut calibration (C5), the flat-histogram sampler (N17 Phase 2).
- **Built, off by default:** F-M2's local-move options; no arm cleared the pass.
- **Built, no consumer yet:** the extension seams F-X1, F-X2, F-X3, F-X6, N-X1, N-X2, N-X3, O-X4.
- **Not built:** E1–E18 ([extensions](bass/extensions-brief.md)); rendering and formatting above maths; the efficiency pass (priced, prefix sums first).

## The measured limitations

These are records this audit read, each where it was measured; none was rerun.

| Limit | Reading | Where |
|---|---|---|
| Friedman's posterior, MCMC | Chains hold their own bands of M; births and deaths accept about 1.5%; R̂ of M 1.5–2.5 at 5,000 + 2,000 sweeps; every D10 group fails | S10 and F-M2's pass, [fitting checklist](bass/port-fitting-checklist.md#reports) |
| Friedman's posterior, the reference | The flat histogram computes it, replicated on three seeds: a1 about 66–70%, ln Ẑ about −841 | [warm-start](bass/flat-histogram-warmstart-brief.md), [loose ends](bass/flat-histogram-loose-ends-brief.md) |
| Friedman under SMC | Replicates disagree 13–66 nats and end in 1–2 families at S10's sizes; band runs open | [fitting SMC](bass/stage8-fitting-smc-brief.md), [band](bass/friedman-smc-band-brief.md) |
| The deterministic profile, F-R1 | No arm clears: predictive R̂ above 1.05 in all 144 fits, least 1.23; ESS about 5 of 8,000 | [its report](bass/deterministic-profile-brief.md) |
| N15's cost | Construction 97–99.8% of a σ² proposal; 1.4 ms at m = 20, k = 39 | [cost half](bass/performance-cost-brief.md#reports) |
| The change at the cap | 54–61% of a sweep at n ≈ 1,000, 79–96% from n = 10,000 | [cost half](bass/performance-cost-brief.md#reports) |
| Spectral gap | Exponential in N on small fixtures, all of it at the count cap | [its brief](bass/spectral-gap-brief.md) |
| Run lengths | R's 9,000 + 1,000 stay provisional; D10 made them a budget | The doc's status section; the statistical half, not run |
| One training case | Refused: every input would be constant, which D01 and D07 refuse. Two cases fit | Probe `OneCase_IsRefusedWithItsReason_AndTwoCasesFit` |

## The gap list, ranked

### Contract gaps

Ranked as the audit found them. G-1 to G-6 are closed; G-7 and G-8 stay open.

| # | Gap | Now | Where |
|---|---|---|---|
| G-1 | **`other`, the residual of basis sensitivity.** P05 names it | **Built,** `0ab1e4d`: `SobolDecomposition.OtherVariance`, `OtherIndex`, `MeanOtherIndex`, and the summary's two fields. The variance less the subsets reported, which agrees with R's where `nind` is null | Matrix F11 |
| G-2 | **Conditioning on one of several output variables** (`func.var`) | **Built,** `4c34948`: `SensitivitySpec.ConditionOn` holds the variables it names at the coordinates and integrates the others. The output side's normalizer stays with the coordinate. The probe that was skipped now runs | Matrix E06, E07 |
| G-3 | **`mcmc.use` on `sobol` and `sobolBasis`** | **Built,** `1364cbb`, `b844f53`: `SensitivitySpec.Draws`, as `BassPredictOptions.Draws`, over a basis fit's joint draws too. Effects evaluate only the draws selected | Matrix E05, F04 |
| G-4 | **`temp.ladder` as a list** | **Built,** `3a83d0a`: `TemperingSpec.InverseTemperatures`, in place of the count and the hot end. A list of the built rungs runs bit-identical to them. A tempered config's fingerprint moves, because its JSON gains the field, so a tempered checkpoint from before this commit stops resuming. An untempered config's holds: its JSON has `"Tempering":null`, and `BassMovePassTests`' pin did not move. A list runs as given by default (the user's ruling, 2026-11-02, `3424912`), and `Adaptive = true` re-spaces it; a tempered config's JSON moved again with that, so a tempered checkpoint from before it stops resuming | Matrix A22 |
| G-5 | **`start.temper`** | **Left out by decision,** the user's ruling: it changes no target, and the self-sizing ladder reads its swap rejections during burn-in (`port-fitting-brief.md:487`) | Matrix A23 |
| G-6 | **`predict.bassBasis`'s `nugget` and `trunc.error`** | **Built,** `26befff`: `BassPredictOptions.BasisParts` takes the noise or the truncation alone, and the covariance follows. D07's full draw stays the default, bit for bit | Matrix D-p09 |
| G-7 | **Statistical error per unit time.** Known | **Open.** Long runs on the user's schedule; `BASS_QUALIFICATION` left unset here | Matrix part 5; [its brief](bass/performance-statistical-brief.md) |
| G-8 | **The sensitivity moment work alone.** Measured only inside Sobol and effects per draw. Known | **Open.** In the statistical half's cost blocks, with G-7 | Matrix part 5; [its brief](bass/performance-statistical-brief.md) |

### Built but untested at the start of this audit

Each is **now covered by a probe that passed**, so none remains in the matrix as untested. They are listed because they were not covered by the work that built them.

| Path | Probe |
|---|---|
| Restart across a basis fit (each component's checkpoint, then `BassComponents.Assemble`) | `Cell_Basis_ReconstructsDecomposesAndRestarts` |
| Categorical-only fits, scalar, direct and through a basis; mixed inputs with a direct response; sensitivity and effects on them | `Cell_Direct_PredictsDecomposesAndRestarts`, `Cell_Basis_ReconstructsDecomposesAndRestarts` |
| A supplied nonorthogonal basis through a fit | `Cell_Basis_*(*, Supplied)` |
| More than 64 category levels, through a fit | `Levels_AboveSixtyFour_DrawAndScoreUnderBothMeasures`, `Levels_AboveSixtyFour_FitPredictMeasureAndRestart`. There is no bitset, so no 64-level seam |
| A wholly intercept-only fit through every view, against its closed form | `InterceptOnly_PosteriorIsItsClosedForm`, `InterceptOnly_EveryViewReadsIt` |
| One draw, through every view | `OneDraw_EveryViewReadsIt` |
| One predictor's Sobol, which R refuses | `OnePredictor_FitPredictAndSensitivity` |
| Two output variables through a fit; degrees 2 and 3 at new inputs | `TwoOutputVariables_PredictAndIntegrate`, `HigherDegree_PredictsDecomposesAndHasEffects` |
| A measure wider than the training range; truncated measures and a zero-probability level through a fit | `WideMeasure_ExtrapolatesTheHinges`, `TruncatedMeasure_SensitivityMatchesQuadrature`, `ZeroProbabilityLevel_SensitivityMatchesEnumeration` |
| Declared `npart.func`, `w1`, `w2`, `h1`, `h2` and interaction caps below the variables', through the config | `DeclaredControls_ReachTheModelAndTheRetainedTerms` |
| Full PCA rank through a fit | `FullPcaRank_FitsPredictsAndDecomposes` |

### Record drift

The main session has fixed these (the brief's Reports). The extension seams' write-up landed on 2026-11-02 ([the section](../../docs/Maths/BayesianRegression.md#bass-the-extension-seams)), and the rows now name the seams' types. The other public types stay as listed, because the rows name entry points. The matrix's part 6 has the lists.

- **The source map** omits `BassObservationModel.cs`, `DesignKernel.cs`, `DesignKernels.cs` (doc line 2138). The **extension seams' write-up** in `docs/Maths` was owed, listed once.
- **The API rows** name no removed thing, but at the audit 71 of 242 public types appeared in no row (60 of 243 on 2026-11-02, with the seams named): 10 seams, 31 a caller reads, 30 readouts and building blocks. `BassComponents.Assemble`, which carries basis restart, is unnamed.
- **Statements found false:** `port-fitting-brief.md:503` (a cached training fit: none is built); `opus-bass-completion-remainder.md:1` ("Everything R's BASS can do is now built"); `BassRunOptions.cs:15` (tempering and SMC "when they join"); `completion-audit-brief.md` ("E1–E15": the checklist holds E18); `port-outputs-brief.md:82` ("condition on one functional coordinate", which G-2 qualifies).
- **The open boxes** N15, N17, F-M2 and F-X4 are not stale: each is open by its own condition (the matrix, 6a).

## The rulings on what the brief left open

The audit left three questions open. The user ruled on them on 2026-11-01, recorded in `port-fitting-brief.md` ("Not carried"), `deterministic-profile-brief.md` ("For a rerun") and `docs/Maths/BayesianRegression.md`.

1. **Are G-4, G-5 and G-6 gaps or recorded choices?** G-4 and G-6 are in the contract and built, with D07's full draw staying the default. G-5 is left out by decision.
2. **Is G-3 in the contract?** Yes, and it is built as `BassPredictOptions.Draws` does it.
3. **Do the per-iteration records count?** No. `swap`, `temp.val`, a `log.post.cold` trace and calibration's hot rungs stay limited, not gaps (matrix B-r10, B-r11, G14): they are working diagnostics that the readouts stand in for.

One more, which the main session raised: **`s2.lower`** (matrix A14) stays out of the contract shape, with nothing to build. A truncated σ² prior joins as a pre-registered arm when F-R1 reruns.

## How this was checked

- **The standard:** each Fulfilled row names a test read for this audit. About 80 cited tests were read whole. For the other 192, a script kept each one's first nine assertions and comments, and those were read; a helper's assertions fall outside that. A script confirms every cited name exists (482 checked), and every `file:line` lies inside its file (164 checked).
- **The re-trace** of every tenth Fulfilled row (18 rows) from R to the test: 17 stood. **`npart.func` did not:** its test built the admission by hand and read the config's value at its default. Widening to its neighbours moved **7 rows** (`maxInt`, `maxInt.func`, `maxInt.cat`, `npart.func`, `h1`/`h2`, `w1`, `w2`); all seven are Fulfilled again through `DeclaredControls_ReachTheModelAndTheRetainedTerms`.
- **Probes:** at the audit, 32 cases in `tests/regression/BassCompletionAuditProbes.cs`, 31 passing and 1 skipped as G-2, none over about 2 s, with nothing under `src/` changed. With the gaps closed, 35 cases, all passing, none skipped.
- **The gaps' tests** each compare against something outside the new path: tensor quadrature and the law of total variance (G-2), the whole analysis's own draws (G-1, G-3), bit identity with the ladder the count and hot end build (G-4), and the residuals and a covariance assembled by hand (G-6). A fresh-eyes review of the five changes then moved basis effects to evaluate only the draws selected, and added a test at conditioned coordinates (`b844f53`).
- **Not checked:** the long runs behind the measured limitations; anything gated by `BASS_QUALIFICATION` or `SEED_SWEEP`; the new options against R runs, since their references are those above. `Maths.LinAlg.Tests` and `Maths.Distributions.Tests` were not run here; the main session ran them at `e1eea5d`, 139 of 139 and 148 of 148, and the gaps touched neither.
