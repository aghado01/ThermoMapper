# BASS completion report

Date: 2026-11-01. Written by the audit in [bass/completion-audit-brief.md](bass/completion-audit-brief.md), against the contract in [the hub](bass-full-port-brief.md#scope-and-completion-contract). The row-by-row detail, with each test and reference, is in [bass/completion-audit-matrix.md](bass/completion-audit-matrix.md).

## Is the contract fulfilled?

**Not yet in full. Every capability P01–P15 exists and its central behaviour is verified, but six small items of R's surface differ from the contract, and two benchmark measures are open.** No statistical defect was found. All 215 rows of the matrix have a verdict; the probes this audit added found no wrong number.

What stands between:

1. **Three things the contract names or R exposes that are missing:** the `other` member of basis sensitivity (P05 names it), conditioning on one of several output variables (P04, `func.var`), and selecting draws for sensitivity (`mcmc.use`).
2. **Three narrower controls:** a declared ladder as a list (`temp.ladder`), the swap start (`start.temper`), and `predict.bassBasis`'s two switches, here one.
3. **The benchmark list's two open measures:** statistical error per unit time ([its brief](bass/performance-statistical-brief.md), not run) and the sensitivity moment work alone.

Five of the six are a small change or a one-line recorded decision; conditioning on one output variable is the largest, small to medium. Rendering and text formatting remain above maths; the contract lists them as such, so they are not gaps.

**My reading, apart from what I found.** On R's default paths the port has each capability, and it corrects the inventory's defects (18 rows in the matrix) with fixtures on both sides. The largest honest limit is not in this list: Friedman's posterior, where the default chains do not agree at practical lengths (below). That limit is measured and stated, and the contract does not promise otherwise.

## The audited revision and the suites

- **Audited revision:** `8f4d2da0081c7f17ecb524d5514a27f994191d02`. R: BASS 1.3.1, `26710146c623a258854162367a939041efedd091`.
- **Suites,** once at the end, one at a time, at below-normal priority with `-m:2` and `xUnit.MaxParallelThreads=4`, `BASS_QUALIFICATION` and `SEED_SWEEP` unset:
  - `Maths.Samplers.Tests`: 284 of 284, 1.3 minutes.
  - `Maths.Regression.Tests`: 1,507, of which 1,506 passed and 1 is skipped on purpose (G-2), 5.1 minutes.
  - The BASS classes of `Maths.Oracle.Tests` (`FullyQualifiedName~Maths.Oracle.Tests.Bass`): 90 of 90, 6.1 minutes.
- **Later commits touching BASS:** this audit's `348ad80` (probes) and `ea6bb46` (matrix), and another session's `cb0400b`, a flat-histogram refusal: a b = 1 rung on one level refuses K instead of overflowing it (`BassFlatHistogramRun.cs`, `Wham.cs`, tests). It touches no row's verdict: the flat histogram is optional and beyond the contract. The other commits since `8f4d2da` are SPC, graphs and persistent homology.
- **The tree they ran on:** the main checkout, not a clean snapshot, so that the R oracle's library was present. It began at `121005e` with `cb0400b`'s five files as another session's uncommitted edits, which it committed during the run; it ended clean in those paths at `8614a95`. So these counts are the audited revision plus the probes and that refusal. Under `src/maths/` they differ from the audited revision by `cb0400b` alone.

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
| P01 | `bass` | Fulfilled, limited | `temp.ladder` a geometric or self-sizing ladder; `start.temper` absent; tempering and log-posterior records aggregated |
| P02 | `bassPCA` | Fulfilled | |
| P03 | `bassBasis` | Fulfilled | |
| P04 | `sobol` | Fulfilled, limited | No draw selection; several output variables condition together |
| P05 | `sobolBasis` | Fulfilled, limited | No `other` member; no draw selection |
| P06 | `calibrate.bassBasis` | Fulfilled | Cold rung's draws only; the adaptation scale is on the block |
| P07 | `predict.bass` | Fulfilled | |
| P08 | `predict.bassBasis` | Fulfilled, limited | One switch for noise and truncation together |
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

| # | Gap | Where | Size |
|---|---|---|---|
| G-1 | **`other`, the residual of basis sensitivity.** P05 names it. `Variance` less the reported mains is it, but no member holds it | Matrix F11; probe `BasisOther_IsTheVarianceLessTheReportedMains` pins the arithmetic | One member and a summary field |
| G-2 | **Conditioning on one of several output variables** (`func.var`). A coordinate names every output variable; R conditions on one and integrates the others. Whole for one output variable | Matrix E06; probe `ConditionOnOneOfTwoOutputVariables_IntegratesTheOther`, `[Fact(Skip)]`. Reproduction: `Coordinates = [[0.5]]` on a two-variable fit refuses, "Coordinate 0 has 1 values for 2 output variables" (`sensitivity/BassSensitivity.cs:391`) | Small to medium: the measure already integrates output variables; the conditioned form takes a partial coordinate |
| G-3 | **`mcmc.use` on `sobol` and `sobolBasis`.** Every retained draw is decomposed. Reachable only as `fit with { Draws = … }` | Matrix E05, F04 | Small: one field, as `BassPredictOptions.Draws` |
| G-4 | **`temp.ladder` as a list.** A geometric ladder of declared count and hot end, adaptive, frozen or self-sizing; not an arbitrary list. Calibration's `Tempering` takes one | Matrix A22; `BassLadder.cs:56` already takes a list, internally | Small: one field and its validation, or a recorded decision |
| G-5 | **`start.temper`.** Swaps start at the first scan. No record names the argument; `stage4-tempering-brief.md:206` left the cadence as a reversible default. Burn-in behaviour only, with no change of target | Matrix A23 | A recorded decision, or one field |
| G-6 | **`predict.bassBasis`'s `nugget` and `trunc.error`, one switch here.** The members to form R's noise-only and truncation-only draws are returned; the draws are not. D07 fixes the full draw and no record names the two-switch form | Matrix D-p09 | A recorded decision, or one option |
| G-7 | **Statistical error per unit time.** Known | Matrix part 5; [its brief](bass/performance-statistical-brief.md) | Long runs on the user's schedule; `BASS_QUALIFICATION` left unset here |
| G-8 | **The sensitivity moment work alone.** Measured only inside Sobol and effects per draw. Known | Matrix part 5; `performance-cost-brief.md:380` | A few timed points in the cost harness |

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

Each is a few lines for the main session; the matrix's part 6 has the lists.

- **The source map** omits `BassObservationModel.cs`, `DesignKernel.cs`, `DesignKernels.cs` (doc line 2138). The **extension seams' write-up** in `docs/Maths` is owed, listed once.
- **The API rows** name no removed thing, but 71 of 242 public types appear in no row: 10 seams, 31 a caller reads, 30 readouts and building blocks. `BassComponents.Assemble`, which carries basis restart, is unnamed.
- **Statements found false:** `port-fitting-brief.md:503` (a cached training fit: none is built); `opus-bass-completion-remainder.md:1` ("Everything R's BASS can do is now built"); `BassRunOptions.cs:15` (tempering and SMC "when they join"); `completion-audit-brief.md` ("E1–E15": the checklist holds E18); `port-outputs-brief.md:82` ("condition on one functional coordinate", which G-2 qualifies).
- **The open boxes** N15, N17, F-M2 and F-X4 are not stale: each is open by its own condition (the matrix, 6a).

## What the brief does not rule on

Three questions would change the conclusion above. I have decided none of them.

1. **Are G-4, G-5 and G-6 gaps or recorded choices?** If you rule that a declared ladder list, `start.temper` and the two switches stay out, each becomes "left out by decision", and the list shortens to G-1, G-2, G-3 and the two measures.
2. **Is G-3 in the contract?** The hub's P04 row does not name draw selection; the inventory's "Sensitivity" paragraph does. I took the inventory as the checklist, as the brief says.
3. **Do the per-iteration records count?** R's `log.post.cold`, `swap` and `temp.val` are aggregates here (matrix B-r10, B-r11), and calibration keeps the cold rung only (G14). I read the contract's "posterior retention" as not asking for them and marked them limited, not gaps.

## How this was checked

- **The standard:** each Fulfilled row names a test read for this audit. About 80 cited tests were read whole. For the other 192, a script kept each one's first nine assertions and comments, and those were read; a helper's assertions fall outside that. A script confirms every cited name exists (482 checked), and every `file:line` lies inside its file (164 checked).
- **The re-trace** of every tenth Fulfilled row (18 rows) from R to the test: 17 stood. **`npart.func` did not:** its test built the admission by hand and read the config's value at its default. Widening to its neighbours moved **7 rows** (`maxInt`, `maxInt.func`, `maxInt.cat`, `npart.func`, `h1`/`h2`, `w1`, `w2`); all seven are Fulfilled again through `DeclaredControls_ReachTheModelAndTheRetainedTerms`.
- **Probes:** 32 cases in `tests/regression/BassCompletionAuditProbes.cs`, 31 passing and 1 skipped as G-2, none over about 2 s. Nothing under `src/` changed.
- **Not checked:** the long runs behind the measured limitations; `Maths.LinAlg.Tests` and `Maths.Distributions.Tests` at this revision; anything gated by `BASS_QUALIFICATION` or `SEED_SWEEP`.
