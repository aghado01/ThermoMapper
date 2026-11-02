# The BASS completion audit: matrix

Audited revision: `8f4d2da0081c7f17ecb524d5514a27f994191d02` (`git rev-parse HEAD` at the start). R: BASS 1.3.1 at `D:/aghado01/packages/BASS`, commit `26710146c623a258854162367a939041efedd091`. The brief is [completion-audit-brief.md](completion-audit-brief.md); the contract is the hub's [Scope and completion contract](../bass-full-port-brief.md#scope-and-completion-contract). The report is [bass-completion-report.md](../bass-completion-report.md).

**Verdicts** (the brief's): **F** Fulfilled; **FL** Fulfilled, limited; **U** Built, untested; **G** Gap; **D** Left out by decision. A verdict of F names a test read for this audit and its reference. Paths are under `src/maths/regression/spline/bass/` for code, and tests under `tests/regression/` unless they say `oracle/`. An `oracle/` test returns early where R is absent, so its verdict rests on a run with R present (the suites in the report).

**Probes** are in `tests/regression/BassCompletionAuditProbes.cs` (the class `BassCompletionAuditProbes`, written "Probe." below).

**The suites at the end,** in the main checkout, below normal, `-m:2`, `xUnit.MaxParallelThreads=4`, `BASS_QUALIFICATION` and `SEED_SWEEP` unset: Samplers 284 of 284; Regression 1,507, with 1,506 passed and the G-2 probe skipped; the BASS oracle classes 90 of 90, R present (`Basis_MatchesRBass` took 4 s per degree). **Later commits touching BASS:** this audit's own, and `cb0400b` (a flat-histogram refusal, another session's), uncommitted when the suites began and committed during them. It changes no verdict below.

**Closed since the audit, 2026-11-02, on the user's rulings.** G-1 to G-4 and G-6 are built, and G-5 is left out by decision; the report's gap list has each commit. Rows moved: A22, D-p09, E05, E06, F04 and F11 to Fulfilled, A23 to Left out by decision, and the roll-up's P04, P05 and P08 to Fulfilled. The notes on A14, B-r10, B-r11 and G14 record the rulings. **Verdicts now: F 192, FL 5, D 18, G 0** (at the audit F 183, FL 12, D 17, G 3). Those rows' lines are at `b844f53`; every other row's are still the audited revision's. **The suites at `b844f53`,** run the same way: `Maths.Samplers.Tests` 284 of 284; `Maths.Regression.Tests` 1,544 of 1,544, none skipped; the BASS classes of `Maths.Oracle.Tests` 90 of 90, R present.

## Part 1: capabilities

### 1a. `bass`: arguments

| # | R argument | C# | Test read | Reference | Verdict | Note |
|---|---|---|---|---|---|---|
| A01 | `xx` | `BassData.cs:18`, `ContinuousPredictor` :72, `CategoricalPredictor` :104 | `BassContractTests.Data_RefusesItsDefects_WithReasons`, `Data_CopiesWhatItIsGiven`; `oracle/BassBasisParityTests.Basis_MatchesRBass` (scaled columns to 10⁻¹⁵) | R parity; refusals by reason | F | Typed schema, by name. R guesses factors and orientation; here neither |
| A02 | `y` scalar or matrix | `ScalarResponse` `BassData.cs:169`, `FunctionalResponse` :202 | `BassContractTests.Data_*`; `BassRunTests.Fit_RunsEveryOption("functional")` | Refusals by reason; structures in F | F | Non-finite refused, as R refuses missing |
| A03 | `maxInt` | `BassConfig.cs:260`, resolved `BassResolved.cs:156` | `BassContractTests.Resolved_MatchesDirectComputation` (cap cut to variables, logged); `BassProposalTests.Prior_ReadsTheConfigsCapsAndOrders`; probe `DeclaredControls_ReachTheModelAndTheRetainedTerms` | Direct computation, enumeration; retained terms | F | Moved by the re-trace |
| A04 | `maxInt.func` | `BassConfig.cs:266`, `BassResolved.cs:158` | `BassContractTests.Resolved_Functional_CountsBothParts`; `BassProposalTests.Prior_ReadsTheConfigsCapsAndOrders`; probe `DeclaredControls_ReachTheModelAndTheRetainedTerms` | Direct computation; retained terms | F | Moved by the re-trace |
| A05 | `maxInt.cat` | `BassConfig.cs:263`, `BassResolved.cs:157` | as A03 | Direct computation; retained terms | F | Moved by the re-trace |
| A06 | `xx.func` | `FunctionalResponse.OutputVariables` `BassData.cs:244` | `oracle/BassBasisParityTests.Basis_MatchesRBass` (two output variables); `BassBasisTests.FunctionalDesign_SeparableForm_MatchesTheExplicitDesign`; probe `TwoOutputVariables_PredictAndIntegrate` | R parity, explicit tensor, hand evaluation | F | Locations need not form a grid; the probe's do not |
| A07 | `degree` | `BassConfig.cs:257` | `oracle/BassBasisParityTests.Basis_MatchesRBass` at 1, 2, 3; `BassBasisTests.CaseFactors_MatchRsConstruction_ByHand`; `oracle/BassSensitivityDefectParityTests.EffectCurves_UseTheDeclaredDegree`; probe `HigherDegree_PredictsDecomposesAndHasEffects` | R parity; hand values; a grid | F | R's validation is commented out; `Validate` refuses < 1. Before the probe, degrees above 1 had no prediction test at new inputs on a fit |
| A08 | `maxBasis` | `BassConfig.cs:269`, `BassResolved.cs:160` | `BassContractTests.Resolved_MatchesDirectComputation` | Direct computation | F | min(1000, N − 1), the hub's corrected bound; R allows one column too many |
| A09 | `npart` | `SupportRule.Minimum` `BassConfig.cs:331`, `BassResolved.cs:174` | `BassContractTests.Resolved_MatchesDirectComputation`; `BassDesignTests.EachRefusalReason_HasItsCase`; `BassBasisTests.Support_CountsEachPart` | Direct computation | F | |
| A10 | `npart.func` | `SupportRule.MinimumFunctional`, `BassResolved.cs:175` | `BassDesignTests.EachRefusalReason_HasItsCase` (its own threshold, the admission built by hand); `BassContractTests.Resolved_Functional_CountsBothParts` (the default); probe `DeclaredControls_ReachTheModelAndTheRetainedTerms` (a declared value through a fit) | Direct computation; retained terms | F | Moved by the re-trace: was untested as a declared control |
| A11 | `nmcmc`, `nburn` | `BassSampling.BurnIn` `BassConfig.cs:392`, `Samples` :395 | `BassRunTests.Fit_RebuildsFromTheSeedTreeByHand`; `BassContractTests.Config_HasTheDecidedDefaults` | By-hand rebuild | F | Kept count is explicit; no fractional quotient |
| A12 | `thin` | `BassSampling.Thin` `BassConfig.cs:398` | `BassRunTests.Thinning_KeepsTheTrajectory`, `Fit_RebuildsFromTheSeedTreeByHand(2)` | Same trajectory | F | |
| A13 | `g1`, `g2` | `Sigma2Prior.InverseGamma`, `Jeffreys`; `BassConfig.cs:36` | `BassTargetTests.DeclaredInverseGamma_ReachesTheTarget`; `BassRunTests.Fit_RunsEveryOption("inverse gamma", "jeffreys")`; `oracle/BassChainParityTests("A", "jeffreys")` | Closed form, enumeration | F | R's default (0, 0) is `Jeffreys`; the default here is D03's scaled prior |
| A14 | `s2.lower` | none | none | | D | Not carried: contract shapes, `port-fitting-brief.md:487` (a truncated σ² prior, which would have to be stated as one). The user's ruling, 2026-11-01: nothing to build; a truncated σ² prior joins F-R1's rerun as a pre-registered arm (`deterministic-profile-brief.md:153`) |
| A15 | `h1`, `h2` | `BassCountPrior` `BassConfig.cs:334`, `BassStructurePrior.cs:178,189` | `BassProposalTests.CountPrior_IsTheNegativeBinomial`; `BassSweepTests.ConstantLikelihood_SamplesThePriorOnF`; probe `DeclaredControls_ReachTheModelAndTheRetainedTerms` | Closed form, enumeration | F | λ integrated out, D02. Moved by the re-trace: the non-default path ran only in the gated F-R1 harness |
| A16 | `a.tau` | `TauPrior.UnitInformationGamma(shape)` | `BassTargetTests.TauPrior_IsTheUnitInformationGamma`; `BassCoefficientPriorTests.TauPosterior_OnAFixedStructure_MatchesQuadrature` | Quadrature, closed form | F | |
| A17 | `b.tau` | rate a·n from the family; no free rate | `BassContractTests.Resolved_MatchesDirectComputation` (0.5, N/2) | | D | D03, **Corrected**: R's 2/n contradicts its own documentation; hub register, doc line 2074 "not R's rate" |
| A18 | `w1` | `OrderWeight` `BassConfig.cs:284` | `BassProposalTests.BirthProbability_ByHand`, `BirthProposal_SumsToOne_ForEveryCountState`; probe `DeclaredControls_ReachTheModelAndTheRetainedTerms` | By hand, enumeration | F | Validated positive, D01's reason. Moved by the re-trace |
| A19 | `w2` | `VariableWeight` `BassConfig.cs:287` | as A18; `oracle/BassTargetParityTests.Marginals_MatchTheOraclesPosterior` | By hand; `dmwnchBass` at R's weights | F | Moved by the re-trace |
| A20 | `beta.prior = 'g'` | `BassCoefficientPrior.G()` `BassConfig.cs:33` | `oracle/BassFitParityTests.Fit_MatchesTheOraclesPosterior`; `oracle/BassChainParityTests` | Enumeration by the oracle | F | Prior corrected (D03); level-centred (D03) |
| A21 | `beta.prior = 'jeffreys'` | `G(TauPrior.Held())` keeps its practical intent; R's target not ported | `oracle/BassTargetParityTests.ExactFits_AreTheOraclesExclusion`; `BassRunTests.Fit_RunsEveryOption("fixed g")`; `BassPriorSamplerTests.FixedG_HoldsTau` | Enumeration | D | D03: no single target in R (F-O6); `port-fitting-brief.md:487` |
| A22 | `temp.ladder` | `TemperingSpec.InverseTemperatures` `BassConfig.cs:523`, taken at `BassLadder.cs:269`; or `Levels`, `HotEnd`, `Adaptive`, `MaxLevels` | `BassTemperingTests.List_OfTheBuiltRungs_RunsBitIdenticalToTheCountAndHotEnd`, `List_RoundTrips_MovesTheFingerprintWithARung_AndReadsAsItsCountAndHotEnd`, `Fit_ReportsTheLaddersReadout` ("declared list", non-adaptive, equals the declared rungs); `BassTemperingResumeTests.Resume_IsBitIdentical`; `oracle/BassTemperingParityTests.Ladder_EachRungMatchesTheOraclesTemperedLaw` | Bit identity with the ladder the count and hot end build; oracle rung laws | F | Closed by G-4, `3a83d0a`. R's list runs as given under `Adaptive = false`. On the default `true` it is a start: the interior rungs re-space, the ends and the count hold (`List_IsCopiedPerChain_SoAnAdaptiveRunLeavesTheDeclarationAlone`) |
| A23 | `start.temper` | none: swaps from the first scan | none | | D | Left out by the user's ruling, 2026-11-01 (`port-fitting-brief.md:487`, "Not carried"): it changes no target, and the self-sizing ladder reads its swap rejections during burn-in |
| A24 | `curr.list` | `BassFit.Checkpoint`, `Bass.Resume` `Bass.cs:38` | `BassRunResumeTests.Resume_IsBitIdentical` (10 options), `Resume_RefusesAnotherRun`; `BassTemperingResumeTests.Resume_IsBitIdentical` | Bit identity with the uninterrupted run | F | Stronger than R: RNG, adaptation, replica map and a data/config fingerprint |
| A25 | `save.yhat` | `BassFit.Predict`, `TrainingFit`: views | `BassPredictTests.AtTheTrainingInputs_ReproducesTheStageFourFits`; `BassPresentationTests.TrainingFit_IsThePredictionAtTheTrainingInputs` | Direct evaluation | D | Contract shapes, `port-fitting-brief.md:487,503`: a view. No cached view is built, nor needed for the contract |
| A26 | `small` | `BassRunOptions.Checkpoint`; `DrawStore` thinning | `BassRunTests.Fit_DoesNotDependOnTheRunsOtherOptions`; `BassPresentationTests.Traces_AreTheDrawStores_WithTheirSweeps` | Same fit | D | Contract shapes, `port-fitting-brief.md:487` |
| A27 | `verbose` | `BassRunOptions.Progress` | `BassRunTests.Progress_CountsEverySweep`; `BassRunResumeTests.Resume_ReportsItsOwnSweeps` | Exact counts | F | |
| A28 | `ret.str` | `BassFit.Config`, `Resolved`, `Transform`; `BassModel.Create` | `BassPresentationTests.Metadata_RoundTripsWhatTheFitHolds`; `BassContractTests.Config_RoundTripsThroughJson_ByValue` | Identity | F | Values, not live closures, as the inventory asks |

### 1b. `bass`: returned elements

| # | R element | C# | Test read | Reference | Verdict | Note |
|---|---|---|---|---|---|---|
| B-r01 | `yhat`, `yhat.mean` | `BassFit.Predict` → `BassPrediction.Latent`, `Mean()` | `BassPredictTests.AtTheTrainingInputs_ReproducesTheStageFourFits` | `BassDrawFit` evaluation | F | On demand |
| B-r02 | `call` | `BassConfig` as JSON, `RunIdentity` | `BassContractTests.Config_RoundTripsThroughJson_ByValue`, `Fingerprint_MovesWithAnySetting_AndAnyDataBit` | Round trip | F | |
| B-r03 | `beta` | `BassDraw.Coefficients` `BassState.cs:56` | `BassSweepTests.FixedStructure_TauAndComposition_MatchQuadratureAndAJointGibbsReference`, `Compose_AlignsIdsAndUncentredCoefficients` | Joint Gibbs, quadrature | F | Uncentred basis, R's |
| B-r04 | `s2` | `BassDraw.Sigma2` | as B-r03; probe `InterceptOnly_PosteriorIsItsClosedForm` | Closed form | F | |
| B-r05 | `lam` | `BassDraw.Lambda` | `BassSweepTests.ConstantLikelihood_SamplesThePriorOnF` (E[λ] against the Gamma mixture) | Enumeration | F | Drawn at retention, D02 |
| B-r06 | `nbasis` | `BassDraw.Structure.TermIds.Length`, `BassTraces.Count` | `BassPresentationTests.Traces_AreTheDrawStores_WithTheirSweeps` | The stores | F | |
| B-r07 | `beta.prec` | `BassDraw.Tau` | `BassCoefficientPriorTests.TauPosterior_OnAFixedStructure_MatchesQuadrature` | Quadrature | F | Corrected prior |
| B-r08 | `degree`, `nmcmc`, `nburn`, `thin`, `p` | `BassMetadata.Config`, `Lengths`, `Data` `presentation/BassSummary.cs:42` | `BassPresentationTests.Metadata_RoundTripsWhatTheFitHolds` | Identity | F | |
| B-r09 | `y` | not held: `BassData` rides the view, checked by fingerprint `BassFit.cs:172` | `BassFitViewsTests.PointwiseLogDensity_ReadsTheFitsOwnData`; `BassPresentationTests.TrainingFit_RefusesOtherData_AndProbabilitiesOutsideTheOpenInterval` | Refusal by fingerprint | D | Contract shapes; reversible default in the doc's status section |
| B-r10 | `log.post.cold` | `BassDiagnostics.LogPosterior` | `BassVerdictTests.LogPosterior_IsTheChainsOwnDensity` (bit identity with the chain's density) | Chain's own density | FL | A convergence reading, not a per-draw trace. Not in `BassChainTrace`. The user's ruling, 2026-11-01: stays limited, not a gap (`port-fitting-brief.md:487`) |
| B-r11 | `swap`, `count.swap`, `count.swap.prop`, `temp.val`, `temp.ladder` | `BassTemperingReadout` `BassDiagnostics.cs:131`: `PairAcceptance`, `PairRejection`, `Ladders`, round trips | `BassTemperingTests.Fit_ReportsTheLaddersReadout`; `BassTemperingResumeTests` | Own consistency; oracle swaps | FL | Per-pair aggregates over the sampling phase. R's per-iteration `swap` and `temp.val` traces are not kept. The user's ruling, 2026-11-01: stays limited, not a gap (`port-fitting-brief.md:487`) |
| B-r12 | `n.models`, `model.lookup` | `BassFit.Terms`, `BassDraw.Structure` | `BassRunTests.Fit_RebuildsFromTheSeedTreeByHand` | By hand | F | Interning by term, not by model |
| B-r13 | `des`, `func`, `cat`, `type`, `cx` | Predictor and response types; `BassDataShape` | `BassPresentationTests.Metadata_RoundTripsWhatTheFitHolds` | Identity | F | |
| B-r14 | `curr.list` | `BassFit.Checkpoint` | A24 | | F | Beside the posterior |
| B-r15 | `knotInd.des`, `signs.des`, `vars.des`, `n.int.des` | `BassTerm.Continuous`, `HingeFactor` `BassTerm.cs:11` | `BassBasisTests.Terms_AreCanonical_AndCompareByValue`; `oracle/BassBasisParityTests` | R parity | F | |
| B-r16 | `maxInt.des`, `maxInt.cat`, `maxInt.func` | `BassResolved` `BassResolved.cs:37-43` | `BassContractTests.Resolved_*` | Direct | F | |
| B-r17 | `pdes`, `pcat`, `pfunc`, `nlevels` | `BassDataShape`, `TrainingLevels` | `BassPresentationTests.Metadata_RoundTripsWhatTheFitHolds` | Identity | F | |
| B-r18 | `xx.des`, `xx.cat`, `xx.func`, `des.basis`, `cat.basis`, `func.basis` | `BassInputs` from the data; not held in the fit | `BassRunTests.Fit_RunsEveryOption`; `BassBasisTests.Inputs_ScaleColumns_AndTakeDistinctScaledValues_AsKnots` | | D | Contract shapes: the fit holds neither data nor working bases |
| B-r19 | `range.des`, `range.func` | `BassInputTransform`, `RangeScale` | `BassContractTests.Transform_MapsTheTrainingRange_OntoTheUnitInterval_WithExactEnds` | Direct | F | |
| B-r20 | `unique.ind.des`, `unique.ind.func` | `BassInputs.ContinuousKnots`, `FunctionalKnots` `BassInputs.cs:105,108` | `BassBasisTests.Inputs_ScaleColumns_AndTakeDistinctScaledValues_AsKnots` | Direct | F | |
| B-r21 | `vars.cat`, `sub.size`, `sub.list`, `n.int.cat` | `BassTerm.Categorical`, `CategoricalFactor` `BassTerm.cs:17` | `BassBasisTests.CategoricalFactors_AreIndicators_AndAnUnseenLevelIsOutsideEverySubset` | Direct | F | |
| B-r22 | `knotInd.func`, `signs.func`, `vars.func`, `n.int.func` | `BassTerm.Functional` | `oracle/BassBasisParityTests` | R parity | F | |
| B-r23 | `data`, `prior`, `funcs` (`ret.str`) | A28 | A28 | | F | |

### 1c. `bassPCA` and `bassBasis`

| # | R item | C# | Test read | Reference | Verdict | Note |
|---|---|---|---|---|---|---|
| C01 | `bassPCA(xx)` | `BassComponents.Fit(.., predictors, ..)` `output/BassComponents.cs:129` | `BassOutputComponentsTests.EachComponent_IsFittedToItsCoefficientResponse_UnderTheSeedTree` | By hand | F | |
| C02 | `bassPCA(y)` | `BassDecomposition.FromPca` `output/BassDecomposition.cs:354` | `BassOutputDecompositionTests.Pca_ReconstructsTheTrainingResponses_InRsConvention` | Reconstruction, sd by hand | F | |
| C03 | `bassPCA(dat)`, `bassBasis(dat)` | `FromBasis` `output/BassDecomposition.cs:453` | `BassOutputDecompositionTests.SuppliedBasis_Nonorthogonal_ReconstructsWithItsResiduals`; probe `Cell_Basis_ReconstructsDecomposesAndRestarts(*, Supplied)` | Hand reconstruction | F | Nonorthogonal |
| C04 | `n.pc` | `BassPcaSpec.Rank` :260 | `BassOutputDecompositionTests.RankRule_TakesTheSmallestRankReachingTheShare` | Direct | F | Capped at ρ, flagged |
| C05 | `perc.var` | `VarianceShare` :257 | as C04; `FullRetention_IsTheNumericalRank_WithZeroInSampleTruncation`; `oracle/BassOutputDefectsParityTests.FullRetention_IsTheNumericalRank` | Direct; R's NaN at 100 | F | R's failure at 100 corrected |
| C06 | `center` | `Centre` :263 | `BassOutputDecompositionTests.LeaveOneOut_UnderAFixedCentre_EqualsExplicitRefits`, `Pca_*` | Direct | F | A vector; zeros for R's `center=F` |
| C07 | `scale` | `ScaleByColumnSd` :266, `Scale` :269 | `Pca_*`; `ConstantColumn_IsCentred_AndScaledByOne`; `Refusals_SayWhy` | Direct | F | |
| C08 | `n.cores` | `BassRunOptions.MaxDegreeOfParallelism` | `BassOutputComponentsTests.Draws_DoNotDependOnParallelism` | Bit identity | F | |
| C09 | `parType` | none | | | D | R's fork and socket are process plumbing; the hub's "no R object system" |
| C10 | `...` | `BassConfig` per component `output/BassComponents.cs:141` | `BassOutputComponentsTests.EachComponent_*`, `UnequalCounts_AreThinnedEvenlyToTheShortest` | By hand | F | |
| C11 | setup's `basis`, `newy` | `Map.Basis`, `Coefficient`, `CoefficientResponse` | `Pca_*` | Identity | F | |
| C12 | `y.m`, `y.s` | `Map.Offsets`, `Scales` | `Pca_*` | By hand | F | |
| C13 | `trunc.error` | `Truncation` (`Values`, `InSample`, `LeaveOneOut`) `output/BassDecomposition.cs:152` | `BassOutputDecompositionTests.Residual_IsObservedMinusReconstructed_InStandardizedUnits`, `LeaveOneOut_EqualsExplicitRefits`; `oracle/BassOutputTruncationParityTests.Residuals_MatchTheOracle` | R parity where right, refits | F | Orientation corrected; D07's leave-one-out default |
| C14 | `ev` | `BassPcaAccount.SingularValues` :280 | `Pca_*` | Direct | F | |
| C15 | `y`, `xx` in `dat` | `BassBasisFit.Predictors`; `BassTrainingFit.Observed` rebuilt as K w + in-sample e `presentation/BassTrainingFit.cs:71` | `BassPresentationBasisTests.TrainingFit_IsThePredictionAtTheTrainingCases_InPhysicalUnits` (`Observed` equals y to 12 places) | The training y | F | The responses are not held; the decomposition represents them exactly (D07) |
| C16 | `bassBasis` → `mod.list` | `BassBasisFit.Components` | `BassOutputComponentsTests.EachComponent_*` | By hand | F | |
| C17 | `bassBasis` → `dat` | `BassBasisFit.Decomposition` | `BassPresentationBasisTests.Metadata_RoundTripsTheDecompositionTheJointDrawsAndTheComponents` | Identity | F | |

### 1d. `predict.bass` and `predict.bassBasis`

| # | R item | C# | Test read | Reference | Verdict | Note |
|---|---|---|---|---|---|---|
| D-p01 | `object`, `newdata` | `BassFit.Predict(data, PredictionCases, ..)` `BassFit.cs:81` | `BassPredictTests.NewCases_EvaluateTheHinges_AndFlagExtrapolation`; probe `Cell_Direct_*` | Hand evaluation | F | By name; extrapolation flagged |
| D-p02 | `newdata.func` | `BassPredictOptions.OutputLocations` `posterior/BassPrediction.cs:30` | `BassPredictTests.NewOutputLocations_ForADirectFunctionalFit`; probe `Cell_Direct_*(*, Direct)` | Hand evaluation | F | |
| D-p03 | `mcmc.use` | `BassPredictOptions.Draws` :17 | `BassPredictTests.DrawOrder_AndSingletonAxes_AreKept`; `oracle/BassOutputDefectsParityTests.DrawOrder_AxesAndNoise_FollowTheRequest` | R's reordering a fixture | F | Caller's order kept; R's corrected |
| D-p04 | `verbose` | none | | | D | No long loop to report; the hub's "formatting above the engine" |
| D-p05 | `nugget` | `Predictive` :20 | `BassPredictTests.Noise_IsIndependentAcrossDrawsCasesAndLocations` | Statistical; R's recycling a fixture | F | Default true, with `Latent` always kept; R's FALSE returns latent |
| D-p06 | return, scalar and functional | `BassPrediction` `posterior/BassPrediction.cs:50`: `Latent`, `Predictive`, draw × case × location | `DrawOrder_AndSingletonAxes_AreKept` | Fixed axes | F | R's `drop` corrected |
| D-p07 | `predict.bassBasis(newdata)` | `BassBasisFit.Predict` `output/BassComponents.cs:63`, `posterior/BassBasisPrediction.cs:147` | `BassPredictBasisTests.Latent_ReconstructsEachComponentsOwnPrediction`; probe `Cell_Basis_*` | Each component's own and hand | F | No new output locations: refused, as the hub allows |
| D-p08 | `mcmc.use` (basis) | `Draws` over `JointDraws` | `BassPredictBasisTests.DrawOrder_SingletonAxes_AndTheSeed` | Identity | F | |
| D-p09 | `trunc.error`, `nugget` (basis) | `BassPredictOptions.BasisParts` `posterior/BassPrediction.cs:40`, `BasisPredictiveParts` :12: both, D07's default; `Noise`, R's defaults; `Truncation`. `Latent` is neither | `BassPredictBasisTests.Parts_TakeEachSourceAlone_AndTheDefaultIsTheFullDrawBitForBit`, `Predictive_DrawsComponentNoisePerCase_AndWholeResidualsOverEveryCase`; `oracle/BassOutputDefectsParityTests.TruncationSampler_DrawsWholeResidualsOverEveryCase`, `ComponentNoise_IsDrawnPerCase` | Recovered from the draws; the covariance by hand; R fixtures | F | Closed by G-6, `26befff`. The default stays the full draw, bit for bit. A direct fit, with no truncation residuals, refuses another value |
| D-p10 | `n.cores`, `parType` (basis) | none: sequential | | | D | Not a result; R's own doc says 1 is often fastest |
| D-p11 | return (basis) | `BassBasisPrediction` :24: draw × case × location | `DrawOrder_SingletonAxes_AndTheSeed`, `OneLocation_KeepsItsAxis` | Fixed axes | F | |

### 1e. `sobol`

| # | R item | C# | Test read | Reference | Verdict | Note |
|---|---|---|---|---|---|---|
| E01 | `bassMod` | `BassFit.Sensitivity`, `BassSensitivity.Of` | `BassSensitivityTests.Indices_MatchPickFreezeMonteCarloOnAFittedDraw` | Pick-freeze Monte Carlo | F | |
| E02 | `prior` continuous | `SensitivitySpec.Measures`, `UniformMeasure`, `MixtureMeasure`, `EmpiricalMeasure` `sensitivity/InputMeasures.cs:20,27,75` | `BassSensitivityMomentTests.TermMoments_MatchTensorQuadratureOfTheBasisColumns`; `oracle/BassMomentParityTests.HingeMoments_MatchQuadrature_UnderEveryMeasure`; probe `TruncatedMeasure_SensitivityMatchesQuadrature` | Quadrature, R moments | F | Weights before or after truncation are declared |
| E03 | `prior` categorical weights | `CategoricalMeasure` :81 | `oracle/BassSensitivityDefectParityTests.CategoryWeights_AreRead`; `oracle/BassMomentParityTests.CategoricalMoments_MatchEnumeration`; probe `ZeroProbabilityLevel_SensitivityMatchesEnumeration` | Enumeration | F | R ignores them: corrected |
| E04 | `prior.func` | `OutputMeasures` `sensitivity/BassSensitivity.cs:54` | `BassSensitivityTests.FunctionalForms_ObeyTheLawOfTotalVariance` (empirical grid); probe `Cell_Direct_*(*, Direct)` (uniform) | Enumeration, grid | F | |
| E05 | `mcmc.use` | `SensitivitySpec.Draws` `sensitivity/SobolDecomposition.cs:60`, `sensitivity/BassSensitivity.cs:384` | `BassSensitivityTests.Draws_SelectAndOrderTheDecomposedDraws`, `Draws_KeepChainOrder_AndTheirWeights` | Each selected draw against the whole analysis's | F | Closed by G-3, `1364cbb`, `b844f53`. Order kept and repeats allowed, as `BassPredictOptions.Draws`. Effects and conditioned coordinates follow; a weighted fit's selected draws keep their weights |
| E06 | `func.var` | `Functional = Condition` `sensitivity/SobolDecomposition.cs:74`, `ConditionOn` :79 | `BassSensitivityTests.ConditioningOnSome_ObeysTheLawOfTotalVariance`, `ConditioningOnEvery_IsTheConditionedForm_InAnyOrder`, `ConditionedVariance_MatchesTensorQuadratureAtEachCoordinate`, `FunctionalForms_ObeyTheLawOfTotalVariance`; probes `ConditionOnOneOfTwoOutputVariables_IntegratesTheOther`, `ConditionOnOne_HandBuiltTerms_CarryTheNormalizerOnce` | Tensor quadrature; the law of total variance | F | Closed by G-2, `4c34948`. `ConditionOn` holds the variables it names and integrates the others; R's one variable is one name. The output side's normalizer stays with the coordinate |
| E07 | `xx.func.var` | `Coordinates` :72 | `BassPresentationSensitivityTests.ConditionedSummary_IsOnePerCoordinate` | Identity; outside flagged | F | A coordinate holds the variables `ConditionOn` names, or every output variable; the default takes their distinct training values |
| E08 | `verbose` | none | | | D | As D-p04 |
| E09 | `getEffects` | `BassFit.Effects`, `BassSensitivity.Effects` `BassFit.cs:103` | `BassSensitivityEffectTests.Effects_AreConditionalExpectationsLessTheirLowerOrders`; `oracle/BassSensitivityDefectParityTests.EffectMoments_UseTheDeclaredMeasure`, `EffectCurves_UseTheDeclaredDegree`; probe `Cell_Direct_*` | Conditional expectations; quadrature | F | Beyond R: categorical, functional, any measure and degree |
| E10 | return `S`, `T` | `MainIndex`, `TotalIndex`, `MainVariance`, `TotalVariance` `sensitivity/SobolDecomposition.cs:109-139` | `BassSensitivityTests.SobolIdentities_HoldOnEveryDraw`, `Indices_MatchPickFreeze*` | Identities, Monte Carlo | F | Totals on the functional path too; R omits them |
| E11 | return `var.tot` | `Variance` :115 | as E10; probes `Cell_*` | Grid | F | |
| E12 | return `f0` | `EffectFunctions.Mean` | `BassSensitivityEffectTests.*`; probe `Cell_Direct_*` | Conditional means | F | |
| E13 | return `ints` | `BassMoments` `sensitivity/BassMoments.cs` | `BassSensitivityMomentTests.PairStructure_SumsToTheCovariance_AndItsTotalsToTheirComponents` | Quadrature | F | A service, not a result field |
| E14 | return `prior` | `SobolAnalysis.Measure` | `BassSensitivityMomentTests.Resolution_RefusesUnknownNamesAndWrongKinds` | Direct | F | |
| E15 | return `effects` | `EffectAnalysis` `sensitivity/EffectFunctions.cs` | E09 | | F | |
| E16 | return `names.ind`, `func` | `Inputs`, `Subsets` :102,106; `SobolAnalysis.Form` | `BassSensitivityTests.SobolIdentities_*` | Identity | F | |
| E17 | functional return `S`, `S.var`, `xx` | `SobolAnalysis.Decompositions`, `Coordinates` `SobolDecomposition.cs:176` | `BassSensitivityTests.ConditionedVariance_*`; `BassPresentationSensitivityTests.ConditionedSummary_*` | Quadrature | F | |

### 1f. `sobolBasis`

| # | R item | C# | Test read | Reference | Verdict | Note |
|---|---|---|---|---|---|---|
| F01 | `mod` | `BassBasisFit.Sensitivity`, `BassBasisSensitivity.Of` :27 | `BassSensitivityBasisTests.CrossComponentTerms_MatchTensorQuadratureOfTheReconstruction`; `OneComponent_IsTheDirectDecompositionScaled`; probe `Cell_Basis_*` | Tensor quadrature | F | |
| F02 | `int.order` | `MaxOrder` | `BassSensitivityBasisTests.OrderAndBudget_BoundTheSubsets` | The totals held | F | |
| F03 | `prior` | `Measures` | `CrossComponentTerms_*`; `oracle/BassSensitivityDefectParityTests.TruncationBounds_AreClipped`; probe `TruncatedMeasure_SensitivityMatchesQuadrature(Pca)` | Quadrature | F | R's truncation error corrected |
| F04 | `mcmc.use` | `SensitivitySpec.Draws` over the joint draws, `sensitivity/BassBasisSensitivity.cs:53`, `sensitivity/BassBasisEffects.cs:36` | `BassSensitivityBasisTests.Draws_SelectTheJointDraws`, `BassSensitivityBasisEffectTests.Draws_SelectTheJointDraws`, `BassSensitivityBasisTests.JointDraws_PairTheComponentsDraws_ThinnedOrResampled` | Each selected draw against the whole analysis's | F | Closed by G-3, `1364cbb`, `b844f53`. R's one draw is `Draws = [d]`. Effects evaluate only the component draws selected, so a selection fits a budget the whole run would exceed |
| F05 | `nind` | `SubsetBudget` refuses before computing | `OrderAndBudget_BoundTheSubsets` | | D | `port-outputs-brief.md:88`: "the port exposes a budget instead" |
| F06 | `n.cores`, `parType`, `verbose` | none | | | D | As D-p10 |
| F07 | `plot` | `BasisSobolSummary` above maths | P11 below | | D | Rendering above maths, the hub's P09–P11 row |
| F08 | return `S`, `S.var` | `AtLocation(l)`, `PerLocation()` | `BassSensitivityBasisTests.Forms_ObeyTheSobolIdentities_AndTheLawOfTotalVariance`; `oracle/BassSensitivityDefectParityTests.DecompositionUnits_AreThePhysicalOnes` | R's indices; physical variances corrected | F | |
| F09 | return `Var.tot` | `AtLocation(l).Variance` | as F08 | Reference in physical units | F | |
| F10 | return `names.ind` | `Subsets`, `Inputs` | `DecompositionUnits_*` | R's names | F | |
| F11 | return `other` | `SobolDecomposition.OtherVariance` `sensitivity/SobolDecomposition.cs:149`, `OtherIndex` :173, `MeanOtherIndex` :183; `SobolSummary.OtherIndex`, `OtherVariance` `presentation/SobolSummary.cs:26` | `BassSensitivityTests.Other_IsTheVarianceOutsideTheReportedSubsets`; probe `BasisOther_IsTheVarianceLessTheReportedMains`; `BassPresentationSensitivityTests.DirectSummary_AgreesWithDrawSample_OverTheDefinedDraws_UnderEqualAndUnequalWeights` | Full decomposition; the draws' own reductions | F | Closed by G-1, `0ab1e4d`. The variance less the subsets reported, which agrees with R's where `nind` is null (F05). 0 on a zero-variance draw, and within the rounding tolerance |
| F12 | return `xx` | none: `Locations`, R's is an artificial 0–1 axis | | | D | Not a quantity; the hub's "no R object system" |
| F13 | cross-component terms | `ComponentCovariance`, `ComponentMainCovariance`, `ComponentTotalCovariance` | `BassSensitivityBasisTests.CrossComponentCovariances_AreKept_AndTheComponentOrderIsFree`, `CrossComponentTerms_*` | Tensor quadrature | F | |

### 1g. `calibrate.bassBasis`

| # | R item | C# | Test read | Reference | Verdict | Note |
|---|---|---|---|---|---|---|
| G01 | `y` | rides the call, `CalibrationModel` | `BassCalibrationConditionalTests.Target_IsTheDenseGaussianWithEveryPriorWrittenOut` | Dense Gaussian | F | |
| G02 | `mod` | `BassBasisEmulator`, `BassFitEmulator` `calibration/BassEmulators.cs` | `BassCalibrationEmulatorTests.BasisFit_PredictsAsTheFit_WithItsNoiseAndResidualsInPhysicalUnits`, `ScalarFit_WithFixedInputs_PredictsAsTheFit`, `DirectFunctionalFit_PredictsAsTheFit_AtDeclaredLocations` | The fit's own prediction | F | Beyond R: direct fits and fixed inputs |
| G03 | `type` | `FunctionalFormulation`, `MultivariateFormulation`, `TruncationTerm` `calibration/CalibrationCovariance.cs:24,45,64` | `oracle/BassCalibrationCovarianceParityTests.Functional_MatchesTheDenseReferences_AtC02sConditionings`, `Multivariate_MatchesItsDenseReference`; `oracle/BassCalibrationDefectsParityTests.BranchByShape_TheDeclaredFormulationHolds` | Dense references | F | Declared, not read from a shape |
| G04 | `sd.est`, `s2.df` | `InverseGammaVariance` `calibration/MeasurementModel.cs:45`, `HalfCauchyDeviation` :51 | `BassCalibrationMeasurementTests.InverseGamma_IsRsKernel_Normalized`, `HalfCauchy_IsNormalized_AndScaleEquivariant`; `oracle/BassCalibrationDefectsParityTests.S2DfSwitch_EachGroupKeepsItsOwnProperPrior` | R's kernel; closed form | F | s2.df = 0 → half-Cauchy, D09 |
| G05 | `s2.ind` | `VarianceGroup.Outputs` :28 | `BassCalibrationMeasurementTests.Groups_AreBoundByOutputIndex`, `Declarations_ThatAreNotAPartitionOrHaveBadPriors_AreRefused`; `BassCalibrationRecoveryTests.WithTwoVarianceGroups_RecoversThetaStar` | Controlled recovery | F | |
| G06 | `meas.error.cor` | `MultivariateFormulation.Correlation` | `oracle/BassCalibrationCovarianceParityTests.Multivariate_MatchesItsDenseReference` | Dense reference | F | |
| G07 | `bounds` | `ContinuousParameter` `calibration/CalibrationParameter.cs:26` | `BassCalibrationSamplerTests.Outputs_ArePhysical_Reproducible_AndThinned`; `oracle/BassCalibrationDefectsParityTests.Coordinates_ThetaIsReportedAsTheEmulatorSawIt` | Bounds held | F | Physical θ returned; R's normalized corrected |
| G08 | `discrep.mean` | `FunctionalDiscrepancy.Mean`, `MultivariateFormulation.DiscrepancyMean` | `BassCalibrationConditionalTests.Target_*`; `BassCalibrationCutSmcTests.AgreesWithTheNestedReference_WithADiscrepancy_TwoGroups_AndALevel` | Dense Gaussian | F | R's functional path drops it |
| G09 | `discrep.mat` | `FunctionalDiscrepancy.Basis`, `Scale` :57; `MultivariateFormulation.DiscrepancyCovariance` | `BassCalibrationConditionalTests.Discrepancy_DrawsFromItsGaussianConditional`; `BassCalibrationRecoveryTests.WithADeclaredDiscrepancy_RecoversThetaStar` | Conditional; recovery | F | Scale declared, D09 |
| G10 | `nmcmc` | `CalibrationSampling.BurnIn`, `Samples`, `Thin` `calibration/CalibrationSampler.cs:23-29` | `Outputs_ArePhysical_Reproducible_AndThinned` | Identity | F | R has no burn-in step |
| G11 | `temperature.ladder` | `CalibrationSampling.Tempering` :40 | `BassCalibrationSamplerTests.Tempered_EveryRungsLawIsItsLikelihoodPathTarget`, `Swap_IsTheLikelihoodPathsRatio`, `ProprietyRule_ALadderToZeroNeedsProperPriors` | Rung laws | F | Inverse temperatures; no default ladder (doc line 2029) |
| G12 | `decor.step.every` | `DecorrelationEvery` :33, `ThetaDecorrelation` :202 | `BassCalibrationSamplerTests.Untempered_TheChainsLawIsTheConditional` | The conditional's law | F | |
| G13 | `verbose` | none | | | D | As D-p04 |
| G14 | return `theta` | `CalibrationDraw.Theta` :364 | G07 | | F | The cold rung's only; R returns every temperature's. Recorded by the user's ruling, 2026-11-01 (`port-fitting-brief.md:487`) |
| G15 | return `s2` | `CalibrationDraw.Variances` | `BassCalibrationConditionalTests.Variances_MoveOnLogSigma2_ToTheirConditional` | Conditional | F | |
| G16 | return `discrep.vars` | `CalibrationDraw.Coefficients` | `BassCalibrationConditionalTests.Discrepancy_DrawsFromItsGaussianConditional`; `BassCalibrationSamplerTests.Outputs_ArePhysical_Reproducible_AndThinned` (the state rebuilt from the draw's coefficients) | Conditional; the likelihood | F | |
| G17 | return `count`, `count.decor`, `tau` | `ConditionalRun`'s move statistics; `ThetaAdaptation.LogScale` | `BassCalibrationSamplerTests.Outputs_*` | | FL | Statistics kept; the adaptation scale is on the block, not in `ConditionalRun` |
| G18 | return `ii` | `CalibrationDraw.EmulatorDraw` | `Outputs_*`; `BassCalibrationNestedTests.OuterDraws_AreASampleOfTheRetainedDraws` | Identity | F | |
| G19 | return `curr`, `dat` | none | | | D | R's working state; the config and call hold the rest. No restart is contracted for calibration |
| G20 | the inference target | `NestedCalibrationSampler`, `CutSmcSampler` | `BassCalibrationNestedTests.Pooled_IsTheEqualMixtureOfTheConverged_Conditionals`; `BassCalibrationCutSmcTests.Pooled_IsTheExactEqualMixture_ByQuadrature`; `TheCutIsKept_NoStageReweightsTheEmulator` | Quadrature | F | D09 |
| G21 | emulator and truncation uncertainty, discrepancy | `CalibrationCovariance` | `BassCalibrationCovarianceTests.Functional_ResidualSet_IsTheLowRankFormOfTheDenseSecondMoment`, `Multivariate_IsTheDenseSum`, `Units_ScaledOutputs_MoveTheLogLikelihoodByAConstant`; `oracle/BassCalibrationDefectsParityTests.Units_TheCovarianceIsInPhysicalUnits` | Dense | F | |

### 1h. P09–P15: the plot, print and summary methods

The contract asks for their data in regression, with rendering above maths. **Rendering and text formatting are remaining work, not gaps.**

| # | R method | Data in regression | Test read | Reference | Verdict | Remaining |
|---|---|---|---|---|---|---|
| H01 | P09 `plot.bass` traces | `BassFit.Traces()` `presentation/BassTraces.cs:18,24` | `BassPresentationTests.Traces_AreTheDrawStores_WithTheirSweeps` | The stores | F | Rendering |
| H02 | P09 fit, interval, residuals | `BassFit.TrainingFit` `presentation/BassTrainingFit.cs:36`: `Observed`, `Mean`, `Quantiles`, `Residuals`, `NoiseSd` | `BassPresentationTests.TrainingFit_IsThePredictionAtTheTrainingInputs`, `TrainingFit_Intervals_AreThePredictiveMixtures_NotRsPerDrawEndpoints`; `oracle/BassOutputDefectsParityTests.Intervals_AreThePredictiveMixturesQuantiles` | The mixture; R's fixture | F | Rendering. An SMC or flat-histogram fit has no traces, by design (`BassFit.cs:113`) |
| H03 | P10 `plot.bassBasis` traces | `BassBasisFit.Traces()` | `BassPresentationBasisTests.Traces_AreEachComponentsOwn` | Identity | F | Rendering |
| H04 | P10 fit, residuals | `BassBasisFit.TrainingFit` `output/BassComponents.cs:97` | `BassPresentationBasisTests.TrainingFit_IsThePredictionAtTheTrainingCases_InPhysicalUnits` | Prediction | F | Rendering. `pred=FALSE` is the caller not asking |
| H05 | P11 `plot.bassSob` scalar | `SobolAnalysisSummary` `presentation/SobolSummary.cs:68` | `BassPresentationSensitivityTests.DirectSummary_AgreesWithDrawSample_OverTheDefinedDraws_UnderEqualAndUnequalWeights`, `DirectSummary_LeavesTheZeroVarianceDrawsOut_OfTheIndicesAlone` | `DrawSample` | F | Rendering |
| H06 | P11 functional | `ConditionedSummary`, `BasisSobolSummary` :87 | `ConditionedSummary_IsOnePerCoordinate`, `BasisSummary_IsPerLocationAndIntegrated` | Per coordinate | F | Rendering |
| H07 | P12 `print.bass` | `BassFit.Metadata()` `presentation/BassSummary.cs:42` | `BassPresentationTests.Metadata_RoundTripsWhatTheFitHolds` | Identity: R's print fields | F | Formatting |
| H08 | P13 `print.bassBasis` | `BassBasisFit.Metadata()` :94 | `BassPresentationBasisTests.Metadata_RoundTripsTheDecompositionTheJointDrawsAndTheComponents` | Identity | F | Formatting |
| H09 | P14 `summary.bass` | `BassFit.Summary()` :78 | `BassPresentationTests.Summary_IsTheDrawsPooledOverTheChains_UnderTheFitsWeights` | The draws; R's range and mean | F | Formatting |
| H10 | P15 `summary.bassBasis` | `BassBasisFit.Summary()` :113 | `BassPresentationBasisTests.Summary_IsEachComponentsPosteriorSummary` | Each component's | F | Formatting |
| H11 | `quants` (P09, P10) | `TrainingFit(probabilities)` | `TrainingFit_RefusesOtherData_AndProbabilitiesOutsideTheOpenInterval` | | F | |

### 1i. The cross table

Each cell is a test read, each comparing against a reference outside the code under test. **Cell_D** is the probe `Cell_Direct_PredictsDecomposesAndRestarts(inputs, response)`; **Cell_B** is `Cell_Basis_ReconstructsDecomposesAndRestarts(inputs, response)`, run for PCA and for a supplied nonorthogonal basis. Their references: prediction at new inputs against the hinge written out by hand, sensitivity's variance and every main variance against enumeration and a midpoint grid, each main effect against the same grid's conditional mean, reconstruction through K by hand, restart bit for bit against the uninterrupted run.

| Representation | Inputs | Prediction at new inputs and coordinates | Sensitivity | Effects | Reconstruction | Restart |
|---|---|---|---|---|---|---|
| Scalar | Continuous | `BassPredictTests.NewCases_*`; Cell_D | `oracle/BassSensitivityDefectParityTests`; Cell_D | `BassSensitivityEffectTests`; Cell_D | Cannot apply: no basis | `BassRunResumeTests("defaults")`; Cell_D |
| Scalar | Categorical | Cell_D | Cell_D (exact enumeration) | Cell_D (exact) | Cannot apply | Cell_D |
| Scalar | Mixed | `BassPredictTests.UnseenLevel_*`; Cell_D | `BassSensitivityTests.SobolIdentities_*`, `Indices_MatchPickFreeze*`; Cell_D | `BassSensitivityEffectTests`; Cell_D | Cannot apply | `BassRunResumeTests("categorical")`; Cell_D |
| Direct | Continuous | `BassPredictTests.NewOutputLocations_*`; Cell_D | `BassSensitivityTests.FunctionalForms_*`, `ConditionedVariance_*`; Cell_D | `BassSensitivityEffectTests`; Cell_D | Cannot apply: no basis; the surface itself | `BassRunResumeTests("functional")`; Cell_D |
| Direct | Categorical | Cell_D | Cell_D | Cell_D | Cannot apply | Cell_D |
| Direct | Mixed | Cell_D | Cell_D | Cell_D | Cannot apply | Cell_D |
| PCA | Continuous | `BassPredictBasisTests.Latent_*`; Cell_B. New output coordinates: **refused**, the hub's extension | `BassSensitivityBasisTests.CrossComponentTerms_*`; Cell_B | `BassSensitivityBasisEffectTests.EachLocationAndTheIntegratedForm_AgreeWithTensorQuadrature`; Cell_B | `BassOutputDecompositionTests.Pca_*`; Cell_B | Cell_B |
| PCA | Categorical | Cell_B | Cell_B (exact) | Cell_B (exact) | Cell_B | Cell_B |
| PCA | Mixed | `BassCalibrationEmulatorTests.BasisFit_WithACategoricalTheta_EnumeratesAsItPredicts`; Cell_B | Cell_B | Cell_B | Cell_B | Cell_B |
| Supplied | Continuous | Cell_B | Cell_B | Cell_B | `SuppliedBasis_Nonorthogonal_*`; Cell_B | Cell_B |
| Supplied | Categorical | Cell_B | Cell_B (exact) | Cell_B (exact) | Cell_B | Cell_B |
| Supplied | Mixed | Cell_B | Cell_B | Cell_B | Cell_B | Cell_B |

Restart across a basis fit has no entry point of its own: each component's checkpoint continues it, and `BassComponents.Assemble` joins them. Before this audit no test exercised that path; Cell_B now does.

### 1j. Capability roll-up

| ID | Verdict | What limits it |
|---|---|---|
| P01 `bass` | FL | B-r10 and B-r11 (aggregated traces, by the user's ruling). A22 is closed; A14 and A23 are left out by decision |
| P02 `bassPCA` | F | |
| P03 `bassBasis` | F | |
| P04 `sobol` | F | E05 and E06 closed |
| P05 `sobolBasis` | F | F04 and F11 closed |
| P06 `calibrate.bassBasis` | F | G17's note |
| P07 `predict.bass` | F | |
| P08 `predict.bassBasis` | F | D-p09 closed |
| P09–P11 | F (data) | Rendering remains |
| P12–P15 | F (data) | Formatting remains |

## Part 2: the baseline

**The qualified existing-engine revision:** `cf37584` (2026-10-28), the last requalification, by the numerics checklist's report "Baseline requalified before stage 4". The chain before it: gate closed at `4c97de2` (B01–B03, B06), B04 and B05 at `fbe6d0d`, requalified at `2efa134`, `3a04601`, `cf37584`. Those suites are recorded runs of that work; this audit read the tests below and reran the three suites the brief names at the audited revision, not the baseline's four. Thirteen commits under `spline/bars`, `spline/baps`, `bgp`, `samplers/ensemble` and `linalg` follow `cf37584` (`git log cf37584..8f4d2da`); N14's verdict and BARS's groups, `DrawSample`'s two quantile fixes (`6d37bab`, `d0a0fd8`), the ladder's resize (`ec565f6`), the seams and the inverse temperature's rename (`f5ba714`) are among them, port changes each with its bit-identity pin or report, and **no requalification has run since**. The hub asks for one only before the first integration tests, so this is a note, not a gap.

### 2a. B01–B06

| # | Batch, findings | C# | Test read | Reference | Verdict | Note |
|---|---|---|---|---|---|---|
| Bb01 | B01 Diagnostics, R09 | `samplers/ensemble/ChainDiagnostics.cs`, `ConvergenceVerdict.cs`, `ChainEnsemble.cs` | `samplers/DiagnosticsPolicyTests.RHat_ChainsFrozenAtDifferentValues_IsFrozen`, `RHat_OneChainOrOneDraw_IsUnavailable`, `RHat_NonFiniteDraw_IsUnavailable`, `RHat_MomentsForm_MatchesHandComputedValue`; `samplers/ConvergenceVerdictTests.FrozenAndUnavailableReadings_Fail`, `OneChain_IsUndetermined` | Hand value; states | F | No invalid diagnostic can stop a run |
| Bb02 | B02 Mutable state, R05 and R07 | `bars/BarsEnsemble.cs` | `KernelOwnershipTests`; `AdaptationAcrossResamplingTests.AdaptiveTau_IsUnaffectedByAResamplerThatKeepsTheWeights` | Identical trajectory | F | |
| Bb03 | B03 BARS contracts, R01, R04, R08, R12, R13, R14 | `bars/*` | `CountLawTests.KnotCount_FollowsTheCountPrior`; `BarsDrawTests.Fit_IsTheMeanOfTheDrawsConditionalMeans_UnderRobustWeights`, `EveryDraw_KeepsKnotsCoefficientsAndSigma2`; `DesignAdmissionTests.NearlyCollinearDesign_IsRejected`, `ColumnScaling_DoesNotChangeAdmission`; `WeightedNormalMarginalTests.Marginal_ZeroResponse_IsRejected_WhileAnInterpolatingFitIsNot`; `SpanCoverageTests.NestedPeaks_CoverageStaysAProbability` | Flat-likelihood count law; dense | F | Baseline repairs, kept apart from the port |
| Bb04 | B04 BAPS, R02 and R03 | `baps/*` | `BapsMarginalEvidenceTests.ProfiledEvidence_Plateaus_SoAFlatPriorIsImproper`; `HierarchicalPSplineConditionalTests.VarianceDraws_MatchTheirFullConditionals`, `StationaryMeans_MatchQuadratureReference` | Quadrature, full conditionals | F | R03's commit `0ea0329` |
| Bb05 | B05 BGP, R06 and R10 | `bgp/*` | `BgpConcurrencyTests.ConcurrentFits_MatchSerialFits_BitForBit`; `BgpSupportTests.EmptySupport_IsRejected_BeforeAnyChainStarts`, `ConstantCoordinate_LeavesPriorAndFitUnchanged` | Bit identity; closed form | F | |
| Bb06 | B06 Numerical, R08 and R11 | `linalg/Cholesky*`, `BorderedCholesky.cs` | `linalg/DenseCholeskyTests.NonPositiveDefinite_StopsAtTheFailingPivot_WithoutFlooring`, `TinyScaleIdentity_IsFactoredExactly`; `linalg/BorderedCholeskyTests.RandomUpdates_TrackFreshFactorizations`, `RejectedProposals_LeaveTheAcceptedFactorUntouched`, `Solves_AndInnovations_HoldAfterRotations` | Fresh factorization | F | `Maths.LinAlg.Tests` is not one of this audit's three suites; its last recorded count is a record, 139/139 by the seams' chip |

### 2b. The inventory's observed defects

The inventory's table, its calibration row split by defect. Each row: the R oracle keeps R's value, the C# asserts the corrected one, as the parity rule says.

| # | R defect | R's value kept in | Corrected value asserted by | Verdict |
|---|---|---|---|---|
| Df01 | Truncation sampler selects columns by the row count | `oracle/BassOutputDefectsParityTests.TruncationSampler_DrawsWholeResidualsOverEveryCase` (`r_error`, `r_cases_drawn`) | Same test (draws every case); `BassPredictBasisTests.Predictive_DrawsComponentNoisePerCase_AndWholeResidualsOverEveryCase` | F |
| Df02 | One component-noise draw per posterior draw, recycled | `oracle/BassOutputDefectsParityTests.ComponentNoise_IsDrawnPerCase` (`r_max_spread_across_cases` 0, correlation 1) | Same test | F |
| Df03 | Truncation residual's sign and units | `oracle/BassOutputDefectsParityTests.Residual_IsObservedMinusReconstructed_AddedInStandardizedUnits` (`r_reconstruction_max_error`, `skewness_r`) | Same test, `BassOutputDecompositionTests.Residual_IsObservedMinusReconstructed_InStandardizedUnits` | F |
| Df04 | `perc.var = 100` fails | `oracle/BassOutputDefectsParityTests.FullRetention_IsTheNumericalRank` | Same; `BassOutputDecompositionTests.FullRetention_IsTheNumericalRank_WithZeroInSampleTruncation` | F |
| Df05 | Category weights ignored by `sobol` | `oracle/BassSensitivityDefectParityTests.CategoryWeights_AreRead` (`r_direct.sub_cnt`, `r_second_moment` 0.5) | Same test against enumeration | F |
| Df06 | Basis Sobol pair moments ignore truncation | `oracle/BassSensitivityDefectParityTests.TruncationBounds_AreClipped` (`r_basis`, five wrong, 0.45 against 0.1105903) | Same test; probe `TruncatedMeasure_SensitivityMatchesQuadrature(Pca)` | F |
| Df07 | Effects hard-code uniform and degree 1 | `oracle/BassSensitivityDefectParityTests.EffectMoments_UseTheDeclaredMeasure` (`r_c1_temp`), `EffectCurves_UseTheDeclaredDegree` (`r_max_abs_vs_degree1`) | Same tests | F |
| Df08 | Basis sensitivity in decomposition units | `oracle/BassSensitivityDefectParityTests.DecompositionUnits_AreThePhysicalOnes` (`r_var_total`, `r_effects`) | Same test | F |
| Df09 | Draws reordered and axes dropped | `oracle/BassOutputDefectsParityTests.DrawOrder_AxesAndNoise_FollowTheRequest` (`r_rows`, `r_dim_*`) | Same test | F |
| Df10 | Constant column scaling differs between training and prediction | `oracle/BassOutputDefectsParityTests.ConstantColumns_AreRefusedAsInputs_AndScaledByOneAsOutputs` (`r_prediction` NaN, Inf) | Same test | F |
| Df11 | Calibration: stale likelihood | `oracle/BassCalibrationDefectsParityTests.StaleLikelihood_TheCachedLikelihoodIsAlwaysTheStatesOwn` (`cached_minus_recomputed`) | Same test | F |
| Df12 | Calibration: normalized θ | `oracle/BassCalibrationDefectsParityTests.Coordinates_ThetaIsReportedAsTheEmulatorSawIt` | Same test | F |
| Df13 | Calibration: branch by shape | `oracle/BassCalibrationDefectsParityTests.BranchByShape_TheDeclaredFormulationHolds` | Same test | F |
| Df14 | Calibration: response-scale mismatch | `oracle/BassCalibrationDefectsParityTests.Units_TheCovarianceIsInPhysicalUnits` | Same test; `oracle/BassCalibrationCovarianceParityTests.Units_ThePhysicalLikelihoodIsRsPhysicalOne` | F |
| Df15 | Calibration: `s2.df` switch for all groups | `oracle/BassCalibrationDefectsParityTests.S2DfSwitch_EachGroupKeepsItsOwnProperPrior` (`r_group2`) | Same test | F |
| Df16 | Calibration: tempered Jacobians and swaps | `oracle/BassCalibrationDefectsParityTests.Tempering_TheKernelsAndTheSwapReachTheLikelihoodPath` | Same test | F |
| Df17 | Fitting targets: fixed thirds, partial collapse, tempering, count normalization | `oracle/BassTemperingParityTests.Swap_RsDensity_IsADefectFixture`; `r/oracles/bass_enumeration_oracle.R` | `BassSweepTests.ConstantLikelihood_SamplesThePriorOnF` (π(1)/π(0) within 5%, not R's ×3); `oracle/BassFitParityTests.Fit_MatchesTheOraclesPosterior`, `oracle/BassChainParityTests.Chain_MatchesTheOraclesPosterior`, `oracle/BassTemperingParityTests.Swap_BalancesOnTheOraclesRungLaws` | F |
| Df18 | Interval quantiles of Gaussian endpoints in `plot.bass` | `oracle/BassOutputDefectsParityTests.Intervals_AreThePredictiveMixturesQuantiles` (`r`) | Same test | F |

## Part 3: the decisions

### 3a. D01–D10

| # | As recorded | Code | Test read | Reference | Verdict |
|---|---|---|---|---|---|
| Dd01 | D01 term measure, ordering factor, knots on distinct values, global F, move probabilities per state, Wallenius, change redraw, categorical subsets | `BassStructurePrior.cs:102,333`, `BassTermProposal.cs`, `BassSweep.cs`, `distributions/Wallenius.cs` | `BassProposalTests.TermPrior_SumsToOne_OverTheTermSpace`, `TermPrior_IsD01sConstant`, `BirthProbability_ByHand`, `DeathRatio_ReversesTheBirthRatio`, `ChangeProposal_IsItsFactorsPriorRedrawn`; `BassSweepTests.ConstantLikelihood_SamplesThePriorOnF`, `Schedules_AtTheBoundaries`; probe `Levels_AboveSixtyFour_DrawAndScoreUnderBothMeasures` | Enumeration, by hand | F |
| Dd02 | D02 global truncation, λ integrated out, λ given M at retention | `BassStructurePrior.cs:178,189`, `BassSweep.cs:181` | `BassProposalTests.CountPrior_IsTheNegativeBinomial`; `BassSweepTests.ConstantLikelihood_SamplesThePriorOnF` (E[λ]); `BassRunDiagnosticsTests.MassAtCap_ReadsTheDraws` | Closed form, enumeration | F |
| Dd03 | D03 τ's family, level centring, σ² scaled, two axes, ridge, proper defaults | `BassConfig.cs:33,36,43,349`, `BassTarget.cs`, `gprior/TauPrior.cs`, `gprior/CollapsedGPrior.cs` | `BassTargetTests.TauPrior_IsTheUnitInformationGamma`, `Shift_MovesOnlyTheZeroCentre`, `Scale_MovesEveryStructureAlike`, `ExactFits_AreRefusedUnderJeffreysOnly`; `BassCoefficientPriorTests.TauPosterior_OnAFixedStructure_MatchesQuadrature`; `BassRidgeTests.Marginal_MatchesTheDenseComputation`; `oracle/BassTargetParityTests.Marginals_MatchTheOraclesPosterior`, `ZeroCentre_MatchesRsStatistics` | Quadrature, oracle enumeration, R's statistics | F |
| Dd04 | D04 full collapse, τ by slice, composition at retention | `BassSweep.cs:181,555`, `samplers/blocks/SliceSampler.cs` | `BassSweepTests.FixedStructure_TauAndComposition_MatchQuadratureAndAJointGibbsReference`; `BassSweepSlotTests.TheSlots_AreD04sOrder`; `BassPointwiseTests.BothContents_MatchDenseReferencesAndRefits` | Quadrature, joint Gibbs | F |
| Dd05 | D05 likelihood path, integrated by default, linear and powered selectable, reference rung, DEO, stepping-stone, SMC on the same sequence | `BassLadder.cs`, `BassTarget.cs`, `BassEvidence.cs`, `BassSmc.cs`, `BassSmcRun.cs` | `BassTargetTests.Forms_AndTheirEndpoints`; `BassTemperingTests.Ladder_EachRungSamplesItsTemperedLaw`, `Evidence_MatchesTheExactLogZ`; `oracle/BassTemperingParityTests.Ladder_EachRungMatchesTheOraclesTemperedLaw`, `Swap_BalancesOnTheOraclesRungLaws`, `Evidence_MatchesTheOraclesLogZ`; `oracle/BassSmcParityTests.Kernel_LeavesEachStagesLawInvariant`, `Fit_NormalizerMatchesTheOraclesLogZ` | Oracle enumeration | F |
| Dd06 | D06 order-free angle gate on the centred Gram, s_min = 10⁻¹⁰, support count, rebuild before refusal, R's rule a fixture | `BassDesign.cs` | `BassDesignTests.TheGate_IsOrderFree_WhereANewestColumnRuleIsNot`, `EachRefusalReason_HasItsCase`, `ANearExactFit_FormsTheCandidatesResidualExplicitly`, `ALongRun_KeepsTheResidualAndTheFactor`; `oracle/BassAdmissionParityTests.ScalarAngles_MatchTheExplicitlyCentredQr`, `ANearConstantIndicator_KeepsItsAngle`, `DirectAngles_MatchTheMaterializedDesignsQr` | QR by R | F |
| Dd07 | D07 units, leave-one-case-out, per-case noise, joint draws, covariance, axes, intervals, mean surface | `output/*`, `posterior/*` | `BassOutputDecompositionTests.LeaveOneOut_EqualsExplicitRefits`, `TwoCases_LeaveOneOutIsTheirDifference`, `TruncationCovariance_MatchesADenseComputation`; `BassPredictBasisTests.Covariance_InPhysicalUnits_MatchesDenseComputations`, `MeanSurface_ReconstructsEachComponentsMeanSurface`; `BassPredictTests.MeanSurface_IsTheGPriorsConditionalMean`; `BassOutputComponentsTests.WeightedComponents_AreResampledWithinEach`, `UnequalCounts_AreThinnedEvenlyToTheShortest`; `BassOutputDensityTests.DiagonalTruncation_IsADenseGaussianInPhysicalUnits`; `oracle/BassOutputTruncationParityTests.Residuals_MatchTheOracle` | Refits, dense, R | F |
| Dd08 | D08 uniform default, empirical option, Student by quadrature, zero variance undefined, effects for the measure and degree | `sensitivity/*`, `distributions/IIntervalMoments.cs` | `BassSensitivityMomentTests.TermMoments_MatchTensorQuadratureOfTheBasisColumns`, `EmpiricalMarginal_IsTheAverageOverTheTrainingValues`; `BassSensitivityTests.ZeroVariance_LeavesTheIndicesUndefined`, `Status_ClampsWithinTheTolerance_AndFlagsBeyondIt`; `oracle/BassMomentParityTests.HingeMoments_MatchQuadrature_UnderEveryMeasure`; probes `Cell_*`, `TruncatedMeasure_*`, `ZeroProbabilityLevel_*` | Quadrature, enumeration | F |
| Dd09 | D09 cut posterior, formulations declared, half-Cauchy, discrepancy scale declared, categorical parameters enumerated, cut SMC | `calibration/*` | `BassCalibrationConditionalTests.Target_IsTheDenseGaussianWithEveryPriorWrittenOut`, `Categorical_EnumeratedConditional_IsBruteForce`, `Categorical_AboveTheCap_SingleSiteUpdatesReachTheEnumeratedLaw`, `RsKernel_IsImproperAndUnitDependent_WhereTheHalfCauchyIsNeither`; `BassCalibrationCutSmcTests.Pooled_IsTheExactEqualMixture_ByQuadrature`, `AgreesWithTheNestedReference_OnBassFits`; `BassCalibrationRecoveryTests.BasisEmulator_RecoversThetaStar` | Dense, quadrature, recovery | F |
| Dd10 | D10 adaptation policy; a verdict on every run, three groups, a declared scope, lengths as a budget | `samplers/ensemble/ConvergenceVerdict.cs`, `BassDiagnostics.cs:69-83`, `BassConfig.cs:412` | `BassVerdictTests.Functionals_AreNamed_AndGroupedByD10`, `Fit_CarriesTheVerdict_AndTheModeSeparation`, `PredictiveScope_StopsNoLaterThanAll_AndReportsEveryGroup`; `samplers/ConvergenceVerdictTests.FrozenAndUnavailableReadings_Fail`; `BassTemperingTests.Ladder_RefusesWhatItCannotDo` (frozen after burn-in); `BassRunResumeTests.Resume_IsBitIdentical` (adaptation replayed) | Generic verdict | F |

### 3b. The defaults, mechanically

Every BASS row of `docs/Maths/BayesianRegression.md`'s Configuration table against the code. **No difference in value.**

| Doc line | Setting | Doc | Code | Same |
|---|---|---|---|---|
| 2016 | Degree / MaxInteraction / MaxCategorical / MaxFunctional | 1 / 3 / 3 / 3 | `BassConfig.cs:257-266` | yes |
| 2017 | MaxBasis | null: min(1000, N − 1) | `BassConfig.cs:269`, `BassResolved.cs:160,169` | yes |
| 2018 | Support | `Count`, min(20, 0.1·n) per part | `BassConfig.cs:272,331`, `BassResolved.cs:358` | yes |
| 2019 | AngleFloor / MaxVarianceInflation | 10⁻¹⁰ / null | `BassConfig.cs:275,278` | yes |
| 2020 | CategoricalSubsets / OrderWeight / VariableWeight | `SizeUniform` / 5 / 5 | `BassConfig.cs:281,284,287` | yes |
| 2021 | Count / Coefficients | (10, 10) / `G()`: Gamma(½, rate N/2) | `BassConfig.cs:29,33,334`; `TauPrior.cs:82` | yes |
| 2022 | Sigma2 / Intercept | `ScaledToResponse(1e-3)` / `Mean` | `BassConfig.cs:36,43` | yes |
| 2023 | Chains / MasterSeed / BurnIn / Samples / Thin | 4 / 0 / 9000 / 1000 / 1 | `BassConfig.cs:385-398` | yes |
| 2024 | RHatTarget / BatchSize / MinEssPerChain / ReportEvidence | 0 / 0 / 100 / true | `BassConfig.cs:403-409,448`; `BassContractTests.Config_HasTheDecidedDefaults` | yes |
| 2025 | StoppingScope | `All` | `BassConfig.cs:412` | yes |
| 2026 | Smc | null; 1,000, 4, 0.99, 5, 0.5, multinomial, integrated, 1,000, false | `BassConfig.cs:549-574` | yes |
| 2027 | Calibration Parameters / Measurement / Covariance | required | `calibration/CalibrationConfig.cs:17-23` | yes |
| 2028 | EnumerationCap | 1,024 | `CalibrationConfig.cs:27` | yes |
| 2029 | Calibration Sampling | 2,000, 2,000, 1, 100, null | `calibration/CalibrationSampler.cs:23-40` | yes |
| 2030 | Reference | null, 4, 4,000, 250, 1.01, 100 | `calibration/NestedCalibration.cs:18-30` | yes |
| 2031 | Cut Smc | null; 1,000, 0.99, 5, 0.5, multinomial, as C3's, Tour, 1,000 | `calibration/CutSmc.cs:17-40` | yes |
| 2032 | Moves | `EqualOverAvailable`; `PriorCancelling(C = 0.4)`, C ∈ (0, ½] | `BassConfig.cs:415,141,647` | yes |
| 2033 | LocalSteps | null; (0.5, 0.02) | `BassConfig.cs:418,657` | yes |
| 2034 | FactorMoves | null; (1, 1, 1) | `BassConfig.cs:422,676` | yes |
| 2035 | MonitoredPoints | up to three at the quartiles, a direct response's at the middle | `BassDiagnostics.cs:100-107` | yes |
| 2036 | FlatHistogram | null; linear with the reference rung, self-sizing, 4 ladders | `BassConfig.cs:445,619-628` | yes |
| 2037 | Tempering | null; likelihood, integrated, self-sizing from 8 within 64, hot end auto, adaptive, DEO, no rung | `BassConfig.cs:483-511`, `BassResolved.cs:287` | yes |
| 2038 | Split | null; 0.9 | `BassConfig.cs:515,683` | yes |
| 2039 | VarianceShare / ScaleByColumnSd / Residual | 0.99 / false / `LeaveOneCaseOut` | `output/BassDecomposition.cs:257,266,272` | yes |
| 2040 | Predict options | all / true / 0 / `Refuse` / null / false | `posterior/BassPrediction.cs:17-34` | yes |
| 2041 | Measures / OutputMeasures / MaxOrder / SubsetBudget | null / null / null / 10,000 | `sensitivity/BassSensitivity.cs:50-62`; the uniform defaults read back by probes `Cell_*` | yes |
| 2042 | Surface / Functional / Coordinates / RoundingTolerance | `Mean` / `Integrate` / null / 10⁻¹⁰ | `sensitivity/BassSensitivity.cs:65-76` | yes |
| 2043 | GridPoints / CellBudget | 100 / 10⁷ | `sensitivity/EffectFunctions.cs` | yes |

`docs/Maths/Samplers.md`'s Configuration table has no row naming BASS. Its rows BASS reads (`ChainEnsemble.MinEssPerChain` = 100; `ParallelTempering` schedule DEO; `LadderEvidence` stepping-stone; `FlatHistogramSpec`) match the code BASS reads. `TemperingSpec.Schedule` is DEO only, by `BassConfig.cs:241`.

### 3c. The register's Corrected entries

| # | Entry | Code | Test read | Verdict |
|---|---|---|---|---|
| Co1 | Move probabilities per state (D01) | `BassSweep.cs` | `BassSweepTests.ConstantLikelihood_SamplesThePriorOnF`, `Schedules_AtTheBoundaries` | F |
| Co2 | Unordered terms, one (M + 1) (D01); KHAOS's double count a fixture | `BassStructurePrior.cs:178` | `BassProposalTests.DeathRatio_ReversesTheBirthRatio`; `BassSweepTests.ConstantLikelihood_SamplesThePriorOnF` | F |
| Co3 | τ's shape fixed (D03; F-O8) | `gprior/TauPrior.cs:82` | `BassTargetTests.TauPrior_IsTheUnitInformationGamma`; `oracle/BassTargetParityTests.Marginals_MatchTheOraclesPosterior` | F |
| Co4 | τ's rate unit information (D03) | `TauPrior.cs` | `BassContractTests.Resolved_MatchesDirectComputation` (rate 20 at N = 40) | F |
| Co5 | Refuse unsupported columns rather than drop them (Q6, D06) | `BassDesign.cs` | `BassDesignTests.EachRefusalReason_HasItsCase`; `BassRidgeTests.Admission_IsTheGForms` | F |

## Part 4: the validation policy

### 4a. The independent checks

| # | Item | Scalar | Direct | Basis | Verdict |
|---|---|---|---|---|---|
| V01 | Boundary proposal ratios | `BassSweepTests.Schedules_AtTheBoundaries`; `BassSweepTests.ConstantLikelihood_SamplesThePriorOnF` (π(1)/π(0)) | same chain; `BassPriorSamplerTests.Draws_MatchThePriorOnF_ForADirectFunctionalFit` | Cannot apply: each component is scalar | F |
| V02 | Enumerable model spaces | `oracle/BassChainParityTests.Chain_MatchesTheOraclesPosterior`; `oracle/BassFitParityTests.Fit_MatchesTheOraclesPosterior`; `BassLocalStepTests.DetailedBalance_Enumerated`, `BassFactorMoveTests.DetailedBalance_Enumerated` | `BassTemperingTests.Evidence_MatchesTheExactLogZ_OnADirectFunctionalResponse`; `BassFlatHistogramTests.Fit_OnFixturesBAndC_ReadsTheOraclesEvidence` | Cannot apply | F |
| V03 | Gaussian conditionals | `BassSweepTests.FixedStructure_TauAndComposition_MatchQuadratureAndAJointGibbsReference`; `BassRidgeDrawTests.FixedStructure_ComposedDraws_MatchAJointGibbsReference` | `BassSweepTests.Compose_AlignsOnDirectFunctionalData`; `BassPointwiseTests.BothContents_MatchDenseReferencesAndRefits` | Components: scalar's; `BassPredictBasisTests.MeanSurface_ReconstructsEachComponentsMeanSurface` | F |
| V04 | Precision sampling | `linalg/BorderedCholeskyTests.Solves_AndInnovations_HoldAfterRotations`; `linalg/DiagonalPlusLowRankTests.Draws_HaveExactlyTheirCovariances`; `BassSweepTests.FixedStructure_TauAndComposition_MatchQuadratureAndAJointGibbsReference` (through the joint Gibbs reference) | same factor | Calibration's: `BassCalibrationConditionalTests.Discrepancy_DrawsFromItsGaussianConditional` | F |
| V05 | Empirical covariance | Cannot apply | Cannot apply | `BassOutputDecompositionTests.TruncationCovariance_MatchesADenseComputation`; `BassPredictBasisTests.Covariance_InPhysicalUnits_MatchesDenseComputations`; `BassCalibrationCovarianceTests.Functional_Diagonal_IsRsCentredSampleVariance`, `Functional_ResidualSet_IsTheLowRankFormOfTheDenseSecondMoment`, `Multivariate_IsTheDenseSum` | F |
| V06 | Functional contractions | Cannot apply | `BassBasisTests.FunctionalDesign_SeparableForm_MatchesTheExplicitDesign`; `BassDesignKernelTests.FunctionalTerms_FactorsWithTheScale_ReproduceBothSides`; `oracle/BassAdmissionParityTests.DirectAngles_MatchTheMaterializedDesignsQr`; `BassOutputProductTests.ADeclaredIdentity_OnADirectFit_IsTheDesignsOwnArithmetic` | Cannot apply | F |
| V07 | Sensitivity identities | `BassSensitivityTests.SobolIdentities_HoldOnEveryDraw`, `Indices_MatchPickFreezeMonteCarloOnAFittedDraw`; probes `Cell_Direct_*` | `BassSensitivityTests.FunctionalForms_ObeyTheLawOfTotalVariance`, `ConditionedVariance_MatchesTensorQuadratureAtEachCoordinate`; probes `Cell_Direct_*` | `BassSensitivityBasisTests.Forms_ObeyTheSobolIdentities_AndTheLawOfTotalVariance`, `CrossComponentTerms_MatchTensorQuadratureOfTheReconstruction`; probes `Cell_Basis_*` | F |
| V08 | Calibration targets | `BassCalibrationEmulatorTests.ScalarFit_WithFixedInputs_PredictsAsTheFit`; `BassCalibrationConditionalTests.Target_IsTheDenseGaussianWithEveryPriorWrittenOut` | `BassCalibrationRecoveryTests.DirectEmulator_RecoversThetaStar`; `BassCalibrationEmulatorTests.DirectFunctionalFit_PredictsAsTheFit_AtDeclaredLocations` | `BassCalibrationEmulatorTests.BasisFit_PredictsAsTheFit_WithItsNoiseAndResidualsInPhysicalUnits`; `BassCalibrationRecoveryTests.BasisEmulator_RecoversThetaStar`; `oracle/BassCalibrationCovarianceParityTests.Functional_MatchesTheDenseReferences_AtC02sConditionings` | F |

### 4b. The SMC tests

| # | Item | Generic (`Maths.Samplers`) | BASS fitting | Cut calibration | Verdict |
|---|---|---|---|---|---|
| V09 | A known initial distribution | `samplers/SmcTests.ConjugateGaussian_RecoversPosteriorAndEvidence` | `BassPriorSamplerTests.Draws_MatchThePriorOnF`; `BassSmcTests.Fit_UnderEachOption_ReachesThePosterior` | `BassCalibrationBridgeTests.DrawPrior_IsThePrior_AndRefusedWhenImproper` | F |
| V10 | Incremental weights | `samplers/SmcTests.WithoutResamplingOrMoves_EvidenceTelescopesToImportanceSampling`, `LargestStep_KeepsTheConditionalEssAtItsTarget_OnCachedIncrements`, `Weights_NormalizeAtExtremeScales` | `oracle/BassSmcParityTests.Kernel_LeavesEachStagesLawInvariant` | `BassCalibrationBridgeTests.Bridge_IsTheGeometricMixtureOfItsTwoConditionals`, `BassCalibrationCutSmcTests.FirstStep_IsTheBisectionOnTheCachedIncrements` | F |
| V11 | Normalizer ratios | `samplers/SmcTests.ConjugateGaussian_RecoversPosteriorAndEvidence`; `samplers/LadderEvidenceTests.SmcNormalizerError_IsNaN_AfterAResamplingItDoesNotDescribe` | `oracle/BassSmcParityTests.Fit_NormalizerMatchesTheOraclesLogZ` | Cannot apply: a cut has no one normalizer | F |
| V12 | Resampling | `samplers/SmcTests.Resampling_OffspringAverageNTimesTheWeights`, `GenealogyVariance_IsUnbiased_UnderMultinomialResamplingAtEveryStage` | `BassOutputComponentsTests.WeightedComponents_AreResampledWithinEach` | `BassCalibrationCutSmcTests.WeightedDraws_TheResampledSetIsToured_AndPooledAtOneOverSPerVisit` | F |
| V13 | State ownership | `samplers/SmcTests.Results_DoNotDependOnParallelism` | `BassSmcTests.Fork_IsAnExactCopy_AndIndependentOfItsParent`, `Run_InParallel_MatchesTheSerialRun` | `BassCalibrationCutSmcTests.Run_IsDeterministicInTheSeed_WhateverTheParallelism` | F |
| V14 | Populations whose modes are known | `samplers/SmcTests.BimodalLikelihood_RecoversTheModesMasses` | `oracle/BassSmcParityTests.Fit_MatchesTheOraclesPosterior`; `BassFlatHistogramTests.Fit_OnFixtureA_MatchesTheEnumeration` | `BassCalibrationCutSmcTests.Pooled_IsTheExactEqualMixture_ByQuadrature` | F |

### 4c. The essential edge cases

| # | Edge case | C# | Test read | Reference | Verdict | Note |
|---|---|---|---|---|---|---|
| Ed01 | Intercept-only models | `BassSweep.cs:181`, `posterior/BassSurface.cs` | probes `InterceptOnly_PosteriorIsItsClosedForm`, `InterceptOnly_EveryViewReadsIt`; `BassPointwiseTests.EmptyStructure_MatchesUnderBothContents` | Closed-form IG; a₀'s centre | F | The probes are new: before them M = 0 had pointwise and mock coverage, no whole fit |
| Ed02 | One predictor | `BassData.cs:18` | probe `OnePredictor_FitPredictAndSensitivity` | Hand evaluation; S = T = 1 | F | R's `sobol` refuses; here it holds the whole variance |
| Ed03 | One response component | `output/BassComponents.cs` | `BassOutputDecompositionTests.OneLocation_DecomposesToOneComponent`; `BassSensitivityBasisTests.OneComponent_IsTheDirectDecompositionScaled`; `BassSensitivityBasisEffectTests.OneComponent_ReproducesItsDirectEffects_ScaledByTheBasis`; `BassPredictBasisTests.OneLocation_KeepsItsAxis`; probe `OneLocation_DirectIsRefused_AScalarFitStandsIn` | Direct decomposition | F | A direct fit on one location is refused, its output variable constant; a scalar fit is that case |
| Ed04 | One case | Refused: `ContinuousPredictor` `BassData.cs:84`, `CategoricalPredictor` :126, `FromPca` | probe `OneCase_IsRefusedWithItsReason_AndTwoCasesFit`; `BassOutputDecompositionTests.Refusals_SayWhy`, `TwoCases_LeaveOneOutIsTheirDifference` | Refusal by reason; a two-case fit runs | FL | **A fit needs two cases**: one case makes every input constant, which D07 and D01 refuse. Prediction at one case is whole: `BassPredictTests.DrawOrder_AndSingletonAxes_AreKept` |
| Ed05 | One draw | `posterior/*` | probe `OneDraw_EveryViewReadsIt`; `BassCalibrationCutSmcTests.OneDraw_IsTheFirstLegAlone`; `BassPredictTests.DrawOrder_AndSingletonAxes_AreKept` | Identity; verdict not converged | F | Probe new, through every view |
| Ed06 | Constant columns | Refused as inputs `BassData.cs:84`; scaled by one as outputs | `BassContractTests.Data_RefusesItsDefects_WithReasons`; `BassOutputDecompositionTests.ConstantColumn_IsCentred_AndScaledByOne`; `oracle/BassOutputDefectsParityTests.ConstantColumns_AreRefusedAsInputs_AndScaledByOneAsOutputs` | R's NaN a fixture | F | D07 |
| Ed07 | Constant responses | Refused by convention `BassResolved.cs:124-132`; `FromPca` | `BassContractTests.Resolved_RefusesWhatTheConventionCannotFit`; `BassTargetTests.ExactFits_AreRefusedUnderJeffreysOnly`; `BassOutputDecompositionTests.Refusals_SayWhy` ("no component") | Refusal by reason; D03 | F | |
| Ed08 | Nearly dependent terms | `BassDesign.cs` | `BassDesignTests.TheGate_IsOrderFree_WhereANewestColumnRuleIsNot`, `EachRefusalReason_HasItsCase`; `oracle/BassAdmissionParityTests.ANearConstantIndicator_KeepsItsAngle`; `linalg/BorderedCholeskyTests.DeathOfANearDuplicate_KeepsItsPartnersInverseDiagonalAccurate`; `BassRidgeTests.BorderedCandidates_MatchFreshDecompositions` | QR by R; fresh factors | F | |
| Ed09 | More than 64 category levels | `CategoricalFactor` `BassTerm.cs:17` holds sorted `int[]`; `BassStructurePrior.cs:298,333` | probes `Levels_AboveSixtyFour_DrawAndScoreUnderBothMeasures` (70, 300), `Levels_AboveSixtyFour_FitPredictMeasureAndRestart` (70) | Closed-form prior; hand evaluation; a grid | F | The representation has no 64-level seam: no bitset exists, so the hub's "bitsets with a representation above 64" is met by not using one. Before the probes the only runs above 64 were the env-gated cost qualification's (L up to 1,024). Membership costs 9–65 ns per case by the cost report |
| Ed10 | Category probabilities of zero | `sensitivity/InputMeasures.cs:286` | `BassSensitivityMomentTests.TermMoments_MatchTensorQuadratureOfTheBasisColumns`; `oracle/BassMomentParityTests.CategoricalMoments_MatchEnumeration`; probe `ZeroProbabilityLevel_SensitivityMatchesEnumeration` | Enumeration | F | Whole mass on one level leaves V = 0: undefined, not 0 |
| Ed11 | Bounded and truncated input measures | `InputMeasures.cs:20,27`; `sensitivity/BassMoments.cs` | `BassSensitivityMomentTests.TermMoments_MatchTensorQuadratureOfTheBasisColumns`; `oracle/BassSensitivityDefectParityTests.TruncationBounds_AreClipped`, `EffectMoments_UseTheDeclaredMeasure`; probe `TruncatedMeasure_SensitivityMatchesQuadrature` | Quadrature; a grid | F | Direct and basis, through a fit. A measure wider than the training range: probe `WideMeasure_ExtrapolatesTheHinges` |
| Ed12 | Full PCA rank | `output/BassDecomposition.cs:354` | `BassOutputDecompositionTests.FullRetention_IsTheNumericalRank_WithZeroInSampleTruncation`, `FullRank_WithFewerLocationsThanCases`; `oracle/BassOutputDefectsParityTests.FullRetention_IsTheNumericalRank`; probe `FullPcaRank_FitsPredictsAndDecomposes` | Rank; reconstruction | F | Probe new, through the fit |
| Ed13 | Rejected proposals leave the state untouched | `BassDesign.cs` candidate workspaces; `BassSweep.cs` | `linalg/BorderedCholeskyTests.RejectedProposals_LeaveTheAcceptedFactorUntouched`; `BassDesignTests.Load_MatchesTheIncrementalPath_AndARefusalLeavesTheDesign`, `RandomMoves_MatchAFreshReference`; `BassMoveOptionsTests.DuplicateAndSupport_AreRefusedAsTheirKind_TheDesignUntouched`, `NullMoves_LeaveTheStateTheDesignAndTheTalliesConsistent`; `BassSweepTests.Design_FollowsTheChain` | Fresh design every 100 sweeps | F | |

## Part 5: the benchmark list

Run no benchmark. Read against [the cost half's report](performance-cost-brief.md#reports), "The hub's list" (`performance-cost-brief.md:367-381`). These verdicts rest on that report by the brief's instruction; its harness is `tests/regression/BassCostSweep.cs` and `BassCostQualification*.cs`, gated by `BASS_QUALIFICATION`, left unset here.

| Axis or measure | Covered by | Left | Reading |
|---|---|---|---|
| Case count | n 200–50,000; n × M_max | | Covered |
| Predictor count | p 2–200 | past p = 50 at a lower M | Covered, limited |
| Active terms | M_max 5–200 | M reached 94 of 200 at n = 1,000 | Covered, limited |
| Interaction order | K 1–5 | the chain holds mostly one-way terms, mean order 1.05 | Covered, limited |
| Response size and rank | direct n_x × m; basis r 1–30 | | Covered |
| Category cardinality | L 2–1,024 | | Covered |
| Allocation and peak memory | bytes per sweep and part, GC, retained, peak | | Covered |
| Factor work | proposal, commit, rebuild, death, change | | Covered |
| Basis and cross-product scans | column and cross products per unit | | Covered |
| Sensitivity moment work | inside Sobol and effects per draw | **`BassMoments.On` alone is not isolated** | Open (known; in [the statistical half](performance-statistical-brief.md)'s cost blocks, waits on a quiet window) |
| Statistical error per unit time | none | **[its brief](performance-statistical-brief.md), not run** | Open (known; waits on the user's machine) |

## Part 6: the records

### 6a. The open boxes

| Box | Why open | Stale? |
|---|---|---|
| N15 (`port-numerics-checklist.md:27`) | Built when real sizes make it binding; priced by the cost half (the capacitance count about 7× below one Jacobi sweep), not built | Not stale: open by its own condition |
| N17 (`port-numerics-checklist.md:29`) | Phases 1 and 2 built; Phase 3 is SPC's | Not stale; BASS's part is done, and the box spans SPC |
| F-M2 (`port-fitting-checklist.md:19`) | Phase 1 built, Phase 2 ran, no arm cleared; the run lengths, the worse modes and D10's treatment of chains that never agree are the user's next calls | Not stale: the calls are open |
| F-X4 (`port-fitting-checklist.md:56`) | Proposal kernels behind one contract; built with KHAOS, an extension | Not stale |
| `extensions-checklist.md:9-26` | E1–E18 | Open by design. This brief's "E1–E15" lags (Fa4) |

### 6b. The API and source-map rows

- **Rows naming what no longer exists:** none. Every capitalized identifier in the BASS rows of "API at a glance" exists under `src/maths`, checked by script.
- **Files on disk the BASS source map does not name** (doc line 2138): `BassObservationModel.cs`, `DesignKernel.cs`, `DesignKernels.cs`, the extension seams.
- **Public types under `spline/bass/` no API or configuration row names:** 71 of 242 at the audit, by script; 60 of 243 on 2026-11-02, with the seams, `BassBasisPrediction` and `BasisPredictiveParts` named.
  - *The seams, owed at the audit and named since 2026-11-02:* `IBassObservationModel`, `WorkingResponse`, `GaussianObservation`, `BassSweepSlot`, `BassSweepSlots`, `IDesignKernel`, `DesignKernels`, `DesignFactor`, `HingeDesignFactor`, `SubsetDesignFactor`.
  - *Entry points and shapes a caller reads:* `BassBasisPrediction`, `BassSensitivity`, `BassBasisSensitivity`, `BassBasisEffects`, `Effect`, `EffectFunctions`, `EffectGrid`, `ProductMeasure`, `InputMeasure`, `MixtureComponent`, `ContinuousMeasure`, `CalibrationSampling`, `ReferenceSampling`, `CovarianceFormulation`, `VariancePrior`, `MeasurementErrorModel`, `CalibrationParameter`, `CalibrationChain`, `BassPredictor`, `BassResponse`, `PredictionColumn`, `RangeScale`, `TrainingLevels`, `PredictorScale`, `MoveSchedule`, `BassProgress`, `BassPcaAccount`, `BassTruncationSummary`, `IntervalSource`, `DecompositionStatus`, `DecomposedSurface`. Eight of them appear in the doc's prose, none in a row.
  - *Readouts and building blocks:* `TermSupport`, `BassChainAccount`, `BassMonitoredValue`, `BassOrderSupport`, `BassSmcSpread`, `BassStart`, `BassSupportTally`, `CaseInputs`, `OutputInputs`, `BassAutoPick`, `BassFactorShares`, `BassStepWindows`, `BassSmcAccount`, `BassPart`, `CalibrationBlock`, `CategoricalBlock`, `DiscrepancyBlock`, `VarianceBlock`, `CalibrationEvaluation`, `ThetaAdaptation`, `ThetaDecorrelation`, `ThetaWalk`, `CutSmcBlockRates`, `CutSmcOrder`, `InnerRun`, `BasisUncertainty`, `ReconstructionDomain`, `BassDataShape`, `SubsetCovariance`, `VariableCovariance`.
- `BassComponents.Assemble` is public and carries basis restart; no row names it.

### 6c. Statements found false along the way

| # | Where | Statement | What holds |
|---|---|---|---|
| Fa1 | `port-fitting-brief.md:503` | "Training fits are a cached view only on request (R's `save.yhat`), stored apart from the draws." | No cache is built: `TrainingFit` and `Predict` compute on each call |
| Fa2 | `opus-bass-completion-remainder.md:1` | "Everything R's BASS can do is now built at the maths layer." | `start.temper`, `mcmc.use` on `sobol`, a declared ladder list, conditioning on one of several output variables and the `other` member are not |
| Fa3 | `BassRunOptions.cs:15` | "so tempering and SMC read it too when they join" | Both joined; the sentence has aged |
| Fa4 | `completion-audit-brief.md` "Not in this brief" | "The extensions (E1–E15)" | The checklist holds E1–E18 |
| Fa5 | `port-outputs-brief.md:82` | "condition on one functional coordinate" | Built as: a coordinate naming every output variable (E06) |

The extension seams' write-up, owed at the audit, landed in `docs/Maths` on 2026-11-02 (6b).

## Part 7: what remains, and the limits

### 7a. Optional backends and extensions

| Item | State | Where |
|---|---|---|
| SMC for fitting, stage 8 C4 | Built, optional; tracked apart from capability completion | `BassSmc.cs`, `BassSmcRun.cs`; `oracle/BassSmcParityTests`, `BassSmcTests` |
| SMC for cut calibration, stage 8 C5 | Built, optional | `calibration/CutSmc*.cs`; `BassCalibrationCutSmcTests` |
| Flat-histogram fitting, N17 Phase 2 | Built, optional, beyond the contract | `BassFlatHistogramRun.cs`; `BassFlatHistogramTests` |
| Local-move options, F-M2 Phase 1 | Built, off by default; no arm cleared the pass | `BassConfig.cs:418,422` |
| Seams F-X1, F-X2, F-X3, F-X6, N-X1, N-X2, N-X3, O-X4 | Built; no consumer | `extensions-preparation-brief.md` Reports |
| E1–E18 | Not built; deferred, parked or candidate | `extensions-brief.md`, `extensions-checklist.md` |
| Rendering and text formatting above maths | Not built; waits on a consumer and viz v2's contracts | The hub's P09–P15 rows |
| Efficiency pass; prefix sums first | Priced, not built | `performance-cost-brief.md` Reports |

### 7b. Measured limitations, each where it was measured

These are records this audit read, not reruns.

| Limit | Reading | Measured in |
|---|---|---|
| Friedman's posterior, MCMC | Chains hold their own bands of M; births and deaths accept about 1.5%; R̂ of M 1.5–2.5 at 5,000 + 2,000 sweeps; every D10 group fails | S10 and F-M2's pass, `port-fitting-checklist.md` Reports; `BassVerdictTests.S10_Fits_FailEveryGroup` reads the verdict off the pass's recorded R̂, not a new fit |
| Friedman's posterior, the reference | Computed by the flat histogram on three seeds; a1 about 66–70%, ln Ẑ about −841 | `flat-histogram-warmstart-brief.md`, `flat-histogram-loose-ends-brief.md` Reports |
| Friedman under SMC | Replicates disagree 13–66 nats and end in 1–2 families at S10's sizes; 5.6 nats low on the linear form at b = 0.10 | `stage8-fitting-smc-brief.md`, `flat-histogram-loose-ends-brief.md` Reports. The band runs are open (`friedman-smc-band-brief.md`) |
| The deterministic profile, F-R1 | No arm clears: predictive R̂ above 1.05 in all 144 fits, least 1.23; ESS about 5 of 8,000 | `deterministic-profile-brief.md` Reports |
| N15's cost | Construction 97–99.8% of a σ² proposal; 1.4 ms at m = 20, k = 39 | `performance-cost-brief.md` Reports |
| The change proposal at the cap | 54–61% of a sweep at n ≈ 1,000, 79–96% from n = 10,000 | `performance-cost-brief.md` Reports |
| Spectral gap | Exponential in N on small fixtures, all of it at the count cap | `spectral-gap-brief.md` |
| Run lengths | R's 9,000 + 1,000 stay provisional; D10 made them a budget | The doc's status section; `performance-statistical-brief.md` not run |
| One training case | Refused | Ed04, a probe |

## The review pass

**Completeness.** Rows checked by hand against R's signatures and return statements and the brief's lists: P01–P15 in 1j and 1h; every formal argument of `bass` (30, in 28 rows), `bassPCA`, `bassBasis`, `predict.bass`, `predict.bassBasis`, `sobol`, `sobolBasis` and `calibrate.bassBasis`, and every returned element; B01–B06; D01–D10; the Corrected entries; each edge case of the hub's list; the benchmark list; the open boxes.

**The re-trace of every tenth Fulfilled row,** from R to the test, ignoring the note: A10, A27, B-r12, B-r23, C11, D-p05, E09, F02, G05, G16, H07, Bb02, Df06, Df16, Dd08, V03, V13, Ed10. The question for each: would this test fail if the capability were removed?

| Row | Result |
|---|---|
| A10 `npart.func` | **Failed.** `EachRefusalReason_HasItsCase` builds `BassAdmission` by hand, and `Resolved_Functional_CountsBothParts` reads the default only. A resolver that ignored the declared `MinimumFunctional` would pass both |
| G16 `discrep.vars` | Stood, on a second test: the row cited a block's conditional. `Outputs_ArePhysical_Reproducible_AndThinned` rebuilds the state from the draw's coefficients, which is what would fail. Citation corrected |
| The other 16 | Stood. Each compares against a reference outside the code, and a broken capability moves the compared value |

**Widening A10 to its neighbours,** the controls a config declares: `maxInt`, `maxInt.func`, `maxInt.cat` (a declared cap below the variables' count was unread), `h1`/`h2` and `w1`, `w2` (the tests read `From` at the defaults; a non-default path ran only in the gated F-R1 harness). `maxBasis`, `npart`, `thin` and `degree` were sound. **Seven rows moved** (A03, A04, A05, A10, A15, A18, A19), through the new probe `DeclaredControls_ReachTheModelAndTheRetainedTerms`, to Fulfilled again on its passing.

**Demotion.** Two verdicts leaned on a record alone and were reworded rather than kept as read evidence: part 5's readings rest on the cost half's report, as the brief instructs, and say so; 7b's limitations are records and say so. `BassVerdictTests.S10_Fits_FailEveryGroup` reads the verdict off the pass's recorded R̂, not a new fit; 7b says so.

**Citations.** A script confirms every cited test name exists (482) and every `file:line` lies inside its file (164). About 80 tests were read whole; for the other 192 a second script kept each one's first nine assertions and comments, and those were read, so a helper's assertions fall outside that reading.
