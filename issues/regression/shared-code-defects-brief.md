# Two defects in shared numerics: the weight reader and the legacy Cholesky factor (chip brief)

Commissioned by the user on 2026-11-13 through the coordinating session. Both defects were found by conformal slices and parked for their owners ([handoff](coordination-handoff-2026-11-11.md#open-defects-in-shared-code)). Each returns wrong numbers with no error, to every caller that reads through it.

**Read first:** [AGENTS.md](../../AGENTS.md), [CONTRIBUTING.md](../../CONTRIBUTING.md), then the two records linked below.

## Slice 1: the weight reader at extreme weights

**The site.** `PredictiveChecks.Rows` (`src/maths/regression/PredictiveChecks.cs`) takes `log(weight / total)`. `Pointwise.CheckWeights` (`src/maths/regression/Pointwise.cs`) admits the weights it then misreads.
- `[1e300, 1e-300]`: the second ratio underflows and its log weight reads −∞. Subtracting logs keeps about −1381.55.
- `[1e308, 1e308]`: the sum overflows, `CheckWeights` accepts it, and both log weights read −∞.

**The record.** Conformal slice 10's report, ["Parked shared weighted-row follow-up" and "Follow-up verification criterion"](archive/bass-gbass-khaos-activespaces/briefs/conformal-layer-brief.md#2026-11-11--slice-10-density-and-cdf-error-diagnostics). The criterion there is this slice's acceptance:
- the first input keeps its finite second log weight; the second returns two weights of 1/2; no NaN and no lost row;
- the normalized log weights, and the weighted predictive, CDF and carried-error readings built on them, agree with an independent high-precision or log-sum-exp reference;
- a finite log-likelihood tilt that makes the first input's tiny row material;
- invariance to a positive rescaling of the weights, and zero-weight controls;
- the shared consumers' selection and the whole assembly pass.

**Its readers.** The predictive checks (two sites in `PredictiveChecks.cs`) and nine conformal sites: `ConformalBayes` (three), `SmoothedPValues` (two), `CrossConformal`, `SplitConformal`, `SplitPitCalibration`, `Tournament`. Search for any reader that normalizes weights on its own; the repair is one shared routine that all of them read.

**The constraint.** Weights in the ordinary range keep their bits, so every pinned hash and fingerprint holds. Scaling by a power of two is exact, which makes this reachable. If no repair keeps the ordinary range bit for bit, stop and report the pins that would move: regenerating a pin is the user's call.

**Gates.** The criterion above; chosen cases and a seeded random sweep over the weights' whole dynamic range, subnormals and values near `double.MaxValue` included; the pinned tier at each commit; the whole regression assembly once at the end.

## Slice 2: `CholeskyDecomposition` on a near-singular covariance

**The site.** `src/maths/linalg/Cholesky.cs`. A pivot at or below 10⁻¹² is set to 10⁻¹², whatever the matrix's scale.

**The record.** Conformal slice 11's report, ["Defects found" and "The shared factor's floor"](archive/bass-gbass-khaos-activespaces/briefs/conformal-layer-brief.md#2026-11-11--slice-11-bgps-drawn-f-the-retained-parameters-content-the-full-step-takes).
- On a smooth kernel's covariance, whose spectrum falls below rounding, the rows under a floored pivot divide their rounding error by 10⁻⁶. Later pivots went to −3.7·10¹⁵⁷ and were floored again; 55 entries of the factor were infinite.
- `LogDet` stayed finite (−683.6). Only `Status.Finite` reported the failure.
- BGP left for `DenseCholesky` (`5db2e40`), and its drawn function carries its own diagonally pivoted factor (`BgpDrawnFunction`).

**Remaining callers.** GMM: `GaussianComponent` and `BhattacharyyaCoefficient` (two sites), in `src/clustering/statistical/gmm/`. Viz-core: `src/viz-core/serializer.cs`. Tests and cost profiles under `tests/linalg/` and `tests/regression/` also construct it.

**Reading first.** Reproduce the failure in a test on that class of covariance. Then measure, on a ladder of condition numbers and scales, where the floored factor's solves and `LogDet` leave a reference's tolerance while every entry is still finite. The overflow is the visible end of the defect; wrong finite numbers come before it. Read how `DenseCholesky` treats the same inputs (`Status`, `FailedAt`, a NaN `LogDet`, `RelativePivot`), and what each remaining caller does with `Status` today.

**Build.** A factor with a non-finite entry is a failed factorization, and the class refuses it the way `DenseCholesky` does: no finite `LogDet`, and no solve, sample, inverse or factor read from it. Each remaining caller reads the outcome and handles it. The two factors then fail one way.

**For the user, reported with the measurements and a recommendation, and left unbuilt.** Each changes GMM's numbers or removes a surface:
- *The floor.* An absolute 10⁻¹² against a floor relative to the matrix's scale, and whether a pivot negative beyond rounding is floored or reported.
- *A pivoted factor in the linear algebra.* `BgpDrawnFunction`'s diagonal pivoting (Harbrecht, Peters and Schneider 2012) as a shared primitive that near-singular callers take.
- *Retiring `CholeskyDecomposition` into `DenseCholesky`.* What is lost (the soft floor, `WriteInverseTo`'s L⁻¹, `Sample`), and which callers, tests and cost profiles move. A removal is its own question.

**Gates.** The reproduction fails before the build and passes after. Every existing test of the class holds bit for bit on inputs the floor never reached. The pinned tier at each commit; the linear algebra, clustering and viz-core test projects that reach the class, and the regression assembly, once at the end.

## Slice 6: the samplers' condition check as a setting

Commissioned by the user on 2026-11-14 through the coordinating session, on the [measurement](#2026-11-14--slice-5s-decision-1-measured-the-condition-estimates-cost-per-sampler-step) (`7dca8ff`). It feeds the BAPS follow-up's [admission input](bass/extensions-brief.md#input-for-the-follow-up-admission-and-singular-systems-2026-11-14).

**The decision.** How often a sampler estimates κ(E) becomes a setting of three values. The default is what slice 5 landed, so no default changes in this slice.
- *Every factorization* (the default): slice 5's behaviour, bit for bit.
- *Screened:* one estimate per run at a reference point, a scalar bound at each step, and the estimate itself at any step the bound does not clear.
- *Reference alone:* one estimate per run at the reference point, and positive definiteness at each step. It is the opt-out for a model already known to be sound.

The measurement's options 2 and 4 are not built. Option 2 costs what *reference alone* costs per step and misses an unidentified design. Option 4's k has no value a caller can choose on principle, and its refusal can arrive late or never.

**The bound.** Each sampler's system is A(c) = Σₖ cₖSₖ. Every Sₖ is positive semidefinite and fixed by the design, and every cₖ > 0 is drawn or fixed at 1: ZᵀZ + λP, ZᵀZ + Σⱼ λⱼPⱼ, and ZᵀZ + Σᵢ λᵢdᵢdᵢᵀ for the locally adaptive penalty. With ρₖ = cₖ / c_ref,ₖ, every Rayleigh quotient of A(c) lies between min ρₖ and max ρₖ times A(c_ref)'s, so

κ₂(A(c)) ≤ κ₂(A(c_ref)) · maxₖ ρₖ / minₖ ρₖ.

- *Its standing.* Derived by the coordinating session on 2026-11-14 and checked against no run. Confirm the form at each sampler site, the hierarchical model's five included. A site whose system does not take the form is reported, not forced into it.
- *From the bound to the gate's reading.* The gate reads a 1-norm estimate on the equilibrated matrix, and the estimate errs low. Carrying the bound to that reading costs three factors: Jacobi scaling against the best diagonal scaling (van der Sluis; Higham, *Accuracy and Stability of Numerical Algorithms*, §7.3), the 1-norm against the 2-norm, and the estimator's shortfall. The screen admits a step when S · κ̂_ref · max ρ / min ρ ≤ `MaxCondition`, where S is a named constant covering the three. Derive S, state the derivation, and let the sweep below test it.
- *The reference point* is fixed by the configuration alone: trace balance, λ₀ = tr ZᵀZ / tr P for each penalty, which the tensor optimizer's box already centres on. State the choice for the per-row multipliers and the hierarchical terms. A resumed run then decides every step as the uninterrupted one did.

**The invariant.** *Screened* decides every step as *every factorization* does. It admits by the bound only what the estimate would admit, and it estimates otherwise, so the two differ in cost alone. *Every factorization* is the reference the tests compare against.

**Build.**
- *The setting,* declared on the sampler's configuration. Propose its name and its values' names in the report, by [CONTRIBUTING's vocabulary](../../CONTRIBUTING.md). It reaches the three `Baps` sites, the five `HierarchicalPSpline` sites and the collapsed scheme's solve (`PenalizedSpline.SolveFactored`, the 2.00 case). A deterministic solve outside a sampler keeps its estimate at every factorization.
- *The refusal's message.* Under the screened rule it states the reference estimate, the ratio and the step's estimate, so a reader sees whether the design or the drawn coefficients failed.
- *Records on the run's result,* in every mode: the reference estimate where the mode takes one, the widest ratio visited, how many estimates ran and the largest estimate seen. They are reductions of what the steps saw and leave the draws untouched. They ride in the chain's state, so a resumed run reports what the uninterrupted one does.
- *The measurement's option 5,* as its own commit, since it changes LinAlg's shared factor: the estimate reuses workspace vectors, and `BandCholesky`'s ‖E‖₁ pass runs only for a factorization whose estimate is asked. If the factor's interface does not take that cleanly, report it and leave it.

**For the user, reported and left unbuilt.**
- *Whether "screened" becomes the default.* Reported with the sweep's slack, the cost table and how often the screen fell through. A default change is the user's.
- *S,* if its derivation and the sweep disagree.
- *Any site outside the form.*

**Gates.**
- *The default mode:* every pinned hash, fingerprint and existing test holds bit for bit.
- *Same chain:* the three modes draw the same chain bit for bit on every admitted run, over the measurement's 19 cases at least.
- *The invariant,* on chosen cases and a seeded random sweep, for each sampler form:
  - designs from well identified to nearly unidentified (y's spread shrinking in the tensor design; x confined toward r distinct values in 1-D);
  - coefficient ratios over the whole range, up to the refusal and past it;
  - no step where the screen admits and the estimate refuses;
  - the slack seen, reported: the largest κ̂(E(c)) / (κ̂_ref · max ρ / min ρ).
- *Refusals:*
  - all three modes refuse an unidentified design before the first draw;
  - *every factorization* and *screened* refuse a coefficient drawn far out at the same step;
  - *reference alone* admits that step, which is its stated loss.
- *Resume:* an interrupted and resumed run matches the uninterrupted one bit for bit in each mode, the records included.
- *Cost:* the measurement's probe (`BapsGateCostProbe`) rerun with the three modes and positive definiteness alone as arms, under the same blocks and load record, with the screen's fall-through count for each case.
- The pinned tier at each commit; the test projects the changed files select classes in; the regression assembly once at the end.

**Files.** `src/maths/regression/spline/baps/**`, `src/maths/linalg/BandCholesky.cs` for option 5, and their tests.

## Slice 7: "screened" as the samplers' default

Decided by the user on 2026-11-14, on [slice 6's report](#2026-11-14--slice-6-the-samplers-condition-check-as-a-setting) and its [verification](#2026-11-14--slice-6-verified-from-clean-archives-coordinating-session).

**The change.**
- `ConditionCheck`'s default becomes `Screened` on `Baps` and `HierarchicalPSpline`. `EveryFactorization` and `ReferenceOnly` stay as named modes, and `EveryFactorization` stays the reference the tests compare against.
- The screen takes the sharper constant, (1 + √ν)/2 in place of √ν (`cond1_unitScale_le_sharp`, verified at `dd57650`).
- Both land together and are verified in one cycle.

**First, before any code moves: checkpoints and pins.** Today only a non-default mode enters a run's fingerprint. Read what the default change does to:
- a checkpoint written under the old default and resumed after the change;
- every pinned hash, fingerprint and golden value a default run reaches;
- any test that reads a default run's refusal text or its `ConditionRecord`.

The shape to find moves no pinned value and leaves no existing checkpoint without a way to resume. The two modes draw the same chain and refuse the same steps, which is what makes such a shape possible. If none meets both, stop and report the candidates with what each one moves: regenerating a pin, or orphaning a checkpoint, is the user's call and its own question.

**What a default run now shows differently.** State each in the report:
- the refusal's text, where a default run refuses;
- `ConditionRecord` on a healthy run: no estimates and no largest estimate, with the reference condition and the widest ratio read in their place;
- one exact reading per run, ν solves, before the first draw.

**Gates.**
- A default run draws the chain it drew before, bit for bit, over the measurement's 19 cases and every pinned draw.
- The invariant sweep again under the sharper constant: no step where the screen admits and the estimate refuses, with the slack table.
- Resume in each mode, and across the change as the first section settles it.
- Every test project the changed files select classes in, named in the report with its count. `Maths.Samplers.Tests` was missed in slice 6's gates.
- The pinned tier at each commit; the regression assembly once at the end.
- No cost probe: the constant moves fall-through only, and slice 6's probe saw none.

**Files.** `src/maths/regression/spline/baps/**` and their tests.

## Slice 8: three problems in the test tooling, found by BPPR's gates

Added by the user on 2026-11-14, from the BPPR chip's [slice 0 report](bass/bppr-brief.md#slice-0-sources-references-and-the-power-spherical-law-2026-11-14), under "Not E36's, found by its gates". Three parts, one commit per concern, in this order. Each starts with a reading of what is true today, and the report states that reading before any change.

**8a. Test runs keep the shared-machine rule.**
- *The problem.* `scripts/test-select.ps1` lowers its own priority while it starts the harness (lines 346 to 350), so the harness is below normal. The BPPR chip read the chain during a regression run: testhost (Normal) ← dotnet (Normal) ← dotnet (Normal) ← the harness's dotnet (BelowNormal). So the processes that do the work run at normal priority, in every test-select run and every `scripts/verify-csharp.ps1` verification so far.
- *First:* reproduce the reading, and find which process leaves below normal and why. The harness starts `dotnet test` from `src/test-harness/` (`FactRunner.cs`, `Program.cs`, `TestDiscovery.cs`). MSBuild has a low-priority switch; that is cited from memory, to be read before it is relied on.
- *The fix to find:* a run's whole tree at below normal, build nodes and test host included, by the tools' own switches where they have them.
- *Keep:* a block that times itself sets its own priority (High, pinned) and never runs below normal; that rule is the user's. Read which runs time blocks and how they set priority, and say in the report how each stands after the fix.
- *Gates:*
  - the chain's priority classes read during a run, before and after, through test-select and through verify-csharp;
  - the thread and processor caps unchanged;
  - one project's selection and counts the same before and after;
  - `TestHarness.Runner.Tests` and every project the changed files select.

**8b. The transport oracle's time against its limit.**
- *The problem.* `tests/oracle/ROracle.cs` gives every R oracle 120 s (line 74). `oracles/transport_oracle.R`, under `DiagramMetricsTransportParityTests.DiagramDistances_MatchT4TransportAndLpAssign` at p = 1 and 2, takes about 115 s alone, so any load fails both. They failed in the clean-archive Oracle runs of slice 6 and of BPPR's slice 0.
- *First:* where the 115 s go: R's start and its package loads, or the transport problem itself. Time each alone.
- *Then,* by what the reading shows: remove time that is not the reference's (a load repeated, a problem solved twice), or give this oracle a stated limit of its own with the reading behind it.
- *Keep:* the fixture, the reference values and the assertion's tolerance. A smaller or easier case is a removal, and the user's question.
- *Gates:* the two cases pass alone and inside a full Oracle run under the load of the day, each time recorded; every other oracle's limit unchanged; Oracle 207 of 207 at commit scope.

**8c. The Distributions case that runs 15 to 22 minutes.**
- *The problem.* `GeneralizedBetaPrimeTests.TheIntegratedDensity_MeetsTheGeneralizedHyperbolicsClosedForm` at (−0.5, 0.09, 1.0, 1.2, 0.5) ran 22 min at commit scope and 15 min in the whole assembly in BPPR's gates. Its comment says its table took 20 s before the integral about its peak, so it should now take less. `Maths.Distributions.Tests` has no `slow-tests.txt`, so every scope runs it.
- *First:* time the case alone.
  - Seconds alone: the cause is contention, of which 8a may be part, and the report says so.
  - Minutes alone: find where the time goes and since which commit. A path that no longer takes the integral about its peak is a defect to fix, with the test as it is.
- *The slow tier.* If a method still reaches the slow tier's threshold after that, write Distributions' list by the tool's own route, `-WriteSlowList` from a profile of the whole assembly ([CONTRIBUTING](../../CONTRIBUTING.md#testing-standards)). Report what it lists, since commit scope then leaves those out.
- *Gates:* the case's time alone, before and after; the assertion and its tolerance unchanged; Distributions whole; and, if code under `src/` moves, every project the changed files select, since `Maths.Distributions` is a hub.

**Times under load.** The machine is shared and loaded. Read a time alone, at least twice and apart, with the load of the moment beside it. No limit is set from one reading.

**Files.** `scripts/test-select.ps1`, `scripts/verify-csharp.ps1`, `src/test-harness/**`, `tests/oracle/ROracle.cs`, `r/oracles/transport_oracle.R`, the code under the generalized hyperbolic's integrated density, `tests/Maths.Distributions.Tests/slow-tests.txt`, and their tests.
- The BPPR chip (`local_4d3c75b2`) writes in `src/maths/distributions/` (`GammaRatio.cs`, `PowerSpherical.cs`, `SpecialFunctions.cs`) and, from its slice 1, in BASS's core.
- Both coordinating sessions run `verify-csharp.ps1` and `test-select.ps1`. Message the second coordinating session before the commit that changes either script, so that no verification is mid-run.

## Slice 9: `LogBeta` through `LogGammaRatio`

Decided by the user on 2026-11-14, as a distributions change of its own. It comes from the BPPR chip's slice 0 follow-up and the user's own check of the current functions.

**The defect.** `SpecialFunctions.LogBeta(a, b)` is `LogGamma(a) + LogGamma(b) − LogGamma(a + b)` (`SpecialFunctions.cs`, line 80). Where one shape is large, two large ln Γ cancel. The user's check, against ln B(x, 1) = −ln x:
- at x = 10⁵ the absolute error falls from 7.75·10⁻¹¹ to 1.78·10⁻¹⁵ through `LogGammaRatio`;
- at x = 10¹⁶ `LogBeta` returns 0, and the ratio form returns −36.841361487904734, which is −`Math.Log(x)`.

**The change.** ln B(a, b) = ln Γ(s) − `LogGammaRatio`(l, s), with l the larger shape and s the smaller (`GammaRatio.cs`, as fixed in `117230e`).

**The order of work,** the user's:
1. *Qualify first:* balanced and unequal shapes against a high-precision reference from outside the engine, far ranges included, by chosen cases and a seeded sweep. State where each form holds its digits. If the ratio form loses digits in some regime, report it with the form that would hold them before building further.
2. *Then run the existing pins,* before any expectation is updated.
3. *Report which pins move,* and for each whether a density value changed or an accept/reject decision changed. BARS's locality kernel reads the beta log density in its Hastings correction (`bars/KnotKernels.cs`, `LogDensity`), so a last-bit change can reach a trajectory. The count of affected pins is not known.
4. *Re-pin where the numerical change explains the difference,* and only there. A pin whose movement the change does not explain is a finding: report it and leave it.

**Read as well.**
- `LogBeta`'s callers under `src/` (`Beta.cs`, `GeneralizedBetaPrime.cs`, `StudentT.cs` and `SpecialFunctions.cs` itself) and what reads them: beta densities and quantiles, the incomplete beta, the conformal layer's Beta levels.
- A run resumed across the change: say in the report what a checkpoint written before it does after it.

**Gates.**
- The qualification's table: the worst error by regime, before and after, against the outside reference.
- Every test project the changed files select, named with its count, since `SpecialFunctions` is a hub. The pinned tier at each commit; the whole Distributions and Regression assemblies at the end.
- The list of pins that moved, each with its cause, and a statement of those that did not.

**Sequence.** After 8a, so that its gate runs keep the machine rule. BPPR's slice 2 does not wait on it.
- The BPPR chip (`local_4d3c75b2`) owns `GammaRatio.cs`; a change inside `LogGammaRatio` is agreed with it by message first.
- The first coordinating session's chips hold pins that read beta laws (the conformal layer). Tell the second coordinating session which of those move before updating any.

**Files.** `src/maths/distributions/SpecialFunctions.cs` (`LogBeta`), its tests, and the pinned expectations that step 4 allows.

## Working rules

- **Order.** Slice 1, then slice 2. They touch different code, so slice 2's reading starts when slice 1 lands. The coordinating session verifies each landed slice from a clean archive ([scripts/verify-csharp.ps1](../../scripts/verify-csharp.ps1)).
- **Files.** This chip owns `PredictiveChecks.cs`, the weight check in `Pointwise.cs`, `Cholesky.cs`, the three caller files, and their tests. Three other chips are at work:
  - the conformal chip (`local_69bae1e5`), on the remaining slices 14 to 18, owns `src/maths/regression/conformal/**` and reads `Rows`; a change inside that folder is agreed with it by message first;
  - the GBASS chip owns BASS's sampler files for the u–s funnel fix;
  - two Lean chips work in `lean/`.
- **Git.** On `main`, one commit per concern with its tests, committed by pathspec, nothing left staged, nothing pushed, no attribution lines.
- **Builds and runs.** `ArtifactsRoot` set inline to `artifacts\chip-defects\`, deleted at the end. Heavy runs at below-normal priority with capped threads; the machine is shared. Tests by [targeted runs](../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`).
- **Review.** Each slice ends with a fresh-eyes pass for correctness and coverage gaps, fixed before the report.
- **Scope.** New beside old. Findings are framed by problem class. A follow-up is reported with the decision it would change and what each outcome would do. `docs/` stays untouched; the report lists "Doc changes owed". Sub-agents need the user's permission.
- **Notes.** An "In progress" entry under Reports holds uncommitted work, running checks and the next step. Probe code and data go to `artifacts/`.
- **Events.** End a turn with one line and the report's link when a slice lands or needs a decision. Results go in the report.

## Reports

Newest first. Each gives the commits, the gates' results, defects found, decisions for the user and "Doc changes owed".

### 2026-11-14 — Slice 8a verified (coordinating session)

- **Commit:** `0dad092`, `scripts/test-select.ps1` alone: `-lowPriority` on the `dotnet test` it starts, and two comments.
- **Read apart from the chip's reading:** a LinAlg run through `scripts/verify-csharp.ps1` on a snapshot of `14a77f9`, with every process in the run's tree sampled each 1.5 s for its priority class.
  - Every process was BelowNormal in every sample: `dotnet test -lowPriority`, the worker node with `/low:true`, csc under both, `vstest.console`, testhost, and their console hosts.
  - 116 of 116, as the chip read before and after.
- **Tests:** this entry covers the diff and the chain. The chip's report names what the changed file selects.
- **Open, for the chip's report:** `scripts/fact-harness.ps1` starts `dotnet test` through `src/test-harness`, a separate entry point. Its chain has not been read under a project that starts a worker node.
- **What follows from it.** Every test-select run and every clean-archive verification before `0dad092` did its building and testing at normal priority. Their results stand; their load on the machine was above what the rule allows.
- **Receipts:** [verify-defects-s8a-0dad092](receipts/verify-defects-s8a-0dad092/), with the sampled chain.

### 2026-11-14 — Slice 8: in progress

**8a, the reading before any change.** Every process under the run, read from `Win32_Process` every 0.5 s, from `2ded631` staged.
- *Where the chain leaves below normal: MSBuild's worker node.*
  - test-select starts `dotnet test` below normal, and its start holds.
  - MSBuild, inside that process, starts its worker node with `/nodemode:1 /nodeReuse:false /low:false`. The build was not asked for low priority, so the node sets itself to Normal.
  - Everything the node starts inherits Normal: csc for the projects it builds, `vstest.console`, and testhost.
- *Through test-select,* LinAlg (two project references), Slice scope: dotnet test BelowNormal; node Normal; csc Normal under the node (BelowNormal in the main node); vstest.console Normal; testhost Normal. 193/193.
- *Through verify-csharp,* LinAlg, Commit scope over `BandCholesky.cs`: the same chain. 116/116.
- *With no worker node,* in `TestHarness.Runner.Tests` (one reference), MSBuild started no node and the whole tree stayed BelowNormal, 19/19. So the defect shows wherever MSBuild starts a node, as it does for the regression project.
- *The fact harness is not in the chain.* test-select starts `dotnet test` itself; `src/test-harness` is the runner behind `scripts/fact-harness.ps1`, a separate entry point.
- *MSBuild's switch, read:* `dotnet msbuild -h` lists `-lowPriority[:True|False]`, short `-low`.
- *Trial in a staged copy,* with `-lowPriority` on test-select's `dotnet test`: the node gets `/low:true`, and the whole tree stays BelowNormal. 193/193.

**Next.** 8a's edit waits for the coordinating session to confirm that no verification is mid-run. 8b's timings are running.

### 2026-11-14 — Slice 7 verified from a clean archive (coordinating session)

- **Run:** `scripts/verify-csharp.ps1` at HEAD `1441600`, the Commit selection over `1441600`'s four files.
- **Result:** regression 228 of 228 (25 of 325 classes and the pinned tier). A sweep of all 27 test projects with `test-select -ListOnly` found no other project selecting a class.
- **Read:** the source diff. The default flips on both constructors, the screen's factor is m·(1 + √ν)/2 as `cond1_unitScale_le_sharp` states it, and the fingerprint rule's code is untouched.
- **Slice 6's sweep, closed.** That sweep covered the 22 projects under folders named `*.Tests`. The five it missed are smoke and analysis programs with no test class, so nothing was skipped.
- **Not run here:** the whole regression assembly, 3,458 of 3,458 on the chip's own run.
- **Slice 6's open question, settled since.** Slice 6's report left open whether a smaller factor holds for every family. It does not.
  - The Lean chip proved that no constant bounds the inverse part over every family (`inverse_part_unbounded`, `fc793df`, verified at `95c3500`), on a dense family.
  - A tridiagonal family grows like √n as well, once max ρ / min ρ passes the pattern's spectral gap: 0.90, 2.07 and 3.09 at n = 256 for ratios of 10², 10⁴ and 10⁶. That one is numeric, rebuilt by the coordinating session in 50-digit arithmetic. The banded case is proved in Lean since (`inverse_part_unbounded_tridiagonal`, verified at `6f6bb27`), on a path family with three nonzeros a column.
  - So the screen's m·(1 + √ν) is the right order, and the gap between the proved factor and the sweep's slack is how far BAPS's own systems sit from the worst case.
- **Receipts:** `artifacts/verify-defects-s7-1441600/`, with the counters, snapshot, changed files and the sweep tracked in [receipts](receipts/verify-defects-s7-1441600/).

### 2026-11-14 — Slice 7: "screened" as the samplers' default

**Commits.**
- `740e081`: the In-progress entry.
- `1441600`: the default and the sharper constant.
- This report.

**Checkpoints and pins, read before any code.** The shape found moves no pinned value and leaves no checkpoint without a way to resume.
- *Checkpoints.* A BAPS checkpoint (`EnsembleCheckpoint`) is an in-memory object; nothing in `src/` writes one to disk. One taken before the change cannot reach a process running after it.
- *Today's fingerprint rule stays.* A mode enters the fingerprint unless it is `EveryFactorization`.
  - A default run's fingerprint now carries `Screened`.
  - A checkpoint taken under the old default is refused by the new default and resumes when `conditionCheck: EveryFactorization` is named.
  - A checkpoint taken under `Screened` resumes under the default, slice 6's included.
  - So each checkpoint resumes only under the mode its records came from. `BapsConditionCheckTests` checks all three cases.
- *Pins.* No test in the pinned tier reaches `Baps` or `HierarchicalPSpline`. The draws a default run makes do not move.
- *Callers.* None outside `spline/baps/`. The existing BAPS tests now run screened and pass as before.

**What a default run now shows differently.**
- *The refusal's text.*
  - A step that is not positive definite reads "not positive definite".
  - A step refused on its estimate reads "its equilibrated condition … past 1e14, where the reference point's is … and the coefficients' ratio to it …".
  - Every factorization keeps the combined text it always gave.
- *`ConditionRecord` on a healthy run.* Estimates 0 and no largest estimate (NaN). In their place: the reference condition, read exactly, and the widest ratio, recorded as before.
- *Before the first draw:* one exact reading per run and system. That is one factorization and ν solves (dense for the hierarchical collapsed precision).
- *The run's fingerprint* carries the mode, as above.

**The sharper constant.**
- *The new constant.* The screen takes (1 + √ν)/2 in place of √ν. Its basis is `cond1_unitScale_le_sharp`, whose column step ‖P‖₁ ≤ λ_max(P)·(1 + √ν)/2 for a positive definite P is attained. The Lean chip proved it in `7f535b3`; the coordinating session verified it at `dd57650` and `6fc199f`.
- *The screen now* clears a step when 2·m·((1 + √ν)/2)·‖E(c_ref)⁻¹‖₁·max ρ / min ρ ≤ 10¹⁴. So S = m·(1 + √ν) / ‖E(c_ref)‖₁.
- *The sweep, rerun.* On no step did the screen admit what the estimate refused.

  | System | Systems | Screen cleared | Fell through, admitted | Refused | Not PD | Slack, admitted | Slack, all PD | S |
  |---|---|---|---|---|---|---|---|---|
  | ZᵀZ + λP | 848 | 291 | 220 | 192 | 145 | 1.00 | 3.49 | 8.46–33.2 |
  | ZᵀZ + Dᵀdiag(λᵢ)D | 816 | 184 | 195 | 193 | 244 | 1.00 | 1.10 | 11.9–25 |
  | ZᵀZ + Σⱼ λⱼPⱼ | 474 | 80 | 91 | 126 | 177 | 1.00 | 15 | 76.1–189 |
  | ZᵀZ/σ² + I/τ_u² | 124 | 72 | 38 | 13 | 1 | 1.30 | 1.30 | 10.7–25.6 |
  | (R/τ_u²)I + P/τ_0² | 124 | 78 | 30 | 5 | 11 | 0.94 | 1.00 | 15.4 |
  | the collapsed precision | 332 | 119 | 104 | 65 | 44 | 1.03 | 1.07 | 47.6–89.5 |
  | ξ²ZᵀZ/σ² + I/σ_η² | 128 | 72 | 42 | 13 | 1 | 1.30 | 1.30 | 10.7–25.6 |
  | R·ZᵀZ/σ² + P/τ_0² | 124 | 41 | 43 | 21 | 19 | 0.99 | 1.41 | 11.9–25 |

  - S now runs from 8.5 to 189, against 12.5 to 335 before. The slack on admitted systems is unchanged, at most 1.3.
  - A few more systems clear: 291 against 283 for ZᵀZ + λP, 184 against 176 for the local form.

**The gates.**
- *A default run draws the chain it drew before.* The 19 cases draw the same chain under `Screened`, the default now, as under `EveryFactorization`, the default before. A new test checks a default run directly against every factorization on `Baps` and `HierarchicalPSpline`.
- *Every pinned draw:* the pinned tier ran in each Commit-scope selection below.
- *`BapsConditionCheckTests`:* 28/28.
- *Resume:*
  - in each mode, as in slice 6;
  - across the change, the three cases above.
- *Commit scope, every test project found on disk (27)*, from `740e081` plus the changed files:
  - `Maths.Regression.Tests`: 25 of 323 classes and the pinned tier, 228/228.
  - `Maths.Oracle.Tests`: 6 affected files, no class selected.
  - No class selected, so nothing ran:
    - `Maths.Samplers.Tests` (0 of 31), `Maths.LinAlg.Tests`, `Maths.Distributions.Tests`, `Maths.Geometry.Tests`, `Maths.Distance.Tests`, `Maths.Topology.Tests`;
    - `VizCore.Tests`, `VizCoreSmoke`;
    - `TDA.Ph.Tests`, `TDA.Mapper.Tests`, `TDA.DimReduction.Tests`;
    - `Clustering.Graphical.SPC.Tests`, `Clustering.Graphical.CumulativeField.Tests`, `Clustering.Geometric.EMeans.Tests`, `Clustering.Geometric.KMeans.Tests`, `Clustering.HdbScanSmoke`;
    - `Spc.BlattSmoke`, `Spc.BlattAnalyze`, `Spc.HyperbolicSmoke`;
    - `Graphs.Primitives.Tests`, `Synthetic.Tests`, `Archivory.Tests`, `RepoAudit.Tests`, `StorageAudit.Tests`, `TestHarness.Runner.Tests`.

- *The regression assembly at the end:* Slice scope, the whole assembly on `1441600`: 3458/3458.
- *No cost probe.* The constant moves fall-through alone, and slice 6's probe saw none.

**Doc changes owed.**
- `docs/` for BAPS:
  - `Screened` as the default;
  - the sharper constant and its Lean theorem;
  - the refusal text and records a default run now shows;
  - the fingerprint rule and how to resume a checkpoint taken under the old default.
- This supersedes slice 6's "the default's fingerprint string is unchanged".

### 2026-11-14 — Slice 6 verified from clean archives (coordinating session)

- **Runs:** `scripts/verify-csharp.ps1`, the Commit selection over each commit's files in every test project where they select a class. All 22 test projects were swept with `test-select -ListOnly` first.
- **`5f7fac5`, at HEAD `5f7fac5`:** regression 232 of 232, the one project that selects (26 of 322 classes).
- **`514dbcd`, at HEAD `514dbcd`:**
  - linear algebra 129 of 129 (12 of 21 classes);
  - Samplers 8 of 8 (1 of 31), a project the chip's own gates did not name;
  - VizCore 12 of 12 (2 of 85);
  - regression 3,116 of 3,116 (315 of 322);
  - Oracle 145 of 147 (28 of 39). The two failures are `DiagramMetricsTransportParityTests` at p = 1 and p = 2, each `transport_oracle.R` timing out at 120 s. That is the standing finding for the TDA owner, and no changed code reaches it.
- **`8b0a7fb`, at HEAD `1f2f23e`:** its diff changes comment lines only; linear algebra 129 of 129.
- **Read:** `SamplerAdmission.cs` whole, the reference construction and the three step sites in `Baps.cs`, and `514dbcd`'s change to `BandCholesky` and `ConditionEstimate`.
  - The fixed ZᵀZ coefficient enters the ratios at each `Baps` site.
  - The screen's test is the chain the Lean module proves (`cond1_unitScale_le`, verified at `94f6b37`).
  - `condition` defaults to true, so every caller outside BAPS keeps its bits.
  - The five hierarchical sites were not read line by line; the sweep's hierarchical test and the chain-by-chain check cover them.
- **Not run here:** the whole regression assembly (the chip's Slice scope, 3,453 of 3,453) and the cost probe, whose table stands on the chip's run.
- **Receipts:** `artifacts/verify-defects-s6a-5f7fac5/`, `-s6b-514dbcd/` and `-s6c-8b0a7fb/`, with the counters, snapshots, changed files and the sweep tracked in [receipts](receipts/).

### 2026-11-14 — Slice 6: the samplers' condition check as a setting

**Commits.**
- `76851ee`: the In-progress entry.
- `5f7fac5`: the setting, the screen and the records.
- `514dbcd`: option 5.
- `8b0a7fb`: the doc-comment fix from the review of `514dbcd`.
- This report.

**The setting.** The name, the values' names and the records' type are proposed by [CONTRIBUTING's vocabulary](../../CONTRIBUTING.md#naming--vocabulary): functional names, with *screen* used in its numerical sense, a cheap test that rules out an expensive one.
- *`ConditionCheck`* is a constructor parameter, `conditionCheck`, on `Baps` and `HierarchicalPSpline`. Its three values:
  - `EveryFactorization` is the default, slice 5's rule at each step.
  - `Screened`.
  - `ReferenceOnly`.
- *Where it reaches:* the three `Baps` sites, the hierarchical model's five sites, and the collapsed scheme's `PenalizedSpline.SolveFactored`, through an internal overload. Deterministic solves keep their estimate at every factorization: `SolveFactored(λ)`, `PenaltySpectrum` and `AnisotropicTensorSpline`.
- *`ConditionRecord`*, one per system, sits on `BapsResult.ConditionRecords` and `HierarchicalResult.ConditionRecords`. It is filled in every mode with:
  - the reference condition (NaN under `EveryFactorization`, which takes none, and +∞ for a reference that is not positive definite);
  - the widest ratio visited;
  - the number of factorizations;
  - the number of estimates;
  - the largest estimate (NaN when none ran).
- *Where the records live.* They ride in each chain's state (`ConditionTally`), so a snapshot carries them.
- *The fingerprint.* A non-default mode enters the run's fingerprint, so a checkpoint resumes only under the mode that started it. The default's fingerprint string is unchanged.
- *The screened refusal* names the step's estimate, the reference point's condition and the coefficients' ratio to it. For example, for λ/λ₀ = 10¹⁴ on a 1-D fixture: "its equilibrated condition 2.39E+14 past 1e14, where the reference point's is 2.72 and the coefficients' ratio to it 1E+14".
- *The every-factorization message* keeps its text exactly.

**The form at each site.** Confirmed at all nine sites, so no site is reported outside the form.
- *Two cases need care.*
  - *ξ²ZᵀZ/σ² + I/σ_η².* The coefficient ξ² can be 0. The bound is then infinite, and the screen hands the step to the estimate.
  - *The collapsed precision.* Its data term R·Zᵀ(σ²I + τ_u²ZZᵀ)⁻¹Z is Σₖ R/(σ² + τ_u²λₖ)·λₖvₖvₖᵀ: ν positive semidefinite terms from ZᵀZ's eigenvectors, which the block already builds. Each term's ratio to the reference, (λ̄ + λₖ)/(σ² + τ_u²λₖ), is monotone in λₖ. So the extremes sit at the smallest positive and the largest eigenvalue, and the ratio costs O(1) per step.
- *What the screen assumes of a penalty.* It assumes every `IBandPenalty` term is positive semidefinite, as a prior precision is. `DifferencePenalty` and `TensorPenalty` are: DᵀD, D_xᵀD_x ⊗ I and I ⊗ D_yᵀD_y.

**The reference point.** It is fixed by the configuration alone: trace balance of each term against its system's leading term. A resumed run therefore decides every step as the uninterrupted one does.
- `Baps`:
  - λ₀ = tr ZᵀZ / tr P;
  - λ₀ⱼ = tr ZᵀZ / tr Pⱼ;
  - the per-row multipliers at their prior mean, ηᵢ = 1, so every local term sits at λ₀.
- `HierarchicalPSpline`:
  - ZᵀZ/σ² + I/τ_u² and ξ²ZᵀZ/σ² + I/σ_η²: I at λ̄ = tr ZᵀZ / ν;
  - (R/τ_u²)I + ΣPⱼ/τ_0j²: each Pⱼ at ν / tr Pⱼ;
  - R·ZᵀZ/σ² + ΣPⱼ/τ_0j²: each Pⱼ at tr ZᵀZ / tr Pⱼ;
  - the collapsed precision: the data term at σ²/τ_u² = λ̄, and each Pⱼ at trace balance against the data term there.

**From the bound to the gate's reading, and S.**
- *The bound* κ₂(A(c)) ≤ κ₂(A(c_ref))·max ρ / min ρ holds as the brief states it.
- *Two steps carry it to the gate's κ₁(E(c)).*
  1. *Jacobi scaling.* diag A(c) is diag A(c_ref) times factors in [min ρ, max ρ]. So λ_min(E(c)) ≥ λ_min(E(c_ref))·min ρ / max ρ, and van der Sluis's factor does not enter.
  2. *The 1-norm.* ‖E(c)‖₁ ≤ m, the most nonzeros in a column, read from the terms' patterns. ‖E(c)⁻¹‖₁ ≤ √ν·‖E(c)⁻¹‖₂. ‖E(c_ref)⁻¹‖₂ ≤ ‖E(c_ref)⁻¹‖₁.
- *The result:* κ₁(E(c)) ≤ m·√ν·(max ρ / min ρ)·‖E(c_ref)⁻¹‖₁.
- *The estimator's shortfall has no bound.* `ConditionEstimate`'s summary says so, and so does Higham's account of the method. A screen built on the reference's estimate would rest on an unproven factor. So the reference is read exactly instead: ‖E(c_ref)⁻¹‖₁ from ν solves, once per run, in both modes that take a reading.
- *This departs from the brief's "one estimate per run".* The cost is ν solves once per run, in place of 5 to 11.
- *What S comes to.* The screen clears a step when 2·m·√ν·‖E(c_ref)⁻¹‖₁·max ρ / min ρ ≤ 10¹⁴, the 2 covering rounding. In the brief's form, S = 2·m·√ν / ‖E(c_ref)‖₁. It depends on ν and the pattern, so it is not one named constant.
- *The screen never refuses on its own.* A reference that cannot clear the screen leaves every step to the estimate. So screened decides every step as every-factorization does by construction, the reference included.
- *What each mode refuses on.*
  - *Every factorization:* the estimate at every step, positive definite and κ̂(E) ≤ 10¹⁴.
  - *Screened:* the same rule at every step. The reference reading only decides where the estimate is skipped.
  - *Reference only:* the exact κ₁(E(c_ref)) at the reference, before the first draw, then positive definiteness at each step. The exact value is at least the estimate. So it can refuse a borderline design that every factorization would admit at that same point: one with κ̂ ≤ 10¹⁴ < κ₁. Every factorization never factors the reference itself, so no run shows that difference unless a chain visits the reference.
- *Settled, not for the user.* The coordinating session checked the chain by hand on 2026-11-14. The derived form and the exact reference reading stand as built.
- *Proved in Lean.* The Lean chip proved the chain in `212752f`: `cond1_unitScale_le` in `lean/Enthymemata/PenalizedAdmission.lean`. The coordinating session verified it at `cea554f`.
  - *Its hypotheses:* A(c_ref) positive definite; m·c_ref,ₖ ≤ cₖ ≤ M·c_ref,ₖ for every k; at most q nonzeros in any column of A(c). All three are satisfiable.
  - *Its conclusion, in exact arithmetic:* κ₁(E(c)) ≤ q·√n·(M/m)·‖E(c_ref)⁻¹‖₁.
  - *Outside the proof:* rounding, which the factor 2 covers, and a coefficient at 0 (ξ²), where the screen hands the step to the estimate.

**The sweep.** `BapsConditionCheckTests`, the four invariant tests. On no step did the screen admit what the estimate refused.
- *What it covers:* each system's form on its own fixtures:
  - designs from well identified to nearly unidentified: x confined to widths from 1 down to 10⁻⁶ in 1-D; y's spread shrinking from 1 to 0 in the 8 × 8 tensor;
  - coefficient ratios from 10⁻¹² to 10²⁰, chosen and seeded at random, past the refusal.
- *The slack against the derived S.* The largest κ̂(E(c)) / (κ_ref·max ρ / min ρ):

  | System | Systems | Screen cleared | Fell through, admitted | Refused | Not PD | Slack, admitted | Slack, all PD | S |
  |---|---|---|---|---|---|---|---|---|
  | ZᵀZ + λP | 848 | 283 | 228 | 192 | 145 | 1.00 | 3.49 | 12.5–58 |
  | ZᵀZ + Dᵀdiag(λᵢ)D | 816 | 176 | 203 | 193 | 244 | 1.00 | 1.10 | 19.8–41.5 |
  | ZᵀZ + Σⱼ λⱼPⱼ | 474 | 79 | 92 | 126 | 177 | 1.00 | 15 | 135–335 |
  | ZᵀZ/σ² + I/τ_u² | 124 | 72 | 38 | 13 | 1 | 1.30 | 1.30 | 17.7–42.5 |
  | (R/τ_u²)I + P/τ_0² | 124 | 78 | 30 | 5 | 11 | 0.94 | 1.00 | 25.6 |
  | the collapsed precision | 332 | 117 | 106 | 65 | 44 | 1.03 | 1.07 | 79.1–149 |
  | ξ²ZᵀZ/σ² + I/σ_η² | 128 | 72 | 42 | 13 | 1 | 1.30 | 1.30 | 17.7–42.5 |
  | R·ZᵀZ/σ² + P/τ_0² | 124 | 38 | 46 | 21 | 19 | 0.99 | 1.41 | 19.8–41.5 |

  - On admitted systems the slack is at most 1.3. Above that it appears only near singularity: up to 15, on tensor systems the gate refuses.
  - S stands 10 to 300 times above the slack seen. The derivation and the sweep disagree by that much. It is recorded below, not put as a decision.
  - "Fell through, admitted" counts the cost of that gap in this sweep: steps the estimate admitted that a tighter S would have cleared. The sweep's ratios reach 10²⁰, far past what the chains visit.

**The gates.**
- *Same chain.* Over the measurement's 19 cases the three modes drew the same chain bit for bit, 20 burn-in steps and 30 draws each:
  - every coefficient;
  - λ, σ², τ² and η, and each component's τⱼ²;
  - the conditional means;
  - every hierarchical draw.
- *The same factorizations.* The records show the three modes factor the same systems and visit the same widest ratio.
- *Refusals.*
  - *An unidentified design.* All three modes refuse it before the first draw: 1-D with one distinct x under Gibbs; the 8 × 8 tensor with fixed y under the anisotropic scheme; the hierarchical collapsed and parameter-expanded sweeps. Reference only refuses at its reference reading, and says so.
  - *A coefficient far out.* At λ/λ₀ = 10¹⁴ the factor is positive definite with κ̂ = 2.4·10¹⁴. Every factorization and screened refuse it, and reference only admits it.
  - *In chains.* The collapsed scheme with a log-step of 12 reaches such λ. Every factorization and screened refuse at the same λ on every seed that refuses. Reference only refuses there only on positive definiteness.
- *Resume.* An interrupted and resumed run matched the uninterrupted one bit for bit in each mode, records included, for `Baps` Gibbs and the centred hierarchical sweep. A screened checkpoint resumed under the default mode is refused.
- *Option 5.* After a factorization the estimate reads the same bits as on a fresh factor, and allocates 0 bytes per call once warm. A factorization with `condition: false` has the same factor, solves and log-determinant bit for bit, and refuses the estimate.
- *Commit and slice scope.*
  - `5f7fac5`: Commit scope 227/227.
  - `514dbcd`: Commit scope from `5f7fac5` plus its files: LinAlg 129/129; Regression 3116/3116, 315 classes; Oracle 145/147, the two failures the transport oracle's R timeout; VizCore 12/12; SPC, CumulativeField, EMeans and KMeans select none. Two of the first readings failed on this chip's own tests and were fixed before the commit: an allocation count taken during the JIT's tier-up, and the far-λ probe asking a skipped pass for its estimate.
  - The regression assembly at the end: Slice scope, the whole assembly on `514dbcd` plus `8b0a7fb`'s file: 3453/3453.

**What the chains saw: R1, on the measurement's fixtures.** In every sampler and both hierarchical sizes, over 50 steps, the screen cleared every step: 0 estimates under screened.
- The reference condition ran from 2.6 to 9.8.
- The widest ratio ran from 6.5 to 2.4·10³.
- The largest estimate, 3.0·10³, is eleven decades below the bound.
- The extensions brief's R1 reading still owes the nearly polynomial truth and the long chains.

**Cost.** The measurement's probe was rerun on `514dbcd` with four arms on the same seeds: the three modes, and positive definiteness alone.
- *How it was run.*
  - Each case first checked that the four arms draw the same chain bit for bit. All 19 did.
  - Then 7 interleaved repetitions of 0.5 s blocks, each arm leading in turn. The blocks ran at High priority with the thread at Highest, pinned to 0xFF, after the measurement's warm-up and settle.
  - 0.9 to 14.8 of 18 logical processors were busy, mostly under 6. The coordinating session held its runs for the duration.
- *The positive-definite arm skips the ‖E‖₁ pass as well.* So the every-factorization column now carries that pass, and runs above `7dca8ff`'s ratios on small systems.
- *Where the files are.* The probe, its stage patch, its run script and the CSVs are in `artifacts/chip-defects/slice6/probe/`.

  | Case | PD alone, µs/step | Every factorization | Screened | Reference only | Screen fell through | Solves per estimate |
  |---|---|---|---|---|---|---|
  | Gibbs ν=12, n/ν=2 | 2.68 | 1.87 | 1.00 | 0.90 | 0/200 | 6.5 |
  | Gibbs ν=12, n/ν=10 | 4.23 | 1.46 | 0.96 | 1.04 | 0/200 | 7.1 |
  | Gibbs ν=12, n/ν=50 | 15.3 | 1.03 | 1.04 | 0.94 | 0/200 | 7.2 |
  | Gibbs ν=48, n/ν=2 | 13.9 | 1.41 | 1.02 | 0.99 | 0/200 | 5.1 |
  | Gibbs ν=48, n/ν=10 | 46.5 | 1.08 | 0.97 | 0.99 | 0/200 | 5.0 |
  | Gibbs ν=48, n/ν=50 | 198 | 1.05 | 1.09 | 1.02 | 0/200 | 5.0 |
  | Gibbs ν=192, n/ν=2 | 130 | 1.17 | 0.94 | 0.94 | 0/200 | 6.7 |
  | Gibbs ν=192, n/ν=10 | 583 | 1.09 | 0.96 | 1.00 | 0/200 | 6.3 |
  | Gibbs ν=192, n/ν=50 | 3,065 | 1.03 | 1.01 | 0.97 | 0/200 | 5.1 |
  | Locally adaptive ν=48 | 44.9 | 1.09 | 0.96 | 1.01 | 0/200 | 5.1 |
  | Collapsed ν=48 | 4.03 | 2.31 | 0.99 | 1.00 | 0/201 | 5.0 |
  | Hierarchical centred ν=24 | 62.6 | 0.99 | 0.95 | 0.96 | 0/400 | 5.8 |
  | Hierarchical collapsed ν=24 | 99.4 | 1.11 | 1.01 | 1.03 | 0/400 | 5.8 |
  | Hierarchical expanded ν=24 | 127 | 1.09 | 1.00 | 1.09 | 0/400 | 6.2 |
  | Hierarchical centred ν=96 | 641 | 1.10 | 1.11 | 1.04 | 0/400 | 6.0 |
  | Hierarchical collapsed ν=96 | 2,634 | 1.21 | 0.98 | 1.05 | 0/400 | 6.4 |
  | Hierarchical expanded ν=96 | 1,594 | 1.02 | 0.93 | 0.99 | 0/400 | 5.8 |
  | Tensor 8 × 8 | 108 | 1.42 | 1.04 | 1.05 | 0/200 | 7.0 |
  | Tensor 16 × 16 | 1,355 | 1.36 | 0.99 | 1.10 | 0/200 | 6.9 |

- *What the table shows.*
  - Screened and reference only cost what positive definiteness alone costs: 0.93 to 1.11, inside the blocks' ranges.
  - Every factorization costs 1.03 to 1.87 on 1-D Gibbs steps, 2.31 in the collapsed scheme, 1.36 to 1.42 on the tensor, and 0.99 to 1.21 hierarchical.
  - The screen fell through on no step, the reference condition ran from 2.7 to 9.8, and the widest ratio reached 1.0·10⁵ (locally adaptive).
- *Allocation, in bytes per step over the positive-definite arm.*
  - *+32 on the band sites under every factorization.* That is the probe's own counting closure, whose display class C# allocates at the method's entry. Without the patch the workspace estimate allocates 0, as `ConditionEstimateTests` checks.
  - *+1,416 in the collapsed scheme.* `SolveFactored` builds a new factor per solve, so its workspace is built per estimate.
  - *+816 and +2,008 for the hierarchical collapsed precision.* Its `DenseCholesky` keeps the allocating estimate.
  - *The screened and reference arms* carry the per-run reference reading spread over a short block: +574 on the 16 × 16 tensor over 395 steps, and +128 at ν = 192, n/ν = 50, over 181 steps.

**Option 5** (`514dbcd`). The factor's interface took it cleanly.
- `BandCholesky.Decompose` and `DecomposeBanded` take `condition: true` by default. False skips the ‖E‖₁ pass, and `EquilibratedCondition` then refuses.
- The estimate runs in a `ConditionWorkspace` the factor keeps: three work vectors and a cached E⁻¹ apply. It computes the same bits.
- The BAPS sampler sites pass the screen's decision, so a step the screen clears skips the pass as well.
- *Not changed here.*
  - `DenseCholesky`, `WideBandCholesky` and `BorderedCholesky` keep their allocating estimate.
  - `PenalizedSpline.SolveFactored` builds a new factor per solve, so the collapsed scheme's estimate still allocates.
  - Holding one factor in the collapsed sampler would take that away. Its accepted solve keeps the factor for the conditional draw, so that change needs its own reading.
- *Reviewed by the coordinating session.* Its one finding: `ConditionWorkspace` had been placed between `ConditionEstimate`'s summary and its declaration. That is fixed in `8b0a7fb`.

**For the user. Reported, left unbuilt.**
1. *Whether screened becomes the default.*
   - *For it.* On the measurement's fixtures a screened step costs 0.93 to 1.11 of the positive-definite arm's, against every-factorization's 0.99 to 2.31, and the screen fell through on no step. It decides as every-factorization does: proved, swept, and checked chain by chain on the 19 cases.
   - *Against it.*
     - Its first call per run costs ν solves, visible only on short runs.
     - The provable S falls through about two decades early, against a run whose ratios reach about 10¹⁴ / (S·κ_ref). The sweep shows that region. No measured chain came near it.
     - The records under screened show fewer estimates, so the largest estimate seen covers only the steps the screen did not clear. The widest ratio is recorded in every mode.

**For the record, not a decision.**
- *S against the sweep.*
  - The derivation gives S = 2·m·√ν / ‖E(c_ref)‖₁, 12.5 to 335 on these systems.
  - The sweep's slack on admitted systems is at most 1.3.
  - So the provable screen falls through one to two and a half decades before a screen tuned to the sweep would. On the measurement's fixtures that gap cost nothing, since no step came near it. `RoundingMargin` is the one named constant; m and √ν are structure.
  - *The Lean chip's own reading agrees.* Script: `lean/Protolemmata/BAPS/penalized-admission-checks.py`, readings PA-7(d).
    - 240 banded families of the three 1-D forms; n from 6 to 29; coefficients over e^±7.
    - κ₁(E(c)) / ((max ρ / min ρ)·κ₁(E(c_ref))) reached at most 1.007, against q·√n up to 38.
    - κ₁(E(c)) reached at most 0.145 of the proved bound.
  - Whether a smaller factor holds for every family is open.
- *Sites outside the form:* none.

**Doc changes owed.**
- `docs/` for BAPS:
  - the `conditionCheck` setting and its three modes;
  - the records on both results;
  - the reference point's choices;
  - the screen's derivation and S;
  - that `ReferenceOnly` admits a drawn λ far out.
- The LinAlg docs:
  - `BandCholesky`'s `condition` flag;
  - the estimate's workspace.
- The extensions brief's input for the follow-up:
  - R1 now rides on every run's records;
  - R2 is answered by `7dca8ff` and this cost table.

### 2026-11-14 — The named-output readouts verified from a clean archive (coordinating session)

- **Run:** `scripts/verify-csharp.ps1` on one clean archive of HEAD `79e858f`, the Commit selection over the seven files of `c543297` in each test project where they select a class.
  - `Maths.Regression.Tests`: 229 of 318 classes, 2,603 of 2,603. `MultiOutputActiveSubspaceTests` is among them.
  - `Maths.Oracle.Tests`: 13 of 39 classes, 57 of 57.
- **The source was read:** nothing found.
  - Scalar and functional fits take the arithmetic they took, and one named output under the identity takes the scalar's.
  - The sampled route reads `BassMultiOutput.Gradient` in bounded draw batches.
  - Pairs within one fit go through `BassMultiOutput.CoActive` alone.
- **What it removes:** `RefuseSeveralOutputs` and its test, which `d9152e6` added, and nothing else of E34's. The three entries now read such a fit, so the first coordinating session's verification of `d9152e6` is superseded by design. Locations, output patterns and a metric that is not K × K are still refused.
- **Receipts:** `artifacts/verify-defects-e8-named-79e858f/`, with the changed-file list, counters and snapshot line tracked in [receipts](receipts/verify-defects-e8-named-79e858f/).

### 2026-11-14 — Slice 5's decision 1, measured: the condition estimate's cost per sampler step

The coordinating session relayed the user's yes on 2026-11-14. This is a measurement only: nothing is built and no default changes. The decision it serves is whether `BapsAdmission`'s estimate stays on every sampler step.

**Method.**
- *The stage.* A `git archive` of HEAD `79e858f` with two probe-only patches, never committed:
  - `BapsAdmission` gains a switch, so one arm admits on positive definiteness alone, and an estimate counter;
  - `ConditionEstimate` counts the solves per call.
- *Where the files are.* The probe, its stage patch, its run script and its CSVs are in `artifacts/shared-code-defects/probe/` (`BapsGateCostProbe.cs`, `baps-gate-cost-{summary,rows}.csv`).
- *Steps per case.*
  - First, the two arms' draws are compared bit for bit: 20 burn-in steps and 30 draws, covering every coefficient, λ, σ², τ², η and component τ². All 19 cases matched.
  - Then the estimates per step and the solves per estimate are counted over 200 steps.
  - Then come 7 interleaved pairs of blocks of about 0.5 s each, with the gate first on even repetitions. The blocks run at process High, thread Highest and affinity 0xFF, after a warm-up of at least 40 runs or 1.5 s and a 300 ms settle.
  - A step's time is the block over the steps (burn-in plus 2 draws, one chain). The table gives the median of 7 and the range.
- *Load.* Over the blocks, 6 to 17 of the 18 logical processors were busy (the last column). Interleaving puts the drift on both arms, and the ranges show it. A ratio near 1.05 or below sits inside them.
- *What the difference is.* Both arms pay the ‖E‖₁ pass, which `BandCholesky` makes before every factorization. So the difference is the estimate's solves and its three work vectors.

| sampler | regime | gate µs/step (range) | PD alone µs/step (range) | ratio | extra bytes/step | estimates/step | solves/estimate (seen) | busy cores |
|---|---|---|---|---|---|---|---|---|
| Baps Gibbs | 1-D, ν = 12, n/ν = 2 | 6.9 (6.3–10.0) | 4.1 (3.5–4.5) | 1.66 | 496 | 1.00 | 6.49 (5–11) | 6.0–14.0 |
| Baps Gibbs | 1-D, ν = 12, n/ν = 10 | 10.1 (9.0–11.5) | 7.0 (6.3–7.3) | 1.45 | 496 | 1.00 | 7.07 (5–10) | 9.1–13.4 |
| Baps Gibbs | 1-D, ν = 12, n/ν = 50 | 24.6 (21.8–28.3) | 20.8 (19.7–24.4) | 1.18 | 496 | 1.00 | 7.17 (6–9) | 8.1–11.3 |
| Baps Gibbs | 1-D, ν = 48, n/ν = 2 | 32.8 (30.1–36.0) | 23.6 (22.8–26.3) | 1.39 | 1360 | 1.00 | 5.10 (5–7) | 7.9–11.5 |
| Baps Gibbs | 1-D, ν = 48, n/ν = 10 | 90.7 (89.3–96.4) | 85.2 (75.1–88.1) | 1.06 | 1360 | 1.00 | 5.04 (5–7) | 9.7–14.0 |
| Baps Gibbs | 1-D, ν = 48, n/ν = 50 | 381.9 (345.3–451.7) | 363.2 (345.0–366.5) | 1.05 | 1360 | 1.00 | 5.04 (5–7) | 9.1–15.4 |
| Baps Gibbs | 1-D, ν = 192, n/ν = 2 | 309.9 (275.7–337.5) | 256.8 (209.2–315.1) | 1.21 | 4816 | 1.00 | 6.72 (5–9) | 7.7–17.0 |
| Baps Gibbs | 1-D, ν = 192, n/ν = 10 | 1209.8 (991.1–1367.1) | 1121.1 (1014.5–1288.8) | 1.08 | 4816 | 1.00 | 6.29 (5–7) | 7.1–16.1 |
| Baps Gibbs | 1-D, ν = 192, n/ν = 50 | 5281.1 (4879.9–5692.0) | 5132.3 (4802.9–5442.4) | 1.03 | 4816 | 1.00 | 5.08 (5–9) | 8.8–10.9 |
| Baps LocallyAdaptiveGibbs | 1-D, ν = 48, n/ν = 10 | 93.1 (89.8–98.3) | 84.2 (77.3–87.2) | 1.11 | 1360 | 1.00 | 5.07 (5–7) | 8.2–12.2 |
| Baps Collapsed | 1-D, ν = 48, n/ν = 10 | 17.5 (16.9–19.1) | 8.7 (7.9–9.6) | 2.00 | 1360 | 1.00 | 5.03 (5–7) | 7.3–12.9 |
| Hierarchical Centred | 1-D, ν = 24, R = 8, n/ν = 5 each | 110.2 (98.2–118.5) | 90.0 (82.8–109.1) | 1.22 | 1568 | 2.00 | 5.82 (5–9) | 6.2–9.8 |
| Hierarchical CollapsedPopulation | 1-D, ν = 24, R = 8, n/ν = 5 each | 165.1 (156.4–174.8) | 141.1 (132.9–160.5) | 1.17 | 1568 | 2.00 | 5.79 (5–9) | 8.6–11.5 |
| Hierarchical ParameterExpanded | 1-D, ν = 24, R = 8, n/ν = 5 each | 231.2 (205.0–256.4) | 215.5 (198.8–240.6) | 1.07 | 1568 | 2.00 | 6.19 (5–10) | 8.9–12.3 |
| Hierarchical Centred | 1-D, ν = 96, R = 8, n/ν = 5 each | 1404.7 (1258.7–1441.4) | 1312.8 (1187.0–1530.8) | 1.07 | 5024 | 2.00 | 6.02 (5–9) | 10.1–14.6 |
| Hierarchical CollapsedPopulation | 1-D, ν = 96, R = 8, n/ν = 5 each | 5433.6 (4945.4–5998.2) | 5146.6 (4955.7–5628.6) | 1.06 | 5024 | 2.00 | 6.42 (5–11) | 11.6–13.7 |
| Hierarchical ParameterExpanded | 1-D, ν = 96, R = 8, n/ν = 5 each | 2772.6 (2663.7–2999.8) | 2833.8 (2611.7–3314.6) | 0.98 | 5024 | 2.00 | 5.78 (5–11) | 9.1–12.0 |
| Baps AnisotropicGibbs | tensor 8x8 | 266.6 (251.0–302.4) | 209.6 (187.8–217.8) | 1.27 | 1744 | 1.00 | 7.03 (5–9) | 7.7–11.9 |
| Baps AnisotropicGibbs | tensor 16x16 | 3053.9 (2555.1–3109.8) | 2492.4 (2135.3–2567.3) | 1.23 | 6352 | 1.00 | 6.91 (5–9) | 7.4–10.4 |

**Readings.**
- *1-D Gibbs.* The ratio falls as the rest of the step grows: from 1.66 at ν = 12, n/ν = 2 to 1.03 at ν = 192, n/ν = 50.
  - Basis-dominated steps pay most: 1.21 to 1.66 at n/ν = 2.
  - Data-dominated steps pay 3% to 18% at n/ν = 50.
- *The collapsed scheme* doubles its step (2.00). Each step is one factored solve and the evidence, and the estimate costs as much as both.
- *Locally adaptive*, at the middle point: 1.11.
- *Hierarchical*, two estimates per step:
  - ν = 24: 1.07 to 1.22;
  - ν = 96: 0.98 to 1.07, inside the ranges.
- *Tensor*, anisotropic: 1.27 at 8×8 and 1.23 at 16×16.
- *Solves per estimate:* 5.0 to 7.2 on average. Between 5 and 11 were seen, never 4.
- *Allocation.* The estimate allocates per call, so the gate arm allocates about 0.5 KB more per step at ν = 12 and about 4.8 KB more at ν = 192. That churn is part of its time.

**Options.** Each is listed with what it costs and what it loses. No default is recommended from this grid.
1. *Keep the estimate on every factorization, as landed.*
   - Costs: the ratios above. That is up to 1.66 on small basis-dominated 1-D steps, 2.0 in the collapsed scheme, and 1.2 to 1.3 on tensor steps, but little on large data-dominated steps.
   - Loses: nothing.
2. *Estimate in the deterministic solves only (`PenalizedSpline`, `AnisotropicTensorSpline`), and admit sampler factorizations on positive definiteness.*
   - Costs: nothing per Gibbs step. The collapsed scheme's step is `PenalizedSpline.SolveFactored`, so it keeps its 2.0 unless that solve is counted as a sampler step.
   - Loses: refusal of a sampler system that stays positive definite with κ(E) above 1e14. That includes an unidentified design (strict alone admitted 88 of 234 in slice 5's probe) and a drawn λ far out (the 1e20 points). Such a draw carries no reliable digits.
3. *Estimate once per sampler, at the trace-balance point λ_ref, and check positive definiteness per step.*
   - The rule depends on the configuration alone, so a resumed run matches bit for bit.
   - Costs: one estimate per run.
   - Loses: the per-step check of a drawn λ far out. Slice 5's probe refused every unidentified design at every box point, its centre included, so the design check itself holds.
4. *Estimate every k-th step*, with k fixed by the configuration and the step index carried in the chain's state.
   - Costs: about 1/k of option 1's.
   - Loses: a refusal arrives up to k − 1 steps late, or never for an excursion that returns between checks.
5. *Alongside any of these:* an estimate that reuses workspace vectors cuts option 1's allocation and part of its time, and checks the same thing. Under options 2 to 4, `BandCholesky`'s ‖E‖₁ pass could also run only when an estimate is asked. Neither was measured.

**What the numbers cover.**
- The 1-D cubic Gibbs sampler on a 3 × 3 grid: ν = 12, 48 and 192, each at n/ν = 2, 10 and 50.
- The locally adaptive and collapsed schemes at one point: ν = 48, n/ν = 10.
- `HierarchicalPSpline`'s three schemes, which reach its five sites, at ν = 24 and 96, with R = 8 and n/ν = 5 each.
- Two tensor layouts.
- Not covered: the locally adaptive and collapsed schemes across sizes, other tensor layouts, the hierarchical sampler with anisotropic smoothing, and a quiet machine.

### 2026-11-14 — E8's readouts read a fit of several named outputs (outside the defects, at the user's word)

*The request.* A message the user wrote for the active-subspaces work reached this session. It asked for three things:
- route both forms of `fit.ActiveSubspace` for a fit of several named outputs;
- keep E34's same-fit pairing and E8's cross-fit pairing;
- refuse output patterns and declared locations.

It also asked for an acceptance test: W_y = aaᵀ against the co-active C_fg of (a, a). The user then let this session continue. The files belong to the E8 and E34 chips, and both heard which files and hunks before the commit. Files held:
- `BassActiveSubspace.cs`, `BassMultiOutputSubspaces.cs`, `ActiveSubspaces.cs`, `BassMultiOutput.cs` and `BassFit.cs`;
- `MultiOutputReadoutTests.cs` and `MultiOutputActiveSubspaceTests.cs`.

This builds what E34's brief records as unbuilt in its entry "E8's readouts on a fit of several outputs" (2026-11-13). The work fills these entries:
- both routes of `BassFit.ActiveSubspace`;
- the one-fit pair through `BassFit.CoActiveSubspace`;
- pairs across fits, with a named combination on either side.

**Commits.**
- `c543297` feat(bass): E8's active subspace and co-active pairs read a fit of several named outputs.

**What changed.**
- *The active subspace.*
  - *The draw.* Each draw is read as its terms and H_ka = B_{a+1,k}, taken from its output coefficients. Row zero is the intercept and has no gradient.
  - *The metric.* W_y is K × K and positive semidefinite, in `MultiOutputResponse.Names` order, with the identity as the default.
  - *The exact route* reads M = HᵀW_yH through the moment service, as a functional fit does.
  - *The sampled route* reads `BassMultiOutput.Gradient` in draw batches of about 2²² partials. Changing H alone would have left this route going through `fit.Gradient`, the scalar readout.
  - *K = 1.* Under the identity, the readout takes the scalar arithmetic (M = bbᵀ), so it is the scalar fit's on the same coefficients, bit for bit.
  - *The Poincaré check.* Its Jacobi ratio reads a scalar fit's coefficients, so it is not read here. The totals, the bounds and the matrix check are read.
- *Refused, each with its reason:*
  - `OutputMetric.Locations`, since named outputs have no location grid;
  - `OutputPatterns` > 0, since E[JᵀW_yJ] in input space stands alone and the output basis from E[J W_x Jᵀ] is outside this brief;
  - a metric that is not K × K. The message names the rows as outputs in Names order.
  - Basis adequacy now refuses a direct fit of several outputs as such; it used to call it scalar.
- *Pairs within one fit.* On such a fit, `BassFit.CoActiveSubspace(data, spec)` calls `BassMultiOutput.CoActive`, so there is one implementation. Draw d pairs with itself and carries w_d once.
  - Its sampled route now reads the gradient in the same bounded draw batches. A draw's pair has the same bits in any batch.
  - `BassMultiOutput.Functional` takes an optional parameter name, which the cross-fit side passes, so its refusal names `Output` or `OtherOutput`.
- *Pairs across fits.* Either side may take a named combination (`Output` or `OtherOutput`, one weight per output), read through H and E34's gradient. E8's positional selection and product weights are unchanged.
- *Removed.* `d9152e6` made these three entries refuse such a fit; that refusal and its test, `E8sReadouts_RefuseSeveralOutputs_NamingE34s`, are gone, because the entries now read such a fit. `BassActiveSubspace.RefuseSeveralOutputs` had no other caller.
- *Unchanged.* Scalar and functional fits compute the same bits. A scalar draw's coefficients now come from its surface's row, which copies the same values.

**Tests.** `MultiOutputActiveSubspaceTests` has 5 tests.
- *`UnderAaT_TheActiveMatrixIsTheCoActivePairOfAWithItself`.* With a = (0.6, −1.3) and W_y = aaᵀ, each draw's C is compared with C_fg of (a, a) on the same draws:
  - exact: 4.4e-16 relative;
  - sampled, against E34's sampled pair on the same 40,000 points: 4.3e-16;
  - sampled against exact: within five reported standard errors (largest |z| 1.68);
  - every draw, read across two draw batches, has the same bits as when read alone.
- *`TwoCombinationsOfOneFit_AreE34sPair_AndEachDrawCarriesItsWeightOnce`.* On both routes the pair is bit for bit, and the weights are w_d (the fixture's are unequal, 1:2:3).
- *`TwoFits_PairByPosition_WithProductWeights_ANamedCombinationOnEitherSide`.* The named fit is paired against a scalar fit that holds the projected surfaces.
  - Exact, on each side: the pairs are the one fit's, bit for bit, with weights w_d²/Σw².
  - Sampled: within 2e-16.
- *`AtOneOutput_UnderTheIdentity_TheReadoutIsTheScalarFits_BitForBit`.* With no metric and with [[1]]: each draw's C, E[C], the Poincaré totals and bounds, and the matrix check's ratios.
- *`NamedOutputs_RefuseLocationsPatternsAndAMetricNotKByK`.*
- `MultiOutputReadoutTests.Fixed` now takes `outputs`, which defaults to 2.

**Gates.** The snapshot is HEAD `8b68d63` plus the 7 files, with whole assemblies run.
- Regression: 3,408 of 3,408.
- Oracle: 212 of 214. The two failures are `DiagramMetricsTransportParityTests`, the transport oracle's R timeout recorded in slice 4. The active-subspace and co-active parity tests pass..

**Doc changes owed.** These supersede the ones in E34's entry of 2026-11-13.
- *`docs/Maths/Regression/ActiveSubspaces.md`*, where the readouts' responses are listed:
  - a fit of several named outputs reads C = E[JᵀW_yJ] with H_ka = B_{a+1,k} and W_y K × K in `MultiOutputResponse.Names` order;
  - it refuses locations and output patterns;
  - its pairs are read within one fit and across fits.
- *`docs/Maths/Regression/MultiOutput.md`'s readout paragraph:* E8's entry points read such a fit, and `BassFit.CoActiveSubspace(data, spec)` is `BassMultiOutput.CoActive`.

**Follow-up.** The output basis from E[J W_x Jᵀ] for named outputs would decide whether `OutputPatterns` stays refused for such a fit. Building it would give a K-dimensional output space's patterns, and the refusal would become a reading.
- *HEAD since the snapshot.* It moved to `94b11ae` with conformal changes (`CvPlus.cs`, `WeightedCrossFitTests.cs`). Those share no file and no code path with this change, so the gates were not rerun.

### 2026-11-14 — Slice 5 verified from a clean archive (coordinating session)

- **Run:** `scripts/verify-csharp.ps1` on one clean archive of HEAD `f1356cb`, the Commit selection over the eleven files of `f128e04` in each test project where they select a class.
  - `Maths.Regression.Tests`: 129 of 317 classes, 1,461 of 1,461.
  - `Maths.LinAlg.Tests`: 10 of 21 classes, 113 of 113.
  - `Maths.Samplers.Tests`: 1 of 31 classes, 8 of 8.
  - `Maths.Oracle.Tests`: 2 of 39 classes, 1 of 3.
- **The two Oracle failures** are `DiagramMetricsTransportParityTests` at p = 1 and p = 2: `transport_oracle.R` timed out after 120 s inside the R call. It is slice 4's finding for the TDA owner; the test reaches neither `BandCholesky` nor BAPS.
- **The source was read:** nothing found in the default's change, `BapsAdmission` or the eleven sites. `BapsAdmissionTests` ran in the regression selection.
- **What this does not cover:** the report's decision 1. The gate throws inside the samplers' steps, and neither its cost per step nor the condition numbers the samplers' drawn systems reach has been measured. The probe behind the 10¹⁴ bound covers the tensor layouts only.
- **Receipts:** `artifacts/verify-defects-s5-f1356cb/`, with the changed-file list, counters and snapshot line tracked in [receipts](receipts/verify-defects-s5-f1356cb/).

### 2026-11-13 — Slice 5: `BandCholesky` strict by default

The user chose on 2026-11-13: strict by default, with the rollout qualified in three directions. Those are a valid small-scale fit that BAPS refuses today, singular-design refusal, and the public λx = λy = 1e20 tensor points. `PivotedCholesky`'s solve is still undecided and untouched.

**The user's correction to slice 4.** Slice 4's decision 1 said that under a strict default the banded marginals would refuse a near-interpolating design the floor carries. That is backwards. BARS's banded g-prior already factors strictly. BAPS refused any factor that reported a floored pivot. So the floor carried no caller. The switch moves BAPS the other way: it can admit a valid fit in small units that the floor refused.

**Commits.**
- `f128e04` fix(linalg, baps): BandCholesky strict by default; BAPS admits a system by its equilibrated condition.

**What changed.**
- `BandCholesky`'s `strict` defaults to true, and `strict: false` keeps the 1e-12 floor. Every factor in `src` is now strict: BARS named it already, and BAPS takes the default.
- `BapsAdmission` (new, internal) admits a factored system when it is positive definite and its estimated κ(E) is at most 1e14. It replaces `Status is { PositiveDefinite: true, FlooredPivots: 0 }` at all eleven BAPS sites:
  - `PenalizedSpline`: 2;
  - `Baps`: 3;
  - `AnisotropicTensorSpline`: 1;
  - `HierarchicalPSpline`: 5. One of these is β₀'s dense collapsed precision, which never floored, so it gains the condition bound.
  - The messages say "not admissible".
- Five places in four test files exercise the floored mode itself, so they now ask for it by name (`strict: false`):
  - `BandCholeskyStrictTests`: 2;
  - `ConditionEstimateTests`: 1;
  - `FactorReferenceTests`: the band kinds of the status helper;
  - `CostProfileKernels`: the floored kernel.

**Why a condition gate: the probe.**
- *Setup.* The probe is in `artifacts/shared-code-defects/probe/` and is not committed. It covers the nine tensor layouts of `issues/linalg/spectral-tensor-factorization-report.md`.
  - *Points per layout.* The box has 26 points: the optimizer's 5×5 grid around each axis's trace balance, plus the nominal (0.7, 2.0). Then λ = 1e12 and λ = 1e20.
  - *Three families.* Identified. Small-scale: the same design × 2⁻²⁴ with λ × 2⁻⁴⁸. Unidentified: y fixed at 0.5, which leaves the penalty's y-linear null space unidentified.
  - That makes 756 systems.
- *Systems admitted by each gate.*

| family, points | today (floored, no floored pivot) | strict alone | strict, relative pivot ≥ 1e-12 | ≥ 1e-10 | BARS's `DesignAdmission` (1e-10, κ ≤ 1e10) | `BapsAdmission` |
|---|---|---|---|---|---|---|
| identified, box (234) | 234 | 234 | 234 | 234 | 234 | 234 |
| identified, λ = 1e12 (9) | 9 | 9 | 9 | 0 | 0 | 9 |
| identified, λ = 1e20 (9) | 1 | 1 | 0 | 0 | 0 | 0 |
| small-scale, box (234) | 0 | 234 | 234 | 234 | 234 | 234 |
| small-scale, λ = 1e12 (9) | 0 | 9 | 9 | 0 | 0 | 9 |
| small-scale, λ = 1e20 (9) | 1 | 1 | 0 | 0 | 0 | 0 |
| unidentified, box (234) | 55 | 88 | 5 | 2 | 0 | 0 |
| unidentified, λ = 1e12 (9) | 0 | 0 | 0 | 0 | 0 | 0 |
| unidentified, λ = 1e20 (9) | 1 | 1 | 0 | 0 | 0 | 0 |

- *κ(E) of the systems that factor strictly:*
  - identified box: 6.4 to 7.8·10⁹;
  - λ = 1e12: 1.3·10¹² to 2.5·10¹³;
  - unidentified box: 1.2·10¹⁶ to 3.1·10¹⁹;
  - the one λ = 1e20 survivor (6×16): 2.7·10¹⁹.
  - The small-scale family reads the same κ as the identified family; the estimate is scale-free.
- *Strict alone loosens singular refusal* (55 → 88 admitted). What refused those systems was BAPS's floored-pivot check, and it is vacuous once nothing floors.
- *Relative-pivot floors leak* (5 and 2 admitted). The 1e-10 floor also refuses the public λ = 1e12 points, and so does BARS's gate. Its 1e10 sits only 1.3× above the box's largest κ.
- *Only the condition separates the families.* The bound of 1e14 is four decades above the box, about 4× above the λ = 1e12 maximum, and two decades below the smallest singular system. At κ·ε ≈ 0.02 a backward-stable solve still keeps about two digits.

**Bits.**
- On a system the old gate admitted, nothing floored. The two modes differ only where the floor fires, so the strict factor is the same bit for bit. The gate adds a read-only estimate that draws no random numbers. Every fit and draw the old gate admitted is therefore unchanged, and the existing pinned tests pass as they were.
- A newly admitted small-scale fit is the same fit in other units: β' = β/c bit for bit.
- *Newly refused systems:* the 55 unidentified box systems, the 6×16 point at 1e20, and any system with κ(E) above 1e14 that never floored. No existing test reaches one.

**Tests.** `tests/regression/BapsAdmissionTests.cs` has 19 cases.
- `AValidFitInSmallUnits_IsAdmitted_AndIsTheSameFit`. The design is scaled by 2⁻²⁴, on `PenalizedSpline` at λ ∈ {1e-3, 1, 1e3} and on the 8×8 tensor across its box. β' = β/c bit for bit, RSS_pen is equal, and log|A| shifts by ν·log c².
- `AnUnidentifiedDesign_IsRefusedAcrossTheBox`, on the nine layouts.
- `AnIdentifiedDesign_IsAdmittedAcrossTheBox_AndAtTheTwelfthPower_ButNotTheTwentieth`, on the nine layouts.
- *On HEAD before the change, 10 of 19 fail:* the small-scale fit, 8 of the 9 unidentified layouts, and 6×16 at 1e20. The receipt is `artifacts/shared-code-defects/receipts/slice5-reproduction-at-HEAD.trx` (local).

**Gates.** The snapshot is HEAD `23ce73d` plus the slice's 11 files, built in Release, with whole assemblies run.
- LinAlg: 190 of 190.
- Samplers: 296 of 296.
- Regression: 3,403 of 3,403.
- While the gates ran, HEAD moved to `01eb63d` (`d9152e6`, E8's readouts). That change touches none of this slice's files and uses neither `BandCholesky` nor BAPS, so the gates were not rerun.

**Fresh-eyes review.**
- Every `BandCholesky` in `src` is strict: BARS names it, and BAPS takes the default behind `BapsAdmission`. No caller lost a floor it relied on.
- `Admits` short-circuits on a failed factor, so the estimate never runs on one. The estimate also refuses a floored factor, and BAPS never builds one.
- Every test that constructs a band factor lies in the three gated assemblies.

**Decisions for the user.**
1. *The estimate's cost in BAPS's Gibbs steps* (not measured). `Admits` estimates κ on every factorization, including every sampler step.
   - *The arithmetic.* Higham's refinement of Hager's method takes 4 to 11 solves with the factor, usually 4 or 5. Each solve costs about 2νp, against the LDLᵀ factor's νp²/2. That makes the estimate about 16/p to 44/p factorizations' worth.
     - For the 1-D cubic P-spline (p = 3), that is 5 to 15 factorizations per step.
     - For the tensor layouts (p in the tens), it is under one.
   - *No cheap way to admit.* A small relative pivot proves a large κ, but a large one proves nothing; the probe's leaks show it. So the estimate is what admits.
   - *The options.*
     - Keep it on every step.
     - Estimate in the deterministic solves (penalized, REML) only, and let the samplers check positive definiteness alone. The cost is that a drawn λ far out would go unchecked.
   - A step-time measurement on a 1-D BAPS fixture is worth taking only if the answer to this question turns on it.
2. *`PivotedCholesky`'s solve* is still open, as in slice 4.

**Doc changes owed.**
- *`docs/Maths/LinearAlgebra.md`:*
  - `BandCholesky` is strict by default, and the floor is opt-in by `strict: false` (around lines 269 and 277, and the floor's table row at line 447).
  - Callers judge a band factor by `EquilibratedCondition`.
  - Slice 4's "`BandCholesky`'s floor no longer mirrors anything" is superseded: the floor is opt-in.
- *`docs/Maths/Regression/BAPS.md`, "The rule" (line 171):* replace it with `BapsAdmission`, which admits a system that is positive definite with κ(E) ≤ 1e14. Add the probe's table and the alternatives it rules out: the floored-pivot rule, the relative-pivot floors, and BARS's 1e10.

### 2026-11-13 — Slice 4 verified from clean archives (coordinating session)

Over the nineteen files of `f448931`, `43bd61f`, `ff5d635`, `e6dc220` and the review fix `78603a2`.

- **The regression project:** `scripts/verify-csharp.ps1` on a clean archive of HEAD `a3ed9c0`, the Commit selection: 308 of 314 classes, 3,049 executed, 3,049 passed, 0 not executed.
- **The other projects,** from a clean archive of HEAD `8db288f`, after a sweep of every test project for the classes the files select (`closure.txt`):
  - `Maths.LinAlg.Tests`, whole: 190 of 190.
  - `TDA.Mapper.Tests`, whole: 28 of 28.
  - `VizCore.Tests`, whole: 703 passed and 2 skipped of 705. The skips are the two Möbius golden drafts.
  - `Maths.Samplers.Tests`, by selection: 8 of 8.
  - `Maths.Oracle.Tests`, by selection with the pinned tier: 146 of 147.
- **The one failure** is `DiagramMetricsTransportParityTests.DiagramDistances_MatchT4TransportAndLpAssign(p: 1)`: `transport_oracle.R` timed out after 120 s, inside the R call. The test names none of the changed types, and its p = 2 case passed. It is the report's finding for the TDA owner, and no part of this slice.
- **The other nine projects** that compile a changed file have no selected class and did not run.
- **The source was read** at `e6dc220`, and the two notes it drew landed in `78603a2`: the scene export's ε relative at every scale, and the remainder pass's cost in `PivotedCholesky`'s doc. Nothing else was found.
- **Receipts:** `artifacts/verify-defects-s4-a3ed9c0/` and `artifacts/verify-defects-s4-projects-8db288f/`, tracked in [receipts](receipts/) with the project script as it ran.

### 2026-11-13 — Slice 4: slice 2's three decisions, built

The user chose on 2026-11-13: the review's recommendations (`discussion/astra-shared-code-cleanup.md`), with this chip's two adjustments. The scene export draws a singular block instead of refusing it, and the pivoted factor's remainder check reads the off-diagonals too. The review's order: the floor first, the pivoted factor apart, the callers migrated and qualified, the legacy class deleted in its own commit.

**Commits.**
- `f448931` fix(linalg): CholeskyDecomposition drops its pivot floor.
- `43bd61f` feat(linalg): PivotedCholesky, the semidefinite factor BGP's drawn function reads.
- `ff5d635` refactor(gmm, viz): a Gaussian's covariance factor is DenseCholesky.
- `e6dc220` chore(linalg): retire CholeskyDecomposition and BGP's own pivoted routine.
- `78603a2` fix(viz): a singular block's ε is relative below unit scale too (after the coordinating session's review).

**1. The floor dropped (`f448931`).**
- *The class is strict.* The first pivot that is not finite and positive stops it, at the row `DenseCholesky` stops at, and an L⁻¹ past `double.MaxValue` fails it. On every input of a ladder the two agree on the failure row, `Status`, factor and `LogDet`. The ladder: 60 SPD matrices at scales 2⁻⁸⁰ to 2³⁹, 16 squared-exponential kernels, 10⁻¹³·I, and indefinite and non-finite inputs. 10⁻¹³·I now factors exactly, where the floor read every variance ten times too large.
- *The earlier class, kept verbatim and pinned*, reads bit for bit on 60 inputs the floor never reached.
- *The callers.*
  - GMM's `UpdateCache` and the Bhattacharyya distance refuse a covariance that is not positive definite. A warm start with a singular covariance is refused, an acceptance case.
  - The scene export draws a singular block as a thin ellipsoid. It adds ε = 10⁻⁸ of its largest variance, and refuses only a block that is still indefinite. The viewer loads all nine entries as a general 3×3 map (`viewer.html`, `buildGaussianGroup`), so a rank-deficient factor would make its normal matrix singular, and the ε keeps it invertible.
  - Rank one at 10²⁰⁰, which overflowed the floored factor, now exports a finite needle. `[[1, 2, 0], [2, 1, 0], [0, 0, 1]]` is refused, naming its layer and component.

**2. The pivoted factor promoted (`43bd61f`).**
- *The contract is narrow and sampling-first:* the rank, the loadings in the matrix's own row order, and `Map` from Rank standard normals to a draw of N(0, LLᵀ).
- *BGP's loop, verbatim.* The loadings match the earlier routine bit for bit (pinned, 40 kernels and low-rank Gram matrices). BGP's drawn function keeps its seeded draws: its tests, including the bit-for-bit candidate, and the conformal tests that read it pass.
- *The remainder is checked.* The earlier routine read `[[1, 2], [2, 1]]` as rank 1, truncating −3, and `[[0, 1], [1, 0]]` as rank 0. A semidefinite remainder whose diagonals are at most the stop has every entry within it, so the check reads the off-diagonals as well as the diagonal.
- *Measured: two tolerances.* BGP's joint posterior covariances leave remainders up to 1.2·10⁻¹¹ of their largest variance. That is rounding of the computed covariance, past the 10⁻¹² stop: the first form, which judged the remainder at the stop, refused two of BGP's tests.
  - The rank still stops at 10⁻¹².
  - Indefiniteness is judged at its own tolerance, √ε (1.5·10⁻⁸) by default, far below the refused cases.
  - A remainder `[[0, δ], [δ, 0]]` at δ = 5·10⁻¹¹ is kept at the default and refused at 10⁻¹².
- *The cost.* The factor costs O(p·k²) at rank k, and the remainder's check O((p − k)²·k), which is above the factor when p is far above k. For BGP both cost less than forming the covariance. The class's doc says so (`78603a2`).
- *The solve waits*, as the review asked. On a singular matrix it must refuse or project a right-hand side outside L's range, and BGP's sampling needs neither.

**3. The callers moved to `DenseCholesky` (`ff5d635`).**
- *`TryWriteInverseTo`* forms A⁻¹ on demand. It ports the retired L⁻¹ arithmetic, so GMM's Σ⁻¹ keeps its bits (pinned against the arithmetic kept verbatim). It refuses an inverse that overflows: `[1e-320]` factors, and its reciprocal does not fit, and a power-of-two bidiagonal L⁻¹ overflows at row 51. GMM's component is then left unusable (NaN caches, an exception) rather than holding a partial precision matrix.
- *`AddSample`* accumulates L·z into the mean in the retired order. GMM's seeded Box–Muller draws match bit for bit (pinned).
- *The exporter* reads `Factor(i, j)`: the same factor, the same bits.
- *Bhattacharyya* reads its mean term as ‖L⁻¹δ‖² through `InverseQuadraticForm`, so its distances move in their last bits. They agree with the explicit-inverse formula within 10⁻¹³ relative.

**4. The legacy class retired (`e6dc220`).**
- *Deleted:* `CholeskyDecomposition` with its own tests, and BGP's own pivoted routine.
- *Moved to `DenseCholesky`:* `FactorReferenceTests`, `DenseCholeskyTests` and `BandCholeskyTests`.
- *The cost profiles* time `DenseCholesky` as production. Their two rows that timed the retired path's `Decompose` are gone; the replica kernels of its two passes stay.

**The gates.**
- *Each commit*, from a snapshot of HEAD plus its files, by selection with the pinned tier:
  - A: linear algebra 90 of 90, VizCore 9 of 9, regression 80 of 80;
  - B: linear algebra 14 of 14, regression 145 of 145;
  - C: linear algebra 93 of 93, VizCore 12 of 12.
- *The end*, from a snapshot with D:
  - `Maths.LinAlg.Tests` 190 of 190;
  - `TDA.Mapper.Tests` 28 of 28;
  - `VizCore.Tests` 703 passed, with the 2 Möbius skips that stood before;
  - `Maths.Regression.Tests` 3,370 of 3,370 (18 min);
  - `Maths.Oracle.Tests` by selection, 145 of 147. The 2 failures are `transport_oracle.R` timing out at 120 s, alone as well as under load. Their test reads `TDA.Ph` and R only, and no TDA code reaches a changed factor, so the timeout stands apart from this slice (below).
- Receipts are in `artifacts/shared-code-defects/receipts/` (`slice4-*`).

**Defects found.**
- *One in the slice's code, found in the coordinating session's review and fixed in `78603a2`.* The export's ε started its largest variance at 1, so below unit scale it was absolute: a singular block with variances near 10⁻⁸ drew as thick as it was long.
  - It is now 10⁻⁸ of the block's largest variance at every scale, which is what the display needs: a flat ellipsoid stays flat.
  - An all-zero block keeps the missing dimensions' absolute 10⁻⁸, a sphere of radius 10⁻⁴.
  - Gates: VizCore 7 of 7 and linear algebra 14 of 14 by selection with the pinned tier.
- My own first remainder bound refused BGP's covariances, as above. One test used an assertion precision the analyzer refuses, and another rested on noise that could move the rank; both were rewritten before any run counted.
- *Finding: the transport oracle no longer fits its timeout.* `DiagramMetricsTransportParityTests` (both p) times out in `transport_oracle.R` at 120 s, even run alone. It has nothing to do with this slice, and is for the TDA owner: either the timeout or the fixture's size (s = 100 slices).
- *Finding: `BandCholesky`'s default floor is the same rule.* It "mirrored" the retired one: an absolute 10⁻¹², on by default, with strict mode a choice. By the same arithmetic it can introduce order-one error at small scale; this slice did not measure it.

**Decisions for the user.**
1. *`BandCholesky`'s default.* If it becomes strict, the banded marginals refuse a near-interpolating design that the floor now carries, and their numbers move only where the floor fires. If it stays floored, it is the one floored factor left, and its docs should say so. Recommendation: strict by default, as here, with the over-knotted design refused by its caller's prior. *Corrected in slice 5:* the floor carried no caller, since BARS factors strictly and BAPS refused a floored pivot; the switch admits a valid small-scale fit instead.
2. *`PivotedCholesky`'s solve*, when a caller needs one: refuse a right-hand side outside the retained range, or project it to the minimum-norm least-squares answer.

**Doc changes owed.**
- *`docs/Maths/LinearAlgebra.md`:*
  - "Dense Cholesky (`CholeskyDecomposition`)" (line 193) and its rows in the API list (line 95), the module table (line 12), the diagram (line 33), the allocation note (line 73) and the storage table (line 81): retired.
  - `DenseCholesky` gains `TryWriteInverseTo` and `AddSample`, and serves the Gaussians.
  - A section for `PivotedCholesky`: its contract, the two tolerances, the remainder rule, its cost (the check's O((p − k)²·k) above the factor's O(p·k²) when p is far above k), and the measurements above.
  - `BandCholesky`'s floor no longer mirrors anything.
- *`docs/Clustering/GaussianMixtures.md`:*
  - line 11: `DenseCholesky`;
  - line 45: `UpdateCache` refuses a covariance not positive definite, or one whose inverse overflows;
  - the warm start (line 89): a singular covariance is refused;
  - "Bhattacharyya distance" (line 162): refusal, and the mean term through the factor.
- *`docs/Maths/Regression/BGP.md`*, the drawn function: it reads `PivotedCholesky`, whose remainder is checked.
- *Slice 2's "Doc changes owed"* for the floor are superseded: there is no floor left to describe.

### 2026-11-13 — Slice 3 verified from a clean archive (coordinating session)

- **Run:** `scripts/verify-csharp.ps1` on a clean archive of HEAD `badc0d9`, the Commit selection over slice 3's fourteen files (`c9f9a84`). All fourteen lie in the regression project's compile closure.
- **Selection:** 237 of 313 classes, with the pinned tier and without the slow tier's 71 methods.
- **Result:** 2,659 executed, 2,659 passed, 0 failed, 0 not executed (29 min). `DrawWeightRangeTests` is among them.
- **The source change was read:** nothing found. Each reader forms a ratio, so the power-of-two scaling cancels and ordinary weights keep their bits. The positivity checks still read the raw weights.
- **The two shared files:** the commit's diff of `SubspacePosterior.cs` and `CoActivePosterior.cs` holds this slice's lines only, and E8 slice 7e's orientation stands in the archive.
- **Not rerun here:** the `Maths.Oracle.Tests` selection, whose R oracles time out under this machine's load. It stands on the chip's run (112 of 116, the four timeouts passing when rerun singly).
- **Receipts:** `artifacts/verify-defects-s3-f1d70d5/` (named for the HEAD when the run was queued; its snapshot line records `badc0d9`), tracked in [receipts](receipts/verify-defects-s3-f1d70d5/).

### 2026-11-13 — Slice 2 verified from clean archives (coordinating session)

- **The standard run:** `scripts/verify-csharp.ps1` on a clean archive of HEAD `78d5077`, the Commit selection over slice 2's six files (`37eb63c`). It selects within the regression test project only: 2 of 308 classes, 79 executed, 79 passed, 0 not executed.
- **The other projects.** Five of the six files lie outside that project, so each test project that reaches them ran whole, from a clean archive of HEAD `ed7df5d`:
  - `Maths.LinAlg.Tests`: 179 of 179, `CholeskyDivergenceTests` among them.
  - `TDA.Mapper.Tests`: 28 of 28.
  - `VizCore.Tests`: 699 passed and 2 skipped of 701, `GaussianFactorRefusalTests` among them. The skips are the two Möbius golden drafts, skipped before this slice.
- **The sweep behind that choice** (`closure-78d5077.txt`): `Cholesky.cs` compiles into 14 test projects.
  - The name graph selects classes in three: linear algebra (7 of 19), VizCore (2 of 85) and regression (2 of 308).
  - Mapper holds three of the files and reaches the component without naming it.
  - `Clustering.Graphical.CumulativeField.Tests` also ran whole in the sweep's first pass: 15 of 15.
  - The other nine have no selected class and did not run.
- **The whole regression assembly** less its slow tier passed with the refusal in place, under slice 1's run (2,922 of 2,922 at `e39924c`).
- **The source change was read:** nothing found. A non-finite entry of L shows in its own row's pivot, and every read is gated on one flag.
- **Receipts:** `artifacts/verify-defects-s2-e429b29/` and `artifacts/verify-defects-s2-projects-ed7df5d/`. The lists, the counters and the one-off project script are tracked in [receipts](receipts/).

### 2026-11-13 — Slice 1 verified from a clean archive (coordinating session)

- **Run:** `scripts/verify-csharp.ps1` on a clean archive of HEAD `e39924c`, the Commit selection over slice 1's three files (`f50753f`).
- **Selection:** 288 of 308 classes of the regression assembly, with the pinned tier and without the slow tier's 71 methods.
- **Result:** 2,922 executed, 2,922 passed, 0 failed, 0 not executed (26 min). `WeightReaderTests` is among them.
- **The source change was read:** ordinary weights take the earlier expression, so their bits hold; nothing was found.
- **Slice 2** was verified later the same day (the entry above). This archive holds its source (`37eb63c`), so the regression assembly passed with the refusal in place.
- **Receipts:** `artifacts/verify-defects-s1-e39924c/`, with the changed-file list, counters and snapshot line tracked in [receipts](receipts/verify-defects-s1-e39924c/).

### 2026-11-13 — Slice 3: BASS's draw weights over their whole range

Added by the user to this chip on 2026-11-13, after the review in `discussion/astra-shared-code-cleanup.md`. It is slice 1's decision for BASS's own readers.

**Commit.** `c9f9a84` fix(regression): BASS's draw-weight readers keep every finite weight's share.

**The readers.** Eleven sites normalized linear weights through a plain sum or summed raw weights times values. All of them read `BassFit.DrawWeights`, an emulator's weights, or a subspace posterior's weights:
- *the selection:* `Draws3.Select`, which prediction, gradients, sensitivity, the active subspaces and the multi-output readouts read;
- *the summaries:* `PosteriorSummary.Over` (the fit's summary, the Sobol summaries) and `SobolDecomposition`'s means;
- *the resamplers:* `NestedCalibrationSampler.OuterDraws` and `BassJointDraws`' resampling;
- *the residual check:* its posterior mean, and `BassDrawFit.MeanOffset`;
- *the projection report:* its mean term count;
- *the subspaces:* `Whitening.Normalized`, `SubspacePosterior.Of` and `CoActivePosterior.Of`.

`DrawSample.Quantile` already reads weights relative to the largest, and `PredictiveChecks.Rows` was repaired in slice 1. The pointwise and candidate densities pass the weights to `Rows`.

**The repair.**
- *One routine, `WeightScale`* (`src/maths/regression/WeightScale.cs`). Each reader takes its weights scaled by 2^−e, the power of two that brings the largest into [1, 2), then does its own arithmetic unchanged.
- *Why ordinary weights keep their bits.* Each reader forms a ratio: a normalized weight, a weighted mean, a cumulative pick against u·Σ. The scaling is exact wherever the scaled weight is normal, and it cancels in the ratio, so ordinary weights give the bits they gave.
- *Why the extremes now read.* A sum of n scaled weights is below 2n, and a weighted sum of values stays below twice the values' own absolute sum.
- *Where the exponent comes from.* The selection takes it over the requested draws, repeats included; the Sobol means over the defined draws; every other reader over all its weights.
- *What stays as it was.* Zero-weight behavior and every refusal message at each entry point. The summary still drops a zero-weight draw from its count and range, testing the raw weight.
- *`PredictiveChecks.Rows`' fallback* now calls the same routine, with the same arithmetic.

**The gates.** Acceptance from the review: both extreme rows reproduce the ordinary row; equal weights resample the same indices at a fixed seed; the selection normalizes over the requested draws, repeats included; zero-weight behavior and the ordinary-range pins hold.
- *The reproduction*, `tests/regression/DrawWeightRangeTests.cs`, with the rows `[1, 1]`, `[1e308, 1e308]` (the sum overflows) and `[5e307, 5e307]` (the sum is finite but the weighted sums overflow). It also uses exact power-of-two multiples (2^±1000) and subnormal weights where their ratios are exact.
  - On HEAD's code, 9 of the 14 tests fail:
    - the selection read `[0, 0]`;
    - the summary's mean read NaN at `[1e308, 1e308]` and ∞ at `[5e307, 5e307]`;
    - every outer calibration draw fell to the last draw;
    - the joint draws resampled other ancestors;
    - the subspace weights normalized to 0;
    - the projection's mean term count read NaN;
    - the residual check threw on a non-finite residual.
  - The 5 that pass on HEAD are the cases that never overflowed, and the pinned comparison.
  - On the build, all 14 pass.
- *What the tests check.*
  - The two summary values read mean 3, SD 1 and median 3 at every row.
  - The selection, the subspace posteriors, the outer draws, the joint draws, and the fit's summary, prediction and sensitivity read the ordinary row's values. So do the residual check and the projection report. Power-of-two multiples match bit for bit; the extreme rows match within 10⁻¹¹ relative.
  - The pinned test compares 2,000 SMC-like weight vectors with the earlier selection and summary, kept verbatim, bit for bit.
- *The whole regression assembly*, from a snapshot of `77f80ba` plus the fourteen files: 3,341 of 3,341 (50 min).
- *`Maths.Oracle.Tests`*, which reaches the readers through 22 classes, by selection with its pinned tier: 112 of 116.
  - The 4 failures are R oracles timing out at 120 s on the loaded machine: `bass_enumeration_oracle.R` three times and `gbass_readouts_oracle.R` once. Each failed inside the R call, before any comparison.
  - Rerun one at a time, all 21 cases of those 4 tests pass (17 min).
- *The merged state.* E8 slice 7e (`199abd7`, the eigenvectors' signs) and conformal slice 15 (`55f4b93`) landed during the run, slice 7e in two of this slice's files. So the commit gate ran again from `048bbc4` plus the fourteen files, by selection (236 of 312 classes) with the pinned tier: 2,646 of 2,646 (12 min).
- Receipts are in `artifacts/shared-code-defects/receipts/` (`slice3-*`).

**Shared files.** E8 slice 7e committed its eigenvector orientation to `SubspacePosterior.cs` and `CoActivePosterior.cs` while this slice's gates ran.
- This slice's three lines in each sit beside that work, and its commit is the working tree, which is that HEAD plus them.
- A first staging, from the snapshot taken before slice 7e landed, would have reverted slice 7e's lines. Its staged diff showed it before any commit, and it was unstaged.

**Defects found.** None in the slice's code. My own test first built `2·1e308/3`, an infinite weight; it was rewritten to divide first before any run counted.

**Decisions for the user.** None for this slice. Slice 2's three decisions stand.

**Doc changes owed.**
- *`docs/Maths/Regression/BASS.md`*, `DrawWeights` (line 699): every view reads any finite weights, over their whole range, by `WeightScale`'s rule. Ordinary weights read as before.
- *`docs/Maths/Regression/Regression.md`*, the weighted draws (line 214): one rule for linear and log weights, `WeightScale` and `Rows`.
- *`docs/Maths/Regression/ActiveSubspaces.md`*, the summaries (line 65) and the posterior (line 105): `SubspacePosterior` and `CoActivePosterior` read their weights over the same range.

### 2026-11-13 — Slice 2: `CholeskyDecomposition` on a near-singular covariance

**Commit.** `37eb63c` fix(linalg): CholeskyDecomposition refuses a factor that is not finite, and its callers read it.

**The reading.** The probe is `artifacts/shared-code-defects/probe/cholesky-floor.cs`, with its tables in `results-before.md` and `results-edges.md` beside it.
- *The reproduction: squared-exponential covariances*, variance 1, on 50 to 401 points of [0, 1], lengths 0.05 to 1, no jitter.
  - 19 of the 20 configurations floor pivots, and 15 of those overflow: 351 to 73,536 entries of the factor are not finite, the first at rows 17 to 57.
  - `LogDet` still reads a finite −1,156 to −11,019 on every one.
  - `DenseCholesky` stops on all 19, at rows 5 to 16.
  - The other 4 floor 40 to 398 pivots and stay finite.
- *Where the floored factor goes wrong while finite.* A = s·Q·diag(κ^(−i/39))·Qᵀ with n = 40, s from 10⁻¹⁶ to 10⁶ and κ from 1 to 10²⁰. The reference is a double-double factor of the stored A.
  - *The floor fires on scale alone.* The pivots lie in [s/κ, s], so the floor fires when s/κ reaches 10⁻¹². At κ = 1, 10⁻¹³·I floors all 40 pivots.
  - *Once a pivot is floored, the outputs leave any tolerance while every entry stays finite.*
    - With 4 to 8 pivots floored, `LogDet` is off by 3.6 to 13 nats and the solves and inverse by 73 to 93 % (s = 1 at κ = 10¹⁴, 10⁻⁶ at 10⁸, 10⁻⁸ at 10⁶).
    - Across every cell where flooring changed a pivot, `LogDet` is off by 3.6 to 1,200 nats and the solves by 46 % up to 2·10¹¹ relative. At s = 10⁻¹², κ = 1, the floor replaces 10⁻¹² by itself, and the error there is 0.
    - The law `Sample` draws from, LLᵀ, is up to 5.9·10⁴ from A in relative terms, among the cells that have a reference.
  - *`DenseCholesky` on the same matrices is accurate to its κ·ε*: 2.4·10⁻¹⁰ nats and 1.4·10⁻¹⁰ at s = 10⁻⁶, κ = 10⁸.
    - It refuses only where the stored A has no factor in double: κ = 10²⁰ at every scale, and κ = 10¹⁸ at most.
    - There the floored factor returns finite numbers for a matrix that is not positive definite. At s = 10⁶, κ = 10²⁰, its LLᵀ is 5.9·10⁴⁴ from A.
  - *Where no pivot is floored, the two factors agree bit for bit.*
- *Two edges.*
  - A finite matrix can overflow a pivot while L stays finite: rank one at 10²⁰⁰ with one entry 2ε off reads `LogDet` 405.3, and its Σ⁻¹ has four non-finite entries.
  - L⁻¹ can overflow while L and every pivot are finite. The tests construct one.
  - `Status.Finite` catches the first edge and not the second.
- *How `DenseCholesky` treats the same inputs.* It stops at the first pivot that is not finite and positive and names it in `FailedAt`. `LogDet` is NaN, and every read, `RelativePivot` included, throws.
- *What the callers did with `Status`: none read it.*
  - `GaussianComponent` wrote its inverse and log-normalizer from any factor.
  - `BhattacharyyaCoefficient` read `LogDet` and the inverse.
  - The scene export wrote `WriteLTo`'s factor as it came: the edge above exported a finite factor whose LLᵀ is not the covariance.

**The build.**
- *The refusal.* `Decompose` fails when a pivot before flooring, an entry of L or an entry of L⁻¹ is not finite.
  - `FailedAt` (new) names the first such row, and `LogDet` is NaN.
  - `Solve`, `WriteInverseTo`, `WriteLTo` and `Sample` throw, and they also throw before any `Decompose`.
  - A floored factor that stays finite is read as before.
- *The callers.*
  - `GaussianComponent.UpdateCache` throws and leaves `CovarianceInverse` and `LogNormalizationFactor` NaN, so its density reads NaN and `Sample` is refused.
  - `BhattacharyyaCoefficient` refuses by name (`cov1`, `cov2`, `a`, `b`, or the averaged covariance).
  - The scene export refuses the component, naming its layer and index.
- *The two factors now fail one way on a factor that is not finite.* They still differ on a floored factor that stays finite (decision 1).

**The gates.**
- *The reproduction.*
  - On HEAD's code, 8 of the 11 new linear-algebra tests fail, every one that asserts a refusal. All 3 caller tests fail too: the component's cache, the distance (which read a finite −1,156 for `LogDet`), and the export of a meaningless factor.
  - The 3 that pass on HEAD hold what the build leaves alone: floored factors that stay finite are read, and the pinned comparison.
  - On the build, all of them pass.
- *Bit for bit.* The class's existing tests hold, `DenseCholeskyTests`' exact agreement included. The new pinned test compares 63 matrices with the earlier class kept verbatim: the factor, `LogDet`, the inverse, a solve and a sample, bit for bit.
  - 60 are random SPD matrices at scales 2⁻⁸⁰ to 2³⁹; the test requires at least 10 of them to floor.
  - The other 3 are two floored kernels and 10⁻¹³·I.
- *Whole projects*, from a snapshot of HEAD plus the six files:
  - `Maths.LinAlg.Tests`: 179 of 179.
  - `TDA.Mapper.Tests`: 28 of 28 (`MapperGMM` reaches the component).
  - `VizCore.Tests`: 699 passed, 2 skipped. The skips are the Möbius golden drafts, skipped before this slice.
  - `Maths.Regression.Tests`: 3,313 of 3,313 (27 min).
  - The clustering projects (SPC, cumulative field) do not reach the class: their selection is empty and they have no pinned tier. Both build.
- Receipts are in `artifacts/shared-code-defects/receipts/` (`slice2-*`).

**Defects found.**
- None in the slice's code.
- *The defect left by design* (decision 1): a floored factor that stays finite is read, and its `LogDet` and solves can carry order-one error. `Status.FlooredPivots` is its only signal, and no caller reads it.
- *Where GMM leans on the floor.*
  - `CovarianceRegularization` may be 0.
  - `InitializeWithParameters` takes covariances as given.
  - By the arithmetic, not measured: the 10⁻⁶ ridge falls below the pivots' rounding at large covariance scales.
  - In each, a singular covariance is floored, stays finite, and is wrong.

**Decisions for the user.** All three are unbuilt.
1. *The floor's rule.*
   - *Measured.* The absolute 10⁻¹² fires on scale alone, and flooring can introduce order-one error. A pivot at the floor is unchanged. Wherever flooring raised a pivot by orders of magnitude, the outputs carry the errors above: at κ = 1, 10⁻¹³·I has every variance multiplied by ten.
     - A floor relative to the largest diagonal would make the class scale-invariant: each κ column of the ladder would read the same at every s. Once it fires, it is just as wrong.
     - Pivots negative beyond rounding (the stored matrices at κ ≥ 10¹⁸) are floored silently today.
   - *Recommendation: no floor.* Report a pivot that is not positive, as `DenseCholesky` does, and let callers regularize. If a floor stays, make it relative to the largest diagonal, and refuse pivots below −n·ε·max A_jj.
     - Corrected 2026-11-13 (review in `discussion/astra-shared-code-cleanup.md`): not every caller regularizes.
     - GMM's EM adds its ridge, but its warm start and a regularization of 0 take covariances as given.
     - The export's ε fills only the display dimensions a layer lacks; it copies a supplied block as given.
     - Under no floor, a singular covariance in any of these is refused. Those refusals are acceptance cases.
   - *What changes.* GMM's numbers move only where the floor fires today, in the cases above. There `UpdateCache` would throw instead of returning a floored fit.
2. *A pivoted factor in the linear algebra.*
   - *Measured.* Diagonal pivoting, stopped at 10⁻¹² of the largest diagonal, factors every kernel of the reproduction at rank 8 to 57. Its ‖A − LLᵀ‖_F/‖A‖_F is at most 1.3·10⁻¹³, on the same kernels where the floored factor overflowed or floored up to 398 pivots.
   - *Recommendation.* Promote `BgpDrawnFunction`'s factor (Harbrecht, Peters and Schneider 2012) into `Maths.LinAlg` as the semidefinite factor, with a low-rank solve and sample. BGP's drawn function would be its first reader. It is the tool for smooth kernels' covariances; GMM does not need it.
3. *Retiring `CholeskyDecomposition` into `DenseCholesky`.* This is a removal, so it is its own question.
   - *What is lost:*
     - the soft floor (decision 1);
     - `WriteInverseTo`'s Σ⁻¹ through L⁻¹, which GMM's `CovarianceInverse`, its Mahalanobis distances and its mode ascent read;
     - `Sample`'s μ + Lz (`DenseCholesky` offers L⁻ᵀz);
     - `WriteLTo` (`DenseCholesky` offers `Factor(i, j)`).
   - *What moves:*
     - the callers: `GaussianComponent`, `BhattacharyyaCoefficient` (both overloads) and the scene export;
     - the tests: `DenseCholeskyTests` (2 comparisons), `FactorReferenceTests` (3 sites), `BandCholeskyTests` (3) and `CholeskyDivergenceTests`;
     - the cost profiles `CostProfileKernels`, `CostProfileLevers` and `CostProfileEngines`, which time the class as BGP's earlier path.
   - *Bits.* `DenseCholesky`'s factor, `LogDet` and `Status` equal this class's bit for bit where no pivot is floored.
     - An inverse writer that ports this class's L⁻¹ path keeps GMM's Σ⁻¹ bit for bit, where one rebuilt from solves would not.
     - It must refuse an inverse that overflows: `[1e-320]` has a finite factor and `LogDet`, but its reciprocal overflows.
     - Reading Bhattacharyya's mean term through `InverseQuadraticForm` would move its distances in their last bits.
   - *Recommendation.* Retire it once decision 1 removes the floor, and give `DenseCholesky` an inverse writer and an Lz sample so GMM loses nothing. Keep it if the floor stays: the floor is then the whole difference between the two.

**Doc changes owed.**
- *`docs/Maths/LinearAlgebra.md`:*
  - the API list (line 95): `CholeskyDecomposition` gains `FailedAt`;
  - "Dense Cholesky (`CholeskyDecomposition`)" (line 193): the refusal, its three triggers, and the reads that throw, before a first `Decompose` too. Also the floored factor that stays finite: it is read, and flooring can introduce order-one error, with the ladder's numbers; and the floor fires on scale alone.
- *`docs/Clustering/GaussianMixtures.md`:*
  - line 45: `UpdateCache` throws on a covariance with no finite factor, leaving Σ⁻¹ and the log-normalizer NaN;
  - the warm start (line 89): a singular warm-start covariance is floored and read, or refused when its factor is not finite;
  - "Bhattacharyya distance" (line 162): it refuses a covariance with no finite factor, by name.
- *The scene export's `cholesky_l`:* no doc in `docs/` describes it.

### 2026-11-13 — Slice 1: the weight reader at extreme weights

**Commit.** `f50753f` fix(regression): the weight reader keeps every finite weight's normalized log.

**What changed.**
- *`PredictiveChecks.Rows`* reads every finite weight to full precision.
  - Where the sum is finite and W_d/ΣW a normal double, it takes ln(W_d/ΣW), the expression it took before. Ordinary weights keep their bits by construction.
  - Where the sum passes `double.MaxValue`, it takes the same ratio of the weights scaled by the power of two that brings the largest into [1, 2). The scaling is exact.
  - Where the ratio is still below 2⁻¹⁰²², it takes ln W_d − e ln 2 − ln(the scaled sum), a difference of logs. These are the ratios the old form read as −∞ or as an imprecise subnormal.
- *`PointwiseLogDensity.CheckWeights`* tests for one positive weight instead of forming a sum that can overflow. It admits and refuses exactly what it did.
- *One shared routine.* WAIC, PSIS-LOO and the conformal layer's nine sites all read `Rows`. No other reader of the pointwise or candidate densities' weights normalizes them, so the conformal folder is unchanged.

**The gates.**
- *Slice 10's criterion*, in `WeightReaderTests`:
  - `[1e300, 1e-300]` reads rows [0, 1] with log weights [0, −1381.551055796427]; `[1e308, 1e308]` reads two weights of exactly ½. No NaN, no lost row.
  - The log weights against an exact reference: the weights summed as integers in units of 2⁻¹⁰⁷⁴, and the ratio taken to at least 62 bits by integer division. The tolerance is (n + 4)·2⁻⁵²·max(1, |reference|).
  - Nine chosen cases, and a seeded sweep of 4,000 vectors: 2 to 41 weights, binary exponents over random bands of [−1074, 1024), with zeros, the smallest subnormals and values within 10⁻³ of `double.MaxValue`. The worst error is 0.15 of the tolerance. 2,694 of the vectors hold a ratio the old form lost.
  - The readings, against log-sum-exp references on the exact log weights: the weighted predictive (`SplitConformal.LogPredictive`, WAIC's per-unit elpd), the CDF scores (`CdfScores.Fixed`), and slice 10's carried CDF error.
  - A finite tilt of about −1381 nats on the dominant draw makes the first input's tiny row material. It moves WAIC's lppd by 0.35 and 0.61 nats, and the carried bound's weights to .62 and .38.
  - Positive rescaling: 4,000 vectors scaled by a power of two. The 1,353 whose weights and ratios all stay normal agree bit for bit, and the rest within the tolerance. Scales that are not powers of two agree within three tolerances, and the carried error holds under 2^±20.
  - Zero weights drop out and leave the others' log weights bit for bit (2,000 vectors).
- *The reproduction.* On HEAD's code 12 of the 18 new tests fail, each extreme input reading −∞. After the repair, 18 of 18 pass.
- *The ordinary range.* The new pinned test compares 4,000 ordinary vectors with the old form, kept verbatim, bit for bit. The `[1, double.Epsilon]` receipts also keep their bits: that ratio is exact, and the difference of logs gives the same ln ε.
- *The whole regression assembly*, pinned tier included, from a `git archive HEAD` snapshot plus the three files: 3,306 passed, 0 failed (29 min). The selection by change was the whole assembly less its slow tier, since `Pointwise.cs` is a hub, so the slice's run covered it.
- Receipts are in `artifacts/shared-code-defects/receipts/` (`slice1-*`).

**Defects found.**
- None in the slice's code.
- *Finding, the same class outside the brief's files.* BASS's draw-weight readers normalize linear weights by a sum formed in linear space: `Draws3.Select` (`spline/bass/posterior/BassPrediction.cs`) and `PosteriorSummary.Over` (`spline/bass/presentation/`). Under `[1e308, 1e308]` the sum is +∞, and `Draws3.Select` returns weights of 0. They read `BassFit.DrawWeights`, not the pointwise densities. The samplers never produce such weights (SMC's and the flat histogram's are self-normalized), but `fit with { DrawWeights = … }` can.

**Decisions for the user.**
- *Whether BASS's draw weights get the same full-range reading.* If yes, a linear-space twin of `Rows`' rule serves `Draws3.Select` and `PosteriorSummary`, and prediction, sensitivity and the summaries read any finite weights. If no, `BassFit.DrawWeights` states that its sum must stay below `double.MaxValue`, and those readers refuse one that does not.

**Doc changes owed.**
- *`docs/Maths/Regression/Regression.md`*, the weighted draws (line 214): `PointwiseLogDensity.Weights` takes any finite, non-negative weights with one positive, over the whole finite range. The shared reader's rule: the ratio's log where it is normal, so ordinary weights read as before; the power-of-two-scaled ratio past `double.MaxValue`; a difference of logs below the normal range.
- *`docs/Maths/Regression/ConformalPrediction.md`*, "Weighted draws work as they are" (line 76): the shared rows read extreme weights to full precision. This closes slice 10's parked weighted-row follow-up.
