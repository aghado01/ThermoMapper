# E21 slice 1 cold-start hinge investigation

Date 2026-11-09. Commissioned by the coordinating chat after the [slice 1 report](ridge-blocks-brief.md#slice-1--drawn-scale-on-one-gaussian-intercept-grouping-2026-11-09). Source baseline `69398e4`. The [diagnostic tests](../../../tests/regression/RandomEffectsColdStartDiagnosticsTests.cs) and this report are separate from the BASS production files.

The trapped pair has accessible, correctly scored exits with a severe proposal barrier. At fixed s = 1.3 and ln τ = −4, the production kernel matches the enumerated transition row, and all 299 admissible structures are reachable. Solving that finite kernel's first-passage equations gives 1.11 million structure steps to reach either one-term linear representation. This explains how a state with posterior probability around 0.0002 can retain a prior-started chain over 210,000 sweeps. The original moving-nuisance chain's first-passage distribution remains unmeasured.

## The pair and its routes

The saved move/block seeds 303/403 trace occupied
`x0−k2 · x1+k0 | x0−k2 · x1−k2` at 418 of 420 checkpoints spaced 500 sweeps apart. With the fixture's knots 0, ½ and 1, its columns are (1−x₀)x₁ and (1−x₀)(1−x₁). Their sum is 1−x₀. The pair retains the response's strong linear x₀ signal, while either interaction alone loses it.

The default palette at the two-term cap selects death and change with probability ½ each. A death picks either term with probability ½, giving a proposal probability of ¼ for each reduced structure. From that reduced structure the reverse birth has probability 1/198, including selection of birth at M = 1. Their ratio is 2/99. Changes redraw one factor's sign and knot while preserving the term's variables and interaction order. Each distinct factor redraw has unconditional proposal probability 1/48 and an equal reverse probability. Birth has zero selection weight at the cap.

At s = 1.3, ln τ = −4, the better death has Δlog target = −10.21730 and log acceptance ratio = −14.11928, including the reverse-proposal correction. Its accepted probability per structure step is 1.84508 × 10⁻⁷. Both deaths together contribute 1.85706 × 10⁻⁷. Changes contribute 7.31610 × 10⁻⁵ to the probability of leaving the exact pair, but keep the chain among two-interaction structures. Many such excursions return to the pair.

The total proposal mass refused by support or numerical admission is 0.2083333. Unchanged redraws account for a further 0.0833333. These are included in the transition's self-loop. The remaining proposals are scored with the full selection and forward/reverse probabilities.

## Exact rows and production controls

The test reuses `BassKernelEnumeration` to enumerate candidates from the move definitions, sums paths to each atom, and obtains each reverse proposal by enumerating its reverse move. Production `Propose` calls then check candidate frequencies, refusals and log ratios. Finally, 150,000 production `ReversibleJumpChain.Step` trials per panel reset to the fixed source before each draw. They exercise the actual selection and acceptance code.

| Control | Enumerated probability of leaving the exact pair | Observed exits / 150,000 | Expected exits |
|---|---:|---:|---:|
| s = 1.3, ln τ = −4, b = 1 | 0.00007334669 | 11 | 11.002 |
| s = 0.05894, ln τ = −3.94330, b = 1 | 0.00023518215 | 30 | 35.277 |
| s = 1.11035, ln τ = −2.62045, b = 1 | 0.00069768022 | 104 | 104.652 |
| s = 1.3, ln τ = −6, b = 1 | 0.00002158428 | 3 | 3.238 |
| s = 1.3, ln τ = −4, b = 0.4 | 0.015074786 | 2,223 | 2,261.218 |
| s = 1.3, ln τ = −4, b = 1, `FactorMoves(1,1,0)` | 0.062536803 | 9,308 | 9,380.520 |
| Grouping removed, ln τ = −4, b = 1 | 0.00023171955 | 31 | 34.758 |

The two controls drawn from the saved trace use its rounded printed coordinates. Across the panels, the largest log proposal-ratio discrepancy was 2.22 × 10⁻¹⁶ and the largest target-difference discrepancy against a strict fresh fit was 2.13 × 10⁻¹⁴. All candidate frequency and acceptance checks passed. The fixed-coordinate posterior probability of the pair ranges from 0.0001116 to 0.0005268 in these controls, consistent with the earlier fully integrated enumeration's small mass.

## Connectivity and passage times

The full transition graph reaches all 299 states from the trapped pair. The largest detailed-balance discrepancy is 3.55 × 10⁻¹⁴ in log flow. First-passage times solve the killed kernel's equations through a reversible, symmetric Dirichlet matrix. The largest residual in the original passage equations is 1.55 × 10⁻⁹ against a right-hand side of one.

| Frozen control at ln τ = −4 | Mean steps to leave the sector with two second-order interactions | Mean steps to reach x₀ or 1−x₀ as a single term |
|---|---:|---:|
| Default, s = 1.3, b = 1 | 795,133 | 1,114,442 |
| Default, s = 1.3, b = 0.4 | 86.27 | 103.51 |
| `FactorMoves(1,1,0)`, s = 1.3, b = 1 | 16.00 | 20.49 |
| Grouping removed, b = 1 | 1,491,680 | 1,907,823 |

One default path goes through a rare death, births x₀ beside the surviving interaction, then deletes that interaction. Its first step has probability 1.84508 × 10⁻⁷. With the existing factor moves enabled, pruning x₁ from either interaction produces a main effect beside the other interaction and retains the same column span. Each of these two routes has forward probability 1/32, reverse probability 1/96 and acceptance one. Their prior and proposal ratios cancel. A subsequent death reaches a one-term linear representation. The no-group control shows that the geometric barrier also occurs in ordinary Gaussian BASS.

## Moving nuisance coordinates

Quadrature over the nuisance law conditional on the exact pair gives the following expected outgoing probabilities. The fixed-s control retains the drawn-scale model's frozen prior reference and integrates τ. The coarse grid uses ln s in [−18, 10], ln τ in [−35, 4] at spacing 0.2. The fine grid expands these to [−22, 12] and [−40, 5] at spacing 0.1. The largest change in log mass, log means or outgoing probabilities is 1.60 × 10⁻⁷.

| Conditional control | Leave the exact pair | Direct death out of the two-interaction sector | Prune |
|---|---:|---:|---:|
| Default, b = 1, s and τ integrated | 0.000119610626 | 3.96874516 × 10⁻⁷ | 0 |
| Default, b = 1, s = 1.3, τ integrated | 0.000094384328 | 5.12736585 × 10⁻⁷ | 0 |
| Default, b = 0.4, s and τ integrated | 0.032878465 | 0.0171607381 | 0 |
| `FactorMoves(1,1,0)`, b = 1, s and τ integrated | 0.0625601121 | 3.96874516 × 10⁻⁷ | 0.0625001084 |

These conditional averages price the available exits as s and τ move. The passage times above hold the nuisance coordinates fixed. They establish the bottleneck and its sensitivity to the proposal palette and temperature, while the original moving-nuisance chain still needs its own first-passage qualification if cold initialization is to be made reliable.

## Receipts and reproduction

The [saved receipts](../../../artifacts/e21-cold-start-qualification/) ([tracked inputs](../receipts/e21-cold-start-qualification/)) contain `diagnostics.trx`, `commit.trx`, both runsettings files and the source baseline. The diagnostic snapshot held `69398e4` plus the new test file. Builds used Release, below-normal priority, `DOTNET_PROCESSOR_COUNT=2`, `MSBUILDDISABLENODEREUSE=1`, `-m:1` and one xUnit worker. With `BASS_QUALIFICATION=1`, all eight diagnostic bodies passed in seven seconds.

```powershell
dotnet test <snapshot>/tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj -c Release -p:ArtifactsRoot=<private-build-root> -m:1 --disable-build-servers --settings <diagnostics.runsettings> --filter FullyQualifiedName~RandomEffectsColdStartDiagnosticsTests --logger 'trx;LogFileName=diagnostics.trx' --results-directory <receipts>
scripts/test-select.ps1 -SourceRoot <snapshot> -Scope Commit -Changed tests/regression/RandomEffectsColdStartDiagnosticsTests.cs -ArtifactsRoot <private-build-root> -Threads 2 -NoBuild
```

The Commit selection ran with the qualification flag unset and passed 44/44. Its eight diagnostic rows returned early and its pinned tests executed. The production files retain the slice 1 implementation during this investigation.
