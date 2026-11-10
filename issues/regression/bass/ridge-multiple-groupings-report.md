# E21 slice 2 — named crossed and nested Gaussian intercepts

Work commissioned by [the ridge-block brief](ridge-blocks-brief.md), slice 2. The public API follows the user's decision to migrate the singular API to named collections. The shared runtime was reconciled against E34's committed handoff `41727ac`; qualification uses a clean archive of that commit with the E21 source overlay. Concurrent conformal work and other sessions' notes are outside that source snapshot.

## Implemented shape

- `BassData.Groupings` holds named categorical columns beside the predictors. `BassConfig.RandomIntercepts` names the columns each block reads. A grouping can serve only as a parent. Missing, repeated and unused names are refused.
- `ObservationCentred` enforces Σ n_g u_g = 0 over a grouping. `WithinParent` names a parent and enforces that sum separately inside each parent level. A child level spanning multiple parents is refused. Singleton children have zero contrast inside their parent.
- A single observation-centred block delegates to the original closed form, including its residual-sum arithmetic and random draws. Several blocks share one scale-free thin factorization of their concatenated constrained indicators. Each scale update factors only the reduced weighted matrix R·S.
- Each drawn scale has its own step in ln s, after τ. The target sums each half-Cauchy log density and log-coordinate Jacobian. Prior samplers draw each declared scale independently before admitting the structure on that reading.
- `BassState.GroupScales` reads the fitted scale vector. `BassDraw.Groups` holds named `BassGroupDraw` records with scale and effects. Prediction, conditional means, pointwise readers and gradients use the collections; retained records resolve by name even when reordered.
- Exact crossed aliases are logged, including aliases with a within-parent child in a mixed layout. The full latent contrast dimension is retained under the reduced observation rank, with no additional alias constraint.

The model and migration examples are in [RandomEffects.md](../../../docs/Maths/RandomEffects.md). These blocks take a scalar Gaussian response. Slopes, moving noise weights and the conformal candidate density remain later slices.

## Numerical review and correction

The initial reduced Gram implementation passed the ordinary dense algebra and law checks. A separate review added a stable analytic two-contrast reference for substantially unequal crossed scales. At `(10⁶, 0.6)` and `(0.9, 10⁶)`, forming R S² Rᵀ and applying S²RᵀK⁻¹ lost the small component: maximum whitening/mean gaps were `7.04284e-6` and `1.091e-5`. The original `1e-9` bound failed in both orders; the failing receipt is retained.

The correction takes a reduced SVD `R S = U diag(d) Vᵀ`, avoiding the squared condition number. The conditional mean applies `S V diag(d/(1+d²)) UᵀQᵀr`. In standardized contrast coordinates the conditional square root is `I + V diag((1+d²)^(-1/2)-1) Vᵀ`. Drawing all q standard coordinates preserves the q−r latent complement under aliases.

| Scales | Whitening/mean gap after correction | Conditional covariance gap | Largest mean/covariance z over 50,000 actual draws |
| --- | ---: | ---: | ---: |
| `(10⁻⁶, 0.6)` | 2.220e-16 | 6.939e-17 | 1.242 |
| `(10⁶, 0.6)` | 4.441e-16 | 1.554e-10 | 1.242 |
| `(0.9, 10⁶)` | 3.331e-16 | 2.055e-10 | 1.348 |

Both deterministic gates retain the `1e-9` bound. The covariance reference inverts the independent two-contrast precision; the ordinary full-case references construct level-space projectors directly from counts and parent membership, without the production contrast bases or thin factors. Ordinary, nested, exact-alias and mixed-alias cases also check every covariance entry and each draw's constraint.

## Laws, recovery and retained state

The finite enumeration uses three admissible structures with substantial posterior mass on each, a held τ and both scales drawn. Gauss–Legendre quadrature integrates each half-Cauchy through its full prior probability coordinate, covering `(0,∞)` without a log-scale cutoff. Refining 64 to 96 nodes per scale changed ln Z by at most `3.553e-15`, and mean log scales by at most `3.544e-4`. The chains meet structure probabilities and both mean log scales within the specified batch-mean errors.

| Layout | Dense ln Z | SMC ln Z ± reported MC error | SMC flags |
| --- | ---: | ---: | ---: |
| Crossed | −18.134866 | −18.120279 ± 0.009880 | 0 |
| Nested | −17.598988 | −17.598812 ± 0.006705 | 0 |

Four independent Geweke rows each use 12,000 prior/successor pairs: crossed and nested layouts, two seeds each, with τ and both scales drawn. Functions cover both log scales, both log effect spreads, M, ln τ and ln σ². The largest absolute paired-difference z is `2.672` against the `4.5` gate.

Recovery uses four datasets per layout, 320 observations and a nonlinear surface, with unbalanced groups. One effect law is two-point and the other shifted exponential. The fitted prior is Gaussian and both scales are unknown. Child means are assigned to their parents before generation under the within-parent convention. Held-out prediction uses 400 cases of observed levels.

| Layout | Effects covered by 90% intervals | Standardized error mean square | Scale intervals covering the declared target | σ intervals covering 0.35 |
| --- | ---: | ---: | ---: | ---: |
| Crossed | 45/52 | 0.989 | 8/8 | 3/4 |
| Nested | 58/64 | 1.066 | 8/8 | 4/4 |

Under these non-Gaussian generators the scale target is explicitly the realized constrained contrast RMS divided by σ. These counts do not establish population-variance coverage under arbitrary misspecification. In every dataset the fitted level is within two posterior standard deviations of the convention's target; held-out grouped RMSE is less than half the fit without blocks, and the population surface satisfies the permitted 10% margin against that fit.

Plain and reference-rung ladder runs resume bit for bit with two scales and both effect arrays retained. Constraints hold to `1e-12` in every checked draw. Reordering retained named records keeps fitted readings and values. Population predictive spread includes both missing groupings' draw-specific scales.

## Migration and E34 reconciliation

Old `RandomIntercept` JSON is rejected. Convert its object to a one-entry `RandomIntercepts` array with `Grouping`, and pass a grouping collection to `BassData`. Read each named retained record's `Effects` and `Scale` in place of the old flat fields. [The conversion example](../../../docs/Maths/RandomEffects.md#migrating-a-saved-one-group-configuration) includes the checkpoint consequence: grouped JSON changes the fingerprint, so converting a config alone cannot resume an old grouped checkpoint under the existing guard. New named checkpoints resume exactly.

The ungrouped fingerprints remain `77a50a29e3922982` and `5663b1f7fddc270e` for their pinned fixtures. All three existing held-scale draw hashes remain unchanged. E34's `LogOutputRates`, `OutputRateArray`, propagation through state transformations and structure proposals, `WithOutputLogRate` and `BassDraw.OutputRates` are retained. Its `EverySamplerCarriesTheOutputDrawsAndResumesExactly` cases pass on the combined source.

## Cost and qualification receipts

At 512 cases, two crossed groupings and q = r, the indicative timings below were measured at below-normal priority with two visible .NET processors. They are cost rows, not portable performance bounds. Scale factorization costs `O(qr²+r³)` and column whitening `O(nr+r²)`.

| Retained rank r | Scale update, µs | Column whitening, µs |
| ---: | ---: | ---: |
| 2 | 2.66 | 4.62 |
| 6 | 8.20 | 9.31 |
| 14 | 26.41 | 50.68 |
| 30 | 216.17 | 67.54 |

The [receipts](../../../artifacts/e21-ridge-s2-qualification/) include the baseline commit, a SHA-256 manifest matching the tested C# files to main, run settings and saved TRXs. Runs use the isolated source snapshot and two visible .NET processors, below-normal priority and disabled node reuse. Qualification and Commit use one xUnit thread; full Slice uses two with qualification off, following slice 1's signed-off protocol. A serial direct build prepares the assembly; the broad tiers use the selector's generated filter and its supported `NoBuild` mode.

| Gate | Result | Receipt |
| --- | --- | --- |
| Combined qualification, `BASS_QUALIFICATION=1` | 66/66 passed | `combined-qualified.trx` |
| Held/drawn slots and second-block histogram refusal | 1/1 passed | `held-drawn-slots.trx` |
| Commit selection plus pins | 2,037/2,037 passed | `Commit.trx`, `commit.runsettings` |
| Current-main test assembly compile after concurrent work landed | Passed, zero errors, at `5ed8bd1` plus E21 | `main-build.log`, `main-build-baseline.txt` |
| Full Slice assembly, qualification off, two xUnit threads | 2,874/2,874 passed | `Slice.trx`, `slice.runsettings` |
| Fresh review | Unequal-scale and mixed-alias defects corrected; deterministic covariance and actual-draw checks passed | `unequal-scales.trx`, `stable-whitening.trx`, `combined-qualified.trx` |

Early failed receipts preserve test-assertion corrections and the actual unequal-scale numerical failure; they are not substituted for the final gates. The prior cold-start investigation remains [its own bounded result](ridge-cold-start-investigation.md).
