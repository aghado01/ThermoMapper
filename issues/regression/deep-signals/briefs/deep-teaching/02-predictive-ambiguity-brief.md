# Predictive consequences of the depth-four measurement completion

Commissioned 2026-11-17 as question 2 of the [Deep Teaching research brief](../deep-teaching-audio-research-brief.md), following [measurement completion](01-depth-four-completion-brief.md).

## Information and loss

Training supplies exact, separately labeled population laws of the stated windows. Deployment supplies the signed length-four history h. The source has the uniform history law in the collision family. The target is the conditional probability of the next sign being positive. Evaluate the history-averaged excess squared probability loss, E_i(t)=Σ_h 2⁻⁴(t(h)−q_i(h))². For a Bernoulli target this is Brier excess when the score is (t−Y)².

## Proof obligations

1. Derive the conditional probability and exact squared prediction distance for two parameters in the stationary family proved in slice 1. Quantify the residual directions under each single-taper augmentation.
2. Prove the sharp minimax excess for two specified observationally identical sources, with the midpoint predictor attaining it. Instantiate the constants on the physical Fourier collision.
3. Characterize the family’s linear predictive functionals invariant under the omitted directions. Include average next-positive probability and prediction given shorter suffixes. Specify deployment history and training information separately.
4. Prove that both added tapers remove the population ambiguity throughout the positive depth-four class. Explain the resulting zero ambiguity term under exact summaries. Finite-data estimation and stability require further results.
5. Report declarations, builds, axioms and Doc changes owed. Formulate concrete subsequent stability, teacher-transfer and sampling questions from the established constants.

## Acceptance

The main result quantifies loss due to the actual stationary Fourier collision and what the completed design repairs. The generic two-point identity supports that result. Retaining ordered training recordings or pairing observations changes the information contract.

## Result and proof report — 2026-11-17

Training receives the exact separate population power laws specified in slice 1. Deployment receives the signed history h=(h₀,h₁,h₂,h₃), uniformly distributed on the 16 histories in this family. Let H_r(h) be the product of the three history signs except position r. The stationary next-positive probability is

\[
q_{a,b}(h)=\tfrac12\{1+aH_3(h)+bH_2(h)-(a+b)H_1(h)\}.
\]

The three characters are orthogonal under the uniform history law. For Δa=a−c and Δb=b−d,
\[
\|q_{a,b}-q_{c,d}\|_{L^2(\pi)}^2
=\frac{\Delta a^2+\Delta b^2+(\Delta a+\Delta b)^2}{4}.
\]
Thus the exact two-source minimax Brier excess is the numerator divided by 16. The midpoint of the two conditional probabilities attains it. Admissibility makes both conditionals and their midpoint strictly between zero and one. The Brier convention here is the single Bernoulli score (t−Y)²; summing scores over both labels would double these constants.

### Cost on the actual measurement collision

For each possible added four-site taper, use its invisible direction v from the [completion report](01-depth-four-completion-brief.md#result-and-proof-report--2026-11-17).

- The fair source and parameter εv have exact two-source minimax excess ε²/8.
- The parameters +εv and −εv have exact two-source minimax excess ε²/2, attained by the fair predictor 1/2.
- For 0<ε<1/2 these are positive stationary sources with the same rectangular Fourier-vector laws and the same whole added taper law.
- For any predictor and any c<1/8, such a pair has worst-source excess greater than c. Conversely, the fair predictor has excess strictly below 1/8 at every admissible parameter in the entire family.

Consequently the minimax value over the family's fiber matching the fair observations is **1/8**, for the rectangular baseline and also after any one four-site indicator taper. This is a supremum over the open positivity region, not an attained maximum at a positive source. It is a lower bound for the broader stationary depth-four class; no matching upper bound over that entire class is claimed. The theorem `observed_excess_lower` includes the actual matching taper and Fourier-vector laws, rather than leaving indistinguishability as a hypothesis.

### Which targets remain visible

For a fixed real function f on histories, set F_f(q)=E[f(H)q(H)], u=E[f(H)(H₃−H₁)] and v=E[f(H)(H₂−H₁)]. Then
\[
F_f(q_{a,b})=F_f(q_{0,0})+(au+bv)/2.
\]
This gives a complete classification of **linear** predictive functionals on this two-parameter family:

| Retained information beyond rectangular laws | Necessary and sufficient condition |
| --- | --- |
| None; or a taper omitting site 0 or 4 | u=v=0 |
| A, omitting site 3 | v=0 |
| B, omitting site 2 | u=0 |
| C, omitting site 1 | u=v |
| Both A and B | Every f is identified |

The whole single-taper histogram is equivalent, on this family, to its coefficient: one direction follows from `taper_law`, the reverse from `taper_maximum`. This justifies using coefficient equality in `one_taper_functionals`. The separate `single_taper_functional` theorem classifies constancy along the selected invisible line only.

The unconditional next-positive probability is 1/2. Averaging over the oldest sign conditional on the last three signs also gives 1/2, despite nontrivial prediction from all four signs. The proof therefore distinguishes a target visible at a shorter suffix from the deeper signed-history target.

With both added tapers, `completed_kernel_tree` gives pointwise equality of the conditionals throughout the stated positive depth-four class. `augmented_prediction` makes the resulting zero ambiguity excess explicit. This removes population ambiguity; it does not prove finite-data accuracy, stability under approximate summaries, or correct teacher transfer outside the model.

### Declarations and verification

All declarations are in `TimeFrequency.ContextTree.DepthFour`, in [DepthFourPrediction.lean](../../../../../lean/Enthymemata/TimeFrequency/DepthFourPrediction.lean), with the source construction in [DepthFourStationary.lean](../../../../../lean/Enthymemata/TimeFrequency/DepthFourStationary.lean).

| Obligation | Declarations |
| --- | --- |
| Physical conditional and valid probability range | `prediction_formula`, `prediction_range`, `midpoint_range`, `brier_excess` |
| Exact prediction distance and two-point minimax | `prediction_distance`, `midpoint_risk`, `unavoidable_excess`, `stationary_minimax` |
| Collision cost and sharp supremum | `one_taper_predictive_cost`, `symmetric_pair_cost`, `sharp_eighth_lower`, `observed_excess_lower`, `fair_predictor_upper` |
| Visible targets | `suffix_prediction`, `functional_formula`, `average_prediction`, `identifiable_linear_functional`, `single_taper_functional`, `one_taper_functionals` |
| Completed measurement design | `augmented_prediction` |

The final combined `lake build` passed **3,751 jobs** with no warnings in the four new modules. This module's 19 public axiom listings contain only `propext`, `Classical.choice` and `Quot.sound`; combined with slice 1, there are 54 checked public listings. `lean/scripts/meta-ci.ps1 -NoBuild` passed with all four modules proof-closed. Lean ran BelowNormal with `LEAN_NUM_THREADS=2`. The proofs remain in Enthymemata pending promotion review.

### Doc changes owed

Update the approach with the training/deployment distinction and the exact ambiguity cost, keeping the supremal family result separate from any full-class upper bound. The [research brief](../deep-teaching-audio-research-brief.md#concrete-follow-up-statements-from-the-depth-four-slices) records concrete, unproved stability, teacher-transfer and sampling questions 3–5. No engine or probe was commissioned.
