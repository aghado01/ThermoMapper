# Depth-four completion by two four-site power laws

Commissioned 2026-11-17 by the user through the research coordinator. This is question 1 of the [Deep Teaching research brief](../deep-teaching-audio-research-brief.md). Complete the mathematical assignment in this brief, then the [predictive consequence](02-predictive-ambiguity-brief.md).

## Observation and model

Use stationary, sign-symmetric binary sources of memory at most four with strictly positive five-word probabilities. The baseline consists of separate population DC-power laws of rectangular widths one through five. Add the separate DC-power laws with indicator tapers at A={0,1,2,4} and B={0,1,3,4}. There is no correspondence between samples of different laws. Use the existing Fourier power and finite source definitions.

## Proof obligations

1. Prove that the augmented observations recover the five-word distribution, its transition kernel and the minimal suffix tree. Derive the two added fourth moments from the actual DC maxima and recover the third from the width-five maximum. Verify the finite inversion.
2. Check the proposed two-parameter family with quartic coefficients a, b, −a−b. Verify strict positivity under |a|+|b|+|a+b|<1, endpoint marginals, a stochastic shift-and-append kernel and its invariant history law. Identify histories witnessing exact memory four whenever a chosen collision direction is nonzero.
3. Prove equality of the full joint rectangular Fourier-power vector law at each separate width. Check cyclic invariance and the induced cancellation on observable fibers.
4. For every four-site indicator taper on the five sites, give a nonzero positive stationary family that also matches its entire DC-power law. This establishes the claimed minimum of two within the stated class, if the completion result holds.
5. Link declaration names, build commands and axiom checks. Record any correction to the conjecture and its consequences. End with Doc changes owed.

## Acceptance

The theorem must connect actual observation operators and stationary histories to kernel and tree recovery. A dimension count is supporting evidence. Historical novelty and extension to general speech remain unassessed.

## Result and proof report — 2026-11-17

**The proposed boundary holds.** Two separate added four-site DC-power laws suffice; one such law does not suffice, even with all rectangular Fourier powers retained jointly within each separate width. The minimum is two within the four-site indicator class on these five sites. This says nothing about arbitrary real tapers or paired measurements.

The positive proof uses four lag moments from rectangular power means, the consecutive quartic from the width-four law, and its translated copy. For a four-site set S, q_S = 8 Pr(Z_S=16) − 1 − Σ_{i<j in S} r_{j−i}. The A and B laws therefore give q_A and q_B. The width-five maximum supplies q_C through 16 Pr(Y_5=25) = 1 + 4r_1 + 3r_2 + 2r_3 + r_4 + 2q + q_A + q_B + q_C. Sign symmetry removes the odd Walsh moments. The finite inversion recovers all 32 word masses.

The algebraic block-equality theorem actually requires neither positivity nor a memory bound: it applies to the stationary, sign-symmetric cylinder laws in `Source`. Its consequence is equality of the four-history conditional ratio and the suffix tree defined from that ratio. Strictly positive five-word masses make the ratio an ordinary conditional probability (`ContextTree.kernel_denominator_pos`). The stated memory-at-most-four assumption is what makes this kernel and tree describe the whole source; it is not a claim that five-word observations determine arbitrary longer-memory sources.

For the lower bound, the proposed family has total mass one, strict positivity under the stated parameter restriction, sign symmetry and both four-word endpoint marginals equal to 1/16. The formal construction gives the shift-and-append transition kernel, its stochasticity, the identity between its stationary five-word joint law and `weight`, and invariance of the uniform history law against every test function. This is a finite-state stationary source construction; an infinite path-space measure is not introduced. Nonzero (a,b) has histories with the same last three signs and different next-positive probabilities, proving exact memory four. The fair kernel has only the empty context; a nonzero family member does not.

The full Fourier-vector cancellation follows from cyclic rotation of the five sites. Each quartic has the same sum on each rotation-invariant observable fiber, and their coefficients sum to zero. Proper contiguous subblocks are uniform by an omitted-endpoint cancellation. Thus every full rectangular Fourier-vector law through width five agrees with the fair source, separately at each width and offset.

A four-site indicator taper omits one site r. Its entire law depends only on the following coefficient:

| Omitted site r | Coefficient retained | A nonzero invisible direction (a,b) |
| --- | --- | --- |
| 0 or 4 | 0 | (ε,0) |
| 1 (C) | −a−b | (ε,−ε) |
| 2 (B) | b | (ε,0) |
| 3 (A) | a | (0,ε) |

Here 0<ε<1/2. Each direction is strictly positive, stationary, and exact depth four; it matches both the fair rectangular laws and the whole added taper law. Every four-element subset of the five sites has the omitted-site form. Conversely, Pr(Z_A=16)=(1+a)/8 and Pr(Z_B=16)=(1+b)/8 recover both family parameters.

### Declarations and verification

All new declarations are in `TimeFrequency.ContextTree.DepthFour`.

| Module | Main declarations |
| --- | --- |
| [DepthFourCompletion.lean](../../../../../lean/Enthymemata/TimeFrequency/DepthFourCompletion.lean) | `taper_is_power`, `tapered_fourth`, `five_top`, `five_inversion`, `lags_identified`, `quartics_identified`, `completed_block`, `completed_kernel_tree` |
| [DepthFourStationary.lean](../../../../../lean/Enthymemata/TimeFrequency/DepthFourStationary.lean) | `positive`, `symmetry`, `endpoint_masses`, `total_mass`, `stochastic`, `kernel_from_block`, `stationary`, `sign_symmetric_kernel`, `proper_laws`, `prediction_formula`, `exact_four` |
| [DepthFourSharpness.lean](../../../../../lean/Enthymemata/TimeFrequency/DepthFourSharpness.lean) | `rotate_power`, `full_vector_law`, `all_fourier_laws`, `taper_law`, `taper_maximum`, `two_parameters_recovered`, `one_taper_insufficient`, `four_site_shape`, `fair_tree`, `collision_tree` |

Native Lean checks passed for the three modules. The final combined `lake build` passed **3,751 jobs** with both slices registered in the TimeFrequency aggregate; all 35 public axiom listings for this slice contain only `propext`, `Classical.choice` and `Quot.sound`. The four new modules have no warnings in the final run. `lean/scripts/meta-ci.ps1 -NoBuild` passed: all four new modules are proof-closed, with no apologies or declared axioms. Lean ran BelowNormal with `LEAN_NUM_THREADS=2`. The modules remain in Enthymemata pending promotion review.

### Doc changes owed

Update the approach's measurement boundary with the exact depth-four design and its restricted minimality class. Keep population identification separate from stability, finite-sample estimation, and general speech claims. The [prediction slice](02-predictive-ambiguity-brief.md) owns the consequences for a causal next-sign target.
