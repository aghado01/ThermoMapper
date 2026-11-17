# Which context-tree memories are identifiable from separate window powers?

Proposed Lean slice, 2026-11-16. Can separate window-power distributions identify which past events affect the next event, even when they cannot recover the signs of individual events? The assignment is to establish or refute a sharp cutoff in memory depth. It informs which measurements to retain before fitting a context tree.

## Experiment and claim

Fix a maximum memory depth D and set L = D+1. Events occur at known ordered sites and have unit amplitudes with signs Sₜ ∈ {−1,+1}. The source is stationary, invariant in law under simultaneous sign reversal, has strictly positive L-word probabilities, and has a transition law of order at most D. This includes variable-length context trees up to depth D. Define the **minimal tree** by taking, for each full history, its shortest suffix on whose set of extensions the next-sign law is constant. Positivity makes every such history observable in the source law.

For each width w from 1 through L, retain the exact distribution of the unnormalized rectangular DC power

$$Y_w=(S_1+\cdots+S_w)^2.$$

Widths are labeled. Within a width we know the entire power distribution. Across widths and successive blocks we retain no pairing or time order. Stationarity allows any fixed start index. This is a controlled source class and measurement family, distinct from the pilot's Hann cells.

**Candidate theorem.** These observations identify the transition kernel, and hence the minimal context tree, throughout this class when D ≤ 3. For every D ≥ 4 there exist a fair IID source and a source of minimal memory depth exactly D with identical observations. For that pair, even the joint vector of *all* rectangular Fourier powers at any one width, delivered separately for each width up to L, has the same law. The cutoff concerns widths through D+1 and this unpaired observation protocol.

## Positive direction — reconstruct the block law

1. Write M_w = E[Y_w], M₀ = 0 and rₖ = E[S₀ Sₖ]. Expanding the square yields M_w = w + 2 Σ_{k=1}^{w−1}(w−k)rₖ. Second differences recover each rₖ for k < L.
2. For D = 3 let q = E[S₀S₁S₂S₃]. The event Y₄ = 16 means all four signs agree. Show that 8 P(Y₄ = 16) = 1 + 3r₁ + 2r₂ + r₃ + q. Thus the *law* of the longer-window power supplies an interaction its mean omits.
3. Expand the finite block law in sign products χ_A(s) = ∏_{i∈A}sᵢ. Odd moments vanish by sign symmetry. Through length four, the constant, the pair moments at their respective lags, and (at length four only) q exhaust the even moments. Use finite inversion to recover every word probability. Ratios of word and prefix probabilities give the kernel, then the suffix rule gives the minimal tree. Include the empty-history case D = 0.

## Negative direction — exhibit two different memory structures

For each L ≥ 5, use zero-based coordinates. Set A = {0,1,2,L−1}, let A′ be its reversal under i ↦ L−1−i, and choose 0 < ε < 1/2. Compare the uniform L-word law with

$$p_\varepsilon(s)=2^{-L}\big[1+\varepsilon(\chi_A(s)-\chi_{A'}(s))\big].$$

1. Prove normalization, strict positivity and sign symmetry. Both products contain both endpoints, so summing over either endpoint leaves the uniform (L−1)-word law. **Realize the witness as a stationary source.** Put the uniform law on length-D histories and use K(b∣h) = 2ᴰ pε(hb) on the shift-and-append chain. Verify stationarity and the claimed block probabilities. For exact memory D, fix the interior signs to +1 except at index 1. Flipping the oldest sign then changes the probability that the next sign is +1 while preserving the last D−1 signs. Fair IID has the single empty context, so the two minimal trees differ.
2. Every proper contiguous subblock has the uniform law. At full width the perturbation is antisymmetric under time reversal. Show that summing it on a fiber of any reversal-invariant observable gives zero. For real rectangular frames, reversal conjugates each Fourier coefficient up to a unit phase, so its entire power vector is reversal-invariant. This proves observational equivalence with fair IID under the stated protocol.
3. **A measurement that sees the difference.** Retain a DC cell with unit weights at the four sites in A and zero elsewhere. Its power reaches 16 when those four signs agree. Derive probability 1/8 under fair IID and (1+ε)/8 under the perturbed law. This one asymmetric taper distinguishes the constructed pair. Uniform recovery of all depth-D trees with an augmented design would be a further theorem.

## Relation to the literature and existing Lean

- Kontoyiannis et al. give exact evidence, prediction and tree selection from an observed discrete sequence. Those results provide the source model. The observations here discard the signs before inference. [Bayesian Context Trees](https://arxiv.org/abs/2007.14900)
- The posterior representations for Bayesian context trees quantify uncertainty about branches given a symbol sequence. If the proposed collision holds, independent samples from the prescribed width-specific power laws have likelihood ratio one between the two witness sources at every sample size. A prior favoring either source would retain those odds. [Posterior Representations for Bayesian Context Trees](https://arxiv.org/abs/2202.02239)
- Dumont proves strong consistency of penalized tree selection, up to hidden-label permutation, under identifiability and regularity assumptions for a variable-length hidden Markov model. Its likelihood uses an ordered sequence of observations emitted independently conditional on the hidden sequence, each from its current state. Identifiability under the specific window operator is the question posed here. [Dumont, model and main theorem](https://arxiv.org/abs/1109.0392)
- Collet, Galves and Leonardi bound the effects of small independent symbol flips and recover finite truncations of a variable-memory tree under their continuity and non-nullness conditions. Here a deterministic power can erase memory even with exactly known observation laws and zero measurement noise. [Random perturbations of stochastic chains](https://arxiv.org/abs/0707.2796)
- Reversal ambiguity also occurs in reconstruction of individual strings from substring compositions. The object here is a *stationary source law* with one population distribution per width. The scoped literature check is not a claim of publication priority for the proposed cutoff. [Acharya et al.](https://arxiv.org/abs/1403.2439)
- Reuse the project's finite [Fourier-power definitions](../../../../lean/Enthymemata/TimeFrequency/TwoSampleRecovery.lean), [finite power measurements](../../../../lean/Enthymemata/TimeFrequency/FinitePowerMeasurements.lean) and [unpaired-view incidence criterion](../../../../lean/Enthymemata/TimeFrequency/UnpairedModels.lean) where they fit. The missing implication in the incidence argument is whether a direction hidden by the measurements can change memory while still defining a positive, stationary source.

## Acceptance and independent checks

Close the slice with the low-depth recovery theorem, the all-D high-depth source construction, and the asymmetric-cell witness, or record the first counterexample to the proposed boundary and the corrected scope. A stationary finite history chain with a verified invariant law suffices to realize the source. Link the named declarations, build and axiom checks in this brief under the [Lean harness conventions](../../../../lean/README.md).

Drafting checks used exact integer enumeration. For each D = 1,…,7, group nonempty even subsets of an L-word by translations, sum χ_A within each group, and sum those columns over each labeled DC-power bin. The resulting ranks were 1, 2, 4, 6, 9, 12, 16 against 1, 2, 4, 8, 16, 32, 64 columns. For D = 4,…,7 the proposed perturbation summed to zero on every within-width autocorrelation-vector fiber, and shifted the asymmetric-cell event by ε/8. These checks support the posed statements. The general Lean proofs remain open.
