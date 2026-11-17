# Deep signals: the approach, what is established, and what is missing

Date: 2026-11-16. Status: **a direction record, written at the user's request** so the why and the what of this workstream survive the session. It assigns no work. Each item under [What is missing](#what-is-missing-in-build-order) gets its own brief when the user takes it up. Drafted by the coordinating session from one day's discussion with the user, two chips' reports, and three source collections.

## Why

The aim is a **native audio foundation model trained on human conversation, without labels, and data-efficiently.** The regression engines cannot be that model: they sample shallow, low-order structure exactly, and audio needs deep composition at a scale they cannot touch. Their part is different and specific:

- **The exact teacher at small scale.** A posterior over a small latent under a known likelihood, where a network would only have a point estimate.
- **The measuring stick.** Evidence in bits, and coverage that holds without trusting the model.
- **The place where missing mathematics gets built soundly first,** with Lean statements as specifications and pre-stated gates on numerical claims.

The user's standing direction for the thread: the data is human conversation; arguments are from first principles, and what a field does or omits is context, never a refutation; the engine comes first and grows by coherent, additive pieces.

## The approach

**Three roles.**
- *Physics, known.* The measurement operators (windows at several lengths, magnitudes) and their noise law.
- *The prior, learned.* A generator from a low-dimensional code to the signal.
- *The solver, with guarantees.* It recovers the code from a lossy view. Voroninski's theorems say when that answer is unique, stable and reachable, for random expansive networks.

**The signal model.**
- *Three geometries in the time–frequency plane:* lines along time (tones), lines along frequency (transients), and filled area (noise). The first two are atoms. The third is described by a smooth spectral surface.
- *Speech is the atom model in physical form:* glottal pulses driving vocal-tract resonances. Between excitations a resonator obeys Lu = 0, and an excitation starts a new piece. That is an L-spline whose knots are the innovations (Schumaker, Chapter 10).
- *Conversation adds structure that is constant by construction:* one vocal tract per speaker, one room and microphone per recording. Prosody is the slow layer; turns, gaps and overlaps are the slowest events, between two streams.

**Representation principles, each argued from first principles.**
- Time resolution should scale with frequency, with a floor: bandwidth Δf(f) ≈ Δf₀ + f/Q. Three arguments give the scaling (the period bound on an envelope, the scale covariance of the wave equation, stability to small warps); harmonic sources give the floor.
- Tokens should follow the rate of innovation, not the bandwidth. Timing precision is cheap as an attribute of an event and expensive as a grid rate.
- A tiling moves information between magnitude and phase. The right one puts what matters in magnitudes or in atom parameters.

**Supervision without labels,** of three kinds:
- *Consistency.* Complex frame coefficients must reproduce themselves under synthesis with a dual frame followed by analysis. The Gram matrix itself is this consistency map for a Parseval frame. For a stationary signal, each window smooths the expected power spectrum by its own spectral kernel. Individual magnitude views can lose the phase needed to predict another view. For example, (1, 1) and (1, −1) have identical length-one magnitudes, while their unit-energy length-two rectangular DC powers are 2 and 0. Speaker identity across that speaker's turns and a stable recording channel are additional modelling assumptions.
- *A teacher.* A solver, or the engine's posterior, wherever the guarantees hold.
- *Synthesis.* The prior's own samples carry exact labels for an encoder.

The recovery theorems do not perform these. They certify the tasks: they say when a lossy view determines the latent, which is when the encoder's target is well defined.

**Conformal-Bayes.**
- A held-back view supplies its own label. A fixed consistency score can be calibrated on held-out exchangeable units to give marginal coverage for the corresponding observable on a new unit.
- Synthetic data also supplies the true latent. Transferring its coverage to real audio by covariate-shift weighting requires common support in the direction of transfer, the same conditional law of the latent given the view, and the target-to-source view density ratio for the exact weighted guarantee. Whether the simulator and real speech meet those assumptions is open. [Tibshirani et al., *Conformal Prediction Under Covariate Shift*](https://arxiv.org/abs/1904.06019)
- Choose the exchangeable unit to match deployment. For inference about a new speaker, split by speaker and account for speakers shared across conversations.
- Turn-taking is a risk-control problem: start speaking before the other has finished at most α of the time.
- Test martingales give change detection with false-alarm control.

## What is established

How far each row was checked is in the last column. "Elaborated" means the coordinating session re-elaborated the statements and printed their axioms; "gates checked" means it diffed the pre-stated gates against their final text and matched headline numbers to receipts; "chip's report" means neither was done.

**Proved in Lean** (no `sorry`, standard axioms). The first nine modules were promoted to `lean/Lemmas/` on 2026-11-16 (`f0cc4e2`). The noise-law and generative-prior modules are in `lean/Enthymemata/`. Paths below are relative to the stage.

| Result | Module | Checked |
|---|---|---|
| Edge products and magnitudes pin a signal up to one phase exactly when the graph is connected; polarization supplies the edge products from magnitudes | `TimeFrequency/PhasePropagation.lean` | Elaborated |
| A unit-energy window of length w reads \|A\|²·w from a tone, \|a\|²/w from an impulse, σ² from white noise | `TimeFrequency/WindowScaling.lean` | Elaborated |
| For symmetric 2×2 H, ‖SᵀHS‖²_F ≥ 2(det S)²\|det H\|, attained at every cell area | `TimeFrequency/FlatnessBound.lean` | Elaborated |
| A regularity inequality on a ball gives linear convergence to the noise level | `TimeFrequency/Contraction.lean` | Elaborated |
| The marginal under noise σ²(I ⊗ R), with the prior on the whitened Gram, is the unit-noise marginal on whitened data plus a determinant term; the Khatri–Rao Gram under Q; the level under Q | `BassWhitenedMarginal.lean` | Elaborated |
| A Kronecker design's Gram and projection factor by axis | `BamrsKroneckerGram.lean` | Elaborated |
| A design's full column rank is equivalent to a greedy ordered matching of sites to supports, given Theorem 4.61 as a premise (met at order 1) | `Spline/DesignRank.lean` | Chip's report |
| For order-2 hat splines with nondecreasing knots and strictly increasing sites, full column rank is equivalent to the greedy support matching. S16 proves the sufficiency premise for this basis, including repeated knots under the right-continuous convention | `Enthymemata/Spline/HatDesignRank.lean` | Clean archive of `e0b70c4` built and passed all ten axiom checks and the taxonomy check. [S16 report](../../../lean/noise-law-statements-brief.md#2026-11-16--s16-order-2-hat-design-rank) |
| The order-2 trigonometric B-spline is positive while its support is under half the carrier's period, and the condition is sharp; the damped oscillator's Green's function gives a one-sided basis | `Spline/TrigonometricBSpline.lean` | Chip's report |
| Free knots approximate any monotone f by steps within (f(b) − f(a))/(2(k + 1)); uniform knots leave some f at ½ for every k | `Spline/FreeKnotSteps.lean` | Chip's report |
| For a unit jointly circular complex Gaussian pair with squared coherence κ, the covariance of its log powers is Li₂(κ). The joint Laplace transform gives the series directly for arbitrary correlation phase. At perfect coherence the log-power variance is π²/6 | `TimeFrequency/CircularGaussianPowers.lean` and its three supporting modules | Clean archive of `2d9494c` built and passed all 26 axiom checks and the taxonomy check. [S10c proof report](../../../lean/noise-law-statements-brief.md#2026-11-16--s10c-gaussian-log-power-covariance) |
| Cov(\|a\|², \|b\|²) = \|E[a·conj b]\|² + \|E[a·b]\|², and Var(xᵀHx) = 2·tr((HΣ)²), each with the Gaussian pairing rule as a named premise | `TimeFrequency/LogPowerCovariance.lean`, `TimeFrequency/QuadraticFormVariance.lean` | Chip's report |
| ∫₀^∞ (1 − (1 − e^{−x})ⁿ) dx = H_n; the rectangular window's turnover is the root of eᵘ = 1 + 2u, between 1.25 and 1.26 | `TimeFrequency/PeakOfExponentials.lean`, `TimeFrequency/WindowTurnover.lean` | Chip's report |
| The response's second moments under the g-prior on the whitened Gram are the covariance the whitened marginal reads | `BassWhitenedMarginal.lean` | Chip's report |
| A range-restricted isometry bounds the squared norm of every generator secant and the compressed objective. At ε < 1 it preserves injectivity on each set where the generator is injective | `GenerativePrior/IsometryCondition.lean` | Elaborated |
| For isotropic raw second moments σ²I, the expected squared norm after a fixed rank-k orthogonal projection is σ²k, against σ²n in the ambient space | `GenerativePrior/IsotropicProjection.lean` | Elaborated |
| The exact ReLU angle map is strictly increasing and 1-Lipschitz on [0, π]. Its iterates decrease to zero from every positive starting angle and stay at zero when started there. A one-step angular error at most δ accumulates to at most dδ after d layers | `GenerativePrior/AngleMap.lean` | Elaborated, all 25 axiom prints checked |
| For depth d ≥ 2 and nonzero codes, the deterministic linear comparison field vanishes exactly at x★ and −ρ_d x★. Its antipodal radius satisfies ρ₂ = 1/π, ρ_d > 0 and ρ_d → 1 | `GenerativePrior/ComparisonField.lean` | Elaborated, all 19 axiom prints checked |
| For two planar unit directions at angle θ ∈ [0, π], the second moment of a standard Gaussian gated on both positive projections is Q(θ). Its cross term and the two output norms give the S17 angle map. The derivation includes the Gaussian radial integral | `GenerativePrior/WeightCondition.lean`, `GenerativePrior/GaussianWeightMoment.lean` | Elaborated both modules, all 17 axiom prints checked |
| For n ≥ 1 central hyperplanes in ℝᵏ, the number of open regions is at most 2·Σ_{j<k} C(n − 1, j), including repeated hyperplanes. Each component of the complement realizes one strict sign pattern and the ReLU layer agrees with a single linear map there. Coordinate hyperplanes attain 2ᵏ regions when n = k | `GenerativePrior/SignPatterns.lean` | Elaborated, all 18 axiom prints checked |

**Measured,** under gates written before their readings. The receipts are under `issues/regression/receipts/` in `latent-recovery/`, `spline-space/` and `generative-prior/`.

| Finding | Checked |
|---|---|
| For unit-energy window cells of Gaussian white noise whose coefficients are jointly circular, the measured log-power correlation follows Li₂(\|G\|²)/Li₂(1). Real-signal cells near DC and Nyquist can have appreciable pseudo-covariance, as the one-bin-from-DC probe demonstrates | Gates checked. [B1 readings](../../latent-recovery-pilot-brief.md) |
| Real signals within about two bins of DC or Nyquist need the pseudo-Gram entry as well | Gates checked |
| Log power against log window length has slope +1, −1, 0 for tone, impulse, noise, exactly | Gates diffed; numbers from the chip's report |
| With the synthesizer itself as the generator, a code of dimension k ≤ 5 is recovered from about 2k to 4k frequency bins at 30–40 dB | Gates checked |
| For a harmonic source the pitch posterior is multimodal with octave traps; tempering and cut SMC recover it in 10 of 10 replicates where the nested sampler's chains stay trapped | Chip's report |
| The built support gate admits no singular spline design and refuses 78% of the exactly full-rank ones that have more sites than columns, nearly all well conditioned. The ordered matching read rank exactly on 23 400 designs | Chip's report |
| A one-sided Green's-function basis fits a sum of decaying resonances with about two knots per resonance, twenty to eighty times fewer than cubic pieces at lower error | Chip's report |
| Noiseless two-layer Gaussian ReLU generator, widths 250 and 600, k = 10, m = 150. Descent with a sign check recovered 30/30 codes at relative error < 10⁻³. One reached the 50,000-update cap with relative error 6.84 × 10⁻¹⁶ | Gate checked, endpoint errors recomputed and replay receipts matched |
| At fixed m = 150 and widths 250/600, the signed recovery curve passed its preregistered comparison with all 15 published counts. At k = 80, 90, 100 it recovered 20, 9 and 1 out of 30, against 17, 7 and 4 in the paper | Recomputed the 450 signed outcomes and all 15 exact comparisons |
| Across matched sign-on and sign-off runs at k = 10, 80 and 100, only k = 80 showed discordant pairs (six on-only, two off-only). Across all 540 signed and unsigned endpoints, all 105 missed recoveries reached the 50,000-update cap at angles of 105°–171° from the negative truth direction | Recomputed all endpoint classes and the three paired tests |
| Holding k = 80, both layer widths and all 30 problems fixed, signed recovery reached 0/30, 4/30, 20/30 and 29/30 at m = 80, 120, 150 and 240. The primary paired comparison gained 29 recoveries with zero reversals (two-sided exact McNemar p = 3.73 × 10⁻⁹). Repetition 27 reached the cap at every m. | Recomputed all 120 relative errors and endpoint classes, all six paired comparisons and the 30 truth/start preservation checks |

## What was corrected

Each of these was claimed in the discussion and then turned by a proof or a reading.

- **"The engines are an exact teacher with no new code."** The inverse problem is sound and the generator is the bottleneck. Through the BASS emulator, recovery fails for every k ≥ 2, with errors 3 to 32 times the true generator's. Neither a functional nor a directional layout emulates a harmonic comb (held-out R² at most 0.34).
- **"A resonance turns over near its decay time, a chirp near 1/√c."** The scaling holds; the constant is the window's. Rectangular: 1.26τ and 1.43/√c. Hann: 2.38τ and 2.09/√c.
- **"Slopes classify content."** In expectation. A single frame's largest bin grows like ln w for noise, so a classifier reads averaged power.
- **"A trigonometric B-spline needs support under one period."** Half the carrier's period at order 2. The degeneracy belongs to the positive local element; the one-sided basis shows nothing at that spacing.
- **"A knot count estimates the number of events."** It does not: a kink costs 0.3 knots on average under a free count, and events read from the count carry a standard deviation of 6.8.
- **"Free knots pay for sparse events."** In theory. At the sizes the engine runs (n = 2000, σ = 10⁻³) the free-knot error reaches the noise floor by 8 to 16 knots, before the rate can show.
- **"The whitener is fixed by the windows."** For noise-like cells away from DC and Nyquist. A cell with a strong deterministic component has log-power variance about 2/SNR.
- **"Measurements scale linearly with k."** Consistent with the pilot, on thin evidence: k ≤ 5, four ratios, the same model generating and inverting, a Gaussian likelihood on log power.

<a id="the-compressed-sensing-guarantees-used-by-analogy-not-yet-tested"></a>

## The compressed-sensing guarantees

All six Lean statements are proved. They cover the angle map and Gaussian gated second moment, the secant bounds, the deterministic comparison field's two zeros, isotropic projection energy and the count of central sign regions. The G1 probes reproduced the decline in signed recovery across the published noiseless depth-two linear experiment, compared the sign check on matched problems and isolated the effect of measurement count at fixed dimension. The guarantees for learned decoders and audio measurements remain open. Statements below are from the papers' LaTeX (`1705.07576v3`, `1812.04176v1`, `1807.04261v1`, `2008.10579v1`, `1805.08855v2`, `2006.07953v2`). [Lean reports](generative-prior-statements-brief.md#round-1-report--2026-11-16) · [G1 linear results](generative-prior-probes-brief.md#g1-report--linear-curve-and-paired-sign-comparison) · [Fixed-dimension results](generative-prior-probes-brief.md#g1-report--fixed-k-measurement-count)

**What they say.**
- *The setting.* A generator G : ℝᵏ → ℝⁿ, d ReLU layers with no bias, and a signal y₀ = G(x₀). Observe Ay₀ with m ≪ n rows and minimize f(x) = ½‖AG(x) − Ay₀‖².
- *The landscape theorem* (Hand and Voroninski). If every layer satisfies the weight distribution condition and A satisfies the range restricted isometry condition, f has a strict descent direction at every nonzero x outside two small balls, one around x₀ and one around −ρ_d·x₀, and zero is a local maximizer. With A = I this says such a network is approximately invertible.
- *The weight distribution condition.* For every pair x, y, the sum Σᵢ 1{wᵢ·x > 0}·1{wᵢ·y > 0}·wᵢwᵢᵀ over a layer's rows is within ε of Q_{x,y} = ((π − θ)/2π)·I + (sin θ/2π)·M, with θ the angle between x and y and M the map exchanging their directions. Q is that sum's expectation for Gaussian rows, so a layer acts on a pair through its angle alone.
- *The range restricted isometry condition.* A preserves inner products between differences of points in G's range, to within ε.
- *When they hold.* Gaussian layers that expand by nᵢ ≥ c·nᵢ₋₁·log nᵢ₋₁, and a Gaussian A with m > c·d·k·log(n₁⋯n_d), satisfy both with high probability. The count is linear in k.
- *The algorithm.* Subgradient descent with a sign check converges to x₀ up to the noise level.
- *Around it.* The same from magnitudes, where a sparsity prior needs k² measurements; noise reduced by about k/n, by an encoder in one pass too; and the lifted form, ‖G(x)G(x)ᵀ − M‖²_F with M a sample covariance.

**Why they matter here.**
- The measurement budget follows the code's dimension, with an efficient algorithm and not only an existence claim.
- Under the papers' hypotheses, descent with the prescribed sign check can supply a recovered code as a teacher target. Transfer to a learned speech decoder and structured measurements is open.
- The two conditions suggest design constraints on the decoder weights and on pairs of secants of its range. A local pullback metric measures infinitesimal changes. The range-restricted isometry condition also compares differences between distant signals.

**What the work has used.**
- The deterministic half of the amplitude-flow convergence, as the contraction result in Lean.
- The mechanism of the three-pattern uniqueness theorem, as phase propagation in Lean.
- The Bayesian analogue on a parametric synthesizer: a code recovered from 2k to 4k bins, and octave traps standing in for the impostor at −ρ_d·x₀.

**What is untested.**
- Depth three, magnitude-only observations and the calibration stack on a ReLU generator remain open. The fixed-k experiment isolates measurement count in one depth-two noiseless setting. It measures the finite solver's recovery frequency. Whether another code exactly reproduces the measurements on an individual failed trial remains undetermined.
- The weight condition has not been measured for anything of ours, such as a layer of directional hinges under its prior, which has sparse directions and thresholds.
- The isometry condition has not been measured for spectrogram rows or magnitudes on any signal manifold, whitened or not.
- The pilot's held-back log-power reconstruction error can be read against k/m. The rank-k projection identity concerns expected squared error in the coordinates carrying white noise. Here the pilot added white noise to the waveform before taking log power.

**Remaining work in the two briefs.**
- S17–S22 complete the required Lean brief. S18 derives the angular integral for Q directly under the product standard Gaussian law. S22 proves the bound on full-dimensional sign regions with no general-position assumption. [Statements and reports](generative-prior-statements-brief.md)
- The numerical brief orders the rest as G1 depth three and magnitude observations, G3 weight-condition errors across layer widths and direction priors, G4 distortion of secants under structured audio measurements, G2 calibration samplers against descent, and G5 analysis of noise in the existing pilot receipts. [Probe definitions and order](generative-prior-probes-brief.md#the-probes)
- *What they would decide:* whether a gradient-informed move on the code is justified, which the benign landscape suggests; how wide a decoder must be; whether our views are isometric enough; and how the engine's inversion scales in k against descent.

### Analytical controls for the generative-prior probes

- **Stationary points (G1).** From the two papers' formulas, the predicted negative multiple at depth two has radius 1/π for linear observations and approximately 0.669152 for magnitudes. In the phaseless formula put θ = arccos(1/π), giving 2 sin(θ)/π + (1 − 2θ/π)/π. The frequency with which random starts reach either neighborhood is a separate empirical question. [Linear landscape](https://arxiv.org/abs/1705.07576) · [Phaseless landscape](https://arxiv.org/abs/2008.10579)
- **Why the negative ray matters (S17, S20, G1).** Write θ₀ = θ, θᵢ₊₁ = g(θᵢ), and cᵢ = cos θᵢ. Then cᵢ₊₁ = (1 − θᵢ/π)cᵢ + sin(θᵢ)/π. At θ₀ = π the first multiplier vanishes, so expanding this recurrence gives exactly the paper's finite sum for ρ_d = cos(gᵈ(π)) for d ≥ 1. Thus convergence of the S17 iterates to zero also gives ρ_d → 1. In the deterministic linear comparison objective, two codes of norms r,t and angle θ have loss (r² + t² − 2rt cos(gᵈ(θ)))/2ᵈ⁺¹. Along the opposite ray its minimum occurs at r = tρ_d. For k ≥ 2 and d ≥ 2 one can rotate away from that ray at fixed radius. This strictly lowers the loss because g is strictly increasing, whereas changing the radius in either direction raises it. The predicted point is therefore a saddle of this comparison objective. The argument does not establish a capture frequency for finite networks. [Comparison objective and ρ_d, §5.1](https://arxiv.org/abs/1812.04176)
- **Finite-network control (G1).** Each bias-free generator is positively homogeneous. With y = AG(x₀) and v = AG(−x₀), the noiseless linear loss at −a x₀ is ½‖av − y‖² for a ≥ 0. When v ≠ 0 its exact minimizing multiplier is max(0, ⟨v,y⟩/‖v‖²). If v = 0 the entire ray has the same loss. This per-network control distinguishes displacement along the ray from a terminal code pointing away from it.
- **A sparse-direction witness (G3).** Take the threshold-free independent direction prior on k ≥ 3 coordinates with active count a uniform on {1, 2, 3}. Scale n independent unit rows by √(k/n), so the expected full Gram is I. For x = y = e₁ the WDC target is I/2. Coordinate 2 contributes to the gated diagonal only if both coordinates are active and the first weight is positive. Conditional on a, that diagonal has expectation k × a(a−1)/(k(k−1)) × 1/(2a) = (a−1)/(2(k−1)). Averaging over a gives 1/(2(k−1)). The population operator-norm gap is therefore at least (k−2)/(2(k−1)), already 3/7 at k = 8. This gives a direction and a quantitative prediction before the width sweep. Thresholded gates depend on the input radius as well. [Engine direction law](../../../../src/maths/regression/spline/bass/DirectionalTermSpace.cs)
- **Measurement units (G4).** Count real measurement rows and report their rank. Test the inner product of two normalized secants, including the same-secant case. A finite search gives a lower bound on the uniform worst case. For real linear rows B, whitening with the positive square root of the pseudoinverse of BBᵀ gives a squared norm equal to the energy projected onto the row space of B. This supplies a check even with redundant rows. The pilot's kept log powers are nonlinear observables and require their own measurement-space interpretation. [Pilot transform](../../receipts/latent-recovery/probe/Synthesizer.cs)

## What is missing, in build order

Stated for any application. The order follows the pilot's main finding: generator adequacy first, then honest noise, then robust inference, then the time–frequency pieces.

**A. Generator adequacy, the present bottleneck.**

| Piece | Why | Specification or evidence | State |
|---|---|---|---|
| Declared derived inputs: a coordinate computed from existing inputs by a declared map, such as ln f − ln f₀ | Where structure lives on a known invariant, the basis should knot along it | The comb is a near-periodic function of f/f₀; dilation covariance says the same | Proposed |
| An L-spline term family: one-sided Green's-function columns of a declared operator, its parameters sampled or profiled | Two knots per resonance in place of tens; the atom model's exact home | S8; the spacing sweep and fits of the spline-space probes | Runs in probe code under BARS's `IBasis`. The pre-stated test did not decide it; the user's call |
| A support rule on supports, with the ordered matching in front of the factorization | The mass gate refuses most full-rank designs and had to be turned off for signed columns | S7; the admission probe | Proposed. Whether the mass gate gives way is a modelling question for the user |
| Births that can climb a flat likelihood | Cubic free knots stayed near zero on 13 of 15 oscillatory signals | The oscillator fits | Proposed; the built tempering ladder is the first thing to try |

**B. Honest noise.**

| Piece | Why | Specification or evidence | State |
|---|---|---|---|
| The output inner product Q in BASS's core | Residuals correlated along the output index are counted as independent today | S4, the Khatri–Rao Gram, the level under Q; the Li₂ table needs no data | Five implementation slices are proposed. Both prior Grams and all three precision-weight patterns are in scope, including a different location profile for each case. [Q extension decisions](../../bass/output-inner-product-brief.md#decisions) |
| An emulator's error covariance across outputs, from held-out residual vectors | A functional emulator's residuals read lag-1 r of 0.70–0.83 while its draws carry σ²·I | The pilot's Stage 1 | Proposed |
| A noise law for log power | Its noise is skewed with known shape; the Gaussian core is biased on it | The dilogarithm results; the Whittle model exists in BARS only | Not posed |

**C. Robust inference.**

| Piece | Why | State |
|---|---|---|
| Tempering or cut SMC as the default where the posterior is multimodal | The nested sampler's chains trap on octaves | Built; a default to declare |
| A move for narrow, curved conditionals | The random walk fails the convergence verdict at k ≥ 3 with a stable answer | Proposed |
| A discrepancy kept out of the parameters' span | The declared one was confounded with gain and skirt: up to 3.3 times the error | Proposed |

**D. Time–frequency proper.**

| Piece | Why | State |
|---|---|---|
| A blur-aware design: a term's factor averaged under each cell's window kernel | Several window lengths become one latent surface seen through known blurs | Argued, not probed. The moment code already averages a factor under a measure |
| A grid-aware product design | On gridded data the Gram factors by axis | S5 proved; BAMRS not built |
| An SMC path that takes arriving data, on a sliding window | Coefficients are real-time given structure; structure needs a sequential path | Argued, not probed |
| Local refinement in time by frequency | A tensor product refines time for the whole band at once | The multiresolution candidate the BAMRS brief sets aside |

**E. Beyond the engine,** once A to C serve.
- A neural generator behind the emulator contract, replacing the BASS emulator without changing the inversion.
- An amortized encoder trained on the engine's posteriors and on the generator's samples.
- Event tokens with continuous attributes in place of a fixed frame grid.
- Conformal guarantees on consistency losses, and risk control on turn-taking.

**TorchSharp as a candidate for the learned generator.** The first concrete use would be training a decoder on batches of signals, then freezing its weights and differentiating the measurement loss with respect to the latent code during recovery. A decoder-specific .NET assembly could own the tensors, optimiser and checkpoints. A fixed checkpoint can supply one draw to `ICalibrationEmulator`, with its prediction in physical units and an explicit emulator-error model. The existing contract supports concurrent single-code predictions. Its `PredictLevels` operation enumerates categorical levels for one continuous code. Training batches and latent gradients belong with the decoder. [Calibration contract](../../../../src/maths/regression/spline/bass/calibration/CalibrationEmulator.cs) · [TorchSharp training example](https://github.com/dotnet/TorchSharp#torchsharp)

The new `deep-signals` project owns the initial audio path inside ThermoMapper. NAudio.Core decodes bounded WAV windows into channel-first TorchSharp tensors. Its CPU tests cover the recorder sample's 16 kHz mono PCM format, stereo IEEE-float audio and a known gradient update. The example recorder, player and waveform display in the workstream's `csharp/` folder remain references for later application work. [Audio project](../../../../src/deep-signals/README.md)

For a first comparison, load identical G1 weights, measurements and initial codes into TorchSharp. Check forward outputs, signed losses and latent gradients in float64 against the C# reference before timing. Measure end-to-end recovery for one code and for the intended batch size, including transfers and peak native memory. Lower precision needs its own stopping tolerance and recovery check. Production adoption would depend on those measurements and a deployment test with the pinned CPU or CUDA LibTorch runtime on the intended host. [G1 probe](../../receipts/generative-prior/probe/Network.cs) · [Native runtime dependencies](https://www.nuget.org/packages/TorchSharp/0.107.0) · [Tensor lifetime management](https://github.com/dotnet/TorchSharp/blob/main/docfx/articles/memory.md)

## Inspecting regression fits before training a network

The [phaseless audio statements sequence](phaseless-audio/README.md) asks what the recorded window powers determine about a waveform, how errors propagate across quiet intervals, and when a population of such measurements determines a signal model.

The current implementation work is the Lean and numerical briefs. The user wants to inspect audio fits, BASS uncertainty and the generative-prior results in that order. The first signal-processing literature pass examines windowed measurements, speech prediction residuals and sparse innovations. [Audio analysis brief](audio-analysis-brief.md)

| Proposed view | What goes into it | What to check on the first example |
|---|---|---|
| 1. Waveform fit and residual | One WAV channel, the regression prediction at the same sample indices, and observed minus fitted amplitude. Use the recording's sample rate and starting frame to put both panels on one time axis. Identify the fit, channel and fitted or held-out interval. | Use a known tone with a second component deliberately omitted from the fit. The residual should contain the omitted component at the original amplitude. Then inspect a short recording, zoom across a transient, and compare its timing in both panels. |
| 2. BASS predictions and uncertainty | Predictions from the existing regression engine on an identified output axis. Use posterior draws and their weights to show the fitted mean and a pointwise interval with its level and target. A predictive interval for a new observation includes the observation-noise model. | Replay the plotted quantiles from the draws and weights. An audio pilot that predicts log power gets a frequency axis and log-power units. A predicted waveform requires a fit whose output is waveform amplitude. |
| 3. Generative-prior recovery | Read the probe receipts by matched trial ID. Show relative latent error, the recovery threshold, and aggregate counts at each measurement setting. Preserve the sign setting and update-cap flag. | Reproduce the signed recovery counts and identify any trial that recovers at one setting and fails at a larger one. |

For a first inspection, a workstream probe can write a reproducible figure and its input readings under `artifacts/`, with the compact replay source and receipt alongside the existing regression receipts. Production WAV handling and audio measurements belong in `src/deep-signals/` (`projects/deep-signals/`, `tests/deep-signals.Tests/`). Existing regression fits and posterior prediction stay with their owners in `src/maths/regression/`.

For permanent charts, the Viz v2 plan assigns producer-specific translation to `Viz.Adapters.*` and drawing to the workbench. A future audio adapter can map sample indices, physical units and fit provenance into `Viz.Contracts`. A chart panel can then share time selections across the observed, fitted and residual series. The current v2 implementation has neutral study and scene contracts. Chart panels are scheduled for phase 4, after the EyeTorus slice. This placement is a proposal for review before creating adapter or renderer projects. [Viz assembly direction](../../../viz/program.md#assembly-direction) · [Migration roadmap](../../../viz/migration-roadmap.md)

## Open questions for the user

- Does the oscillatory term family earn its own brief?
- Does the mass gate give way to "one site per support, in order", or stay as a declared prior beside it?
- For the Q extension, decisions 3–7 and the default prior Gram for a newly declared nonidentity Q remain open. [Decision brief](../../bass/output-inner-product-brief.md#decisions)
- Promotion of the noise-law modules; the first nine are promoted.
- Which reading leads, fitting the spectrogram's surface or fitting the waveform with atoms? For a linear transform with full column rank, Gram whitening preserves the waveform's squared residual norm. A fit to magnitudes or log powers still has its own observation model.
- Two source notes from the numerical chip: a suspected misprint at the opening of Schumaker § 4.8, unverified, and a conversion defect at Lemma 10.16.

## Records

- **Briefs and their reports,** at their paths on 2026-11-16 (the user is moving them under this folder): `issues/lean/time-frequency-statements-brief.md`, `issues/lean/spline-space-statements-brief.md`, `issues/lean/noise-law-statements-brief.md`, `issues/regression/latent-recovery-pilot-brief.md`, `issues/regression/spline-space-probes-brief.md`. Each report ends in exposition written to be spliced into docs.
- **Reading list:** `reading-list.md` in this workstream's folder.
- **Sources:** Voroninski's papers as LaTeX under `codex-scientiae/supellex-private/gauntlet/voroninski`; Schumaker as markdown under `graveyard/pdf2md-prototype/schumaker`, with `errata.md` for the book's misprints.
- **Probes and receipts:** `artifacts/latent-recovery/`, `artifacts/spline-space/`, and the matching paths under `issues/regression/receipts/`.
