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
- *Consistency.* An overcomplete representation must reproduce itself under its atoms' Gram matrix; views at different window lengths are related by a known blur; speaker and channel are constant across turns.
- *A teacher.* A solver, or the engine's posterior, wherever the guarantees hold.
- *Synthesis.* The prior's own samples carry exact labels for an encoder.

The recovery theorems do not perform these. They certify the tasks: they say when a lossy view determines the latent, which is when the encoder's target is well defined.

**Conformal-Bayes.**
- A held-back view is its own label, so calibration data is free and every consistency loss can carry a finite-sample coverage statement.
- Coverage is on observables. A latent gets it only on synthetic data, and weighted conformal carries that to real audio.
- The exchangeable unit is the conversation or the speaker.
- Turn-taking is a risk-control problem: start speaking before the other has finished at most α of the time.
- Test martingales give change detection with false-alarm control.

## What is established

How far each row was checked is in the last column. "Elaborated" means the coordinating session re-elaborated the statements and printed their axioms; "gates checked" means it diffed the pre-stated gates against their final text and matched headline numbers to receipts; "chip's report" means neither was done.

**Proved in Lean** (no `sorry`, standard axioms). The first nine modules were promoted to `lean/Lemmas/` on 2026-11-16 (`f0cc4e2`); the noise-law modules are in `lean/Enthymemata/TimeFrequency/`. Paths below are relative to the stage.

| Result | Module | Checked |
|---|---|---|
| Edge products and magnitudes pin a signal up to one phase exactly when the graph is connected; polarization supplies the edge products from magnitudes | `TimeFrequency/PhasePropagation.lean` | Elaborated |
| A unit-energy window of length w reads \|A\|²·w from a tone, \|a\|²/w from an impulse, σ² from white noise | `TimeFrequency/WindowScaling.lean` | Elaborated |
| For symmetric 2×2 H, ‖SᵀHS‖²_F ≥ 2(det S)²\|det H\|, attained at every cell area | `TimeFrequency/FlatnessBound.lean` | Elaborated |
| A regularity inequality on a ball gives linear convergence to the noise level | `TimeFrequency/Contraction.lean` | Elaborated |
| The marginal under noise σ²(I ⊗ R), with the prior on the whitened Gram, is the unit-noise marginal on whitened data plus a determinant term; the Khatri–Rao Gram under Q; the level under Q | `BassWhitenedMarginal.lean` | Elaborated |
| A Kronecker design's Gram and projection factor by axis | `BamrsKroneckerGram.lean` | Elaborated |
| A design's full column rank is equivalent to a greedy ordered matching of sites to supports, given Theorem 4.61 as a premise (met at order 1) | `Spline/DesignRank.lean` | Chip's report |
| The order-2 trigonometric B-spline is positive while its support is under half the carrier's period, and the condition is sharp; the damped oscillator's Green's function gives a one-sided basis | `Spline/TrigonometricBSpline.lean` | Chip's report |
| Free knots approximate any monotone f by steps within (f(b) − f(a))/(2(k + 1)); uniform knots leave some f at ½ for every k | `Spline/FreeKnotSteps.lean` | Chip's report |
| The covariance of log powers is Σκⁿ/n², from the Laguerre coefficients −1/n of the logarithm, with the bivariate expansion as a named premise | `TimeFrequency/LogPowerCovariance.lean` | Chip's report |
| Cov(\|a\|², \|b\|²) = \|E[a·conj b]\|² + \|E[a·b]\|², and Var(xᵀHx) = 2·tr((HΣ)²), each with the Gaussian pairing rule as a named premise | `TimeFrequency/LogPowerCovariance.lean`, `TimeFrequency/QuadraticFormVariance.lean` | Chip's report |
| ∫₀^∞ (1 − (1 − e^{−x})ⁿ) dx = H_n; the rectangular window's turnover is the root of eᵘ = 1 + 2u, between 1.25 and 1.26 | `TimeFrequency/PeakOfExponentials.lean`, `TimeFrequency/WindowTurnover.lean` | Chip's report |
| The response's second moments under the g-prior on the whitened Gram are the covariance the whitened marginal reads | `BassWhitenedMarginal.lean` | Chip's report |

**Measured,** under gates written before their readings (`issues/regression/receipts/latent-recovery/` and `receipts/spline-space/`).

| Finding | Checked |
|---|---|
| The covariance of log powers of circular Gaussians with squared coherence κ is Li₂(κ). The correlation of log power between any two cells of white noise is Li₂(\|G\|²)/Li₂(1) from the windows' Gram entry alone | Gates checked |
| Real signals within about two bins of DC or Nyquist need the pseudo-Gram entry as well | Gates checked |
| Log power against log window length has slope +1, −1, 0 for tone, impulse, noise, exactly | Gates diffed; numbers from the chip's report |
| With the synthesizer itself as the generator, a code of dimension k ≤ 5 is recovered from about 2k to 4k frequency bins at 30–40 dB | Gates checked |
| For a harmonic source the pitch posterior is multimodal with octave traps; tempering and cut SMC recover it in 10 of 10 replicates where the nested sampler's chains stay trapped | Chip's report |
| The built support gate admits no singular spline design and refuses 78% of the exactly full-rank ones that have more sites than columns, nearly all well conditioned. The ordered matching read rank exactly on 23 400 designs | Chip's report |
| A one-sided Green's-function basis fits a sum of decaying resonances with about two knots per resonance, twenty to eighty times fewer than cubic pieces at lower error | Chip's report |

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

## The compressed-sensing guarantees: used by analogy, not yet tested

These theorems are the certificate the whole approach leans on, and so far the work has touched them only at the edges. Statements are from the papers' LaTeX (`1705.07576v3`, `1812.04176v1`, `1807.04261v1`, `2008.10579v1`, `1805.08855v2`, `2006.07953v2`).

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
- Local methods suffice, which is what makes a solver a trustworthy teacher.
- The two conditions are design constraints: a decoder that expands with well-spread weights, and views that are near-isometric on the signal manifold. The second is the global form of the whitened metric pulled back through the decoder.

**What the work has used.**
- The deterministic half of the amplitude-flow convergence, as the contraction result in Lean.
- The mechanism of the three-pattern uniqueness theorem, as phase propagation in Lean.
- The Bayesian analogue on a parametric synthesizer: a code recovered from 2k to 4k bins, and octave traps standing in for the impostor at −ρ_d·x₀.

**What is untested.**
- No ReLU generator has been run. The papers' regime has not been reproduced, and the engine's inversion has not been tried on their object.
- The weight condition has not been measured for anything of ours, such as a layer of directional hinges under its prior, which has sparse directions and thresholds.
- The isometry condition has not been measured for spectrogram rows or magnitudes on any signal manifold, whitened or not.
- The k/n noise reduction has not been read from the pilot's receipts.

**Posable next.**
- *In Lean:* Q_{x,y} as the Gaussian expectation; the angle map it induces, cos θ′ = ((π − θ)·cos θ + sin θ)/π, with θ′ < θ on (0, π]; that the isometry condition bounds the compressed objective's inner products by the uncompressed ones; the one-layer landscape in expectation.
- *Numerically:* success against m/k for the papers' algorithm on random expansive networks; the same inverse problem through the calibration stack with G behind the emulator contract, reading threshold, coverage and the impostor; the weight condition's residual against width, for Gaussian layers and for directional-hinge layers; the isometry constant for structured views on the synthesizer's range.
- *What they would decide:* whether a gradient-informed move on the code is justified, which the benign landscape suggests; how wide a decoder must be; whether our views are isometric enough; and how the engine's inversion scales in k against descent.

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
| The output inner product Q in BASS's core | Residuals correlated along the output index are counted as independent today | S4, the Khatri–Rao Gram, the level under Q; the Li₂ table needs no data | Prepared, not built: `issues/regression/bass/output-inner-product-brief.md`, seven decisions and five slices, awaiting the user's approval |
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

## Open questions for the user

- Does the oscillatory term family earn its own brief?
- Does the mass gate give way to "one site per support, in order", or stay as a declared prior beside it?
- The Q extension's seven decisions, each at a recommendation in its preparation brief.
- Promotion of the noise-law modules; the first nine are promoted.
- Theorem 4.61 at order 2, left unattempted, which would make the rank result unconditional for degree-1 splines.
- Which reading leads: fitting the spectrogram's surface, or fitting the waveform with atoms. Whitening the first collapses it into the second.
- Two source notes from the numerical chip: a suspected misprint at the opening of Schumaker § 4.8, unverified, and a conversion defect at Lemma 10.16.

## Records

- **Briefs and their reports,** at their paths on 2026-11-16 (the user is moving them under this folder): `issues/lean/time-frequency-statements-brief.md`, `issues/lean/spline-space-statements-brief.md`, `issues/lean/noise-law-statements-brief.md`, `issues/regression/latent-recovery-pilot-brief.md`, `issues/regression/spline-space-probes-brief.md`. Each report ends in exposition written to be spliced into docs.
- **Reading list:** `reading-list.md` in this workstream's folder.
- **Sources:** Voroninski's papers as LaTeX under `codex-scientiae/supellex-private/gauntlet/voroninski`; Schumaker as markdown under `graveyard/pdf2md-prototype/schumaker`, with `errata.md` for the book's misprints.
- **Probes and receipts:** `artifacts/latent-recovery/`, `artifacts/spline-space/`, and the matching paths under `issues/regression/receipts/`.
