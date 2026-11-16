# Latent recovery from lossy magnitude views: a numerical pilot

Date: 2026-11-16. Status: **commissioned by the user on 2026-11-16** ("one opus chip to do the numerical research"). Drafted by the coordinating session. Its sibling is the Lean chip's [statements brief](../lean/time-frequency-statements-brief.md).

## Why

A discussion with the user on 2026-11-16 explored unsupervised training of an audio model with the regression engines as teacher and measuring stick. The picture rests on one claim: a lossy magnitude view of a signal determines a low-dimensional latent code, with measurements in proportion to the code's dimension k. Voroninski and co-authors prove versions of it for random expansive networks. The engine's calibration stack is the Bayesian form of the same problem at small scale: an emulator as the generator, a posterior over the code, a discrepancy for misfit. This pilot runs that form on built code, before any deep model exists, to learn where it holds and which primitive is missing first.

This is exploration on probes. The engine stays as it is: a gap found here is a primitive to propose in the report, stated for any application, and built later under its own brief.

## Sources

- **The papers,** LaTeX under `D:\aghado01\codex-scientiae\supellex-private\gauntlet\voroninski` (each folder a versioned arXiv id; `inventory.jsonl` names the entry `.tex`):
  - `1812.04176v1`: recovery from compressive measurements with sample count linear in k; its introduction's one-parameter example is the pilot's starting case.
  - `1807.04261v1` and `2008.10579v1`: the same from magnitudes, where a sparsity prior needs k² measurements.
  - `1805.08855v2`: noise energy reduced by about k/n, by optimization over the range and by an encoder in one pass.
  - `1705.07576v3`: the landscape's one impostor, a negative multiple of the truth.
- **The engine,** under `src/maths/regression/spline/bass/`: `BassData.cs` (functional and scalar responses), `calibration/` (the emulator contract, `BassEmulators`, the conditional sampler, cut SMC, the measurement and covariance models), `BassResidualDependence.cs`, the tempering and SMC runs; `src/maths/regression/subspaces/RankPosterior.cs`; `src/maths/linalg/Fft.cs`; `src/synthetic/`.
- **Receipts:** [receipts/README.md](receipts/README.md) for what is tracked and what stays under `artifacts/`.

## Part A: the pilot

1. **A synthesizer outside the estimator.** A code of dimension k maps to a waveform and then to a magnitude spectrum on a fixed frame and window, with additive noise at a stated level. Ground truth is the code that made each signal. Record the frame, the window, the sampling rate, the noise and every range.
2. **An emulator of it.** A BASS fit with the code as inputs and log magnitude over frequency as the response: the generator, with a posterior. Report its held-out accuracy and the residual-dependence check across frequency.
3. **The inversion.** The calibration sampler recovers the code from one observed spectrum of which only m bins are kept. Vary m/k from 1 upward, with the discrepancy off and on.
4. **The readings,** over R replicates:
   - code recovery error against m/k;
   - coverage of the code's 90% intervals;
   - modality: chains' agreement, and whether tempering or cut SMC changes the answer;
   - reconstruction error at the bins held back, which is the consistency-across-views signal;
   - reconstruction error against the noise level and k/m.

**Stages,** each with its own gate:
- *Stage 1, resonances.* The code sets the centre frequency, bandwidth and gain of one or two decaying resonances, k from 1 to about 5. The spectrum is smooth, so the emulator should serve and the posterior should have one mode.
- *Stage 2, a harmonic source.* The code includes a fundamental f₀ driving a harmonic comb through the resonances. Two difficulties are expected and are part of the result: the comb lies along diagonals in (log f₀, log f), which axis-aligned terms follow poorly, and the octave is a structured impostor. Try the layouts the engine offers (a functional response; a scalar response with frequency as a fixed input and a directional basis) and report which serves.

**What each outcome points to:**

| Outcome | Next primitive |
|---|---|
| Recovery holds near a small multiple of k, with honest coverage | Correlated noise along the output index (the Q the working-response contract declares) |
| Multimodal posterior | Moves aimed at the known impostor, or tempering as the default for this use |
| Emulator misfit dominates | The discrepancy term, or the basis |
| Cost stops it near k = 5 | Reduction first: activity scores, rotate-and-refit |

**Gates come first.** Before the first reading of each stage, write its gate into the report: sizes, seeds, tolerances and the rerun rule for a marginal failure. A gate changed after a failure is recorded as a change.

## Part B: two sweeps

- **B1, the covariance of log powers.** For jointly circular complex Gaussian (a, b) with squared coherence κ, the coordinating session claims Cov(log|a|², log|b|²) = Li₂(κ), the dilogarithm. It was stated from memory, and only the endpoints were checked (0 at κ = 0, π²/6 at κ = 1). Test it by seeded Monte Carlo over a grid of κ with stated error. Then test the form it implies for spectrograms of white noise: the correlation of log power between cells at offsets in time and frequency, and across window lengths, against the value predicted from the windows' Gram entries alone.
  - *Decision:* whether a fixed whitener can be tabulated from the windows, as the correlated-noise extension would need.
- **B2, window scaling.** With energy-normalized windows of length w, confirm that peak power scales as w for a tone, 1/w for an impulse and w⁰ for white noise; that a decaying resonance turns over near its decay time; and that a linear chirp of rate c turns over near w = 1/√c. Use the normalization the Lean chip's statement S2 fixes, and tell that chip if a reading contradicts a posed statement.
  - *Decision:* whether the window-length axis classifies content as claimed, before any engine fit is built on it.

## Scope

- **Probes live under `artifacts/latent-recovery/`:** code, data, logs. Small settings, probe sources and result summaries are tracked at the matching path under `issues/regression/receipts/`.
- **The engine is read and called as built.** `src/`, `tests/` and `docs/` stay as they are.
- **Names say what a thing does.** The papers' illustrative examples stay in prose.
- **Every run is seeded and reproducible** from the tracked settings on a pinned commit.

## Done

- Stage 1 complete against its gate. Stage 2 run and reported as it fell, including a layout that failed.
- B1 and B2 complete with their decisions answered.
- A report appended to this brief: the gates as pre-stated, the readings, what each outcome table row came to, the primitives proposed with what each would change, and what was left unread or unrun. A **Doc changes owed** heading. An **In progress** entry kept current while work is open.

## Working rules

- Commits land on `main` by pathspec, each one informative about its work, with no co-author trailer. The index is shared with other sessions, so nothing stays staged.
- Builds set `ArtifactsRoot` inline on each command. Heavy runs go at below-normal priority with capped threads; the machine is shared.
- Each step ends with a fresh-eyes review of its own output, fixed before it is reported.
- Messages to the coordinating session or the Lean chip are terse events; results travel in the report.

## Report (the numerical chip, from 2026-11-16)

### In progress

- 2026-11-16: probe project at `artifacts/latent-recovery/probe/` against a source snapshot of `be979dd` (`git archive` into `artifacts/latent-recovery/source/`), so other sessions' uncommitted edits never reach a reading. B1 and B2 gates written below; their runs next. Stage 1's gate follows once the synthesizer is fixed and before any emulator is fitted.

### Gates as pre-stated

Each gate was written here before the first reading of its stage or sweep. A gate changed after a failure is recorded below it as a change, with the reason.

**B1, the covariance of log powers.** Written before its first reading.

- *Pairs.* κ ∈ {0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1}; a = z₁, b = √κ·e^{iφ}·z₁ + √(1 − κ)·z₂ with z₁, z₂ independent circular standard complex normals (E|z|² = 1) and φ = 1 rad. 4×10⁶ pairs per κ in 200 batches of 2×10⁴; the estimate is the mean of the batches' unbiased sample covariances of (ln|a|², ln|b|²), its standard error the batches' sd/√200. Seeds: `SeedTree.Derive(1116001, 34)`, child i for κᵢ, child 15 + i its rerun seed, children 30 and 31 the complex and real noise of the cells below, 32 and 33 their rerun seeds.
- *Pass:* |z| ≤ 4 at every κ, z = (estimate − Li₂(κ))/SE. *Marginal:* exactly one κ with 4 < |z| ≤ 5 and none above: rerun that κ on its rerun seed, pass if the rerun reads |z| ≤ 4. *Fail:* any |z| > 5, two or more marginal κ, or a failed rerun. Recorded beside each row: z against the linear alternative κ·π²/6.
- *Spectrogram cells.* Periodic Hann windows of length w ∈ {64, 128, 256} at unit energy (Σg² = 1); a 1024-sample signal; the reference cell centred at sample 512 at f = 0.25 cycles/sample; comparison cells at Δt ∈ {0, w/8, w/4, 3w/8, w/2, 3w/4, w} samples × Δf ∈ {0, 0.5, 1, 2}/w (the pair Δt = Δf = 0 left out); cross-length pairs (w₁, w₂) ∈ {(64, 128), (128, 256), (64, 256)} with the reference at w₁, the other cell at w₂, Δt ∈ {0, w₁/2}, Δf ∈ {0, 1/w₂}. Unit-variance complex circular and real white noise, 2×10⁵ realizations per kind in 100 batches (batch seeds `SeedTree.Derive` of the kind's seed). The estimate is the pooled sample correlation of the two cells' ln|X|², its SE the batches' correlations' sd/10. The prediction is Li₂(κ)/Li₂(1) with κ = |Σₙ conj(a₁ₙ)·a₂ₙ|², the two unit-energy atoms' Gram entry. For real noise only, a near-DC family as well: reference at f = 1/w, Δf ∈ {−1, 0, 0.5, 1, 2}/w.
- *Pass:* |z| ≤ 4 for every complex-noise pair and for every real-noise pair in the mid-band and cross-length families. *Marginal,* per noise kind: at most two pairs in (4, 5] and none above: rerun that kind on its rerun seed, pass if those pairs read |z| ≤ 4. The real-noise near-DC family is a reading, not a pass condition: it is expected to depart where the pseudo-Gram entry Σ conj(a₁)·conj(a₂) is not small, and its departures are recorded against that entry.
- *Decision:* a fixed whitener can be tabulated from the windows if the complex family and the real mid-band and cross-length families pass. The near-DC reading says where the pseudo-Gram entry must join the Gram entry in the table.

**B2, window scaling.** Written before its first reading. S2's normalization throughout: Σ g² = 1, power |Σ g_j·x_j·e^{−2πi f j}|².

- *Windows and lengths.* Rectangular and periodic Hann; even w on the quarter-octave grid 2·round(2·2^{j/4}), j = 0…40 (4 to 4096, repeats dropped). Slopes are least-squares slopes of ln P on ln w.
- *Tone,* x_j = e^{2πi f₀ j}, f₀ = 0.1234567: the DFT at f₀ equals |Σg|² within 10⁻⁹ relative at every w; the peak of a 16-times zero-padded FFT does not exceed it by more than 10⁻⁹ relative; slope over the whole grid 1 ± 0.01.
- *Impulse,* amplitude 1, best placement in the window: P = max g²; slope −1 ± 0.01.
- *White noise,* σ² = 1, complex circular and real: 4000 frames per (window, w, kind); P is the mean over frames and over the bins of the FFT zero-padded to the next power of two. Pass: |P − 1| ≤ 4·SE at every w (SE from the frames' means), and slope 0 ± 0.01. Marginal: one w per (window, kind) in (4, 5] SE: rerun that w on its rerun seed (seeds `SeedTree.Derive(1116002, 4·lengths)`, the second half the reruns). Recorded as a reading, not a pass condition: the mean of each frame's largest bin, which grows with the number of bins and is not w⁰.
- *Decaying resonance,* x_n = e^{−n/τ}·e^{2πi f₀ n} for n ≥ 0, τ ∈ {64, 128, 256} samples (amplitude decay time). g·envelope is non-negative, so the peak over f sits at f₀ and equals (Σ_j g_j·e^{−(s+j)/τ})², maximized over the window's start s. Pass, for each window: the turnover w* (the grid's argmax) has w*/τ ∈ [0.5, 4] at each τ; w*/τ across the three τ spans at most one grid step (ratio ≤ 1.25); slope over w ∈ [4, τ/8] ≥ 0.85; slope over w ∈ [8τ, 4096] ≤ −0.9. Recorded beside: the rectangular window's closed form at onset alignment, (Σ_{j<w} e^{−j/τ})²/w, whose continuum argmax is w ≈ 1.2564·τ.
- *Linear chirp,* x_j = e^{iπ c j²}, c ∈ {10⁻⁴, 2.5×10⁻⁴, 6.25×10⁻⁴}, the window centred at j = 0, w ≤ 0.8/c so the instantaneous frequency stays in band; P is the 16-times padded FFT's peak. Pass, for each window: w*·√c ∈ [0.5, 2] at each c; w*·√c across the three c spans at most one grid step; slope over w ≤ 1/(4√c) ≥ 0.9; slope over w ≥ 8/√c ≤ −0.8.
- *Decision:* the window-length axis classifies content as claimed if tone, impulse and noise pass and both turnovers sit at their scales in proportion to them. These sweeps are deterministic but for the noise, so the rerun rule applies to the noise alone.
