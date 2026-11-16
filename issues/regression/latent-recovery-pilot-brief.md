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

- 2026-11-16, stopped at a usage limit. Probe at `artifacts/latent-recovery/probe/` against a `be979dd` snapshot (`artifacts/latent-recovery/source/`); build `artifacts/latent-recovery/build/`. **Done:** B1 and B2 run and read (below); the Stage 1 and Stage 2 gates committed before their first readings. **Running when stopped:** a background loop running `s1-run` for k = 1…5 and then `s1-gate` into `artifacts/latent-recovery/results/stage1/` (logs in `logs/s1-k*.log`, `logs/s1-gate.log`). k = 1–3 had finished and k = 4 was mid-grid. If the loop died, rerun `s1-run 4`, `s1-run 5` and `s1-gate`; the seeds make the rows reproducible. **Not yet written:** Stage 1's readings. **Next:** run `s2-run 2`, `s2-run 3` and `s2-gate` (into `results/stage2/`, after Stage 1, since the machine is loaded); then the outcome table, the primitives, the unread and unrun items, **Doc changes owed**, and receipts for stage1/stage2 under `receipts/latent-recovery/results/`.
- *Unread, provisional notes from k = 1–4, not readings against the gate:*
  - **Fits:** held-out R² was 0.995, 0.988, 0.986 and 0.919 for k = 1–4 (RMSE 0.15, 0.25, 0.28 and 0.57 nepers); each fit took 3–12 s. Residuals are correlated along ln f, with lag-1 r 0.70–0.81 against −0.016 for independent noise.
  - **k = 1:** R1 and R2 pass at both low noise levels, and the coverage is conservative (30/30 in most cells). The posterior σ is dominated by the IG(2, 0.3) prior at small m.
  - **Declared discrepancy:** at ratio ≤ 4 it fails to converge in 50–73% of replicates, and its error is 3–7 times the error without it.
  - **Reconstruction:** at the held-back bins it is set by emulator error, 0.16 against 0.01 for the reference at the lowest noise.

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

**Stage 1, resonances.** Written before any emulator was fitted. The synthesizer's levels were read first (`s1-peek`, no estimator involved), to set the floor and the noise.

- *Synthesizer* (`ResonanceSynthesizer`): fs = 8 kHz; a 512-sample frame with the onset at sample 128; x(t) = Σ_r G_r·e^{−πB_r t}·sin(2πF_r t) for t ≥ 0; a unit-energy periodic Hann window; 64 analysis frequencies log-spaced from 150 Hz to 3.6 kHz; the feature y = ln(|X(f)|² + ε), ε = 10⁻⁴. Code coordinates in order ln F₁ ∈ [ln 300, ln 1200], ln B₁ ∈ [ln 60, ln 240], ln G₁ ∈ [ln 0.5, ln 2], ln F₂ ∈ [ln 1400, ln 3000], ln B₂ ∈ [ln 80, ln 320]; the first k are the code, the rest fixed (B₁ = 120 Hz, G₁ = 1, the second resonance absent below k = 4, B₂ = 160 Hz at k = 4, G₂ = 0.7 at k ≥ 4). Noise: real white Gaussian on the waveform, σ²ₙ ∈ {10⁻⁶, 10⁻⁵, 10⁻⁴} per sample, which is about 40, 30 and 20 dB SNR over the frame and 0.01, 0.1 and 1 times ε per bin. The peek: the clean feature runs from the floor (ln ε = −9.2, where the far tails sit) up to about −1.5…0.1; noisy minus clean has rms about 0.03, 0.09 and 0.32 nepers at the three levels. The likelihood is Gaussian in y, so the noise model is misspecified by construction (the log of a noisy power), most at the highest level.
- *Emulator,* per k: BASS on a direct functional response, the k code coordinates as inputs and ln f as the output variable at the 64 frequencies; a Latin hypercube of n = 80 + 80k runs; MaxInteraction 3, MaxBasis 80; 2 chains, burn-in 8000, 500 kept per chain thinned by 4. Held out: 100 codes uniform on the box, RMSE and R² of the posterior-mean prediction against the clean feature; the fit's lag-1 residual dependence along ln f.
- *Inversion:* the nested reference sampler (the cut posterior) with 8 outer emulator draws; 3 chains per inner run from prior draws; burn-in 1000; R's decorrelation step every 100 scans; batches of 250 up to 2000 kept per chain; R̂ target 1.05 and ESS ≥ 50 per chain. One variance group, σ² ~ IG(ν = 2, ŝ = 0.3). Discrepancy off everywhere; on (D = [1, standardized ln f], s_v = 1) at the middle noise level with the BASS emulator. Unconverged inner runs are kept in the readings and counted.
- *Reference inversion:* the same calibration with the synthesizer itself as the emulator (one draw, no emulator noise), at every discrepancy-off cell. It separates what the BASS emulator costs from what the m bins and the noise allow.
- *Replicates:* R = 30 per k. θ* uniform on the box; one noise realization per level; one random order of the 64 bins, whose first m = ratio·k are kept, ratio ∈ {1, 2, 4, 8}, the rest held back. Shared across ratios, emulators and modes, so their comparisons are paired.
- *Modality subset:* k ∈ {2, 5}, ratio ∈ {1, 4}, the middle noise level, the first 10 replicates, both emulators: the nested sampler on the ladder (1, 0.6, 0.35, 0.2, 0.1, 0.05), and cut SMC with 400 particles and 8 outer draws.
- *Seeds:* `SeedTree.Derive(SeedTree.Derive(2026111601, 5)[k − 1], 16)`: design, test, fit, truth, noise, subsets, calibration, modality; children 8–15 the same roles for a rerun's fresh replicates (the design, test set and fit are reused).
- *Readings per calibration:* posterior mean and 5% and 95% quantiles per coordinate in unit coordinates; the recovery error, the RMS over coordinates of |mean − θ*| in unit coordinates; 90% coverage; converged inner runs; the posterior σ; reconstruction RMSE at the held-back bins against the clean and against the noisy feature; RMSE at the kept bins against the noisy feature; wall time.
- *R1, recovery holds near a small multiple of k:* at the two lower noise levels, with the BASS emulator and the discrepancy off, for every k from 1 to 5 there is a ratio ≤ 4 at which the median recovery error is ≤ 0.05, and every larger ratio meets it too.
- *R2, honest coverage:* in each cell R1 counts, the pooled 90% coverage over the R·k intervals lies in the two-sided 99% binomial acceptance band for p = 0.9 and n = R·k. Intervals within one replicate are not independent; the band is the stated tolerance all the same.
- *R3, modality* (a reading for the outcome table): a cell is flagged when more than 20% of its replicates have an unconverged inner run, or, in the modality subset, when tempering or cut SMC moves a coordinate's posterior mean by more than 2 posterior sd (of the untempered run) in more than 2 of the 10 replicates. Ratio 1 is expected to be multimodal (one bin cannot place a peak on one side); the outcome row reads ratios ≥ 2.
- *R4, emulator misfit* (a reading): in a cell, the BASS median error exceeds twice the reference's plus 0.01, or the BASS coverage falls below the band while the reference's lies in it.
- *R5, cost* (a reading): the fit's and the median calibration's wall time per k. Cost stops it near k = 5 if the k = 5 fit exceeds 2 h or the k = 5 grid cannot finish in 8 h at 6 threads.
- *Marginal failure and rerun:* R1 with a median error in (0.05, 0.06], or R2 with coverage outside the 99% band but inside the 99.9% band: rerun that k's affected cells once on the rerun seeds; the criterion passes only if the rerun passes outright, and both readings are reported.
- *Outcome rows:* row 1 (recovery near a small multiple of k with honest coverage) if R1 and R2 pass; row 2 (multimodal) if R3 flags a cell at ratio ≥ 2; row 3 (emulator misfit dominates) if R4 fires in cells where the reference meets R1's threshold; row 4 (cost) by R5.

**Stage 2, a harmonic source.** Written before any Stage 2 emulator was fitted; the synthesizer's levels were read first (`s2-peek`). Run and reported as it falls.

- *Synthesizer* (`HarmonicSynthesizer`): fs = 8 kHz; a 1024-sample frame of a steady source x_n = Σ_h A_h·sin(2π h f₀ n/fs), harmonics up to 3.8 kHz, A_h = 1/√(1 + (2(h f₀ − F)/B)²), one resonance; a unit-energy periodic Hann window; 128 analysis frequencies linearly spaced from 100 Hz to 2 kHz (14.96 Hz apart, against the window's main-lobe half-width of 15.6 Hz, so each harmonic is sampled at about two points); y = ln(|X|² + 10⁻²). Code: ln f₀ ∈ [ln 100, ln 400], ln F ∈ [ln 400, ln 1600], ln B ∈ [ln 100, ln 400], k ∈ {2, 3} (B = 200 Hz at k = 2). The f₀ box spans two octaves, so for half of it an octave impostor lies inside the box. Noise σ²ₙ = 10⁻³, about 25–30 dB SNR; the peek reads noisy minus clean rms 0.11–0.16 nepers with mean +0.07…+0.11 (the floor's bias); the clean feature runs from the floor (−4.6) to about 5.
- *Layouts,* on one design: (a) the direct functional response, ln f the output variable; (b) a scalar response with ln f an input fixed per observation and a directional basis (`DirectionalBasis` at its defaults: K = 1 directional hinges, independent directions; MaxInteraction 3). A direction in (ln f₀, ln f) can follow the comb's diagonals; an axis-aligned hinge product cannot.
- *Training:* a Latin hypercube of n = 100 + 100k runs (300 and 400; N = 128·n rows); MaxBasis 120; 2 chains, burn-in 8000, 500 kept per chain thinned by 4. Held out: 100 uniform codes, RMSE and R², and the lag-1 residual correlation along the bins (held out for both layouts; the engine's training check for (a)).
- *Inversion:* Stage 1's sampler settings and prior, discrepancy off; the reference inversion through the synthesizer; R = 20 replicates; ratio ∈ {2, 4, 8, 16} (m up to 48 of 128); the modality subset at k = 2, ratios {4, 16}, 10 replicates, all three emulators, tempered and cut SMC as in Stage 1.
- *Seeds:* `SeedTree.Derive(SeedTree.Derive(2026111602, 5)[k − 1], 16)` in Stage 1's roles.
- *Readings:* Stage 1's per cell, and for f₀ the posterior mass within a semitone of f₀*, of 2f₀* and of f₀*/2, the replicates whose mass at an impostor exceeds 0.5, and f₀'s 90% coverage.
- *Serves,* stated in advance: a layout serves when its held-out R² ≥ 0.9 and, at some ratio ≤ 8, its median recovery error is ≤ 0.05 with the 90% coverage inside Stage 1's 99% band. Nothing else is pass or fail and there is no rerun; a layout that does not serve is a result.

### Readings

**B1: passes. The dilogarithm claim holds, and a whitener can be tabulated from the windows' Gram entries away from DC and Nyquist.** Run on `be979dd`'s snapshot; rows in `receipts/latent-recovery/results/sweeps/`.

- *Pairs:* across the 15 κ, |z| ≤ 1.56 against Li₂(κ) (SE 0.7–1.9 ×10⁻³). Against the linear alternative κ·π²/6 the same estimates sit 7 to 311 SE away at every interior κ. The endpoints read 0.0008 ± 0.0009 at κ = 0 and 1.6450 ± 0.0019 at κ = 1 (π²/6 = 1.6449).
- *White-noise cells* (2×10⁵ realizations per kind): the complex family's 93 pairs read max |z| = 2.84, χ²/df 0.97 within lengths and 1.58 across lengths; the largest |estimate − prediction| is 0.0053 in correlation. The real family's mid-band and cross-length pairs (93) read max |z| = 3.26, largest gap 0.0059. So the log-power correlation of any two cells of white noise, within one window length or across lengths and at any time-frequency offset, is Li₂(|G₁₂|²)/Li₂(1), G₁₂ the unit-energy atoms' Gram entry.
- *Near DC, real noise* (the reading, not gated): 26 of 102 pairs depart beyond 4 SE, by up to 0.076 in correlation (z up to 38), all of them predicted too low. The worst pairs have the largest pseudo-Gram entry |Σ conj(a₁)·conj(a₂)| (0.67, 0.62, 0.48 for the eight worst), with the reference's own pseudo entry 1/6 at one bin from DC. A real signal's coefficients are not circular there, and the table needs the pseudo entry beside the Gram entry. The same holds by symmetry one bin from Nyquist; it was not run.
- *Decision:* yes. A fixed whitener for the log-power noise of a stationary white source is a table of Li₂(|G|²) over the windows' Gram entries, and its correlation part needs no data. Two qualifications travel with it: real signals within about two bins of DC or Nyquist need the pseudo-Gram entry too, and the table is the noise of the log power of noise, not of a noisy signal (a cell with a strong deterministic component has a smaller log-power variance, 2/SNR in place of π²/6).

**B2: fails its gate as pre-stated, on two Hann tolerances. The scaling claims hold; the turnover constants are the window's.** Rows and a post-gate supplement in `receipts/latent-recovery/results/sweeps/`.

- *Tone, impulse, noise, at S2's normalization:* all pass. The DFT at f₀ equals |Σg|² to 6×10⁻¹⁶ relative, and the padded grid's peak never exceeds it; the slopes read 1.0000 (tone) and −1.0000 (impulse) for both windows, and w·P for the impulse at w = 4096 reads 1.0000 (rect) and 2.6667 = 8/3 (Hann). Noise power is σ² at every w (max |z| 2.76, slopes within ±0.0006). No reading contradicts S2, so no message went to the Lean chip.
- *A frame's peak bin is not w⁰ for noise.* The mean of each frame's largest bin grows with slope 0.19–0.22 over w = 4…4096, to 8.1–8.9 times σ² at w = 4096: the maximum of about w exponential variables, near ln w. A classifier on the window-length axis has to read an expected or averaged power, never a single frame's peak.
- *Resonance:* the turnover sits in proportion to τ to three digits across τ = 64, 128, 256. Rect: refined w* = 1.255·τ, the closed form's 1.2564·τ. Hann: 2.376·τ. Rect passes every condition. Hann fails one: beyond the turnover its slope approaches −1 slowly, −0.90, −0.85 and −0.77 over [8τ, 4096] at τ = 64, 128, 256 against the gate's ≤ −0.9. The supplement's local slopes for Hann at τ = 64 read −0.78 over w/τ ∈ [8, 16), −0.92 over [16, 32) and −0.98 over [32, 64): the limit is −1, reached later than the gate assumed. The gate's tolerance rested on a wrong estimate of the window's curvature correction, which is of order (πτ/w)² with a large constant.
- *Chirp:* the turnover sits in proportion to 1/√c across the three rates (grid spread 1.13 rect, 1.11 Hann; refined constants 1.42–1.43 rect and 2.09–2.10 Hann). Slopes read 1.000 below and −1.00 ± 0.03 above. Hann fails the band: w*·√c = 2.16, 2.02, 2.25 on the grid against [0.5, 2].
- *The gate stands as written,* its failure recorded here; nothing was rerun, since these readings are deterministic.
- *Decision:* the window-length axis classifies content as claimed, with one change to the claim. Slopes +1, −1 and 0 separate tones, impulses and noise exactly in expectation. The turnover of a decaying resonance and of a linear chirp scales with τ and with 1/√c, but the constant belongs to the window: about 1.26τ and 1.43/√c for the rectangular window, 2.38τ and 2.09/√c for Hann, in its nominal length. Read the "near its decay time" and "near w = 1/√c" claims as "at a window-specific multiple of". A classifier has to carry that table and read averaged power.
