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
