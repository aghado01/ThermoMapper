Both chips stopped on the usage limit, which resets at 5am Pacific. Before stopping, the Lean chip committed proofs for all six first-round statements with a report, and the numerical chip committed its gates, a passing log-power covariance sweep, a window-scaling sweep that failed two Hann-window tolerances, and a hand-over note. I have not reviewed either against its diff yet.

## The picture, with conversation as the data

- **The atom model is speech production.** Voiced speech is glottal pulses driving vocal-tract resonances. That is sparse excitation through a slowly varying linear operator, which is the L-spline of Schumaker's Chapter 10 in physical form.
- **Conversation adds slow latents that are constant by construction.** One speaker has one vocal tract for the whole recording, and one room and microphone colour everything. Constancy across turns is free supervision.
- **Prosody is the slow layer** and turn-taking the slowest: gaps, overlaps and backchannels are sparse events between two streams.
- **Nothing is labeled, and nothing needs to be.** Linguistic content is what remains once speaker, channel and prosody are factored out.

## Where conformal-Bayes fits

- **Self-supervision makes calibration data free.** Conformal's usual bottleneck is a labeled calibration set. Here the held-back view is the label, so every consistency loss can carry a finite-sample coverage statement at no labeling cost.
- **Coverage is on observables, not latents.** It certifies the prediction of a held-back view. A latent code gets coverage only on synthetic data where it is known. Weighted conformal is then the tool for carrying that to real audio.
- **The exchangeable unit is the conversation or speaker.** Frames and turns inside one are dependent. Your group cells are the right shape for this.
- **Turn-taking is a risk-control problem.** "Start speaking before the other person has finished at most α of the time" is a distribution-free guarantee on an interaction decision.
- **Test martingales give change detection with false-alarm control.** A speaker change or a new data source is a break in exchangeability.
- **Set size is an assumption-free measuring stick.** A better representation gives smaller valid sets, which complements evidence in bits.

Your collection already holds the core for this (risk control, beyond exchangeability, covariate shift, functional bands, martingales). What it lacks is listed below.

## Schumaker and wavelets

You are not hallucinating, but it is not in the chapters. Chapters 1–13 never mention wavelets. The third edition's Supplement (PDF `.016`) has a "Wavelets" passage covering spline wavelets, tensor-product spline wavelets and L-spline wavelets, and the New References (PDF `.018`) list about twenty wavelet items. Both are unconverted, and both are worth converting: the Supplement is the bridge and the references are a reading map.

## Already held, worth faithful markdown

| Source | Why |
|---|---|
| Schumaker Supplement and New References | The spline-to-wavelet bridge, including non-uniform-knot and L-spline wavelets |
| Ghosal and van der Vaart, Appendix E | Splines and wavelets stated for Bayesian nonparametrics |
| "Multivariate Lévy Adaptive B-Spline Regression" (2108.11863) | Reversible jump over a continuous dictionary, the closest relative of the atom model |
| "Adaptive B-spline knot selection using multi-resolution basis set" | Local refinement with free knots |
| The three Bayesian context-tree papers | Exact Bayesian inference for discrete sequences, relevant once events become tokens |
| "Conformal Prediction Bands for Multivariate Functional Data" | Bands over curves and surfaces |

## Worth hunting down

These are from memory, so check the titles.

| Topic | Source |
|---|---|
| Sparse innovations, L-splines and wavelets together | Unser and Tafti, *An Introduction to Sparse Stochastic Processes* |
| Time–frequency foundations | ~~Mallat, *A Wavelet Tour of Signal Processing*~~; Gröchenig, *Foundations of Time-Frequency Analysis* |
| Splines to wavelets | Chui, *An Introduction to Wavelets*; Lyche and Schumaker, "L-spline wavelets" |
| Rate of innovation | Vetterli, Marziliano and Blu, "Sampling signals with finite rate of innovation" |
| Spectrogram phase retrieval | Jaganathan, Eldar and Hassibi on STFT phase retrieval; Grohs and Rathmair, "Stable Gabor phase retrieval and spectral clustering" |
| Learning from lossy views alone | Bora, Price and Dimakis, "AmbientGAN"; Tachella, Chen and Davies on sensing theorems for unsupervised learning |
| Denoising without clean targets | Stein (1981) on the unbiased risk estimate; Efron, "Tweedie's formula and selection bias" |
| Speech production | Fant, *Acoustic Theory of Speech Production*; Makhoul, "Linear prediction: a tutorial review" |
| Prosody as a driven system | Fujisaki and Hirose (1984); Prom-on, Xu and Thipakorn (2009) |
| Turn-taking | Sacks, Schegloff and Jefferson (1974); Levinson and Torreira (2015) |
| Event streams | Daley and Vere-Jones, *An Introduction to the Theory of Point Processes*; Hawkes (1971) |
| Conformal for this data | Lee, Barber and Willett on hierarchical data; Gibbs and Candès on adaptive conformal inference |

Two of these stand out:

- **Grohs and Rathmair** tie the stability of phase retrieval to how well connected the spectrogram is. That would link the Lean chip's phase-propagation statement to your clustering work.
- **Tachella, Chen and Davies** give conditions for when a signal model can be learned from incomplete measurements alone. That is the step none of the Voroninski papers takes.


Lee, Barber & Willett — Distribution-free inference with hierarchical data (2306.06342): Not present (you have Yonghoon Lee's 2409.13990v5 Batch Predictive Inference and Barber's 2605.29200v5 / 2202.13415v5 / 1904.06019v3 / 1905.02928v3, but not 2306.06342)
Gibbs & Candès — Adaptive Conformal Inference Under Distribution Shift (2106.00170) and Conformal Inference for Online Prediction with Arbitrary Distribution Shifts (2208.08401): Not present (you have Candès's 1904.06019v3, 1905.02928v3, 1905.03222v1, 2202.13415v5, and 2409.19712v2, but neither of the Gibbs & Candès papers)



Ran list PDFs directory structure
Ran cross-check reading list against all inventories
Ran cross-check reading list against all PDFs and gauntlet inventories
Ran verify ambiguous SplineWavelets and other titles

Almost the entire [reading-list.md](file:///d:/aghado01/ThermoMapper/issues/regression/deep-signals/reading-list.md) is already in your `PDFs` or `gauntlet` inventories (including Lee, Barber & Willett and Gibbs & Candès, which are now in `gauntlet/conformal-bayes`).

Here is the exact gap list:
Measurements"* (ICLR 2018 / OpenReview `Hy7fDog0b`) — you have the third-party reproduction `1810.10108v1` in `gauntlet/signal-processing`, not the original paper.

### Partial / Chapter-only
* **Unser & Tafti**, *An Introduction to Sparse Stochastic Processes* — only [UnserTafti-preview.pdf](file:///d:/aghado01/codex-scientiae/supellex-private/PDFs/other/UnserTafti-preview.pdf) in `PDFs/other`.
* **Chui**, *An Introduction to Wavelets* — only one chapter ([b978-0-08-052084-1.50013-2.pdf](file:///d:/aghado01/codex-scientiae/supellex-private/PDFs/SplineWavelets/b978-0-08-052084-1.50013-2.pdf), *"Spline Wavelets over R, Z, R/NZ, and Z/NZ"*) in `PDFs/SplineWavelets`.
* **Fant**, *Acoustic Theory of Speech Production* — you have two related chapters/articles in `PDFs/SplineWavelets` ([978-94-011-4657-9_3.pdf](file:///d:/aghado01/codex-scientiae/supellex-private/PDFs/SplineWavelets/978-94-011-4657-9_3.pdf) and [9783110873429.13.pdf](file:///d:/aghado01/codex-scientiae/supellex-private/PDFs/SplineWavelets/9783110873429.13.pdf)), not the 1960 monograph itself.

---

### Already in your collection
* **`gauntlet/conformal-bayes`**: `2106.01792v1` (Diquigiovanni et al.), `2306.06342v4` (Lee, Barber & Willett), `2106.00170v3` & `2208.08401v3` (Gibbs & Candès).
* **`gauntlet/signal-processing`**: `1508.02820v3` (Jaganathan, Eldar & Hassibi), `1706.04374v1` & `1903.01104v1` (Grohs & Rathmair), `2203.12513v2` & `2201.12151v4` (Tachella, Chen & Davies).
* **`gauntlet/bayesian` & `arboreal-bayes`**: `2108.11863v4` (Multivariate Lévy Adaptive B-Spline), `2007.14900v3`, `2106.03023v3`, `2211.02676v4` (Bayesian Context Trees).
* **`PDFs/Schumaker` & `PDFs/GhosalVanDerVaart2017` & `PDFs/Mallat` & `PDFs/textbooks`**: Full chapter sets including `.016`/`.018` for Schumaker, Ghosal & van der Vaart, Mallat, and `Daley2003.pdf` (Daley & Vere-Jones).
* **`PDFs/Bayesian` & `PDFs/other` & `PDFs/SplineWavelets`**: `0740817X.2012.726758.pdf` (Yuan et al. 2013), `2240405.pdf` (Stein 1981), `jasa.2011.tm11181.pdf` (Efron 2011), `proc.1975.9792.pdf` (Makhoul 1975), `tsp.2002.1003065.pdf` (Vetterli, Marziliano & Blu 2002).
* **`PDFs/voice-contours`**: `ast.5.233.pdf` (Fujisaki & Hirose 1984), `1.3037222.pdf` (Prom-on, Xu & Thipakorn 2009), `lan.1974.0010.pdf` (Sacks, Schegloff & Jefferson 1974), `fpsyg.2015.00731.pdf` (Levinson & Torreira 2015), `hawkes1971.pdf` (Hawkes 1971).
* **Bayesian context-tree (4th paper)**: you have `2007.14900v3`, `2106.03023v3`, and `2211.02676v4` in `gauntlet` (plus `Kontoyiannis2022.pdf` and `Papageorgiou2026.pdf` in `PDFs/Bayesian`), `2202.02239` (*Posterior Representations for Bayesian Context Trees*).
