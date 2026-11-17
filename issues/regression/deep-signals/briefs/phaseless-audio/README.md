# Phaseless audio measurements — statements sequence

Posed 2026-11-16 for the deep-signals statements session. Start with slice 1. Each linked brief owns its proof obligations and report. The [audio approach](../approach-brief.md) calls for training from several magnitude views of an utterance. Three different objects could be identified by such data: a waveform up to sign, a set of admissible waveforms, or a probability law on that set.

| Slice | Question | Dependencies |
| --- | --- | --- |
| [1. Two-sample frames and support](01-two-sample-recovery-brief.md) | Which real sample products do actual Fourier powers determine? | Existing [phase propagation](../../../../../lean/Lemmas/TimeFrequency/PhasePropagation.lean) |
| [2. Relative signs across a weighted cut](02-weighted-cut-brief.md) | How small a change in edge products can accompany a large waveform change? | Uses the observation in 1 for the audio corollary |
| [3. Finite windows and exact recovery](03-finite-window-identifiability-brief.md) | Which relative signs do the retained Hann cells see? Which additional ambiguities survive? | Extends 1 |
| [4. Stability in the measured powers](04-measured-stability-brief.md) | What would let a cut or a restricted stability constant certify the sampled measurements? | Uses 2 and the operator specified in 3 |
| [5. Paired views and waveform laws](05-paired-model-identification-brief.md) | When does the joint observation law determine the law of a waveform modulo sign? | Uses 1 for an exact witness. The general result is conditional on the operator and signal class. |
| [6. Unpaired views and model identification](06-unpaired-model-identification-brief.md) | What can the separate image sets or marginal laws determine? | Reuses the witness from 5 |

The [conformal audio sequence](../conformal-audio/README.md) uses slice 1's gap witness to bound the size of a valid prediction set and to test synthetic-to-real latent coverage.

**Sources.** The local collection is `codex-scientiae/supellex-private/gauntlet/signal-processing`. The versioned directories contain PDFs and TeX. The following entry points identify the arguments used here.

- [Jaganathan, Eldar and Hassibi](https://arxiv.org/abs/1508.02820v3), `1508.02820v3-tex/ProblemSetup.tex` and `Uniqueness.tex`, theorem labelled `STFTUN`. The authors use a specified finite boundary convention and resolve local Fourier ambiguities generically before propagating phase through overlaps.
- [Grohs and Rathmair](https://arxiv.org/abs/1706.04374v1), `1706.04374v1-tex/phase_poincare_may15.tex`, theorem labelled `thm:main`. Its stability bound uses continuous Gaussian-window magnitudes, a weighted Cheeger constant and a data norm including derivatives. Slices 2 and 4 pose finite questions with explicit data norms.
- [Tachella, Chen and Davies](https://arxiv.org/abs/2203.12513v2), `2203.12513v2-tex/arxiv.tex`, theorem labelled `theo: multiple op`. For signal sets satisfying their dimension hypothesis, they identify the set from its separate linear images for almost every operator family when the measurements per operator exceed the model dimension plus the ambient dimension divided by the number of operators. Slices 5 and 6 ask distinct questions for fixed nonlinear powers.

**Completion.** Use the [Lean harness conventions](../../../../../lean/README.md). For each slice, append the precise statement proved or refuted, a witness that exercises its hypotheses, declaration and source paths, and the verification command and outcome. Mark conditional results with the outstanding premise. Record any necessary corrections to the proposed statements in their own brief. A counterexample to a proposed transfer is a result. Record any resulting documentation work under **Doc changes owed**. These are mathematical assignments. The numerical probes already have their own session.
