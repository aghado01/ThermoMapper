# Slice 1 — two-sample Fourier recovery on the nonzero support

**Question.** For finite real audio, exactly which samples share one recoverable sign when the measurements consist of one-sample magnitudes and both powers of each two-sample frame? This supplies actual Fourier measurements to the [existing phase-propagation theorem](../../../../../lean/Lemmas/TimeFrequency/PhasePropagation.lean). The [sequence](README.md) gives the paper references and reporting convention.

**Measurements.** Let the samples be indexed from 0 to N − 1. Take every complete adjacent pair with hop one, unit window coefficients, unnormalized two-point DFT, no wrapping or padding. Retain the one-sample magnitudes and the DC and Nyquist powers of each pair. State the normalization factor separately if the proof also covers unit-energy windows.

1. Prove that the power difference on edge (t,t+1) is (xₜ+xₜ₊₁)² − (xₜ−xₜ₊₁)² = 4xₜxₜ₊₁. Equal measured powers give equal edge products. The singleton view gives equal magnitudes, including at the endpoints.
2. Prove the full ambiguity statement in both directions. Two waveforms have the same views exactly when they agree on the zero set and differ by one sign per connected component of the induced path on nonzero samples. Restrict phase propagation to this support, or reuse its edge step and show why the relative unit phase is real. Derive global recovery up to sign for one nonempty component and handle the all-zero waveform and a single sample explicitly.
3. Give an exact witness with samples (a,0,b) and (a,0,−b), ab ≠ 0. Verify equality of every measured cell and failure of a global sign. Then add the DC cell of a three-sample rectangular frame and compute its power difference. This demonstrates what a longer window contributes at this gap.

**Done when** the declarations connect the specified powers to the componentwise equivalence and the gap witness compiles. This is a deterministic real finite result. Jaganathan et al.'s generic uniqueness theorem concerns a different STFT measurement family.
