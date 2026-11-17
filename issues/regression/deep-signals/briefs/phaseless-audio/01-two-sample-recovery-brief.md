# Slice 1 — two-sample Fourier recovery on the nonzero support

**Question.** For finite real audio, exactly which samples share one recoverable sign when the measurements consist of one-sample magnitudes and both powers of each two-sample frame? This supplies actual Fourier measurements to the [existing phase-propagation theorem](../../../../../lean/Lemmas/TimeFrequency/PhasePropagation.lean). The [sequence](README.md) gives the paper references and reporting convention.

**Measurements.** Let the samples be indexed from 0 to N − 1. Take every complete adjacent pair with hop one, unit window coefficients, unnormalized two-point DFT, no wrapping or padding. Retain the one-sample magnitudes and the DC and Nyquist powers of each pair. State the normalization factor separately if the proof also covers unit-energy windows.

1. Prove that the power difference on edge (t,t+1) is (xₜ+xₜ₊₁)² − (xₜ−xₜ₊₁)² = 4xₜxₜ₊₁. Equal measured powers give equal edge products. The singleton view gives equal magnitudes, including at the endpoints.
2. Prove the full ambiguity statement in both directions. Two waveforms have the same views exactly when they agree on the zero set and differ by one sign per connected component of the induced path on nonzero samples. Restrict phase propagation to this support, or reuse its edge step and show why the relative unit phase is real. Derive global recovery up to sign for one nonempty component and handle the all-zero waveform and a single sample explicitly.
3. Give an exact witness with samples (a,0,b) and (a,0,−b), ab ≠ 0. Verify equality of every measured cell and failure of a global sign. Then add the DC cell of a three-sample rectangular frame and compute its power difference. This demonstrates what a longer window contributes at this gap.

**Done when** the declarations connect the specified powers to the componentwise equivalence and the gap witness compiles. This is a deterministic real finite result. Jaganathan et al.'s generic uniqueness theorem concerns a different STFT measurement family.

## Verified result — 2026-11-16

**Proof revision `18ef361`.** In [TwoSampleRecovery.lean](../../../../../lean/Enthymemata/TimeFrequency/TwoSampleRecovery.lean), `frameDFT` sums the complex Fourier kernel with negative exponent and no normalization. `framePower` takes its squared complex modulus. The coefficients of the window are both one. `twoSample_powers` evaluates the actual bins as (a+b)² and (a−b)². Subtraction gives 4ab. `twoSampleViews_iff_products` identifies equality of the entire data with equality of singleton magnitudes and all adjacent products. The indices of a measured pair satisfy v = u+1 as natural numbers. Both endpoints have singleton observations.

`twoSampleViews_iff_componentSigns` proves both directions for every finite real waveform. The support graph is the ordinary path induced on indices where x is nonzero. At a nonzero sample, the real ratio yᵢ/xᵢ has absolute value one, so it equals +1 or −1. The existing complex `phase_step` applied to real samples makes the ratio constant along each support edge. It therefore defines one sign on each connected component. Conversely any choice of those signs preserves singleton magnitudes, products on support edges, and zero products at a gap. `global_sign_of_connected_support` gives y = x or y = −x when the support is connected and nonempty. `zero_signal_views` and `single_sample_views` handle the all-zero and length-one cases explicitly.

`gap_witness` verifies every measured cell for (a,0,b) and (a,0,−b). For ab ≠ 0 neither global sign works. Adding the rectangular three-sample DC power changes the first waveform’s power relative to the second by exactly 4ab. The posed statements hold as written under the specified normalization.

### Verification

The targeted module, its import module and the aggregate built successfully. Direct elaboration audited `framePower_dc`, `twoSample_powers`, `twoSample_power_difference`, `twoSampleViews_iff_products`, `twoSampleViews_iff_componentSigns`, `global_sign_of_connected_support`, `zero_signal_views`, `single_sample_views` and `gap_witness`. All nine printed exactly `propext`, `Classical.choice` and `Quot.sound`. The taxonomy gate passed, recording five definitions and nine theorems with zero proof holes and zero declared axioms. Replay from `lean/` at below-normal priority with `LEAN_NUM_THREADS=2`.

```powershell
lake build +Enthymemata.TimeFrequency.TwoSampleRecovery +Enthymemata
lake env lean Enthymemata/TimeFrequency/TwoSampleRecovery.lean
pwsh -NoProfile -File scripts/meta-ci.ps1 -NoBuild
```

### Doc changes owed

Record the factor four with the frame convention, the support-component ambiguity, and the three-sample DC witness beside any audio claim that uses adjacent-pair recovery.
