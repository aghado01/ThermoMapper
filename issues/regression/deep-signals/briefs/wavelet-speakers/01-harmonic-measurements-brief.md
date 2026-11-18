# 1. Which wavelet scales resolve a voiced harmonic balance?

Proposed first slice of the [speaker-signature study](../wavelet-speaker-signatures-brief.md). Settle whether an extra analysis scale can distinguish a change in pitch from a change in the fraction of energy carried by the second harmonic.

**Signal and observation.** Work on a circular block of \(N=24\) real samples. For pitch bin \(k\in\{3,4\}\), form two real cosines at bins \(k,2k\), with arbitrary phases. Give them amplitudes \(\sqrt{2(1-t)}\) and \(\sqrt{2t}\) for \(0\le t\le1\). The total mean square is one. Here \(t\) is a harmonic energy fraction. Treat it as one observable acoustic attribute. At \(t=1\), \(k\) means the synthesis parameter even though its fundamental component vanishes.

For \(b=1,2,4\), take the length-\(2b\) Haar filter with its first \(b\) taps \(+1/\sqrt{2b}\) and its last \(b\) taps \(-1/\sqrt{2b}\). Convolve circularly at **every** translation, square, then average over the block. Let \(v_j\) collect the squared Fourier gains of the three filters at bin \(j\). Prove from the finite sums that the measured vector is

\[
F_k(t)=(1-t)v_k+t v_{2k}.
\]

The cross term cancels by full-block orthogonality. State exactly where normalization, full translation averaging and circular boundaries enter. For validation, the first two gains at angular frequency \(\omega\) should be \(1-\cos\omega\) and \((1+\cos\omega)(1-\cos 2\omega)\).

**Claims to settle.**

1. Compute \(v_3,v_6,v_4,v_8\) exactly. One numerical drafting check found intersecting segments in the first two coordinates at \(F_3(t)=F_4(u)\) *after projection*, with \(0.630<t<0.632\) and \(0.238<u<0.240\). Prove existence with interior, unequal fractions and render both witnesses as actual signals.
2. In three coordinates, set \(a=v_6-v_3,\ c=v_8-v_4,\ w=v_4-v_3\). Prove \(w\cdot(a\times c)\ne0\), and show both segments have nonzero direction. This would identify \((k,t)\) from the measured energies despite arbitrary phases. The distance between their supporting lines is \(|w\cdot(a\times c)|/\|a\times c\|\). Derive an explicit positive lower bound for use in slice 2. This is a lower bound on segment distance.
3. Include a one-bin STFT or full-block Fourier control with its actual window and retained bins specified. Show how a resolved harmonic-power observation distinguishes the *interior collision pair*. This pins the loss of information to the selected features.

**Closure.** Supply the filter definitions, finite-waveform identities, collision, and separating constant. If a proposed comparison fails under these conventions, give the exact response and corrected condition. The test concerns these scales and this signal family. Relate the necessary orthogonality facts to the [finite measurement work](../phaseless-audio/03-finite-window-identifiability-brief.md).


## Slice 1 result and proof report — 2026-11-17

Two Haar scales give identical measurements for different interior harmonic fractions on this 24-sample generator. Adding the half-width-four filter identifies the synthesis bin and the harmonic fraction for every phase pair, including at the fraction endpoints. At \(t=1\), the recovered bin indexes the synthesis family although the fundamental has zero amplitude.

For sample \(i=0,\ldots,23\), the waveform is
\[
 x_i=\sqrt{2(1-t)}\cos(2\pi ki/24+\phi)
      +\sqrt{2t}\cos(2\pi(2k)i/24+\psi).
\]
The samples form \(\mathbb Z/24\mathbb Z\). The mean square of this waveform is one. Convolution uses delays \(d\) and \(d+b\) modulo 24, for \(0\le d<b\), with the second tap negated. Dividing the resulting difference sum by \(\sqrt{2b}\) gives unit filter energy. Squaring and averaging over all 24 translations cancels the cross term by finite Fourier orthogonality. The phase drops out and each measured coordinate equals \((1-t)v_k+t v_{2k}\).

In the table \(r=\sqrt{2}\). The columns are squared filter gains, in order of half-width.

| Bin | \(b=1\) | \(b=2\) | \(b=4\) |
| --- | --- | --- | --- |
| 3 | \(1-r/2\) | \(1+r/2\) | \(2+r\) |
| 6 | \(1\) | \(2\) | \(0\) |
| 4 | \(1/2\) | \(9/4\) | \(9/8\) |
| 8 | \(3/2\) | \(3/4\) | \(3/8\) |

The first two response polynomials are \(1-z\) and \((1+z)(1-z^2)\), with \(z=e^{-i\omega}\). After tap normalization their squared gains are \(1-\cos\omega\) and \((1+\cos\omega)(1-\cos 2\omega)\). The third polynomial is \((1+z)(1+z^2)(1-z^4)\). These identities come from the finite tap sums. The earlier finite-measurement work represents power as the squared modulus of a specified Fourier row. Here its rectangular row supplies the control, and character orthogonality evaluates the effect of averaging translated Haar powers. [Finite waveform measurements](../../../../../lean/Enthymemata/TimeFrequency/FinitePowerMeasurements.lean)

**Collision witnesses.** Use the displayed waveform at \((k,t)=(3,(3+r)/7)\) and \((4,(9-4r)/14)\). Setting all phases to zero gives two concrete signals. The theorem permits independent arbitrary phases on both signals. Their first two measured coordinates agree although their fractions lie in the disjoint intervals \((0.630,0.632)\) and \((0.238,0.240)\).

**Margin for Slice 2.** With \(a=v_6-v_3,\ c=v_8-v_4,\ w=v_4-v_3,\ n=a\times c\), Lean proves
\[
 w\cdot n=(2+r)/8,\qquad n\cdot n<256(w\cdot n)^2.
\]
Every difference between points on the two lines has the same projection onto \(n\). Cauchy–Schwarz gives the Euclidean separation \( \|F_4(u)-F_3(t)\|_2>1/16 \) for **all real** \(t,u\). Slice 2 can therefore take \(h=1/16\) in these normalized energy coordinates. Both directions are nonzero, so equal measurements within one pitch force equal fractions. Across pitches, the separation excludes equality.

**Fourier control.** A rectangular window of 24 ones retains only positive-frequency bin 3. With power normalized as \(|\mathrm{DFT}(x)_3|^2/24^2\), the bin reads \((1-t)/2\) on the bin-3 synthesis family and zero on the bin-4 family. The first collision witness has \(t<1\), so this control strictly distinguishes the pair for every phase choice.

**Lean declarations.** The aggregate reaches the following modules in namespace `TimeFrequency.Haar`.

| Module | Main results |
| --- | --- |
| [HaarHarmonics](../../../../../lean/Enthymemata/TimeFrequency/HaarHarmonics.lean) | `tone_cos`, `tone_inner`, `signal_energy` |
| [HaarMeasurements](../../../../../lean/Enthymemata/TimeFrequency/HaarMeasurements.lean) | `filtered_energy`, `gain_formulas` |
| [HaarGains](../../../../../lean/Enthymemata/TimeFrequency/HaarGains.lean) | `bin_gains`, `measured_signal` |
| [HaarGeometry](../../../../../lean/Enthymemata/TimeFrequency/HaarGeometry.lean) | `two_scale_collision`, `normal_certificate`, `cross_pitch_gap`, `three_scale_identification` |
| [HaarFourierControl](../../../../../lean/Enthymemata/TimeFrequency/HaarFourierControl.lean) | `fourier_control_eq`, `bin_three_signal`, `fourier_distinguishes_collision` |

**Verification.** Full `lake build` passed (3829 jobs, BelowNormal process priority, `LEAN_NUM_THREADS=2`). All 20 public theorem axiom listings use only subsets of `[propext, Classical.choice, Quot.sound]`. The five new modules emitted no warnings. The taxonomy gate `lean/scripts/meta-ci.ps1 -NoBuild` passed with zero proof holes and zero postulated axioms in all five modules. They remain in Enthymemata pending semantic review for promotion.

## Doc changes owed

In the next documentation batch, record the two-scale collision and the three-scale margin alongside the generator, averaging rule and feature units. The identified acoustic quantity is the harmonic fraction in a controlled block.
