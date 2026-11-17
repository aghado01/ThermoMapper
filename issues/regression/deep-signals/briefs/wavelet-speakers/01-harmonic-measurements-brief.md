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
