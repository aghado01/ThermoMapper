# 2. Separate families of utterances under pitch, articulation and channel changes

Proposed follow-up to [the harmonic measurements](01-harmonic-measurements-brief.md). Give an acoustic meaning to separation of two voice families before fitting an encoder.

**Experiment.** Fix a comparable context \(u\), such as a controlled sustained vowel. Model source \(s\)'s harmonic fraction by the closed interval \(I_{s,u}=[\theta_{s,u}-r_{s,u},\theta_{s,u}+r_{s,u}]\subset(0,1)\). Both pitches and arbitrary phases are allowed. Overall gain is positive. Normalize each filter energy by raw block mean square. The clean measured family is the union of \(F_3(I_{s,u})\) and \(F_4(I_{s,u})\). This generator allows changes within one voice as well as differences between voices.

**Proposed boundary.** Let \(g=\mathrm{dist}(I_{1,u},I_{2,u})\), \(\lambda=\min_k\|v_{2k}-v_k\|\), and let \(h\) be the positive cross-pitch line-distance bound from slice 1. For the clean three-scale features prove a between-family margin of at least \(\min(\lambda g,h)\). For the same pitch a target difference \(|t-t'|\) has measured distance \(\|v_{2k}-v_k\||t-t'|\). At fixed matched context, prove that disjoint compact intervals suffice, and that overlapping intervals admit identical waveforms from the two sources. Choose narrow, disjoint intervals around the slice 1 collision parameters to show that removing a scale can destroy the margin.

**Noise and inverse use.** Given feature error \(\eta\) in the Euclidean norm after normalization, prove correct nearest-family assignment if \(2\eta\) is strictly less than the certified margin. For an off-range recording and approximate solver, state the error budget contributed by measurement noise, mismatch to the synthesis family, and fit residual. Tie at least one bound to the waveform, for example by propagating an additive block RMS error through convolution, squared energy and normalization with a positive clean RMS floor. Compare the resulting restricted separation with the [fixed-generator isometry](../../../../../lean/Enthymemata/GenerativePrior/IsometryCondition.lean) and [measured residual bounds](../phaseless-audio/04-measured-stability-brief.md).

**Channel obstruction.** At a fixed pitch let the channel apply independent positive amplitude gains \(c_1,c_2\) to the two harmonics. The observed fraction is \(c_2^2t/(c_1^2(1-t)+c_2^2t)\). Prove a collision between two distinct interior fractions under allowable channel choices, including the normalized wavelet observations. For a bounded log gain ratio, derive the remaining margin in log odds from the shift \(2\log(c_2/c_1)\). Give separate assumptions for a shared channel and independently varying channels.

**Closure.** Report a usable noise tolerance for one explicit pair of intervals and a physical collision for an excluded condition. Make clear which parameters are recovered for a single turn and which ambiguities require additional turns or pairing.

## Slice 2 result and proof report

For the 24-sample generator at a fixed context, one gain-normalized observation identifies the pitch and the harmonic fraction reaching the microphone. With only overall gain variation and known, disjoint source intervals, it also identifies which of the two families produced the block. Any positive overall gain and any two phases give the same features at a given pitch and fraction. If the intervals overlap, choosing the same fraction, pitch, phases and gain produces an identical raw block in both families. Associating unlabeled observations across different contexts requires additional information, such as the linked turns studied in Slice 3.

### Family distance and assignment

The normalized image of the synthesis family with fractions in \(I\) is exactly \(F_3(I)\cup F_4(I)\). Write \(s_p=\|v_{2p}-v_p\|_2\) and \(\lambda=\min(s_3,s_4)\). On either pitch line, \(\|F_p(t)-F_p(u)\|_2=s_p|t-u|\). For nonempty closed intervals \(I=[a,b]\) and \(J=[c,d]\), their attained distance is \(g=\max(0,c-b,a-d)\). Combining the same-pitch identity with the Slice 1 cross-pitch certificate gives the lower bound
\[
M=\min(\lambda g,1/16).
\]
Disjoint intervals have \(M>0\). A closest point in the union of the two compact feature families exists. With normalized Euclidean feature error \(\eta\) satisfying \(2\eta<M\), every minimizer belongs to the true family. [Family image, attained gap and nearest assignment](../../../../../lean/Enthymemata/TimeFrequency/HaarFamilies.lean)

For **\(I=[0.238,0.240]\) and \(J=[0.630,0.632]\)**, \(g=0.390\) and \(\lambda\ge 1/2\), so the certified margin is \(1/16\). Feature error strictly below \(1/32\) suffices. Dropping the third measurement makes the Slice 1 witnesses in these very intervals agree on both retained coordinates, even with independently chosen phases. [Two-scale witnesses](../../../../../lean/Enthymemata/TimeFrequency/HaarGeometry.lean)

### Error before and after measurement

Let \(R=\operatorname{RMS}(x)>0\) for the clean block and \(e=\operatorname{RMS}(x-y)<R\) for its perturbation. The observed block \(y\) then has positive RMS. For Haar half-width \(b\), the filtered RMS is at most \(\sqrt{2b}\) times the input RMS. Applying this bound to both blocks and their difference bounds the error of the normalized squared band energy by \(8b\,e/R\). Combining half-widths \(1,2,4\) gives
\[
\|N(x)-N(y)\|_2\le 8\sqrt{21}\,e/R.
\]
Thus **\(e/R<1/(256\sqrt{21})\)**, about **0.0852% of clean block RMS**, suffices for correct assignment in the example. An absolute error limit follows from a known floor \(R\ge R_0>0\) by requiring \(e<R_0/(256\sqrt{21})\). These are sufficient tolerances. [Convolution, normalization and raw waveform tolerance](../../../../../lean/Enthymemata/TimeFrequency/HaarFamilyNoise.lean)

For an off-family recording, let \(m\) bound the feature distance from a clean model point to the noiseless recording, \(\eta\) the distance from that recording to the observation, and \(r\) the fitted residual against the observation. The distance from the clean model point to the fit is at most \(m+\eta+r\). If the fitted residual is within \(\epsilon\) of the best residual, this becomes \(2(m+\eta)+\epsilon\). A model-to-model bound \(E<1/16\) forces the fitted pitch to agree and gives \(|t-\hat t|\le E/s_p\). A bound \(E<M\) excludes the competing family.

The fixed-generator isometry theorem assumes control of secant inner products by a linear map and obtains waveform-range injectivity. Here the normalized energy map discards phase and gain. The proved inverse control concerns pitch and harmonic balance. The earlier measured-residual argument converts a residual into parameter error using a lower gain in the same measurement norm. Here \(s_p\) supplies that gain on each pitch line, the cross-pitch bound excludes the other line, and \(m\) accounts for model mismatch. [Fixed-generator isometry](../../../../../lean/Enthymemata/GenerativePrior/IsometryCondition.lean) · [Measured-residual bounds](../phaseless-audio/04-measured-stability-brief.md)

### Channel collision and shared-channel comparisons

Positive amplitude gains \(c_1,c_2\) give raw power \(P=c_1^2(1-t)+c_2^2t\) and observed fraction \(u=c_2^2t/P\). The resulting samples equal \(\sqrt P\) times the generator at \(u\) with the original pitch and phases. In particular, source fraction **\(1/5\)** through gains **\((1,2)\)** and source fraction **\(4/5\)** through **\((2,1)\)** have the same power \(8/5\), the same observed fraction \(1/2\), and identical raw samples for any common pitch and phases. A single such recording cannot determine the source fraction.

Writing \(\ell(t)=\log(t/(1-t))\), the observed log odds are \(\ell(t)+2\log(c_2/c_1)\). A common unknown gain ratio across two recordings preserves their log-odds difference exactly. Absolute source fractions require calibration or a known reference. If the ratios vary independently and each satisfies \(|\log(c_2/c_1)|\le B\), a clean log-odds gap \(G\) leaves observed separation at least \(\max(0,G-4B)\). A positive certificate in these units requires \(G>4B\). [Waveform identity, collision and channel bounds](../../../../../lean/Enthymemata/TimeFrequency/HaarChannels.lean)

### Lean declarations and verification

The time-frequency aggregate imports all three modules in the `TimeFrequency.Haar` namespace. The proof reuses the explicit pitch directions from Slice 1. Direct Mathlib imports cover Euclidean three-space and block RMS, minimization on compact sets, the finite-sum Cauchy–Schwarz inequality, and the real logarithm.

| Module | Main public declarations |
| --- | --- |
| [HaarFamilies](../../../../../lean/Enthymemata/TimeFrequency/HaarFamilies.lean) | `normalized_family`, `overlap_waveform`, `same_pitch_distance`, `intervalGap_attained`, `family_separation`, `nearest_assignment`, `off_range_fit`, `approximate_fit`, `fitted_parameters`, `collision_family_margin`, `collision_family_two_scales` |
| [HaarFamilyNoise](../../../../../lean/Enthymemata/TimeFrequency/HaarFamilyNoise.lean) | `filtered_energy_le`, `normalized_band_change`, `normalized_waveform_error`, `collision_family_waveform_tolerance` |
| [HaarChannels](../../../../../lean/Enthymemata/TimeFrequency/HaarChannels.lean) | `channel_as_signal`, `normalized_channel`, `channel_collision`, `logOdds_observed`, `shared_channel_difference`, `independent_channel_gap` |

Full `lake build` passed (3832 jobs, BelowNormal priority, `LEAN_NUM_THREADS=2`). The three new modules compiled without warnings. All 38 public theorem axiom listings across those modules, plus the newly exported `directions` identity, use only `propext`, `Classical.choice` and `Quot.sound`. The taxonomy gate `lean/scripts/meta-ci.ps1 -NoBuild` passed and marked all three modules proof-closed in Enthymemata.

## Doc changes owed

In the next audio documentation batch, give the interval example and its strict relative RMS tolerance, and state whether the harmonic gain ratio is calibrated, shared, or independently varying when interpreting a source balance.
