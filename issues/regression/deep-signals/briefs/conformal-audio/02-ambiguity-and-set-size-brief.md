# Slice 2 — the size of a valid set when short windows lose a sign

**Question.** At fixed observable, units, nominal coverage and sampling protocol, how much precision can an extra view buy? Use [slice 1](01-conversation-band-brief.md) for the meaning of marginal coverage. The following lower bound applies to every predictor with that coverage.

1. Give the real waveforms (1,0,1) and (1,0,−1) equal probability. Use `gap_witness` and the rectangular DC convention to verify that all singleton magnitudes and both powers of every adjacent two-sample frame agree, while the three-sample DC targets are 4 and 0. Retain unnormalized powers in both experiments.
2. A prediction set C may depend on the short views, independent calibration data and independent algorithmic randomness. Prove that the new fair sign remains independent of the two membership indicators. By expanding the two cases, show P(Y ∈ C) = E[card(C ∩ {0,4})]/2. Hence 1−α marginal coverage requires at least 2(1−α) accepted target values on average. At 95% this is 1.9. Use counting measure on the two target values, even when C is an interval.
3. Show sharpness. At 95%, return both values with probability 0.9 and a singleton otherwise, using independent randomness. Check size and coverage exactly. Then keep the same target and give the predictor one other observable that distinguishes the signs. For example, verify the unnormalized frequency-one power of the three-point DFT is 1 on (1,0,1) and 3 on (1,0,−1). The DC target can then be predicted by a singleton with probability one. It remains held back.

**Acceptance.** Verify the lower bound, its attainment and the singleton control with the extra measured bin. State the independence premise and the reference measure. This is a two-point sharpness result. The separate `hoff_cb1_expectedSize_le` theorem uses its own hypothesis comparing coverage across orderings.

[Finite DFT and gap proof](../../../../../lean/Enthymemata/TimeFrequency/TwoSampleRecovery.lean) · [Existing expected-size theorem](../../../../../lean/Enthymemata/ConformalValidity.lean)
