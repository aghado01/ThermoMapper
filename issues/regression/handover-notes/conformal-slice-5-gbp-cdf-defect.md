# GBP mixing with a drift: the tabled CDF near x = 0 (found by the conformal layer's slice 5, 2026-11-08)

For the noise laws' owner (`src/maths/distributions/NormalVarianceMeanMixture.cs`). Not fixed by the conformal chip.

- **What.** `GbpMixing.Noise(s², m)` (a `NormalGeneralizedBetaPrime`) reads ln F about 1e-3 off near x = 0 when m ≠ 0. Its density there is right to 12 digits.
- **Reference.** GBP(1, 1, 1.5), m/s = −0.2, so the same law at every σ². ln F from ∫ Φ((x − m v)/(s√v)) g(v) dv, by 260k-point Simpson in u = ln v:

| x / s | C# ln F | Reference ln F |
|---|---|---|
| −0.1 | −0.674213 | −0.674379 |
| −0.001 | −0.550488 | −0.551580 |
| 0 | −0.549285 | −0.550321 |
| 0.001 | −0.548083 | −0.549064 |
| 0.1 | −0.441803 | −0.441820 |
| 0.5 | −0.207449 | −0.207454 |

- **Its local slope.** F(0.001) − F(−0.001) reads 0.001389 against 2·0.001·p(0) = 0.001453, so the slope is about 4% low.
- **Without a drift** F(0) is ½ exactly. Only a = 1 was checked against the integral; a = 0.5 and a = 2 were not. The GIG mixtures are within their stated 1e-5.
- **Reproduce.** `EngineCandidateTests.TheLogCdf_UnderGbpMixingWithADrift_IsItsDensitysIntegral` (eef6f47) is skipped until the fix. Unskip it to check: its assertions are the correct references.
- **Effect.** BARS's GBP noise reads its predictive distribution's log CDF from this table. Its density, and so its conformal sets, are unaffected.
