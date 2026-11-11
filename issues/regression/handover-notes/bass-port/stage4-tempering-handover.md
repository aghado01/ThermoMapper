# Stage-4 tempering (S8, chip D): hand-over for S10 and the main session (2026-10-28)

From the chip running `issues/regression/bass/stage4-tempering-brief.md`. Full report: that brief, "Reports". Commits: `58fdeb2` (tempered run), `8c1bb66` (evidence), `fb22f73` + `f2d7149` (review pass 1), `aaf46b6` (pass 2), `25cb47f` (report). Regression 756/756, Samplers 146/146, BASS oracle 35/35.

- **Entry points.** `BassSampling.Tempering` (`TemperingSpec`, null = untempered); `BassLadder` (public, one per chain); `BassFit.Evidence` (`BassEvidence`, `BassEvidencePriors`); `BassDiagnostics.Tempering` (`BassTemperingReadout`).
- **Departure from S6.** A swapped state's design travels with it (an exact copy), instead of a whole reload; only reference draws are loaded whole. Internal `BassStructureMoves.Holds/DesignFor/Adopt`.
- **For S10 / the default-levels decision.** At the default 8 levels, the N = 175 functional fixture (`BassRunTests.FunctionalData`, maxBasis 3) had adjacent rejections 0.96–1.00; the stepping-stone chains spread 122–183 against errors of 2–5. At 16–32 levels both forms agree near 180–186. The reported Monte Carlo error under-covers when rungs barely overlap. Recommendation in the report: levels auto-picked from the burn-in's barrier (`TemperingLadder.RecommendedLevels`).
- **Not mine, edited:** `tests/oracle/BassTargetParityTests.cs` (`FixtureB` private → internal).
