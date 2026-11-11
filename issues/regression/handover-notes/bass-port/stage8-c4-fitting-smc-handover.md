# Stage 8, C4 (fitting SMC): hand-over, 2026-10-29

Full report: `issues/regression/bass/stage8-fitting-smc-brief.md`, under Reports (commit `9036a38`).

- Commits on `main`, none pushed: `d0a0fd8` (DrawSample = Akinshin's weighted type 7), `9e8c349` (particles + kernel), `b0978d5` (SmcSpec on BassSampling.Smc, BassSmcFit, BassSmcReadout, evidence, PointwiseLogDensity.Weights, stage budget + bit-identical resume), `73352f0` (qualification behind BASS_QUALIFICATION), `9036a38` (report).
- Suites: Samplers 206/206, Regression 1169/1169, Maths.Oracle.Tests.Bass* 71/71. The 8 BassBitIdentityTests hashes hold.
- Shared code touched outside the BASS files: `src/maths/regression/Pointwise.cs` (CRLF kept) and `PredictiveChecks.cs` (optional weights, null = old path), `BassModel.Create` (reads SmcSpec.Form), `BassMovePassTests`' data hash re-pinned (config JSON gained Smc).
- Qualification: mixed and functional replicates find the best MCMC chains' representation at S10's sizes; Friedman's do not (−897…−915 against −892.4), and replicates disagree 13–66 nats. Every Friedman replicate ends in 1–2 families, so the genealogy error degenerates to √ln 2 = 0.83 and the between-run error is up to 21× larger. No default recommended.
- For the main session: docs (DrawSample's rule, weighted checks, the SMC fit); BassSummary/BassTraces are not SMC-aware (stage 9's files); TemperToPosterior throws at MaxStages, so the adapter runs its loop itself.
