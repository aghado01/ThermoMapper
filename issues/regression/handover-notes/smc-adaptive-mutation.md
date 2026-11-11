# SMC adaptive mutation count (chip, 2026-11-08): done

Landed on `main`, nothing pushed:
- `8ec1bc8`: core `SmcOptions.AdaptiveMutation` and `SmcStageReadout`.
- `9261c9f`: `BassSmcReplicate.StageReadout`.
- `ebe5757`: oracle law test.
- `7d0869a`: measurement harness, tables and report.
- `11436c4`: `Samplers.md`.
- `7dee840`: `SmcSpec.AdaptiveMutation` wiring.
- `03b6d13`: `BayesianRegression.md` rows.
- `a6b641b`, `7abceac`, `02453db`: report updates.

The report is `issues/regression/bass/smc-adaptive-mutation-report.md`; `artifacts/chip-smc-mutation/` has been deleted.

**Decision:** adaptive stays an option, and the fixed count stays `SmcSpec`'s default.
- On GBASS t(5) and skewed GH at a budget of 5, the adaptive count is the fixed count bit for bit. At a budget of 20, its intervals cross 1.
- On KHAOS every arm is trapped. Adaptive's lower work ratio there is lower cost at an equal bias.
- The rule's letter would have made adaptive the default on KHAOS. I did not apply it, the report says why, and it is the user's call.

**For KHAOS slice 2 (important):** at N 500 with up to 20 sweeps per stage, the Ishigami raw-ridge SMC traps most replicates at near-empty structures.
- The trapped runs read M ≈ 0.1–0.8, σ² ≈ 12, and ln Ẑ about 60 nats low, with genealogy error under 1. Only 0–2 of 16 replicates reach M ≈ 10.
- Don't rely on KHAOS SMC fits or evidence until this has a remedy, a refusal, or a flag.
- `StageReadout.LeastMoved` reads 0.02–0.06 there.
