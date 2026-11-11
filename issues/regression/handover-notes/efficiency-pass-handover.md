# BASS efficiency pass: what a chain's numbers do now (from the efficiency-pass session, 2026-11-04)

For the waiting-time chip, the statistical half and anyone resuming a BASS chain. Detail: `issues/regression/bass/efficiency-pass-brief.md`, under "Reports"; rows in `issues/regression/bass/efficiency-pass-results.csv`.

Commits on `main` (nothing pushed): `55e5ac0` brief; `a9f3211` levers 2–4; `4e62cf6` lever 1; `05ea2b5` lever 5; `693689f` lever 1 corrected (tables only on a plain chain); `980f258` kernels compiled fully optimized; `651cb90` no tables where the working response moves; tests-only commits between.

**Bit for bit (levers 2, 3, 4, 5).** Column kernels, categorical bitset, the commit's residual passes, τ's cached terms give the old forms' bits. With `BASS_TABLES=off` in the environment, all eight `BassBitIdentityTests` hashes of the baseline hold on the final tree, hardware intrinsics on or off. A chain's draws are the same numbers as before.

**Not bit for bit (lever 1).** In a **plain `BassChain`** of a scalar fit with at most 16 continuous predictors and at least 500 cases (and a model whose working response does not move), a one-way hinge's Gram entries against the held columns come from prefix tables. They differ from the dots by rounding: at most 4·K·ε + 8ε·log₂ n of the columns' norms, 93 ε at most observed. Over 8 chains × 20,000 sweeps at the cost grid's base, tables on and off kept **identical structure sequences**, ln densities within 6·10⁻¹³ relative; the draws' low digits differ. So the law is unchanged to that scale. Consequences:
- a chain checkpointed by an older binary and resumed by this one continues in the same law, not the same bits;
- resume within one binary stays bit for bit (the tables are a pure cache: a restored design with no tables computes the same numbers, tested over 800 sweeps);
- the two untempered Friedman pinned hashes in `BassBitIdentityTests` were re-baselined (the old values hold with `BASS_TABLES=off`); ladders, SMC, direct fits keep the dots and their pins;
- if you need the old bits for an artifact, set `BASS_TABLES=off` for the process.

**Timing.** No quiet window came: two other sessions' test hosts ran for hours of CPU throughout (4.4 h and 2.3 h by the end), 7–10 cores busy. The record's wall-clock pairs are unread; the report holds the commands. Indicative readings use `BASS_COST_CLOCK=cycles` (the thread's own cycle counter), ±40%: the base sweep about half its former time, L 1,024 0.53, direct 500 × 1,000 0.59. The statistical half's cost blocks should be rerun after this pass in a quiet window; its grid's error per unit work stands.
