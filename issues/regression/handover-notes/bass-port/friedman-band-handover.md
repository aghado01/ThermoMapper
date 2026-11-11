# Friedman's band and modes (R1–R4): hand-over, 2026-10-29

Full report: `issues/regression/bass/friedman-band-brief.md`, under Reports ("R1–R4 measured").

- Commits on `main`, none pushed: `5673978` (harness `tests/regression/BassFriedmanBandTests.cs` + `BassMixing.Fixture.Basin`), `f4577a2` (report). No production change.
- R1: ℓ₁ at fixed b is bimodal at b ≈ 0.23–0.28 on all four chains (BC > 5/9, D > 2), and the split holds within each M (M = 6, 7 lumps 120–190 nats apart); M itself slides unimodally 8 → 6. A second split at b ≈ 0.12–0.15 is between M classes (6 vs 4–5). Weights rest on ESS 2–19 per chain. Λ̂: band 2.0–2.4 (self-sized), 2.7–3.0 (dense).
- R2: ln Ẑ replicate sd 8.4 / 7.1 / 6.5 at 1k / 4k / 16k particles; pooled −875.8 / −874.4 / −861.8 (still rising). Families collapse inside the band at every N (1–10 left at b = 0.32). The 16k best replicate holds a basin (4 × x1·x2 + 2 × x3 + x4 + x5) no MCMC chain visits. Basin log masses disagree across replicates by up to 17 nats.
- R3: 20 sweeps at 1k: ln Ẑ −867.0, sd 7.2; families at b = 0.32 3–10; per-sweep acceptance unchanged.
- R4: angle floor refuses 0% everywhere (no 1e-14/1e-18 reruns); support refuses 15.6% of births; halving the minimum to 10 leaves R̂ of M at 2.13.
- A concurrent session began uncommitted edits to BASS production/test files at 21:06, after this chip's last build (20:52); no reading depends on them, none are in these commits.
- `artifacts\chip-friedman\` deleted.
