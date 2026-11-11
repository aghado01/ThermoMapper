# Flat-histogram loose ends (chip fh4): hand-over, 2026-10-31

Full report: `issues/regression/bass/flat-histogram-loose-ends-brief.md`, under Reports (`a712cf8`).

- Commits on `main`, none pushed:
  - `9d15376`: pins of today's toy runs and gap fits.
  - `72788a2` + `0454bf8`, item 1: `ExtendTopUntil` and `WindowTopMargin`, off by default.
  - `25e5124`, item 2: sweeps folded into both hosts' `ComparativeSweep`.
  - `f093987`, item 3: test before its run.
  - `0ec16d3`, item 4: `TiedNodes` (engine result and checkpoint, plus the BASS pilot readout), `Wham.NodeDraws` in the adapter, and `SamcLearning.Infeasibility` naming the refusal's segment.
- Item 1, seed 8017, one run each:
  - **Today's window** reproduces the qualification: −841.36, a1 65.9%, learning 29% inert, ESS 34,256, 47 round trips.
  - **Extended** (until 50k): the top climbed 26.6 nats to −751.65 in 4,000 iterations, into a sparse tail. Learning was 0.1% inert, but production made 3 round trips (two chains none), with ESS 7,462. ln Ẑ −841.45, a1 69.6%.
  - **Margin** (3 sd, unheld): top −767.73, start tied, 1.1% inert, 24 round trips, ESS 14,543, ln Ẑ −841.11, a1 65.5%. a1 + x6·x8 parts across chains (p < 0.001).
  - The open variant: a quantile or capped extension rule.
- Item 2: the folded claims' rows are bit for bit at `SEED_SWEEP` = 1. Two warm-start claims fail at their first chain-axis seed (ratio 1.254 double well, 1.188 BASS). This predates the chip.
- Item 3, SMC from the prior at b = 0.10, against the flat histogram:
  - **Linear form:** 5.6 nats low at 5 sweeps (families 188–217). At 20 sweeps 3.9 low; one replicate came 1.0 short on 69 families.
  - **Integrated form (D05):** 0.4 nats low at both sweep counts (families 442–855); within 0.01 at b ≤ 0.07.
  - The integrated reference is the flat histogram's draws reweighted by that form's tempered marginal.
- Suites at `0ec16d3`: Samplers 282, Regression 1,409, BASS oracle classes 90, seed 8017's pin, all passing; `BassBitIdentityTests` hold.
- A samplers compile with an empty artifacts root (after a failed edit script) wrote bin/obj under `projects/` at 21:07. It was removed.
- `artifacts\chip-fh4\` is deleted. docs/Maths and the checklists are untouched. No default recommended.
