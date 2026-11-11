# F-M2 Phase 2: the pre-registered pass, hand-over (2026-10-29)

The full report is prepended under "Reports" in `issues/regression/bass/local-move-pass-brief.md` (repo `D:\aghado01\ThermoMapper`).

- Commits on `main`, none pushed:
  - `6f1e99a`: the harness `BassMovePassTests` and `Friedman(n)`, committed before the timed run;
  - `e45cdc0`: `local-move-pass-results.csv` (160 fits plus 4 controls) and a minimum edge-case check;
  - `60f4637`: the report.
- The verdict: under the rule as written, no arm clears. The scores:

  | A1 | A2 | A3 | A4 | A5 | A6 | A7 |
  | --- | --- | --- | --- | --- | --- | --- |
  | 0.929 | 1.031 | 0.998 | 0.864 | 0.984 | 0.078 | 0.067 |

  The defaults stay D01's. An independent re-score of the CSV agrees.
- Context: the untempered fits' minimum ESS sits at the bulk-ESS floor.
  - The median is 4.96, and 119 of 120 fits read below 10.
  - M's R̂ is about 2 to 2.5, so the untempered ratios mostly compare times.
  - Tempering lifts the ESS 1.8–3.8× and brings M's R̂ to ≤ 1.55, with σ covered 40 of 40. It costs 24–43× the time, and no fit completed a round trip.
- Times: Release, about 7–9 µs per chain sweep at n ≤ 500, against Phase 1's smoke reading of 35–62. Those earlier runs' build configuration is unrecorded.
- The run-length sizing "from the winner's ESS per sweep" has no winner. Where it starts from is the main session's call.
