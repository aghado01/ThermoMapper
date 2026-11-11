# Stage-4 oracle modes chip: findings for the main session (2026-10-28)

From the chip running `issues/regression/bass/stage4-oracle-modes-brief.md`. The full report, with every emitted field and the checks' numbers, is under that brief's "Reports". Commits: `6acc2ef` (A), `e76da3c` (B), `bf608b5` (C), `6582edb` (D), `4a563f7` (E), `4d7ebdd` (report). `r/scripts/bass_stage4_checks.R` prints 39 PASS; stage-1 checks are unchanged.

- **S6.** Fixture A (N = 27, M_max = 3, 2,263 states) is the chain fixture, and C (a three-level factor, 1,949 states) the mixed one. The posteriors have exact ties: x and 1 − x with the intercept, and a factor main effect on D and one on its complement. Compare tied states' visits summed.
- **S8.** On the integrated form, the slope at b = 0 is infinite under Zellner–Siow, because it carries −(k/2)E[1/τ]; the oracle emits null there. Even at fixed g it is about −N/(2c), −13,471 on fixture A. A ladder on that form needs geometric rungs near 0. Under flat the form is undefined for b ≤ (1 − 2c)/N.
- **S4/S5.** R's `s2.rate` cancels: compare against `s2_rate_direct`.
- **A decision for the user.** D01 leaves the mixed (J, L) order prior open, and R's has none (its zero-order part scores `lbmcmp = 0`). The oracle reads it as uniform per part on {0..cap}. If the user decides otherwise, the oracle changes two lines (`cont_lp`, `cat_lp`).
