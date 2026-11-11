# The bracket's T_ps anchor: a transcription error (2026-11-01, stopping round 2)

Full record: `issues/spc-BARS/stopping-theory-primer.md`, §13.14, claim C44.

What other sessions should know:
- **The error.** `src/graphs/models/potts/BwdPottsCriticalEstimate.cs` computes (e^(-1/2)/4)·ln(1 + √q). BWD 1997 eq. (4.2) has the log in the denominator: e^(-1/2)/(4 ln(1 + √q)). At q = 20 that is 0.258 against 0.089, and the code's value rises with q, which the comparison inequalities forbid.
- **Fixed 2026-11-01, at the user's call.** See git log for `fix(potts): T_ps divides by the logarithm`. The multipliers [0.05, 4] stand, re-read against measurements, so the q-only bracket at q = 20 is now [0.0045, 0.357]. `auto` schedules under `MeanEdgeDistance` move 1.06 e-folds colder. Use `BwdPottsCriticalEstimate.Bracket`.
- **It bears on the graph compiler work:** any earlier reasoning about where a graph's temperatures sit against T_ps, or the heuristic's "K̂× too cold" against it, used the old value, off by 2.89 at q = 20.
- **A usable closed form, checked:** a dense cluster's melt sits near s/(d·K_d(q)), K_d(q) = ln[(q − 2)/((q − 1)^(1 − 2/d) − 1)], with s the mean node strength; about 0.093 at s = 0.6, q = 20. Measured melts run about 1.2× hotter.
