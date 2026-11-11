# Performance qualification, cost half: chip hand-over

Date: 2026-10-31. Full report: `issues/regression/bass/performance-cost-brief.md`, "Reports" ("The cost half measured"). Per-point table: `performance-cost-results.csv` beside it.

Commits on `main`, nothing pushed:
- `fa39a1d`: the harness, written before its timed runs.
- `c4caee4`, `b3ee63d`, `e216f8e`, `b484fca`, `abe8a7c`, `ca729db`: its revisions.
- `f7d54f7`: the report and the table.
- `50650e6`: repairs `fa39a1d`, which swept in the warm-start chip's staged brief. That chip committed its report itself (`3990c63`).

No production file changed.

Headlines:
- **At the cap the change proposal is the sweep:** 54–61% at n ≈ 1,000 and 79–96% from n = 10,000. A change is a column evaluation (2–2.5 ns per element pass) plus its cross products (0.3 ns per multiply-add in cache, 0.9 beyond it). Base sweep: 22 µs.
- **The composition is no longer dominant.** `9ed600a`'s column-mean cache holds it at about 2 µs per draw, which supersedes M1's 34–94%.
- **`PairwiseSum.Dot` already runs at the machine's vector rate.** Columns run 3–7× slower than that rate; categorical membership costs 9–65 ns per case.
- **Levers at the base:**
  - prefix sums for a one-way hinge: 16% of a sweep (cross products alone), 35% with the column skipped;
  - the column kernels: 25%;
  - the composition and τ: 18%;
  - the commit's residual: 4%, rising to 56% in direct fits;
  - GC: at most 6.5% anywhere;
  - the rebuild: 0.1% at most.
- **N15:** the construction is 97–99.8% of a σ² proposal. The capacitance count is about 7× below one Jacobi sweep.
- **Machine:** Core Ultra 5 125H (hybrid). Two quiet grid passes differ by a median 1.29× (IQR 0.98–1.51); shares hold to about 5 points.

R context: `truncdist` 1.0-2 (with `evd`) is now pinned in `r/`'s renv (`91155ea`). R BASS ran as written in a gated quiet window: 691, 825 and 525 µs per sweep (`5237b84`). Its chains are identical to those of the earlier inverse-CDF stand-in.

Artifacts root `artifacts/chip-perf/` deleted. Scratch rows were kept in the session's scratchpad only.
