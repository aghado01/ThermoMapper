# IRLS scale probe (chip-irls, 2026-11-08)

Exploratory probe behind `issues/geometry/irls-scale-brief.md` (ThermoMapper), commits 61a7fac (spread unit) and the
optimality check at the nearest data point that follows it. The committed tests are `tests/geometry/IrlsScaleTests.cs`;
this folder keeps the probe and its raw rows for reference.

- `Program.cs` / `IrlsProbe.csproj`: sweeps the median and mean over spreads 1e-12..1e6 on R^d (d 2, 5, 20; centres at
  0, 100·s, and 10; warm starts medoid, mean, first point), the sphere S^2 and S^9 and Gr(2,10), Gr(3,20) (spreads
  1e-1..1e-10), collinear and duplicated-median cases, and MoMPCA rescaled by 1e-12..1e6. Truth from outside the solver
  (plain weighted Weiszfeld at scale 1; exp_{p0}(s·m_u) on curved manifolds). It referenced a `git archive` snapshot of
  the repo with the chip's files copied in (`..\work\projects\Maths.Geometry`).
- `PinProbe.cs`: prints the pinned cases of `DataScale_ReproducesTheSolverInDataUnits_BitForBit`, built once against
  afc814d (`Snap=base`) and once against the change (`Snap=work`, `IrlsScale.Data`); the outputs agreed bit for bit,
  except the seven-frame Grassmann median after the optimality check, which became frame 0, the true median.
- `probe-rows-*.csv`: columns kind, space, dim, seed, wts, centre, warm, method (spread | data | ambient-*), s, err
  (relative to s, in the truth's u-units), iters (Exp calls under spread, Distance calls/(n+1) under data), floor.
  These rows are from the spread unit before the optimality check: `special` and the mean/first warm starts in
  `euclid` still show the vertex thrash (200 iterations) that the check removes. `agg.py <csv>` aggregates the curved
  rows.

Headline readings (the report in the brief has the full list):
- Data units: medians run to 200 iterations below s = 1e-4 (1e-3 on S^2) and miss by 2% to 11% of the spread.
- Spread unit: flat errors identical at every s (median <= 8e-8, mean <= 1e-14 relative); curved medians ~1e-8 down
  to 1e-4, then ~1e-12/s (the manifolds' Exp drops tangent steps below 1e-12), unbounded iterations from ~1e-7 (Log's
  arccos drops principal angles below ~1e-8).
- MoMPCA: equivariant while its alpha clamp [1e-6, 2-1e-6] does not bind (1e-3, 1e2 on a balanced fixture: 1e-14);
  the clamp binds at 1e-6 and 1e3 and moves the mean ~1% and the subspace ~0.01 rad.
