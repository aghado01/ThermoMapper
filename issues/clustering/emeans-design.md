# e-means: design note

Source: Bessa, Mazaré, Szilvasy, Touvron, Picard, Douze, Jégou, "e-means: outperforming k-means (in one-tenth the
epochs)", arXiv:2611.09501v1 (`supellex-private/staging/2611.09501v1`). Module: `src/clustering/geometric/emeans`,
project `Clustering.Geometric.EMeans`, tests `tests/clustering/geometric/emeans`.

## What is built

- `EMeans.Cluster(data, config, initialCentroids?)`: Algorithm 1 with the cosine step, dead-cluster revival, the
  slow EMA, and the three centroid views (EMA, slow EMA, last-epoch mean). The result keeps the decayed statistics
  `(Counts, Sums, Scatter)` and the slow pair; labels are a final nearest-centroid pass over the selected view.
- `EMeansScalingLaw`: the batch and step laws of Table 2, with `Extrapolations` naming each bound of the fitting
  grid a problem crosses. `Auto` uses the laws and logs the extrapolation in the manifest.
- `KMeansPlusPlusSeeding.Seed` in the k-means module: D² seeding on a pinned `Xoshiro256PlusPlus` stream, used by
  `EMeansInit.KMeansPlusPlus`. `KMeansPlusPlus.Cluster` still seeds from `System.Random` through its own private
  copy; moving it over changes every existing caller's stream and is left for a decision.

## Departures from the paper's text

1. **Trailing partial batch is trained on.** The paper counts `M = T·⌊n/B⌋` updates, which drops up to `B − 1`
   points per epoch and leaves them without a last-epoch assignment. Here `M = T·⌈n/B⌉`. At the paper's sizes the
   remainder is negligible; at a few hundred points it is not.
2. **A cluster with no mass holds its position.** The step law caps `α₀` at 1, so `γ = 0` on the first update, and
   a cluster that receives no point in that batch has `N = 0, S = 0` and a centroid of `0/0`. The laws reach the cap
   often (`n = 10⁶, k = 1024, T = 4` gives 1.26 before the cap). `Decay` carries such a cluster at a floor mass of
   `1e-200` at its previous position. Test: `WithoutRevival_AClusterWithNoMassHoldsItsPosition`.
3. **A split is offset along a random direction by a fraction of the donor's own radius.** The paper scales the
   donor's coordinates by `1 ± δ` in a checkerboard. That offset is zero for a centroid at the origin (both children
   coincide and the revived cluster stays empty), it grows with the centroid's distance from the origin rather than
   with the cluster's size, and a checkerboard direction is orthogonal to data that varies only along duplicated
   coordinates. Here the offset is `δ·r·u`: `r` the donor's root-mean-square radius from a decayed scatter kept
   beside `(N, S)`, `u` a random unit direction. Tests: `Revival_SplitsADonorWhoseCentroidSitsAtTheOrigin`,
   `Revival_SplitsATightDonorFarFromTheOriginWithinItsOwnRadius`, `Run_IsEquivariantUnderTranslationAndScale`.
   A first version sized the offset by the data's overall spread; the seeded sweep showed it stranding children
   outside a tight clump that another centroid already covered.
4. **A split in the last epoch restarts the last-epoch records of both slots.** Records made before the split
   describe the donor's whole cluster, and would pull its last-epoch mean back towards the undivided centroid.

Readings of ambiguous text: the step schedule follows the equation (`t` from 0, `π·t/M`), so the last update sits
just above the floor; "nearest power of two" is taken in log space; donors are drawn among live clusters only.

## Measured below the fitting grid (2026-11-14)

Probe `artifacts/probes/emeans-small-scale` (code, `results.md`, `results-eighth.md`, `runs.csv`; local, not
tracked). Sixteen cells, all below the grid: iris (k=3), landsat (k=6, 24), isolet (k=26), and balanced mixtures,
Zipf-sized mixtures and uniform cubes at `(n, d, k)` = (2000, 2, 8), (2000, 16, 8), (10000, 8, 32),
(20000, 32, 64). Five seeds per cell, one data draw per synthetic cell. Baseline: Lloyd from k-means++ seeds until
no label changes (cap 100 updates; median 4 to 100 updates to converge). e-means: `Auto` batch and step, default
revival, last-epoch mean. Metric: training distortion; a pass is one epoch or one Lloyd update.

- **From the same k-means++ seeds, `Auto` does not lose to Lloyd.** After 3 passes e-means is below Lloyd on 79 of
  80 seed-cells. After 10 epochs its median is within 0.2% of Lloyd at convergence, or below it, in all 16 cells.
  After 25 it is at or below Lloyd at convergence, to within 0.01%, in all 16 (strictly below on 64 of 80
  seed-cells; most of the rest are ties at the same optimum). At 3 epochs it is level with Lloyd at 25 updates:
  better in 7 cells, worse in 8, tied in 1, all by under 1% except two cells in its favour.
- **What the laws resolve to here.** The step law sits at its cap (`α₀ = 1`) in 14 of 16 cells at `T = 10`, and the
  batch law gives `B ≈ n/8`, about 32 points per cluster per batch. Against a 60-point grid of explicit `(B, α₀)`
  at `T = 10` from random points (median of 2 to 3 seeds, so the grid's best is optimistic), `Auto` is within 0.3%
  of the best in 11 cells, 0.6 to 2.8% behind in 3, and 10.8% and 18.4% behind in two mixture cells where smaller,
  noisier batches reach better optima across a whole block of the grid. That block costs 0.2 to 6% on the uniform
  cube of the same shape, and no single `(B ≤ 128, α₀)` ends at or below `Auto` in more than 11 of 16 cells. Full
  batch is the worst batch size in 14 of 16 cells.
- **The paper's divided-by-8 batch law, run directly, shows no consistent gain** over the undivided law here: it
  wins on 0 to 5 of 5 seeds depending on the cell, in both directions.
- **Random-point initialization is the weak part at small `k`.** From random points e-means stays 7.8% above Lloyd
  on the well-separated mixture at every budget (ARI 0.83 against 0.99) and ends below Lloyd at convergence on
  only 43 of 80 seed-cells at 25 epochs. It also wins one cell by 4.8% (landsat k=6, ARI 0.53 against 0.38): at
  these sizes the optimum the seeding lands in decides more than the optimizer does. Revival does not help; it
  fired in 3 of 16 cells.

Reading: `Auto` stays as the batch and step default. The default `Init` is the open choice: it is still the
paper's `RandomPoints`, and at these sizes k-means++ seeding costs one pass and is the start from which e-means
never ended above Lloyd.

## Not done

- Assignment is exact, single-threaded and scalar. No approximate assignment (the paper's IVF), no balanced variant
  (present in the paper's source only as commented-out text).
- No Mapper adapter and no REPL surface.
- Soft assignments (mixtures) and anything SPC-side are discussion only.
