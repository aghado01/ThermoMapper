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

## Not done

- No measurement of the laws below their fitting grid (`n < 10⁵`, `k < 256`), which is where ThermoMapper's
  datasets sit. The paper says the laws are not optimal on small, low-dimensional data. Whether `Auto` deserves to
  be the default there is open until measured against Lloyd from k-means++.
- Assignment is exact, single-threaded and scalar. No approximate assignment (the paper's IVF), no balanced variant
  (present in the paper's source only as commented-out text).
- No Mapper adapter and no REPL surface.
- Soft assignments (mixtures) and anything SPC-side are discussion only.
