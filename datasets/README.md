# Datasets

The centralized canonical-dataset store. One home for the real-world benchmark
assets, the reference implementations they're compared against, and the prep
scripts that produced them — reachable the same way from tests, an interactive
REPL, and demo staging.

## How to reach them

- **C# (tests / REPL):** `UserRepl.Datasets.Path("iris.csv")` /
  `Datasets.Path("reference/spc_n/data")`. Resolves this directory by walking up
  to the `datasets/` folder holding the `iris.csv` anchor — no bin-depth-pinned
  `../../../../../` paths.
- **CLI:** `spc --dataset-file datasets/iris.csv` (the CLI takes a path; it does
  not yet resolve canonical names — that's a candidate follow-up).
- **Python adapters:** `notebook/mvp/datasets.py` (loaders for non-CSV formats →
  CSV the CLI can ingest).

## Real-world benchmarks

- `iris.csv` — the standard 150×4 Fisher Iris (four unnormalized numeric features:
  sepal/petal length & width + species label). UCI ML Repository.

- `landsat.csv` — Statlog (Landsat Satellite), UCI id 146 (`sat.trn` + `sat.tst`
  concatenated, space→comma). 6435 rows × 36 raw integer features (4 spectral
  bands × the 3×3 pixel neighbourhood) + a class label last. Classes
  `{1,2,3,4,5,7}` (class 6, "mixture", absent in the source); sizes
  1533/703/1358/626/707/1508. The BWD1996 §G peripheral-capture oracle (cluster
  density decreasing toward the perimeter). Full set — used as-is, the
  mean-edge-distance bandwidth handles the raw scale.

- `isolet.csv.gz` — ISOLET (Isolated Letter Speech Recognition), UCI id 54
  (`isolet1+2+3+4.data` + `isolet5.data`, the `.Z` members `uncompress`-ed), full
  7797 rows × 617 continuous features in [-1,1] + a letter label `1..26` last.
  Gzipped (~10 MB vs 35 MB raw). **Deferred — no active test.** The Domany1999
  §3.2 letter-hierarchy oracle (Fig 3) needs dimensionality reduction, but a PCA
  front-end was built and characterized (2026-06-15) and it does **not** clear the
  wall: unsupervised PCA tops out at **~14/26 letters** (raw 617-dim ~12), ~19%
  covered at ~0.78 purity — the acoustic confusables ({B,D}, {M,N}, the E-set)
  stay fused because PCA keeps high-*variance*, not high-*discriminative*, axes
  (the curse of dimensionality CRQ2018 sidesteps with *supervised* feature
  selection — which would violate validation independence). Reaching the published
  ~20+/3%-unclassified needs a supervised/discriminative or learned-metric
  front-end. Full experiment + data table + unblock paths:
  [`.discussion/issues/spc-parity/isolet-pca-wall.md`](../.discussion/issues/spc-parity/isolet-pca-wall.md).
  Kept prepared for when that front-end lands.

`prep/fetch_uci.py` is the provenance + regenerator for the two UCI sets
(download, decompress the `.Z` members, build the CSVs).

## Reference implementations (`reference/`)

Third-party impls whose own outputs serve as a "second ground truth" alongside
true labels. Adapted in Python (`notebook/mvp/datasets.py`) into CSV the CLI can
ingest — the CLI itself only knows synthetic-generator names and CSV feature
matrices.

- `reference/spc_n/` — Fernando Chaure's C/Haskell port of Domany's original SPC
  code. A 63-point toy that ships with both the input distances (`data` — sparse
  k-NN `#NAME/#DISTANCES` headers + a dense `(mclmatrix …)`) AND the clustering
  Domany's algorithm produces (`model`: per-(point, temperature) assignments +
  the susceptibility curve). Loader: `datasets.py::load_spc_n_example()`.
  (`embedded_features.csv` is its MDS-embedded output, not an input.)

**Adding a reference dataset:** drop the source under `reference/<project>/`, add
a parser to `notebook/mvp/datasets.py` returning `{features|distance_matrix,
labels?, source_path}`, embed distances via `embed_via_mds(...)` if needed, and
emit CSV with `save_features_csv(...)`.

## Synthetic canonical sets (code, not files)

Generated on the fly in `src/synthetic/` — no file assets, indexed here for
discoverability: `euclidean/` (Bwd1995Toy, BlattHierarchy, BlattThreeCluster,
EyeTorusToy, TwoMoons, SpatialBlobs, CrescentEllipsoid, …) and `manifolds/`
(HyperbolicEyeTorus, GaussianManifold, Simplex, …). Reached via
`--dataset <Name>` on the CLI.
