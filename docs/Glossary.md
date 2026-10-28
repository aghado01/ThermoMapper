# Glossary

Terms as ThermoMapper's code uses them, in alphabetical order. The last column names the doc that treats each term in depth. Identifiers are given exactly as they appear in the code.

## A–C

| Term | Definition | Home |
|---|---|---|
| **Abstention** | A resolver returning `Assignment.Unassigned` (−1) for a point instead of forcing a label; generalizes HDBSCAN noise | [ClusterResolution](./ClusterResolution.md) |
| **`AccumulationSpec`** | Declarative choice of which optional statistics a Swendsen–Wang run collects (`None`, `Currencies`, or individual flags) | [SpcEngine](./SpcEngine.md) |
| **Accumulator** | The sufficient statistics of one Swendsen–Wang run at one temperature and replica (moments, optional per-edge counts, per-node sums), plus spins and RNG state for resuming | [SpcEngine](./SpcEngine.md) |
| **Affinities** | The universal per-edge currency: bond activity or edge survival $G_{ij}\in[0,1]$ at one temperature, produced by every sampler | [Concepts](./Concepts.md#per-edge-currencies) |
| **Alexander duality** | Correspondence between $H_{p-1}$ of a complex embedded in $\mathbb R^p$ and $\tilde H_0$ of its dual graph, used for codimension-one zigzags | [ZigzagPersistence](./ZigzagPersistence.md) |
| **Alignments** | Per-edge spin agreement $\langle\delta(\sigma_i,\sigma_j)\rangle$; produced by Swendsen–Wang only | [Concepts](./Concepts.md#per-edge-currencies) |
| **Assignment** | The crisp clustering result: labels dense in `[0, Count)` or −1, with `Coverage` reported separately | [Concepts](./Concepts.md#resolution-carriers) |
| **Axis-alignment law** | A landscape walk runs only when `Dendrogram.CostAxis` equals `Landscape.Axis` (ordinal string comparison) | [ClusterResolution](./ClusterResolution.md) |
| **Bandwidth (δ)** | Kernel length scale turning distances into couplings; δ = 0 on a descriptor means "estimate automatically" | [GraphConstruction](./GraphConstruction.md#bandwidth-estimation) |
| **BAPS** | Bayesian adaptive penalized splines: coefficients integrated out (marginal-evidence mode, REML) or drawn from Gaussian conditionals (Gibbs modes, hierarchical fits); smoothing strength λ sampled or selected by REML | [BayesianRegression](./BayesianRegression.md) |
| **Bar / Barcode** | A persistence interval (birth, death, dimension, optional payloads, end types) / the list of bars with an axis label; the output of every persistence engine | [Concepts](./Concepts.md#barcodes) |
| **Barcode forest** | Level-annotated forest used by the graph zigzag H0 engines to date merges, splits and departures | [ZigzagPersistence](./ZigzagPersistence.md) |
| **BARS** | Bayesian adaptive regression splines: reversible-jump MCMC over knot number and positions, with coefficients marginalized | [BayesianRegression](./BayesianRegression.md) |
| **BGP** | Bayesian Gaussian-process regression that samples a kernel bandwidth under an empirical-Bayes prior | [BayesianRegression](./BayesianRegression.md) |
| **Blatt hierarchy** | SPC hierarchical partitioning after Blatt–Wiseman–Domany, producing a `PartitionHierarchy` across temperatures | [SpcAnalysis](./SpcAnalysis.md) |
| **Bottom-k** | The k smallest eigenpairs of a symmetric matrix or Laplacian, returned in ascending order | [LinearAlgebra](./LinearAlgebra.md) |
| **Canonical slot** | The CSR slot of an undirected edge in the smaller endpoint's row (`Targets[e] > i`); the only slot whose per-edge values are meaningful | [Concepts](./Concepts.md#the-csr-graph) |
| **Changepoint (exact)** | Piecewise-constant regression computed exactly by dynamic programming | [BayesianRegression](./BayesianRegression.md) |
| **Chosen temperature ($T_{\text{clus}}$)** | Midpoint of the superparamagnetic plateau (the χ peak when no high-T cliff is found), at which the sweep mints its currencies | [SpcEngine](./SpcEngine.md#sweeps) |
| **Clearing** | Persistence optimization that skips reducing columns already known to be paired (H0 via union–find, then pivots of the dimension below) | [PersistentHomology](./PersistentHomology.md) |
| **CoMembership** | Per-edge same-cluster frequency $\langle n_{ij}\rangle$ over Swendsen–Wang draws; the default SPC cut field | [Concepts](./Concepts.md#per-edge-currencies) |
| **Condensed tree** | HDBSCAN-style simplification of a merge tree that keeps clusters of at least `minClusterSize` and measures their stability in λ | [ClusterResolution](./ClusterResolution.md) |
| **Coning** | Turning a zigzag into one standard filtration by adding an apex ω and a cone for every deletion (`FastZigzag`) | [ZigzagPersistence](./ZigzagPersistence.md) |
| **Core distance** | Distance from a point to its `minPts`-th nearest *other* point | [Hdbscan](./Hdbscan.md) |
| **Cost axis** | Free-text `Dendrogram.CostAxis` naming what merge heights mean, e.g. `"temperature"`, `"neg_coupling"` | [Concepts](./Concepts.md#resolution-carriers) |
| **Coupling (J)** | Similarity edge weight (larger = closer) produced by a kernel from a distance; the Potts interaction strength | [Concepts](./Concepts.md#edge-weights-distance-coupling-affinity) |
| **Coupling fidelity** | How curvature enters the kernel: `GeodesicLinear` (flat kernel on geodesic distance) or `Intrinsic` (manifold heat kernel) | [GraphConstruction](./GraphConstruction.md) |
| **Cover** | Overlapping bins over a Mapper lens's range: `Cover.Uniform`, `Cover.Balanced` (quantile) or `ProductCover` | [Mapper](./Mapper.md) |
| **CSR graph** | `CsrGraph`: symmetric compressed-sparse-row adjacency storing each undirected edge in both endpoint rows | [Concepts](./Concepts.md#the-csr-graph) |
| **Cumulative field** | Deterministic graphical clusterer that evaluates edge-survival affinities in closed form from a cumulative-energy ladder, without Monte Carlo | [CumulativeField](./CumulativeField.md) |
| **Currency** | A shared data type passed between modules; specifically, a per-edge field at one temperature (Affinities, Alignments, CoMembership) | [Concepts](./Concepts.md) |

## D–H

| Term | Definition | Home |
|---|---|---|
| **Dendrogram** | Index-based binary merge tree: leaves `0..N−1`, merge i is node `N+i`, each with a height and size | [Concepts](./Concepts.md#resolution-carriers) |
| **`DiagnosticsLog`** | Info/Warning messages recording every automatic decision of a graph build, returned with the result | [GraphConstruction](./GraphConstruction.md) |
| **Diagram distance** | Wasserstein-type distance between barcodes: exact (Hungarian), sliced, or entropic (Sinkhorn) | [PersistentHomology](./PersistentHomology.md) |
| **Dominated vertex** | A vertex whose maximal simplices share another vertex; removing it is a strong collapse | [ZigzagPersistence](./ZigzagPersistence.md) |
| **Dynamic MSF** | Minimum spanning forest maintained under insertions and deletions (Euler-tour tree + link-cut tree), used by the fast graph H1 engine | [ZigzagPersistence](./ZigzagPersistence.md) |
| **Effective resistance** | Resistance distance between nodes computed from Laplacian eigenpairs; an alternative filtration edge value | [SpectralGraphs](./SpectralGraphs.md) |
| **Elder rule** | When components merge, the older one survives and the younger one's bar ends | [PersistentHomology](./PersistentHomology.md) |
| **Ensemble shell** | `ChainEnsemble`: runs C independent chains with derived seeds, burn-in, sampling rounds and R̂ checks around any `IChain` | [Samplers](./Samplers.md) |
| **EOM (excess of mass)** | Cluster selection maximizing summed stability over a tree, used by HDBSCAN condensation and by landscape walks | [ClusterResolution](./ClusterResolution.md) |
| **Essential bar** | A bar with `Death = +∞`, i.e. a class that never dies in the filtration | [Concepts](./Concepts.md#barcodes) |
| **Euler-tour tree** | Tree stored as its Euler tour in an implicit treap, supporting link, cut and connectivity in $O(\log n)$ expected time | [ZigzagPersistence](./ZigzagPersistence.md) |
| **Eye skeleton** | Coordinate-free description of the "eye" dataset (stroke arcs carrying tubes, optional pupil), realized in flat or hyperbolic space | [SyntheticDatasets](./SyntheticDatasets.md) |
| **Fiedler vector** | Eigenvector of the second-smallest Laplacian eigenvalue; also a Mapper graph lens | [SpectralGraphs](./SpectralGraphs.md) |
| **Filtration** | Total order of simplices, faces first, with non-decreasing values; positions in it are simplex indices | [PersistentHomology](./PersistentHomology.md) |
| **FK bond probability** | Fortuin–Kasteleyn $p_e=1-e^{-J_e/T}$: the chance an aligned edge freezes in a Swendsen–Wang draw | [SpcEngine](./SpcEngine.md) |
| **Flag (graph-Rips) filtration** | Filtration of a graph's clique complex: edges born at their distance, triangles at their longest edge | [PersistentHomology](./PersistentHomology.md) |
| **Frame** | SPC: one `Accumulator` (one temperature and replica). Mapper: one nerve at one sweep parameter | [SpcAnalysis](./SpcAnalysis.md), [Mapper](./Mapper.md) |
| **Friends-of-friends cut** | Threshold a per-edge field at θ and take connected components | [SpcAnalysis](./SpcAnalysis.md) |
| **Functional** | Fixed-width scalar reduction of each MCMC draw (e.g. knot count) on which convergence diagnostics are computed | [Samplers](./Samplers.md) |
| **Generalized rank** | Rank of the limit-to-colimit map over a window of a zigzag; Möbius inversion of these ranks gives the barcode | [ZigzagPersistence](./ZigzagPersistence.md) |
| **`GraphCompiler`** | Static five-stage graph-construction engine: topology, filter, repair, refinement, scaling/projection | [GraphConstruction](./GraphConstruction.md) |
| **Grassmann median** | Geometric median of subspaces on the Grassmannian, used to aggregate distributed SPRED blocks | [Spred](./Spred.md) |
| **Groups** | Soft clustering result: an N×K membership field, resolved by `Argmax()` | [Concepts](./Concepts.md#resolution-carriers) |
| **HDT** | Holm–de Lichtenberg–Thorup dynamic connectivity: levelled spanning forests with $O(\log^2 n)$ amortized updates | [ZigzagPersistence](./ZigzagPersistence.md) |
| **Heat-kernel calibration** | Choosing an intrinsic kernel's bandwidth on a curved space by matching the second moment of the k-th-neighbor sample | [RiemannianGeometry](./RiemannianGeometry.md) |
| **Hubness** | Skew of kNN in-degree in high dimensions; graph construction's first diagnostic gate | [GraphConstruction](./GraphConstruction.md) |

## I–O

| Term | Definition | Home |
|---|---|---|
| **Incarnation** | One continuous presence of a cell in a zigzag; a re-added id starts a new incarnation | [ZigzagPersistence](./ZigzagPersistence.md) |
| **`IntervalEnd`** | `Closed` or `Open` endpoint of a bar; set by zigzag engines from the creating and killing arrows | [Concepts](./Concepts.md#barcodes) |
| **Involuted persistence** | Čufar–Virk method: cohomology finds pairs, a restricted homology reduction supplies cycle representatives | [PersistentHomology](./PersistentHomology.md) |
| **IRLS** | Iteratively reweighted least squares for robust location and scatter on a Riemannian manifold (`Irls.Solve`) | [RiemannianGeometry](./RiemannianGeometry.md) |
| **k-th sample** | Per-node distances to the k-th nearest neighbor, carried unchanged through graph construction for bandwidth estimation (unless mutual-kNN filtering uses `MutualBandwidthSource.MutualKth`) | [GraphConstruction](./GraphConstruction.md) |
| **λ (persistence)** | Inverse merge height, $\lambda=1/h$, used by condensation to measure cluster stability | [ClusterResolution](./ClusterResolution.md) |
| **Landscape** | Per-node scalar field sampled on an ascending grid over a named axis (e.g. mean cluster size vs temperature) | [Concepts](./Concepts.md#node-landscapes) |
| **Landscape walk** | Cluster selection that integrates a landscape along a dendrogram whose cost axis matches | [ClusterResolution](./ClusterResolution.md) |
| **Lens (filter)** | Map from points or nodes to $\mathbb R$ or $\mathbb R^D$ that Mapper covers | [Mapper](./Mapper.md) |
| **Level interval** | Zigzag interval $[b,d]$ of complexes $K_b..K_d$ in which a class is alive; reported as `Bar(b−1, d)` | [ZigzagPersistence](./ZigzagPersistence.md) |
| **Link-cut tree** | Splay-based dynamic tree supporting path-maximum queries in $O(\log n)$ amortized time | [ZigzagPersistence](./ZigzagPersistence.md) |
| **LMP** | Local mutual proximity: replaces each coupling by a product of local ranks at its endpoints, rescuing edges a global bandwidth would crush | [GraphConstruction](./GraphConstruction.md) |
| **LOBPCG** | Locally optimal block preconditioned conjugate gradient: the iterative bottom-k eigensolver | [LinearAlgebra](./LinearAlgebra.md) |
| **MAD consistency factor** | Constant c in the default bandwidth $\delta=c\cdot\mathrm{MAD}$, depending on kernel and geometry | [GraphConstruction](./GraphConstruction.md#bandwidth-estimation) |
| **Metric registry** | `MetricRegistry`: parses metric spec strings and bridges them to struct metric types via a visitor | [DistanceMetrics](./DistanceMetrics.md) |
| **Mirror slot** | The reverse half-edge of a CSR slot | [Concepts](./Concepts.md#the-csr-graph) |
| **Möbius addition** | The Poincaré-ball group operation $c\oplus y$ used for translations and exponential maps | [SyntheticDatasets](./SyntheticDatasets.md) |
| **MoMPCA** | Median-of-means PCA: robust principal subspaces from block estimates | [RiemannianGeometry](./RiemannianGeometry.md) |
| **Mutual reachability distance** | $\max\{\kappa(u),\kappa(v),d(u,v)\}$ with κ the core distance; HDBSCAN's MST metric | [Hdbscan](./Hdbscan.md) |
| **`NeighborSelection`** | Directed kNN result (neighbors, nearest and k-th distances) from graph construction's first stage | [Concepts](./Concepts.md#neighbor-selections-and-eigenpairs) |
| **Nerve** | Mapper's output graph: nodes are (bin, cluster) pairs, edges join nodes sharing points, weighted by overlap count | [Mapper](./Mapper.md) |
| **Node matching** | `NerveDiff`'s best-overlap mapping from nodes of one Mapper frame to the next | [Mapper](./Mapper.md) |

## P–S

| Term | Definition | Home |
|---|---|---|
| **Participation ratio** | $(\sum\lambda_i)^2/\sum\lambda_i^2$ of a node's PCA eigenvalues; Mapper's intrinsic-dimension estimate | [Mapper](./Mapper.md) |
| **`PartitionHierarchy` (T-stack)** | One partition per temperature, cold to hot; nests when each hotter cluster lies inside one colder cluster | [SpcAnalysis](./SpcAnalysis.md) |
| **Pathology interrupt** | Configured threshold that makes a graph build throw `GraphPathologyException` | [GraphConstruction](./GraphConstruction.md) |
| **Periphery completion** | Policies that fill unassigned points without relabelling assigned ones | [ClusterResolution](./ClusterResolution.md) |
| **Persistence** | Bar length `Death − Birth` | [PersistentHomology](./PersistentHomology.md) |
| **Poincaré ball** | Unit-ball model of hyperbolic space (curvature −1) used by hyperbolic metrics, kernels and datasets | [RiemannianGeometry](./RiemannianGeometry.md) |
| **Potts model** | q-state spin model on a coupling graph whose superparamagnetic phase reveals clusters | [SpcEngine](./SpcEngine.md) |
| **Protect set** | Edges (typically H1 loop edges) that local mutual proximity must leave untouched; supplied by a `ProtectedEdgeSource` | [PersistentHomology](./PersistentHomology.md#cycle-representatives-and-the-h1-protect-set) |
| **R̂ (R-hat)** | Between/within-chain convergence diagnostic computed on functionals | [Samplers](./Samplers.md) |
| **Replica** | Independent Swendsen–Wang run at the same temperature with its own seed | [SpcEngine](./SpcEngine.md) |
| **Reversible jump (RJ-MCMC)** | MCMC that moves between models of different dimension (e.g. knot counts) | [Samplers](./Samplers.md) |
| **Rich sweep** | Sweep whose frames carry per-edge and per-node arrays, required by the hierarchical resolvers | [SpcAnalysis](./SpcAnalysis.md) |
| **`SeedTree`** | Deterministic derivation of child seeds by index from one master seed (SplitMix64) | [ProbabilityAndInformation](./ProbabilityAndInformation.md) |
| **Single linkage** | Merge tree from a maximum spanning tree of couplings; cost axis `"neg_coupling"` | [ClusterResolution](./ClusterResolution.md) |
| **Slot** | Index into a CSR graph's `Targets`/`Weights`: one directed half-edge | [Concepts](./Concepts.md#the-csr-graph) |
| **SP plateau** | Superparamagnetic temperature band between the FM→SP and SP→PM transitions | [SpcAnalysis](./SpcAnalysis.md) |
| **SPC** | Superparamagnetic clustering: clusters read from Potts-model correlations in the superparamagnetic phase | [SpcEngine](./SpcEngine.md) |
| **SPCX** | Versioned binary checkpoint format for Swendsen–Wang accumulators, with presence bytes for optional arrays | [SpcEngine](./SpcEngine.md) |
| **Spectral noise field** | Random field with power-law spectrum (e.g. pink noise), used for textured backgrounds | [Samplers](./Samplers.md) |
| **SPRED** | Persistence-preserving linear dimensionality reduction by subspace annealing | [Spred](./Spred.md) |
| **Strong collapse** | Repeated removal of dominated vertices, preserving homotopy type; yields a complex's core | [ZigzagPersistence](./ZigzagPersistence.md) |
| **Subspace annealing** | Simulated annealing over orthonormal k-frames (`SubspaceAnnealer`) | [RiemannianGeometry](./RiemannianGeometry.md) |
| **Superparamagnetic phase** | Intermediate-temperature regime where dense regions align internally but not with each other | [SpcEngine](./SpcEngine.md) |
| **Survivor** | Zigzag class alive at the end of the sequence, reported with finite `Death = m` and an open end | [ZigzagPersistence](./ZigzagPersistence.md) |
| **Susceptibility (χ)** | Fluctuation observable of a sweep profile whose peaks mark transitions (FK, reduced, or magnetic variants) | [SpcAnalysis](./SpcAnalysis.md) |
| **Swendsen–Wang** | Cluster Monte Carlo update: freeze bonds on aligned edges, find clusters, recolour each cluster | [SpcEngine](./SpcEngine.md) |

## T–Z

| Term | Definition | Home |
|---|---|---|
| **Tempering** | Parallel tempering: a ladder of chains at inverse temperatures β that swap states to cross barriers | [Samplers](./Samplers.md) |
| **Thermal dendrogram** | Tree built from per-edge decoupling temperatures; cost axis `"temperature"`, heights non-increasing | [ClusterResolution](./ClusterResolution.md) |
| **Thermal EOM** | SPC resolver that walks a thermal dendrogram built from per-edge δ̄(T) curves | [SpcAnalysis](./SpcAnalysis.md) |
| **Track 1 / Track 2** | SPC hierarchical resolvers: excess of mass over a nested T-stack (`HierarchyEom`) / lineage persistence (`LineagePersistence`) | [SpcAnalysis](./SpcAnalysis.md) |
| **Tracking barcode** | H0 barcode from a Mapper parameter sweep; a component-tracking summary, not a filtration's persistence module | [Mapper](./Mapper.md) |
| **Wasserstein distance ($W_p$)** | Optimal-matching distance between persistence diagrams, with unmatched points sent to the diagonal | [PersistentHomology](./PersistentHomology.md) |
| **xoshiro256++** | The library's pseudo-random generator (`Xoshiro256PlusPlus`), seeded via SplitMix64 and checkpointable | [ProbabilityAndInformation](./ProbabilityAndInformation.md) |
| **Zigzag filtration** | Sequence of single-cell additions and deletions (`ZigzagFiltration`), with faces referenced by global cell id | [ZigzagPersistence](./ZigzagPersistence.md) |

**See also:** [Concepts](./Concepts.md) · [Architecture](./Architecture.md) · [DesignPrinciples](./DesignPrinciples.md)
