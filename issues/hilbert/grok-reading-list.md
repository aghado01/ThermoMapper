**Curated arXiv reading list** (versioned IDs preferred). Grouped by theme and ordered roughly from foundational → more specialized / recent. These are the papers most relevant to the ThermoMapper threads we’ve been discussing (magnetic/sheaf Laplacians, product structures, efficient spectral methods, persistent sheaves, and the continuum complex-analysis side).

### Persistent Sheaves & Sheaf / Persistent Laplacians

- **1808.01513v2** — Hansen & Ghrist, _Toward a Spectral Theory of Cellular Sheaves_
  Foundational paper introducing the sheaf Laplacian and spectral sheaf theory.
- **2112.10906v4** — Wei & Wei, _Persistent sheaf Laplacians_
  Core reference for persistent sheaf Laplacians on point clouds / cellular sheaves.
- **2012.02808v3** — Mémoli, Wan & Wang, _Persistent Laplacians: properties, algorithms and implications_
  Best single reference for the theory, algorithms, Schur-complement view, and stability of ordinary persistent Laplacians.
- **2602.14846v1** — Wang & Wei, _Multi-dimensional Persistent Sheaf Laplacians for Image Analysis_
  Recent multi-dimensional / multi-parameter extension.
- **2509.20220** (check latest version) — Wolf, Fan & Monod, _Generalized Persistent Laplacians and their Spectral Properties_
  Unifying operator-theoretic framework that covers discrete and continuum cases.

### Magnetic Laplacians (discrete & continuum)

- **2102.11391** — Zhang et al., _MagNet: A Magnetic Laplacian based Graph Neural Network_ (or the closely related MagNet papers)
  Clean modern treatment of the magnetic (complex Hermitian) Laplacian for directed graphs.
- **2205.13459** — related Hermitian / Sign-Magnetic Laplacian constructions.
- **2005.08080** — good discrete magnetic weighted Laplacian development and spectral preorder material.
- Continuum side (complex analysis / magnetic Schrödinger operators): look at Aharonov–Bohm and magnetic Laplacian papers such as **2306.08910** and related works on magnetic Laplacians on domains or Riemann surfaces.

### Product Complexes, Product Manifolds & Künneth-type results

- **1910.05656** — Gakhar & Perea, _Künneth Formulae in Persistent Homology_
  Persistent Künneth theorems for categorical and tensor products of filtered spaces.
- **2303.10495** — signal processing / Hodge Laplacian on product spaces of simplicial and cellular complexes (Cartesian products).
- **2505.18844v3** — Park & You, _Geometric medians on product manifolds_
  Directly relevant to your geometric-median / Karcher-scatter implementation on products (shows the mean separates while the median couples).

### Graph Laplacians + Efficient Computation (Chebyshev, Krylov, matrix-free, wavelets)

- **0912.3848** (Hammond–Vandergheynst–Gribonval) — _Wavelets on graphs via spectral graph theory_
  The classic SGWT paper; introduces the Chebyshev approximation that avoids full eigendecomposition.
- **1509.04537** — Lanczos-based spectral graph filtering (comparison with Chebyshev).
- Papers on LOBPCG / matrix-free eigensolvers for graph Laplacians (Knyazev line and follow-ups) and block-Krylov / Chebyshev methods for matrix functions of the Laplacian (heat kernels, filters, wavelets).
- Additional useful computational papers: works on polynomial / rational approximation of matrix functions of sparse Laplacians and on fast graph Fourier transforms for structured graphs.

### Quick starter subset (highest leverage for your current docket)

1. 1808.01513v2 (sheaf Laplacian foundations)
2. 2112.10906v4 (persistent sheaf Laplacians)
3. 2012.02808v3 (persistent Laplacians theory + algorithms)
4. 0912.3848 (spectral graph wavelets + Chebyshev)
5. 2102.11391 (magnetic Laplacian)
6. 1910.05656 (persistent Künneth / products)
7. 2505.18844v3 (geometric medians on products)
8. 2602.14846v1 or 2509.20220 (recent multi-dimensional / generalized persistent sheaf / Laplacian work)

These give a solid spine connecting complex phases (magnetic), sheaf-theoretic and persistent spectral methods, product structures (both continuous and discrete), and the matrix-free / polynomial approximation techniques needed to make the expensive operators practical.

If you want a tighter shortlist focused on one sub-area (e.g., only computational methods, or only the sheaf + magnetic intersection), or PDFs / notes on any of these, just say.

**Yes.** Here is an extension of the reading list that bridges the spectral / sheaf / magnetic / product side with your **statistical–geometry** machinery (Karcher / Fréchet means, geometric medians, scatter on product manifolds, intrinsic estimators) and with the **graph-compiler** side (turning data into graphs / complexes, learning or optimizing the graph/Laplacian structure itself).

### Statistical Geometry ↔ Persistence / Spectra / Sheaves

- **1605.08912** — _A Riemannian Framework for Statistical Analysis of Topological Persistence Diagrams_
  Models persistence diagrams as densities on a Hilbert sphere; gives closed-form Fréchet-type means, geodesics, and classification. Directly links your manifold statistics toolkit to TDA summaries.
- **2410.17214** — _Fréchet Means in Infinite Dimensions_
  Modern asymptotic theory for Fréchet means in general metric spaces (including spaces of persistence diagrams and shape spaces). Useful for theoretical backing of your Karcher estimators when the objects are themselves diagrams or spectra.
- **2204.09155v3** — _Approximating Persistent Homology for Large Datasets_ (Cao & Monod)
  Multiple-subsampling + Fréchet means / mean persistence measures for large data. Practical statistical geometry on the output of (sheaf) persistent homology.
- **2207.03943** (and related) — work on uniqueness and geometry of Fréchet means of persistence diagrams / flat groupings.

### Product Manifolds + Robust Location (your existing implementation)

- **2505.18844v3** — Park & You, _Geometric medians on product manifolds_
  Already highlighted; the clean theoretical companion to your product-manifold geometric-median / Weiszfeld code (mean separates, median couples; robustness and algorithms).

### Graph Learning / Topology Inference / “Compiler” side

These sit at the interface between statistical estimation and the construction of the graph (or sheaf / connection Laplacian) that your compiler would emit:

- **2510.11245** — _Learning the structure of Connection Graphs_
  Learns consistent connection graphs (sheaf-like / orthogonal-valued edge maps) via spectral constraints + Riemannian optimization. Very close in spirit to learning a sheaf Laplacian or magnetic-type structure from data.
- **2506.10628** — _Leveraging Low-rank Factorizations of Conditional Correlation Matrices in Graph Learning_
  Riemannian optimization on structured positive-semidefinite matrices for Laplacian-constrained graph learning. Practical for large-scale “compile a graph from observations.”
- Papers on **Laplacian-constrained Gaussian graphical models** and **graph learning under Laplacian constraints** (many recent works; the Riemannian / low-rank factorizations are the computationally interesting ones).
- **2402.08105** — _Learning Cartesian Product Graphs with Laplacian Constraints_
  Explicitly learns product-graph structure under Laplacian constraints — ties product complexes / product manifolds to the graph-construction problem.

### Riemannian / Manifold-valued Graph Signal Processing

- **2401.14381** — Graph neural network layers for features living on Riemannian manifolds; diffusion driven by a manifold-valued graph Laplacian that pushes features toward weighted Fréchet means of neighbors.
  Conceptual bridge between your intrinsic geometric estimators and discrete diffusion / heat-kernel methods on graphs.
- Related continuum-limit and manifold-learning papers that study when a graph Laplacian converges to a Laplace–Beltrami operator under non-Euclidean metrics (e.g. **2503.16187** and follow-ups).

### Spectral Sparsification & Efficient Graph Construction

- Recent deterministic / almost-linear-time spectral sparsification results (e.g. **2608.13910** and the Spielman–Srivastava lineage).
  Useful if the graph compiler needs to produce sparse, spectrally faithful graphs at scale before feeding them to LOBPCG / Chebyshev / sheaf-Laplacian pipelines.

### Short prioritized add-ons for your docket

1. **2505.18844v3** (product geometric medians — already yours)
2. **1605.08912** (Riemannian stats on persistence diagrams)
3. **2510.11245** (learning connection / sheaf-like graphs)
4. **2402.08105** (learning Cartesian product graphs)
5. **2401.14381** (manifold-valued graph diffusion ↔ Fréchet means)
6. **2204.09155v3** or **2410.17214** (statistical geometry of persistent homology summaries)

These papers let the statistical–geometry layer (means, medians, scatter on products) talk directly to the objects your spectral / sheaf / magnetic engines produce, and to the graph-construction (“compiler”) stage that turns raw data into the operators those engines act on.

If “graph compiler” refers to a more specific component (code generation, scheduling of sparse kernels, IR for graph algorithms, etc.), point me at the relevant files or name and I can refine the list further.
