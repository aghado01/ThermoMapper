**Yes, but only in structured cases or via approximations.** Exact FFTs diagonalize graph (and certain sheaf) Laplacians precisely when the underlying complex has the right symmetry; for the general proximity graphs, Rips filtrations, or arbitrary cellular sheaves that appear in ThermoMapper-style pipelines, they do not. There are useful accelerations and work-arounds that fit the project’s existing spectral and complex-analytic machinery.

### Exact case: circulant / lattice / product structure

When the graph (or the 1-skeleton of the complex) is a cycle, a rectangular lattice with periodic boundaries, a torus grid, or more generally a Cayley graph of an abelian group (or Cartesian product of such), the combinatorial Laplacian \(L = D - A\) is a circulant or multilevel block-circulant matrix. These are diagonalized by the discrete Fourier transform:

\[
L = F^\* \Lambda F
\]

(with the usual unitary DFT matrix \(F\)). Eigenvalues are closed-form (e.g., for the cycle \(C_n\): \(2 - 2\cos(2\pi k/n)\); for a 2-D grid product the sums of the 1-D spectra). Matvecs, eigendecompositions, heat-kernel evaluations, and powers of \(L\) all become \(O(N\log N)\) via FFT.

Your codebase already has the primitive: `Maths.LinAlg.Fft` implements radix-2 Cooley–Tukey (1-D and separable N-D) over power-of-two lengths, used for spectral noise fields and convolution. Extending it to the Laplacian spectrum on a regular grid is therefore immediate and exact.

The same holds for the magnetic Laplacian you already implement (`MagneticLaplacianOperator`). On a lattice with a constant magnetic flux (or a translation-invariant connection), the magnetic Laplacian remains diagonalized by a Fourier (or magnetic Bloch) basis; the phases \(e^{i 2\pi q\Theta}\) simply shift the eigenvalues. This sits cleanly inside the “complex-analytic TDA thread” already flagged in the magnetic operator remarks.

### Sheaf and persistent-sheaf setting

A cellular sheaf Laplacian \(L*{\mathcal{F}} = \delta^\*\delta\) (degree-0) reduces to the ordinary graph Laplacian when stalks are \(\mathbb{R}\) and restrictions are identities. When stalks are higher-dimensional or restrictions are non-trivial linear maps, \(L*{\mathcal{F}}\) is a block matrix whose blocks are built from the restriction maps. On a regular lattice (or product complex) with _translation-invariant_ restriction maps, the whole operator remains circulant / block-circulant and is still FFT-diagonalizable. The sheaf Fourier transform is then just the ordinary Fourier transform acting fiberwise on the stalks.

Persistent sheaf Laplacians (Wei–Wei and subsequent multi-dimensional variants) further enlarge the matrix by considering pairs of filtration values \((X_a,X_b)\). If the underlying filtration is itself a product filtration on a grid (or can be approximated by one), the same Fourier structure can be exploited at each scale; otherwise the cost is the same as ordinary dense spectral work on the enlarged complex.

### General graphs and point-cloud data (the expensive regime)

ThermoMapper’s typical workflow—point cloud \(\to\) proximity graph (kNN / \(\varepsilon\)-ball / LMP) \(\to\) dense Laplacian materialization (`GraphLaplacian.BuildDense` / column-major) \(\to\) dense eigen or LOBPCG—produces an unstructured sparse or dense matrix. No exact FFT exists. Construction and storage already scale poorly (\(O(n^2)\) for dense affinities), and full eigendecomposition is cubic.

Practical accelerations that stay compatible with the existing stack:

- **Matrix-free Krylov + fast matvecs.** For dense RBF / Gaussian affinities the kernel matrix (hence the Laplacian) admits fast multipole / nonequispaced FFT (NFFT) matvecs in near-linear time. The 2018 NFFT-meets-Krylov line of work embeds exactly into Lanczos / LOBPCG; you already have LOBPCG and a matrix-free magnetic operator, so the pattern is familiar.
- **Approximate graph Fourier transforms.** Greedy sparse orthogonal factorizations (modified Jacobi / Givens rotations, FAUST-style multilayer sparse transforms) produce an approximate eigenbasis that can be applied in \(O(n^\alpha)\) with \(1<\alpha<2\). These are useful when you only need a few low-frequency modes or a fast approximate spectral embedding.
- **Structure exploitation / multiscale.** If the data live near a manifold that admits a good chart or a hierarchical decomposition (octrees, neural fields, etc.), one can transfer the problem to a regular grid (or a sequence of grids) where FFT is exact, then pull the spectrum back. Heat-semigroup engines already discussed in the Hilbert issues are natural consumers of such spectra.
- **Sparse direct methods or iterative solvers for the linear systems that arise in sheaf cohomology / harmonic projection**, rather than full eigendecomposition, when only the kernel or a few modes are required.

### Project-specific fit

- The magnetic Laplacian already injects the complex structure you want for directed / oriented cycles and Aharonov–Bohm-type discriminants; on any graph that is close to a lattice it can ride the FFT.
- Persistent homology / zigzag / nerve pipelines in `TDA.Ph` and Mapper produce filtrations whose 1-skeleta are still proximity graphs; the same caveats apply, but grid or product filtrations (synthetic eye-torus, hyperbolic hierarchies, etc.) are FFT-friendly.
- Heat-kernel / semigroup methods (Hilbert notes) become \(O(N\log N)\) per evaluation once the spectrum is known via FFT.

**Bottom line.** For any experiment whose underlying space is (or can be faithfully discretized as) a regular lattice, torus, or abelian Cayley graph, use the existing FFT directly on the Laplacian (or magnetic / sheaf Laplacian). For the generic point-cloud / arbitrary sheaf case the answer is “approximate or matrix-free,” not a pure FFT. The cheapest high-leverage next step is probably an NFFT or fast-multipole matvec for the dense kernel Laplacian, feeding the LOBPCG you already have, while keeping exact FFT for the structured synthetic benchmarks and continuum limits that appear in the complex-analytic and heat-semigroup threads.

**Wavelet methods are a strong fit**—often stronger than pure FFT for the irregular proximity graphs and sheaf Laplacians that dominate ThermoMapper workloads—because they deliver multiscale, localized analysis while remaining matrix-free.

### Spectral graph wavelets (the main practical route)

The spectral graph wavelet transform (SGWT; Hammond–Vandergheynst–Gribonval) constructs wavelets by spectral multipliers of the Laplacian:

\[
\psi\_{s,n} = g(sL)\,\delta_n
\]

where \(g\) is a band-pass kernel on the spectrum of \(L\) (or the normalized / magnetic / sheaf Laplacian). The resulting family is localized in both the vertex domain and the “frequency” domain defined by the Laplacian spectrum.

**The key computational win**: you never need the full eigendecomposition. Approximate \(g(sL)\) by a Chebyshev (or other) polynomial of moderate degree \(K\). Applying the operator then reduces to a short sequence of sparse (or block-sparse) matrix–vector products with \(L\). Complexity is essentially \(O(K\cdot|E|)\) per scale / per signal—linear in the number of edges for the sparse kNN / LMP graphs that the project already builds. This is exactly the same matrix-function action primitive already sketched in the heat-semigroup notes (`issues/hilbert/claude-heat-semigroup-engines.md`), where spectral graph wavelets are listed among the natural compositions of the Chebyshev / Krylov / LOBPCG toolkit.

On regular lattices or tori the construction recovers classical discrete wavelets (and can be realized with FFTs), so the structured cases you already handle remain efficient. On general point-cloud graphs the Chebyshev route stays cheap while the pure Fourier basis would be cubic.

### Why this helps persistent sheaves and complex analysis

- **Sheaf Laplacians**. Replace the ordinary graph Laplacian by the degree-0 sheaf Laplacian \(L*{\mathcal{F}}=\delta^\*\delta\). The same spectral-multiplier + Chebyshev pipeline works; the only requirement is a fast matvec with \(L*{\mathcal{F}}\). For cellular sheaves whose restriction maps are local, that matvec remains sparse or block-sparse. Persistent sheaf Laplacians (pairs of filtration values) simply give a family of such operators; you can extract wavelet signatures at each scale of the filtration without ever forming the dense matrices that `GraphLaplacian.BuildDense` currently materializes.
- **Magnetic / complex setting**. The magnetic Laplacian you already implement (`MagneticLaplacianOperator`) is a Hermitian operator. Spectral wavelets go through unchanged (or with complex-valued kernels). This sits inside the same “complex-analytic TDA thread” already flagged for the magnetic operator—localized phase-sensitive multiscale features rather than global Fourier modes.
- **Multiscale by design**. Wavelets are inherently multi-resolution; they mesh naturally with the filtration / persistence viewpoint and with the heat-semigroup / diffusion observables already under discussion in the Hilbert notes (HKS fields, heat content, local dimension, etc.).

### Other wavelet families worth keeping in view

- **Diffusion wavelets** (Coifman–Maggioni). Built by dyadic powers of a diffusion operator (heat kernel or random-walk matrix). Again matrix-free once you can apply powers of the operator; they produce orthonormal multiscale bases and have a clear continuum limit. They sit especially close to the heat-semigroup engines already being analyzed.
- **Lifting-scheme / spatial wavelets**. Constructed directly on the graph by prediction/update steps (no spectrum required). Useful when you want extremely local, adaptive bases or when the graph changes over a filtration.
- **Tight frames / filter banks** on graphs. Polynomial approximations of trigonometric or sinusoidal spectral filters give low-order Chebyshev expansions with rigorous error bounds and perfect reconstruction.

### Practical placement in the existing stack

Everything needed is already present or sketched:

- Sparse (or matrix-free) Laplacian operators → Chebyshev recurrence.
- LOBPCG for the low-frequency end of the spectrum when needed.
- FFT for the structured / grid cases.
- Heat-kernel actions (Chebyshev shared basis across a time grid) that the Hilbert notes already treat as a first-class primitive.

A minimal next step is therefore a thin spectral-wavelet layer on top of the existing `ILinearOperator` + Chebyshev action machinery: choose a kernel \(g\), approximate at a discrete set of scales, and obtain either a multiscale feature map on vertices or a sparse coefficient representation of a signal (or of a sheaf cochain). For persistent sheaves one simply repeats the construction across the filtration, obtaining a multiscale signature that encodes both geometry and the sheaf data.

In short: wavelets give you the multiscale localization that pure Fourier lacks, the linear-cost evaluation that dense eigendecompositions lack, and a clean extension path from ordinary graphs through magnetic Laplacians to persistent sheaves—all while reusing the matrix-function and LOBPCG infrastructure already in the repo.

**Wavelet methods are a strong fit**—often stronger than pure FFT for the irregular proximity graphs and sheaf Laplacians that dominate ThermoMapper workloads—because they deliver multiscale, localized analysis while remaining matrix-free.

### Spectral graph wavelets (the main practical route)

The spectral graph wavelet transform (SGWT; Hammond–Vandergheynst–Gribonval) constructs wavelets by spectral multipliers of the Laplacian:

\[
\psi\_{s,n} = g(sL)\,\delta_n
\]

where \(g\) is a band-pass kernel on the spectrum of \(L\) (or the normalized / magnetic / sheaf Laplacian). The resulting family is localized in both the vertex domain and the “frequency” domain defined by the Laplacian spectrum.

**The key computational win**: you never need the full eigendecomposition. Approximate \(g(sL)\) by a Chebyshev (or other) polynomial of moderate degree \(K\). Applying the operator then reduces to a short sequence of sparse (or block-sparse) matrix–vector products with \(L\). Complexity is essentially \(O(K\cdot|E|)\) per scale / per signal—linear in the number of edges for the sparse kNN / LMP graphs that the project already builds. This is exactly the same matrix-function action primitive already sketched in the heat-semigroup notes (`issues/hilbert/claude-heat-semigroup-engines.md`), where spectral graph wavelets are listed among the natural compositions of the Chebyshev / Krylov / LOBPCG toolkit.

On regular lattices or tori the construction recovers classical discrete wavelets (and can be realized with FFTs), so the structured cases you already handle remain efficient. On general point-cloud graphs the Chebyshev route stays cheap while the pure Fourier basis would be cubic.

### Why this helps persistent sheaves and complex analysis

- **Sheaf Laplacians**. Replace the ordinary graph Laplacian by the degree-0 sheaf Laplacian \(L*{\mathcal{F}}=\delta^\*\delta\). The same spectral-multiplier + Chebyshev pipeline works; the only requirement is a fast matvec with \(L*{\mathcal{F}}\). For cellular sheaves whose restriction maps are local, that matvec remains sparse or block-sparse. Persistent sheaf Laplacians (pairs of filtration values) simply give a family of such operators; you can extract wavelet signatures at each scale of the filtration without ever forming the dense matrices that `GraphLaplacian.BuildDense` currently materializes.
- **Magnetic / complex setting**. The magnetic Laplacian you already implement (`MagneticLaplacianOperator`) is a Hermitian operator. Spectral wavelets go through unchanged (or with complex-valued kernels). This sits inside the same “complex-analytic TDA thread” already flagged for the magnetic operator—localized phase-sensitive multiscale features rather than global Fourier modes.
- **Multiscale by design**. Wavelets are inherently multi-resolution; they mesh naturally with the filtration / persistence viewpoint and with the heat-semigroup / diffusion observables already under discussion in the Hilbert notes (HKS fields, heat content, local dimension, etc.).

### Other wavelet families worth keeping in view

- **Diffusion wavelets** (Coifman–Maggioni). Built by dyadic powers of a diffusion operator (heat kernel or random-walk matrix). Again matrix-free once you can apply powers of the operator; they produce orthonormal multiscale bases and have a clear continuum limit. They sit especially close to the heat-semigroup engines already being analyzed.
- **Lifting-scheme / spatial wavelets**. Constructed directly on the graph by prediction/update steps (no spectrum required). Useful when you want extremely local, adaptive bases or when the graph changes over a filtration.
- **Tight frames / filter banks** on graphs. Polynomial approximations of trigonometric or sinusoidal spectral filters give low-order Chebyshev expansions with rigorous error bounds and perfect reconstruction.

### Practical placement in the existing stack

Everything needed is already present or sketched:

- Sparse (or matrix-free) Laplacian operators → Chebyshev recurrence.
- LOBPCG for the low-frequency end of the spectrum when needed.
- FFT for the structured / grid cases.
- Heat-kernel actions (Chebyshev shared basis across a time grid) that the Hilbert notes already treat as a first-class primitive.

A minimal next step is therefore a thin spectral-wavelet layer on top of the existing `ILinearOperator` + Chebyshev action machinery: choose a kernel \(g\), approximate at a discrete set of scales, and obtain either a multiscale feature map on vertices or a sparse coefficient representation of a signal (or of a sheaf cochain). For persistent sheaves one simply repeats the construction across the filtration, obtaining a multiscale signature that encodes both geometry and the sheaf data.

In short: wavelets give you the multiscale localization that pure Fourier lacks, the linear-cost evaluation that dense eigendecompositions lack, and a clean extension path from ordinary graphs through magnetic Laplacians to persistent sheaves—all while reusing the matrix-function and LOBPCG infrastructure already in the repo.

**Product structures give some of the cleanest algebraic and computational leverage available** for the spectral, heat, wavelet, and sheaf machinery already under discussion in ThermoMapper.

### Product manifolds (continuous spectral geometry)

On a Riemannian product \(M\times N\) with the product metric the Laplace–Beltrami operator separates completely:

\[
\Delta\_{M\times N}=\Delta_M\otimes I+I\otimes\Delta_N.
\]

Eigenvalues are therefore pure sums \(\lambda_i^M+\lambda_j^N\) and eigenfunctions are products of the factor eigenfunctions. The heat kernel factors as an ordinary product of the two heat kernels. The same holds for higher-order operators, for the Hodge Laplacian on forms, and (with the appropriate product structures) for the magnetic Laplacian or other twisted operators.

This is the continuum limit one expects when discrete product graphs converge to product manifolds. It also supplies exact model problems on which one can test numerical approximations (Chebyshev heat/wavelet actions, LOBPCG, etc.) before moving to irregular point-cloud graphs.

### Product graphs and product complexes

For the combinatorial (or normalized) Laplacian the **Cartesian product** \(G\square H\) is the most useful:

\[
L\_{G\square H}=L_G\otimes I+I\otimes L_H
\]

(Kronecker sum). Consequently the spectrum is again the sumset of the factor spectra and the eigenvectors are products. When the factors themselves admit fast transforms (paths, cycles, grids, trees with recursive structure, etc.) the whole product does as well—FFT on each factor, or separable Chebyshev expansions, or tensorized LOBPCG.

Other graph products (tensor/categorical, strong, lexicographic) have more complicated Laplacian formulae and generally lose the clean additive spectrum unless the factors are regular. For simplicial complexes there are analogous product constructions; the resulting chain complexes and Laplacians inherit Künneth-type exact sequences and spectral additivity under suitable conditions.

In the ThermoMapper stack this is immediately actionable: any time the underlying space (or a cover, or a parameter domain) factors, the expensive dense `GraphLaplacian.BuildDense` path and the matrix-free LOBPCG/Chebyshev paths both become dramatically cheaper. The same Kronecker structure turns the heat-semigroup and spectral-wavelet engines already sketched in the Hilbert notes into separable operations.

### Persistent homology and Künneth

The classical Künneth theorem relates the homology of a product to the homology of the factors (with Tor terms). In the filtered setting one obtains Künneth formulae for persistent homology, though the presence of the filtration makes the statements more delicate; multiparameter persistence appears naturally once each factor carries its own filtration parameter.

### Persistent sheaves and product complexes

Cellular sheaves behave well under products. One can form:

- the external product of two sheaves \(\mathcal{F}\boxtimes\mathcal{G}\) on the product complex,
- the pull-back of a sheaf along a projection,
- sheaves defined directly on a product filtration.

Because the restriction maps of a product sheaf are built from the restriction maps of the factors, the resulting sheaf Laplacian again inherits a Kronecker-sum structure when the underlying complex is a Cartesian product. Persistent sheaf Laplacians (or multi-dimensional persistent sheaf Laplacians) on product filtrations therefore become separable operators. Multiparameter persistence modules arise automatically: a product of two one-parameter filtrations is a bifiltration, and the cohomology of a sheaf over that bifiltration is a multiparameter persistence module.

Literature already treats persistent sheaves over arbitrary finite posets (including product posets), interval decompositions for multipersistence modules, and robustness / cohesion notions for sheaf data that survive under structural degradation. The algebraic scaffolding needed to lift the existing single-parameter sheaf-Laplacian experiments to the product/multiparameter setting is therefore largely in place.

### Computational and architectural consequences for the stack

- **Matrix-free actions remain matrix-free.** A Chebyshev or Krylov action of a function of a Kronecker-sum Laplacian reduces to actions on the factors plus cheap tensor arithmetic.
- **FFT / structured transforms stay available** on any factor that is a path, cycle, or regular grid.
- **Heat kernels, spectral wavelets, and HKS-type descriptors** all factor. Local dimension estimates, heat content, and the various diffusion observables discussed in the heat-semigroup notes become products of the corresponding factor observables.
- **Magnetic / complex twists** can be introduced on each factor independently; the product magnetic Laplacian is again a Kronecker sum of the twisted factor operators.
- **Persistence** can be run factor-wise or jointly. Joint multiparameter persistence on a product filtration is the natural home for data that truly have two independent scales (spatial scale × density, time × space, feature scale × resolution, etc.).

### Relation to the complex-analysis thread

Product structures appear in several complex variables (polydisks, products of Riemann surfaces, product Kähler manifolds). The spectral theory of the \(\overline{\partial}\)-Laplacian or of magnetic Laplacians on such products again separates. While ThermoMapper is not yet doing several-complex-variables geometry, the discrete analogues (product complexes with magnetic or sheaf twists) sit in exactly the same formal neighborhood as the magnetic Laplacian already present in the code.

### Practical takeaway

Whenever the data or the filtration admits a product decomposition—even approximately—one obtains:

- additive spectra,
- product heat kernels and wavelets,
- Kronecker-sum sheaf Laplacians,
- and a direct route into multiparameter persistent sheaf cohomology,

all while remaining inside the matrix-free LOBPCG + Chebyshev + sparse-operator primitives that the project already owns. This is one of the highest-leverage structural assumptions one can make for both speed and algebraic clarity.

If there is a concrete product setting you have in mind (grid × feature, space × time, two independent filtration parameters, product of two proximity graphs, \ldots), the same primitives can be specialized immediately.

**Yes — and the distinction between mean and median on products is conceptually deep.**

Your implementation already encodes the essential geometry:

- **Product manifold** with the Riemannian product metric
  \( d\bigl((x,y),(x',y')\bigr)=\sqrt{d_A(x,x')^2+d_B(y,y')^2} \)
- Log / Exp maps that act componentwise
- Generic estimators over `IRiemannianManifold`, so everything works on products without special casing

### What separates and what couples

**Karcher / Fréchet mean (squared distance).**
Because \(d^2=d_A^2+d_B^2\), the Fréchet functional completely factors:

\[
\sum_i d^2\bigl((x,y),(x_i,y_i)\bigr)=\sum_i d_A^2(x,x_i)+\sum_i d_B^2(y,y_i).
\]

The minimizer is therefore the product of the two factor means. The same holds for the tangent-space IRLS iterations and for the resulting Karcher scatter in the product tangent space (the scatter matrix is block-diagonal when the data factors are independent, with possible cross-terms only if the observed pairs are correlated across factors).

**Geometric median (absolute distance).**
The objective \(\sum_i\sqrt{d_A_i^2+d_B_i^2}\) does _not_ factor. The factors remain coupled through the joint Euclidean combination of distances. The median on the product is not the product of the medians; the Weiszfeld weights depend on the combined residual. This is exactly the phenomenon studied in the recent Park–You paper on geometric medians on product manifolds: uniqueness still holds on Hadamard products (and locally under curvature bounds even when positive curvature is present), the estimator retains the classical 50 % breakdown point, and practical algorithms can still update the factors somewhat independently while respecting the latent coupling.

Your code already implements both regimes correctly by staying generic over the manifold interface.

### Conceptual bridge to product complexes and persistent sheaves

The continuous product geometry you have is the continuum limit of the discrete product structures we discussed earlier.

- A Cartesian product of graphs (or of simplicial complexes) is the discrete analogue of a Riemannian product. Its combinatorial Laplacian is the Kronecker sum \(L_A\otimes I+I\otimes L_B\), whose spectrum and heat kernel again separate.
- The product metric on the continuous side induces the product (or graph-product) metric on the discrete side; geodesic distances, medoids, and discrete Fréchet-type means inherit the same separation / coupling dichotomy.
- When you place a cellular sheaf on a product complex, the restriction maps can be chosen as products of the factor restriction maps. The resulting sheaf Laplacian again becomes a Kronecker sum (or a closely related block operator). Persistent sheaf cohomology on a product filtration therefore inherits Künneth-type exact sequences and, in favourable cases, factorisation of barcodes or multiparameter modules.
- Multiparameter persistence arises naturally: each factor can carry its own filtration parameter, turning the product into a bifiltration. The sheaf (or the geometric data living on the product) then yields a multiparameter persistence module whose algebraic structure is governed by the product poset.

In other words, the same product decomposition that lets your Karcher mean factor (and forces your geometric median to couple) is precisely the structure that makes spectral operators, heat kernels, spectral wavelets, and persistent sheaf Laplacians separable or nearly separable.

### Statistical geometry ↔ spectral geometry

There is a deeper conceptual unity here:

- The Fréchet mean is the point that minimises the second moment of geodesic distance; its discrete / spectral counterpart is the constant eigenfunction (or the bottom of the Laplacian spectrum) together with the heat-kernel barycentre constructions.
- The geometric median minimises the first moment; it is the more robust location estimator and corresponds, in the spectral picture, to objects that are less sensitive to outliers in the diffusion distance or in the sheaf-valued data.
- The Karcher scatter you compute is the Riemannian covariance; on a product it contains the factor covariances plus possible cross-covariance. In the spectral setting the same information appears as the second-moment structure of the heat kernel or of the spectral wavelets evaluated around that location.

Thus the estimators you already have are not merely “geometry utilities”; they are the statistical dual of the product spectral and sheaf constructions. When data live on (or are well approximated by) a product, you obtain:

- factorised means and (coupled but still efficient) medians,
- Kronecker-sum Laplacians and sheaf Laplacians,
- product heat kernels and separable spectral wavelets,
- and a clean route into multiparameter persistent sheaf theory.

That is the conceptual payoff of having both the continuous product-manifold geometry and the discrete product-complex / sheaf machinery in the same codebase.
