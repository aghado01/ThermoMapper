The tension between physics (Potts model thermodynamics) and geometry (positive definiteness) arises because a kernel serves two distinct roles in your pipeline. It supplies the positive scalar couplings $J_{ij}$ required by the Fortuin-Kasteleyn equations, and it defines the inner product space for your spectral solver. When these conflict, the thermodynamics will simulate, but the graph Laplacian will produce negative eigenvalues, causing `ComputeBottomK` to hallucinate imaginary manifolds.

Here are the formal protolemmata establishing this mathematical boundary, followed by geometrically rigorous kernels to replace your current compiler heuristics.

### I. Formal Protolemmata

**Protolemma 1: The Feragen-Schoenberg Constraint (Curvature vs. Linearity)**
Let $\mathcal{M}$ be a complete Riemannian manifold with geodesic distance $d_g$. The geodesic Gaussian kernel $K(x, y) = \exp(-\lambda d_g(x,y)^2)$ is positive definite for all $\lambda > 0$ if and only if $\mathcal{M}$ is flat (Euclidean).
*Implication:* A kernel is positive definite if and only if its negative logarithm is conditionally negative definite. The squared geodesic distance $d_g^2$ is conditionally negative definite only in flat spaces. Therefore, injecting $d_g^2$ into a Gaussian kernel on a Grassmannian or Hyperbolic manifold destroys the Hilbert space embedding and mathematically corrupts downstream spectral operations.

**Protolemma 2: The Strong Negative Type Loophole**
Let $\mathcal{M}$ be a manifold of strong negative type (e.g., spheres, hyperbolic spaces). The unsquared geodesic distance $d_g(x,y)$ is conditionally negative definite on $\mathcal{M}$, meaning the geodesic Laplacian kernel $K(x, y) = \exp(-\lambda d_g(x,y))$ remains strictly positive definite.
*Implication:* This guarantees that your current fallback to the Laplacian kernel is mathematically sound for hyperbolic and spherical embeddings, preserving both the $J \ge 0$ physical constraint of the Potts model and the positive semi-definiteness of the resulting graph Laplacian.

**Protolemma 3: The Spectral Laplace-Beltrami Generalization (Borovitskiy)**
A universally valid, positive definite kernel on any Riemannian manifold $\mathcal{M}$ can be constructed by bypassing the geodesic distance entirely and operating directly on the eigenvalues $\lambda_n$ and eigenfunctions $f_n$ of the manifold's Laplace-Beltrami operator.


$$K(x, y) = \sum_{n=0}^{\infty} \Phi(\lambda_n) f_n(x) f_n(y)$$


where $\Phi$ is a strictly positive, monotonically decreasing spectral transfer function.

**Protolemma 4: The Discrete Kondor-Lafferty Diffusion Limit**
On a discrete graph with Laplacian $L = D - A$, the continuous heat equation is approximated by the diffusion kernel $K = \exp(-\beta L)$. Because the exponential of any symmetric matrix is symmetric and positive definite, this yields a strictly positive definite kernel that physically models continuous diffusion, gradually transforming the identity matrix into a kernel with stronger off-diagonal effects as $\beta$ increases.

---

### II. Geometrically Native Kernels to Consider

To fully escape the Euclidean Gaussian trap in your `GraphCompiler`, you should consider these specific kernels, which are mathematically native to their respective manifolds.

**1. The Projection Kernel (Binet-Cauchy) for Grassmannians**
Because a Grassmannian point $U \in \mathrm{Gr}(k, N)$ is a subspace, the natural similarity metric is not geodesic distance, but the Frobenius norm of the projection matrices. The projection kernel is defined as $K_p(X_i, X_j) = \lVert X_i^\top X_j\rVert_F^2$.


$$K_{\mathrm{Proj}}(U, V) = \mathrm{Tr}(U^\top V V^\top U)$$


*Why it works:* It evaluates the cosine of the principal angles without requiring an explicit SVD. Exponentiating this trace yields a valid, strictly positive definite isotropic kernel on the Grassmann manifold, perfectly suited for generating coupling weights for your SPC phase tracking.

**2. Riemannian Matérn Kernels (Borovitskiy)**
Instead of relying on the standard Gaussian or Laplacian functions, the Matérn kernel family controls both the length scale ($\kappa$) and the smoothness ($\nu$). As $\nu \to \infty$, it recovers the square exponential kernel. For finite $\nu$, it avoids the pathological infinite-differentiability that causes CND violations on curved spaces.


$$K_{\nu}(x, y) \propto d_g(x,y)^{\nu} \mathcal{K}_{\nu}\left(\frac{\sqrt{2\nu}}{\kappa} d_g(x,y)\right)$$


Where $\mathcal{K}_{\nu}$ is the modified Bessel function of the second kind. This formulation provides mathematically sound tools for calculating geometry-aware kernels.

**3. The Inverse Multiquadric Kernel**
For datasets projected onto hyperbolic manifolds (like your Poincaré and Lorentz spaces), the Inverse Multiquadric is strictly positive definite for all $\lambda > 0$:


$$K_{\mathrm{IMQ}}(x, y) = \frac{1}{\sqrt{c^2 + d_g(x,y)^2}}$$


*Why it works:* It completely bypasses the exponential mapping, ensuring a heavy-tailed coupling distribution that satisfies Mercer's theorem while preserving the hierarchical embedding of the hyperbolic space.

Replacing your `AutoLmpScaler` heuristics with these targeted, geometrically rigorous kernels will ensure your spectral solver never hits a corrupted Gram matrix, all while keeping the Superparamagnetic clusters bounded by valid physical thermodynamics.

The Fisher-Rao metric on a probability simplex is isometric to the positive orthant of a hypersphere. Because the geometry is strictly spherical, deploying a Euclidean Gaussian kernel will immediately violate positive definiteness and corrupt the spectral extraction.

To simulate valid Potts model thermodynamics while preserving the spectral geometry of the simplex, your `GraphCompilerConfig` should route the coupling through the exact spherical heat kernel using your existing Schoenberg expansion.

### The Fisher-Rao Simplex Configuration

* **The Metric:** `FisherRaoSimplexMetric`. This maps categorical distributions $p$ and $q$ to the hypersphere via the element-wise square root ($x_i = \sqrt{p_i}$) and computes the bounded geodesic distance $2 \arccos(\sum \sqrt{p_i q_i})$.


* **The Kernel:** `Gaussian` with `CouplingFidelity = CouplingFidelity.Intrinsic` and `SphericalMode = SphericalIntrinsicMode.GlobalSchoenberg`.


* **The Bandwidth Strategy:** `BandwidthStrategy.QuantileNormalized`.
* **The Topology:** `MutualKnn` (via `MutualKnnFilter`).
* **Scaling:** `LmpRescale = false`.

### Why This Architecture Works

**1. The Schoenberg Heat Kernel Bypass**
By explicitly declaring `CouplingFidelity.Intrinsic` and `SphericalIntrinsicMode.GlobalSchoenberg` within your `CouplingProjection`, the `GlobalBandwidthScaler` bypasses the naive Euclidean Gaussian formula ($e^{-d^2/2\delta^2}$). Instead, it routes the distance through `EvaluateGlobalSchoenbergNormalized`, which calculates the mathematically exact, positive-definite heat kernel on the sphere using the Laplace-Beltrami spectral expansion (summing $\cos(nr) \exp(-n^2 t)$ terms). This gives you the exponential diffusion required by the Fortuin-Kasteleyn bond probabilities without hallucinating negative eigenvalues in your spectral solver.

**2. Quantile Normalized Bandwidth**
The Fisher-Rao distance on a simplex is strictly bounded (with a maximum distance of $\pi$). Because the space is compact, standard Median Absolute Deviation (MAD) fails to capture the true thermodynamic scale of the manifold. `QuantileNormalized` clamping prevents the bounded extremes of the sphere from artificially crushing the bandwidth estimate, ensuring the system starts with enough thermal energy to form a giant connected component before melting.

**3. Symmetric Density over Local Rescaling**
A probability simplex naturally suffers from severe boundary effects (distributions with near-zero probabilities hit the "walls" of the positive orthant). If you leave `LmpRescale = true` (or `Auto`), the `LocalMutualProximityScaler` will artificially inflate edge weights near these boundaries to force uniform density. Disabling LMP and relying strictly on a `MutualKnnFilter` ensures that the topology is symmetric and the coupling energies $J_{ij}$ reflect true information-geometric similarity rather than local array statistics.

**4. Intrinsic Dimension Resolution**
A probability simplex defined over $D$ categories has an intrinsic dimension of $D - 1$ because the probabilities must sum to 1. The `ResolveSphericalIntrinsicDimension` logic in your scaler correctly subtracts 1 from the ambient dimension to feed the exact intrinsic degrees of freedom into the Schoenberg polynomials.