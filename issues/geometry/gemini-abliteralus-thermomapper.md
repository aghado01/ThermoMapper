
Here is a detailed mathematical and architectural comparison between **`abliteralus` (public repo)**, your **private research work**, and **`thermomapper`**, along with concrete mathematical primitives that can be transferred into ThermoMapper.

---

### 1. Capability Landscape Matrix

| Domain | `abliteralus` (Public) | Your Private Work (`abliteralus/private`) | `thermomapper` (`src/maths`, `tda`) |
| :--- | :--- | :--- | :--- |
| **Grassmannian Subspaces** | [`grassmann.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/grassmann.py): Gauge-invariant $Gr(k, d)$, Björck–Golub `atan2(sin, cos)` angles, chordal/projection distance across different ranks, Log/Exp, Karcher mean, diameter bounds. | [`guard_study/grassman.py`](file:///d:/aipithicus/abliteralus/private/src/guard_study/grassman.py): Multi-layer subspace trajectories through depth (lattice $\to$ seam $\to$ plateau), causal-map projection gaps. | [`GrassmannManifold.cs`](file:///D:/aghado01/ThermoMapper/src/maths/geometry/manifolds/GrassmannManifold.cs), [`SubspaceAnnealer.cs`](file:///D:/aghado01/ThermoMapper/src/maths/geometry/dim-reduction/SubspaceAnnealer.cs), [`MoMPCA.cs`](file:///D:/aghado01/ThermoMapper/src/maths/geometry/dim-reduction/MoMPCA.cs): $Gr(r, n)$ with thin SVD, Log/Exp, geodesic annealing, product-manifold median of means. |
| **Curved & Non-Euclidean Manifolds** | [`riemannian_manifold.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/riemannian_manifold.py): Empirical sectional curvature $K$ via Gauss–Bonnet angle excess, geodesic vs Euclidean ratio, 2nd-order Jacobi field curvature correction. | Direction bundle observation contracts, local curvature checks on activation manifolds. | [`IRiemannianManifold.cs`](file:///D:/aghado01/ThermoMapper/src/maths/geometry/IRiemannianManifold.cs): Stiefel, Lorentz, Poincaré ball, Spherical, Scaled, Product manifolds. IRLS Karcher mean, Weiszfeld median, scatter matrices. |
| **Conic & Polyhedral Geometry** | [`concept_geometry.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/concept_geometry.py): Polyhedral concept cones, minimal enclosing cone solid angle, Direction Specificity Index (DSI). | Subspace orientation vs condition-labeled axes (e.g. Miller paper discussions). | None (ThermoMapper operates primarily on linear subspaces or Riemannian manifolds, not polyhedral cones). |
| **Spectral & Random Matrix Theory** | [`spectral_certification.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/spectral_certification.py): BBP phase transition (Baik-Ben Arous-Péché), Marchenko–Pastur noise floor, anisotropic noise separation. | Finite-difference eigenvalue tracking across layer depths. | LinAlg (`DenseEigen`, `ThinSvd`), participation ratio in [`IntrinsicDimension.cs`](file:///D:/aghado01/ThermoMapper/src/tda/mapper/Diagnostics/IntrinsicDimension.cs), but no formal RMT/BBP detector. |
| **Optimal Transport & Information Geometry** | [`wasserstein_optimal.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/wasserstein_optimal.py), [`leace.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/leace.py), [`whitened_svd.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/whitened_svd.py): $W_2^2$ Bures metric decomposition, generalized Rayleigh quotient, FLD within-class covariance. | Contrastive directions with empirical scale calibration (projection std, projection gap). | [`Wasserstein1.cs`](file:///D:/aghado01/ThermoMapper/src/maths/distance/geodesic/Wasserstein1.cs), [`FisherRaoSimplex.cs`](file:///D:/aghado01/ThermoMapper/src/maths/distance/geodesic/FisherRaoSimplex.cs), [`FisherRaoHalfPlane.cs`](file:///D:/aghado01/ThermoMapper/src/maths/distance/geodesic/FisherRaoHalfPlane.cs), Shannon & JS divergences. |
| **Dynamics & Differential Interventions** | SVD component selection, layer ranges. | [`layer_selection.py`](file:///d:/aipithicus/abliteralus/private/src/guard_study/layer_selection.py): $n$-th order finite difference $\Delta^n$ for seam detection; [`interventions.py`](file:///d:/aipithicus/abliteralus/private/src/guard_study/interventions.py): Orthogonal sham controls. | TDA Zigzag filtrations ([`GraphZigzag.cs`](file:///D:/aghado01/ThermoMapper/src/tda/ph/GraphZigzag.cs), [`RuVineyard.cs`](file:///D:/aghado01/ThermoMapper/src/tda/ph/RuVineyard.cs)), but no differential discrete gradient/seam operators. |

---

### 2. Deep Mathematical Comparison

#### A. Grassmannian Geometry: Numerical Stability & Variable Rank
- **The Angle Catastrophe at $\theta \approx 0$**:
  - In `thermomapper` ([`GrassmannManifold.Distance`](file:///D:/aghado01/ThermoMapper/src/maths/geometry/manifolds/GrassmannManifold.cs#L47-L59)), distance is computed via $\arccos(\sigma_i(Y^T Z))$. For small angles $\theta \approx 0$, $\cos \theta \approx 1 - \theta^2/2$. Floating-point cancellation in $\arccos(1 - \epsilon)$ wipes out half the precision bits.
  - In `abliteralus` ([`grassmann.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/grassmann.py#L140-L159)), principal angles use the Björck–Golub (1973) / Knyazev–Argentati (2002) stable formulation:
    $$\cos \theta_i = \sigma_i(Y^T Z), \quad \sin \theta_i = \sigma_{k-i+1}(Y - (Y Z^T) Z), \quad \theta_i = \operatorname{atan2}(\sin \theta_i, \cos \theta_i)$$
    This preserves machine precision even when two subspaces are nearly coincident (crucial for Karcher iteration convergence and fine manifold tracking).
- **Comparing Subspaces of Different Ranks**:
  - In `thermomapper`, `GrassmannManifold` assumes fixed rank $r$. It cannot compare a 1D direction with a 2D plane (a frequent occurrence in neural decoding, e.g., the LPFC working memory analysis, or layer-by-layer rank changes).
  - `abliteralus` provides the chordal/projection metric:
    $$d_{\text{proj}}(Y, Z) = \frac{\|P_Y - P_Z\|_F}{\sqrt{2}} = \sqrt{\frac{k_y + k_z}{2} - \sum \cos^2 \theta_i}$$
    This is well-defined across unequal ranks (e.g., distance between a line and a plane containing it is $1/\sqrt{2}$, not 0).
- **Normalized Grassmann Diameter**:
  - `abliteralus` implements the maximum geodesic distance bound:
    $$D_{\max}(k, d) = \sqrt{\min(k, d - k)} \cdot \frac{\pi}{2}$$
    which allows normalizing Grassmannian distances into an intrinsic $[0, 1]$ scale.

#### B. Intrinsic Curvature & Geodesic Projection
- **Curvature Estimation on Empirical Clouds**:
  - `thermomapper` has Riemannian structures (Lorentz, Poincaré, Spheres, Grassmann, Stiefel) with known analytic metrics, but has no mechanism to measure whether an empirical data cloud (e.g., neural activations or latent embeddings) is actually flat or curved.
  - `abliteralus` ([`riemannian_manifold.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/riemannian_manifold.py)) introduces two empirical curvature proxies:
    1. **Gauss–Bonnet Angle Excess**: Discrete triangles on local tangent spaces measuring $K \approx (\sum \theta_i - \pi) / \text{Area}$.
    2. **Geodesic-to-Euclidean Ratio**: Comparing Euclidean chord length $\|a - b\|$ to piecewise geodesic lengths projected along local kNN PCA tangent planes.
- **Second-Order Geodesic Projection**:
  - Standard projection removes a linear component: $x_{\text{proj}} = x - (x^T r) r$.
  - When the underlying manifold has sectional curvature $K$, `abliteralus` applies a Jacobi field second-order correction:
    $$x_{\text{geodesic}} = x_{\text{linear}} - \frac{K}{6} \|x^T r\|^2 v_{\perp}$$
    where $v_{\perp}$ is the principal normal direction of maximum variance.

#### C. Polyhedral / Conic Geometry vs. Subspaces
- In `abliteralus` ([`concept_geometry.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/concept_geometry.py)), concepts are modeled as **polyhedral cones** rather than subspaces:
  - **Direction Specificity Index (DSI)**:
    $$\text{DSI}(d_i) = 1 - \frac{1}{|C|-1}\sum_{j \neq i} |\cos(d_i, d_j)|$$
    Quantifies whether a geometric axis is category-specific ($\text{DSI} \to 1$) or universal/shared ($\text{DSI} \to 0$).
  - **Minimal Enclosing Cone Solid Angle & Effective Dimensionality**: Formulated via SVD participation ratios over cone rays.

#### D. Random Matrix Theory (RMT) & Noise Separation
- In `thermomapper`, [`IntrinsicDimension.cs`](file:///D:/aghado01/ThermoMapper/src/tda/mapper/Diagnostics/IntrinsicDimension.cs) computes the participation ratio $(\sum \lambda_i)^2 / \sum \lambda_i^2$.
- In `abliteralus` ([`spectral_certification.py`](file:///d:/aipithicus/abliteralus/abliteralus/analysis/spectral_certification.py)), spectral structure is tested against the **Marchenko–Pastur** upper edge $\lambda_+ = \sigma^2 (1 + \sqrt{\gamma})^2$ and the **BBP phase transition threshold** $\theta_{\text{BBP}} = \sigma^2(1 + \sqrt{\gamma})$.
  - This mathematically distinguishes true low-rank geometric signal from anisotropic sampling noise.

#### E. Phase Transitions & Finite Differences (Private Work)
- In `abliteralus/private` ([`layer_selection.py`](file:///d:/aipithicus/abliteralus/private/src/guard_study/layer_selection.py)), layer ranking is performed via $n$-th finite differences:
  $$\Delta^1 E_L = E_L - E_{L-1}$$
  Locates the **seam** (where a representation *emerges* or transitions) rather than the **plateau** (where it has accumulated and saturated).
- In [`interventions.py`](file:///d:/aipithicus/abliteralus/private/src/guard_study/interventions.py), **orthogonal sham controls** construct deterministic null-space vectors to verify that geometric effects are not artifacts of generic norm perturbation.

---

### 3. Ideas to Transfer into ThermoMapper as Mathematical Primitives

Here are 5 high-impact mathematical primitives ready to be implemented in ThermoMapper (`Maths.Geometry` or `Maths.LinAlg`):

#### 1. Numerically Stable & Universal-Rank Grassmann Primitives
Add to `Maths.Geometry.GrassmannManifold` (or a static `Grassmann` helper):
- **`PrincipalAngles(ReadOnlySpan<double> y, ReadOnlySpan<double> z)`**:
  Implement the Björck–Golub `atan2(sin, cos)` with the orthogonal residual $Y - P_Z Y$, replacing naive `Acos` to eliminate catastrophic precision loss at small angles.
- **`ProjectionDistance(ReadOnlySpan<double> y, int ky, ReadOnlySpan<double> z, int kz, int n)`**:
  The rank-agnostic chordal distance $\|P_Y - P_Z\|_F / \sqrt{2}$, enabling comparison and clustering between subspaces of different dimensions (e.g. comparing 1D axes with 2D/3D planes).
- **`MaxGeodesicDistance(int r, int n)`**:
  Returns $\sqrt{\min(r, n - r)} \cdot \pi / 2$ for relative diameter scaling.

#### 2. Empirical Sectional Curvature & Geodesic Deviation
Add a primitive under `Maths.Geometry.Estimators.Intrinsic.SectionalCurvature`:
- **`EstimateCurvature(double[][] points, int kNearest = 10)`**:
  Sample discrete triangles on local kNN neighborhoods and compute Gauss–Bonnet angular excess $\Delta \theta / \text{Area}$.
- **`GeodesicEuclideanRatio(double[][] points, int samplePairs = 50)`**:
  Compare straight-line Euclidean distance with local piecewise tangent-projected paths. Returns a diagnostic ratio $\rho \ge 1.0$ (with $\rho \approx 1$ signifying a flat submanifold and $\rho > 1$ quantifying curvature).

#### 3. Bures-Wasserstein Metric on Covariances ($W_2(\mathcal{N}_1, \mathcal{N}_2)$)
Add to `Maths.Distance.Geodesic`:
- **`BuresWassersteinDistance`**:
  $$W_2^2((\mu_1, \Sigma_1), (\mu_2, \Sigma_2)) = \|\mu_1 - \mu_2\|_2^2 + \operatorname{tr}\left(\Sigma_1 + \Sigma_2 - 2(\Sigma_1^{1/2} \Sigma_2 \Sigma_1^{1/2})^{1/2}\right)$$
  ThermoMapper already has $W_1$ and Fisher–Rao; adding Bures–Wasserstein gives the exact geometric distance between Gaussian distributions on the positive semi-definite cone $\mathcal{S}_{++}^n$.

#### 4. Random Matrix Theory / Marchenko–Pastur Spectral Filter
Add to `Maths.LinAlg.SpectralRmt`:
- **`BbpPhaseTransition(double[] eigenvalues, int nSamples, int ambientDim)`**:
  Calculates aspect ratio $\gamma = d/n$, estimates noise variance $\sigma^2$ from the bulk median, computes the Marchenko–Pastur boundary $\lambda_+ = \sigma^2 (1 + \sqrt{\gamma})^2$, and flags supercritical eigenvalue spikes (signals) vs subcritical bulk (noise).
  - Can be plugged directly into `MoMPCA`, `SubspaceAnnealer`, or `Mapper` diagnostics.

#### 5. Discrete Finite-Difference Seam Detector & Orthogonal Controls
Add to `Maths.LinAlg` or `TDA.Ph`:
- **`DifferentialSeamFilter(ReadOnlySpan<double> series, int order = 1)`**:
  Finite-difference operator $\Delta^n$ for discovering sudden geometric phase transitions / onset points along filtration steps, time steps, or network layers.
- **`OrthogonalShamControl(ReadOnlySpan<double> direction, int seed)`**:
  Deterministic generation of null-space vectors orthogonal to a candidate geometric direction for hypothesis falsification.