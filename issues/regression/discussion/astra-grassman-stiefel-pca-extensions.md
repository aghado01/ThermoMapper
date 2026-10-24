**Yes—this is a useful direction for functional BASS and the wider regression engine.** I would treat PCA as one way to construct a response representation, with Grassmann, Stiefel and whitening addressing different parts of that representation.

| Concept | What it controls | Potential value |
|---|---|---|
| **Grassmann** | The retained subspace, regardless of which axes represent it | Robust aggregation of subspaces, comparing fits, tracking changing response structure |
| **Stiefel** | A particular orthonormal frame within the response space | Learning component directions suited to regression; representing uncertainty in those directions |
| **Whitening** | The scale and correlation of coordinates | Comparable coefficient scales, numerical conditioning, or accounting for known measurement noise |

PCA already selects a point on a Grassmann manifold. Keeping its individual component vectors gives a Stiefel frame. Merely changing the optimization machinery preserves ordinary PCA if the objective remains the same; the statistical gains come from changing the objective, metric or uncertainty model. [Edelman, Arias and Smith](https://math.mit.edu/~edelman/publications/geometry_of_algorithms.pdf)

There is an interesting detail in **BASS’s existing implementation**. It computes
\[
Y_{\mathrm{scaled}}^\top=UDV^\top,
\qquad K=U_rD_r,\qquad W=V_r^\top.
\]
The component regressions receive \(W\), while the singular values live in the reconstruction basis \(K\). Thus the training scores already have equal energy and are mutually orthogonal—essentially whitened scores up to the sample-size normalization when mean-centred. [BASS source](D:/aghado01/packages/BASS/R/bassPCA.R:88)

Consequently, “add whitened SVD” could mean two different changes:

- **Normalize the retained scores:** largely a refinement of BASS’s existing convention. Preserve the inverse transformation, prior scales and predictive covariance.
- **Whiten measurement noise before choosing the subspace:** potentially a substantial improvement when different output locations have different noise levels. This selects structure relative to noise, rather than raw response variance.

Whitening alone supplies neither robustness nor independence of prediction errors. Very small estimated variances also make inverse scaling unstable, so rank truncation or explicit regularization matters. [Kessy, Lewin and Strimmer](https://arxiv.org/html/1512.00809)

**The important regression-specific distinction is that the axes can matter even when the subspace is unchanged.** BASS fits separate adaptive spline models to component scores. Rotating those scores can turn several simple component functions into complicated mixtures. Independent component priors and independently selected spline structures generally do not preserve the model under that rotation.

That suggests three useful extensions, in increasing order of scope:

1. **A robust, fixed response basis.** Estimate the centre and subspace robustly, freeze them, and run the existing component regressions. Your [MoMPCA](D:/aghado01/ThermoMapper/src/maths/geometry/dim-reduction/MoMPCA.cs:12) already aggregates centres and Grassmann subspaces. This connects directly to our median discussion and the [robust-basis extension already recorded in the brief](D:/aghado01/ThermoMapper/issues/regression/bass-full-port-brief.md:333).

2. **A basis chosen for predictive usefulness.** Ordinary PCA retains directions of greatest variation. A large-variance direction might be mostly unpredictable noise, while a smaller direction might contain a clean response to the predictors. A regression-aware objective could select a subspace—and potentially rotate its axes—to make the component functions easier to predict. Your generic subspace optimizer is relevant, but a Grassmann objective must give the same answer for every frame spanning the same subspace. An axis-dependent objective needs frame-aware treatment.

3. **A Bayesian basis with uncertainty.** Treat the component frame as unknown and propagate its uncertainty into reconstructed predictions. This matters particularly when the data weakly identify the components. There is concrete precedent: **FAST Bayesian FPCA** represents eigenfunctions through an orthonormal spline basis and samples an orthonormal coefficient matrix using a polar-decomposition parameterization. Adapting that to BASS would be a new joint model, beyond fixed-basis preprocessing. [Sartini et al.](https://arxiv.org/abs/2412.11340)

The current code provides useful starting points, but I found prerequisites worth separating from those opportunities:

- [PCA whitening](D:/aghado01/ThermoMapper/src/maths/linalg/PCA.cs:152) rescales `Components` in place. Those vectors cease to be an orthonormal frame. A shared representation should retain the **frame and coordinate scaling separately**.
- [MoMPCA’s reported variances](D:/aghado01/ThermoMapper/src/maths/geometry/dim-reduction/MoMPCA.cs:102) are calculated using ordinary squared deviations after robust subspace aggregation. They can still be contaminated and should not automatically become robust whitening scales.
- **The current [Stiefel distance](D:/aghado01/ThermoMapper/src/maths/geometry/manifolds/StiefelManifold.cs:45) needs correction.** Its formula gives zero for distinct frames within the same subspace—even a unit vector and its negative. That defeats the distinction we would need Stiefel to preserve.

My preferred first expansion is **robust fixed bases with explicit score scaling**, followed by predictive basis selection if benchmarks justify it. LinAlg should own decompositions and numerical transforms; Geometry should own manifold operations and generic optimization; Regression should own the response metric, basis-selection objective, priors and reconstruction. Geometry need not import Regression: the regression objective can be supplied through the existing delegate boundary.

If you also mean PCA of the **predictors**, rotations create another capability: oblique spline surfaces in transformed inputs. That changes BASS’s axis-aligned model family and original-variable sensitivity interpretation, so I would evaluate it separately from response-basis work.

**Yes. Across the existing tools, the strongest opportunities are learned input geometry for BGP, better coefficient coordinates for BAPS, and shared functional-response models across engines.** They have different statistical implications:

| Area | Useful extension | Effect |
|---|---|---|
| BGP | Learn a low-dimensional predictive input subspace | New kernel/model |
| BAPS | Orthogonalize the design and smoothing penalty | Can preserve the existing model |
| Hierarchical BAPS | Learn a few shared patterns of variation between curves | New hierarchical model |
| BARS | Fit splines along learned input projections | New regression capability |
| All continuous-response engines | Shared output basis and reconstruction | Functional/multivariate extension |
| Changepoints | Detect shared changes in several response components | New segment model; reusable DP recurrence |

**BGP is the clearest use of Grassmann geometry.** It already accepts multivariate predictors, but its [current kernel](D:/aghado01/ThermoMapper/src/maths/regression/bgp/GpKernel.cs:5) uses ambient Euclidean distance. A projected kernel could instead use

\[
k(x,x')=\exp\!\left[-\frac{\|Q^\top(x-x')\|^2}{2t}\right],
\qquad Q^\top Q=I.
\]

The engine would learn which combinations of inputs predict the response. That can retain a low-variance but highly predictive direction that input PCA would discard. This approach has direct precedent in GP regression with learned orthogonal projections. [Bilionis, Tripathy and Gonzalez](https://arxiv.org/abs/1602.04550)

With an isotropic kernel in the projected space, only the subspace matters: **Grassmann is the appropriate geometry**. With different length scales along its axes, the orientation also matters.

Your present BGP already pursues intrinsic-dimension adaptation through its bandwidth prior, without requiring an explicit linear projection. A projected GP would therefore be a complementary mode. Its distance cache and bandwidth-prior interpretation would need updating. Reducing input dimension also does **not** remove the dense \(n\times n\) factorization cost.

**For BAPS, the first opportunity is penalty-aware orthogonalization.** Ordinary PCA orders directions by data variance. A *Demmler–Reinsch* representation instead chooses coefficient coordinates that simplify both the design geometry and the smoothing penalty. It distinguishes unpenalized modes from increasingly rough modes—a much better match to spline regression. [Bach and Klein](https://arxiv.org/abs/2109.04288)

That could make repeated global-smoothing calculations, effective degrees of freedom and coefficient draws cheaper. There is already a starting point: [PenalizedSpline](D:/aghado01/ThermoMapper/src/maths/regression/spline/baps/PenalizedSpline.cs:19) computes a generalized penalty spectrum.

The qualifications matter:

- Retaining all modes and transforming the prior correctly can preserve the model; truncating modes changes it.
- The penalty’s null space needs explicit treatment rather than an inverse square root.
- A reusable diagonalization is most attractive for a fixed design and fixed penalty shape. Locally adaptive penalties change that shape.
- Dense transformed coordinates can lose the advantages of your existing banded implementation. This warrants a measured alternative, not a blanket replacement.

**Hierarchical BAPS offers a more substantial Stiefel application.** Its [current replicate model](D:/aghado01/ThermoMapper/src/maths/regression/spline/baps/HierarchicalPSpline.cs:89) uses independent coefficient deviations around the population curve. An extension could learn

\[
\beta_r=\beta_0+A z_r+\epsilon_r,
\]

where the columns of \(A\) describe shared patterns of curve variation and \(z_r\) describes each replicate’s expression of those patterns.

For example, replicates might differ mainly in peak height, peak position and shoulder width. Learning those shared patterns could improve borrowing of information when individual curves are sparsely observed. Here, Stiefel constraints help identify the loading frame—but orthogonality should use the appropriate **function-space metric**, rather than depend arbitrarily on the spline coefficient representation. This would connect hierarchical smoothing and functional PCA within one statistical model.

**BARS could become a projected-input spline engine.** A single-index model,
\[
y=f(w^\top x)+\epsilon,
\]
would reuse the existing one-dimensional adaptive spline machinery while learning a direction through multivariate predictors. Multiple such components lead toward projection-pursuit regression. That is a meaningful capability alongside BASS’s coordinate-based interactions.

For ordinary one-dimensional BARS, however, Grassmann machinery adds little directly. QR or coefficient whitening can improve numerical coordinates, but it must preserve the g-prior and should not casually displace efficient banded solves or local proposal updates.

**The output-basis machinery should eventually be shared across regression engines.** A fixed functional-response decomposition could feed component regressions to BARS, BAPS, BGP or BASS, then reconstruct predictions and uncertainty consistently. That belongs in a shared regression response layer when those consumers are implemented. It has a different role from the current spline [`IBasis`](D:/aghado01/ThermoMapper/src/maths/regression/spline/IBasis.cs), which constructs the input-side design from knots.

For [changepoints](D:/aghado01/ThermoMapper/src/maths/regression/changepoint/ExactChangepoint.cs:23), a fixed response basis could support shared change boundaries across several components. The dynamic-programming recurrence remains reusable if the joint evidence still factorizes by segment; the current scalar segment-summary contract would need expansion.

I would prioritize **shared functional-response preparation, BAPS spectral-coordinate benchmarks, and a projected BGP experiment**. Learned hierarchical factors and projected BARS are larger model additions. Poisson counts and Whittle periodograms need their own treatment: arbitrary response rotations do not preserve those likelihoods, although suitable latent Gaussian or coefficient blocks can still be transformed.

**Yes—this line of work gives us a concrete bridge between BASS, active subspaces, and your Grassmann ideas. I would put posterior active-subspace analysis ahead of jointly learning projection directions during fitting.**

There is an identifier distinction:

- **2306.01911v1** is [*Generalized Bayesian MARS: Tools for Emulating Stochastic Computer Models*](https://arxiv.org/pdf/2306.01911v1).
- The directly relevant companion is **2307.11241v1**, [*Discovering Active Subspaces for High-Dimensional Computer Models*](https://arxiv.org/pdf/2307.11241v1), by overlapping authors.

The companion paper explicitly uses **BASS**. Its central calculation is

\[
C_\rho=\mathbb E_{x\sim\rho}\!\left[\nabla f(x)\nabla f(x)^\top\right].
\]

The leading eigenvectors identify input directions along which the fitted response changes most, averaged over the chosen input distribution. For tensor-product representations—including MARS—the authors reduce this calculation to one-dimensional integrals. They also propagate BASS posterior uncertainty and cache integrals across structural changes. [Active-subspace paper, §§3 and supplementary §11](https://arxiv.org/pdf/2307.11241v1)

That is a useful distinction from our PCA discussion: **PCA describes variation in the data; active subspaces describe variation in the response as inputs change.** Grassmann geometry then gives us a way to compare the resulting subspaces across posterior draws.

For this repository, my proposed placement would be:

| Capability | Reuse or addition | Ownership |
|---|---|---|
| Integrals of basis products | Reuse planned Sobol machinery | BASS sensitivity |
| Integrals involving basis derivatives | Add beside existing moments | BASS sensitivity |
| Assemble \(C_\rho\), eigenvalues and activity summaries | New posterior analysis | Regression |
| Symmetric eigendecomposition | Reuse numerical infrastructure | LinAlg |
| Compare subspaces across posterior draws | Reuse qualified Grassmann operations | Geometry |
| Posterior weights | Consume existing MCMC or future SMC results | Regression consumes Samplers |

This fits the [existing brief’s boundaries](/D:/aghado01/ThermoMapper/issues/regression/bass-full-port-brief.md:77): fitting remains independent of sensitivity analysis. **Neither SMC nor manifold sampling is a prerequisite.**

There is concrete reuse underneath that placement. The R implementation already computes pairwise basis moments for Sobol analysis in [sobol_normTmix.R](/D:/aghado01/packages/BASS/R/sobol_normTmix.R:817). My proposed C# extension would add derivative–basis and derivative–derivative moments to that family, retaining sparse term structure and reusable buffers.

Several design choices would make this more useful across the engine:

- **Posterior subspace uncertainty.** Compare projectors \(QQ^\top\) or principal angles across draws. Averaging eigenvector coordinates can cancel equivalent directions because of arbitrary signs and rotations. Weak eigenvalue separation should appear as uncertainty.
- **Collapsed coefficient uncertainty.** A posterior mean of \(C_\rho\) requires coefficient second moments, not merely the outer product of conditional mean coefficients. That calculation would still be distinct from obtaining the full posterior distribution of subspaces.
- **Functional responses.** My proposed extension is
  \[
  C_\rho=\mathbb E[J_f(x)^\top W_yJ_f(x)],
  \]
  where \(W_y\) specifies output units, quadrature weights, or another meaningful response metric. This connects active subspaces to the output-basis work without privileging arbitrary PCA scaling.
- **Other engines.** Tensor BAPS is a natural second consumer. BGP could use derivative-based analysis and eventually use an estimated subspace to initialize a projected model. Ordinary gradient analysis would need a different treatment for categorical predictors and discontinuous changepoint models.

Whitening also has a precise role here: the paper’s correlated-Gaussian route fits the surrogate in whitened coordinates and transforms the gradient matrix back by congruence. It is therefore more consequential than changing an eigensolver setting. [§3.3.2](https://arxiv.org/pdf/2307.11241v1)

The **GBMARS paper you supplied** matters separately for our likelihood and prior decisions. It extends MARS through normal mean–variance mixtures while retaining conditional coefficient marginalization. Crucially, it explains that a coefficient prior based on sampled observation weights complicates those weights’ full conditionals, and recommends a ridge prior. That supports our concern about weight-dependent priors, but **does not force us to adopt ridge**; it belongs in D03 alongside the intercept and prior-scaling decisions. It also explicitly retains an approximate count-intensity update caused by Poisson truncation—something our existing prerequisite work already requires us to resolve explicitly. [GBMARS, §§3–3.2](https://arxiv.org/pdf/2306.01911v1)

I would treat analytic active-subspace analysis as a scoped extension beside the Sobol stage, with polynomial/hinge reference cases, numerical integration checks, and gradient-accuracy checks before using its directions to reduce another model. No files changed.

# ive also been looking at 2405.13342v2 and https://doi.org/10.1093/biomtc/ujag121 , at least the supplemental materials and repositories for the latter

**Both are pertinent. They point toward two distinct extensions: a spectral GP capability for curved input domains, and adaptive tensor-product splines alongside the BASS work.**

I checked the GP paper’s requested version, your local revised ABMRS supplement, and the repositories. The three core ABMRS source files I compared match current upstream. My assessment of the Biometrics method is based on that supplement and code, rather than its paywalled main text.

| Approach | What it adapts | Most relevant connection |
|---|---|---|
| Active-subspace analysis | Response-sensitive input directions | BASS/BAPS sensitivity and projected regression |
| **2405.13342v2: FLGP** | Covariance to the geometry of the input domain | BGP and graph spectral infrastructure |
| **ujag121: ABMRS** | Knot locations and model complexity | BARS, tensor BAPS, shared spline machinery |

**FLGP gives a substantive reason to expand the GP architecture.** It constructs a sparse relationship between observations and representative points, then uses truncated SVD to obtain a heat-kernel approximation:

\[
K_t \approx \Phi\,\operatorname{diag}(e^{-t\lambda_j})\,\Phi^\top .
\]

This lets covariance follow the domain’s geometry. A linear active subspace and a curved domain describe different structures; either could be useful independently. [Requested FLGP paper](https://arxiv.org/html/2405.13342v2)

Your current [GP kernel contract](/D:/aghado01/ThermoMapper/src/maths/regression/bgp/GpKernel.cs:5) takes an ambient squared distance and bandwidth. A supplied spectral representation would therefore need an additional GP path. It would also allow a low-rank Gaussian solve; the authors’ implementation already switches to a smaller coefficient-space system when the number of observations exceeds the retained rank. [Prediction implementation](https://github.com/junhuihe2000/FLGP/blob/8fc1c1204f77e978016c651ed4b3a0466d377b56/src/Predict.cpp#L40-L74)

Three qualifications matter for our design:

- **The linear-time claim holds other dimensions fixed.** Their stated cost still depends on representative-point count, neighbourhood size, ambient dimension, spectral rank, and labelled sample count. We should benchmark these separately. [Complexity analysis, §3.3](https://arxiv.org/html/2405.13342v2)
- **Prediction currently participates in geometry construction.** The implementation combines training and prediction inputs before constructing the spectrum. A reusable C# fitted model needs an explicit choice between that transductive workflow and extending a frozen representation to previously unseen points. [Spectrum construction](https://github.com/junhuihe2000/FLGP/blob/8fc1c1204f77e978016c651ed4b3a0466d377b56/src/Spectrum.cpp#L48-L75)
- **Its fitting routine optimizes hyperparameters.** Sampling diffusion time and noise would be a further statistical choice; conditional predictive uncertainty does not include all uncertainty in those optimized quantities. [Fitting implementation](https://github.com/junhuihe2000/FLGP/blob/8fc1c1204f77e978016c651ed4b3a0466d377b56/src/Fit.cpp#L42-L79)

The dependency boundary can remain clean: `graphs/spectral` produces geometry-derived numerical data; an orchestration layer passes that data into Regression. LinAlg owns decomposition and solves; Regression owns the GP prior, likelihood and posterior. **Maths.Regression would not import Graphs.** Likewise, using the existing clustering implementation to select representative points would happen above that boundary, with the points supplied as data.

**ABMRS is especially relevant to extending adaptive splines beyond scalar-input BARS.** Its implementation constructs full tensor products of coordinate-wise B-spline bases. That differs materially from BASS’s selected collection of hinge-product terms: adding one coordinate knot can add an entire slice of tensor coefficients. [Tensor construction](/D:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/src/utils.cpp:297)

My proposed reuse would be B-spline evaluation, tensor design, Gaussian coefficient calculations, proposal infrastructure and posterior summaries. The adaptive tensor model would retain its own structural state and prior. This could be particularly useful for low-dimensional surfaces with interpretable knot locations, while the full BASS port retains its broader predictor-selection and interaction capabilities.

Your current [EbicPrior](/D:/aghado01/ThermoMapper/src/maths/regression/spline/bars/ComplexityPriors.cs:83) already incorporates the relevant count-law distinction:

\[
\pi(k,\xi)\propto {C\choose k}^{-\gamma}
\quad\Longrightarrow\quad
\pi(k)\propto {C\choose k}^{1-\gamma}.
\]

So that part of the paper is already reflected in your preparation work.

The supplement and repositories also supply useful **references to qualify**, rather than outputs to reproduce indiscriminately:

- **The Gaussian derivation and sampler differ.** Appendix A derives the collapsed unit-information \(g\)-prior likelihood. The bundled multivariate sampler instead uses an \(m^{-\Delta\nu/2}\) and residual-scale ratio, and draws coefficients using a fitted noise scale. Matching the derivation requires explicit coefficient-prior, intercept and variance conventions. [Appendix](/D:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/web_appendix_extracted.txt:21), [sampler](/D:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/src/mebars.cpp:184)
- **Numerical handling changes the effective fit.** The utility drops low-support basis columns and adds a fixed diagonal ridge. Our port should specify those choices statistically and numerically. [Implementation](/D:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/src/utils.cpp:19)
- **The upstream peak example is a useful prototype.** It refits least squares for sampled knot configurations and compares predictions in a small neighbourhood of knot intersections. A C# posterior peak analysis should also propagate coefficient uncertainty and define boundary and nonsmooth extrema. [Peak example](https://github.com/junhuihe2000/exampleABMRS/blob/840ce21831b8db4901c52be0c1fda1d5343900db/scripts/surface_peak.r#L94-L120)

I would keep three separately scoped additions on the roadmap: **BASS active-subspace analysis, spectral BGP, and adaptive tensor splines**. Shared numerical improvements can serve all three, but their statistical models remain distinct. Both reference packages already use compiled C++/Eigen, so the C# opportunity lies chiefly in better contracts, numerical algorithms, caching and uncertainty handling—not an automatic speedup from changing languages.

No files changed.