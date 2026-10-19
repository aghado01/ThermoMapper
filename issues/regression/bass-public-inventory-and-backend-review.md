# BASS: public capability inventory and C# backend review

Date: 2026-10-18. Status: source investigation and implementation planning; no engine changes.

BASS would add adaptive multivariate regression, mixed continuous/categorical interactions, two approaches to functional responses, posterior prediction, analytic sensitivity analysis, and emulator calibration. The port should cover that complete capability set. Predictor dimension, response dimension, and the dimension of coordinates indexing a response are separate concepts; none is limited to two.

The most valuable implementation strategy is to retain BASS's statistical representation and its useful algebraic structure, while rebuilding the numerical and sampling contracts. Incremental factorization is an important improvement, but it is only one part of the work. Prediction uncertainty, sensitivity measures, calibration, and validation also need attention.

## Evidence and scope

- BASS checkout: `D:/aghado01/packages/BASS`, commit `26710146c623a258854162367a939041efedd091`, clean when inspected. [DESCRIPTION](D:/aghado01/packages/BASS/DESCRIPTION:1) identifies version 1.3.1, dated 2023-07-01, GPL-3.
- ThermoMapper checkout: commit `c82c6368b151d28d14dfa1f2076de633efd125b1`. Existing unrelated changes were left alone.
- Reviewed the [NAMESPACE](D:/aghado01/packages/BASS/NAMESPACE:3), all public help entries, the 19 R source files, test scripts, examples, [NEWS](D:/aghado01/packages/BASS/NEWS.md:1), and relevant parts of the bundled [JSS paper](D:/aghado01/packages/BASS/vignettes/jss_paper.pdf). Public behavior below follows executable source when documentation differs.
- Compared with [the adaptation notes](D:/aghado01/ThermoMapper/issues/gemini-bass-adaptation-notes.md) and [BASS 2.0 notes](D:/aghado01/ThermoMapper/issues/gemini-bass-2.0.md).
- Ran five small, isolated R helper checks, described below. These did not fit models, install packages, or modify BASS. No posterior validation, performance benchmark, full R suite, or C# tests were run for this report.

BASS's BMARS ancestry dates to 1998, but the package paper is from 2020 and this checkout is dated 2023. Historical lineage does not establish why a particular compromise was made. There is direct evidence for some packaging compromises: NEWS records replacing GSL with `hypergeo` to remove a system dependency, replacing a BiasedUrn calculation, removing a deprecated wavelet dependency, and removing JSON model import/export. Those are stronger explanations than attributing every limitation to old computers. [Publication](https://www.jstatsoft.org/article/view/v094i08), [local history](D:/aghado01/packages/BASS/NEWS.md:16).

## Complete public entry-point inventory

The namespace exposes **six named functions and nine S3 methods**. The IDs below are coverage items for the full port; they are not proposed C# type names.

| ID | Public entry point | Capability and actual contract |
|---|---|---|
| P01 | [`bass`](D:/aghado01/packages/BASS/R/bass.R:54) | Fit a Bayesian adaptive spline model with scalar or directly indexed functional output. Continuous, categorical, and mixed predictors; birth/death/change sampling; hierarchical priors; optional tempering; retained posterior draws and optional training predictions. |
| P02 | [`bassPCA`](D:/aghado01/packages/BASS/R/bassPCA.R:36) | Construct a response PCA decomposition, or accept an already prepared `dat`, then fit one BASS model per retained response coefficient. Returns `bassBasis`, not a separate prediction class. |
| P03 | [`bassBasis`](D:/aghado01/packages/BASS/R/bassPCA.R:126) | Fit independent coefficient regressions under a supplied output basis. The basis need not be orthogonal; independence of coefficient models is an explicit modeling assumption. Optional parallel execution across components. |
| P04 | [`sobol`](D:/aghado01/packages/BASS/R/sobol_normTmix.R:41) | Draw-wise analytic Sobol decomposition for `bass`, including mixed inputs and direct functional responses. Global scalar indices or indices varying with one selected response coordinate; optional main/two-way effect functions for continuous scalar models. Some advertised options have implementation gaps below. |
| P05 | [`sobolBasis`](D:/aghado01/packages/BASS/R/bassPCA.R:307) | Sensitivity of reconstructed basis/PCA output, including cross-component terms. Select interaction order and largest indices; retain an `other` variance component. Executable path uses **one retained draw index**, defaulting to the last, despite plural wording elsewhere in its help. |
| P06 | [`calibrate.bassBasis`](D:/aghado01/packages/BASS/R/calib.R:26) | Infer bounded numeric simulator parameters from an observed output vector using a `bassBasis` emulator. Includes measurement uncertainty, truncation uncertainty, emulator residual uncertainty, discrepancy, adaptive proposals, decorrelation proposals, and tempering. Two covariance formulations. |
| P07 | [`predict.bass`](D:/aghado01/packages/BASS/R/predict.R:34) | Evaluate selected posterior draws at new predictor values and, for direct functional fits, new output-coordinate locations. Defaults to latent mean draws; `nugget=TRUE` adds residual noise. Actual scalar layout is draw × case, contrary to the help's transposed description. R drops singleton dimensions. |
| P08 | [`predict.bassBasis`](D:/aghado01/packages/BASS/R/bassPCA.R:183) | Reconstruct output from component predictions, restore response scaling/centering, optionally add empirical truncation error. Returns draw × case × output position. Defaults to `nugget=TRUE`, unlike P07. Uses the supplied basis locations; has no `newdata.func` equivalent. |
| P09 | [`plot.bass`](D:/aghado01/packages/BASS/R/plot.R:25) | Basis-count and residual-variance traces; training predicted-versus-observed plot, optional intervals, and residual histogram with Gaussian reference. Prediction panels depend on saved training predictions. |
| P10 | [`plot.bassBasis`](D:/aghado01/packages/BASS/R/plot.R:80) | Component count/variance traces; reconstructed training predictions and residual histogram. `pred=FALSE` suppresses prediction panels. |
| P11 | [`plot.bassSob`](D:/aghado01/packages/BASS/R/plot.R:126) | Scalar index/total-index boxplots, or cumulative functional sensitivity and variance curves. Used by both sensitivity entry points. |
| P12 | [`print.bass`](D:/aghado01/packages/BASS/R/generics.R:15) | Call, predictor/category counts, sample size, and functional-coordinate/grid sizes. |
| P13 | [`print.bassBasis`](D:/aghado01/packages/BASS/R/generics.R:47) | Predictor/category counts, sample size, and number of output-basis components. |
| P14 | [`summary.bass`](D:/aghado01/packages/BASS/R/generics.R:97) | Basic fit metadata plus range of retained term counts and posterior mean residual variance. Not a coefficient-significance or convergence report. |
| P15 | [`summary.bassBasis`](D:/aghado01/packages/BASS/R/generics.R:72) | Essentially the same basic metadata as its print method. |

Plotting parity means preserving the useful diagnostics and their underlying data. In ThermoMapper, regression should expose those data; rendering belongs above `maths`, using the existing visualization layer when that integration is undertaken.

### Fitting controls that must remain representable

These are the full formal controls of `bass`, grouped by purpose. Defaults are from source, not assumed from the help.

| Controls | Behavior/default | Port implication |
|---|---|---|
| `xx`, `y`, `xx.func` | Predictor matrix/data frame; scalar response or case × output-location matrix; coordinates indexing output locations. Factors identify categorical inputs. Missing `xx`/`y` values are rejected. | Typed predictor schema, explicit response shape and coordinate metadata. Validate all nonfinite values and coordinate data too. Avoid guessing orientation by transposing mismatched inputs. |
| `maxInt`, `maxInt.cat`, `maxInt.func` | Each defaults to 3, capped by its group's variable count. These are **separate group caps**, not one total interaction cap. In mixed models a group may contribute zero factors. | Preserve these model-space controls. A single hard-coded two-dimensional basis or total-order-three restriction would lose functionality. |
| `degree` | Defaults to 1; powered hinge factors permit higher polynomial degrees. Validation is commented out. | Support the intended positive integer degrees with explicit validation and qualified sensitivity support. |
| `maxBasis` | Defaults to 1000; clamped to number of scalar response entries. Excludes intercept. | Capacity limit and statistical term-count prior are distinct. Full-rank designs actually require `M+1 <= N`, including intercept. |
| `npart`, `npart.func` | Defaults to `min(20, 0.1 × group sample count)`; reject candidate factors supported on too few observations/locations. Continuous/category products share the case-side support check. | Preserve support constraints and include their precise definition in the target/proposal specification. |
| `nmcmc`, `nburn`, `thin` | 10000, 9000, 1. | Record actual retained count; define edge cases and divisibility. Do not size results from an unchecked fractional quotient. |
| `g1`, `g2`, `s2.lower` | Inverse-gamma shape/scale for residual variance: 0, 0, lower bound 0. Default is an improper prior kernel. | Stable gamma/inverse-gamma primitives and truncation; posterior propriety and heated-target validity require checking. |
| `h1`, `h2` | Gamma shape/rate for Poisson term-count intensity, both 10. | Preserve the complexity hierarchy, with explicit normalization at the term-count cap. |
| `a.tau`, `b.tau`, `beta.prior` | Coefficient-precision hyperprior: shape 0.5, default rate **2/N**; coefficient options `g` and `jeffreys`. Help says 1/N for the rate; paper and code use 2/N. | Implement distinct, derived targets for both advertised coefficient-prior modes. Validate option names. |
| `w1`, `w2` | Both 5; nominal interaction-order and variable proposal weights, added to current-model counts. | These tune proposals; they are not sparse Dirichlet model priors. |
| `temp.ladder`, `start.temper` | Default single temperature 1; default swap start `min(1000, ceiling(nburn/2))` when multiple replicas. R takes temperatures; ThermoMapper takes inverse temperatures. | Explicit ladder convention and validation; derive hot transitions and swaps from one target definition. |
| `curr.list` | Resume from per-replica working states with the same setup. | Starting from saved state is supported; exact replay is not established. Add data/config identity, iteration/adaptation state, replica mapping, and RNG state to a real checkpoint. |
| `save.yhat`, `small` | Defaults TRUE/FALSE. First controls saved training draws and their mean; second omits restart states and some working basis matrices. | Separate portable posterior storage, optional predictions, diagnostics, and restart state. `small` does not remove the main preallocations during fitting. |
| `verbose`, `ret.str` | Progress printing; optional return of internal data/prior/function structures. | Structured progress and inspectable model metadata. Live R closures are not a portable result contract. |

Sources: [validation and input handling](D:/aghado01/packages/BASS/R/bass.R:60), [priors](D:/aghado01/packages/BASS/R/bass.R:257), [storage](D:/aghado01/packages/BASS/R/bass.R:379), [returned object](D:/aghado01/packages/BASS/R/bass.R:570).

### Three response workflows

Let there be `n` cases, `p` predictors, `m` response locations/components, and `d` coordinates describing a response location.

| Workflow | Representation | What it permits | Assumption/limit to preserve or explicitly revise |
|---|---|---|---|
| Scalar | `X: n × p`, `y: n` | A surface over any supported number of numeric/categorical predictors. | Gaussian, independent, common-variance residuals in the fitting model. |
| Direct functional | `X: n × p`, `T: m × d`, `Y: n × m` | Learn `f(x,t)` jointly; predict at new `x` and new `t`. Each term factors into a predictor part and an output-coordinate part. | All cases share the same complete list of response locations. Those locations need not be a Cartesian lattice. Missing/ragged observations require another likelihood/design path. |
| PCA/supplied output basis | `K: m × r`, coefficient responses `r × n`, plus output offsets/scales | Learn independent coefficient surfaces, reconstruct arbitrary vector/curve/field outputs. A supplied basis can be nonorthogonal. | The basis is fixed during regression; its estimation uncertainty is not integrated out. Prediction at new output locations requires a separately supplied basis evaluator, absent from the R API. |

`bassPCA` adds `n.pc` or `perc.var=99`, `center=TRUE`, `scale=FALSE`, optional prepared `dat`, and `n.cores`/`parType`. Center/scale can also be supplied vectors. `bassBasis` expects the prepared basis/coefficient decomposition; it does not solve for coefficients under an arbitrary basis for the caller. Its `...` forwards fitting controls to each BASS component. See [setup](D:/aghado01/packages/BASS/R/bassPCA.R:46) and [supplied-basis contract](D:/aghado01/packages/BASS/R/bassPCA.R:109).

Prepared `dat` includes `xx`, `n.pc`, `basis`, `newy`, `y.m`, `y.s`, and optional truncation residuals; training plots also need `y`. The source stores truncation residuals as **output location × training case**, although help text says the reverse.

### Sensitivity, calibration, and operational coverage

**Sensitivity.** `sobol` accepts retained draw indices, input measures (`prior`, `prior.func`), an optional conditioning coordinate `func.var`, its evaluation grid `xx.func.var`, and `getEffects`. Its algebra assumes an independent product measure over inputs. Uniform and truncated normal/Student mixtures are represented; default continuous ranges come from training data. Mixture truncation renormalizes the whole mixture, so the specification must distinguish weights before truncation from weights conditional on truncation. The Student `sd` argument is used as a location-scale parameter, not converted from a marginal standard deviation.

For a scalar decomposition it returns main/interactions `S`, totals `T`, and additional variance/mean/effect data. These decompose the fitted latent surface under the input measure; they do not add observation noise to the variance denominator. With direct functional output, either integrate over output coordinates as additional inputs or condition on **one** coordinate to obtain varying indices. The functional result contains `S`, `S.var`, labels, and grid; it does not expose the scalar path's `T` result. `getEffects` is disabled for categorical or direct functional fits. Its implemented grids are 100 points per axis, with main and two-way effects only; degree and nonuniform-measure limitations are discussed below.

`sobolBasis` takes required `int.order`, optional `nind`, one draw index, input measures, parallel controls, and `plot`. It returns output-location-dependent `S`, `S.var`, `Var.tot`, labels including `other`, and an artificial 0–1 plotting coordinate. It includes covariance between different fitted coefficient surfaces under the input measure; independent coefficient fitting does **not** let us discard these cross terms. Its current public API does not compute a full posterior distribution of sensitivity indices in one call. [Direct sensitivity](D:/aghado01/packages/BASS/R/sobol_normTmix.R:194), [basis sensitivity](D:/aghado01/packages/BASS/R/bassPCA.R:418).

**Calibration.** All of `y`, `mod`, `type`, `sd.est`, `s2.df`, `s2.ind`, `meas.error.cor`, `bounds`, `discrep.mean`, and `discrep.mat` are formal inputs. Optional controls are `nmcmc=10000`, `temperature.ladder=1.05^(0:30)`, `decor.step.every=100`, and verbosity. It calibrates one observed output vector, with all bounded numeric emulator inputs treated as unknown. There is no general categorical-parameter proposal, fixed-design/unknown-parameter partition, or repeated-experiment interface.

| Calibration formulation | Intended covariance/discrepancy behavior |
|---|---|
| Type 1, dense functional output | Diagonal empirical truncation covariance; grouped independent measurement-error variances; low-rank emulator residual covariance; discrepancy represented by supplied basis `D` with Gaussian coefficient updates. Measurement correlation and fixed discrepancy mean are not used in this branch. |
| Type 2, shorter multivariate output | Full empirical truncation covariance, supplied measurement correlation scaled by unknown variance groups, fixed full discrepancy covariance and fixed discrepancy mean, plus emulator residual covariance. |

The implementation actually chooses its main branch from `ncol(discrep.mat) < length(y)`, while `type` chooses the truncation covariance. C# should use an explicit formulation and validate matching dimensions.

Calibration returns `theta`, `s2`, `discrep.vars` for every iteration and temperature, plus acceptance/decorrelation counts, proposal scales, emulator draw indices `ii`, and final working state/data. It performs no public burn-in/thinning step. Returned `theta` remains in **normalized [0,1] coordinates**: conversion to supplied bounds happens only inside emulator evaluation. Variance priors use `a=s2.df/2`, `b=a*sd.est^2`; if **any** `s2.df` is zero, a different kernel `-log(1+s2)` is selected for **all** groups. These details need explicit contracts, rather than inherited behavior. [Calibration setup](D:/aghado01/packages/BASS/R/calib.R:42), [return](D:/aghado01/packages/BASS/R/calib.R:466).

**Diagnostics and storage.** Fit results carry coefficients, variance/intensity/precision traces, term counts, unique structure descriptions and draw-to-structure lookup, input transforms, training responses, optional fitted draws, and tempering records. They do not provide modern multi-chain convergence diagnostics. Print/summary methods are lightweight displays. R session RNG state is external to `curr.list`; there is no exported seed or serialization API. Fork/socket options apply to component fitting/prediction and basis sensitivity; the fitting tempering loop itself is sequential.

### Things that are not active public BASS capabilities

- `rjmcmc_des_ridge.R` is unreachable from the six-type dispatch in `bass`. Its filename does not establish a directional ridge-spline feature: it still constructs ordinary axis-aligned candidate hinges. Do not count it as supported geometric regression.
- `checkMono`/`checkMono2` have no active callers. There is no public monotonicity constraint.
- `predict_fast` and its basis wrapper are unexported shortcuts, not separate parity requirements. The scalar shortcut disables predictor scaling and selects the continuous-only multiplication path.
- `plot_prior`/`sample.prior` and integration S3 methods are internal helpers. Some examples use `BASS:::` to access them; that does not make them namespace exports.
- JSON import/export was removed. The removed wavelet example does not remove the ability to supply a wavelet basis through `bassBasis`.
- Classification, Poisson/count likelihoods, robust/heteroscedastic fitting, correlated scalar residual fitting, correlated-input Sobol decompositions, learned projection directions, and missing-response imputation are extensions, not existing public baseline functionality.

Sources: [dispatch](D:/aghado01/packages/BASS/R/bass.R:169), [unused ridge path](D:/aghado01/packages/BASS/R/rjmcmc_des_ridge.R:11), [monotonicity helpers](D:/aghado01/packages/BASS/R/rjmcmc_all.R:279), [fast prediction](D:/aghado01/packages/BASS/R/predict.R:117).

## How the backend works, and what to retain

### Adaptive terms and six specialized fitting paths

A model is an intercept plus a variable number of products of normalized hinge factors and categorical membership indicators:

\[
f(x,t)=a_0+\sum_{j=1}^{M} a_j B_j(x,t),\qquad
B_j=\prod_{v\in J_j}\frac{[s_{jv}(x_v-k_{jv})]_+^{q}}{c_{jv}}
       \prod_{v\in C_j}1(x_v\in A_{jv}),
\]

with analogous continuous factors for output coordinates in direct functional fits. Normalization bounds a term by one on the normalized training domain. A variable appears at most once within a term. Different terms can reuse it.

The R code specializes birth/death/change and prediction for continuous, categorical, mixed, continuous-functional, categorical-functional, and mixed-functional data. This avoids some runtime branching but duplicates algebra and acceptance logic across six files. A C# engine can share the statistical target and move accounting, while retaining specialized evaluation/design kernels for scalar and separable functional observations. This is a concrete consumer for more general design contracts; it does not require turning the existing ordered 1D `KnotConfig` into a universal object.

Retain normalized factors and compact term descriptions. Evaluate a proposed term into a private reusable buffer, accumulate support and cross-products, and commit it only after acceptance. Encode categories once with a durable dictionary. Use `ulong` for small category sets and a multiword bitset or another representation above 64 levels; a single word is not a full categorical implementation. [Basis primitives](D:/aghado01/packages/BASS/R/basis_funcs.R:13), [candidate construction](D:/aghado01/packages/BASS/R/rjmcmc_all.R:11).

### Sufficient statistics and the Gaussian target

Write `B` for the complete design, including intercept, `k=M+1`, `N` for the number of scalar response entries, `G=BᵀB`, `h=Bᵀy`, and `Q=hᵀG⁻¹h`. Under the g-prior

\[
a\mid\sigma^2,\tau,B\sim N(0,\sigma^2G^{-1}/\tau),
\]

integrating out `a` gives

\[
\log p(y\mid B,\sigma^2,\tau)
=-\frac N2\log(2\pi\sigma^2)
+\frac k2\log\frac{\tau}{1+\tau}
-\frac{y^Ty-Q/(1+\tau)}{2\sigma^2}.
\]

This is a useful independent specification for the g-prior path. It is not the existing BARS weighted-normal target with a different basis. The conditional coefficient mean is `G⁻¹h/(1+τ)` and covariance is `σ²G⁻¹/(1+τ)`.

R already changes only the necessary entries of `G` and `h` when a term changes. However, `getQf` fully factors the candidate Gram matrix for every numerically evaluated move; accepted moves also materialize a triangular inverse. Preserve the sufficient statistics, replace the repeated full factorization with updates, and sample coefficients by a triangular solve. In R, despite the name `R.inv.t`, accepted-move helpers store `R⁻¹`, which is the correct orientation for multiplying a normal vector when `G=RᵀR`. Do not infer a transpose defect from that field name. [Proposal algebra](D:/aghado01/packages/BASS/R/rjmcmc_des.R:18), [factorization](D:/aghado01/packages/BASS/R/mcmc_funcs.R:113), [accepted-state factor](D:/aghado01/packages/BASS/R/rjmcmc_all.R:93).

The special `jeffreys` option needs its own derivation, including dimension-dependent normalization. Likewise, specify whether the term-count hierarchy means a conditionally normalized truncated Poisson with an independent Gamma prior on its intensity, or a globally truncated joint prior. The former has a cap-normalizer depending on intensity, so its intensity conditional is not simply Gamma. Source uses the simple Gamma update. The distinction matters when the cap is active.

### Incremental factorization and numerical admission

For a proposed new design column `b`, with current factorization `G=LLᵀ`, use

\[
G'=\begin{bmatrix}G&u\\u^T&v\end{bmatrix},\quad
u=B^Tb,\ v=b^Tb,\quad
Lz=u,\quad \delta=v-z^Tz.
\]

The new Cholesky border is `zᵀ` and `sqrt(δ)`. Factor maintenance costs `O(k²)` instead of `O(k³)`. This is a **bordered extension**, not a same-size rank-one update. For death, delete the appropriate factor column from an upper-triangular representation and restore triangularity with row Givens rotations; a lower-triangular representation needs the transposed procedure. A changed term replaces one Gram row and column, a symmetric change of rank at most two. A staged delete/insert is a simpler initial implementation than a delicate direct update/downdate pair.

The full proposal is not `O(k²)` regardless of data size: scalar cross-products still cost `O(Nk)`, plus factor evaluation. Functional contractions have their own cost below. Benchmark factor work separately from data scans.

R rejects Cholesky failures and uses a diagonal-ratio threshold of 1000 as a conditioning heuristic. A fixed `δ <= 1e-12` replacement is not generally scale-aware. Use finite checks, relative admission criteria, diagnostics, and a full-factorization/rebuild fallback for suspicious accumulated error. Pivoted QR or SVD can diagnose rank and offer a stable reference; they cannot make the original singular g-prior mathematically defined. Adding a ridge changes the prior and must be an explicit model choice. Floating-point admission rules themselves need a documented relationship to the intended model space.

For `G=RᵀR`, precision-normal innovation is `R⁻¹z`, obtained by solving `Ru=z`; covariance-normal sampling `Lz` is a different operation. Avoid full inverse formation. Also avoid cancellation in `yᵀy-Q/(1+τ)` and unstable tail-CDF subtraction. Centering responses silently is not automatically model-preserving because the intercept also has a prior; any response transformation must carry the prior/offset semantics with it.

### Functional separability is already exploited

In a direct functional fit, let `A` contain the case-side basis factors and `H` the output-location factors. A full design column would be `A_j ⊗ H_j`, but R avoids constructing the full `nm × k` design:

\[
G=(A^TA)\circ(H^TH),\qquad h_j=A_j^TYH_j.
\]

This is the key reusable idea for large functional outputs. It reduces design storage to approximately `(n+m)k`, while the observed response still costs `nm`. A new Gram border costs approximately `(n+m)k`; contracting the new term with dense `Y` generally costs `nm`. It is incorrect to describe the whole functional proposal as independent of response-grid size. Mixed categorical/continuous case factors are multiplied before this contraction. [Implemented contractions](D:/aghado01/packages/BASS/R/rjmcmc_des_func.R:24), [mixed counterpart](D:/aghado01/packages/BASS/R/rjmcmc_des_cat_func.R:27).

### Proposal construction and mixing

Birth chooses interaction sizes, distinct variables, signs, observed-value knots, and category subsets. Death chooses a term; change refreshes a sign/knot or categorical subset of an existing factor. The local change path does not switch to another predictor. Interaction and variable weights are based on the **current active model**, with nominal counts `w1`/`w2`; death removes counts. This is state-dependent proposal construction, not adaptation from an ever-growing history.

The proposal is not wholly uniform. There is an important exception: the one-variable branch chooses its variable uniformly, whereas higher-order selections use the weights. Sorted weighted sampling without replacement requires the probability of an unordered subset. `dmwnchBass` calculates that through inclusion–exclusion over subsets, with exponential dependence on interaction order. This is a real reason to examine high-order proposal cost. Merely multiplying probabilities for one selected ordering would be wrong. Preserve exact proposal accounting, or introduce an explicit auxiliary ordering with a correctly derived target. [Selection](D:/aghado01/packages/BASS/R/rjmcmc_all.R:14), [subset probability](D:/aghado01/packages/BASS/R/mcmc_funcs.R:12).

Safe model-preserving extensions include additional variable-replacement or coordinated moves and residual-informed proposal mixtures with a full-support exploratory component. For a candidate scalar column, the gain in the fitted quadratic form is

\[
\Delta Q=\frac{(b^Tr_{\mathrm{OLS}})^2}{\delta}.
\]

This is a more relevant score than raw predictor correlations. Marginal `Xᵀr` screening can miss a pure interaction such as `y=x₁x₂` under centered independent inputs. Any informed proposal still needs its actual reverse probability and normalizers. Locally balanced proposals supply useful principles; they do not automatically make a random candidate-pool scheme valid. [Zanella](https://arxiv.org/abs/1711.07424).

### Sampling correctness must precede sampler optimization

Three source patterns need a derivation and independent checks before R fitting output becomes the posterior oracle:

1. **Boundary move selection.** At zero terms the sampler forces birth; at the cap it chooses only death/change. Acceptance code still receives fixed one-third move probabilities. The actual forward and reverse selection probabilities must appear in the acceptance ratio. [Dispatch](D:/aghado01/packages/BASS/R/mcmc_update.R:16), [fixed probabilities](D:/aghado01/packages/BASS/R/bass.R:250).
2. **Partial collapse.** Coefficients are sampled, variance is then sampled with coefficients marginalized out, and coefficient precision uses the earlier coefficient draw together with the new variance. Retained coefficients and variance also come from this ordering. Reordering collapsed and uncollapsed updates requires a stationary-distribution argument; it is not an innocuous implementation detail. [Updates](D:/aghado01/packages/BASS/R/mcmc_update.R:35), [van Dyk and Jiao](https://arxiv.org/abs/1309.3217).
3. **Tempering.** Birth/death code multiplies proposal/prior-correction terms together with the likelihood contribution by inverse temperature. Standard MH for a specified powered target leaves proposal ratios outside that power. Swap density calculation also disables basis-prior blocks. Marginalizing coefficients and tempering a joint density do not generally commute. Derive hot within-replica transitions, hyperparameter updates, and swaps from the same chosen density. [Birth ratio](D:/aghado01/packages/BASS/R/rjmcmc_des.R:34), [swap density](D:/aghado01/packages/BASS/R/mcmc_funcs.R:64).

These are source-grounded correctness concerns, not measurements of posterior bias in this report. A transparent baseline could update the collapsed structural/variance state and then regenerate coefficients in a justified order. Another option is a fully coefficient-marginal target, with precision updated by MH or slice sampling and coefficients generated for retained predictions. The latter cannot keep the existing coefficient-dependent precision Gibbs formula unchanged.

Retain tempering as a full capability, but do not assume it guarantees exploration. Independent cold-chain replicates, occupancy/round-trip diagnostics, move rejection reasons, and invariant prediction readouts are more informative than acceptance alone. Comparing coefficient index 4 across changing basis structures is not a meaningful convergence diagnostic.

### PCA and supplied-basis reconstruction

R performs an actual SVD of the transposed, centered/scaled response matrix:

\[
Y_{\mathrm{scaled}}^T=UDV^T,\qquad K=U_rD_r,\quad W=V_r^T.
\]

It models rows of `W`, not ordinary PCA scores. Reconstruction is `K W`, followed by output scaling and centering. Existing ThermoMapper PCA uses covariance/Gram eigendecomposition and reports eigenvalues scaled by `1/n`. An adapter would therefore use singular values `sqrt(n × eigenvalue)`, scale loadings by them, and divide ordinary scores by them. Passing existing scores unchanged would alter coefficient-response scale, priors, and nugget interpretation. [R decomposition](D:/aghado01/packages/BASS/R/bassPCA.R:88), [C# PCA](D:/aghado01/ThermoMapper/src/maths/linalg/PCA.cs:11).

A stable thin SVD is a strong shared-library improvement; a truncated/randomized decomposition can be an explicit option with a residual-error criterion. Covariance eigendecomposition is mathematically related but squares conditioning and is not numerically equivalent to direct SVD. Preserve explicit rank, zero-variance-column handling, component sign conventions for comparisons, and reconstruction residuals. Equal singular values make componentwise comparisons ambiguous; compare reconstructed quantities or subspaces where appropriate.

### Analytic Sobol machinery

For a product input measure, each term's mean and pairwise product moment factor into one-dimensional hinge moments and categorical subset probabilities. The code caches these quantities per unique basis structure, forms quadratic forms in posterior coefficients, and subtracts lower-order contributions to obtain ANOVA components. This is worth preserving. “Analytic” describes integration conditional on a sampled spline model, not exact posterior inference.

Direct sensitivity discovers interactions present in the fitted terms. Basis sensitivity enumerates all variable subsets through `int.order`, whose count is `sum choose(p,j)`, and loops over pairs of response components and pairs of spline terms. These can dominate runtime even after matrix-update improvements. `nind` filters results after computation; it is not a computational budget on subset enumeration. [Structure reuse](D:/aghado01/packages/BASS/R/sobol_normTmix.R:238), [subset enumeration](D:/aghado01/packages/BASS/R/bassPCA.R:429), [cross-component integration](D:/aghado01/packages/BASS/R/bassPCA.R:624).

C# should share one tested moment service between direct and basis sensitivity, cache only useful moments under a memory budget, exploit symmetry and sparse variable supports, and batch coefficient quadratic forms. Categorical first and second moments should be sums of supplied category probabilities, including intersection probabilities. Avoid ratio formulations that create `0/0` for zero-probability factors.

Uniform polynomial moments are comparatively simple. Normal/Student mixtures require stable CDF differences and truncated first/second moments; the source restricts normal/Student first moments to degree one and uses `hypergeo` for Student antiderivatives. Student formulas contain factors involving `df-1`; finite-interval moments need a well-defined treatment at such special values even when the untruncated distribution lacks a moment. A moment recurrence or bounded one-dimensional quadrature fallback, split at knots, may be more maintainable than a literal hypergeometric port. That is a numerical method choice requiring error bounds and parity tests, not permission to silently change the input distribution. [Moment implementations](D:/aghado01/packages/BASS/R/sobol_normTmix.R:934).

Define what happens at zero total variance: sensitivity is undefined, not evidence of an influential or unimportant variable. Stabilize small variance differences and distinguish rounding error from substantive negative variance. Current debug `browser()` calls are unsuitable library failure behavior.

### Calibration is a separate inference subsystem

Calibration combines a draw of emulator mean predictions with covariance contributions from coefficient-model nuggets, estimated basis truncation residuals, measurement error, and discrepancy. It chooses a new emulator posterior index uniformly each iteration, shared across temperatures; updates log measurement variances by MH; updates discrepancy coefficients by Gaussian sampling in the basis formulation; proposes a block of bounded parameters using an adapting covariance; periodically makes coordinatewise uniform independence proposals; and performs replica swaps. [Iteration](D:/aghado01/packages/BASS/R/calib.R:227).

The functional formulation already uses the Woodbury identity and determinant lemma. However, it materializes diagonal and full inverse matrices. For covariance `C=D+U Uᵀ`, keep the diagonal of `D` as a vector and factor `I+UᵀD⁻¹U`. Compute `C⁻¹v` and `log|C|` through this small factor, without an `m × m` inverse. Setup then costs approximately `O(mr²+r³)` with `O(mr+r²)` storage; applying the solve costs `O(mr+r²)`. The dense full-covariance formulation remains distinct and generally costs `O(m³)` to factor. [Existing Woodbury helpers](D:/aghado01/packages/BASS/R/calib.R:119), [functional covariance](D:/aghado01/packages/BASS/R/calib.R:183).

Important calibration-specific issues:

- After drawing discrepancy coefficients, the code changes `curr$discrep` without refreshing `curr$llik` before the parameter acceptance test. Candidate and current likelihoods can then use different discrepancy states. This is a concrete cache-consistency defect. [Discrepancy update](D:/aghado01/packages/BASS/R/calib.R:299), [parameter acceptance](D:/aghado01/packages/BASS/R/calib.R:350).
- Variance-proposal Jacobians are multiplied by inverse temperature, while swap ratios omit those terms and discrepancy-prior terms. Specify the base measure and full heated target before reusing these formulas. [Variance MH](D:/aghado01/packages/BASS/R/calib.R:277), [swaps](D:/aghado01/packages/BASS/R/calib.R:425).
- With nonunit response scaling, predictions are restored to physical units, but covariance assembly uses unscaled `basis` and truncation residual covariance. All uncertainty terms must be expressed in the same response units. [Prediction reconstruction](D:/aghado01/packages/BASS/R/bassPCA.R:224), [calibration covariance](D:/aghado01/packages/BASS/R/calib.R:151).
- Unconditionally refreshing an emulator draw and taking a few calibration transitions does not, by itself, establish exact sampling from a cut posterior. The desired target must be written down: plug-in calibration, fully joint calibration, and a modular cut posterior are different inferential objects. For a cut distribution `p(e|training) p(theta,noise,discrepancy|y,e)`, a nested conditional reference sampler is a useful baseline; approximation error from finite inner runs must be assessed.

There is newer research directly relevant to this last issue. Mathews, Gopalan, Gattiker, Smith, and **BASS author Devin Francom** develop SMC for cut-posterior computation with theoretical error control. Hutchings, Rumsey, Bingham, and Huerta study emulating conditional posteriors to reduce repeated conditional-sampling cost. These are candidates for a better calibration backend after the target and reference sampler are established, not proven speedups for this engine. [Cut SMC, 2024 revision](https://arxiv.org/abs/2406.07555), [conditional-posterior emulation, 2024](https://arxiv.org/abs/2410.19028).

## Observed contract defects and qualification gaps

“Probe” means a small helper call reproduced the behavior; “source” means directly traced code; “derivation” means a statistical correctness question still needs a full target/invariance check.

| Finding | Evidence | C# requirement |
|---|---|---|
| Empirical truncation sampler selects column indices using the **row** count. | **Probe:** a 3 × 2 residual matrix with 100 requested samples raises `subscript out of bounds`. If there are fewer rows than columns, it samples only a subset of training residuals. [Helper](D:/aghado01/packages/BASS/R/bassPCA.R:254). | Use named dimensions and sample training-case residual vectors. Preserve their within-output correlation. |
| Basis prediction recycles one component-noise draw per posterior draw across every prediction case. | **Probe:** with a zero-mean mock predictor, the maximum difference between noise columns for two cases is exactly 0. [Noise addition](D:/aghado01/packages/BASS/R/bassPCA.R:230). | Draw per case and posterior draw under the intended independent residual model; output correlation then comes through the basis. |
| Truncation residual sign and units are inconsistent with predictive reconstruction. | **Source:** residual is reconstructed minus observed standardized output, then added after restoring output units. [Residual](D:/aghado01/packages/BASS/R/bassPCA.R:100), [addition](D:/aghado01/packages/BASS/R/bassPCA.R:218). | Define residual orientation and units once; apply the output scale to sampled residuals. Sign matters for an asymmetric empirical distribution. |
| `perc.var=100` fails. | **Probe:** `bassPCAsetup` on a small nondegenerate matrix produces `NA/NaN argument`; selection uses a strict `>` comparison. Zero total variance and excessive requested rank also need handling. [Selection](D:/aghado01/packages/BASS/R/bassPCA.R:92). | Rank-aware threshold selection, explicit full-rank behavior, and constant-response handling. |
| Categorical weights advertised by `sobol` do not enter its category moment formulas. | **Source + helper probe:** subset/intersection moments use counts divided by level count. A category with probability 0.9 in a two-level example yields second moment 0.5. Basis sensitivity likewise uses uniform category moments. [Direct moments](D:/aghado01/packages/BASS/R/sobol_normTmix.R:791), [intersection](D:/aghado01/packages/BASS/R/sobol_normTmix.R:1221), [basis moments](D:/aghado01/packages/BASS/R/bassPCA.R:727). | Implement the documented weighted-category intent explicitly and validate against enumeration. Numeric weight-vector input also needs a real parser contract. |
| Basis Sobol pair moments do not respect supplied truncation bounds. | **Probe:** for positive hinges at 0.2 and 0.4 under Uniform(0.25,0.75), `C22Basis` gives 0.45 versus numerical integration approximately 0.1105903. It integrates to 0/1 or knots rather than clipping to the measure's interval. [Helper](D:/aghado01/packages/BASS/R/bassPCA.R:801). | Share interval-aware moments with the direct path. |
| Effect functions and returned mean are less general than the sensitivity API suggests. | **Source:** the direct path uses uniform expressions for `f0`/effect integrals and hard-coded degree-one basis evaluation; only main/two-way effect grids are built. [Effects](D:/aghado01/packages/BASS/R/sobol_normTmix.R:379). | Either implement effects for the declared measure/degree or reject unsupported combinations explicitly. |
| Basis sensitivity returns variances in decomposition units. | **Source:** uses `dat$basis` without `y.s`; per-location normalized indices can cancel this scaling, absolute variance cannot. Ranking indices by average variance can change under nonuniform output scales. [Assembly](D:/aghado01/packages/BASS/R/bassPCA.R:418). | Keep dimensionless indices and output-unit variance contributions distinct. |
| Prediction can reorder requested draws and erase singleton axes. | **Source:** results are appended by unique model group and then `drop`ped; arbitrary nonmonotone `mcmc.use` can lose caller order. Noise variance indexing still follows the supplied indices. [Storage](D:/aghado01/packages/BASS/R/predict.R:90). | Scatter cached calculations back to original draw positions; fixed result axes, including one draw/case/location. |
| Training and prediction use different constant-column scaling behavior. | **Source:** `scale_range` handles a zero span, `scale_range_mat` divides by it. [Training helper](D:/aghado01/packages/BASS/R/misc.R:35), [prediction helper](D:/aghado01/packages/BASS/R/predict.R:11). | One fitted transform reused everywhere; explicit unknown-category and out-of-domain policy. |
| Calibration state/coordinate/covariance inconsistencies. | **Source:** stale likelihood after discrepancy draw, normalized returned parameters, formulation selected by shape, response-scale mismatch; detailed above. | Typed formulation, consistent units, physical-coordinate results plus transform metadata, correct cache invalidation. |
| Fitting and calibration target consistency is not established by shipped tests. | **Derivation:** boundary move probabilities, partial collapse, tempering, count-prior normalization, modular calibration. | Independent algebra and statistical validation; preserve intended capabilities rather than accidental transitions. |

Other necessary edge cases include an intercept-only posterior, one response component, one predictor/case/draw, nearly dependent terms, extreme input measures, all-zero variance, rank-deficient output decompositions, and rejected proposals leaving current state unchanged. Several R paths assume dimensions survive automatic simplification or enter an interactive debugger on numerical anomalies. These are qualification requirements, not evidence that every such input has already been tested here.

The scalar diagnostic plot also forms intervals by taking quantiles of draw-specific Gaussian interval endpoints. That is not generally the quantile of the posterior predictive mixture. Preserve the diagnostic purpose, but calculate interval summaries from a specified predictive distribution. [Plot interval code](D:/aghado01/packages/BASS/R/plot.R:41).

## C# reuse and renovation map

The existing structure is compatible with BASS as another engine under `Maths.Regression`, with BASS-specific code naturally beside the other spline engines. The [assembly project](D:/aghado01/ThermoMapper/projects/Maths.Regression/Maths.Regression.csproj:8) already includes that subtree and references Samplers, Distributions, LinAlg, and Rng. There is no need for another assembly solely to introduce this engine.

| Existing area | Useful reuse | Required change or limit |
|---|---|---|
| [RJMCMC proposal/chain](D:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/ReversibleJumpChain.cs:105) | Generic state/target/move separation; proposal/Jacobian/selection accounting; move telemetry. | BASS supplies its own mixed-term state and moves. Refresh cached target values when conditional updates change hyperparameters. Keep candidate factor workspaces private until acceptance. |
| [Ensemble `IChain`](D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/IChain.cs:1) | A hybrid BASS sweep can present as a chain; ensemble execution and readout infrastructure remain useful. | Samples have changing term structures. Readouts should include invariant predictions, count, noise, variable inclusion, and sensitivity summaries. |
| [Parallel tempering](D:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/ParallelTempering.cs:18) | Ladder, independent RNG streams, swap infrastructure. | Currently constructs concrete RJ chains with shared moves/target/start and advances them sequentially. Hybrid sweeps need a replica factory/state-exchange contract and isolated mutable workspaces. Validate finite ordered inverse temperatures with cold value 1. |
| [Chain diagnostics](D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainDiagnostics.cs:17) | Existing trace/summary infrastructure. | Add rank-normalized split/folded R-hat and suitable bulk/tail ESS where traces are available. The [early-stop path](D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainEnsemble.cs:131) maps nonfinite R-hat to zero; undefined diagnostics must not mean convergence. |
| [MatrixOps](D:/aghado01/ThermoMapper/src/maths/linalg/MatrixOps.cs:37) | SIMD dot products, column layout, Gram calculations. | Add BASS factor evaluation and functional contractions without materializing the full tensor design. Prefer hardware dispatch with scalar fallback, not mandatory AVX2. |
| [Dense Cholesky](D:/aghado01/ThermoMapper/src/maths/linalg/Cholesky.cs:58) | Dense solves and log determinant are relevant. | Silent pivot flooring and eager inverse preparation are unsuitable as the strict posterior factor contract. Introduce explicit factor status, solves into buffers, bordered insertion/deletion, rebuild diagnostics, and precision-normal sampling. Existing covariance callers may still need explicit regularization/inverses; do not silently change their semantics. |
| [Band Cholesky/LDLT](D:/aghado01/ThermoMapper/src/maths/linalg/BandCholesky.cs:1) | Existing precision-sampling and factorization patterns; useful reference calculations. | BASS's active-term Gram is generally dense. Local hinge support does not imply a narrow band. Do not force BASS into a B-spline banded layout. |
| [PCA](D:/aghado01/ThermoMapper/src/maths/linalg/PCA.cs:11) | Projection/result conventions and small-case reference calculations. | Response scaling adapter plus stable SVD/rank handling. Existing covariance/Gram eigensolver is not a direct SVD or a partial-spectrum implementation. |
| [Distributions](D:/aghado01/ThermoMapper/src/maths/distributions/Gamma.cs:15), [RNG](D:/aghado01/ThermoMapper/src/maths/rng/Xoshiro256PlusPlus.cs:105) | Normal/Gamma/Beta sampling, log-gamma, seeded streams, export/restore of RNG state. | R Gamma parameters use rate in relevant calls; C# Gamma uses scale. Add stable inverse-gamma/truncation and input-measure CDF/moment functionality. No complete Student/truncated-mixture moment service exists here. |
| [BARS weighted-normal model](D:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WeightedNormalModel.cs:32) | Normal-equation and solve patterns. | Its coefficient prior/target is different. Reusing the likelihood class unchanged would not implement BASS. |
| [Spline `IBasis`](D:/aghado01/ThermoMapper/src/maths/regression/spline/IBasis.cs:10), [BAPS tensor design](D:/aghado01/ThermoMapper/src/maths/regression/spline/baps/TensorProductDesign.cs:19) | Evidence for existing spline and tensor conventions. | `IBasis` is 1D/ordered-knot oriented; BAPS tensor design is bivariate fixed-grid spline machinery. Neither is a multivariate adaptive hinge term set. |
| [BGP](D:/aghado01/ThermoMapper/src/maths/regression/bgp/GpRegression.cs:36), [BAPS sampler](D:/aghado01/ThermoMapper/src/maths/regression/spline/baps/BapsSampler.cs:79) | Existing multivariate predictor handling and conditional sampling patterns. | BGP already accepts multiple predictors; the whole regression engine does not need a blanket “1D to 2D” conversion. BAPS's current sampler expects its own difference penalty. |
| [R oracle bridge](D:/aghado01/ThermoMapper/tests/oracle/ROracle.cs:1) | Existing process/fixture infrastructure and oracle conventions. | Add pinned BASS deterministic fixtures and independent statistical oracles. No BASS oracle is present in the inspected oracle set. |

Use owned immutable posterior records and separately owned per-chain numeric workspaces. Pooling is useful for scratch buffers but dangerous if retained samples still point into reused arrays. A checkpoint must contain actual RNG and sampler state, not just an initial seed. Parallelism across response components or independent cold chains is straightforward; avoid nesting component × replica × BLAS parallelism without a shared resource budget.

The R package's large dense operations already call compiled linear algebra. A language switch alone is not a speedup argument. Measure allocation, factor work, data scans, moment integration, and effective posterior samples per second separately. [R linear-algebra documentation](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Linear-algebra). For diagnostics, nonfinite R-hat should remain unavailable rather than favorable. [Reference behavior](https://mc-stan.org/posterior/reference/rhat.html).

## Reconciliation with the modernization notes

| Proposed idea | Assessment for the full port |
|---|---|
| Bordered Cholesky, Givens deletion, speculative buffers | Adopt as core engineering work after an independently checked full-factorization baseline. Include drift monitoring and rebuilds. They improve factor maintenance; data-dependent proposal cost remains. |
| Branchless SIMD hinges | Useful implementation candidate. The actual R positive-part helper already uses arithmetic `(abs(x)+x)/2`; the notes' scalar branch and 50% misprediction story does not describe this source. Measure vectorization/layout effects without an invented branch-rate premise. |
| Category bitmasks | Good fast path with a scalable representation above 64 categories and explicit label mapping. Preserve arbitrary supported cardinality. |
| Column-major streaming and low allocation | Good fit for repeated term/cross-product calculations, but R arrays are already column-major. Benefit comes from ownership, reuse, fused kernels, and avoiding copies—not discovering column-major storage. |
| Residual-informed / locally balanced proposals | Worth designing as optional kernels under the same target. Exact forward/reverse probabilities and exploratory support are mandatory; screening must not eliminate interaction-only signals. |
| Sparse Dirichlet variable prior | A distinct prior extension, motivated by sparse Bayesian tree work. Current BASS proposal weights do not implement it. Without-replacement variable sets introduce normalizers, so tree split-count conjugacy cannot simply be copied. Positive Dirichlet weights are not exact spike-and-slab selection. [Linero](https://www.tandfonline.com/doi/full/10.1080/01621459.2016.1264957). |
| Intrinsic-dimension gating of interaction order | Not a sound universal restriction. A one-dimensional circle embedded as `(cos θ, sin θ)` can support a response involving the product of both coordinates. Use dimension estimates to inform proposals or optional representations, not silently remove valid interactions. |
| Projection/Grassmann/Stiefel or directional ridge terms | Interesting extensions after full baseline parity. A span alone does not determine axes for axis-aligned hinges; global rotation changes every design column and may destroy product-measure Sobol factorization. The existing [subspace annealer](D:/aghado01/ThermoMapper/src/maths/geometry/dim-reduction/SubspaceAnnealer.cs:38) is an optimizer, not automatically a posterior transition. |
| Prefix-sum candidate scoring | Useful for selected one-dimensional knot sweeps. Residuals and interaction contexts change, so prefix statistics and proposal normalization may need rebuilding. This is not automatically an O(1) proposal mechanism. |
| Fixed speedups, acceptance rates, reuse percentages | Not established by the inspected code or tests. Replace them with benchmark questions stratified by case count, predictor count, active terms, interaction order, response size/rank, and category cardinality. |

This separates implementation improvements that preserve a specified model from new priors, new basis families, and approximate inference. All can be useful; their validation obligations differ.

## Validation and implementation sequence for complete coverage

The shipped [test-f1 through test-f5 tests](D:/aghado01/packages/BASS/tests/testthat/test-f1.R:1) read saved fitted models and compare prediction/Sobol output to saved values at `1e-12`. They cover mixed categorical/continuous inputs, continuous examples, direct functional output, PCA output, and several input measures. The PCA prediction test explicitly disables nugget and truncation error. They do not exercise fresh fitting, prove transition invariance, or test calibration. Example-generation scripts are not a substitute for those tests.

| Stage | Work and completion evidence | Coverage |
|---|---|---|
| 1. Specify targets and contracts | Write the g-prior and Jeffreys targets, count/basis prior measures, support/rank policy, result dimensions/units, checkpoint contract, and calibration target. Pin/export source fixtures early. Check tiny enumerable model spaces and analytic Gaussian conditionals. | Statistical foundation for all P01–P15. |
| 2. Shared numerics | Strict full factor/solve reference; then bordered insert/delete/change with randomized sequence comparisons to full rebuilds, rejection immutability, near-collinearity tests, precision-draw covariance tests, stable distributions, low-rank covariance operators. | Fitting, calibration, uncertainty. |
| 3. Full mixed predictor and direct-functional fitting | One coherent target/move layer with all six data modes, all advertised priors/controls, direct functional contractions checked against explicit tiny tensor designs. Statistical comparisons across independent chains, not identical R/C# random trajectories. | P01 plus fit-side parts of P07/P09/P12/P14. |
| 4. Hybrid tempering, storage, and diagnostics | Derive hot targets; verify boundary proposal ratios and swaps; isolate replica workspaces; snapshot/restore RNG and state; retain cold-chain structures/draws; invariant diagnostics and reproducibility across execution schedules where promised. | Complete P01 operation and reusable sampler improvements. |
| 5. Supplied basis, PCA, and prediction | Validate reconstruction/rank/scaling, arbitrary fixed bases, prediction shapes/order, new direct-functional coordinates, independent case noise, and empirical truncation residuals. Test component parallelism and one-component cases. | P02, P03, P07, P08, P10, P13, P15. |
| 6. Sensitivity and effect reconstruction | Shared moments tested against quadrature and categorical enumeration; independent-input variance identities; all direct/PCA modes; restricted ranges and mixtures; zero variance; cross-component terms; requested draw ordering. | P04, P05, P11 and intended documented features currently defective in R. |
| 7. Calibration and diagnostic presentation | Both covariance formulations against dense references; cache/unit/coordinate consistency; declared modular target against a conditional reference; adaptive/tempered inference checks; reusable summary and plot data. | P06 and remaining P09–P15 presentation coverage. |

All seven stages belong to the full port. Ordering them does not redefine the goal as a scalar-only or two-predictor subset. Begin correctness fixtures with the first stage rather than postponing oracle construction until the implementation is finished.

For performance, report both kernel timings and end-to-end effective samples per second, with allocation/peak-memory figures. For statistical validation, compare predictive distributions, term/variable inclusion, variance/precision marginals, sensitivity uncertainty, and calibration recovery on controlled examples. The R implementation is valuable evidence and a deterministic formula oracle where verified; its accidental bugs and unqualified stochastic transitions should not become C# acceptance criteria.

## Backend source map

| Source group | Responsibility |
|---|---|
| `bass.R` | Input preparation, transforms, priors, six-way dispatch, initialization, replica loop, retained storage, result assembly. |
| `basis_funcs.R` | Positive part, normalization, continuous and categorical term evaluation. |
| `rjmcmc_all.R` | Candidate construction, weighted selection, term insertion/deletion/replacement, working-state updates; unused monotonicity helpers. |
| `rjmcmc_des.R`, `rjmcmc_cat.R`, `rjmcmc_des_cat.R`, `rjmcmc_des_func.R`, `rjmcmc_cat_func.R`, `rjmcmc_des_cat_func.R` | Six active birth/death/change implementations. |
| `rjmcmc_des_ridge.R` | Unselected alternative/dead path; not a public ridge-direction feature. |
| `mcmc_update.R`, `mcmc_funcs.R`, `misc.R` | Gibbs/structural sweep, proposal probabilities, log density, factorization, tempered distribution helpers and fitted-value helpers. |
| `predict.R` | Input transforms, structure-grouped prediction, six evaluation paths, unexported shortcut. |
| `bassPCA.R` | Response decomposition, supplied-basis fits, reconstruction/noise, basis sensitivity and cross-component moments. |
| `sobol_normTmix.R` | Direct sensitivity, optional effects, one-dimensional mixture moments and categorical moments. |
| `calib.R` | Active calibration function through line 471; the remaining roughly 800 lines are commented historical implementations. |
| `generics.R`, `plot.R` | Public printing, summaries, and diagnostic plotting. |

The next design decision should be the explicit statistical targets and data/result contracts, informed by this complete coverage list. Those determine which shared-library renovations are necessary and prevent a fast implementation from silently becoming a different model.
