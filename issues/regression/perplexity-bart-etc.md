BART—Bayesian additive regression trees—is the general-purpose Bayesian regression model underneath Bayesian causal forests. It represents an unknown function as a sum of many small trees, puts shrinkage priors on those trees, and samples a posterior over the entire ensemble. Chipman, George, and McCulloch introduced it in 2010. [arxiv](https://arxiv.org/html/0806.3286v2)

The key distinction from BCF: BART supplies a flexible model for a function; BCF uses that machinery in a particular parameterization designed for causal-effect estimation. [projecteuclid](https://projecteuclid.org/journals/bayesian-analysis/volume-15/issue-3/Bayesian-Regression-Tree-Models-for-Causal-Inference--Regularization-Confounding/10.1214/19-BA1195.pdf)

## 1. The mathematical model

For continuous outcomes, standard BART assumes

\[
Y_i=f(X_i)+\epsilon_i,
\qquad
\epsilon_i\overset{\mathrm{iid}}{\sim}\mathcal N(0,\sigma^2),
\]

with

\[
f(x)=\sum_{j=1}^{m}g(x;T_j,M_j).
\]

Here:

- \(T_j\) is a binary regression-tree structure, including its split variables and thresholds.
- \(M_j=\{\mu_{j1},\ldots,\mu_{jb_j}\}\) contains its leaf values.
- \(g(x;T_j,M_j)\) routes \(x\) to a leaf and returns that leaf’s value.

The model therefore learns both a partition-based representation and the coefficients associated with it. [arxiv](https://arxiv.org/html/0806.3286v2)

A useful mathematical interpretation is an adaptive indicator expansion:

\[
f(x)=
\sum_{j=1}^{m}\sum_{\ell=1}^{b_j}
\mu_{j\ell}\,
\mathbf 1\{x\in R_{j\ell}\}.
\]

Unlike a fixed spline basis, the regions \(R_{j\ell}\) are themselves unknown. Posterior inference explores alternative bases as well as alternative coefficients. The original paper describes this as regression using “dimensionally adaptive random basis elements.” [rob-mcculloch](https://www.rob-mcculloch.org/2021_dt/webpage/BART.pdf)

## 2. The prior makes the ensemble work

Without constraints, a sum of trees is massively redundant: many ensembles can represent the same function. BART makes this useful through **shrinkage** rather than trying to identify a unique tree decomposition. [arxiv](https://arxiv.org/html/0806.3286v2)

### Small trees

The usual tree prior gives a node at depth \(d\) probability

\[
P(\text{split at depth }d)
=
\alpha(1+d)^{-\beta}.
\]

Common defaults are \(\alpha=0.95\) and \(\beta=2\), making deep trees increasingly unlikely. Small trees favor low-order interactions, while the ensemble can still express complicated functions. [gatsby.ucl.ac](https://www.gatsby.ucl.ac.uk/~balaji/pgbart_aistats15.pdf)

### Small leaf contributions

Leaf values receive centered Gaussian priors:

\[
\mu_{j\ell}\sim\mathcal N(0,\sigma_\mu^2).
\]

The prior scale decreases with the number of trees, conventionally as

\[
\sigma_\mu\propto \frac{1}{k\sqrt m}.
\]

Thus, each tree contributes only a small amount, while the total prior scale remains roughly stable as \(m\) changes. In the original response-scaling convention, \(\sigma_\mu=0.5/(k\sqrt m)\), with \(k=2\) a recommended default. [arxiv](https://arxiv.org/html/0806.3286v2)

That is the central construction: individually weak trees, collectively flexible regression.

## 3. How fitting works: Bayesian backfitting

BART is not fitted by greedily appending trees until the training loss stops improving. It usually fixes the number of trees and repeatedly updates them through MCMC. [arxiv](https://arxiv.org/html/0806.3286v2)

For tree \(j\), form the partial residuals

\[
R_{ij}
=
Y_i-\sum_{k\ne j}g(X_i;T_k,M_k).
\]

Then:

1. Update tree \(j\)’s structure conditional on these residuals and the current noise variance.
2. Sample its leaf values.
3. Move to the next tree.
4. Update the noise variance and repeat.

Typical tree proposals grow, prune, change, or swap parts of a tree. Gaussian leaf priors allow leaf parameters to be integrated out when evaluating structure proposals, making the updates more tractable. [gatsby.ucl.ac](https://www.gatsby.ucl.ac.uk/~balaji/pgbart_aistats15.pdf)

After burn-in, each retained ensemble gives a function draw

\[
f^{(s)}(x)=\sum_j g(x;T_j^{(s)},M_j^{(s)}).
\]

You can average these for a point estimate or retain them for posterior inference. The individual trees are auxiliary components; the principal object of inference is their sum. [arxiv](https://arxiv.org/html/0806.3286v2)

## 4. Relationship to boosting

BART resembles boosting because both build predictions by adding weak trees. The important difference is how they construct and regularize the ensemble. [arxiv](https://arxiv.org/html/0806.3286v2)

| Aspect | BART | Classical gradient boosting |
|---|---|---|
| Fitting principle | Posterior sampling under a likelihood and priors | Sequential optimization of a loss |
| Tree updates | Revisit trees conditional on the rest of the ensemble | Typically append trees to the current fit |
| Weak-learner control | Priors on tree structure and leaf magnitudes | Tree restrictions and shrinkage/learning rate |
| Main inferential output | Posterior draws of the regression function | A fitted prediction function |

The original BART paper explicitly makes this comparison: both use residual-based fitting, but BART uses Bayesian backfitting over a fixed ensemble rather than a forward boosting path. [arxiv](https://arxiv.org/html/0806.3286v2)

## 5. What uncertainty does it provide?

BART distinguishes uncertainty about the regression function from uncertainty about a future observation:

\[
\underbrace{f^{(s)}(x_*)}_{\text{conditional-mean draw}}
\qquad\text{versus}\qquad
\underbrace{Y_*^{(s)}
=
f^{(s)}(x_*)+\epsilon_*^{(s)}}_{\text{posterior-predictive draw}}.
\]

The second includes observation noise and therefore answers a different question. Both follow from posterior draws of the ensemble and noise variance. [arxiv](https://arxiv.org/html/0806.3286v2)

These are model-based uncertainty statements. They depend on the likelihood, priors, and adequacy of posterior computation—not just on calling the estimator Bayesian.

For binary outcomes, a classic extension uses a probit link,

\[
P(Y=1\mid x)=\Phi(f(x)),
\]

while retaining the sum-of-trees representation for the latent regression function. [cran.ms.unimelb.edu](https://cran.ms.unimelb.edu.au/web/packages/BayesTree/BayesTree.pdf)

## 6. How BART becomes BCF

You can use ordinary BART for causal estimation by fitting

\[
Y=f(X,Z)+\epsilon
\]

and, for each posterior draw, contrasting

\[
\tau^{(s)}(x)
=
f^{(s)}(x,1)-f^{(s)}(x,0).
\]

But this places the prior on the outcome surface \(f\); the prior on its treatment contrast is indirect. That is the issue BCF addresses. [projecteuclid](https://projecteuclid.org/journals/bayesian-analysis/volume-15/issue-3/Bayesian-Regression-Tree-Models-for-Causal-Inference--Regularization-Confounding/10.1214/19-BA1195.pdf)

BCF instead writes

\[
Y
=
\mu(X,\widehat e(X))
+
Z\,\tau(X)
+
\epsilon,
\]

and gives \(\mu\) and \(\tau\) separate BART priors. This allows complex baseline prognosis and comparatively simple effect heterogeneity, while propensity-score inclusion helps mitigate regularization-induced confounding. [projecteuclid](https://projecteuclid.org/journals/bayesian-analysis/volume-15/issue-3/Bayesian-Regression-Tree-Models-for-Causal-Inference--Regularization-Confounding/10.1214/19-BA1195.pdf)

So the conceptual hierarchy is:

- BART: a Bayesian prior and inference procedure for an unknown regression function.
- BCF: a causal model that uses separate BART components to regularize prognosis and treatment effects differently. [arxiv](https://arxiv.org/html/0806.3286v2)

Yes—BART looks like a natural complement to your spline work, and your existing sampler architecture is already close to what it needs. I inspected `aghado01/ThermoMapper`, including the regression layout, BARS likelihood/conditional-fit contracts, `IRjMove<TState>`, and `BlockScan<TState>`.

My recommendation is: implement BART as a sibling regression family, reuse your sampling infrastructure, and build BCF as a structured composition of regression components—not as another independent sampler. The key architectural distinction is **composition**.

## Where it fits in your engine

Your current layout already separates the right kinds of responsibilities:

- `regression/spline` contains basis representations and the BARS/BAPS/BASS model families.
- `samplers/rjmcmc` owns proposal bookkeeping, including reverse moves and move-selection probabilities.
- `samplers/blocks` owns systematic scans, while model-specific blocks retain their sufficient statistics and scratch state.
- BARS separates structure uncertainty from conditional coefficient inference through `IObservationModel` and `IConditionalFit`.

That is a good starting point for BART because its standard algorithm is a partially collapsed Metropolis-within-Gibbs sampler: update one tree’s structure with its leaf values integrated out, then sample those leaves conditional on the selected structure and the rest of the ensemble. [par.nsf](https://par.nsf.gov/servlets/purl/10181031)

I would provisionally organize the addition like this:

```text
maths/regression/
    spline/                  existing
    tree/
        TreeState
        TreePrior
        SplitRule
        LeafStatistics
        GaussianLeafConditional
        TreeStructureMoves
    bart/
        BartState
        BartPrior
        TreeBackfitBlock
        NoiseVarianceBlock
        BartReadout
    causal/
        BcfState
        BcfPrior
        PrognosticForestBlock
        TreatmentForestBlock
        TreatmentEffectReadout
```

Those are proposed responsibilities, not a claim that you need exactly this many classes.

## 1. BART complements splines at the representation level

Both approaches infer a function using an uncertain, adaptive representation:

\[
\begin{aligned}
\text{splines:}\quad
f(x)&=\sum_k \beta_k B_k(x;\xi),\\
\text{BART:}\quad
f(x)&=\sum_j\sum_{\ell}
\mu_{j\ell}\mathbf 1\{x\in R_{j\ell}(T_j)\}.
\end{aligned}
\]

For splines, the structure includes knots and basis configuration. For BART, it includes tree partitions. In both cases, structure determines a basis, and coefficients determine the function conditional on that structure. BART explicitly treats the trees as adaptive basis elements. [ar5iv.labs.arxiv](https://ar5iv.labs.arxiv.org/html/0806.3286)

That suggests sharing concepts such as:

- posterior function draws;
- prediction and function readouts;
- conditional coefficient uncertainty;
- structure-move diagnostics;
- chain execution and reproducibility.

But I would not force trees into the spline basis machinery merely because they can be written as indicator columns.

For BART’s Gaussian leaf updates, you mostly need per-leaf sufficient statistics such as

\[
n_\ell,\qquad
S_\ell=\sum_{i\in\ell}r_i,
\qquad
Q_\ell=\sum_{i\in\ell}r_i^2.
\]

A dense design-matrix route would obscure the local computations that make tree updates efficient. The standard algorithm exploits integrated leaf likelihoods rather than repeatedly fitting a generic dense regression. [par.nsf](https://par.nsf.gov/servlets/purl/10181031)

So: share statistical contracts where they actually align; preserve representation-specific numerical machinery.

## 2. Your block-scan abstraction is the strongest immediate match

Your `BlockScan<TState>` already says that the model owns its blocks, sufficient statistics, factors, and scratch space. That is almost exactly the division I would want for Bayesian backfitting.

For tree \(j\), maintain the total fitted vector \(F\), and compute

\[
r^{(j)}=y-F+g_j.
\]

Then:

1. Propose a new tree structure using the residuals.
2. Accept or reject using the collapsed leaf likelihood, tree prior, and proposal ratio.
3. Sample leaf values under the resulting structure.
4. Update \(F\) by replacing the old tree contribution.
5. After the forest sweep, update the shared noise variance.

This is the standard BART backfitting structure. [par.nsf](https://par.nsf.gov/servlets/purl/10181031)

In your architecture, I would make one tree update—or possibly an entire forest sweep—an `IBlock<BartState>`. Which granularity is preferable depends on how much per-tree scheduling and diagnostics you want.

One concrete implementation concern: your scan snapshots copy model state. If that state contains cached predictions, leaf memberships, or mutable trees, the copy operation must make those caches independent or reconstruct them. Otherwise reproducible snapshots can quietly become aliasing bugs.

## 3. Reuse RJ bookkeeping, but do not force full RJ coefficient maps

Your `IRjMove<TState>` makes a useful separation:

- the engine handles the target ratio and normalized move-selection ratio;
- the move handles its auxiliary proposal terms and Jacobian;
- reverse moves are explicit.

Grow/prune proposals fit that bookkeeping well.

However, standard BART does not require you to propose a dimension-matching transformation of all leaf coefficients. Instead, it integrates the updated tree’s coefficients out for the structure decision, then draws them afterward. The resulting algorithm is partially collapsed. [par.nsf](https://par.nsf.gov/servlets/purl/10181031)

That gives you two possible integration paths:

| Path | My assessment |
|---|---|
| Reuse the existing RJ machinery on a collapsed tree-structure state | Reasonable if its target and state contracts make conditional updates straightforward |
| Write a specialized tree-structure MH step inside a backfitting block | Probably the simplest first implementation |
| Build full leaf-parameter split/merge bijections immediately | Unnecessary complexity for canonical Gaussian BART |

The important correctness condition is the order: select the structure against the leaf-integrated target, then refresh leaf values conditional on that structure. Do not treat the collapsed structure step as an arbitrary interchangeable update while retaining stale leaf parameters. [par.nsf](https://par.nsf.gov/servlets/purl/10181031)

## 4. Your BARS likelihood interface is related—but not directly reusable

This is the main boundary I noticed in the code.

Your current BARS `IObservationModel` defines a marginal likelihood with coefficients and scale integrated out. It also supports conditional draws and readouts given the knot configuration.

Canonical BART’s tree update integrates out one tree’s leaf coefficients while conditioning on:

- the other trees;
- the current shared noise variance;
- the resulting partial residuals. [par.nsf](https://par.nsf.gov/servlets/purl/10181031)

These are different marginalization contracts:

\[
\underbrace{
p(y\mid \xi)
}_{\text{your BARS marginal contract}}
\qquad\text{versus}\qquad
\underbrace{
p(r^{(j)}\mid T_j,\sigma^2)
}_{\text{BART tree-update marginal}}.
\]

I would not silently stretch `IObservationModel` to mean both.

Instead, my proposed first version would give trees a small conditional-leaf contract:

```text
LogIntegratedLeafLikelihood(stats, noiseVariance, leafPrior)
SampleLeaf(stats, noiseVariance, leafPrior, rng)
```

Once BART works, you can decide whether there is a genuinely useful common abstraction above this and your spline conditional fits. That is safer than generalizing the spline interface in anticipation.

## 5. BCF should sit above the forest components

BCF introduces a structured predictor:

\[
y_i
=
\mu(x_i,\widehat e_i)
+
z_i\tau(x_i)
+
\epsilon_i.
\]

The statistical distinction is separate priors for prognosis and treatment-effect heterogeneity, plus propensity-score inclusion in the prognostic component—not a fundamentally new tree sampler. [arxiv](https://arxiv.org/pdf/1706.09523v3.pdf)

I would therefore build BCF around two forest components sharing a likelihood and residual state:

```text
BcfState
    PrognosticForest
    TreatmentForest
    NoiseVariance
    FittedValues
```

There is one useful generalization worth designing for early: a component multiplier.

Write a tree contribution as

\[
a_i\,g(x_i).
\]

For Gaussian noise and a leaf prior
\(\mu_\ell\sim\mathcal N(0,s_\mu^2)\), the conditional posterior uses

\[
V_\ell^{-1}
=
s_\mu^{-2}
+
\sigma^{-2}\sum_{i\in\ell}a_i^2,
\]

\[
m_\ell
=
V_\ell\,\sigma^{-2}
\sum_{i\in\ell}a_i r_i.
\]

Here \(r_i\) excludes the tree being updated. These follow directly by completing the square in the Gaussian likelihood.

This one proposed abstraction handles:

- ordinary BART: \(a_i=1\);
- the BCF prognostic forest: \(a_i=1\);
- the canonical BCF treatment forest: \(a_i=z_i\).

It also prevents an easy mistake: controls with \(z_i=0\) do not supply direct likelihood information about treatment-tree leaf values, even if their covariates are routed through that tree.

## 6. What this opens up alongside your splines

There are two distinct opportunities.

### Complementary model families

First, keep BART and splines separate and compare posterior predictive behavior. That gives your engine different structural priors without requiring a hybrid model.

The meaningful comparison is not merely “which fits better?” but which representation supports the readout you care about. Standard BART uses hard partitions and produces piecewise-constant function draws; that is quite different from a smooth spline representation. [ar5iv.labs.arxiv](https://ar5iv.labs.arxiv.org/html/0806.3286)

For your extrema and curve-based readouts, this matters: an interface that accepts any posterior function does not imply every function family gives the same meaning to derivatives, isolated peaks, or peak counts.

### A later additive hybrid

A possible later extension is

\[
y_i=s(t_i)+h(x_i)+\epsilon_i,
\]

where a spline captures a smooth coordinate-dependent component and BART captures residual interactions.

I would treat this as a new model, not just two interchangeable fitters. If both components can explain the same variation, their decomposition becomes prior-dependent. You would need explicit centering or other identifiability conventions, and then test component recovery—not only total prediction.

Similarly, replacing a BCF forest with a spline component is an interesting design space, but it no longer reproduces canonical BCF’s prior.

## A practical implementation order

I would proceed in four stages:

1. Gaussian single-tree inference  
   Verify leaf conditionals, integrated likelihoods, and grow/prune proposal ratios against tiny enumeratable examples.

2. Gaussian BART  
   Add shrinkage priors, cached backfitting residuals, shared variance updates, posterior prediction, and snapshot/restart tests.

3. Canonical BCF  
   Add component multipliers, separate forest priors, propensity inputs, and treatment-effect readouts. Validate using simulated data with known effects and targeted selection—the setting that motivated BCF. [arxiv](https://arxiv.org/pdf/1706.09523v3.pdf)

4. Cross-family composition  
   Only then extract shared additive-component machinery or experiment with spline–tree hybrids.

The strongest reuse opportunity is your existing separation between model-owned updates and sampler-owned execution. 

Yes—there is a substantial literature, and the interesting part is that it branches along several largely independent axes: smoothness, variable selection, leaf models, observation distributions, causal structure, and inference algorithms. It is much richer than a succession of differently named BART implementations. [par.nsf](https://par.nsf.gov/servlets/purl/10181031)

For your engine, I would think of this as a family of Bayesian adaptive function representations connected by **modularity**: change the partition prior, change what lives inside a leaf, change the likelihood, or change how components are assembled.

## First: it is a family resemblance, not one straight genealogy

Your discovery path—BARS → BASS/BAPS → trees—makes conceptual sense, but it is not necessarily the historical derivation.

For example, BASS describes itself as closely related to Bayesian multivariate adaptive regression splines, or Bayesian MARS. It supports scalar and functional responses, uses RJMCMC, and can use parallel tempering. That makes it a neighboring adaptive-basis lineage rather than simply “BARS extended to surfaces.” [github](https://github.com/lanl/BASS)

On the tree side, the main reference points are Bayesian CART, BART, and then extensions of their priors, representations, and likelihoods. Reviews explicitly treat Bayesian CART and BART together. [par.nsf](https://par.nsf.gov/servlets/purl/10181031)

A useful conceptual map—not a literal historical family tree—is:

```text
Bayesian adaptive function estimation
│
├─ Adaptive spline representations
│    ├─ Free-knot curve models
│    └─ Adaptive multivariate spline surfaces
│
├─ Bayesian partition models
│    ├─ Bayesian CART
│    └─ Treed local models / treed Gaussian processes
│
└─ Additive tree ensembles
     ├─ BART
     ├─ Sparse / variable-selecting forests
     ├─ Soft and smooth forests
     ├─ Richer terminal-node models
     ├─ Non-Gaussian and distributional models
     ├─ Causal decompositions
     └─ Alternative inference algorithms
```

The last six branches are where a lot of the interesting design space lives. [par.nsf](https://par.nsf.gov/servlets/purl/10181031)

## 1. Soft BART: learn smooth functions without abandoning trees

This is probably the most immediately relevant branch for your spline interests.

Ordinary BART sends an observation deterministically left or right at each split. SBART replaces that with soft routing: an observation contributes to multiple leaves through smooth gating weights. Linero and Yang’s 2018 paper, “Bayesian Regression Tree Ensembles that Adapt to Smoothness and Sparsity,” develops this construction and provides posterior-concentration results adapting to unknown smoothness and sparsity under its assumptions. [academic.oup](https://academic.oup.com/jrsssb/article/80/5/1087/7048381)

Instead of

\[
f(x)=\sum_{j,\ell}\mu_{j\ell}\mathbf 1\{x\in R_{j\ell}\},
\]

you have something like

\[
f(x)=\sum_{j,\ell}\mu_{j\ell}\phi_{j\ell}(x),
\]

where \(\phi_{j\ell}(x)\) is a product of soft branch probabilities along a path. The tree now defines an adaptive smooth basis rather than a hard partition. [academic.oup](https://academic.oup.com/jrsssb/article/80/5/1087/7048381)

That is a particularly nice bridge to your spline work: both become adaptive bases with locally organized support or influence.

There is a computational price. Hard-tree leaf memberships are disjoint, whereas soft leaf weights overlap. SBART therefore needs more involved coefficient marginalization and more leaf–observation calculations; it is not merely swapping a comparison operator for a sigmoid. [academic.oup](https://academic.oup.com/jrsssb/article/80/5/1087/7048381)

Implementation to inspect: `SoftBart`, whose vignette explains the change from hard to soft regression trees. [cran.rstudio](https://cran.rstudio.com/web//packages/SoftBart/vignettes/SoftBart-Vig.pdf)

## 2. Sparse forests: put a prior on which coordinates matter

Another branch changes the split-variable prior rather than the tree’s output.

DART-style constructions use a sparse Dirichlet prior on variable-selection probabilities. Instead of assigning every predictor the same fixed probability of being chosen at a split, the model learns a shared preference over predictors. SBART also incorporates sparsity-inducing split-variable priors, using conditionally conjugate Dirichlet updates. [academic.oup](https://academic.oup.com/jrsssb/article/80/5/1087/7048381)

Schematically,

\[
s\sim\operatorname{Dirichlet}(\eta/p,\ldots,\eta/p),
\qquad
v_b\mid s\sim\operatorname{Categorical}(s),
\]

where \(v_b\) is the coordinate selected at branch \(b\).

The architectural significance is that the split prior becomes shared state across the ensemble, rather than a local detail inside each tree.

There is also a direct causal descendant: “Shrinkage Bayesian Causal Forests for Heterogeneous Treatment Effects Estimation” develops sparsity-inducing BCF. That asks which covariates matter for heterogeneous effects, not just which predict outcomes. [tandfonline](https://www.tandfonline.com/doi/abs/10.1080/10618600.2022.2067549)

For your engine, these are separate knobs:

- shrink tree depth;
- shrink leaf amplitudes;
- shrink the set of active coordinates;
- shrink treatment-effect heterogeneity.

Those are different structural preferences, not interchangeable implementations of “regularization.”

## 3. Treed Gaussian processes and GP-BART: change what lives in a region

A leaf does not have to contain a constant.

The treed-GP literature combines tree-based partitioning with local Gaussian-process modeling. Related partitioned-process work uses binary-tree structure to capture regional nonstationarity and improve computational tractability. [rdrr](https://rdrr.io/cran/hetGP/f/inst/doc/hetGP_vignette.pdf)

GP-BART takes the ensemble route: it extends BART by putting Gaussian-process priors on terminal-node predictions rather than using only scalar leaf contributions. [arxiv](https://arxiv.org/html/2204.02112v4)

These are related but distinct constructions:

| Construction | Main idea |
|---|---|
| Treed GP | Partition the domain and use local GP models |
| GP-BART | Sum trees whose terminal-node contributions have GP structure |
| Soft BART | Keep coefficient-based leaves but make routing smooth |

The attraction is that partitioning can handle regime changes while local models handle within-region structure. But local GPs and soft routing do not automatically imply the same global continuity properties; their boundaries and priors matter. [academic.oup](https://academic.oup.com/jrsssb/article/80/5/1087/7048381)

This is where your existing GP, spline, and tree components could eventually meet—but I would regard each combination as a statistical model requiring its own prior and inference design.

## 4. Distributional extensions: stop assuming all the complexity is in the mean

Canonical Gaussian BART puts its flexibility in \(f(x)\) while retaining a simple residual distribution.

DPMBART changes that by combining the tree ensemble with a Dirichlet-process mixture model for the errors:

\[
Y=f(X)+\epsilon,
\qquad
\epsilon\sim G,
\]

with a Bayesian nonparametric mixture construction for the error distribution. The paper’s motivation is precisely to relax the restrictive IID normal-error assumption while keeping BART’s flexible regression surface. [arxiv](https://arxiv.org/html/1807.00068v2)

Other work extends soft forests to quantile regression and skew-normal or skew-\(t\) observation models, including clustered survey settings. These models target different aspects of the response distribution rather than just its conditional mean. [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC13007972/)

For your architecture, this suggests separating:

1. the predictor representation;
2. the observation model;
3. any latent-variable augmentation;
4. the conditional-update machinery.

But it also reinforces the caveat from the previous answer: “supports a likelihood interface” does not mean every likelihood supports the same collapsed tree update.

## 5. Survival BART: model an entire time-to-event distribution

BART has also become a practical survival-analysis tool.

One approach models a discrete-time hazard,

\[
h(t,x)=P(T=t\mid T\ge t,X=x),
\]

using a tree-based probability model over time and covariates. Survival probabilities then follow from

\[
S(t\mid x)=\prod_{u\le t}\bigl(1-h(u,x)\bigr).
\]

The `BART-Survival` documentation describes precisely this discrete-time construction, including the augmented training and prediction datasets used to fit hazards and recover survival curves. [cdcgov.github](https://cdcgov.github.io/BART-Survival/build/html/methods.html)

This is not just “run a forest on survival times”: the observation representation handles risk sets and event status. The R `BART` package also provides models for continuous, binary, categorical, and time-to-event outcomes. [rsparapa.r-universe](https://rsparapa.r-universe.dev/BART)

For your interests in scientific pipelines, this is an example of a substantial application built from a reusable function prior plus a carefully constructed observation model.

## 6. Causal forests: change the decomposition, not only the likelihood

BCF is the most relevant example here:

\[
Y=\mu(X,\widehat e(X))+Z\tau(X)+\epsilon.
\]

The model separates prognosis from treatment response so that their priors can differ. Its motivation includes regularization-induced confounding and the problems that arise when treatment selecttion tracks baseline prognosis. [arxiv](https://arxiv.org/pdf/1706.09523v3.pdf)

[User note: interesting direction for studying Bias in machine learning and AI?]

The broader lesson is more general than treatment effects: if the scientific question concerns a particular contrast or component, place the prior on a parameterization that exposes it.

Sparse BCF then adds another layer: regularize which variables drive the effect surface. That is a composition of causal structure and sparse function learning. [tandfonline](https://www.tandfonline.com/doi/abs/10.1080/10618600.2022.2067549)

So BCF is not merely another output distribution for BART. It occupies a different axis in the design space.

## 7. Computational descendants: change how the forest is explored

There is a parallel literature on inference itself.

XBART—“Accelerated Bayesian Additive Regression Trees”—replaces the familiar small local tree edits with a grow-from-root strategy. Rather than making only incremental changes to an existing tree, it constructs new trees through stochastic split decisions. [dblp.uni-trier](https://dblp.uni-trier.de/rec/conf/aistats/HeYH19.html)

That matters for your sampler work because it separates two questions:

- What function prior and likelihood define the statistical model?
- What transition mechanism or approximation makes inference practical?

An accelerated algorithm should not automatically be treated as a drop-in exact posterior kernel merely because it returns an ensemble of sampled trees. The target and invariance properties need to be established for the particular algorithm.

## What I would explore first for ThermoMapper

Given your spline work, my research order would be:

1. SBART  
   The closest conceptual bridge from smooth adaptive bases to tree ensembles. Read Linero and Yang’s 2018 paper first. [academic.oup](https://academic.oup.com/jrsssb/article/80/5/1087/7048381)

2. Sparse split priors / DART  
   A relatively clean extension of ensemble-level prior state, and useful groundwork for sparse treatment-effect models. [academic.oup](https://academic.oup.com/jrsssb/article/80/5/1087/7048381)

3. Treed GPs and GP-BART  
   Particularly interesting alongside your existing GP regression work: region adaptation and local smoothness become separate modeling decisions. [arxiv](https://arxiv.org/html/2204.02112v4)

4. BCF and shrinkage BCF  
   Examples of scientifically meaningful component decompositions rather than just richer predictive models. [tandfonline](https://www.tandfonline.com/doi/abs/10.1080/10618600.2022.2067549)

5. DPMBART and survival BART  
   Good tests of whether your function representation and observation-model abstractions are genuinely independent. [arxiv](https://arxiv.org/html/1807.00068v2)

The unifying design question is not “how many B-prefixed regressors can I add?” It is: which aspects of an unknown function should be random—its basis, partition, smoothness, active coordinates, local model, distribution, or scientific decomposition? This literature contains concrete constructions for each of those choices, which makes it a useful source of both models and architectural counterexamples.


you have spline and GP implementations on hand; the proposed addition is the adaptive tree representation and its inference machinery

## The interesting question is the coupling

With your existing engine, there are three distinct experiments:

- A treed GP: a partition determines which observations belong to each local GP. The new work is partition inference and managing the local fits.
- GP-BART: an additive ensemble whose terminal regions carry GP contributions. The new work includes conditional residual updates and managing multiple interacting local GP states.
- An additive GP-plus-BART model:
  \[
  y_i=g_{\mathrm{GP}}(x_i)+f_{\mathrm{BART}}(x_i)+\epsilon_i.
  \]
  This keeps the two representations separate and alternates their conditional updates.

Those are proposed integration directions, not claims that your current GP API already supports them unchanged.

## What I should inspect next

The relevant question is **conditionality**: can your GP engine be fitted or updated against changing partial residuals while conditioning on a shared noise state, or does its current inference contract integrate those quantities out?

That is the same boundary we identified between your BARS marginal-likelihood interface and BART’s backfitting updates. Having both engines implemented does not yet establish that their marginalization conventions compose correctly.

[User note: It doesn't matter if the current engine has these limitation because I own the project and can develop it]

So the useful next architectural investigation is your existing GP engine’s conditioning, marginalization, and caching contracts—not a generic explanation of why GPs might be useful.

They are genuinely competitive on prediction—not merely “less accurate models that happen to give Bayesian intervals”—but the evidence does not support declaring BART the general winner over modern boosting. The strongest direct results are for accelerated descendants such as XBART against XGBoost; ordinary MCMC BART’s clearest disadvantage is computational cost. [proceedings.mlr](https://proceedings.mlr.press/v89/he19a/he19a.pdf)

For your engine, the distinctive value is **inference**: a competitive nonlinear function model that can participate in a larger Bayesian calculation, rather than simply produce a prediction.

## What the head-to-head evidence actually says

### Ordinary BART: a serious predictor, with a computational handicap

The original BART paper included a predictive “bake-off” across 42 datasets. It established BART as a competitive general-purpose prediction method, not just a framework for uncertainty quantification. But that is historical evidence: a 2010 comparison cannot establish superiority over today’s tuned XGBoost, LightGBM, or CatBoost. [asu.elsevierpure](https://asu.elsevierpure.com/en/publications/bart-bayesian-additive-regression-trees/)

The computational problem is real. Standard BART repeatedly revisits trees using local Metropolis–Hastings proposals and performs posterior sampling, rather than stopping after constructing one optimized ensemble. The XBART literature explicitly identifies much longer fitting times relative to random forests and XGBoost as a barrier to BART’s adoption. [arxiv](https://arxiv.org/pdf/2002.03375.pdf)

So the defensible assessment is:

- Predictive accuracy: competitive, rather than inherently inferior.
- Training efficiency: a major weakness of conventional BART MCMC.
- Posterior inference: an explicit objective of the algorithm, not an incidental extra. [arxiv](https://arxiv.org/pdf/2002.03375.pdf)

### XBART: direct evidence that the gap is not intrinsic to Bayesian tree regularization

The 2019 XBART paper provides a particularly informative comparison. Across its simulated data-generating processes and sample sizes, the authors reported that XBART was “31% more accurate” than cross-validated XGBoost and typically faster. Untuned default XGBoost was faster than XBART, but substantially less accurate in those experiments. [proceedings.mlr](https://proceedings.mlr.press/v89/he19a/he19a.pdf)

That result needs its boundaries attached:

- These were the paper’s chosen simulations—not a universal tabular-data ranking.
- “Faster than cross-validated XGBoost” includes the cost of tuning XGBoost.
- It does not establish superiority over equally budgeted, modern versions of all three mainstream boosting libraries.

The subsequent expanded XBART work added real-data comparisons and continued to find strong predictive performance. It also showed that initializing ordinary BART MCMC with XBART could improve credible-interval coverage and reduce total runtime. [arxiv](https://arxiv.org/pdf/2002.03375.pdf)

The substantive conclusion is that Bayesian-style shrinkage and stochastic forest search can be competitive on both speed and point prediction. You do not necessarily have to accept a large predictive penalty to obtain a Bayesian-compatible representation. [proceedings.mlr](https://proceedings.mlr.press/v89/he19a/he19a.pdf)

## There are several different competitions

| What you are evaluating | What the evidence supports |
|---|---|
| Out-of-sample point prediction | BART/XBART can compete strongly; there is no established universal winner here.  [asu.elsevierpure](https://asu.elsevierpure.com/en/publications/bart-bayesian-additive-regression-trees/) |
| Fitting speed | Conventional BART MCMC is disadvantaged; accelerated forest algorithms materially change the comparison.  [arxiv](https://arxiv.org/pdf/2002.03375.pdf) |
| Time including hyperparameter search | XBART’s defaults performed particularly well against tuned XGBoost in its published experiments.  [proceedings.mlr](https://proceedings.mlr.press/v89/he19a/he19a.pdf) |
| Posterior uncertainty about the regression function | This is native to BART’s model and inference procedure.  [asu.elsevierpure](https://asu.elsevierpure.com/en/publications/bart-bayesian-additive-regression-trees/) |
| Causal-effect estimation | BCF targets a different problem from ordinary predictive boosting through its prognosis/effect decomposition.  [arxiv](https://arxiv.org/pdf/1706.09523v3.pdf) |
| Composition inside a probabilistic model | BART components have been explicitly developed for arbitrary probabilistic modeling, rather than only standalone prediction.  [ar5iv.labs.arxiv](https://ar5iv.labs.arxiv.org/html/2206.03619) |

I would not collapse those into one leaderboard.

## Where the descendants matter

### SBART: when smoothness is useful signal

Soft BART changes the function class, not just the optimizer. It replaces hard routing with smooth gating and has theory for adapting to unknown smoothness and sparsity. That makes it particularly relevant when your target is a smooth response surface rather than a collection of abrupt thresholds. [academic.oup](https://academic.oup.com/jrsssb/article/80/5/1087/7048381)

But “more suitable smoothness prior” is not the same as “always beats boosting.” The accelerated-soft-BART literature itself reports differences across noise regimes and between BART, XBART, and SBART-like constructions. The ranking depends on the data-generating structure. [arxiv](https://arxiv.org/html/2310.13975v1)

For your spline work, that is an interesting comparison: not Bayesian versus non-Bayesian, but which structural prior spends the available data most effectively.

### BCF: judge it on effects, not merely outcome RMSE

BCF’s central argument is that excellent outcome prediction does not guarantee an accurate treatment contrast. It parameterizes baseline prognosis and treatment response separately so they can be regularized differently, with propensity-score inclusion addressing regularization-induced confounding. [arxiv](https://arxiv.org/pdf/1706.09523v3.pdf)

A meaningful comparison is therefore BCF against a complete causal estimator using boosted nuisance models—not BCF against raw XGBoost trained to predict \(Y\).

For example, an experiment should separately evaluate:

- outcome prediction;
- average-effect error;
- conditional-effect error when ground truth is available;
- interval coverage;
- sensitivity to overlap and targeted selection.

That is a proposed evaluation design, not a claim that BCF necessarily wins each metric.

## The uncertainty advantage needs a qualification

BART gives posterior draws of the latent regression function, so quantities such as

\[
\Delta^{(s)}=f^{(s)}(x_a)-f^{(s)}(x_b)
\]

can be computed using paired evaluations from the same draw. The original formulation explicitly supports posterior inference for the function and its predictor effects. [asu.elsevierpure](https://asu.elsevierpure.com/en/publications/bart-bayesian-additive-regression-trees/)

That is useful for your purposes because you care about downstream functionals, not only pointwise predictions.

But posterior sampling is not a certificate of calibrated uncertainty. The XBART follow-up’s finding that a different initialization improves BART credible-interval coverage is a concrete reminder that computation and exploration affect inferential quality. [arxiv](https://arxiv.org/pdf/2002.03375.pdf)

Likewise, I would not automatically equate variation across accelerated stochastic forests with samples from the exact canonical BART posterior. The original XBART paper calls its method a heuristic; its predictive success and its inferential interpretation are separate questions. [proceedings.mlr](https://proceedings.mlr.press/v89/he19a/he19a.pdf)

## What I would conclude for ThermoMapper

I would not add BART on the premise that you are implementing a superior replacement for CatBoost or LightGBM. The sources here do not establish that broad claim, particularly for current, equally tuned implementations.

I would add it because it offers:

1. A genuinely competitive adaptive function prior, supported by both historical BART comparisons and newer XBART results. [arxiv](https://arxiv.org/pdf/2002.03375.pdf)
2. Posterior function draws that fit the purpose of a Bayesian regression engine. [asu.elsevierpure](https://asu.elsevierpure.com/en/publications/bart-bayesian-additive-regression-trees/)
3. Structured extensions such as BCF and SBART that change the scientific model, not just the training procedure. [arxiv](https://arxiv.org/pdf/1706.09523v3.pdf)

For benchmarking your implementation, I would report two tracks: equal wall-clock budget for prediction, and inferential quality for posterior functionals. Include tuning time in the first, and coverage plus effective sampling performance in the second. That would reveal whether your implementation is a useful predictor, a useful posterior engine, or both—without asking one RMSE number to answer all three questions.

