**Yes—the next phase should be an integrated Bayesian inference framework, not simply another set of regression implementations.** After KHAOS, GBMARS, and adaptive Bayesian multivariate regression splines, I would make the existing regression families share a clear contract for **posterior targets, conditional updates, predictive distributions, diagnostics, and evidence**.

The strongest starting point is already in ThermoMapper: its sampling target distinguishes the prior from temperature-dependent likelihood evaluations, including collapsed models. That is an important architectural foundation to preserve. [Sampling target](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/src/maths/samplers/targets/SamplingTarget.cs)

I’m treating “branch” as the next architectural phase rather than assuming you specifically mean tree-based regression. Also, I don’t have the titles of “all those papers” in this conversation, so I won’t pretend to have matched your newly acquired collection against the code.

## Where to concentrate

The existing documentation describes three substantial building blocks: model-neutral sampling, reusable conditional-update primitives, and shared probabilistic distributions. [Samplers](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/docs/Maths/Samplers.md) [Probability primitives](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/docs/Maths/ProbabilityAndInformation.md)

**My recommendation is to integrate around mathematical capabilities—not around algorithm names.** A regression method should expose what it can evaluate, integrate, sample conditionally, and predict; the inference machinery should consume those capabilities without knowing whether the model uses knots, hinges, trees, or kernels.

| Layer | Proposed responsibility | Avoid |
|---|---|---|
| Distributions | Probability laws, stable tails, support, parameterizations, sampling | Model-specific posterior-update logic |
| Samplers | Transitions, populations, adaptation, diagnostics, execution state | Regression-specific design matrices and basis semantics |
| Regression | Likelihoods, priors, structures, sufficient statistics, predictive laws | Reimplementing generic sampling machinery |
| Assessment | Predictive evaluation, calibration, model comparison | Treating one convergence number as a correctness certificate |

These are proposed boundaries, not a claim that every current implementation already follows them.

## Preserve the tempering distinction

This is the most important technical point in the current target contract: **tempering a likelihood and then integrating nuisance parameters is generally different from powering an already integrated likelihood**. ThermoMapper’s interface explicitly accommodates that distinction. [Target contract](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/src/maths/samplers/targets/SamplingTarget.cs)

For structure \(s\) and nuisance parameters \(\eta\), the collapsed likelihood path is

\[
\ell_b(s)
=
\log \int
p(\eta\mid s)\,
p(y\mid s,\eta)^b\,d\eta,
\]

not generally

\[
b\log \int
p(\eta\mid s)\,
p(y\mid s,\eta)\,d\eta.
\]

That distinction is already reflected in the source’s definition of `LogLikelihood(b)`. [Target contract](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/src/maths/samplers/targets/SamplingTarget.cs)

For your upcoming spline additions, I would make this a mandatory design review:

- **Explicit target:** Which coefficients, scales, latent variables, and structures belong to the sampled state?
- **Collapsed target:** Which variables have been integrated out, under exactly which priors?
- **Path:** Is the implementation evaluating the integrated tempered likelihood or a different, deliberately chosen bridge?
- **Reference:** Does the zero-temperature endpoint exist and have a known normalizer?
- **Constants:** Are likelihood normalizers retained when computing evidence?

**Do not force the explicit and collapsed representations into identical update code.** Instead, require them to agree on posterior quantities where both representations describe the same model.

## Make capabilities explicit

The current interface allows models to refuse unsupported temperatures or likelihood slopes. [Target contract](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/src/maths/samplers/targets/SamplingTarget.cs)

For the next phase, I would supplement that interface with explicit capability metadata. The purpose would be to reject unsupported inference plans *before* starting a long run.

| Proposed capability | Why expose it |
|---|---|
| Supports arbitrary likelihood temperature | Determines whether a proposed tempering or annealing path is usable |
| Supports the prior endpoint | Prevents assuming that a run can start at \(b=0\) |
| Supplies absolute normalization | Separates posterior sampling from evidence estimation |
| Supplies temperature derivatives | Enables compatible path-sampling calculations |
| Supplies conditional draws | Enables model-aware blocked inference |
| Supplies gradients | Enables optional gradient-based continuous updates |
| Supplies pointwise predictive quantities | Enables compatible predictive assessment |
| Supplies exact prior draws | Enables prior initialization and suitable reference constructions |

I would keep these as narrow interfaces or descriptors. A giant `IBayesianRegression` interface that requires every model to provide everything would make the framework less extensible.

There is also an important naming issue for GBMARS: **“generalized” does not automatically mean generalized Bayes in the loss-based sense.** If you later add loss-based updating, give its learning rate a separate name from the sampler’s inverse temperature.

## Build a regression composition layer

For the post-spline branch, my preferred organizing principle would be:

> **Structure + likelihood + coefficient prior + hyperprior + inference plan + predictive readout.**

But I would make composition **compatibility-checked**, not unrestricted. A likelihood, prior, and basis representation should only advertise combinations for which the relevant mathematics is implemented.

For example, the framework could distinguish:

- A model with an analytic coefficient marginal
- A model with exact Gaussian conditional draws
- A model requiring latent-variable augmentation
- A model requiring non-conjugate coefficient updates
- A model whose structure and coefficients must move jointly

The literature investigation should ask **which combinations preserve useful structure**, rather than asking how many named regression algorithms can be added.

A particularly useful deliverable from your KHAOS/GBMARS/adaptive-spline work would be a small implementation matrix:

| Component | Questions to settle |
|---|---|
| Basis structure | What changes dimension, and what remains fixed? |
| Coefficients | Conditional draw, analytic integration, or numerical update? |
| Scale parameters | Which parameterization and hyperprior? |
| Structural prior | On ordered objects, labelled objects, or equivalence classes? |
| Prediction | Conditional mean only, or full posterior predictive draws? |
| Tempering | Which representation defines the path? |
| Validation | What independent calculation can verify it? |

That matrix would tell you which abstractions are genuinely shared and which should remain model-specific.

## Samplers: improve composition first

The documented sampler layer already includes reversible-jump machinery, blocked updates, tempering, SMC-related infrastructure, and flat-histogram methods. [Sampler architecture](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/docs/Maths/Samplers.md)

I would therefore **not make “add another general sampler” the first milestone**. First make it straightforward to assemble, inspect, and validate a model-specific inference plan.

For a structurally adaptive regression model, a proposed sweep might comprise:

1. A structural transition
2. A coefficient update, unless coefficients are collapsed
3. Observation-scale or latent-weight updates
4. Shrinkage and smoothing updates
5. Optional representation-changing updates
6. Predictive and diagnostic readouts

Each block should declare the state it modifies, the target it preserves, its adaptation policy, and its checkpoint state.

Only then would I prioritize additional continuous kernels. My proposed order is:

- Exploit exact conditional updates where available
- Improve parameterization and blocking
- Add a gradient contract where a real model benefits
- Benchmark additional kernels against existing alternatives

That keeps “state of the art” tied to demonstrated performance on ThermoMapper’s models rather than the presence of fashionable sampler names.

## Distributions: target inference needs

The distribution layer already emphasizes named parameterizations and log-space tail calculations. [Distribution design](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/docs/Maths/ProbabilityAndInformation.md)

My priority would be **correctness and consistency at inference-critical boundaries**, followed by additions justified by concrete regression clients.

| Proposed workstream | What to investigate |
|---|---|
| Support semantics | Endpoint densities, invalid parameters, infinities, and degenerate cases |
| Tail arithmetic | Narrow intervals, extreme truncation, underflow, and complementary tails |
| Sampling consistency | Agreement between generated draws and the declared probability law |
| Constrained variables | Truncated draws and transformations, when required by a model |
| Shrinkage priors | Conditional representations needed by the chosen regression priors |
| Discrete laws | Structural counts, allocations, or augmented likelihoods |
| Multivariate operations | Shared Gaussian precision-factor operations without duplicating linear algebra |

There is one concrete documentation issue worth fixing before expanding this surface: the probability guide lists quantile APIs but later states that no distributions have quantiles. [Probability guide](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/docs/Maths/ProbabilityAndInformation.md)

I would reconcile that inventory against source and tests rather than using the guide alone to decide which primitives are missing.

## Literature: organize by decisions

For the next literature pass, I would use the following research agenda. This is a proposed reading structure—not a claim that I have verified the contents of your downloaded papers.

| Literature strand | Question to extract | Implementation consequence |
|---|---|---|
| Adaptive spline regression | What is the actual posterior over structure and coefficients? | State representation and structural prior |
| Bayesian backfitting | Which conditional decomposition makes fitting practical? | Block design and sufficient-statistic caches |
| Transdimensional inference | How are forward/reverse proposals and reference measures defined? | Structural transition correctness |
| Shrinkage and smoothing | Which priors control complexity without pathological geometry? | Hyperparameterization and conditional updates |
| Robust regression | Which latent representations are computationally useful? | Observation-level augmentation |
| Distributional regression | Which response parameters vary with predictors? | More than one predictor component |
| Model assessment | Which predictive estimand is being evaluated? | Pointwise likelihood and validation-unit contracts |
| Multimodal computation | What failure remains after ordinary local updates? | Evidence for tempering, SMC, or flat-histogram use |

Your current sampler documentation already provides entry points into several relevant strands, including Lang–Brezger P-splines, Neal’s slice sampling, nonreversible parallel tempering, and stepping-stone evidence estimation. [Literature connections](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/docs/Maths/Samplers.md)

For each new paper, I would extract **the target, the update, the assumptions, the computational shortcut, and the validation evidence**. Those five items are more actionable than a prose synopsis.

## Define the validation standard

The existing documentation distinguishes monitored functionals, convergence verdicts, and several evidence-error limitations. [Diagnostics and evidence](https://github.com/aghado01/ThermoMapper/blob/440d8a7da8bae5a31233eaa2acf193d36c0edb44/docs/Maths/Samplers.md)

The next step should be a common qualification protocol across regression families:

- **Target checks:** Compare small cases against enumeration, quadrature, or analytic calculations.
- **Representation checks:** Compare explicit and collapsed implementations on common posterior estimands.
- **Transition checks:** Validate structural probabilities and conditional laws independently of end-to-end predictive accuracy.
- **Predictive checks:** Distinguish latent-function uncertainty from future-observation uncertainty.
- **Evidence checks:** Require a valid reference, complete constants, and independent repeated runs.
- **Execution checks:** Test checkpoint continuation, parallel reproducibility, and separation of readout randomness.
- **Performance checks:** Compare at matched wall time and matched inferential accuracy, including adaptation cost.

I would explicitly separate four verdicts: **target correctness, computational adequacy, predictive adequacy, and evidence qualification**. Passing one should not silently certify the others.

## Recommended next branch

My proposed sequence after the spline implementations is:

1. **Write a model-capability inventory.** Use the completed spline families to discover the real shared seams.
2. **Specify inference-plan contracts.** Establish block invariance, tempering semantics, normalization requirements, and restart behavior.
3. **Standardize posterior outputs.** Make predictions, pointwise quantities, weighted draws, and uncertainty semantics comparable.
4. **Close distribution gaps.** Add only primitives demanded by qualified model implementations.
5. **Qualify two contrasting models.** Choose one heavily collapsed model and one explicit blocked model.
6. **Add the next model family.** Use it to test whether the integration genuinely generalizes.

**If your intended next branch is specifically Bayesian tree ensembles, put them at step six—not at step one.** They would then serve as a strong test of whether the regression engine is truly basis-agnostic, rather than becoming a second independent inference stack.

**Yes—Bayesian trees deserve their own research branch, and AI evaluation gives that branch a much sharper purpose than “add another nonlinear regressor.”** My recommendation is to aim at **hierarchical evaluation models that estimate where systems fail, how reliably judges measure those failures, and which evaluations are worth running next**.

I haven’t been able to inspect the contents of `opus-bayesian-evals.md` or refresh ThermoMapper’s GitHub source in this turn, so the architecture below is a proposal—not a code audit or a review of the attachment.

## Separate the research axes

The interesting BART literature is not one succession of increasingly superior algorithms. Different developments change **routing geometry, leaf expressiveness, likelihood support, or posterior computation**; these distinctions matter when deciding what should be reusable in ThermoMapper. [Oblique BART](https://arxiv.org/html/2411.08849v1) [ridgeBART](https://arxiv.org/html/2411.07984v4)

| Research direction | What changes | Proposed AI-evaluation application |
|---|---|---|
| Classical BART | A regularized additive ensemble of hard-routed, constant-leaf trees. [Oblique BART, background](https://arxiv.org/html/2411.08849v1) | A reference model for nonlinear failure probabilities and interactions |
| SoftBART | Smooth combinations of leaf outputs through soft routing. [ridgeBART, discussion](https://arxiv.org/html/2411.07984v4) | Smooth capability or calibration surfaces |
| Generalized BART | Nonconjugate likelihood support through reversible-jump updates. [Linero](https://arxiv.org/abs/2202.09924) | Outcomes beyond Gaussian scores: investigate ordinal judgments, counts, and latency models |
| ridgeBART | Ridge-function expansions inside leaves. [Yee, Ghosh, Deshpande](https://arxiv.org/html/2411.07984v4) | Locally smooth performance surfaces with regime changes |
| obliqueBART | Splits along sparse combinations of continuous features. [Nguyen, Yee, Deshpande](https://arxiv.org/html/2411.08849v1) | Failure boundaries depending on combinations of task properties |
| RoBART | A learned rotation shared by the internal nodes of each tree. [Heo, Jeong](https://arxiv.org/abs/2610.10214) | An experimental alternative to independently oblique splits |
| GP-based BART extensions | Functional, Gaussian-process leaf outputs. [ridgeBART, background](https://arxiv.org/html/2411.07984v4) | Structured smoothness over time, budgets, or controlled perturbations |
| Treed-GP models outside additive BART | Partition-based modeling with GP components appears in the broader Bayesian-tree literature. [Oblique BART, background](https://arxiv.org/html/2411.08849v1) | A contrasting architecture for local regimes and uncertainty |

**Those application mappings are hypotheses to test, not published demonstrations of superiority for AI evaluation.** In particular, “newer,” “more expressive,” and “better calibrated” should remain separate claims.

### Three distinctions worth preserving

**Soft routing is not the same as smooth leaves.** ridgeBART retains hard partitions while changing the functions represented inside leaves; SoftBART changes how leaf contributions are combined. A smooth leaf does not, by itself, impose continuity across a hard partition boundary. [ridgeBART](https://arxiv.org/html/2411.07984v4)

**Oblique splits are not the same as tree-specific rotations.** obliqueBART permits node-level hyperplane rules, whereas RoBART shares a rotation across a tree; RoBART was submitted on October 7, 2026, so it belongs on the experimental reading list rather than being treated as an established implementation default. [Oblique BART](https://arxiv.org/html/2411.08849v1) [RoBART](https://arxiv.org/abs/2610.10214)

**Generalized BART is not automatically loss-based generalized Bayes.** Linero’s paper expands the likelihoods available to BART through nonconjugate posterior computation; it is particularly relevant to your inference-engine ambitions. [Generalized BART](https://arxiv.org/abs/2202.09924)

## AI evaluation as measurement

The application I would prioritize is **Bayesian evaluation as a measurement system**, not merely prediction of benchmark scores.

I would distinguish three objects throughout the API:

- **System behavior:** What a model or agent actually does under a specified task and execution policy
- **Measurement:** What tests, human raters, or automated judges report about that behavior
- **Decision:** Whether the evidence supports deployment, escalation, further testing, or a claimed improvement

That separation suggests several distinct products.

| Application | Proposed estimand | Role for Bayesian trees | Essential qualification |
|---|---|---|---|
| Failure-surface estimation | Failure probability conditional on task and configuration | Learn nonlinear interactions | A conditional model is not automatically a deployment-risk estimate |
| Comparative evaluation | Performance difference under a specified task distribution | Model heterogeneous advantages | Preserve pairing and shared task effects |
| Judge calibration | Relationship between judge reports and an anchored criterion | Model nonlinear judge errors | Agreement alone does not establish validity |
| Evaluation allocation | Which next observation improves a decision | Supply uncertainty and predicted outcomes | Preserve sampling-design information |
| Reliability envelopes | Conditions under which a requirement is satisfied | Estimate a conditional response surface | Do not certify unsupported regions |
| Intervention analysis | Effect of changing a system component | Model heterogeneous responses | Identification comes from the design, not the tree |

The last distinction is important: **trees can help represent an effect without making that effect identifiable**. I would reserve causal language for evaluations with an explicit intervention and a defensible design.

### An illustrative model

For repeated pass/fail evaluations, one proposed starting point is

\[
Y_{mir}\mid p_{mi}\sim\operatorname{Bernoulli}(p_{mi}),
\]

\[
\Phi^{-1}(p_{mi})
=
\alpha_m-d_i+h_m(x_i)+u_{g(i)}.
\]

Here, \(m\) indexes systems, \(i\) items, \(r\) repeated runs, \(x_i\) observed task properties, and \(g(i)\) a task family. The intended roles are system level \(\alpha_m\), residual item difficulty \(d_i\), a tree-based system–task interaction \(h_m\), and a family effect \(u_{g(i)}\).

This is a **proposed model**, not a uniquely identified decomposition. I would impose centering constraints and hierarchical shrinkage, and explicitly test whether item effects, family effects, and flexible tree terms are distinguishable with the available data.

For automated judges, I would add a measurement layer rather than silently identifying the judge’s output with \(Y\). For preferences, I would use a separate pairwise or ordinal observation model—including ties where appropriate—rather than forcing preferences into a pass/fail interpretation.

## The strongest application opportunities

### Comparative failure maps

My first prototype would answer:

> **Where does system A improve over system B, where does it regress, and how much does that matter under the intended deployment distribution?**

Define the conditional contrast

\[
\Delta(x)=p_A(x)-p_B(x),
\]

and a deployment-weighted contrast

\[
\Delta_{\mathrm{deploy}}
=
\mathbb{E}_{X\sim P_{\mathrm{deploy}}}[\Delta(X)].
\]

The important product would be the **posterior over these contrasts**, with the task distribution stated explicitly—not simply two independently fitted performance surfaces.

I would begin with interpretable task properties and controlled perturbations. Embeddings could be additional predictors, but I would not let an embedding-defined cluster substitute for an explanation of a failure mechanism.

### Judge auditing

A second prototype would investigate where automated judges disagree with a defined human or executable criterion.

The proposed model should distinguish:

- System failure
- Judge false positives and false negatives
- Position or presentation effects
- Task-family-specific judge behavior
- Disagreement in the criterion itself

I would require an anchored subset of observations before attempting to infer judge accuracy. Without that anchor, a flexible model may explain disagreement without establishing which judgment is correct.

### Reliability and cost surfaces

Your spline work offers a particularly interesting bridge: **smooth variation along controlled experimental axes, combined with tree-discovered interactions elsewhere**.

For example, I would investigate performance as a function of inference budget, context length, or perturbation strength, conditional on task family and system configuration. ridgeBART and GP-leaf extensions provide relevant modeling ideas, although their representational and computational trade-offs differ. [ridgeBART](https://arxiv.org/html/2411.07984v4)

I would not impose smoothness everywhere. A tool permission change or parser failure may call for a different representation from a gradually changing compute–performance relationship.

### Adaptive evaluation

This is where the regression engine, sampler, and evaluation application could become especially coherent:

> Fit a posterior → identify decision-relevant uncertainty → choose evaluations → update → reassess the decision.

But **active learning and valid population inference are different objectives**. Active statistical inference directly studies label allocation for inferential accuracy, rather than only predictive training. [arxiv](https://arxiv.org/abs/2403.03208)

There is also a valuable warning in the newer literature: unreliable uncertainty estimates can make active sampling worse than uniform sampling, motivating robust combinations of active and uniform allocation. [arxiv](https://arxiv.org/abs/2511.08991)

For ThermoMapper, I would therefore require:

- A declared target population or finite evaluation pool
- Logged selection probabilities where the estimator requires them
- Randomized exploration rather than exclusively deterministic acquisition
- Separation of discovery data from confirmation where appropriate
- Evaluation of the allocation policy, not just the fitted predictor

**The posterior should guide collection; it should not silently define the population being measured.**

## What belongs in ThermoMapper

I would organize the tree branch around **partition rules, leaf models, and inference strategies**, while keeping AI-evaluation semantics outside the generic mathematics layer.

| Proposed component | Responsibility | Boundary to enforce |
|---|---|---|
| Tree topology | Parent–child structure and structural mutations | No observation likelihood |
| Routing rule | Axis-aligned, categorical, oblique, or soft routing | No assumption that all rules have the same computational dependencies |
| Leaf model | Constant, linear, ridge-feature, or other outputs | Declare conditional and marginal calculations explicitly |
| Ensemble state | Contributions, caches, hyperparameters, retained draws | Do not equate tree labels with scientific meaning |
| Tree transition | Structural and parameter updates | Account for the actual forward and reverse proposals |
| Observation model | Response likelihood and optional augmentation | No hard-coded evaluation-domain semantics |
| Evaluation layer | Systems, tasks, judges, runs, selection design, estimands | Do not contaminate generic sampler interfaces |

These are conceptual responsibilities, not suggested names for existing files.

### Do not over-generalize routing

The distinction between hard and soft routing is architectural, not cosmetic. Soft routing combines leaf contributions, so a leafwise factorization available to a hard-routed model should not be assumed to survive unchanged. [ridgeBART](https://arxiv.org/html/2411.07984v4)

Similarly, the published obliqueBART construction calculates admissible projection ranges using linear programs. Replacing that construction with empirical cutpoint bounds would need a fresh examination of the prior and proposal—not merely an optimization pass. [Oblique BART](https://arxiv.org/html/2411.08849v1)

I would start with specialized, correct implementations and extract common interfaces only after two genuinely different variants work.

### Reuse the right algebra

ridgeBART is an especially useful bridge to your regression engine: conditional on the inner ridge features, its outer leaf weights have Gaussian updates, and those weights can be integrated out during part of the structural update. [ridgeBART](https://arxiv.org/html/2411.07984v4)

Its reported linear-in-data scaling assumes a fixed small leaf expansion; the fuller update cost also depends on the number of ridge features and leaves. That makes it a promising reuse opportunity for Gaussian precision algebra, not a reason to assume every richer leaf remains cheap. [ridgeBART](https://arxiv.org/html/2411.07984v4)

My proposed dependency direction is therefore:

**Shared numerical primitives → qualified leaf calculations → tree transitions → evaluation models.**

### Temper the declared model

For collapsed tree updates, retain the distinction discussed previously:

\[
\log\int p(\eta\mid T)\,p(y\mid T,\eta)^b\,d\eta
\]

is generally not

\[
b\log\int p(\eta\mid T)\,p(y\mid T,\eta)\,d\eta.
\]

That is a mathematical design requirement for the proposed implementation. I would also distinguish **integrating the leaves of one tree conditional on the others** from integrating the complete ensemble: the former does not automatically supply the latter.

## What “state of the art” requires

I would define the research target through a qualification matrix, not a claim that one descendant wins universally.

| Dimension | Proposed test |
|---|---|
| Predictive performance | Held-out log score and task-appropriate loss |
| Uncertainty | Calibration and interval behavior for the actual estimand |
| Structural exploration | Independent chains and function-level agreement |
| Generalization | Held-out task families, templates, and later system versions |
| Judge robustness | Anchored judge audits and sensitivity analyses |
| Evaluation efficiency | Decision accuracy at matched total evaluation cost |
| Reproducibility | Restart behavior, random-stream control, and traceable data selection |
| Decision validity | Explicit assumptions for any claimed guarantee |

Posterior contraction results in ridgeBART and RoBART are valuable theoretical evidence, but they are **not, by themselves, demonstrations of finite-sample credible-interval coverage for an AI-evaluation application**. [ridgeBART](https://arxiv.org/html/2411.07984v4) [RoBART](https://arxiv.org/abs/2610.10214)

Where a formal risk guarantee is required, I would investigate a separate calibration layer. Conformal risk control is relevant to expected monotone-loss control, but its assumptions and target must be checked; it is not an automatic guarantee for arbitrary adaptive evaluations or distribution shift. [arxiv](https://arxiv.org/abs/2208.02814)

## A concrete first programme

My recommended first application is **paired comparative evaluation with adaptive human auditing**. It exercises the important scientific and engineering questions without requiring a universal evaluation platform.

1. **Specify the estimand.** Choose two systems, a task population, a human-defined criterion, and a practically meaningful difference.
2. **Establish a reference implementation.** Fit hard-routed BART with a response model appropriate to the criterion.
3. **Add hierarchy and repetition.** Represent shared items, task families, system configurations, and repeated runs.
4. **Add a judge measurement model.** Use an anchored audit subset and retain uncertainty about judge errors.
5. **Compare one smooth alternative.** Choose SoftBART or ridgeBART because the experimental axes justify it—not because both exist.
6. **Compare one geometric alternative.** Test oblique rules where combinations of task properties plausibly define boundaries.
7. **Introduce adaptive collection.** Compare against uniform or stratified allocation using the same total evaluation budget.
8. **Confirm independently.** Reassess important discovered regressions on fresh tasks or a reserved confirmation set.

For the immediate reading order, I would put **generalized BART first for inference architecture, ridgeBART second for reuse of your regression algebra, and obliqueBART third for partition geometry**; SoftBART supplies the contrasting routing design, while RoBART remains an experimental comparator. [Generalized BART](https://arxiv.org/abs/2202.09924) [ridgeBART](https://arxiv.org/html/2411.07984v4) [Oblique BART](https://arxiv.org/html/2411.08849v1) [RoBART](https://arxiv.org/abs/2610.10214)

**The opportunity I would pursue is not “BART for benchmark scores.” It is a Bayesian evaluation engine that connects heterogeneous behavior, imperfect measurement, and cost-aware evidence collection.** Trees would be one powerful response-surface family inside that system—not the definition of 

**Yes—both belong on the shortlist, and my previous map underweighted them.** DPMBART expands the **error distribution**, while survival BART expands the **observation process to censored time-to-event outcomes**; NFT BART provides a particularly relevant connection between these directions. [pubmed.ncbi.nlm.nih](https://pubmed.ncbi.nlm.nih.gov/36932826/)

For your AI-evaluation interests, I would investigate **survival BART for budget-dependent agent reliability**, and **DPMBART as reusable distributional machinery**, rather than treating either as another tree-splitting variant.

## Where they fit

These methods address different limitations, so they should not be ranked as interchangeable successors.

| Method | What it changes | Important restriction |
|---|---|---|
| **DPMBART** | Replaces Gaussian errors with a Dirichlet-process mixture of Gaussian components.  [arxiv](https://arxiv.org/pdf/1807.00068.pdf) | The original construction uses a shared residual distribution—not an arbitrary covariate-dependent distribution.  [arxiv](https://arxiv.org/pdf/1807.00068.pdf) |
| **Discrete-time survival BART** | Models conditional event probabilities over time using probit BART. [Original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4899272/) | Expands each subject into time-indexed observations, increasing computational cost. [Original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4899272/) |
| **AFT BART with DPM errors** | Models log event time with a BART location function and flexible residual distribution. [NFT paper](https://academic.oup.com/biometrics/article/79/4/3023/7587539) | Covariates shift log time without changing its residual scale—the accelerated-failure-time restriction. [NFT paper](https://academic.oup.com/biometrics/article/79/4/3023/7587539) |
| **NFT BART** | Combines a BART location function, heteroskedastic BART scale function, and DPM errors.  [pubmed.ncbi.nlm.nih](https://pubmed.ncbi.nlm.nih.gov/36932826/) | More flexible than AFT, but still a structured location–scale model, not unrestricted conditional-density estimation. [NFT paper](https://academic.oup.com/biometrics/article/79/4/3023/7587539) |

## DPMBART: distributional flexibility

The original DPMBART construction can be written as

\[
Y_i=f(x_i)+\epsilon_i,\qquad
\epsilon_i\mid\mu_i,\sigma_i^2\sim N(\mu_i,\sigma_i^2),
\]

\[
(\mu_i,\sigma_i)\mid G\sim G,\qquad
G\sim\operatorname{DP}(\alpha,G_0).
\]

The shared random mixing distribution allows the residual law to depart from a single Gaussian, including asymmetry or multiple modes. [arxiv](https://arxiv.org/pdf/1807.00068.pdf)

**The subtlety is that observation-specific mixture allocations are not the same thing as covariate-driven heteroskedasticity.** Conditional on the shared mixture, the original model does not make its mixing weights or component parameters functions of \(x\). [arxiv](https://arxiv.org/pdf/1807.00068.pdf)

That matters for evaluation: a proposed application is modeling an irregular continuous outcome distribution, but learning that “hard tasks have greater variability” requires additional covariate-dependent structure. Also, I would **not interpret inferred mixture components as discovered failure mechanisms**; the NFT literature explicitly cautions that these components serve density estimation rather than necessarily representing meaningful populations. [NFT paper](https://academic.oup.com/biometrics/article/79/4/3023/7587539)

Computationally, DPMBART is attractive for your integration goal: its sampler alternates weighted Gaussian BART updates, residual-mixture updates, and a concentration-parameter update. That suggests a reusable residual-mixture subsystem rather than a monolithic `DPMBART` implementation. [arxiv](https://arxiv.org/pdf/1807.00068.pdf)

## Survival BART: evaluation trajectories

The original nonparametric survival BART models the probability of an event at a time point **conditional on no earlier event**, using time and covariates as predictors. Survival is recovered by multiplying the conditional non-event probabilities, and the formulation does not impose proportional hazards. [Original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4899272/)

My proposed AI-evaluation applications would include:

- **Time or tool calls to verified completion**
- **Budget consumed before an unsafe action**
- **Attempts required to find a counterexample**
- **Time to recovery after an injected fault**

This changes the question from “did it pass?” to **“what is the probability of reaching a specified outcome within a specified resource budget?”**

But the event definition is crucial. The competing-risks BART literature distinguishes an event that prevents another event from occurring from censoring that merely ends observation.  For an agent evaluation, I would therefore distinguish: [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC6954340/)

| Proposed observation | Proposed treatment |
|---|---|
| Run remains active when a preassigned observation budget expires | Candidate right-censoring |
| Agent irreversibly terminates unsuccessfully | Competing terminal outcome |
| Success occurs, ending the run | Success event |
| Evaluator stops because the trajectory looks hopeless | Potentially informative stopping requiring explicit treatment |

**Do not quietly code every unsuccessful run as censored.** That would change what the model is being asked to estimate.

## NFT BART connects them

NFT BART is especially relevant because it combines flexible log-time location, covariate-dependent scale, and a DPM residual law. Its paper also addresses centering and scaling constraints needed to separate these components. [NFT paper](https://academic.oup.com/biometrics/article/79/4/3023/7587539) [pubmed.ncbi.nlm.nih](https://pubmed.ncbi.nlm.nih.gov/36932826/)

For your application, the proposed distinction would be:

- The **location function** models typical completion or failure timing
- The **scale function** models how variable timing becomes across conditions
- The **mixture residual** captures departures from a simple standardized Gaussian shape

However, **NFT does not mean assumption-free**: after location and scale adjustment, it still shares an error distribution across covariates. Its published performance evidence supports investigating it—not declaring it universally superior. [NFT paper](https://academic.oup.com/biometrics/article/79/4/3023/7587539)

## Concrete ThermoMapper implications

There is a verified implementation issue here: your current `Truncated<T>` supports normalized densities, probabilities, and interval moments, but **does not implement sampling**. [Truncated.cs](https://github.com/aghado01/ThermoMapper/blob/HEAD/src/maths/distributions/Truncated.cs)

That becomes relevant immediately: original probit survival BART uses truncated-normal latent draws, while AFT/NFT approaches use truncated draws to augment censored event times. [Original survival paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4899272/)[NFT paper](https://academic.oup.com/biometrics/article/79/4/3023/7587539)

I would organize the work around four shared capabilities:

1. **Tail-stable truncated-normal sampling**, qualified separately from truncated-density evaluation.
2. **Weighted Gaussian tree updates**, shared with residual-mixture models.
3. **DPM inference state and transitions**, distinct from a fixed mixture distribution.
4. **Censoring-aware likelihoods and predictive readouts**, including survival curves and event-specific probabilities.

Your distribution documentation already lists log-survival and log-interval probabilities, so the numerical probability layer offers a useful starting point; that is not yet equivalent to a qualified survival implementation. [Probability primitives](https://github.com/aghado01/ThermoMapper/blob/HEAD/docs/Maths/ProbabilityAndInformation.md)
