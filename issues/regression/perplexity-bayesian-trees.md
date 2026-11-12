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

**The opportunity I would pursue is not “BART for benchmark scores.” It is a Bayesian evaluation engine that connects heterogeneous behavior, imperfect measurement, and cost-aware evidence collection.** Trees would be one powerful response-surface family inside that system—not the definition of the system itself.