
## The primitives

1. **Exact readouts are roots and integrals of piecewise polynomials.** BARS's peak and span readouts are one case of something more general. Any feature defined by a polynomial equation in x, f, f′ and f″ has an exact posterior per draw. So do features defined between two fits, since f − g on merged knots is still piecewise polynomial. That covers:
   - peaks, onsets and inflection points (f″ = 0);
   - time above a threshold and area under the curve;
   - the point where two curves cross;
   - break-even points;
   - any condition that clears denominators into a polynomial. For example, revenue p·D(p) peaks where elasticity equals −1, which is the polynomial condition D + pD′ = 0.
2. **Variance decomposition with uncertainty.** Sobol indices and active subspaces work on any tensor-product representation, including a surrogate fitted to a black box.
3. **A noise family.** GBMARS gives heavy tails, skew and quantiles.
4. **Calibration and comparison.** Conformal sets over retained draws, LOO, evidence and proper scores sit on top of every engine. That makes "does the conclusion survive a change of prior" a routine check.
5. **Curves as data.** Functional outputs, elastic alignment, whole-curve bands, and replicates as exchangeable units.
6. **Sequential design aimed at features.** You sample where the feature is uncertain, rather than where the function is. Bayesian optimization is one special case.

**Causal inference.**
  - Heterogeneous treatment effects as a surface over covariates get Sobol indices that say which covariates moderate the effect. Bayesian causal forests (BART) are the comparison point.
  - The point where one treatment overtakes another is a crossing readout.
  - Regression discontinuity and kink designs are knots at a known or unknown cutoff.
**Interpretable ML.** Any trained model is a deterministic simulator, so emulating it gives analytic global importance with interactions, plus active subspaces. Functional ANOVA on random forests (Hutter, Hoos and Leyton-Brown 2014) is the incumbent. The caveat: Sobol under a product input measure is an interventional question. On correlated tabular features it extrapolates off the data, which is the same problem interventional SHAP has.

# Cross-pollination worth exploiting
- Step or ramp? Fit the changepoint engine and BARS to the same data and compare them by LOO. That answers whether an apparent discontinuity is a jump or a steep smooth rise. The question comes up in regression discontinuity, regime shifts and dose thresholds, and few toolkits can pose it.
- One readout layer, several engines, deliberate disagreement    - .
    BASS’s degree-1 maxima sit at knot    - s.
    ABMRS gives smooth extre    - ma.
    BGP needs a grid.
- When peak posteriors disagree across engines, the disagreement shows which part of the answer comes from the p- rior.
- Three independent Sobol computations. KHAOS’s sums of squared coefficients, BASS’s hinge integrals and pick-freeze Monte Carlo compute the same quantity three ways. Agreement among them is a correctness test for-  each.
- Conformal as engine selection. All engines meet the same coverage guarantee, so interval width at fixed coverage is a fair efficiency comparison.
- The thermodynamic layer applies to everything. The tempering ladder, free-energy curve and real log canonical threshold came from SPC but apply to every engine. The threshold measures complexity in a way that compares across very different priors.

`typed_evals` does group-aware splits, refuses splits too small to fit, switches to Venn–Abers below 2,000 labels, reports held-out Brier and log loss, and fingerprints the model and rubric. What it calibrates is a one-dimensional monotone map from judge score to P(human pass), one metric at a time, used as point predictions. The gaps sit around that core.

## typed_evals: where it stops, and what fits

**1. The dataset-level number is biased and has no interval.** This is the most important gap, and the fix needs none of your engines.
- `pass_rate` counts items whose calibrated score clears the threshold. That is a classifier's positive rate, biased whenever the judge errs. Their own benchmark shows it: at the default cutoff, calibration flagged fewer answers.
- `mean_score` estimates the human pass rate only if the fitted curve holds on the new data. It carries neither the curve's uncertainty nor any check that it holds.
- The benchmark's intervals were computed outside the library, with the curve held fixed.

Prediction-powered inference fixes this (Angelopoulos, Bates, Fannjiang, Jordan and Zrnic, *Science* 2023). Take the mean of g(s) over all items, plus the mean of (y − g(s)) over a labeled sample independent of the fit. The interval is valid however good or bad the judge is, and it gets tighter as the judge gets better. PPI++ tunes the weight on the judge term so the result is never worse than using labels alone. It's a pure-Python change.

**2. Venn–Abers intervals are computed, then thrown away.** `predict` collapses the (p0, p1) pair to a point. The docstring is right that the pair isn't a confidence interval, but its width is the small-sample warning, and nothing downstream sees it.

**3. Calibration uses the score alone.** The same judge score means different things across task types, agents and response lengths. This is where your stack fits:
- *Smooth curves with a posterior.* BARS with your monotonization gives a smooth monotone calibration curve with a posterior. Isotonic gives a step function that is unstable at small n.
- *Binary labels ride your GBMARS seam.* Use Pólya–Gamma augmentation: with ω ~ PG(1, ψ) and κ = y − ½, the logit model is Gaussian in ψ, with working response κ/ω and precision ω. That's a working response plus weights, the same seam you're building for E1. You'd need a PG sampler in `Maths.Distributions`.
- *Covariates show where the judge is unreliable.* With covariates, the calibration map becomes a surface g(s, x). Its Sobol indices then tell you which covariates the judge's reliability depends on.
- *Pooling.* Partial pooling across metrics or domains lets small domains borrow strength.

**4. Thresholds are point choices.**
- The benchmark picks cutoffs by F1 on validation. The point where the calibrated curve crosses a target is a span readout, so it gets a posterior.
- The tool guard's 0.9 is a fixed number. Conformal risk control (Angelopoulos and coauthors) sets the threshold from labeled safe and unsafe calls. It guarantees, in finite samples, that the expected fraction of unsafe calls let through stays at or below α. A guard should be able to make that kind of statement.

**5. Drift goes undetected, as their docs admit.** A small fresh labeled sample each period gives a stream of residuals y − g(s). The exact changepoint posterior on that stream says when calibration broke, for example after a provider silently changed weights behind the same version string. Shifts in the raw score distribution are a cheaper, label-free first alarm.

**6. Nothing chooses which items to label.** The library only does a random hold-out. Active statistical inference (Zrnic and Candès, 2024) chooses which items get human labels to minimize the PPI variance. Your feature-targeted design idea applies too: spend labels where the threshold crossing is uncertain.

## Agent evaluation more broadly

The field looks ripe because most results are point estimates on small task sets, with the factors that matter left uncontrolled.

- **Factorial capability surfaces.** Cross model × scaffold × horizon × complexity, as in your benchmark proposal, and fit the surface with BASS using categorical inputs. The model × scaffold interaction index puts a number on the confound everyone acknowledges but nobody measures.
- **Time horizons.** METR's 50% time horizon comes from a logistic fit of success against log task length. A monotone BARS fit makes the 50% crossing a span readout without assuming the logistic shape. Comparing the changepoint engine with BARS then asks whether there's a cliff rather than a slope.
- **Repeated stochastic runs.** Repeated runs of one task are exchangeable units, so E16–E18 apply directly. BayesBag's paper is literally about reproducible model selection, which is the leaderboard problem: rankings that flip under resampling.
- **Trajectories as curves.** Treat cumulative progress against steps as a curve per run.
  - Elastic alignment separates an agent's pace from its path.
  - Changepoints locate where it stalls or loops.
  - Whole-curve conformal bands bound a new run's progress curve.
  - This is also the process-level scoring idea from your proposal, with uncertainty attached.
- **Auditing the judge.** At temperature zero a judge is a deterministic simulator. Perturb answer order, length, formatting and rubric wording in a factorial design, emulate the judge's score, and Sobol indices quantify position, length and format bias, including how those biases interact.
- **Cost against performance.** Fit success against tokens or dollars. The derivative readout is the marginal return on inference compute. The threshold readout is the budget needed to reach a target success rate.

## Where to start

For typed_evals itself, adding PPI to the pass rate is a small PR with immediate value, and an easy way into the project. Covariate-aware Bayesian calibration is where your stack offers something no one else has. It depends on the PG route and a Python binding, so it follows naturally once GBMARS's seam lands. For your own research line, the factorial capability surface with an oracle that needs no judge is still the stronger piece. It also has the judge audit as a natural companion.

Sources:
- [typed_evals on GitHub](https://github.com/TrustifAI/typed_evals)
- [typed-evals on PyPI](https://pypi.org/project/typed-evals/)

For an evals project, I'd build the statistics as a spine rather than add-ons. Evan Miller's *Adding Error Bars to Evals* (arXiv 2411.00640) covers the basics. Beyond it:

1. **Declare the estimand before the run.** Fix which population of tasks, which quantity (pass rate, mean score, a difference between systems, a curve feature), and at what level. Most evals never state the population, so their numbers have nothing to generalize to.
2. **Declare units and groups.** Say what is exchangeable (tasks, runs, source documents) and what clusters (same template, same source). Intervals, conformal sets and PPI all rest on this. typed_evals' `group_id` is the start of it.
3. **Pair every comparison.** The same tasks run under two systems give paired differences, which usually shrink the interval a lot. Most leaderboards compare unpaired means.
4. **Design the allocation.** With T tasks and R runs per task, the variance of a mean is σ²_task/T + σ²_run/(TR), so extra repeats only shrink the second term. Under a budget, the classic two-stage sampling optimum is R* = √((σ²_run/σ²_task)·(c_task/c_run)). Authored agent tasks are expensive and runs are comparatively cheap, which pushes toward repeats. Cheap tasks push toward more tasks and fewer repeats. A pilot gives the two variance components.
5. **Correct for multiplicity.** Many metrics times many systems means control the false discovery rate on the claims you report.
6. **Monitor.** Changepoint posteriors on scores or residuals over time catch silent changes in models or judges.

Take outcomes y for agent a, task t, run r, application k, rubric item j. Agents and tasks are crossed. Applications sit inside tasks but recur across them, rubric items nest inside application × task, and runs nest inside agent × task. Most of the statistics follows from writing that down. Here's how I'd think about it, assuming several agent systems are compared on the same tasks; collaborating agents get a note at the end.

**1. Aggregation changes the answer.** A task that needs every application to succeed is conjunctive, while averaging rubric items gives partial credit. The two can rank agents differently, so declare which is the estimand and report both. Failures within a run are positively correlated: one misread fact poisons everything downstream. So the product of per-application success rates misstates task success, and the gap between them measures the correlation.

**2. A model that keeps the structure.** Roughly:

logit P(y = 1) = θ_ak − d_tkj + u_atr

- **θ_ak is agent a's skill on application k.** Each agent gets a vector, which is its capability profile.
- **d_tkj is the difficulty of rubric item j,** nested within application and task.
- **u_atr is a run-level shared shock.** Its variance is the "went off the rails" correlation from point 1.

This is multidimensional item response theory with crossed random effects. Partial pooling also handles most of the multiplicity across agents × applications × items: shrinkage does what corrections would, as in Gelman, Hill and Yajima's argument for why hierarchical models rarely need them.

**3. Dependence between applications is directional.** In cross-application tasks, application B often needs a fact from application A. Two consequences:
- *Scoring rule.* Decide whether B is scored against ground truth or against what A actually returned. The conditional score isolates B's skill; the absolute score measures the task. Record both.
- *Root-cause attribution.* Model success along the dependency graph: P(B ok | A ok) and P(B ok | A failed), a continuation-ratio structure. Attribute each failure to the first failing node, so a downstream failure isn't counted as a separate error. That's the root-cause grouping you practiced, made into an estimator.

**4. Design the application combinations.** If tasks are composed from applications, you can ask whether a pair is harder together than its parts predict. That integration penalty, an interaction term for each pair, is estimable only if the pairs co-occur in a balanced way. A balanced incomplete block design does it. With 6 applications and 3 per task, a (6, 3, 2) design uses 10 task templates: each application appears in 5 and each pair in exactly 2. Tasks are the expensive unit, so the allocation formula from before applies at the task level. Cluster every interval at the task, never at the rubric item.

**5. If rubrics are judged by an LLM**, the judges' errors are correlated across applications, especially when one model judges all of them. For conjunctive success you need joint calibration, not per-metric calibration. So label whole tasks, all applications at once, rather than sampling rubric items independently. typed_evals' sparse per-metric labels can't support the conjunctive estimand.

**6. Redundant rubric items.** Items within one application's rubric often measure the same thing, so a long rubric outweighs a short one. Item discrimination from the model in point 2, or a check on the effective number of items, tells you whether weights need fixing.

**Where your engines come in.** The core here is a generalized linear mixed model, not a curve fitter. Your engines enter on the continuous design variables:
- Success against the number of applications in a task shows whether decay follows the independence prediction p^k or falls off a cliff.
- Horizon and complexity can be crossed with the application factors.
- BASS with categorical inputs gives a nonparametric check on the model's variance decomposition, through the share of agent × application interaction.
- Binary outcomes go through the Pólya–Gamma route.

**If "multi-agent" means agents collaborating,** say one sub-agent per application, the attribution question becomes which agent's failure. A factorial swap design answers it: replace each sub-agent with a reference agent in turn, and decompose the effect on task success across the swaps.

That's the right reading: an eval becomes a stochastic computer experiment once tasks are points in a space you design. The user simulators I covered last time are a separate track. Most of what you pasted holds up. One claim needs correcting, two parts of it conflict with each other, and there are a few leads it doesn't have.

## One claim to correct: success on every application at once

BASS-PCA's joint draws describe your uncertainty about the mean success rates (μ₁(x), …, μ_K(x)). The chance that one run succeeds on every application, P(all succeed | x), is something else. It depends on how the outcomes vary together *within a single run*, and a model of the means says nothing about that, however good its joint posterior is. There are two ways to get it:

- **Treat "everything succeeded" as its own output.** It's cheap, but you lose the per-application information.
- **A multivariate probit.** Each run gets latent scores z ~ N(μ(x), R) across applications, and Y_k = 1{z_k > 0}. On the seams:
  - F-X1's working response, with the latent z as its Values (Albert–Chib). This needs only unit weights, so the current core takes it. Logit through Pólya–Gamma needs the weighted Gram, which is GBMARS's extension.
  - F-X6 with a dense Q = R⁻¹ across outputs.
  - A sampler for R. Chib & Greenberg (1998): R has to be a correlation matrix to be identified.

  Today the design refuses Q ≠ I and nothing samples R, so this is an extension, not a configuration. A coregionalized GP has the same gap unless its noise is also correlated across outputs.

## Where the paste conflicts with itself

**The difference estimator against adaptive designs.** The difference estimator is unbiased only because the n runs are a random sample. Replication-or-exploration and contour designs put runs where the emulator points. Once they do, the average residual no longer estimates the population's. Two ways to have both:

- **An audit stream.** A fixed share of runs is drawn at random and used for the design-based estimates. The emulator fits on all runs. It's simple, but the audit runs do less to locate boundaries.
- **Known sampling probabilities.** Draw each configuration with a probability π_t that the emulator sets, kept above a floor, and weight its residual by 1/π_t. This is active statistical inference (Zrnic & Candès 2024). Every run then serves both purposes, but the variance grows as the floor drops. The variance-optimal π_t is proportional to the root-mean-square residual expected at x, which is the emulator's predictive sd (noise plus its own uncertainty). So a heteroscedastic emulator supplies the sampling plan.

That trade-off belongs in docs/Maths with both options kept.

**Tune λ instead of fixing it at 1.** With a poor emulator, Var(y − f̂) can exceed Var(y), and the difference estimator does worse than a plain mean. The regression estimator fixes this. Survey sampling calls it GREG; the ML version is PPI++ (Angelopoulos, Duchi & Zrnic 2023). It uses λ·f̂ with λ* = Cov(y, f̂)/((1 + n/N)·Var(f̂)). λ = 1 gives the difference estimator, λ = 0 the plain mean, and the tuned λ is asymptotically never worse than either. It's still a few lines.

**Cross-fitting is unbiased only approximately.** Suppose the n runs are drawn from a finite pool of N configurations. Then each fold's residuals estimate the mean residual over units the fit didn't use, so cross-fitting carries a bias of order n/N. It becomes exact when:
- run configurations are independent draws from the generator, and
- the first term is the emulator's integral under the generator's measure.

That integral is what BassMoments computes under a product measure (the Monte Carlo tier covers other measures). So the integral and the correction become one estimator rather than two pieces.

## One estimator for both kinds of cost

Running an agent is expensive. Scoring a run is a second, separate cost: a judge is cheap, a human is expensive. The two stack:

μ̂ = ∫ f̂ dμ + (1/n) Σ_runs (judge − f̂) + (1/m) Σ_{human-scored ⊂ runs} (human − judge)

- Each correction term averages over a random subsample of the level above it.
- The telescoping sum is unbiased for the mean of the human scores.
- This is multilevel/multifidelity Monte Carlo (Giles 2008; Peherstorfer, Willcox & Gunzburger 2016, and their 2018 SIAM Review survey). It sets N, n and m from the cost of each level and the correlation between neighbouring levels.
- With a tuned λ at each level, this becomes the approximate-control-variate framework (Gorodetsky et al. 2020).
- A fixed judge needs no cross-fitting. A judge calibrated on the human labels does.

Because it's general UQ machinery, it belongs beside the engine, taking predictions from any emulator (BASS, BGP, KHAOS, an IRT fit). It isn't an evals feature.

## Smaller points

- **The judge's error belongs to each run, not each configuration.** The judge and the human score the same transcript. So the discrepancy is δ(x, run, agent), not just δ(x). Judges favour some agents' styles, including their own model family's (Panickssery et al. 2024). That agent term is what biases comparisons between agents, and a δ(x) averaged over agents can hide it.
- **Six applications as fixed levels or as a sample.** As a Sobol share over these six applications, the number is well defined. As a variance component for applications in general, six levels give about five degrees of freedom, whatever the run count. Say which one you're reporting.
- **Reweighting needs deployment tasks placed in the knob space.** Real tasks don't come labelled with horizon or complexity, so estimating those is a measurement problem; METR uses human completion time. Then report the Kish effective sample size of the weighted correction term. That one number tells you whether the design covered the deployment mix.
- **Design-based guarantees cover only some readouts.** They cover linear quantities like means and reweighted means, and, through PPI, regression coefficients and quantiles. Sobol shares, boundaries and active subspaces depend on the emulator being right. They're checked by held-out prediction and posterior predictive checks, not by unbiasedness.

## Leads the paste doesn't have

- **Psychophysics, for capability boundaries.** A 50%-success horizon is a psychometric threshold. Adaptive psychophysics has placed trials to estimate thresholds and slopes for forty years:
  - QUEST (Watson & Pelli 1983) and the Psi method (Kontsevich & Tyler 1999).
  - QUEST+ (Watson 2017) for several stimulus dimensions.
  - qCSF (Lesmes et al. 2010), which estimates a whole threshold curve along a second dimension. That's the direct analogue of a boundary over horizon × complexity.

  The psychometric function also has a lapse rate: success tops out below 1 even on easy trials. That's the "agent derails on easy tasks" mode, kept separate from the slope. Contour designs from computer experiments usually don't model lapses.
- **Differential item functioning (DIF).** A task feature that hurts one agent more than another is what psychometrics calls DIF, and it has its own detection literature. The BASS interaction index between agent identity and a task knob measures the size of that effect.
- **BASS's categorical factors share strength across applications in a third way.**
  - LVGP places each application at a point in a latent space.
  - Multidimensional IRT gives each one a latent vector.
  - BASS's CategoricalFactor groups applications into subsets within each basis function.

  How often two applications share a subset across the posterior, weighted by each term's variance, would give an application similarity matrix straight from BASS. That's an untested readout idea.

## What's built and what isn't

- **Built:** the integral under a product measure (BassMoments), categorical inputs, the working-response seam (the probit form fits it now; the logit form waits on weights), and cut posteriors (D09, stage 8).
- **Missing:**
  - the estimator layer (tuned λ, multilevel, known sampling probabilities);
  - the multivariate probit (dense Q plus a sampler for R);
  - a heteroscedastic noise model to set the sampling probabilities.

I'd still build the estimator layer first, but in the tuned, multilevel form rather than the fixed difference estimator. It's small, stays valid whatever the emulator does, and every later emulator plugs into it.


Yes, once tasks have features you control. When a task is a point in a design space (which applications, horizon, complexity, perturbations) rather than a one-off artifact, an eval becomes a stochastic computer experiment. Runs are expensive, you choose the design, and you never see the response surface whole. Several pieces from last time turn out to be emulation results under other names.

**The benchmark score is an integral of the emulator.** A benchmark's headline number is the mean of the success surface f(x) over the task distribution. With an emulator, computing it is Bayesian quadrature (O'Hagan 1991). It's also f_∅, the constant term of the Sobol decomposition, with the task distribution as the input measure. For BASS's hinge products it's exact per draw, which is the job of E7's moment service; for KHAOS it's the constant coefficient. Two consequences:
- *Reweighting.* Change the input measure and integrate again, and you get the score the agent would earn on a different task mix, such as a deployment mix. That only holds where the design covered the space. Outside it, the number comes from the prior.
- *Validity without trusting the emulator.* Suppose you generate N task configurations cheaply but run only a random n of them. Then (1/N)Σ f̂(x_t) + (1/n)Σ (y_t − f̂(x_t)) is unbiased for the mean over all N, whatever f̂ is. The condition is that f̂ is fitted on separate runs or by cross-fitting.
  - This is the difference estimator from model-assisted survey sampling (Särndal, Swensson and Wretman), and PPI has the same form.
  - The design supplies validity and the emulator supplies efficiency.
  - Agent evals suit it because generating a task is cheap and running one isn't.

**IRT is an emulator with a categorical input.**
- *Task features.* A random effect on task identity is a GP whose kernel on task labels is σ² when the labels match and zero otherwise. Give tasks features and the random effect becomes a function of them. Explanatory IRT (Fischer's LLTM; De Boeck and Wilson) does this linearly, and an emulator does it nonparametrically. The task variance left over after fitting the features measures how much your generator's knobs control difficulty. If it's large, something other than the knobs makes tasks hard.
- *Latent skill spaces.* Multidimensional IRT's logit, θ_a·α_t − d_t, is bilinear in latent agent and task vectors. That's the function class of a GP with a product of linear kernels on latent embeddings.
  - Latent-variable GPs for qualitative inputs (Zhang, Tao, Chen and Apley, *Technometrics* 2020) generalize this nonlinearly by placing each category at a point in a low-dimensional space.
  - Applied to applications, that map shows which ones call on similar skills.
  - Under a separable kernel over agent × application, the application correlation matrix is the correlation of last time's skill vector θ_ak across applications.
- tinyBenchmarks (Polo et al. 2024) already uses IRT to predict full-benchmark scores from a small subset, which is emulation by another name.

**Variance components are Sobol indices.** In a balanced crossed design, the random-effects variance components are the Sobol partial variances under the design's uniform measure, up to finite-sample corrections. So the mixed model's agent × application share and BASS's interaction index estimate the same quantity, which is why the BASS fit works as a check. The catch is that Sobol indices need a declared input measure, and a block design isn't a product measure over application subsets. In the (6, 3, 2) design only 10 of the 20 possible triples appear, so pairwise interactions are estimable but three-way interactions are aliased with each template's own effect.

**Tasks against repeats is the stochastic-emulator design problem.**
- *Replication versus exploration.* Binois, Huang, Gramacy and Ludkovski (*Replication or exploration?*, *Technometrics* 2019) treat it as a sequential design choice with noise that varies across the space. Last time's R* is the special case with constant noise and no smoothness across tasks.
- *Capability boundaries.* If the target is where success crosses 50% over horizon and complexity, that's a contour of the surface. Contour-estimation designs (Ranjan, Bingham and Michailidis 2008) put new runs near it. That's your feature-targeted design, with an existing literature behind it.

**Judges against humans is a multi-fidelity problem.**
- *The model.* Kennedy and O'Hagan's autoregressive model (2000), y_human(x) = ρ·y_judge(x) + δ(x), treats the judge as the cheap code and humans as the expensive one. The discrepancy δ can vary across tasks.
- *Relation to PPI.* PPI handles the same pair for a single estimand without a model. The multi-fidelity model is what tells you where in the task space the judge goes wrong.
- *Fitting jointly.* If you fit judge calibration and agent ability together, agent outcomes feed back into the judge's calibration, and a misspecified ability model then distorts it. D09's cut posterior blocks that feedback, and E14's semi-modular inference tunes how much of it to allow.
- *Benchmark against deployment.* The same structure applies, with δ(x) as the validity gap. It's estimable only where you observe both.

**Where it doesn't help.**
- *No task features.* With hand-authored tasks that have no features, the emulator reduces to the random-effects model and there's no smoothness to exploit.
- *Agents.* The same goes for agents, unless you describe them by features (model size, scaffold settings) or by latent embeddings.
- *Conjunctive success.* It needs a joint emulator across applications, because marginal fits can't give the probability that everything succeeds. BASS-PCA's joint draws or a coregionalized GP can.

Of these, I'd build the difference estimator first. It's a few lines of code, it stays valid whatever the emulator does, and it makes any emulator you already have reduce the number of runs you need.


R is the correlation matrix of the latent scores across applications within one run, the same R as in the seams brief's Q = R⁻¹. Everything would be in C#. Calling it a "sampler for R" was a bad choice in a repo where R also means the R oracles.

In plain terms, a multivariate probit gives each run a hidden score per application. The run succeeds on application k when its score is positive. The matrix R says how those hidden scores move together within a run. That's what's needed for "succeeds on all of them" or "fails on two at once". Building it in C# needs three pieces:

- **Draw the latent scores.** Draw them from a truncated multivariate normal, conditional on the observed successes and failures. They go in as F-X1's working response, Values = z.
- **Fit the mean.** Fit μ(x) with the existing collapsed core, using F-X6 with a dense Q = R⁻¹. Today the design accepts only Q = I.
- **Draw R itself.** This is the hard part, because R has to be a correlation matrix to be identified. Chib & Greenberg did it with Metropolis–Hastings steps. Parameter-expanded data augmentation is easier (Liu & Wu 1999; Lawrence et al. 2008): draw an unconstrained covariance with a conjugate inverse-Wishart step, then rescale it to a correlation matrix. That's one conjugate draw and a rescaling.

You only need it if you want readouts about how outcomes combine within a run, with the per-application models kept. If the only joint quantity you care about is "everything succeeded", model that indicator directly as one more binary output. That needs nothing new.

Sorry, my notation was ambiguous. R there means the matrix, not the R language. It's the correlation matrix of the latent scores across applications within one run, the same R as in the seams brief's Q = R⁻¹. Everything would be in C#. Calling it a "sampler for R" was a bad choice in a repo where R also means the R oracles.

In plain terms, a multivariate probit gives each run a hidden score per application. The run succeeds on application k when its score is positive. The matrix R says how those hidden scores move together within a run. That's what's needed for "succeeds on all of them" or "fails on two at once". Building it in C# needs three pieces:

- **Draw the latent scores.** Draw them from a truncated multivariate normal, conditional on the observed successes and failures. They go in as F-X1's working response, Values = z.
- **Fit the mean.** Fit μ(x) with the existing collapsed core, using F-X6 with a dense Q = R⁻¹. Today the design accepts only Q = I.
- **Draw R itself.** This is the hard part, because R has to be a correlation matrix to be identified. Chib & Greenberg did it with Metropolis–Hastings steps. Parameter-expanded data augmentation is easier (Liu & Wu 1999; Lawrence et al. 2008): draw an unconstrained covariance with a conjugate inverse-Wishart step, then rescale it to a correlation matrix. That's one conjugate draw and a rescaling.

You only need it if you want readouts about how outcomes combine within a run, with the per-application models kept. If the only joint quantity you care about is "everything succeeded", model that indicator directly as one more binary output. That needs nothing new.

Right. I drifted into designing the eval application, the same drift I made with SPC. I've added that to the memory note. Taken as engine questions, the eval ideas mostly land on seams the Francom lineage already put in place. They also point to a few pieces that are actually new.

| Idea from the evals thread | Engine axis | Seam | State |
|---|---|---|---|
| Probit / ordinal / Pólya–Gamma (logit) outcomes | Observation model | F-X1 working response | Probit fits now (unit weights); Pólya–Gamma and GH need the weighted Gram |
| Heavy tails shared across a run | Observation model | GBMARS latent scale, drawn per declared group instead of per observation | New indexing on a planned piece |
| Correlated outputs, factor-structured R | Output axis | F-X6, Q = (ΛΛᵀ + D)⁻¹ via Woodbury | New inner product + a covariance sampler |
| Random effects on many-level categoricals | Basis + prior | Prior-precision builder in the collapsed core | New, but KHAOS forces it anyway |
| Learned positions for categorical levels (LVGP / IRT) | Basis | F-X3 per-part contract | Lead |
| Multi-fidelity, y_hi = ρ·f_lo + δ | Calibration | F-X1 offset + scalar ρ + stage-8 cut | Almost just configuration |
| Guarded estimators (GREG / multifidelity Monte Carlo / ACV) | Functionals | BassMoments integral + predictions + cross-fit hooks | A new layer beside the engine |

What this says about the engine:

**1. F-X1 is a data-augmentation seam, not just a noise seam.** Several observation models make the likelihood Gaussian once you condition on an extra latent variable:
- GBMARS's GH mixture
- Albert–Chib probit truncation
- Pólya–Gamma weights
- ordinal cutpoints
- per-group scales

All of them hand the core (values, weights). The weighted Gram is what unlocks most of them, and it's part of GBMARS. That makes it the most useful single piece this thread surfaced. GBMARS's core then splits naturally: first the weighted Gram, an engine change to centring and the marginal; then the GIG/GH sampler, which is model-specific. This adds detail to the KHAOS → E8 → GBMARS order; it doesn't reorder it.

**2. KHAOS and random effects need the same seam.** KHAOS's modified g-prior already needs the collapsed core to accept a prior precision other than (1/g)·ΨᵀΨ, namely Σ⁻¹ = G∘ΨᵀΨ. If that becomes a declared prior-precision builder, ridge blocks are a second kind of prior it can build. Random-effect intercepts on task identity or on any grouping variable are ridge blocks. So are the scale-sensitive priors from the earlier mixed-dictionary discussion. The KHAOS work then pays off for every grouped-data problem, not just polynomial chaos.

**3. F-X6's inner-product family gets a third member.** It currently has the identity and the tridiagonal Markov residual. A low-rank-plus-diagonal form covers multi-response simulators, multi-task emulation and correlated functional residuals. Its sampler is general: factor loadings with parameter expansion (Ghosh & Dunson 2009), or a correlation matrix via inverse-Wishart then rescale. By the samplers-apart-from-models rule, it belongs in Samplers.

**4. Multi-fidelity comes almost free.** The offset was built for calibration discrepancy, y − f(x, θ). The same offset carries Kennedy–O'Hagan's AR1 multi-fidelity model: fit δ against y_hi − ρ·f_lo, and draw ρ as a scalar. The stage-8 cut decides whether the high-fidelity data feed back into f_lo. I'm not aware of a multi-fidelity BASS in the Francom papers I've read, so this is a lead the lineage itself doesn't have.

**5. The estimator layer is UQ-native, not specific to evals.** Multifidelity Monte Carlo and ACV come from the Sandia/LANL computer-experiments community, the same one as duqling. What it needs from the engine:
- a per-draw integral under a declared measure (BassMoments);
- predictions at points;
- a way to fit on folds.

The closed-form functionals the lineage insists on make the cross-fit estimator exactly unbiased for the superpopulation target. That's a concrete reason the lineage's emphasis on product measures pays off beyond Sobol indices.

**The one real gap:** a categorical level placed at a learned continuous position, with a hinge on that position. This is LVGP's idea expressed as a BASS factor. It stays separable, so under a categorical measure the moments are a finite sum. But it isn't conjugate, and it's identified only up to rotation and translation. That gives it a statistical cost, not just a build cost. I'd record it as a lead rather than a seam.