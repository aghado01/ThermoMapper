# BASS port: calibration and SMC backends

A segment of the [full port brief](../bass-full-port-brief.md), which holds the scope and completion contract, the dependency boundaries, the stage table and the decision register. This segment covers stages 7 and 8: calibration against an explicit target with its reference sampler, then the fitting and cut-calibration SMC adapters, with the stage-1 decision D09. D05's target sequence, which the SMC adapters follow, is drafted in the [fitting segment](port-fitting-brief.md). Its checklist, where the reports go as work lands, is [port-calibration-checklist.md](port-calibration-checklist.md).

The decision drafts here are proposals for the user's design pass, not decisions. The register, with each decision's status, and the table of sources the drafts were drawn from are in the hub's [decisions section](../bass-full-port-brief.md#statistical-decisions-to-settle).

## SMC integration

Sequential Monte Carlo belongs beside existing general sampling algorithms. It is not a replacement for required BASS MCMC/tempering capability, and it does not define the statistical model.

The existing-engine baseline and the shared primitives used by its adapter must pass their qualification gates first. SMC does not repair an improper final target, a missing prior-density factor or an invalid conditional mutation. Reuse the ownership and cache-consistency checks for particle cloning and target-stage changes, with independent offspring RNG streams and workspaces.

| Consumer | General machinery in `samplers/smc` | Consumer-owned work |
|---|---|---|
| BASS fitting | Weighted population, resampling, stage advancement, mutation execution, ancestry and normalizer-ratio estimates | Initial distribution, target sequence, structural/hyperparameter moves, posterior readout. |
| BASS cut calibration | Reuse the population machinery across conditional distributions | Emulator-draw selection, conditional calibration targets, bridges, and outer averaging. |
| Future SPC population annealing | Reuse population mechanics | Potts energy, SW updates, physical temperature selection, state transfer, and thermal measurements. |

The generic module needs log-weight normalization; resampling methods; a population driver over supplied target distributions; optional adaptive intermediate stages; mutation execution with private workspaces; particle diagnostics; and checkpoint state. Its contracts must support density evaluation, initialization, advancing a particle, and safe state transfer/forking. Exact type names should be chosen after mapping these operations to existing contracts.

Keep weight ESS, ancestry/family diversity, and population diagnostics separate from independent-chain R-hat and autocorrelation ESS. Equal weights after resampling do not establish independent particles. A normalizer ratio becomes Bayesian model evidence only when the model and starting normalization justify that interpretation.

Step 3.2 replaced the whole-posterior `TemperedTarget` with declared paths, `TemperingPath.Powered` and `TemperingPath.Likelihood` on [`ISamplingTarget`](../../../src/maths/samplers/targets/SamplingTarget.cs), and each engine declares its own (Q5). A likelihood bridge, a powered collapsed posterior, and a marginalized powered joint model are different paths. BASS must declare its path (D05); replica swaps and all mutation blocks must preserve the corresponding stage target, and swaps evaluate the cross-target densities rather than assume a simple power form.

**Exact reference rung (tracked as OQ2 in the preparation brief).** The same machinery gives parallel tempering a β = 0 rung drawn exactly from a normalized reference, which turns each annealed restart into an independent draw, the setting of Syed et al.'s round-trip analysis. BARS has one on both paths: the knot prior on the likelihood path (with a proper variance prior for the Normal model), and k ~ Poisson(1) with uniform ordered locations on the powered path. Build it here, then rerun the preparation brief's OQ1 tempering comparison.

BASS's default improper variance prior cannot initialize a prior-to-posterior SMC run. Select a proper-prior mode or a specified normalized reference distribution with appropriate support. A proper final posterior alone does not make the prior a valid initial particle distribution or establish a marginal likelihood.

For cut calibration, the intended form to investigate is

\[
p_{\mathrm{cut}}(e,\theta\mid D_{\mathrm{sim}},y)
=p(e\mid D_{\mathrm{sim}})\,p(\theta\mid y,e),
\]

where `theta` includes calibration, noise, and discrepancy parameters as applicable. Establish a nested conditional reference first, with finite inner-sampling error assessed. Then use SMC to transport conditional populations while preserving the upstream emulator-draw distribution. Ordinary joint SMC reweights the emulator and therefore solves a different problem. Variable-dimensional emulator draws should be handled through their induced conditional densities, not by blindly interpolating incompatible term arrays. See [Mathews et al.](https://arxiv.org/html/2406.07555v2).

## Decision draft: D09

### D09 (draft): calibration target and covariance formulations

Sources: R's [calibration](D:/aghado01/packages/BASS/R/calib.R:42) and its documentation (a "modular Bayesian calibration model"), the inventory's calibration section, Mathews et al. (arXiv 2406.07555), and Francom et al.'s elastic calibration (arXiv 2305.08834), which pairs BMARS emulators with modular calibration.
- **The target is the cut posterior,** p_cut(e, θ | D_sim, y) = p(e | D_sim) · p(θ | y, e).
  - e is a retained emulator draw: structure, coefficients and component noise variances.
  - θ collects the calibration inputs, the measurement variances and the discrepancy.
  - The observations do not update the emulator. This matches R's design and the lineage's practice.
  - *Not adopted:* plug-in calibration on a single emulator point estimate, which understates uncertainty; and fully joint calibration, a different inferential object that the cut is meant to avoid under misspecification.
- **The model, in physical units (D07):** y = y_m + y_s ∘ K w_e(θ) + δ + ε.
- **The covariance formulation is a declared type, with validated dimensions.** R picks its main branch by whether ncol(`discrep.mat`) < length(y), while `type` picks the truncation covariance (defect fixture).
  - *Functional* (R's type 1): Σ = diag(truncation variance + σ²_{g(j)}) + (y_s∘K)·diag(σ²_{e,k})·(y_s∘K)ᵀ, with the discrepancy δ = D·v on a supplied basis and v drawn from its Gaussian conditional. R fixes v's prior precision at 1, so the port declares v's prior and its units. Solves and log-determinants go through the factor of I + UᵀD⁻¹U, the stage-2 diagonal-plus-low-rank path, without forming m × m inverses.
  - *Multivariate* (R's type 2): Σ = Σ_trunc + diag(σ_g)·R·diag(σ_g) + Σ_δ + (y_s∘K)·diag(σ²_{e,k})·(y_s∘K)ᵀ, with a supplied correlation R over variance groups, a fixed discrepancy covariance Σ_δ and a fixed mean shift. It uses the strict dense factor.
- **Variance groups.**
  - Each group declares its prior. The default is R's IG(ν_g/2, (ν_g/2)·ŝ_g²) from `sd.est` and `s2.df`.
  - R switches every group to a −log(1 + σ²) kernel when any `s2.df` is 0. The port makes that choice per group.
  - Proposals on log σ² carry their Jacobian.
- **Tempering** follows D05's likelihood path. R multiplies the variance proposals' Jacobians by the inverse temperature and leaves terms out of its swaps (defect fixture).
- **Parameter coordinates.** The bounded inputs are sampled in normalized [0, 1] coordinates, and results are reported in physical coordinates with the transform recorded. R returns the normalized values.
- **Categorical calibration parameters,** an extension beyond R's API, which has no categorical proposal. Francom et al. 2019 sample them jointly from their discrete full conditional: they enumerate every combination (162 in their study), evaluate the conditional posterior at each with the current emulator draw, and precompute each term's categorical part so that a combination costs one product with the continuous part. The port does the same while the combination count allows, and otherwise uses a Metropolis step on one categorical variable at a time.
- **Cache consistency.** Any change to δ, σ² or the emulator draw invalidates the cached likelihood before the next acceptance test. R's stale likelihood after a discrepancy draw is a defect fixture.
- **Reference sampler (stage 7): nested conditional.**
  - For s = 1, …, S, draw e_s from the retained emulator draws.
  - Run a chain on p(θ | y, e_s) to stationarity, with its own diagnostics, and keep its draws.
  - Measure the error from finite inner runs by varying their length.
  
  R refreshes e uniformly every iteration and takes a few θ steps, which does not equilibrate under each e, so it is not this sampler.
- **Cut SMC (stage 8)**, following Mathews et al.:
  - carry particles through μ_s = p(θ | y, e_s) over independent e_s, reweighting by the ratio of conditional densities, then resampling and mutating with the inner kernel;
  - order the e_s so that successive conditionals stay close (their ordering by a travelling-salesman path, here by distances between the draws' predictions);
  - bridge successive conditionals through their induced densities, p(θ | y, e_{s−1})^{1−t}·p(θ | y, e_s)^t, because emulator draws are discrete and their term arrays cannot be interpolated. Mathews et al.'s linear tempering interpolates a continuous cut parameter and does not apply directly.
- **Outputs:** θ in physical coordinates, the variances, the discrepancy and the emulator index for each draw, with burn-in, thinning and diagnostics applied explicitly.
- **Independent checks:**
  - the dense covariance against the low-rank path;
  - units: scaled and unscaled y give the same physical posterior;
  - recovery of known θ on synthetic simulators;
  - the nested reference's error against inner run length;
  - the R defect fixtures above.

## Room for extensions

Prep work in this segment, so that the [extensions](extensions-brief.md) plug in without reworking the port:
- **The measurement-error model is a declared component:** Gaussian with variance groups now. Robust measurement models then plug in through GBMARS's seam; the derivations review's §7, item 4, names the calibration measurement model as that seam's largest benefit.
- **Emulator draws behind an interface.** The cut posterior needs only draws from p(e | D_sim) and their predictions, so a non-BASS emulator (BGP, ABMRS) can implement it.
- **Cut SMC's ordering and bridging** are written against that interface, not against BASS's term arrays.
