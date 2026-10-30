# The coefficient prior's two axes, and the ridge (chip brief)

Date: 2026-10-29. Builds the ridge priors D03 decided as options and stage 4 deferred (its choice 1). The prior is factored into a covariance form and a scale treatment, as the user decided on 2026-10-29. The decision is D03's *(2026-10-29)* bullet in [port-fitting-brief.md](port-fitting-brief.md); read it first. This also builds the fitting checklist's F-X5, the coefficient priors as a declared set.
- *What it builds:*
  - `BassCoefficientPrior(Form, Scale)`, replacing `SampledG` and `FixedG`;
  - the ridge form, with its marginal, tempered forms, draws and evidence;
  - the scale family in the shared `Maths.Regression.GPrior` namespace: the unit-information Gamma, hyper-g, hyper-g/n, and a held value.
- *Not in scope:*
  - F-R1's pre-registered comparison, which follows this build;
  - the robust prior and empirical Bayes;
  - BARS taking the family;
  - any default change: the default stays the g-form with the unit-information Gamma at a = ½, and draws bit-identically;
  - `docs/Maths` and the checklists, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then, in the fitting brief:
- D03 whole (the g-form, τ's family, propriety, the evidence, the other coefficient priors and the decision bullet);
- D04 (the collapsed and conditional update order) and D05 (the tempered forms);
- D06 (admission: a ridge is never a numerical repair);
- the contract shapes (`BassCoefficientPrior`, "the ridge priors join when built");
- the deterministic-simulator regime note.

Then `issues/regression/argon-friedman-etc.md`'s F-R1 section, with the Opus chat's corrections: the spectral formula, and the column scale.

In the code:
- *The shared seam:* `src/maths/regression/gprior/CollapsedGPrior.cs` (`FromCentredLeastSquares`, `CentredMarginalAt`, `GPriorMarginal`, `GPriorConditional`) and `Sigma2Prior.cs`, whose record style the scale family follows.
- *BASS:*
  - `BassConfig.cs` (`BassCoefficientPrior`, `SampledG`, `FixedG`, validation);
  - `BassResolved.cs` (the resolution of τ's prior);
  - `BassTarget.cs` (`TauPrior`, `FixedTau`, `SamplesTau`, `LogTauDensity`, the marginal and the tempered forms);
  - `BassDesign.cs` (the bordered Cholesky of the centred Gram, D06's gates);
  - `BassSweep.cs` (`BassTauBlock`, the composition of a retained draw);
  - `BassPriorSampler.cs` (τ's draw);
  - `BassEvidence.cs` (`BassEvidencePriors`);
  - `BassPointwise.cs` (the collapsed content);
  - `BassSmc.cs` and `BassLadder.cs`, which read τ through the target.
- *The oracles:* `r/oracles/bass_enumeration_oracle.R` (its `posterior` mode's `coef` key), `r/oracles/bass_tau_prior_oracle.R`, and `r/scripts/bass_stage1_checks.R`.

## The design, decided

**The configuration.** `BassCoefficientPrior` becomes a record of two closed sets, each with a `Kind` discriminator in `BassConfig`'s JSON style.
- *Form:*
  - `G`: σ²G_c⁻¹/τ, today's;
  - `Ridge`: σ²I/τ on the centred columns scaled to unit root mean square.
- *Scale:*
  - `Sampled(prior)`, with the prior from the shared family below;
  - `Held(τ)`. Under the g-form, null means unit information, 1/N, today's `FixedG`. Under the ridge, the value is required (D03: the scale is part of the model).
- `SampledG` and `FixedG` are deleted outright, with no shims, and every caller and test moves to the new shape. Fingerprints move; say so in the report.

**The scale family** (shared, `Maths.Regression.GPrior`, beside `Sigma2Prior`). Each member supplies its log density on u = ln τ with the Jacobian, an exact draw, whether it is proper, and its parameters. Each is defined on the form's unit-information precision τ₀ (1/N for the g-form, 1 for the normalized ridge) where its source defines it so:
- *the unit-information Gamma:* τ ~ Gamma(a, rate a/τ₀), whose mean is τ₀. a = ½ is Zellner–Siow, today's default (D03);
- *hyper-g* (Liang, Paulo, Molina, Clyde and Berger 2008): on g = 1/τ, p(g) = ((a − 2)/2)(1 + g)^{−a/2}, proper for a > 2 only;
- *hyper-g/n:* the same law on τ₀/τ. On the normalized ridge, where τ₀ = 1, the two coincide. State each mapping in its docstring.

**The ridge's marginal.** In the raw centred coordinates, the normalized ridge is the generalized ridge τD, with D = diag(G_c)/N.
- The posterior precision at b is bG_c + τD.
- The collapsed marginal carries ½ln|τD| − ½ln|bG_c + τD| and the quadratic b²·ỹᵀZ_c(bG_c + τD)⁻¹Z_cᵀỹ.
- σ² is integrated with its declared prior, as the g-form's is.
- Write it so the marginal does not depend on the coordinates the prior is stated in, and check that.
- *The spectral form* (the Opus chat's correction of argon's): with D^{−½}G_cD^{−½} = VΛVᵀ and w = VᵀD^{−½}Z_cᵀỹ, the quadratic is Σ_j w_j²/(bλ_j + τ), entering the exponent multiplied by b², and the log-determinant ratio is Σ_j ln(τ/(bλ_j + τ)). Evaluating it is O(M) once the decomposition is cached.

**The fast path.**
- Births, deaths and changes need the marginal at the current τ. Either maintain a second bordered factor of G_c + τD beside the design's factor of G_c (O(M²) per proposal), or decompose the candidate's normalized Gram (O(M³)). Choose by measured cost, and state the choice and the timings.
- τ's slice step and the tempered rungs read the spectral form, decomposing once per structure.
- The g-form's path is untouched.

**What the ridge changes elsewhere.**
- *The coefficient draws:* β | σ², τ ~ N((G_c + τD)⁻¹Z_cᵀỹ, σ²(G_c + τD)⁻¹), mapped back to the uncentred basis as the g-form's are.
- *σ²'s conditional,* with the ridge's quadratic.
- *The pointwise collapsed content:* the ridge's Student-t predictive per structural draw.
- *The evidence:* `BassEvidencePriors` records the form, its normalization and the scale prior.
- *The prior sampler and the reference rung:* τ drawn from the scale prior.
- *Admission:* D06's gates apply as under the g-form, so F is unchanged.

**Propriety.** A b = 0 rung, SMC from the prior, and the evidence each need a proper scale prior: refuse hyper-g with a ≤ 2 there, as D05 refuses the improper options. A held τ is proper.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **The Friedman measurement chip may still be running.** It changes no production file, but it reads BASS's code. Rebuild before you read its numbers, and expect its harness to compile against your changes. If its build breaks because of a file you changed, fix your file; never edit its harness.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-ridge\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-ridge\` when you finish.
- **Your files:**
  - in `src/maths/regression/spline/bass/`: `BassConfig.cs`, `BassResolved.cs`, `BassTarget.cs`, `BassDesign.cs`, `BassSweep.cs`, `BassPriorSampler.cs`, `BassEvidence.cs`, `BassPointwise.cs`, `BassModel.cs`, and the files that construct `SampledG` or `FixedG`;
  - `src/maths/regression/gprior/`: a new file for the scale family, and additions to `CollapsedGPrior.cs` that leave BARS bit-identical;
  - the R oracles above and `bass_stage1_checks.R`, with every existing check still passing;
  - the tests, and new ones;
  - this brief's Reports section.
- **Commit your own files only, by explicit path,** one commit per step, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked notes under `issues/regression/` stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs, C# and R are LF.
- **Bit-identity:** the default configuration and every g-form configuration draw exactly as before (`BassBitIdentityTests`' eight hashes), and BARS is bit-identical.

## The steps

Each step lands with its checks and ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next step. Name the passes and their fixes in the report.

1. **The shape.**
   - *Build:* the two-axis record, the scale family and validation, with today's configurations as G with the unit-information Gamma, and G with a held τ.
   - *Check:* bit-identity; the JSON round trip; every refusal, including hyper-g's a ≤ 2 where a proper prior is needed and a ridge without a held value where one is required.
2. **The scale family.**
   - *Build:* each member's density on u, its draw and its propriety.
   - *Check:* each density's normalization by quadrature on u; draws against the law within 5 SE; τ's posterior on a fixed structure against quadrature under both forms (D03's quadrature check). The R oracle `bass_tau_prior_oracle.R` takes the new members where that helps.
3. **The ridge's marginal and its fast path.**
   - *Build:* the marginal, its tempered forms and the fast path.
   - *Check:*
     - the marginal against a dense computation to rounding, at b = 1 and at two tempered b;
     - its independence of the coordinates;
     - the spectral form against the dense one;
     - the bordered (or decomposed) proposals against fresh evaluations over randomized birth, death and change sequences.
   - The R enumeration oracle gains the ridge in its `posterior` mode: the exact posterior over F on fixture A, held and sampled τ (the sampled one by quadrature over u).
4. **The draws and the rest.**
   - *Build:* the coefficient and σ² draws, the pointwise content, the prior sampler and the reference rung, the evidence's record, and SMC and tempering under the ridge.
   - *Check:*
     - the chain's structure posterior against the oracle's exact one, within 5 SE;
     - the composed draws against a joint-Gibbs reference on a fixed structure, as `BassSweepTests` does for the g-form;
     - the evidence against the oracle's exact ln Z under the ridge on fixtures B and C, from a ladder with the reference rung, and from SMC;
     - a direct functional fit and a mixed one under the ridge, end to end, every view running.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - the g-form bit-identical;
  - the ridge's marginal independent of the coordinates;
  - b entering the quadratic as b²;
  - the scale prior's Jacobian on u;
  - propriety enforced where a proper prior is needed;
  - the coefficients back in the uncentred basis;
  - BARS untouched.
- *Coverage:*
  - both forms under every scale member;
  - held and sampled τ;
  - the tempered forms and SMC;
  - D03's improper σ² options under the ridge;
  - a structure at M = 1;
  - near-collinear columns (D06 refuses them under both forms).

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the API and the configuration;
- the fast path's choice and timings;
- each check with its tolerance;
- the fingerprints and tests that moved;
- the suites with their counts (Regression, Samplers, the BASS oracle classes, R's stage-1 checks, run from `r/`);
- the review passes and their fixes;
- anything the main session must know.

## Reports

Newest first.
