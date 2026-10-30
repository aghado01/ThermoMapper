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

### 2026-10-29: the two axes and the ridge, built (chip report)

**What landed** (on `main`, nothing pushed):
- `79e0edb`, step 1, the shape: `BassCoefficientPrior(Form, Scale)`, the shared τ family and validation; `SampledG` and `FixedG` deleted.
- `90ce809`, step 2: the family's checks, and hyper-g in the τ-prior oracle.
- `e225c22`, step 3: the ridge's marginal, its fast path, the design layer and the R oracle's `form=ridge`. The draws, the collapsed pointwise content and the conditional means route through the form's block here, so that no commit draws the ridge's coefficients from the g-form.
- `0f92af0`, step 4: the checks of the draws, the pointwise content, the prior sampler, the evidence, SMC and tempering under the ridge. Their code landed in step 3; step 4 adds no production change.

**The API and the configuration.**
- *The record.* `BassCoefficientPrior(CoefficientForm Form, TauPrior Scale)`, with `BassCoefficientPrior.G(scale = Zellner–Siow)` and `BassCoefficientPrior.Ridge(scale)`. The form is a closed set of records: `GForm` (Kind `"G"`) and `RidgeForm` (Kind `"Ridge"`). The default serializes as `{"Form":{"Kind":"G"},"Scale":{"Kind":"UnitInformationGamma","Shape":0.5,"Tau":null}}`.
- *A shape call to confirm.* The scale axis is the shared family itself: `TauPrior`, in `Maths.Regression.GPrior` beside `Sigma2Prior` and in its record style (`Kind`, `Shape`, `Tau`). Its kinds are `UnitInformationGamma`, `HyperG`, `HyperGOverN` and `Held`.
  - The decision lists "a held value" among the family's members. So the brief's `Sampled(prior)` is every kind but `Held`, and `Held(τ)` is the family's own held member.
  - A wrapper with a `Sampled` variant would have held a second `Held` inside the family.
- *Resolution.* `TauPrior.For(n₀)` resolves the prior on the form's unit information τ₀ = 1/n₀ (n₀ = N for the g-form, 1 for the ridge) to a `ResolvedTauPrior`. It gives `LogDensity(u)` with the Jacobian, `SampleLogTau(rng)`, `HeldTau` and `IsProper`. The Gamma's rate is computed as a·n₀, so the g-form's numbers are bit for bit today's.
- *Hyper-g's domain.* Hyper-g and hyper-g/n are constructible for a > 1 and proper for a > 2. Between 1 and 2 the density is stated up to its constant and has no draw.
  - Validation refuses an improper member under the flat intercept, whose intercept-only structure keeps no column inside the prior, so τ has no posterior there. It also refuses one with SMC and with a reference rung.
  - Such a run reports no evidence, and a chain's start takes ln τ₀ in place of a prior draw.
  - Below 1 even a structure with a column inside the prior has no posterior for τ, so the constructor refuses it.
- *The ridge's intercept.* Under level centring and the zero centre the intercept is a unit-RMS column inside the ridge: α | σ² ~ N(level, σ²/τ). So k counts it, as under the g-form, and both forms give it unit information at τ₀. Under Liang's flat intercept it stays flat.
- *What moved in BASS:*
  - `BassResolved.Scale` replaces `TauPrior` and `FixedTau`.
  - `BassTarget` gains `Scale`, `HeldTau` (was `FixedTau`), `CoefficientForm` and `CoefficientBlock(design)`. The last is the (τ, b) → conditional map that the composition, the collapsed pointwise content and the conditional means read.
  - `BassEvidencePriors(Coefficients, Scale, …)` replaces `(Coefficients, TauPrior, FixedTau, …)`. It carries the declared form and scale, and the resolved scale with its τ₀ (the ridge's normalization).
  - `BassDesign(…, ridge)`, `BassDesign.Spectrum`, `DesignProposal.Ridge` and `BassStructureFit.Ridge` are new.
- *A held ridge's hot end.* It auto-picks min(τ/N, 1/N): the b at which the held τ halves a unit-RMS column's shrinkage τ/(τ + bN). Under the g-form it is min(τ, 1/N), as before.

**The marginal** (shared, `CollapsedGPrior.cs`).
- `RidgeSpectrum.Of(G_c, Z_cᵀỹ, N)` decomposes the unit-diagonal scaling G̃/N of the normalized Gram by cyclic Jacobi (`DenseEigen`, absolute tolerance 10⁻¹⁵), so λ = N·μ. It keeps w = VᵀD^(−½)c and qⱼ = wⱼ²/λⱼ.
- `RidgeColumns.At(τ, b)` gives `RidgeTerms`: ½Σ ln(τ/(bλⱼ + τ)), the shrinkage Σ qⱼτ/(bλⱼ + τ), and their slopes in b.
- `CollapsedGPrior.CentredRidgeMarginalAt` adds the intercept, a direction with λ₀ = N and energy Nȳ² under the zero centre, and returns a `GPriorMarginal`.
- S = b(‖r‖² + shrinkage), a sum of non-negative terms. This is the corrected spectral form, the exponent's reduction b²·Σ wⱼ²/(bλⱼ + τ), rewritten without cancellation.
- `CollapsedRidgePrior.FromCentred(spectrum, ‖r‖², ȳ, centre).At(τ, b)` gives the `GPriorConditional`. Its precision bG_c + τD is an `ISpdFactor` read through the eigenvectors, so the draws, `LeaveOutLogDensity` and `ConditionalVariance` are the shared ones.
- BARS is untouched. `CollapsedGPrior`'s and `GPriorMarginal`'s members are unchanged, but for two doc lines on `GPriorMarginal`; the log helpers moved into an internal `LogFunctions`; BARS's hash holds.

**The fast path: bordered on the current spectrum** (chosen by measured cost).
- A design under the ridge keeps its current structure's spectrum, decomposed lazily once per structure and carried in its snapshots.
- A birth, a death or a change is evaluated through that spectrum's eigenvectors (`RidgeSpectrum.Birth`, `Death`, `Change`): the parent's precision is Schur-reduced by the removed column, then bordered by the new one, and S keeps its non-negative form.
  - This costs O(M²) to start, one rotation of the border, and O(M) per (τ, b). It needs no factor of its own and no refactorization when τ moves.
  - It is the brief's bordered option, without a second factor to hold at the current τ.
- Where the design's strict factorization decided admission, within a factor of two of s_min, the candidate is decomposed fresh, and an acceptance keeps that decomposition.
- Measured (scratch benchmark, N = 500 hinge-like columns, per proposal, warm):

  | M | fresh decomposition | bordered birth | death | change | the design's border products (N·M) |
  |---|---|---|---|---|---|
  | 5 | 8 μs | 1.0 μs | 0.7 μs | 1.4 μs | 6 μs |
  | 10 | 14 μs | 1.7 μs | 1.4 μs | 2.8 μs | 4 μs |
  | 20 | 168 μs | 3.9 μs | 1.9 μs | 5.6 μs | 7 μs |
  | 40 | 932 μs | 11 μs | 3.8 μs | 17 μs | 12 μs |
  | 80 | 5,600 μs | 3.6 μs | 1.6 μs | 6.1 μs | 19 μs |

  The bordered path stays at or below the design's own per-proposal cost. The fresh one passes it from M ≈ 10 and grows as M³. Under the bordered path an acceptance pays one decomposition; Friedman's acceptance rate is 1–5%.
- End to end on S10's Friedman (N = 500, mean M ≈ 10–12, one chain, 3,500 sweeps, warm): the g-form took 0.25 s and the ridge 0.31 s.

**Checks, with their tolerances.**
- *The shape* (`BassContractTests`):
  - the JSON round trip of ridge, hyper-g and held configs, and the exact JSON of the default and of a held ridge;
  - every refusal: the family's own parameters, also on read-back from JSON; a missing form or scale; a held ridge without τ; an improper hyper-g under the flat intercept, with SMC and with a reference rung, while the proper members pass there.
- *Bit-identity.* `BassBitIdentityTests`' eight hashes held after every step, and so did sixteen more g-form options hashed before step 1. Those cover:
  - sampled τ at shape 0.25, held τ₀ and held 0.05;
  - the flat and zero centres, and Jeffreys;
  - functional data;
  - tempering (reference, linear, powered, self-sizing) and SMC (integrated, and linear with the flat intercept);
  - over draws, accounts, ladders, the evidence and the collapsed pointwise values.
- *The family* (`TauPriorTests`):
  - each member's density normalizes on u to 1e-12, and matches the law on τ (or on g′) with its Jacobian to 1e-12;
  - the closed forms hold: E[τ] = τ₀ and E[ln τ] = ψ(a) − ln(a·n₀) to 1e-11, and hyper-g's Beta(1, a/2 − 1) moments to 1e-11;
  - draws meet the law within 5 SE (40,000 draws, three quantiles and a smooth statistic);
  - hyper-g's Bayes factor from the C# density meets Liang et al.'s eq. (17) to 1e-9 (seen 1e-14);
  - an improper hyper-g's tail falls at the improper rate and it refuses a draw;
  - hyper-g/n is hyper-g on the ridge.
  - R's `bass_tau_prior_oracle.R` takes hyper-g and hyper-g/n, and the stage-1 checks read hyper-g's quadrature Bayes factors against eq. (17), at 4e-13.
- *τ's posterior on a fixed structure* (`BassCoefficientPriorTests`): the slice step against quadrature under both forms. The members are the Gamma at a = ½ and 2, hyper-g, hyper-g/n (the g-form), and improper hyper-g at 1.5. Each checks E[u], P(u ≤ median) and E[τ/(τ + τ₀)] within 5 batch SE over 60,000 sweeps.
- *The marginal* (`BassRidgeTests`):
  - against a dense factorization of the explicit design, under every centre × {scaled, Jeffreys}, at M = 0, 1 and 3, τ ∈ {10⁻⁴, 0.7, 40} and b ∈ {1, 0.3, 0.02}: S and the determinant term to 1e-11, ν exact;
  - independent of the columns' scale (×10⁻³, 7, 10⁴) to 1e-11;
  - on orthogonal columns it is the g-form at τ/N (S, the determinant term, the log marginal and its slope) to 1e-12.
- *Bordered against fresh.* On 600-step randomized birth, death and change sequences (the grid fixture, and 160 close knots of nearly collinear hinges, at s_min 10⁻¹⁰ and 2·10⁻³), every admitted candidate's terms at five (τ, b), b = 0 included, meet a fresh decomposition within 1e-8.
  - Seen: worst 6e-13 on the grid and 4e-10 at close knots near 10⁻¹⁰.
  - 417–526 candidates went bordered and 1–40 fresh.
- *Admission is the g-form's.* The same proposals give bit-identical accounts under both forms: 157 and 191 shared refusals.
- *The R oracle* (`bass_enumeration_oracle.R`, `posterior`, new key `form=ridge`). Each state's design is scaled to unit RMS and its spectrum taken by `eigen()`, an independent route. The oracle's own σ² quadrature of the explicit Gaussian integral agrees to ≤ 2e-13. A held ridge needs `tau=`.
  - `BassTargetParityTests`, fixture A, 299 states: ‖r‖², Q, the prior and each state's τ-integrated marginal match under mean, flat, zero, Jeffreys, flat + Jeffreys, held (2.5) and shape 2, the marginal to ≤ 7e-14, and E[τ | y] to 1e-8.
  - Fixture B's rungs, sampled, held (0.4) and flat + Jeffreys, integrated and linear: log Z_b to 1e-8 and E[slope] to 1e-7.
- *The draws* (`BassRidgeDrawTests`): the composed τ, ln τ, σ², σ⁴, the coefficients in the uncentred basis, a₁² and a₁a₂ on a fixed structure meet a joint-Gibbs reference within 5 SE. The reference samples a | σ², τ, then σ², then τ under τD₁ (100,000 composed draws against 400,000 Gibbs sweeps). For example, E[τ] is 1.4999 against 1.4942 ± 0.0041, and E[σ²] 0.27313 against 0.27349 ± 0.00015.
- *The pointwise content* (`BassPointwiseTests`): both contents against dense references and explicit leave-unit-out refits, under the ridge's precision τ·diag(XᵀX)/N, to 1e-10 retained and 1e-8 collapsed. The cases are scalar, categorical and functional × mean, zero and flat × scaled and Jeffreys × held and sampled.
- *Prior draws* (`BassCoefficientPriorTests`):
  - the sampler's τ meets its law within 5 SE under both forms (Gamma at ½ and 4, hyper-g, hyper-g/n), uncorrelated with M;
  - a held τ is held, and an improper hyper-g starts at τ₀.
- *Against the oracle's exact posterior, under the ridge:*
  - `BassChainParityTests` on fixtures A (sampled and held) and C: the count law within 5 SE + 0.003, the inclusions within 5 SE + 0.005, and E[τ] within 5 SE (0.49929 against 0.49917 ± 0.0017);
  - `BassTemperingParityTests`: each rung of a reference ladder (count, inclusions and E[τ]), and the self-sized ladder's evidence within 4 SE + 0.01. That is −27.893 ± 0.078 against −27.950 on B, and −30.855 ± 0.048 against −30.899 on C, and the record carries the ridge and τ₀ = 1;
  - `BassSmcParityTests`: the normalizer within 4 SE + 0.01, per replicate and pooled: −27.942 ± 0.028 against −27.950 (B), −30.891 ± 0.020 against −30.899 (C), and −27.228 ± 0.029 against −27.219 (B, held).
- *End to end* (`BassRidgeDrawTests`): fits under the ridge each run every view. The fits are functional, mixed, held, hyper-g, improper hyper-g, tempered with the reference rung, SMC, and flat + Jeffreys, with the g-form under hyper-g/n beside them. The views are the summary, traces (refused for SMC as designed), metadata, both pointwise contents with PSIS-LOO and WAIC (retained only), prediction with the mean surface, the training fit, sensitivity, effects and the residual check. The evidence exists exactly where it should.
- *Resume:* bit-identical under the ridge (`BassRunResumeTests`, scalar and functional), and the self-sizing ladder restored across its resizes (`BassTemperingResumeTests`).

**Fingerprints and tests that moved.**
- Every config fingerprint moves with the coefficient prior's serialization. `BassMovePassTests.S10Friedman`, a data hash over `Bass.Fingerprint(new BassConfig(), data)`, went from `a2946a493a466d54` to `2f8091b7c4adc06a`, the data unchanged.
- Tests moved to the new shape:
  - every `SampledG`/`FixedG` construction;
  - `FixedTau` → `HeldTau`, and `Resolved.TauPrior`/`FixedTau` → `Resolved.Scale`;
  - the contract tests' JSON strings and refusals;
  - `BassTemperingSizingTests`' hot-end reasons, which now name the b at which a held τ halves the shrinkage;
  - the evidence-priors record in `BassTemperingTests`;
  - `BassTemperingResumeTests`' state comparison, which reads a fit's ridge terms by value.

**Suites** (Release, `artifacts/chip-ridge`): Regression 1,319 passed, 0 failed (7.9 min; 16 of them the uncommitted before/after hash cases, since removed, and the sibling session's then-untracked Friedman reference file compiled in); Samplers 206 of 206; the BASS oracle classes (`Maths.Oracle.Tests`, every `Bass*` class) 90 of 90 (9.1 min); R's stage-1 checks from `r/`, exit 0, 74 checks pass, none FAIL.

**Review passes and their fixes.**
- *Step 1:*
  - the hot-end reasons and a dozen docstrings still said "fixed g" or "the g mode", and were reworded to "held" and "sampled";
  - the evidence's gates (`ReportsEvidence`, SMC's) gained τ's propriety;
  - `BassEvidencePriors.Of` refuses an improper τ prior;
  - the pinned data hash was moved, stated as such.
- *Step 2:*
  - the quadrature checks of the indicator statistics jumped by one node's weight, a test-side error; the CDF is now integrated on grids that end at the quantile;
  - the C# check against Liang's closed form was added, so the density is tied to an independent formula, not to quadrature alone.
- *Step 3:*
  - `EigenFast`'s absolute stopping threshold is unreachable on a Gram whose diagonal is N, so the spectrum decomposes the unit-diagonal G̃/N and scales back;
  - `GPriorMarginal`'s docs now state the ridge's S and determinant term;
  - the admission-equality test was added (F unchanged), which nothing had covered;
  - the parity test's τ nodes now centre on each form's τ₀.
- *Step 4:*
  - the end-to-end test misused WAIC on the collapsed content and traces on SMC, both of which refuse by design; fixed in the test;
  - resume under the ridge was added, MCMC and the self-sizing ladder, and the ladder comparison now reads the ridge terms by value;
  - end-to-end fits for improper hyper-g and the g-form's hyper-g/n were added.
- *Final, against the brief's list:* the g-form is bit-identical; the marginal is coordinate-free (the scale test, the dense route, R's route); b enters the exponent's quadratic as b² (dense checks at b = 0.3 and 0.02, and the oracle's rungs); the Jacobian is on u; propriety is enforced at the b = 0 rung, SMC and the evidence; the coefficients are uncentred (joint Gibbs); BARS is untouched. Coverage: both forms under every member, held and sampled, the tempered forms and SMC, the improper σ² options, M = 1, and near-collinear columns.

**For the main session.**
- `docs/Maths` and the checklists (F-X5) are yours. The API above is what they describe.
- The scale axis's shape (the family itself, `Held` its held member) is the one call made beyond the brief's letter; see above.
- Hyper-g between 1 and 2 is improper and accepted only where the posterior stays proper. Say if you would rather refuse a ≤ 2 outright.
- The ridge's qⱼ = wⱼ²/λⱼ divides by the normalized eigenvalues. Near s_min its digits follow the Gram's conditioning, as the g-form's Q does through its factor: bordered against fresh, 4e-10 at close knots near 10⁻¹⁰.
- F-R1 can run from here, with `BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma())` and `Ridge(TauPrior.Held(τ))` beside the g-form. The g-form also takes hyper-g and hyper-g/n now.
- While I worked, a sibling session's new `tests/regression/BassFriedmanReferenceTests.cs` briefly failed to compile (`Math.Log1p`). It blocked the shared test build until that session fixed it. I did not touch it.
