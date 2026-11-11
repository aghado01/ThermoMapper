# Active subspaces on the BASS core: E8's build (chip brief)

Date: 2026-11-08. Status: slices 0–6 complete; slice 7's plan written and waiting for the user's call. Commissioned by the user on 2026-11-08 ("kick off e8"), after KHAOS's last slice; slice 5, the Grassmann discussion's items the user settled, commissioned 2026-11-08 ("ill take as recommended"). Every slice landed and was verified by the coordinating session on 2026-11-08: slice 0 (`15fb495`, `9e3028d`, `aa62879`), slice 1 (`95294f8`, `9dfcc6d`, `1385b86`), slice 2 (`8c7db77`, `3cdc2e8`), slice 3 (`84fe230`, `1049055`), slice 4 (`20aeb06`, `bc72014`) and slice 5 (`def59f3`, `9d2ff23`, `95d2900`, `a54a537`). The follow-ups, commissioned 2026-11-10: slice 6 (`276af70`, `e9d8827`, `30df59e`, `9a669b6`, report `3c26660`) verified by the coordinating session on 2026-11-10; slice 7's plan (below, "Slice 7's plan") waits for the user's decisions; E36's brief follows its verification. Written by the main session from [E8's entry](extensions-brief.md#active-subspace-analysis), [gap G3](lineage-gaps-register.md#g3-sobol-indices-and-active-subspaces-have-different-standing), [the lineage review's separability question](lineage-review.md#separability-is-what-the-closed-forms-rest-on), the [Astra notes](../discussion/astra-grassman-stiefel-pca-extensions.md), and the [Grassmann discussion](../discussion/grassmann-stiefel-subspaces-brief.md)'s report so far. Reports go under "Reports" below, newest first.

## Why

An active subspace is the span of the directions in input space along which a function changes most, averaged over an input measure ρ: the top eigenvectors of C_ρ = E_ρ[∇f ∇fᵀ]. Sobol indices rank inputs and their subsets; C finds linear combinations of inputs, which Sobol cannot. Computed per posterior draw, the subspace comes with its uncertainty.

It serves any smooth surrogate of a costly code or system: reducing a calibration's parameter space, choosing directions to refit along, comparing two models' responses (co-active subspaces), and reading which output patterns move with the inputs. The build order puts it after KHAOS by the user's call of 2026-11-07, and E34 follows it.

## Sources

- **Rumsey et al., *Discovering Active Subspaces for High-Dimensional Computer Models*, arXiv 2307.11241.** The method on BASS: C under tensor-product representations through one-dimensional integrals, BASS's posterior propagated, integrals cached across structure changes. **Held** (v1) in `D:\aghado01\codex-scientiae\supellex-private\gauntlet\bayesian\2307.11241v1`, added by the user on 2026-11-08; slice 0 read it.
- **Co-active subspaces, arXiv 2311.18146,** held in `D:\aghado01\codex-scientiae\supellex-private\gauntlet\bayesian`.
- **Not held:** 1801.07922 (Zahm et al., vector-valued active subspaces), 1907.11572 (Wycoff, Binois and Wild, C for GPs), 2210.09181 (BPPR).
- **The oracle: `coactivity`,** pinned at `knrumsey/coactivity@35e4573` in `r/renv.lock` (`0a35e5d`), with `lhs` 1.3.0 and `zipfR` 0.6-70, which its code calls and its DESCRIPTION does not declare. It holds `C_bass`, `Cfg_bass`, `C_bass_mc`, `E_bass` and the bassPCA and bassfunc converters, and succeeds `knrumsey/concordance`, which announces itself deprecated (CRAN's `concordance` is an unrelated package). Read-only clones of both sit under `D:\aghado01\packages\bayesian-related`. Slice 0's oracle mode, `r/oracles/active_subspace_oracle.R`, runs `C_bass` on given draws and confirms its defects (Reports). E8 does not depend on it for correctness: it consumes what is built, and its internal checks are strong (below).

## What is already built

- **The moment service** (`BassMoments`): term-pair moments under a `ProductMeasure`, cached, with each factor's value and derivative (`FactorDerivative`). Its consumers are Sobol (`BassSensitivity.Of`), effects, and the basis fit's sensitivity (`BassBasisSensitivity`).
- **The measures** (`InputMeasures.cs`): `UniformMeasure`, `MixtureMeasure` (normal and Student components), `EmpiricalMeasure`, `CategoricalMeasure`. Each is one input's marginal, so every measure built today is a product measure. A dependent measure enters as a joint sample of points; G4's joint mixture of products is E25's.
- **The gradient readout** (`BassFit.Gradient` → `BassGradient`, KHAOS slice 4, `8f57487`): ∂f/∂z per draw, case and continuous input, in the input's own units, for hinge and Legendre fits, read from each draw's sampled coefficients. Averaged over any measure's points, it gives a C for any fit.
- **Coefficient readouts** (`BassSensitivity.FromCoefficients`): exact Sobol indices and the Jacobi derivative measure for Legendre fits under the product uniform law (Lean SG-4(b), `faca090`).
- **Draw weights** (`BassFit.DrawWeights`) from SMC and the flat histogram.
- **Geometry and linear algebra:** `GrassmannManifold` (principal angles by thin SVD), `GeometricMedian`, `ThinSvd`, the symmetric eigensolver. **`StiefelManifold` is not to be used:** its distance and log map ignore the frame (the Grassmann discussion's report). E8 needs no frames.
- **Lean:** SG-2 (√V_u is 1-Lipschitz in L²(μ); C has no such bound), SG-5 (`top_eigenspace_ne_subspace_mean`: E[C]'s top eigenspace is not the draws' mean subspace, 54° apart in a three-draw example), SG-4(b) (the Legendre bound total(v) ≤ ½·E[x_v(1−x_v)(∂_v f)²]).

## Decisions

These are the main session's, taken in commissioning. Each goes on the page with its alternatives and trade-offs. The user may reverse any of them.

1. **Two routes, one readout, chosen per query.** This settles the lineage review's open question for E8.
   - *Exact:* when the measure is a `ProductMeasure` and every term of the fit is separable, C per draw comes from the moment service with derivative factors.
   - *Sampled:* otherwise, C per draw is the average of the gradient readout over the measure's points.
   - *The readout says which route ran,* and for the sampled one its Monte Carlo error. The user may force the sampled route on a product measure, for a check.
   - *Alternatives:* sampled always (simple, never exact); exact only, refusing the rest (excludes dependent inputs and non-separable kernels).
2. **The input metric is declared.** For x = Lz, C_z = LᵀC_xL, so the subspace moves under any non-orthogonal change of units.
   - *Default:* each continuous input scaled to its box, the declared domain or the training range, as the fit's transform already does.
   - *Options:* the inputs' own units; or whitening by the measure's covariance, the generalized problem C v = λ Σ_ρ⁻¹ v, which is affine-equivariant and the natural choice under a dependent measure.
   - *Under a product uniform law* the box scaling equals whitening up to a constant. Every readout states its metric.
3. **Categorical inputs take no part in C;** it has no derivative along them. The readout names them, and their effect stays Sobol's. A fit's categorical terms still enter through each continuous input's gradient within its level.
4. **Two posterior summaries, reported apart** (G3, SG-5).
   - *E[C] and its eigendecomposition:* the eigenvalues, the top-k vectors, each draw's eigengap.
   - *The subspaces' projection mean:* the top-k eigenvectors of P̄ = Σ w_s W_s W_sᵀ. It is the Bayes action under chordal loss, and its eigenvectors in order give the nested estimates for every k (the flag mean, cited from memory).
   - *P̄'s eigenvalues are reported as what they are,* mean squared cosines. A value of 0.5 arises both when every draw sits at 45° to a direction and when half the draws contain it exactly, so the readout also gives each draw's principal angles to the mean.
   - *Outlying draws are reported, not down-weighted.* A minority mode is real posterior mass, so a robust median is at most a labelled option.
5. **Numerics.**
   - *Each draw's subspace* comes from the thin SVD of its N × p (square-root-weighted) gradient matrix, not from forming C, which squares the gaps.
   - *On the sampled route, the draws share the measure's points,* so their spread is the posterior's and not Monte Carlo's.
   - *The arccos floor* (about 1.5·10⁻⁸ rad) is below any posterior spread that matters. A sine-based angle is needed only where angles are differenced. Slice 5 moves the angles to LinAlg's stable routine, by the user's settlement of 2026-11-08.
6. **Per-draw consistency checks, reported, not enforced.** For a uniform input on [a, b], total(v) ≤ ((b − a)²/π²)·C_vv (Poincaré). Under a Legendre fit there is also SG-4(b)'s weighted bound with ½. A violation beyond Monte Carlo error is a defect.

## The slices

One commit per concern; tests with the code they test; each slice ends with its gates and a fresh-eyes review before the next starts.

| Slice | Delivers | Its own gates |
| --- | --- | --- |
| 0. Sources and references | The paper read (once supplied); whether the authors' code exists; reference functions with known C: ridge functions g(wᵀx), sums of ridges, polynomials, Legendre expansions; C by quadrature for each | References against quadrature to 10⁻¹²; each function's C in closed form where it has one |
| 1. Derivative moments | ∫φ′ψ, ∫φψ and ∫φ′ψ′ per factor pair beside the value moments, hinge and Legendre, under every `ProductMeasure` kind the value moments take; cached as the pairs are; C per draw on the exact route | Against brute-force quadrature per factor pair and seeded sweeps, far ranges included; C on fixed draws against the gradient readout averaged over a fine grid; SG-4(b)'s identity exactly on Legendre draws |
| 2. The readout | `BassFit.ActiveSubspace` (or the module's name for it) with its spec: measure, metric, k, route; per-draw C and subspaces, E[C], the projection mean, P̄'s spectrum, the principal angles and gaps, the Poincaré check, all under draw weights | Recovery at growing N on functions from outside the model (G3), with E[C] and the mean subspace apart; SG-5's three-draw example reproduced; the exact and sampled routes agree within Monte Carlo error; LᵀCL under a declared change of metric; weighted draws against replicated ones |
| 3. Functional outputs and co-active pairs | C = E[Jᵀ W_y J] with a declared output metric W_y (F-X6's Q), for direct and basis fits; C_fg = E[∇f ∇gᵀ] for two fits on the same inputs, and for two outputs of one | The basis fit's C against the direct computation on its reconstructed outputs; C_fg's symmetry and its reduction to C at f = g; 2311.18146's examples where they can be run |
| 4. Dependent inputs | The sampled route on a joint sample of points, with the whitened metric; a measure that says it is a joint sample | Against a product measure the route already gives exactly; equivariance under affine input maps with whitening; a dependent example where the box metric and whitening differ, both reported |
| 5. The Grassmann discussion's settled items | Per-draw k by a declared rule, opt-in; the output side G = E[J W_x Jᵀ]; the matrix Poincaré check; basis adequacy beside a direct-mode fit; the angles through LinAlg's routine (below) | Without a rule, the fixed-k readout bit for bit as before, apart from the angles; sums of r ridges recover r; tr P̄ against the weighted mean of k_s; the rank summary unbeaten by random projections of every rank, on seeded sweeps; tr G = tr C under the same metrics at W_y = I; a basis fit's G against the direct computation on its reconstructed outputs; the matrix check with equality on cos(πx), failing a deliberately mis-scaled D; adequacy zero on outputs inside the span, the known share outside it; slice 2's angles within arccos's error of the new ones above 10⁻⁶ rad |

### Slice 5: the Grassmann discussion's settled items

The user settled the [discussion](../discussion/grassmann-stiefel-subspaces-brief.md)'s open items on 2026-11-08 ("ill take as recommended"; "ill go with per-draw k per your recommendation"). Slice 3 had closed without the output side, so they form a slice of their own. The discussion's report has the derivations.

- **Per-draw k, opt-in.** The spec takes a declared rule for each draw's k_s and has no default: a share of tr C_s, or the largest relative eigengap up to a declared k_max.
  - *Without a rule,* the readout is unchanged but for the angles' new routine.
  - *With one,* it reports each draw's k_s and their weighted distribution beside the fixed-k readout. P̄ = Σ w_s P_s keeps its formula over mixed ranks, and tr P̄ is the mean rank. A rank-deficient draw (decision 3) enters at its own rank.
  - *A draw's distance to the mean subspace across ranks* is the projection distance, d² = |k_s − k|/2 + Σ sin²θ_i over the min(k_s, k) angles: a rank term plus a containment defect (Ye and Lim 2016, cited from memory). Both terms are reported.
  - *A rank summary, derived in the discussion and checked here.* Over projections of every rank, Σ_s w_s ‖P − P_s‖²_F = tr P + E[k_s] − 2 tr(P P̄). It is least when P keeps P̄'s eigenvectors whose eigenvalue exceeds ½. For coordinate subspaces that is the median probability model (Barbieri and Berger 2004, cited from memory). It is reported beside the mean rank, not used as a default k.
- **The output side.** G = E[J W_x Jᵀ] for direct and basis fits, in the declared input metric W_x, with the output metric W_y entering by its root as slice 3 takes it. It comes from the thin SVD of the stacked, metric-weighted Jacobians (m × Np), without forming G when m is large. The readout gives G's spectrum and top eigenvectors under draw weights: the output patterns that move with the inputs, which the Why names and slices 0–4 did not deliver.
- **The matrix Poincaré check.** Under a product law, Efron–Stein and each input's one-variable Poincaré inequality give Cov_ρ(f) ≼ E[J D Jᵀ], with D the diagonal of the inputs' Poincaré constants: (b − a)²/π² for a uniform input, so Cov_ρ(f) ≼ G/π² under the box metric.
  - Per draw, the readout gives the largest generalized eigenvalue of Cov_ρ(f) against E[J D Jᵀ], at most 1 up to Monte Carlo error.
  - It is decision 6's check in matrix form: reported, not enforced, and for the marginals decision 6 covers. A joint sample has no such bound and is not checked.
- **Basis adequacy.** For a basis fit beside a direct-mode fit of the same cases, it reports the share of the direct fit's G outside the basis's span, tr((I − P_K) G)/tr G, and the containment defect of G's top-k eigenspace in span(K).
  - A basis fit's own G lies inside span(K) by construction, so it needs the direct fit.
  - It asks D07's question: did truncation drop structure that moves with the inputs? [E37](extensions-brief.md#predictability-of-the-truncation-residual-e37) asks it from the truncation residuals, without a direct fit.
- **The angles through LinAlg.** `SubspacePosterior.PrincipalAngles` delegates to LinAlg's `PrincipalAngles`, and the mixed-rank distances read it.
  - That routine is the stable one (atan2 of cosines and residual sines) that the IRLS scale session is landing on main, under the user's settlement of 2026-11-08.
  - If it has not landed when the slice reaches this step, the rest proceeds and this step waits. Regression keeps no copy of it.

### The follow-ups, commissioned (2026-11-10)

Commissioned by the user on 2026-11-10 ("lets send the E8 follow-ups to E8 slice"). Already done elsewhere: per-draw k and the angles (slice 5), `khaos_helpers.R`'s copies, and test-select's name reading. `DenseEigen`'s relative threshold stays the user's LinAlg decision.

**Slice 6, the readout follow-ups.** Build these:
1. *A basis fit paired with another fit* (slice 3's follow-up): one overload on each class, in the terms-and-H form.
2. *The co-active projection mean,* by slice 3's rule: top-q by |λ| per pair, P̄ over the pairs, and the concordance's sign beside it.
3. *Output locations beyond the training ones* for a direct fit: `OutputMetric` gains its locations.
4. *The matrix Poincaré check's cost.* Cache the term-pair Grams (Σ_B, K_D) per structure. It is a pure cache, sized by a constant and never by hit history. Measure the check against the readout. Whether it defaults on again comes back to the user with those numbers.
5. *Poincaré constants for log-concave laws,* so the matrix check covers truncated normals: D takes σ² per such input.
6. *Activity scores.* C's diagonal is reported always, as the paper's sum over every eigenpair. Constantine's truncation at k is reported beside it when a rank or rank rule is given.
7. *Legendre moments at high degree:* expand in the Legendre basis under a uniform law, where the integrals are exact rationals, so that high-degree KHAOS fits keep 10⁻¹⁴. The KHAOS chip is idle. Check `git diff` on the moment service first.

**Slice 7, rotate-and-refit: plan first.** It changes the model, so write the plan as a section here before building, with the decisions for the user. Then wait for the user's call. The plan covers:
- the frame: the ordered eigenvectors of P̄ or of E[C], with no Stiefel code unless a simplicity rotation is wanted;
- the rotated inputs' box and law, and how BASS's knot prior and the moment service read them;
- what the refit reports beside the unrotated fit (evidence, recovery);
- what a user declares.

**E36 (BPPR) and the directional hinge: a build brief for the user.** Draft `issues/regression/bass/bppr-brief.md` in the shape of this brief: decisions, slices and gates. Draw it from [E36's entry](extensions-brief.md#projection-pursuit-ridges-on-the-bass-core-bppr-e36), the Grassmann discussion's settlement of its direction prior, and the read-only `gcbeck/bppr` clone (GPL-3.0, ideas only; validate against the paper as if no repo existed). The directional hinge enters as a term type behind `IDesignKernel` (F-X3). Nothing is built until the user approves the brief.

**For all three:**
- *Order:* slice 6, then slice 7's plan, then E36's brief. Each is reported and verified before the next.
- *Other writers:* the GBASS chip (E21 slice 4, the noise and block files), the E34 chip (its slice 5, the PCA mode) and the conformal chip (`ConformalBayes.cs`). On shared files (`BassConfig.cs`, `BassModel.cs`, `BassResolved.cs`, `BassSweep.cs`, the moment service), message the other chip before editing and commit each such edit by itself.
- *Docs:* chips leave `docs/` alone (AGENTS.md, 2026-11-10). This replaces the `docs/Maths` gate below and the shared-docs line under Working rules: list the doc changes owed in each report.
- *Messages* are terse, and no session reads another's transcript.

### Slice 7's plan: rotate-and-refit (2026-11-10), for the user's call

Nothing is built until the user decides the items marked **Decision** below. The plan follows the brief's four headings, then proposes slices and gates.

**What the source does** (2307.11241, read again for this plan):
- *§4.1.3 projects and refits.* It fits Bayesian MARS on the data projected onto the first d active directions, d = 1 to 6, for g(x) = f(x·t₁, x·t₂), with f = x₁² + x₁x₂ + x₂³/9, t₁ = (12, 9, 6, 1, 1, 1)/3 and t₂ = (1, 1, 1, 15, 6, 6)/3. Each refit is reported as the percentage of g's sd it accounts for: 91.06% at d = 1, 98.37% at d = 2, 99.45% at d = 6.
- *§4.3 refits on the projected data.* On 240 nuclear-data inputs it whitens by the prior's covariance and refits on the projected data. That refit runs in 38 s against 89 s, uses M = 44 basis functions against 223, and has an order of magnitude smaller RMSE.
- *The conclusion names the obstacle and the next step.* Iterating on the projected data is complicated because "the measure corresponding to the projected variables is a high dimensional composition of the measure for x and dependence among the projected variables will be introduced". It suggests a projected Bayesian MARS, which is E36's ground.
- *It never states which summary of the posterior's C it takes the directions from.*

**The frame.** The frame is the ordered eigenvectors W of a summary in the analysis's own metric. The active variables are u = W₁ᵀL⁻¹(x_box − c), with L and c the analysis's map (for the box metric, u = W₁ᵀx_box; for whitening, the paper's Σ^(−1/2)(x − μ)).
- BASS ignores an input's sign and the inputs' order, so W₁ fully specifies the refit, and no Stiefel code is needed (the Grassmann discussion, settled).
- **Decision 1, the summary.**
  - *E[C]'s eigenvectors (recommended).* For any frame W, the posterior mean of each draw's inactive energy Σ_{j>d} w_jᵀC_sw_j is Σ_{j>d} w_jᵀE[C]w_j, and E[C]'s trailing eigenvectors minimize it. That energy is what bounds a ridge approximation's error by Poincaré's inequality (Constantine's bound, cited from memory).
  - *P̄'s eigenvectors,* the Bayes action under chordal loss, each draw counted once. The two agree when the draws agree; under SG-5's disagreement they part by up to 54°.
- **Decision 2, a rotation within the span for simplicity** (varimax's role). Not recommended. It needs Stiefel code, and BASS's axis alignment already rewards a frame whose directions are E[C]'s eigenvectors.

**The inputs the refit takes.** **Decision 3:**
- *(a) Projected (recommended), as the paper does:* the d active variables, plus the categorical predictors unchanged. Variation along W₂ goes into the noise.
- *(b) Rotated:* all p rotated variables, so nothing is dropped and only the axes move.
- *(c) Augmented:* the p native inputs plus the d active variables, so the sampler chooses between axis-aligned and oblique terms, and the native model is a special case.

**The rotated inputs' box and law, and how BASS reads them.**
- *The box.* The input box maps onto a zonotope in u, whose bounding box is exact: [Σᵢ min(0, aᵢⱼ), Σᵢ max(0, aᵢⱼ)] per variable for A = L⁻ᵀW₁, shifted by c's image.
  - **Decision 4:** that bounding box as each active variable's declared domain (`InputDomain`) (recommended), so every image of the input box lies inside the refit's box and E8 readouts of the refit by a joint sample are never refused. *Alternative:* BASS's default, the projected training range, which a pushforward point can leave.
- *The knot prior* reads knots at the training values of each input. The active variables' training values are the projected cases, almost all distinct, so every case is a knot candidate. Nothing in the prior changes.
- *The law.* ρ_u is the pushforward of the analysis's law, and in general not a product (the paper's obstacle). So the moment service's exact route does not hold for the refit under ρ_u.
  - The refit's readouts under ρ_u take E8's sampled route on a pushforward joint sample: N points from the original law, mapped. That route exists (slice 4).
  - Readouts in x follow by the chain rule: ∇ₓf = A∇ᵤf̃, so C_x = A C_u Aᵀ, of rank at most d.
  - A product law declared on u is allowed, but it is another law, and the readout says so.

**What the refit reports beside the native fit.**
- *Size and fit:* M (the basis functions in use) and the run's time for both fits. With a test set, held-out RMSE and the percentage of sd accounted for, the paper's measure.
- *Predictive criteria:* PSIS-LOO's elpd for both fits (`PredictiveChecks.PsisLoo`).
- *Evidence:* ln p(y) for both, when both runs are tempered to the exact reference rung (`BassFit.Evidence`). It is opt-in, as the ladder is.
- *Consistency:* the refit's own C mapped back to x against the frame (principal angles), and the frame's posterior spread (each draw's angles to the mean), which the refit does not carry.
- *The caveat, in every report:* the frame was estimated from the same data, so in-sample criteria favour the refit.
  - **Decision 5:** cross-fitting as an opt-in (recommended). Over K folds, the frame comes from a native fit on the other folds and the refit is scored on the fold held out, at K + 1 native fits and K refits. *Alternative:* in-sample criteria with the caveat alone.
- *Frame uncertainty is a plug-in,* as the PCA basis's is (settled). Propagating it means refitting per frame draw, or E36's joint model.

**What a user declares:** the analysis the frame comes from (decision 1's summary), d, the inputs' form (decision 3), the refit's configuration (the native fit's by default), and optionally the folds and the evidence ladder. A name for the types, after the paper's "projected Bayesian MARS": `BassProjection` for the frame and the map, and `BassProjectedFit` for the refit, which predicts, differentiates and reads E8 in the original coordinates.

**Decided by the user on 2026-11-10 ("Approve as recommended"), relayed by the main session:**
1. The frame is E[C]'s eigenvectors.
2. No simplicity rotation. The report gives the eigengaps within the active block, where near-equal eigenvalues make the eigenvectors arbitrary.
3. The refit takes the d active variables (a). Form (b) is (a) at d = p, so it comes with it. (c), the native inputs plus the active ones, stays a later option.
4. Each active variable's domain is the zonotope's bounding box. A gate is added: how the ridge forms and any knot-spacing rule read that scaling when, in high p, the box is far wider than the projected data's spread.
5. Cross-fitting is an opt-in. The caveat about reusing the data stays in every report, the evidence included.

**Proposed slices, after the user's call:**

| Slice | Delivers | Gates |
| --- | --- | --- |
| 7a. The frame and the projected data | `BassProjection` from an analysis: the summary, d, the map; projected `BassData` and cases with the bounding box as their domains; the pushforward joint sample | the bounding box against the input box's vertices, every projected vertex inside, one on each face; the pushforward's mean and covariance against AᵀE[x] and AᵀΣA; refusals (another fit's inputs, d outside 1 to p) |
| 7b. The projected fit | `BassProjectedFit`: fit, predict and gradient in x by the chain rule; E8's readout in u on the pushforward sample, and in x by A | the gradient against central differences of the prediction; C_x in span(W₁); a ridge g(wᵀx) recovered at d = 1 with fewer basis functions and lower held-out error than the native fit; the paper's g at d = 1 to 6, against its 91.06% and 98.37% at d = 1 and 2 (exploration, not validation) |
| 7c. The comparison | The report above, with cross-fitting opt-in | the criteria against their direct computation; cross-fitting's folds disjoint and its frames from the other folds only; a function with no low-dimensional structure, where the refit does not win |

**The cost row:** native against projected fit time and M, and the pushforward readout per draw.

## Gates common to every slice

- **Fitting untouched.** E8 reads fits, so no draw moves: the pinned tier at each commit, the whole assembly at the slice's end ([targeted runs](../../../CONTRIBUTING.md#testing-standards)).
- **Truth from outside the model.** Recovery functions are not drawn from BASS's prior. Their subspaces are known in closed form.
- **Primitives.** Chosen cases and a seeded random sweep, far ranges included.
- **The cost row.** Derivative moments per pair, C per draw on each route, the readout at a few thousand draws. Cycle clock, marked indicative.
- **`docs/Maths`.** The decisions above land on the sensitivity page, or a new page, with their alternatives and implications. New writing takes the module's symbols, and crossings go in `issues/notation-ledger.md`.

## Not in this brief

- **Rotate-and-refit, BPPR and a directional hinge.** They change the model, not a readout. They wait for the user. Rotate-and-refit can take its frame from the ordered eigenvectors of P̄ or E[C], which needs no Stiefel code; Stiefel enters only if a rotation within the span is optimized for a simplicity criterion, as varimax is (the Grassmann discussion, settled 2026-11-08). BPPR's direction prior is still open there.
- **A robust median of subspaces as the default,** per decision 4.
- **C for BGP** in closed form (Wycoff et al., not held). BGP's gradients are a later consumer.
- **G4's joint mixture of products** and E25's readouts (ALE, Shapley).
- **The response envelope and E34.** E34's co-active pairs across outputs use slice 3 once E34 exists.
- **`StiefelManifold`'s repair,** which waits for a frame consumer.

## Working rules

- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../../CLAUDE.md)). Leave `basin-hold-validation.csv` and `path-profile-validation.csv` out of every commit.
- Other sessions may be writing in the checkout; the Grassmann discussion edits only its own brief. Before editing a file, check `git diff -- <file>` for changes you did not make; build from a snapshot (`git archive HEAD` plus your files); commit only compiling states. Shared docs (`docs/Maths/BayesianRegression.md`, the notation ledger, `extensions-brief.md`) are edited and committed in one step.
- Builds set `ArtifactsRoot` inline into `artifacts\chip-e8\`, deleted at the end.
- Heavy runs at below-normal priority with capped threads, MSBuild node reuse off; the machine is shared. Tests run by [targeted runs](../../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`, `-SourceRoot` at your snapshot): the selected classes per edit, with the pinned tier per commit, the whole assembly once at the slice's end.
- **Sources:** read a paper before stating what it says, and mark anything cited from memory. When a source is not on disk, ask the main session; it asks the user, who may add it. Fetch nothing yourself.
- **The Grassmann discussion** (`../discussion/grassmann-stiefel-subspaces-brief.md`) is the user's. Read its report before slice 2, and take its conclusions where the user has settled them. Raise a conflict with this brief's decisions with the main session; don't resolve it on your own.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.
- A follow-up goes in a report with the decision it would change and what each outcome would do.
- After each slice, report to the main session and wait for its verification note before starting the next.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; decisions taken beyond this brief; what was left unbuilt and why; follow-ups with the decision each would change.

## Reports

### Slice 6: the readout follow-ups (2026-11-10)

**Landed:**
- `276af70`: items 2 to 6, and the core that pairs a basis fit with another fit.
- `e9d8827`: item 7.
- `30df59e`: `BassFit`'s overload, committed alone through the index, since the file held the weighted-fits session's uncommitted hunk, which stays in the working tree.
- `9a669b6`: `BassBasisFit`'s overload and the pair gates.

**What was built:**
1. *Pairs with basis fits.*
   - `BassFit.CoActiveSubspace(data, BassBasisFit other, spec)` pairs a scalar or direct functional fit with a basis fit.
   - `BassBasisFit.CoActiveSubspace(BassBasisFit other, spec)` pairs two basis fits.
   - Both run through one internal core over two sides, each a direct fit or a basis fit, in the terms-and-H form: draws by position, product weights. The direct pairing of slice 3 now takes the same core, with the same calls in the same order.
2. *The co-active projection mean.* `CoActivePosterior.ProjectionMean(q)` returns a `CoActiveProjection`:
   - each pair's top-q co-active subspace by |λ|, and P̄ over the pairs with its eigensystem;
   - each pair's q selected eigenvalues, signed, and its angles to the mean subspace;
   - along each mean direction u, the directional concordance uᵀE[V]u/√(uᵀE[C_f]u·uᵀE[C_g]u), in [−1, 1] by Cauchy–Schwarz.
3. *Output locations.*
   - `OutputMetric.Locations` takes columns by output-variable name, as prediction's `OutputLocations` does. A direct fit's C, output patterns and checks are then read there.
   - Refused: locations outside the output variables' boxes, and locations asked of a basis fit or of basis adequacy.
   - The sampled route reads the gradient there through the internal `BassGradients.Of`, which gains the locations.
4. *The matrix check's Grams cached.* Σ_B and K_D are kept per structure within a readout, keyed by the terms' identities.
   - It is a pure cache: a constant budget of 2²² values, filled in draw order, never evicted.
   - The check stays off by default (below).
5. *Poincaré constants for truncated normals.* A mixture of one normal component, N(μ, σ²) truncated to [a, b], takes min(σ², (b − a)²/π²), in both the per-input check and the matrix check.
   - σ² is Bakry–Émery's bound for a law with potential curvature 1/σ²; (b − a)²/π² is the weighted Payne–Weinberger bound for a log-concave law on an interval. Both are cited from memory and checked numerically below.
6. *Activity scores.* `ActiveSubspaceAnalysis.Activity` (`ActivityScores`), read when asked from the posterior and the ranks:
   - C's diagonal per draw and its weighted mean, always;
   - the truncation at k, per draw and for E[C]'s own eigensystem;
   - under a rank rule, the truncation at each draw's k_s.
7. *Legendre moments under a uniform law.* `LegendreDesignFactor` skips the monomial expansion:
   - on [0, 1], closed forms;
   - on another interval, the Gauss–Legendre rule exact for the integrand's degree, each factor by its recurrence;
   - its error estimate is 0. Other laws keep the expansion.

**Gates:**

| Gate | Reading |
| --- | --- |
| Pairs: C_f, C_g, C_fg per pair against tensor quadrature of the two gradients (the direct fit's readout; the basis fit's reconstructed from its components'), split at the knots | direct with basis 1.3·10⁻¹⁴; two bases 1.6·10⁻¹⁴ |
| Pairs: the routes, direct with basis, 8 pairs at 8,192 points | within 2.53 standard errors |
| A basis fit paired with itself against its two outputs read together | every pair's matrices the same bits |
| The projection mean on 20 pairs with known co-active subspaces and signs (p = 5, groups of 70% and 30%) | P̄'s eigenvalues (1, 0.7, 0.3, 0, 0) to 2·10⁻¹⁵; the mean subspace to 6·10⁻¹⁶ rad; the concordances 1 and −0.908383, as the closed form |
| Opposed pairs, g = −f; weights against replicates | P̄ equals the C_f draws' projection mean to 1.1·10⁻¹⁶; every concordance −1; weighted against replicated 1.1·10⁻¹⁶ |
| The training locations declared as the metric's locations | C and the patterns the same bits, on both routes |
| Nine new locations, one of them a training location, against tensor quadrature of the gradient readout there; W = I and a full W | 4.1·10⁻¹⁵ and 6.8·10⁻¹⁵; tr E[G] = tr E[C] to 2.7·10⁻¹⁵; the training location's partials the same bits |
| New locations, the routes | within 0.97 standard errors |
| The Gram cache: every draw's ratio read alone against among every draw | the same bits; scalar 60 draws over 3 structures, functional 60 over 5 |
| Truncated normal at σ = 0.08 and 2 of the box on [0.1, 0.9]: Bound/C_vv | the constant taken, σ² and (b − a)²/π², to 2·10⁻¹⁶; no violation in either check |
| The Gaussian equality: a hinge linear on the support under N(0.5, 0.05²) truncated at eight sds | Var/(σ²E[φ′²]) = 1 − 8.1·10⁻¹⁴; at 0.95σ², 1.053 (fails) |
| The cited bounds on 300 seeded truncated normals and random smooth functions; their two ends | worst 0.99935; linear under a wide truncation 1.00000000, cos under a nearly flat normal 0.99999968 |
| Laws with no constant here: two normals, a Student t, the empirical marginal | not applicable; the per-input ratio NaN |
| Activity scores on known spectra, k = 1 to 4, with a trace-share rule | to 8.7·10⁻¹⁵; the mean is E[C]'s diagonal bit for bit through the readout |
| Legendre on [0, 1] to degree 40 | the means and pairs exactly; the slopes 1.3·10⁻¹⁶ from the closed forms, 4.3·10⁻¹⁴ from 64-node quadrature |
| Legendre on [0.1, 0.85], [−0.2, 1.3] and [0.4, 0.45] to degree 30 | within 3.1·10⁻¹³ of a wider rule; within 1.3·10⁻¹⁴ of the expansion to degree six, which at degree 30 is 42% off on [0.1, 0.85] |
| C of three Legendre terms to degree 25 against a 30 × 30 tensor rule of the recurrence | within 10⁻¹⁴ of C's largest entry; the service's error estimate 0 |
| Fitting untouched | the whole regression assembly 3,049 of 3,049 on `9a669b6` |

**The matrix check's cost (item 4), for the user's decision on its default.** The machine was loaded throughout, so the readings are on the thread-cycle clock and indicative (identical binaries differ up to 2×). Per draw at 3,000 draws, the exact readout with the check and without, µs:

| Build | Scalar | Functional | Basis |
| --- | --- | --- | --- |
| Cached, run 1 | 19.5 / 13.4 | 58.5 / 31.0 | 72.9 / 19.9 |
| Cached, run 2 | 31.8 / 25.4 | 49.2 / 44.1 | 110.2 / 24.2 |
| Before the cache (`cb4849b`), interleaved | 44.5 / 20.5; 36.8 / 27.7 | 46.6 / 56.5; 53.4 / 54.0 | 135.9 / 34.9; 142.9 / 34.6 |
| Cached, interleaved | 31.3 / 22.7; 15.9 / 11.8 | 65.4 / 52.2; 56.6 / 24.9 | 135.8 / 29.6; 69.1 / 18.7 |

- On the scalar fit the check's added cost fell from +24 and +9 µs to +8.6 and +4.1 µs.
- The functional readings are below the noise: two of the pre-cache pairs read the check as costing less than nothing.
- On the basis fit it stays about +50 to +108 µs either way.
- So the Grams were the scalar check's cost. The functional and basis checks' cost is the per-draw m × m generalized eigenproblem over the output locations (15 for the basis fixture).
- With the cache the check still costs 1.3 to 2 times the readout on scalar and functional fits, and 3 to 5 times on the basis fit. It therefore stays off by default, and the default is the user's call with these numbers.

**Decisions beyond the brief:**
1. *Pairs: the direct fit with a basis fit on `BassFit`, two basis fits on `BassBasisFit`.* A basis fit first with a direct fit is the first with the roles swapped: V and the concordances are the same, and C_f and C_g trade places. *Alternative:* that ordering on the basis class, leaving two basis fits unpaired.
2. *q is the call's,* `ProjectionMean(q)`, as `CoActivityScores(q)` takes it, so the spec is unchanged. *The signs beside the mean:* each pair's selected eigenvalues, and the directional concordance along each mean direction. *Alternative:* each pair's concordance alone, which the posterior already holds.
3. *Locations outside the output variables' boxes are refused,* as a joint sample's points are. Prediction allows them with a flag. The co-active output functionals stay over the training locations.
4. *The cache lives for one readout,* so nothing is kept between calls. Its budget is counted in values, and its keys are term identities, so equal terms from other instances only miss.
5. *The constant is built for one normal component alone.* A mixture of several is not log-concave in general, a Student t is not log-concave, and an atom has no Poincaré inequality. Slice 5's gate that pinned one normal component as not applicable now uses a mixture of two normals.
6. *Activity scores are a property of the analysis,* computed from the posterior and the ranks, so the record's shape is unchanged. The truncation at k is always given, since every readout has a k.
7. *The brief's "exact rationals"* are read as closed forms, with the standardization's square roots, on [0, 1], and as an exact rule on other intervals. The change moves the low bits of every Legendre moment under a uniform law: Sobol indices, effects and C alike.

**Left unbuilt:**
- A public `outputLocations` on `BassFit.Gradient`. The internal path carries it; the public file was in another session's hands.
- Co-active output functionals at new locations.
- Constants for several components, Student t and empirical laws.
- The matrix check's reduction to K × K for a basis fit.

**Follow-ups:**
- *The matrix check's default.* If on, every readout pays the costs above; if off, it stays opt-in as now.
- *A K × K matrix check for basis fits.* A basis fit's outputs are G·w, so the generalized eigenvalue of (AΣAᵀ, AKAᵀ) with A = RG of full column rank K is that of the K × K pair (Σ, K). If built, the basis check costs about what the scalar's does, which could change the default question. If not, it stays at 3 to 5 times the readout.
- *A public `outputLocations` on `BassFit.Gradient`.* If wanted, one optional parameter, which the internal path already takes. If not, the readout alone reads gradients at new locations.

**Doc changes owed:**
- `docs/Maths/Regression/ActiveSubspaces.md`: a slice-6 section with the seven items, their alternatives and the cost table above; the matrix check's paragraph (the truncated-normal constant, the cache); `OutputMetric`'s locations.
- `docs/Maths/Regression/PolynomialChaos.md:41`: moments under a uniform law no longer go through the expansion, on [0, 1] or elsewhere.

### Slice 5: the Grassmann discussion's settled items (2026-11-08)

**Landed** on `main`, nothing pushed:
- `def59f3`, per-draw ranks, the output side, the matrix Poincaré check, basis adequacy, the angles through LinAlg, the gates and the cost row.
- The docs commit carrying this report: the active-subspace page's slice-5 section and the notation ledger's entry.

The slice was commissioned by the Grassmann discussion's session (`c5859ad`), from the user's settlement of 2026-11-08 ("ill take as recommended … ill go with per-draw k per your recommendation"). The derivations are that brief's synthesis and its "settled by the user" entry (`07c60ec`). LinAlg's `PrincipalAngles` landed in `0b4ed4e` before the slice reached the delegation step, so nothing waited.

**What was built:**
- **Per-draw k:** `ActiveSubspaceSpec.Rank` takes a `RankRule` (trace share, or the gap against λ₁ up to k_max), with no default. It reads `ActiveSubspaceAnalysis.Ranks`, a `RankPosterior` with:
  - each k_s and their distribution;
  - P̄ over mixed ranks, with tr P̄ = E[k_s];
  - the rank summary (P̄'s eigenvalues above ½) and its loss, with `Loss` for any basis;
  - each draw's rank term and containment defect.
  `DrawSubspace` now keeps all its eigenvectors (`Vectors`) for the ranks to read.
- **The output side:** `ActiveSubspaceSpec.OutputPatterns` reads G = E[J W_x Jᵀ] for direct and basis fits under the output metric's symmetric root, on both routes, as an `OutputPosterior`. The exact route reads it through `BassMoments.GradientGram`; the sampled route through stacked Jacobians. Neither forms G.
- **The matrix Poincaré check:** `ActiveSubspaceAnalysis.MatrixPoincare`, from the moment route (`BassMoments.CovarianceGram`, `GradientGram` under D), when every input is continuous and uniform. Opt-in by `ActiveSubspaceSpec.MatrixPoincare`, off by default, with `MatrixPoincareStatus` saying whether it was read, not asked, or not applicable (`a54a537`, at the verification's request). The fit-agnostic primitive is `MatrixPoincare.Ratio`.
- **Basis adequacy:** `BassBasisFit.Adequacy(direct, data, spec)`, which needs the same predictors and locations. The fit-agnostic primitive is `BasisAdequacy.Of`.
- **The angles:** `SubspacePosterior.PrincipalAngles` delegates to `Maths.LinAlg.PrincipalAngles.Between` and accepts bases of different ranks; an empty basis gives no angles.

**Gates:**

| Gate | Result |
| --- | --- |
| The brief's: without a rule, the fixed-k readout bit for bit as before, apart from the angles. With and without ranks and patterns asked, both routes | every draw's C, eigenvalues, basis and angles, E[C], P̄ and their eigensystems identical |
| The rules on chosen spectra (trailing zeros, k = p, negative rounding); refusals | as stated |
| The brief's: sums of r ridges recover r. BASS fits of a ridge on [0, 1]³ and two on [0, 1]⁴, N = 900, trace share 0.95 | summary rank 1 and 2, every draw at r; at 0.99, one more; the gap rule, 1 for both |
| The brief's: tr P̄ against the weighted mean of k_s; the rank summary unbeaten by random projections of every rank. A seeded sweep of 150 posteriors, with the summary's loss and each draw's distance against direct Frobenius norms | 1.8·10⁻¹⁵; never beaten (margin ≥ 0); 1.2·10⁻¹⁴ and 7.6·10⁻¹⁵ |
| The brief's: tr G = tr C under the same metrics at W_y = I. Here every metric, W_y = I and full, both routes, direct and basis fits, 96 draws | 5.5·10⁻¹⁴ |
| The brief's: a basis fit's G against the direct computation on its reconstructed outputs; a direct fit's against its own readout's; E[G]'s spectrum | 7.0·10⁻¹⁵; a basis fit's G of rank at most K |
| The brief's: the matrix check with equality on cos(πx), failing a deliberately mis-scaled D | (cos πx₁, cos πx₂): 1 to 10⁻¹⁵; D halved: 2; a mixed field 0.94; mass outside the bound's range: ∞ |
| The matrix check on 280 fitted draws, asked (scalar, direct under a full W_y, basis); by default; asked where it does not apply | no violation, worst 0.867; not asked; not applicable under a categorical input, a normal law, a sample |
| The brief's: adequacy zero inside the span, the known share outside. Chosen factors | 0; b²/(a² + b²) to 10⁻¹⁵; containment 0 and 1 as constructed |
| Adequacy on fits: a(x)φ₁ + b(x)φ₂ + noise, a direct fit beside PCA bases of rank 1 and 2, against the truth's share for each basis's own span | 0.2251 against 0.2251 (span(φ₁) itself 0.2296); 0.0011 against 0.0001 |
| The brief's: slice 2's angles within arccos's error of the new ones above 10⁻⁶ rad. 300 random pairs of mixed ranks; rotations by 10⁻¹² to 10⁻⁶ rad | within 0.11 of arccos's own error bound; the small angles to 2·10⁻¹⁶ relative, where arccos returns 0 below 10⁻⁸ |
| Slices 0 to 4's E8 classes, rerun | unchanged |
| Fitting untouched | whole regression assembly 2,565 of 2,565 on `15eb151` plus the slice; the E8 oracle classes 10/10 |

**The cost row** (thread cycles, indicative), per draw at 3,000 draws:
- the exact readout with the matrix check and without it: scalar 20.4 and 8.6 µs, functional 35.9 and 22.9 µs, basis 59.2 and 19.1 µs (a third reading: 13.4 and 7.1, 28.9 and 15.6, 58.2 and 31.6; a fourth, on the wall clock with the final code: 16.6 and 12.9, 28.4 and 14.1, 65.2 and 30.0); per-draw ranks add less than the readings resolve (0.3 µs in the first row, −2.8 µs in the fourth);
- a functional fit's exact readout, 25 µs, and 43 µs with two output patterns;
- the sampled output side's own cost, 2.9 ms at 4,096 points;
- basis adequacy, 0.16 ms per direct draw at 12 locations.

The first cost row read the scalar readout at 8.8 µs with the check and had no baseline without it; measured beside one, the check multiplies the exact readout by 1.3 to 3.1 over four readings, the added cost never small against it, so it is off by default everywhere, as the verification asked for where it doubles the cost.

**Decisions beyond the brief,** each with its alternatives:
1. *The gap rule measures each gap against λ₁, over k from 1 to min(k_max, p), with λ_{p+1} = 0.* The first form tried, (λ_k − λ_{k+1})/λ_k, read a trailing zero eigenvalue (an input no term holds) as a gap of 1. It took rank 3 on two ridges over four inputs, and any ratio or log-scale gap shares that failure. Against λ₁ the rule favours the first gap under a geometrically decaying spectrum: two ridges with λ₂ = 0.37λ₁ read as rank 1. *Alternatives:* a log-ratio gap with a floor, which then needs a floor; only the trace share, which needs a share.
2. *The ridge gate reads the trace share at 0.95.* At 0.99 each fit's off-span gradient energy (a few percent, as slice 2's recovery found) adds a rank. Both are printed.
3. *A draw whose C is zero has k_s = 0,* the empty subspace, its distance its rank term.
4. *G's output metric enters by its symmetric root,* so patterns stay over the locations. *Alternative:* slice 3's eigen-root rows. The trace equals tr C either way.
5. *Output patterns are opt-in by count, and so is the matrix check,* off by default because it about doubles the exact readout on scalar, functional and basis fits alike: the user's standing rule that a payload which makes a path expensive becomes opt-in for that path and keeps its output. Slice 2's per-input check stays on. Both checks read the moment route on either route, so they carry no Monte Carlo error.
6. *The matrix check applies only when every input is continuous and uniform.* A categorical input's variance has no gradient bound, and a Poincaré constant for other laws is not built. *Alternative:* known constants for log-concave laws (a truncated normal's is at most σ², cited from memory), a follow-up if wanted.
7. *Adequacy takes its pattern count from `OutputPatterns`,* 1 when 0, and needs the same predictors and the same number of locations.

**Left unbuilt:** Poincaré constants beyond the uniform law; patterns at locations beyond the training ones (slice 3's follow-up, still the user's).

**Follow-ups:**
- *Poincaré constants for log-concave laws* would extend the matrix check to truncated normals. If wanted, D takes σ² per such input. If not, the check stays with uniform laws.
- *The matrix check's cost.* Its term-pair Grams (Σ_B, K_D) depend only on a draw's terms, which repeat across draws, so caching them per structure would leave it a sandwich and a small eigenproblem per draw. If that is built and measured small against the readout, the check could default on again; if not, it stays opt-in.

### Slice 4: dependent inputs (2026-11-08)

**Landed** on `main`, nothing pushed:
- `20aeb06`, the joint-sample law, whitening by its full covariance, the basis fits' sampled route, the unit-scale eigendecompositions, the gates and the cost row.
- The docs commit carrying this report: the active-subspace page's slice-4 section, decisions 1 and 2 marked built, and the notation ledger's entry.

2307.11241 was read again for its dependent laws. §3.3.2 whitens a multivariate normal and refits the emulator on z. §3.3.3 approximates a dependent law by a mixture of products, as in §4.1.1's triangle with L uniform boxes and the supplement's normal mixture fitted to calibration draws. It has no joint-sample route. The triangle is this slice's dependent gate, with its C reached here from a fit and a sample.

**What was built:**
- **`JointSample`** (points as prediction cases, optional weights) on **`ActiveSubspaceSpec.Sample`** and **`CoActiveSpec.Sample`**. Under one:
  - `Auto` takes the sampled route at the sample's points, every draw sharing them; the exact route is refused, and so are `Measures` beside a sample;
  - points must lie in the boxes; weights are normalized, with at least two positive;
  - the analyses' `Measure` and `Poincare` are null, and `Points` is the sample's size;
  - scalar, direct functional and basis fits all take it, and so do the three co-active overloads.
- **Weighted standard errors:** Σw²(v − m)²/(1 − Σw²) per entry, slices 2 and 3's formula at equal weights. They are streamed, so no draw holds N matrices.
- **`InputMetric.Whitened`** is L = Σ_ρ^{1/2} by the symmetric root. Under a product law that is diag(σ), unchanged; under a sample it is the full covariance. **`InputMetric.Standardized`** names the diagonal form under either law. `Map` is a full L; congruences take its diagonal entries alone when it is diagonal, in slice 2's product order, so the exact route's C keeps its bits under a product law (the eigensystems move by the rescaling below).
- **`Maths.Regression.Subspaces.Whitening`**, fit-agnostic: `Moments` (weighted mean and covariance, two passes, the sample as its own law) and `SymmetricRoot` (refused below a 10⁻¹² eigenvalue ratio).
- **A basis fit's sampled route,** deferred from slice 3. It reads each component's gradient readout at the shared points, with the output metric's root taken in component space, GᵀW_yG, for C and for two outputs.
- **`ScaledEigen`:** every E8 eigendecomposition (each draw's C, P̄, V, W_y's root, whitening's root) now runs at unit scale with a 10⁻¹⁴ threshold. `DenseEigen` stops rotating below an absolute 10⁻¹², which left the root of a small-entry covariance 1.1·10⁻⁶ from squaring back. The equivariance sweep found it at 3.8·10⁻³ in QᵀQ − I; after the change, 6.8·10⁻¹⁵ and 1.8·10⁻¹¹.
- Slice 3's two refusal tests no longer assert that a basis fit's sampled route is refused.

**Gates:**

| Gate | Result |
| --- | --- |
| The brief's: against a product law the route already gives exactly. A joint sample drawn from it (truncated normal, uniform, uniform levels) on a scalar fit with a categorical input, 16,384 points, against the exact route | every entry of 12 draws within 2.45 standard errors |
| The same for a direct functional fit under a full W_y; a basis fit on its sampled route by the law's points and by a sample; two scalar fits' pairs; a basis fit's two outputs both ways | within 1.73; 0.30 and 1.14; 1.81; 1.86 and 1.54 standard errors; the basis checks without a violation |
| The brief's: equivariance under affine input maps with whitening. A seeded sweep of 200 (p to 5, 50 to 400 points, weighted in half; reference gradient fields; random maps, scales 10⁻³ to 10³) | whitened eigenvalues 4.3·10⁻¹³; QᵀQ − I 1.8·10⁻¹¹; C_u′ against QᵀC_uQ 6.2·10⁻¹⁴; the box spectrum moves by up to 1.7·10⁸ |
| Whitening in the readout: the map against Denman–Beavers' root of the weighted sample covariance; each draw's C against L̂ᵀC_boxL̂; standardizing against diag(σ̂); weights against repeats | 9.6·10⁻¹⁵; 1.6·10⁻¹⁴; 1.1·10⁻¹⁴ |
| The whitened spectrum in the inputs' own units (boxes [−3, 5] and [10, 12], a correlated sample) against L_rawᵀC_rawL_raw | 4.9·10⁻¹⁵; box and raw differ by up to 303% |
| The primitives: roots of covariances to condition 10⁸ and scales 10⁻⁶ to 10⁶; weighted moments against repeats; refusals | 6.8·10⁻¹⁵, exactly symmetric; 6.3·10⁻¹⁸ |
| The brief's: a dependent example where the box metric and whitening differ, both reported. The band \\|x₁ − x₂\\| < 0.2, correlation 0.91, f = 16(x₁ − x₂)² + (x₁ + x₂)²: truths by quadrature; a BASS fit on 500 Latin hypercube points read under 20,000 points of the law | box leads with x₁ − x₂, whitened with x₁ + x₂; the readouts 0.036 and 0.051 rad from their truths, 1.48 rad apart |
| 2307.11241's triangle law: its C_ρ₂ against quadrature; a BASS fit read under 20,000 points | closed form to 10⁻¹⁴; E[C] within 0.38% (box) and 0.36% (whitened); leading direction 3·10⁻⁴ rad from the truth's; the two metrics' leading functionals 1.7·10⁻⁴ rad apart, C being near rank one |
| The refusals: a sample beside Measures, the exact route, points outside the boxes, bad weights, a missing column, whitening a degenerate sample, standardizing a constant input | each refused; the box metric reads the degenerate sample |
| Slices 1 to 3's E8 classes, rerun (the readout restructured, the eigendecompositions rescaled) | unchanged at their tolerances |
| Fitting untouched | whole regression assembly 2,530 of 2,531 on `e82e8b5` plus the slice, the GBP test skipped as at `e82e8b5`; the E8 oracle classes 10/10 |

**The cost row** (thread cycles, indicative): under a joint sample, 0.53 ms per draw at 4,096 points and 1.8 ms at 16,384, whitened alike. The basis fit's sampled route takes 0.92 ms per joint draw at 4,096 points, and whitening at 16,384 points 0.23 ms for 2 inputs and 2.4 ms for 10. Slice 2's sampled readout read 0.68 ms in the same run.

**A finding on where the metric matters.** An exactly rank-k C's top-k subspace is the same set of input functionals in every invertible metric: LᵀCL has range LᵀG, and L⁻ᵀ maps it back to G. A metric decides only among the directions a C of higher rank leaves. So the paper's triangle, near rank one, cannot separate the metrics, and the band law, of full rank, separates them by 90°. The page records it.

**Decisions beyond the brief,** each with its alternatives:
1. *The symmetric root* for whitening. *Alternatives:* a Cholesky factor, which depends on the inputs' order; PCA whitening, which loses the input axes. All three give the same eigenvalues.
2. *The sample is its own law,* its covariance divided by the weights' sum. *Alternative:* N − 1, a population estimate, which the sample is not here.
3. *Points outside the boxes are refused,* the region a product law's readout covers. *Alternatives:* extrapolating with a count reported; truncating silently, which changes the law.
4. *`Standardized` is added; `Whitened` keeps its name,* with its product-law values unchanged and the full covariance under a sample. *Alternative:* renaming the diagonal form and moving `Whitened`'s meaning, which would change slice 2's readout.
5. *`Measure` and `Poincare` are nullable* in the analyses, null under a sample. *Alternative:* NaN-filled checks, which would read as checks that ran.
6. *A joint sample may carry categorical columns,* a dependent law over the levels too; C still takes no part in them.
7. *Unit-scale eigendecompositions in E8's code* rather than a change to `DenseEigen`, whose other consumers would move.

**Left unbuilt:** an exact route under a dependent law, by G4's joint mixture of products or 2307.11241's whitening transform with a refit (not in this brief); the slice-3 follow-ups the user holds.

**Follow-ups:**
- *`DenseEigen`'s absolute threshold.* It reads 10⁻¹² against the matrix's entries, not their scale, so any consumer that decomposes a small-entry matrix gets eigenvectors good only to that threshold over the eigengap. E8 works around it. If LinAlg adopts a relative criterion (Demmel and Veselić's, |a_pq| ≤ tol·√|a_pp a_qq|, cited from memory), `ScaledEigen` retires and every consumer gains; the decision is LinAlg's owner's, with a bit-for-bit check of what moves. If not, E8 keeps its wrapper and the hazard stays for others.
- *Rotate-and-refit, BPPR and a directional hinge* wait for the user after this slice, per the brief.

### Slice 3: functional outputs and co-active pairs (2026-11-08)

**Landed** on `main`, nothing pushed:
- `84fe230`, the readouts, the posterior, the gates, the oracle's sections and the cost row.
- The docs commit carrying this report: the active-subspace page's slice-3 section, four rows in BayesianRegression's corrections, the notation ledger's entry and the R README's oracle row.

The papers and the oracle were read as slice 0 read them. 2311.18146 v2 was read from its TeX source; `coactivity` at `35e4573` was read and run for every comparison below; `concordance` at `b34e7c6` was read for the paper's script and not run. Neither paper defines a functional output's C: 2307.11241 is scalar, and 2311.18146 pairs scalar functions, so C = E[JᵀW_yJ] follows the brief and F-X6. Zahm et al.'s vector-valued form (1801.07922) is not held, and nothing here cites it.

**What was built:**
- **Every readout reads a draw as its terms and H,** each term's coefficient at each output location: one row for a scalar fit; c_a·h_a(t_ℓ) for a direct functional fit, its output-side factors and normalizer at the training output locations; G_ℓk·c_ka for a basis fit over its joint draws, G = diag(y_s)K. Each quantity below is then a quadratic form in the term pairs' gradient moments, summed by the new **`BassMoments.GradientCross`** for two term lists and a weight matrix.
- **Functional C:** `ActiveSubspaceSpec.OutputMetric`, an **`OutputMetric`** W_y over the training output locations (`OutputMetric.Diagonal` for weights), the identity by default. C = Σ M_ab E[∇B_a∇B_bᵀ], M = HᵀW_yH. W_y is refused unless m × m, finite, symmetric and positive semidefinite, on either route. `BassFit.ActiveSubspace` now reads a direct functional fit; the sampled route maps each point's Jacobian by W_y's root. **`BassBasisFit.ActiveSubspace(spec)`** reads a basis fit on the exact route; its sampled route is refused until slice 4.
- **The checks** take their totals as Σ M_ab·T_ab(v), the scalar readout included (M = bbᵀ), in place of slice 2's call to the Sobol readout.
- **Co-active pairs:** **`BassFit.CoActiveSubspace(data, other, otherData, spec)`** for two fits on the same inputs at one degree, scalar or functional, draws paired by position with product weights; **`BassFit.CoActiveSubspace(data, spec)`** and **`BassBasisFit.CoActiveSubspace(spec)`** for two outputs of one fit, each draw paired with itself. **`CoActiveSpec`** carries the law, metric, route, points, seed, both draw selections, and the output functionals `Output` and `OtherOutput` (a_f, a_g) for functional sides. **`CoActiveAnalysis`** returns the pairs' posterior and, on the sampled route, every entry's standard error (**`CoActiveError`**).
- **`Maths.Regression.Subspaces.CoActivePosterior`**, fit-agnostic beside `SubspacePosterior`. It holds E[C_f], E[C_g], E[C_fg] and V, with V's eigenpairs ordered by |λ|; the contributions, summing to the concordance of the means; each pair's concordance; and the signed and unsigned co-activity scores. **`CoActivePair`** holds one pair's matrices, traces, concordance and discordance.
- **The oracle's sections** (`active_subspace_oracle.R`): `Cfg_bass` on draw pairs as written, with slice 0's corrections patched into a copy of `Cfg_bass_univariate`, and the quadrature reference; `Cfg_bass`'s default and dispatch; a functional draw read per location through `bassfunc2bass_fixed_t` and corrected `C_bass_univariate`, beside the quadrature of the draw folded at t; the functional converters' readings. Row 3's last column now says what a draw is for, so draws sent only within pairs skip the single-draw section and the run stays inside the oracle's two minutes.

**Gates:**

| Gate | Result |
| --- | --- |
| The brief's: a basis fit's C (2 components, 15 locations) against the direct computation on its reconstructed outputs, J = G·J_w by each component's gradient readout, Σ w·JᵀW_yJ by tensor quadrature split at the knots; identity and full W_y | 6.9·10⁻¹⁶, 9.4·10⁻¹⁶ |
| A direct functional fit's C (7 locations) against quadrature of the gradient readout's Jacobian: identity, trapezoid weights, full W_y | 4.6·10⁻¹⁵, 2.0·10⁻¹⁵, 1.6·10⁻¹⁵ |
| W_y's readings: e_ℓe_ℓᵀ against that location's quadrature; linearity; no W_y as the identity | 5.4·10⁻¹⁶ |
| The direct fit's exact and sampled routes under a full W_y, 8,192 points, 12 draws | within 1.22 standard errors |
| The brief's: C_fg's symmetry, two fits swapped; its reduction to C at f = g, a fit paired with itself draw for draw, concordance 1 | 3.8·10⁻¹⁶; 2.0·10⁻¹⁶ |
| Two fits against quadrature of both gradient readouts; exact against sampled (C_f, C_g, C_fg, 10 pairs) | 6.1·10⁻¹⁵; within 2.12 standard errors |
| Two outputs of a direct fit: swapped; a_f = a_g against C under W_y = a_fa_fᵀ; quadrature; the two-fit route with the fit paired with itself; the routes | 2.4·10⁻¹⁶; 6.0·10⁻¹⁶; 8.1·10⁻¹⁵; equal; within 1.82 standard errors |
| A functional fit through a_fᵀy paired with a scalar fit on its inputs, against quadrature | 5.5·10⁻¹⁵ |
| Two outputs of a basis fit: a_f = a_g against C under a_fa_fᵀ; quadrature of the reconstructed a_fᵀJ, a_gᵀJ | 7.8·10⁻¹⁶; 1.1·10⁻¹⁵ |
| The brief's: 2311.18146's examples. The simple polynomial example (§3.1.1) in closed form at β = 0, ½, 3, −12: C₁, C₂, C₁₂ against the paper's matrices and quadrature, the traces, the concordance, π summing to it, the \|λ\| order, the scores | 10⁻¹⁴; every formula holds |
| Its quoted readings | 0.944 (β = ½), 0.551 (β = 3), w₁⁽¹⁾ and w₁⁽²⁾ hold; three slips, below |
| Its estimate (§3.2): BASS fits of f₁, f₂ at β = 3 from n = 200 and 1,000 Latin hypercube points, E[C₁₂]'s Frobenius distance from the closed form | 0.061 and 0.015, 0.87% and 0.22% of ‖C₁₂‖_F; the paper's single fits through `concordance`: 0.050 and 0.0086; the concordance of the means 0.548 and 0.551 (truth 0.551) |
| Oracle (`CoActiveParityTests`): two engine fits' C_fg against `Cfg_bass` on the same 6 draw pairs, and the mode's quadrature | uniform law: as written, 1.2·10⁻¹⁵; a truncated normal on x0: corrected, 1.5·10⁻¹⁵, and 27% from it as written; quadrature 3.3·10⁻¹⁴, 2.3·10⁻¹⁴ |
| Oracle: a direct functional fit's C at a training location (W_y = e_ℓe_ℓᵀ) against `bassfunc2bass_fixed_t` then corrected `C_bass_univariate`, and the quadrature of the draw folded at t_ℓ; 4 draws at 4 of 9 locations | 3.5·10⁻¹⁵; 2.8·10⁻¹⁴ |
| Oracle: `Cfg_bass`'s default and dispatch; the converters' readings and stops | as read, below |
| Slice 1's and slice 0's oracle classes, rerun on the extended mode | 7/7, unchanged; with the co-active class, 10/10 |
| The checks on functional draws (direct and basis, identity and full W_y) and the sampled route | no violation; worst Poincaré ratio 0.85 |
| Slice 2's gates, rerun (the readout's totals now come through M) | unchanged: routes within 1.30, 1.14, 1.30 standard errors; ratios 0.81, 0.21, 0.38 |
| Fitting untouched | whole regression assembly 2,484/2,484 on `e70afaf` plus the slice; the E8 oracle classes 10/10 |

The paper's application (§4, PBX 9501) reads data from its authors' repository (github.com/knrumsey/PBX9501), which is not under `packages` or `codex-scientiae`, so it was not run.

**The cost row** (thread cycles, indicative), per draw at 3,000 draws: a direct functional fit's C over 7 locations, 39 µs (identity) and 44 µs (full W_y); a basis fit's over 15 locations and 2 components, 43 µs; two outputs, 48 µs; per pair, two scalar fits, 15 µs exact and 0.93 ms sampled at 4,096 points. Slice 2's readout in the same run: 10 µs exact, 0.98 ms sampled. Slice 3's first cost run read the sampled route at 2.6 ms: it allocated each point's Jacobian, its W_y-mapped copy and a Gram matrix. It now reuses one buffer and streams the sums, as the co-active pair does. Against slice 2's build (`git archive 3cdc2e8`), interleaved twice under one load, the sampled readout reads the same (1.33 and 1.19 ms against 1.27 and 0.93) and the exact one about a quarter slower (34.7 and 32.5 µs against 27.5 and 25.8). Its checks' totals now run over every term pair through M, where slice 2 read them from the Sobol readout. That is about 7 µs a draw, and it stays. The oracle's reference gradient now takes its factor products vectorized, bit for bit the same; the functional section ran past the runner's two minutes without it.

**The paper's slips,** each against its own definitions; the tests assert the definitions' values and pin each quoted one to what it reads:
- *−0.131 "for β = −12"* is the concordance at β = −15; at −12 it is −0.1095.
- *w₁⁽¹²⁾ ≈ (0.907, 0.422) at β = ½* is V's leading eigenvector at β = 0.4; at ½ it is (0.8971, 0.4418).
- *π⁽¹²⁾ = (0.9518, −0.0077) at β = ½* are the unsymmetrized C₁₂'s eigenvalues over √(t₁t₂). V's, which its eq. (contrib) defines, give (0.9684, −0.0242). Both pairs sum to the concordance, since tr V = tr C₁₂.

Its script (`concordance`'s `inst/CoASM/simple_poly_concordance.R`, read) symmetrizes as the text says, so the quoted values did not come from it as it stands.

**`coactivity`'s readings,** each run or counted in the pinned body:
- *`Cfg_bass_univariate`* shares `C_bass_univariate`'s integrals and prior handling, so its pairs carry slice 0's second-moment and weight defects and read degree one only. With both patched in, it meets the engine.
- *Defaults:* `Cfg_bass` without `mcmc.use` reads one pair, `min(length(s2), length(s2′))`: the shorter chain's last draw with the same index of the other (run).
- *List dispatch:* its third branch tests `class(mod1)` where mod2 is meant, so a list in mod2 alone is refused as not a model (run). Its recursive calls pass `func.use[i], pca_method, verbose` positionally into `func.use, func.true, pca_method` (read; a scalar model ignores them).
- *Functional models* are read per location only, with no output metric. `bassfunc2bass_fixed_t` folds the hinges of variable `pdes + pfunc` alone, the last functional variable (counted). It stops on a one-draw model, since `beta[, 1:(1 + maxBasis)]` drops a one-row matrix to a vector (run). It also stops when the largest model loses no term at t: it trims beta to max(nbasis) + 1 columns, then tests column max(nbasis after merging) + 2 (run, "subscript out of bounds"). With one basis function, `knotInds[i, 1:1, ]` drops to a vector before the duplicate search (read).
- *Basis models:* `bassPCA2bass_fixed_t` and `get_u` place the basis rows at `seq(0, 1, length.out = nt)`, equispaced whatever the output locations are (counted: three and one). `func.true`, documented as the output locations, appears nowhere in the body (counted). `C_bass_pca_v2` defaults to nt points strictly between that grid's, so its default C is read between the training locations; for nt = 15 only ½ falls on the grid (counted and evaluated). `C_bass_pca_v1` and `Cfg_bass_pca_v1` return TRUE (run). `basslc2bass` takes each model's structure at `mcmc.use[mm]` and its coefficients at row `mm` (lines 145 and 231), so a subset `mcmc.use` pairs one draw's terms with another's coefficients (read). `bassPCA2bass_fixed_t` calls it without `mcmc.use`, so the per-location reading is untouched.

The rows in BayesianRegression's corrections take the paper's slips, `Cfg_bass`'s inheritance, the basis converters and `bassfunc2bass_fixed_t`. The defaults and the dispatch are readings, as `C_bass`'s default was in slice 0.

**Decisions beyond the brief,** each with its alternatives:
1. *A functional side is read through an output functional,* f = a_fᵀy, which holds both coactivity's per-location reading (a = e_ℓ) and integrated outputs. Two fits take one each, so two functional fits, or a functional and a scalar fit, pair as two scalar fits do. *Alternative:* per-location matrices only, as coactivity returns them.
2. *Two fits pair draws by position with product weights,* which samples the posteriors' product. coactivity pairs its `mcmc.use` rows the same way. *Alternatives:* every pair of draws, the exact product expectation at quadratic cost; random pairing, which adds noise.
3. *W_y is declared and checked on both routes,* the identity by default (F-X6's Q). *Alternatives:* diagonal weights only, which loses the outputs' correlations; a W_y built from the output grid, which would choose the units for the user.
4. *The co-active sampled route carries standard errors* for every entry of C_f, C_g and C_fg, as slice 2's readout does for C.
5. *The co-active posterior is E[C_f], E[C_g], E[C_fg] and V's eigensystem,* with each pair's concordance beside the concordance of the means. *Not built:* a projection mean over the pairs' co-active subspaces, decision 4's second summary carried over. V's eigenvalues are signed, so a top-q subspace per pair needs the |λ| order to hold across pairs. The question is open below.
6. *The readout reads at the training output locations,* as the gradient readout and the conditioned Sobol form do; a basis fit exists only there.
7. *Notation:* W_y for the output metric, the brief's own symbol, since the page already writes W_s for a draw's basis; a_f and a_g for the output functionals, since b is a draw's coefficients. The splits are logged in the notation ledger.

**Left unbuilt:** a basis fit paired with another fit; a direct fit's output locations beyond the training ones; the basis fits' sampled route (slice 4); the co-active projection mean (decision 5 above); the PBX application (data not on disk).

**Follow-ups:**
- *A basis fit paired with another fit.* Each is a new H in the same terms-and-H form: the pair's draws by position, the basis side over its joint draws. If wanted, it lands as one overload on each class. If not, two outputs of one basis fit and every direct pairing stand.
- *The co-active projection mean.* If the user wants decision 4's second summary for pairs, it needs a rule for q and for signed eigenvalues: top-q by |λ| per pair, P̄ over the pairs, with the concordance's sign reported beside it. If not, E[C_fg]'s eigensystem and the per-pair concordances are the summary.
- *Output locations beyond the training ones for a direct fit.* Its output-side factors are hinges in the output variables, defined across their box, so H can be read at any location. The decision is whether W_y may span new locations. If yes, `OutputMetric` gains its locations; if no, the training locations stand.

### Slice 2: the readout (2026-11-08)

**Landed** on `main`, nothing pushed:
- `8c7db77`, the readout, the summary, the gates and the cost row.
- The docs commit carrying this report: the active-subspace page's slice-2 section, with decisions 4 to 6 marked built.

The Grassmann discussion's report was read first. Its settled readings are taken: P̄'s eigenvalues as mean squared cosines, the projection mean as the Bayes action under chordal loss, and its eigenvectors in order as the nested estimates for every k. Its proposals, not yet settled by the user, are treated as follows. Item 1, a per-draw k, was asked of the coordinating session on 2026-11-08 and is not built; the readout's shape allows it. Item 2, a sine-based angle primitive in LinAlg, is not built either; decision 5's arccos stands. The proposed home for the code, posterior summaries in `Maths.Regression` with no Geometry import, is taken (decision 1 below).

**What was built:**
- **`BassFit.ActiveSubspace(data, spec)`** with **`ActiveSubspaceSpec`**: the laws, the metric (`InputMetric`: box by default, raw, whitened), the dimension k (1 by default), the route (`Auto`, `Exact`, `Sampled`), the sampled route's points (4,096 by default) and seed, and the draws. It returns **`ActiveSubspaceAnalysis`**: the inputs, the excluded categorical predictors, the route that ran, the metric and its map L, the law, the posterior, the sampled route's per-entry standard errors, and the checks.
- **The routes.**
  - *Exact:* the moment service, per draw.
  - *Sampled:* N points drawn once from the product law and shared by every draw, the gradient readout in batches of 64 draws, and each draw's subspace from the thin SVD of its gradient factor F with C = FᵀF.
- **The metric:** C = LᵀC_boxL with L diagonal under a product law: 1 for the box, 1/width for raw units, the law's sd from its interval moments for whitening.
- **`Maths.Regression.Subspaces.SubspacePosterior`**, engine-neutral and needing only LinAlg: E[C] and its eigensystem; P̄, its spectrum and its eigenvectors; each draw's C, eigenvalues, top-k basis and gap; each draw's principal angles to the mean subspace; and the angles between E[C]'s subspace and the mean subspace. `DrawSubspace` reads a draw's subspace from C (`FromC`) or from a factor (`FromFactor`).
- **`PoincareCheck`** (decision 6): totals from the moment route under the same law and draws.
  - *Poincaré* for each input under a uniform law;
  - *SG-4(b)* when every term is a Legendre product and every law uniform on its box;
  - the violations, counted against rounding, or against three standard errors on the sampled route.

**Gates:**

| Gate | Result |
| --- | --- |
| SG-5's three-draw example: the Lean's rational matrices exactly; the summaries' leading directions against the closed forms | 37.99° (rational) and 54.29° (1.2 rad) apart, each to 10⁻¹² |
| Weighted draws against replicated ones: the summary; a fit with `DrawWeights` against `Draws` repeated | 1.9·10⁻¹⁴; 2.1·10⁻¹⁵ |
| Exact against sampled, 8,192 points, on a hinge fit with a categorical input, under the uniform law, a mixture and x0's empirical marginal: 180 entries each | worst 1.30, 1.14 and 1.30 standard errors |
| LᵀC_boxL under raw units and whitening (truncated normal and uniform laws) on a fit whose box is [−3, 5] | 5.4·10⁻¹⁵ |
| Decision 6 on hinge (degree 2) and Legendre (degree 5) draws; SG-4(b)'s equality on a degree-one Legendre fit | no violation, worst ratios 0.81 (Poincaré, hinge), 0.21 (Poincaré, Legendre), 0.38 (SG-4(b)); equality to 3.3·10⁻¹⁶ on 88 draw-inputs |
| Recovery at N = 100, 300, 900 from outside the model: sin(2.5wᵀx) on [0, 1]³, and two sin ridges on [0, 1]⁴ with one inert input; E[C]'s subspace, the mean subspace and E[C]'s relative error | ridge: 0.046, 0.016, 0.010 rad (E[C]'s error 8.3%, 3.6%, 2.4%); two ridges: 0.026, 0.024, 0.0105 rad (6.0%, 5.8%, 2.0%); the mean subspace within 0.001 rad of E[C]'s; within the Davis–Kahan bound at every N; the falls registered in `ComparativeSweep` |
| Fitting untouched | whole regression assembly 2,461/2,461 on `3063af9` plus the slice |

**The cost row** (thread cycles, indicative): the exact readout at 3,000 draws, 15 µs per draw with its checks and 3.9 µs for the summary alone; the sampled readout, 0.63 ms per draw at 4,096 points. Slice 1's row reran beside it: C per draw warm, 5.7 to 6.7 µs.

**Decisions beyond the brief**, each with its alternative:
1. *The posterior over subspaces is engine-neutral,* in `Maths.Regression.Subspaces`, parallel to `Maths.Regression.Conformal`, on LinAlg alone, as the Grassmann discussion proposes. Co-active pairs (slice 3), BGP and BPPR read it too. Alternative: under BASS's sensitivity folder, closer to its one consumer today.
2. *Defaults:* k = 1; `Auto` takes the exact route, which every product law allows today; 4,096 points; seed 1. Alternative: no default k, which would refuse until one is declared.
3. *P̄ takes every draw, rank-deficient ones included,* each with its gap reported. Where a draw's C has rank below k, its k-subspace holds arbitrary directions. Alternative: leave those draws out of P̄, which drops posterior mass; per-draw k (item 1) is the principled answer.
4. *The checks' totals come from the moment route even on the sampled route,* so only the bound carries Monte Carlo error, and the slack is three of its standard errors.
5. *The recovery gate's bound is Davis–Kahan's* (Yu, Wang and Samworth's variant, cited from memory), from E[C]'s own error and the truth's gap, not a fixed angle. Its rank-two function changed from a polynomial of two forms, whose relative gap of 0.0045 makes any subspace gate vacuous, to two ridges with a gap of 0.37.

**Left unbuilt:** per-draw k (asked); the Karcher mean, the median and their scatter, which need Geometry and an orchestration layer; the LinAlg angle primitive (the user's call); functional outputs (slice 3); joint samples and full whitening (slice 4).

**Follow-ups:**
- *Per-draw k* waits for the user to settle it in the Grassmann discussion (the coordinating session, 2026-11-08). The decision is whether the readout reports k_s by a declared rule beside the fixed k. If yes, P̄ over mixed ranks keeps its formula and its trace becomes the mean rank; if no, the gaps stay the signal. The shape it would join:
  - each `DrawSubspace` already holds all p eigenvalues, so a rule reads k_s from them with no second pass;
  - `DrawSubspace.FromC` and `FromFactor` already take a dimension per call;
  - `SubspacePosterior.Of` sums W_sW_sᵀ draw by draw, so mixed ranks change only its shared-k check and the angles to the mean, which become min(k_s, k̄) angles plus a containment defect;
  - `ActiveSubspaceSpec` gains a rule beside `Dimension`.
- *The LinAlg principal-angle primitive.* If the user adopts it, `SubspacePosterior.PrincipalAngles` delegates to it. If not, arccos stands; it is accurate to about 1.5·10⁻⁸ rad, below any posterior spread the gates met (the tightest, 0.01 rad).

### Slice 1: derivative moments and C on the exact route (2026-11-08)

**Landed** on `main`, nothing pushed:
- `95294f8`, slice 0's two `coactivity` corrections recorded as matched expressions, as the coordinating session asked. `replace_lang`, `patch_body` and `patch` move into `_common.R`; `khaos_helpers.R` keeps identical copies, whose removal is a cleanup for later.
- `9dfcc6d`, the derivative moments, C, the gates, the oracle parity and the cost row.
- The docs commit carrying this report: `docs/Maths/ActiveSubspaces.md` (new), a pointer from BayesianRegression.md's sensitivity section, and a notation-ledger entry.

**What was built:**
- **Factor level:** `DesignFactor.SlopeMean` (E[φ′]) and `SlopePair` (E[φ′ψ], E[φ′ψ′]) for hinge and Legendre factors, under the variable's own `ContinuousMeasure`. A categorical factor refuses.
  - *Hinges:* the slopes come through `PiecewisePolynomial.HingeProduct`'s new slope flags. The value path keeps its bits.
  - *Legendre:* the slopes differentiate the expansion about each piece's centre.
  - *The kink convention:* under an empirical marginal, an atom on a degree-one knot takes the derivative's mean of one-sided limits, s/2, as the readout does (F-X3, G3). The pieces, closed on the right, would give it the left limit.
- **Term level:** `BassMoments.GradientPair` gives E[∂ᵢB_a·∂ⱼB_b] over each term's continuous predictors (`TermPairGradientMoments`). `BassMoments.GradientOuterProduct(terms, coefficients)` gives C over `GradientVariables`, the continuous predictors, in transformed units.
  - Slope means, slope pairs and term pairs are cached as the value moments are.
  - Categorical and output factors enter as values (decision 3).

**Gates:**

| Gate | Result |
| --- | --- |
| Hinge factor pairs, degrees one to three, every sign pair, knots in, at the ends of and outside the support, against one-variable quadrature of the factors' own values and derivatives; uniform, normal–Student mixture, and empirical with atoms on the knots | ≤ 2.4·10⁻¹⁵ of the integrand's size |
| Legendre factor pairs to degree ten, same laws, and on [0, 1] against the closed forms of slice 0's reference | ≤ 4·10⁻¹⁴ to degree six, ≤ 7·10⁻¹² to ten (the expansion's cancellation, as the value moments document) |
| Seeded sweep, 400 hinge pairs: knots in [−1, 2], mixture components up to 40 scales outside the interval, knots on empirical atoms | ≤ 2.7·10⁻¹⁴ |
| Term pairs against a tensor quadrature of `BassBasis.CaseDerivative`, categorical factors included, degrees one to three | ≤ 5.1·10⁻¹⁴ |
| C on fitted draws against the gradient readout averaged over the law (rules exact for the integrand): hinge degrees one and two with a categorical input, 60 draws, under the uniform law, a mixture and x0's empirical marginal; a Legendre fit, 80 draws, uniform and empirical | ≤ 8.1·10⁻¹⁴ of C's largest entry |
| SG-4(b) through the slope moments under the Jacobi weight 6x(1 − x): factor pairs; E[x_v(1 − x_v)(∂_v f)²] on the Legendre fit's draws against `PolynomialChaosSensitivity.Jacobi` | ≤ 6.9·10⁻¹⁴ to degree six (2.8·10⁻¹¹ to ten); 5·10⁻¹⁴ on the draws |
| Oracle: an engine fit's draws (degree one, two inputs, up to 9 terms) against `C_bass` on the same draws | uniform law: as written, 8.6·10⁻¹⁶; truncated normal: corrected, 1.1·10⁻¹⁵, and 6.7% from it as written |
| Fitting untouched | whole regression assembly 2,449/2,449 on `95294f8` plus the slice, pinned tier included |

**The cost row** (`ActiveSubspaceCostTests`, behind `BASS_QUALIFICATION`, thread cycles, indicative):

| Law | Factor pair, value moments | Factor pair, slope moments | Term pair, cold | C per draw, cold | C per draw, warm |
| --- | --- | --- | --- | --- | --- |
| Uniform | 2.9 µs | 2.9 µs | 5.0 µs | 21 µs | 7.0 µs |
| Normal–Student mixture | 81 µs | 15 µs | 10 µs | 44 µs | 9.2 µs |
| Empirical, 60 atoms | 1.9 µs | 2.0 µs | 4.8 µs | 23 µs | 8.4 µs |

The value moments under the mixture pay for the Student component's quadrature twice, once for the cross and once for the covariance.

**Decisions beyond the brief**, each with its alternative:
1. *The slope moments take the variable's `ContinuousMeasure`,* not the `ProductMeasure` the value moments take, because only a continuous factor has a slope. This is what lets SG-4(b)'s gate pass the Jacobi weight as a law. Alternative: the `ProductMeasure`, for one shape across the factor API, at the cost of that gate.
2. *Names:*
   - `GradientOuterProduct` for C, the field's term for E[∇f∇fᵀ], where "gradient" alone already names the pointwise readout `BassGradient`;
   - `TermPairGradientMoments` beside `TermPairMoments`;
   - "slope" for the derivative moments, beside the pointwise `Derivative`.
   The splits are in the notation ledger.
3. *The empirical kink correction* matches the readout's convention exactly. Alternative: leave the left limit, and the exact and sampled routes would then disagree on empirical laws by half each knot atom's slope.
4. *A functional fit's output factors* enter C as values over whatever the measure covers, as the value moments treat them. The output metric W_y is slice 3's.

**Left unbuilt:** slice 2's readout (`BassFit.ActiveSubspace`, the sampled route, metrics, summaries). Directional hinges and other non-separable kernels, which have no factors, are still the sampled route's.

**Follow-ups:**
- *Legendre moments at high degree.* The monomial expansion about each piece's centre loses about five digits by degree ten, for the value moments as for the slopes. The decision is whether to expand in the Legendre basis itself under a uniform law, where the integrals are exact rationals. If yes, KHAOS fits of high degree keep 10⁻¹⁴; if not, the documented loss stands and the gates stay banded by degree.
- *`khaos_helpers.R`'s copies* of the patch helpers. The decision is the cleanup: remove them now that `_common.R` holds them. Nothing is lost, since they are identical.

### Slice 0: sources, references and the oracle (2026-11-08)

**Landed** on `main`, nothing pushed:
- `15fb495`, the reference functions and their gates (`tests/regression/ActiveSubspaceReferences.cs`, `ActiveSubspaceReferenceTests.cs`).
- `9e3028d`, the oracle mode and its defect tests (`r/oracles/active_subspace_oracle.R`, `tests/oracle/ActiveSubspaceDefectParityTests.cs`).
- The docs commit carrying this report: `r/README.md`'s row for the mode, five rows in [Corrections to the sources](../../../docs/Maths/BayesianRegression.md#corrections-to-the-sources), and the Sources lines above.

**The paper, read** (v1, main text and supplement). What E8 takes from it:
- C_ij through I₁ = ∫h′_{m₁}h_{m₂}ρ, I₂ = ∫h_{m₁}h_{m₂}ρ and I₃ = ∫h′_{m₁}h′_{m₂}ρ per input (its eq. Cij). That is the moment service's factor pairs with derivative factors, slice 1's work.
- Its MARS case is degree one only: integrands of degree ≤ 2, truncated moments r = 0, 1, 2.
- Theorem 1, C = AᵀC̃A for z = Ax + b, is decision 2's congruence.
- For a correlated Gaussian it refits the emulator on whitened inputs (§3.3.2), a model change. Decision 1's sampled route and decision 2's whitened metric reach dependent laws without a refit.
- Mixtures of products give C = Σ w_l C_l (§3.3.3, G4's form).
- Its activity scores Σ_{j=1}^{p} λ_j w_ij² are C's diagonal exactly (Constantine's truncate the sum at r).
- Posterior propagation: one C per draw, recomputed only when the structure changes.

**Slips in the paper's text, struck for the code** (`coactivity` has each right): the supplement's eq. (Cijv2) swaps the i = j and i ≠ j cases; eq. (I2mars) writes t_{im_j} for t_{im₂}; the Beta moment's exponent is written j for r.

**The oracle, read and run.** Every comparison below read `coactivity` at `35e4573`. `concordance` at `b34e7c6` was read, not run.
- **Parity:** on degree-one draws under uniform laws, `C_bass` meets the mode's reference (each draw's own gradient by tensor Gauss–Legendre split at its knots) to 2.0·10⁻¹⁵ at two inputs and 1.6·10⁻¹⁴ at three.
- **Confirmed, and also in `concordance` by reading:**
  - *No degree is read.* A degree-two draw gets the degree-one C, 13% from its reference.
  - *The truncated normal's second moment flips the lower end's term* (Z₂ = Δ₁ − Aφ(A) − Bφ(B), in the code and in the paper's ξ_TNorm). ξ(2) is off by up to 59%, and patched it meets `integrate` to 10⁻¹².
  - *A normal mixture's components are normalized twice:* weights w/cc, then each component's truncated moment over z again. With the second moment patched, `C_bass` meets the law Σ(w_l/cc)φ_l/z_l to 10⁻¹⁵; with the weights w z/cc as well, the stated law to 2.3·10⁻¹⁵. One component per input inflates C by Π 1/z_k: 1.97 on the fixture, checked against this side's normal CDF.
  - *The gamma prior:* `build_prior`'s gamma reaches `C_bass` as NA, and a hand-set Gamma(2, 3) gives C at its mass on the box (0.80), not renormalized.
  - *maxInt = 1 stops* ("dim(X) must have a positive length"), as the paper's supplement says of `concordance`.
- **Confirmed by reading, not run:** a fit's categorical factors are never read, so a level indicator counts as one.
- **Only in `concordance`, by reading:** its `XI_FUNC_TNORM_SVW` has no component index, so a mixture's per-component moments come back as a vector. `coactivity` added the index.
- **Readings:** without `mcmc.use`, `C_bass` reads the last draw only; `use_native_scale` defaults to false where its documentation says true.
- Each correction is the engine's already or lands in slice 1: the truncated-normal and mixture moments through `Truncated` and `MixtureMeasure`, the degree through the design kernels, the categorical part by decision 3. Nothing is ported.

**The references' gates:**

| Gate | Result |
| --- | --- |
| Closed form against quadrature, on chosen cases: ridges (sin, exp, damped wave, two in four inputs, three in five), two ridges under a correlated Gaussian, a quadratic form, a cubic ridge, an interacting polynomial, a Legendre expansion, and the paper's three cases | ≤ 2.7·10⁻¹⁵ of C's largest entry |
| The paper's benchmark x₁² + x₁x₂ + x₂³/9 on [0, 1]² against its stated C = [[120, 50], [50, 21]]/45 | 8·10⁻¹⁷; first direction (0.9229, 0.3850), the paper's (0.923, 0.385); inert inputs pad it with exact zeros |
| The paper's composite f(t₁ᵀx, t₂ᵀx) on [0, 1]⁶ | range in span(t₁, t₂) to 3.6·10⁻¹⁵; √(λ₂/λ₁) = 0.1099 exact, where the paper's table has 0.125 from a fit |
| The supplement's f₃ = x₁² + x₁x₂x₃ + x₂³/9, which the paper estimates by 10⁶ Monte Carlo draws | exact: [[35/18, 41/72, 5/12], [·, 17/90, 1/8], [·, ·, 1/9]] |
| Gradients against fourth-order central differences | ≤ 4.7·10⁻¹¹ |
| ∫ψ′_mψ_n and ∫ψ′_mψ′_n against Gauss–Legendre, m, n ≤ 12 | 8·10⁻¹⁶ |
| Quadratic form ½yᵀAy against A(ddᵀ + H²/3)A | 2·10⁻¹⁶ |
| Rank one on the direction: a cubic ridge's monomial C, and a tanh ridge's quadrature C | 3·10⁻¹⁶, 2·10⁻¹⁶; the cubic against its quadrature held as a profile 8·10⁻¹⁶ |
| Seeded sweep, 400 cases: one to three inputs, laws centred up to 10⁴ from the origin, widths 10⁻³ to 10² | within 10⁻¹² plus the quadrature's argument rounding at the far end (worst 0.23 of that tolerance); where that rounding is below 10⁻¹³, ≤ 8.6·10⁻¹⁵ |
| Fitting untouched | whole regression assembly 2,429/2,429 on `0a35e5d` plus the slice; oracle commit scope 5/5 |

The cost row has nothing to time in slice 0; it starts with slice 1's derivative moments.

**Decisions beyond the brief**, each with its alternative:
1. *The references are test-side* (`tests/regression`), as the oracles are. Alternative: `src/synthetic`, if an application later wants them as benchmark functions.
2. *A correlated Gaussian beside the box,* with ridge sums in closed form through its MGF exp(½μᵀΣμ). This gives slice 4's dependent example an exact C, where the box metric and whitening differ. Alternative: Monte Carlo there.
3. *Polynomials of several linear forms* (`OfLinearForms`), so slice 2's recovery has rank-r truth from outside the model, the paper's composite among it. Alternative: ridge sums only.
4. *The sweep's tolerance adds a conditioning term,* because quadrature in raw units rounds wᵀx − c at 10⁴. Alternative: keep far laws wide, which would hide the far-range regime.
5. *The mode patches `coactivity`* (`coactivity_corrected`: the second moment and the weights, each matched once in the deparsed body), the way `khaos_corrected` does. The patched code is the reference the engine meets under normal laws; the degree has no patch, since its correction is structural.

**Left unbuilt:**
- `QuadratureRule` has no breakpoints and no atoms. Slice 1's C# hinge quadrature needs both; the R mode already splits at knots.
- The categorical defect was not run, because `C_bass` never reads those fields.
- The Beta prior was not read further: the engine has no Beta measure.
- The paper's triangle law ρ = 2·1(0 < x₂ < x₁ < 1), whose C = [[1710, 741], [741, 322]]/540 was checked by hand, is not built. It is slice 4's dependent example.

**Follow-ups:**
- *`scripts/test-select.ps1` reads a record's positional parameter as a type.* `RidgeColumns? Ridge` in `BassDesign.cs:74` made a test type named `Ridge` select 126 classes instead of 1; the type is now `RidgeTerm`. The decision is whether to fix the selector's name reading, which would cut over-selection wherever a parameter shares a type's name. If fixed, edits like this one select only their own classes. If not, new types should avoid names that are parameter names elsewhere.
- *Activity scores in slice 2's readout.* The decision is whether to report C's diagonal as the activity score, which is the paper's sum over every eigenpair, or to add Constantine's truncation at k. The diagonal costs nothing; the truncation is a second number that depends on k.
