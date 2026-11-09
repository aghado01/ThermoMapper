# Active subspaces on the BASS core: E8's build (chip brief)

Date: 2026-11-08. Status: commissioned by the user on 2026-11-08 ("kick off e8"), after KHAOS's last slice; slices 0 and 1 verified by the coordinating session; slice 2 landed (`8c7db77`) and waits for its verification. Written by the main session from [E8's entry](extensions-brief.md#active-subspace-analysis), [gap G3](lineage-gaps-register.md#g3-sobol-indices-and-active-subspaces-have-different-standing), [the lineage review's separability question](lineage-review.md#separability-is-what-the-closed-forms-rest-on), the [Astra notes](../discussion/astra-grassman-stiefel-pca-extensions.md), and the [Grassmann discussion](../discussion/grassmann-stiefel-subspaces-brief.md)'s report so far. Reports go under "Reports" below, newest first.

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
   - *The arccos floor* (about 1.5·10⁻⁸ rad) is below any posterior spread that matters. A sine-based angle is needed only where angles are differenced.
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

## Gates common to every slice

- **Fitting untouched.** E8 reads fits, so no draw moves: the pinned tier at each commit, the whole assembly at the slice's end ([targeted runs](../../../CONTRIBUTING.md#testing-standards)).
- **Truth from outside the model.** Recovery functions are not drawn from BASS's prior. Their subspaces are known in closed form.
- **Primitives.** Chosen cases and a seeded random sweep, far ranges included.
- **The cost row.** Derivative moments per pair, C per draw on each route, the readout at a few thousand draws. Cycle clock, marked indicative.
- **`docs/Maths`.** The decisions above land on the sensitivity page, or a new page, with their alternatives and implications. New writing takes the module's symbols, and crossings go in `issues/notation-ledger.md`.

## Not in this brief

- **Rotate-and-refit, BPPR and a directional hinge.** They change the model, not a readout. They wait for the user after slice 4, and for the Grassmann discussion: rotate-and-refit's frame is a Stiefel choice.
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
