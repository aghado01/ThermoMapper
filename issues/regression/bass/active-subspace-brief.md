# Active subspaces on the BASS core: E8's build (chip brief)

Date: 2026-11-08. Status: commissioned by the user on 2026-11-08 ("kick off e8"), after KHAOS's last slice; slice 0 starts. Written by the main session from [E8's entry](extensions-brief.md#active-subspace-analysis), [gap G3](lineage-gaps-register.md#g3-sobol-indices-and-active-subspaces-have-different-standing), [the lineage review's separability question](lineage-review.md#separability-is-what-the-closed-forms-rest-on), the [Astra notes](../discussion/astra-grassman-stiefel-pca-extensions.md), and the [Grassmann discussion](../discussion/grassmann-stiefel-subspaces-brief.md)'s report so far. Reports go under "Reports" below, newest first.

## Why

An active subspace is the span of the directions in input space along which a function changes most, averaged over an input measure ρ: the top eigenvectors of C_ρ = E_ρ[∇f ∇fᵀ]. Sobol indices rank inputs and their subsets; C finds linear combinations of inputs, which Sobol cannot. Computed per posterior draw, the subspace comes with its uncertainty.

It serves any smooth surrogate of a costly code or system: reducing a calibration's parameter space, choosing directions to refit along, comparing two models' responses (co-active subspaces), and reading which output patterns move with the inputs. The build order puts it after KHAOS by the user's call of 2026-11-07, and E34 follows it.

## Sources

- **Rumsey et al., *Discovering Active Subspaces for High-Dimensional Computer Models*, arXiv 2307.11241.** The method on BASS: C under tensor-product representations through one-dimensional integrals, BASS's posterior propagated, integrals cached across structure changes. **Not on disk.** Slice 0 asks for it (see the working rules).
- **Co-active subspaces, arXiv 2311.18146,** held in `D:\aghado01\codex-scientiae\supellex-private\gauntlet\bayesian`.
- **Not held:** 1801.07922 (Zahm et al., vector-valued active subspaces), 1907.11572 (Wycoff, Binois and Wild, C for GPs), 2210.09181 (BPPR).
- **The oracle.** None is pinned. The authors' code may be public; slice 0 finds out, and a download waits for the user's go-ahead through the main session. E8 does not depend on one: it consumes what is built, and its internal checks are strong (below).

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
