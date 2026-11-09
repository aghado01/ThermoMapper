# Grassmann and Stiefel geometry for active subspaces and response bases (discussion brief)

Date: 2026-11-08. Status: opened for the user's own exploration; a discussion, not a build. Written by the main session at the user's request ("why dont we chip off so i can dig into it"), from a discussion that followed E8's reminder. Reports go under "Report" below, newest first.

## Purpose

Work through where subspace geometry earns its place in the regression engine: posterior active subspaces (E8), BASS's PCA response basis (D07, `BassDecomposition`), learned projection directions (BPPR, rotate-and-refit), and output reduction for several outputs (E34). The user leads; the session's job is to be a critical partner, pressure-testing the ideas below, scouting the established terms for each idea across fields, and reading sources before restating them. The engine serves many applications, so frame everything by problem class.

## Read first

- `issues/regression/discussion/astra-grassman-stiefel-pca-extensions.md`: the earlier discussion (Grassmann for subspaces, Stiefel for frames, whitening; the PCA-basis observations; the active-subspace placement table).
- E8's entry, "Active-subspace analysis", in `issues/regression/bass/extensions-brief.md`, and E34's ("Several outputs on one basis, with Σ integrated").
- The lineage review's open question on separability: `issues/regression/bass/lineage-review.md`, "Separability is what the closed forms rest on".
- The code: `src/maths/geometry/manifolds/GrassmannManifold.cs`, `StiefelManifold.cs`, `src/maths/geometry/dim-reduction/MoMPCA.cs`, `src/maths/regression/spline/bass/output/BassDecomposition.cs`.
- What KHAOS slice 4 just built: `BassFit.Gradient` (`BassGradient`, pointwise ∂f/∂z per draw) and `BassSensitivity.FromCoefficients` (`8f57487`), and the Lean SG-4(b) identities (`faca090`, `lean/Enthymemata/BassIntegrals.lean`).

## The main session's synthesis, to pressure-test

1. **The posterior over active subspaces lives on Gr(k, p).** Each draw's top-k eigenspace of C = E[∇f∇fᵀ] is a Grassmann point.
   - The projection (extrinsic) mean has a closed form: the top-k eigenvectors of P̄ = Σ w_s W_s W_sᵀ, with SMC or flat-histogram weights where draws carry them.
   - P̄'s eigenvalues lie in [0, 1]. Read as a stability spectrum, they show how many draws hold each direction.
   - The eigenspace of E[C] and the mean of the draws' subspaces differ when the draws disagree, and that difference is itself a readout.
   - By Davis–Kahan, a subspace's uncertainty scales with ‖ΔC‖ over the eigengap, so a small gap means an unidentified subspace. Where the order of directions matters, the flag manifold is the right space, not Gr(k, p).
   - `MoMPCA`'s Grassmann geometric median (after Kisung You), which already down-weights weak eigengaps, is a ready robust summary.
2. **BASS's PCA basis lives on Stiefel, because its components are fitted separately.** Rotating the scores changes how simple the component functions are, so the frame is part of the model, not only its span. Today the frame is a plug-in; its uncertainty enters only through the truncation residuals. A Bayesian frame is a joint model, as in FAST Bayesian FPCA's polar parameterization.
3. **Learned directions are frames too.** In projection pursuit (BPPR), each direction lives on a sphere, identified only up to sign under hinges, since (w, t, s) and (−w, −t, −s) give the same function. BASS is axis-aligned, so rotate-and-refit's choice of frame inside the active subspace changes the refit's sparsity: a Stiefel or flag choice, not a Grassmann one. Candidate samplers and priors, cited from memory and not yet read: matrix Bingham–von Mises–Fisher; Jauch, Hoff and Dunson's polar expansion and Cayley transform; Byrne and Girolami's geodesic Monte Carlo.
4. **A two-sided picture for functional outputs (the main session's reading, to test).** With J(x) the Jacobian, outputs by inputs:
   - E[Jᵀ W_y J] is the functional active subspace, on the input side.
   - E[J W_x Jᵀ] gives the output patterns that move with the inputs. PCA's Cov(f) gives the patterns that vary most.
   - A Poincaré bound links them one way only: low gradient energy implies low variance, not the converse. So PCA truncation can drop a small-amplitude mode that varies fast with the inputs.
   - The fit cannot see a truncated mode. D07's leave-one-case-out truncation residuals can. If they are predictable from x, the truncation discarded input-driven structure, and the basis predictive, which adds them as independent of x, is locally miscalibrated.
   - The field's term to check is Cook's response envelopes (material against immaterial variation, estimated on a Grassmannian), with a Bayesian version (Khare, Pal and Su 2017, from memory). E34 is multivariate linear regression on BASS's basis, so an envelope on E34 would be a predictor-aware output reduction, the alternative to PCA preprocessing. Also nearby: reduced-rank regression, and Zahm et al.'s vector-valued active subspaces (arXiv 1801.07922, unread).
5. **Cost.** E8's first cut needs only SVDs (projector averages and principal angles). Manifold sampling enters only with learned directions or a Bayesian frame.

## Code facts found on 2026-11-08

- **`StiefelManifold.Distance` is still wrong,** as the Astra notes flagged. It sums the principal angles' squares and ‖(I − YYᵀ)Z‖², so two frames of one subspace (v and −v, for instance) read 0, and the movement out of the subspace is counted twice. It needs a fix before anything that depends on the frame uses it.
- **`GrassmannManifold.Distance` runs on a thin SVD** (`e585b18`), so the Gram-squaring trap the Argon notes raised is gone. It reads angles by arccos, which loses near-identical subspaces; a sine-based route (Björck–Golub; Knyazev–Argentati) keeps them, which matters for well-identified posteriors and for SPC's adjacent temperatures.

## Questions the user may want to take up

- Which posterior summaries E8 should report, and how they read under weights and under a dependent input measure (the gradient readout averaged over the measure's points).
- Whether the truncation-residual check belongs with the basis fit now, and what decision its outcome would change.
- Whether an envelope on E34, or a regression-aware basis, should become an extension entry, and how it relates to the Astra notes' "basis chosen for predictive usefulness".
- Whether the Stiefel distance fix and the sine-based Grassmann angles should be done now as a small geometry fix.

## Working rules

- **A discussion first.** Code, entries in `extensions-brief.md` and checklist changes wait for the user's go-ahead. Notes and readings go in this file's Report.
- **Sources:** read a paper before stating what it says; mark anything cited from memory as such. Sources live under `D:\aghado01\codex-scientiae` (`supellex-private/PDFs`, `gauntlet/`, with `inventory.jsonl`).
- **Field terms:** when hunting a name for an idea, survey the established terms across fields before coining one.
- **Shared checkout:** other sessions are working in it (KHAOS slice 3 at the time of writing). Touch only this file unless the user asks for more. Commit by pathspec (`git commit -- <paths>`), leave nothing staged, push nothing, and leave `issues/regression/bass/basin-hold-validation.csv` and `path-profile-validation.csv` out of every commit.

## Report

### 2026-11-08, the IRLS chip's outcome, and what it adds to the angle-primitive question

The [IRLS chip](../../geometry/irls-scale-brief.md#report) landed `61a7fac` and `e70afaf`.
- **What the solver now does.** It reads its constants in units of the data's spread, a weighted lower median of the radii, recomputed each iteration. The old absolute semantics stay reachable as `IrlsScale.Data`, pinned bit for bit. A median that is a data point is now returned exactly.
- **The arccos item now weighs more.** In spread units the median still stalled below s = 1e-4, because `Distance`'s arccos error (about ε/2θ²) entered the radii, so the solver now reads its radii from ‖Log‖. The sphere has the same floor as the Grassmannian. The Exp maps' absolute 1e-12 guard bounds the median's relative accuracy at about 1e-12/s. The cut-locus fallback depends on `Distance` being right near π/2.
- **This strengthens the case for the stable principal-angle routine** (synthesis item 2, and the chip's finding 4). Its placement in LinAlg still waits for the user.
- **Outside this discussion, and also waiting for the user:**
  - `MoMPCA`'s α̂ clamp breaks its equivariance beyond about 1e3 either way, so the 1e-6 case was not tested.
  - Weiszfeld's slow rate near a nearly optimal data point.
  - `EuclideanMedian` stalls when warm-started on a data point.
  - GMM's robust initialization has no test.

### 2026-11-08, where the subspace code should live (proposed)

The user is weighing a specialized implementation under Regression, since Geometry was not meant to be imported by Bayesian code.
- **Regression does not import Geometry today.** `Maths.Regression` references Samplers, Distributions, LinAlg and Rng.
- **E8's default path is linear algebra only.** Each draw's subspace is a thin SVD of its gradient matrix, E[C] and P̄ (mixed ranks included) need a symmetric eigendecomposition, principal angles are two SVDs and an atan2, and the rank and containment split is a residual's norm. That is about a hundred lines on `ThinSvd` and `DenseEigen`, with no manifold interface and no IRLS. Rotate-and-refit's frame, a FAST-style polar frame and BPPR's power spherical proposal are also linear algebra or direct sampling.
- **Geometry is needed only for the intrinsic options:** the Karcher mean, the median and their scatter, which the synthesis ranks as options, and whose solver needs the scale fix below.
- **Proposed split** (the ownership the Astra notes recorded):
  - The stable principal-angle primitive goes in LinAlg beside `ThinSvd`. `GrassmannManifold.Distance` and `LogMap` call it, and so does E8, so there is one implementation and Regression never imports Geometry.
  - The posterior summaries go in Regression: P̄, the posterior over k, the readouts, and the Bayes-action reading.
  - The manifold estimators stay in Geometry. Should E8 ever offer them, an orchestration layer above both modules passes the draws across.
  - A manifold proposal for a future sampler, such as an envelope's Γ step or a Givens or Cayley move, belongs with that sampler, not with IRLS.

### 2026-11-08, probe: the IRLS median and mean on tight clusters

The probe ran at `d571dff`, built from a `git archive` snapshot. Its code is outside the repo, copied to the shared notes folder's `grassmann-probe/`.
- **Setup.** Gr(2, 10) and Gr(3, 20), three seeds each, with 21 points exp_{W0}(s·u_i) for Gaussian horizontal u_i and spreads s from 1e-2 to 1e-9.
- **Truth from outside the solver.** The median is exp_{W0}(s·m_u), where m_u is the Euclidean median of the u_i, found by plain Weiszfeld at scale 1. The mean is exp_{W0}(s·ū). Both are exact to O(s³), so the 1e-2 column's 3e-6 for Grassmann is that approximation.
- **Comparisons.** A Euclidean control with points c + s·u_i, whose truth c + s·m_u is exact. A probe-only sine-based Grassmann beside the struct, with atan2 angles and Edelman–Arias–Smith's atan-form log.
- **Warm start.** At each manifold's medoid, as `MoMPCA` does.

The worst error over the six cases, as a fraction of s, with iteration counts where the solver hit `MaxIterations`:

| | 1e-2 | 1e-4 | 1e-6 | 1e-7 | 1e-8 | 1e-9 |
| --- | --- | --- | --- | --- | --- | --- |
| Median, defaults, Euclidean | 2e-7 | 5e-2 (200 it.) | 5e-2 (200) | 5e-2 (200) | 5e-2 (200) | 5e-2 (200) |
| Median, defaults, Grassmann (arccos) | 3e-6 | 5e-2 (200) | 5e-2 (200) | 5e-2 (200) | 0.92 | 1.1 |
| Median, defaults, Grassmann (sine) | 3e-6 | 5e-2 (200) | 5e-2 (200) | 5e-2 (200) | 5e-2 (200) | 5e-2 (200) |
| Median, Weiszfeld only, Euclidean and sine | 3e-6 | 6e-6 | 3e-4 | 3e-3 | 4e-2 | 4e-2 |
| Median, Weiszfeld only, arccos | 3e-6 | 2e-5 | 2e-3 | 2e-2 | 0.92 | 1.1 |
| Karcher mean, arccos | 3e-6 | 9e-10 | 7e-10 | 7e-9 | 0.91 | 1.1 |
| Karcher mean, sine | 3e-6 | 9e-10 | 6e-10 | 7e-9 | 7e-8 | 5e-7 |

- **The solver's absolute constants are the larger defect, and they act on every manifold.** `IrlsOptions.Default` switches the median to subgradient steps of size η₀/√k = 1/√k, absolute, whenever the iterate is within 1e-4 of a data point.
  - Below a spread of 1e-4, the Euclidean median included, it runs to `MaxIterations` and returns an estimate about 5% of the spread away from the median.
  - With Weiszfeld only, the absolute tolerance of 1e-8 stops the iteration early as the spread approaches it. So the solver is not scale-equivariant.
  - Callers: `MoMPCA`'s four medians (its Euclidean factor is in the data's units) and `DistributedSpred`. `EuclideanMedian`, GMM's, has its own loop with the same absolute 1e-8 tolerance; it was not probed.
  - The fix makes the threshold, the step and the tolerance relative to the data's scale, for instance to the warm start's median distance to the data.
- **The arccos floor is real but smaller.**
  - `Distance` loses accuracy as ε/θ²: a relative error of 1e-5 at 1e-5 rad, 4e-4 at 1e-6 and 5e-2 at 1e-7. At 1e-9 it reads up to ten times the true distance.
  - `LogMap` holds to 1e-7, because its norm comes from the residual scaled by θ/sin θ, which barely depends on θ's error. From 1e-8 it reads zero.
  - So the Karcher mean is accurate down to a spread of 1e-7 and below that returns its warm start. The sine reference holds at every spread, within 5e-7 of s at 1e-9.
- **The decision.**
  - Moving `Distance` and `LogMap` to the residual route is a fix below a spread of about 1e-7, and below 1e-5 for distances read directly. Above that it is only tidying.
  - The solver's scale constants are a fix for every consumer at spreads below 1e-4, whatever the manifold.
  - Neither touches P̄, E8's default.

### 2026-11-08, the manifold module and what it changes

The module: `IRiemannianManifold` (`Distance`, `LogMap`, `ExpMap`, `AddScaled`, `Norm`, static `IsFlat`); the manifolds; the wrappers `ScaledManifold<T>` (a constant metric scale, with Log and Exp unchanged) and `ProductManifold<TA, TB>`; `Irls.Solve<TManifold, TLoss>`, which dispatches on three questions (a closed-form L² Karcher step or not, flat or curved, Weiszfeld or a projected subgradient); the losses `L2Loss` and `L1Loss`, the only two built, although the comments name Huber and Tukey; `GeometricMean` and `GeometricMedian`, each with its scatter and consistency factors; and `SubspaceAnnealer`, simulated annealing on Gr(k, d) with a black-box objective. `MoMPCA` composes them as a median on Scaled(ℝᵖ, α) × Scaled(Gr, 2 − α).
- **E8's summaries are three losses on one interface.**
  - P̄ is the chordal one, in closed form.
  - `GeometricMean` on `GrassmannManifold(p, k)` gives the Karcher mean, for squared geodesic loss, and `GeometricMedian` gives the median, for geodesic loss. Both take the draw weights as their `weights`.
  - The intrinsic two need every draw inside the injectivity radius of the iterate, all principal angles below π/2, which a multimodal posterior breaks. That is one more reason for P̄ as the default and the other two as options.
- **The Karcher scatter is a readout that is already built.** At the mean it is the tangent covariance of the draws' subspaces, principal geodesic analysis (Fletcher et al. 2004, cited from memory). Its top eigenvectors are the coupled ways the subspace moves across draws, such as direction i tilting toward complement direction j, which P̄'s spectrum does not show. It is D × D on the n·k buffer, of rank at most k(p − k).
- **Mixed ranks enter through the flat embedding, not through a new manifold.** One instance has a fixed `Dimension`, and no geodesic joins subspaces of different ranks. Each P_s as a vector on `EuclideanVectorManifold(p²)` does admit mixed ranks:
  - the flat mean is P̄;
  - the flat median is a robust extrinsic summary;
  - the flat distance divided by √2 is the projection distance.

  That shrinks synthesis item 2 to the stable principal-angle routine, with `Distance` and `LogMap` moved onto it, plus the split readout (rank term and containment defect).
- **The solver's defaults sit at the arccos floor, which corrects the orientation's "costs nothing".** `IrlsOptions.Default` has Tolerance = 1e-8 and Epsilon = 1e-10, while `GrassmannManifold` resolves angles only to about 1.5e-8. A Grassmann mean or median therefore runs its convergence test and its coincidence check at or below the resolution of the angles it reads.
  - This does not touch P̄, E8's default, but it does touch `MoMPCA`, `DistributedSpred` and SPRED, which run the median today.
  - The decision it changes: whether the residual-route `Distance` and `LogMap` are a fix or only tidying. A seeded test decides it, a Grassmann median on tight clusters with spreads from 1e-6 down to 1e-9, checking termination and accuracy against a known centre.
- **Products and wrappers.** E8 could summarize (subspace, log eigenvalues) jointly on Gr(k, p) × ℝᵏ with MoMPCA's α calibration, and co-active pairs on Gr × Gr. The Astra notes' robust fixed PCA basis for D07 already composes: `MoMPCA`'s centre and basis enter through `BassDecomposition.FromBasis`, with the truncation residuals supplied.
- **`SubspaceAnnealer` covers span objectives only.** It fits a predictor-aware output basis posed as an objective on the span. Separate fits' sparsity depends on the frame, so that objective cannot be posed on it, and rotate-and-refit's frame choice cannot either.

### 2026-11-08, BPPR read: its directions are sparse points on the sphere, not lines

Read: arXiv 2210.09181v1 (`gauntlet/bayesian/2210.09181v1`, main text §2.3 to §2.4 and §2.7) and the manual of the R package BayesPPR 0.2.0 (`PDFs/Bayesian/BayesPPR.pdf`).
- **The ridge is one-sided, so a direction's sign carries information.** Each ridge is β′b(xᵀθ | t), a modified natural spline whose first basis function is (u − t₀)₊, so the ridge is zero below t₀ (eq. spline_expansion). The prior probability that t₀ falls inside the data's projections, which makes the ridge local, is p₀ = 2/3 (the package's `prob_relu`). θ and −θ therefore place the ridge at opposite ends of the projection, and the prior on θ is uniform on the whole sphere, not on a half-sphere.
  - This corrects synthesis point 3 and the orientation's "each direction is a line": the (w, t, s) ↔ (−w, −t, −s) symmetry needs a ridge family closed under reflection, such as a hinge with a sign parameter, which BPPR's is not.
  - Only global ridges, those with t₀ below every projection, are close to sign-symmetric.
- **The directions are not a frame, which confirms the orientation.** Each θ_m has an independent uniform prior on the unit sphere, with at most A = min(3, p) nonzero coordinates, a set 𝒥_m drawn as BMARS draws variables. No orthogonality is imposed. A draw is a variable-size collection of sparse points on coordinate subspheres, and M moves by reversible jump.
- **The sampler is plain spherical Metropolis–Hastings.** A change step proposes the active coordinates from a power spherical distribution (De Cao and Aziz, cited in the paper) centred on the current direction, with κ = 1000 in the paper and `scale_proj_dir_prop` ≈ 0.002 in the package. Neither matrix Bingham–von Mises–Fisher nor geodesic Monte Carlo is needed: synthesis point 3's candidate samplers belong to frames, not to BPPR. `SphericalManifold` covers its geometry.
- **The direction prior has the input-metric question that E8's decision 2 settles.** The paper standardizes the inputs to mean 0 and variance 1, so that the uniform prior gives every input an equal chance. That makes the prior depend on the input metric, just as C's eigenspace does. Under dependent inputs, standardizing and whitening differ.
- **The multivariate mode is D07's separate-fit setting, without the truncation residual.** `bppr_pca` fits each principal component independently, with `prop_var` = 0.99 by default, and the paper drops the truncation error r_i "under the assumption that D⁻ is large enough". D07's leave-one-case-out residuals already go further.
- **E8 on a BPPR fit.** A ridge is not separable, so C runs on E8's sampled route. C = Σ_{m,l} θ_m θ_lᵀ E[g_m′ g_l′] lies in span{θ_m}, of rank at most M. Comparing draws' directions directly would meet label switching across ridges as well as a varying M, so E8's subspace summaries are the natural readout of a BPPR posterior.

### 2026-11-08, synthesis: what is worth doing (proposed, not yet settled by the user)

E8's build brief (`e31d7dd`) has already taken the orientation's corrections: the declared input metric (its decision 2), both posterior summaries with P̄'s eigenvalues read as mean squared cosines (4), the SVD of the gradient matrix and shared points (5), and `StiefelManifold` set aside. E34's entry took the removal of the frame and the envelope's scope (`6e5e53e`). What remains, in order of worth:

1. **The active dimension as a posterior quantity (E8 slice 2).** Each draw gets its own k by a declared rule, either an eigengap ratio or a share of tr C, and the readout reports the distribution of k_s beside a fixed k.
   - The cheapest form needs no new geometry: apply the rule per draw and give the histogram.
   - With mixed ranks, P̄ = Σ w_s P_s stays defined and its trace is the mean rank. Its ordered eigenvectors are the flag mean for subspaces of mixed dimensions (Draper et al. 2014, cited from memory).
2. **Comparing subspaces across ranks, in Geometry, with `GrassmannManifold` moved onto the same routine.**
   - The routine: principal angles from the cosines and the residual's sines, combined by atan2; the projection distance split into its rank term and its containment defect; and the containment alone.
   - `Distance` and `LogMap` move to the residual route. Nothing needs to lose the arccos floor today, but once the routine exists the move costs nothing.
   - Its consumers in the plan: item 1's mixed ranks, item 3's basis adequacy, and co-active comparisons between fits whose active dimensions differ. It lands before slice 2.
3. **The output side, in slice 3.** G = E[J W_x Jᵀ], computed from the thin SVD of the stacked, metric-weighted Jacobians (m × Np) without forming the m × m matrix. E8's brief names "reading which output patterns move with the inputs" in its Why, but no slice delivers it. Two readouts follow:
   - *A matrix Poincaré check (derived here).* Under a product law, Efron–Stein and the one-variable Poincaré inequality give Cov_ρ(f) ≼ E[J D Jᵀ], where D is the diagonal of the inputs' Poincaré constants: (b − a)²/π² for a uniform input, so Cov_ρ(f) ≼ G/π² under the box metric. It is the matrix form of decision 6's check, and it makes "low gradient energy implies low variance" exact, per draw.
   - *Basis adequacy.* The containment defect of G's top-k eigenspace in the PCA basis's span measures how much input-driven output structure the truncation dropped. It needs G from a direct-mode fit, because a basis fit's own G lies inside span(K) by construction.
4. **A predictor-aware output basis needs no new decomposition.** G's top eigenvectors enter through `BassDecomposition.FromBasis`, with the truncation residuals supplied. Whether that basis beats PCA is an experiment for problems with functional outputs, gated by item 3's adequacy readout.
5. **The truncation-residual check (D07).** It asks whether the leave-one-case-out residuals are predictable from x, against a permutation of the cases as the null. The decision it changes is D07's rank rule:
   - If it fires on real cases, a rank rule by predictability joins variance share. That rule keeps components while their fits carry signal.
   - If it does not fire, variance share stands.
   - Gavish and Donoho's threshold (cited from memory) is a third rule, for noisy outputs only.

Later, and only under the stated condition:
- **Rotate-and-refit, after slice 4,** can take its frame from the ordered eigenvectors of P̄ or E[C]. BASS ignores an input's sign and the inputs' order, so that frame fully specifies the refit, and no Stiefel code is needed. Stiefel enters only if a rotation is optimized for a simplicity criterion, the role varimax plays in factor analysis. The refit treats the subspace as fixed, the same plug-in the PCA basis is.
- **A Bayesian PCA frame** (FAST), only if the stability ratio shows an unsettled basis on real cases. E34's joint fit is the cheaper way to remove the frame.
- **Projected BGP and BPPR,** later engines. BGP with an isotropic projected kernel needs only Grassmann geometry. BPPR's directions are sparse points on the sphere, with a meaningful sign (the BPPR reading above).

Not worth doing now:
- the Stiefel rewrite, which has no consumer;
- flag-manifold code, since the flag mean is P̄'s eigenvectors;
- Davis–Kahan bounds as a readout, since the draws give the angles directly;
- manifold samplers, which are needed only with a Bayesian frame or learned directions;
- a robust median as the default summary.

### 2026-11-08, what the AbliteraLUS comparison note offers this brief

The note (`issues/geometry/gemini-abliteralus-thermomapper.md`, a comparison Gemini wrote) proposes five ports. Its Grassmann claims were checked against `abliteralus/analysis/grassmann.py`.
- **Sine-based angles.** `principal_angles` takes the cosines from svdvals(Y Zᵀ) and the sines from the residual of Y off span(Z), with bases stored as rows, and pairs them by atan2. That code's `log_map` and `karcher_mean` still use acos, so the note's claim that the stable angles matter for Karcher convergence does not describe it. A port here has to cover `LogMap` as well, since `GeometricMedian` iterates on it. Edelman, Arias and Smith's form, U atan(Σ) Vᵀ from the SVD of (I − YYᵀ)Z(YᵀZ)⁻¹, reads the residual directly (cited from memory).
- **The projection distance across ranks applies to E8.** It is the projection mean's own geometry, so each draw can keep its own k, chosen by its own gap, and both P̄ and the distance stay defined. For k_y ≤ k_z it splits exactly: d² = (k_z − k_y)/2 + ‖(I − P_Z)Y‖²_F, a rank term plus a containment defect.
  - Reported apart, the two terms say whether the draws differ in how many directions they hold or in which directions.
  - The implemented form, (k_y + k_z)/2 − Σcos², has the same floor as arccos. The split form computes the defect from the residual and so avoids it.
  - Ye and Lim (2016) is the established reference for distances between subspaces of different dimensions (cited from memory).
- **The diameter bound.** √min(k, p − k)·π/2 is correct for the geodesic distance. The projection distance has a different maximum, so normalizing a distance needs that metric's own diameter.
- **Marchenko–Pastur and BBP bear on D07's rank, not on E8.** For noisy outputs with m ≫ n, the leading sample eigenvector's angle to the population one does not vanish (Paul 2007, cited from memory), which quantifies the PCA basis's Grassmann error. The note's noise-floor rule, σ² estimated from the bulk median, is Gavish and Donoho's optimal hard threshold (cited from memory). Neither applies to a deterministic simulator, whose outputs have no noise bulk, and both still rank components by noise rather than by how predictable they are.
- **Not for this brief:**
  - Bures–Wasserstein: a PSD barycenter of the draws' C is possible, but no decision rests on it.
  - Empirical curvature, cones and DSI, and the seam detector.
  - The orthogonal sham control is a negative control. E8's counterpart is the response's variation along random directions in the complement; the truncation-residual test's counterpart is a permutation of the cases.

### 2026-11-08, orientation pass: the code read, and corrections raised before discussion

Checked against the code:
- **`StiefelManifold` ignores the frame in `LogMap` as well as in `Distance`.** `LogMap` forms log(YᵀZ) as U diag(log σ) Vᵀ, so for Z = YQ with Q orthogonal every σ is 1 and log_Y(YQ) = 0. `ExpMap` re-orthonormalizes Y + Δ, a projection retraction, not the geodesic its comment names. Nothing in `src/` or `tests/` uses the struct, so nothing downstream is wrong today. The repair is a rewrite with its first tests, made when a frame consumer arrives.
- **`MoMPCA`'s variances no longer have the problem the Astra notes raised.** They are now the eigenvalues of the geometric median of the blocks' projected covariances (`ComputeMoM`, step 6).
- **`MoMPCA`'s down-weighting does not carry over to E8.** Its α̂ balances the subspace factor against the mean factor of the product metric, and a posterior over subspaces has no mean factor. What carries over is `GeometricMedian` on `GrassmannManifold`, which down-weights draws far from the bulk, not draws with weak gaps. For a posterior those outlying draws are real mass, a minority mode, so they should be reported rather than resisted.
- **The arccos floor.** With σ near 1, arccos resolves θ only down to √(2ε) ≈ 1.5e-8, with relative error about ε/θ². A posterior's spread loses nothing unless its draws agree to within about 1e-6 rad. The floor matters where small angles are divided or differenced, as between SPC's adjacent temperatures.

Corrections to the synthesis, raised in the orientation and not yet discussed:
- **P̄'s eigenvalues are mean squared cosines, not counts.** λ_j = Σ w_s ‖W_sᵀu_j‖², so a value of 0.5 arises both when every draw sits at 45° to u_j and when half the draws contain u_j exactly. Their sum, tr P̄, is k.
- **The projection mean is the Bayes action under chordal loss.** E‖P − P_s‖²_F = 2k − 2 tr(P P̄), which the top-k eigenspace of P̄ minimizes. Its eigenvectors in order are the flag mean (Draper et al. 2014, cited from memory), so the estimates for every k come nested from one eigendecomposition.
- **BPPR's directions do not form a frame.** Projection-pursuit ridges carry no orthogonality, and the hinge's sign symmetry makes each direction a line. M ridges therefore sit on (ℝP^{p−1})^M, up to relabelling. Stiefel geometry enters only where the rotation is orthogonal, as in rotate-and-refit.
- **The separate fits identify the PCA frame only up to signed permutation.** Flipping a component's sign or relabelling the components leaves their joint posterior unchanged: the coefficient priors are symmetric, and the structure prior never sees the response. The frame's space is therefore St(r, m) modulo the hyperoctahedral group, locally Stiefel. A distance between two bases must minimize over that group, and `BassDecomposition`'s largest-entry sign convention can flip between nearby bases.
- **E34 makes the frame irrelevant.** Its equivariance check (Y → YL with Ψ₀ → LᵀΨ₀L) means that fitting the PCA scores jointly, with shared inclusion, depends only on their span. The frame matters only through the components' separate structures, so Stiefel geometry is the price of per-component sparsity.
- **The input metric fixes C's eigenspace.** For x = Lz, C_z = Lᵀ C_x L, so the active subspace moves under any non-orthogonal reparametrization, and E8 has to declare its input metric. Under a product uniform law, the box scaling equals whitening by the law's covariance up to a constant. Under a dependent law the two differ, and whitening, the generalized problem C_x v = λ Σ_ρ⁻¹ v, is the affine-equivariant choice.
- **Per-draw C under Monte Carlo.** Each draw's eigenspace should come from the thin SVD of its N × p gradient matrix, not from forming C, which squares the gaps. The draws should share the measure's points, so that their spread is the posterior's and not the Monte Carlo error's. `BassGradient` already reads each draw's sampled coefficients, not its conditional means.
- **An envelope on E34 suits only small K.** It needs Σ over the outputs, which E34 keeps only for small, unordered K. On PCA scores it reduces only within the retained span, so it cannot recover a mode PCA dropped. For a functional output, the predictor-aware alternative is an output basis from E[J W_x Jᵀ], computed from a direct-mode pilot fit whose gradient readout already exists (a reading here). Given Γ, the envelope seems to keep E34's collapse: Γᵀ Y is a u-output E34 and Γ₀ᵀ Y is noise alone (the structure of Khare, Pal and Su, cited from memory).

Terms and sources to check, all cited from memory: Baptista, Marzouk and Zahm's gradient-based data and parameter dimension reduction, a likely established home for the output-side matrix; varimax and quartimax (Kaiser 1958), for choosing a frame within a span by simplicity; MatchAlign (Poworoznek, Ferrari and Dunson), for aligning frames across draws; the flag mean; Knyazev and Argentati (2002); Cook, Li and Chiaromonte (2010); Khare, Pal and Su (2017). Held locally, in `supellex-private/gauntlet/bayesian`: 2311.18146 (co-active subspaces) and 2412.11340 (FAST Bayesian FPCA). Not held: 2307.11241, 1801.07922, 2210.09181 and 1907.11572.
