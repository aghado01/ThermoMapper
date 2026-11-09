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
