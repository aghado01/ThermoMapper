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
