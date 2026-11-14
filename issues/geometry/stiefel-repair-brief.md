# StiefelManifold's repair (brief and report)

Date: 2026-11-13. Status: done (`6a55e2f`). Commissioned by the user in the Grassmann discussion's session: "lets do StiefelManifold's repair, which need not wait for a frame that needs it. i think this is in the e8 brief". E8's brief had deferred it until a consumer optimized over frames. Reports go under "Report" below, newest first.

## Why

`StiefelManifold` (`src/maths/geometry/manifolds/StiefelManifold.cs`) promised orientation-sensitive canonical geometry (Edelman, Arias and Smith 1998), but its code did not deliver it.
- **Distance** read only the column spans, so d(Y, YQ) = 0 for every orthogonal Q.
- **Log** was an approximation.
- **Exp** was the Gram–Schmidt retraction MGS(Y + Δ).
- **Norm** was the Frobenius norm, which is the embedded metric, while the docstring named the canonical one.

The type had no consumers, and no tests.

## Decisions

1. **The canonical metric**, ⟨Δ, Δ⟩ = tr(Δᵀ(I − ½YYᵀ)Δ). It is the quotient of O(n)'s bi-invariant metric by O(n − r): frames are read as rotations, and turning a frame within its span weighs half what moving the span does.
   - *Why this one:* the type and `docs/Maths/RiemannianGeometry.md` already promised it. Its Exp needs only exponentials of skew matrices, and its Log has a standard algorithm (Zimmermann 2017, arXiv 1604.05054v5, read after the build; see decision 3).
   - *The alternative, the embedded (Euclidean) metric tr(ΔᵀΔ):* Absil, Mahony and Sepulchre 2008, Example 5.4.2 (read) give its geodesic. That geodesic needs a general matrix exponential of a non-skew 2r × 2r matrix, which LinAlg lacks, and its Log needs a shooting method. A one-parameter family of metrics joins the two (Hüper, Markina and Silva Leite 2021, cited from memory). Neither is built.
   - *Implication:* `Norm` changed from Frobenius to the canonical norm. No caller existed.
2. **Exp in closed form:** Y(t) = [Y Q] exp(t·[A, −Bᵀ; B, 0]) [I; 0]. Here A = skew(YᵀΔ), and QB is the thin SVD of (I − YYᵀ)Δ over the directions above 10⁻¹³‖Δ‖_F, of which there are at most n − r. Δ is read through its tangent part.
3. **Log by Zimmermann's iteration, with a completion the size of the normal part.**
   - [YᵀZ; QᵀZ] is completed to an (r + k) × (r + k) rotation, where k is the rank of Z's part off Y's span (sines above 10⁻¹⁴, at most n − r). So n < 2r works, a pure turn within the span (k = 0) is log(YᵀZ) directly, and k = 1 needs no iteration because C is then 1 × 1 and skew.
   - The completion is Procrustes-aligned so that its lower block is symmetric positive semidefinite, which starts V near the identity for nearby frames. Its determinant is then made +1 by reflecting along the lower block's least singular direction.
   - The stopping rule is ‖C‖_F ≤ 10⁻¹⁵(1 + ‖log V‖_F), or ‖C‖_F no longer falling once it is below 10⁻¹¹(1 + ‖log V‖_F).
   - **Against the paper,** read after the build when the user added it (Zimmermann 2017, arXiv 1604.05054v5):
     - *Agrees:* §3's algorithm step for step. That is M = YᵀZ, the thin QR of Z − YM, the orthogonal completion, the same Procrustes step (Y₀ ← Y₀R₂D₂ᵀ, his remedy for a Gram–Schmidt completion), and the updates log, exp(−C) and V ← V·diag(I, Φ). §2's canonical metric and Exp, after EAS §2.4, agree too.
     - *His guarantee:* the convergence theorem of §4 proves linear convergence, ‖C‖₂ at least halving each step, for ‖Y − Z‖₂ < 0.0912. His §5.3 experiments converged wherever ‖log V₀‖₂ < π, and never failed below 0.9π for p > 10.
     - *His failure, removed here:* in St(4,2) his algorithm starts failing near distance π/2, because the completion can give V₀ with det −1, so an eigenvalue at −1 and no principal logarithm. His appendix "A critical special case" gives an explicit pair and notes that a sign flip of the completion would work. He leaves choosing a good V₀ open. The det +1 completion here does that, and the reduced completion makes his pair a k = 1 case. The probe's St(4,2) row converged at every length to 2.7, where his starts failing near π/2.
     - *Extensions not in the paper:* the rank-k completion, which his 2p completion lacks; his thin QR of a rank-deficient normal part also leaves columns not guaranteed orthogonal to Y. Then the det +1 completion, the Frobenius stopping rule with its floor, and the eigensystem-based log and exp. He computes those with MATLAB's `logm` and `expm`, and suggests Schur forms.
     - *His other note:* the convergence domain is related to the injectivity radius without being the same. He cites Rentmeesters (2013) for St(4,2)'s injectivity radius, estimated at about 0.71π ≈ 2.23. The probe's 40 random directions per length at 2.4 and 2.7 still returned the geodesic Exp had followed, which fits a cut distance that depends on the direction.
4. **Refusal, not a short log.** Log and Distance throw `InvalidOperationException` in three cases: the iteration does not converge, a rotation block comes within 10⁻⁸ rad of a half-turn, or Z = YQ with det Q = −1.
   - The Grassmannian and the sphere instead return a short log at the cut locus, and IRLS compensates by checking it against `Distance`.
   - Here `Distance` is the log's own norm. A short log would understate it silently.
5. **Two LinAlg primitives in `MatrixFunctions`.**
   - `SkewExponential`: exp(Ω) = cos √X + (sin √X/√X)·Ω, with X = ΩᵀΩ.
   - `RotationLogarithm` / `TryRotationLogarithm`: L = Σ (θ_q/sin θ_q)·v_q v_qᵀ·K over the eigenvectors of S = (V + Vᵀ)/2, where K = (V − Vᵀ)/2 and each angle is θ_q = atan2(‖K v_q‖, λ_q).
   - Both read symmetric eigensystems (Jacobi, element-wise tolerance 10⁻¹⁵). A first form took θ = arccos λ. Near a half-turn that amplifies the cosine's rounding by 1/sin θ, and it left 3·10⁻¹⁰ at π − θ = 10⁻³. The atan2 form with the tighter tolerance leaves 4·10⁻¹².

## Report

### 2026-11-13: done (`6a55e2f`)

**Landed.**
- `StiefelManifold`: the canonical Distance, Log, Exp and Norm.
- `MatrixFunctions.SkewExponential`, `RotationLogarithm` and `TryRotationLogarithm`. The existing members are unchanged, byte for byte.
- Tests: `tests/geometry/StiefelManifoldTests.cs` (15 cases) and `tests/linalg/RotationFunctionsTests.cs` (8 cases).

**Gates.** Each reference originates outside the manifold's own code.

| Check | Result |
|---|---|
| Rotations in known planes (closed form), 300 seeded trials, sizes 2–8, angles to π − 10⁻³ | exp to 2·10⁻¹⁴; log worst 4·10⁻¹² |
| r = 1 against the sphere, θ from 10⁻¹⁰ to 3 | distance to 10⁻¹⁵·max(1, 1/(π − θ)); log to 10⁻¹⁴ |
| r = n against O(n), two planes at 2.2 and −0.7 | distance and log to 10⁻¹³ |
| A turn by φ within the span | distance φ to 10⁻¹⁵, Grassmann distance 0 |
| Grassmann's distance | a lower bound on 20 random pairs, and equal on horizontal geodesics to 10⁻¹² |
| The first-variation identity ∇½d(·, Z)² = −log_Y Z in the canonical metric | to 10⁻⁸ by central differences |
| Zimmermann's critical pair in St(4,2) (his appendix), at t = 1, π/2 and 2 | Exp reproduces his Ũ(π/2) to 10⁻¹⁵; the log returns tΔ and the distance t to 10⁻¹⁴ (added in the follow-up commit) |
| Seeded round trip: 8 shapes (n ≤ 8, r ≤ 4) × lengths 10⁻¹⁰ to 2.5 | log(exp Δ) − Δ, length, and symmetry to 10⁻¹²; orthonormality to 10⁻¹⁴ |
| Invariance under left O(n) and right O(r) | to 10⁻¹³ |
| IRLS on a symmetric cloud of 8 frames | the mean within 10⁻⁹ and the median within 10⁻⁸ of the centre |
| Whole assemblies (`-Scope Slice`, from a `git archive HEAD` snapshot) | LinAlg 186/186, Geometry 106/106 |

**The Log's domain** (probe, `artifacts/stiefel-probe/`; 40 pairs per shape and canonical length, shapes (3,1), (4,2), (5,2), (3,3), (6,3), (8,4)):
- Every pair converged up to length 2.7, with round-trip error at most 2.2·10⁻¹⁴.
- At 2.9, 4 of 40 refused for (4,2) and 1 of 40 for (5,2); the rest converged, to 3·10⁻¹³.
- At 3.05 some pairs refused, and some returned a shorter geodesic than the one Exp had followed, as expected past the cut locus.
- Wall clock, indicative on a shared machine: under 1 ms per log up to length 2, and up to about 6 ms near 2.9.

**Not run.** The `Maths.Regression.Tests` selection, which reaches 225 classes by name through `MatrixFunctions`. Its existing members are byte-identical, so no existing path changes.

**Doc changes owed.**
- `docs/Maths/RiemannianGeometry.md`:
  - The "Stiefel (partial)" paragraph and the limitations bullet "`StiefelManifold` is partial" give way to the canonical metric, Exp, Zimmermann's Log, the domain and the refusals.
  - The "Grassmann SVD" paragraph's "Stiefel keeps the 10⁻¹² rule throughout" no longer holds.
  - The manifold table's Stiefel entries need the same update.
- `docs/Maths/LinearAlgebra.md`: `MatrixFunctions.SkewExponential` and `RotationLogarithm`, with the method, the atan2 angles, the accuracy near a half-turn and the refusal.

**Follow-ups, each with the decision it would change.**
- *The embedded metric as an option.* It changes which metric frames are averaged under. It needs a general matrix exponential (Padé scaling and squaring) in LinAlg, and a shooting Log. It waits for a consumer that wants the embedded metric.
- *Zimmermann 2017, resolved.* The user added arXiv 1604.05054v5 on 2026-11-13: R. Zimmermann, "A matrix-algebraic algorithm for the Riemannian logarithm on the Stiefel manifold under the canonical metric", *SIAM J. Matrix Anal. Appl.* 38(2), 2017. It is read against the build under decision 3, and its critical pair is now a test.
  - Its proven radius, ‖Y − Z‖₂ < 0.0912, lies far inside the measured range, so the docstring states both.
  - arXiv 1603.06757, added first and checked, is a different Zimmermann. Hernando, Igual and Quintana-Ortí (2017) speed up the Brouwer–Zimmermann algorithm for the minimum distance of a random linear code over GF(2), so it bears on nothing here.
