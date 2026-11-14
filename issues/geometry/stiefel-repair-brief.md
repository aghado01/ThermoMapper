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
   - *Why this one:* the type and `docs/Maths/RiemannianGeometry.md` already promised it. Its Exp needs only exponentials of skew matrices, and its Log has a standard algorithm (Zimmermann 2017, cited from memory, not on disk).
   - *The alternative, the embedded (Euclidean) metric tr(ΔᵀΔ):* Absil, Mahony and Sepulchre 2008, Example 5.4.2 (read) give its geodesic. That geodesic needs a general matrix exponential of a non-skew 2r × 2r matrix, which LinAlg lacks, and its Log needs a shooting method. A one-parameter family of metrics joins the two (Hüper, Markina and Silva Leite 2021, cited from memory). Neither is built.
   - *Implication:* `Norm` changed from Frobenius to the canonical norm. No caller existed.
2. **Exp in closed form:** Y(t) = [Y Q] exp(t·[A, −Bᵀ; B, 0]) [I; 0]. Here A = skew(YᵀΔ), and QB is the thin SVD of (I − YYᵀ)Δ over the directions above 10⁻¹³‖Δ‖_F, of which there are at most n − r. Δ is read through its tangent part.
3. **Log by Zimmermann's iteration, with a completion the size of the normal part.**
   - [YᵀZ; QᵀZ] is completed to an (r + k) × (r + k) rotation, where k is the rank of Z's part off Y's span (sines above 10⁻¹⁴, at most n − r). So n < 2r works, and a pure turn within the span (k = 0) is log(YᵀZ) directly.
   - The completion is Procrustes-aligned so that its lower block is symmetric positive semidefinite, which starts V near the identity for nearby frames. Its determinant is then made +1 by reflecting along the lower block's least singular direction.
   - The stopping rule is ‖C‖_F ≤ 10⁻¹⁵(1 + ‖log V‖_F), or ‖C‖_F no longer falling once it is below 10⁻¹¹(1 + ‖log V‖_F).
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
- *Zimmermann 2017 is not on disk.* Its algorithm is checked here by properties, not against the paper. If the user adds it, the convergence radius it proves can replace the probe's measured range in the docstring.
