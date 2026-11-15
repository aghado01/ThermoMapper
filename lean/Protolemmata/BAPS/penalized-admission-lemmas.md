# A penalized system's admission: identification free of the coefficients, the bound, the chain to the gate's κ₁

BAPS factors A(c) = Σₖ cₖSₖ at every sampler step, and a gate (`BapsAdmission`) admits a factor that is positive definite with an estimated κ₁ of its unit-diagonal scaling E(c) at most 10¹⁴. Slice 6 of the [defects brief](../../../issues/regression/shared-code-defects-brief.md#slice-6-the-samplers-condition-check-as-a-setting) builds a screen on a bound relative to one reference point, and the [BAPS follow-up's admission input](../../../issues/regression/bass/extensions-brief.md#input-for-the-follow-up-admission-and-singular-systems-2026-11-14) states the two failures behind the gate. These are the statements both rest on. They were commissioned on 2026-11-14, as targets 5 to 7 of the [chip brief](../BASS/active-subspace-lean-chip-brief.md#added-targets-a-penalized-systems-admission-2026-11-14). The coordinating session derived them that day and checked them against no run. The numeric checks are [`penalized-admission-checks.py`](penalized-admission-checks.py).

**Status (2026-11-14).** All three are proved, with no apologies, in `lean/Enthymemata/PenalizedAdmission.lean`, namespace `PenalizedAdmission`. Each main result prints only `propext`, `Classical.choice` and `Quot.sound`. Target 7 is the restated chain (`90fa6ce`). No statement proved false. Two follow-ups landed the same day (`7f535b3`): the Rayleigh extremes are the extreme eigenvalues (PA-6), and the column step's constant is (1 + √n)/2, sharp, which makes the composed factor q·(1 + √n)/2 (PA-7(c2′), (d′)). The user then asked for the conjecture they left open to be tried, and it proved false (finding 6, `fc793df`). The checks pass 27 of 27. Each target's own status line follows.

**What posing changed** (findings before any Lean):
1. **PA-7: the first statement's chain held, and the same pieces gave a constant √q smaller.**
   - *‖E(c)‖₁ ≤ q directly,* from the unit diagonal and entries of at most 1 in size. The first chain went through ‖E‖₁ ≤ √q·‖E‖₂ ≤ √q·λ_max(E), and λ_max(E) can itself reach q.
   - *No scaling theorem is needed.* The diagonal of A(c) is A(c′)'s times factors in [m, M], so λ_min(E(c)) ≥ (m/M)·λ_min(E(c′)) directly, where the first chain went through van der Sluis's κ₂(B) ≤ q·κ₂(FBF).
   - The defects chip reached the same chain on its own (its slice 6 entry, `76851ee`). The brief restated target 7 to it on 2026-11-14 (`90fa6ce`), with the Jacobi-scaling bound off the path. PA-7 below is the restated target.
   - The first statement, κ₁(E(c)) ≤ q·√(qn)·(M/m)·κ₁(E(c′)), follows from the restated one, since ‖E(c′)‖₁ ≥ 1 and q ≥ 1.
2. **What the banded families use of the proved factor (a reading).** The checks ran 240 families built from a cubic B-spline Gram and second-difference penalties, with n from 6 to 29 and q up to 7. The largest κ₁(E(c)) / ((M/m)·κ₁(E(c′))) was 1.007, against a proved q·√n of up to 38. (d)'s own bound was met to at most 0.145 of itself. Whether a smaller factor holds for every family is open. The sweep measures the slack.
3. **Van der Sluis's q in practice (a reading), off the path now.** On three banded unit-diagonal matrices, a local search over diagonal scalings found κ₂(B) / min_F κ₂(FBF) ≈ 1.000, 1.013 and 1.105, against q = 3, 5 and 7.
4. **Follow-up: the column step's √n is (1 + √n)/2, and that constant is sharp.**
   - For P symmetric positive semidefinite with λ_max(P) ≤ μ, a column u = Peⱼ has ‖u‖² ≤ μ·uⱼ. So the rest of its absolute sum is at most √(n − 1)·√(uⱼ(μ − uⱼ)), and ‖u‖₁ ≤ μ(1 + √n)/2.
   - A rank-one μvvᵀ attains it. At n = k² that is v = (k + 1, 1, …, 1), and at n = 4 it is a Lean theorem.
   - For n ≥ 4, E⁻¹ also attains it in the limit, for a unit-diagonal E. Take unit vectors xᵢ with Σᵢ vᵢxᵢ = 0: their Gram matrix has a unit diagonal and v in its kernel.
   - So the step admits no better constant, even on the matrices the gate sees. Against √n it gains a factor approaching 2.
5. **Follow-up: whether a smaller factor holds for every family (searches, readings).**
   - *The factor q is reached:* c₁J + c₂I gives κ₁(E(c)) / ((M/m)·‖E(c′)⁻¹‖₁) → n = q.
   - *The inverse part's reading stays near 1.* T = ‖E(c)⁻¹‖₁ / ((M/m)·‖E(c′)⁻¹‖₁) read 1.00, 1.18, 1.20, 1.23, 1.00 and 1.05 at n = 3, 4, 6, 8, 12 and 16, while its proved bound (1 + √n)/2 grows from 1.37 to 2.50. R/q read at most 1.20.
   - *The search covers every family.* Every family gives m·A(c′) ≼ A(c) ≼ M·A(c′), and every pair satisfying it is a two-term family: S₁ = (A − mA′)/(M − m), S₂ = (MA′ − A)/(M − m). So the search runs over pairs.
   - *Each reading is a lower bound* from a local search; the higher dimensions search worst.
   - *The conjecture they suggested,* that the inverse part is bounded by a constant, is false (finding 6). The searches missed a family of measure near zero.
6. **Follow-up: the inverse part is unbounded, so no factor of smaller order holds for every family.**
   - *The family.* S₁ = I − vvᵀ and S₂ = vvᵀ for a unit v, with c′ = (1, 1) and c = (M, m). Then A(c′) = I, so ‖E(c′)⁻¹‖₁ = 1, and A(c) = M·I − (M − m)·vvᵀ.
   - *Closed form,* for v ≥ 0, with ε = m/M and Δᵢ = 1 − (1 − ε)vᵢ²: column j of E(c)⁻¹ gives T = ε·Δⱼ + (1 − ε)·vⱼ√Δⱼ·Σᵢ vᵢ√Δᵢ. At v = (a, t, …, t) with a = 1/√3, as ε → 0, T → 2/9 + (2/(3√3))·√(n − 5/3), about 0.385·√n. The checks match the closed form to 10⁻⁸ from n = 16 to 1024.
   - *So the inverse part's supremum lies between about 0.385·√n and (1 + √n)/2.* The proved order is the true order.
   - *Lean.* `inverse_part_unbounded`: for every C there is a family meeting every hypothesis of (d′), with m = 1 and M = 2, whose inverse part exceeds C.
     - It goes through a dual bound, |yⱼ| ≤ ‖E⁻¹‖₁·maxₖ|(Ey)ₖ| (`abs_le_norm1_inv_mul`), which lower-bounds ‖E⁻¹‖₁ without computing the inverse.
     - The instance is rational: v = (3/5, 4/(5r), …, 4/(5r)) on r² + 1 coordinates and c = (2, 1). It has ‖E(c)⁻¹‖₁ ≥ 6r/25 (`norm1_inv_spike`).
   - *The whole chain.* A fixed dense q-block beside the family gives R/q growing like 0.035·√n: 0.25, 0.42 and 0.80 at n = 32, 128 and 512, with q = n/2. So no C·q bound holds for every family either, though R stays far below q·(1 + √n)/2.
   - *Why the searches missed it.* The family sits where A(c′) is diagonal and the change is rank one along a spiked v.
   - *Banded patterns too.* A local search over banded two-term families read T between 1.01 and 1.54 for q = 3 and 5 up to n = 64. A construction shows the same growth at q = 3:
     - S₁ = L, a path-graph Laplacian with heavy edges at one node, and S₂ = 2·diag(L) − L, its signless Laplacian. Both are tridiagonal and positive semidefinite.
     - Then A(c′) = 2·diag(L) and E(c′) = I, while E(c) = (1 − δ)·Ľ + δ·I. Here Ľ is the normalized Laplacian, whose null vector is spiked over the node and its two neighbours.
     - With M/m = 10¹², T reads 0.99, 1.70, 3.12 and 5.95 at n = 16, 64, 256 and 1024, about 0.19·√n.
     - At n = 256 it reads 0.90, 2.07 and 3.09 at M/m = 10², 10⁴ and 10⁶. The growth needs M/m large against the path's spectral gap, about n².
     - So √n cannot leave S on a pattern argument alone, even at q = 3. In practice M/m against the pattern's spectral gap limits it.
     - *Lean, the target* (the user asked for it). `inverse_part_unbounded_tridiagonal`: for every C, some family whose every Sₖ is zero beyond the first off-diagonals meets every hypothesis of (d′), with m = 1, and has T > C.
       - The instance is a path on N + 2 nodes. Its edges have unit weight except the last, (N, N + 1), which has weight N. S₁ = L and S₂ = |L|, the signless Laplacian, with c′ = (1, 1) and c = (ρ, 1), where ρ = 2(N + 1)² + 1. So A(c′) = 2·diag(L) and E(c′) = I.
       - The certificate is w_j = 1 − s·j(j + 1)/2 for j ≤ N and w_{N+1} = ((ρ − 1)/(ρ + 1))·w_N, where s makes (Aw)_N = 0. Then (Aw)_{N+1} = 0, |Aw| ≤ 10 and w_N ≥ 1/2.
       - The dual bound at node N gives ‖E(c)⁻¹‖₁ ≥ (ρ + 1)·√(N + 1)/20 (`norm1_inv_pathFamily`). So T > C once √(N + 1) > 20C, and the instance takes N = r².
       - The checks verify the certificate in exact arithmetic at N = 1, 4, 16 and 64.


**Symbols.** n the order. Sₖ are n × n real symmetric positive semidefinite, for k in a nonempty finite index set. c and c′ are coefficient vectors with every entry positive, c′ the reference point (c_ref in the defects brief), and A(c) = Σₖ cₖSₖ. ρₖ = cₖ / c′ₖ, m = minₖ ρₖ and M = maxₖ ρₖ. λ_min and λ_max are the extreme eigenvalues of a symmetric matrix, and κ₂ = λ_max / λ_min where it is positive definite. ‖A‖₁ is the largest column absolute sum, and κ₁(A) = ‖A‖₁‖A⁻¹‖₁. E(c) = D^{-1/2}A(c)D^{-1/2} with D = diag A(c), the unit-diagonal scaling. q bounds the nonzeros in any row (equally, any column) of every A(c).

**The smallest instance**, used throughout: S₁ = diag(1, 0) and S₂ = diag(0, 1). Each is singular, they share no null vector but 0, and A(c) = diag(c₁, c₂).

---

## PA-5 — identification does not read the coefficients

**Status:** proved (`form_comb_eq_zero_iff`, `comb_mulVec_eq_zero_iff`, `comb_posDef_iff`, `comb_posDef_iff_of_pos`, `penalized_det_eq_zero_iff`). Checks PA-5(i) to (iii) pass.

**Claims.** For c > 0 entrywise:
- **(a)** xᵀA(c)x = 0 exactly when Sₖx = 0 for every k, so ker A(c) = ⋂ₖ ker Sₖ at every such c.
- **(b)** A(c) is positive definite at one c > 0 exactly when it is at every c > 0, and exactly when the Sₖ share no null vector but 0.
- **(c) BAPS's instance.** ZᵀZ + Σⱼ λⱼPⱼ, with every λⱼ > 0 and every Pⱼ positive semidefinite, is singular exactly when some x ≠ 0 has Zx = 0 and Pⱼx = 0 for every j.

**Proof.**
- xᵀA(c)x = Σₖ cₖ·xᵀSₖx is a sum of nonnegative terms with positive weights, so it is 0 exactly when every xᵀSₖx is.
- For S positive semidefinite, xᵀSx = 0 forces Sx = 0, by Cauchy–Schwarz for the form: (yᵀSx)² ≤ (yᵀSy)(xᵀSx), at y = Sx.
- A(c) is itself positive semidefinite, so A(c)x = 0 exactly when xᵀA(c)x = 0.
- (b) follows from (a), since the kernel's description does not name c.
- (c) is (a) with S₀ = ZᵀZ (c₀ = 1) and xᵀZᵀZx = ‖Zx‖².

**Hypotheses, and what each does.**
- *Every cₖ > 0.* With cₖ = 0, Sₖ drops out and the kernel can grow. The hierarchical model's ξ² can reach 0 (the defects chip's slice 6 entry). There the claim is about the remaining terms.
- *Every Sₖ positive semidefinite.* With an indefinite term the form can vanish off the common kernel. For example, diag(1, −1) + diag(0, 1) is singular at c = (1, 1) and not at c = (1, 2).

**Checks.** Exact rationals. Nine families: the smallest instance, plus random integer Sₖ = GₖᵀGₖ of rank below n, four with a common null vector and four without. Each at five random rational c, with ker A(c) computed by row reduction. For (c), hat-function designs with n = 7 and three penalty forms. Data at one point gives a singular system. Two points identify a single second-difference penalty but not a block-split pair, and the checks find each of these at every c tried.

**Decision.** Whether one check of the design at construction stands on a proved footing (the follow-up's direction 1). Under (b), a design that passes once passes at every c > 0. A per-step positive-definiteness check can then catch only a loss in the computed factor, never a change in the system's identification.

**Lean (proved).** The statements are over Mathlib's `Matrix.PosSemidef` and `Matrix.PosDef`, through `posSemidef_iff_dotProduct_mulVec` and `posDef_iff_dotProduct_mulVec`. `comb S c` is A(c) = Σₖ cₖ • Sₖ.
- (a): `form_comb_eq_zero_iff` and `comb_mulVec_eq_zero_iff`.
- (b): `comb_posDef_iff` and `comb_posDef_iff_of_pos`.
- (c): `penalized_det_eq_zero_iff`, for a rectangular Z. It goes through a family indexed by `Option J`, whose extra term is ZᵀZ with coefficient 1.
- Cauchy–Schwarz for a positive semidefinite form is `form_sq_le`, from `discrim_le_zero`. The proof's second step is `mulVec_eq_zero_of_form`. Mathlib's `PosSemidef.dotProduct_mulVec_zero_iff` states it too, but its module imports the C*-algebra order.

**Prior art.** Textbook linear algebra; no source consulted.

## PA-6 — the bound

**Status:** proved (`form_comb_bound`, `lamMin_comb_ge`, `lamMax_comb_le`, `cond2_comb_le`, `cond2_diagPair`, `bracket_one`). Checks PA-6(i) to (iv) pass.

**Claims.** With ρₖ = cₖ / c′ₖ, m = minₖ ρₖ and M = maxₖ ρₖ:
- **(a)** m·xᵀA(c′)x ≤ xᵀA(c)x ≤ M·xᵀA(c′)x for every x.
- **(b)** λ_min(A(c)) ≥ m·λ_min(A(c′)) and λ_max(A(c)) ≤ M·λ_max(A(c′)).
- **(c)** Where A(c′) is positive definite, κ₂(A(c)) ≤ (M/m)·κ₂(A(c′)).
- **(d) Sharpness.** M/m is attained by the smallest instance with c′ = (1, 1): κ₂(A(c′)) = 1 and κ₂(A(c)) = M/m. So no smaller constant serves every family.
- **(e)** A coefficient held at 1 in both, as ZᵀZ's is, enters with ratio 1, so m ≤ 1 ≤ M.

**Proof.** xᵀA(c)x = Σₖ ρₖ·(c′ₖ·xᵀSₖx), and each bracket is nonnegative. That gives (a). Then (b) by the Rayleigh characterization of the extreme eigenvalues, and (c) by dividing.

**Hypotheses, and what each does.**
- *Positive semidefinite terms* make each bracket nonnegative. With an indefinite term the ratios can fall outside [m, M].
- *(c) needs A(c′) positive definite.* By PA-5(b), A(c) then is too.

**Checks.**
- (a): exact rationals, 160 draws over eight random families, with rank-deficient Sₖ of both kinds.
- (b) and (c): 300 seeded float families, n from 4 to 15, with up to five terms and coefficients log-uniform over e^±7. The largest readings against a bound of 1 were 0.99969 (λ_min), 0.99972 (λ_max) and 0.99912 (κ₂).
- (d): exact, at three c.

**Decision.** Slice 6's screen rests on it. Sharpness says the window of c that the solve carries widens only with another form of the solve (the follow-up's direction 2), never with a better inequality.

**Lean (proved), through Rayleigh bounds.**
- *Definitions.* `lamMin A` and `lamMax A` are the infimum and supremum of xᵀAx / xᵀx over x ≠ 0, and `cond2 A` is their ratio.
  - For a symmetric A these are the extreme eigenvalues, so `cond2` is κ₂. The bridge is `lamMin_eq_iInf_eigenvalues` and `lamMax_eq_iSup_eigenvalues`, through Mathlib's `Matrix.IsHermitian.spectral_theorem` (`form_eq_sum_eigenvalues`, `eigenvector_facts`). It landed as a follow-up (`7f535b3`). The results above do not use it.
- *The ratios enter as a bracket:* m·c′ₖ ≤ cₖ ≤ M·c′ₖ for every k. With c′ > 0 that is m ≤ ρₖ ≤ M, and min ρ and max ρ are the tightest bracket.
- (a): `form_comb_bound`.
- (b): `lamMin_comb_ge` and `lamMax_comb_le`.
- (c): `cond2_comb_le`. It uses `lamMin_pos`: a positive definite matrix has a positive least Rayleigh quotient, proved through its inverse with no spectral theorem.
- (d): `cond2_diagPair`. `diagPair_posSemidef` and `comb_diagPair_posDef` witness that the instance meets the hypotheses: each term is singular, yet A(c) is positive definite.
- (e): `bracket_one`.

**Prior art.** Elementary. The defects brief's slice 6 states (a) to (c) as derived, with no source.

## PA-7 — from the bound to what the gate reads

Restated on 2026-11-14 (`90fa6ce`) from the defects chip's chain; finding 1 says how it differs from the first statement.

**Status:** proved (`diag_comb_bound`, `lamMin_unitScale_ge`, `norm1_le_card`, `norm1_le_sqrt_card`, `mulVec_sq_le_norm1`, `cond1_unitScale_le`). The follow-up proves (c2′) and (d′) below (`norm1_le_half_one_add_sqrt`, `norm1_columnWitness`, `cond1_unitScale_le_sharp`). Checks PA-7(a) to (d), (c2′) and (d′) pass.

**Claims.** With A(c′) positive definite (so A(c) is too, PA-5(b)), ρ against c′, and q bounding the nonzeros in any column of A(c):
- **(a) The diagonal.** Every (Sₖ)ᵢᵢ ≥ 0, so m·A(c′)ᵢᵢ ≤ A(c)ᵢᵢ ≤ M·A(c′)ᵢᵢ.
- **(b) The scaled matrix's smallest eigenvalue.** λ_min(E(c)) ≥ (m/M)·λ_min(E(c′)).
- **(c) Norms.**
  - **(c1)** ‖E(c)‖₁ ≤ q, from a unit diagonal, entries of at most 1 in absolute value, and at most q nonzeros in a column.
  - **(c2)** For any n × n matrix, ‖A‖₁ ≤ √n·‖A‖₂. A column of ones attains it.
  - **(c3)** For a symmetric matrix, ‖A‖₂ ≤ ‖A‖₁.
- **(d) Composed.** κ₁(E(c)) ≤ q·√n·(M/m)·‖E(c′)⁻¹‖₁.

**Proof.**
- (a) is PA-6(a) at x = eᵢ.
- (b): write xᵀE(c)x = wᵀA(c)w with wᵢ = xᵢ/√A(c)ᵢᵢ. PA-6(a) gives wᵀA(c)w ≥ m·wᵀA(c′)w = m·vᵀE(c′)v, with vᵢ = xᵢ·√(A(c′)ᵢᵢ / A(c)ᵢᵢ). By (a), ‖v‖² ≥ ‖x‖²/M, so xᵀE(c)x ≥ (m/M)·λ_min(E(c′))·‖x‖².
- (c1): E(c) is positive semidefinite with a unit diagonal, so Eᵢⱼ² ≤ EᵢᵢEⱼⱼ = 1, by Cauchy–Schwarz for the form. A column's absolute sum is then at most its count of nonzeros, which is A(c)'s.
- (c2): column j's absolute sum is at most √n·‖Aeⱼ‖, by Cauchy–Schwarz.
- (c3): Schur's test. (Ax)ᵢ² ≤ (Σⱼ|aᵢⱼ|)(Σⱼ|aᵢⱼ|xⱼ²), and summing over i gives ‖Ax‖² ≤ ‖A‖∞‖A‖₁‖x‖², where ‖A‖∞ = ‖A‖₁ for a symmetric A.
- (d): κ₁(E(c)) = ‖E(c)‖₁·‖E(c)⁻¹‖₁.
  - ‖E(c)‖₁ ≤ q, by (c1).
  - ‖E(c)⁻¹‖₂ = 1/λ_min(E(c)), so (c2) gives ‖E(c)⁻¹‖₁ ≤ √n/λ_min(E(c)).
  - (b) gives λ_min(E(c)) ≥ (m/M)·λ_min(E(c′)).
  - λ_min(E(c′)) = 1/‖E(c′)⁻¹‖₂ ≥ 1/‖E(c′)⁻¹‖₁, by (c3) for the symmetric E(c′)⁻¹.

**Hypotheses, and what each does.**
- *A(c′) positive definite* makes every diagonal entry positive, so E(c) and E(c′) are defined, and makes λ_min(E(c′)) > 0.
- *Every cₖ > 0.* The hierarchical model's ξ² can reach 0. There m = 0 and the bound is infinite, and the screen hands that step to the estimate.
- *q counts A(c)'s nonzeros.* For c > 0 the pattern of A(c) is the union of the Sₖ's patterns, barring cancellation, so one count from the design serves every c.

**Checks.**
- (a) to (d): 240 banded families, 80 of each of BAPS's three penalty forms, with n from 6 to 29, a cubic B-spline Gram and second-difference penalties, and coefficients log-uniform over e^±7. No piece was violated.
  - The largest κ₁(E(c)) as a share of (d)'s bound was 0.145.
  - Finding 2 is the share of (M/m)·κ₁(E(c′)).
- (c2) and (c3): 300 random matrices per line. A column of ones gives 3.000 at n = 9.

**Outside these hypotheses:** a coefficient at 0, handed to the estimate by the screen.

**Outside Lean: rounding.**
- The screen multiplies the bound by 2 for it, and the defects chip's sweep tests that.
- The per-step estimate is a lower bound on κ₁ by its construction, which is the direction the screen needs.
- The reference's ‖E(c′)⁻¹‖₁ is read exactly, from n solves, so no estimator's shortfall enters.

**Decision.** Whether the screen admits only what the estimate admits, by proof in exact arithmetic, with rounding alone left to the sweep. The defects brief gets a line it can take, and the defects chip a message when the chain is settled.

**Lean (proved).** Elementary: finite sums, Cauchy–Schwarz for a positive semidefinite form and for a finite sum, and Mathlib's matrix inverse.
- *Definitions.* `unitScale A` is E = D^{-1/2}AD^{-1/2}. `norm1` is the largest column absolute sum, and `cond1 A = norm1 A · norm1 A⁻¹`.
- (a): `diag_comb_bound`.
- (b): `lamMin_unitScale_ge`.
- (c1): `norm1_le_card`, through `abs_entry_le_one`.
- (c2): `norm1_le_sqrt_card`.
- (c3): `mulVec_sq_le_norm1`.
- (d): `cond1_unitScale_le`.
  - It uses `inv_norm1_le_lamMin` for λ_min(E(c′)) ≥ 1/‖E(c′)⁻¹‖₁ > 0.
  - Its q bounds the nonzeros in every column of A(c), and E(c) has no more.
- λ_min is the infimum of the Rayleigh quotient, as in PA-6.
- ‖A‖₂ is the least β with ‖Ax‖ ≤ β‖x‖ for every x, so the norm lines are stated for such a β, in squares. (c2) reads: if ‖Ax‖ ≤ β‖x‖ for every x, then ‖A‖₁ ≤ √n·β. (c3) reads: for a symmetric A, ‖Ax‖ ≤ ‖A‖₁‖x‖ for every x. No operator norm is defined.

**Prior art.** The chain is the defects chip's (`76851ee`). Its pieces are textbook; the norm lines are standard comparisons. (c3) is Schur's test.

**Follow-up (2026-11-14): the column step sharpened.**
- **(c2′)** For P symmetric positive semidefinite with zᵀPz ≤ μ‖z‖² for every z: ‖P‖₁ ≤ μ·(1 + √n)/2. Attained: P = vvᵀ/12 with v = (3, 1, 1, 1) has λ_max = 1 and ‖P‖₁ = 3/2 = (1 + √4)/2.
- **(d′)** With the hypotheses of (d): κ₁(E(c)) ≤ q·((1 + √n)/2)·(M/m)·‖E(c′)⁻¹‖₁. (d) follows from it, since (1 + √n)/2 ≤ √n.
- *Proof of (c2′).*
  - For u = Peⱼ: ‖u‖² ≤ μ·uⱼ (`sq_norm_le_form`), so Σ_{i≠j} uᵢ² ≤ μuⱼ − uⱼ², and Cauchy–Schwarz gives Σ_{i≠j}|uᵢ| ≤ √(n−1)·√(μuⱼ − uⱼ²) = s·r.
  - With x = uⱼ − μ/2: x² + r² = μ²/4, and (x + sr)² ≤ (1 + s²)(x² + r²) = nμ²/4.
  - So uⱼ + sr ≤ μ/2 + √n·μ/2.
- *Proof of (d′).* As (d), with (c2′) for P = E(c)⁻¹ and μ = 1/λ_min(E(c)) (`form_inv_le`) in place of (c2).
- *Checks:* 400 random PSD matrices, none above the constant. The rank-one witnesses at n = 4, 9, 16 and 25 are exact. The unit-diagonal Gram construction reads (1 + √n)/2 to five places at n = 4, 9, 16, 25 and 49. On 120 banded families, κ₁(E(c)) reached at most 0.240 of (d′)'s bound.
- *Lean:* `norm1_le_half_one_add_sqrt`, `form_inv_le`, `norm1_columnWitness` and `cond1_unitScale_le_sharp`. `cond1_unitScale_le` is now its corollary, with its statement unchanged.
- *Decision:* the proved part of S becomes q·(1 + √n) with the defects chip's factor 2 for rounding, in place of 2·q·√n.
- *The disproof* (finding 6): no constant bounds the inverse part over every family, so the order of (d′) cannot be improved without a further hypothesis. `abs_le_norm1_inv_mul`, `unitScale_mulVec`, `spikeVec`, `spikeFamily`, `spikeVec_dot`, `spikeFamily_posSemidef`, `comb_spike_one`, `comb_spike_two`, `norm1_inv_spike` and `inverse_part_unbounded` are in the module. A sparse pattern does not remove it either: a tridiagonal family shows the same growth. What a bound free of n would need is a hypothesis on M/m against the pattern's spectral gap (finding 6).

**Off the path: Jacobi scaling.** The first statement's route. Recorded here and not queued: no decision rests on it now.
- *The statement:* for B positive definite with unit diagonal and at most q nonzeros in any row, and F any positive diagonal, κ₂(B) ≤ q·κ₂(FBF).
- *The sketch:*
  - |Bᵢⱼ| ≤ 1, so λ_max(B) ≤ q;
  - λ_max(FBF) ≥ maxᵢ fᵢ²;
  - λ_min(B) ≥ λ_min(FBF) / maxᵢ fᵢ².
- *Checks:* 400 seeded banded B with F log-uniform over e^±7. Every step held, and the largest κ₂(B) / (q·κ₂(FBF)) was 0.064. Finding 3 is the local-search reading.
- **van der Sluis, "Condition numbers and equilibration of matrices", Numer. Math. 14:14–23 (1969). Read** (the PDF on the shelf is filed as `PDFs/other/VanderSluis1963.pdf`).
  - Theorem 4.1: for P Hermitian positive definite with all diagonal elements equal, κ₂(P) ≤ n·min over nonsingular diagonal D of κ₂(DᴴPD).
  - Theorem 4.3: if moreover P has at most q nonzero elements in any row, the factor is q. That is this statement, with the unit diagonal as the equal diagonal.
  - His Theorem 4.2 is a √q lemma of another kind: ‖A‖₂ ≤ √q·maxᵢ ‖rowᵢ(A)‖₂ when A has at most q nonzeros in any column.
- **Bradley and Murray, arXiv 1110.2805. Read in its LaTeX source.** It cites Theorems 4.1 and 4.3 with the same numbers and statements.
- **Higham, *Accuracy and Stability of Numerical Algorithms* (2002), Corollary 7.6. Not on the shelf; read at second hand.** NAG's routine document for F07FFF (DPOEQU), on the shelf as `PDFs/other/LapackF07.pdf`, says this. The scaling sⱼ = 1/√aⱼⱼ gives B = SAS a unit diagonal, which puts κ₂(B) within a factor n of the smallest over all diagonal scalings, and it cites Corollary 7.6 for that. That is van der Sluis's 4.1.
