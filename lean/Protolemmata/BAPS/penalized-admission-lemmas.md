# A penalized system's admission: identification free of the coefficients, the bound, the chain to the gate's κ₁

BAPS factors A(c) = Σₖ cₖSₖ at every sampler step, and a gate (`BapsAdmission`) admits a factor that is positive definite with an estimated κ₁ of its unit-diagonal scaling E(c) at most 10¹⁴. Slice 6 of the [defects brief](../../../issues/regression/shared-code-defects-brief.md#slice-6-the-samplers-condition-check-as-a-setting) builds a screen on a bound relative to one reference point, and the [BAPS follow-up's admission input](../../../issues/regression/bass/extensions-brief.md#input-for-the-follow-up-admission-and-singular-systems-2026-11-14) states the two failures behind the gate. These are the statements both rest on. They were commissioned on 2026-11-14, as targets 5 to 7 of the [chip brief](../BASS/active-subspace-lean-chip-brief.md#added-targets-a-penalized-systems-admission-2026-11-14). The coordinating session derived them that day and checked them against no run. The numeric checks are [`penalized-admission-checks.py`](penalized-admission-checks.py).

**Status (2026-11-14).** Posed and checked, 17 of 17 pass. No Lean yet. Each target's own status line follows.

**What posing changed** (findings before any Lean):
1. **PA-7(c): the chain holds, and the same pieces give a constant √q smaller.** Two of the posed chain's steps spend more than they must.
   - *‖E(c)‖₁ ≤ q directly.* E has a unit diagonal and entries of at most 1 in size, so each column's absolute sum is at most its count of nonzeros. The posed chain goes through ‖E‖₁ ≤ √q·‖E‖₂ ≤ √q·λ_max(E), and λ_max(E) can itself reach q.
   - *The Jacobi step is not needed.* The diagonal of A(c) is A(c_ref)'s times factors in [m, M], which is PA-6(a) at x = eᵢ. So λ_min(E(c)) ≥ (m/M)·λ_min(E(c_ref)) directly, without passing through κ₂(FBF).
   - *The result:* κ₁(E(c)) ≤ q·√n·(M/m)·‖E(c_ref)⁻¹‖₁ ≤ q·√n·(M/m)·κ₁(E(c_ref)). The last step holds because ‖E(c_ref)‖₁ ≥ 1.
   - The defects chip reached the first form on its own (its slice 6 entry, with m for q and ν for n). Its screen reads ‖E(c_ref)⁻¹‖₁ exactly, so that is the form the screen rests on.
   - The posed q·√(qn) follows from both forms, and all three are posed below.
2. **What the banded families use of the proved factor (a reading).** The checks ran 240 families built from a cubic B-spline Gram and second-difference penalties, with n from 6 to 29 and q up to 7. The largest κ₁(E(c)) / ((M/m)·κ₁(E(c_ref))) was 1.007, against a proved q·√n of up to 38. Whether a smaller factor holds for every family is open. The sweep measures the slack.
3. **Van der Sluis's q in practice (a reading).** On three banded unit-diagonal matrices, a local search over diagonal scalings found κ₂(B) / min_F κ₂(FBF) ≈ 1.000, 1.013 and 1.105, against q = 3, 5 and 7.

**Symbols.** n the order. Sₖ are n × n real symmetric positive semidefinite, for k in a nonempty finite index set. c, c′ and c_ref are coefficient vectors with every entry positive, and A(c) = Σₖ cₖSₖ. ρₖ = cₖ / c′ₖ, m = minₖ ρₖ and M = maxₖ ρₖ. λ_min and λ_max are the extreme eigenvalues of a symmetric matrix, and κ₂ = λ_max / λ_min where it is positive definite. ‖A‖₁ is the largest column absolute sum, and κ₁(A) = ‖A‖₁‖A⁻¹‖₁. E(c) = D^{-1/2}A(c)D^{-1/2} with D = diag A(c), the unit-diagonal scaling. q bounds the nonzeros in any row (equally, any column) of every A(c).

**The smallest instance**, used throughout: S₁ = diag(1, 0) and S₂ = diag(0, 1). Each is singular, they share no null vector but 0, and A(c) = diag(c₁, c₂).

---

## PA-5 — identification does not read the coefficients

**Status:** posed and checked (PA-5(i) to (iii) pass).

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

**Lean (planned).** Statements over Mathlib's `Matrix.PosSemidef` and `Matrix.PosDef`, through `posSemidef_iff_dotProduct_mulVec` and `posDef_iff_dotProduct_mulVec`. Cauchy–Schwarz for a positive semidefinite form is proved from `discrim_le_zero`. Mathlib's `PosSemidef.dotProduct_mulVec_zero_iff` states the second step too, but its module imports the C*-algebra order.

**Prior art.** Textbook linear algebra; no source consulted.

## PA-6 — the bound

**Status:** posed and checked (PA-6(i) to (iv) pass).

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

**Lean (planned).** (a) is the core. The eigenvalue lines go through Rayleigh bounds: λ_min and λ_max are defined as the infimum and supremum of xᵀAx / xᵀx over x ≠ 0. A bridge to Mathlib's `Matrix.IsHermitian.eigenvalues` is to be added if `Mathlib.Analysis.Matrix.Spectrum` serves, and the report says which landed.

**Prior art.** Elementary. The defects brief's slice 6 states (a) to (c) as derived, with no source.

## PA-7 — from the bound to what the gate reads

**Status:** posed and checked (PA-7(i) to (iii) pass). Finding 1 sharpens (c).

**Claims.**
- **(a) Jacobi scaling.** Let B be positive definite with unit diagonal and at most q nonzeros in any row, and F any positive diagonal. Then κ₂(B) ≤ q·κ₂(FBF). The sketch:
  - |Bᵢⱼ| ≤ 1 (a 2 × 2 principal minor), so every row's absolute sum is at most q, and λ_max(B) ≤ q;
  - λ_max(FBF) ≥ maxᵢ fᵢ², from the diagonal entry fᵢ²;
  - λ_min(B) ≥ λ_min(FBF) / maxᵢ fᵢ², by writing x = Fy.
- **(b) Norms.**
  - For a symmetric matrix, ‖A‖₂ ≤ ‖A‖₁.
  - With at most q nonzeros in a column, ‖A‖₁ ≤ √q·‖A‖₂, by Cauchy–Schwarz on the column's support and ‖Aeⱼ‖ ≤ ‖A‖₂.
  - For any n × n matrix, ‖A‖₁ ≤ √n·‖A‖₂. A column of ones attains it.
- **(c) Composed.** With A(c_ref) positive definite and q bounding the nonzeros in a row of every A(c):
  - **(c1)** κ₁(E(c)) ≤ q·√n·(M/m)·‖E(c_ref)⁻¹‖₁, with ρ taken against c_ref. This is finding 1's form and the defects chip's.
  - **(c2)** κ₁(E(c)) ≤ q·√n·(M/m)·κ₁(E(c_ref)), from (c1) since ‖E(c_ref)‖₁ ≥ 1.
  - **(c3)** κ₁(E(c)) ≤ q·√(qn)·(M/m)·κ₁(E(c_ref)), as posed. It follows from (c2), and also from (a), (b) and PA-6 along the posed chain.

**Proof of (c1).**
- ‖E(c)‖₁ ≤ q, by the unit diagonal and |Eᵢⱼ| ≤ 1.
- ‖E(c)⁻¹‖₁ ≤ √n·maxⱼ ‖E(c)⁻¹eⱼ‖ ≤ √n / λ_min(E(c)).
- Write xᵀE(c)x = yᵀÃ(c)y, where yᵢ = xᵢ·√(A(c_ref)ᵢᵢ / A(c)ᵢᵢ) and Ã(c) = D_ref^{-1/2}A(c)D_ref^{-1/2}.
  - PA-6(a) at eᵢ gives A(c)ᵢᵢ ≤ M·A(c_ref)ᵢᵢ, so ‖y‖² ≥ ‖x‖²/M.
  - PA-6(a) at D_ref^{-1/2}y gives yᵀÃ(c)y ≥ m·yᵀE(c_ref)y ≥ m·λ_min(E(c_ref))‖y‖².
  - So λ_min(E(c)) ≥ (m/M)·λ_min(E(c_ref)).
- λ_min(E(c_ref)) ≥ 1/‖E(c_ref)⁻¹‖₂ ≥ 1/‖E(c_ref)⁻¹‖₁, by (b) for the symmetric E(c_ref)⁻¹.
- The congruence by D_ref^{-1/2} keeps each Sₖ positive semidefinite. That is how PA-6 enters, in the brief's words. In the proof it is PA-6(a) evaluated at a scaled vector.

**Hypotheses, and what each does.**
- *Unit diagonal and positive definiteness* give |Bᵢⱼ| ≤ 1. Without the unit diagonal, (a)'s first step fails.
- *q bounds every A(c)'s pattern.* For c > 0 the pattern of A(c) is the union of the Sₖ's patterns, barring cancellation, so one q serves every c.

**Checks.**
- (a): 400 seeded unit-diagonal banded B, n from 3 to 13 and half-bandwidth 0 to 3, with F log-uniform over e^±7. Each sketch step was checked, and the largest κ₂(B) / (q·κ₂(FBF)) was 0.064. Finding 3 is the local-search reading.
- (b): 300 random matrices per line. A column of ones gives 3.000 at n = 9.
- (c): 240 banded families, 80 of each of BAPS's three penalty forms, with no violation of (c1), (c2) or (c3). The largest κ₁(E(c)) as a share of each bound was 0.145 for (c1), 0.054 for (c2) and 0.021 for (c3). Each piece of (c1)'s proof was checked too.

**Outside Lean.** The gate reads an estimate of κ₁ that errs low, with no proved bound on its shortfall. That factor stays with the defects chip's sweep. Its slice 6 entry reads ‖E(c_ref)⁻¹‖₁ exactly, from n solves once per run, which removes the shortfall at the reference.

**Decision.** Which part of S is proved and which is measured. The proved part is (c1)'s q·√n against ‖E(c_ref)⁻¹‖₁, or q·√n against κ₁(E(c_ref)), and the factor 2 for rounding is the defects chip's. The defects brief gets a line it can take, and the defects chip a message when the chain's constant is settled.

**Lean (planned).** Elementary: finite sums, Cauchy–Schwarz for a positive semidefinite form and for a finite sum, and Mathlib's matrix inverse. ‖·‖₂ enters only through λ_max and λ_min of symmetric matrices, so no operator norm is defined.

**Prior art.**
- **van der Sluis, "Condition numbers and equilibration of matrices", Numer. Math. 14:14–23 (1969). Read** (the PDF on the shelf is filed as `PDFs/other/VanderSluis1963.pdf`).
  - Theorem 4.1: for P Hermitian positive definite with all diagonal elements equal, κ₂(P) ≤ n·min over nonsingular diagonal D of κ₂(DᴴPD).
  - Theorem 4.3: if moreover P has at most q nonzero elements in any row, the factor is q. That is (a), with the unit diagonal as the equal diagonal.
  - His Theorem 4.2 is a different √q lemma: ‖A‖₂ ≤ √q·maxᵢ ‖rowᵢ(A)‖₂ when A has at most q nonzeros in any column. (b)'s √q line is not his.
- **Bradley and Murray, arXiv 1110.2805. Read in its LaTeX source.** It cites Theorems 4.1 and 4.3 with the same numbers and statements.
- **Higham, *Accuracy and Stability of Numerical Algorithms* (2002), Corollary 7.6. Recalled, not read;** the book is not on the shelf. From memory, the factor n for symmetric positive definite A scaled by diag(aᵢᵢ^{-1/2}), van der Sluis's 4.1 in Higham's §7.3.
