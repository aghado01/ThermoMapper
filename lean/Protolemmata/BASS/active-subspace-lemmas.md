# Active-subspace geometry — P̄ over mixed ranks, the matrix Poincaré check, E36's direction moves, rotate-and-refit

These are the statements that E8's readouts ([active-subspace-brief.md](../../../issues/regression/bass/active-subspace-brief.md), slices 5 and 7) and E36's sampler ([bppr-brief.md](../../../issues/regression/bass/bppr-brief.md), decision 5) rest on. They were commissioned on 2026-11-13 by the [chip brief](active-subspace-lean-chip-brief.md); before that the four existed only as names (handoff of 2026-11-10, item 8). The numeric checks are [`active-subspace-checks.py`](active-subspace-checks.py).

**Status (2026-11-13).** All four are proved, with no apologies, in `lean/Enthymemata/ActiveSubspace/`: `ProjectionMean` (AS-1), `MatrixPoincare` (AS-2), `DirectionMoves` (AS-3) and `Refit` (AS-4). The barrel is `Enthymemata/ActiveSubspace.lean` and the namespace `ActiveSubspace`. Each main result prints only `propext`, `Classical.choice` and `Quot.sound`. The checks pass 44 of 44. The extension (below) closed the two hypotheses the first round left: the uniform's one-variable inequality, now proved on the continuous product law (AS-2(d)), and the reversible-jump moves in measure form (AS-3(c)). The truncated normal's constant still rests on the checks. Each target's own status line follows.

**What posing changed** (findings before any Lean):
1. **AS-1(c): the cross-rank distance is not Ye and Lim's.** ½‖P_s − P‖²_F = |k_s − k|/2 + Σ sin²θᵢ holds (checked to 5·10⁻¹⁵). But Ye and Lim's chordal metric on subspaces of all dimensions is |k − l| + Σ sin²θᵢ, with the rank term at weight 1. That is their Table 3, and their §7 shows it equals the symmetric directional distance. The two agree when the ranks are equal.
   - The rank summary's loss is the Frobenius one, so slice 5's readout is consistent with itself. Only the attribution changes.
   - Under Ye and Lim's metric the expected loss would be E[max(k, k_s)] − tr(P P̄). Its minimizer is not the ½ rule.
2. **AS-1(b): Barbieri and Berger prove optimality for a different loss.** Their median probability model keeps the variables of inclusion probability ≥ ½, and they prove it optimal for predictive squared error, in normal linear models under conditions ("often", in their abstract). On coordinate subspaces the Frobenius loss is the Hamming loss, and its Bayes action is the ½ rule with no conditions. The name fits the rule; the reason it is optimal here is not their theorem.
3. **AS-3: at a = 1 the change move never flips a sign.** On S⁰ = {±μ} the power spherical law puts (1 − 1)^κ = 0 on −μ, so a one-coordinate direction keeps its sign and only t₀ moves.
   - Balance holds, since the proposal is the identity on θ there.
   - The paper's normalizer has Γ(β) with β = (a − 1)/2 = 0 at a = 1, so a sampler must treat a = 1 apart.
   - A one-coordinate ridge changes sign only by death and rebirth.
4. **BPPR's printed birth ratio has a slip at a\* = 1.** Its eq. (mhb) divides by Wallenius' W_{a\*}(𝒥 | υ) at every a\*, but its Algorithm 2 draws 𝒥 uniformly from {1, …, p} when a\* = 1. There the proposal term is 1/p, and it cancels π(𝒥 | a\* = 1) = 1/p.
5. **AS-4(b): "the subspace moves under any other" is true only of a truncation.**
   - Under every invertible L, the exact inactive subspace maps by L⁻¹: ker(LᵀCL) = L⁻¹ ker C. So for an exact k-index function (rank C = k), the leading k is equivariant whatever the metric.
   - What moves under a non-conformal L is the leading k of a C of rank above k.
   - "Orthogonal" also widens to conformal (L = cQ).
6. **AS-4(a): the frame claim is stronger than posed.** A linear map R of u carries the tensor-product hinge class into itself only when each row of R has at most one nonzero entry. For orthogonal R that means R is a signed permutation, which BASS already ignores. So the refit depends on the frame at every frame, and not only where the active block's eigenvalues are near-equal.
7. **AS-3(b), found in the proof: the count read twice changes the prior.** A birth ratio with a further 1/(M + 1) still balances, but against the count weight λ^M/(M!)², not Poisson's λ^M/M! (`double_count_ratio`). So the slip is silent: the chain still runs and passes reciprocity, while targeting a lighter count prior. It shows only against an independent reference, which is what AS-3c's balance gap measures (0.34, against the Poisson target).

**Which reading of AS-4 the record supports.** The frame. The name entered the record on 2026-11-10 (`391f052`, the handoff's item 8). At that point the record held slice 7's plan, which says BASS's axis alignment makes the frame matter (decision 2, after the Grassmann discussion's synthesis point 3), and the user's decision 2, which asks for the eigengaps within the active block. The scale reading's material came later: decision 4's gate in 7a and the training-range rule in 7d, both on 2026-11-11. Both readings are posed below and both hold.

**Symbols.** p inputs; P_s, P orthogonal projections (symmetric and idempotent) of ranks k_s, k; w_s ≥ 0 with Σ w_s = 1; P̄ = Σ w_s P_s with eigenpairs (μᵢ, vᵢ); ‖·‖_F the Frobenius norm. For AS-2: f the m outputs, J = ∂f/∂x (m × p), D diagonal. For AS-3: 𝒥 is the active index set (to keep J for the Jacobian), a = |𝒥|, and t₀ the knot. For AS-4: u = W₁ᵀx_box are the active variables, L a change of input coordinates, B a design, T a column map, and W, W′ the declared boxes' widths beside a training width.

---

## AS-1 — P̄ over mixed ranks

**Status (2026-11-13).** Proved in `ActiveSubspace/ProjectionMean.lean`, no apologies; checks AS-1a to AS-1c pass. The statement stands as posed. The Lean form of (c) holds at every k.

**Claims.** Each P_s is an orthogonal projection with k_s = tr P_s (its rank), the weights are nonnegative and sum to 1, and P̄ = Σ w_s P_s.
- **(a)** P̄ is symmetric, 0 ≼ P̄ ≼ I, and tr P̄ = Σ w_s k_s.
- **(b)** For every orthogonal projection P, of any rank: Σ_s w_s ‖P − P_s‖²_F = tr P + Σ w_s k_s − 2 tr(P P̄).
  - Writing P̄ = Σ μᵢ vᵢvᵢᵀ: tr P − 2 tr(P P̄) = Σᵢ (1 − 2μᵢ)·vᵢᵀPvᵢ, with 0 ≤ vᵢᵀPvᵢ ≤ 1. So the loss is at least Σ w_s k_s + Σᵢ min(0, 1 − 2μᵢ).
  - That bound is attained by the projection onto span{vᵢ : i ∈ S} for any S with {μ > ½} ⊆ S ⊆ {μ ≥ ½}. So the rank summary is least over projections of every rank, and an eigenvalue of exactly ½ is indifferent. No Ky Fan argument is needed: the diagonal bound does the work for every rank at once.
- **(c)** Let k ≤ k_s (the distance is symmetric). Let U be an orthonormal basis of P's range, and cos²θᵢ the k eigenvalues of UᵀP_sU, which lie in [0, 1]. Then ½‖P_s − P‖²_F = (k_s − k)/2 + Σ_{i≤k} sin²θᵢ.

**Hypotheses, and what each does.**
- *Symmetry of the projections* gives ‖P − Q‖²_F = tr P + tr Q − 2 tr(PQ). An oblique idempotent breaks it.
- *Σ w_s = 1* puts tr P in (b) with coefficient 1, and gives I − P̄ = Σ w_s(I − P_s) for the upper order.
- *w_s ≥ 0* gives the two orders.

**Decision.** The readout's rank summary (P̄'s eigenvalues above ½) and tr P̄ = E[k_s] stand as documented. The cross-rank distance stands as ½‖·‖²_F, with its citation corrected (finding 1).

**Prior art.**
- *Ye and Lim, SIAM J. Matrix Anal. Appl. 2016 (arXiv 1407.0900v3), read in its LaTeX source.*
  - Their Table 2 gives the equal-rank chordal distance (Σ sin²θᵢ)^{1/2} = ‖AAᵀ − BBᵀ‖_F/√2.
  - Their Theorem 12 makes the distance from a k-plane to the Schubert variety of l-planes Σ_{i≤min(k,l)} over the same per-angle terms, with no rank term. That is this readout's containment defect.
  - Their Table 3 makes the metric on all dimensions |k − l| + Σ sin²θᵢ. §7 shows it equals the symmetric directional distance max(k, l) − Σ_{ij}(aᵢᵀb_j)².
- *Barbieri and Berger, Ann. Statist. 2004 (arXiv math/0406464), from the abstract:* the median probability model, inclusion probability ≥ ½, is often optimal for predictive squared error in normal linear models.
- *The flag mean for mixed dimensions (Draper et al. 2014),* cited in the Grassmann discussion: from memory, unread.

**Checks** (exact rationals for the algebra, on 200 seeded posteriors of mixed ranks):
- (a): symmetry, idempotence, tr P_s = k_s and tr P̄ = Σ w_s k_s exactly; P̄'s eigenvalues within [−3·10⁻¹⁶, 1 + 10⁻¹⁵].
- (b): the identity exactly, at a projection of every rank. The summary was never beaten, by any eigen-subset or by 9,000 random projections (least margin −1.8·10⁻¹⁵). The ½ case reads ½ either way. On coordinate subspaces P̄ is the diagonal of inclusion probabilities, and the Frobenius loss is the Hamming loss.
- (c): to 4.7·10⁻¹⁵ on 300 pairs of mixed ranks. Ye and Lim's metric² exceeds the half Frobenius by |k_s − k|/2, to the same precision.

**Lean** (`Enthymemata/ActiveSubspace/ProjectionMean.lean`). A projection is `IsOrthProj` (Pᵀ = P, P² = P), and `frobSq` is the entry sum of squares.
- (a): `IsOrthProj.posSemidef`, `projMean_posSemidef`, `one_sub_projMean_posSemidef`, `trace_projMean`. `IsOrthProj.diag_mem` gives the diagonal in [0, 1], which (b) uses.
- (b): `projLoss_eq`, then `trace_sub_two_trace_mul_eq` (tr Q − 2 tr(QP̄) = Σᵢ (1 − 2μᵢ)(UᵀQU)ᵢᵢ) and `trace_sub_two_trace_mul_ge` (the bound, over projections of every rank). `summaryProj_isOrthProj` and `summaryProj_trace_sub` show it is attained for any S between {μ > ½} and {μ ≥ ½}. `summaryLoss_le` puts them together.
- (c): `half_frobSq_cross_rank` holds at every k, with cos²θᵢ the eigenvalues of UᵀP_sU, and `compress_eigenvalues_mem` puts those in [0, 1]. Taking U as the basis of the smaller rank gives the |k_s − k|/2 form; the norm is symmetric.
- *Not formalized:* that tr P_s is the rank, which for a projection is its count of unit eigenvalues.

---

## AS-2 — the matrix Poincaré check

**Status (2026-11-13).** (a), (b), the matrix form and (c)'s algebra are proved in `ActiveSubspace/MatrixPoincare.lean`, no apologies, reusing `efron_stein` (`Enthymemata/BassSensitivity.lean`, SG-1(d)). The one-variable inequality is a hypothesis here, as posed, in its integrated form (corrected after review, below). AS-2(d) in the extension proves it for the uniform and proves the matrix check on the continuous product of uniforms (`matrix_poincare_pi`). AS-2(e) proves the truncated normal's σ² (Bakry–Émery on an interval). Checks AS-2a to AS-2e pass.

**Claims.** ρ is a product law on the inputs, f has m outputs, and D = diag(c₁, …, c_p).
- **(a) The reduction.** Cov_ρ(f) ≼ E_ρ[J D Jᵀ] holds exactly when Var(vᵀf) ≤ Σᵢ cᵢ E[(∂ᵢ vᵀf)²] for every v ∈ ℝᵐ. This holds because vᵀCov(f)v = Var(vᵀf) and vᵀE[JDJᵀ]v = Σᵢ cᵢ E[(Σ_j v_j J_ji)²], and ∂ᵢ(vᵀf) = Σ_j v_j J_ji.
- **(b) The product step.** Efron–Stein gives Var g ≤ Σᵢ E[Varᵢ g]. If each coordinate satisfies Varᵢ h ≤ cᵢ Eᵢ[(∂ᵢh)²] for every h, with the other coordinates held fixed, then Var g ≤ Σᵢ cᵢ E[(∂ᵢ g)²]. Each coordinate keeps its own constant on D's diagonal, where a scalar bound would take max cᵢ.
- **(c) The constants.**
  - For the uniform on [a, b], c = (b − a)²/π², with equality at cos(π(x − a)/(b − a)) (Wirtinger's inequality, with the mean taken out).
  - In the box metric J_box = J·diag(bᵢ − aᵢ), so E[J D Jᵀ] = G/π² and Cov_ρ(f) ≼ G/π². The largest generalized eigenvalue is then at most 1.
  - For one normal component N(μ, σ²) truncated to [a, b], c = min(σ², (b − a)²/π²) is an upper bound on the constant: Bakry–Émery for curvature 1/σ², and Payne–Weinberger for a log-concave law on an interval.

**Hypotheses, and what each does.**
- *Independence* is used only in (b), through Efron–Stein, which fails under dependence. A joint sample has no such bound, and E8 does not check one.
- *Each coordinate's one-variable inequality* needs f absolutely continuous in that coordinate, with ∂ᵢf in L²: hinges of degree ≥ 1 and Legendre terms qualify.
- *Mathlib carries no Poincaré or Wirtinger inequality.* A search of the pinned Mathlib for Wirtinger, Poincaré inequality and Efron finds none. So in Lean the one-variable inequality is a hypothesis, stated over an abstract energy, and the finite core is (a) and (b). The constants rest on the checks below.

**Decision.** A reading above 1 beyond rounding marks a defect for every marginal that decision 6 and slice 6 cover:
- the uniform, whose constant is exact;
- one truncated normal component, whose constant is a valid upper bound. A reading there sits further below 1, and is less sharp.

It does not mark one for the laws E8 already reports as not applicable: mixtures of several normals, Student t, empirical marginals, and joint samples. For those, no constant is given, and a bimodal law's true constant exceeds both bounds (AS-2e: 248 against 0.0036 and 0.10). The moment route carries no Monte Carlo error, so "beyond rounding" is the threshold.

**Prior art.**
- *Efron and Stein (1981),* and the generator form already proved here (SG-1(d)).
- *Wirtinger's inequality* for the uniform constant (classical, from memory).
- *Bakry and Émery (1985), and Payne and Weinberger (1960)* in its weighted one-dimensional form: from memory, as slice 6 recorded, and checked below.
- *The matrix form* was derived in the Grassmann discussion (synthesis item 3). No source for it was found.

**Checks.**
- (a): both quadratic forms exactly on 100 finite laws.
- (b): Efron–Stein exactly on 400 product laws of 1 to 4 coordinates, with equality on every additive g and strict inequality on 227 of the rest.
- (c), the uniform: equality at the cosine to 10⁻¹⁴. The largest ratio on 300 random smooth functions is 0.954. A Galerkin supremum over degree 30 reads 1 to 10⁻¹⁴.
- (d), the matrix check: (cos πx₁, cos πx₂) reads 1 to 8·10⁻¹⁵, and 2 with D halved. On 60 random three-output fields on random boxes the worst reading is 0.991. The diagonal D reads 1.005 to 5.5 times the scalar max(D)'s reading, so it is the tighter bound.
- (e), truncated normals: the Galerkin constant over min(σ², (b − a)²/π²) is at most 1.000000 on 200 laws. Narrow normals inside their interval sit at equality, as Bakry–Émery's sharp case predicts.
  - A first version solved B⁻¹A directly and read 8.6 there. It came from ill-conditioning, and whitening B removes it. The script records why.

**Lean** (`Enthymemata/ActiveSubspace/MatrixPoincare.lean`), on a finite product law, with `BASS.meanOf` and `BASS.coordMean`.
- (a): `loewner_iff`, `covForm_eq_var`, `energyForm_eq`, through `quadForm_eq` and `meanOf_sq_lin`.
- (b): `meanOf_mul_sub_coordMean` (Efron–Stein's summand is E[Varᵢ g]); `var_le_of_coord_poincare`, from `BASS.efron_stein` and the hypothesis `hcoord : ∀ i x, condVar μ i g x ≤ c i * e i x`; and `coordwise_le_scalar`.
- *The matrix check:* `matrix_poincare`, with the one-variable inequality assumed in its integrated form for every linear combination vᵀf: Varᵢ(vᵀf) ≤ cᵢ Eᵢ[(Σₐ vₐJₐᵢ)²], both sides means over coordinate i (`coordMean`), so neither reads xᵢ. The first version bounded Varᵢ by the squared partial at x itself, pointwise in xᵢ. That form is true but stronger than the one-variable inequality, which cannot discharge it: on a symmetric box, g = xᵢ² has zero partial at 0 and nonzero variance. The coordinating session's review found it on 2026-11-13, and it is corrected.
- (c): `energy_box`, for any r; at r = 1/π² it is G/π².
- *Not formalized:* the constants themselves, the uniform's (b − a)²/π² and the truncated normal's bound. They rest on the checks.

---

## AS-3 — E36's direction moves

**Status (2026-11-13).** Proved in `ActiveSubspace/DirectionMoves.lean`, no apologies, reusing `mhRatio` and `mh_detailed_balance` (`Enthymemata/BassStructureChain.lean`, SC-1). Checks AS-3a to AS-3d pass. One statement changed: the count read twice changes the target, and still balances (finding 7).

**Claims.**
- **(a) The change move.** The power spherical density on S^{a−1} is q(x | μ) = N(κ, a)⁻¹(1 + μᵀx)^κ, with a normalizer that does not read μ. So on one sphere q(θ′ | θ) = q(θ | θ′).
  - The move keeps a and 𝒥, proposes θ′_𝒥 from q(· | θ_𝒥), and draws t₀′ from its prior given θ′.
  - With θ's prior uniform on the sphere of 𝒥, its acceptance ratio is the collapsed likelihood ratio alone, p(y | B′, τ)/p(y | B, τ). No proposal term and no prior term remain: t₀'s prior depends on θ through the quantile bounds (L, U), and it cancels because the forward draw reads θ′ and the reverse reads θ.
- **(b) Birth and death.** A birth draws a, then 𝒥 given a, from proposals q(a), q(𝒥 | a) that read the current state only; θ from its prior on S^{a−1}_𝒥; and t₀ from its prior given θ. It inserts the ridge at a uniform position of the M + 1. Its ratio is

  α_birth = [p(y | B⁺, τ)/p(y | B, τ)] · [λ/(M + 1)] · [P_D(M + 1)/P_B(M)] · [π(a)π(𝒥 | a) / (q(a)q(𝒥 | a))],

  and the death's is its reciprocal, read from the state with M + 1 ridges.
  - *The count's factors.* The Poisson ratio π(M + 1)/π(M) = λ/(M + 1) remains. The insertion position's 1/(M + 1) cancels the death's uniform choice of ridge. On multisets, the multinomial coefficient does the same.
  - *What to leave out.* A further 1/(M + 1) reads the count twice (KHAOS's slip, in the port hub's register). The chain still balances, but against the count weight λ^M/(M!)² in place of Poisson's (finding 7).
  - *Dimension matching.* The new ridge's parameters are the proposal's draw itself, so the map (state, u) ↦ (state with u) has Jacobian 1. The densities of θ and t₀ enter as plain ratios, and cancel, provided prior and proposal are densities against one base measure: counting on (a, 𝒥), surface measure on S^{a−1}_𝒥, Lebesgue on t₀. A parameterization by angles would bring its Jacobian into both and cancel the same way, but only if both are written in it.
  - *A cap M_max* refuses births at the cap. P_D and P_B are the move-type probabilities at the two states, and both are ⅓ in the paper.
  - *The paper's display* agrees with this, except at a\* = 1 (finding 4). There π(a)π(𝒥 | a) = 1/(A·binom(p, a)), with A = min(3, p).

**Hypotheses, and what each does.**
- *t₀ is drawn from its prior given the proposed direction.* Drawing it given the old one breaks balance (AS-3b: gap 2·10⁻²).
- *The normalizer does not read μ.* That is the sphere's rotation invariance. A grid whose row normalizer reads θ breaks symmetry and balance (AS-3b, uneven circle).
- *The knots above t₀ are a function of (θ, t₀).* Then they carry no density.
- *β and σ² are integrated out,* as the paper's ratio does. A ratio that kept β would match K more dimensions per ridge.

**Decision: the ratios BPPR's first build slice codes.**
- The change move accepts on the likelihood ratio alone.
- The birth accepts on the form above, with BASS's own variable proposal for q(𝒥 | a) if decision 5 takes it.
- The death takes the reciprocal.
- a = 1 is handled apart in the sampler.

**Prior art.**
- *De Cao and Aziz (arXiv 2006.04437v2), read in its TeX:* the density, the normalizer 2^{α+β}π^β Γ(α)/Γ(α + β) with α = (a − 1)/2 + κ and β = (a − 1)/2, the Beta marginal, and the Householder sampler. Its symmetry about μ is stated there; the symmetry in the pair (μ, x) is this item's.
- *Collins, Francom and Rumsey (arXiv 2210.09181v1), §2.4, read in its TeX:* the birth's eq. (mhb), the death as its inverse, and the change move. The change's ratio is in their supplement, which is not on disk.
- *Green (1995)* for reversible jump: from memory.
- *SC-1 here:* `mhRatio_mul_swap` and `mh_detailed_balance`.

**Checks.**
- (a):
  - the normalizer integrates to one for a = 2 to 7 and κ = 0 to 1000 (log gap 2·10⁻¹³);
  - the sampler's cosine mean is within 1.35 standard errors;
  - the symmetry is exact.
- (b), the change move on 9 equally spaced directions with a θ-dependent t₀ prior: balance to 3·10⁻¹⁸. t₀ drawn given the old direction gives a gap of 2·10⁻²; an uneven grid gives 10⁻³.
- (c), birth and death exactly on multisets of 4 ridge types up to 3 ridges:
  - balance and stationarity both exact (0);
  - the same ratio on ordered tuples with uniform positions also exactly balanced;
  - an extra 1/(M + 1) gives a gap of 0.34 against the Poisson target.
- (d), a continuous toy with M ∈ {0, 1}, θ on S¹ and t₀ on [0, 1]. Births come from the prior, and the change move uses the power spherical at κ = 20. Over 400,000 steps it reads P(M = 1) = 0.838 against the exact 0.837, and E[cos θ | M = 1] = 0.530 against 0.534.

**Lean** (`Enthymemata/ActiveSubspace/DirectionMoves.lean`), reusing SC-1.
- (a): `powerSpherical_symm`, `mhRatio_of_symm`, `change_ratio` (any symmetric kernel, a uniform direction prior, the knot from its prior given the proposed direction), `change_detailed_balance`, `change_ratio_powerSpherical` (the same on `Metric.sphere 0 1` with the power spherical kernel), and `powerSpherical_antipode`.
- (b), on ordered tuples with insertion and deletion at uniform positions (`insertAt`, `tupleTarget` with any count weight w):
  - `birth_green_ratio`: the ratio is L′/L · w(M+1)/w(M) · P_D/P_B · π_r/q_r;
  - `poisson_count_ratio`: λ/(M + 1);
  - `double_count_ratio` (finding 7);
  - `prior_draw_cancel`;
  - `birth_death_balance`, through `mh_detailed_balance` on a two-point space.

*Not formalized:* the measure-theoretic form, on S^{a−1} × ℝ with densities. Mathlib's support for reversible jump and for sphere measure is not enough to state it short of a project. The finite statements carry the algebra of the ratio, and AS-3d checks the continuous form.

---

## AS-4 — rotate-and-refit's non-invariance

**Status (2026-11-13).** Both readings proved in `ActiveSubspace/Refit.lean`, no apologies; checks AS-4a to AS-4d pass. (viii) is proved on its example, with the leading eigenvectors stated (`IsLeading`); the general non-conformal case stays with the check.

### The frame (the record's reading)

**Claims.**
- **(i)** For a ≠ 0 and b ≠ 0, the ridge hinge (a u₁ + b u₂ − t)₊ is not Σ_{m<M} φ_m(u₁)ψ_m(u₂) for any M and any functions.
  - On the grid u₁ = −i/a, u₂ = (t + j + ½)/b, for i, j = 0 to M, its values are (j − i + ½)₊: an upper triangular M + 1 by M + 1 matrix with diagonal ½. It has rank M + 1, where any M-term separable sum has rank at most M.
- **(ii)** Every tensor-product hinge sum in u (BASS's class: an intercept plus terms that are products of hinge factors ±(u_l − t)₊^q, at most one per variable) is such a finite sum across the split {j} | rest, for every j.
- **(iii)** So a linear map R of u carries the class into itself only if each row of R has at most one nonzero entry. If row i has R_ij ≠ 0 and R_ik ≠ 0 with j ≠ k, the one-term hinge (u_i − t)₊ composed with R is not in the class. Restricted to the (u_j, u_k) plane it is (i)'s ridge.
  - For invertible R, such an R is a monomial matrix; for orthogonal R, a signed permutation.
  - Monomial maps do carry the class to itself. Hinges take both signs, a positive scale moves into the knot and the coefficient, and knots at data values move with the data. BASS already ignores those maps.
- **(iv)** So the refit's model class, and with it its prior over functions of x, depends on the frame W₁R within the active span unless R is a signed permutation. That holds at every frame, not only where eigenvalues tie.

**Decision: the eigengap report within the active block is the right warning, for what it warns of.** Where the active block's eigenvalues tie, R is free within the tied block, and (iv) says the refit then depends on an arbitrary choice. A wide gap says only that the frame is stable (Davis–Kahan). It does not say the frame suits the hinge class: E[C]'s eigenvectors are one frame among many that the class tells apart.

### The scale

**Claims.**
- **(v)** For x = Lz, ∇_z f = Lᵀ∇_x f, so C_z = LᵀC_xL.
- **(vi)** If L = cQ with Q orthogonal, the eigenvectors map by Qᵀ with eigenvalues times c². Every leading-k subspace is carried.
- **(vii)** For every invertible L, ker C_z = L⁻¹ ker C_x. So for rank C = k, the leading k of C_z is Lᵀ range C_x whatever L is (finding 5).
- **(viii)** If L is not conformal (its singular values are not all equal), some full-rank C has its leading direction moved. C = diag(1, ½) with L = diag(1, 2) gives LᵀCL = diag(1, 2), and the leading direction moves from x₁ to x₂. In general, C = U diag(c) Uᵀ is built on L's left singular vectors with c_j > c_i but σᵢ²cᵢ > σ_j²c_j.
- **(ix) The projected fit's rules.**
  - The g-form's hat B(BᵀB)⁻¹Bᵀ, and with it R² and the collapsed marginal, is invariant under B ↦ BT for every invertible T, and a column's scale is a diagonal T.
  - The count support reads counts.
  - The raw ridge's marginal moves with a column's scale, and the column mass Σᵢ B(xᵢ) scales with the column.
- **(x) The training-range rule.** For a factor in the box's units, with peak taken over the training range, the evaluated hinge (((x − t)/W)₊)^q / ((U − t)/W)^q equals ((x − t)₊)^q/(U − t)^q, which does not read W; downward factors are the same with t − L. A product of factors is invariant factor by factor. κ = Π ((U − t)/W)^q moves by Π (W′/W)^q.

**Decision: the training range as BASS's default (the handoff's item 2) rests on a proved invariance.**
- Under the rule, every column is the same function of the raw inputs under any declared box (x). So the column mass, the raw ridge, the g-form, and the fitted function all do not read the box.
- κ alone reads it, as 7d's gate measured: by (W/w)^q.
- What a default change would move is existing fits under column mass and the raw ridge, as the handoff says. The g-form with count support does not move under either rule (ix).

**Prior art.**
- *The rank argument* is the standard one: a separable sum is a low-rank kernel. It was not looked up.
- *Constantine's book (2015)* states the equivariance of C under orthogonal maps and the metric's role (from memory).
- *The paper's conclusion (2307.11241, §5)* names the projected law's dependence as its obstacle, which is a different non-invariance: of the law, not of the class.
- *Rotation-invariance failures of additive and tensor bases* are folklore in the MARS and projection-pursuit literature: from memory, no source read.

**Checks.**
- (i): exactly, on 200 random a, b ≠ 0 and t with n ≤ 7: triangular, diagonal ½, nonzero determinant.
- (ii): 100 random M-term hinge sums have grid rank ≤ M.
- (iii): the hinge in u₁ after a rotation by φ has rank 1 at 0°, 90°, 180° and 270°, and 7 to 27 otherwise, on a 30 × 30 grid.
- (vi) to (viii): orthogonal L carries eigenpairs to 2·10⁻¹⁴; ker C maps by L⁻¹ to 7·10⁻¹⁵; a rank-k C's leading k is carried to 10⁻¹⁵ under random L; and the SVD construction moved the leading direction in 200 of 200 cases.
- (ix): the hat exactly on 60 designs; the raw ridge's log marginal moves by 1.17 with a column scaled by 4.
- (x): exactly on 300 products of 1 to 3 factors, degrees 1 to 3, inside and beyond the training range. κ's ratio is Π (W′/W)^q exactly.

**Lean** (`Enthymemata/ActiveSubspace/Refit.lean`).
- *The frame:*
  - (i): `ridgeGrid_det`, `ridgeHinge_not_separable`;
  - (ii)–(iii): `IsHingeFactor`, `IsHingeSum`, `hinge_isHingeSum`, `hinge_comp_not_hingeSum` and `hingeSum_not_closed`;
  - the converse: `hingeFactor_neg`, `hingeFactor_scale` and `IsHingeFactor.monomial` at the factor level, and `IsHingeSum.comp_monomial` for the class. Every permutation of the variables with signs and positive scales keeps it, the factors' multiples moving into the coefficients.
- *The scale:*
  - (v): `outer_pullback`, `activeMatrix_pullback`;
  - (vi): `pullback_eigen_of_orthogonal`, `conformal_pullback`;
  - (vii): `pullback_ker`, `pullback_range`;
  - (viii): `IsLeading` (an eigenvector whose eigenvalue bounds every Rayleigh quotient), `leading_moves`, `leading_C_coord` and `leading_pullback_coord`. C's leading eigenvectors are the multiples of e₁, LᵀCL's those of e₂, and Lᵀe₁ is not leading for LᵀCL. The general non-conformal L stays with the check;
  - (ix): `hat_invariant` (it needs T invertible only), `rawRidge_moves`;
  - (x): `trainingRange_factor`, `trainingRange_hinge`, `kappa_moves`, `kappa_ratio`.

---

## Extension — the open items closed (commissioned by the user, 2026-11-13)

Four of the report's open items, taken in order of size. The numeric checks already cover each one: AS-1a for (d), AS-4b's SVD construction for (viii) in general, AS-3d for the continuous moves, and AS-2c to AS-2e for the constants.

**Status (2026-11-13).** All four are proved, with no apologies: AS-1(d), AS-4(viii) in general, AS-3(c), and AS-2(d) in all three parts. Every declaration prints only `propext`, `Classical.choice` and `Quot.sound`.

### AS-1(d) — a projection's trace is its rank

**Claim.** For Pᵀ = P and P² = P, rank P = tr P. Each eigenvalue satisfies μ² = μ, so it is 0 or 1, and tr P = Σμ counts the units, which is the rank.

**Lean (proved,** `ProjectionMean`**).** `IsOrthProj.eigenvalues_mem` (0 or 1) and `IsOrthProj.rank_eq_trace`, through `rank_eq_card_non_zero_eigs` and `trace_eq_sum_eigenvalues`. k_s = tr P_s is now a theorem.

### AS-4(viii) in general — every non-conformal L moves some leading direction

**Claim.** If LLᵀ is not a multiple of I, there are a positive semidefinite C and a v with `IsLeading C v` and ¬`IsLeading (LᵀCL) (Lᵀv)`.

**Construction, without SVD.** LLᵀ has eigenpairs with σᵢ² > σⱼ², and uᵢ, uⱼ orthonormal. Take C = cᵢuᵢuᵢᵀ + cⱼuⱼuⱼᵀ with cⱼ = 1 and cᵢ = (σᵢ² + σⱼ²)/(2σᵢ²), so cⱼ > cᵢ while σᵢ²cᵢ > σⱼ²cⱼ.
- uⱼ leads C, by Bessel for the orthonormal pair.
- Lᵀuⱼ is an eigenvector of LᵀCL with eigenvalue σⱼ²cⱼ.
- The Rayleigh quotient at Lᵀuᵢ is σᵢ²cᵢ, which is larger. So Lᵀuⱼ does not lead.
- No invertibility is needed: Lᵀuⱼ = 0 cannot lead.

**Lean (proved,** `Refit`**).** `leading_moves_of_not_conformal`, with the helpers `outer_mulVec`, `bessel_pair` and `eigvec_dot`. Its hypothesis is ∀ c, LLᵀ ≠ c·I. It reaches the real spectral theorem through `ProjectionMean`'s `eq_eigU`.

### AS-3(c) — the moves in measure form

**Claims.** On a measurable space with a σ-finite base measure ν, with the target π and the proposal q as densities against ν:
- **(i) Metropolis–Hastings is reversible in measure.** For measurable A and B, the flux from A to B equals the flux from B to A: ∫_A∫_B π(x)q(x, y)min(1, r(x, y)) dν(y)dν(x) is symmetric. The integrand is min(π(x)q(x, y), π(y)q(y, x)) for nonnegative π and q, and Tonelli swaps the order. With the sphere of the active coordinates and its surface measure as ν, this is the change move. The hypothesis that carries the weight is that prior and proposal are densities against the same ν.
- **(ii) Birth and death in measure.** For each position i, insertion (r, c) ↦ insertNth i r c carries ν ⊗ ν^M to ν^{M+1} (`measurePreserving_piFinSuccAbove`). That is the dimension matching, with Jacobian 1 as a measure fact. With π_M = w(M)Π pr(cⱼ) L_M(c) and q densities against ν, the flux from A ⊆ R^M to B ⊆ R^{M+1} through births at i equals the flux back through deaths at i.

**Lean (proved,** `DirectionMoves`**).** `mh_flux_symm`, for any s-finite ν. `birth_death_flux` per position, for a σ-finite ν, with `birthRatio`, `birthFlux`, `deathFlux` and the pointwise `birthFlux_eq_deathFlux`. The death flux's measurability is derived (`measurable_deathFlux`) from the measurability of pr, q, L_M and L_{M+1}, which `birth_death_flux` takes as hypotheses.

### AS-2(d) — the matrix check on continuous laws

**Claims.**
- **(i) Wirtinger on an interval.** For h continuous on [a, b] with derivative h′ (right derivatives at the ends do) and h′ square-integrable, Var_U(h) ≤ ((b − a)/π)² E_U[h′²] under the uniform U on [a, b].
  - Proof: even-extend h to period 2(b − a). The extension takes equal values at the ends, so `fourierCoeffOn_of_hasDeriv_right` has no boundary term, and |ĉₙ| ≤ (b − a)/(π|n|)·|ĉₙ(h′)|. Parseval (`hasSum_sq_fourierCoeffOn`) sums it, with the mean the n = 0 coefficient.
- **(ii) Efron–Stein on a product of probability measures,** for bounded measurable f: Var f ≤ Σᵢ E[Varᵢ f]. The route is `CommutingMeans` instantiated on bounded measurable functions under `Measure.pi`, with Eᵢ integrating out coordinate i. Its integral identity comes from `measurePreserving_piFinSuccAbove`, which reaches SG-1's continuous case as well.
- **(iii) The glue.** For f with C¹ sections along each coordinate, under a product of uniforms: Cov(f) ≼ E[J D Jᵀ] with D = diag((bᵢ − aᵢ)²/π²). That is (ii), then (i) on each section, then AS-2(a)'s reduction.

**Lean.**
- *(i), proved* in `ActiveSubspace/ContinuousPoincare`:
  - `fourierCoeffOn_of_hasDerivAt_off_countable` is the new piece. It is Fourier integration by parts with a derivative off a countable set, which a hinge's kinks need: Mathlib's form asks a one-sided derivative at every point.
  - `wirtinger_periodic` is built on it.
  - `wirtinger_interval` holds for h continuous on [a, b], with a derivative off a countable set that is measurable and bounded.
- *(ii), proved* in `ActiveSubspace/ProductVariance`:
  - `sq_integral_le`, and the two-block step `var_prod_le` (the law of total variance, then Jensen on the inner average, then Fubini);
  - `varM_comp_equiv` and `integral_pi_succ`;
  - `efron_stein_pi`, by induction over coordinates, splitting off coordinate 0 with `measurePreserving_piFinSuccAbove`. It holds for bounded measurable F under any finite product of probability measures.
- *(iii), proved* in `ContinuousPoincare`:
  - `uniformIcc`, with `integral_uniformIcc` and `varM_uniformIcc`;
  - `integral_integral_update`, which integrates out coordinate i;
  - `poincare_pi`, the scalar form;
  - `covM_form` and `energyM_form`;
  - `matrix_poincare_pi`: Cov(f) ≼ E[J diag(((bᵢ − aᵢ)/π)²) Jᵀ] under the product of uniforms.

  The hypotheses that remain are the fit's regularity only: bounded measurable outputs and partials, with sections along each coordinate continuous and differentiable off a countable set.
  - *The kink set may move with the other coordinates* (`s i x`), so a hinge in rotated coordinates (a projected refit read in x) is covered. The coordinating session's review found that the first form, with `s i` fixed, covered only axis-aligned kinks.
  - *Every hypothesis is asked only on the box* in `matrix_poincare_box`. It reads f and J through `clampBox`, which changes nothing under the law because the box has full measure (`ae_mem_box`). So a fit unbounded off the box, as every hinge is, is covered. The same review found that `matrix_poincare_pi` asks its bounds on all of ℝⁿ, which a hinge does not meet without that step.
- *The truncated normal's σ²* is proved in AS-2(e) below. The (b − a)²/π² half of its min, Payne–Weinberger, stays with the checks.

### AS-2(e) — Bakry–Émery on an interval

**Claims.** V on ℝ has derivatives V′ and V″ everywhere, with V″ continuous and V″ ≥ κ > 0 on [a, b]. h is continuous on [a, b], with a measurable derivative h′ off a countable set, bounded on [a, b].
- **(i)** If m is the e^{−V}-weighted mean, ∫ₐᵇ (h − m)e^{−V} = 0, then ∫ₐᵇ (h − m)²e^{−V} ≤ κ⁻¹ ∫ₐᵇ h′²e^{−V}. Divided by the normalizer: Var_π h ≤ κ⁻¹ E_π[h′²] for π ∝ e^{−V} on [a, b].
- **(ii)** For N(μ, σ²) truncated to [a, b], V = (x − μ)²/(2σ²) and κ = σ⁻², so the constant is σ² on every interval.

**Proof.** Let g = h − m and W(x) = ∫ₐˣ g e^{−V}. W vanishes at a, and at b by the mean's definition, so h needs no boundary condition.
- By parts, kinks allowed: ∫ g²e^{−V} = −∫ W h′.
- Ψ = V′W²e^{V} vanishes at both ends, and Ψ′ = V″W²e^{V} + 2V′Wg + V′²W²e^{V}. So ∫ g²e^{−V} − ∫ V″W²e^{V} = ∫ (g + V′We^{V})²e^{−V} ≥ 0.
- Pointwise, −Wh′ ≤ ½V″W²e^{V} + ½h′²e^{−V}/V″.

Together, ∫ g²e^{−V} ≤ ∫ h′²e^{−V}/V″ ≤ κ⁻¹ ∫ h′²e^{−V}. The middle bound is Brascamp–Lieb's form on the interval. It is not stated separately.

**Lean (proved,** `WeightedPoincare`**).** `poincare_weighted` is (i) and `poincare_truncNormal` is (ii). Both derivative identities go through `integral_eq_of_hasDerivAt_off_countable_of_le`, as AS-2(d)(i)'s step by parts did.

**Prior art.** Bakry and Émery (1985), from memory. The route through W and Ψ was derived here, with no source read.
