import Enthymemata.BgpCovering
import Enthymemata.BgpGaussianKL
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.MeasureTheory.Constructions.Pi

/-!
# BGP-4, UV-6 — posterior contraction with an unknown noise variance (enthymema)

The composed theorem of `Protolemmata/BARS/bgp-unknown-variance-lemmas.md`: under TWCD2025's
framework inequalities (59)–(61) for the GP prior and the proper prior [P1], with
`σ₀ ∈ (σ_min, σ_max)`, the posterior puts vanishing mass outside
`{‖f − f*‖ₙ + |σ − σ₀| ≤ M ε̄ₙ}`.

The statement is assembled from:
* the prior-mass condition, from UV-1 (`klDiv_le_and_evariance_le_of_slab`) and UV-2
  (`map_smul_closedBall`, `log_inv_mul_radius_le`), with the [P1] amplitude window and noise slab;
* the sieve and its remaining mass (UV-3, `prior_compl_sieve_le`);
* entropy (UV-4, with UV-4a's `isCover_image2_smul` and `exists_isCover_Icc`);
* the tests (UV-5, `exists_tests`), the forward boundary;
* the general theorem, `gvdv07_contraction`, a cited axiom (Ghosal & van der Vaart 2007,
  Theorem 4), in the style of `Stability.lean`.

**What is cited, and how.** TWCD2025's (59)–(61) and Remark A.1 enter as hypotheses of
`posterior_contraction_unknown_variance`: they are facts about the paper's GP prior and sieve,
which this file does not construct. Two of them carry the protolemma's "to confirm against the
paper" checks and are named for it: `h59` (UV-2: the (59) bound uniform over the amplitude box,
via Remark A.1) and `hRadius` (UV-4: the sieve's sup-norm radius `Rₙ` grows at most
polynomially). The remaining-mass exponent `C₂` of (60) must exceed the prior-mass exponent `C₁`
of (59) by `144/σ₀²` (`hC₂`); see the protolemma's correction under UV-3.

**The forward boundary.** UV-5's tests (`exists_tests`) are the one open obligation, so the
ledger reads `apologizing(1)`.
-/

open MeasureTheory ProbabilityTheory InformationTheory Filter Metric Set
open scoped ENNReal NNReal Topology

namespace BGP

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

section Model

/-! ### The model: fixed design, Gaussian noise with an unknown variance

`ev n i` is evaluation at the `i`-th of `n` design points, `g ↦ g(xᵢ)`, a continuous linear
functional on the function space `E`; with `‖ev n i‖ ≤ 1` the sup norm dominates the empirical
norm. A parameter is `θ = (f, σ)`. -/

variable (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ)

/-- The law of the data at `θ = (f, σ)`: independent `Yᵢ ~ N(f(xᵢ), σ²)`. -/
noncomputable def dataLaw (n : ℕ) (θ : E × ℝ) : Measure (Fin n → ℝ) :=
  Measure.pi fun i => gaussianReal (ev n i θ.1) (θ.2 ^ 2).toNNReal

/-- The likelihood of the data `y` at `θ`: the product of the Gaussian densities. -/
noncomputable def likelihood (n : ℕ) (θ : E × ℝ) (y : Fin n → ℝ) : ℝ≥0∞ :=
  ∏ i, gaussianPDF (ev n i θ.1) (θ.2 ^ 2).toNNReal (y i)

/-- The KL neighbourhood `B*ₙ(θ₀, ε; 2)` of Ghosal & van der Vaart (2007): the summed KL
divergences and summed log-ratio variances are at most `n ε²`. Both are taken in `ℝ≥0∞`, so a
parameter with an infinite divergence is excluded rather than counted as `0`. -/
def klBall (n : ℕ) (θ₀ : E × ℝ) (ε : ℝ) : Set (E × ℝ) :=
  {θ | ∑ i, klDiv (gaussianReal (ev n i θ₀.1) (θ₀.2 ^ 2).toNNReal)
          (gaussianReal (ev n i θ.1) (θ.2 ^ 2).toNNReal) ≤ ENNReal.ofReal (n * ε ^ 2) ∧
       ∑ i, evariance (llr (gaussianReal (ev n i θ₀.1) (θ₀.2 ^ 2).toNNReal)
          (gaussianReal (ev n i θ.1) (θ.2 ^ 2).toNNReal))
          (gaussianReal (ev n i θ₀.1) (θ₀.2 ^ 2).toNNReal) ≤ ENNReal.ofReal (n * ε ^ 2)}

/-- The empirical norm `‖g‖ₙ = ((1/n) Σᵢ g(xᵢ)²)^{1/2}`. -/
noncomputable def empiricalNorm (n : ℕ) (g : E) : ℝ :=
  √((∑ i, ev n i g ^ 2) / n)

/-- The metric the theorem reports: `‖f − f′‖ₙ + |σ − σ′|`. -/
noncomputable def paramDist (n : ℕ) (θ θ' : E × ℝ) : ℝ :=
  empiricalNorm ev n (θ.1 - θ'.1) + |θ.2 - θ'.2|

/-- The empirical norm is a Euclidean norm, scaled by `1/√n`. -/
theorem empiricalNorm_eq (n : ℕ) (g : E) :
    empiricalNorm ev n g
      = ‖(WithLp.toLp 2 fun i => ev n i g : EuclideanSpace ℝ (Fin n))‖ / √n := by
  rw [empiricalNorm, EuclideanSpace.norm_eq, Real.sqrt_div' _ (Nat.cast_nonneg n)]
  simp [Real.norm_eq_abs, sq_abs]

theorem empiricalNorm_add_le (n : ℕ) (g h : E) :
    empiricalNorm ev n (g + h) ≤ empiricalNorm ev n g + empiricalNorm ev n h := by
  rw [empiricalNorm_eq, empiricalNorm_eq, empiricalNorm_eq, ← add_div]
  refine div_le_div_of_nonneg_right ?_ (Real.sqrt_nonneg _)
  have : (WithLp.toLp 2 fun i => ev n i (g + h) : EuclideanSpace ℝ (Fin n))
      = WithLp.toLp 2 (fun i => ev n i g) + WithLp.toLp 2 (fun i => ev n i h) := by
    rw [← WithLp.toLp_add]
    congr 1
    funext i
    simp
  rw [this]
  exact norm_add_le _ _

theorem empiricalNorm_neg (n : ℕ) (g : E) : empiricalNorm ev n (-g) = empiricalNorm ev n g := by
  simp [empiricalNorm]

theorem empiricalNorm_nonneg (n : ℕ) (g : E) : 0 ≤ empiricalNorm ev n g := Real.sqrt_nonneg _

/-- With `‖ev n i‖ ≤ 1`, the sup norm dominates the empirical norm. -/
theorem empiricalNorm_le_norm (hev : ∀ n i, ‖ev n i‖ ≤ 1) (n : ℕ) (g : E) :
    empiricalNorm ev n g ≤ ‖g‖ := by
  rw [empiricalNorm, Real.sqrt_le_left (norm_nonneg g)]
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  rw [div_le_iff₀ hn']
  have hi : ∀ i, ev n i g ^ 2 ≤ ‖g‖ ^ 2 := by
    intro i
    have h1 : |ev n i g| ≤ ‖g‖ := by
      have := (ev n i).le_opNorm g
      rw [← Real.norm_eq_abs]
      nlinarith [hev n i, norm_nonneg g, norm_nonneg (ev n i)]
    have := pow_le_pow_left₀ (abs_nonneg _) h1 2
    rwa [sq_abs] at this
  calc ∑ i, ev n i g ^ 2 ≤ ∑ _i : Fin n, ‖g‖ ^ 2 := Finset.sum_le_sum fun i _ => hi i
    _ = ‖g‖ ^ 2 * n := by simp [mul_comm]

theorem continuous_empiricalNorm (n : ℕ) : Continuous (empiricalNorm ev n) := by
  unfold empiricalNorm
  exact ((continuous_finsetSum _ fun i _ => (ev n i).continuous.pow 2).div_const _).sqrt

theorem paramDist_comm (n : ℕ) (θ θ' : E × ℝ) : paramDist ev n θ θ' = paramDist ev n θ' θ := by
  rw [paramDist, paramDist, ← neg_sub θ'.1, empiricalNorm_neg, abs_sub_comm]

theorem paramDist_triangle (n : ℕ) (θ θ' θ'' : E × ℝ) :
    paramDist ev n θ θ'' ≤ paramDist ev n θ θ' + paramDist ev n θ' θ'' := by
  unfold paramDist
  have h1 : θ.1 - θ''.1 = (θ.1 - θ'.1) + (θ'.1 - θ''.1) := by abel
  rw [h1]
  have h2 := empiricalNorm_add_le ev n (θ.1 - θ'.1) (θ'.1 - θ''.1)
  have h3 := abs_sub_le θ.2 θ'.2 θ''.2
  linarith

theorem continuous_paramDist (n : ℕ) (θ' : E × ℝ) : Continuous fun θ => paramDist ev n θ θ' := by
  unfold paramDist
  exact ((continuous_empiricalNorm ev n).comp (continuous_fst.sub continuous_const)).add
    (continuous_snd.sub continuous_const).abs

/-- `paramDist` is dominated by the sup norm and the noise gap. -/
theorem paramDist_le (hev : ∀ n i, ‖ev n i‖ ≤ 1) (n : ℕ) (θ θ' : E × ℝ) :
    paramDist ev n θ θ' ≤ ‖θ.1 - θ'.1‖ + |θ.2 - θ'.2| := by
  unfold paramDist
  linarith [empiricalNorm_le_norm ev hev n (θ.1 - θ'.1)]

variable [MeasurableSpace E]

/-- The posterior mass of `A` given the data `y`, by Bayes' formula. -/
noncomputable def posterior (prior : Measure (E × ℝ)) (n : ℕ) (y : Fin n → ℝ)
    (A : Set (E × ℝ)) : ℝ≥0∞ :=
  (∫⁻ θ in A, likelihood ev n θ y ∂prior) / ∫⁻ θ, likelihood ev n θ y ∂prior

end Model

variable [MeasurableSpace E] [BorelSpace E]

/-! ### The general theorem (cited) -/

/-- **AXIOM — the general contraction theorem (Ghosal & van der Vaart 2007, Theorem 4).**
For independent, non-identically distributed observations with a dominated model, a prior that
puts mass at least `exp(−C n εₙ²)` on (a measurable part of) the KL neighbourhood
`B*ₙ(θ₀, εₙ; 2)`, a sieve `Θₙ` whose complement has mass at most `exp(−(C + 4) n εₙ²)`, a cover
of the sieve by `exp(n ε̄ₙ²)` balls of radius `ξ ε̄ₙ/2` in a semimetric `d`, and exponentially
powerful local tests for `d` give a posterior that contracts at rate `M ε̄ₙ` in `d`, in mean under
the truth.

Cited, not formalized. The form taken here, as TWCD2025 and van der Vaart & van Zanten use it:
* **two rates**, `εₙ ≤ ε̄ₙ`: `εₙ` for the prior mass and the remaining mass, `ε̄ₙ` for the entropy
  and the conclusion;
* **a fixed constant `M`**, from the prior-mass and remaining-mass form of the conditions;
* **the testing hypothesis** in place of the Hellinger metric. Theorem 4 is the case where `d` is
  the root average squared Hellinger distance, whose tests exist by Le Cam and Birgé; UV-5 builds
  the tests for `‖f − f′‖ₙ + |σ − σ′|` directly.

The hypotheses are stated on the strong side, so the axiom stays a consequence of the cited
theorem: tests at every radius `r > 0`, a global cover of the sieve, measurable sieves and
neighbourhood parts, and `n εₙ² → ∞`. Whether the printed Theorem 4 states the fixed-`M` form or
the "every `Mₙ → ∞`" form is a citation detail to confirm at review (protolemma, UV-6). -/
axiom gvdv07_contraction
    (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (prior : Measure (E × ℝ)) [IsProbabilityMeasure prior]
    (hprior_pos : ∀ᵐ θ ∂prior, 0 < θ.2) (θ₀ : E × ℝ) (hθ₀ : 0 < θ₀.2)
    (d : ℕ → E × ℝ → E × ℝ → ℝ) (hd_comm : ∀ n θ θ', d n θ θ' = d n θ' θ)
    (hd_triangle : ∀ n θ θ' θ'', d n θ θ'' ≤ d n θ θ' + d n θ' θ'')
    (hd_meas : ∀ n θ', Measurable fun θ => d n θ θ')
    (ε εbar : ℕ → ℝ) (hε : ∀ᶠ n in atTop, 0 < ε n ∧ ε n ≤ εbar n)
    (hnε : Tendsto (fun n : ℕ => n * ε n ^ 2) atTop atTop)
    (Θ : ℕ → Set (E × ℝ)) (hΘ : ∀ n, MeasurableSet (Θ n)) (C K ξ : ℝ) (hK : 0 < K) (hξ : 0 < ξ)
    (htests : ∀ᶠ n in atTop, ∀ r > 0, ∀ θ₁ ∈ Θ n, r < d n θ₁ θ₀ →
      ∃ φ : (Fin n → ℝ) → ℝ, Measurable φ ∧ (∀ y, 0 ≤ φ y ∧ φ y ≤ 1) ∧
        ∫ y, φ y ∂dataLaw ev n θ₀ ≤ Real.exp (-(K * n * r ^ 2)) ∧
        ∀ θ ∈ Θ n, d n θ θ₁ ≤ ξ * r →
          ∫ y, (1 - φ y) ∂dataLaw ev n θ ≤ Real.exp (-(K * n * r ^ 2)))
    (hentropy : ∀ᶠ n in atTop, ∃ N : Finset (E × ℝ),
      (N.card : ℝ) ≤ Real.exp (n * εbar n ^ 2) ∧
        Θ n ⊆ ⋃ θ' ∈ N, {θ | d n θ θ' ≤ ξ * εbar n / 2})
    (hmass : ∀ᶠ n in atTop, ∃ B ⊆ klBall ev n θ₀ (ε n), MeasurableSet B ∧
      ENNReal.ofReal (Real.exp (-(C * n * ε n ^ 2))) ≤ prior B)
    (hrem : ∀ᶠ n in atTop,
      prior (Θ n)ᶜ ≤ ENNReal.ofReal (Real.exp (-((C + 4) * n * ε n ^ 2)))) :
    ∃ M : ℝ, Tendsto
      (fun n => ∫⁻ y, posterior ev prior n y {θ | M * εbar n ≤ d n θ θ₀} ∂dataLaw ev n θ₀)
      atTop (𝓝 0)

/-! ### UV-5 — the tests (the forward boundary) -/

/-- **UV-5 (OBLIGATION).** Exponentially powerful local tests for `‖f − f′‖ₙ + |σ − σ′|` on a
compact noise range: if `θ₁` is at distance more than `r` from the truth, one test separates the
truth from the whole `ξ r`-ball around `θ₁`, with both errors at most `exp(−K n r²)`. When
`‖f₁ − f*‖ₙ > r/2`, a likelihood-ratio test for the Gaussian means; when `|σ₁ − σ₀| > r/2`, a χ²
test of the residual scale around `f₁`. The constants depend on `σ_min` and `σ_max`. -/
theorem exists_tests (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) {σmin σmax : ℝ} (hσmin : 0 < σmin)
    (f₀ : E) {σ₀ : ℝ} (hσ₀ : σ₀ ∈ Icc σmin σmax) :
    ∃ K ξ : ℝ, 0 < K ∧ 0 < ξ ∧ ∀ n : ℕ, ∀ r > 0,
      ∀ θ₁ ∈ (univ ×ˢ Icc σmin σmax : Set (E × ℝ)), r < paramDist ev n θ₁ (f₀, σ₀) →
        ∃ φ : (Fin n → ℝ) → ℝ, Measurable φ ∧ (∀ y, 0 ≤ φ y ∧ φ y ≤ 1) ∧
          ∫ y, φ y ∂dataLaw ev n (f₀, σ₀) ≤ Real.exp (-(K * n * r ^ 2)) ∧
          ∀ θ ∈ (univ ×ˢ Icc σmin σmax : Set (E × ℝ)), paramDist ev n θ θ₁ ≤ ξ * r →
            ∫ y, (1 - φ y) ∂dataLaw ev n θ ≤ Real.exp (-(K * n * r ^ 2)) := by
  sorry -- [UV-5] the tests for ‖f − f′‖ₙ + |σ − σ′|: Gaussian LR tests and a χ² scale test

/-! ### Rates: the paper's radius shape `C n^(−a) (log n)^k` -/

/-- The radius tends to `0`: `(log n)^k = o(n^a)`. -/
theorem tendsto_priorMassRadius_zero (C a k : ℝ) (ha : 0 < a) :
    Tendsto (fun n : ℕ => priorMassRadius C a k n) atTop (𝓝 0) := by
  have h := (isLittleO_log_rpow_rpow_atTop k ha).tendsto_div_nhds_zero
  have h' := (h.comp tendsto_natCast_atTop_atTop).const_mul C
  rw [mul_zero] at h'
  refine h'.congr' (Filter.eventually_atTop.2 ⟨1, fun n hn => ?_⟩)
  have hx : (0 : ℝ) < n := by exact_mod_cast hn
  simp only [Function.comp, priorMassRadius, Real.rpow_neg hx.le]
  field_simp

/-- `n εₙ² → ∞`: eventually `n εₙ² ≥ log(1/εₙ)` (UV-2b), and `log(1/εₙ) → ∞`. -/
theorem tendsto_n_priorMassRadius_sq (C a k : ℝ) (hC : 0 < C) (ha : 0 < a) (ha' : a < 1 / 2)
    (hk : 0 ≤ k) :
    Tendsto (fun n : ℕ => n * priorMassRadius C a k n ^ 2) atTop atTop := by
  have h0 : Tendsto (fun n : ℕ => priorMassRadius C a k n) atTop (𝓝[>] 0) :=
    tendsto_nhdsWithin_iff.2 ⟨tendsto_priorMassRadius_zero C a k ha,
      eventually_priorMassRadius_pos C a k hC⟩
  have hlog : Tendsto (fun n : ℕ => Real.log (1 / priorMassRadius C a k n)) atTop atTop := by
    simp_rw [one_div]
    exact Real.tendsto_log_atTop.comp h0.inv_tendsto_nhdsGT_zero
  exact tendsto_atTop_mono' atTop (log_inv_radius_le_n_radius_sq C a k hC ha ha' hk) hlog

/-- `log n` is negligible against `n εₙ²`: eventually `log n ≤ δ n εₙ²` for every `δ > 0`. UV-2b
at `C = 1, k = 0` gives `a log n ≤ c n·n^(−2a)`, and `n εₙ² ≥ C² n·n^(−2a)` once `log n ≥ 1`. -/
theorem log_le_mul_n_priorMassRadius_sq (C a k δ : ℝ) (hC : 0 < C) (ha : 0 < a)
    (ha' : a < 1 / 2) (hk : 0 ≤ k) (hδ : 0 < δ) :
    ∀ᶠ n : ℕ in atTop, Real.log n ≤ δ * (n * priorMassRadius C a k n ^ 2) := by
  filter_upwards [log_inv_mul_radius_le 1 a 0 1 (δ * a * C ^ 2) one_pos ha ha' le_rfl one_pos
      (by positivity), Filter.eventually_ge_atTop ⌈Real.exp 1⌉₊] with n h hn
  have hxe : Real.exp 1 ≤ (n : ℝ) := (Nat.le_ceil _).trans (by exact_mod_cast hn)
  have hx0 : (0 : ℝ) < n := lt_of_lt_of_le (Real.exp_pos 1) hxe
  have hL1 : 1 ≤ Real.log n := by
    rw [← Real.log_exp 1]
    exact Real.log_le_log (Real.exp_pos 1) hxe
  have hρ : priorMassRadius 1 a 0 n = (n : ℝ) ^ (-a) := by
    simp [priorMassRadius]
  rw [hρ, one_mul, one_div, ← Real.rpow_neg hx0.le, neg_neg, Real.log_rpow hx0] at h
  have hLk : 1 ≤ Real.log n ^ k := Real.one_le_rpow hL1 hk
  have hsq : C ^ 2 * (n * ((n : ℝ) ^ (-a)) ^ 2) ≤ n * priorMassRadius C a k n ^ 2 := by
    have hnn : 0 ≤ (n : ℝ) * ((n : ℝ) ^ (-a)) ^ 2 := by positivity
    have : priorMassRadius C a k n ^ 2 = C ^ 2 * ((n : ℝ) ^ (-a)) ^ 2 * (Real.log n ^ k) ^ 2 := by
      unfold priorMassRadius
      ring
    rw [this]
    have h1 : 1 ≤ (Real.log n ^ k) ^ 2 := by nlinarith
    calc C ^ 2 * (n * ((n : ℝ) ^ (-a)) ^ 2) = C ^ 2 * (n * ((n : ℝ) ^ (-a)) ^ 2) * 1 := by ring
      _ ≤ C ^ 2 * (n * ((n : ℝ) ^ (-a)) ^ 2) * (Real.log n ^ k) ^ 2 :=
        mul_le_mul_of_nonneg_left h1 (by positivity)
      _ = n * (C ^ 2 * ((n : ℝ) ^ (-a)) ^ 2 * (Real.log n ^ k) ^ 2) := by ring
  have hmain : a * Real.log n ≤ a * (δ * (n * priorMassRadius C a k n ^ 2)) := by
    calc a * Real.log n ≤ δ * a * C ^ 2 * n * ((n : ℝ) ^ (-a)) ^ 2 := h
      _ = (δ * a) * (C ^ 2 * (n * ((n : ℝ) ^ (-a)) ^ 2)) := by ring
      _ ≤ (δ * a) * (n * priorMassRadius C a k n ^ 2) :=
        mul_le_mul_of_nonneg_left hsq (by positivity)
      _ = a * (δ * (n * priorMassRadius C a k n ^ 2)) := by ring
  exact le_of_mul_le_mul_left hmain ha

/-- Two radii of the same polynomial order: `C₀ n^(−a)(log n)^k ≤ (C₀/Cb) · Cb n^(−a)(log n)^kb`
once `log n ≥ 1`, for `k ≤ kb`. -/
theorem priorMassRadius_le_mul (C₀ Cb a k kb : ℝ) (hC₀ : 0 < C₀) (hCb : 0 < Cb) (hkk : k ≤ kb) :
    ∀ᶠ n : ℕ in atTop, priorMassRadius C₀ a k n ≤ C₀ / Cb * priorMassRadius Cb a kb n := by
  filter_upwards [Filter.eventually_ge_atTop ⌈Real.exp 1⌉₊] with n hn
  have hxe : Real.exp 1 ≤ (n : ℝ) := (Nat.le_ceil _).trans (by exact_mod_cast hn)
  have hx0 : (0 : ℝ) < n := lt_of_lt_of_le (Real.exp_pos 1) hxe
  have hL1 : 1 ≤ Real.log n := by
    rw [← Real.log_exp 1]
    exact Real.log_le_log (Real.exp_pos 1) hxe
  have hpow : Real.log n ^ k ≤ Real.log n ^ kb := Real.rpow_le_rpow_of_exponent_le hL1 hkk
  have hxa : 0 < (n : ℝ) ^ (-a) := Real.rpow_pos_of_pos hx0 _
  have e : C₀ / Cb * (Cb * (n : ℝ) ^ (-a) * Real.log n ^ kb)
      = C₀ * (n : ℝ) ^ (-a) * Real.log n ^ kb := by
    field_simp
  unfold priorMassRadius
  rw [e]
  exact mul_le_mul_of_nonneg_left hpow (by positivity)

/-! ### [P1]: the amplitude and noise priors -/

/-- The map `(s, (w, η)) ↦ (s • w, s √η)`: amplitude times unit-amplitude GP draw, and the noise
standard deviation `σ = s √η`. -/
noncomputable def amplitudeNoiseMap (p : ℝ × (E × ℝ)) : E × ℝ := (p.1 • p.2.1, p.1 * √p.2.2)

theorem measurable_amplitudeNoiseMap : Measurable (amplitudeNoiseMap (E := E)) := by
  apply Continuous.measurable
  unfold amplitudeNoiseMap
  exact (continuous_fst.smul (continuous_fst.comp continuous_snd)).prodMk
    (continuous_fst.mul (Real.continuous_sqrt.comp (continuous_snd.comp continuous_snd)))

/-- **The proper prior [P1].** `f = s • W` and `σ = s √η`, with `s ~ πs`, `W ~ μW` (the GP law of
`Wᵗ`, the bandwidth integrated out) and `η ~ πη` independent. -/
noncomputable def amplitudeNoisePrior (πs : Measure ℝ) (μW : Measure E) (πη : Measure ℝ) :
    Measure (E × ℝ) :=
  (πs.prod (μW.prod πη)).map amplitudeNoiseMap

theorem isProbabilityMeasure_amplitudeNoisePrior (πs : Measure ℝ) (μW : Measure E)
    (πη : Measure ℝ) [IsProbabilityMeasure πs] [IsProbabilityMeasure μW]
    [IsProbabilityMeasure πη] : IsProbabilityMeasure (amplitudeNoisePrior πs μW πη) :=
  Measure.isProbabilityMeasure_map measurable_amplitudeNoiseMap.aemeasurable

/-- **UV-2, step 1: the amplitude window.** If `σ₀ ∈ (s_min √η_lo, s_max √η_hi)`, there is a
nondegenerate window `[s₁, s₂] ⊆ [s_min, s_max]` on which `σ₀/s` stays in `[√η_lo, √η_hi)`, with
room above: `s₂ √η_lo ≤ σ₀ < s₁ √η_hi`. -/
theorem exists_amplitude_window {smin smax ηlo ηhi σ₀ : ℝ} (hsmin : 0 < smin) (hs : smin < smax)
    (hηlo : 0 < ηlo) (hη : ηlo < ηhi) (hσ₀ : smin * √ηlo < σ₀ ∧ σ₀ < smax * √ηhi) :
    ∃ s₁ s₂ : ℝ, smin ≤ s₁ ∧ s₁ < s₂ ∧ s₂ ≤ smax ∧ s₂ * √ηlo ≤ σ₀ ∧ σ₀ < s₁ * √ηhi := by
  have hlo : 0 < √ηlo := Real.sqrt_pos.2 hηlo
  have hhi : √ηlo < √ηhi := Real.sqrt_lt_sqrt hηlo.le hη
  have hσpos : 0 < σ₀ := lt_of_le_of_lt (by positivity) hσ₀.1
  set ℓ := max smin (σ₀ / √ηhi) with hℓ
  set u := min smax (σ₀ / √ηlo) with hu
  have hℓu : ℓ < u := by
    refine max_lt (lt_min hs ?_) (lt_min ?_ (div_lt_div_of_pos_left hσpos hlo hhi))
    · rw [lt_div_iff₀ hlo]
      exact hσ₀.1
    · rw [div_lt_iff₀ (hlo.trans hhi)]
      exact hσ₀.2
  have h1 : smin ≤ ℓ := le_max_left _ _
  have h2 : σ₀ / √ηhi ≤ ℓ := le_max_right _ _
  have h3 : u ≤ smax := min_le_left _ _
  have h4 : u ≤ σ₀ / √ηlo := min_le_right _ _
  refine ⟨(2 * ℓ + u) / 3, (ℓ + 2 * u) / 3, by linarith, by linarith, by linarith, ?_, ?_⟩
  · rw [← le_div_iff₀ hlo]
    linarith
  · rw [← div_lt_iff₀ (hlo.trans hhi)]
    linarith

/-- **UV-2, step 3: the noise slab.** Under a log-uniform `η` on `[η_lo, η_hi]` (here: at least the
log-uniform mass on every subinterval), `σ = s √η` lands within `ρ` of `σ₀` with probability at
least `ρ/(σ₀ log(η_hi/η_lo))`, as long as `σ₀/s ≥ √η_lo`, `ρ ≤ σ₀` and `σ₀ + ρ ≤ s √η_hi`: the
η-interval `[(σ₀/s)², ((σ₀+ρ)/s)²]` has log-length `2 log(1 + ρ/σ₀) ≥ ρ/σ₀`. -/
theorem noise_slab_mass (πη : Measure ℝ) {ηlo ηhi : ℝ} (hηlo : 0 < ηlo) (hη : ηlo < ηhi)
    (hπη : ∀ a b, ηlo ≤ a → a ≤ b → b ≤ ηhi →
      ENNReal.ofReal (Real.log (b / a) / Real.log (ηhi / ηlo)) ≤ πη (Icc a b))
    {s σ₀ ρ : ℝ} (hs : 0 < s) (hσ₀ : 0 < σ₀) (hρ : 0 < ρ) (hρσ : ρ ≤ σ₀)
    (hlo : s * √ηlo ≤ σ₀) (hhi : σ₀ + ρ ≤ s * √ηhi) :
    ENNReal.ofReal (ρ / (σ₀ * Real.log (ηhi / ηlo))) ≤ πη {η | |s * √η - σ₀| ≤ ρ} := by
  have hL : 0 < Real.log (ηhi / ηlo) := Real.log_pos ((one_lt_div hηlo).2 hη)
  have hq0 : 0 ≤ σ₀ / s := by positivity
  have hq1 : 0 ≤ (σ₀ + ρ) / s := by positivity
  have ha' : ηlo ≤ (σ₀ / s) ^ 2 := by
    have h : √ηlo ≤ σ₀ / s := by
      rw [le_div_iff₀ hs, mul_comm]
      exact hlo
    have := pow_le_pow_left₀ (Real.sqrt_nonneg _) h 2
    rwa [Real.sq_sqrt hηlo.le] at this
  have hab : (σ₀ / s) ^ 2 ≤ ((σ₀ + ρ) / s) ^ 2 :=
    pow_le_pow_left₀ hq0 (div_le_div_of_nonneg_right (by linarith) hs.le) 2
  have hb' : ((σ₀ + ρ) / s) ^ 2 ≤ ηhi := by
    have h : (σ₀ + ρ) / s ≤ √ηhi := by
      rw [div_le_iff₀ hs, mul_comm]
      exact hhi
    have := pow_le_pow_left₀ hq1 h 2
    rwa [Real.sq_sqrt (hηlo.trans hη).le] at this
  have hsub : Icc ((σ₀ / s) ^ 2) (((σ₀ + ρ) / s) ^ 2) ⊆ {η | |s * √η - σ₀| ≤ ρ} := by
    rintro η ⟨h1, h2⟩
    have e1 : σ₀ / s ≤ √η := by
      have := Real.sqrt_le_sqrt h1
      rwa [Real.sqrt_sq hq0] at this
    have e2 : √η ≤ (σ₀ + ρ) / s := by
      have := Real.sqrt_le_sqrt h2
      rwa [Real.sqrt_sq hq1] at this
    rw [div_le_iff₀ hs] at e1
    rw [le_div_iff₀ hs] at e2
    change |s * √η - σ₀| ≤ ρ
    rw [abs_le]
    constructor <;> nlinarith
  have hlog : ρ / σ₀ ≤ Real.log (((σ₀ + ρ) / s) ^ 2 / (σ₀ / s) ^ 2) := by
    have e : ((σ₀ + ρ) / s) ^ 2 / (σ₀ / s) ^ 2 = ((σ₀ + ρ) / σ₀) ^ 2 := by
      field_simp
    have hpos : 0 < (σ₀ + ρ) / σ₀ := by positivity
    have h1 := Real.one_sub_inv_le_log_of_pos hpos
    have h2 : 1 - ((σ₀ + ρ) / σ₀)⁻¹ = ρ / (σ₀ + ρ) := by
      field_simp
      ring
    have h3 : ρ / σ₀ ≤ 2 * (ρ / (σ₀ + ρ)) := by
      rw [div_le_iff₀ hσ₀]
      have h4 : 2 * (ρ / (σ₀ + ρ)) * σ₀ = ρ * (2 * σ₀ / (σ₀ + ρ)) := by ring
      rw [h4]
      have h5 : 1 ≤ 2 * σ₀ / (σ₀ + ρ) := by
        rw [le_div_iff₀ (by linarith)]
        linarith
      nlinarith
    rw [e, Real.log_pow]
    push_cast
    linarith
  calc ENNReal.ofReal (ρ / (σ₀ * Real.log (ηhi / ηlo)))
      ≤ ENNReal.ofReal (Real.log (((σ₀ + ρ) / s) ^ 2 / (σ₀ / s) ^ 2) / Real.log (ηhi / ηlo)) := by
        apply ENNReal.ofReal_le_ofReal
        rw [← div_div]
        exact div_le_div_of_nonneg_right hlog hL.le
    _ ≤ πη (Icc ((σ₀ / s) ^ 2) (((σ₀ + ρ) / s) ^ 2)) := hπη _ _ ha' hab hb'
    _ ≤ πη {η | |s * √η - σ₀| ≤ ρ} := measure_mono hsub

/-- **UV-2, steps 2 and 4 (one `n`).** The prior mass of the slab `‖f − f*‖∞ ≤ ρ`, `|σ − σ₀| ≤ ρ`
is at least `πs([s₁, s₂]) · m · m′`, given lower bounds `m` on the GP mass of the rescaled ball
`closedBall (s⁻¹ • f*) (ρ/s)` (the authors' Remark A.1, through UV-2a) and `m′` on the noise slab,
uniformly over the window. By independence the slab's section at `s` is a rectangle. -/
theorem prior_slab_mass (πs : Measure ℝ) (μW : Measure E) (πη : Measure ℝ) [SFinite πs]
    [SFinite μW] [SFinite πη] {f₀ : E} {σ₀ ρ s₁ s₂ m m' : ℝ} (hs₁ : 0 < s₁)
    (hW : ∀ s ∈ Icc s₁ s₂, ENNReal.ofReal m ≤ μW (closedBall (s⁻¹ • f₀) (ρ / s)))
    (hη : ∀ s ∈ Icc s₁ s₂, ENNReal.ofReal m' ≤ πη {η | |s * √η - σ₀| ≤ ρ}) :
    πs (Icc s₁ s₂) * (ENNReal.ofReal m * ENNReal.ofReal m')
      ≤ amplitudeNoisePrior πs μW πη (closedBall f₀ ρ ×ˢ closedBall σ₀ ρ) := by
  have hslab : MeasurableSet (closedBall f₀ ρ ×ˢ closedBall σ₀ ρ) :=
    measurableSet_closedBall.prod measurableSet_closedBall
  rw [amplitudeNoisePrior, Measure.map_apply measurable_amplitudeNoiseMap hslab,
    Measure.prod_apply (measurable_amplitudeNoiseMap hslab), mul_comm,
    ← lintegral_indicator_const measurableSet_Icc]
  refine lintegral_mono fun s => ?_
  by_cases hs : s ∈ Icc s₁ s₂
  · rw [indicator_of_mem hs]
    have hs0 : 0 < s := hs₁.trans_le hs.1
    have hsec : Prod.mk s ⁻¹' (amplitudeNoiseMap ⁻¹' (closedBall f₀ ρ ×ˢ closedBall σ₀ ρ))
        = ((fun w => s • w) ⁻¹' closedBall f₀ ρ) ×ˢ {η | |s * √η - σ₀| ≤ ρ} := by
      ext ⟨w, η⟩
      simp [amplitudeNoiseMap, Real.dist_eq]
    have hball : μW ((fun w => s • w) ⁻¹' closedBall f₀ ρ)
        = μW (closedBall (s⁻¹ • f₀) (ρ / s)) := by
      rw [← map_smul_closedBall μW hs0 f₀ ρ,
        Measure.map_apply (continuous_const_smul s).measurable measurableSet_closedBall]
    rw [hsec, Measure.prod_prod, hball]
    exact mul_le_mul' (hW s hs) (hη s hs)
  · rw [indicator_of_notMem hs]
    exact zero_le

/-! ### UV-1 at the parameter level -/

omit [MeasurableSpace E] [BorelSpace E] in
/-- **UV-1 for `n` observations.** The slab `‖f − f*‖∞ ≤ σ₀ε/6`, `|σ − σ₀| ≤ σ₀ε/6`, with
`0 ≤ ε ≤ 3`, lies in the KL neighbourhood `B*ₙ((f*, σ₀), ε; 2)`: each design point sees
`|f(xᵢ) − f*(xᵢ)| ≤ ‖f − f*‖∞`, so `klDiv_le_and_evariance_le_of_slab` applies termwise. -/
theorem slab_subset_klBall (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (hev : ∀ n i, ‖ev n i‖ ≤ 1)
    (n : ℕ) {f₀ : E} {σ₀ ε : ℝ} (h0 : 0 < σ₀) (hε0 : 0 ≤ ε) (hε3 : ε ≤ 3) :
    closedBall f₀ (σ₀ * ε / 6) ×ˢ closedBall σ₀ (σ₀ * ε / 6) ⊆ klBall ev n (f₀, σ₀) ε := by
  rintro ⟨f, σ⟩ ⟨hf, hσ⟩
  rw [mem_closedBall, dist_eq_norm] at hf
  rw [mem_closedBall, Real.dist_eq] at hσ
  have hpt : ∀ i, |ev n i f - ev n i f₀| ≤ σ₀ * ε / 6 := by
    intro i
    rw [← map_sub, ← Real.norm_eq_abs]
    calc ‖ev n i (f - f₀)‖ ≤ ‖ev n i‖ * ‖f - f₀‖ := (ev n i).le_opNorm _
      _ ≤ 1 * ‖f - f₀‖ := mul_le_mul_of_nonneg_right (hev n i) (norm_nonneg _)
      _ ≤ σ₀ * ε / 6 := by rw [one_mul]; exact hf
  have hsum : ∑ _i : Fin n, ENNReal.ofReal (ε ^ 2) = ENNReal.ofReal (n * ε ^ 2) := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
      ENNReal.ofReal_mul (Nat.cast_nonneg _), ENNReal.ofReal_natCast]
  refine ⟨?_, ?_⟩
  · rw [← hsum]
    exact Finset.sum_le_sum fun i _ =>
      (klDiv_le_and_evariance_le_of_slab h0 hε0 hε3 (hpt i) hσ).1
  · rw [← hsum]
    exact Finset.sum_le_sum fun i _ =>
      (klDiv_le_and_evariance_le_of_slab h0 hε0 hε3 (hpt i) hσ).2

/-! ### UV-3 — the sieve and its remaining mass -/

/-- The sieve `B′ₙ = {s • g | s ∈ [s_min, s_max], g ∈ Bₙ} × [σ_min, σ_max]`. -/
def sieve (B : Set E) (smin smax σmin σmax : ℝ) : Set (E × ℝ) :=
  image2 (· • ·) (Icc smin smax) B ×ˢ Icc σmin σmax

/-- **UV-3 (60′).** Outside the (closed) sieve, the prior has no more mass than the GP law has
outside `Bₙ`: `f = s • W` lies in `s • Bₙ` whenever `W ∈ Bₙ`, and `s` and `σ = s √η` stay in their
compact ranges almost surely. -/
theorem prior_compl_sieve_le (πs : Measure ℝ) (μW : Measure E) (πη : Measure ℝ)
    [IsProbabilityMeasure πs] [IsProbabilityMeasure μW] [IsProbabilityMeasure πη]
    {smin smax ηlo ηhi : ℝ} (hsmin : 0 ≤ smin)
    (hπs : πs (Icc smin smax)ᶜ = 0) (hπη : πη (Icc ηlo ηhi)ᶜ = 0) (B : Set E) :
    amplitudeNoisePrior πs μW πη (closure (sieve B smin smax (smin * √ηlo) (smax * √ηhi)))ᶜ
      ≤ μW Bᶜ := by
  rw [amplitudeNoisePrior, Measure.map_apply measurable_amplitudeNoiseMap
    isClosed_closure.isOpen_compl.measurableSet]
  have hsub : amplitudeNoiseMap ⁻¹' (closure (sieve B smin smax (smin * √ηlo) (smax * √ηhi)))ᶜ
      ⊆ ((Icc smin smax)ᶜ ×ˢ univ ∪ univ ×ˢ (Bᶜ ×ˢ univ))
        ∪ univ ×ˢ (univ ×ˢ (Icc ηlo ηhi)ᶜ) := by
    rintro ⟨s, w, η⟩ hp
    simp only [mem_preimage, mem_compl_iff] at hp
    by_contra hcon
    simp only [mem_union, mem_prod, mem_compl_iff, mem_univ, and_true, true_and, not_or,
      not_not] at hcon
    obtain ⟨⟨hs, hw⟩, hη⟩ := hcon
    apply hp
    refine subset_closure ⟨mem_image2_of_mem hs hw, ?_, ?_⟩
    · exact mul_le_mul hs.1 (Real.sqrt_le_sqrt hη.1) (Real.sqrt_nonneg _) (hsmin.trans hs.1)
    · exact mul_le_mul hs.2 (Real.sqrt_le_sqrt hη.2) (Real.sqrt_nonneg _)
        (hsmin.trans (hs.1.trans hs.2))
  calc (πs.prod (μW.prod πη)) (amplitudeNoiseMap ⁻¹'
        (closure (sieve B smin smax (smin * √ηlo) (smax * √ηhi)))ᶜ)
      ≤ (πs.prod (μW.prod πη)) (((Icc smin smax)ᶜ ×ˢ univ ∪ univ ×ˢ (Bᶜ ×ˢ univ))
          ∪ univ ×ˢ (univ ×ˢ (Icc ηlo ηhi)ᶜ)) :=
        measure_mono hsub
    _ ≤ (πs.prod (μW.prod πη)) ((Icc smin smax)ᶜ ×ˢ univ)
          + (πs.prod (μW.prod πη)) (univ ×ˢ (Bᶜ ×ˢ univ))
          + (πs.prod (μW.prod πη)) (univ ×ˢ (univ ×ˢ (Icc ηlo ηhi)ᶜ)) :=
        (measure_union_le _ _).trans (add_le_add_left (measure_union_le _ _) _)
    _ = μW Bᶜ := by
        simp only [Measure.prod_prod, hπs, hπη, measure_univ, mul_zero, one_mul,
          mul_one, zero_add, add_zero]

/-! ### UV-4 — covering the sieve -/

omit [MeasurableSpace E] [BorelSpace E] in
open scoped Classical in
/-- **UV-4 (61′), the cover.** Covers `N₁ ⊆ B` at radius `ε₁` (so `‖g′‖ ≤ R`), `N₂` of
`[s_min, s_max]` at `ε₂` and `N₃` of `[σ_min, σ_max]` at `ε₃` give a cover of the closed sieve by
`paramDist`-balls of radius `r ≥ s_max ε₁ + R ε₂ + ε₃`, centred at `(s′ • g′, σ′)`. UV-4a covers the
scaled family; the empirical norm is dominated by the sup norm. -/
theorem closure_sieve_subset_cover (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ)
    (hev : ∀ n i, ‖ev n i‖ ≤ 1) (n : ℕ) {B : Set E} {smin smax σmin σmax R r : ℝ}
    (hsmin : 0 ≤ smin) (hsmax : 0 ≤ smax) (hR0 : 0 ≤ R) {ε₁ ε₂ ε₃ : ℝ≥0}
    (N₁ : Finset E) (hN₁B : (N₁ : Set E) ⊆ B) (hR : ∀ g ∈ B, ‖g‖ ≤ R)
    (hN₁ : IsCover ε₁ B N₁) (N₂ : Finset ℝ) (hN₂ : IsCover ε₂ (Icc smin smax) N₂)
    (N₃ : Finset ℝ) (hN₃ : IsCover ε₃ (Icc σmin σmax) N₃) (hr : smax * ε₁ + R * ε₂ + ε₃ ≤ r) :
    closure (sieve B smin smax σmin σmax)
      ⊆ ⋃ θ' ∈ ((N₂ ×ˢ N₁) ×ˢ N₃).image (fun q => (q.1.1 • q.1.2, q.2)),
          {θ | paramDist ev n θ θ' ≤ r} := by
  have hclosed : IsClosed (⋃ θ' ∈ ((N₂ ×ˢ N₁) ×ˢ N₃).image (fun q => (q.1.1 • q.1.2, q.2)),
      {θ | paramDist ev n θ θ' ≤ r}) :=
    isClosed_biUnion_finset fun θ' _ => isClosed_le (continuous_paramDist ev n θ') continuous_const
  refine closure_minimal ?_ hclosed
  have hcov := isCover_image2_smul hsmin (fun g hg => hR g (hN₁B hg))
    (Real.le_coe_toNNReal (smax * ε₁ + R * ε₂)) hN₁ hN₂
  rw [isCover_iff_subset_iUnion_closedBall] at hcov hN₃
  rintro ⟨x, σ⟩ ⟨hx, hσ⟩
  obtain ⟨x', hx'N, hxx'⟩ := mem_iUnion₂.1 (hcov hx)
  obtain ⟨σ', hσ'N, hσσ'⟩ := mem_iUnion₂.1 (hN₃ hσ)
  obtain ⟨s', hs', g', hg', rfl⟩ := hx'N
  rw [mem_closedBall, dist_eq_norm] at hxx'
  rw [mem_closedBall, Real.dist_eq] at hσσ'
  refine mem_iUnion₂.2 ⟨(s' • g', σ'), ?_, ?_⟩
  · exact Finset.mem_image.2 ⟨((s', g'), σ'), Finset.mem_product.2
      ⟨Finset.mem_product.2 ⟨Finset.mem_coe.1 hs', Finset.mem_coe.1 hg'⟩, Finset.mem_coe.1 hσ'N⟩,
      rfl⟩
  · change paramDist ev n (x, σ) (s' • g', σ') ≤ r
    have h12 : (((smax * ε₁ + R * ε₂).toNNReal : ℝ≥0) : ℝ) = smax * ε₁ + R * ε₂ :=
      Real.coe_toNNReal _ (by positivity)
    rw [h12] at hxx'
    calc paramDist ev n (x, σ) (s' • g', σ') ≤ ‖x - s' • g'‖ + |σ - σ'| :=
          paramDist_le ev hev n _ _
      _ ≤ (smax * ε₁ + R * ε₂) + ε₃ := add_le_add hxx' hσσ'
      _ ≤ r := hr

/-- The arithmetic of the count in `eventually_sieve_cover`: three slack factors, each at most
`exp(X/3)`, and `L² ≥ K₀ + 1`. -/
theorem card_le_exp_of_slack {N₁ N₂ N₃ A₂ A₃ K₀ L t ε X : ℝ} (hε : 0 < ε) (ht : 0 < t)
    (hX : 0 ≤ X) (hN₁0 : 0 ≤ N₁) (hN₃0 : 0 ≤ N₃) (hN₁ : N₁ ≤ Real.exp (K₀ * X))
    (hN₂ : N₂ ≤ (|A₂| + 2) * (t / ε)) (hN₃ : N₃ ≤ (|A₃| + 2) * (1 / ε))
    (hA : (|A₂| + 2) * (|A₃| + 2) ≤ Real.exp (X / 3)) (hP : t ≤ Real.exp (X / 3))
    (hE : (1 / ε) ^ 2 ≤ Real.exp (X / 3)) (hL : K₀ + 1 ≤ L ^ 2) :
    N₂ * N₁ * N₃ ≤ Real.exp (L ^ 2 * X) := by
  have h2 : 0 ≤ (|A₂| + 2) * (t / ε) := by positivity
  calc N₂ * N₁ * N₃
      ≤ ((|A₂| + 2) * (t / ε)) * Real.exp (K₀ * X) * ((|A₃| + 2) * (1 / ε)) := by
        gcongr
    _ = ((|A₂| + 2) * (|A₃| + 2)) * t * (1 / ε) ^ 2 * Real.exp (K₀ * X) := by ring
    _ ≤ Real.exp (X / 3) * Real.exp (X / 3) * Real.exp (X / 3) * Real.exp (K₀ * X) := by
        gcongr
    _ = Real.exp ((K₀ + 1) * X) := by
        rw [← Real.exp_add, ← Real.exp_add, ← Real.exp_add]
        congr 1
        ring
    _ ≤ Real.exp (L ^ 2 * X) := Real.exp_le_exp.2 (mul_le_mul_of_nonneg_right hL hX)

omit [MeasurableSpace E] [BorelSpace E] in
open scoped Classical in
/-- **UV-4 (61′), the count.** With TWCD2025's (61) for `Bₙ` at radius `(ξ/(12 s_max)) ε̄ₙ` and a
sup-norm radius `Rₙ ≤ n^p`, the closed sieve is covered by at most `exp(n (L ε̄ₙ)²)` balls of
`paramDist`-radius `ξ L ε̄ₙ/2`, once `L ≥ 1 + |K₀|`. The interval counts are `O(n^p/ε̄ₙ)` and
`O(1/ε̄ₙ)`, and `log n`, `log(1/ε̄ₙ)` are `o(n ε̄ₙ²)`. -/
theorem eventually_sieve_cover (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (hev : ∀ n i, ‖ev n i‖ ≤ 1)
    (B : ℕ → Set E) {smin smax σmin σmax : ℝ} (hsmin : 0 ≤ smin) (hsmax : 0 < smax)
    {Cb a kb ξ L K₀ : ℝ} (hCb : 0 < Cb) (ha : 0 < a) (ha' : a < 1 / 2) (hkb : 0 ≤ kb)
    (hξ : 0 < ξ) (hL : 1 + |K₀| ≤ L)
    (h61 : ∀ᶠ n : ℕ in atTop, ∃ N : ℕ,
      externalCoveringNumber (ξ / (12 * smax) * priorMassRadius Cb a kb n).toNNReal (B n) ≤ N ∧
        Real.log N ≤ K₀ * n * priorMassRadius Cb a kb n ^ 2)
    {p : ℕ} (hR : ∀ᶠ n : ℕ in atTop, ∀ g ∈ B n, ‖g‖ ≤ (n : ℝ) ^ p) :
    ∀ᶠ n : ℕ in atTop, ∃ N : Finset (E × ℝ),
      (N.card : ℝ) ≤ Real.exp (n * (L * priorMassRadius Cb a kb n) ^ 2) ∧
        closure (sieve (B n) smin smax σmin σmax)
          ⊆ ⋃ θ' ∈ N, {θ | paramDist ev n θ θ' ≤ ξ * (L * priorMassRadius Cb a kb n) / 2} := by
  have hL1 : 1 ≤ L := by linarith [abs_nonneg K₀]
  have hLpos : 0 < L := by linarith
  have hL2 : K₀ + 1 ≤ L ^ 2 := by
    have : L ≤ L ^ 2 := by nlinarith
    linarith [le_abs_self K₀]
  have hceil : ∀ x : ℝ, (⌈x⌉₊ : ℝ) ≤ |x| + 1 := fun x => by
    rcases le_or_gt x 0 with hx | hx
    · rw [Nat.ceil_eq_zero.2 hx, Nat.cast_zero]
      positivity
    · have := Nat.ceil_lt_add_one hx.le
      rw [abs_of_pos hx]
      linarith
  have hX := tendsto_n_priorMassRadius_sq Cb a kb hCb ha ha' hkb
  filter_upwards [h61, hR, eventually_priorMassRadius_pos Cb a kb hCb,
    (tendsto_priorMassRadius_zero Cb a kb ha).eventually (gt_mem_nhds one_pos),
    Filter.eventually_ge_atTop 1,
    hX.eventually_ge_atTop (3 * Real.log ((|3 * (smax - smin) / (ξ * L)| + 2)
      * (|3 * (σmax - σmin) / (ξ * L)| + 2))),
    log_le_mul_n_priorMassRadius_sq Cb a kb (1 / (3 * ((p : ℝ) + 1))) hCb ha ha' hkb
      (by positivity),
    log_inv_mul_radius_le Cb a kb 1 (1 / 6) hCb ha ha' hkb one_pos (by norm_num)]
    with n h61n hRn hε hε1 hn1 hXc hlogn hlogε
  clear h61 hR hX
  generalize priorMassRadius Cb a kb n = ε at h61n hε hε1 hXc hlogn hlogε ⊢
  have hn : (1 : ℝ) ≤ n := by exact_mod_cast hn1
  have hnp : (1 : ℝ) ≤ (n : ℝ) ^ p := one_le_pow₀ hn
  have hnp0 : (0 : ℝ) < (n : ℝ) ^ p := by positivity
  have hX0 : 0 ≤ (n : ℝ) * ε ^ 2 := by positivity
  -- the cover of `Bₙ`, internal, from (61)
  obtain ⟨N₀, hN₀, hlogN₀⟩ := h61n
  have hr₁c : (((ξ / (12 * smax) * ε).toNNReal : ℝ≥0) : ℝ) = ξ / (12 * smax) * ε :=
    Real.coe_toNNReal _ (by positivity)
  have hint : coveringNumber (2 * (ξ / (12 * smax) * ε).toNNReal) (B n) ≤ N₀ :=
    (coveringNumber_two_mul_le_externalCoveringNumber _ (B n)).trans hN₀
  obtain ⟨S₁, hS₁B, hS₁fin, hS₁, hS₁card⟩ :=
    exists_set_encard_eq_coveringNumber (ne_top_of_le_ne_top (ENat.coe_ne_top N₀) hint)
  have hN₁coe : (hS₁fin.toFinset : Set E) = S₁ := hS₁fin.coe_toFinset
  have hN₁card : (hS₁fin.toFinset.card : ℝ) ≤ N₀ := by
    have h : (hS₁fin.toFinset.card : ℕ∞) ≤ N₀ := by
      rw [← hS₁fin.encard_eq_coe_toFinset_card, hS₁card]
      exact hint
    exact_mod_cast h
  -- the interval covers
  have hr₂c : (((ξ * L * ε / (6 * (n : ℝ) ^ p)).toNNReal : ℝ≥0) : ℝ)
      = ξ * L * ε / (6 * (n : ℝ) ^ p) := Real.coe_toNNReal _ (by positivity)
  have hr₂0 : 0 < (ξ * L * ε / (6 * (n : ℝ) ^ p)).toNNReal := by
    rw [← NNReal.coe_pos, hr₂c]
    positivity
  have hr₃c : (((ξ * L * ε / 6).toNNReal : ℝ≥0) : ℝ) = ξ * L * ε / 6 :=
    Real.coe_toNNReal _ (by positivity)
  have hr₃0 : 0 < (ξ * L * ε / 6).toNNReal := by
    rw [← NNReal.coe_pos, hr₃c]
    positivity
  obtain ⟨N₂, -, hN₂, hN₂card⟩ := exists_isCover_Icc smin smax hr₂0
  obtain ⟨N₃, -, hN₃, hN₃card⟩ := exists_isCover_Icc σmin σmax hr₃0
  refine ⟨((N₂ ×ˢ hS₁fin.toFinset) ×ˢ N₃).image (fun q => (q.1.1 • q.1.2, q.2)), ?_, ?_⟩
  · -- the count
    have hcard : ((((N₂ ×ˢ hS₁fin.toFinset) ×ˢ N₃).image
        (fun q => (q.1.1 • q.1.2, q.2))).card : ℝ)
        ≤ N₂.card * hS₁fin.toFinset.card * N₃.card := by
      have := Finset.card_image_le (s := (N₂ ×ˢ hS₁fin.toFinset) ×ˢ N₃)
        (f := fun q => (q.1.1 • q.1.2, q.2))
      rw [Finset.card_product, Finset.card_product] at this
      exact_mod_cast this
    have hN₂' : (N₂.card : ℝ) ≤ (|3 * (smax - smin) / (ξ * L)| + 2) * ((n : ℝ) ^ p / ε) := by
      have h1 : (N₂.card : ℝ) ≤ ⌈(smax - smin) / (2 * (ξ * L * ε / (6 * (n : ℝ) ^ p)))⌉₊ + 1 := by
        rw [← hr₂c]
        exact_mod_cast hN₂card
      have h2 : (smax - smin) / (2 * (ξ * L * ε / (6 * (n : ℝ) ^ p)))
          = 3 * (smax - smin) / (ξ * L) * ((n : ℝ) ^ p / ε) := by
        field_simp
        ring
      have h3 : 1 ≤ (n : ℝ) ^ p / ε := by
        rw [le_div_iff₀ hε]
        linarith
      have h4 := hceil ((smax - smin) / (2 * (ξ * L * ε / (6 * (n : ℝ) ^ p))))
      rw [h2, abs_mul, abs_of_pos (by positivity : (0 : ℝ) < (n : ℝ) ^ p / ε)] at h4
      rw [h2] at h1
      nlinarith [abs_nonneg (3 * (smax - smin) / (ξ * L))]
    have hN₃' : (N₃.card : ℝ) ≤ (|3 * (σmax - σmin) / (ξ * L)| + 2) * (1 / ε) := by
      have h1 : (N₃.card : ℝ) ≤ ⌈(σmax - σmin) / (2 * (ξ * L * ε / 6))⌉₊ + 1 := by
        rw [← hr₃c]
        exact_mod_cast hN₃card
      have h2 : (σmax - σmin) / (2 * (ξ * L * ε / 6)) = 3 * (σmax - σmin) / (ξ * L) * (1 / ε) := by
        field_simp
        ring
      have h3 : 1 ≤ 1 / ε := by
        rw [le_div_iff₀ hε]
        linarith
      have h4 := hceil ((σmax - σmin) / (2 * (ξ * L * ε / 6)))
      rw [h2, abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 1 / ε)] at h4
      rw [h2] at h1
      nlinarith [abs_nonneg (3 * (σmax - σmin) / (ξ * L))]
    have hN₁' : (hS₁fin.toFinset.card : ℝ) ≤ Real.exp (K₀ * (n * ε ^ 2)) := by
      refine hN₁card.trans ?_
      rcases Nat.eq_zero_or_pos N₀ with h0 | h0
      · rw [h0, Nat.cast_zero]
        exact (Real.exp_pos _).le
      · rw [← Real.exp_log (Nat.cast_pos.2 h0), ← mul_assoc]
        exact Real.exp_le_exp.2 hlogN₀
    have hA : (|3 * (smax - smin) / (ξ * L)| + 2) * (|3 * (σmax - σmin) / (ξ * L)| + 2)
        ≤ Real.exp (n * ε ^ 2 / 3) := by
      have hpos :
          0 < (|3 * (smax - smin) / (ξ * L)| + 2) * (|3 * (σmax - σmin) / (ξ * L)| + 2) := by
        positivity
      rw [← Real.log_le_iff_le_exp hpos]
      linarith
    have hP : (n : ℝ) ^ p ≤ Real.exp (n * ε ^ 2 / 3) := by
      rw [← Real.log_le_iff_le_exp hnp0, Real.log_pow]
      have hlog0 : 0 ≤ Real.log n := Real.log_nonneg hn
      have h1 : ((p : ℝ) + 1) * Real.log n ≤ n * ε ^ 2 / 3 := by
        have := mul_le_mul_of_nonneg_left hlogn (by positivity : (0 : ℝ) ≤ (p : ℝ) + 1)
        have e : ((p : ℝ) + 1) * (1 / (3 * ((p : ℝ) + 1)) * (n * ε ^ 2)) = n * ε ^ 2 / 3 := by
          field_simp
        rwa [e] at this
      linarith
    have hE : (1 / ε) ^ 2 ≤ Real.exp (n * ε ^ 2 / 3) := by
      have hεsq : 0 < (1 / ε) ^ 2 := by positivity
      rw [← Real.log_le_iff_le_exp hεsq, Real.log_pow]
      rw [one_mul] at hlogε
      have : ((2 : ℕ) : ℝ) = 2 := by norm_num
      rw [this]
      linarith
    have := card_le_exp_of_slack hε hnp0 hX0 (Nat.cast_nonneg _) (Nat.cast_nonneg _) hN₁' hN₂'
      hN₃' hA hP hE hL2
    calc ((((N₂ ×ˢ hS₁fin.toFinset) ×ˢ N₃).image (fun q => (q.1.1 • q.1.2, q.2))).card : ℝ)
        ≤ N₂.card * hS₁fin.toFinset.card * N₃.card := hcard
      _ ≤ Real.exp (L ^ 2 * (n * ε ^ 2)) := this
      _ = Real.exp (n * (L * ε) ^ 2) := by
          congr 1
          ring
  · -- the cover
    have hS₁B' : (hS₁fin.toFinset : Set E) ⊆ B n := by rw [hN₁coe]; exact hS₁B
    have hS₁' : IsCover (2 * (ξ / (12 * smax) * ε).toNNReal) (B n)
        (hS₁fin.toFinset : Set E) := by
      rw [hN₁coe]
      exact hS₁
    refine closure_sieve_subset_cover ev hev n hsmin hsmax.le hnp0.le hS₁fin.toFinset hS₁B' hRn
      hS₁' N₂ hN₂ N₃ hN₃ ?_
    have h2r₁ : (((2 * (ξ / (12 * smax) * ε).toNNReal : ℝ≥0)) : ℝ)
        = 2 * (ξ / (12 * smax) * ε) := by
      rw [NNReal.coe_mul, hr₁c]
      norm_num
    rw [h2r₁, hr₂c, hr₃c]
    have e1 : smax * (2 * (ξ / (12 * smax) * ε)) = ξ * ε / 6 := by
      field_simp
      ring
    have e2 : (n : ℝ) ^ p * (ξ * L * ε / (6 * (n : ℝ) ^ p)) = ξ * L * ε / 6 := by
      field_simp
    rw [e1, e2]
    have := mul_le_mul_of_nonneg_left hL1 (mul_pos hξ hε).le
    linarith

/-- **UV-2 (59′), eventually.** Under (59) uniform over the amplitude box (`h59`), the [P1] prior
puts mass at least `exp(−(C₁ + δ) n εₙ²)` on the slab of radius `εₙ`, which lies in the KL
neighbourhood of radius `6εₙ/σ₀` (UV-1): the window `[s₁, s₂]` and the noise slab contribute
`πs([s₁, s₂]) · εₙ/(σ₀ log(η_hi/η_lo))`, absorbed into `exp(−δ n εₙ²)` by UV-2b′. -/
theorem eventually_prior_mass (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (hev : ∀ n i, ‖ev n i‖ ≤ 1)
    (f₀ : E) {σ₀ : ℝ} (μW : Measure E) [IsProbabilityMeasure μW]
    (πs : Measure ℝ) [IsProbabilityMeasure πs] {smin smax : ℝ} (hsmin : 0 < smin)
    (hs : smin < smax) (hπs_pos : ∀ a b, smin ≤ a → a < b → b ≤ smax → 0 < πs (Icc a b))
    (πη : Measure ℝ) [IsProbabilityMeasure πη] {ηlo ηhi : ℝ} (hηlo : 0 < ηlo) (hη : ηlo < ηhi)
    (hπη : ∀ a b, ηlo ≤ a → a ≤ b → b ≤ ηhi →
      ENNReal.ofReal (Real.log (b / a) / Real.log (ηhi / ηlo)) ≤ πη (Icc a b))
    (hσ₀ : smin * √ηlo < σ₀ ∧ σ₀ < smax * √ηhi)
    {C₀ a k C₁ δ : ℝ} (hC₀ : 0 < C₀) (ha : 0 < a) (ha' : a < 1 / 2) (hk : 0 ≤ k) (hδ : 0 < δ)
    (h59 : ∀ᶠ n : ℕ in atTop, ∀ s ∈ Icc smin smax,
      ENNReal.ofReal (Real.exp (-(C₁ * n * priorMassRadius C₀ a k n ^ 2)))
        ≤ μW (closedBall (s⁻¹ • f₀) (priorMassRadius C₀ a k n / s))) :
    ∀ᶠ n : ℕ in atTop, ∃ S ⊆ klBall ev n (f₀, σ₀) (priorMassRadius C₀ a k n / (σ₀ / 6)),
      MeasurableSet S ∧
      ENNReal.ofReal (Real.exp (-((C₁ + δ) * (σ₀ / 6) ^ 2 * n
          * (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2)))
        ≤ amplitudeNoisePrior πs μW πη S := by
  have hσ₀pos : 0 < σ₀ := lt_of_le_of_lt (by positivity) hσ₀.1
  obtain ⟨s₁, s₂, hs₁, hs₁₂, hs₂, hs₂lo, hs₁hi⟩ := exists_amplitude_window hsmin hs hηlo hη hσ₀
  have hs₁pos : 0 < s₁ := hsmin.trans_le hs₁
  have hS₀ : 0 < πs (Icc s₁ s₂) := hπs_pos s₁ s₂ hs₁ hs₁₂ hs₂
  have hS₀top : πs (Icc s₁ s₂) ≠ ⊤ := measure_ne_top _ _
  set t := (πs (Icc s₁ s₂)).toReal with ht
  have htpos : 0 < t := ENNReal.toReal_pos hS₀.ne' hS₀top
  set Lη := Real.log (ηhi / ηlo) with hLη
  have hLη0 : 0 < Lη := Real.log_pos ((one_lt_div hηlo).2 hη)
  set c₁ : ℝ := t / (σ₀ * Lη) with hc₁
  have hc₁0 : 0 < c₁ := by positivity
  set m : ℝ := min (σ₀ / 2) (s₁ * √ηhi - σ₀) with hm
  have hm0 : 0 < m := lt_min (by positivity) (by linarith)
  filter_upwards [h59, eventually_priorMassRadius_pos C₀ a k hC₀,
    (tendsto_priorMassRadius_zero C₀ a k ha).eventually (gt_mem_nhds hm0),
    log_inv_mul_radius_le C₀ a k c₁ δ hC₀ ha ha' hk hc₁0 hδ] with n h59n hρ hρm habs
  set ρ := priorMassRadius C₀ a k n with hρ_def
  have hρσ : ρ ≤ σ₀ / 2 := hρm.le.trans (min_le_left _ _)
  have hρgap : ρ ≤ s₁ * √ηhi - σ₀ := hρm.le.trans (min_le_right _ _)
  refine ⟨closedBall f₀ ρ ×ˢ closedBall σ₀ ρ, ?_,
    measurableSet_closedBall.prod measurableSet_closedBall, ?_⟩
  · have e : σ₀ * (ρ / (σ₀ / 6)) / 6 = ρ := by field_simp
    have h := slab_subset_klBall ev hev n (f₀ := f₀) hσ₀pos (ε := ρ / (σ₀ / 6))
      (by positivity) (by rw [div_le_iff₀ (by positivity)]; linarith)
    rwa [e] at h
  · have hW : ∀ s ∈ Icc s₁ s₂, ENNReal.ofReal (Real.exp (-(C₁ * n * ρ ^ 2)))
        ≤ μW (closedBall (s⁻¹ • f₀) (ρ / s)) :=
      fun s hs' => h59n s ⟨hs₁.trans hs'.1, hs'.2.trans hs₂⟩
    have hηs : ∀ s ∈ Icc s₁ s₂, ENNReal.ofReal (ρ / (σ₀ * Lη))
        ≤ πη {η | |s * √η - σ₀| ≤ ρ} := by
      intro s hs'
      have hsqlo : 0 ≤ √ηlo := Real.sqrt_nonneg _
      have hsqhi : 0 ≤ √ηhi := Real.sqrt_nonneg _
      refine noise_slab_mass πη hηlo hη hπη (hs₁pos.trans_le hs'.1) hσ₀pos hρ (by linarith) ?_ ?_
      · exact (mul_le_mul_of_nonneg_right hs'.2 hsqlo).trans hs₂lo
      · have := mul_le_mul_of_nonneg_right hs'.1 hsqhi
        linarith
    refine le_trans ?_ (prior_slab_mass πs μW πη hs₁pos hW hηs)
    rw [← ENNReal.ofReal_toReal hS₀top, ← ht, ← ENNReal.ofReal_mul (Real.exp_pos _).le,
      ← ENNReal.ofReal_mul htpos.le]
    apply ENNReal.ofReal_le_ofReal
    -- exp(−δ n ρ²) ≤ c₁ ρ, from UV-2b′
    have hcρ : 0 < c₁ * ρ := mul_pos hc₁0 hρ
    have h1 : 1 / (c₁ * ρ) ≤ Real.exp (δ * n * ρ ^ 2) :=
      (Real.log_le_iff_le_exp (one_div_pos.2 hcρ)).1 habs
    have h2 : Real.exp (-(δ * n * ρ ^ 2)) ≤ c₁ * ρ := by
      rw [Real.exp_neg]
      rw [one_div] at h1
      exact (inv_le_comm₀ (Real.exp_pos _) hcρ).2 h1
    have e : (C₁ + δ) * (σ₀ / 6) ^ 2 * n * (ρ / (σ₀ / 6)) ^ 2
        = C₁ * n * ρ ^ 2 + δ * n * ρ ^ 2 := by
      field_simp
    rw [e, neg_add, Real.exp_add]
    calc Real.exp (-(C₁ * n * ρ ^ 2)) * Real.exp (-(δ * n * ρ ^ 2))
        ≤ Real.exp (-(C₁ * n * ρ ^ 2)) * (c₁ * ρ) :=
          mul_le_mul_of_nonneg_left h2 (Real.exp_pos _).le
      _ = t * (Real.exp (-(C₁ * n * ρ ^ 2)) * (ρ / (σ₀ * Lη))) := by
          rw [hc₁]
          ring

/-! ### UV-6 — the composed theorem -/

/-- **UV-6 (the extended theorem).** Fixed design with `‖ev n i‖ ≤ 1` (the sup norm dominates
`‖·‖ₙ`), truth `(f*, σ₀)`, and the proper prior [P1]:
* `s ~ πs`, supported in `[s_min, s_max]` with positive mass on every subinterval (a positive
  density);
* `W ~ μW`, the GP law of `Wᵗ` with the empirical-Bayes bandwidth integrated out;
* `η ~ πη`, log-uniform on `[η_lo, η_hi]` (at least the log-uniform mass on subintervals);
* the one requirement on the truth, `σ₀ ∈ (s_min √η_lo, s_max √η_hi) = (σ_min, σ_max)`.

TWCD2025's framework enters through hypotheses, at the paper's radii
`εₙ = C₀ n^(−a)(log n)^k` and `ε̄ₙ = C̄ n^(−a)(log n)^k̄` (`C̄, k̄` are `Cb, kb`; `k ≤ k̄`,
`0 < a < 1/2`):
* `h59`: (59) for `f*/s`, uniformly over the box, via Remark A.1 (**UV-2, to confirm**);
* `h60`: (60), the sieve `Bₙ`'s remaining GP mass, with exponent `C₂ > C₁ + 144/σ₀²`;
* `h61`: (61), the sieve's sup-norm entropy at any fixed multiple of `ε̄ₙ`;
* `hRadius`: the sieve's sup-norm radius `Rₙ ≤ n^p` (**UV-4, to confirm**).

Then for some `M`, the expected posterior mass of `{‖f − f*‖ₙ + |σ − σ₀| ≥ M ε̄ₙ}` tends to `0`:
the known-σ rate. The proof is `gvdv07_contraction` fed by UV-1 + UV-2 (`eventually_prior_mass`),
UV-3 (`prior_compl_sieve_le`), UV-4 + UV-4a (`eventually_sieve_cover`) and UV-5 (`exists_tests`,
the open obligation). -/
theorem posterior_contraction_unknown_variance
    (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (hev : ∀ n i, ‖ev n i‖ ≤ 1) (f₀ : E) {σ₀ : ℝ}
    (μW : Measure E) [IsProbabilityMeasure μW]
    (πs : Measure ℝ) [IsProbabilityMeasure πs] {smin smax : ℝ} (hsmin : 0 < smin)
    (hs : smin < smax) (hπs_supp : πs (Icc smin smax)ᶜ = 0)
    (hπs_pos : ∀ a b, smin ≤ a → a < b → b ≤ smax → 0 < πs (Icc a b))
    (πη : Measure ℝ) [IsProbabilityMeasure πη] {ηlo ηhi : ℝ} (hηlo : 0 < ηlo) (hη : ηlo < ηhi)
    (hπη_supp : πη (Icc ηlo ηhi)ᶜ = 0)
    (hπη : ∀ a b, ηlo ≤ a → a ≤ b → b ≤ ηhi →
      ENNReal.ofReal (Real.log (b / a) / Real.log (ηhi / ηlo)) ≤ πη (Icc a b))
    (hσ₀ : smin * √ηlo < σ₀ ∧ σ₀ < smax * √ηhi)
    {C₀ Cb a k kb : ℝ} (hC₀ : 0 < C₀) (hCb : 0 < Cb) (ha : 0 < a) (ha' : a < 1 / 2)
    (hk : 0 ≤ k) (hkk : k ≤ kb) {C₁ C₂ : ℝ} (hC₂ : C₁ + 144 / σ₀ ^ 2 < C₂)
    (h59 : ∀ᶠ n : ℕ in atTop, ∀ s ∈ Icc smin smax,
      ENNReal.ofReal (Real.exp (-(C₁ * n * priorMassRadius C₀ a k n ^ 2)))
        ≤ μW (closedBall (s⁻¹ • f₀) (priorMassRadius C₀ a k n / s)))
    (B : ℕ → Set E)
    (h60 : ∀ᶠ n : ℕ in atTop,
      μW (B n)ᶜ ≤ ENNReal.ofReal (Real.exp (-(C₂ * n * priorMassRadius C₀ a k n ^ 2))))
    (h61 : ∀ c > 0, ∃ K₀ : ℝ, ∀ᶠ n : ℕ in atTop, ∃ N : ℕ,
      externalCoveringNumber (c * priorMassRadius Cb a kb n).toNNReal (B n) ≤ N ∧
        Real.log N ≤ K₀ * n * priorMassRadius Cb a kb n ^ 2)
    (hRadius : ∃ p : ℕ, ∀ᶠ n : ℕ in atTop, ∀ g ∈ B n, ‖g‖ ≤ (n : ℝ) ^ p) :
    ∃ M : ℝ, Tendsto (fun n => ∫⁻ y, posterior ev (amplitudeNoisePrior πs μW πη) n y
        {θ | M * priorMassRadius Cb a kb n ≤ paramDist ev n θ (f₀, σ₀)}
        ∂dataLaw ev n (f₀, σ₀)) atTop (𝓝 0) := by
  haveI := isProbabilityMeasure_amplitudeNoisePrior πs μW πη
  have hkb : 0 ≤ kb := hk.trans hkk
  have hσ₀pos : 0 < σ₀ := lt_of_le_of_lt (by positivity) hσ₀.1
  have hsmax : 0 < smax := hsmin.trans hs
  have hc : 0 < σ₀ / 6 := by positivity
  -- UV-5: the tests, on the compact noise range
  have hσmin : 0 < smin * √ηlo := mul_pos hsmin (Real.sqrt_pos.2 hηlo)
  obtain ⟨K, ξ, hK, hξ, htests⟩ := exists_tests ev hσmin f₀ (σ₀ := σ₀) ⟨hσ₀.1.le, hσ₀.2.le⟩
  -- (61) at the radius the cover needs, and the sieve's radius
  obtain ⟨K₀, h61'⟩ := h61 (ξ / (12 * smax)) (by positivity)
  obtain ⟨p, hR⟩ := hRadius
  set L : ℝ := 1 + |K₀| + C₀ / (Cb * (σ₀ / 6)) with hL
  have hL1 : 1 + |K₀| ≤ L := by
    have : 0 ≤ C₀ / (Cb * (σ₀ / 6)) := by positivity
    linarith
  set δ : ℝ := (C₂ - C₁ - 144 / σ₀ ^ 2) / 2 with hδ_def
  have hδ : 0 < δ := by rw [hδ_def]; linarith
  -- the sieve lies in the compact noise range
  have hsieve : ∀ n, closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi))
      ⊆ (univ ×ˢ Icc (smin * √ηlo) (smax * √ηhi) : Set (E × ℝ)) := fun n =>
    closure_minimal (prod_mono (subset_univ _) subset_rfl) (isClosed_univ.prod isClosed_Icc)
  -- σ > 0 almost surely under the prior
  have hpos : ∀ᵐ θ ∂amplitudeNoisePrior πs μW πη, 0 < θ.2 := by
    rw [amplitudeNoisePrior, ae_map_iff measurable_amplitudeNoiseMap.aemeasurable
      (measurableSet_lt measurable_const measurable_snd), ae_iff]
    refine measure_mono_null (t := (Icc smin smax)ᶜ ×ˢ univ ∪ univ ×ˢ (univ ×ˢ (Icc ηlo ηhi)ᶜ))
      ?_ (measure_union_null ?_ ?_)
    · rintro ⟨s, w, η⟩ hp
      simp only [mem_setOf_eq, amplitudeNoiseMap, not_lt] at hp
      by_contra hcon
      simp only [mem_union, mem_prod, mem_compl_iff, mem_univ, and_true, true_and, not_or,
        not_not] at hcon
      obtain ⟨hs', hη'⟩ := hcon
      have : 0 < s * √η := mul_pos (hsmin.trans_le hs'.1) (Real.sqrt_pos.2 (hηlo.trans_le hη'.1))
      linarith
    · simp [Measure.prod_prod, hπs_supp]
    · simp [Measure.prod_prod, hπη_supp]
  -- the rates
  have hrate : ∀ᶠ n : ℕ in atTop, 0 < priorMassRadius C₀ a k n / (σ₀ / 6) ∧
      priorMassRadius C₀ a k n / (σ₀ / 6) ≤ L * priorMassRadius Cb a kb n := by
    filter_upwards [eventually_priorMassRadius_pos C₀ a k hC₀,
      priorMassRadius_le_mul C₀ Cb a k kb hC₀ hCb hkk,
      eventually_priorMassRadius_pos Cb a kb hCb] with n h1 h2 h3
    refine ⟨div_pos h1 hc, ?_⟩
    rw [div_le_iff₀ hc]
    have hLc : C₀ / (Cb * (σ₀ / 6)) ≤ L := by linarith [abs_nonneg K₀]
    calc priorMassRadius C₀ a k n ≤ C₀ / Cb * priorMassRadius Cb a kb n := h2
      _ = C₀ / (Cb * (σ₀ / 6)) * priorMassRadius Cb a kb n * (σ₀ / 6) := by field_simp
      _ ≤ L * priorMassRadius Cb a kb n * (σ₀ / 6) :=
          mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hLc h3.le) hc.le
  have hnε : Tendsto (fun n : ℕ => (n : ℝ) * (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2)
      atTop atTop := by
    refine ((tendsto_n_priorMassRadius_sq C₀ a k hC₀ ha ha' hk).const_mul_atTop
      (inv_pos.2 (pow_pos hc 2))).congr fun n => ?_
    field_simp
  -- UV-5 on the sieve
  have htests' : ∀ᶠ n : ℕ in atTop, ∀ r > 0,
      ∀ θ₁ ∈ closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)),
      r < paramDist ev n θ₁ (f₀, σ₀) →
      ∃ φ : (Fin n → ℝ) → ℝ, Measurable φ ∧ (∀ y, 0 ≤ φ y ∧ φ y ≤ 1) ∧
        ∫ y, φ y ∂dataLaw ev n (f₀, σ₀) ≤ Real.exp (-(K * n * r ^ 2)) ∧
        ∀ θ ∈ closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)),
          paramDist ev n θ θ₁ ≤ ξ * r →
          ∫ y, (1 - φ y) ∂dataLaw ev n θ ≤ Real.exp (-(K * n * r ^ 2)) :=
    Eventually.of_forall fun n r hr θ₁ hθ₁ hd => by
      obtain ⟨φ, h1, h2, h3, h4⟩ := htests n r hr θ₁ (hsieve n hθ₁) hd
      exact ⟨φ, h1, h2, h3, fun θ hθ hθθ₁ => h4 θ (hsieve n hθ) hθθ₁⟩
  -- UV-3: the remaining mass, in the exponent the general theorem asks for
  have hrem : ∀ᶠ n : ℕ in atTop, amplitudeNoisePrior πs μW πη
      (closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)))ᶜ
      ≤ ENNReal.ofReal (Real.exp (-(((C₁ + δ) * (σ₀ / 6) ^ 2 + 4) * n
          * (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2))) := by
    filter_upwards [h60] with n h60n
    refine (prior_compl_sieve_le πs μW πη hsmin.le hπs_supp hπη_supp (B n)).trans
      (h60n.trans ?_)
    apply ENNReal.ofReal_le_ofReal
    apply Real.exp_le_exp.2
    have e : ((C₁ + δ) * (σ₀ / 6) ^ 2 + 4) * n * (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2
        = (C₁ + δ + 144 / σ₀ ^ 2) * (n * priorMassRadius C₀ a k n ^ 2) := by
      field_simp
      ring
    rw [e]
    have hcoef : C₁ + δ + 144 / σ₀ ^ 2 ≤ C₂ := by rw [hδ_def]; linarith
    have hnρ : 0 ≤ (n : ℝ) * priorMassRadius C₀ a k n ^ 2 := by positivity
    have := mul_le_mul_of_nonneg_right hcoef hnρ
    linarith
  -- the general theorem
  obtain ⟨M, hM⟩ := gvdv07_contraction ev (amplitudeNoisePrior πs μW πη) hpos (f₀, σ₀) hσ₀pos
    (paramDist ev) (paramDist_comm ev) (paramDist_triangle ev)
    (fun n θ' => (continuous_paramDist ev n θ').measurable)
    (fun n => priorMassRadius C₀ a k n / (σ₀ / 6)) (fun n => L * priorMassRadius Cb a kb n)
    hrate hnε (fun n => closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)))
    (fun n => isClosed_closure.measurableSet) ((C₁ + δ) * (σ₀ / 6) ^ 2) K ξ hK hξ htests'
    (eventually_sieve_cover ev hev B hsmin.le hsmax hCb ha ha' hkb hξ hL1 h61' hR)
    (eventually_prior_mass ev hev f₀ μW πs hsmin hs hπs_pos πη hηlo hη hπη hσ₀ hC₀ ha ha' hk hδ
      h59)
    hrem
  refine ⟨M * L, hM.congr fun n => ?_⟩
  simp only [mul_assoc]

end BGP
