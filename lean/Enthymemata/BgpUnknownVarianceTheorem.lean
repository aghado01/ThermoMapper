import Enthymemata.BgpCovering
import Enthymemata.BgpGaussianKL
import Enthymemata.BgpGaussianTails
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
* the tests (UV-5, `exists_tests`), proved from the Gaussian tails of `BgpGaussianTails.lean`;
* the general theorem, as a cited axiom in the style of `Stability.lean`, in two forms.

**The general theorem, two ways.**
* `gvdv07_theorem1` and `gvdv07_lemma1` are Ghosal & van der Vaart (2007), Theorem 1 and Lemma 1,
  as printed, for independent observations. They give
  `posterior_contraction_unknown_variance_gvdv07`:
  contraction at `Mₙ ε̄ₙ` for **every** `Mₙ → ∞`.
* `fixed_constant_contraction` is the fixed-`M`, two-rate form that TWCD2025 (Theorem 3.1),
  van der Vaart & van Zanten (2008, Theorem 3.3; 2009, Section 3) and the protolemma state. Its
  source is Ghosal & van der Vaart (2001, Theorem 2.1), the two-rate version of Ghosal, Ghosh &
  van der Vaart (2000, Theorem 2.1). For independent non-identically distributed observations,
  with GvdV07's tests, it is not printed as a theorem. It gives
  `posterior_contraction_unknown_variance`: contraction at `M ε̄ₙ` for one fixed `M`.

**What is cited from TWCD2025, and how.** The paper's (59)–(61) and Remark A.1 enter as hypotheses:
they are facts about the paper's GP prior and sieve, which this file does not construct. Two carry
the protolemma's "to confirm against the paper" checks, both confirmed against Appendix A.1:
* `h59`: (59) uniform over the amplitude box. Part I uses `f*` only through (A2)'s constants
  `(ν₁, ν₂)`, and (A2) is monotone in them, so `f*/s` satisfies it with `(ν₁/s_min, ν₂/s_min²)`
  for every `s` in the box.
* `hRadius`: the sieve's sup-norm radius is at most `Nₙ (rₙ/δₙ)^{D/2} + εₙ`, polynomial in `n`,
  since the kernel has `h(0) = 1`.

The remaining-mass exponent `C₂` must exceed `C₁` by a margin (`hC₂`; see the protolemma under
UV-3). The paper's sieve is tunable, which `posterior_contraction_unknown_variance_of_tunable_sieve`
records.

**The ledger.** Every proof here is complete, so the ledger reads PROOF-CLOSED. The axioms keep
the module's results conditional on the cited theorems; `#print axioms` names them per theorem.
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

/-- **AXIOM — the general contraction theorem, fixed-constant two-rate form.**
For independent, non-identically distributed observations with a dominated model, a prior that
puts mass at least `exp(−C n εₙ²)` on (a measurable part of) the KL neighbourhood
`B*ₙ(θ₀, εₙ; 2)`, a sieve `Θₙ` whose complement has mass at most `exp(−(C + 4) n εₙ²)`, a cover
of the sieve by `exp(n ε̄ₙ²)` balls of radius `ξ ε̄ₙ/2` in a semimetric `d`, and exponentially
powerful local tests for `d` give a posterior that contracts at rate `M ε̄ₙ` in `d`, in mean under
the truth.

Cited, not formalized. Checked against the sources (protolemma, UV-6):
* **The printed Theorem 4 of Ghosal & van der Vaart (2007)** has one rate and concludes for every
  `Mₙ → ∞`. Its tests may be for any semimetric with Lemma 2's conclusion. That form is
  `gvdv07_theorem1` with `gvdv07_lemma1`.
* **The fixed-`M` two-rate form here** is Ghosal & van der Vaart (2001, Theorem 2.1), the two-rate
  version of Ghosal, Ghosh & van der Vaart (2000, Theorem 2.1). That is the form van der Vaart &
  van Zanten (2009, Section 3: "a slight generalization of the main result in [GGvdV00],
  formulated as Theorem 2.1 in [GvdV01]") and TWCD2025 (Appendix A.1) invoke. Van der Vaart &
  van Zanten (2008, Theorem 3.3) state it for fixed-design GP regression with an unknown scale on
  `[a, b]`, with the tests left as in GvdV07.

The hypotheses are stated on the strong side, so the axiom stays a consequence of that argument:
* the tests (GvdV07's (2.2)) at every radius `r > 0`;
* a global cover of the sieve;
* measurable sieves and neighbourhood parts;
* `n εₙ² → ∞`. -/
axiom fixed_constant_contraction
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

/-! ### The general theorem as printed: Ghosal & van der Vaart (2007), Theorem 1 and Lemma 1 -/

/-- **AXIOM — Ghosal & van der Vaart (2007), Theorem 1**, for independent observations, with
`eₙ = dₙ`, `k = 2`, the parameter set `S` and the KL neighbourhood `B*ₙ(θ₀, εₙ; 2)` in its
independent-observations form (`klBall`). The hypotheses are the theorem's:
* `Θₙ ⊆ S`, the sieves inside the parameter set;
* (2.2): for every `r > 0` and `θ₁ ∈ S` with `dₙ(θ₁, θ₀) > r`, a test with
  `P₀φ ≤ e^{−Knr²}` and `Pθ(1 − φ) ≤ e^{−Knr²}` on the `ξr`-ball around `θ₁`;
* (2.4): the local entropy `log N(ξr/2, {θ ∈ Θₙ : dₙ(θ, θ₀) < r}, dₙ) ≤ nεₙ²` for `r > εₙ`, here
  by closed balls of radius `ξr/4`;
* (2.5): `Πₙ(θ ∈ Θₙ : jεₙ < dₙ ≤ 2jεₙ) ≤ e^{Knεₙ²j²/2} Πₙ(B*ₙ(θ₀, εₙ; 2))` for every sufficiently
  large `j`;
* `εₙ → 0` and `(nεₙ²)⁻¹ = O(1)`.

The conclusion is the printed one: for every `Mₙ → ∞`, the expected posterior mass of
`{θ ∈ Θₙ : dₙ(θ, θ₀) ≥ Mₙεₙ}` tends to `0`. Hypotheses are on the strong side: tests for every `n`,
closed balls at half the radius, measurable sieves and neighbourhood parts, eventually in `n`. -/
axiom gvdv07_theorem1
    (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (prior : Measure (E × ℝ)) [IsProbabilityMeasure prior]
    (S : Set (E × ℝ)) (hS : prior Sᶜ = 0) (hS_pos : ∀ θ ∈ S, 0 < θ.2) (θ₀ : E × ℝ) (hθ₀ : θ₀ ∈ S)
    (d : ℕ → E × ℝ → E × ℝ → ℝ) (hd_comm : ∀ n θ θ', d n θ θ' = d n θ' θ)
    (hd_triangle : ∀ n θ θ' θ'', d n θ θ'' ≤ d n θ θ' + d n θ' θ'')
    (hd_meas : ∀ n θ', Measurable fun θ => d n θ θ')
    (ε : ℕ → ℝ) (hε : ∀ᶠ n in atTop, 0 < ε n) (hε0 : Tendsto ε atTop (𝓝 0))
    (hnε : ∃ c > 0, ∀ᶠ n in atTop, c ≤ n * ε n ^ 2)
    (Θ : ℕ → Set (E × ℝ)) (hΘ : ∀ n, MeasurableSet (Θ n)) (hΘS : ∀ n, Θ n ⊆ S)
    (K ξ : ℝ) (hK : 0 < K) (hξ : 0 < ξ)
    (htests : ∀ n : ℕ, ∀ r > 0, ∀ θ₁ ∈ S, r < d n θ₁ θ₀ →
      ∃ φ : (Fin n → ℝ) → ℝ, Measurable φ ∧ (∀ y, 0 ≤ φ y ∧ φ y ≤ 1) ∧
        ∫ y, φ y ∂dataLaw ev n θ₀ ≤ Real.exp (-(K * n * r ^ 2)) ∧
        ∀ θ ∈ S, d n θ θ₁ ≤ ξ * r →
          ∫ y, (1 - φ y) ∂dataLaw ev n θ ≤ Real.exp (-(K * n * r ^ 2)))
    (hentropy : ∀ᶠ n in atTop, ∀ r > ε n, ∃ N : Finset (E × ℝ),
      (N.card : ℝ) ≤ Real.exp (n * ε n ^ 2) ∧
        {θ | θ ∈ Θ n ∧ d n θ θ₀ < r} ⊆ ⋃ θ' ∈ N, {θ | d n θ θ' ≤ ξ * r / 4})
    (hshell : ∃ j₀ : ℕ, ∀ᶠ n in atTop, ∃ B ⊆ klBall ev n θ₀ (ε n), MeasurableSet B ∧
      ∀ j : ℕ, j₀ ≤ j →
        prior {θ | θ ∈ Θ n ∧ j * ε n < d n θ θ₀ ∧ d n θ θ₀ ≤ 2 * j * ε n}
          ≤ ENNReal.ofReal (Real.exp (K * n * ε n ^ 2 * j ^ 2 / 2)) * prior B) :
    ∀ M : ℕ → ℝ, Tendsto M atTop atTop →
      Tendsto (fun n => ∫⁻ y, posterior ev prior n y {θ | θ ∈ Θ n ∧ M n * ε n ≤ d n θ θ₀}
        ∂dataLaw ev n θ₀) atTop (𝓝 0)

/-- **AXIOM — Ghosal & van der Vaart (2007), Lemma 1**, for independent observations and `k = 2`:
if `Πₙ(Θ ∖ Θₙ)/Πₙ(B*ₙ(θ₀, εₙ; 2)) = o(e^{−2nεₙ²})`, the posterior mass outside the sieve tends to
`0` in mean. The denominator is taken on a measurable part of the neighbourhood with positive mass,
which only strengthens the hypothesis. -/
axiom gvdv07_lemma1
    (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (prior : Measure (E × ℝ)) [IsProbabilityMeasure prior]
    (S : Set (E × ℝ)) (hS : prior Sᶜ = 0) (hS_pos : ∀ θ ∈ S, 0 < θ.2) (θ₀ : E × ℝ) (hθ₀ : θ₀ ∈ S)
    (ε : ℕ → ℝ) (hε : ∀ᶠ n in atTop, 0 < ε n) (hnε : ∃ c > 0, ∀ᶠ n in atTop, c ≤ n * ε n ^ 2)
    (Θ : ℕ → Set (E × ℝ)) (hΘ : ∀ n, MeasurableSet (Θ n)) (B : ℕ → Set (E × ℝ))
    (hB : ∀ᶠ n in atTop, B n ⊆ klBall ev n θ₀ (ε n) ∧ MeasurableSet (B n) ∧ 0 < prior (B n))
    (hratio : Tendsto (fun n => prior (Θ n)ᶜ / prior (B n)
      * ENNReal.ofReal (Real.exp (2 * n * ε n ^ 2))) atTop (𝓝 0)) :
    Tendsto (fun n => ∫⁻ y, posterior ev prior n y (Θ n)ᶜ ∂dataLaw ev n θ₀) atTop (𝓝 0)

/-! ### UV-5 — the tests -/

omit [MeasurableSpace E] [BorelSpace E] in
instance isProbabilityMeasure_dataLaw (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (n : ℕ) (θ : E × ℝ) :
    IsProbabilityMeasure (dataLaw ev n θ) := by
  unfold dataLaw
  infer_instance

omit [MeasurableSpace E] [BorelSpace E] in
theorem empiricalNorm_sq (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (n : ℕ) (g : E) :
    empiricalNorm ev n g ^ 2 = (∑ i, ev n i g ^ 2) / n := by
  rw [empiricalNorm, Real.sq_sqrt (div_nonneg (Finset.sum_nonneg fun i _ => sq_nonneg _)
    (Nat.cast_nonneg n))]

/-- A rejection set's mass bounds the test `1_R`'s first error. -/
theorem integral_indicator_le {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    {R : Set Ω} (hR : MeasurableSet R) {b : ℝ} (hb : 0 ≤ b) (h : P R ≤ ENNReal.ofReal b) :
    ∫ y, R.indicator (1 : Ω → ℝ) y ∂P ≤ b := by
  rw [integral_indicator_one hR]
  exact ENNReal.toReal_le_of_le_ofReal hb h

/-- The acceptance set's mass bounds the test `1_R`'s second error. -/
theorem integral_one_sub_indicator_le {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    [IsProbabilityMeasure P] {R : Set Ω} (hR : MeasurableSet R) {b : ℝ} (hb : 0 ≤ b)
    (h : P Rᶜ ≤ ENNReal.ofReal b) :
    ∫ y, (1 - R.indicator (1 : Ω → ℝ) y) ∂P ≤ b := by
  have hpt : (fun y => 1 - R.indicator (1 : Ω → ℝ) y) = Rᶜ.indicator (1 : Ω → ℝ) := by
    funext y
    by_cases hy : y ∈ R
    · simp [hy]
    · simp [hy]
  rw [hpt]
  exact integral_indicator_le P hR.compl hb h

/-- `exp(−X) ≤ exp(−Y)` in `ℝ≥0∞` from `Y ≤ X`. -/
theorem ofReal_exp_neg_le {X Y : ℝ} (h : Y ≤ X) :
    ENNReal.ofReal (Real.exp (-X)) ≤ ENNReal.ofReal (Real.exp (-Y)) :=
  ENNReal.ofReal_le_ofReal (Real.exp_le_exp.2 (neg_le_neg h))

/-- The tests' exponent constant, `σ_min²/(1024 σ_max⁴)`. -/
noncomputable def testConst (σmin σmax : ℝ) : ℝ := σmin ^ 2 / (1024 * σmax ^ 4)

theorem testConst_pos {σmin σmax : ℝ} (hσmin : 0 < σmin) (hσmax : 0 < σmax) :
    0 < testConst σmin σmax := by
  unfold testConst
  positivity

theorem testConst_le {σmin σmax : ℝ} (hσmin : 0 ≤ σmin) (hle : σmin ≤ σmax) (hσmax : 0 < σmax) :
    testConst σmin σmax ≤ 1 / (1024 * σmax ^ 2) := by
  unfold testConst
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  have h := mul_le_mul_of_nonneg_right (pow_le_pow_left₀ hσmin hle 2)
    (by positivity : (0 : ℝ) ≤ 1024 * σmax ^ 2)
  linarith

theorem testConst_mul {σmin σmax : ℝ} (hσmax : 0 < σmax) :
    testConst σmin σmax * σmax ^ 4 = σmin ^ 2 / 1024 := by
  unfold testConst
  rw [div_mul_eq_mul_div, div_eq_div_iff (by positivity) (by norm_num)]
  ring

omit [MeasurableSpace E] [BorelSpace E] in
/-- An empirical-norm bound `‖g‖ₙ ≤ b` bounds the sum of squares by `n b²`. -/
theorem sum_sq_le_of_empiricalNorm_le (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) {n : ℕ} (hn : 0 < n)
    {g : E} {b : ℝ} (h : empiricalNorm ev n g ≤ b) : ∑ i, ev n i g ^ 2 ≤ n * b ^ 2 := by
  have h2 := pow_le_pow_left₀ (empiricalNorm_nonneg ev n g) h 2
  rw [empiricalNorm_sq, div_le_iff₀ (by exact_mod_cast hn)] at h2
  linarith

omit [MeasurableSpace E] [BorelSpace E] in
/-- **UV-5, the mean case.** When `‖f₁ − f*‖ₙ > r/2`, reject when the likelihood-ratio statistic
`Σ δᵢ(Yᵢ − f*(xᵢ))`, `δ = f₁ − f*`, exceeds `‖δ‖²/2`. Under an alternative with
`‖f − f₁‖ₙ ≤ r/8` its mean is at least `3‖δ‖²/4`, by Cauchy–Schwarz. -/
theorem tests_mean_case (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) {σmin σmax σ₀ r : ℝ}
    (hσmin : 0 < σmin) (hσ₀ : σ₀ ∈ Icc σmin σmax) (hr : 0 < r) {n : ℕ} (hn : 0 < n)
    (f₀ f₁ : E) (hA : r / 2 < empiricalNorm ev n (f₁ - f₀)) :
    ∃ R : Set (Fin n → ℝ), MeasurableSet R ∧
      dataLaw ev n (f₀, σ₀) R
        ≤ ENNReal.ofReal (Real.exp (-(testConst σmin σmax * n * r ^ 2))) ∧
      ∀ f σ, σ ∈ Icc σmin σmax → empiricalNorm ev n (f - f₁) ≤ r / 8 →
        dataLaw ev n (f, σ) Rᶜ
          ≤ ENNReal.ofReal (Real.exp (-(testConst σmin σmax * n * r ^ 2))) := by
  obtain ⟨hσ₀lo, hσ₀hi⟩ := hσ₀
  have hσ₀pos : 0 < σ₀ := hσmin.trans_le hσ₀lo
  have hσmax : 0 < σmax := hσ₀pos.trans_le hσ₀hi
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have hK1 := testConst_le hσmin.le (hσ₀lo.trans hσ₀hi) hσmax
  obtain ⟨K, hK⟩ : ∃ K, testConst σmin σmax = K := ⟨_, rfl⟩
  rw [hK] at hK1 ⊢
  have hK0 : 0 < K := hK ▸ testConst_pos hσmin hσmax
  obtain ⟨D, hD⟩ : ∃ D, ∑ i, ev n i (f₁ - f₀) ^ 2 = D := ⟨_, rfl⟩
  have hDgt : n * r ^ 2 < 4 * D := by
    have h1 : (r / 2) ^ 2 < empiricalNorm ev n (f₁ - f₀) ^ 2 :=
      pow_lt_pow_left₀ hA (by positivity) two_ne_zero
    rw [empiricalNorm_sq, hD, lt_div_iff₀ hnR] at h1
    linarith
  have hDpos : 0 < D := by
    have := mul_pos hnR (pow_pos hr 2)
    linarith
  -- the common bound K n r² ≤ D/(32 σmax²)
  have hKD : K * n * r ^ 2 ≤ D / (32 * σmax ^ 2) := by
    rw [le_div_iff₀ (by positivity)]
    have h1 : K * (32 * σmax ^ 2) ≤ 1 / 32 := by
      have := mul_le_mul_of_nonneg_right hK1 (by positivity : (0 : ℝ) ≤ 32 * σmax ^ 2)
      rw [show 1 / (1024 * σmax ^ 2) * (32 * σmax ^ 2) = 1 / 32 by
        rw [div_mul_eq_mul_div, one_mul, div_eq_div_iff (by positivity) (by norm_num)]
        ring] at this
      exact this
    have h2 := mul_le_mul_of_nonneg_right h1 (by positivity : (0 : ℝ) ≤ n * r ^ 2)
    linarith
  have hR : MeasurableSet {y : Fin n → ℝ | D / 2 ≤ ∑ i, ev n i (f₁ - f₀) * (y i - ev n i f₀)} :=
    measurableSet_le measurable_const (Finset.measurable_sum _ fun i _ =>
      measurable_const.mul ((measurable_pi_apply i).sub measurable_const))
  refine ⟨_, hR, ?_, ?_⟩
  · -- the truth: P₀(T ≥ D/2) ≤ exp(−D/(8σ₀²))
    have h := pi_gaussian_linear_tail (fun i => ev n i f₀) (fun i => ev n i (f₁ - f₀))
      (v := σ₀ ^ 2) (t := D / 2) (by positivity) (by simp only [hD]; exact hDpos) (by linarith)
    simp only [hD] at h
    refine h.trans (ofReal_exp_neg_le ?_)
    have e : (D / 2) ^ 2 / (2 * σ₀ ^ 2 * D) = D / (8 * σ₀ ^ 2) := by
      rw [div_eq_div_iff (by positivity) (by positivity)]
      ring
    have h3 : D / (32 * σmax ^ 2) ≤ D / (8 * σ₀ ^ 2) := by
      apply div_le_div_of_nonneg_left hDpos.le (by positivity)
      have := pow_le_pow_left₀ hσ₀pos.le hσ₀hi 2
      linarith [sq_nonneg σmax]
    rw [e]
    linarith
  · -- an alternative within r/8 of f₁: P(T < D/2) ≤ exp(−D/(32σ²))
    intro f σ hσ hf1
    have hσpos : 0 < σ := hσmin.trans_le hσ.1
    -- the mean of T under (f, σ) is at least 3D/4
    obtain ⟨mθ, hmθ⟩ : ∃ mθ, ∑ i, ev n i (f₁ - f₀) * (ev n i f - ev n i f₀) = mθ := ⟨_, rfl⟩
    have hmθ_ge : 3 * D / 4 ≤ mθ := by
      have hsplit : mθ = D + ∑ i, ev n i (f₁ - f₀) * ev n i (f - f₁) := by
        rw [← hmθ, ← hD, ← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun i _ => ?_
        simp only [map_sub]
        ring
      have hcs := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (fun i => ev n i (f₁ - f₀))
        (fun i => ev n i (f - f₁))
      rw [hD] at hcs
      have hUf := sum_sq_le_of_empiricalNorm_le ev hn hf1
      obtain ⟨s, hs⟩ : ∃ s, ∑ i, ev n i (f₁ - f₀) * ev n i (f - f₁) = s := ⟨_, rfl⟩
      rw [hs] at hcs hsplit
      have h1 := mul_le_mul_of_nonneg_left hUf hDpos.le
      have h2 := mul_lt_mul_of_pos_left hDgt hDpos
      have hs2 : s ^ 2 < (D / 4) ^ 2 := by nlinarith
      have := (abs_lt.1 (abs_lt_of_sq_lt_sq hs2 (by positivity))).1
      linarith
    have hsub : {y : Fin n → ℝ | D / 2 ≤ ∑ i, ev n i (f₁ - f₀) * (y i - ev n i f₀)}ᶜ
        ⊆ {y | mθ - D / 2 ≤ ∑ i, (fun i => -ev n i (f₁ - f₀)) i * (y i - ev n i f)} := by
      intro y hy
      simp only [Set.mem_compl_iff, Set.mem_setOf_eq, not_le] at hy ⊢
      have e : ∑ i, -ev n i (f₁ - f₀) * (y i - ev n i f)
          = mθ - ∑ i, ev n i (f₁ - f₀) * (y i - ev n i f₀) := by
        rw [← hmθ, ← Finset.sum_sub_distrib]
        exact Finset.sum_congr rfl fun i _ => by ring
      rw [e]
      linarith
    have hQ : ∑ i, (-ev n i (f₁ - f₀)) ^ 2 = D := by
      rw [← hD]
      exact Finset.sum_congr rfl fun i _ => by ring
    have h := pi_gaussian_linear_tail (fun i => ev n i f) (fun i => -ev n i (f₁ - f₀))
      (v := σ ^ 2) (t := mθ - D / 2) (by positivity) (by simp only [hQ]; exact hDpos)
      (by linarith)
    simp only [hQ] at h
    refine (measure_mono hsub).trans (h.trans (ofReal_exp_neg_le ?_))
    have ht : D / 4 ≤ mθ - D / 2 := by linarith
    have h1 : D / (32 * σ ^ 2) ≤ (mθ - D / 2) ^ 2 / (2 * σ ^ 2 * D) := by
      rw [div_le_div_iff₀ (by positivity) (by positivity)]
      have h4 := mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (by positivity) ht 2)
        (by positivity : (0 : ℝ) ≤ 32 * σ ^ 2)
      linarith
    have h4 : D / (32 * σmax ^ 2) ≤ D / (32 * σ ^ 2) := by
      apply div_le_div_of_nonneg_left hDpos.le (by positivity)
      have := pow_le_pow_left₀ hσpos.le hσ.2 2
      linarith
    linarith

omit [MeasurableSpace E] [BorelSpace E] in
/-- **UV-5, the scale-up case.** When `σ₁ > σ₀ + r/2`, reject when
`Σ(Yᵢ − f*(xᵢ))² ≥ nτ`, `τ = σ₀² + 3σ₀r/8`. Under the truth the sum is central; under an
alternative with `|σ − σ₁| ≤ r/8` only its lower tail matters, where the means drop out. -/
theorem tests_scale_up_case (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) {σmin σmax σ₀ σ₁ r : ℝ}
    (hσmin : 0 < σmin) (hσ₀ : σ₀ ∈ Icc σmin σmax) (hσ₁max : σ₁ ≤ σmax) (hr : 0 < r) {n : ℕ}
    (hn : 0 < n) (f₀ : E) (hB1 : σ₀ + r / 2 < σ₁) :
    ∃ R : Set (Fin n → ℝ), MeasurableSet R ∧
      dataLaw ev n (f₀, σ₀) R
        ≤ ENNReal.ofReal (Real.exp (-(testConst σmin σmax * n * r ^ 2))) ∧
      ∀ f σ, σ ∈ Icc σmin σmax → |σ - σ₁| ≤ r / 8 →
        dataLaw ev n (f, σ) Rᶜ
          ≤ ENNReal.ofReal (Real.exp (-(testConst σmin σmax * n * r ^ 2))) := by
  obtain ⟨hσ₀lo, hσ₀hi⟩ := hσ₀
  have hσ₀pos : 0 < σ₀ := hσmin.trans_le hσ₀lo
  have hσmax : 0 < σmax := hσ₀pos.trans_le hσ₀hi
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have hK1 := testConst_le hσmin.le (hσ₀lo.trans hσ₀hi) hσmax
  have hKe := testConst_mul (σmin := σmin) hσmax
  obtain ⟨K, hK⟩ : ∃ K, testConst σmin σmax = K := ⟨_, rfl⟩
  rw [hK] at hK1 hKe ⊢
  have hK0 : 0 < K := hK ▸ testConst_pos hσmin hσmax
  have hrmax : r ≤ 2 * σmax := by linarith
  obtain ⟨τ, hτ⟩ : ∃ τ, σ₀ ^ 2 + 3 / 8 * σ₀ * r = τ := ⟨_, rfl⟩
  have hτpos : 0 < τ := by rw [← hτ]; positivity
  have hR : MeasurableSet {y : Fin n → ℝ | n * τ ≤ ∑ i, (y i - ev n i f₀) ^ 2} :=
    measurableSet_le measurable_const (Finset.measurable_sum _ fun i _ =>
      ((measurable_pi_apply i).sub measurable_const).pow_const 2)
  refine ⟨_, hR, ?_, ?_⟩
  · have hM : ∑ i, (ev n i f₀ - ev n i f₀) ^ 2 = 0 := by simp
    have hlt : n * σ₀ ^ 2 + ∑ i, (ev n i f₀ - ev n i f₀) ^ 2 < n * τ := by
      rw [hM, add_zero, ← hτ]
      have := mul_pos hnR (by positivity : (0 : ℝ) < 3 / 8 * σ₀ * r)
      linarith
    have h := pi_gaussian_sq_upper_tail (fun i => ev n i f₀) (fun i => ev n i f₀)
      (v := σ₀ ^ 2) (τ := τ) (by positivity) hlt
    simp only [hM, add_zero] at h
    refine h.trans (ofReal_exp_neg_le ?_)
    have hτmax : τ ≤ 2 * σmax ^ 2 := by
      rw [← hτ]
      have h1 : σ₀ ^ 2 ≤ σmax ^ 2 := pow_le_pow_left₀ hσ₀pos.le hσ₀hi 2
      have h2 : σ₀ * r ≤ σmax * (2 * σmax) := mul_le_mul hσ₀hi hrmax hr.le hσmax.le
      nlinarith
    have h3 : K * τ ≤ 9 / 1024 := by
      have := mul_le_mul hK1 hτmax hτpos.le (by positivity)
      rw [show 1 / (1024 * σmax ^ 2) * (2 * σmax ^ 2) = 1 / 512 by
        rw [div_mul_eq_mul_div, one_mul, div_eq_div_iff (by positivity) (by norm_num)]
        ring] at this
      linarith
    have hτσ : n * τ - n * σ₀ ^ 2 = n * (3 / 8 * σ₀ * r) := by rw [← hτ]; ring
    rw [hτσ, le_div_iff₀ (by positivity)]
    have h4 := mul_le_mul_of_nonneg_left h3
      (by positivity : (0 : ℝ) ≤ 16 * σ₀ ^ 2 * (n : ℝ) ^ 2 * r ^ 2)
    linarith
  · intro f σ hσ hs1
    have hσpos : 0 < σ := hσmin.trans_le hσ.1
    have hσlow : σ₀ + 3 / 8 * r < σ := by
      have := (abs_le.1 hs1).1
      linarith
    have hsq' : (σ₀ + 3 / 8 * r) ^ 2 < σ ^ 2 := pow_lt_pow_left₀ hσlow (by positivity) two_ne_zero
    have hgap : 3 / 8 * σ₀ * r ≤ σ ^ 2 - τ := by
      rw [← hτ]
      nlinarith [sq_nonneg r]
    have hτlt : τ < σ ^ 2 := by
      have := mul_pos hσ₀pos hr
      linarith
    have hsub : {y : Fin n → ℝ | n * τ ≤ ∑ i, (y i - ev n i f₀) ^ 2}ᶜ
        ⊆ {y | ∑ i, (y i - ev n i f₀) ^ 2 ≤ n * τ} := by
      intro y hy
      simp only [Set.mem_compl_iff, Set.mem_setOf_eq, not_le] at hy ⊢
      exact hy.le
    have h := pi_gaussian_sq_lower_tail (fun i => ev n i f) (fun i => ev n i f₀)
      (v := σ ^ 2) (τ := τ) (by positivity) hτlt
    refine (measure_mono hsub).trans (h.trans (ofReal_exp_neg_le ?_))
    rw [le_div_iff₀ (by positivity)]
    have hσ4 : (σ ^ 2) ^ 2 ≤ σmax ^ 4 := by
      have := pow_le_pow_left₀ (by positivity) (pow_le_pow_left₀ hσpos.le hσ.2 2) 2
      rw [show (σmax ^ 2) ^ 2 = σmax ^ 4 by ring] at this
      exact this
    have hsq : (3 / 8 * σ₀ * r) ^ 2 ≤ (σ ^ 2 - τ) ^ 2 := pow_le_pow_left₀ (by positivity) hgap 2
    have hσ₀2 : σmin ^ 2 ≤ σ₀ ^ 2 := pow_le_pow_left₀ hσmin.le hσ₀lo 2
    have h1 := mul_le_mul_of_nonneg_left hσ4 (by positivity : (0 : ℝ) ≤ 8 * K * n * r ^ 2)
    have hKe2 : K * n * r ^ 2 * σmax ^ 4 = σmin ^ 2 / 1024 * (n * r ^ 2) := by
      rw [show K * n * r ^ 2 * σmax ^ 4 = (K * σmax ^ 4) * (n * r ^ 2) by ring, hKe]
    have h2 := mul_le_mul_of_nonneg_right hσ₀2 (by positivity : (0 : ℝ) ≤ n * r ^ 2)
    have h3 := mul_le_mul_of_nonneg_left hsq hnR.le
    have h0 : 0 ≤ σ₀ ^ 2 * (n * r ^ 2) := by positivity
    linarith

omit [MeasurableSpace E] [BorelSpace E] in
/-- **UV-5, the scale-down case.** When `σ₁ < σ₀ − r/2`, reject when
`Σ(Yᵢ − f₁(xᵢ))² ≤ nτ`, `τ = σ₀² − σ₀r/8`. Under the truth only the lower tail matters; under an
alternative with `‖f − f₁‖ₙ ≤ r/8` and `|σ − σ₁| ≤ r/8` the means are within `r/8` of `f₁`, and
`r < 2σ₀` keeps their shift below the scale gap. -/
theorem tests_scale_down_case (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) {σmin σmax σ₀ σ₁ r : ℝ}
    (hσmin : 0 < σmin) (hσ₀ : σ₀ ∈ Icc σmin σmax) (hσ₁min : σmin ≤ σ₁) (hr : 0 < r) {n : ℕ}
    (hn : 0 < n) (f₀ f₁ : E) (hB2 : σ₁ < σ₀ - r / 2) :
    ∃ R : Set (Fin n → ℝ), MeasurableSet R ∧
      dataLaw ev n (f₀, σ₀) R
        ≤ ENNReal.ofReal (Real.exp (-(testConst σmin σmax * n * r ^ 2))) ∧
      ∀ f σ, σ ∈ Icc σmin σmax → empiricalNorm ev n (f - f₁) ≤ r / 8 → |σ - σ₁| ≤ r / 8 →
        dataLaw ev n (f, σ) Rᶜ
          ≤ ENNReal.ofReal (Real.exp (-(testConst σmin σmax * n * r ^ 2))) := by
  obtain ⟨hσ₀lo, hσ₀hi⟩ := hσ₀
  have hσ₀pos : 0 < σ₀ := hσmin.trans_le hσ₀lo
  have hσmax : 0 < σmax := hσ₀pos.trans_le hσ₀hi
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have hK1 := testConst_le hσmin.le (hσ₀lo.trans hσ₀hi) hσmax
  obtain ⟨K, hK⟩ : ∃ K, testConst σmin σmax = K := ⟨_, rfl⟩
  rw [hK] at hK1 ⊢
  have hK0 : 0 < K := hK ▸ testConst_pos hσmin hσmax
  have hr2 : r < 2 * σ₀ := by linarith
  obtain ⟨τ, hτ⟩ : ∃ τ, σ₀ ^ 2 - σ₀ * r / 8 = τ := ⟨_, rfl⟩
  have hσr : σ₀ * r < 2 * σ₀ ^ 2 := by nlinarith
  have hτpos : 0 < τ := by
    rw [← hτ]
    linarith [pow_pos hσ₀pos 2]
  have hτle : τ ≤ σ₀ ^ 2 := by
    rw [← hτ]
    have := mul_pos hσ₀pos hr
    linarith
  have hR : MeasurableSet {y : Fin n → ℝ | ∑ i, (y i - ev n i f₁) ^ 2 ≤ n * τ} :=
    measurableSet_le (Finset.measurable_sum _ fun i _ =>
      ((measurable_pi_apply i).sub measurable_const).pow_const 2) measurable_const
  refine ⟨_, hR, ?_, ?_⟩
  · have hτlt : τ < σ₀ ^ 2 := by
      rw [← hτ]
      have := mul_pos hσ₀pos hr
      linarith
    have h := pi_gaussian_sq_lower_tail (fun i => ev n i f₀) (fun i => ev n i f₁)
      (v := σ₀ ^ 2) (τ := τ) (by positivity) hτlt
    refine h.trans (ofReal_exp_neg_le ?_)
    have hτσ : σ₀ ^ 2 - τ = σ₀ * r / 8 := by rw [← hτ]; ring
    rw [hτσ, le_div_iff₀ (by positivity)]
    have h3 : K * σ₀ ^ 2 ≤ 1 / 1024 := by
      have hσ2 : σ₀ ^ 2 ≤ σmax ^ 2 := pow_le_pow_left₀ hσ₀pos.le hσ₀hi 2
      have := mul_le_mul hK1 hσ2 (by positivity) (by positivity)
      rw [show 1 / (1024 * σmax ^ 2) * σmax ^ 2 = 1 / 1024 by
        rw [div_mul_eq_mul_div, one_mul, div_eq_div_iff (by positivity) (by norm_num)]
        ring] at this
      exact this
    have h4 := mul_le_mul_of_nonneg_right h3
      (by positivity : (0 : ℝ) ≤ 8 * n * r ^ 2 * σ₀ ^ 2)
    have h5 : (0 : ℝ) ≤ n * r ^ 2 * σ₀ ^ 2 := by positivity
    linarith
  · intro f σ hσ hf1 hs1
    have hσpos : 0 < σ := hσmin.trans_le hσ.1
    have hσup : σ < σ₀ - 3 / 8 * r := by
      have := (abs_le.1 hs1).2
      linarith
    obtain ⟨U, hUdef⟩ : ∃ U, ∑ i, (ev n i f - ev n i f₁) ^ 2 = U := ⟨_, rfl⟩
    have hUle : U ≤ n * (r / 8) ^ 2 := by
      rw [← hUdef]
      have := sum_sq_le_of_empiricalNorm_le ev hn hf1
      simp only [map_sub] at this
      exact this
    have hσsq : σ ^ 2 ≤ (σ₀ - 3 / 8 * r) ^ 2 := pow_le_pow_left₀ hσpos.le hσup.le 2
    have hgap : 5 / 16 * σ₀ * r * n ≤ n * τ - (n * σ ^ 2 + U) := by
      rw [← hτ]
      have h1 : r * r < 2 * σ₀ * r := mul_lt_mul_of_pos_right hr2 hr
      have h2 : (n : ℝ) * σ ^ 2 ≤ n * (σ₀ - 3 / 8 * r) ^ 2 :=
        mul_le_mul_of_nonneg_left hσsq hnR.le
      have h3 : (n : ℝ) * (r * r) ≤ n * (2 * σ₀ * r) := mul_le_mul_of_nonneg_left h1.le hnR.le
      linarith
    have hlt : n * σ ^ 2 + U < n * τ := by
      have : 0 < 5 / 16 * σ₀ * r * n := by positivity
      linarith
    have hsub : {y : Fin n → ℝ | ∑ i, (y i - ev n i f₁) ^ 2 ≤ n * τ}ᶜ
        ⊆ {y | n * τ ≤ ∑ i, (y i - ev n i f₁) ^ 2} := by
      intro y hy
      simp only [Set.mem_compl_iff, Set.mem_setOf_eq, not_le] at hy ⊢
      exact hy.le
    have h := pi_gaussian_sq_upper_tail (fun i => ev n i f) (fun i => ev n i f₁)
      (v := σ ^ 2) (τ := τ) (by positivity) (by simp only [hUdef]; exact hlt)
    simp only [hUdef] at h
    refine (measure_mono hsub).trans (h.trans (ofReal_exp_neg_le ?_))
    rw [le_div_iff₀ (by positivity)]
    have hsq : (5 / 16 * σ₀ * r * n) ^ 2 ≤ (n * τ - (n * σ ^ 2 + U)) ^ 2 :=
      pow_le_pow_left₀ (by positivity) hgap 2
    have hσ2 : σ ^ 2 ≤ σmax ^ 2 := pow_le_pow_left₀ hσpos.le hσ.2 2
    have hden := mul_le_mul hσ2 (mul_le_mul_of_nonneg_left hτle hnR.le) (by positivity)
      (by positivity)
    have hK2 : K * σmax ^ 2 ≤ 1 / 1024 := by
      have := mul_le_mul_of_nonneg_right hK1 (by positivity : (0 : ℝ) ≤ σmax ^ 2)
      rw [show 1 / (1024 * σmax ^ 2) * σmax ^ 2 = 1 / 1024 by
        rw [div_mul_eq_mul_div, one_mul, div_eq_div_iff (by positivity) (by norm_num)]
        ring] at this
      exact this
    have h1 := mul_le_mul_of_nonneg_left hden (by positivity : (0 : ℝ) ≤ 16 * K * n * r ^ 2)
    have h2 := mul_le_mul_of_nonneg_right hK2
      (by positivity : (0 : ℝ) ≤ 16 * (n : ℝ) ^ 2 * r ^ 2 * σ₀ ^ 2)
    have h3 : (0 : ℝ) ≤ (n : ℝ) ^ 2 * r ^ 2 * σ₀ ^ 2 := by positivity
    linarith

omit [MeasurableSpace E] [BorelSpace E] in
/-- **UV-5.** Exponentially powerful local tests for `‖f − f′‖ₙ + |σ − σ′|` on a compact noise
range: if `θ₁` is at distance more than `r` from the truth, one test separates the truth from the
whole `r/8`-ball around `θ₁`, with both errors at most `exp(−K n r²)`,
`K = σ_min²/(1024 σ_max⁴)` (`testConst`). The test is `1_R` for the rejection set of the case
that applies: `tests_mean_case` when `‖f₁ − f*‖ₙ > r/2`, else `tests_scale_up_case` or
`tests_scale_down_case` as `σ₁ > σ₀ + r/2` or `σ₁ < σ₀ − r/2`. Each scale test is centred where
its means do the least harm, which keeps the gap linear in `r` on the whole compact range. -/
theorem exists_tests (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) {σmin σmax : ℝ} (hσmin : 0 < σmin)
    (f₀ : E) {σ₀ : ℝ} (hσ₀ : σ₀ ∈ Icc σmin σmax) :
    ∃ K ξ : ℝ, 0 < K ∧ 0 < ξ ∧ ∀ n : ℕ, ∀ r > 0,
      ∀ θ₁ ∈ (univ ×ˢ Icc σmin σmax : Set (E × ℝ)), r < paramDist ev n θ₁ (f₀, σ₀) →
        ∃ φ : (Fin n → ℝ) → ℝ, Measurable φ ∧ (∀ y, 0 ≤ φ y ∧ φ y ≤ 1) ∧
          ∫ y, φ y ∂dataLaw ev n (f₀, σ₀) ≤ Real.exp (-(K * n * r ^ 2)) ∧
          ∀ θ ∈ (univ ×ˢ Icc σmin σmax : Set (E × ℝ)), paramDist ev n θ θ₁ ≤ ξ * r →
            ∫ y, (1 - φ y) ∂dataLaw ev n θ ≤ Real.exp (-(K * n * r ^ 2)) := by
  have hσmax : 0 < σmax := hσmin.trans_le (hσ₀.1.trans hσ₀.2)
  refine ⟨testConst σmin σmax, 1 / 8, testConst_pos hσmin hσmax, by norm_num, ?_⟩
  intro n r hr θ₁ hθ₁ hd
  obtain ⟨f₁, σ₁⟩ := θ₁
  have hσ₁ : σ₁ ∈ Icc σmin σmax := hθ₁.2
  have hd' : r < empiricalNorm ev n (f₁ - f₀) + |σ₁ - σ₀| := hd
  have hexp_pos : 0 ≤ Real.exp (-(testConst σmin σmax * n * r ^ 2)) := (Real.exp_pos _).le
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · refine ⟨fun _ => 0, measurable_const, fun _ => ⟨le_rfl, zero_le_one⟩, ?_, ?_⟩
    · simp
    · intro θ _ _
      simp
  -- a rejection set with both bounds gives the test 1_R
  have key : ∀ R : Set (Fin n → ℝ), MeasurableSet R →
      dataLaw ev n (f₀, σ₀) R
        ≤ ENNReal.ofReal (Real.exp (-(testConst σmin σmax * n * r ^ 2))) →
      (∀ f σ, σ ∈ Icc σmin σmax → paramDist ev n (f, σ) (f₁, σ₁) ≤ 1 / 8 * r →
        dataLaw ev n (f, σ) Rᶜ
          ≤ ENNReal.ofReal (Real.exp (-(testConst σmin σmax * n * r ^ 2)))) →
      ∃ φ : (Fin n → ℝ) → ℝ, Measurable φ ∧ (∀ y, 0 ≤ φ y ∧ φ y ≤ 1) ∧
        ∫ y, φ y ∂dataLaw ev n (f₀, σ₀) ≤ Real.exp (-(testConst σmin σmax * n * r ^ 2)) ∧
        ∀ θ ∈ (univ ×ˢ Icc σmin σmax : Set (E × ℝ)), paramDist ev n θ (f₁, σ₁) ≤ 1 / 8 * r →
          ∫ y, (1 - φ y) ∂dataLaw ev n θ
            ≤ Real.exp (-(testConst σmin σmax * n * r ^ 2)) := by
    intro R hR h0 h1
    refine ⟨R.indicator 1, measurable_one.indicator hR, fun y => ?_,
      integral_indicator_le _ hR hexp_pos h0, ?_⟩
    · by_cases hy : y ∈ R
      · simp [hy]
      · simp [hy]
    · rintro ⟨f, σ⟩ hθ hθb
      exact integral_one_sub_indicator_le _ hR hexp_pos (h1 f σ hθ.2 hθb)
  have hball : ∀ f σ, paramDist ev n (f, σ) (f₁, σ₁) ≤ 1 / 8 * r →
      empiricalNorm ev n (f - f₁) ≤ r / 8 ∧ |σ - σ₁| ≤ r / 8 := by
    intro f σ h
    have h' : empiricalNorm ev n (f - f₁) + |σ - σ₁| ≤ 1 / 8 * r := h
    have := empiricalNorm_nonneg ev n (f - f₁)
    have := abs_nonneg (σ - σ₁)
    constructor <;> linarith
  by_cases hA : r / 2 < empiricalNorm ev n (f₁ - f₀)
  · obtain ⟨R, hR, h0, h1⟩ := tests_mean_case ev hσmin hσ₀ hr hn f₀ f₁ hA
    exact key R hR h0 fun f σ hσ hb => h1 f σ hσ (hball f σ hb).1
  · have hσgap : r / 2 < |σ₁ - σ₀| := by linarith [not_lt.mp hA]
    rcases lt_abs.1 hσgap with hB1 | hB2
    · obtain ⟨R, hR, h0, h1⟩ := tests_scale_up_case ev hσmin hσ₀ hσ₁.2 hr hn f₀ (by linarith)
      exact key R hR h0 fun f σ hσ hb => h1 f σ hσ (hball f σ hb).2
    · obtain ⟨R, hR, h0, h1⟩ :=
        tests_scale_down_case ev hσmin hσ₀ hσ₁.1 hr hn f₀ f₁ (by linarith)
      exact key R hR h0 fun f σ hσ hb => h1 f σ hσ (hball f σ hb).1 (hball f σ hb).2

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
    ∀ᶠ n : ℕ in atTop,
      closedBall f₀ (priorMassRadius C₀ a k n) ×ˢ closedBall σ₀ (priorMassRadius C₀ a k n)
          ⊆ klBall ev n (f₀, σ₀) (priorMassRadius C₀ a k n / (σ₀ / 6)) ∧
      ENNReal.ofReal (Real.exp (-((C₁ + δ) * (σ₀ / 6) ^ 2 * n
          * (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2)))
        ≤ amplitudeNoisePrior πs μW πη
          (closedBall f₀ (priorMassRadius C₀ a k n) ×ˢ closedBall σ₀ (priorMassRadius C₀ a k n))
      := by
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
  refine ⟨?_, ?_⟩
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
* `h59`: (59) for `f*/s`, uniformly over the box (UV-2; confirmed against Part I);
* `h60`: (60), the sieve `Bₙ`'s remaining GP mass, with exponent `C₂ > C₁ + 144/σ₀²`;
* `h61`: (61), the sieve's sup-norm entropy at any fixed multiple of `ε̄ₙ`;
* `hRadius`: the sieve's sup-norm radius `Rₙ ≤ n^p` (UV-4; confirmed against Parts II–III).

Then for some `M`, the expected posterior mass of `{‖f − f*‖ₙ + |σ − σ₀| ≥ M ε̄ₙ}` tends to `0`:
the known-σ rate. The proof is `fixed_constant_contraction` fed by UV-1 + UV-2
(`eventually_prior_mass`), UV-3 (`prior_compl_sieve_le`), UV-4 + UV-4a (`eventually_sieve_cover`)
and UV-5 (`exists_tests`). `posterior_contraction_unknown_variance_gvdv07` is the same composition
on the printed general theorem. -/
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
  obtain ⟨M, hM⟩ := fixed_constant_contraction ev (amplitudeNoisePrior πs μW πη) hpos (f₀, σ₀)
    hσ₀pos
    (paramDist ev) (paramDist_comm ev) (paramDist_triangle ev)
    (fun n θ' => (continuous_paramDist ev n θ').measurable)
    (fun n => priorMassRadius C₀ a k n / (σ₀ / 6)) (fun n => L * priorMassRadius Cb a kb n)
    hrate hnε (fun n => closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)))
    (fun n => isClosed_closure.measurableSet) ((C₁ + δ) * (σ₀ / 6) ^ 2) K ξ hK hξ htests'
    (eventually_sieve_cover ev hev B hsmin.le hsmax hCb ha ha' hkb hξ hL1 h61' hR)
    ((eventually_prior_mass ev hev f₀ μW πs hsmin hs hπs_pos πη hηlo hη hπη hσ₀ hC₀ ha ha' hk hδ
      h59).mono fun n h => ⟨_, h.1, measurableSet_closedBall.prod measurableSet_closedBall, h.2⟩)
    hrem
  refine ⟨M * L, hM.congr fun n => ?_⟩
  simp only [mul_assoc]

/-! ### Plumbing for the composition from the printed theorem -/

/-- The likelihood is jointly measurable in the data and the parameter. -/
theorem measurable_likelihood (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (n : ℕ) :
    Measurable fun p : (Fin n → ℝ) × (E × ℝ) => likelihood ev n p.2 p.1 := by
  unfold likelihood
  refine Finset.measurable_prod _ fun i _ => ?_
  exact measurable_uncurry_gaussianPDF.comp
    (((ev n i).continuous.measurable.comp (measurable_fst.comp measurable_snd)).prodMk
      ((measurable_real_toNNReal.comp ((measurable_snd.comp measurable_snd).pow_const 2)).prodMk
        ((measurable_pi_apply i).comp measurable_fst)))

/-- The posterior mass of a set is a measurable function of the data. -/
theorem measurable_posterior (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (prior : Measure (E × ℝ))
    [IsFiniteMeasure prior] (n : ℕ) (A : Set (E × ℝ)) :
    Measurable fun y => posterior ev prior n y A := by
  have h1 : Measurable fun y => ∫⁻ θ in A, likelihood ev n θ y ∂prior :=
    (measurable_likelihood ev n).lintegral_prod_right' (ν := prior.restrict A)
  have h2 : Measurable fun y => ∫⁻ θ, likelihood ev n θ y ∂prior :=
    (measurable_likelihood ev n).lintegral_prod_right' (ν := prior)
  unfold posterior
  simp only [div_eq_mul_inv]
  exact h1.mul h2.inv

omit [BorelSpace E] in
/-- The posterior is subadditive across a sieve and its complement. -/
theorem posterior_le_add (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (prior : Measure (E × ℝ)) (n : ℕ)
    (y : Fin n → ℝ) (A Θ : Set (E × ℝ)) :
    posterior ev prior n y A ≤ posterior ev prior n y (A ∩ Θ) + posterior ev prior n y Θᶜ := by
  unfold posterior
  rw [← ENNReal.add_div]
  apply ENNReal.div_le_div_right
  calc ∫⁻ θ in A, likelihood ev n θ y ∂prior
      ≤ ∫⁻ θ in (A ∩ Θ) ∪ Θᶜ, likelihood ev n θ y ∂prior := by
        refine lintegral_mono_set fun θ hθ => ?_
        by_cases h : θ ∈ Θ
        · exact Or.inl ⟨hθ, h⟩
        · exact Or.inr h
    _ ≤ _ := lintegral_union_le _ _ _

omit [MeasurableSpace E] [BorelSpace E] in
/-- The KL neighbourhood grows with its radius. -/
theorem klBall_mono (ev : (n : ℕ) → Fin n → E →L[ℝ] ℝ) (n : ℕ) (θ₀ : E × ℝ) {ε ε' : ℝ}
    (h0 : 0 ≤ ε) (h : ε ≤ ε') : klBall ev n θ₀ ε ⊆ klBall ev n θ₀ ε' := by
  rintro θ ⟨h1, h2⟩
  have hle : ENNReal.ofReal (n * ε ^ 2) ≤ ENNReal.ofReal (n * ε' ^ 2) :=
    ENNReal.ofReal_le_ofReal (mul_le_mul_of_nonneg_left (pow_le_pow_left₀ h0 h 2)
      (Nat.cast_nonneg n))
  exact ⟨h1.trans hle, h2.trans hle⟩

/-- Under [P1], the noise standard deviation stays in `[σ_min, σ_max]` almost surely. -/
theorem prior_compl_noiseRange (πs : Measure ℝ) (μW : Measure E) (πη : Measure ℝ)
    [IsProbabilityMeasure πs] [IsProbabilityMeasure μW] [IsProbabilityMeasure πη]
    {smin smax ηlo ηhi : ℝ} (hsmin : 0 ≤ smin) (hπs : πs (Icc smin smax)ᶜ = 0)
    (hπη : πη (Icc ηlo ηhi)ᶜ = 0) :
    amplitudeNoisePrior πs μW πη (univ ×ˢ Icc (smin * √ηlo) (smax * √ηhi))ᶜ = 0 := by
  rw [amplitudeNoisePrior, Measure.map_apply measurable_amplitudeNoiseMap
    (MeasurableSet.univ.prod measurableSet_Icc).compl]
  refine measure_mono_null (t := (Icc smin smax)ᶜ ×ˢ univ ∪ univ ×ˢ (univ ×ˢ (Icc ηlo ηhi)ᶜ))
    ?_ (measure_union_null (by simp [Measure.prod_prod, hπs]) (by simp [Measure.prod_prod, hπη]))
  rintro ⟨s, w, η⟩ hp
  simp only [mem_preimage, mem_compl_iff, mem_prod, mem_univ, true_and, amplitudeNoiseMap]
    at hp
  by_contra hcon
  simp only [mem_union, mem_prod, mem_compl_iff, mem_univ, and_true, true_and, not_or,
    not_not] at hcon
  obtain ⟨hs, hη⟩ := hcon
  exact hp ⟨mul_le_mul hs.1 (Real.sqrt_le_sqrt hη.1) (Real.sqrt_nonneg _) (hsmin.trans hs.1),
    mul_le_mul hs.2 (Real.sqrt_le_sqrt hη.2) (Real.sqrt_nonneg _) (hsmin.trans (hs.1.trans hs.2))⟩

/-- The shells condition (2.5) from the prior mass alone: `P(A) ≤ 1 ≤ e^{K n e′² j²/2} e^{−C n e²}`
once `j ≥ ⌈2|C|/K⌉ + 1` and `0 ≤ e ≤ e′`. -/
theorem measure_le_exp_mul_of_mass {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    [IsProbabilityMeasure P] (A B : Set Ω) {K C e e' : ℝ} {n j : ℕ} (hK : 0 < K)
    (hj : ⌈2 * |C| / K⌉₊ + 1 ≤ j) (he : 0 ≤ e) (hee : e ≤ e')
    (hB : ENNReal.ofReal (Real.exp (-(C * n * e ^ 2))) ≤ P B) :
    P A ≤ ENNReal.ofReal (Real.exp (K * n * e' ^ 2 * j ^ 2 / 2)) * P B := by
  have hj1 : ((⌈2 * |C| / K⌉₊ + 1 : ℕ) : ℝ) ≤ j := by exact_mod_cast hj
  push_cast at hj1
  have h2C : 2 * |C| ≤ ⌈2 * |C| / K⌉₊ * K := (div_le_iff₀ hK).1 (Nat.le_ceil _)
  have hKj : 2 * |C| + K ≤ j * K := by
    have := mul_le_mul_of_nonneg_right hj1 hK.le
    linarith
  have hj0 : (1 : ℝ) ≤ j := by
    have := Nat.cast_nonneg (α := ℝ) ⌈2 * |C| / K⌉₊
    linarith
  have hjj : (j : ℝ) ≤ j ^ 2 := by nlinarith
  have hCj : |C| ≤ K * j ^ 2 / 2 := by
    have := mul_le_mul_of_nonneg_left hjj hK.le
    linarith [abs_nonneg C]
  have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have he2 : e ^ 2 ≤ e' ^ 2 := pow_le_pow_left₀ he hee 2
  have hexp : 0 ≤ K * n * e' ^ 2 * j ^ 2 / 2 + -(C * n * e ^ 2) := by
    have h1 : C * (n * e ^ 2) ≤ |C| * (n * e' ^ 2) :=
      (mul_le_mul_of_nonneg_right (le_abs_self C) (by positivity)).trans
        (mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left he2 hn0) (abs_nonneg C))
    have h2 := mul_le_mul_of_nonneg_right hCj (by positivity : (0 : ℝ) ≤ n * e' ^ 2)
    linarith
  calc P A ≤ 1 := prob_le_one
    _ ≤ ENNReal.ofReal (Real.exp (K * n * e' ^ 2 * j ^ 2 / 2))
        * ENNReal.ofReal (Real.exp (-(C * n * e ^ 2))) := by
      rw [← ENNReal.ofReal_mul (Real.exp_pos _).le, ← Real.exp_add]
      exact ENNReal.one_le_ofReal.2 (Real.one_le_exp hexp)
    _ ≤ _ := mul_le_mul_right hB _

/-- The ratio in Lemma 1 from exponential bounds: `a/b · e^w ≤ e^z` when `a ≤ e^u`, `e^v ≤ b` and
`u − v + w ≤ z`. -/
theorem div_mul_ofReal_exp_le {a b : ℝ≥0∞} {u v w z : ℝ} (ha : a ≤ ENNReal.ofReal (Real.exp u))
    (hb : ENNReal.ofReal (Real.exp v) ≤ b) (h : u - v + w ≤ z) :
    a / b * ENNReal.ofReal (Real.exp w) ≤ ENNReal.ofReal (Real.exp z) := by
  calc a / b * ENNReal.ofReal (Real.exp w)
      ≤ ENNReal.ofReal (Real.exp u) / ENNReal.ofReal (Real.exp v)
        * ENNReal.ofReal (Real.exp w) :=
        mul_le_mul' (ENNReal.div_le_div ha hb) le_rfl
    _ = ENNReal.ofReal (Real.exp (u - v + w)) := by
        rw [← ENNReal.ofReal_div_of_pos (Real.exp_pos _),
          ← ENNReal.ofReal_mul (div_nonneg (Real.exp_pos _).le (Real.exp_pos _).le),
          ← Real.exp_sub, ← Real.exp_add]
    _ ≤ ENNReal.ofReal (Real.exp z) := ENNReal.ofReal_le_ofReal (Real.exp_le_exp.2 h)

/-- **UV-6 from the printed general theorem.** The hypotheses of
`posterior_contraction_unknown_variance`, with the remaining-mass margin `C₂ > C₁ + 72/σ₀²` that
Lemma 1 of Ghosal & van der Vaart (2007) asks for, give the printed conclusion: for **every**
`Mₙ → ∞`, the expected posterior mass of `{‖f − f*‖ₙ + |σ − σ₀| ≥ Mₙ ε̄ₙ}` tends to `0`.

The proof applies `gvdv07_theorem1` on the sieve at the entropy rate `L ε̄ₙ` and
`gvdv07_lemma1` off the sieve at the prior-mass rate `6εₙ/σ₀`, and adds the two
(`posterior_le_add`). This is the two-rate composition: Theorem 1 needs only the entropy and the
shells at `ε̄ₙ`, and Lemma 1 only the prior mass at `εₙ`. -/
theorem posterior_contraction_unknown_variance_gvdv07
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
    (hk : 0 ≤ k) (hkk : k ≤ kb) {C₁ C₂ : ℝ} (hC₂ : C₁ + 72 / σ₀ ^ 2 < C₂)
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
    ∀ M : ℕ → ℝ, Tendsto M atTop atTop →
      Tendsto (fun n => ∫⁻ y, posterior ev (amplitudeNoisePrior πs μW πη) n y
        {θ | M n * priorMassRadius Cb a kb n ≤ paramDist ev n θ (f₀, σ₀)}
        ∂dataLaw ev n (f₀, σ₀)) atTop (𝓝 0) := by
  intro M hM
  haveI := isProbabilityMeasure_amplitudeNoisePrior πs μW πη
  have hkb : 0 ≤ kb := hk.trans hkk
  have hσ₀pos : 0 < σ₀ := lt_of_le_of_lt (by positivity) hσ₀.1
  have hsmax : 0 < smax := hsmin.trans hs
  have hc : 0 < σ₀ / 6 := by positivity
  have hσmin : 0 < smin * √ηlo := mul_pos hsmin (Real.sqrt_pos.2 hηlo)
  obtain ⟨K, ξ, hK, hξ, htests⟩ := exists_tests ev hσmin f₀ (σ₀ := σ₀) ⟨hσ₀.1.le, hσ₀.2.le⟩
  obtain ⟨K₀, h61'⟩ := h61 (ξ / 2 / (12 * smax)) (by positivity)
  obtain ⟨p, hR⟩ := hRadius
  obtain ⟨L, hL⟩ : ∃ L : ℝ, 1 + |K₀| + C₀ / (Cb * (σ₀ / 6)) = L := ⟨_, rfl⟩
  have hL1 : 1 + |K₀| ≤ L := by
    have : 0 ≤ C₀ / (Cb * (σ₀ / 6)) := by positivity
    rw [← hL]
    linarith
  have hLpos : 0 < L := by linarith [abs_nonneg K₀]
  obtain ⟨δ, hδ_def⟩ : ∃ δ : ℝ, (C₂ - C₁ - 72 / σ₀ ^ 2) / 2 = δ := ⟨_, rfl⟩
  have hδ : 0 < δ := by rw [← hδ_def]; linarith
  have hC2 : C₂ = C₁ + 72 / σ₀ ^ 2 + 2 * δ := by rw [← hδ_def]; ring
  -- the parameter set, the sieves, the rates
  set S : Set (E × ℝ) := univ ×ˢ Icc (smin * √ηlo) (smax * √ηhi) with hS_def
  have hS : amplitudeNoisePrior πs μW πη Sᶜ = 0 :=
    prior_compl_noiseRange πs μW πη hsmin.le hπs_supp hπη_supp
  have hS_pos : ∀ θ ∈ S, 0 < θ.2 := fun θ hθ => hσmin.trans_le hθ.2.1
  have hθ₀ : (f₀, σ₀) ∈ S := ⟨mem_univ _, hσ₀.1.le, hσ₀.2.le⟩
  have hρpos := eventually_priorMassRadius_pos C₀ a k hC₀
  have hrate : ∀ᶠ n : ℕ in atTop, 0 < priorMassRadius C₀ a k n / (σ₀ / 6) ∧
      priorMassRadius C₀ a k n / (σ₀ / 6) ≤ L * priorMassRadius Cb a kb n := by
    filter_upwards [hρpos, priorMassRadius_le_mul C₀ Cb a k kb hC₀ hCb hkk,
      eventually_priorMassRadius_pos Cb a kb hCb] with n h1 h2 h3
    refine ⟨div_pos h1 hc, ?_⟩
    rw [div_le_iff₀ hc]
    have hLc : C₀ / (Cb * (σ₀ / 6)) ≤ L := by rw [← hL]; linarith [abs_nonneg K₀]
    calc priorMassRadius C₀ a k n ≤ C₀ / Cb * priorMassRadius Cb a kb n := h2
      _ = C₀ / (Cb * (σ₀ / 6)) * priorMassRadius Cb a kb n * (σ₀ / 6) := by field_simp
      _ ≤ L * priorMassRadius Cb a kb n * (σ₀ / 6) :=
          mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hLc h3.le) hc.le
  have hnρ := tendsto_n_priorMassRadius_sq C₀ a k hC₀ ha ha' hk
  have hnεb := tendsto_n_priorMassRadius_sq Cb a kb hCb ha ha' hkb
  have hsieve : ∀ n, closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)) ⊆ S :=
    fun n => closure_minimal (prod_mono (subset_univ _) subset_rfl)
      (isClosed_univ.prod isClosed_Icc)
  -- Theorem 1, on the sieve, at the rate L ε̄ₙ
  have hT1 := gvdv07_theorem1 ev (amplitudeNoisePrior πs μW πη) S hS hS_pos (f₀, σ₀) hθ₀
    (paramDist ev) (paramDist_comm ev) (paramDist_triangle ev)
    (fun n θ' => (continuous_paramDist ev n θ').measurable)
    (fun n => L * priorMassRadius Cb a kb n)
    ((eventually_priorMassRadius_pos Cb a kb hCb).mono fun n h => mul_pos hLpos h)
    (by simpa using (tendsto_priorMassRadius_zero Cb a kb ha).const_mul L)
    ⟨1, one_pos, (hnεb.eventually_ge_atTop 1).mono fun n h => by
      have hL2 : 1 ≤ L ^ 2 := one_le_pow₀ (by linarith [abs_nonneg K₀])
      have := mul_le_mul hL2 h zero_le_one (by positivity)
      calc (1 : ℝ) = 1 * 1 := by ring
        _ ≤ L ^ 2 * (n * priorMassRadius Cb a kb n ^ 2) := this
        _ = n * (L * priorMassRadius Cb a kb n) ^ 2 := by ring⟩
    (fun n => closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)))
    (fun n => isClosed_closure.measurableSet) hsieve K ξ hK hξ htests
    ?hentropy ?hshell
  case hentropy =>
    filter_upwards [eventually_sieve_cover ev hev B (σmin := smin * √ηlo) (σmax := smax * √ηhi)
      hsmin.le hsmax hCb ha ha' hkb (by positivity : 0 < ξ / 2) hL1 h61' hR]
      with n ⟨N, hcard, hcov⟩
    intro r hr
    refine ⟨N, hcard, fun θ hθ => ?_⟩
    obtain ⟨θ', hθ'N, hθθ'⟩ := mem_iUnion₂.1 (hcov hθ.1)
    have hθθ'' : paramDist ev n θ θ' ≤ ξ / 2 * (L * priorMassRadius Cb a kb n) / 2 := hθθ'
    have := mul_lt_mul_of_pos_left hr hξ
    exact mem_iUnion₂.2 ⟨θ', hθ'N,
      (le_trans hθθ'' (by linarith) : paramDist ev n θ θ' ≤ ξ * r / 4)⟩
  case hshell =>
    obtain ⟨C, hC⟩ : ∃ C : ℝ, (C₁ + δ) * (σ₀ / 6) ^ 2 = C := ⟨_, rfl⟩
    refine ⟨⌈2 * |C| / K⌉₊ + 1, ?_⟩
    filter_upwards [eventually_prior_mass ev hev f₀ μW πs hsmin hs hπs_pos πη hηlo hη hπη hσ₀
      hC₀ ha ha' hk hδ h59, hrate] with n ⟨hBsub, hBmass⟩ ⟨hpos, hle⟩
    refine ⟨_, hBsub.trans (klBall_mono ev n (f₀, σ₀) hpos.le hle),
      measurableSet_closedBall.prod measurableSet_closedBall, fun j hj => ?_⟩
    rw [hC] at hBmass
    exact measure_le_exp_mul_of_mass _ _ _ hK hj hpos.le hle hBmass
  -- Lemma 1, off the sieve, at the rate 6εₙ/σ₀
  have hL1' := gvdv07_lemma1 ev (amplitudeNoisePrior πs μW πη) S hS hS_pos (f₀, σ₀) hθ₀
    (fun n => priorMassRadius C₀ a k n / (σ₀ / 6)) (hρpos.mono fun n h => div_pos h hc)
    ⟨1, one_pos, (((hnρ.const_mul_atTop (inv_pos.2 (pow_pos hc 2))).eventually_ge_atTop 1)).mono
      fun n h => by
        have e : ((σ₀ / 6) ^ 2)⁻¹ * (n * priorMassRadius C₀ a k n ^ 2)
            = n * (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2 := by rw [div_pow]; ring
        rw [← e]
        exact h⟩
    (fun n => closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)))
    (fun n => isClosed_closure.measurableSet)
    (fun n => closedBall f₀ (priorMassRadius C₀ a k n) ×ˢ closedBall σ₀ (priorMassRadius C₀ a k n))
    ?hB ?hratio
  case hB =>
    filter_upwards [eventually_prior_mass ev hev f₀ μW πs hsmin hs hπs_pos πη hηlo hη hπη hσ₀
      hC₀ ha ha' hk hδ h59] with n ⟨hBsub, hBmass⟩
    exact ⟨hBsub, measurableSet_closedBall.prod measurableSet_closedBall,
      lt_of_lt_of_le (ENNReal.ofReal_pos.2 (Real.exp_pos _)) hBmass⟩
  case hratio =>
    have hlim : Tendsto (fun n : ℕ => ENNReal.ofReal
        (Real.exp (-(δ * (n * priorMassRadius C₀ a k n ^ 2))))) atTop (𝓝 0) := by
      rw [← ENNReal.ofReal_zero]
      exact ENNReal.tendsto_ofReal (Real.tendsto_exp_atBot.comp
        (tendsto_neg_atTop_atBot.comp (hnρ.const_mul_atTop hδ)))
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hlim
      (Eventually.of_forall fun _ => zero_le) ?_
    filter_upwards [eventually_prior_mass ev hev f₀ μW πs hsmin hs hπs_pos πη hηlo hη hπη hσ₀
      hC₀ ha ha' hk hδ h59, h60] with n ⟨_, hBmass⟩ h60n
    have hrem := (prior_compl_sieve_le πs μW πη hsmin.le hπs_supp hπη_supp (B n)).trans h60n
    refine div_mul_ofReal_exp_le hrem hBmass ?_
    have hc0 : (σ₀ / 6) ^ 2 ≠ 0 := pow_ne_zero 2 hc.ne'
    have hcan : (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2 * (σ₀ / 6) ^ 2
        = priorMassRadius C₀ a k n ^ 2 := by rw [div_pow, div_mul_cancel₀ _ hc0]
    have hcan2 : (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2
        = 36 * priorMassRadius C₀ a k n ^ 2 / σ₀ ^ 2 := by
      rw [div_pow, div_eq_div_iff hc0 (pow_ne_zero 2 hσ₀pos.ne')]
      ring
    have e1 : (C₁ + δ) * (σ₀ / 6) ^ 2 * n * (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2
        = (C₁ + δ) * (n * priorMassRadius C₀ a k n ^ 2) := by
      rw [show (C₁ + δ) * (σ₀ / 6) ^ 2 * n * (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2
          = (C₁ + δ) * n * ((priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2 * (σ₀ / 6) ^ 2) by ring,
        hcan]
      ring
    have e2 : 2 * n * (priorMassRadius C₀ a k n / (σ₀ / 6)) ^ 2
        = 72 / σ₀ ^ 2 * (n * priorMassRadius C₀ a k n ^ 2) := by
      rw [hcan2]
      ring
    rw [e1, e2, hC2]
    exact le_of_eq (by ring)
  -- add the two
  have hM' : Tendsto (fun n => M n / L) atTop atTop := hM.atTop_div_const hLpos
  have h1 := hT1 (fun n => M n / L) hM'
  have h12 := h1.add hL1'
  rw [add_zero] at h12
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds h12
    (fun _ => zero_le) (fun n => ?_)
  have hset : {θ | θ ∈ closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)) ∧
        M n / L * (L * priorMassRadius Cb a kb n) ≤ paramDist ev n θ (f₀, σ₀)}
      = {θ | M n * priorMassRadius Cb a kb n ≤ paramDist ev n θ (f₀, σ₀)}
        ∩ closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)) := by
    have e : M n / L * (L * priorMassRadius Cb a kb n) = M n * priorMassRadius Cb a kb n := by
      field_simp
    ext θ
    simp only [mem_setOf_eq, mem_inter_iff, e]
    exact and_comm
  dsimp only
  rw [hset]
  calc ∫⁻ y, posterior ev (amplitudeNoisePrior πs μW πη) n y
        {θ | M n * priorMassRadius Cb a kb n ≤ paramDist ev n θ (f₀, σ₀)} ∂dataLaw ev n (f₀, σ₀)
      ≤ ∫⁻ y, (posterior ev (amplitudeNoisePrior πs μW πη) n y
          ({θ | M n * priorMassRadius Cb a kb n ≤ paramDist ev n θ (f₀, σ₀)}
            ∩ closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)))
          + posterior ev (amplitudeNoisePrior πs μW πη) n y
            (closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)))ᶜ)
          ∂dataLaw ev n (f₀, σ₀) :=
        lintegral_mono fun y => posterior_le_add ev _ n y _ _
    _ = _ := lintegral_add_right _ (measurable_posterior ev _ n _)

/-- **UV-6 with the paper's tunable sieve.** TWCD2025's construction (Appendix A.1, Part II, with
the constants of (80) scaled) gives a sieve for every remaining-mass exponent `C₂`, with (61) at
the same rate and the same polynomial radius. In that form the margin `hC₂` needs no checking:
pick `C₂ = C₁ + 144/σ₀² + 1`. -/
theorem posterior_contraction_unknown_variance_of_tunable_sieve
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
    (hk : 0 ≤ k) (hkk : k ≤ kb) {C₁ : ℝ}
    (h59 : ∀ᶠ n : ℕ in atTop, ∀ s ∈ Icc smin smax,
      ENNReal.ofReal (Real.exp (-(C₁ * n * priorMassRadius C₀ a k n ^ 2)))
        ≤ μW (closedBall (s⁻¹ • f₀) (priorMassRadius C₀ a k n / s)))
    (hsieve : ∀ C₂ : ℝ, ∃ B : ℕ → Set E,
      (∀ᶠ n : ℕ in atTop,
        μW (B n)ᶜ ≤ ENNReal.ofReal (Real.exp (-(C₂ * n * priorMassRadius C₀ a k n ^ 2)))) ∧
      (∀ c > 0, ∃ K₀ : ℝ, ∀ᶠ n : ℕ in atTop, ∃ N : ℕ,
        externalCoveringNumber (c * priorMassRadius Cb a kb n).toNNReal (B n) ≤ N ∧
          Real.log N ≤ K₀ * n * priorMassRadius Cb a kb n ^ 2) ∧
      (∃ p : ℕ, ∀ᶠ n : ℕ in atTop, ∀ g ∈ B n, ‖g‖ ≤ (n : ℝ) ^ p)) :
    ∃ M : ℝ, Tendsto (fun n => ∫⁻ y, posterior ev (amplitudeNoisePrior πs μW πη) n y
        {θ | M * priorMassRadius Cb a kb n ≤ paramDist ev n θ (f₀, σ₀)}
        ∂dataLaw ev n (f₀, σ₀)) atTop (𝓝 0) := by
  obtain ⟨B, h60, h61, hRadius⟩ := hsieve (C₁ + 144 / σ₀ ^ 2 + 1)
  exact posterior_contraction_unknown_variance ev hev f₀ μW πs hsmin hs hπs_supp hπs_pos πη hηlo
    hη hπη_supp hπη hσ₀ hC₀ hCb ha ha' hk hkk (C₂ := C₁ + 144 / σ₀ ^ 2 + 1) (by linarith) h59 B
    h60 h61 hRadius

end BGP
