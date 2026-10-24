import Mathlib.Analysis.Normed.Module.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.MeasureTheory.Constructions.BorelSpace.Basic
import Mathlib.MeasureTheory.Measure.Map
import Mathlib.Topology.MetricSpace.Pseudo.Lemmas
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# BGP-4 — bookkeeping lemmas for contraction with an unknown noise variance

The load-bearing elementary pieces of `Protolemmata/BARS/bgp-unknown-variance-lemmas.md`, which extends
TWCD2025's contraction theorem to an unknown noise variance:

* **UV-1a** `gaussian_kl_scale_part_le` — the scale part of the Gaussian KL divergence,
  `−log r + (r² − 1)/2` with `r = σ₀/σ`, is at most `(5/2)(r − 1)²` for `r ≥ 1/2`.
* **UV-1b** `gaussian_v_scale_part_le` — the scale part of the log-ratio variance, `(r² − 1)²/2`, is at most
  `(9/2)(r − 1)²` on `[1/2, 2]`.
* **UV-2a** `map_smul_closedBall` — a ball's prior mass under an amplitude rescaling `w ↦ s • w` is the mass of the
  rescaled ball: the device the authors' Remark A.1 uses, carried to `f*/s`.
* **UV-2b** `log_inv_radius_le_n_radius_sq` — the noise interval's prior mass is absorbed: for the paper's radius
  shape `εₙ = C n^(−a) (log n)^k` with `0 < a < 1/2` and `k ≥ 0`, `log(1/εₙ) ≤ n εₙ²` eventually.
-/

open MeasureTheory Filter

namespace BGP

/-! ### UV-1: the Gaussian KL divergence and log-ratio variance, scale parts -/

/-- **UV-1a.** The scale part of `KL(N(μ₀, σ₀²) ‖ N(μ, σ²))`, written in `r = σ₀/σ`, is locally quadratic:
`−log r + (r² − 1)/2 ≤ (5/2)(r − 1)²` for every `r ≥ 1/2`. From `1 − 1/r ≤ log r`, the gap is
`(r − 1)²(2r − 1)/r ≥ 0`. (The protolemma states it on `[1/2, 2]`; the upper end is not needed.) -/
theorem gaussian_kl_scale_part_le (r : ℝ) (hr : 1 / 2 ≤ r) :
    -Real.log r + (r ^ 2 - 1) / 2 ≤ (5 / 2) * (r - 1) ^ 2 := by
  have hr0 : 0 < r := by linarith
  have hr0' : r ≠ 0 := hr0.ne'
  have hlog : 1 - r⁻¹ ≤ Real.log r := Real.one_sub_inv_le_log_of_pos hr0
  have gap : (5 / 2) * (r - 1) ^ 2 - (r⁻¹ - 1 + (r ^ 2 - 1) / 2) = (r - 1) ^ 2 * (2 * r - 1) / r := by
    field_simp
    ring
  have hgap : 0 ≤ (r - 1) ^ 2 * (2 * r - 1) / r :=
    div_nonneg (mul_nonneg (sq_nonneg _) (by linarith)) hr0.le
  linarith

/-- **UV-1b.** The scale part of the log-ratio variance, `(r² − 1)²/2`, is at most `(9/2)(r − 1)²` on `[1/2, 2]`,
since `(r² − 1)² = (r − 1)²(r + 1)²` and `(r + 1)² ≤ 9` there. -/
theorem gaussian_v_scale_part_le (r : ℝ) (hr : r ∈ Set.Icc (1 / 2 : ℝ) 2) :
    (r ^ 2 - 1) ^ 2 / 2 ≤ (9 / 2) * (r - 1) ^ 2 := by
  obtain ⟨hlo, hhi⟩ := hr
  have h1 : (r + 1) ^ 2 ≤ 9 := by
    nlinarith [mul_nonneg (sub_nonneg.2 hhi) (by linarith : (0 : ℝ) ≤ r + 4)]
  have h2 : (r ^ 2 - 1) ^ 2 = (r - 1) ^ 2 * (r + 1) ^ 2 := by ring
  have h3 : (r - 1) ^ 2 * (r + 1) ^ 2 ≤ (r - 1) ^ 2 * 9 := mul_le_mul_of_nonneg_left h1 (sq_nonneg _)
  rw [h2]
  linarith

/-! ### UV-2: prior mass -/

/-- **UV-2a.** Under an amplitude rescaling `w ↦ s • w` with `s > 0`, a closed ball's mass is the mass of the
rescaled ball: `μ.map (s • ·) (closedBall f ε) = μ (closedBall (s⁻¹ • f) (ε / s))`. The preimage of
`closedBall f ε` is `closedBall (s⁻¹ • f) (ε / s)`, because `s • w − f = s • (w − s⁻¹ • f)`. -/
theorem map_smul_closedBall {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [MeasurableSpace E]
    [BorelSpace E] (μ : Measure E) {s : ℝ} (hs : 0 < s) (f : E) (ε : ℝ) :
    (μ.map (fun w => s • w)) (Metric.closedBall f ε) = μ (Metric.closedBall (s⁻¹ • f) (ε / s)) := by
  rw [Measure.map_apply (continuous_const_smul s).measurable Metric.isClosed_closedBall.measurableSet]
  congr 1
  ext w
  simp only [Set.mem_preimage, Metric.mem_closedBall, dist_eq_norm]
  have hsw : s • w - f = s • (w - s⁻¹ • f) := by
    rw [smul_sub, smul_smul, mul_inv_cancel₀ hs.ne', one_smul]
  rw [hsw, norm_smul, Real.norm_of_nonneg hs.le, le_div_iff₀ hs, mul_comm]

/-- The paper's prior-mass radius shape, `εₙ = C n^(−a) (log n)^k` (TWCD2025, (62)). -/
noncomputable def priorMassRadius (C a k : ℝ) (n : ℕ) : ℝ :=
  C * (n : ℝ) ^ (-a) * Real.log n ^ k

/-- **UV-2b.** The noise interval's prior mass `≍ εₙ` is absorbed into `exp(−n εₙ²)`: for `0 < a < 1/2` and
`k ≥ 0`, `log(1/εₙ) ≤ n εₙ²` eventually. With `b = (1 − 2a)/2`, `y = nᵇ` and `L = log n ≥ 1`:
`log(1/εₙ) ≤ −log C + a·L`, `L ≤ y/b` (from `log y ≤ y − 1`), and `n εₙ² ≥ C² y²`; the quadratic wins once
`y ≥ max 1 ((|log C| + a/b)/C²)`. -/
theorem log_inv_radius_le_n_radius_sq (C a k : ℝ) (hC : 0 < C) (ha : 0 < a) (ha' : a < 1 / 2) (hk : 0 ≤ k) :
    ∀ᶠ n : ℕ in atTop, Real.log (1 / priorMassRadius C a k n) ≤ n * priorMassRadius C a k n ^ 2 := by
  set b : ℝ := (1 - 2 * a) / 2 with hb_def
  have hb : 0 < b := by rw [hb_def]; linarith
  set M : ℝ := max 1 ((|Real.log C| + a / b) / C ^ 2) with hM_def
  have hM1 : 1 ≤ M := le_max_left _ _
  refine Filter.eventually_atTop.2 ⟨max ⌈Real.exp 1⌉₊ ⌈M ^ (1 / b)⌉₊, fun n hn => ?_⟩
  have hn1 : ⌈Real.exp 1⌉₊ ≤ n := le_of_max_le_left hn
  have hn2 : ⌈M ^ (1 / b)⌉₊ ≤ n := le_of_max_le_right hn
  set x : ℝ := (n : ℝ) with hx_def
  have hxe : Real.exp 1 ≤ x := (Nat.le_ceil _).trans (by rw [hx_def]; exact_mod_cast hn1)
  have hxM : M ^ (1 / b) ≤ x := (Nat.le_ceil _).trans (by rw [hx_def]; exact_mod_cast hn2)
  have hx0 : 0 < x := lt_of_lt_of_le (Real.exp_pos 1) hxe
  set L : ℝ := Real.log x with hL_def
  have hL1 : 1 ≤ L := by
    rw [hL_def, ← Real.log_exp 1]
    exact Real.log_le_log (Real.exp_pos 1) hxe
  have hL0 : 0 < L := lt_of_lt_of_le one_pos hL1
  have hxa : 0 < x ^ (-a) := Real.rpow_pos_of_pos hx0 _
  have hLk : 1 ≤ L ^ k := Real.one_le_rpow hL1 hk
  have hLk0 : 0 < L ^ k := lt_of_lt_of_le one_pos hLk
  have hε : priorMassRadius C a k n = C * x ^ (-a) * L ^ k := rfl
  -- log(1/εₙ) ≤ −log C + a·L
  have hlog : Real.log (1 / priorMassRadius C a k n) ≤ -Real.log C + a * L := by
    rw [hε, one_div, Real.log_inv, Real.log_mul (mul_pos hC hxa).ne' hLk0.ne', Real.log_mul hC.ne' hxa.ne',
      Real.log_rpow hx0, Real.log_rpow hL0]
    have : 0 ≤ k * Real.log L := mul_nonneg hk (Real.log_nonneg hL1)
    linarith
  -- y = xᵇ is past the threshold M
  set y : ℝ := x ^ b with hy_def
  have hy0 : 0 < y := Real.rpow_pos_of_pos hx0 _
  have hyM : M ≤ y := by
    have h1 : (M ^ (1 / b)) ^ b ≤ x ^ b := Real.rpow_le_rpow (Real.rpow_nonneg (by linarith) _) hxM hb.le
    have h2 : (M ^ (1 / b)) ^ b = M := by
      rw [← Real.rpow_mul (by linarith), one_div, inv_mul_cancel₀ hb.ne', Real.rpow_one]
    rw [h2] at h1
    exact h1
  have hy1 : 1 ≤ y := hM1.trans hyM
  -- L ≤ y/b, from b·L = log y ≤ y − 1
  have hLy : L ≤ y / b := by
    have h1 : Real.log y = b * L := by rw [hy_def, Real.log_rpow hx0]
    have h2 : Real.log y ≤ y - 1 := Real.log_le_sub_one_of_pos hy0
    rw [le_div_iff₀ hb]
    linarith
  -- −log C + a·L ≤ (|log C| + a/b)·y
  have hstep : -Real.log C + a * L ≤ (|Real.log C| + a / b) * y := by
    have h1 : -Real.log C ≤ |Real.log C| * y := by
      nlinarith [neg_le_abs (Real.log C), abs_nonneg (Real.log C)]
    have h2 : a * L ≤ a / b * y := by
      calc a * L ≤ a * (y / b) := mul_le_mul_of_nonneg_left hLy ha.le
        _ = a / b * y := by ring
    linarith
  -- (|log C| + a/b)·y ≤ C²y²
  have hquad : (|Real.log C| + a / b) * y ≤ C ^ 2 * y ^ 2 := by
    have hC2 : 0 < C ^ 2 := pow_pos hC 2
    have hM2 : (|Real.log C| + a / b) / C ^ 2 ≤ y := (le_max_right _ _).trans hyM
    rw [div_le_iff₀ hC2] at hM2
    nlinarith
  -- C²y² ≤ n εₙ²
  have hsq : C ^ 2 * y ^ 2 ≤ x * priorMassRadius C a k n ^ 2 := by
    have hxy : x * (x ^ (-a)) ^ 2 = y ^ 2 := by
      have e1 : x * (x ^ (-a)) ^ 2 = x ^ (1 + (-a + -a)) := by
        rw [Real.rpow_add hx0, Real.rpow_add hx0, Real.rpow_one]
        ring
      have e2 : y ^ 2 = x ^ (b + b) := by rw [hy_def, pow_two, ← Real.rpow_add hx0]
      rw [e1, e2]
      congr 1
      rw [hb_def]
      ring
    have hLk2 : 1 ≤ (L ^ k) ^ 2 := by nlinarith
    have hnn : 0 ≤ C ^ 2 * (x * (x ^ (-a)) ^ 2) := mul_nonneg (sq_nonneg C) (mul_nonneg hx0.le (sq_nonneg _))
    calc C ^ 2 * y ^ 2 = C ^ 2 * (x * (x ^ (-a)) ^ 2) * 1 := by rw [hxy, mul_one]
      _ ≤ C ^ 2 * (x * (x ^ (-a)) ^ 2) * (L ^ k) ^ 2 := mul_le_mul_of_nonneg_left hLk2 hnn
      _ = x * priorMassRadius C a k n ^ 2 := by rw [hε]; ring
  calc Real.log (1 / priorMassRadius C a k n) ≤ -Real.log C + a * L := hlog
    _ ≤ (|Real.log C| + a / b) * y := hstep
    _ ≤ C ^ 2 * y ^ 2 := hquad
    _ ≤ x * priorMassRadius C a k n ^ 2 := hsq

end BGP
