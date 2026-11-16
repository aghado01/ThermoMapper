import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
import Mathlib.Analysis.SpecialFunctions.Trigonometric.InverseDeriv
import Mathlib.Dynamics.FixedPoints.Topology
import Mathlib.Topology.MetricSpace.Lipschitz
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Topology.Order.MonotoneConvergence
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# The deterministic ReLU angle map (S17)

Hand and Voroninski, equation `defn-g`. For a starting angle `θ`, the comparison angle after
`n` layers is `reluAngle^[n] θ`. Angles are in radians.
-/

open Real Set Function Filter
open scoped Topology

namespace GenerativePrior

noncomputable def reluAngleCos (θ : ℝ) : ℝ := ((π - θ) * cos θ + sin θ) / π

noncomputable def reluAngle (θ : ℝ) : ℝ := arccos (reluAngleCos θ)

theorem hasDerivAt_reluAngleCos (θ : ℝ) :
    HasDerivAt reluAngleCos (-((π - θ) / π) * sin θ) θ := by
  have h := ((((hasDerivAt_id θ).const_sub π).mul (hasDerivAt_cos θ)).add
    (hasDerivAt_sin θ)).div_const π
  exact h.congr_deriv (by simp only [id_eq]; ring)

theorem continuous_reluAngle : Continuous reluAngle :=
  continuous_arccos.comp (continuous_iff_continuousAt.mpr fun x ↦
    (hasDerivAt_reluAngleCos x).continuousAt)

@[simp] theorem reluAngleCos_zero : reluAngleCos 0 = 1 := by simp [reluAngleCos, pi_ne_zero]

@[simp] theorem reluAngleCos_pi : reluAngleCos π = 0 := by simp [reluAngleCos]

theorem strictAntiOn_reluAngleCos : StrictAntiOn reluAngleCos (Icc 0 π) := by
  apply strictAntiOn_of_deriv_neg (convex_Icc ..)
    (fun x _ ↦ (hasDerivAt_reluAngleCos x).continuousAt.continuousWithinAt)
  intro x hx
  rw [interior_Icc] at hx
  rw [(hasDerivAt_reluAngleCos x).deriv]
  exact mul_neg_of_neg_of_pos (neg_lt_zero.mpr (div_pos (sub_pos.mpr hx.2) pi_pos))
    (sin_pos_of_pos_of_lt_pi hx.1 hx.2)

theorem reluAngleCos_mem {θ : ℝ} (hθ : θ ∈ Icc 0 π) : reluAngleCos θ ∈ Icc 0 1 := by
  have anti := strictAntiOn_reluAngleCos.antitoneOn
  exact ⟨by simpa using anti hθ (right_mem_Icc.mpr pi_pos.le) hθ.2,
    by simpa using anti (left_mem_Icc.mpr pi_pos.le) hθ hθ.1⟩

theorem cos_reluAngle {θ : ℝ} (hθ : θ ∈ Icc 0 π) : cos (reluAngle θ) = reluAngleCos θ :=
  cos_arccos (by linarith [(reluAngleCos_mem hθ).1]) (reluAngleCos_mem hθ).2

@[simp] theorem reluAngle_zero : reluAngle 0 = 0 := by simp [reluAngle]

@[simp] theorem reluAngle_pi : reluAngle π = π / 2 := by simp [reluAngle]

theorem reluAngle_mem {θ : ℝ} (hθ : θ ∈ Icc 0 π) : reluAngle θ ∈ Icc 0 (π / 2) :=
  ⟨arccos_nonneg _, arccos_le_pi_div_two.mpr (reluAngleCos_mem hθ).1⟩

theorem reluAngle_mapsTo : MapsTo reluAngle (Icc 0 π) (Icc 0 π) :=
  fun _ h ↦ ⟨(reluAngle_mem h).1, (reluAngle_mem h).2.trans (by linarith [pi_pos])⟩

theorem strictMonoOn_reluAngle : StrictMonoOn reluAngle (Icc 0 π) := by
  intro x hx y hy hxy
  exact arccos_lt_arccos (by linarith [(reluAngleCos_mem hy).1])
    (strictAntiOn_reluAngleCos hx hy hxy) (reluAngleCos_mem hx).2

theorem reluAngle_image : reluAngle '' Icc 0 π = Icc 0 (π / 2) := by
  simpa using continuous_reluAngle.continuousOn.image_Icc_of_monotoneOn pi_pos.le
    strictMonoOn_reluAngle.monotoneOn

theorem mul_cos_lt_sin {θ : ℝ} (hθ : θ ∈ Ioc 0 π) : θ * cos θ < sin θ := by
  have inc : StrictMonoOn (fun x : ℝ ↦ sin x - x * cos x) (Icc 0 π) := by
    apply strictMonoOn_of_deriv_pos (convex_Icc ..)
      (continuous_sin.sub (continuous_id.mul continuous_cos)).continuousOn
    intro x hx
    rw [interior_Icc] at hx
    have hd := (hasDerivAt_sin x).sub ((hasDerivAt_id x).mul (hasDerivAt_cos x))
    rw [hd.deriv]
    dsimp only [id_eq]
    linarith only [mul_pos hx.1 (sin_pos_of_pos_of_lt_pi hx.1 hx.2)]
  have h := inc (left_mem_Icc.mpr pi_pos.le) ⟨hθ.1.le, hθ.2⟩ hθ.1
  simpa using h

theorem reluAngle_lt {θ : ℝ} (hθ : θ ∈ Ioc 0 π) : reluAngle θ < θ := by
  have hcos : cos θ < reluAngleCos θ := by
    apply (lt_div_iff₀ pi_pos).mpr
    linarith [mul_cos_lt_sin hθ]
  simpa only [arccos_cos hθ.1.le hθ.2, reluAngle] using
    arccos_lt_arccos (neg_one_le_cos θ) hcos (reluAngleCos_mem ⟨hθ.1.le, hθ.2⟩).2

/-- The interior derivative. The closed-interval bound below also covers the endpoint at zero. -/
theorem hasDerivAt_reluAngle {θ : ℝ} (hθ : θ ∈ Ioo 0 π) :
    HasDerivAt reluAngle (((π - θ) / π * sin θ) / √(1 - reluAngleCos θ ^ 2)) θ := by
  have h₁ : reluAngleCos θ < 1 := by
    simpa using strictAntiOn_reluAngleCos (left_mem_Icc.mpr pi_pos.le) (Ioo_subset_Icc_self hθ) hθ.1
  have hn := (reluAngleCos_mem (Ioo_subset_Icc_self hθ)).1
  have hd := (hasDerivAt_arccos (by linarith : reluAngleCos θ ≠ -1) h₁.ne).comp θ
    (hasDerivAt_reluAngleCos θ)
  exact hd.congr_deriv (by ring)

theorem deriv_reluAngle_mem {θ : ℝ} (hθ : θ ∈ Ioo 0 π) : deriv reluAngle θ ∈ Icc 0 1 := by
  let a := (π - θ) / π
  let b := sin θ / π
  have hs := sin_pos_of_pos_of_lt_pi hθ.1 hθ.2
  have ha : 0 < a := div_pos (sub_pos.mpr hθ.2) pi_pos
  have hb : 0 ≤ b := (div_pos hs pi_pos).le
  have hab : a + b ≤ 1 := by
    dsimp [a, b]
    rw [← add_div, div_le_one pi_pos]
    linarith [sin_le hθ.1.le]
  have hcos : reluAngleCos θ = a * cos θ + b := by dsimp [reluAngleCos, a, b]; ring
  have htrig := congrArg (fun t : ℝ ↦ a ^ 2 * t) (sin_sq_add_cos_sq θ)
  -- (a cos θ + b)² + (a sin θ)² ≤ (a + b)² ≤ 1.
  have hsq : (a * sin θ) ^ 2 ≤ 1 - reluAngleCos θ ^ 2 := by
    rw [hcos]
    nlinarith only [htrig, mul_nonneg (mul_nonneg ha.le hb) (sub_nonneg.mpr (cos_le_one θ)),
      pow_le_pow_left₀ (add_nonneg ha.le hb) hab 2]
  have hden := sqrt_pos.mpr ((pow_pos (mul_pos ha hs) 2).trans_le hsq)
  rw [(hasDerivAt_reluAngle hθ).deriv]
  exact ⟨(div_pos (mul_pos ha hs) hden).le, (div_le_one hden).mpr (le_sqrt_of_sq_le hsq)⟩

/-- The deterministic step does not amplify angular error, including between endpoints. -/
theorem lipschitzOnWith_reluAngle : LipschitzOnWith 1 reluAngle (Icc 0 π) := by
  have ordered (x y : ℝ) (hx : x ∈ Icc 0 π) (hy : y ∈ Icc 0 π) (hxy : x ≤ y) :
      dist (reluAngle y) (reluAngle x) ≤ dist y x := by
    rw [Real.dist_eq, Real.dist_eq, abs_of_nonneg (sub_nonneg.mpr
      (strictMonoOn_reluAngle.monotoneOn hx hy hxy)), abs_of_nonneg (sub_nonneg.mpr hxy)]
    have hd : DifferentiableOn ℝ reluAngle (Ioo 0 π) :=
      fun _ hz ↦ (hasDerivAt_reluAngle hz).differentiableAt.differentiableWithinAt
    have hb := (convex_Icc (0 : ℝ) π).image_sub_le_mul_sub_of_deriv_le
      continuous_reluAngle.continuousOn (by rwa [interior_Icc])
      (fun z hz ↦ (deriv_reluAngle_mem (show z ∈ Ioo 0 π by simpa only [interior_Icc] using hz)).2)
      x hx y hy hxy
    simpa using hb
  apply LipschitzOnWith.of_dist_le_mul
  intro x hx y hy
  rcases le_total x y with h | h
  · rw [dist_comm (reluAngle x), dist_comm x]
    simpa using ordered x y hx hy h
  · simpa using ordered y x hy hx h

/-- After any positive number of layers the angle is in the first quadrant. -/
theorem reluAngle_iterate_succ_mem {θ : ℝ} (hθ : θ ∈ Icc 0 π) (n : ℕ) :
    reluAngle^[n + 1] θ ∈ Icc 0 (π / 2) := by
  rw [iterate_succ_apply']
  exact reluAngle_mem (reluAngle_mapsTo.iterate n hθ)

theorem reluAngle_iterate_antitone {θ : ℝ} (hθ : θ ∈ Icc 0 π) :
    Antitone fun n : ℕ ↦ reluAngle^[n] θ := by
  apply antitone_nat_of_succ_le
  intro n
  rw [iterate_succ_apply']
  have hn := reluAngle_mapsTo.iterate n hθ
  rcases hn.1.eq_or_lt with hz | hp
  · rw [← hz, reluAngle_zero]
  · exact (reluAngle_lt ⟨hp, hn.2⟩).le

theorem reluAngle_iterate_strictAnti {θ : ℝ} (hθ : θ ∈ Ioc 0 π) :
    StrictAnti fun n : ℕ ↦ reluAngle^[n] θ := by
  have maps : MapsTo reluAngle (Ioc 0 π) (Ioc 0 π) := by
    intro x hx
    have hxc : x ∈ Icc 0 π := ⟨hx.1.le, hx.2⟩
    exact ⟨by simpa using strictMonoOn_reluAngle (left_mem_Icc.mpr pi_pos.le) hxc hx.1,
      (reluAngle_mapsTo hxc).2⟩
  apply strictAnti_nat_of_succ_lt
  intro n
  rw [iterate_succ_apply']
  exact reluAngle_lt (maps.iterate n hθ)

theorem reluAngle_iterate_tendsto_zero {θ : ℝ} (hθ : θ ∈ Icc 0 π) :
    Tendsto (fun n : ℕ ↦ reluAngle^[n] θ) atTop (𝓝 0) := by
  have hb : BddBelow (range fun n : ℕ ↦ reluAngle^[n] θ) :=
    ⟨0, by rintro _ ⟨n, rfl⟩; exact (reluAngle_mapsTo.iterate n hθ).1⟩
  let L := ⨅ n : ℕ, reluAngle^[n] θ
  have hL : Tendsto (fun n : ℕ ↦ reluAngle^[n] θ) atTop (𝓝 L) :=
    tendsto_atTop_ciInf (reluAngle_iterate_antitone hθ) hb
  have hmem : L ∈ Icc 0 π := isClosed_Icc.mem_of_tendsto hL
    (Eventually.of_forall fun n ↦ reluAngle_mapsTo.iterate n hθ)
  have hfix : reluAngle L = L := isFixedPt_of_tendsto_iterate hL continuous_reluAngle.continuousAt
  have hzero : L = 0 := by
    rcases hmem.1.eq_or_lt with hz | hp
    · exact hz.symm
    · exact ((reluAngle_lt ⟨hp, hmem.2⟩).ne hfix).elim
  rwa [hzero] at hL

/-- The antipodal start attains the two-layer angular bound. -/
theorem reluAngle_iterate_two_pi : reluAngle^[2] π = arccos (1 / π) := by
  rw [show reluAngle^[2] π = reluAngle (reluAngle π) by rfl, reluAngle_pi]
  simp [reluAngle, reluAngleCos]

theorem reluAngle_iterate_le_arccos_inv_pi {θ : ℝ} (hθ : θ ∈ Icc 0 π) {n : ℕ} (hn : 2 ≤ n) :
    reluAngle^[n] θ ≤ arccos (1 / π) := by
  apply (reluAngle_iterate_antitone hθ hn).trans
  rw [← reluAngle_iterate_two_pi]
  exact strictMonoOn_reluAngle.monotoneOn (reluAngle_mapsTo hθ)
    (reluAngle_mapsTo (right_mem_Icc.mpr pi_pos.le))
    (strictMonoOn_reluAngle.monotoneOn hθ (right_mem_Icc.mpr pi_pos.le) hθ.2)

/-- Errors in a sequence of measured angles accumulate by at most one step tolerance per layer. -/
theorem reluAngle_iterate_error {u : ℕ → ℝ} {δ : ℝ} (hu : ∀ n, u n ∈ Icc 0 π)
    (hstep : ∀ n, |u (n + 1) - reluAngle (u n)| ≤ δ) (n : ℕ) :
    |u n - reluAngle^[n] (u 0)| ≤ n * δ := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [iterate_succ_apply']
    have hLip := lipschitzOnWith_reluAngle.dist_le_mul _ (hu n) _
      (reluAngle_mapsTo.iterate n (hu 0))
    simp only [NNReal.coe_one, one_mul, Real.dist_eq] at hLip
    calc
      _ ≤ |u (n + 1) - reluAngle (u n)| +
          |reluAngle (u n) - reluAngle (reluAngle^[n] (u 0))| := abs_sub_le ..
      _ ≤ δ + (n : ℝ) * δ := add_le_add (hstep n) (hLip.trans ih)
      _ = _ := by push_cast; ring

#print axioms hasDerivAt_reluAngleCos
#print axioms continuous_reluAngle
#print axioms reluAngleCos_zero
#print axioms reluAngleCos_pi
#print axioms strictAntiOn_reluAngleCos
#print axioms reluAngleCos_mem
#print axioms cos_reluAngle
#print axioms reluAngle_zero
#print axioms reluAngle_pi
#print axioms reluAngle_mem
#print axioms reluAngle_mapsTo
#print axioms strictMonoOn_reluAngle
#print axioms reluAngle_image
#print axioms mul_cos_lt_sin
#print axioms reluAngle_lt
#print axioms hasDerivAt_reluAngle
#print axioms deriv_reluAngle_mem
#print axioms lipschitzOnWith_reluAngle
#print axioms reluAngle_iterate_succ_mem
#print axioms reluAngle_iterate_antitone
#print axioms reluAngle_iterate_strictAnti
#print axioms reluAngle_iterate_tendsto_zero
#print axioms reluAngle_iterate_two_pi
#print axioms reluAngle_iterate_le_arccos_inv_pi
#print axioms reluAngle_iterate_error

end GenerativePrior
