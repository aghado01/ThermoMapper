import Enthymemata.TimeFrequency.HaarFamilies
import Mathlib.Analysis.SpecialFunctions.Log.Basic

/-! # Harmonic amplitude channels and the ambiguity in source balance -/

namespace TimeFrequency.Haar
noncomputable section

def channelPower (c₁ c₂ t : ℝ) := c₁^2*(1-t) + c₂^2*t
def observedFraction (c₁ c₂ t : ℝ) := c₂^2*t / channelPower c₁ c₂ t

/-- The same pitch and phases, with each harmonic amplitude scaled separately. -/
def channelSignal (p : Pitch) (t φ ψ c₁ c₂ : ℝ) (i : Site) :=
  tone (phasor (c₁ * Real.sqrt (2*(1-t))) φ) (bin p) i +
  tone (phasor (c₂ * Real.sqrt (2*t)) ψ) (2*bin p) i

theorem channelPower_pos {c₁ c₂ t : ℝ} (h₁ : 0 < c₁) (h₂ : 0 < c₂) (ht : t ∈ Set.Ioo 0 1) :
    0 < channelPower c₁ c₂ t :=
  add_pos (mul_pos (pow_pos h₁ _) (sub_pos.mpr ht.2)) (mul_pos (pow_pos h₂ _) ht.1)

theorem observedFraction_mem {c₁ c₂ t : ℝ} (h₁ : 0 < c₁) (h₂ : 0 < c₂)
    (ht : t ∈ Set.Ioo 0 1) : observedFraction c₁ c₂ t ∈ Set.Ioo 0 1 := by
  have hM := channelPower_pos h₁ h₂ ht
  exact ⟨div_pos (mul_pos (pow_pos h₂ _) ht.1) hM, (div_lt_one hM).mpr
    (lt_add_of_pos_left _ (mul_pos (pow_pos h₁ _) (sub_pos.mpr ht.2)))⟩

private theorem tone_scale (a r φ : ℝ) (j i : Site) :
    tone (phasor (a*r) φ) j i = a * tone (phasor r φ) j i := by
  simp [tone, phasor, mul_assoc]

theorem channel_as_signal (p : Pitch) (φ ψ : ℝ) {t c₁ c₂ : ℝ} (ht : t ∈ Set.Ioo 0 1)
    (h₁ : 0 < c₁) (h₂ : 0 < c₂) :
    channelSignal p t φ ψ c₁ c₂ =
      Real.sqrt (channelPower c₁ c₂ t) • signal p (observedFraction c₁ c₂ t) φ ψ := by
  let M := channelPower c₁ c₂ t
  let u := observedFraction c₁ c₂ t
  have hM : 0 < M := channelPower_pos h₁ h₂ ht
  have hMu : M*u = c₂^2*t := by
    change M*(c₂^2*t/M) = _
    field_simp
  have hMc : M*(1-u) = c₁^2*(1-t) := by
    rw [mul_sub, mul_one, hMu]
    dsimp [M,channelPower]
    ring
  have amp (c v w : ℝ) (hc : 0 < c) (he : M*w = c^2*v) :
      c * Real.sqrt (2*v) = Real.sqrt M * Real.sqrt (2*w) := by
    rw [← Real.sqrt_sq hc.le, ← Real.sqrt_mul (sq_nonneg c), ← Real.sqrt_mul hM.le]
    congr 1
    linarith only [he]
  have ha := amp c₁ (1-t) (1-u) h₁ hMc
  have hb := amp c₂ t u h₂ hMu
  funext i
  change _ + _ = Real.sqrt M * (tone _ _ _ + tone _ _ _)
  rw [ha, hb, tone_scale, tone_scale, mul_add]

theorem normalized_channel (p : Pitch) (φ ψ : ℝ) {t c₁ c₂ : ℝ} (ht : t ∈ Set.Ioo 0 1)
    (h₁ : 0 < c₁) (h₂ : 0 < c₂) :
    normalized (channelSignal p t φ ψ c₁ c₂) = features p (observedFraction c₁ c₂ t) := by
  rw [channel_as_signal p φ ψ ht h₁ h₂,
    normalized_smul (Real.sqrt_pos.2 (channelPower_pos h₁ h₂ ht)).ne']
  exact normalized_signal _ (Set.Ioo_subset_Icc_self (observedFraction_mem h₁ h₂ ht)) _ _

/-- Two distinct interior source fractions produce the very same recorded samples. -/
theorem channel_collision (p : Pitch) (φ ψ : ℝ) :
    channelSignal p (1/5) φ ψ 1 2 = channelSignal p (4/5) φ ψ 2 1 ∧
    normalized (channelSignal p (1/5) φ ψ 1 2) = features p (1/2) ∧
    normalized (channelSignal p (4/5) φ ψ 2 1) = features p (1/2) := by
  have a : (1/5 : ℝ) ∈ Set.Ioo 0 1 := ⟨by norm_num,by norm_num⟩
  have b : (4/5 : ℝ) ∈ Set.Ioo 0 1 := ⟨by norm_num,by norm_num⟩
  constructor
  · rw [channel_as_signal p φ ψ a (by norm_num) (by norm_num),
      channel_as_signal p φ ψ b (by norm_num) (by norm_num)]
    norm_num [observedFraction, channelPower]
  · rw [normalized_channel p φ ψ a (by norm_num) (by norm_num),
      normalized_channel p φ ψ b (by norm_num) (by norm_num)]
    norm_num [observedFraction, channelPower]

def logOdds (t : ℝ) := Real.log (t/(1-t))

theorem logOdds_observed {t c₁ c₂ : ℝ} (ht : t ∈ Set.Ioo 0 1) (h₁ : 0 < c₁) (h₂ : 0 < c₂) :
    logOdds (observedFraction c₁ c₂ t) = logOdds t + 2 * Real.log (c₂/c₁) := by
  have hc := sub_pos.mpr ht.2
  have hM := channelPower_pos h₁ h₂ ht
  have elim : 1 - observedFraction c₁ c₂ t = c₁^2*(1-t)/channelPower c₁ c₂ t := by
    unfold observedFraction
    field_simp
    unfold channelPower
    ring
  unfold logOdds
  rw [elim]
  have odds : observedFraction c₁ c₂ t / (c₁^2*(1-t)/channelPower c₁ c₂ t) =
      (t/(1-t)) * (c₂/c₁)^2 := by
    unfold observedFraction
    field_simp
  rw [odds, Real.log_mul (div_pos ht.1 hc).ne' (pow_ne_zero _ (div_pos h₂ h₁).ne'), Real.log_pow]
  norm_num

theorem shared_channel_difference {t u c₁ c₂ : ℝ} (ht : t ∈ Set.Ioo 0 1)
    (hu : u ∈ Set.Ioo 0 1) (h₁ : 0 < c₁) (h₂ : 0 < c₂) :
    logOdds (observedFraction c₁ c₂ u) - logOdds (observedFraction c₁ c₂ t) =
      logOdds u - logOdds t := by
  rw [logOdds_observed ht h₁ h₂, logOdds_observed hu h₁ h₂]
  ring

theorem independent_channel_gap {t u c₁ c₂ d₁ d₂ B : ℝ} (ht : t ∈ Set.Ioo 0 1)
    (hu : u ∈ Set.Ioo 0 1) (hc₁ : 0 < c₁) (hc₂ : 0 < c₂) (hd₁ : 0 < d₁) (hd₂ : 0 < d₂)
    (hc : |Real.log (c₂ / c₁)| ≤ B) (hd : |Real.log (d₂ / d₁)| ≤ B) :
    |logOdds u - logOdds t| - 4*B ≤
      |logOdds (observedFraction d₁ d₂ u) - logOdds (observedFraction c₁ c₂ t)| := by
  rw [logOdds_observed ht hc₁ hc₂, logOdds_observed hu hd₁ hd₂]
  obtain ⟨hc',hc''⟩ := abs_le.mp hc
  obtain ⟨hd',hd''⟩ := abs_le.mp hd
  cases le_total (logOdds t) (logOdds u) with
  | inl h =>
    rw [abs_of_nonneg (sub_nonneg.mpr h)]
    linarith [le_abs_self ((logOdds u + 2*Real.log (d₂/d₁)) - (logOdds t + 2*Real.log (c₂/c₁)))]
  | inr h =>
    rw [abs_of_nonpos (sub_nonpos.mpr h)]
    linarith [neg_le_abs ((logOdds u + 2*Real.log (d₂/d₁)) - (logOdds t + 2*Real.log (c₂/c₁)))]

#print axioms channelPower_pos
#print axioms observedFraction_mem
#print axioms channel_as_signal
#print axioms normalized_channel
#print axioms channel_collision
#print axioms logOdds_observed
#print axioms shared_channel_difference
#print axioms independent_channel_gap

end
end TimeFrequency.Haar
