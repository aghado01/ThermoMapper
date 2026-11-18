import Enthymemata.TimeFrequency.HaarGains

/-! # Two-scale collision and three-scale separation on the harmonic family -/

namespace TimeFrequency.Haar

noncomputable section
abbrev Point := Fin 3 → ℝ

def difference (x y : Point) : Point := fun i ↦ x i - y i
def dot (x y : Point) := x 0*y 0+x 1*y 1+x 2*y 2
def cross (x y : Point) : Point :=
  ![x 1*y 2-x 2*y 1,x 2*y 0-x 0*y 2,x 0*y 1-x 1*y 0]
def length3 (x : Point) := Real.sqrt (dot x x)
def distance3 (x y : Point) := length3 (difference x y)

def direction (p : Pitch) := difference (gains (2*bin p)) (gains (bin p))
def offset := difference (gains 4) (gains 3)
def normal := cross (direction .three) (direction .four)

private theorem dot_self (x : Point) : 0 ≤ dot x x := by
  unfold dot
  exact add_nonneg (add_nonneg (mul_self_nonneg _) (mul_self_nonneg _)) (mul_self_nonneg _)

private theorem cauchy (x y : Point) : (dot x y)^2 ≤ dot x x * dot y y := by
  unfold dot
  nlinarith [sq_nonneg (x 0*y 1-x 1*y 0), sq_nonneg (x 0*y 2-x 2*y 0),
    sq_nonneg (x 1*y 2-x 2*y 1)]

theorem directions :
    direction .three = ![Real.sqrt 2/2,1-Real.sqrt 2/2,-(2+Real.sqrt 2)] ∧
    direction .four = ![1,-3/2,-3/4] ∧
    offset = ![Real.sqrt 2/2-1/2,5/4-Real.sqrt 2/2,-7/8-Real.sqrt 2] := by
  obtain ⟨h3,h6,h4,h8⟩ := bin_gains
  have a : 2*bin .three = 6 := by decide
  have b : 2*bin .four = 8 := by decide
  simp only [direction, a, b, offset]
  simp only [bin, h3, h4, h6, h8]
  refine ⟨?_,?_,?_⟩ <;> funext i <;> fin_cases i <;>
    norm_num [difference, Matrix.cons_val_two] <;> ring

theorem directions_nonzero : direction .three ≠ 0 ∧ direction .four ≠ 0 := by
  obtain ⟨ha,hc,_⟩ := directions
  constructor <;> intro h
  · have hh := congrFun (ha.symm.trans h) 0
    norm_num at hh
  · have hh := congrFun (hc.symm.trans h) 0
    norm_num at hh

theorem normal_certificate :
    dot offset normal = (2+Real.sqrt 2)/8 ∧ dot normal normal < 256*(dot offset normal)^2 := by
  obtain ⟨ha,hc,hw⟩ := directions
  have hr := Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 2)
  have hp := Real.sqrt_nonneg 2
  have hn : normal = ![-(15/4+9*Real.sqrt 2/8),-(2+5*Real.sqrt 2/8),
      -(1+Real.sqrt 2/4)] := by
    rw [normal, ha, hc]
    funext i
    fin_cases i <;> norm_num [cross, Matrix.cons_val_two] <;> ring
  have hd : dot offset normal = (2+Real.sqrt 2)/8 := by
    rw [hw,hn]
    norm_num [dot, Matrix.cons_val_two]
    nlinarith
  refine ⟨hd,?_⟩
  rw [hd,hn]
  norm_num [dot, Matrix.cons_val_two]
  nlinarith

theorem normal_offset (t u : ℝ) :
    dot (difference (features .four u) (features .three t)) normal = dot offset normal := by
  have hw : offset = difference (gains (bin .four)) (gains (bin .three)) := rfl
  rw [hw]
  simp [dot, normal, cross, features, direction, difference, Matrix.cons_val_two]
  ring

theorem cross_pitch_gap (t u : ℝ) :
    1/16 < distance3 (features .four u) (features .three t) := by
  let x := difference (features .four u) (features .three t)
  have hc := cauchy x normal
  -- The projection onto the normal is constant along both supporting lines.
  have hx : dot x normal = dot offset normal := normal_offset t u
  rw [hx] at hc
  have hh := normal_certificate.2
  have hn := dot_self normal
  have hd : (1/16:ℝ)^2 < dot x x := by
    by_contra! h
    have hb := mul_le_mul_of_nonneg_right h hn
    nlinarith
  change 1/16 < Real.sqrt (dot x x)
  nlinarith [Real.sq_sqrt (dot_self x), Real.sqrt_nonneg (dot x x)]

def collisionT : ℝ := (3+Real.sqrt 2)/7
def collisionU : ℝ := (9-4*Real.sqrt 2)/14

theorem collision_ranges :
    630/1000 < collisionT ∧ collisionT < 632/1000 ∧
    238/1000 < collisionU ∧ collisionU < 240/1000 := by
  have h := Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 2)
  have hp := Real.sqrt_nonneg 2
  have hlo : (1414/1000:ℝ) < Real.sqrt 2 := by nlinarith
  have hhi : Real.sqrt 2 < (1415/1000:ℝ) := by nlinarith
  dsimp [collisionT,collisionU]
  exact ⟨by linarith,by linarith,by linarith,by linarith⟩

theorem two_scale_collision (φ ψ α β : ℝ) :
    ∀ i : Fin 2, measured (signal .three collisionT φ ψ) i.castSucc =
      measured (signal .four collisionU α β) i.castSucc := by
  obtain ⟨h1,h2,h3,h4⟩ := collision_ranges
  have ht : collisionT ∈ Set.Icc 0 1 := ⟨by linarith, by linarith⟩
  have hu : collisionU ∈ Set.Icc 0 1 := ⟨by linarith, by linarith⟩
  rw [measured_signal _ ht, measured_signal _ hu]
  obtain ⟨g3,g6,g4,g8⟩ := bin_gains
  have a : 2*bin .three = 6 := by decide
  have b : 2*bin .four = 8 := by decide
  intro i
  simp only [features, a, b]
  fin_cases i <;> simp [bin, g3, g4, g6, g8, collisionT, collisionU] <;>
    nlinarith [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 2)]

private theorem balance_identified (p : Pitch) {t u : ℝ}
    (h : features p t = features p u) : t = u := by
  by_contra hne
  have hd : direction p = 0 := by
    funext i
    have hi := congrFun h i
    change (1-t)*_+t*_ = (1-u)*_+u*_ at hi
    change gains _ i - gains _ i = 0
    apply (mul_eq_zero.mp (show (t-u)*(gains (2*bin p) i-gains (bin p) i) = 0 by
      linarith only [hi])).resolve_left (sub_ne_zero.mpr hne)
  cases p
  · exact directions_nonzero.1 hd
  · exact directions_nonzero.2 hd

theorem three_scale_identification (p q : Pitch) {t u : ℝ}
    (ht : t ∈ Set.Icc 0 1) (hu : u ∈ Set.Icc 0 1) (φ ψ α β : ℝ)
    (h : measured (signal p t φ ψ) = measured (signal q u α β)) : p = q ∧ t = u := by
  rw [measured_signal _ ht, measured_signal _ hu] at h
  have sep (x y : ℝ) : features .four y ≠ features .three x := by
    intro he
    have hb := cross_pitch_gap x y
    rw [he] at hb
    norm_num [distance3, length3, difference, dot] at hb
  cases p <;> cases q
  · exact ⟨rfl, balance_identified _ h⟩
  · exact False.elim (sep _ _ h.symm)
  · exact False.elim (sep _ _ h)
  · exact ⟨rfl, balance_identified _ h⟩

#print axioms normal_offset
#print axioms directions
#print axioms directions_nonzero
#print axioms normal_certificate
#print axioms cross_pitch_gap
#print axioms collision_ranges
#print axioms two_scale_collision
#print axioms three_scale_identification

end
end TimeFrequency.Haar
