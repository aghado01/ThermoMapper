import Enthymemata.TimeFrequency.HaarGeometry
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Topology.Order.Compact

/-! # Gain-normalized families at a fixed articulation context -/

namespace TimeFrequency.Haar
noncomputable section

abbrev Block := Site → ℝ

def normalized (x : Block) : Point := fun j ↦ measured x j / meanSquare x

theorem meanSquare_smul (a : ℝ) (x : Block) : meanSquare (a • x) = a^2 * meanSquare x := by
  simp [meanSquare, mul_pow, ← Finset.mul_sum, mul_div_assoc]

theorem filtered_smul (b : ℕ) (a : ℝ) (x : Block) : filtered b (a • x) = a • filtered b x := by
  funext i
  simp [filtered, ← mul_sub, ← Finset.mul_sum, mul_div_assoc]

theorem normalized_smul {a : ℝ} (ha : a ≠ 0) (x : Block) : normalized (a • x) = normalized x := by
  funext j
  fin_cases j <;> simp [normalized, measured, filtered_smul, meanSquare_smul,
    mul_div_mul_left _ _ (pow_ne_zero 2 ha)]

theorem normalized_signal (p : Pitch) {t : ℝ} (ht : t ∈ Set.Icc 0 1) (φ ψ : ℝ) :
    normalized (signal p t φ ψ) = features p t := by
  funext j
  simp [normalized, signal_energy p ht, measured_signal p ht]

def rawFamily (I : Set ℝ) : Set Block :=
  {x | ∃ (p : Pitch) (t : ℝ), t ∈ I ∧ ∃ (a φ ψ : ℝ), 0 < a ∧ x = a • signal p t φ ψ}

def featureFamily (I : Set ℝ) : Set Point := features .three '' I ∪ features .four '' I

theorem normalized_family (I : Set ℝ) (hI : I ⊆ Set.Icc 0 1) :
    normalized '' rawFamily I = featureFamily I := by
  ext z
  constructor
  · rintro ⟨_, ⟨p,t,ht,a,φ,ψ,ha,rfl⟩,rfl⟩
    rw [normalized_smul ha.ne', normalized_signal _ (hI ht)]
    cases p
    · exact Or.inl ⟨_,ht,rfl⟩
    · exact Or.inr ⟨_,ht,rfl⟩
  · rintro (⟨t,ht,rfl⟩ | ⟨t,ht,rfl⟩) <;>
      exact ⟨signal _ t 0 0, ⟨_,t,ht,1,0,0,one_pos,(one_smul _ _).symm⟩,
        normalized_signal _ (hI ht) 0 0⟩

theorem overlap_waveform {I J : Set ℝ} {t : ℝ} (hi : t ∈ I) (hj : t ∈ J)
    (p : Pitch) (φ ψ : ℝ) {a : ℝ} (ha : 0 < a) :
    a • signal p t φ ψ ∈ rawFamily I ∩ rawFamily J :=
  ⟨⟨p,t,hi,a,φ,ψ,ha,rfl⟩, ⟨p,t,hj,a,φ,ψ,ha,rfl⟩⟩

private theorem distance3_euclidean (x y : Point) :
    distance3 x y = dist (WithLp.toLp 2 x : EuclideanSpace ℝ (Fin 3)) (WithLp.toLp 2 y) := by
  rw [EuclideanSpace.dist_eq]
  simp [distance3, length3, dot, difference, Fin.sum_univ_succ, Real.dist_eq, pow_two, add_assoc]

theorem distance3_comm (x y : Point) : distance3 x y = distance3 y x := by
  rw [distance3_euclidean, dist_comm, ← distance3_euclidean]

theorem distance3_triangle (x y z : Point) : distance3 x z ≤ distance3 x y + distance3 y z := by
  simpa only [distance3_euclidean] using dist_triangle
    (WithLp.toLp 2 x : EuclideanSpace ℝ (Fin 3)) (WithLp.toLp 2 y) (WithLp.toLp 2 z)

theorem same_pitch_distance (p : Pitch) (t u : ℝ) :
    distance3 (features p t) (features p u) = length3 (direction p) * |t-u| := by
  have h : dot (difference (features p t) (features p u))
      (difference (features p t) (features p u)) = (t-u)^2 * dot (direction p) (direction p) := by
    simp only [dot, difference, features, direction]
    ring
  unfold distance3 length3
  rw [h, Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq_eq_abs, mul_comm]

def minSlope := min (length3 (direction .three)) (length3 (direction .four))

theorem minSlope_lower : 1 / 2 ≤ minSlope := by
  obtain ⟨h,k,_⟩ := directions
  unfold minSlope length3
  rw [h,k, le_min_iff]
  constructor <;> apply Real.le_sqrt_of_sq_le <;> norm_num [dot, Matrix.cons_val_two]
  nlinarith [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 2), Real.sqrt_nonneg 2]

/-- Distance between closed intervals with endpoints in increasing order. -/
def intervalGap (a b c d : ℝ) := max 0 (max (c-b) (a-d))

theorem intervalGap_le {a b c d t u : ℝ} (ht : t ∈ Set.Icc a b) (hu : u ∈ Set.Icc c d) :
    intervalGap a b c d ≤ |t-u| := by
  exact max_le (abs_nonneg _) (max_le (by linarith [ht.2, hu.1, neg_le_abs (t-u)])
    (by linarith [ht.1, hu.2, le_abs_self (t-u)]))

theorem intervalGap_attained {a b c d : ℝ} (hab : a ≤ b) (hcd : c ≤ d) :
    ∃ t ∈ Set.Icc a b, ∃ u ∈ Set.Icc c d, |t-u| = intervalGap a b c d := by
  by_cases h : b < c
  · exact ⟨b,⟨hab,le_rfl⟩,c,⟨le_rfl,hcd⟩, by
      rw [intervalGap, max_eq_left (by linarith : a-d ≤ c-b), max_eq_right (by linarith),
        abs_of_neg (by linarith)]
      ring⟩
  by_cases k : d < a
  · exact ⟨a,⟨le_rfl,hab⟩,d,⟨hcd,le_rfl⟩, by
      rw [intervalGap,max_eq_right (by linarith : c-b ≤ a-d),max_eq_right (by linarith),
        abs_of_pos (by linarith)]⟩
  refine ⟨max a c, ⟨le_max_left .., max_le hab (le_of_not_gt h)⟩, max a c,
    ⟨le_max_right .., max_le (le_of_not_gt k) hcd⟩, ?_⟩
  rw [sub_self, abs_zero, intervalGap, max_eq_left (max_le (by linarith) (by linarith))]

theorem disjoint_gap_pos {a b c d : ℝ} (hab : a ≤ b) (hcd : c ≤ d)
    (h : Disjoint (Set.Icc a b) (Set.Icc c d)) : 0 < intervalGap a b c d := by
  obtain ⟨t,ht,u,hu,he⟩ := intervalGap_attained hab hcd
  rw [← he, abs_pos, sub_ne_zero]
  rintro rfl
  exact Set.disjoint_left.mp h ht hu

def familyMargin (a b c d : ℝ) := min (minSlope * intervalGap a b c d) (1/16)

theorem familyMargin_pos {a b c d : ℝ} (hab : a ≤ b) (hcd : c ≤ d)
    (h : Disjoint (Set.Icc a b) (Set.Icc c d)) : 0 < familyMargin a b c d :=
  lt_min (mul_pos (by linarith [minSlope_lower]) (disjoint_gap_pos hab hcd h)) (by norm_num)

theorem family_separation {a b c d : ℝ} {x y : Point} (hx : x ∈ featureFamily (Set.Icc a b))
    (hy : y ∈ featureFamily (Set.Icc c d)) : familyMargin a b c d ≤ distance3 x y := by
  have same (p : Pitch) {t u : ℝ} (ht : t ∈ Set.Icc a b) (hu : u ∈ Set.Icc c d) :
      familyMargin a b c d ≤ distance3 (features p t) (features p u) := by
    have hs : minSlope ≤ length3 (direction p) := by
      cases p
      · exact min_le_left ..
      · exact min_le_right ..
    rw [same_pitch_distance]
    exact (min_le_left ..).trans (mul_le_mul hs (intervalGap_le ht hu)
      (le_max_left ..) (by linarith [minSlope_lower]))
  rcases hx with ⟨t,ht,rfl⟩ | ⟨t,ht,rfl⟩ <;> rcases hy with ⟨u,hu,rfl⟩ | ⟨u,hu,rfl⟩
  · exact same _ ht hu
  · rw [distance3_comm]
    exact (min_le_right ..).trans (cross_pitch_gap t u).le
  · exact (min_le_right ..).trans (cross_pitch_gap u t).le
  · exact same _ ht hu

/-- The true witness beats every candidate from the competing family. -/
theorem nearest_family {a b c d η : ℝ} {x z : Point}
    (hx : x ∈ featureFamily (Set.Icc a b)) (herr : distance3 z x ≤ η)
    (hη : 2 * η < familyMargin a b c d) {y : Point} (hy : y ∈ featureFamily (Set.Icc c d)) :
    distance3 z x < distance3 z y := by
  have h := distance3_triangle x z y
  rw [distance3_comm x z] at h
  linarith [family_separation hx hy]

theorem featureFamily_compact (a b : ℝ) : IsCompact (featureFamily (Set.Icc a b)) := by
  have h (p : Pitch) : Continuous (features p) := by
    apply continuous_pi
    intro i
    unfold features
    fun_prop
  exact (isCompact_Icc.image (h .three)).union (isCompact_Icc.image (h .four))

/-- A nearest fit over both families exists, and every minimizer belongs to the true family. -/
theorem nearest_assignment {a b c d η : ℝ} {x z : Point}
    (hx : x ∈ featureFamily (Set.Icc a b)) (herr : distance3 z x ≤ η)
    (hη : 2 * η < familyMargin a b c d) :
    let S := featureFamily (Set.Icc a b) ∪ featureFamily (Set.Icc c d)
    (∃ y ∈ S, IsMinOn (distance3 z) S y) ∧
      ∀ y ∈ S, IsMinOn (distance3 z) S y → y ∈ featureFamily (Set.Icc a b) := by
  dsimp only
  have hcont : Continuous (distance3 z) := by
    unfold distance3 length3 dot difference
    fun_prop
  refine ⟨((featureFamily_compact a b).union (featureFamily_compact c d)).exists_isMinOn
    ⟨x, Or.inl hx⟩ hcont.continuousOn, ?_⟩
  intro y hy hmin
  rcases hy with h | h
  · exact h
  · exact False.elim ((nearest_family hx herr hη h).not_ge (isMinOn_iff.mp hmin x (Or.inl hx)))

/-- A clean family proxy, the clean recording, the observation and the fitted point. -/
theorem off_range_fit {x w z y : Point} {m η r : ℝ} (hm : distance3 x w ≤ m)
    (hη : distance3 w z ≤ η) (hr : distance3 z y ≤ r) : distance3 x y ≤ m + η + r := by
  linarith [distance3_triangle x w z, distance3_triangle x z y]

theorem approximate_fit {x w z y : Point} {m η ε : ℝ} (hm : distance3 x w ≤ m)
    (hη : distance3 w z ≤ η) (hfit : distance3 z y ≤ distance3 z x + ε) :
    distance3 x y ≤ 2*(m+η)+ε := by
  rw [distance3_comm z x] at hfit
  linarith [distance3_triangle x w z, distance3_triangle x z y]

/-- A sub-gap fit has the correct pitch and a fraction error controlled by the pitch slope. -/
theorem fitted_parameters (p q : Pitch) (t u e : ℝ)
    (hfit : distance3 (features p t) (features q u) ≤ e) (he : e < 1 / 16) :
    p = q ∧ length3 (direction p) * |t-u| ≤ e := by
  have hp : p = q := by
    cases p <;> cases q
    · rfl
    · rw [distance3_comm] at hfit
      linarith [cross_pitch_gap t u]
    · linarith [cross_pitch_gap u t]
    · rfl
  subst q
  exact ⟨rfl, same_pitch_distance p t u ▸ hfit⟩

theorem collision_family_margin :
    familyMargin (238/1000) (240/1000) (630/1000) (632/1000) = 1/16 := by
  norm_num [familyMargin, intervalGap, min_eq_right (show (1:ℝ)/16 ≤ minSlope*(39/100) by
    linarith [minSlope_lower])]

theorem collision_family_two_scales (φ ψ α β : ℝ) :
    ∃ x ∈ rawFamily (Set.Icc (238/1000) (240/1000)),
    ∃ y ∈ rawFamily (Set.Icc (630/1000) (632/1000)),
      ∀ i : Fin 2, normalized x i.castSucc = normalized y i.castSucc := by
  obtain ⟨ht,hT,hu,hU⟩ := collision_ranges
  refine ⟨signal .four collisionU α β, ⟨.four,_,⟨hu.le,hU.le⟩,1,_,_,one_pos,(one_smul _ _).symm⟩,
    signal .three collisionT φ ψ, ⟨.three,_,⟨ht.le,hT.le⟩,1,_,_,one_pos,(one_smul _ _).symm⟩, ?_⟩
  intro i
  simp only [normalized, signal_energy _ (show collisionT ∈ Set.Icc 0 1 from
    ⟨by linarith, by linarith⟩), signal_energy _ (show collisionU ∈ Set.Icc 0 1 from
    ⟨by linarith, by linarith⟩), div_one]
  exact (two_scale_collision φ ψ α β i).symm

#print axioms meanSquare_smul
#print axioms filtered_smul
#print axioms normalized_smul
#print axioms normalized_signal
#print axioms normalized_family
#print axioms overlap_waveform
#print axioms distance3_comm
#print axioms distance3_triangle
#print axioms same_pitch_distance
#print axioms minSlope_lower
#print axioms intervalGap_le
#print axioms intervalGap_attained
#print axioms disjoint_gap_pos
#print axioms familyMargin_pos
#print axioms family_separation
#print axioms featureFamily_compact
#print axioms nearest_assignment
#print axioms fitted_parameters
#print axioms nearest_family
#print axioms off_range_fit
#print axioms approximate_fit
#print axioms collision_family_margin
#print axioms collision_family_two_scales

end
end TimeFrequency.Haar
