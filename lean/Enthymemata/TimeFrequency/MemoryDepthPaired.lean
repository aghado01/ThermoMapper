import Enthymemata.TimeFrequency.MemoryDepthObservations

/-! # Paired adjacent powers recover a sign-symmetric block experiment -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset

def advance (b e : Bool) : Bool := if e then b else !b

def relative {d : ℕ} (x : Word (d + 1)) : Word d :=
  fun i ↦ decide (x i.castSucc = x i.succ)

def liftPath : (d : ℕ) → Bool → Word d → Word (d + 1)
  | 0, b, _ => fun _ ↦ b
  | d + 1, b, e => Fin.cons b (liftPath d (advance b (e 0)) (Fin.tail e))

@[simp] theorem liftPath_zero (d : ℕ) (b : Bool) (e : Word d) : liftPath d b e 0 = b := by
  cases d <;> rfl

private theorem advance_relative (a b : Bool) : advance a (decide (a = b)) = b := by
  cases a <;> cases b <;> decide

private theorem relative_advance (b e : Bool) : decide (b = advance b e) = e := by
  cases b <;> cases e <;> decide

theorem relative_liftPath (d : ℕ) (b : Bool) (e : Word d) :
    relative (liftPath d b e) = e := by
  induction d generalizing b with
  | zero => exact Subsingleton.elim _ _
  | succ d ih =>
    funext i
    refine Fin.cases ?_ (fun j ↦ ?_) i
    · simp [relative, liftPath, relative_advance]
    · exact congrFun (ih (advance b (e 0)) (Fin.tail e)) j

theorem liftPath_relative {d : ℕ} (x : Word (d + 1)) :
    liftPath d (x 0) (relative x) = x := by
  induction d with
  | zero => funext i; have hi : i = 0 := Fin.ext (by omega); subst i; rfl
  | succ d ih =>
    have ht : Fin.tail (relative x) = relative (Fin.tail x) := rfl
    have ha : advance (x 0) (relative x 0) = x (Fin.succ 0) := advance_relative _ _
    funext i
    refine Fin.cases ?_ (fun j ↦ ?_) i
    · rfl
    · simpa only [liftPath, Fin.cons_succ, ha, ht, Fin.tail] using
        congrFun (ih (Fin.tail x)) j

theorem liftPath_flip (d : ℕ) (b : Bool) (e : Word d) :
    liftPath d (!b) e = flip (liftPath d b e) := by
  induction d generalizing b with
  | zero => rfl
  | succ d ih =>
    have ha : advance (!b) (e 0) = !(advance b (e 0)) := by
      cases b <;> cases e 0 <;> rfl
    funext i
    refine Fin.cases ?_ (fun j ↦ ?_) i
    · rfl
    · simpa only [liftPath, Fin.cons_succ, ha, flip] using congrFun (ih _ (Fin.tail e)) j

def pathEquiv (d : ℕ) : Bool × Word d ≃ Word (d + 1) where
  toFun z := liftPath d z.1 z.2
  invFun x := (x 0, relative x)
  left_inv z := by simp [relative_liftPath]
  right_inv x := liftPath_relative x

def adjacentPower {d : ℕ} (x : Word (d + 1)) (i : Fin d) : ℤ :=
  (sign (x i.castSucc) + sign (x i.succ)) ^ 2

theorem adjacent_power_relative {d : ℕ} (x : Word (d + 1)) (i : Fin d) :
    adjacentPower x i = if relative x i then 4 else 0 := by
  unfold adjacentPower relative sign
  cases x i.castSucc <;> cases x i.succ <;> decide

theorem paired_reconstruction {d : ℕ} (x : Word (d + 1)) :
    liftPath d (x 0) (fun i ↦ decide (adjacentPower x i = 4)) = x := by
  have he : (fun i ↦ decide (adjacentPower x i = 4)) = relative x := by
    funext i
    rw [adjacent_power_relative]
    cases relative x i <;> decide
  rw [he, liftPath_relative]

noncomputable def relativeLaw {d : ℕ} (P : Law (Word (d + 1))) : Law (Word d) where
  mass e := P.mass (liftPath d false e) + P.mass (liftPath d true e)
  positive e := add_pos (P.positive _) (P.positive _)
  total := by
    have he := (pathEquiv d).sum_comp P.mass
    rw [Fintype.sum_prod_type] at he
    simpa only [pathEquiv, Equiv.coe_fn_mk, Fintype.sum_bool, sum_add_distrib,
      P.total, add_comm] using he

theorem relative_bins {d : ℕ} (P : Law (Word (d + 1))) (e : Word d) :
    bin P.mass relative e = (relativeLaw P).mass e := by
  classical
  unfold bin
  rw [← (pathEquiv d).sum_comp (fun x ↦ if relative x = e then P.mass x else 0)]
  simp [pathEquiv, relative_liftPath, Fintype.sum_prod_type, relativeLaw, add_comm]

private theorem half_parameter : (1 / 2 : ℝ) ∈ Set.Ioo 0 1 := by
  constructor <;> norm_num

theorem fair_lift_factorization {d : ℕ} (P : Law (Word (d + 1)))
    (hs : ∀ x, P.mass (flip x) = P.mass x) (b : Bool) (e : Word d) :
    (product (bernoulli (1 / 2) half_parameter) (relativeLaw P)).mass
        (b, e) = P.mass (liftPath d b e) := by
  have he : P.mass (liftPath d true e) = P.mass (liftPath d false e) := by
    rw [show true = !false from rfl, liftPath_flip, hs]
  cases b <;> simp [product, bernoulli, bitMass, relativeLaw, he] <;> ring

theorem randomized_lift_law {d : ℕ} (P : Law (Word (d + 1)))
    (hs : ∀ x, P.mass (flip x) = P.mass x) (x : Word (d + 1)) :
    bin (product (bernoulli (1 / 2) half_parameter) (relativeLaw P)).mass
        (fun z ↦ liftPath d z.1 z.2) x = P.mass x := by
  classical
  unfold bin
  simp_rw [fair_lift_factorization P hs]
  change (∑ z : Bool × Word d, if pathEquiv d z = x then P.mass (pathEquiv d z) else 0) = _
  rw [(pathEquiv d).sum_comp (fun y ↦ if y = x then P.mass y else 0)]
  simp

theorem randomized_lift_records {d : ℕ} (P : Law (Word (d + 1)))
    (hs : ∀ x, P.mass (flip x) = P.mass x) (n : ℕ) (x : Fin n → Word (d + 1)) :
    bin (iid n (product (bernoulli (1 / 2) half_parameter)
        (relativeLaw P))).mass (fun records i ↦ liftPath d (records i).1 (records i).2) x =
      (iid n P).mass x := by
  rw [iid_observation_bins n _
    (fun (_ : Fin n) (z : Bool × Word d) ↦ liftPath d z.1 z.2)]
  simp only [randomized_lift_law P hs, iid]

def historyEquiv (d : ℕ) : Word d × Bool ≃ Word (d + 1) :=
  (Equiv.prodComm _ _).trans (Fin.snocEquiv (fun _ : Fin (d + 1) ↦ Bool))

noncomputable def historyResponseLaw {d : ℕ} (P : Law (Word (d + 1))) :
    Law (Word d × Bool) where
  mass z := P.mass (historyEquiv d z)
  positive _ := P.positive _
  total := ((historyEquiv d).sum_comp P.mass).trans P.total

theorem history_response_bins {d : ℕ} (P : Law (Word (d + 1))) (z : Word d × Bool) :
    bin P.mass (historyEquiv d).symm z = (historyResponseLaw P).mass z := by
  classical
  unfold bin
  rw [← (historyEquiv d).sum_comp
    (fun x ↦ if (historyEquiv d).symm x = z then P.mass x else 0)]
  simp [historyResponseLaw]

theorem history_response_records {d : ℕ} (P : Law (Word (d + 1)))
    (n : ℕ) (z : Fin n → Word d × Bool) :
    bin (iid n P).mass (fun records i ↦ (historyEquiv d).symm (records i)) z =
      (iid n (historyResponseLaw P)).mass z := by
  rw [iid_observation_bins]
  simp only [history_response_bins, iid]

theorem paired_record_law {d : ℕ} (P : Law (Word (d + 1))) (n : ℕ)
    (e : Fin n → Word d) :
    bin (iid n P).mass
      (fun records i j ↦ decide (adjacentPower (records i) j = 4)) e =
        (iid n (relativeLaw P)).mass e := by
  have he : (fun records : Fin n → Word (d + 1) ↦
      fun i j ↦ decide (adjacentPower (records i) j = 4)) =
      (fun records i ↦ relative (records i)) := by
    funext records i j
    rw [adjacent_power_relative]
    cases relative (records i) j <;> decide
  rw [he, iid_observation_bins]
  simp only [relative_bins, iid]

#print axioms paired_record_law

#print axioms randomized_lift_records
#print axioms history_response_bins
#print axioms history_response_records

#print axioms liftPath_zero
#print axioms relative_liftPath
#print axioms liftPath_relative
#print axioms liftPath_flip
#print axioms adjacent_power_relative
#print axioms paired_reconstruction
#print axioms relative_bins
#print axioms fair_lift_factorization
#print axioms randomized_lift_law

end TimeFrequency.ContextTree.DepthEvidence
