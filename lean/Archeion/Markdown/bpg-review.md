The proofs are mathematically sound, well-documented, and cleanly structured. Before promoting this into the main library, there is **one substantive downstream API trap** (in **UV-2b**) and **two minor Lean/linter items** worth addressing.

---

### 1. Substantive Catch: The `Real.log` Junk-Value Trap in `UV-2b`

Look at the signature of `log_inv_radius_le_n_radius_sq`:

```lean
theorem log_inv_radius_le_n_radius_sq (C a k : ℝ) (hC : 0 < C) (ha : 0 < a) (ha' : a < 1 / 2) (hk : 0 ≤ k) :
    ∀ᶠ n : ℕ in atTop, Real.log (1 / priorMassRadius C a k n) ≤ n * priorMassRadius C a k n ^ 2

```

In Mathlib, `Real.log x = 0` for `x ≤ 0`. Because of this junk value, the inequality `Real.log (1 / εₙ) ≤ n * εₙ ^ 2` is **trivially satisfied whenever `εₙ ≤ 0**` (since the LHS is `0` and the RHS is non-negative) and therefore **carries zero positivity information about `εₙ` to the caller**.

When the main contraction theorem consumes `UV-2b` to lower-bound the noise slab's prior mass via


$$\Pi(\sigma \in I_n) \gtrsim \varepsilon_n \ge \exp(-n \varepsilon_n^2),$$


rewriting `Real.log (1 / εₙ) ≤ n * εₙ ^ 2` into `Real.exp (-(n * εₙ ^ 2)) ≤ εₙ` requires `0 < priorMassRadius C a k n`. Right now, you prove `0 < x ^ (-a)` and `0 < L ^ k` *inside* `log_inv_radius_le_n_radius_sq`, and that positivity proof gets discarded when the theorem closes—forcing the caller to re-prove `n ≥ ⌈Real.exp 1⌉₊ → 0 < priorMassRadius C a k n` from scratch.

**Two ways to make `UV-2b` plug-and-play for the caller:**

* **Option A (Expose positivity & the exp form directly):**
Add a companion lemma (or conjunction) so the caller never has to unfold `priorMassRadius` or `Real.exp 1`:
```lean
theorem eventually_priorMassRadius_pos (C a k : ℝ) (hC : 0 < C) :
    ∀ᶠ n : ℕ in atTop, 0 < priorMassRadius C a k n := by
  refine Filter.eventually_atTop.2 ⟨⌈Real.exp 1⌉₊, fun n hn => ?_⟩
  have hxe : Real.exp 1 ≤ (n : ℝ) := (Nat.le_ceil _).trans (by exact_mod_cast hn)
  have hx0 : 0 < (n : ℝ) := lt_of_lt_of_le (Real.exp_pos 1) hxe
  have hL0 : 0 < Real.log (n : ℝ) := by
    have : 1 ≤ Real.log (n : ℝ) := by
      rw [← Real.log_exp 1]; exact Real.log_le_log (Real.exp_pos 1) hxe
    linarith
  exact mul_pos (mul_pos hC (Real.rpow_pos_of_pos hx0 _)) (Real.rpow_pos_of_pos hL0 _)

theorem exp_neg_n_radius_sq_le_radius (C a k : ℝ) (hC : 0 < C) (ha : 0 < a) (ha' : a < 1 / 2) (hk : 0 ≤ k) :
    ∀ᶠ n : ℕ in atTop, Real.exp (-(n * priorMassRadius C a k n ^ 2)) ≤ priorMassRadius C a k n := by
  filter_upwards [eventually_priorMassRadius_pos C a k hC,
    log_inv_radius_le_n_radius_sq C a k hC ha ha' hk] with n hpos hlog
  rw [← inv_inv (priorMassRadius C a k n), Real.exp_neg, inv_le_inv₀ (Real.exp_pos _) hpos,
    ← Real.log_le_iff_le_exp (inv_pos.2 hpos), ← one_div]
  exact hlog

```


* **Option B (Constant multipliers for the slab density):**
In the contraction proof, the $\sigma$-interval has width $c_1 \varepsilon_n$ (from the prior density lower bound on $[\sigma_0, \sigma_0(1 + c\varepsilon_n)]$) while the exponent budget is $c_2 n \varepsilon_n^2$. Notice that applying your current lemma to $C' = \min(c_1, \sqrt{c_2}) \cdot C$ already handles this cleanly by monotonicity of $\log$, *provided* you have `eventually_priorMassRadius_pos` available.

---

### 2. Minor Lean / Linter Cleanups

1. **Unused variable warning on `hr0'` in `gaussian_kl_scale_part_le` (Line 34):**
```lean
have hr0 : 0 < r := by linarith
have hr0' : r ≠ 0 := hr0.ne'

```


In current Mathlib, `field_simp` uses `positivity` as its discharger, which sees `hr0 : 0 < r` directly and never touches `hr0'`. If your project runs with `linter.unusedVariables` enabled, `hr0'` will trigger a warning. You can safely delete `have hr0' : r ≠ 0 := hr0.ne'`.
2. **Optional tightening in `gaussian_v_scale_part_le` (Lines 49–50):**
```lean
have h1 : (r + 1) ^ 2 ≤ 9 := by
  nlinarith [mul_nonneg (sub_nonneg.2 hhi) (by linarith : (0 : ℝ) ≤ r + 4)]

```


This works fine as written (`linarith` uses `hlo` implicitly to get `0 ≤ r + 4`), though `nlinarith` doesn't need the `r + 4` hint if given `r + 1 ≥ 0` and `2 - r ≥ 0`, or you can leave it as-is since it is fast and explicit.

---

### 3. Quick Check on Scope vs. the Main Contraction File

One architectural question to verify against your `Protolemmata/BARS/bgp-unknown-variance-lemmas.md` spec before locking the file:

* **UV-1a / UV-1b** bound the scale parts purely in terms of $r = \sigma_0 / \sigma$. Does the main contraction theorem expect a small bridge lemma in this file converting $\vert{}\sigma - \sigma_0\vert{} \le \sigma_0 \varepsilon_n$ (or $\vert{}\sigma^2 - \sigma_0^2\vert{}$) into $r \in [1/2, 2]$ and $(r - 1)^2 \le 4 (\sigma - \sigma_0)^2 / \sigma_0^2$, or is the KL neighborhood parameterized directly in $r = \sigma_0 / \sigma$?

If the neighborhood is parameterized directly in $r$ (or that bridge lives in the main file), then dropping `hr0'` and adding `eventually_priorMassRadius_pos` (and/or `exp_neg_n_radius_sq_le_radius`) makes this 100% ready for promotion.