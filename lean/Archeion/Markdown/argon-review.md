The five modules form a complete, proof-closed formalization of the unknown-variance contraction extension (`UV-1` through `UV-6`) modulo the three cited general contraction axioms in the final file.

---

### 1. Axiom Boundary & Statistical Soundness (`BgpContraction.lean`

)

When a formalization is conditional on cited axioms (`fixed_constant_contraction`, `gvdv07_theorem1`, `gvdv07_lemma1`), the primary soundness risk is a mismatch between the concrete definitions (`dataLaw`, `likelihood`, `posterior`, `klBall`) and the hypotheses of the axioms.

* **Model-to-Likelihood Coupling**: Neither `fixed_constant_contraction` nor `gvdv07_theorem1`/`gvdv07_lemma1` takes an abstract dominated experiment; instead, they are specialized directly to `dataLaw ev n` and `posterior ev prior n` (which integrates `likelihood ev n`). This is sound because `hprior_pos` (`∀ᵐ θ ∂prior, 0 < θ.2`) and `hS_pos` (`∀ θ ∈ S, 0 < θ.2`) force `(θ.2 ^ 2).toNNReal ≠ 0` on the support of the prior, ensuring `likelihood ev n θ` is genuinely the Radon–Nikodym density of `dataLaw ev n θ` with respect to Lebesgue product measure.


* **Dead Hypothesis / Unused `hsieve` in `posterior_contraction_unknown_variance_gvdv07**`:
On line 699 of the fourth file, you prove:


```lean
have hsieve : ∀ n, closure (sieve (B n) smin smax (smin * √ηlo) (smax * √ηhi)) ⊆ S := ...

```


However, `hsieve` is **never used** anywhere in `posterior_contraction_unknown_variance_gvdv07` (unlike in `posterior_contraction_unknown_variance`, where it restricts `htests` to the sieve).


* *Why it still holds mathematically*: In `gvdv07_theorem1`, `htests` only provides tests for centres $\theta_1 \in S$, whereas `hentropy` asks for an external cover $N \subseteq E \times \mathbb{R}$ of $\{ \theta \in \Theta_n \mid d_n(\theta, \theta_0) < r \}$ at radius $\xi r / 4$ without requiring $\Theta_n \subseteq S$. Because `hS : prior Sᶜ = 0` is a hypothesis of `gvdv07_theorem1`, one can intersect $\Theta_n$ with $S$ up to a prior-null set and then pass from the external $\xi r / 4$-cover to an internal $\xi r / 2$-cover with centres inside $\Theta_n \cap S \subseteq S$.


* *Fix*: Either remove the unused `hsieve` from `posterior_contraction_unknown_variance_gvdv07` or add `(hΘS : ∀ n, Θ n ⊆ S)` to `gvdv07_theorem1` and pass `hsieve`.




* **Rate Margin Accounting (`144 / σ₀²` vs `72 / σ₀²`)**:
Since the slab radius $\rho = \epsilon_n$ maps to a KL-ball radius of $\epsilon_n' = \frac{6}{\sigma_0}\epsilon_n$ via `slab_subset_klBall` (where $6^2 = 36 > 18 + 16 = 34$ absorbs both KL and $V_2$ scale/location constants in `klDiv_le_and_evariance_le_of_slab`), the factor of $(\sigma_0 / 6)^{-2} = 36 / \sigma_0^2$ multiplies the excess exponent: $4 \times 36 / \sigma_0^2 = 144 / \sigma_0^2$ in `posterior_contraction_unknown_variance` (for the $C + 4$ remainder condition) and $2 \times 36 / \sigma_0^2 = 72 / \sigma_0^2$ in `posterior_contraction_unknown_variance_gvdv07` (for the $e^{2 n (\epsilon_n')^2}$ factor in GvdV07 Lemma 1 at $k = 2$). Both constants are tight and consistent.



---

### 2. Mathematical & Proof Architecture Highlights

* Asymmetric Centring in the Scale Tests (`tests_scale_up_case` vs `tests_scale_down_case`):
Centring the $\chi^2$ statistic at the true mean $f_0$ when $\sigma_1 > \sigma_0 + r/2$ and at the alternative centre $f_1$ when $\sigma_1 < \sigma_0 - r/2$ is the cleanest part of the test construction:


* In `tests_scale_up_case`, the null error is an upper tail at $f_0$ (zero non-centrality) and the alternative error is a lower tail where `pi_gaussian_sq_lower_tail` drops the non-centrality $\sum (f(x_i) - f_0(x_i))^2 \ge 0$ unconditionally.


* In `tests_scale_down_case`, centring at $f_1$ makes the null error a lower tail (again dropping $\sum (f_0(x_i) - f_1(x_i))^2 \ge 0$ for free) while bounding the alternative's upper-tail non-centrality by $U \le n(r/8)^2$ on the $\xi r$-ball around $f_1$. Had it been centred at $f_0$, the non-centrality could be arbitrarily large when $\Vert{}f_1 - f_0\Vert{}_n \le r/2$ wasn't assumed or would be $16\times$ larger ($n(r/2)^2$), eating the scale gap $\frac{3}{8}\sigma_0 r$.




* Elementary Non-Central $\chi^2$ Upper Tail (`BgpGaussianTails.lean`):
In `pi_gaussian_sq_upper_tail`, using `one_div_sqrt_one_sub_le_exp` ($(1 - x)^{-1/2} \le \exp(\frac{x/2}{1-x})$) puts the determinant factor onto the exact same $(1 - 2\kappa v)$ denominator as the non-centrality term $\frac{\kappa(m_i - c_i)^2}{1 - 2\kappa v}$, collapsing the product MGF bound into $\exp(\frac{\kappa A}{1 - x})$ before applying $\frac{1}{1-x} \le 1 + 2x$ on $[0, 1/4]$.


* Internal-to-External Cover Conversion (`eventually_sieve_cover` & `BgpCovering.lean`):
Because `isCover_image2_smul` requires a norm bound `∀ g ∈ N₁, ‖g‖ ≤ R` on the cover centres $N_1$, `eventually_sieve_cover` converts the external covering number hypothesis `h61` at radius $\frac{\xi}{12 s_{\max}}\bar{\epsilon}_n$ into an *internal* cover $S_1 \subseteq B_n$ at radius $2 \cdot \frac{\xi}{12 s_{\max}}\bar{\epsilon}_n = \frac{\xi}{6 s_{\max}}\bar{\epsilon}_n$ via `coveringNumber_two_mul_le_externalCoveringNumber`. That guarantees $S_1 \subseteq B_n$, so `hRn : ∀ g ∈ B n, ‖g‖ ≤ n^p` applies directly to the centres without an extra $+\epsilon_1$ norm slack.


* `ℝ≥0∞`-Valued `posterior` (`BgpContraction.lean`):
Keeping `posterior` in `ℝ≥0∞` and stating contraction via `∫⁻ y, posterior ... ∂dataLaw` lets `posterior_le_add` use `ENNReal.add_div` and `ENNReal.div_le_div_right` unconditionally, avoiding a.e.-positivity and finiteness proofs for the marginal likelihood denominator $\int \text{likelihood}\,d\text{prior}$.



---

### 3. Lean 4 / Mathlib Craft Suggestions

1. Replace `obtain ⟨w, hW⟩ : ∃ w, expr = w := ⟨_, rfl⟩` with `generalize`:
Across `BgpGaussianTails.lean` (lines 56, 68–69, 179, 181, 214–225, 297, 316) and `BgpContraction.lean` (lines 281–284, 369, 439, 667–673), you use the idiom:


```lean
obtain ⟨w, hW⟩ : ∃ w, 1 - 2 * a * v = w := ⟨_, rfl⟩
rw [hW]

```


to introduce an opaque variable `w` that `field_simp` and `ring` will not unfold (unlike `set w := ...`). You can write this more concisely using Lean's `generalize` tactic:
```lean
generalize hW : 1 - 2 * a * v = w

```


which introduces `w`, binds `hW : 1 - 2 * a * v = w`, and rewrites the target in one step (and `generalize hW : 1 - 2 * a * v = w at hG ⊢` rewrites hypotheses simultaneously).
2. Quadratic vs. Quartic Moment Helper in `BgpGaussianKL.lean`:
In `gaussianLogRatio_eq`, `integral_gaussianLogRatio`, and `integrable_gaussianLogRatio`, the quadratic log-ratio is padded with `+ 0 * (x - μ₀) ^ 3 + 0 * (x - μ₀) ^ 4` to reuse `integral_quartic_gaussianReal_zero` and `integrable_quartic_gaussianReal_zero`. This works cleanly as a low-LOC shortcut; if `gaussianLogRatio_eq` is meant to be read externally as a normal-form lemma, a 4-line `integral_quadratic_gaussianReal_zero` wrapper (specializing `a₃ = 0, a₄ = 0`) would keep the statement of `gaussianLogRatio_eq` free of the dummy cubic and quartic terms.


3. Standalone API in `BgpCovering.lean`:
The cardinal-level lemmas `coveringNumber_image2_smul_le`, `coveringNumber_image2_smul_le_half`, `externalCoveringNumber_image2_smul_le`, and `coveringNumber_Icc_le` in `BgpCovering.lean` are not called in `BgpContraction.lean` (which works directly with explicit `Finset` covers via `isCover_image2_smul` and `exists_isCover_Icc`). They are worth keeping as the clean `ENat`-level covering number API for `UV-4a`, as they match the prose statement of the lemma.