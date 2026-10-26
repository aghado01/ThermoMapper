# BGP-4 — contraction with an unknown noise variance

Status: protolemma, 2026-10-23. Written at the user's request to complete an item TWCD2025 leave open.
Companion to [`bgp-twcd-lemmas.md`](bgp-twcd-lemmas.md) (BGP-1 to BGP-3).

TWCD2025's contraction theorems (Theorems 3.1–3.3) assume the noise standard deviation σ is known. Their
discussion section (§6, "Observation variance") says the theory could be extended to an unknown σ by putting a
prior on it, following van der Vaart & van Zanten, but they leave it there. Their Remark A.1 tracks how the rate's
constant depends on σ/‖f*‖ when σ is known. The implementation has since gone further than the theory:
`GpNoiseModel.UnknownScale` (commit `e8bb73c`) infers both the signal scale and the noise. This note states the
extension and splits it into pieces that are cited, pieces that need a proof, and pieces worth formalizing.

**Citation boundary.** The following are given and cited, not re-proved:
- **From TWCD2025:** Theorems 3.1–3.3, Assumptions 1–2 ((A1)–(A3)), the three framework inequalities (59)–(61)
  of Appendix A.1, and Remark A.1.
- **The general theorem** for independent, non-identically distributed observations (Ghosal & van der Vaart 2007,
  Theorem 4), which the paper's own framework rests on.
- **Hellinger tests** for independent observations (Le Cam; Birgé).

**Citation check (2026-10-25), recorded rather than applied silently.** The citation boundary is unchanged.
- **The printed GvdV07 Theorem 4 does not state the fixed-constant conclusion.** It has a single rate εₙ and concludes
  for every Mₙ → ∞. After it, the paper remarks that the Hellinger dₙ "can be replaced by any other distance dₙ that
  satisfies (3.2)–(3.3) and for which the conclusion of Lemma 2 holds". The composition with UV-5's tests rests on
  that remark, through Theorem 1 and Lemma 1.
- **UV-6's fixed constant C is the fixed-M, two-rate form.** That is Ghosal & van der Vaart (2001, Theorem 2.1), the
  two-rate version of Ghosal, Ghosh & van der Vaart (2000, Theorem 2.1). van der Vaart & van Zanten (2009, §3) name it
  as such, and TWCD2025's Appendix A.1 and Theorem 3.1 rely on it.
- **For independent, non-identically distributed observations the fixed-M form is not printed as a theorem.**
  van der Vaart & van Zanten (2008, Theorem 3.3) state it for fixed-design GP regression with an unknown σ on [a, b],
  with the tests left to the reader.

The Lean module takes both routes (UV-6). The Hellinger tests are not needed: UV-5 builds tests for
‖f − f′‖ₙ + |σ − σ′| directly.

**Prior art.** van der Vaart & van Zanten (2008, Theorem 3.3) is the unknown-σ result for a single GP prior W, with a
σ prior on [a, b] whose density is bounded below. Their proof says only that the known-σ argument "can be extended to
the present case of an unknown scale". BGP-4 is the extension to TWCD2025's hierarchical bandwidth prior and to the
implementation's amplitude and noise-ratio parameterization, with that step written out (UV-5).

**The textbook treatment.** Gaussian regression with an unknown error variance is treated by Ghosal & van der Vaart
(2017, *Fundamentals of Nonparametric Bayesian Inference*). What follows is that argument specialized to TWCD2025's
hierarchical GP prior and to the implementation's parameterization.

The Lean blocks are sketches with scoped imports: statements with proof strategies, not compiled. "Lemma-tier"
means complete up to citable `sorry`s. "Enthymeme" means a forward boundary not yet discharged.

---

## The setting

**Data (fixed design, as in Theorem 3.1).** Yᵢ = f*(Xᵢ) + wᵢ, with wᵢ ~ N(0, σ₀²) independent, and σ₀ > 0 unknown.

**The implementation's prior** (`GpNoiseModel.UnknownScale`):
- the bandwidth t has the empirical-Bayes prior p(t), which satisfies (A3) (Prop. 4.4);
- given t and s, the regression function is f = s·Wᵗ, with Wᵗ ~ GP(0, hₜ) of unit amplitude;
- the noise variance is σ² = s²η;
- s² has the improper prior 1/s², and log η is uniform on [log η_lo, log η_hi].

**The theory's prior [P1]** is the proper restriction:
- s has a continuous density, positive on [s_min, s_max] and zero outside. For example, 1/s² restricted to the box,
  renormalized;
- η is as implemented, so σ lies in the compact range [σ_min, σ_max] = [s_min √η_lo, s_max √η_hi];
- t, s and η are independent a priori;
- **the one requirement on the truth:** σ₀ ∈ (σ_min, σ_max). There is no "true s": s is a prior scale for f, and
  any interior value serves.

**Why compact support for s.** It is what the standard sieve argument needs (UV-3, UV-4). A heavy-tailed amplitude
prior is open (UV-8). The noise variance alone could take inverse-gamma tails, because its entropy is counted on a
log scale (UV-4). The amplitude cannot, because f's entropy is counted in sup norm, where s multiplies.

---

## UV-1 — KL divergence and variance for one Gaussian observation (Lean target)

For one observation with mean μ₀ and sd σ₀ under the truth, and μ, σ under the model, write δ = μ − μ₀ and
r = σ₀/σ:

```text
K(μ,σ) := KL(N(μ₀,σ₀²) ‖ N(μ,σ²)) = −log r + (r² − 1)/2 + δ² r²/(2σ₀²)
V(μ,σ) := Var_{N(μ₀,σ₀²)}[log p₀/p]   = (r² − 1)²/2 + δ² r⁴/σ₀²
```

**Derivation.** Write y = μ₀ + σ₀z with z ~ N(0, 1). The log ratio is
log(σ/σ₀) + (σ₀z − δ)²/(2σ²) − z²/2. Take its mean and variance, using E z = E z³ = 0 and Var z² = 2.

**Local quadratic bound.** For r ∈ [1/2, 2] and |δ| ≤ σ₀,
- K ≤ (5/2)(r − 1)² + 2δ²/σ₀², because g(r) = −log r + (r² − 1)/2 has g(1) = g′(1) = 0 and g″(r) = 1 + 1/r² ≤ 5 there;
- V ≤ (9/2)(r − 1)² + 16δ²/σ₀².

So with σ₀ fixed, the sup-norm and |σ − σ₀| neighbourhood of radius cε lies inside the KL neighbourhood
{K ≤ ε², V ≤ ε²}, averaged over the n observations. That is the inclusion the prior-mass step (UV-2) needs.

**Tier: Lemma.** **Discharged 2026-10-23** in `Enthymemata/BgpUnknownVariance.lean`:
- `BGP.gaussian_kl_scale_part_le` holds for every r ≥ 1/2, since the upper end is not needed. The proof skips
  Taylor: from 1 − 1/r ≤ log r, the gap is (r − 1)²(2r − 1)/r ≥ 0.
- `BGP.gaussian_v_scale_part_le` holds on [1/2, 2].
- `BGP.scale_ratio_bridge` (UV-1d, added after review): |σ − σ₀| ≤ σ₀/2 puts r = σ₀/σ in [1/2, 2] with
  (r − 1)² ≤ 4(σ − σ₀)²/σ₀². The neighbourhood is stated in |σ − σ₀|, so this lands it where UV-1a and UV-1b apply.

**UV-1c, the closed forms: discharged 2026-10-24** in `Enthymemata/BgpGaussianKL.lean` (PROOF-CLOSED), in the
`(μ₀, σ₀, μ, σ)` parameterization above, with variances entered as `(σ ^ 2).toNNReal`:
- `BGP.klDiv_gaussianReal`: K against Mathlib's `InformationTheory.klDiv`, as `ENNReal.ofReal` of the closed form.
- `BGP.variance_llr_gaussianReal`: V against Mathlib's `variance` of the log-likelihood ratio `llr`;
  `BGP.evariance_llr_gaussianReal` is the same identity in `ℝ≥0∞`.
- The pinned Mathlib has no Gaussian fourth moment. `BGP.moments_gaussianReal_zero` reads E y, E y², E y³, E y⁴ of
  N(0, v) (0, v, 0, 3v²) off the moment-generating function exp(v t²/2), via `iteratedDeriv_mgf_zero`. No gap
  remains.
- `BGP.klDiv_le_and_evariance_le_of_slab` is the inclusion above, assembled from UV-1a, UV-1b, UV-1c and UV-1d: the slab
  |μ − μ₀| ≤ σ₀ε/6, |σ − σ₀| ≤ σ₀ε/6 with ε ≤ 3 gives K ≤ ε²/3 and V ≤ 17ε²/18. So c = σ₀/6.

Notes:
- The bound is elementary real analysis. It is the best Lean target here.

```lean
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.Calculus.Taylor

/-- UV-1a. The variance part of the Gaussian KL is locally quadratic in r = σ₀/σ.
    Strategy: g(r) = −log r + (r² − 1)/2 has g 1 = 0, g' 1 = 0 and g'' r = 1 + 1/r² ≤ 5 on [1/2, 2];
    Taylor with the Lagrange remainder. -/
theorem gaussian_kl_scale_part_le (r : ℝ) (hr : r ∈ Set.Icc (1 / 2 : ℝ) 2) :
    -Real.log r + (r ^ 2 - 1) / 2 ≤ (5 / 2) * (r - 1) ^ 2 := by
  sorry

/-- UV-1b. The variance of the log ratio's scale part: ((r² − 1)²)/2 ≤ (9/2)(r − 1)² on [1/2, 2],
    since (r + 1)² ≤ 9 there. Pure algebra. -/
theorem gaussian_v_scale_part_le (r : ℝ) (hr : r ∈ Set.Icc (1 / 2 : ℝ) 2) :
    (r ^ 2 - 1) ^ 2 / 2 ≤ (9 / 2) * (r - 1) ^ 2 := by
  sorry
```

---

## UV-2 — prior mass of a KL neighbourhood (59′)

**Claim.** For the paper's εₙ ((62): εₙ = C̄₁ n^{−s/(2s+ϱ)} (log n)^{k₁}) and n large,

```text
Π( ‖f − f*‖∞ ≤ εₙ , |σ − σ₀| ≤ εₙ ) ≥ exp(−C nεₙ²)
```

By UV-1, this gives the KL-neighbourhood mass the general theorem needs.

**Proof sketch.**
1. **Pick an interior s\*.** Choose s* with σ₀/s* in the interior of [√η_lo, √η_hi], and a compact neighbourhood
   S₀ ⊂ (s_min, s_max) of it on which σ₀/s stays interior.
2. **Rescaling, from the authors' own Remark A.1.** For s ∈ S₀,
   Π(‖s·Wᵗ − f*‖∞ ≤ ε | t, s) = Π(‖Wᵗ − f*/s‖∞ ≤ ε/s | t). Now f*/s satisfies (A2) with constants (ν₁/s, ν₂/s²).
   Remark A.1 records that the paper's Part I bound for (59) depends on f* only through (ν₁, ν₂), and it depends on
   the radius polynomially. So the (59) bound holds for f*/s at radius εₙ/s with constants uniform over the compact
   S₀, which is bounded away from 0.
3. **The noise interval.** Given s ∈ S₀, σ = s√η has a density bounded below by some c > 0 near σ₀, because log η is
   uniform with σ₀/s interior. So Π(|σ − σ₀| ≤ εₙ | s) ≥ c εₙ.
4. **Combine.** By the independence of t, s and η, the joint mass is at least
   π_s(S₀) · c εₙ · exp(−C′ nεₙ²). Since log(1/εₙ) = O(log n) = o(nεₙ²) (here nεₙ² ≍ n^{ϱ/(2s+ϱ)}(log n)^{2k₁}),
   the factor c εₙ costs at most one more nεₙ² in the exponent.

**Tier: Lemma**, up to the citable (59) bound and Remark A.1.
- **To confirm against the paper:** that Part I's constants are continuous, or at least locally bounded, in
  (ν₁, ν₂) and in the radius scaling. Remark A.1 says so in words; the proof's constants C̄₁ and C̄₂ should be read
  once to make sure.
  - **Confirmed 2026-10-25, by a stronger route than continuity**, against TWCD2025 Appendix A.1, Part I (eqs. 64–74).
  - Part I uses f* only through (A2): some q_t ∈ H̃_t with ‖q_t − f*‖∞ ≤ ν₁t^{s/2} and ‖q_t‖² ≤ ν₂t^{−ϱ/2}.
    Here K₃ = ν₂ + C₄, C = 1/ν₁, and K₄, C̄′₁ and C̄₁ are explicit in them.
  - (A2) is monotone in (ν₁, ν₂). f*/s satisfies it with (ν₁/s, ν₂/s²) via q_t/s, so for every s in the box it
    satisfies it with the single pair (ν₁/s_min, ν₂/s_min²).
  - Part I run with that pair therefore gives (59) for every f*/s at once: the same rate constant and the same
    threshold in n.
  - It holds at radius εₙ/s once C₀ ≥ s_max·C̄₁(ν₁/s_min, ν₂/s_min²), with C₁ ≤ 1 in the units of `h59`.
- **Lean-worthy pieces:** the scaling identity (step 2) and the absorption of log(1/εₙ) (step 4). **Both discharged
  2026-10-23** in `Enthymemata/BgpUnknownVariance.lean`:
  - `BGP.map_smul_closedBall`;
  - `BGP.log_inv_radius_le_n_radius_sq`, for k ≥ 0, which covers the paper's k₁ = (1 + D)/(2 + ϱ/s) > 0. Mathlib's
    `Pow.Asymptotics` module is not in the local cache, so the proof is elementary instead:
    - log n ≤ nᵇ/b, from log y ≤ y − 1;
    - an explicit threshold N = max(⌈e⌉, ⌈M^(1/b)⌉), with b = (1 − 2a)/2.
  - companions, added after review. `Real.log` is 0 at nonpositive arguments, so UV-2b alone holds vacuously
    wherever εₙ ≤ 0 and hands the caller no positivity:
    - `BGP.eventually_priorMassRadius_pos`: εₙ > 0 eventually;
    - `BGP.exp_neg_n_radius_sq_le_radius`: exp(−nεₙ²) ≤ εₙ eventually, the form the prior-mass step consumes;
    - `BGP.log_inv_mul_radius_le`: log(1/(c₁εₙ)) ≤ c₂nεₙ² eventually, for the slab width c₁εₙ and the exponent
      budget c₂nεₙ². It applies UV-2b at c·C with c = min(c₁, 1, c₂).
- **(59′) assembled 2026-10-24** in `Enthymemata/BgpUnknownVarianceTheorem.lean`, with no gap beyond the cited (59):
  - step 1 is `BGP.exists_amplitude_window`: a window [s₁, s₂] with s₂√η_lo ≤ σ₀ < s₁√η_hi;
  - step 3 is `BGP.noise_slab_mass`: mass at least εₙ/(σ₀ log(η_hi/η_lo)), since the η-interval has log-length
    2 log(1 + εₙ/σ₀) ≥ εₙ/σ₀;
  - steps 2 and 4 are `BGP.prior_slab_mass` (by independence the slab's section at s is a rectangle, and UV-2a turns
    its W-side into the rescaled ball) and `BGP.eventually_prior_mass`. The exponent is C₁ + δ for any δ > 0, so the
    slab factor costs an arbitrarily small part of nεₙ², not a whole one.
  - The "to confirm" item above is the named hypothesis `h59`: (59) uniform over the whole box [s_min, s_max].
  - Pow.Asymptotics is now cached; `BGP.tendsto_priorMassRadius_zero` uses it.

```lean
import Mathlib.MeasureTheory.Measure.Map
import Mathlib.Analysis.Normed.Module.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics

open MeasureTheory Filter

/-- UV-2a. Prior mass of a ball under an amplitude rescaling is the mass of the rescaled ball.
    Strategy: the preimage of `closedBall f ε` under `w ↦ s • w` is `closedBall (s⁻¹ • f) (ε / s)`. -/
theorem map_smul_closedBall {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [MeasurableSpace E]
    [BorelSpace E] (μ : Measure E) {s : ℝ} (hs : 0 < s) (f : E) (ε : ℝ) :
    (μ.map (fun w => s • w)) (Metric.closedBall f ε) = μ (Metric.closedBall (s⁻¹ • f) (ε / s)) := by
  sorry

/-- UV-2b. The noise interval's mass is absorbed: for εₙ = C n^(−a) (log n)^k with 0 < a < 1/2,
    log(1/εₙ) ≤ n εₙ² eventually. Strategy: n εₙ² = C² n^(1−2a) (log n)^(2k) dominates a·log n. -/
theorem log_inv_eps_le_n_eps_sq (C a k : ℝ) (hC : 0 < C) (ha : 0 < a) (ha' : a < 1 / 2) :
    ∀ᶠ n : ℕ in atTop,
      Real.log (1 / (C * (n : ℝ) ^ (-a) * Real.log n ^ k)) ≤ n * (C * (n : ℝ) ^ (-a) * Real.log n ^ k) ^ 2 := by
  sorry
```

---

## UV-3 — sieve and remaining mass (60′)

**The sieve.** B′ₙ = {s·g : s ∈ [s_min, s_max], g ∈ Bₙ} × [σ_min, σ_max], with Bₙ the paper's sieve.

**Claim.** Π((f, σ) ∉ B′ₙ) ≤ exp(−4nεₙ²).

**Proof.** f = s·Wᵗ lies in s·Bₙ exactly when Wᵗ ∈ Bₙ, so Π(f ∉ B′ₙ) ≤ Π(Wᵗ ∉ Bₙ). That is the paper's (60)
unchanged. σ always lies in its compact range. Nothing new is needed.

**Tier: Lemma**, citing (60). This is where compact support for s pays: with an unbounded s, the sieve would need
an s-cutoff sₙ with Π(s > sₙ) ≤ exp(−4nεₙ²). For polynomial tails sₙ grows exponentially in nεₙ², and UV-4 then
fails (UV-8).

**Discharged 2026-10-24** in `Enthymemata/BgpUnknownVarianceTheorem.lean`: `BGP.prior_compl_sieve_le` bounds the
prior mass outside the closure of B′ₙ by the GP mass outside Bₙ. The closure is closed, so measurable, which the
general theorem asks of its sieve.

**Correction (2026-10-24, found while composing UV-6): the exponent 4 is not enough by itself.** The general theorem
compares the remaining mass with the prior mass of the KL neighbourhood, not with a fixed exponent. GvdV07's condition
asks Π(Θ ∖ Θₙ)/Π(B*ₙ(εₙ)) = o(e^{−2nεₙ²}). The fixed-constant form asks Π(Θ ∖ Θₙ) ≤ e^{−(C+4)nεₙ²} against
Π(B*ₙ(εₙ)) ≥ e^{−Cnεₙ²}. Both count in units of the KL radius.
- **(59′)'s C is not 1.** It comes from (59) for f*/s at radius εₙ/s, with s-dependent constants (Remark A.1), plus
  the slab factor.
- **The KL radius is the slab radius over c** (UV-1; c = σ₀/6 in the Lean inclusion).
- **So the composition needs a margin.** Write (59) as e^{−C₁nεₙ²} and (60) as e^{−C₂nεₙ²}, both in units of the slab
  radius. The composition needs C₂ > C₁ + 4/c² = C₁ + 144/σ₀². This is the Lean statement's `hC₂`. A sharper UV-1
  constant lowers the margin, but a σ₀-dependent margin is inherent.
- **The paper's (60) at a fixed 4 suffices only if that margin happens to hold.** The sieve is a proof device and can
  be enlarged. Borell's inequality gives sieves of the form Mₙℍ₁ + εₙ𝔹₁ at any remaining-mass exponent, at the
  cost of the entropy constant, as in van der Vaart & van Zanten (2008, Theorem 2.1).
- **To confirm against the paper:** that TWCD2025's B_{Nₙ,rₙ,δₙ,εₙ} gives (60) at any exponent, with (61) keeping its
  order. (60′) should then read Π((f, σ) ∉ B′ₙ) ≤ Π(Wᵗ ∉ Bₙ) ≤ exp(−C₂nεₙ²), with C₂ chosen for the margin.
  - **Confirmed 2026-10-25** against Appendix A.1, Parts II–III.
  - The exponent 4 comes only from the free constants in eq. (80): rₙ^ϱ = (8/K₂)nεₙ² and
    Nₙ² = max{32, 128C₄/K₂}nεₙ²(log)^{1+D}.
  - The bound (81), P[f^t ∉ Bₙ] ≲ rₙ^{…}e^{−K₂rₙ^ϱ} + e^{−Nₙ²/8}, gives exp(−C₂nεₙ²) once rₙ^ϱ and Nₙ² are scaled
    up by C₂/4.
  - Part III's entropy (86) then grows only through the constant C̄₂ of ε̄ₙ, which the paper chooses large, so (61)
    keeps its order.
  - van der Vaart & van Zanten (2008, Theorem 2.1) state the same tunability for their sieves ("any C > 1").
  - In Lean this is `BGP.posterior_contraction_unknown_variance_of_tunable_sieve`.

---

## UV-4 — entropy (61′)

**Claim.**

```text
log N(ε̄ₙ, B′ₙ, d) ≤ log N(ε̄ₙ/(2 s_max), Bₙ, ‖·‖∞) + log N(ε̄ₙ/(2Rₙ), [s_min, s_max], |·|)
                   + log N(c ε̄ₙ, [σ_min, σ_max], |·|) ≤ C n ε̄ₙ²
```

Here Rₙ = sup_{g∈Bₙ} ‖g‖∞, and d is the metric the tests use (UV-5). On the compact σ range, d is Lipschitz in
‖f − f′‖∞ + |σ − σ′|.

**Proof sketch.**
- **Scaled sets.** s·g − s′·g′ = s(g − g′) + (s − s′)g′, so a cover of Bₙ at radius ε/(2s_max) and a cover of the
  s-interval at radius ε/(2Rₙ) cover the scaled set.
- **The first term.** The paper's (61) is polylogarithmic in 1/radius, so a constant factor in the radius changes
  it by a constant factor.
- **The two interval terms** are O(log(Rₙ/ε̄ₙ)) and O(log(1/ε̄ₙ)): O(log n) when Rₙ grows at most polynomially in
  n, and log n = o(nε̄ₙ²).

**Tier: Lemma**, up to the citable (61).
- **To confirm against the paper:** that Rₙ grows at most polynomially. Sieves of the van der Vaart–van Zanten form
  Mₙ·ℍ₁ + εₙ·𝔹₁ have Rₙ ≲ Mₙ ≍ √n εₙ, and the paper's B_{Nₙ,rₙ,δₙ,εₙ} (proof of Theorem 3.1) should be checked.
  - **Confirmed 2026-10-25** against eq. (75): Bₙ = (Nₙ(rₙ/δₙ)^{D/2}H̃¹_{rₙ⁻²} + εₙ𝔹₁) ∪ ⋃_{t>δₙ⁻²}(NₙH̃¹_t + εₙ𝔹₁).
  - The kernel has h_t(x, x) = h(0) = 1, so the RKHS unit balls have sup norm at most 1, and Rₙ ≤ Nₙ(rₙ/δₙ)^{D/2} + εₙ.
  - With eq. (80) and δₙ = εₙ/(2√Dτ_hNₙ), the ratio rₙ/δₙ = C(Nₙ/εₙ)(nεₙ²)^{1/ϱ} is polynomial in n, and so is Nₙ.
- **The heavier-tailed variant.** With inverse-gamma tails on σ² instead of a compact range, cover σ on the log
  scale: the location-scale Hellinger distance is controlled by |log σ − log σ′|. An exponentially wide range then
  costs only log(log σ⁺ / ε̄ₙ) = O(log n). This is the textbook route for an unknown variance.
- **Lean-worthy:** the covering bound for a scaled family and for an interval (metric-space covering numbers;
  check the Mathlib covering-number API at our pinned version).

**UV-4a: discharged 2026-10-24** in `Enthymemata/BgpCovering.lean` (PROOF-CLOSED), over any seminormed space on ℝ, so
the sup-norm space fits:
- **The cover.** `BGP.isCover_image2_smul` covers the scaled family at any radius ε ≥ b ε₁ + R ε₂. It needs 0 ≤ a,
  covers N₁ of B at radius ε₁ and N₂ of [a, b] at radius ε₂, and N₁'s points bounded by R.
- **Covering numbers.** `BGP.coveringNumber_image2_smul_le` gives N(ε) ≤ N(ε₂, [a, b]) · N(ε₁, B), with the brief's
  radii ε/(2R) and ε/(2b) in `…_le_half`. Mathlib's internal covering numbers are used, since the product of internal
  covers lies in the scaled family. `BGP.externalCoveringNumber_image2_smul_le` is the external form, through
  Mathlib's coveringNumber (2δ) ≤ externalCoveringNumber δ.
- **The interval.** `BGP.exists_isCover_Icc` and `BGP.coveringNumber_Icc_le` give N(ε, [a, b]) ≤ ⌈(b − a)/(2ε)⌉ + 1
  for ε > 0.

**(61′) assembled** in `Enthymemata/BgpUnknownVarianceTheorem.lean`:
- `BGP.closure_sieve_subset_cover` covers the closed sieve in ‖f − f′‖ₙ + |σ − σ′|, which the sup norm dominates
  when the design evaluations have norm at most 1.
- `BGP.eventually_sieve_cover` counts the cover: at most exp(n(Lε̄ₙ)²) balls of radius ξLε̄ₙ/2.
- **(61) enters as a hypothesis**, at any fixed multiple of ε̄ₙ with a constant depending on the multiple. That is
  the "constant factor in the radius" remark above. It is taken for external covering numbers.
- **Rₙ ≤ nᵖ is the named hypothesis `hRadius`, not a `sorry`.** It is a property of the paper's sieve, which Lean does
  not construct. A `sorry` there would assert it for an arbitrary sieve, where it is false.

---

## UV-5 — tests and the metric

- **The general theorem.** Ghosal & van der Vaart (2007, Theorem 4) gives contraction for independent, non-identically
  distributed observations. Its metric is the root average squared Hellinger distance dₙ between the observations'
  laws N(f(Xᵢ), σ²) and N(f*(Xᵢ), σ₀²). It needs KL-neighbourhood mass (UV-1, UV-2), a sieve (UV-3) and entropy in
  dₙ (UV-4, via the Lipschitz bound on the compact σ range).
- **From dₙ to the quantities we report.** Contraction in dₙ must be turned into contraction of ‖f − f*‖ₙ and
  |σ − σ₀|. For normal location–scale laws with σ in a compact range, H² ≥ c·min(1, δ²/σ₀² + (r − 1)²). The averaging
  over i needs care, since the average of the minimum is not the minimum of the average. The clean route is the one
  the known-σ regression literature uses: build the tests directly for ‖f − f′‖ₙ + |σ − σ′|, with likelihood-ratio
  tests for Gaussian means and a χ² test for the scale, as in the textbook treatment.
- **Tier: Enthymeme** until that step is written out. It is standard, but it is the one place the extension is more
  than bookkeeping. **Not a Lean target now.**
- **Proved 2026-10-25:** `BGP.exists_tests` in `Enthymemata/BgpUnknownVarianceTheorem.lean`. The Gaussian tail
  bounds are in `Enthymemata/BgpGaussianTails.lean`, which is PROOF-CLOSED.
  - **What the tests do.** They are for ‖f − f′‖ₙ + |σ − σ′| on the compact noise range E × [σ_min, σ_max] and work at
    every radius r > 0. Each separates the truth from the r/8-ball around any θ₁ at distance more than r, with both
    errors at most exp(−Knr²), where K = σ_min²/(1024σ_max⁴) and ξ = 1/8.
  - **‖f₁ − f*‖ₙ > r/2:** the likelihood-ratio statistic Σδᵢ(Yᵢ − f*(xᵢ)), with δ = f₁ − f*, against ‖δ‖²/2. Cauchy–
    Schwarz keeps the alternatives' mean above 3‖δ‖²/4.
  - **σ₁ > σ₀ + r/2:** reject when Σ(Yᵢ − f*(xᵢ))² ≥ nτ, with τ = σ₀² + 3σ₀r/8. Under the truth the sum is central.
    The alternatives need only its lower tail, where the means drop out.
  - **σ₁ < σ₀ − r/2:** reject when Σ(Yᵢ − f₁(xᵢ))² ≤ nτ, with τ = σ₀² − σ₀r/8. The truth needs only the lower tail,
    and the alternatives' means are within r/8 of f₁. Here r < 2σ₀.
  - **The centring is the point.** Each scale test is centred where its means do the least harm. A single centring at
    f₁ leaves a mean shift of order r under the truth, which can cancel the scale gap when r is comparable to σ₀.
    The sketch originally centred both scale tests at f₁, and that is corrected here.
  - **The tails** come from completing the square under `gaussianReal`, independence under `Measure.pi`, and
    Chernoff:
    - the linear tail exp(−t²/(2vΣcᵢ²));
    - the χ² upper tail exp(−(nτ − A)²/(16v·nτ)), with A = nv + Σ(mᵢ − cᵢ)²;
    - the χ² lower tail exp(−n(v − τ)²/(8v²)).

---

## UV-6 — the extended theorem

**Statement (to be proved from UV-1 to UV-5).** Under Assumptions 1–2 with factors (ϱ, s), and the proper prior
[P1] with σ₀ ∈ (σ_min, σ_max):

```text
Π( ‖f − f*‖ₙ + |σ − σ₀| > C ε̄ₙ | data ) → 0 in probability,   ε̄ₙ ≲ n^{−s/(2s+ϱ)} (log n)^{D+1}
```

- **The rate is the known-σ rate.** C depends on σ₀, on the s and η ranges, and (through Remark A.1 applied to
  f*/s*) on ‖f*‖/s*.
- **Random design (the analogue of Theorem 3.3)** follows with the truncation f_M, by the paper's own random-design
  argument (its Lemma A.1).
- **The posterior-mean estimator (the analogue of Theorem 3.2)** needs the paper's moment argument redone with σ
  integrated out: enthymeme.
- **What it says about σ.** The posterior of σ contracts at rate ε̄ₙ. That is slower than √n, and UV-7 is the
  sharper statement.

**Tier:** a theorem composed of cited and proved pieces. As an enthymema it would carry `sorry`s exactly at UV-5 and
the two "to confirm" checks.

**2026-10-25: no `sorry` left.** UV-5 is proved, and `Enthymemata/BgpUnknownVarianceTheorem.lean` has no `sorry`.
It rests on cited axioms, in two routes:
- **`BGP.posterior_contraction_unknown_variance_gvdv07`, on GvdV07 as printed.**
  - `BGP.gvdv07_theorem1` (Theorem 1) applies on the sieve at the entropy rate Lε̄ₙ. `BGP.gvdv07_lemma1` (Lemma 1)
    applies off the sieve at the prior-mass rate 6εₙ/σ₀. The two are added in Lean (`BGP.posterior_le_add`).
  - The conclusion is the printed one: for every Mₙ → ∞, the expected posterior mass of
    {‖f − f*‖ₙ + |σ − σ₀| ≥ Mₙε̄ₙ} tends to 0.
  - Lemma 1 needs the smaller margin C₂ > C₁ + 72/σ₀². The shells condition (2.5) follows from the prior mass, since
    εₙ ≤ ε̄ₙ.
  - 2026-10-26: `BGP.gvdv07_theorem1` takes `Θₙ ⊆ S`, the printed Θₙ ⊂ Θ. The closed sieves lie in the noise range
    `S`, so the hypothesis comes free. Without it the axiom leaned on an unchecked reduction to Θₙ ∩ S. An outside
    review (`Archeion/Markdown/argon-review.md`) found the proof of `Θₙ ⊆ S` unused.
- **`BGP.posterior_contraction_unknown_variance`, the protolemma's fixed-C statement,** on the axiom
  `BGP.fixed_constant_contraction`. This is the fixed-M two-rate form of the citation check above, under its
  corrected citation.
- **`BGP.posterior_contraction_unknown_variance_of_tunable_sieve`** is the fixed-C statement with (60)–(61) taken in
  the paper's tunable form, which removes the margin hypothesis.

`#print axioms` for each theorem should show only its own cited axioms and Lean's standard ones. The record of
2026-10-24 follows:

**Enthymema 2026-10-24:** `BGP.posterior_contraction_unknown_variance` in
`Enthymemata/BgpUnknownVarianceTheorem.lean`, `apologizing(1)`. `#print axioms` shows `sorryAx` (only through
`exists_tests`) and the cited `BGP.gvdv07_contraction`, since renamed `BGP.fixed_constant_contraction`. Every
supporting lemma uses only Lean's standard axioms.
- **Conclusion.** For some M, the expected posterior mass of {‖f − f*‖ₙ + |σ − σ₀| ≥ M ε̄ₙ} under the truth tends to 0.
  This is the general theorem's form; for a mass in [0, 1] it is equivalent to convergence in probability. The rate
  is ε̄ₙ = C̄ n^{−a}(log n)^{k̄}, the shape of (62), with the same a as εₙ and k ≤ k̄: here k₁ = (1+D)/(2+ϱ/s) ≤ D + 1.
- **The design.** Evaluations at the design points are continuous linear functionals of norm at most 1, so the sup
  norm dominates ‖·‖ₙ.
- **TWCD2025's (59)–(61) are hypotheses**, being facts about the paper's GP prior and sieve:
  - `h59` is (59) for f*/s at radius εₙ/s, uniform over [s_min, s_max]: UV-2's check, confirmed (UV-2).
  - `h60` is (60), at exponent C₂ with the margin `hC₂` (see the correction under UV-3).
  - `h61` is (61) (see UV-4).
  - `hRadius` is UV-4's check, confirmed (UV-4).

  The two checks stay hypotheses, now confirmed against the paper, for the reason given under UV-4.
- **[P1] enters by the properties the proof uses**, which [P1] has:
  - the joint law is the image of πs ⊗ (μW ⊗ πη) under (s, w, η) ↦ (s • w, s√η);
  - πs lives on [s_min, s_max] and gives positive mass to every subinterval, as a positive density does;
  - πη lives on [η_lo, η_hi] and gives every subinterval [a, b] at least log(b/a)/log(η_hi/η_lo), which the
    log-uniform law gives with equality.

  The window [s₁, s₂] of UV-2 step 1 is built from σ₀ ∈ (σ_min, σ_max) (`BGP.exists_amplitude_window`). The noise
  slab of step 3 has mass at least εₙ/(σ₀ log(η_hi/η_lo)) (`BGP.noise_slab_mass`).
- **The cited axiom** `BGP.gvdv07_contraction` takes GvdV07 Theorem 4 in the form TWCD2025 and van der Vaart &
  van Zanten use: two rates (εₙ ≤ ε̄ₙ), a fixed constant M, and a testing hypothesis for a general semimetric in
  place of the Hellinger one. Its hypotheses are stated on the strong side: measurable sieves and neighbourhood parts,
  tests at every radius, a global cover, and nεₙ² → ∞. **To confirm at review:** that the printed Theorem 4, or the
  paper's general testing theorem, states the fixed-M form rather than only "every Mₙ → ∞". The fixed-M form follows
  from the same proof under the prior-mass/remaining-mass conditions, as in Ghosal, Ghosh & van der Vaart (2000,
  Theorem 2.1).
  - **Resolved 2026-10-25:** it states only the "every Mₙ → ∞" form (see the citation check at the top). The axiom is
    renamed `BGP.fixed_constant_contraction` and re-cited. The printed route is the separate pair
    `BGP.gvdv07_theorem1` and `BGP.gvdv07_lemma1`.

---

## UV-7 — the σ Bernstein–von Mises pointer

de Jonge & van Zanten (2013), *Semiparametric Bernstein–von Mises for the error standard deviation* (Electron.
J. Statist. 7), give conditions for GP regression under which σ's marginal posterior is asymptotically normal. It is
centred at an efficient estimator, with variance σ₀²/(2n). Credible intervals for σ are then √n-consistent
frequentist confidence intervals.

- **Open:** whether its conditions hold for TWCD2025's empirical-Bayes bandwidth prior, and with an amplitude.
  **Enthymeme.**
- **A cheap numerical check before any proof:** the frequentist coverage of `BgpResult.NoiseVarianceInterval` over
  replicated simulated data sets should approach 95% as n grows. The fact harness can run that.

---

## UV-8 — the bridge to the implementation, and what stays open

- **Improper 1/s² against the compact box (bridged, no code change).** [P1] with s's density proportional to 1/s²
  on [s_min, s_max] gives exactly the implementation's posterior restricted to the box and renormalized. So a run's
  own draws measure the gap: the posterior mass of s² outside a chosen box. When that mass is negligible, the
  theorem's prior and ours give posteriors that cannot be told apart. Reporting it would be a small readout,
  computed from the retained draws.
- **Heavy-tailed amplitude priors (open).** UV-3 and UV-4 fail as stated. A different sieve, or a result on
  amplitude-rescaled GP priors, would be needed. **Enthymeme.**
- **The flat prior on the constant mean (a separate delta).** The implementation integrates a constant mean out
  under a flat prior (commit `f24d9d7`); the paper's GP is centred. A constant with a N(0, v) prior keeps the prior
  Gaussian (kernel hₜ + v), and the flat prior is its limit. This is a candidate BGP-5, not part of this note.
- **The empirical-Bayes bandwidth prior.** Proposition 4.4 concerns t alone and carries over unchanged.
- **Empirical-Bayes σ (the paper's other open direction).** A maximum-marginal-likelihood plug-in for σ would fall
  under Rousseau & Szabó (2017, Ann. Statist. 45) rather than under this argument.

---

## Lean-worthy pieces

| Piece | Statement | Difficulty | Status |
|---|---|---|---|
| UV-1a | −log r + (r² − 1)/2 ≤ (5/2)(r − 1)² for r ≥ 1/2 | Easy: from 1 − 1/r ≤ log r | **Proved** (`BGP.gaussian_kl_scale_part_le`) |
| UV-1b | (r² − 1)²/2 ≤ (9/2)(r − 1)² on [1/2, 2] | Easy: algebra | **Proved** (`BGP.gaussian_v_scale_part_le`) |
| UV-1c | Gaussian KL and log-ratio variance closed forms | Moderate: Gaussian moments in Mathlib | **Proved** (`BGP.klDiv_gaussianReal`, `BGP.variance_llr_gaussianReal`, `BGP.evariance_llr_gaussianReal`; fourth moment from the MGF, `BGP.moments_gaussianReal_zero`; slab inclusion `BGP.klDiv_le_and_evariance_le_of_slab`) |
| UV-1d | \|σ − σ₀\| ≤ σ₀/2 ⇒ r ∈ [1/2, 2], (r − 1)² ≤ 4(σ − σ₀)²/σ₀² | Easy: algebra | **Proved** (`BGP.scale_ratio_bridge`) |
| UV-2a | Ball mass under an amplitude rescaling | Easy: preimage of a ball under scaling | **Proved** (`BGP.map_smul_closedBall`) |
| UV-2b | log(1/εₙ) ≤ nεₙ² eventually, for the paper's εₙ (k ≥ 0) | Easy: elementary asymptotics | **Proved** (`BGP.log_inv_radius_le_n_radius_sq`) |
| UV-2b′ | εₙ > 0, exp(−nεₙ²) ≤ εₙ, and log(1/(c₁εₙ)) ≤ c₂nεₙ², each eventually | Easy: from UV-2b | **Proved** (`BGP.eventually_priorMassRadius_pos`, `BGP.exp_neg_n_radius_sq_le_radius`, `BGP.log_inv_mul_radius_le`) |
| UV-4a | Covering a scaled family and an interval | Moderate: the covering-number API | **Proved** (`BGP.isCover_image2_smul`, `BGP.coveringNumber_image2_smul_le` and `…_le_half`, `BGP.externalCoveringNumber_image2_smul_le`, `BGP.exists_isCover_Icc`, `BGP.coveringNumber_Icc_le`) |
| UV-5 | Tests for ‖f − f′‖ₙ + \|σ − σ′\| on the compact noise range | Moderate: Gaussian tails under `Measure.pi` | **Proved** (`BGP.exists_tests`, by `BGP.tests_mean_case`, `BGP.tests_scale_up_case`, `BGP.tests_scale_down_case`; tails `BGP.pi_gaussian_linear_tail`, `BGP.pi_gaussian_sq_upper_tail`, `BGP.pi_gaussian_sq_lower_tail`) |
| UV-6 | The composed theorem | Heavy; cites GvdV07 and TWCD2025 | **No `sorry`**, conditional on cited axioms. `BGP.posterior_contraction_unknown_variance_gvdv07` rests on GvdV07 as printed (`BGP.gvdv07_theorem1`, `BGP.gvdv07_lemma1`; every Mₙ → ∞). `BGP.posterior_contraction_unknown_variance` gives the fixed C (`BGP.fixed_constant_contraction`; GvdV01 Thm 2.1 form). `…_of_tunable_sieve` is the paper's tunable sieve. `h59` and `hRadius` are confirmed against the paper; UV-3 is proved (`BGP.prior_compl_sieve_le`) |

**Where the pieces live.**
- **`Enthymemata/BgpUnknownVariance.lean`** compiles with scoped imports and no `sorry`. meta-CI's ledger marks it
  PROOF-CLOSED, which makes it a candidate for semantic review and promotion to `Lemmas/`; the names keep their `BGP`
  namespace on promotion. The 2026-10-23 review (`Archeion/Markdown/bpg-review.md`) asked for the UV-2b companions and
  the UV-1d bridge; both have landed, and the module compiles with no warnings under the lakefile's linter options.
- **`Enthymemata/BgpGaussianKL.lean` (UV-1c)** is PROOF-CLOSED, and so is **`Enthymemata/BgpCovering.lean` (UV-4a)**.
  Both are candidates for the same review. `BgpGaussianKL` imports `BgpUnknownVariance` for UV-1a, UV-1b and UV-1d, so
  promoting it means promoting that module too.
- **`Enthymemata/BgpUnknownVarianceTheorem.lean` (UV-6)** holds the composed statements and their assembly: the model,
  [P1], UV-2's prior-mass steps, the sieve and its cover, the rate asymptotics, the tests (UV-5) and the cited
  axioms.
  - Since 2026-10-25 it has no `sorry`, so meta-CI's ledger reads PROOF-CLOSED.
  - Its theorems are conditional on the cited axioms (`#print axioms`), unlike the three bookkeeping modules.
  - Whether a module that rests on cited axioms moves to a verified stage, as `Stability.lean`'s would, is the
    review's call.
- **`Enthymemata/BgpGaussianTails.lean` (UV-5's tails)** is PROOF-CLOSED and axiom-free: completing the square, the
  product moment-generating function, Chernoff, and the three tails.

UV-1a, UV-1b, UV-2a and UV-2b are the boring load-bearing claims the rigor harness exists for. They are short, they
carry the whole bookkeeping, and an error in any of them would quietly break the extension.
