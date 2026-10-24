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

Notes:
- The bound is elementary real analysis. It is the best Lean target here.
- The closed forms need Gaussian moments (Mathlib's `gaussianReal`); worth formalizing once that API is in reach at
  our pinned version.

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
- **Lean-worthy pieces:** the scaling identity (step 2) and the absorption of log(1/εₙ) (step 4). **Both discharged
  2026-10-23** in `Enthymemata/BgpUnknownVariance.lean`:
  - `BGP.map_smul_closedBall`;
  - `BGP.log_inv_radius_le_n_radius_sq`, for k ≥ 0, which covers the paper's k₁ = (1 + D)/(2 + ϱ/s) > 0. Mathlib's
    `Pow.Asymptotics` module is not in the local cache, so the proof is elementary instead:
    - log n ≤ nᵇ/b, from log y ≤ y − 1;
    - an explicit threshold N = max(⌈e⌉, ⌈M^(1/b)⌉), with b = (1 − 2a)/2.

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
- **The heavier-tailed variant.** With inverse-gamma tails on σ² instead of a compact range, cover σ on the log
  scale: the location-scale Hellinger distance is controlled by |log σ − log σ′|. An exponentially wide range then
  costs only log(log σ⁺ / ε̄ₙ) = O(log n). This is the textbook route for an unknown variance.
- **Lean-worthy:** the covering bound for a scaled family and for an interval (metric-space covering numbers;
  check the Mathlib covering-number API at our pinned version).

```lean
import Mathlib.Topology.MetricSpace.Bounded

/-- UV-4a. A cover of `B` at radius ε/(2b) and of `[a, b]` at radius ε/(2R) (R bounding `B`) cover the scaled
    family `{s • g | s ∈ [a, b], g ∈ B}` at radius ε. Strategy: s•g − s'•g' = s•(g − g') + (s − s')•g'. -/
theorem scaled_family_cover_card {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] :
    -- coveringNumber ε (scaled family) ≤ coveringNumber (ε/(2b)) B * coveringNumber (ε/(2R)) (Icc a b)
    sorry := by
  sorry
```

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
| UV-1c | Gaussian KL and log-ratio variance closed forms | Moderate: Gaussian moments in Mathlib | Open; Lemma once the API is in reach |
| UV-2a | Ball mass under an amplitude rescaling | Easy: preimage of a ball under scaling | **Proved** (`BGP.map_smul_closedBall`) |
| UV-2b | log(1/εₙ) ≤ nεₙ² eventually, for the paper's εₙ (k ≥ 0) | Easy: elementary asymptotics | **Proved** (`BGP.log_inv_radius_le_n_radius_sq`) |
| UV-4a | Covering a scaled family and an interval | Moderate: the covering-number API | Open |
| UV-6 | The composed theorem | Heavy; cites GvdV07 and TWCD2025 | Open; an enthymema with `sorry` at UV-5 |

**Where the proved pieces live.** `Enthymemata/BgpUnknownVariance.lean` compiles with scoped imports and no
`sorry`. meta-CI's ledger marks it PROOF-CLOSED, which makes it a candidate for semantic review and promotion to
`Lemmas/`; the names keep their `BGP` namespace on promotion.

UV-1a, UV-1b, UV-2a and UV-2b are the boring load-bearing claims the rigor harness exists for. They are short, they
carry the whole bookkeeping, and an error in any of them would quietly break the extension.
