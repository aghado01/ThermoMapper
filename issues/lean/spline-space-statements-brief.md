# Spline-space statements from Schumaker: the Lean chip's follow-up

Date: 2026-11-16. Status: **posed at the user's request on 2026-11-16** ("pose the follow up lean and numerical work as appropriate based on schumaker with corrections as needed"). Drafted by the coordinating session. It follows the [time–frequency statements brief](time-frequency-statements-brief.md) and starts once that brief's done criteria are met. Its sibling is the numerical [spline-space probes brief](../regression/spline-space-probes-brief.md).

## Why

Schumaker's *Spline Functions: Basic Theory* holds the classical theorems behind three things the first round argued informally: when a spline design has full column rank, why knots spaced against a carrier's period degenerate, and why free knots pay for sparse events. Posing small, exact instances turns each into something the engine can lean on: an admission check, an admission rule for an oscillatory term type, and a statement of what a knot count measures.

## Sources

- **The text,** as markdown: `D:\aghado01\graveyard\pdf2md-prototype\schumaker`, files `schumaker_003_converted.md` to `schumaker_015_converted.md` for Chapters 1 to 13. Each theorem is its own third-level heading. The chapter PDFs are at `D:\aghado01\codex-scientiae\supellex-private\PDFs\Schumaker\`.
- **The errata:** `errata.md` in the same folder lists misprints the user verified against the printed book. Take the corrected forms. Those known on 2026-11-16: Lemma 10.16's last sentence sums to k where d is meant; the remark after Theorem 10.42 drops the adjoint's star; the proof of Theorem 7.11 swaps its bounds on t; the proof of (7.33) drops a factor 2ᵐ on ε. Report any new one to the user for that file.
- **Where to read:** § 4.8 (Theorems 4.61 and 4.63), § 7.4 (Theorems 7.10 to 7.12), § 10.3 (Definition 10.14 to Theorem 10.17), § 10.8 (Theorem 10.42 to 10.46).

## The statements

Numbering continues the first brief's. Each is the coordinating session's posing from a reading of the sections above; sharpen as Lean requires and say what changed.

**S7. Ordered supports and full column rank** (after Theorem 4.61).
- *S7a, zero block.* An n×n matrix over a field with r rows and s columns on which every entry is zero, r + s > n, has determinant zero.
- *S7b, necessity.* Sites t₁ < … < tₙ; functions φ₁, …, φₙ with φⱼ(t) = 0 unless ℓⱼ < t < uⱼ, the ends ℓ and u nondecreasing in j. If det(φⱼ(tᵢ)) ≠ 0, then ℓᵢ < tᵢ < uᵢ for every i. The book's supports close on the left at a full-multiplicity knot (4.133); take the open case first.
- *S7c, designs.* For N ≥ n sorted sites, if the N×n matrix (φⱼ(tᵢ)) has full column rank, some strictly increasing selection i₁ < … < iₙ has ℓₖ < t_{iₖ} < uₖ for every k.
- *S7d, the greedy check.* Taking, for k = 1, …, n in turn, the first unused site after the previous pick that lies in the k-th support succeeds exactly when such a selection exists.
- *S7e, sufficiency, cited.* For B-splines the selection makes the collocation matrix nonsingular, so the design has full column rank. State it with Theorem 4.61 as a named premise, the way `Enthymemata/BgpUnknownVarianceTheorem.lean` carries its cited results. Proving 4.61 itself needs the book's zero counts and is outside this brief.
- *Why it matters.* Together they say a spline design's rank can be decided from the knots and the sorted sites alone, before any factorization: data inside every column's support is not enough, since one site in an overlap cannot serve two columns.

**S8. The order-2 trigonometric B-spline and its spacing condition** (after § 10.3 and § 10.8).
- For y₀ < y₁ < y₂ with y₂ − y₀ < 2π, let Q(x) = sin((x − y₀)/2) / (sin((y₁ − y₀)/2)·sin((y₂ − y₀)/2)) on [y₀, y₁), Q(x) = sin((y₂ − x)/2) / (sin((y₂ − y₁)/2)·sin((y₂ − y₀)/2)) on [y₁, y₂), and 0 elsewhere. This is (10.103) applied to (10.99).
- *To show:* each piece lies in span{cos(x/2), sin(x/2)}, the null space of D² + ¼; Q is continuous on ℝ; Q > 0 on (y₀, y₂).
- *Sharpness.* When y₂ − y₀ = 2π, the single function sin((x − y₀)/2) on [y₀, y₂], zero outside, is continuous, has its pieces in that span and needs no interior knot. So a local-support element exists on two knots, which the corrected Lemma 10.16 excludes when its Tchebycheff hypothesis holds; the hypothesis is necessary.
- *Scaled form.* With carrier ω the pieces lie in the null space of D² + ω² and the condition reads y₂ − y₀ < π/ω, half the carrier's period.
- *The damped oscillator.* For γ ≥ 0 and ω > 0, g(x) = e^{−γx}·sin(ωx)/ω for x ≥ 0 and 0 for x < 0 is continuous, has g(0) = 0 and right derivative 1 at 0, and satisfies g″ + 2γg′ + (γ² + ω²)g = 0 for x > 0. As the instance of Theorem 10.15 at order 2 with simple knots x₁ < … < x_k: the two null-space functions together with the g(· − xᵢ) are linearly independent, and they span the continuous functions on [a, b] whose pieces lie in the null space.
- *Why it matters.* It is the admission rule an oscillatory term type would need, and the exact sense in which a decaying resonance started at a knot is a spline.

**S9. Free knots against fixed knots, for piecewise constants** (the σ = 1 case of Theorem 7.12).
- *Free.* For f monotone on [a, b] and k ≥ 0, some step function with at most k breakpoints in (a, b) has sup |f − s| ≤ (f(b) − f(a)) / (2(k + 1)).
- *Fixed.* Some monotone f with f(b) − f(a) = 1 has sup |f − s| ≥ ½ for every k and every step function s on the k uniform breakpoints a + j(b − a)/(k + 1).
- *Optional.* The free-knot bound for functions of bounded variation, with the total variation in place of f(b) − f(a) and a constant the proof supplies.
- *Why it matters.* With free knots the rate depends on how much variation there is, not on where it sits. It is the simplest exact form of the claim that knots follow innovations.

## Scope and order

- **Order:** S9 and S7a to S7d first, since they need no spline theory. Then S8. S7e last.
- **A witness for every statement,** and `#print axioms` for every closed result, as in the first brief.
- **Placement** beside the regression-family modules, named by what each module states.
- **Left to the numerical chip:** how the engine's built admission compares with S7's check, and every rate measured on data.

## Done

- S7a to S7d and S9 proved without `sorry`, each with its witness.
- S8 with every part that closes proved and the rest compiling as statements; S7e compiling with its named premise.
- A report appended to this brief in the first brief's form.

## Working rules

Those of the [first brief](time-frequency-statements-brief.md#working-rules).

## Report

### 2026-11-16 — S7 to S9: all proved, none false as posed

**Landed:** `8a3cd9d` (S7, S9) and the commit that adds this report (S8).

**Per statement.** Modules sit in `lean/Enthymemata/Spline/`, namespace `Spline`.

| | Module | Status | Changes from the posed form, and why |
|---|---|---|---|
| S7a | `DesignRank.lean` | proved | **Field** relaxed to any commutative ring (`det_eq_zero_of_zero_block`). |
| S7b | `DesignRank.lean` | proved | **Sites** need only be nondecreasing. **Supports** may close on the left at ℓ j when `closed j` holds, which covers (4.133)'s full-multiplicity case (`inSupport_of_det_ne_zero`). That needs one structural hypothesis: a support closing at ℓ μ forces every earlier support with the same left end to close. B-splines satisfy it, since equal y_i = y_μ with μ > i makes y_i = … = y_{i+m−1}. **The open case as posed** is `mem_Ioo_of_det_ne_zero`. |
| S7c | `DesignRank.lean` | proved | **Full column rank** is stated as injectivity of c ↦ Mc. **N ≥ n** follows and is not assumed. The route is `exists_strictMono_det_ne_zero`: n independent rows, sorted, give a nonsingular square submatrix. |
| S7d | `DesignRank.lean` | proved | **The greedy** is correct for any decidable membership predicate, not only for interval supports (`exists_isSelection_iff_greedy`). It is stated on ℕ indices. `isSelection_iff_strictMono` translates to Fin-indexed strictly increasing selections. |
| S7e | `DesignRank.lean` | proved, with premise | **The premise** `h461` is Theorem 4.61's sufficiency half, carried as a hypothesis (`fullRank_of_selection`). `fullRank_iff_greedy` closes the loop: full column rank ⟺ greedy succeeds. **Satisfiable:** `h461_indicator` proves the premise at order 1. **Sites** must be strictly increasing, as 4.61 requires. |
| S8 | `TrigonometricBSpline.lean` | proved | **The carrier** is a general ω, with `trigBSpline_half` the posed ω = ½. **The definition** is written with continuous branches; `trigBSpline_eq_left`, `_right` and `_zero` give the posed half-open pieces. **The null space** claim is split: pieces lie in span{cos ωx, sin ωx}, and that span solves f″ + ω²f = 0. That nothing else solves it (ODE uniqueness) is not formalized. **Sharpness** adds `collocation_det`: the span's Tchebycheff determinant at s, t is sin(ω(t − s)), zero at t − s = π/ω. **The Green's function** is built on `nullDamped`, the null space of D² + 2γD + (γ² + ω²). **Theorem 10.15 at order 2** is proved in full (`spline_span`, `spline_indep`); its pieces are taken on closed intervals, so continuity is implied, not assumed. |
| S9 | `FreeKnotSteps.lean` | proved | **"k breakpoints"** is stated as at most k + 1 interval pieces (`PiecewiseConst`), each interval open or closed at either end. A monotone f's jumps need that: half-open pieces cannot follow a jump that belongs to the left. **Free:** the constant is exact, (f b − f a)/(2(k + 1)) (`exists_piecewiseConst_monotone`). **Fixed:** the step function is constrained only on the open pieces, which makes the lower bound strongest (`exists_monotone_far_from_uniform`). One f for every k needs a jump at an irrational relative position, a + (b − a)(√2 − 1); a rational one is a breakpoint for some k. **The optional BV form** is not done. |

**Declarations to print axioms for:**
- *S7a, S7b:* `det_eq_zero_of_zero_block`, `inSupport_of_det_ne_zero`, `mem_Ioo_of_det_ne_zero`.
- *S7c:* `exists_strictMono_det_ne_zero`, `exists_selection_of_fullRank`.
- *S7d:* `exists_isSelection_iff_greedy`, `isSelection_iff_strictMono`.
- *S7e:* `fullRank_of_selection`, `fullRank_iff_greedy`, `h461_indicator`.
- *S8:* `trig_ode`, `continuous_trigBSpline`, `trigBSpline_pos`, `trigBSpline_half`, `trigBSpline_eq_left`, `trigBSpline_eq_right`, `trigBSpline_eq_zero`, `collocation_det`, `twoKnotBump_spec`, `nullDamped_ode`, `green_spec`, `green_ode`, `spline_span`, `spline_indep`.
- *S9:* `exists_piecewiseConst_monotone`, `exists_monotone_far_from_uniform`.
- Each prints only `propext`, `Classical.choice` and `Quot.sound`. The two greedy results print only `propext` and `Quot.sound`.

**Witnesses:**
- *S7:* order-1 splines on knots 0, 1, 2 at sites ½, 3/2. The chain runs through `h461_indicator`, so the cited premise is met, not assumed.
- *S8:* the posed spline on knots 0, 1, 2 (2 < 2π); the two-knot bump at ω = ½ on [0, 2π]; the Green's function at γ = 1, ω = 2; independence on knots 0, 1, 2.
- *S9:* f = id on [0, 1] with k = 1, and the fixed-knot f on [0, 1].

**Gates:**
- `lean` on each file during drafting.
- `lake build Enthymemata.Spline` with the linter set: no warnings after fixes.
- `scripts/meta-ci.ps1`: PASS with its build gate, all three modules `PROOF-CLOSED`.
- Cache: one scoped fetch, `lake exe cache get Mathlib.NumberTheory.Real.Irrational` (3 files). D: had 22 GB free at the end.

**Errata for the user.** No new misprint in the book was found in §§ 4.8, 7.4, 10.3 or 10.8. The In-progress entry of `8a3cd9d` said the conversion drops Theorem 7.10's display. That was wrong: (7.30) is in `schumaker_009_converted.md`, and I had read past it.

**Follow-ups:**
- *Promotion review* of the three modules, all `PROOF-CLOSED`.
- *Theorem 4.61 itself,* which would discharge `h461`. It needs the book's zero counts for splines (Theorem 4.56), outside this brief.
- *S9's optional BV form.* The quantize-the-variation construction gives TV/(k + 1) directly, and TV/(2(k + 1)) with midranges on each piece.
- *S8's null space as the full solution set* (ODE uniqueness), if a statement needs "exactly".

#### Methods

**S7a, the zero block.**
- Write det M = Σ_σ sgn σ·Π_i M_{σ(i), i}.
- If some permutation avoided the block, σ⁻¹ would map the r block rows injectively into the n − s columns outside it, so r ≤ n − s.
- Since r + s > n, every product has a factor in the block.

**S7b, necessity.**
- *Too far left.* Suppose t_i ≤ ℓ_i, and t_i is not the closed end. For ν ≤ i ≤ μ, t_ν ≤ t_i ≤ ℓ_i ≤ ℓ_μ, so φ_μ(t_ν) = 0. The only escape is t_ν = ℓ_μ with support μ closed. The chain then collapses to ℓ_i = ℓ_μ = t_i, and the closure hypothesis makes support i closed, a contradiction. That is a zero block of (i + 1) rows by (n − i) columns.
- *Too far right.* Suppose t_i ≥ u_i. Then ν ≥ i ≥ μ gives t_ν ≥ t_i ≥ u_i ≥ u_μ, a zero block of (n − i) rows by (i + 1) columns.
- *Either way* r + s = n + 1, and S7a applies.

**S7c, designs.**
- *The rows span.* Injectivity of c ↦ Mc means no nonzero functional ⟨·, c⟩ kills every row. A proper span would admit one (`Submodule.exists_le_ker_of_lt_top`), so the rows span ℝⁿ.
- *Square submatrix.* A linearly independent subfamily that spans is a basis, so it has n members. Sorting their indices gives an order embedding e, and the rows of M at e are linearly independent, so that submatrix is nonsingular.
- *Selection.* S7b on the sites t ∘ e, still sorted, gives the selection.

**S7d, the greedy.**
- *Dominance,* by induction on k: if e is any selection, the greedy's k-th pick p_k exists and p_k ≤ e_k.
  - Column 0: e_0 witnesses a site in the first support, and the greedy takes the first one.
  - Column k + 1: e_{k+1} > e_k ≥ p_k witnesses a site after p_k.
- *Converse.* When every pick exists, the picks form a selection by construction.
- Nothing uses intervals, so the check is exact for any membership pattern.

**S7e, sufficiency with 4.61 as premise.**
- On a selection e, the premise makes the square submatrix S of rows e nonsingular.
- Mc = 0 restricts to Sc = 0, so c = 0.
- With S7c, S7d and `isSelection_iff_strictMono`: rank full ⟺ selection exists ⟺ greedy succeeds.
- *The order-1 premise.* For indicators of (y_j, y_{j+1}) at sites in their own supports, a site in support k lies in no other. So the collocation matrix is the identity.

**S8, the trigonometric B-spline.**
- *Span.* sin(ω(x − c)) = cos(ωc)·sin ωx − sin(ωc)·cos ωx. Differentiation maps (A, B) to (ωB, −ωA), and twice to −ω²(A, B), which is the equation.
- *Continuity* is three `Continuous.if_le` joins.
  - At y₀ and y₂ the pieces vanish, since sin 0 = 0.
  - At y₁ both equal 1/sin(ω(y₂ − y₀)). That needs only sin(ω(y₁ − y₀)) ≠ 0 and sin(ω(y₂ − y₁)) ≠ 0.
- *Positivity.* Each numerator's argument lies in (0, ω(y₂ − y₀)) ⊂ (0, π), and so does each denominator's.
- *Sharpness.*
  - The order-2 Tchebycheff determinant is det[[cos ωs, sin ωs], [cos ωt, sin ωt]] = sin(ω(t − s)). Haar holds on intervals shorter than π/ω and fails at π/ω.
  - There sin(ω(x − y₀)), cut to [y₀, y₂], vanishes at both ends. It is the continuous two-knot element that Lemma 10.16, with its corrected Σ_{i=1}^{d}, rules out under the ET hypothesis.

**S8, the damped oscillator.**
- *The null space.* Differentiation maps e^{−γx}(A cos ωx + B sin ωx) to the same form with (A, B) ↦ (−γA + ωB, −ωA − γB). Applied twice, A₂ + 2γA₁ + (γ² + ω²)A = 0 and likewise for B. The whole null space solves the equation, and the Green's function is its member (0, 1/ω) cut at 0.
- *Its value at 0.* The right derivative at 0 is the first coefficient after one step, −γ·0 + ω·(1/ω) = 1. For x > 0 the function agrees with the smooth member near x, so its derivatives do too (`EventuallyEq.deriv`).

**S8, Theorem 10.15 at order 2.**
- *Shift.* On [τ, ∞), green(· − τ) is in the null space, by angle subtraction with the factor e^{γτ}.
- *Vanishing.* A member vanishing at τ is a multiple of green(· − τ) there. Write ωt = ω(t − τ) + ωτ; the cos(ω(t − τ)) coefficient is the value at τ, which is zero.
- *Span,* by induction over pieces. On [τ_{k+1}, ∞) the current representation is one member of the null space. The next piece agrees with it at τ_{k+1}, so their difference vanishes there and is c·green(· − τ_{k+1}). That term is zero left of τ_{k+1}, so the earlier pieces are untouched.
- *Independence.*
  - On the first piece only the null-space part survives. At two points δ apart, with 0 < ωδ < π, the collocation determinant sin(ωδ) ≠ 0 forces A = B = 0.
  - Then left to right: on [τ_i, τ_{i+1}], earlier coefficients are zero and later Green's terms vanish. Evaluating at τ_i + δ, where green(δ) > 0, forces d_i = 0.

**S9, free against fixed.**
- *Free.* Cut [f(a), f(b)] into k + 1 bands of height h = (f(b) − f(a))/(k + 1). The band index q(x) = min(⌊(f(x) − f(a))/h⌋, k) is monotone in x, so each band's preimage is an interval, possibly open or closed at either end. Setting s = f(a) + (q + ½)h gives \|f − s\| ≤ h/2, because qh ≤ f(x) − f(a) ≤ (q + 1)h. The cap at k covers f(x) = f(b).
- *Fixed.* With r = √2 − 1, the position ρ = (k + 1)r is irrational, so it lies strictly between ⌊ρ⌋ and ⌊ρ⌋ + 1. The jump at a + (b − a)r is then interior to an open uniform piece. A constant v on that piece is at distance max(\|v\|, \|1 − v\|) ≥ ½ from one side of the jump.
