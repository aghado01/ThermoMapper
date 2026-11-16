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
