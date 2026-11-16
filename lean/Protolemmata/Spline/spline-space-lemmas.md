# Spline-space lemmas after Schumaker

Candidate statements for [the spline-space brief](../../../issues/lean/spline-space-statements-brief.md), as sharpened for Lean. The coordinating session posed S7–S9 on 2026-11-16 from Schumaker §§ 4.8, 7.4, 10.3 and 10.8, with the errata in `D:\aghado01\graveyard\pdf2md-prototype\schumaker\errata.md`. The brief's report carries each change from the posed form and the proofs' methods.

## S7. Ordered supports and full column rank (after Theorem 4.61)

**Status:** proved, `Enthymemata/Spline/DesignRank.lean`. S7e is conditional on Theorem 4.61's sufficiency half, carried as the premise `h461` and proved at order 1 (`h461_indicator`).

- **S7a.** A zero block of r rows and s columns with r + s > n gives det = 0, over any commutative ring.
- **S7b.** A nonsingular collocation matrix puts site i in support i.
  - Supports may close on the left at a full-multiplicity knot, given that closing propagates to earlier supports with the same left end.
  - Sites need only be nondecreasing.
- **S7c.** Full column rank (c ↦ Mc injective) gives an increasing selection of sites, one in each support.
- **S7d.** The first-fit greedy succeeds iff a selection exists, for any decidable membership predicate.
- **S7e.** Given h461 and strictly increasing sites: full column rank ⟺ a selection exists ⟺ the greedy succeeds.

## S8. The order-2 trigonometric B-spline (after § 10.3 and § 10.8)

**Status:** proved, `Enthymemata/Spline/TrigonometricBSpline.lean`.

- **The B-spline** with carrier ω: (10.103) on (10.99). It is continuous on ℝ and positive on (y₀, y₂) when ω(y₂ − y₀) < π; at ω = ½ the condition reads y₂ − y₀ < 2π.
- **Pieces** lie in span{cos ωx, sin ωx}, which solves f″ + ω²f = 0.
- **Sharpness at ω(y₂ − y₀) = π.** The two-knot bump exists, and the Tchebycheff determinant sin(ω(t − s)) vanishes.
- **The damped oscillator.**
  - The null space of D² + 2γD + (γ² + ω²) is e^{−γx}(A cos ωx + B sin ωx).
  - The Green's function is continuous, zero at 0, has right derivative 1, and solves the equation for x > 0.
  - Theorem 10.15 at order 2, simple knots: span and independence.

## S9. Free against fixed knots for piecewise constants (σ = 1 of Theorem 7.12)

**Status:** proved, `Enthymemata/Spline/FreeKnotSteps.lean`.

- **Free.** A monotone f is within (f(b) − f(a))/(2(k + 1)) of a piecewise constant with k + 1 interval pieces.
- **Fixed.** One monotone f with unit rise, its jump at a + (b − a)(√2 − 1), stays ½ away from every step function on the k uniform breakpoints, for every k.
- **Open:** the optional bounded-variation form.
