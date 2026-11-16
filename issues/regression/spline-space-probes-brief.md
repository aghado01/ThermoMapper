# Spline-space probes from Schumaker: the numerical chip's follow-up

Date: 2026-11-16. Status: **posed at the user's request on 2026-11-16** ("pose the follow up lean and numerical work as appropriate based on schumaker with corrections as needed"). Drafted by the coordinating session. It follows the [latent-recovery pilot](latent-recovery-pilot-brief.md) and starts once that brief's done criteria are met. Its sibling is the Lean chip's [spline-space statements brief](../lean/spline-space-statements-brief.md).

## Why

Schumaker's *Spline Functions: Basic Theory* gives exact theorems for three questions the engine answers today by convention or not at all: when a spline design has full column rank, what free knots buy on a curve with sparse events, and whether pieces that oscillate would serve a ringing signal better than polynomial ones. Each probe below measures one of them against the built engine and names the decision it informs. The engine stays as it is; a gap comes back as a primitive proposed in the report.

## Sources

- **The text,** as markdown: `D:\aghado01\graveyard\pdf2md-prototype\schumaker`, files `schumaker_003_converted.md` to `schumaker_015_converted.md` for Chapters 1 to 13. Each theorem is its own third-level heading; most chapters carry much HTML markup, so read stripped. The chapter PDFs are at `D:\aghado01\codex-scientiae\supellex-private\PDFs\Schumaker\`.
- **The errata:** `errata.md` in the same folder lists misprints the user verified against the printed book. Take the corrected forms, and report any new one to the user for that file.
- **Where to read:** § 4.8 (Theorem 4.61) for P1; Chapter 7, with § 7.6 on saturation classes, for P2; § 9.9, § 10.2, § 10.3 and § 10.8 for P3.
- **The engine:** `src/maths/regression/spline/` (`BSpline`, the bases), `spline/bars/` (the free-knot sampler, its design admission and support rule), `spline/baps/TensorProductDesign.cs`.

## The probes

**P1. Design admission against the exact rank condition.** Theorem 4.61: a B-spline collocation matrix is nonsingular exactly when each site lies in the support of its own B-spline. For a least-squares design this becomes an ordered matching of sites to supports, which a single pass over sorted sites decides (the Lean brief's S7).
- Over a seeded sweep of knot configurations and site sets, compare three verdicts on each design: the matching test, the numerical rank of the design, and the engine's own admission as built.
- Then product designs on scattered sites in two dimensions: whether a matching on each axis is necessary, and whether it is sufficient, for the product design's rank. Keep the counterexamples.
- *Decisions:* whether the built rule refuses designs that are exactly full rank, or admits ones that are not; whether the matching test could stand in front of the factorization; what the scattered-site product case needs that the one-dimensional theorem does not give.

**P2. What free knots buy on sparse events.** Theorems 7.10 to 7.12: with k free knots the error is governed by a variational modulus, and order σ holds whenever the (σ−1)th derivative has bounded variation.
- On synthetic curves that are smooth except at a few events (kinks, jumps in a derivative) at known places, measure error against knot count for least squares on uniform knots and for the free-knot sampler as built.
- Compare where the posterior puts its knots with the density the theory calls for. Read § 7.6 for that density before stating it.
- *Decisions:* whether the measured rates separate as the theorems say at the sizes the engine runs; whether posterior knot density tracks the curve's variation closely enough for a knot count to be read as an estimate of how many events a signal holds.

**P3. Pieces that oscillate.** § 10.3: an L-spline's pieces lie in the null space of a differential operator, and its one-sided basis is the operator's Green's function placed at each knot. § 10.8: local trigonometric B-splines exist only while a basis function's support stays short against the carrier's period.
- In probe code, build the one-sided basis for the damped oscillator, g(x) = e^{−γx}·sin(ωx)/ω started at each candidate knot, with its two null-space functions.
- On synthetic sums of decaying resonances with known onsets, compare the knots needed for a given error against a polynomial free-knot fit of the same signal.
- Sweep knot spacing against the carrier's period and record the design's conditioning. The Lean brief's S8 predicts degeneracy when two knots sit half a carrier period apart at order 2; say whether the numbers agree.
- *Decisions:* whether an oscillatory term type saves enough knots to earn its own brief; whether it would need an admission rule on spacing, or should work on a band after mixing down.

## Scope

- **Probes live under `artifacts/spline-space/`,** with small settings, probe sources and result summaries tracked at the matching path under `issues/regression/receipts/`.
- **The engine is read and called as built.** `src/`, `tests/` and `docs/` stay as they are.
- **Ground truth comes from outside the estimator:** the synthetic curve's own events, a rank computed independently of the engine's admission.
- **Gates come first.** Before each probe's first reading, write its gate into the report: sizes, seeds, tolerances and the rerun rule for a marginal failure.

## Done

- P1, P2 and P3 each run against its gate, with its decisions answered or the reason one could not be.
- A report appended to this brief in the pilot brief's form, with the primitives proposed and what each would change.

## Working rules

Those of the [pilot brief](latent-recovery-pilot-brief.md#working-rules).
