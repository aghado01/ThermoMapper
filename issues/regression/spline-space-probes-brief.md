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

## Report (the numerical chip, from 2026-11-16)

### In progress

- 2026-11-16: started after the pilot closed (`9128238`). The probe is at `artifacts/spline-space/probe/`, built against the same `be979dd` snapshot as the pilot (`artifacts/latent-recovery/source/`). The gates for P1, P2 and P3 are below, written before any probe code ran. Runs next, P1 first.

### Gates as pre-stated

Each gate was written here before the first reading of its probe. A gate changed after a failure is recorded below it as a change, with the reason. Schumaker is read in the user's conversion, with the corrected forms from `errata.md`.

**P1, design admission against the exact rank condition.**

- *Bases.* Clamped B-splines of degree d ∈ {1, 2, 3} (order m = d + 1) on [0, 1], interior knots ξ₁ < … < ξ_k, k ∈ {0, …, 12}, n = k + m columns. The probe evaluates them with its own Cox–de Boor recursion and checks it against the engine's `SplineBasis.Design` (max |difference| recorded; above 10⁻¹² the comparison stops). Each support σ_j follows (4.133) and Theorem 4.63: open (y_j, y_{j+m}) when y_j < y_{j+m−1}, closed on the left otherwise, and closed at x = 1 for the last B-spline.
- *Designs.* Knot laws: uniform order statistics, or clustered (half the knots inside one window of width 0.02 placed uniformly). Site counts N ∈ {n − 1, n, n + 1, n + 3, 2n}. Site laws: uniform on [0, 1]; one to three clusters of width 0.05; uniform with a random gap of width 0.2 removed. 20 draws per (degree, k, N, knot law, site law): 23 400 designs. Seed `SeedTree.Derive(2026111701, 4)`: [designs, products, reruns, spare].
- *Three verdicts per design.* M: the greedy ordered matching (the Lean chip's S7d) on the sorted sites and the supports. R: the singular values of the probe's design by one-sided Jacobi, with ρ = σ_min/σ_max recorded. E: the engine's `SplineTarget.Evaluate` on the design (cubic, quadratic or linear `SplineBasis`, `WeightedNormalModel`, `UniformPrior(12)`, the default support gate), read as its `DesignRejection`.
- *Pass, the theorem's two halves on these designs:*
  - (i) every design with M false has ρ ≤ 10⁻¹³, singular to rounding;
  - (ii) the designs with M true are read, not gated: exact full rank cannot be read numerically, so the reading is ρ's distribution and the count below 10⁻¹².
  - *Rerun rule for (i):* a violation with ρ ≤ 10⁻¹¹ is rechecked on the design with every entry outside its support set to exactly zero. If the recheck clears it, it is reported as rounding; anything above 10⁻¹¹, or not cleared, is a counterexample to report.
- *Engine readings:* the cross-table of M against E's reason, which answers the decisions: designs refused although M is true (by reason), and designs admitted although M is false.
- *Product designs in two dimensions:* cubic on both axes, k_x, k_y ∈ {0, …, 4}, n = n_x·n_y columns; scattered sites, N ∈ {n, 2n, 4n}; site laws uniform in the square, on the diagonal y = x, near it (|y − x| ≤ 0.02), and in clusters; 20 draws per cell (3000 designs). Verdicts:
  - per-axis matching on the distinct projected sites (M_x, M_y);
  - a bipartite matching of sites to product supports (augmenting paths; a site serves a column when it lies in the product of the two supports);
  - the numerical rank.

  Read four counts: full rank with M_x or M_y false (necessity of per-axis matching), M_x and M_y true but rank-deficient (sufficiency), full rank without a bipartite matching (its necessity), and a bipartite matching but rank-deficient (its sufficiency). Rank-deficient here means ρ ≤ 10⁻¹³. The counterexamples are kept.

**P2, what free knots buy on sparse events.**

- *Theory, read in § 7.4 to § 7.6 before stating it.* With free knots, order m holds on BV^{m−1} (Theorem 7.12) and on functions smooth except at finitely many points (Theorem 7.20, knots placed at the singularities). A balanced partition equalizes the local error, so on the smooth parts the knot density is ∝ |D^m f|^σ with σ = 1/(m + 1/p): σ = 2/9 for cubic least squares (m = 4, p = 2). A kink, a jump in f′, is a spline knot of multiplicity m − 1 = 3: the theory puts three knots at it. With uniform knots a kink costs the order: the L₂ error falls as h^{3/2}.
- *Curves* on [0, 1]: background b(x) = sin(2πx) + 0.5·cos(5πx), plus E kinks Σ a_e·|x − x_e|, E ∈ {0, 1, 2, 4, 8}, a_e = ±U[0.5, 1], x_e uniform on [0.1, 0.9] with gaps ≥ 0.05; one more family with E = 2 jumps in f″ (a_e·(x − x_e)₊²).
- *Data:* n = 2000 equispaced points, Gaussian noise σ = 10⁻³; 5 replicates per curve family, fresh events and noise each.
- *Fits:*
  - Least squares on k ∈ {2, 4, 8, 16, 32, 64} uniform interior knots, through the engine's `SplineBasis` and the probe's solve.
  - BARS as built (cubic, `WeightedNormalModel`, its default kernel), once with a free count (`UniformPrior(80)`) and once per k ∈ {2, 4, 8, 16, 32} with the count held near k by a steep probe-side prior (log p = −50·|j − k|).
  - Each fit uses 2 chains, burn-in 2000, 1000 kept.
- *Readings:*
  - RMS error against the clean curve on 4001 grid points: for least squares per k; for BARS, the posterior-mean fit per run, at its posterior mean knot count.
  - The slope of ln error on ln k over k ∈ [4, 32], by curve family and method.
  - The pooled posterior knots: the share within ±0.02 of an event against that window's share of [0, 1]; knots per event; and, away from events, the histogram against the density ∝ |b⁽⁴⁾|^{2/9}, as a total-variation distance.
  - Under the free count: the posterior mean k against E, its fitted slope and residual sd.
- *Pass, rate separation:* on the kink families (E ≥ 1), the free-knot slope is steeper than the uniform slope by at least 1.0 (the theory's gap is 2.5: −4 against −1.5). *Marginal:* a gap in [0.8, 1.0): rerun the family on the rerun seeds, pass only if the rerun reaches 1.0. The density, knots-per-event and k-against-E rows are readings.
- *Seeds:* `SeedTree.Derive(2026111702, 8)`: curves, noise, least squares (unused), BARS free, BARS held, reruns, spare, spare.

**P3, pieces that oscillate.**

- *Basis, in probe code:* the damped oscillator's one-sided basis on [0, 1]. Two null-space functions, e^{−γx}·cos(ωx) and e^{−γx}·sin(ωx), and g(x − ξᵢ) at each interior knot, g(x) = e^{−γx}·sin(ωx)/ω for x ≥ 0 and 0 before (S8). It implements the engine's `IBasis`, so BARS's free-knot sampler runs on it as built. The support gate is turned off (`MinColumnMass` 0), because its column-mass rule assumes nonnegative columns and these are signed. The coefficient prior is centred at zero for both bases (`GPriorCentre.Zero`).
- *Signals:* n = 1000 equispaced points; ω = 2π·20 (20 cycles over the domain); γ = 8; R ∈ {1, 2, 4} resonances A_r·g(x − t_r)·ω, onsets uniform on [0.05, 0.8] with gaps ≥ 0.05, A_r = ±U[0.5, 1.5]; noise σ = 0.01; 5 replicates per R. Mismatch variants at R = 2: the data's ω × 1.05, and γ × 1.5, the basis keeping the nominal values.
- *Fits:*
  - BARS cubic (`UniformPrior(150)`, the default support gate).
  - BARS on the oscillator basis (`UniformPrior(20)`, gate off).
  - Least squares on the oscillator basis with knots at the true onsets.

  Each BARS run uses 2 chains, burn-in 2000, 1000 kept. Readings: RMS error of the posterior mean against the clean signal, the posterior mean knot count, and wall time.
- *Decision threshold, stated in advance:* an oscillatory term type earns its own brief when, on the matched signals, BARS on it reaches an error at most 1.5 times the cubic fit's with at most a quarter of its knots, and under each mismatch stays within twice the cubic fit's error.
- *Spacing sweep:* knots uniformly spaced at h with hω/π ∈ {0.05, 0.10, …, 1.20}, plus 0.48–0.52 in steps of 0.005, on 4000 sites, γ ∈ {0, 8}. Three readings per spacing:
  - the one-sided basis's equilibrated Gram condition number (Jacobi eigenvalues);
  - the order-2 local element Q of S8, scaled by e^{−γ(x−y₀)}: its sup norm and whether it takes negative values;
  - the equilibrated condition number of the local elements written unnormalized (S8's numerators).

  *S8 agrees* if Q's sup norm diverges at y₂ − y₀ = π/ω (at hω/π = 0.5, at least 10⁶ times its value at 0.25) and Q takes negative values for y₂ − y₀ beyond π/ω. The conditioning rows are readings.
- *Seeds:* `SeedTree.Derive(2026111703, 4)`: signals, noise, BARS, spare.
