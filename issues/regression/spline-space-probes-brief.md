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

- Closed on 2026-11-16. P1, P2 and P3 each ran against its gate and are read below; nothing is left running. The probe is at `artifacts/spline-space/probe/`, built against the pilot's `be979dd` snapshot (`artifacts/latent-recovery/source/`); its sources and summaries are tracked in `receipts/spline-space/`. Post-gate readings (`p3-uniform`, `p2-diagnostic`) are labelled where they are read and changed no gate. Verdicts: P1 passes; P2 fails; P3's sweep agrees with S8, and its knot threshold is not met as written.

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

### Readings

**P1: the theorem's necessity half holds on every design. The engine admits nothing singular, but its support gate refuses most exactly full-rank designs, nearly all of them well conditioned.** Rows and summary in `receipts/spline-space/results/p1-summary.md`.

- *The probe's basis matches the engine's:* the largest |difference| between the probe's Cox–de Boor design and `SplineBasis.Design` over all 23 400 designs is 5.6×10⁻¹⁶.
- *(i) passes.* The matching fails on 19 083 designs, and every one is singular to rounding: max ρ = 2.2×10⁻¹⁷, so the recheck was never needed. The greedy matching reads full column rank exactly in this sweep.
- *(ii), read:* 4317 designs match. Their ρ has median 0.034 and 1st percentile 5×10⁻⁸. Six fall below 10⁻¹⁰ and two below 10⁻¹², all cubic with a site barely inside a support: exactly full rank, numerically singular. A matching cannot stand in for the condition check.
- *The engine against the matching:*

  | matching | engine verdict | designs | with N > n |
  |---|---|---|---|
  | false | ThinSupport | 19 083 | 10 531 |
  | true | admitted | 787 | 787 |
  | true | ThinSupport | 3530 | 2722 |

  The engine admitted no unmatched design; its support gate caught every one before a factorization. Of the 3509 matched designs with N > n, it refused 2722 (78%), all at the support gate (column mass Σᵢ Bⱼ(xᵢ) ≥ 1). 2669 of those have ρ ≥ 10⁻⁶ and 49 have ρ between 10⁻¹⁰ and 10⁻⁶. The refusals spread over every degree and both knot laws, and over every site law (uniform 1167, clustered 418, gapped 1137).
- *Decisions:*
  - The built rule refuses designs that are exactly full rank, and well conditioned ones at that: a column needs a full point's worth of mass, where the rank needs only one site per support, in order. It admits none that is not full rank.
  - The matching test can stand in front of the factorization as an exact, O(N + n) rank test, but not in place of the condition check behind it.
  - Whether the mass gate should give way to the matching is a modelling question (the gate is a declared prior restriction, not a rank test), raised below as a primitive and not decided here.
- *Products, 6000 scattered designs, cubic on both axes:*
  - Per-axis matching is necessary (no full-rank design fails it) and far from sufficient (3040 of the 4211 designs that match on both axes are rank-deficient).
  - A bipartite matching of sites to product supports is also necessary (no full-rank design lacks one).
  - It is sufficient on every uniform design (681 of 681), but not on sites along a line: 374 of 374 exact-diagonal designs with a matching are rank-deficient, as are 3 near the diagonal and 1 clustered.

  The kept counterexamples are the diagonal ones: on the curve y = x the product space restricts to univariate splines of higher degree with merged knots, so its rank there falls below n_x·n_y whatever the matching says. The scattered case needs what the one-dimensional theorem cannot give: a genericity condition on the sites' joint position, beyond any matching.

**P3, the spacing sweep: S8 agrees as posed, but its degeneracy is the local element's, not the design's.** Rows in `receipts/spline-space/results/p3-sweep.csv`.

- *S8's element:* Q's sup norm is 1.0 at hω/π = 0.25 and 8.2×10¹⁵ at 0.5, where y₂ − y₀ = π/ω (7.4×10¹⁵ with γ = 8). Q takes negative values at every spacing beyond, from 0.505 to 1.2, with or without damping. That is the gate's agreement: the normalized, positive local element exists exactly while y₂ − y₀ < π/ω, half the carrier's period. A second divergence sits at hω/π = 1 (sup 3×10³¹), where one interval spans half a period.
- *The design's conditioning has no feature at 0.5:*
  - The one-sided basis's equilibrated condition number falls smoothly as the knots spread: 7.6×10⁵ at hω/π = 0.05, 3.1×10³ at 0.5 and 1.3×10³ at 0.7 undamped; with γ = 8, 4.9×10⁵, 2.1×10² and 51. Its one bump is at hω/π = 1 (4.2×10³ undamped, 91 damped), where a knot's g and its neighbour's cancel beyond the second knot.
  - S8's elements written unnormalized (its numerators) stay at κ ≈ 3–4.5 through 0.5 and rise only toward 1 (680 at 1.0).
- *Decision, for the spacing question:* an oscillatory term built on the one-sided basis needs no admission rule on spacing, only the usual guard against crowded knots (its conditioning grows as h → 0). A term built on local, positive elements needs the rule y₂ − y₀ < π/ω, and it must avoid h = π/ω. The knot-count question waits on the BARS runs below.

**P3, the fits: the oscillator basis fits a sum of decaying resonances with about two knots per resonance at the noise level; cubic BARS as built stays near zero knots on 13 of 15 matched signals. The pre-stated threshold is not met as written, because its knot clause compares against a cubic fit that collapsed.** Rows in `receipts/spline-space/results/p3-summary.md` and `p3-uniform-cubic.csv`.

- *Matched signals* (5 replicates each; medians):

  | R | oscillator k̄ | oscillator rms | least squares at the true onsets | knots within ±0.01 of an onset | cubic k̄ | cubic rms |
  |---|---|---|---|---|---|---|
  | 1 | 2.1 | 9.4×10⁻⁴ | 5.8×10⁻⁴ | 0.91 | 0.1 | 0.17 |
  | 2 | 3.6 | 1.1×10⁻³ | 6.8×10⁻⁴ | 0.92 | 0.1 | 0.26 |
  | 4 | 7.9 | 1.6×10⁻³ | 7.9×10⁻⁴ | 0.94 | 0.0 | 0.34 |

  The oscillator posterior puts its knots at the onsets and comes within a factor of two of least squares with the onsets known.
- *Cubic BARS* starts at 0–4 knots (its overdispersed start) and stays there: one knot of a cubic buys nothing against a 20-cycle oscillation, so the one-knot birth moves meet a flat likelihood and are refused. Two of the fifteen matched-signal runs climbed, both at R = 1: to 17 knots (rms 0.021) and to 44 (rms 0.0038). So did three of the five ω-mismatch runs, to 34–61 knots (rms 0.0043–0.070). When the chain does climb, free knots beat uniform ones: 44 free knots reach 0.0038, below uniform least squares at 160. Against the oscillator at R = 1, that is twenty times the knots for four times the error.
- *Post-gate, apart from the search:* cubic least squares on uniform knots for the same signals (`p3-uniform`) reads, median over the matched signals: rms 0.26, 0.13, 0.020, 0.0097 and 0.0062 at 10, 40, 80, 120 and 160 knots. Its approximation needs about 80 knots for 0.02 and 160 for 0.006, still four to seven times the oscillator's error with 2–8 knots.
- *Mismatch* (R = 2, the basis at the nominal γ and ω):
  - With ω × 1.05 the oscillator fit uses 20 knots, the count prior's cap, for an rms of 6.8×10⁻³ (cubic BARS 0.070).
  - With γ × 1.5 it uses 19.3, for 3.5×10⁻³ (cubic BARS 0.185, at about zero knots).
  - Least squares at the true onsets is poor under both, 0.097 and 0.035: one knot per onset no longer spans the signal, and the free knots buy the shortfall back at up to ten times the count.
- *The threshold as stated:* the error clause passes, with errors two orders below cubic BARS's in every matched cell, and so does the mismatch clause (6.8×10⁻³ against 0.070, 3.5×10⁻³ against 0.185). The knot clause (at most a quarter of the cubic fit's knots) cannot be met against a cubic fit near zero knots. The gate stands as written and is recorded as not met, for a reason that belongs to the cubic sampler's search. The post-gate least-squares reading puts the saving at twenty to eighty times in knots, at lower error.
- *Decisions:*
  - On this evidence an oscillatory term type saves knots by more than an order of magnitude when its operator matches the signal. Whether it earns its own brief is the user's call, since the pre-stated test did not decide it.
  - It needs no admission rule on spacing in its one-sided form (the sweep). Its weakness is the operator's parameters: a 5% carrier error spends the knot budget. That argues either for working on a band after mixing down, where the carrier is known, or for making γ and ω part of the term (sampled or profiled), not fixed.

**P2: fails its gate on every kink family, not marginally. At these sizes the free-knot error reaches a floor set by the noise before the theory's rate can show, and the knot count is no event counter.** Rows and summary in `receipts/spline-space/results/p2-summary.md` and `p2-diagnostic.csv`.

- *Slopes of ln rms on ln k over k ∈ [4, 32],* medians over 5 replicates:

  | family | uniform LS | BARS held | gap | verdict |
  |---|---|---|---|---|
  | kink, E = 1 | −1.48 | −1.07 | −0.41 | fail |
  | kink, E = 2 | −1.24 | −0.91 | −0.33 | fail |
  | kink, E = 4 | −1.15 | −0.73 | −0.42 | fail |
  | kink, E = 8 | −1.24 | −0.71 | −0.53 | fail |
  | kink, E = 0 | −2.03 | −1.07 | −0.95 | reading |
  | f″ jumps, E = 2 | −2.00 | −1.12 | −0.88 | reading |

  Uniform knots follow the theory: −1.2 to −1.5 with kinks (h^{3/2}), −2 on the smooth and f″-jump curves at these k. The held BARS fits flatten instead. Their posterior-mean error stops at 4×10⁻⁴ to 1.2×10⁻³ by k = 32, a few tenths of the noise level σ = 10⁻³, so every gap is negative.
- *At equal k, free knots win only where kinks dominate:*
  - At k = 16, BARS held reads 0.9–2.0×10⁻³ against uniform's 1.3–3.4×10⁻³ on the kink families.
  - At k = 32 the two are within 25% of each other, E = 1 already favouring uniform (5.4×10⁻⁴ against 4.7×10⁻⁴).
  - On the smooth curve uniform wins outright: 1.3×10⁻⁴ against 4.4×10⁻⁴ at k = 32.
  - The held prior did not hold the small counts: at k = 2 and 4 the posterior sat at 4.5–7.0 knots, and those cells read at their actual k̄.
- *Post-gate diagnostic* (`p2-diagnostic`: replicate 0, k = 32, least squares refitted at every 50th retained draw's knots):
  - On the smooth curve, least squares at BARS's own knots reads a median of 5.6×10⁻⁴ (best 3.2×10⁻⁴), no better than the posterior mean's 4.9×10⁻⁴, and both trail uniform knots (1.3×10⁻⁴). The knots, not the averaging, set the error. At an error below the noise, the data cannot tell placements apart, and the posterior wanders among them.
  - On the four-kink curve, BARS beats uniform (6.3×10⁻⁴ against 1.5×10⁻³) and places 2.3 knots within ±0.02 of each kink, against the theory's 3.
- *Free count* (`UniformPrior(80)`):
  - The posterior holds 8.7–11.7 knots for 0 to 8 kinks, rms 1.1–2.3×10⁻³.
  - Knots crowd toward events at 1.7–2.4 times the windows' share of the domain (0.096 against 0.040 at E = 1; 0.557 against 0.320 at E = 8), but only 0.8–1.0 knots per kink and 0.5 per f″ jump.
  - Away from events the knot histogram is far from the density ∝ |b⁽⁴⁾|^{2/9}: total-variation distance 0.61–0.72.
  - k̄ against E has slope 0.29 knots per event and residual sd 2.0 knots, so E read from k̄ carries an sd of 6.8 events.
- *Decisions:*
  - The measured rates do not separate as the theorems say at the sizes the engine runs (n = 2000, σ = 10⁻³, k ≤ 32). The order-m rate needs errors well above the noise; here the free-knot posterior hits the noise floor by k ≈ 8–16, where uniform knots are still kink-limited.
  - Posterior knot density does not track the curve's variation closely enough for a knot count to estimate how many events a signal holds. Under the free count, a kink costs 0.3 knots on average.

### What each decision came to

| Probe | Decision | Answer |
|---|---|---|
| P1 | Does the built rule refuse exactly full-rank designs, or admit singular ones? | It refuses 78% of the full-rank designs with N > n, nearly all well conditioned, at the column-mass gate. It admits no singular design |
| P1 | Could the matching test stand in front of the factorization? | Yes, as an exact rank test in O(N + n), with the condition check kept behind it |
| P1 | What do scattered product designs need beyond one dimension? | A bipartite matching of sites to product supports is necessary; neither it nor per-axis matching is sufficient on sites in degenerate position (a line). A genericity condition is the missing part, and the numerical rank stays the decider |
| P2 | Do the rates separate at the engine's sizes? | No: the free-knot error floors near the noise before the order-m rate shows |
| P2 | Can a knot count be read as an event count? | No (sd of 6.8 events) |
| P3 | Does an oscillatory term save enough knots? | It saves twenty to eighty times in knots at lower error when its operator matches. The pre-stated test did not decide it, because its knot clause compared against a cubic BARS fit that collapsed. The user's call |
| P3 | Does it need a spacing rule, or a band after mixing down? | No spacing rule in its one-sided form. It is sensitive to its operator's parameters, so: a band after mixing down, or γ and ω inside the term |
| P3 | Does S8's half-period degeneracy match the numbers? | Yes for the positive local element; the design's conditioning shows nothing at that spacing |

### Primitives proposed

Stated for any application and not built.

1. **A support rule on supports, not mass, with the ordered matching in front of the factorization.**
   - *The gap:* the column-mass gate refuses most full-rank designs.
   - *The proposal:* the S7 matching (Theorem 4.61) decides a B-spline design's exact rank from the knots and the sorted sites in one pass, before any design is built. A declared support rule can then be stated as "one site per support, in order", and the mass rule kept only where it is wanted as a prior.
   - *Why it reaches further:* the same rule is the one a signed basis needs, since the mass gate had to be turned off for the oscillator basis. Products would take the bipartite matching as their necessary pre-check, the rank still deciding.
   - *What it would change:* the free-knot sampler's support would stop excluding well-conditioned configurations, and admission would cost less.
2. **An L-spline term family:** the one-sided Green's-function columns of a declared linear differential operator (here D² + 2γD + γ² + ω²), with the operator's parameters sampled or profiled rather than fixed. It runs under BARS's `IBasis` today in probe code. The engine lacks a basis with continuous parameters inside the sampler and a support rule for signed columns (primitive 1). It would let a free-knot fit spend one or two knots per resonance instead of tens.
3. **Births that can climb a flat likelihood.**
   - *The gap:* cubic BARS stayed near zero knots on 13 of 15 oscillatory signals, because a single new knot buys nothing there.
   - *The proposal:* block births and deaths (several knots at once), or a start scaled to the data rather than the fixed `StartDispersion` of 4. The ensemble's built tempering ladder may already serve and is the first thing to try.
   - *What it would change:* free-knot fits of any signal whose first knots carry no likelihood gain.

### Left unread or unrun

- *P1:* only the Normal model's admission (Poisson and Whittle factor weighted matrices); BAPS's `TensorProductDesign` was not called (the product designs were built in probe code); knot multiplicities beyond the clamped ends, which BARS's open simplex excludes; degrees above 3.
- *P2:* noise levels other than 10⁻³, though the rate-against-floor question turns on noise; tempered BARS; window widths other than ±0.02.
- *P3:* cubic BARS with a larger start or tempering; S8's positive local basis under BARS; γ and ω sampled; orders above 2.
- *The source:*
  - Lemma 10.16's statement is damaged in the user's conversion: the "<" of "τ₁ < … < τ_d" and the sentence around it are lost to markup (`schumaker_012_converted.md`, near line 359). That is a conversion defect, beside the book's own misprint already in `errata.md`.
  - A suspected new misprint, unverified: § 4.8's opening sentence reads "y₁ < y₂ < ⋯ < y_{n+m} … with y_i < y_{i+m}", where the second condition and (4.133)'s case y_i = y_{i+m−1} both need ≤ in the first. The PDF's text layer cannot tell ≤ from <, and no page renderer is installed here, so it is reported to the user for checking against the printed page, not added to `errata.md`.

### Doc changes owed

- `docs/Maths/Regression/BARS.md`:
  - The support gate (column mass ≥ 1) refuses most exactly full-rank designs (78% of the full-rank designs with N > n here, nearly all well conditioned) and admits no singular one.
  - Theorem 4.61 (Schoenberg–Whitney) is the exact rank condition, decidable in O(N + n).
  - The gate assumes nonnegative columns, so a custom `IBasis` with signed columns needs it off.
  - Single-knot births can stall on signals whose first knots buy nothing (oscillatory ones).
  - At noise-level errors, a free-knot posterior's count is not an event count.
- `docs/Maths/Regression/BAPS.md`, or wherever scattered tensor designs are documented: per-axis Schoenberg–Whitney is necessary and not sufficient, and a site-to-support matching is necessary and not sufficient in degenerate site geometry.
