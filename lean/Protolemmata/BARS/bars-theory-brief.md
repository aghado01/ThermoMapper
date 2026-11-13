# BARS: the theory to formalize

Status: protolemma brief, 2026-11-13. Nothing here is compiled.

This brief gathers every mathematical claim the BARS engine rests on, says where each one stands in the Lean estate, and adds the ones that are missing. It covers one method. BAPS and BAMRS have their own theory and are out of scope, so that an application can lean on one set of guarantees.

BARS is a general curve engine. Nothing below assumes what the curve measures.

## How to read it

Each statement carries one status, in the vocabulary the revisions note asked for:

| Status | Meaning |
|---|---|
| PROVED | Verified in Lean today |
| DRAFTED | Proof written in `BARS.v2.lean.md`, not yet elaborating |
| OBLIGATION | Statement exists or is given here; proof owed |
| AXIOM | Cited theorem, to be stated and not proved |
| NEW | Not in the estate before this brief |
| CITED | Left as a citation or arithmetic; no Lean target |

In quoted Lean, `p` is a polynomial, as the estate writes it. In prose, p is the spline degree, k the knot count, ξ the knots, ν = k + p + 1 the basis dimension, m the number of observations, θ the relative prominence and ρ the drop fraction, as [BARS.md](../../../docs/Maths/Regression/BARS.md) writes them.

## What exists

| File | Holds | Status |
|---|---|---|
| [`Enthymemata/BARS.lean`](../../Enthymemata/BARS.lean) | First statements of MP-1, MP-2, SP-1, SP-2, SP-3 | Six open proofs; superseded by the draft below |
| [`BARS.v2.lean.md`](BARS.v2.lean.md) | Interval-free MP-1, SP-2 with its dual, bracket width, SP-3 split, MP-2′ witness | DRAFTED. Blocked by one renamed lemma and two measure rewrites |
| [`spc-bars-readout-obligations.md`](spc-bars-readout-obligations.md) | O1 argmax, O2 globalization, O3 policy-bearing readout, a six-step discharge order | OBLIGATION |
| [`Fable5-BARS-revisions.md`](Fable5-BARS-revisions.md) | Audit of the first file; `count_of_mean_unconstrained`, `peakCount_unstable`; merge tree as the carrier | OBLIGATION |
| [`bars-multipeak-lemmas.md`](bars-multipeak-lemmas.md), [`bars-span-lemmas.md`](bars-span-lemmas.md) | MP-1 to MP-4, SP-1 to SP-3 in prose | Source notes |
| [`Enthymemata/Stability.lean`](../../Enthymemata/Stability.lean) | `diagram_stability` (AXIOM), `confidence_pushforward` (PROVED), `landscape_stable` | One open proof |
| [`Enthymemata/BassStructureChain.lean`](../../Enthymemata/BassStructureChain.lean) | `mhRatio_mul_swap`, `not_isMhPair_of_mul_ne_one`, `mh_detailed_balance` | No open proof; generic, usable here |

The renamed lemma is `Polynomial.finite_setOf_isRoot` (`Mathlib/Algebra/Polynomial/Roots.lean`).

## What this pass found

These are corrections to the estate or the code's stated behaviour. Each is taken up in its layer below.

1. **The coded prominence is not persistence.** It differs for a peak with no higher ground on one side, and it jumps at the domain boundary (B3). The claim "the thresholded count is stable" does not hold for the coded definition as it stands.
2. **A fold does not change the gated count.** MP-3 lists pair annihilation as a jump of `SignificantPeakCount`. It is a jump of the raw count only (B5).
3. **Extrema of f″ on a cubic sit on knots.** Which feature orders belong to the curve depends on the degree (B0).
4. **`Stability.lean` reads sublevel sets; the readout reads peaks.** The file also lacks the tameness hypothesis and a domain window (D1).
5. **The argmax witness in `fable-BARS-lemma.md` does not state an inequality**, as O1 already records. A two-parabola witness replaces it (C3).
6. **No derivative evaluation exists in code.** [BARS.md](../../../docs/Maths/Regression/BARS.md) says derivatives use the conditional mean "by extension"; `BSpline` has basis values only. Layer A is the theory that capability needs.
7. **The root solvers use absolute tolerances.** The invariance statements of B8 hold for the mathematics and not for the code until the numerical contract of O2 item 5 is scale-relative.

## Layer 0: the carrier

Every later layer is stated on one carrier, so it is defined once.

- **A curve** is a finite ascending break set t₀ < … < t_n and one real polynomial per span, of degree at most p, with a stated order of agreement at each interior break: the pieces and their first r derivatives agree there.
- **A spline of degree p with simple knots** is the case r = p − 1.
- **Why pieces and not B-spline coefficients.** Mathlib has no B-splines (searched the local checkout). Pieces are closed under differentiation, integration and the difference of two draws, and every readout statement needs only pieces. The B-spline form is a bridge, taken last (A4 to A8), and it is what ties Lean to `BSpline.EvaluateBasis`.
- **A draw** is (ξ, β, σ²), and the posterior is a measure on draws. Layers B and D need one curve; layer C needs the measure.

## Layer A: spline calculus

All NEW. These are the guarantees about derivatives, stated so the readout can use them.

| | Statement | Route |
|---|---|---|
| A1 | A curve with agreement order r is r times continuously differentiable | Pieces are polynomials; agreement at breaks glues one-sided derivatives |
| A2 | The derivative of a curve of degree p and order r is a curve of degree p − 1 and order r − 1 on the same breaks | Differentiate each piece (`Polynomial.derivative`); `natDegree_derivative_lt` |
| A2′ | The p-th derivative is constant on each span and the (p + 1)-th vanishes there | `iterate_derivative_eq_zero` |
| A3 | An antiderivative of a curve is a curve of degree p + 1 and order r + 1 | Integrate each piece; fix constants by continuity |
| A4 | In B-spline form the derivative's coefficients are β′ⱼ = p(βⱼ₊₁ − βⱼ)/(tⱼ₊ₚ₊₁ − tⱼ₊₁) on the knot vector without its two end knots | CITED to de Boor; bridge module |
| A5 | The r-fold knot-aware difference of β is zero exactly when the curve is a polynomial of degree below r | From A4 and linear independence of B-splines |
| A6 | B-splines are nonnegative and sum to one, so min β ≤ f ≤ max β and the sup-norm of f is at most max \|β\| | Bridge; gives a sup-norm band from a coefficient band |
| A7 | For a cubic with simple knots, the values of f″ at the breaks are its degree-1 coefficients. So f is convex on a stretch exactly when those values are nonnegative, and the number of inflections equals the number of sign changes in that sequence | Elementary: f″ is continuous and linear on each span |
| A8 | The design has full column rank exactly under the Schoenberg–Whitney interlacing condition | CITED |

- **Checked numerically this session**, on a clamped cubic basis with 25 even knots: A4 against a central difference (error 5.5e-8, the difference's own accuracy), and A5 for a line at order 2 (2e-13) and x² at order 3 (1e-11). The plain difference of the coefficients is not zero for a line on that basis (h/3 at the first two and last two positions), so A5 needs the knot-aware form.
- **A8 and the support gate.** Column mass at least 1 gives each basis function a data point under it. That does not give the distinct matching Schoenberg–Whitney needs, so the gate does not imply full rank; the conditioning gate (`DesignAdmission`) is what refuses a deficient design.
- **Sign-change bounds in general.** For any degree, f has at most as many sign changes as β (variation diminishing, CITED). A7 is the one case where the bound is an equality, and it needs no citation.

## Layer B: features of one curve

### B0. Which feature orders belong to the curve

| Feature | Zero set | Honest when |
|---|---|---|
| Level crossings | f = ℓ | p ≥ 1 |
| Peaks, valleys | f′ = 0 | p ≥ 2 |
| Inflections, steepest and flattest points | f″ = 0 | p ≥ 3 |
| Extrema of curvature f″ | f‴ = 0 | p ≥ 4 |

NEW. A feature defined by zeros of the j-th derivative needs that derivative continuous and not piecewise constant: p ≥ j + 1. On a cubic, f‴ is piecewise constant, so a maximum of f″ is always a knot or an end. Reading one as a feature of the curve needs degree 4 or more. `SplineBasis` takes any degree.

### B1. Candidate sets and enumeration

| | Statement | Status |
|---|---|---|
| MP-1 | `local_max_is_critical`: a local maximiser of a polynomial is a derivative root, with no interval hypothesis | DRAFTED |
| | `finite_local_max` | DRAFTED (name fix) |
| O1 | `argmax_in_closed_form_set`: a maximiser on a closed span is an end or a derivative root; no non-constant hypothesis | OBLIGATION. Route: `IsMaxOn.isLocalMax` on the interior, then MP-1 |
| O2.1 | The candidate set (ends, breaks, per-span interior derivative roots) is finite when no span is flat | OBLIGATION |
| O2.2 | Every local maximum of a continuous curve lies in it | OBLIGATION |
| O2.3 | Between consecutive candidates the curve is strictly monotone | OBLIGATION. `strictMonoOn_of_deriv_pos`, `strictAntiOn_of_deriv_neg` |
| O2.4 | An interior candidate is a strict local maximum exactly when it exceeds both neighbours | OBLIGATION |
| O2.5 | The enumeration agrees with `CriticalPoints`, under a stated numerical contract | OBLIGATION |
| B1.6 | A span of degree p has at most p − 1 interior critical points and p − 2 inflections | NEW. A nonzero polynomial has at most as many roots as its degree |
| B1.7 | The same five results hold for f′ and f″, with A2 supplying the carrier | NEW, by instantiation |

A flat span is a stated degeneracy with its own result, as O2 says.

### B2. Level crossings and widths

| | Statement | Status |
|---|---|---|
| SP-1 | `level_set_finite` | DRAFTED (name fix) |
| SP-2 | `unique_crossing_in_monotone_bracket` and its antitone dual | DRAFTED, complete |
| | `superlevel_of_peak_bracket`: on a rise-then-fall bracket the superlevel set is the interval between the two crossings | DRAFTED, one open proof |
| | `peak_width_eq_bracket_roots` | DRAFTED, measure rewrite open |
| B2.5 | With ℓ = h − ρ(h − c) and 0 < ρ < 1, c < ℓ < h. Scanning outward from the peak, the first candidate at or below ℓ and its inward neighbour form a monotone bracket straddling ℓ | NEW. From O2.3; this is the step `SignificantPeakSpans` relies on before it calls the crossing solver |
| B2.6 | If no candidate on a side is at or below ℓ, the curve stays above ℓ to the domain end: the span is clipped there, and no crossing exists | NEW |

`levelWidth` keeps the name the draft gives it: occupancy above a level, a width only on a single-peak bracket.

### B3. Prominence and the gate

This is the definition the code runs, read from `SignificantPeaks` and `ValleyFloor`. It has to be pinned before O3 can be stated.

- **As coded.** Candidates x₀ < … < x_{n−1} with values fᵢ. For a strict local maximum i, Lᵢ is the least value over the run of candidates to its left that are not higher than fᵢ, stopping before the first higher one or at the end; Rᵢ likewise on the right. At i = 0, Lᵢ = −∞; at i = n − 1, Rᵢ = −∞. Prominence is fᵢ − max(Lᵢ, Rᵢ), and the peak is significant when that is at least θ(max f − min f), with a positive range.
- **Against persistence.** Under the elder rule on superlevel sets, a peak dies at the higher of the cols toward higher ground, and a side with no higher ground imposes nothing. The coded value uses that side's minimum instead. So:
  - with higher ground on both sides, coded prominence equals persistence;
  - with higher ground on one side only, it is at most persistence;
  - for the global maximum, it is fᵢ − max(min left, min right), where persistence is unbounded or the full range.
- **At the boundary it jumps.** A curve that rises to its maximum at the right end has prominence fᵢ − Lᵢ. Move the maximum inward by a perturbation of size δ, so f(1) = max − δ: the right side now contributes f(1), and the prominence is δ. An arbitrarily small change takes the gated count from 1 to 0.
- **Statements.**

| | Statement | Status |
|---|---|---|
| B3.1 | The coded prominence, as a function of the candidate list | NEW definition |
| B3.2 | It equals elder-rule persistence when both sides have higher ground, and is at most it otherwise | NEW |
| B3.3 | The boundary witness above | NEW, one witness |
| B3.4 | Prominence is unchanged by increasing maps of x and scales with positive affine maps of y; the gate is unchanged by both | NEW |
| B3.5 | The global maximum can fail the gate | NEW, one witness: an interior maximum whose curve barely descends to one end |

A search of `PeakTests.cs` for boundary, edge and endpoint finds nothing, so no test pins this behaviour. Which definition is the specification is a decision (below). O3's instruction stands either way: prominence comes from the carrier-generic merge-tree foundation, not from a second definition local to BARS.

### B4. Inflections

All NEW.

| | Statement |
|---|---|
| B4.1 | At an inflection, \|f′\| has a local maximum when f′·f‴ < 0 (a steepest point) and a local minimum when f′·f‴ > 0 (a flattest point: a shoulder) |
| B4.2 | On a cubic span with f′(u) = Au² + Bu + C and Δ = B² − 4AC, the kind is the sign of Δ: f′(u*)·c₃ = −Δ/12 |
| B4.3 | On a strictly monotone stretch between two consecutive critical points, with f′ having only simple critical points there, steepest points outnumber flattest points by one |
| B4.4 | The inflection count of a cubic is the sign-change count of A7 |

B4.3 is Rolle and counting on f′, which vanishes at both ends of the stretch (`exists_deriv_eq_zero`).

### B5. The cubic normal form and the fold

All NEW; algebra on one polynomial.

| | Statement |
|---|---|
| B5.1 | A cubic with c₃ ≠ 0 is f(u*) + m(u − u*) + c₃(u − u*)³ about u* = −B/(2A), with m = −Δ/(4A). It is point-symmetric about its inflection |
| B5.2 | For Δ > 0 the critical points are u* ± √Δ/(2\|A\|) and their values differ by Δ^(3/2)/(6A²) |
| B5.3 | For Δ < 0 there is no critical point and the least \|f′\| is \|Δ\|/(4\|A\|) |
| B5.4 | A pair born at Δ = 0 on a monotone stretch has prominence equal to the height difference of B5.2, which is 0 at birth. So for θ > 0 the gated count is unchanged across the fold |

- B5.1 is the depressed cubic of Cardano's method, so the crossing solver and the feature reading share one normal form.
- B5.4 corrects MP-3's event (a). One feature then has two posterior probabilities: that the pair exists, and that its peak passes the gate.
- MP-3's `sigCount_locally_constant` stays an OBLIGATION, with the regular set now: no prominence on the threshold, no tie between candidate values, and the boundary case of B3. `catastrophe_set_null` stays CITED.

### B6. Features built from the above

These need no new proof beyond being well defined. Each is NEW as a definition.

- **A level between two feature heights**, with its unique crossing on a monotone stretch (SP-2): half-rise points, rise widths, widths at any fraction, asymmetry ratios.
- **A tangent at a feature point and the intersection of two tangents.** Well defined when the slopes differ; at a steepest point of a strictly monotone stretch the slope is nonzero.
- **An integral of f·xʲ between two feature points** (A3): area, centroid, second-moment width.

### B7. Whole-curve sums

NEW.

| | Statement |
|---|---|
| B7.1 | Total variation is the sum of \|fᵢ₊₁ − fᵢ\| over consecutive candidates (from O2.3) |
| B7.2 | On a cubic span of width h with end values a, b of f″: ∫(f″)² = h(a² + ab + b²)/3 |
| B7.3 | ∫\|f″\| = h(\|a\| + \|b\|)/2 when a and b share a sign, and h(a² + b²)/(2(\|a\| + \|b\|)) otherwise |

### B8. Invariance

NEW. One equivariance lemma for g(x) = a·f((x − x₀)/s) + b with a, s > 0, and its corollaries:

| Change of axes | What is preserved |
|---|---|
| Any increasing map of x | Critical values, prominence, the gate, total variation |
| Affine in x and in y separately | Also the candidate set (mapped), widths (scaled by s), inflections and their kinds, the sign of Δ |
| One shared unit | Also geometric curvature and arc length |

These are the statements that keep the engine free of any application's units. They fail in code wherever a tolerance is absolute.

### B9. Counts

| | Statement | Status |
|---|---|---|
| | `count_of_mean_unconstrained (m n k)`: any three counts are attainable by two curves and their mean | OBLIGATION (revisions note) |
| | `peakCount_unstable`: a sup-norm perturbation of any size can add any number of peaks | OBLIGATION |
| | The gated count is sandwiched between gated counts at θ ± 2δ under a δ perturbation | OBLIGATION, downstream of D1, and subject to B3 |

## Layer C: pooling over draws

| | Statement | Status |
|---|---|---|
| C1 | A readout affine in the curve commutes with the posterior mean, and its conditional mean given the knots is exact | NEW. Linearity of the integral, tower property |
| C1′ | A readout that commutes with every two-point average is midpoint-affine | NEW. The converse that justifies "all others use the draw" |
| C2 | Var f(x) = Var(conditional mean) + E(conditional variance) | NEW as an instance. Mathlib has it: `integral_condVar_add_variance_condExp`. This is `FitSd` (`BarsEnsemble.cs:560`) |
| C3 | Argmax, count and width do not commute with averaging | See below |
| C4 | Intensity and coverage | See below |
| C5 | A count of one in every draw does not imply a unimodal location posterior | OBLIGATION (MP-2′, measure form) |
| C6 | f′(x) and f″(x) are affine in β, so C1 and C2 apply with derivative design rows | NEW. `ConditionalVariance` already takes any rows (`BarsEnsemble.cs:974`) |

- **C3, the witnesses.**
  - Argmax, NEW: f = −a(x − u)², g = −b(x − v)². The mean peaks at (au + bv)/(a + b), which is not (u + v)/2 when a ≠ b and u ≠ v.
  - Count: MP-2 (`mp2F`, `mp2G`) and MP-2′ (tilted corners), both DRAFTED with routes.
  - Width: `levelWidth_not_affine` (DRAFTED) and `fwhm_of_mean_ne_mean_of_fwhm` (DRAFTED, route given). The branch-relative form with ℓ = h − ρ(h − c) is an OBLIGATION, restricted to the carrier (ledger step 5).
- **C4, as coded** (`PeakCounts`, `CoveredGridPoints`).
  - λ at a grid cell is the mean number of significant peaks whose location is nearest that cell. Summed over the grid it equals the mean count exactly, since each peak lands in one cell. NEW, a finite identity.
  - λ is additive over disjoint sets (`peak_intensity_additive`, OBLIGATION, MP-4).
  - π at a grid point is the probability that the point lies in some significant span, counted once per draw, so 0 ≤ π ≤ 1. NEW.
  - Binning moves a location by at most half a cell. NEW; it is the only approximation in the pooled fields.
  - A fixed-K vector of peaks has no canonical labelling when the count varies (`fixedK_readout_ill_posed`, OBLIGATION).
- **`PeakModes`** reads the pooled fields, and nested modes share one span. It is a summary of λ and π, with no per-draw statement to prove.

## Layer D: stability

| | Statement | Status |
|---|---|---|
| D1 | `diagram_stability`: the diagram moves no further in bottleneck distance than the curve in sup-norm | AXIOM, to be restated |
| D2 | `stability_event_subset`, `confidence_pushforward` | PROVED |
| D3 | The gated-count sandwich of B9 | OBLIGATION |
| D4 | A band on f gives no band on f′; features of f′ need a band on f′ itself | NEW, one witness (a small fast oscillation) |
| D5 | `landscape_stable` | Open; an AXIOM by the estate's own policy until the landscape is defined |

- **D1 needs three repairs.**
  - A tameness hypothesis, which B1's finiteness results discharge for a curve.
  - The sup-norm over the domain window, not the line.
  - A statement of which sets it reads. The docstring says sublevel; peaks are superlevel features, so the readout applies it to −f.
- **D4 and C6 together** are how derivative features get a guarantee: the band on f′ comes from the posterior of f′, through C6.
- **A bound on f′ from f** exists inside one spline space with knot gaps bounded below (a Markov-type inequality, CITED). BARS has no gap rule, so the brief does not lean on it.

## Layer E: the model and the sampler

The obligations ledger left the sampler out. This is the map, taken from [BARS.md](../../../docs/Maths/Regression/BARS.md); I did not re-read the move code in this pass.

| | Statement | Status |
|---|---|---|
| E1 | Given k, the knots are the order statistics of k uniforms: density k! on the ordered simplex. The target is log p(k) + log k! + the marginal, so p(k) is the count law | NEW |
| E2 | The Normal marginal's scale satisfies a = ‖ỹ − Zβ̂‖² + bᵀG⁻¹b/(m + 1), so a > 0 unless ỹ = 0 | NEW, linear algebra |
| E3 | The birth and death ratios multiply to one: −log(k + 1) − log q(ξ*∣ξ) and +log(k + 1) + log q(ξ*∣ξ) | NEW; an instance of `mhRatio_mul_swap` |
| E4 | Each move pair satisfies detailed balance | NEW; instances of `mh_detailed_balance` |
| E5 | The DMGK schedule gives b_k p(k) = d_{k+1} p(k + 1) | NEW; r·min(1, 1/r) = min(1, r) |
| E6 | Setting the target to −∞ off an admissible set keeps each pair reversible for the restricted target | NEW |
| E7 | Under the support gate the count law is p(k)s(k)/Z, with s(k) the chance that k uniform knots pass | NEW; `CountLawTests` checks it by sampling |
| E8 | A mixture of reversible kernels with a common target is reversible for it | NEW |

CITED, no Lean target: the Gaussian integral behind the marginal, the swap acceptance and the stepping-stone identity of tempering, ergodicity and mixing, posterior contraction, and the approximation rate of splines.

## Order of work

This extends the ledger's six steps and keeps their order.

1. Land the revised per-polynomial statements with the name fix, O1 and the two-parabola witness.
2. Build the carrier with A1 to A3, and prove O2 on f, f′ and f″.
3. Prove B4, B5 and B7, which need only the carrier.
4. Build the merge-tree foundation, settle B3, and state O3 with B2.5 and B2.6.
5. Restrict the non-commutation witnesses to the carrier, and add the measure carrier for C1 to C6.
6. Repair D1, then derive the sandwich.
7. Build the B-spline bridge (A4 to A8) and state the sampler instances (E).

## Decisions owed

- **Which prominence is the specification** for a peak with no higher ground on one side, and at the boundary: the coded one, with its instability stated, or persistence.
- **The boundary policy**, interior only or boundary included, which O3 already names as an axis.
- **Whether the B-spline bridge is in scope.** Without it Lean certifies curves; with it Lean certifies what `EvaluateBasis` returns.
- **The degree to run** when curvature extrema are wanted.

## Reading list

Citations are from memory; confirm editions before acquiring.

| Source | Needed for | Held |
|---|---|---|
| de Boor, *A Practical Guide to Splines* | A4, A5, A6, A8, variation diminishing, knot insertion | Not on the textbook shelf |
| Schumaker, *Spline Functions: Basic Theory* | Zeros of splines, sign-change counts, the Markov-type inequality | Not on the shelf |
| Piegl and Tiller, *The NURBS Book* | The derivative algorithm beside A2.2, which the code already cites | Not on the shelf |
| Cohen-Steiner, Edelsbrunner and Harer (2007), stability of persistence diagrams | D1's exact hypotheses | Not checked |
| Edelsbrunner and Harer, *Computational Topology* | Elder rule and extended persistence, for the B3 decision | Not on the shelf |
| Dey and Wang, *Computational Topology for Data Analysis* | Merge trees, stability | On the shelf (`DeyWang2021`) |
| Fasy et al. (2014), confidence sets for persistence diagrams | D2's prior art | Listed to acquire in [literature-followup.md](literature-followup.md) |
| Green (1995); DiMatteo, Genovese and Kass (2001) | E3 to E5 | Held, per that note; not re-checked |
| Denison, Mallick and Smith (1998) | The uniform birth | Not checked |
| Tierney (1998), Metropolis–Hastings kernels on general state spaces | E4, E6, E8 in measure-theoretic form | Not checked |
| Liang et al. (2008); Kass and Wasserman (1995) | E2 and the unit-information prior | Liang held; Kass and Wasserman not checked |
| Gasser and Kneip (1995); Kneip and Gasser (1992) | Prior art for a pooled intensity of feature locations | Not checked |

## Not verified

- Nothing was compiled and `lake build` was not run. Mathlib names were checked for existence in the local checkout only.
- B3 is from reading `SplineExtrema`; no test was run.
- The numerical checks of A4 and A5 used a standalone script, not the project's code.
- [BARS.md](../../../docs/Maths/Regression/BARS.md) was read in full to line 320 and by headings after that. `SplineTarget`, `KnotMoves` and `KnotBirthDensity` were not re-read in this pass.
- `fable-BARS-multipeak-lemmas.md` was read by headings. Its sufficiency premise (BARS-S) is not covered here.
- The three `bgp-*` notes in this folder are not BARS and were not read.
