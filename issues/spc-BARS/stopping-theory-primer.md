# Stopping conditions for SPC × BARS under dyadic refinement — a theory primer

Date: 2026-10-24. Status: a briefing primer for a separate exploration thread. It records no decisions. The engine-side
stopping strategy (OQ7, OQ8) stays parked until the SPC readout design is complete (bass-preparation-brief.md, "Parked
decisions").

**Sources.** The primer consolidates three discussions and adds one synthesis of its own (§5.3):
- [opus-ESS-bars-spc-etc.md](discussions/opus-ESS-bars-spc-etc.md): the dyadic scheme, the two regimes, and the
  unsampled-feature problem;
- [opus-fisher-zeros-thermo-lengths.md](../../lean/Protolemmata/opus-fisher-zeros-thermo-lengths.md): feature
  concentration, conformal Bayes, Fisher zeros, thermodynamic length, and the a priori point counts;
- [argon-fisher-notes.md](../../lean/Protolemmata/argon-fisher-notes.md): a gate-by-gate stopping protocol and a
  critique of the note above.

Each claim carries its status:
- **checked here:** verified while writing this primer;
- **derived, unverified:** derived in conversation, not yet checked;
- **disputed:** the sources disagree;
- **withdrawn:** retracted in the sources.

---

## 1. The setting

**Two regimes doing unrelated jobs** (from the dyadic-scheme discussion):
- **Regime A, location.** Swendsen–Wang runs on a uniform, nested temperature grid. BARS is refitted on the union of
  every point sampled so far. Its job is the joint posterior over the whole thermal curve, cheaply, to find the
  regions of interest.
- **Regime B, accuracy.** Dense, expensive SW at those regions, with richer observables. Clustering happens only here.
- **"Densify" means different things.** In A it is grid refinement; in B, replica and equilibration depth. An argument
  that slides between the two is wrong wherever it lands.

**The dyadic schedule.**
- The bracket is normalized to [0, 1]. §3 argues the normalization should be in ln T.
- It starts at S₀ = {0, 1}. Each volley adds the complementary midpoints, and BARS refits on the union Uᵣ: never
  composing earlier fits, and never discarding a point.
- With N independent replicates per point, the spacing after v volleys is h_v = 2^(−v), uniform at every volley.
- **Why refitting on the superset matters.** Successive fits are directly comparable, with no compounded
  approximation error, and early knot configurations don't constrain later ones. So agreement across volleys is a real
  convergence signal, not a fit agreeing with itself.

**The bracket** is expansion-only on a log axis, keeps and remaps every prior observation, and extrapolates nothing.
That forecloses the degenerate direction of narrowing until everything looks resolved.

**Replicates** give the observation variance per point empirically and without a model. That is the natural input to
BARS's per-point weights.
- **The caveat:** replicates measure precision, not accuracy.
- **Why it is mostly moot here:** SW is the standard cure for critical slowing down, so its residual bias is small and
  measured rather than modelled.

**The closed-form readout is what licenses the sparse grid.** Peak and span locations are read analytically from each
posterior draw, from knots, endpoints and derivative roots. So feature location is limited by the posterior, not by
the grid spacing.

---

## 2. What a stop has to certify

**The core problem.** A feature that falls entirely between grid points is invisible to every internal diagnostic:
- the neighbouring points read as baseline;
- BARS fits a smooth curve through them, with no knots and no span to flag;
- the posterior is confidently featureless.

Dyadic refinement finds the feature eventually, since h halves each volley. But **the stopping rule may fire first,
because everything visible looks fine.** No criterion computed from the samples alone can close this gap.

**So a stop certifies "no features narrower than w_min",** not "no unresolved features." That claim is honest,
falsifiable, and calibrated by measurement.

**The pieces that compose** (the Fisher-zeros note):
- **Concentration** of the features that define regions of interest, computed per draw as the readout already does:
  - P(K = k̂) ≥ 1 − α for the significant-peak count;
  - each peak location's credible width below the window the dense phase will sample;
  - the same for the span endpoints.
- **Resolution**, the clause above: h_v small against w_min, so each span contains several grid points.
- **Calibration**, by conformal Bayes (Fong & Holmes 2021): full conformal intervals from the existing posterior draws,
  by importance-weighting each draw with the left-out point. That needs exactly the per-draw pointwise log densities
  that `BarsEnsemble.PointwiseLogDensity` now provides (OQ17, built 2026-10-24).
  - **What it gives:** a check that the posterior's confidence is calibrated, not a stopping criterion.
  - **Its limit:** it covers χ at a temperature, not derived features such as peak locations.
- **The learning curve** across volleys, as a report and not as the criterion: concentration needn't be monotone, since a
  later volley can reveal a second peak.

**Two knobs, one symptom** (the dyadic-scheme discussion). A wide posterior has two causes with opposite remedies:
- the grid is too coarse;
- the runs at existing points are too short.

Replicate spread separates them. Large spread at measured points means the measurements are inadequate; small spread
with a wide posterior between points means the grid is.

---

## 3. Thread I — a priori resolution from the graph: the discovery count

Status: **derived, unverified** (the Fisher-zeros note), except where marked.

**The giant-spin switch.** Treat a cluster as one rigid spin, switching between two states with energy gap ΔE and
entropy gap ΔS = ln g, where g is the ratio of state counts. Then:
- **The switch temperature:** T* = ΔE/ΔS.
- **Its 25–75% width:** 2 ln 3/ΔE in β, which is about 2.2/ΔS in ln T whatever ΔE is. **Checked here.**
- **The nearest Fisher zero:** Im β₁ = π/ΔE, at the zeros of 1 + g e^(−βΔE). **Checked here.**
- **Its thermodynamic length.** With p the high-state probability, σ_E dβ = −dp/√(p(1−p)), so traversing p from p₁ to p₂
  costs 2(arcsin√p₂ − arcsin√p₁), whatever ΔE and ΔS. **Checked here, with a correction to the source note.**
  - **The note's "exactly π"** is the whole real β-line, p from 0 to 1, which includes negative temperatures.
  - **At physical temperatures** p runs from 0 (cold) to g/(1+g) (T → ∞). The full length is 2 arcsin√(g/(1+g)): 2.69 at
    g = q − 1 = 19, approaching π as g grows.

**SPC's two kinds of feature:**
- **A split** (A ∪ B breaks into A and B; J_AB is the total coupling between them) has ΔS = ln(q − 1). So its
  T* = J_AB/ln(q − 1), and its width in ln T is about **w_q = 2.2/ln(q − 1)**, roughly 0.75 at q = 20, whatever the
  cluster sizes.
  - No bump in χ can be narrower than w_q.
  - A cluster that lives less than about w_q in ln T shows up as a shoulder, not a bump.
- **A melt** (a cluster of n_c spins going free) has ΔS ≈ (n_c − 1) ln q, so it is sharp. But it is a step at the
  plateau's hot edge, and a coarse grid can only narrow a step down to one interval; it can't lose it.

**The count.**
- **The narrowest feature.** w_min = min over significant clusters of max(lifetime, w_q). A lifetime is
  ln(T_breakup/T_separation), read from PKWang's closed-form ladder (`PKWang.EstimateBracket` already uses that
  ladder's ends to set the window). In ln T, the constants that are hard to pin down cancel: ln 2 against ln(q − 1),
  the 1/K̂ gauge, and multiplicative mean-field shifts.
- **The sampling constant is geometric.** Spacing h ≤ w/2 guarantees at least two samples above half-maximum on any bump
  of half-max width w.
- **On a grid evenly spaced in ln T** over a window of E e-folds:
  n_min = 1 + ⌈2E/w_min⌉ ≤ 1 + ⌈2E·ln(q − 1)/2.2⌉.
  That is 14 points for E = 4.6 (a 100× window) at q = 20, for any graph.
- **Under the dyadic schedule** from S₀ = {0, 1}, the floor is v* = ⌈log₂(2E/w_min)⌉ volleys, which is v* = 4
  (17 points) for that case.
- **Evenly spaced in T instead,** the cold end forces about 260 points for the same window. So normalize ln T, and take
  geometric midpoints.

**Idealizations to pressure-test:**
- **Clusters aren't rigid.** Real clusters peel apart gradually. That widens features, so the count stays an upper bound.
- **k-way simultaneous splits** have ΔS = (k − 1) ln(q − 1), so they are k − 1 times narrower.
- **Mean-field shifts don't fully cancel.** They cancel in lifetimes, but not in the window E, so pad E.
- **A size floor for "significant" clusters is required.** Without one, the weakest cut isolates one loose point and E
  blows up.

---

## 4. Thread II — analyticity and the dense count

Status: **derived, unverified** (the Fisher-zeros note); its operational use is **disputed** (the argon notes).

- **Why the curve is smooth.** For N spins, Z(β) is a finite sum of exponentials, so its zeros (Fisher 1965) lie off the
  real axis. χ is analytic in a strip whose half-width is set by the nearest zero, Im β₁.
- **The spacing.** Cauchy estimates bound every derivative there, roughly |χ⁽ⁿ⁾| ≲ n!·M/(Im β₁)ⁿ. Fed into the cubic
  spline's h⁴ error bound, the spacing for accuracy ε is h ≈ a·(16ε/5)^(1/4), with a the strip width.
- **The dense count.** For a melt, a ≈ π/((n_c − 1) ln q) in ln T. A 200-spin cluster at ε = 1% then needs about 2,000
  points across the window, but only near the melt's edges. That is regime B's job.
- **Finite-size scaling.** Im β₁ ~ L^(−1/ν) on lattices (Itzykson, Pearson & Zuber 1983). Data graphs have no known
  exponents, which is an argument for estimating from data rather than prescribing from theory.

**The dispute** (argon notes):
- **The claim:** Im β₁ can't be estimated from real-axis Monte Carlo before the peak is already resolved, since it needs
  analytic continuation. And χ⁽⁴⁾ is a sixth-order cross-cumulant, too noisy to estimate in practice.
- **The response here:** both points hold for estimating a from sampled data. Neither touches Thread I's a priori route,
  which takes switch widths from the graph alone, before any sampling. The Fisher-zero bound remains the explanation of
  why refinement converges, and the dense count's order of magnitude.

---

## 5. Thread III — information geometry: lengths and divergences

**5.1 Thermodynamic length** (Salamon & Berry 1983; Crooks 2007).
- The Gibbs family's Fisher information in β is Var(E), so ds = σ_E dβ.
- The length L = ∫σ_E dβ counts statistically distinguishable temperatures along the path.
- **In practice:** neighbouring temperatures need overlapping energy histograms, Δβ ≲ c/σ_E. That is the replica
  exchange and Ferrenberg–Swendsen reweighting condition.
- **The count it implies** grows like √N, because it promises the whole ensemble at every temperature. So it is the
  wrong count for discovery (Thread I), and the right tool for placing regime B's points.

**5.2 The Jeffreys divergence across an interval** (argon notes; the identity is standard, **checked here**).
- For the Gibbs family, KL(P_{β₁} ‖ P_{β₂}) + KL(P_{β₂} ‖ P_{β₁}) = (β₂ − β₁)(⟨E⟩_{β₁} − ⟨E⟩_{β₂}) =: D_J. It needs only
  the endpoints' mean energies, which every SW run gives with replicate error bars.
- **The argon notes' gate:** D_J ≤ c_J² with c_J ∈ [1, 2], plus a quadrature-closure ratio for curvature.
- **Their bounded-excursion form:** since d⟨A⟩/dβ = −Cov(A, E), Cauchy–Schwarz gives
  |Δ⟨A⟩| ≤ √(∫Var(A) dβ) · √(Δ⟨E⟩).
- **The gap in that form:** as written it needs ∫Var(A) across the unsampled gap, which is circular. 5.3 closes the gap
  for bounded observables.

**5.3 A synthesis: an endpoint-only certificate against hidden switches** (new here). Status: the inequalities are
**checked here**; the corollary rests on one approximation, flagged below.

1. **Length is bounded by the divergence.** Cauchy–Schwarz on L = ∫σ_E dβ over [β₁, β₂]:
   L² ≤ (β₂ − β₁)·∫σ_E² dβ = (β₂ − β₁)(⟨E⟩_{β₁} − ⟨E⟩_{β₂}) = D_J.
   Equality holds exactly when σ_E is constant across the interval.
2. **The excursion of any bounded mean is bounded.** For an observable A with values in [a, b], Var(A) ≤ (b − a)²/4.
   Applying the argon inequality on every sub-interval [β₁, β], whose two factors are monotone in β, gives
   sup over β ∈ [β₁, β₂] of |⟨A⟩_β − ⟨A⟩_{β₁}| ≤ ((b − a)/2)·√D_J.
   Nothing unsampled enters.
3. **The corollary.** A switch is hidden in an interval when both endpoints sit on its plateaus: p ≤ θ at one end and
   p ≥ 1 − θ at the other, so neither endpoint sees it. Across the interval the switch then traverses at least
   ℓ(θ) = 2(arcsin√(1 − θ) − arcsin√θ) of length (§3). So L ≥ ℓ(θ), and therefore D_J ≥ ℓ(θ)².
   - **The values:** ℓ(0.05)² ≈ 5.0, ℓ(0.10)² ≈ 3.4, and ℓ(0.25)² = (π/3)² ≈ 1.1.
   - **The approximation:** L ≥ ℓ(θ) assumes the system's Var(E) is at least the switching cluster's own. That holds
     exactly when the cluster's energy is uncorrelated with the rest, which is the giant-spin picture's assumption. For
     coupled clusters it needs checking.
   - **The certificate:** an interval with D_J < ℓ(θ)² cannot hide a rigid switch traversed beyond its θ-plateaus. A
     switch that isn't hidden shows at the endpoints, as a difference in the switching probability.

**What it gives.** A feature-free, per-interval certificate, from endpoint mean energies only. It complements Thread I's
a priori count, since one is computed after sampling and the other before. It places the argon gate (D_J ≤ c_J² ≤ 4):
that gate excludes a switch hidden beyond its 5% plateaus (5.0), but not one hidden beyond its 10% plateaus (3.4).

**What it doesn't give.** A bound on χ's excursion directly. χ is a variance, so the route is through bounded means such
as ⟨m⟩ and ⟨m²⟩, with a factor of N/T that makes it loose for large N. Normalized observables (cluster fractions, label
entropy) suit it better.

---

## 6. Thread IV — sequential evidence across volleys

- **Each volley's midpoints are a free holdout** (argon notes). The new midpoints are space-filling and out of sample for
  the previous fit. Standardized surprises z = (χ̄_SW − χ̂_{v−1})/s_{v−1} at the midpoints test the previous fit's
  predictive band:
  - **under-resolved,** surprises are large wherever a new point lands on an unresolved feature;
  - **resolved,** they behave like N(0, 1).
- **The argon rule:** wait for two consecutive quiet volleys, to reduce quarter-point aliasing.
- **Why that rule needs care** (a suggestion of this primer):
  - **It is an optional-stopping decision** over an unbounded sequence of volleys. Fixed z-thresholds lose their error
    rate under optional stopping.
  - **A fix with a guarantee:** anytime-valid tests (e-processes, test martingales; Ramdas, Grünwald, Vovk & Shafer 2023)
    keep the error rate whenever the rule fires.
  - **What remains open:** whether an e-process over midpoint surprises has useful power here.
- **The aliasing limit is unchanged.** No surprise statistic sees a feature that no point has touched. Thread I's floor,
  or 5.3's certificate, is what bounds that risk.

---

## 7. Claims register

| # | Claim | Source | Status | How to settle |
|---|---|---|---|---|
| C1 | A switch's width is 2 ln 3/ΔE in β, about 2.2/ΔS in ln T | Fisher-zeros note | Checked here | Lean (L2) |
| C2 | Its nearest Fisher zero is at Im β₁ = π/ΔE | Fisher-zeros note | Checked here | Lean (L2) |
| C3 | Its thermodynamic length is π for any ΔE and ΔS | Fisher-zeros note | Corrected here: π over the whole β-line; 2(arcsin√p₂ − arcsin√p₁) for a traversal, 2.69 for the physical range at q = 20 | Lean (L2) |
| C4 | A split's width is w_q = 2.2/ln(q − 1), independent of cluster size | Fisher-zeros note | Partly refuted by E1 (§12.9): FkReduced respects it; variance observables at melts and at the toy's background-mediated split run 3–9× narrower | done (E1) |
| C5 | n_min ≤ 1 + ⌈2E·ln(q − 1)/2.2⌉ on a grid evenly spaced in ln T | Fisher-zeros note | Partly checked (§12.9–12.10): the width floor fails for variance observables (E1); BARS registers bumps at 0.6–2.1 points per FWHM, under the geometric 2 (E4) | E3 |
| C6 | The dense spacing is h ≈ a·(16ε/5)^(1/4) | Fisher-zeros note | Derived, unverified | Numeric, on a known analytic curve |
| C7 | Im β₁ can't be estimated before the peak is resolved | Argon notes | Disputed in scope (§4) | Literature on zeros from Monte Carlo data |
| C8 | D_J is the symmetrized KL divergence between endpoint Gibbs measures | Argon notes | Checked here (standard) | Lean (L1) |
| C9 | L² ≤ D_J | This primer | Checked here | Lean (L1) |
| C10 | A bounded mean's excursion is at most ((b − a)/2)√D_J | This primer | Checked here; tightened in §12.2 (C15) | Lean (L1) |
| C11 | A rigid switch hidden beyond its θ-plateaus forces D_J ≥ ℓ(θ)² | This primer | Checked here, unconditional (§12.1): holds for coupled clusters; θ-plateaus restated against (q − 1)/q | Lean (L1); E2 on toys |
| C12 | Two quiet volleys reduce aliasing | Argon notes | Heuristic, tested (§12.11, §12.14): on the n₀ = 10 schedule it stops once before a bump registers (iris's minor melt) | anytime-valid reformulation (§6); more curves |
| C13 | The global grid's union is non-uniform after refinement | Early dyadic discussion | Withdrawn (the union stays uniform) | none needed |
| C14 | Critical slowing down makes regime-A noise heteroscedastic | Early dyadic discussion | Withdrawn (SW suppresses it; replicates measure the rest) | none needed |
| C15 | A [0, 1] observable's Fisher angle moves at most L ≤ √D_J, so its excursion is at most (b − a)·sin(min(√D_J, π)/2), sharply | §12.1–12.2 | Checked here (proof; enumeration) | Lean (L1) |
| C16 | D_J is blind to size; its false-alarm floor is about 2.8ρ(Δ ln T)² at q = 20 | §12.3 | Confirmed by E2 (§12.9): D_J scales with N, and reaches 127 on intervals with no cluster change | done (E2) |
| C17 | Pair connectivities, so S₂, are monotone in β; ΔS₂ bounds every merge inside an interval | Argon notes (seed); §12.4 | Checked here (Grimmett 3.21; enumeration; parity dump; all 64 steps of the extended run) | Lean (L4) |
| C18 | Superadditive cluster functionals are monotone; F_n counts significant merges | §12.4 | Checked here (proof; enumeration) | Lean (L4) |
| C19 | S₂ at the bracket ends bounds the merge mass outside the bracket | §12.5 | Checked here (a corollary of C17); usefulness open | E5, against the graph's components |
| C20 | S₂ bounds FkReduced: inside an interval, sup \|Δχ_red\| ≤ N·ΔS₂ | §12.7 | Checked here (proof; enumeration) | Lean (L4) |
| C21 | No monotone cluster functional ignores shedding: merge(n, n) ≤ n·max join(a) over a ∈ [n, 2n) | §12.7 | Checked here (proof) | none needed |
| C22 | The toy dump's bracket clips its cold edge | User observation; §12.7 | Checked here (§12.8), and confirmed by the extended run (§12.9): the cold side of the stripe split lies below the old T_lo; not a miscomputed bracket, and not components | done |
| C23 | A spanning tree's S₂ is a closed-form lower bound on S₂(T) | §12.8 | Checked here (comparison inequality; all 48 grid points on both datasets) | Lean (L4) |
| C24 | A split through a sparse bridging population frees those spins too, so it is sharper than a bare split, towards a melt | §12.9 | Derived from E1, unverified | Toys that vary the bridge density |
| C25 | Per-interval ΔS₂ tests lose power at fine volleys; a test that descends the dyadic tree keeps it | §12.9, §12.14, §12.17 | Designed (§12.17: hereditary quiet, α spent over levels, narrowest-over-threshold localization); no quiet certification contradicted by an independent sweep | More datasets; replicate allocation for indeterminate blocks |
| C26 | BARS registers a bump at about one point per FWHM, roughly placed; its prominence gate, not the grid, drops weak bumps | §12.10, §12.14 | Observed (E4 on both schedules, equal weights): 0.2–2.4 points per FWHM, the low end from lucky lattice phases | E3; weighted fits |
| C27 | Registration depends on the lattice's phase as well as its density, so a minimum should plan on the worst phase | §12.14 | Observed (n₀ = 2 against n₀ = 10 lattices) | Shifted lattices at fixed n₀ |
| C28 | From about 145 points the default BARS fits sit at the replicate noise floor, so tighter fits need replicates, not settings | §12.18 | Observed (against an independent sweep) | More replicates on one curve |
| C29 | Replicate weights through the robust-weight seam over-shrink, because the unweighted g-prior borrows the weighted σ² | §12.18 | Observed harm; mechanism a hypothesis | A weighted entry point with its own prior scale |

---

## 8. Lean-worthy pieces

These are the elementary, load-bearing claims the rigor harness exists for, as candidates for protolemmata under
`lean/Protolemmata/`.

- **L1, the finite exponential-family certifier.** For a finite state space and P_β ∝ e^(−βE):
  - ∂_β⟨E⟩ = −Var(E);
  - ∂_β⟨A⟩ = −Cov(A, E), with |Cov| ≤ σ_A σ_E;
  - C8 to C10, by Cauchy–Schwarz.

  Finite sums and calculus in Mathlib; the argon notes call this the load-bearing theorem, and it is.
- **L2, the two-level switch.** C1 to C3. Closed forms on one explicit function, including the traversal length
  ∫dp/√(p(1 − p)) = 2(arcsin√p₂ − arcsin√p₁).
- **L3, the dense-count bound.** Cauchy estimates feeding a cubic spline error bound: C6 with explicit constants.
  Heavier, but Mathlib has Cauchy estimates.

---

## 9. Experiments that need no new engine work

- **E1. The narrowest bump.** In the existing dense χ(T) sweeps (the Domany parity battery and the BWD toy), measure
  the narrowest distinct bump in ln T and compare it with w_q ≈ 0.75 at q = 20 (C4).
- **E2. D_J against known switches.** On the same sweeps, compute D_J for every adjacent pair of a coarse dyadic
  subsample, and check that every interval hiding a known split has D_J ≥ ℓ(θ)² (C11). Also record how large D_J gets on
  featureless intervals, the certificate's false-alarm side.
- **E3. c_BARS.** Measure how many points per feature width BARS needs to register a bump: synthetic bumps of known
  width at known spacing and replicate counts, in the style of `BirthBenchmark`. This calibrates C5's sampling constant
  against BARS rather than geometry.
- **E4. The miss rate.** On BWD datasets with known feature inventories, measure the volley at which each feature is
  first registered, against its width. That turns the unsampled-feature risk into an operating curve, and tests C12.

---

## 10. Reading list

These are cited from memory; confirm the details before relying on them.
- **Partition-function zeros:** M. E. Fisher (1965), *The nature of critical points*; T. D. Lee & C. N. Yang (1952),
  Phys. Rev. 87.
- **Finite-size scaling of zeros:** C. Itzykson, R. B. Pearson & J.-B. Zuber (1983), Nucl. Phys. B 220; W. Janke &
  R. Kenna (2001), J. Stat. Phys. 102.
- **Thermodynamic length:** P. Salamon & R. S. Berry (1983), Phys. Rev. Lett. 51; G. E. Crooks (2007), Phys. Rev.
  Lett. 99.
- **Histogram reweighting:** A. M. Ferrenberg & R. H. Swendsen (1988, 1989), Phys. Rev. Lett. 61 and 63.
- **The sampler and the method:** R. H. Swendsen & J.-S. Wang (1987), Phys. Rev. Lett. 58; M. Blatt, S. Wiseman &
  E. Domany (1996), Phys. Rev. Lett. 76, and (1997), Neural Computation 9.
- **Conformal Bayes:** E. Fong & C. Holmes (2021), *Conformal Bayesian computation*, NeurIPS.
- **Anytime-valid inference:** A. Ramdas, P. Grünwald, V. Vovk & G. Shafer (2023), *Game-theoretic statistics and safe
  anytime-valid inference*, Statistical Science 38.

---

## 11. Open questions for the exploring thread

1. **Coupled clusters.** Does C11 survive them? The question is whether the system's energy variance is at least each
   cluster's, or what replaces that inequality.
2. **What should the certificate's observable be?** Energy (which D_J uses), or a normalized cluster observable, for
   which 5.3's excursion bound is tighter?
3. **A priori or a posteriori?** Should the a priori floor (Thread I) set the minimum number of volleys before any stop,
   with 5.3's certificate checked per interval after each volley? Or does one make the other redundant?
4. **Volley-to-volley evidence.** Is there an anytime-valid formulation (§6) with useful power, given the replicate
   noise?
5. **The bracket.** How do the window E and its padding interact with the expansion-only bracket loop? And does
   normalizing ln T change that loop's identifiability limit?
6. **What should regime B inherit** from regime A's final grid: the D_J map, the σ_E map, the regions' posterior widths?
   And should regime-B points ever feed back into the curve, which would need the per-point weights to pool
   observations of different provenance?

---

## 12. Findings, 2026-10-24 (exploration thread)

A first pass on §5.3. Nothing here is a decision; §12.5 lists what it changes in §11, for discussion.

**How it was checked.**
- **Exact enumeration on a 9-spin toy:** two 4-cliques (J = 1) joined by a weak cut (four edges of 0.15), plus a loose
  point (one edge of 0.03, so it detaches about 3 e-folds colder than the split). Spins are enumerated at q = 3, the
  random-cluster measure at q = 20. The script is a scratch file, not committed.
- **The Sep-16 parity profile dump,** `artifacts/parity-profiles.tsv`: the toy (N = 3200) and iris (k = 10), 48 log steps
  over [0.05, 4]·T_ps, so Δ ln T = 0.093. It carries the FK channels but not MeanEnergy.

### 12.1 C11 holds for coupled clusters (§11 Q1)

- **Where the assumption came from.** 5.3 compared the system's thermodynamic length with the switching cluster's own,
  which needs the cluster's energy uncorrelated with the rest. Comparing through an observable removes it.
- **The bound.** Take any observable A with values in [0, 1] and mean μ(β). Bhatia–Davis gives Var(A) ≤ μ(1 − μ), so
  |μ′| = |Cov(A, E)| ≤ √(μ(1 − μ))·σ_E. Integrating:
  2|arcsin√μ(β₂) − arcsin√μ(β₁)| ≤ L(β₁, β₂) ≤ √D_J.
  - The left side is the Fisher–Rao angle between Bernoulli(μ₁) and Bernoulli(μ₂).
  - So this is data processing for Fisher information: no coarse-graining of the Gibbs path is longer than the path.
- **The corollary, now unconditional.** For a split, take A = 1{A's and B's majority labels differ}, a function of the
  whole coupled system. Its traversal forces L ≥ ℓ whatever the coupling. **Checked here:** on the coupled cliques,
  three traversals give arcsine length ≤ L ≤ √D_J (1.03 ≤ 2.33 ≤ 3.01; 1.32 ≤ 3.26 ≤ 4.66; 0.43 ≤ 0.90 ≤ 0.97).
- **A correction to the plateaus.** For a split, that indicator's hot plateau is (q − 1)/q, not 1, so θ > 1/q.
  - At q = 20 the "5% plateaus" case is the T → ∞ limit, so ℓ(0.05)² ≈ 5.0 is never reached at finite T.
  - 5.3's placement of the argon gate should be read with θ ∈ (1/q, ½).

### 12.2 The sharp excursion bound (C10 tightened)

- **The bound.** From the arcsine form, for A with values in [a, b]:
  sup over β ∈ [β₁, β₂] of |⟨A⟩_β − ⟨A⟩_{β₁}| ≤ (b − a)·sin(min(L, π)/2) ≤ (b − a)·sin(min(√D_J, π)/2).
- **It is sharp.** The two-level switch traces a Fisher–Rao geodesic and attains the bound on its own indicator:
  sin(ℓ(θ)/2) = 1 − 2θ exactly.
- **It never exceeds the trivial bound** b − a, unlike C10's linear form.
- **Checked here:** 40 random observables × 300 random intervals, with worst ratios 0.75 (angle against L) and 0.61
  (excursion against the sine bound). The two-level length over β ∈ [0, 20] at g = 19 is 2.6902, matching C3's
  corrected closed form.

### 12.3 What D_J cannot do: it is blind to size

- **Every two-level event costs the same length** (C3), whatever the sizes on either side. A significant split and a
  loose point detaching both cost about 2.7 at q = 20.
  - Toy, in one-e-fold windows around each event, each about 2/3 traversed: D_J = 2.31 for the 4|4 split and 1.90 for
    the loose point.
- **So the false-alarm floor counts events, not mass.**
  - An event at T* releases a latent heat of about (1 − 1/q)·T*·ln(q − 1). An interval traversing events at density ρ
    per e-fold then has D_J ≈ (1 − 1/q)·ln(q − 1)·ρ·(Δ ln T)², about 2.8ρ(Δ ln T)² at q = 20.
  - On a graph with many loose points ρ is of order N/E, so a uniform gate D_J ≤ c² needs spacing of order c/√(2.8ρ).
    That is 5.1's √N count, reached from the event side.
  - **Derived, unverified at scale;** E2's false-alarm side measures it.
- **The θ-guarantee is vacuous for rigid splits after the first volleys.** At q = 20 a 10–90% traversal spans 1.93
  e-folds in ln T; the 25–75% width is 0.78, matching C1's linearized 0.75. Once h < 1.93, no rigid split can hide at
  θ = 0.1. The guarantee keeps its force for melts and k-way splits, which are narrower.
- **And D_J ≤ 4 passes intervals carrying most of a split.** The toy's split window passes at 2.31 with 2/3 of the
  traversal inside. That isn't a miss, since the split shows at the endpoints, but D_J can't tell it from background.
  Something else has to read it.

### 12.4 A monotone, size-weighted certificate

The seed is the argon notes' closing remark: the unsubtracted FK mean is monotone, and "acts as an anti-aliasing
certifier". The monotonicity chain itself is already written up in
[fable-potts-griffiths-lemmas.md](../../lean/Protolemmata/fable-potts-griffiths-lemmas.md) (Grimmett Thm 3.21 with the
Edwards–Sokal two-point identity).

- **The fact.** For ferromagnetic Potts (every J_e ≥ 0) with q ≥ 1 on a finite graph, every pair connectivity
  τ_ij(β) = φ_β(i ↔ j) is non-decreasing in β. That is the random-cluster comparison inequality, with the inhomogeneous
  p_e = 1 − e^(−βJ_e). The mean of every increasing bond functional is monotone the same way.
- **The observable is already recorded.** S₂ = ⟨Σ_k c_k²⟩/N² = (1/N²)·Σ_ij τ_ij is SusceptibilityFkCluster/N, the
  default χ (`SusceptibilityKind.FkCluster`).
- **The certificate.** Every τ_ij moves one way, so nothing inside an interval can exceed its endpoints. For disjoint
  vertex sets A and B, the pair-averaged connectivity change at any β in [β₁, β₂] satisfies
  Δτ̄_AB ≤ N²·ΔS₂/(2|A||B|). It needs endpoint means only: no Cauchy–Schwarz, no rigid clusters, no coupling
  assumption.
- **No event hides; an event can only be unlocalized.** An interval either carries a significant ΔS₂ or certifiably
  contains no significant merge. For this observable, §2's unsampled-feature problem becomes a localization problem.
- **The weighting is by size.** A merge of clusters of sizes a and b adds 2ab/N². Toy: ΔS₂ = 0.268 for the split against
  0.131 for the loose point (D_J: 2.31 against 1.90).
- **It weights by a product, not a minimum.** One point leaving a cluster of size a adds 2a/N², which beats a merge of
  two n-clusters once a > n².
  - **A size-floor-aware version.** Any superadditive g makes Σ_k g(c_k) increasing in the bonds. With g(c) = (c − n)₊,
    a merge of two clusters of size ≥ n adds exactly n, and one point joining adds 1.
  - So F_n = ⟨Σ_k (c_k − n)₊⟩/n counts significant merges, with shedding weighted 1/n. Toy (n = 3): 0.68 for the split
    against 0.22.
  - The SW kernel doesn't accumulate F_n today. It would be one more accumulator: engine work, parked with OQ7 and OQ8.
- **Checked here.**
  - Exact random-cluster enumeration at q = 20 over 601 temperatures: S₂, F₃, τ_AB and ⟨E⟩ are monotone throughout.
  - Real SW data, the parity dump: S₂ has no violations on the toy and two on iris, the larger 0.0025 (noise level).
- **A preview of the brackets:** ΔS₂ per interval, cold to hot, on stride subsamples of the 48-point grid (stride 16 is
  Δ ln T ≈ 1.49, stride 4 is 0.37):
  - toy: 0.460 0.180 0.000 at stride 16; 0.371 0.021 0.017 0.051 0.104 0.072 0.004 0 … at stride 4;
  - iris: 0.151 0.334 0.001 at stride 16; 0.018 0.020 0.034 0.078 0.099 0.195 0.035 0.005 0.001 0 … at stride 4.

### 12.5 What it changes in §11 (for discussion; nothing decided)

- **Q1** is answered: C11 holds for coupled clusters (12.1).
- **Q2:** the observable that suits a certificate is a monotone cluster functional, S₂ today or F_n with one
  accumulator. Energy is its size-blind cousin.
- **Q3:** with a monotone certificate an unflagged interval certifies itself, so a volley floor isn't needed for safety.
  Thread I's count becomes a budget prior. The stop would then certify that every flagged interval is narrower than the
  regime-B window, or resolved by the BARS readout.
- **Q4:** each interval's decision is a one-sided test on ΔS₂ against a floor such as 2n²/N², with replicate error bars.
  ΔS₂ telescopes over the dyadic tree: a parent's increment is the sum of its children's. How those tests compose across
  volleys is still open.
- **Q5:** the bracket gets an endpoint coverage check.
  - S₂(0) − S₂(T_lo), with S₂(0) = Σ n_comp²/N² over the graph's connected components, bounds the merge mass colder than
    the bracket. S₂(T_hi) − 1/N bounds what is hotter.
  - In the dump, both datasets sit at 1/N on the hot side. On the cold side the toy reads 0.640 and iris 0.493, and the
    toy's S₂ is still falling steeply at T_lo (0.17 in its first grid step).
  - Whether that is the graph's component floor or merges colder than the bracket is a direct check against the graph.

### 12.6 Experiments and Lean, refined

- **E1** runs on the dump as it is. Δ ln T = 0.093 is about 8 points per w_q, enough to measure the FWHM in ln T of the
  FkReduced and MagnetizationVariance bumps.
- **E2** needs the MeanEnergy channel, which every SweepProfile carries and the dump omits. Stride subsamples (16, 8, 4,
  2) of the 48-point grid are nested, so they stand in for volleys.
- **E5, new: the monotone certificate against known inventories** (toy: three stripes and residue; iris: setosa, then
  versicolor/virginica at stage 2). Score each stride's flagged intervals against the inventory, with ΔS₂'s noise floor
  taken from replicate spread.
- **L1 gains one lemma:** the Fisher-angle bound of 12.1 and 12.2 (Bhatia–Davis, then the arcsine substitution). It is
  elementary.
- **L4, new: monotone cluster functionals.**
  - A superadditive g makes Σ_k g(c_k) increasing under adding bonds: finite combinatorics.
  - Its mean is then monotone in β, either by the comparison inequality or by
    d/dβ φ(f) = Σ_e p_e′·Cov(f, ω_e)/(p_e(1 − p_e)) ≥ 0 under positive association.
  - It shares its one axiom with the potts-griffiths chain.

### 12.7 Continued: the toy's bracket, sufficiency, and weighting

- **The toy dump is suspect** (C22). The user saw the toy's cold edge still rising across several observables, and
  suspects the bracket, possibly through the graph compiler.
  - S₂ agrees that the edge is clipped. The toy's FK structure sits in the cold 1.5 e-folds of a 4.4 e-fold window:
    S₂ ≈ 1/N by T ≈ 0.1, both FkReduced peaks lie within 0.75 e-folds of T_lo, and S₂ falls 0.17 in the first grid step.
    Iris, by comparison, has its events mid-window.
  - Until the bracket is checked, the toy numbers in 12.4 and 12.5 are illustrations, not evidence.
- **Two explanations, and they separate.**
  - Either the couplings sit on a different scale from the one the q-only bracket assumes, which shifts every event
    cold, or the residue's weak edges let single points detach across the cold edge.
  - A clipped split is one large step in S₂; shedding is many small ones.
  - The per-edge cut temperatures T*_e = J_e/ln(q − 1) on the toy graph, read against [T_lo, T_hi], separate them
    without SW. So does the size histogram at T_lo: mostly singletons, or a second large cluster.
- **S₂ witnesses FkReduced** (C20).
  - Under a monotone coupling across an interval (Strassen, for the finite chain β₁ ≤ β ≤ β₂), clusters only merge.
  - The largest cluster's square then grows by at most 1.5 times Σc²'s growth. The worst case is the largest piece
    merging with an equal one.
  - So Σc² − c_max² moves within [−½, 1] times Σc²'s growth, and inside the interval sup |Δχ_red| ≤ N·ΔS₂.
  - FkCluster is S₂ itself. **Checked here:** enumeration on six random weighted graphs, worst ratio 0.44 against the
    bound 1.
- **The magnetization curves aren't witnessed this way.** N_max follows the colourings as well as the clusters. What it
  shows beyond S₂ is colour fluctuation, including the documented cold-edge artefact of weakly coupled spins: fluctuation
  without merging, which S₂ leaves unflagged.
- **The weighting is a choice of currency, and the certificate holds for each.** Every superadditive g gives a monotone
  ⟨Σ_k g(c_k)⟩ (C18):
  - g = c² counts pairs: ΔS₂ is the fraction of point pairs disconnected, a Rand-index currency dominated by big
    clusters;
  - g = c ln c is the partition entropy, so the expected FK partition entropy is monotone in T (checked by enumeration).
    It weights a merge by about the smaller side times 1 + ln(larger/smaller);
  - g = (c − n)₊ is a thresholded count.
  - On the toy (N = 3200), product weighting makes one residue point leaving a 2,500-point cluster outweigh a merge of two
    30-point clusters, by about 3 to 1. The concern about product weighting is right for residue-heavy data.
- **No monotone functional can ignore shedding** (C21).
  - Telescoping g(2n) − g(n) over single-point joins gives merge(n, n) ≤ n·max over a ∈ [n, 2n) of join(a). So a
    significant merge weighs at most about n shed points, whatever g is.
  - On residue-heavy data, every certificate's per-interval total is swamped by shedding once an interval sheds more than
    about n points.
  - What tells a split from shedding is concentration: a split is a step about w_q wide, while shedding spreads over the
    residue's cut temperatures. That is a question of shape across intervals, which is BARS's job, not of per-interval
    totals.
- **What is recorded today.** S₂ and ⟨c_max²⟩ exactly, through the FkCluster and FkReduced channels. The SW kernel's size
  histogram is log-binned (cluster counts in bins [2^k, 2^(k+1))), so it bounds the other functionals without giving
  them exactly.

### 12.8 The toy's cold edge, diagnosed

A scratch probe (not committed) rebuilt the parity battery's graphs: k = 10, mutual KNN, MST repair, a Gaussian kernel
at the mean edge distance, no LMP rescaling. It read the per-edge cut temperatures T*_e = J_e/ln(q − 1) against the
bracket [0.0129, 1.031].

- **The bracket isn't miscomputed; its recipe can't see this graph's weakest significant cut.** T_lo = 0.05·T_ps is what
  the q-only formula gives.
- **The stripes are joined only through the sparse background in the two gaps between them.** The graph has no
  stripe–stripe edges.
  - The background's bonds are weak because the one global bandwidth (0.162, the mean edge distance) is set by the dense
    stripes, and background edges are two to four times longer.
  - Background–background edges have a median T* of 0.005, an e-fold below T_lo, and 66% of them lie below it. Stripe
    edges have a median of 0.030, with 1% below.
- **So the bracket cuts through the toy's defining event.**
  - Single linkage puts the two stripe joins 0.2 and 0.3 e-folds below T_lo, each through a background–background
    bottleneck.
  - SW's collective binding lifts the split just inside the bracket: S₂ falls from 0.64 at T_lo to the three-stripe
    plateau of 0.25 by T = 0.0155, two grid steps up.
  - The split's cold side, with the stripes joined, lies colder than the bracket. That is the edge seen still rising.
- **Below that, the background keeps ordering for several e-folds** (5% of its bonds sit 6.6 e-folds below T_lo). It is
  a shedding tail, and matches the parity test's own note that Var(m)'s maximum sits on the sparse background's cold
  ordering transition.
- **Not components.** Before repair the toy has 13 components with a 3,179-point giant, an S₂ floor of 0.987. Its 12 MST
  repair edges don't matter.
- **Iris is clean, and shows a trap for the coverage check.**
  - Setosa is joined to the rest by a single MST repair edge with T* ≈ 2×10⁻⁶, 8.8 e-folds below T_lo.
  - Iris's S₂(T_lo) = 0.493 sits against a pre-repair floor of 0.504, so its bracket covers its real structure.
  - A coverage check (C19) has to measure against the pre-repair components, or it chases the repair edge.
- **A graph-only lower bound on S₂** (C23).
  - Any subgraph's connectivities are at most the full graph's: dropping an edge sets its p_e to 0, and the comparison
    inequality applies.
  - A spanning tree's random-cluster measure is a product measure with p̃_e = p_e/(p_e + q(1 − p_e)). So S₂_tree(T) is
    closed-form, O(N) per temperature, and S₂(T) ≥ S₂_tree(T) at every T.
  - **Checked here:** the maximum spanning tree's S₂ lies below the measured S₂ at all 48 grid points, on both datasets.
  - It is loose where binding is collective (the toy at T_lo: 0.14 against 0.64), since a tree breaks dense clusters at
    single bonds.
  - It still bounds the cold side from the graph alone. With the stripes apart, S₂ can't exceed 0.375, so
    P(stripes apart) ≤ (1 − S₂)/0.625. The toy's tree S₂ is 0.78 one e-fold below T_lo and 0.89 two e-folds below, so
    P(stripes apart) ≤ 0.36 and ≤ 0.18 there; the real join is warmer.
- **What this means for the experiments.** The toy's interior is sound, from the three-stripe plateau through the melt.
  E1–E5 on the toy need the cold end extended by about an e-fold (11 more steps at the current spacing) to hold the
  split's cold side. The background's tail is shedding, and can stay outside.
- **The graph-compiler link.** The background's couplings follow from one global bandwidth on data with two densities,
  under mutual KNN. LMP rescaling, off in the parity config, is the codebase's existing fix for density inhomogeneity.
  Whether parity wants it is a question for the graph-compiler work, not for this thread.

### 12.9 E1, E2 and E5 on an extended grid

**The run.**
- **The grid.** A scratch harness (not committed) reran the toy and iris sweeps on the parity grid extended 11 steps
  colder and 6 hotter at the same spacing: 65 temperatures from 0.0046 to 1.80, 5.97 e-folds, Δ ln T = 0.093. With 64
  intervals, the stride-2^(6−v) subsamples are exact dyadic volleys v = 1 to 6.
- **The budgets.** The main budget is 8 replicas × (1,000 burn-in, 2,000 cycles), 57 s for the toy. A cheap rerun at the
  parity test's budget, 4 replicas × (200, 400), gives the noise comparison.

**Two checks first.**
- **Monotonicity (C17) holds on real SW data.** At the main budget, S₂ falls and ⟨E⟩ rises at all 64 steps on both
  datasets, with no reversal beyond 2 SE. At the cheap budget, one reversal on iris exceeds 2 SE.
- **The colder bracket holds the stripe split** (C22).
  - The toy's S₂ is 0.857 at the new T_lo, against 0.564 at the old. The size histogram there shows one giant plus
    about 36 singletons, 21 pairs and a few small clusters per draw: the stripes are joined, and only the background's
    shedding tail lies outside.
  - Iris reaches 0.549, near its floor once setosa's repair edge (T* ≈ 2×10⁻⁶) is excluded: (50² + 100²)/150² = 0.556.
  - Iris's old T_lo missed about 0.06 of S₂ mass, consistent with small components joining through repair edges at
    T* = 0.006–0.011.

**E1, the narrowest bump** (C4).
- **The yardstick.** A rigid split's variance bump p(1 − p) has a full width at half maximum (FWHM) of 1.38 in ln T at
  q = 20. C4's 0.75 is the switch's 25–75% width, not a bump width.
- **The measurements:**
  - **FkReduced:** one bump on the toy, FWHM 1.57 at T = 0.025, which is the three stripes' lifetime. It respects the
    floor. Iris shows no distinct FkReduced bump.
  - **MagnetizationVariance:** 0.41 on the toy at T = 0.014, the stripe split itself; 0.24 on iris at T = 0.091.
  - **Heat capacity:** 0.31 on the toy at the stripes' melt (T = 0.11); 0.23 and 0.15 on iris (T = 0.10 and 0.16).
- **So C4's floor fails for variance-type observables,** by 3 to 9 times. Two mechanisms are at work, both outside the
  rigid giant-spin model:
  - **Melts are sharp** (ΔS ≈ (n_c − 1) ln q), as §3 says. But a melt is a step only in S₂. In the variance
    observables it is a narrow spike, and a spike can fall between grid points.
  - **The toy's split runs through the background in the gaps** (C24). It frees the bridging spins as well as cutting
    bonds, so its entropy jump exceeds ln(q − 1) and the split narrows towards a melt. The width implies an entropy jump
    of roughly 9 nats, against 2.9 for a bare split.
- **Bump widths need the main budget.** At the cheap budget noise makes spurious narrow bumps: down to 0.06 in the toy's
  heat capacity, and six bumps of 0.09–0.36 in iris's MagnetizationVariance.
- **For Thread I:** at the a priori count's spacing (17 points over 4.6 e-folds, h ≈ 0.29), iris's melt spikes
  (0.15–0.24) can fall between points in the variance observables. S₂ still shows each one as a step.

**E2, D_J** (C16).
- **The argon gate fails almost everywhere on the toy.** Even at v = 6 (h = 0.093), 28 of 64 intervals have D_J > 4,
  with a maximum of 155. Holding D_J ≤ 4 everywhere would need a spacing near 0.015 e-folds, about 400 points.
- **D_J is extensive.** At v = 6 its maximum is 155 on the toy and 6.5 on iris: a ratio of 24, against an N ratio of 21.
- **Its false alarms are where the theory puts them.** Of the toy's 28 intervals above 4 at v = 6, 9 have ΔS₂ below the
  floor:
  - on the three-stripe plateau (T = 0.023–0.033), D_J ≈ 6 per interval as points shed from the stripes. In
    2.8ρ(Δ ln T)² terms, that is about 250 events per e-fold;
  - after the melt (T = 0.12–0.19), D_J falls from 127 to 6 while S₂ is already near 1/N. That is short-range order
    dissolving, with no cluster mass at all.
- **The quadrature ratio R** (the argon notes' closure check) runs from 0.24 to 37 at the coarse volleys and sits near 1
  only on featureless intervals. Its departures follow σ_E²'s curvature at the events, so it doesn't isolate hidden
  spikes here.

**E5, the monotone certificate** (C17, C19).
- **At the main budget it recovers the toy's inventory.** v = 6 flags three runs:
  - the background's shedding tail at the cold edge, [0.0046, 0.0081] (ΔS₂ = 0.043);
  - the stripe split, [0.0089, 0.0225] (0.545);
  - the stripes' melts, [0.039, 0.121] (0.227).
  The quiet stretch between them, [0.0225, 0.039], is the three-stripe window that BWD cluster in.
- **Iris** flags [0.0272, 0.175] (0.444), where the versicolor/virginica separation and the melts overlap in ln T, plus
  [0.0081, 0.0089] (0.008), one small component's join through a repair edge.
- **Localization across volleys.** The toy's flagged share of intervals goes 2/2, 3/4, 5/8, 9/16, 18/32, 28/64. It
  converges to the share of ln T the events occupy (about 2.6 of 5.97 e-folds), not to zero: a broad event stays
  flagged at every spacing.
- **Noise sets the floor's cost.** At the main budget the median SE(ΔS₂) is 0.0003–0.001 on the toy and 0.0006–0.0015
  on iris, against floors of 0.005 and 0.0057. At the cheap budget iris's SE reaches 0.008 at v = 1, above its floor.
- **The per-interval test loses power at fine volleys** (C25).
  - An event spread over several intervals contributes ΔS₂ ∝ h to each, while SE stays fixed. At the cheap budget the
    toy's v = 6 flags fragment into 17 intervals in four runs, against 28 in three at the main budget.
  - ΔS₂ telescopes over the dyadic tree, so a test that descends from flagged parents keeps the parents' power. That
    bears on Q4.

### 12.10 E4: when BARS first registers each bump

**The run.**
- **The fits.** BARS (cubic, Poisson(5) knot prior, 4 chains × 4,000 draws after 2,000 burn-in, the default 10%
  prominence gate) ran on each volley's subsample of the main-budget curves of §12.9: v = 2 to 6, 5 to 65 points.
  v = 1 has 3 points, too few for a cubic.
- **Equal weights.** The run API takes no per-point weights, so every point weighs the same (W ≡ I).
- **v = 2 can't be fitted.** At 5 points even the zero-knot cubic fails the support gate, so BARS refuses to start
  (since 6bc5134 it says so; before, it returned NaN curves). The first usable fit is v = 3 (9 points).
- **The fits below are on 6bc5134.** Before it, a chain whose overdispersed start held too many knots could stay there,
  and that touched the 9-point fits. Rerun on the fix, every fit from v = 4 up is bit-identical and only v = 3 moved.
- **Registration.** A bump counts as registered when a peak mode with posterior mass ≥ 0.5 falls inside the dense bump's
  FWHM. Every registration here reached mass ≥ 0.9 at once.

| Bump | FWHM | First volley | h there | Points per FWHM | Location error | Interval width |
|---|---|---|---|---|---|---|
| toy FkReduced, T = 0.025 (the three-stripe lifetime) | 1.57 | 3 | 0.75 | 2.1 | +0.17 | 0.72 |
| toy MagnetizationVariance, T = 0.014 (the split) | 0.41 | 4 | 0.37 | 1.1 | +0.13 | 0.60 |
| toy heat capacity, T = 0.110 (the melt) | 0.31 | 4 | 0.37 | 0.8 | +0.16 | 0.27 |
| iris MagnetizationVariance, T = 0.091 | 0.24 | 5 | 0.19 | 1.3 | −0.06 | 0.36 |
| iris heat capacity, T = 0.100 | 0.23 | 4 | 0.37 | 0.6 | +0.03 | 0.24 |
| iris heat capacity, T = 0.160 | 0.15 | never | — | 1.6 at v = 6 | — | — |

Widths, errors and interval widths are in e-folds of ln T.

- **BARS registers a bump at about one point per FWHM** (0.6–2.1). That is below Thread I's geometric two points per
  FWHM (C5), so on these data the geometric constant is conservative. The closed-form readout does what §1 says: it
  places a peak between grid points.
- **At first registration the location is rough:** errors of 0.03–0.17 e-folds, with intervals 0.24–0.72 wide. By v = 5
  the intervals run 0.09–0.36 wide.
- **The location intervals cover the dense peak** (to within its own grid quantization, ±h/2) in 12 of 13 fits from
  first registration on. The exception is iris's heat capacity at v = 5, which misses by about one grid step.
- **An under-resolved volley spreads its mass rather than misplacing it.** At v = 3 (w/h = 0.55) the toy's
  MagnetizationVariance splits its peak mass over three modes between T = 0.017 and 0.026, none above 0.43. A single
  mode of mass 0.57 half an e-fold off the split, seen before 6bc5134, came from a chain stuck at its start.
- **The one miss comes from the prominence gate, not the grid.**
  - Iris's bump at T = 0.160 rises about 5% of the curve's range above its higher neighbour, under BARS's 10% gate, and
    is one point wide even at v = 6.
  - It is plausibly setosa's melt: setosa's edges run hotter (median T* 0.037, against 0.026).
  - S₂'s flagged run at v = 6, [0.027, 0.175], covers it. The certificate doesn't let it hide; the BARS readout just
    doesn't count it. That is the division of labour of §12.7.
- **No false bumps.** No mode with mass ≥ 0.5 falls outside the dense features. Iris's FkReduced, which has no dense
  bump, never gets a mode above 0.2.
- **Not yet tested:**
  - C12's two-quiet-volleys rule, which needs each fit's predictive band at the next volley's midpoints;
  - an operating curve proper, which needs more than six bumps (E3's synthetic bumps, or more datasets);
  - the replicate-variance weights of §1, which need a weighted entry point.

### 12.11 C12: midpoint surprises and the two-quiet-volleys rule

**The run.**
- **The surprises.** Each volley's fit from §12.10 was read on the full grid. At the next volley's new midpoints,
  z = (ȳ − f̂)/√(sd_f² + se²), with ȳ the replicate mean, f̂ and sd_f the fit's posterior mean and sd, and se the
  midpoint's replicate standard error. A second form uses BARS's own noise variance σ̂² in place of se².
- **The rule.** A volley is quiet when max |z| over its midpoints stays under the Bonferroni threshold
  Φ⁻¹(1 − 0.05/2m), as the argon notes have it (2.73, 2.96 and 3.16 for 8, 16 and 32 midpoints). The rule stops at the
  second quiet volley in a row.
- **The first usable fit is v = 3,** since BARS can't start at 5 points. So surprises start at volley 4, and the rule
  can fire at volley 5 at the earliest.
- **Rerun on 6bc5134.** Only the v = 3 fits changed. Volley 4's surprises moved by at most 0.01, and every outcome below
  is the same.

| Curve | Bumps registered (E4) | max \|z\| at volley 4 / 5 / 6 | Rule fires (replicate z) | Rule fires (model z) |
|---|---|---|---|---|
| toy FkReduced | v = 3 | 1.2 / 2.1 / 2.1 | 5 | 5 |
| toy MagnetizationVariance | v = 4 | **13.7** / 1.0 / 1.4 | 6 | 6 |
| toy heat capacity | v = 4 | **6.5** / 2.9 / 1.9 | 6 | never by 6 |
| iris FkReduced | none | **2.9** / 1.8 / 2.4 | 6 | 5 |
| iris MagnetizationVariance | v = 5 | 1.5 / 2.4 / 2.1 | 5 | 5 |
| iris heat capacity | v = 4; the T = 0.16 bump never | 1.8 / **5.9** / 2.3 | never by 6 | never by 6 |

Bold marks a volley that wasn't quiet.

- **No stop came before the last registration** on these six curves.
- **The surprises find features when a midpoint lands on them.**
  - The toy's MagnetizationVariance volley 4 peaks at T = 0.0141, the split's own peak.
  - Its heat capacity peaks at T = 0.133, on the melt's shoulder.
  - Iris's heat capacity volley 5 peaks at T = 0.110, where the v = 4 fit had the melt registered but its shape
    under-resolved.
- **But a quiet volley can mean a wide band, not a resolved curve.**
  - Iris's volley 4 is quiet on MagnetizationVariance and heat capacity, though neither bump had registered in the
    v = 3 fit.
  - Iris's MagnetizationVariance is quiet at volleys 4 and 5, so the rule fires at 5: the very volley its bump first
    registers, with no surprise ever flagging it. One volley less of resolution and the stop would have been premature.
  - This is the power question §6 left open, now observed: the holdout's power falls where the posterior band is wide.
- **The minor melt doesn't move the rule either way.** Iris's heat capacity never fires because of the main melt's
  shoulder at volley 5. The unregistered T = 0.16 bump draws no surprise of its own.
- **The two noise forms mostly agree.** BARS's homoscedastic σ̂² is smaller than the replicate variance at the heat
  capacity peaks, which run noisiest. So the model form is louder there (8.3 against 5.9 on iris) and quieter on flat
  stretches.
- **Coverage.** At quiet volleys, 88–100% of midpoints fall within ±1.96, against a nominal 95%, on 8 to 32 midpoints.
- **Where that leaves C12.** The rule is a reasonable heuristic here, not a guarantee, and its failure mode is the one
  §6 predicted: optional stopping over low-power tests. The S₂ certificate doesn't share it.
  - The flagged runs of §12.9 hold every event bump of both datasets, the unregistered one included, whatever the
    surprises say.
  - The toy's FkReduced bump is a lifetime, not an event. Its peak sits on the quiet three-stripe plateau, between the
    flagged split and melts that bound it, so the flags hold 2 of the toy's 3 bump peaks (and all 3 of iris's).

The experiments are charted in a dashboard, https://claude.ai/artifact/4a51qGP11a5JH2Qr5V4g3n (private; built from the
scratch data, not from the repo).

### 12.12 Discussion: a minimum, then one afterburn volley

The user's working practice, recorded as a heuristic rather than a rule: set a minimum number of volleys, then run one
more so the fit can mature.
- **What E4 says about it.**
  - The last first registration comes at volley 4 on the toy and volley 5 on iris, so an afterburn stops at 5 and 6.
  - One afterburn volley narrows the location intervals 2–7×: toy MagnetizationVariance 0.60 → 0.09 e-folds, heat
    capacity 0.27 → 0.09, FkReduced 0.72 → 0.36; iris MagnetizationVariance 0.36 → 0.12, heat capacity 0.24 → 0.12.
- **Why it fits the geometry.** BARS registers a bump at about one point per FWHM (C26), so one more volley gives about
  two, which is Thread I's geometric constant, h ≤ w/2. The heuristic amounts to that rule, once the minimum is set from
  the narrowest bump that matters.
- **Its caution.** Tighter isn't always right. Iris's heat-capacity afterburn fit (v = 5) is the one location interval
  that missed the dense peak, by about one grid step.
- **What stays open.**
  - Setting the minimum before sampling: E1 puts the narrowest bumps at melts and background-mediated splits (0.15–0.41
    e-folds), not at the rigid split's 1.38, so a priori the minimum needs melt widths.
  - How much of the afterburn's gain the BARS settings could deliver a volley earlier.

### 12.13 What the stopping work asks of SW retention

The SW retention brief (`issues/spc/sw-retention-brief.md`, items S1–S6) predates this thread. Read against it:
- **S6, replica spread in `SweepCurves`, comes first.** The certificate's per-interval test runs on SE(ΔS₂) from the
  replicate spread. The harness had to compute that from the raw runs.
- **S1 with S3, the per-sweep records,** would serve the certificate twice:
  - within-run error bars, so a single-replica sweep could carry the certificate's SE;
  - Ferrenberg–Swendsen reweighting, to predict a midpoint's S₂ or C_v from its endpoints before it is run.
  The second reaches only so far at this size. Reweighting needs overlapping energy histograms, which is the D_J ≲ 4
  regime, and E2 found D_J > 4 on most toy intervals even at h = 0.093.
- **S2, the higher moments,** would put error bars on the heat-capacity peaks of E1 and E4.
- **Not in the brief: an exact cluster-size histogram.**
  - The log-binned histogram bounds Σ_k g(c_k) without giving it.
  - Counts per exact size, summed over sweeps, are a sufficient statistic for the mean of every cluster functional
    Σ_k g(c_k): S₂, F_n at any size floor, the partition entropy.
  - It costs N integers per run (13 KB for the toy), far less than S1, and it lets the size floor and the currency of
    §12.7 be chosen after the run.
  - S5's partitions would give the same means, at gigabytes.
  - Added to the brief as S7 (de7bb49).

### 12.14 The experiments on the user's schedule (n₀ = 10)

**The correction.**
- **What ran before.** §12.9–12.11 ran bisection from S₀ = {0, 1}: n_v = 2^v + 1, so 3, 5, … 65 points.
- **The user's design.** It starts from n₀ points, 10 in practice, so n_v = (n₀ − 1)·2^v + 1: 10, 19, 37, 73, 145, 289.
- **Comparing the two.** Match them by spacing, not by volley index.

**The run.**
- **The grid.** 289 temperatures (9·2⁵ + 1) over the same 5.97 e-folds, so volleys 0–5 are exact subsamples, with h from
  0.663 to 0.021 e-folds.
- **The budget.** 8 replicas × (1,000 burn-in, 2,000 cycles), 251 s for the toy. BARS, at its defaults on 6bc5134,
  fitted every volley.

**The results.**
- **Monotonicity holds at all 288 steps.** No S₂ reversal exceeds 2 SE on either dataset; one ⟨E⟩ reversal does on
  iris.
- **E1 needs a noise-aware bump finder at this density.**
  - Point to point, heat capacity and all of iris's curves move by about one replicate SE.
  - Topographic prominence (at least 10% of the range and 4 SE) recovers the 65-point widths: toy 1.58, 0.43, 0.30;
    iris 0.29, 0.16, 0.09.
  - It adds a one-point bump at the toy's melt (T = 0.1034), at the noise level.
- **E2.** The toy still has 16 intervals over D_J = 4 at v5 (h = 0.021, maximum 7.9); iris clears the gate from v4 on.
- **E5, and C25 at full strength.**
  - The per-interval test flags one run over the active region at v0–v2, then localizes at v3 (toy: 24 of 72 intervals,
    in 5 runs).
  - Then it fragments. At v5 the toy flags 13 of 288 intervals in 8 runs, and the test misses bump peaks it held at
    coarser volleys: the toy's melt at v4–v5, iris's minor melt at v4–v5.
  - A first-cut descending test keeps a child flagged when its parent is flagged and the child's own ΔS₂ exceeds 2 SE,
    without the floor. It holds every bump peak at every volley except iris's main melt at v5, flagging 84 and 90 of 288
    intervals at v5.

**E4, first registration** (mass ≥ 0.5 inside the dense FWHM):

| Bump | FWHM | First volley | Points per FWHM there |
|---|---|---|---|
| toy FkReduced, the three-stripe lifetime | 1.58 | v0 (10 points) | 2.4 |
| toy MagnetizationVariance, the split | 0.43 | v1 (19) | 1.3 |
| toy heat capacity, the melt | 0.30 | v0 (10) | 0.4 |
| toy heat capacity, the one-point bump | 0.02 | never | — |
| iris MagnetizationVariance | 0.29 | v2 (37) | 1.8 |
| iris heat capacity, the main melt | 0.16 | v0 (10) | 0.2 |
| iris heat capacity, the minor melt | 0.09 | v5 (289) | 4.5 |

- **Registration depends on the lattice's phase as well as its density** (C27). The low values (0.2, 0.4 points per
  FWHM) are lucky landings: the 10-point grid happens to put a point on each melt's shoulder. A minimum set from one
  lattice can be optimistic. The worst case across phases is the one to plan on, about 1–2 points per FWHM here.
- **The minimum on this schedule** is volley 1 (19 points) on the toy and volley 2 (37) on iris, or volley 5 to include
  iris's minor melt. An afterburn then stops at 2 and 3 (or 6).
- **Iris's minor melt registers here.** It never registered on the 65-point lattice; here it registers at v5, at 4.5
  points per FWHM.
- **Location intervals turn overconfident at fine volleys.**
  - From v3 on, intervals 0.00–0.04 e-folds wide exclude the dense argmax for the toy's melt and both
    MagnetizationVariance peaks.
  - The dense argmax is itself noisy at this spacing, so this isn't yet a calibration verdict. That needs an independent
    dense sweep, at another seed, as the truth.
- **The default BARS stops mixing on heat capacity.**
  - From 37 points up, the heat-capacity fits reach knot R̂ of 1.08–1.40 and knot ESS of 8–44. On iris their fit R̂
    runs 1.45–2.95, and the knot count climbs to 11.6 as the fit chases noise at the melt peaks.
  - MagnetizationVariance's peak R̂ reaches 1.1–1.2 at 73–145 points.
  - Equal weights on sharp, noisy peaks is the likely cause. Replicate weights and tempering are what a settings sweep
    should try first.
- **C12 on this schedule.**
  - The two-quiet rule stops at v2–v4: toy FkReduced v2, MagnetizationVariance v3, heat capacity v2; iris FkReduced v2,
    MagnetizationVariance v3, heat capacity v4.
  - Every stop comes after the curve's main bump registers. But iris's heat capacity stops at v4, before its minor melt
    registers at v5: the first premature stop.
  - Fine volleys turn loud again after the stops (toy MagnetizationVariance 4.7 at v5; heat capacity 4.3 and 3.9 at
    v4–v5), as the fits mismodel the noisy peaks.

### 12.15 Discussion: which observable for which job

The user's position: MagnetizationVariance is the detection observable, since it peaks where it should and
unambiguously. FkReduced gives a smoother bump, which is worse for detection but perhaps better for extracting features.
Heat capacity is out of scope. Pressure-tested against the 289-point cluster-size histograms:
- **MagnetizationVariance peaks on the giant's breakup.**
  - Toy: at T = 0.0141 the giant (≥ 2,048 points) is present in 22% of draws and the 1,024–2,047 bin in 48%: the
    stripe split, caught in the act.
  - Iris: at T = 0.0858 the versicolor-plus-virginica cluster (64–127 points) is present in half the draws.
- **First caveat: it is giant-dominated.** It weights an event by how far it moves the largest colour group. The toy's
  stripe melts (T ≈ 0.09–0.11) sit at 0.8–4.5% of its peak, so it is blind to the lower levels of a hierarchy. Those
  show as S₂ steps.
- **Second caveat: colour collisions raise its baseline.** Separated clusters share a colour with probability 1/q, so
  N_max still fluctuates on a plateau. On iris it stays at 65% of its peak across the whole setosa-apart plateau
  (T = 0.03–0.06), so the peak stands only about 35% above it.
- **FkReduced measures a lifetime, so its span is the clustering window.** On the toy, its bump (FWHM 1.58) is the three
  stripes' lifetime, bounded by the split and the melts, which suits regime B's windows. But it needs the second
  cluster to be born by a split. Iris's setosa is a separate component from T = 0 (FkReduced ≈ 50²/150 ≈ 16.7 at the
  cold end), so iris's FkReduced only declines and has no bump to read.
- **Neither is covered by C20 the same way.** S₂ bounds FkReduced's excursions (C20), but not MagnetizationVariance's,
  which follows the colourings as well. In practice its peaks sit inside S₂-flagged runs on both datasets.
- **So the jobs split three ways:**
  - MagnetizationVariance locates the giant's transitions;
  - FkReduced gives spans where a second cluster is born by a split;
  - S₂'s flags cover every event, the lower hierarchy levels included.

### 12.16 Why iris reads differently, and how BWD and Domany read it

From the 289-point cluster-size histograms (cluster counts per draw by size bin):
- **The stage-1 plateau starts at about zero temperature.**
  - From T = 0.0046 to about 0.05, every draw holds exactly one 32–63-point cluster (setosa) and one 64–127-point
    cluster (versicolor plus virginica).
  - Setosa splits off through the single MST repair edge at T* ≈ 2×10⁻⁶ (§12.8). So within any sensible bracket there is
    no ordered phase and no cold peak: the superparamagnetic plateau is already there.
  - FkReduced therefore reads as a plateau (setosa²/N ≈ 16.7), then a cliff, with no rise. MagnetizationVariance sits at
    62–71% of its peak across the plateau, from setosa's colour collisions with the larger cluster.
- **The three-species stage is too short to make a plateau.**
  - The versicolor-plus-virginica cluster starts breaking near T = 0.05; at T = 0.0858 it is present in half the draws,
    which is the MagnetizationVariance peak.
  - By T = 0.10 it is gone, and its pieces are already falling below 32 points. By T = 0.12 the system has melted.
  - The three species coexist for roughly 0.2–0.3 e-folds, under a split's own width w_q ≈ 0.75. That is Thread I's
    shoulder case: a cluster that lives less than w_q makes no bump of its own, and here its split overlaps the melt.
- **How the papers read it,** as recorded in the parity battery (`BwdParityTests.cs`):
  - **BWD 1996 read stage 1 off the plateau,** not off a peak: setosa clean, versicolor and virginica merged, 125 correct
    and 25 unclassified. SPC clusters from the spin–spin correlations at a temperature inside a plateau. A peak only
    marks the passage between two plateaus, and the plateau's hot cliff bounds where to read it.
  - **Domany 1999 (§3.1) got the species with a second stage:** SPC rerun on the merged cluster alone, about 2
    unclassified.
  - **The parity battery** reproduces both with FkReduced's plateau-and-cliff landmarks and with the T-stack dendrogram
    (CutToK(2) the cold 50/100 stage, CutToK(3) the brief hot three-way).
- **What it means for regime A's readout.**
  - Peaks alone miss iris's structure. The readout has to read plateaus and cliffs as well. S₂ carries every event as a
    step: the steepest points of S₂, the roots of its second derivative, are a closed-form readout the spline already
    supports. The current readout only reads maxima.
  - That extension would be engine work, parked with the readout design.
- **What it means for regime B** (a hypothesis, unverified).
  - A cluster that lives less than w_q can't be resolved by sampling it more densely: its split overlaps the melt in the
    physics, whatever the grid.
  - Domany's second stage changes the system instead. A global rescale of the couplings wouldn't separate the split from
    the melt, since it moves both temperatures by the same factor. The rebuilt graph's own bandwidth changes the
    contrast between near and far edges nonlinearly, and that can.
  - If so, regime B for short-lived clusters is a re-graphed subset, not a denser sweep.

### 12.17 The descending test (Q4), designed and checked

**The design.**
- **Each interval is tested once.** Its drop ΔS₂ and Welch standard error are fixed once both endpoints exist, and new
  volleys only add children. So the only multiplicity is the growing tree.
  - Levels spend α_v = 6α/(π²(v + 1)²), split evenly over the level's intervals.
  - The union bound then holds at whatever volley the schedule stops. That is the anytime validity §6 asked for,
    without an e-process.
- **Three states,** with t quantiles at the Welch degrees of freedom (8 replicates per temperature):
  - **Q, certified quiet:** ΔS₂ + t·SE < floor. It is hereditary: S₂ is monotone, so every sub-interval is quiet in
    truth, and descendants inherit Q untested.
  - **A, active:** ΔS₂ − t·SE > floor.
  - **I, indeterminate:** neither.
- **Localization.** Events are the maximal runs of non-Q intervals at the current volley. Each is localized
  narrowest-over-threshold, after Baranowski, Chen & Fryzlewicz's NOT (cited from memory): at the narrowest dyadic
  interval, or adjacent pair, inside it whose lower bound clears the floor. An event spread over many fine intervals is
  then reported at the coarser width where its mass still concentrates, which ends C25's fragmentation.
- **What an indeterminate block means.** A block with no active interval is a measurement verdict, not a resolution
  one: it needs more replicates there, not more volleys (the two knobs of §2).

**Checked on the 10-point schedule at α = 0.05,** against an independent sweep at another seed (4242, same grid and
budget):
- **Validity.** No quiet certification on the main sweep (281 on the toy, 265 on iris, over all volleys) is contradicted
  by the independent sweep's active bound. The two sweeps agree on 85–100% of interval states, volley by volley.
- **Quiet share.** The certified share of ln T grows from 44% to 52% on the toy, and from 33% to 51% on iris, over
  v0–v5.
- **Localization on the toy.**
  - v0–v3 report one block over the active region. From v4 it splits into the split region [0.0060, 0.0253], localized
    to 0.08 e-folds, and the melt region [0.0304, 0.1100], localized to 0.04.
  - The certified-quiet gap between them, [0.0253, 0.0304], is the three-stripe plateau. The clustering window falls
    out of the certificate as a certified-quiet run between active blocks.
- **Localization on iris.**
  - The main block [0.0304, 0.1411] localizes to 0.04 e-folds by v5.
  - The stage-1 plateau stays largely indeterminate. Setosa and its partner shed loose points continuously, about 0.1
    of S₂ over two e-folds, which at 8 replicates keeps many fine intervals' bands across the floor. More replicates
    there would certify it.
- **Inventory, against the independent sweep's bump peaks.**
  - Iris: 5 of 5 lie inside non-quiet intervals at v5.
  - Toy: 1 of 3. The MagnetizationVariance split peak is inside. The FkReduced lifetime peak sits on the certified
    plateau, as a lifetime should (§12.11). The heat-capacity melt peak (T = 0.117) sits where S₂ has stopped moving,
    because by then the stripes have broken into fragments below the 160-point floor.
  - The certificate certifies the size floor, not every χ bump.

### 12.18 The BARS settings sweep, scored against the independent sweep

**The run.**
- **What varied.** One knob at a time from BARS's defaults, on every volley of the 10-point schedule (10 to 289 points),
  for FkReduced, MagnetizationVariance and heat capacity on both datasets: 10 settings × 6 volleys × 3 curves × 2
  datasets.
- **The truth** is the independent sweep (seed 4242), which no fit saw.
  - Error is the RMSE against its curve, as a share of its range.
  - Registration and coverage are scored against its bumps (by topographic prominence) and their peaks, within one dense
    grid step.
- **Replicate weights** entered through a fixed-weights stand-in for the robust weight block: 1/SE², with SE clamped at a
  fifth of its median, normalized to mean 1.
- **One combination can't run.** The engine refuses tempering together with a weight resampler (or with adaptive τ), so
  weights-plus-tempering doesn't run.

| Setting | Median error against the default | Bump registrations (of 48 bump-volleys) | Peak intervals covering the truth | Fits with fit R̂ > 1.1 (of 36) | Mean knots | Seconds per fit |
|---|---|---|---|---|---|---|
| Defaults | 1.00 | 28 | 79% | 14 | 5.4 | 1.5 |
| Replicate weights | 7.37 | 3 | 100% (3 intervals) | 0 | 3.4 | 1.5 |
| Support gate 0.5 | 1.00 | 26 | 77% | 13 | 5.9 | 1.3 |
| Support gate 0.25 | 1.00 | 27 | 78% | 12 | 6.0 | 1.3 |
| Poisson(10) knot prior | 0.95 | 29 | 86% | 13 | 6.2 | 1.4 |
| Uniform knot prior, up to 20 | 1.01 | 25 | 76% | 14 | 5.6 | 1.4 |
| Log response | 1.57 | 17 | 88% | 7 | 4.5 | 1.3 |
| Quadratic basis | 0.97 | 24 | 71% | 18 | 5.5 | 1.1 |
| Adaptive kernel | 1.00 | 26 | 85% | 13 | 5.5 | 1.3 |
| Tempering, 4 rungs | 1.00 | 28 | 79% | 10 | 5.4 | 2.0 |

- **From about 145 points, the default fits are already at the measurement noise** (C28).
  - The two independent dense sweeps disagree with each other by 0.4–2.7% of the range: toy FkReduced 0.83%,
    MagnetizationVariance 1.35%, heat capacity 2.74%; iris 0.39%, 2.64%, 2.00%.
  - The default fits at v4 (145 points) sit at 0.7–1.9%. That is at or below the floor for MagnetizationVariance and
    heat capacity, since smoothing beats the raw points, and slightly above it for FkReduced.
  - A tighter fit there needs more replicates, not other settings.
- **Settings matter at the coarse volleys.** At v1 (19 points) the default errors run 5.7–8.7%, against floors of
  0.4–2.7%. Poisson(10) and the quadratic basis recover some of the gap:
  - toy MagnetizationVariance 0.080 → 0.042 (quadratic, though with fit R̂ 1.11);
  - iris FkReduced 0.063 → 0.051 (Poisson(10)).
- **Poisson(10) is the one consistent small gain.** Its median error is 0.95 of the default's, with the most
  registrations (29) and better coverage (86%). With tempering and with the support gate at 0.5, it also registers
  iris's MagnetizationVariance peak a volley earlier, at v1 instead of v2.
- **Replicate weights, as BARS takes them now, do harm** (C29; the mechanism is a hypothesis).
  - The coefficient prior stays the unweighted g-prior, β ~ N(μ₀, mσ²(ZᵀZ)⁻¹), while σ² becomes the noise of a
    unit-weight point.
  - With weights spread over orders of magnitude, the flat, precisely measured stretches set σ² small. The prior then
    shrinks the curve toward flat, peaks included.
  - §1's replicate weights need a prior that doesn't borrow the weighted noise scale, or a different normalization,
    before they can be used.
- **The log response is a trade.** It gives up held-out error (1.57×) and registrations (17, since the 10% prominence
  gate reads differently on a log scale). It gains the best location error (0.018 e-folds), the best coverage (88%) and
  the best convergence (7 fits over R̂ 1.1, against 14).
- **Tempering helps mixing modestly:** 10 fits over R̂ 1.1 against 14, at 1.3× the time. Iris's heat-capacity fits stay
  unconverged under every setting (fit R̂ up to 2.95).
- **The support gate doesn't bind** at these spacings: at 0.5 or 0.25 the errors are unchanged.

This is one pre-scored pass on two datasets: exploration, not a basis for defaults.

### 12.19 The stopping rule as it stands (a synthesis; nothing decided)

**What survived from the original threads.**
- **Thermal intervals stayed, with a different observable.** The rule still works interval by interval on the
  temperature grid, from endpoint means only. The size-blind energy quantity (D_J, the thermodynamic length bound) gave
  way to pair connectivity S₂, which is monotone and weights events by size (§12.3–12.4).
- **Fisher zeros stayed as explanation and budget, not as a test.** They explain why refinement converges and set the
  switch widths behind Thread I's a priori count. They can't be estimated from real-axis data in time to stop on (C7).
  E1 showed the rigid-split widths they imply don't bound melts or background-mediated splits.

**Before volley 0.**
- Set the bracket so the window holds every significant event. Check its cold edge against the spanning tree's
  closed-form S₂ lower bound (C23), and against the graph's pre-repair components (C19).
- Take a minimum volley count from the narrowest feature that matters, planned on the worst lattice phase (C27).

**After each volley.**
1. **Certify (§12.17).** Classify every interval as quiet, active or open by its drop in S₂ against the size floor,
   spending α over the levels.
   - Quiet passes down the tree: nothing of significant size can hide inside a quiet interval.
   - Active runs are events, each localized at the narrowest interval that still clears the floor.
   - An open interval needs more replicates there, not more volleys.
2. **Read (§12.15–12.16).** BARS on the detection curve (MagnetizationVariance) for the giant's transitions, and
   FkReduced or S₂ for spans and cliffs. Plateaus come out as the certified-quiet runs between active blocks.
3. **Stop when all of these hold:**
   - every interval is quiet, or lies inside an active block localized within the regime-B window;
   - no open interval remains, after adding replicates where the test asks for them;
   - every active block contains something the readout registered: a peak, a plateau edge or a cliff;
   - one afterburn volley has run past the last registration (§12.12).

**What a stop certifies.** No merge of two clusters at or above the size floor hides in any quiet interval, with a
family-wise error at most α at whichever volley the rule stops. It does not certify every χ bump: lifetime peaks sit on
quiet plateaus by construction, and melts of fragments below the floor go unflagged.

**What sits beside the rule.** The two-quiet-volleys surprise rule stays as a soft check, since it stopped once too
early (C12). The engine side (OQ7, OQ8) stays parked with the readout design.

### 12.20 Iris across bandwidths: does a second axis open the species split?

**The run.**
- **The slices.** Iris on its fixed edge set (k = 10 mutual KNN with MST repair), with the Gaussian bandwidth scaled by
  λ ∈ {0.5, 0.7, 1, 1.4, 2}. Since J = e^(−d²/2a²)/K̂, the scaled couplings are exactly J_λ = (K̂·J)^(1/λ²)/K̂; here
  K̂ = 7.04, and the strongest edge has K̂·J = 1.
- **The sweep.** 145 temperatures from 10⁻⁴ to 2, 8 replicas × 2,000 cycles, with per-edge co-membership. The
  accumulator counts each undirected edge on one slot only (the j > i one), so a mean over both slots reads half the
  true connectivity.
- **The species window** is where the versicolor–virginica cut edges are mostly broken (mean co-membership < ½) while
  both species hold (> ½).

| Bandwidth | Cut breaks at T | Versicolor, virginica melt at T | Setosa melts at T | Species window | Most 32–63-point clusters per draw |
|---|---|---|---|---|---|
| ×0.5 | 0.020 | 0.030, 0.021 | 0.079 | 0 | 1.83 |
| ×0.7 | 0.058 | 0.067, 0.057 | 0.120 | 0 | 1.77 |
| ×1 | 0.100 | 0.109, 0.099 | 0.152 | 0 | 1.65 |
| ×1.4 | 0.131 | 0.138, 0.130 | 0.167 | 0 | 1.37 |
| ×2 | 0.153 | 0.159, 0.152 | 0.177 | 0 | 1.13 |

- **No bandwidth opens a species window.** At every scale the cut breaks within 0.05 e-folds of the species' own melts,
  against a split's width of about 0.75. Three ~50-point clusters never coexist: the 32–63-point bin peaks at 1.1–1.8
  per draw.
- **What the bandwidth axis does do is separate clusters of different tightness.** Setosa's edges are short (median
  d²/2a² = 0.26), so its melt pulls away from the versicolor–virginica event as the bandwidth shrinks: 1.37 e-folds
  apart at ×0.5, 0.15 at ×2.
- **Why it can't separate these two.**
  - The cut edges are as long as the species' own: median d²/2a² of 0.71 for the 34 cut edges, against 0.53
    (versicolor) and 0.68 (virginica). Scaling the bandwidth scales all of them nearly alike.
  - The cut is also thick. 27 of the 100 points (13 versicolor, 14 virginica) have an edge across, and the cut carries
    17–21% of each species' internal coupling.
  - On this graph, versicolor and virginica are one elongated cluster with a boundary, not two clusters with a bridge.
- **Where a second axis would have to act instead.**
  - On the edge set: the neighbour count K, a density-aware scaler, or the feature space itself.
  - BWD 1996's sparser K ≥ 5 graph found a three-species stage. Plain mutual 5-NN iris has 11 components, so its
    MST-repaired versicolor–virginica link is thin by construction.
  - That points a grid search at the graph's construction parameters, scored by the cut's co-membership against the
    species' own, rather than at the bandwidth alone.

**The neighbour count K** (mutual KNN with MST repair, the default bandwidth rule and 1/K̂; same sweep):

| K | Edges (cut) | Pre-repair components | Cut breaks at T | Versicolor, virginica melt at T | Gap to the first melt | Species window |
|---|---|---|---|---|---|---|
| 5 | 262 (13) | 16 | 0.088 | 0.113, 0.102 | +0.15 | 0.14 |
| 7 | 353 (19) | 8 | 0.095 | 0.113, 0.098 | +0.04 | 0.07 |
| 10 | 528 (34) | 6 | 0.100 | 0.109, 0.099 | −0.01 | 0 |
| 15 | 824 (52) | 4 | 0.100 | 0.104, 0.097 | −0.03 | 0 |

- **Sparser graphs open only a sliver.** The cut stays at 7–10% of the two species' edges at every K, and only one cut
  edge comes from the MST repair. At the cut's half point virginica's own edges hold 0.51–0.61 of the time: virginica
  is the weak side.
- **Edge-averaged windows understate the cores, though.** The mean over all internal edges includes the boundary points.
  BWD 1996's three-way stage left 25 points unclassified, which is to say the cores separate after the boundary sheds.
  - The size bins show that kind of stage: about 3.2 clusters of 16 or more points per draw at K = 5 and T ≈ 0.06, and
    about 2.7 at K = 10 and T ≈ 0.097, with many points shed into small clusters.
  - Whether those are species-pure cores needs BWD's own partition rule (co-membership above ½ on edges, then
    components) and a purity count; see the purity check below.

**BWD's partition rule through the three-core stage.** Per-edge co-membership, pooled over 8 replicas, at 33
temperatures from 0.02 to 0.16. Edges above ½ are kept, components taken. A cluster of 20 or more points counts as
real, residue is the points in clusters under 5, and "correct" is the sum of the real clusters' majority counts.
- **K = 5 has a three-species stage,** from T = 0.0465 to 0.0687, about 0.39 e-folds.
  - At its onset: setosa 43 (pure), versicolor 38 (pure), virginica 36 (33 virginica and 3 versicolor). That is 114
    correct, with residue 28 and 5 more points in small clusters.
  - BWD 1996 read 45/40/38 with 25 unclassified, so this graph reproduces their stage closely.
  - Through the stage the virginica core shrinks to 21 points while versicolor holds at 34–38.
  - Its MagnetizationVariance peak (T = 0.0455) sits at the stage's onset.
- **K = 10 never reaches it.** The merged cluster sheds from 90 to 67 points (T = 0.02 to 0.089), then loses virginica
  to fragments: 47 points (38 versicolor, 9 virginica) at T = 0.095. Virginica melts before it separates as a cluster.
- **So the neighbour count is the second axis that matters for iris, read at the level of cores.**
  - The split happens between cores once the boundary points have shed, which is why the edge-averaged window above
    understated it.
  - On a sparser graph, virginica's core holds together long enough to come apart from versicolor's.
  - Bandwidth and kernel shape don't reach this: they reweight the same edges by length.
  - The stage is still narrower than a split's width w_q ≈ 0.75. So in χ it shows as a shoulder rather than its own
    plateau, and regime A would read it from MagnetizationVariance's peak with a narrow window after it.
- **The parity battery already reads iris at K = 5.** Its Domany 1999 iris tests (capture, the two-stage hierarchy and
  lineage persistence) build k = 5 graphs. Only the BWD 1996 flat test and the profile dump use k = 10, which this
  thread followed until here.

### 12.21 Lineages across temperature

Final assignments come from tracing clusters across temperatures, and the tree already holds the machinery under
`src/clustering/graphical/spc/partitions/hierarchical/`:
- **`DenseTStack`** cuts BWD's partition at every grid temperature, from pooled co-membership. It is schedule-agnostic
  and accepts non-uniform grids.
- **`HierarchyEom`** checks that the stack nests, restores nesting when sampling noise breaks it, bridges the stack to a
  dendrogram and selects clusters by excess of mass.
- **`LineagePersistence`** tracks each cluster lineage across the stack, keeping a cluster that erodes at its periphery
  distinct from one that splits. It scores persistence as levels survived × decidedness.

**How this thread's results sit under it.**
- **Monotone connectivity makes the stack a merge tree.** Every pair connectivity τ_ij is monotone in T (C17), so the
  edge sets kept by BWD's threshold are nested in T, in expectation, and so are their components.
  - The T-stack is then a dendrogram, and a lineage's lifetime is the persistence of a one-parameter filtration.
  - Nesting breaks only through sampling noise near the threshold, which is what `HierarchyEom`'s nesting gate catches,
    or through peripheral capture's argmax step, which isn't monotone and so belongs to the zigzag engine.
- **The certificate localizes the tree's nodes.**
  - Active blocks of §12.17 are where splits and melts happen, the tree's nodes. Certified-quiet runs are its edges,
    the plateaus.
  - So regime A finds the nodes coarsely, regime B samples densely across each node with co-membership, and the T-stack
    resolvers trace the lineages through those windows to the assignments.
  - On iris at K = 5 that is one node near T = 0.046 (the species split) and the melt sequence after about 0.069.
- **One interaction to settle with regime B's design** (engine side, parked).
  - Persistence measured in levels survived depends on the grid once sampling is uneven: a lineage inside a densely
    sampled window survives many more levels than an equally long-lived one elsewhere.
  - The code's reason for levels over spans (an absolute span over-rewards the hot end of a log grid, a log span the
    cold end of a linear one) holds for uniform grids. SPC × BARS's grids won't be uniform.
  - A lifetime in ln T, or in thermodynamic length, is a property of the lineage rather than of the sampling, and would
    keep persistence comparable across a sparse-plus-dense schedule.

### 12.22 The T-stack resolvers on iris at K = 5

**The run.**
- **The sweep.** One rich sweep (co-membership and the cluster-size landscape) on 129 log-spaced temperatures from 0.005
  to 0.3, with 8 replicas × 2,000 cycles.
- **The resolvers.** `HierarchyEom` (θ = ½, minimum cluster 5) and `LineagePersistence` (θ = ½, minimum 5, no window),
  each run on subsamples of that same sweep.
- **Species counts** are written setosa/versicolor/virginica.

| Grid design (levels) | `HierarchyEom` clusters | Unassigned | `LineagePersistence` selection |
|---|---|---|---|
| Uniform (129, 65 or 33) | 92 [0/46/46], 49 [49/0/0] | 9 | the merged 87–74 [0/43/44], two 5-point setosa remnants |
| Dense in [0.04, 0.08], sparse outside (49) | 92, 49 | 9 | merged 74, a 5-point virginica core, setosa remnants |
| Sparse in [0.04, 0.08], dense outside (113) | 92, 49 | 9 | as uniform |
| Cold edge at 0.0129, the parity bracket's (99) | 88 [0/44/44], 48 [48/0/0] | 14 | as uniform |
| Cold edge at 0.025 (78) | 48 [48/0/0], 37 [0/37/0], 21 [0/0/21], 14 [0/3/11], 5 [0/0/5] | 25 | one 5-point setosa lineage (145 unassigned) |
| Cold edge at 0.035 (68) | 37, 21, 14, 13, 9, 5, 5 (setosa in pieces) | 46 | one 5-point setosa lineage |

- **Neither resolver selects the species on its own.** Persistence and excess of mass both favour the long-lived coarse
  stage.
  - The merged versicolor-plus-virginica cluster lives about 2.2 e-folds on this grid, and in truth back to setosa's
    split near 2×10⁻⁶. Its species children live 0.45–0.65 e-folds: the versicolor core (31–33 points) about 0.45–0.48,
    the virginica core (15–17) about 0.5–0.65.
  - The species lineages are all in `LineagePersistence`'s list but never selected.
- **The grid's cold edge, a bracket choice, decides `HierarchyEom`'s answer.** Starting the grid at 0.025 cuts most of
  the merged cluster's cold levels. Its children then win on excess of mass: 48/37/21 with 25 unassigned, close to BWD
  1996's 125 correct and 25 unclassified. Starting at 0.035 cuts into setosa's own stable stage, and setosa fragments.
- **Density inside the window moves `LineagePersistence`'s ranking.** Densifying [0.04, 0.08] lifts the species
  lineages' scores (the versicolor core from 3.8 at 4 levels to 15.2 at 16) and changes the selection. Their ln-spans
  stay put across every design (versicolor 0.45–0.48, merged 2.18–2.24).
  - Scoring in levels makes the final clustering a function of the sampling design, both its density and its extent.
    Lifetimes in ln T are design-invariant, as §12.21 expected.
- **`LineagePersistence` is fragile off its tested grid.** With the cold edge at 0.025 or 0.035 it selects a single
  5-point lineage and leaves 145 points unassigned. It also represents setosa by a 5-point remnant: its representative
  is the lineage's least-decided, most eroded level. The parity test passes on its own grid and window.
- **What it means for the final assignment.** On iris, persistence picks the coarse stage: two clusters, the species
  merged. The species are a short-lived refinement of it.
  - Reporting three species means choosing a level, say a leaf-style selection, or a scale of interest. That is a
    modelling choice, not something a stability criterion settles. Here it was the bracket that made the choice,
    silently.
  - Design-invariant lifetimes would make that choice explicit, but they wouldn't make it for you. In ln T the merged
    cluster outlives its children by far.

### 12.23 A Mahalanobis metric for iris (static check)

**The question** (the user's): GMMs with full covariance give the best iris results the user knows of, so would SPC on a
Mahalanobis metric do better?
- **The reduction.** With a fixed Σ, Mahalanobis distance is Euclidean distance after x → L⁻¹x (Σ = LLᵀ). That makes it
  a change of the edge set and of the edge lengths, unlike the bandwidth or the kernel's shape.
- **The check.** Mutual 5-NN geometry, without the MST repair, for five metrics:

| Metric | Cut edges | P(cut edge longer than a versicolor edge, a virginica edge) | Cut's coupling share (against versicolor, virginica) | 1-NN accuracy, versicolor against virginica |
|---|---|---|---|---|
| Euclidean | 12 | 0.74, 0.56 | 0.117, 0.134 | 0.94 |
| Standardized features | 15 | 0.81, 0.58 | 0.133, 0.169 | 0.92 |
| Global Mahalanobis (all data) | 31 | 0.72, 0.58 | 0.353, 0.466 | 0.87 |
| Within-core Mahalanobis (SPC's cores, no labels) | 10 | 0.80, 0.60 | 0.090, 0.117 | 0.95 |
| Within-class Mahalanobis (true labels, oracle) | 10 | 0.81, 0.55 | 0.095, 0.122 | 0.95 |

- **The global covariance hurts.** Its large directions are the ones the species differ along, so whitening by it shrinks
  exactly the separating directions: the cut grows to 31 edges and its coupling share roughly triples.
- **The within-cluster covariance helps, modestly.** The cut's coupling share falls by about 23% against versicolor and
  13% against virginica.
- **Estimated from SPC's own cores, it matches the label oracle.** The cores are BWD's rule at the start of K = 5's
  three-species stage: 43, 38 and 36 points. So an unsupervised alternation is viable: SPC, then cores, then their
  pooled covariance, then SPC again.
- **It doesn't remove the overlap.** A cut edge is still longer than a virginica edge only about 55–60% of the time.
  Versicolor and virginica overlap in feature space under any linear metric. Connectivity clustering needs a gap or a
  density dip, and here there is none to find.
- **Why a GMM does better here.** It models each class's covariance and assigns the overlap probabilistically, a
  parametric model of the overlap rather than a connectivity criterion.
  - That argues for the user's hybrid: SPC for the cores, where connectivity is unambiguous, and the GMM to fill in the
    overlap.
  - Per-cluster covariances inside SPC would make the metric move with the clustering: the adaptive scheme the first
    dyadic discussion warned degenerates. An alternation run a fixed number of rounds keeps clear of that.
- **Not yet run:** SPC itself on the within-core-whitened features, scored by BWD's rule as in §12.20.

### 12.24 Is versicolor–virginica a projection artifact?

**The question** (the user's intuition): two manifolds can look overlapped in one projection and come apart in another,
so might a non-Euclidean distance or another embedding separate these two species?
- **There is no hidden dimension.** Fisher's iris has exactly four measurements per flower, and every run here used all
  four. The petal measurements carry the separation:
  - Fisher ratios of 4.37 (petal width) and 3.24 (petal length), against 0.65 and 0.21 for the sepals;
  - even so, 37–38% of the two species fall in the range they share on each petal feature.
- **Every unsupervised axis of the merged cloud has one mode.** Along each, with errors at the best single threshold (of
  100 flowers):
  - its main axis (84% of the variance): 12 errors;
  - log-size, the first principal component of the log measurements: 6 errors;
  - the log-shape axes: 32 and 41 errors.
- **The species nearly separate along one direction.**
  - Fisher's discriminant, which needs labels, misplaces only 2 of the 100. Even there the density has only a shallow
    dip, bottoming out at 75% of the lower mode.
  - Learned from SPC's own cores instead (38 points and 36, without labels), the direction lies within about 12° of
    the label oracle's (cos 0.979). It misplaces 3, with a dip at 73% of the lower mode.
- **So it is not a hidden-manifold case.** A nonlinear embedding (geodesic, diffusion or hyperbolic) rearranges the same
  neighbourhoods; it can't open a gap that the four measurements don't hold. The separation that exists is a thin
  margin along one direction, with only a mild density dip, which a connectivity method sees as one cluster and a margin
  or a mixture model reads well.
- **The biology fits** (cited from memory, to confirm): Anderson (1936), whose measurements these are, proposed
  *I. versicolor* as an allopolyploid derived from *I. virginica* and *I. setosa*. That would make it morphologically
  intermediate, on a size continuum with virginica, rather than a separate manifold.
- **It supports the user's hybrid quantitatively.** SPC's cores already carry nearly all the information needed to
  learn the separating direction: 3 errors against the oracle's 2. SPC finds the cores, and a discriminant or a GMM
  trained on them assigns the rest.

**The SPC → GMM hybrid, as the user has run it in spike sorting.** A GMM wins on iris because it is given K and fits
full covariances to a maximum of the likelihood. SPC supplies K, and labels from points near dense cores. Two
refinements from iris:
- **K = 3 needs the right level of the hierarchy.** Persistence and excess of mass chose 2 on iris unless the bracket cut
  the merged cluster (§12.22). The three-species stage lasts about 0.4 e-folds at K = 5 and is absent at K = 10. Spike
  units are usually better separated, which makes K easier there.
- **Core covariances are truncated.** Against each full species, the cores' covariance volumes are 0.29× (setosa),
  0.22× (versicolor) and 0.10× (virginica), with traces of 0.68×, 0.70× and 0.41×.
  - Fixed as class models, the narrowest core, virginica's, would lose its boundary points.
  - Used as EM's starting point on all the points instead, the covariances can grow back to the species' spread.

### 12.25 SPC on its own: a first pass over q, K and the coupling

The user's focus: what generalized SPC can do without the hybrid, through choices in modelling the data and building
the graph.

**The run.**
- **The grid.** Iris over q ∈ {20, 100, 1000} × K ∈ {5, 7} × four couplings on the mutual-KNN edge set with MST repair.
  193 temperatures from 10⁻⁴ to 3, 8 replicas × 2,000 cycles, one seed.
- **The couplings:**
  - the default Gaussian;
  - the compiler's LMP rescaling;
  - shared nearest neighbours, J = max(|N(i) ∩ N(j)|/K, 0.05)/K̂;
  - mutual reachability, a Gaussian on max(core distance of i, core distance of j, d).
- **The scores,** fixed before running:
  - the width of the three-species stage under BWD's rule (§12.20), where at least three clusters of 20 or more points
    have distinct majority species;
  - its best correct count and residue;
  - whether `HierarchyEom` selects three species-led clusters on the whole grid.

**Two cells can't be read as run.**
- **q = 1000 isn't equilibrated at the cold end.** S₂ at T = 10⁻⁴ is 0.22–0.36 where at least the 0.556 setosa floor
  belongs. With 1,000 colours a random start almost never aligns neighbours, so SW orders very slowly at low
  temperature. It needs an ordered start before it says anything.
- **LMP fragments iris.** Its strongest coupling (0.69) exceeds the Gaussian's (0.29), yet S₂ at T = 10⁻⁴ is only 0.09
  (K = 5) and 0.15 (K = 7): many edges are driven to near-zero coupling, leaving about ten pieces at any temperature.

**The readable cells:**

| K | Coupling | q | Stage (e-folds) | Best reading under BWD's rule | Correct | Residue | `HierarchyEom` |
|---|---|---|---|---|---|---|---|
| 5 | Gaussian | 20 | 0.43 | 43 [43/0/0], 38 [0/38/0], 36 [0/3/33] | 114 | 28 | two clusters |
| 5 | Gaussian | 100 | 0.38 | 40, 38, 35 [0/3/32] | 110 | 32 | two clusters |
| 5 | Shared neighbours | 20 | 0.75 | 67 [0/18/49], 34 [34/0/0], 33 [0/32/1] | 115 | 1 | **three:** 67 [0/18/49], 34 [34/0/0], 33 [0/32/1] |
| 5 | Shared neighbours | 100 | 0.64 | as at q = 20 | 115 | 1 | two clusters |
| 5 | Mutual reachability | 20 | 0.21 | 41, 36 [0/3/33], 35 [0/35/0] | 109 | 38 | two clusters |
| 5 | Mutual reachability | 100 | 0.21 | 40, 34, 21 | 95 | 41 | two clusters |
| 7 | Gaussian | 20, 100 | 0 | — | — | — | two clusters |
| 7 | Shared neighbours | 20 | 0.16 | 56 [0/43/13], 34 [0/0/34], 25 [25/0/0] | 102 | 17 | two clusters |
| 7 | Mutual reachability | 20, 100 | 0 | — | — | — | two clusters |

- **q doesn't widen the stage** between 20 and 100: 0.43 to 0.38 for the Gaussian, 0.75 to 0.64 for shared neighbours.
  - The rigid-split width w_q shrinks with q, but the stage's length is set by the ratio of the split's temperature to
    the melt's, and both shift together.
  - So the hope of a wider stage at larger q, drawn from Thread I's w_q, fails between 20 and 100. q = 1000 is untested
    until it is equilibrated.
- **Shared-nearest-neighbour couplings are the one lever that moved the readings.** At K = 5 and q = 20:
  - the stage widens to 0.75 e-folds, about a split's own width w_q;
  - `HierarchyEom` selects three species-led clusters, the only configuration here where it does;
  - the residue falls to 1.
  - The cost is in assignment: 18 versicolor sit in the virginica cluster, and setosa splits, with 34 in the real
    cluster and the rest in a smaller one.
- **Mutual reachability narrows the stage** (0.21) and sheds more points. Iris's density dip is too shallow for core
  distances to sharpen it.
- **K = 7 mostly closes the stage again;** only shared neighbours keep a sliver (0.16).
- **One pass, one seed: exploration, not a basis for defaults.**
- **Natural next steps:**
  - shared-neighbour variants: the neighbour count separate from the graph's K, the floor, Jaccard normalization;
  - q = 1000 from an ordered start;
  - the readout's θ and peripheral capture, applied to the shared-neighbour stage.

### 12.26 Per-cluster covariances across temperature

**The idea** (the user's): bring the GMM's insight inside SPC as an observable. Compute each cluster's covariance at
each temperature and follow how it evolves.
- **The hypothesis.** Shedding from a cluster's outer edge shrinks its variance in every direction. Shedding from the
  saddle between two hidden lobes empties the middle, so variance along the axis joining them holds up while the others
  shrink, until the split.
- **The data.** BWD's partition at each of 33 temperatures (K = 5, Gaussian couplings, the per-edge co-membership of
  §12.20), with 4×4 covariances of each cluster of 20 or more points. The label-based separating direction is used only
  to check.

**What the merged versicolor-plus-virginica lineage does before it splits** (T = 0.020 → 0.044, 90 → 74 points):
- **Variance along the eventual split axis holds while the rest shrinks.**
  - Total variance falls by 32%: the top eigenvalue 0.906 → 0.590, the third 0.068 → 0.040.
  - Variance along the separating direction stays near 0.23–0.26, so its share of the total rises from 0.23 to 0.30.
  - At the split (T = 0.0465), the children's shares drop to 0.07 (versicolor core) and 0.16 (virginica core).
- **No single snapshot shows it.** The top axis stays unimodal (bimodality coefficient 0.41–0.44, under the 0.555
  threshold), and it isn't the split axis (cos 0.45–0.53).
- **Without labels, the direction is only loosely identified.** The direction where variance holds up best between two
  temperatures solves C(T₂)v = λC(T₁)v.
  - Its top ratio barely leads the second (1.07–1.08 against 1.02–1.06), while one direction loses 50–67%.
  - So the trajectory narrows the split axis to a plane rather than a line. The leading vector's agreement with the
    label direction wanders with the window (cos 0.44–0.93).
  - Projected on it, the 100 versicolor and virginica flowers give 11–29 errors at the best threshold, against 3 for
    the direction learned from the cores after the split. The margin is thin enough that a 20° tilt costs a lot.
- **What it is good for.**
  - A split precursor: a lineage whose variance share along some direction rises as it sheds is likely to split. That
    could tell regime B where to sample densely for hierarchy nodes.
  - Each lineage's shape becomes a feature tracked across T, which suits the process-visualization ambition.
  - For assignment, the children's covariances after the split, the cores, still give the sharp direction.

**The kernel's shape can't help here** (a static check on the K = 10 edge lengths, kernels matched in scale):
- A cut edge is longer than a virginica edge 51% of the time, and longer than a versicolor edge 70% of the time.
- A kernel is a monotone map from distance to coupling, so it can only reweight by length. The cut's share of coupling
  barely moves across kernels: 0.169 and 0.212 (against versicolor's and virginica's interiors) for the Gaussian,
  0.168 and 0.205 for the Laplace, 0.170 and 0.207 for the Cauchy, and 0.167 and 0.217 for a sharper e^(−r⁴/4).
- Kernels matter when cut edges are systematically longer, where a sharper kernel buys contrast. They also matter when
  long internal edges melt a diffuse cluster early, where a heavier tail keeps it whole. Iris's versicolor–virginica cut
  is neither: the information that separates them isn't in the edge lengths at this scale.
- **What could carry it:** density, through a density-aware rescaling such as LMP or mutual-reachability distances,
  which lengthen edges across a density dip; or the features themselves.
