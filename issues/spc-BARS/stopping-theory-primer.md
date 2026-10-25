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
| C12 | Two quiet volleys reduce aliasing | Argon notes | Heuristic | E4; anytime-valid reformulation (§6) |
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
| C25 | Per-interval ΔS₂ tests lose power at fine volleys; a test that descends the dyadic tree keeps it | §12.9 | Observed (E5, cheap budget); the remedy is derived, unverified | Q4 design |
| C26 | BARS registers a bump at about one point per FWHM, roughly placed; its prominence gate, not the grid, drops weak bumps | §12.10 | Observed (E4, six bumps, equal weights) | E3; weighted fits |

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
- **Registration.** A bump counts as registered when a peak mode with posterior mass ≥ 0.5 falls inside the dense bump's
  FWHM. Every registration here reached mass ≥ 0.9 at once.

| Bump | FWHM | First volley | h there | Points per FWHM | Location error | Interval width |
|---|---|---|---|---|---|---|
| toy FkReduced, T = 0.025 (the three-stripe lifetime) | 1.57 | 3 | 0.75 | 2.1 | +0.17 | 0.69 |
| toy MagnetizationVariance, T = 0.014 (the split) | 0.41 | 4 | 0.37 | 1.1 | +0.13 | 0.60 |
| toy heat capacity, T = 0.110 (the melt) | 0.31 | 4 | 0.37 | 0.8 | +0.16 | 0.27 |
| iris MagnetizationVariance, T = 0.091 | 0.24 | 5 | 0.19 | 1.3 | −0.06 | 0.36 |
| iris heat capacity, T = 0.100 | 0.23 | 4 | 0.37 | 0.6 | +0.03 | 0.24 |
| iris heat capacity, T = 0.160 | 0.15 | never | — | 1.6 at v = 6 | — | — |

Widths, errors and interval widths are in e-folds of ln T.

- **BARS registers a bump at about one point per FWHM** (0.6–2.1). That is below Thread I's geometric two points per
  FWHM (C5), so on these data the geometric constant is conservative. The closed-form readout does what §1 says: it
  places a peak between grid points.
- **At first registration the location is rough:** errors of 0.03–0.17 e-folds, with intervals 0.24–0.69 wide. By v = 5
  the intervals run 0.09–0.36 wide.
- **The location intervals cover the dense peak** (to within its own grid quantization, ±h/2) in 12 of 13 fits from
  first registration on. The exception is iris's heat capacity at v = 5, which misses by about one grid step.
- **An under-resolved volley can put a mode in the wrong place.** At v = 3 the toy's MagnetizationVariance carries a mode
  of mass 0.57 at T = 0.024, half an e-fold off the split: enough to report, and misplaced. Its w/h there is 0.55.
- **The one miss comes from the prominence gate, not the grid.**
  - Iris's bump at T = 0.160 rises about 5% of the curve's range above its higher neighbour, under BARS's 10% gate, and
    is one point wide even at v = 6.
  - It is plausibly setosa's melt: setosa's edges run hotter (median T* 0.037, against 0.026).
  - S₂'s flagged run at v = 6, [0.027, 0.175], covers it. The certificate doesn't let it hide; the BARS readout just
    doesn't count it. That is the division of labour of §12.7.
- **No false bumps.** Apart from the misplaced mode above, no mode with mass ≥ 0.5 falls outside the dense features.
  Iris's FkReduced, which has no dense bump, never gets a mode above 0.2.
- **Not yet tested:**
  - C12's two-quiet-volleys rule, which needs each fit's predictive band at the next volley's midpoints;
  - an operating curve proper, which needs more than six bumps (E3's synthetic bumps, or more datasets);
  - the replicate-variance weights of §1, which need a weighted entry point.
