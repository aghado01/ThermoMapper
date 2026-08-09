# The confidence pushforward — BARS curve posterior → certified barcode

The **capstone** of the SPC-with-BARS chain. The seam between the **BARS posterior over the thermal
curve** and the **persistence barcode**, stated so the deep statistical/topological content is *cited*
and the composition is a *provable obligation* — same methodology as
[bifiltration-bridge-lemmas.md](bifiltration-bridge-lemmas.md) (quarantine the cited theorems, discharge
the plumbing).

Two earlier files prove the ends; this one joins them:
- **bifiltration-bridge** proves the *estimand* has structure: spin agreement is monotone in β
  (`spin_agreement_mono`) ⇒ the (β,θ) affinity filtration is **monotone + functorial** (T1/T2).
- **fable-BARS-lemma** + **bars-span** prove the *per-draw reduction is exact*: peak (`argmax_in_closed_form_set`)
  and span (`unique_crossing_in_monotone_bracket`, SP-2) are closed-form roots ⇒ **zero optimizer slop**.

What's left is the link that makes the guarantee *flow downstream*: a sup-norm credible band on the
curve pushes forward to a **bottleneck confidence ball on the diagram**, with no inflation. This is
"the posterior over T\* is an exact pushforward of the curve posterior" — lifted from the scalar peak
to the whole barcode.

## The chain, and where each link is grounded

| # | Statement | Status | Source |
|---|-----------|--------|--------|
| L0 | estimand monotone ⇒ filtration functorial | **proved** (mod Grimmett) | bifiltration-bridge (this dir) |
| L1 | per-draw curve → per-draw diagram, **zero optimizer slop** | **proved** (mod spline-def) | fable-BARS-lemma + bars-span (this dir) |
| L2 | **stability**: `d_b(Dgm f, Dgm g) ≤ ‖f − g‖_∞` | **AXIOM** | REF-PH §5.2 (stability for functions); **BCKL2010** §2.2 + §4 (statistical bottleneck + minimax rate); **PNV20XX** (the nerve readout's twin — persistent nerve lemma ⇒ `d_b`(nerve, space) bounded under ε-good cover) |
| L3 | sup-norm credible band ⇒ bottleneck confidence ball | **proved** (this file) | measure-monotonicity over L2 |
| L4 | the pushforward is a *consistent* posterior | **cited context** | **WRD2025** §4 (PH-posterior consistency + rate + misspec feature-recovery); **MNO2019** (PD = Poisson PP, conjugate posterior) / **MMO2019** (its frequentist twin — KDE + near-diagonal noise-likelihood + PHD = point-process intensity) |
| L5 | summarize the diagram-posterior (average the draws) | **OBLIGATION** | **FH2024** §4 (landscape stability) + §3.2 (stats); **MMO2019** §4 (KDE *density* on diagram space + a mean-absolute-*bottleneck*-deviation dispersion — keeps the random cardinality a landscape collapses); REF-PH §7.3 |
| L6 | K observables ⇒ multiparameter (metric, not bars) | **extension** | **REF-MPH** §9.3 (P-module robustness), §9.2.1 (no complete barcode — *guard*); **SGW2025** (multigraded Betti, computable) |
| ⊕ | certified feature **is** the phase transition | **target** | **MR2026** §4.2 (main theorem; persistent entropy); **DLST2026** (the dynamical-systems twin — phase transition = bifurcation, read as a Conley-Morse barcode via zigzag) |

Honest `sorry` count: **one cited** (L2 — the stability theorem; the random band's probability lives in
the BARS posterior, not formalized here), **one built** (L5 — the landscape map, once
`PersistenceLandscape` lands in `TDA.Ph`). L3 is fully discharged — see [Stability.lean](../Enthymemata/Stability.lean).

> **Corpus note.** **MMO2019** and **SGL2022** live in `compendia/ph` (moved from `intersections` this
> session — both PH-proper: statistical PH, and applied-PH-for-phase-transitions). **Acronym guard:**
> MMO2019's dispersion statistic is *mean absolute **bottleneck** deviation* — do **not** abbreviate it
> "MAD" here. "MAD" is taken: the **Median Absolute Deviation** in the `graphs/**` bandwidth-estimation
> stack — an *upstream* robust-scale estimator that sets the metric the filtration runs on (the L0 /
> DBK2023 leg). Same robust-deviation philosophy, opposite ends of the pipeline; keep them lexically apart.

## Notes for the proof pass

- **L3 is the whole point and it's a three-line measure argument.** `{‖f − f̂‖_∞ ≤ ε} ⊆ {d_b(Dgm f, Dgm f̂) ≤ ε}`
  by L2, so `μ(ball) ≥ μ(band) ≥ 1−α` by `measure_mono`. Event inclusion + monotone measure. *All* the depth
  is in L2 (cited); the bridge between honest-curve-confidence and honest-barcode-confidence is plumbing.
- **Zero-slack is L1's gift, and it is load-bearing.** Without the *exact* per-draw diagram, the ε would
  inflate by the optimizer error and the bottleneck ball would be larger than the credible band. Argmax/span
  exactness is precisely what makes the band-ε *equal* the bottleneck-ε. Exactness (L1) ∘ Lipschitz (L2) =
  honest confidence (L3) — neither half suffices alone.
- **Why BARS-only — the disjointed sweep can't even state L2.** Stability needs a *function* f to perturb;
  the classic per-T sweep has points, no f, no sup-norm, no L2, no pushforward. **WRD2025** is the prior art
  that *does* construct the functional object — it places a Bayesian model directly on the PH likelihood and
  proves posterior consistency (§4); BARS-on-f(T) is the thermal-axis instance, and this file is the stability
  layer that inference routes through. The disjointed sweep is outside the theorem's hypotheses by construction.
- **Two readouts, two stability legs.** L2 as stated is the *function/scalar* readout's leg — perturb the curve
  `f`, the diagram moves no further (sup-norm; BCKL2010). The **nerve/Mapper** readout needs its own, and
  **PNV20XX** is it: the persistent nerve lemma gives `Pers(Nerve U) = Pers(W)` functorially for a growing
  **good cover**, and its ε-good generalization **bounds the bottleneck** nerve↔space — so the nerve-over-T
  persistence faithfully reflects the data, with a quantified error when the cover is imperfect. That is the
  ThermoMapper-layer-B prerequisite (the nerve sees the real topology) and the nerve-path twin of DBK2023 (the
  metric path's faithfulness): same chain shape — a Lipschitz bound from construction to truth — applied to the
  nerve instead of the curve. *(PNV20XX now lives in `compendia/mapper` — the Mapper readout's home; cited by key.)*
- **The diagram-space caveat (L5).** You cannot average diagrams (not a vector space), so "posterior over
  barcodes" is reported as "mean ± band over **landscapes**." **FH2024**'s landscapes are the Banach-valued,
  stable, averageable surrogate — and FH2024 is *also* your zigzag readout (§2.3.2) and your bifurcation target
  (§6.2, Hopf), so one summary functor serves both the monotone (β,θ) readout and the Mapper-nerve zigzag. **MMO2019** is
  the *complementary* surrogate: a KDE *density* on diagram space that keeps the random **cardinality** a
  landscape discards, with a mean-absolute-*bottleneck*-deviation dispersion *in the chain's own metric* — landscape-mean for a
  Banach CLT, KDE-density when the *count* of bars is itself the signal (a born-at-transition bar is exactly
  that). And its upper/lower split *is* L4's noise likelihood: long bars tracked individually, near-diagonal
  modeled collectively — "persistence launders noise" stated as a model, not a threshold.
- **Multi-observable (L6) downgrades the metric, not the program.** K monotone observables ⇒ K-parameter
  filtration; **REF-MPH** §9.3 supplies the interleaving/matching-distance robustness, but §9.2.1 ("Barcodes?")
  is the honest guard — *no complete barcode invariant*, so the confidence set lives in the interleaving metric
  and you report computable surrogates (**SGW2025** multigraded Betti; FH2024 multiparameter landscape §2.2.2).
  State that where it bites: a non-monotone observable is a zigzag axis instead (Z5 machinery).
- **The target (MR2026), and its dynamical twin (DLST2026).** The apparatus certifies a *feature*; **MR2026** §4.2
  is the theorem that the feature *is* the phase transition (persistent entropy as the scalar detector). It closes
  the loop to χ(T): the bar you have now certified is the transition χ's peak nominally flags — but with a
  confidence set instead of a heuristic. χ(T) was the first-order parity guarantor; the certified bar is the thing
  it was a proxy for. **DLST2026** (Dey — the Z5 author) is the *same event* read by dynamical systems: a
  **Conley-Morse persistence barcode** from the Conley index over a poset, interval-decomposed via gentle algebras
  and **computed by adapting zigzag** — a bifurcation (the dynamical name for a phase transition) rendered as a
  barcode on the *same machinery* (Z5). Two detectors of one event — persistent entropy (scalar, thermodynamic)
  and the Conley-Morse barcode (structural, dynamical) — a **triangulation** candidate (independent signals
  agreeing), not a redundancy.

## Candidate extensions (ideation — the backbone / slice / operator thread)

Fishable claims from the ordered-backbone discussion; same genre split as fable-BARS-lemma (a positive result
+ one-witness teeth), grouped by which vocabulary coordinate they move. The **slice** trio (SL) keeps the chain
1-parameter; the **zigzag** pair (ZZ) extends it to the non-monotone *structural* readout the Z5 engine actually
computes; the **nerve** pair (NV) carries it to the Mapper construction; the **operator** pair (OP) pushes it past
the barcode — one a surprise in the honest direction; the **triangulation** pair (TR) asks when two detectors
agreeing actually buys confidence. (Loosely: SL/ZZ move the Domain order-type, NV the construction, OP the Degree,
TR the combination — a locating aside, not a hierarchy to enforce.)

- **SL-1 — pullback keeps the barcode** *(positive; the whole warped-curve program rests here)*. A monotone
  path `γ : [0,1] → P` composed with a `P`-indexed module `M` is a 1-parameter module `M∘γ` ⇒ decomposes
  into intervals ⇒ proper barcode. Certifies "a warped slice keeps the barcode exactly where the full
  product-poset (REF-MPH §9.2.1) loses it." **Fishable:** `∘`-functoriality **proved** (trivial); "barcode exists" **cited** — and the cite is now
  *constructive*, not a bare structure theorem. **AL2026** (ph) gives a closed **multiplicity formula** `d_M(V_I)`
  for each interval `I` in a *poset* module via ranks of structure-map matrices (the 1-D barcode generalized),
  plus the **essential cover**: an order-preserving `ζ : Z → P` that *preserves `I`'s multiplicity under
  restriction*, and when `Z` is Dynkin type A (the zigzag/fence poset = a path) the multiplicity reads off the
  **filtration level** directly. That is the warped-curve move as a theorem — the slice recovers a target
  feature's multiplicity, and AL2026 says *which* slice (the essential cover of `I`). **DLST2026** (ph) supplies
  the "barcode exists" half for exactly our setting: the Conley-index poset module is **interval-decomposable
  via gentle algebras**, computed by adapting **zigzag** persistence — a concrete recent grounding stronger than
  the generic Crawley–Boevey pfd structure theorem (REF-PH §3) for the zigzag/poset case the Z5 engine lives in.
- **SL-2 — off-path invisibility** *(one-witness; the cost, with teeth)*. ∃ a 2-parameter `M` and monotone
  `γ` with an interval summand of `M` absent from `Pers(M∘γ)`. Certifies the slice is load-bearing — an
  arbitrary `γ` gives an arbitrary answer. **Fishable:** explicit witness (a feature that opens off-path);
  `argmax_expectation_noncommute` genre. **Empirical witness — SGL2022** (ph): the Nematic XY model
  (Ising + BKT) needs *two* filtrations — one slice catches one transition, neither catches both — i.e. SL-2
  in the wild, and their "design the filtration, not unsupervised" is the off-path cost confirmed on a real
  spin system. So SL-2 graduates from constructed-witness to observed-phenomenon. **AL2026 is the constructive
  dual:** the cost is real for an *arbitrary* `γ`, but the essential cover is the *non-arbitrary* slice — the one
  that provably keeps a *chosen* interval's multiplicity. So SGL2022's empirical "design the filtration, don't go
  unsupervised" gets a theory: design `γ` = the essential cover of the feature you mean to certify, and SL-2's
  teeth bite only the slices that aren't it.
- **SL-3 — slice ≠ reparam** *(one-witness; bounds the bridge's T3)*. ∃ monotone `γ₁,γ₂` through one
  landscape, not reparametrizations, with `Pers(M∘γ₁) ≇ Pers(M∘γ₂)`. The honest dual of
  `barcode_reparam_invariant`: *re-timing one axis is free; choosing the slice is not.* **Fishable:** one-witness.
- **ZZ-1 — zigzag stability for the structural readout** *(positive; nearly an omission, not an extension)*.
  The chain certifies the **scalar** readout — the sublevel barcode of an observable `f(T)` (`χ`, `b₁(T)`),
  monotone, where function-stability (L2) applies directly. But a scalar is a **0-form reduction** (Degree
  grammar): `b₁(T)` is the *count* of bars alive, not the bars. The **native SPC-over-T object is the structural
  zigzag** the **Z5 engine** computes — clusters merge *and* split as `T` varies, and that non-monotonicity is
  *why Z5 exists*. Certifying the scalar certifies a projection; the structural barcode is uncertified. Closing it
  needs **zigzag algebraic stability** (interleaving bounds bottleneck for zigzag modules) in place of
  function-stability — and that leg is **already in the corpus**: **FH2024 §4** (block-extension functor +
  interleaving distance for zigzag modules — and FH2024 is *already* the L5 summary) and **DS2026 §3.1**
  (ZZ-GRIL stability over a quasi-zigzag *bifiltration* — the formal twin of the warped-curve-over-two-axes
  object). **Fishable:** decomposition **cited/built** (CDSM2009 + DLST2026 + AL2026 Dynkin-A gentle algebras;
  Z5 computes it); zigzag stability **cited** (FH2024 §4 / DS2026 §3.1 — Botnan–Lesnick optional, for the
  cleanest general constant). Bridge-shaped — and less of a gap than first flagged.
- **ZZ-2 — the structural readout is more fragile than the scalar** *(one-witness; the SW-jitter worry as a
  theorem)*. ∃ a curve perturbation small in sup-norm that **reorders a merge/split** and moves the *zigzag*
  barcode by `> ε`, while the *scalar* (sublevel) barcode moves `≤ ε`. Certifies the structural readout does
  **not** inherit the scalar's clean L2 — zigzag stability is genuinely weaker, sensitive to event *order*. This
  is "SW's MC jitter manufactures spurious merge/split events" stated precisely. **Mitigation (positive
  companion):** the BARS posterior over `f(T)` (and adaptive / joint-posterior sampling) *averages out* the
  reordering — the certified object is the **posterior zigzag**, not a single jittery draw. Fragility is real
  per-draw; the posterior is the answer — the **L5 landscape-mean argument re-applied to the zigzag**. **Fishable:**
  one-witness for the fragility; the mitigation reuses L5.
- **NV-1 — nerve faithfulness** *(positive; the Mapper readout's L2, promoted from a note)*. **PNV20XX** (now
  `compendia/mapper`): the persistent nerve lemma + ε-good bottleneck bound controls `d_b`(nerve-over-T,
  space-over-T) — the Mapper-nerve sequence faithfully tracks the data's topology. The **nerve construction**
  (cover → nerve) is a stable Transform in the same sense the **sublevel construction** (function → filtration) is.
  **Fishable:** **cited** (PNV20XX); the discharge is showing the SPC cover meets ε-goodness — a graphs/clustering
  obligation, not a topology one.
- **NV-2 — the cover is *random*** *(one-witness; the Mapper-side jitter, twin of ZZ-2)*. The SPC cover is built
  by a **stochastic** clusterer (SW draws), so it is a **random cover**: each realization can be ε-good yet the
  nerve persistence varies across draws. Certifies NV-1's faithfulness is *per-realization*, not automatic for the
  ensemble. **Mitigation:** the *persistent* nerve over `T` + the draw-ensemble posterior — same posterior-averaging
  resolution as ZZ-2. **Fishable:** one-witness (each cover-member ε-good, ensemble nerve unstable).
- **OP-2 — harmonic spectrum = barcode** *(positive; the QW2024 precondition)*. `dim ker Δ_k = b_k`
  (Eckmann / discrete Hodge) ⇒ the persistent Laplacian *contains* the barcode (harmonic) + geometry
  (non-harmonic). **Fishable:** classical (**cite** Eckmann); persistent version is a **build** obligation
  in `TDA.Ph` (T2 tier).
- **OP-1 — non-harmonic spectrum is *not* stable** *(one-witness; the operator extension's boundary)*. ∃ a
  small filtration change with a discontinuous non-harmonic persistent-Laplacian eigenvalue jump ⇒
  `confidence_pushforward` does **not** extend to the geometric spectrum in general; only the harmonic block
  does. **Corrects the earlier "delicate via labeling" gloss** — it genuinely fails (operator dimension
  changes along the filtration; known instability), not a turnkey Weyl. So the geometry rides *only* under an
  added **local spectral-gap hypothesis** — the honest guard to name, à la `hp : p.derivative ≠ 0`.
  **Fishable:** minimal few-simplex witness.
- **TR-1 — independence is the load-bearing axis; coincidence is not the enemy** *(positive-with-a-condition;
  sharpens the DLST2026 ⊕-note)*. Two detectors both peaking at a backbone location **is** genuine corroboration
  **provided they are independent** — coincidence of *independent* signals is evidence, full stop, not circularity.
  The condition lives on **independence** (distinct error sources / the detectors don't factor through one
  filtration), not on the coincidence. Persistent-entropy (MR2026, thermodynamic) and the Conley-Morse barcode
  (DLST2026, dynamical) corroborate to the extent their errors are independent; the *richest* form is
  **relationship-consistency** (agreement across the feature's whole structure, not one location), but independent
  peak-coincidence already counts as real evidence. **Fishable:** joint confidence beats marginal **iff** the
  detectors' confidence events aren't nested — and the positive direction is the *default*, not the exception.
- **TR-2 — the failure is *dependence*, not coincidence** *(one-witness; the teeth, and it audits our own chain)*.
  The degenerate case is two detectors *forced* to agree because one is a **reduction** of the other — `b₁(T)` is a
  0-form face of the H1 barcode (Degree grammar), so their coincidence is mechanical and carries no new information.
  This sharpens the standing rule *triangulation = relationship-consistency of independent signals*: the rider "not
  coincident peaks" means **coincidence alone doesn't *establish* triangulation until independence is checked** — it
  does **not** mean coincident peaks are suspect (independent ones are precisely the evidence). The guard the chain
  must pass: confirm the two detectors don't factor through a common filtration *before* counting their agreement.
  **Fishable:** one-witness (a reduction posing as an independent detector).

**Sequencing.** SL-1 first (cleanest fish; certifies the whole warped-curve move). SL-2 / SL-3 are the
one-witness costs that keep it honest (SL-3 graduates next to T3). **ZZ-1 is the highest-value *new* fish — nearly
an omission, not an extension:** the chain proves the scalar reduction, but the structural zigzag (Z5's output) is
the readout that's actually built, so ZZ-1 closes the gap between what's *proved* and what's *computed*. ZZ-2 / NV-2
are the paired jitter-teeth, and their shared mitigation (posterior-averaging — L5 re-applied) is one argument used
twice. OP-2 then OP-1: OP-2 says the Laplacian contains the barcode; OP-1 is the surprise — the geometry does
**not** come for free, the formal teeth behind "far rung, earns its keep only with a gap hypothesis." TR last —
epistemics guard, cheap, and it audits the chain's own corroboration claims.
