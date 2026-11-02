# Lineage gaps: register

Where the Francom lineage and the literature around it stop short, and what closing each gap would take with this engine. The engine is built for diverse applications, so gaps are stated by problem class. Nothing here is scheduled: an entry starts on the user's go-ahead. Lineage findings that land keep their record in the port hub's [register](../bass-full-port-brief.md#what-the-lineage-taught); extension work keeps its entry in the [extensions brief](extensions-brief.md); the map of the lineage itself is the [lineage review](lineage-review.md).

**How to read it.**
- *Literature status* comes from web searches on 2026-11-01, then from the [chip](../../../lean/Protolemmata/BASS/lean-leads-chip-brief.md) that read most references the same day. Unless a paper is marked as read, it is known from its abstract or a search snippet. "Not found" means those searches did not turn it up, not that it does not exist. The [reading list](#reading-list) says which is which.
- *Statements* are posed in [`lean/Protolemmata/BASS/`](../../../lean/Protolemmata/BASS/), and [`lineage-gap-checks.py`](../../../lean/Protolemmata/BASS/lineage-gap-checks.py) checks them numerically. The self-contained ones are proved in `lean/Enthymemata/` (`BassStructureChain`, `GuardedEstimator`, `BassSensitivity`, `BassIntegrals`, `BassContractionCover`), awaiting semantic review. Those files are lemmas about the mathematics. No theorem about BASS's mixing or contraction is proved.
- *Measurements* cited below are this engine's own, with their reports linked.

| Gap | What the literature leaves | What this engine holds | Kind of work | Statements |
|---|---|---|---|---|
| [G1](#g1-does-the-sampler-reach-the-posterior) | Whether Bayesian MARS's sampler reaches its posterior | F-R1's failures, the flat histogram's a1, enumerated kernels | Measurement; the proof-walk is done | [SC-0, SC-2, SC-4, SC-5](../../../lean/Protolemmata/BASS/structure-chain-lemmas.md) |
| [G2](#g2-posterior-contraction-for-bass) | Any contraction result for Bayesian MARS | The BGP proof chain, in Lean | Proof-walk, then new lemmas | [The walk](../../../lean/Protolemmata/BASS/contraction-walk.md) |
| [G3](#g3-sobol-indices-and-active-subspaces-have-different-standing) | What a posterior's active subspace is worth | The moment service, term derivatives | A lemma, a recovery check, a caveat on E8 | [SG-2, SG-5](../../../lean/Protolemmata/BASS/sensitivity-generator-lemmas.md) |
| [G4](#g4-one-moment-service-across-dependent-inputs) | Closed forms under dependent inputs | `BassMoments`, `IIntervalMoments` | A new measure | [SG-3](../../../lean/Protolemmata/BASS/sensitivity-generator-lemmas.md) |
| [G5](#g5-sensitivity-as-a-markov-generator) | One frame for Sobol, DGSM and active subspaces; a smoothness-matched PCE prior | The mixing layer, graph Laplacians, F-X3's per-part contract | Formalization, a readout, a prior option | [SG-1, SG-4](../../../lean/Protolemmata/BASS/sensitivity-generator-lemmas.md) |
| [G6](#g6-orthogonal-discrepancy-for-cut-calibration) | Cut calibration with a discrepancy | D09's cut posterior, stage 7's D·v discrepancy | A derivation, then a calibration option | Not posed yet |

Two other things sit beside these: a [hypothesis](#h1-the-steep-band-as-a-transition-in-sample-size) to falsify, and [what was checked and is not open](#checked-and-not-open).

---

## G1: Does the sampler reach the posterior?

**What the literature leaves.** The lineage judges fits by prediction and coverage. BART has mixing theory, read: an exponential lower bound for one tree with grow and prune (Ronen et al., arXiv 2210.09352), and for Bayesian CART polynomial upper bounds on connected signals with a superpolynomial lower bound on isolated deep ones (Kim & Roćková, arXiv 2306.00126). Both bottlenecks use the tree's hierarchy. Nothing was found for Bayesian MARS.

**What this engine holds.**
- All 144 of F-R1's untempered fits fail the predictive R̂ ≤ 1.05 gate, at medians 2.14 and 1.74 ([report](deterministic-profile-brief.md)).
- On S10's Friedman fixture (N = 500) the flat histogram finds a1, 4 × x1·x2 + 2 × x3 + x4 + x5, at 66–68%, with ln Ẑ(1) of −841.36 ± 0.17, −841.55 ± 0.07 and −841.14 ± 0.32 over three seeds. That is about 45 nats above the reference arms, whose basins hold at most 1.3% ([docs](../../../docs/Maths/BayesianRegression.md); [report](flat-histogram-warmstart-brief.md)). These are exploration runs, one per seed, on one fixture.
- F-M2 found local modes tens of nats apart, every term load-bearing, and births lost to an Occam cost near ½ ln N ([measurements](local-move-measurements-brief.md)).
- The exact gap on enumerated kernels falls exponentially in N, with all of the growth at modes at the count cap ([report](spectral-gap-brief.md)). Held τ, small capped fixtures.
- The best routes found between Friedman's real modes are not superset routes, and on one route's reading are far too slow to have been taken ([report](path-profile-brief.md#reports)).

**How to close it.**
1. *Exact gaps on enumerated kernels: measured* ([report](spectral-gap-brief.md)). At fixed τ the structure chain is a finite reversible chain (SC-0), so 1 − λ₂ can be computed outright. The fixtures are small, so this shows scaling rather than Friedman's barrier.
   - *The gap.* −ln γ reads linear in N in all 17 replicate cells. Fifteen read so by the registered doubling rule; the two M ≤ 4 cells have four rungs a factor of 4 apart, which that rule cannot read, so theirs come from the fit. From N = 144 the bottleneck's mode is at M_max on all 76 rungs, and through M ≤ 3 no mode with a free slot passes 13.8 nats at any N.
   - *The falsifier that came up.* At M ≤ 4 inside F, depth at one free slot read linear. M ≤ 5 followed that state, with the outcome rule written first: at N = 18,432 it is 2,349.5 nats at 0 free, 801.9 at 1 free and 21.2 at 2 free, gaining about ln 4 per quadrupling at 2 free, as SC-5 gives for k = 2. So it was the cap closing a two-term exchange.
   - *Moves.* The fixed-M family raises γ 33- to 210-fold on the base fixtures. That is read only at N ≤ 24 and the base caps, where none of the three bottleneck modes has the two free slots a two-term exchange needs. So it may be the cap again, and "exchanging terms helps an uncapped fit" is unread.
   - *R̂, against the exact kernel at held τ, S10's lengths.* Of 364 failures, 334 are not mixed, 27 between and 3 mixed. Of 104 passes, 34 are not mixed, 49 between and 21 mixed, so a pass says little. The failures' reliability leans on prevalence: only 24 of 468 runs are mixed, all at N ≤ 36, and 3 of those 24 fail, as do 4 of 12 at ten times the lengths. No false-alarm rate is read past N = 36.
   - *The flat histogram.* With a window declared from the exact density of states, it holds basin TV ≤ 0.021 at every rung, down to γ = 4.9·10⁻¹³¹. That window comes from the enumeration. The pilot route refused 13 of 20 runs, and on the functional fixture misread the top basin's mass (0.58 and 0.89 against 0.718) on 7 and 11 round trips. The error bar is too small on 6 of 27 runs.
2. *The walk of Ronen et al. onto BASS is done, and it changes the expectation* ([SC-4](../../../lean/Protolemmata/BASS/structure-chain-lemmas.md)). The first draft here expected linear barrier growth in N, from load-bearing deaths. Their proof's cut is one edge because a root split changes only through the trivial tree. BASS has no root, and a state can hold a term and its substitute together. So the path S → S ∪ T → T loses at most the Occam term (`collapsedLogMarginal_superset_ge`), not N times a variance. F-M2's measurements fit this: births near −½ ln N, deaths at a median of −234 nats. **The expectation is now polynomial in N, of degree growing with the number of terms exchanged at once,** and superpolynomial only if that number grows with N. That can still be very slow: k = 4 costs about e^{−20} at N = 500. *The exact-gap sweep supports this without testing it where it matters.* With a free slot for every term exchanged, SC-5 bounds the depth along the path, so no ladder on this grid can read linear growth at such a mode; growth where slots look free can only mean the exchange needs more terms than that. So these ladders measure k, with the proposal cost on top. At a closed cap the bound does not apply, and that is where the sweep's growth sits. The free-slot result is close to forced on ladder R: on a fixed grid a useful term gains in proportion to N, so a mode with a slot free and a useful term left does not stay a mode. Friedman's weak directions, where the gain is comparable to ½ ln N, are unread, k = 2 rests on one mode, one data seed and one surface, and k ≥ 3 is unread. The same walk suggests a move exchanging several terms at once, as Kim & Roćková's twiggy proposals do, recorded as a finding.
   - *The path profile, measured* ([report](path-profile-brief.md#reports); Friedman, N = 500, one flat-histogram run, one data seed). The best routes found between modes are not superset routes. They lose a load-bearing term before its replacement pays, and reach 15–26 nats below the chains' states, against 54–114 for the superset route. At their lowest state they hold at most 2 terms beyond a chain's, and at most 4 beyond a1's. k counted by additions moves between exemplars, because swaps between near twins come and go with the knots. Net k, the terms held at the lowest state beyond the start's, stays 2 on every A–a1 row, so read net k. SC-5 bounds the log marginal, not π: in π an unpaid term also costs its measure for as long as it is held, 12.9 nats with one factor and 21.3 with two at N = 500, before a ±25-knot class's credit of at most ln 51 or ln 2,601. So the 3.1 nats first recorded here is only the Occam part, the low end of SC-4's range.
   - *The waiting time on those routes* reads e^{28}–e^{45} sweeps from the chains' states against a 2,000-sweep phase, which by the brief's rule is consistent with the chains staying put. At the costliest edge, 11–19 of 28–45 nats is proposal. This is a weak reading. It is a series over one route, so parallel routes shorten it: crediting every addition's class to all later edges gives e^{11}–e^{35}. Chain 1's state sits 20 nats under its phase's median, which flatters its figures. And it has not been checked against an exact hitting time or a chain's first passage.
   - *k against N* reads flat on carried modes, partly by construction, since carrying an exemplar fixes its term count. At N = 2000, A and B leave for basins with more terms, which is SC-4's "M_n grows with N" and remains unread. So the polynomial expectation is neither broken nor confirmed.
   - *Candidate readings, none run.* (b) The fixed-M family at wider caps, which removes the cap from its 33- to 210-fold gain. (c) The flat histogram's error against its round trips: 27 runs with known truth are already in [`spectral-gap-flat.csv`](spectral-gap-flat.csv), and a1 rests on the same method with 6–21 round trips per chain. (d) The waiting time against exact hitting times on enumerated kernels, and against a long untempered chain's first passage to a1, at about 40,000 sweeps a second a chain. (e) Routes between modes found at N = 2000, not carried.
   - *Not a test.* A ladder that holds a weak term's gain near ln N, with its coefficient shrinking as √(ln N / N), was proposed and withdrawn. Its death costs about ln N by construction, whichever mechanism holds.
   - *A positive theorem to attempt.* The superset route is a canonical path, so an analogue of Yang, Wainwright & Jordan's argument might turn SC-5 from one path's depth into a lower bound on the gap, over a restricted set with the cap at least twice the state size. It would also name where the bound fails: substitutes that do not pay. SC-4 records why their result does not carry as it stands: it rests on asymptotic unimodality, which near-duplicate knots and F-M2's modes do not give. So this has to get past that, and the paper is not on disk. `BassStructureChain` holds the finite-chain Dirichlet machinery it would use. The path profile adds two conditions. D06 closes the superset at a hinge and its mirror, one linear column with the intercept, so the argument needs a redraw edge or a rule excluding such pairs. And the bound has to hold for π, where each term's measure counts, not only for the log marginal.
3. *Persistence* (Woodard, Schmidler & Huber 2009) as the reading of the steep ladder band. It is γ(A, φ) = min{1, φ[A]/π[A]}, and measuring it needs a basin's share per level of a flat-histogram run's ln g, which `ln g` alone does not hold. SC-4 sets out the steps. *The Potts question is resolved:* the toy's fitted N³ at q = 3 is pre-asymptotic. The barrier there is 0.00113 nats per site, about e^{2.2} across the sweep, and an exponential at that rate overtakes N³ only near N ≈ 27,000. The theorem itself is Bhatnagar & Randall at q = 3, α = α_c, widened to q ≥ 3, α ≥ α_c by Woodard et al., whose proofs of the propositions are in Woodard's thesis.

**Would land in:** the mixing layer shared across engines; F-M2's record. *Unlikely to help:* locally balanced proposals, an idea that went in here and came out. F-M2 measured every term load-bearing, and its use-weighted births gained only ×1.07, so there is little sign of an attractive neighbour to scan for. That is a reading of those measurements, not a test of these proposals.

---

## G2: Posterior contraction for BASS

**What the literature leaves.** Nothing found for Bayesian MARS. *MARS via LASSO* (Ki, Fang & Guntuboyina, *AoS* 2024, read) gives n^{−4/5} up to log factors whose power depends on the interaction order, for least squares under a variation constraint. Its constants depend on d, so it is a fixed-d result, and "logarithmic in dimension" does not mean logarithmic in p. It states no M-term approximation rate.

**What this engine holds.** [PosteriorContraction](../../../docs/Maths/PosteriorContraction.md): Gaussian KL, tests and covers in Lean, with the general theorem cited as axioms.

**How to close it.** [The walk](../../../lean/Protolemmata/BASS/contraction-walk.md), with its reading of each step, checked against Ghosal & van der Vaart (2017), Chapter 8.
- The Gaussian pieces transport, with one small lemma: the KL ball's inclusion is packaged through a sup-norm slab, and BASS's step works in ‖·‖ₙ.
- Heavy-tailed τ costs nothing in entropy, unlike the GP case. **Proved:** `exists_dictionary_local_cover` covers the union of subspaces at radius ξr/4 with a count free of r, of any coefficient bound and of any range for σ.
- Theorem 8.19(b) as printed asks only for local entropy, and the book names the global cover as a stronger condition. So the fixed-constant axiom can be restated locally and stay a consequence of the book. Proposed in the walk, not made.
- Prior mass needs an M-term rate. **Proved:** `exists_mTerm_approx` gives V²/M in ‖·‖ₙ for a truth in *MARS via LASSO*'s class, so ε_n² ≈ V·(log(pN)/N)^{1/2} for every interaction order, with constants free of p. A fast rate matching n^{−4/5} holds by hand at K = 1 and is open at K ≥ 2.
- "Fixed g" has to mean a sequence of order N: a constant g makes the coefficient mass cost order N.
- The structure prior's normalization is no obstacle: it lies between 1 and ((h₂ + 1)/h₂)^{h₁}, under one nat, free of N and the design.
- The work left sits in the data-dependent priors. σ²'s heavy upper tail looks removable by one test rather than a sieve; the level-centred intercept is as before.

Start at g_n of order N, compact σ and the Theorem-1 route.

**Would land in:** docs/Maths beside the BGP result; `lean/`. **Pays off in:** G3's Sobol reading, which inherits any L²(μ) rate.

---

## G3: Sobol indices and active subspaces have different standing

**What the literature leaves.** Read: Yoo & Ghosal (arXiv 1411.6716) get derivative rates for tensor-product B-splines, on fixed quasi-uniform knots with the number of splines set by a known smoothness, so without adaptation. Odin, Bachoc & Lagnoux (arXiv 2403.03540) recover a subspace for a prior that holds it as a parameter, with the intrinsic dimension known, through a detectability assumption from L² and not through gradients. Nothing was found for gradient outer products under an adaptive hinge posterior.

**The finding (SG-2, SG-5).**
- √V_u is 1-Lipschitz in L²(μ), so a rate for f passes straight to Sobol (**proved:** `sqrt_compVar_sub_le_form`). C = E[∇f∇fᵀ] has no such bound: a small triangle wave with many teeth is small in L² and large in gradient.
- E[C]'s top eigenspace weights draws by their eigenvalues, so it is not the mean of the draws' subspaces (**proved** on explicit matrices: `top_eigenspace_ne_subspace_mean`). They sit 54° apart in a three-draw example.
- Under `EmpiricalMeasure` the atoms sit on knots, so C reads F-X3's kink convention.

**How to close it.** Before E8's directions reduce anything, add a recovery check for C on functions where it is known, at growing N. Offer E[C] and the subspace mean as separate readouts, and say which measure a C was computed under.

**Would land in:** [E8](extensions-brief.md#active-subspace-analysis).

---

## G4: One moment service across dependent inputs

**What the literature leaves.** Shapley effects are the usual reading under dependence (Iooss & Prieur, arXiv 1707.01334), with analytic forms in Gaussian linear cases and otherwise a double Monte Carlo loop. Nothing was found on mixtures of products here.

**The finding (SG-3).** A joint mixture of products keeps every pair integral closed. One component is `ProductMeasure`; N components with delta marginals are a joint empirical sample. So this replaces a split between a closed-form tier and a Monte Carlo tier with one formula, and latent-class inputs, product-kernel smoothing and possibly Bernstein copulas sit between the two.

**Limits.** Hoeffding does not apply under dependence, so what is exact is the moments and C. Shapley drops from a double loop to a single one, not to a closed form. Separated components weaken the Poincaré constant behind active-subspace bounds. Cost is linear in the components, and the delta case keeps its cheaper direct path.

**Would land in:** `InputMeasures`, `BassMoments`, `SensitivitySpec`; its extension entry is [E25](extensions-brief.md#sensitivity-under-dependent-inputs). It also widens the [guarded estimators](../../../lean/Protolemmata/BASS/guarded-estimator-lemmas.md)' exact case.

---

## G5: Sensitivity as a Markov generator

**What the literature leaves.** Read, and narrower than first posed. Roustant et al. (arXiv 1906.09883) state the commuting-projector expansion that is SG-1(b), build chaos bases from each variable's Poincaré operator, and have the weighted eigen-identity that is SG-4(b). What they do not have is a resampling generator over subsets, or a coefficient prior from the eigenvalues. Σ_v total(v) = Σ_u |u|·V_u is, from memory, mean dimension and total influence, so most of this is restatement. The readout (variance by order from U_ρ's factors) and the prior are the likely new parts. R2D2 on polynomial chaos (Bürkner et al., arXiv 2204.06043, read) is exchangeable across polynomials, and degree enters only through truncation.

**The finding (SG-1, SG-4).** Proved in `BassSensitivity` and `BassIntegrals`, over any commuting family of self-adjoint idempotents for SG-1.
- Hoeffding components are eigenfunctions of resampling dynamics, with eigenvalue −|u|. So Σ_v total(v) = Σ_u |u|·V_u is a Dirichlet energy and Efron–Stein a spectral gap of 1.
- Variance by interaction order comes from one polynomial per term pair, from factors `BassMoments` already caches.
- DGSM and active subspaces are the same construction under the Langevin generator, and a generator per variable covers categorical and ordinal inputs.
- Under dependence the spectrum is lost, but the Dirichlet form remains, with Gibbs dynamics' gap measuring how far totals undercount.
- Shifted Legendre is the Jacobi generator's eigenbasis, so a weighted DGSM comes straight from KHAOS's coefficients, and a prior can use the eigenvalues: (5,1) and (3,3) both score 13 under KHAOS's penalty, and 32 and 24 here.

**Connections.** This is the same language as the MCMC family's spectral gaps and the graph Laplacians, which is why SC-2 and SG-1 share a Dirichlet form.

**How to close it.** Build the variance-by-order readout, and test the Efron–Stein identity against `SobolDecomposition` as a property check. The prior is an option beside KHAOS's, to earn its place predictively.

**Would land in:** the sensitivity stage; [E5](extensions-brief.md#khaoss-modified-g-prior).

---

## G6: Orthogonal discrepancy for cut calibration

**What the literature leaves.** Gopalan & Francom's note says the cut approach still needs hybridizing with discrepancy modelling. Plumlee (2017) keeps θ identified by making the discrepancy orthogonal to the emulator's gradient, through GP kernels; Xie & Xu's projected calibration avoids the GP, and arXiv 2605.30492 treats constraints too. Nothing was found on a basis discrepancy with a linear constraint.

**The candidate.** Stage 7's discrepancy is δ = D·v with v Gaussian. Orthogonality to ∂_θ f at θ is a linear constraint on v, C(θ)·v = 0, whose entries come from the moment service with `FactorDerivative` when D's columns are piecewise polynomial, or a finite sum on a grid. Conditioning a Gaussian on a linear constraint stays conjugate, so the collapsed likelihood becomes a low-rank update at each θ.

**Open, and not small.** The constraint moves with θ, and a hinge emulator's ∂_θ f is piecewise constant, so the constrained prior jumps as θ crosses a knot. What that does to θ's posterior, and to the cut's modular reading, has not been derived. So this one is a derivation to attempt before anything is posed for Lean.

**Would land in:** D09's options; E14's neighbourhood.

---

## H1: The steep band as a transition in sample size

**Hypothesis, untested.** Under the integrated form, a power posterior at b has the likelihood scale and the Occam term of a full posterior on bN cases. So Friedman's steep band at b ≈ 0.16–0.32 ([readings](friedman-band-brief.md)) predicts that fits on nested subsamples of about 80–160 of the 500 cases switch basin. If they do, the band is a Bayesian phase transition in sample size, in singular learning theory's sense, and the ladder reads where more data changes the model's structure.

**Falsifying it.** Flat-histogram fits on nested subsamples, read by basin. The linear form's Occam term scales differently, so the prediction is for the integrated form only.

---

## Checked and not open

- **Multi-fidelity BASS exists in the lineage's orbit.** Tucker, Collins, Huerta & Brown (arXiv 2609.27174, read to the end of its introduction) use BASS emulators for multi-fidelity mapping and fusion, with elastic alignment. An AR1 model through F-X1's offset with a stage-8 cut may still differ; compare after reading it.
- **Closed-form Sobol and Shapley for BART** are in Statistica Sinica SS-2024-0318, so that is not a gap on the tree side.
- **A variance-share prior on polynomial chaos** is R2D2 (Bürkner et al.). The open piece in G5 is the smoothness-matched one.

## Manuscript readings to confirm

KHAOS's printed birth, death and Mutate1 ratios invert some proposal factors, beyond the double count the hub records ([SC-1](../../../lean/Protolemmata/BASS/structure-chain-lemmas.md)). That is a reading of the manuscript only; the `khaos` package's source and the supplement have not been read. A KHAOS port derives its ratios and tests their reciprocity either way, under F-X4's contract.

## Reading list

**On disk, read 2026-11-01 for the [chip](../../../lean/Protolemmata/BASS/lean-leads-chip-brief.md)'s questions** (the statements and proofs the walks use, not every section).
- In `gauntlet/bayesian`: Ronen et al. (2210.09352v1); Kim & Roćková (2306.00126v1); Odin, Bachoc & Lagnoux (2403.03540v2); Roustant et al. (1906.09883v1; 1612.03689v1 by search and its introduction); Bürkner et al. (2204.06043v2).
- In `gauntlet/stats`: Ki, Fang & Guntuboyina, *MARS via LASSO* (2111.11694v5).
- In `PDFs/Bayesian`: Yoo & Ghosal (`Ghosal2015-1411.6716v3.pdf`); Woodard, Schmidler & Huber (`Woodard2009.pdf`).
- Ghosal & van der Vaart (2017), Chapter 8, Sections 8.3–8.3.2.

**On disk, unread.** Tucker et al. (2609.27174v1, read to the end of its introduction); Plumlee (`Plumlee2016.pdf`).

**Not found on disk by name or arXiv id.** Yang, Wainwright & Jordan (2016), known here through Kim & Roćková; Liu & Owen (2006) on mean dimension; Kuo et al. (2010), cited by Roustant et al.; Temlyakov (2018) and Dũng, Temlyakov & Ullrich (2018) on mixed smoothness, cited by *MARS via LASSO*. Bhatnagar & Randall (2004); Woodard et al.'s AAP companion, if `Woodard2009.pdf` is the EJP paper; Shen & Ghosal (1403.0625); Xie & Xu (1803.01231); arXiv 2605.30492; Lamboni et al. (1202.0943); arXiv 2412.04918 and 2605.30979; Iooss & Prieur (1707.01334); O'Donnell's *Analysis of Boolean Functions*; Sancetta & Satchell (2004) on Bernstein copulas; Zrnic & Candès (2024); Angelopoulos, Duchi & Zrnic (2023); Giles (2008); Peherstorfer, Willcox & Gunzburger (2016, 2018).
