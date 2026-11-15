# Conformal validity — prediction sets valid whatever the model does

Pins the contract of the engine's conformal layer, one layer over retained draws and pointwise densities for every model with a pointwise density ([E16](../../../issues/regression/bass/extensions-brief.md#conformal-bayes-the-consumer)). The layer owns coverage and the score owns width. These are the finite statements of that split. They are the sibling of [GuardedEstimator](../BASS/guarded-estimator-lemmas.md): there, means valid whatever the emulator does; here, prediction sets valid whatever the model does. Written from the [chip brief](conformal-chip-brief.md).

**Status (2026-11-15).** Proved in `lean/Enthymemata/ConformalValidity.lean`, namespace `Conformal`, with no apologies (PROOF-CLOSED, 668 declarations counting E35's M-1 and M-2 and slices 9, 10, 12, 13, 14, 15, 16 and 18, slice 14's follow-ups 1, 2 and 3, slice 18's follow-up, RoBAS-Full's intervals, and E35's CUSUM and Shiryaev–Roberts guarantees at the end of this page; `#print axioms` shows `propext`, `Classical.choice`, `Quot.sound` only):
- CB-1: `card_rankLe_le`, `coverage_orderings`, `coverage_orderings_eq`, `coverage_strict_le`, `coverage_orderings_le`, `rank_uniform`, `coverage_tieBreak`, `coverage_given`, `coverage_law`, `coverage_stratum`;
- CB-2: `split_coverage`, `exchangeable_given_training`, `cqr_ordConnected`;
- CB-4: `sum_jointU`, `post_addOneIn`, `integratedScore_equivariant`, `latentScore_equivariant`, `latentScore_eq`, `sum_jointU_addOneIn`, `integratedScore_le_latentScore`, `accepts_latent_of_mixed`, `mixed_counterexample`;
- CB-3, with an error bound per unit: `accepts_or_fragile`, `accepts_shift_of_accepts_close`, `acceptsInflated_of_accepts`, `fragileEstimated_of_fragile`, `fragileEstimated_mono`, `coverage_close`, `coverage_inflated`, `card_fragile_le_estimated`.

[`conformal-checks.py`](conformal-checks.py) holds the exact-arithmetic checks, run before the Lean. Every assertion holds. Changed from the brief: CB-1's upper bound is an equality at the grid levels, and Hoff's strict rule fails under ties; CB-3's band is one-sided and tight, its bounds are per unit, and it gives a third option, the inflated rule. Nothing proved false. CB-3's prior art was read in full on 2026-11-06; the inflated rule is known, and the margin was not found. Open: the law forms on continuous spaces. Closed on 2026-11-10, at the end of this page: CB-3's margin under an ε that holds only with high probability (slice 10), and the Hoff finding under CB-4 (slice 12).

**Conventions, fixed here and used throughout.**
- A score is a **conformity** score: higher conforms better. The engine's pointwise log density is one, and a monotone map of every unit's score changes nothing.
- `C` is the ranked set: every unit, or the test unit's stratum. It has `m = |C|` units, which is n + 1 for n data units and the candidate.
- A unit's rank `rankLe C s u` counts the units of `C` scoring at most its score, itself included.
- **The accept rule at level k** accepts u when `k < rankLe C s u`. This is Fong and Holmes's π(y) > k/(n + 1), with ties counted for the candidate.
- The strict rule accepts u when at least k units score strictly below it. This is Hoff's "c_{n+1} above the k-th order statistic", which he states without ties.
- k is the level. E16's "steps of 1/(k + 1)" count replicates, which here are steps of 1/m.

**What the finite statements cover.** Exchangeability is written as the average over the orderings of a fixed data set (`permsOn C`, the permutations moving only units of `C`). An exchangeable law on a finite space is a mixture of those averages (`sum_law_average`). The law forms need a finite data space. On a continuous one, the same argument runs by conditioning on the multiset, and that is not formalized. The scores are any function of the data set, so a "fit" includes whatever a sampler does with its seeds.

---

## CB-1 — the rank guarantee

**Claim.** Take a data set z of units and a score equivariant on C, so reordering the units of C reorders their scores (`Equivariant`). Over the orderings of z, the test unit is accepted at level k with frequency:
- at least 1 − k/m, whatever the ties (`coverage_orderings`);
- exactly 1 − k/m when the scores are distinct on C (`coverage_orderings_eq`);
- at most 1 − k/m plus the share of units whose ties straddle the threshold, rankLt < k < rankLe (`coverage_orderings_le`).

The strict rule covers at most 1 − k/m (`coverage_strict_le`). With distinct scores the test unit holds each rank 1, …, m on 1/m of the orderings (`rank_uniform`).

**The counting core** (`card_rankLe_le`): whatever the scores, at most k units of C have rank at most k. Take the highest-scoring of them, and every one of them counts toward its rank. Its dual gives the strict rule's bound. With no ties the two rules agree (`accepts_iff_acceptsStrict`).

**The forms.**
- *Random tie-breaking* (`coverage_tieBreak`): a label carried with each unit, distinct within C, breaks ties lexicographically. The score stays equivariant and has no ties, so coverage is exactly 1 − k/m. These are the smoothed p-values of the literature.
- *Under a law* (`coverage_law`): `W` carries anything else the score reads, with any law. Given each w the units of C are exchangeable. `coverage_given` is the same guarantee given w.
- *Per stratum* (`coverage_stratum`): the law and the score are invariant under every permutation that keeps units in their strata. Ranked within the test unit's stratum, coverage there is at least 1 − k/m_stratum. The general `coverage_law` with C the stratum needs less: only permutations inside that stratum. That is Fong and Holmes's form, where the score may read the other strata in any way.

**Changed from the brief.** The brief's "at most 1 − k/(n + 1) + 1/(n + 1)" is Lei et al.'s bound for a level α off the grid, with k = ⌊α(n + 1)⌋. At a grid level, with distinct or tie-broken scores, coverage is exactly 1 − k/(n + 1). With ties, the liberal rule can reach 1: all scores tied accepts every candidate. Hoff's strict rule is not valid under ties: all scores tied refuses every candidate. So the tie rule is not cosmetic.

**Checks.**
- The counting core holds on 2,004 (scores, k) cases with ties, up to m = 5.
- Over every ordering of 156 (score, multiset) pairs (three equivariant scores, alphabet 3, m ≤ 5), liberal ≥ 1 − k/m ≥ strict, with equality for distinct and tie-broken scores.
- A score that reads positions falls 3/4 below the bound.
- Under an iid law and a Pólya urn at m = 4, k = 1, coverage is 0.85 to 0.91 against 3/4. The position-reading score covers 0.40 and 0.33.
- *Per stratum,* a stratum of 2 beside one of 3 with a different law: ranked within the stratum, 0.851 against 1/2. Ranked across all five units at the same nominal level, 0.472 against 0.6.
- *Small strata* at α = 0.1: k = ⌊α m⌋ is 0 until m = 10, and the set is everything; then 9/10, 10/11, …, 18/19, 9/10, 19/21.

**Decision: the layer's tie rule, the C# tests' band, and small strata.**
- *Tie rule:* count ties for the candidate (rank ≤, the candidate included). It is valid whatever the ties. Offer random tie-breaking where an exactly uniform p-value is wanted (a calibrated PIT, the conformal predictive distribution). Never the strict rule.
- *C# tests on enumerable fixtures* (exact enumeration of orderings):
  - with distinct scores, assert coverage = 1 − k/m exactly;
  - with ties, assert 1 − k/m ≤ coverage ≤ 1 − k/m + (straddling units)/m, and the strict rule ≤ 1 − k/m;
  - with tie-breaking labels, assert exactly 1 − k/m.
- *Small strata:* a stratum of m units attains only 1 − k/m. At level α the layer takes k = ⌊α m⌋, and coverage exceeds 1 − α by less than 1/m. Below m = ⌈1/α⌉ the set is the whole space, and the layer reports that, not a set. Ranking is always within the stratum.

**Prior art.**
- Read: Fong and Holmes (§2.1 and their Proposition for groups), citing Vovk, Gammerman and Shafer (2005, ch. 8.7) for validity and Lei et al. (2018, Thm 1) for the 1/(n + 1) upper bound with continuous scores. Hoff 2023, §4.1: the strict rule, "if … there are no ties".
- Smoothed p-values are standard (Vovk, Gammerman and Shafer 2005), cited from memory and not re-read.

---

## CB-2 — the split form holds for any fit

**Claim.** Let `W` be the training split, with anything else the fit reads (the sampler's seeds, an unmixed chain). Given w, let the calibration and test units be exchangeable. Then for any map `fit` from W to a score, the score `fit w (z i)` is equivariant (`splitScore_equivariant`, by `rfl`), and coverage is at least 1 − k/m (`split_coverage`), given w too (`coverage_given`). The hypothesis follows from exchangeability of all the units, training ones included (`exchangeable_given_training`).

**The instance.** CQR's conformity, min(y − q_lo(x), q_hi(x) − y), is the negative of Romano, Patterson and Candès's score. Its accept set in y is an interval (`cqr_ordConnected`). The checks find it to be [q_lo + c₍ₖ₎, q_hi − c₍ₖ₎], with c₍ₖ₎ the k-th smallest calibration conformity: 3,072 sets on a grid, every one.

**Checks.** Four fits, under an iid law and a Pólya urn: a fitted centre, an "unmixed sampler" that reads one training unit, a fixed wrong centre, and a seeded random lookup. With two training units and three calibration and test units, coverage at k = 1 is 0.83 to 0.89 against 2/3, and the worst given any training split is 0.80. CQR with quantile fits that miss covers 27/32 against 3/4.

**Decision: the split form as the guarantee that does not lean on the sampler, cited as a theorem.** `split_coverage` has no hypothesis on the fit. A chain that never mixed, an approximate posterior, or a point estimate all give a valid split set, and only the width suffers. The layer cites it for any request that must not lean on the sampler, and states it given the training split.

**Prior art.**
- Read: split conformal (Lei et al. 2018, as Fong and Holmes summarize it), and their note that hyperparameters chosen on Z₁…Zₙ break full conformal but not split.
- From the abstract: CQR (Romano, Patterson and Candès, arXiv 1905.03222).

---

## CB-4 — latents keep the score symmetric only if every unit treats them alike

**Setting** (finite). θ ∈ Θ with prior π, each unit's latent v ∈ V with law q(v | θ), the density f(z | θ, v), and the integrated density g(z | θ) = Σᵥ q f (`integrated`).

**(a) Add-one-in.** The joint posterior's θ-marginal is the posterior with each latent integrated (`sum_jointU`). The candidate enters as one more factor (`postU_addOneIn`), and normalized it gives Fong and Holmes's self-normalized weights (`post_addOneIn`).

**(b) The integrated score** sᵢ = Σ_θ p(θ | all units) g(zᵢ | θ) is equivariant (`integratedScore_equivariant`).

**(c) The latent-posterior score** sᵢ = E[f(zᵢ | θ, vᵢ)], under the joint posterior of θ and all n + 1 latents, is equivariant too (`latentScore_equivariant`). Its closed form is Σ_θ p(θ | z) Σᵥ q f²/g (`latentScore_eq`). From the n-unit draws (`sum_jointU_addOneIn`), the joint with the candidate is the n-unit joint times q(v′ | θ) f(z_new | θ, v′). That factor is the add-one-in weight g(z_new | θ) times the candidate latent's conditional p(v′ | θ, z_new). So the data units read their retained latents with the weight, and the candidate needs a latent drawn from its conditional for each draw.

**(d) The mixed score breaks symmetry, and it undercovers.** Score the data units by (c) and the candidate by g.
- Unit by unit, (c) ≥ (b) (`integratedScore_le_latentScore`). The gap is the latent's conditional variance of f over g, zero only where f does not vary with the latent.
- So the mixed score accepts no candidate that (c) refuses (`accepts_latent_of_mixed`).
- `mixed_counterexample` has one θ, a fair latent bit, y equal to the bit with probability 3/4, and two units. The data unit scores 5/8 and the candidate 1/2, whatever the data, so the candidate is never accepted at k = 1. Coverage is 0 against 1/2. Both symmetric scores cover 1 there, all tied.
- In the checks' model (two θ, a latent that moves y strongly, x and y), at m = 4, under the model's own law and under an iid law outside it:

| score | k = 1 (bound 0.75) | k = 2 (bound 0.50) |
|---|---|---|
| (b) | 0.808, 0.812 | 0.594, 0.578 |
| (c) | 0.800, 0.804 | 0.604, 0.585 |
| (d) | 0.532, 0.510 | 0.263, 0.224 |
| (d), data latents not reweighted | 0.540, 0.515 | 0.251, 0.202 |
| reversed (data integrated, candidate by its latent) | 0.892, 0.895 | 0.709, 0.728 |

  The reversed mix over-covers, so it is valid but loses the exact band.

**Decision: which density the conformal layer reads for every latent-noise model** (Student t, quantile, Normal-Wald, generalized beta prime, Pólya–Gamma later). Every unit, the candidate included, is read the same way.
- *(b), the default:* the density with every observation-level latent integrated. This is OQ17's scope rule and GBASS's pointwise content now. Its weight is the same add-one-in g, and it needs no per-candidate latent draw.
- *(c), by request:* the latent-posterior density for every unit, with the candidate's latent drawn from p(v | θ, y) for each draw and each candidate.
- *The contract:* the layer refuses a pointwise matrix whose latents are conditioned on, the retained content of a latent-noise model. The retained draws keep their latents, but the conformal reader does not read them for the data units alone.

**Finding (Hoff, out of scope, recorded).** Hoff's Bayes-optimal score is cᵢ = p(zᵢ | every other unit, the candidate included), a leave-one-out predictive. Fong and Holmes's add-one-in score is the in-sample p(zᵢ | all units). Both are equivariant, so both cover. They rank differently: on 12 of 216 data sets of three units in the checks' model, and in Hoff's own setting with no covariate, y = (2, 2, 2, 0) puts the lone 0 lowest in-sample and highest by Hoff's. So E16's line that Hoff's score "is the add-one-in step above", and hence that the step's score is the Bayes-optimal one, is not right as written. From draws, Hoff's cᵢ needs a weight of g(z_new | θ)/g(zᵢ | θ) for each data unit (add-one-in and leave-one-out at once, with PSIS's instability), and its c_{n+1} is the plain n-unit predictive. The bullet sits outside "What the conformal literature adds", and it was corrected there on 2026-11-06 at the user's call.

**Prior art.**
- Read: Fong and Holmes, §2.2, the score's invariance from π(θ | Z₁…Zₙ₊₁) ∝ π(θ) ∏ f, and their limitation that MC samples must not "disrupt exchangeability". Hoff 2023, §4.1, for cᵢ.
- Nothing on mixing latent readings across units was found. No search went beyond those texts.

---

## CB-3 — what an imperfect score costs

**Bounds.** Each unit j carries its own error bound εⱼ, its Monte Carlo error for instance. With every εⱼ equal to ε, the shifts below are 2ε and 4ε, which is how the brief posed it. This was generalized on 2026-11-06, following Ndiaye's per-unit τᵢ.

**Claim, the deterministic part.** Let each ŝⱼ be within εⱼ of an exact score sⱼ on C. Then:
- ŝ refuses a unit u that s accepts only when u is **fragile** (`accepts_or_fragile`). Fragile means accepted by s, but with no more than k units j of C at least εⱼ + ε_u below it, itself counted (`Fragile`).
- ŝ accepts a unit that s refuses only when more than k units have sⱼ ≤ s_u + εⱼ + ε_u, so that u is within reach of acceptance (`accepts_shift_of_accepts_close`, Ndiaye's deflated side).
- So a decision changes only inside a band around the threshold. With uniform ε and distinct scores, the loss side of that band is M = {u : s₍ₖ₊₁₎ ≤ s_u < s₍ₖ₎ + 2ε}: the accepted units within 2ε of the highest refused score.

**Claim, over orderings.** Let S be equivariant, and let the estimate ŝ_σ for each ordering σ be anything within its bounds (Monte Carlo error, an asymmetric computation).
- *The margin* (`coverage_close`): with bounds that travel with the units, coverage is at least 1 − (k + |M|)/m, M the exact score's fragile units. This is the brief's identity P(test ∈ band) = E|band|/(n + 1), with the one-sided band.
- *The inflated rule* (`coverage_inflated`, from `acceptsInflated_of_accepts`): accept u when more than k units have ŝⱼ − εⱼ ≤ ŝ_u + ε_u. It keeps 1 − k/m. Its bounds may differ by ordering and by position, since the candidate's error differs from a data unit's.

**The engine's reading.** |M| is at most the count of the band read from any ordering's estimates and bounds alone (`card_fragile_le_estimated`, from `fragileEstimated_of_fragile`). That band holds the units the inflated rule accepts with at most k units 2εⱼ + 2ε_u below them. It only widens as the bounds grow (`fragileEstimated_mono`), so any dominating bounds will do. With uniform ε and distinct estimates it is ŝ₍ₖ₊₁₎ − 2ε ≤ ŝ_u < ŝ₍ₖ₎ + 4ε, the brief's "within 4ε" on the high side. So the layer can report 1 − (k + m̂)/m from what it holds.

**Where per-unit bounds pay.** Monte Carlo error attaches to a position, the candidate's or a data unit's, more than to a unit.
- The inflated rule takes bounds by position, so it gains directly: the candidate's larger error no longer widens every comparison.
- The margin needs bounds carried by the units. Take each unit's bound as its worst error over the orderings: every unit is the candidate in some ordering, so it inherits the candidate's error. In the checks' conformal-Bayes case that bound equals the uniform one.
- So the margin gains only where units differ in themselves, a unit in a tail for instance.

**Changed from the brief.** The brief's band is two-sided, |s_test − s₍ₖ₎| ≤ 2ε. The coverage loss needs only the side above the threshold, and the bound is tight: an adversary within the bounds, the test unit down and the rest up, refuses exactly the units the bound counts. The bounds are per unit. The k-th order statistic moving at most ε is checked for uniform ε, not formalized; the count form replaces it.

**Checks.**
- The deterministic part holds on 150,435 cases: exhaustive on grids (uniform ε to m = 4, per-unit ε to m = 3) and over a seeded random sweep. A loss lies only in M, the inflated rule keeps every exact acceptance, M lies inside the estimated band, and with uniform ε order statistics move at most ε.
- With uniform or per-unit bounds, the adversary attains the bound for every unit.
- Over orderings, the inflated rule holds with bounds by unit and by position. With bounds carried by the units, the margin holds for an adversarial and a noisy estimate, the adversary attains it with distinct scores, and every ordering's estimated band counts at least |M|.
- Conformal Bayes by add-one-in weights over T = 20 and 200 draws on a Beta–Bernoulli grid: ε is 0.007 to 0.073, and every bound holds, per unit and by position too. Here the draws' decisions match the exact ones, and each unit's worst error equals the uniform ε.

**Decision: whether the full-conformal layer reports its margin, and when it falls back to the split form.** One lemma gives these policies.
- *Report the margin:* coverage ≥ 1 − (k + m̂)/m, with m̂ the estimated band's count, read with bounds carried by the units and taken as the maximum over the candidates on the grid.
- *Inflate by each position's bound* (Ndiaye's stable sets): keep 1 − k/m at the cost of width. The width follows each estimate's own error, the candidate's included.
- *Fall back to the split form (CB-2):* when m̂/m is too large for the requested level, or when the bounds cannot be trusted.
- *Correct by tournament* (Bhattacharyya, Zhang and Barber, under prior art): full data and 1 − 2α whatever the draws, at the cost of n leave-one-out chains.

The bounds are sure in these statements. A Monte Carlo bound holds only with some probability, and its failure probability adds to the miscoverage. The bounds come from the sampler's own error estimate for each add-one-in score, so they are only as good as the mixing. A chain that fails its diagnostics sends the layer to the split form or the tournament.

**Prior art** (read on 2026-11-06 at their theorems, from `gauntlet/conformal-bayes` and `PDFs/Bayesian/ndiaye22a.pdf`).
- *The inflated rule is known.*
  - Ndiaye, *Stable conformal prediction sets* (ICML 2022, PMLR 162; Prop. 3.2 to 3.4). Per-unit sure stability bounds τᵢ sandwich the exact p-value between a deflated and an inflated one. The inflated set contains the exact one and covers 1 − α. `coverage_inflated` is its case with a uniform ε, and `accepts_shift_of_accepts_close` is its deflated side.
  - Razafindrakoto, Celisse and Lacaille, *Approximate full conformal prediction in an RKHS* (arXiv 2601.13102, the approximate region's definition). The same construction, with a slack τᵢ(y) for each unit and candidate.
  - Bhattacharyya, Zhang and Barber, *Approximating full conformal prediction: … the tournament correction* (arXiv 2605.29200, appendix, "Coverage of the approximate method under stronger stability"). They inflate by 2ε, where each unit's error, the test unit's included, is within ε with probability 1 − ν, and they get 1 − α − 3√ν. Their Example 3 is Fong and Holmes's add-one-in score from draws.
- *The uninflated rule's margin was not found.* None of the five bounds an uninflated approximation's coverage by a fragile count. Bhattacharyya et al. write that such approximations "typically do not come with general finite-sample coverage guarantees". Abad Martinez et al. (AAAI 2023, arXiv 2202.01315) bound their influence-function scores' error by O(1/N), but leave finite-sample validity open. `coverage_close` answers that for any sure score error.
- *A full-data guarantee that does not need mixing* (Bhattacharyya et al., the tournament correction, Thm "Validity" and its randomized form).
  - The candidate is compared with each data unit i, using a score fitted without both of them and symmetric in the rest.
  - With no assumption on the approximation, coverage is at least 1 − 2α.
  - Under a weaker stability, coverage of the 2ε-inflated set is at least 1 − α − 3√ν. The stability needed is that moving one data unit's role leaves the candidate's score within ε.
  - For conformal Bayes this takes posterior draws given each leave-one-out set: n chains, each needing to treat its data symmetrically but not to mix.
  - Their shortcut, which rejection-samples from the n-unit draws, is valid only for exact posterior draws and a likelihood bounded below.
- Barber, Candès, Ramdas and Tibshirani, *Conformal prediction beyond exchangeability* (AoS 2023, arXiv 2202.13415, "Nonexchangeable full conformal prediction" with §3.2).
  - An algorithm that treats one position specially is made valid by swapping the test unit with a random unit before fitting. The gap term vanishes under exchangeability.
  - For add-one-in that means draws given every unit but a random one, the candidate included: a refit per candidate, which the step exists to avoid.

---

## E35 — conformal test martingales: exchangeability checked online

Posed and proved on 2026-11-08 for [E35](../../../issues/regression/bass/extensions-brief.md#conformal-test-martingales-exchangeability-checked-online-e35), before any C#. As with CB-1 to CB-4, the point is the contract: the guarantee is the product's, and its hypotheses are easy to break silently.

**Status (2026-11-08).** Both statements are proved in `lean/Enthymemata/ConformalValidity.lean`, with no apologies: M-1's finite core and its law form, and M-2, the link and Ville. Two parts are not formalized, and both are checked exactly in [`conformal-checks.py`](conformal-checks.py): M-1's lift from ranks to continuous p-values, and the e-value at an unbounded stopping time. CUSUM and Shiryaev–Roberts are posed, not proved.

**Conventions.** The scores are conformity scores, as E35's entry writes them. The nth smoothed p-value is pₙ = (#{i ≤ n : sᵢ < sₙ} + θₙ·#{i ≤ n : sᵢ = sₙ})/n, so a small p means a nonconforming newest unit. The Lean indexes the units newest first: `v 0` is the newest, and `v ∘ Fin.succ` are the units before it.

### M-1 — the smoothed p-values are i.i.d. uniform

**Claim.** Let the units' sequence be exchangeable, the θₙ i.i.d. uniform and independent of it, and the scores at each time a symmetric function of the units so far. Then p₁, p₂, … are i.i.d. uniform on [0, 1]. This is Vovk, Nouretdinov and Gammerman 2003, Theorem 1, and Vovk et al. 2021, Proposition 1.

**The finite core** (`seqRank_uniform`). Write L for the units strictly below the newest and T for those tied with it, itself included. The rank part of pₙ is Rₙ = L + ⌈θₙT⌉, whose tie part is uniform on 1, …, T (`tieRankProb`). Over every ordering of a fixed data set, each rank sequence that n units can produce has total probability 1, so 1/n! under a uniform ordering. So the ranks are independent, and the newest of m units has a rank uniform on 1, …, m.

**The proof's two steps.**
- *Which unit is newest.* Every ordering splits into the newest unit p and an ordering e of the rest (`Equiv.Perm.decomposeFin`). Equivariance moves the newest unit's rank onto p's rank in the fixed data set, and the earlier ranks are those of a fixed (n − 1)-unit set in the order e, so the induction applies.
- *The tie class covering a rank.* For each rank r, exactly one tie class straddles it, and its T members each hold it with probability 1/T, so the per-unit probabilities sum to 1 (`sum_tieRankProb`).

**Under a law** (`seqRank_law`): for any law on data sets invariant under reordering, each valid rank sequence has probability 1/n!.

**The lift, not formalized.** θT = (J − 1) + V with J = ⌈θT⌉ uniform on 1, …, T and V uniform on (0, 1], independent of J. So pₙ = (Rₙ − 1 + Vₙ)/n, with the Vₙ i.i.d. uniform and independent of the ranks, and the pₙ are i.i.d. uniform. The checks compute the joint CDF exactly from θ: 11,400 values of P(p₁ ≤ x₁, …, p_N ≤ x_N) equal x₁⋯x_N, ties included, for both kinds of score.

**The premises, each explicit in the statement, and what breaks each in the checks.**
- *Every earlier unit is ranked, calibration and test alike:* `univ` at every time in `seqRankProb`. Two tests ranked against one fixed calibration set each have uniform ranks, but they are dependent: P(both lowest) = 1/10 against 1/16 (three calibration units).
- *Ties broken by a fresh uniform draw at each time:* `tieRankProb`, multiplied across times. Counting ties for the unit, CB-1's rule, makes an all-tied set accept with certainty, so the ranks are not uniform.
- *The scores at each time a symmetric function of the units so far:* `Equivariant` for every n, which covers a fixed score (inductive) and scores all recomputed from the bag (full). The full form keeping each score from its own time breaks it, because in-sample scores drift.
- *Exchangeability of the whole monitored sequence:* the sum over every ordering, or a law invariant under reordering. In the inductive form the training set proper fits the score and is not in the sequence (CB-2's conditioning).

**Checks.**
- The finite core holds exactly on 68 (data multiset, score) cases up to N = 4, ties included, with inductive and full scores.
- Each premise was broken in turn, with the results above.
- The joint CDF of the lift is exact.

**Prior art, read.** Vovk, Nouretdinov and Gammerman 2003: Theorem 1, with its appendix proof by reverse σ-algebras Gₙ. The finite core here is that proof's inner step, conditioning on the bag at time n, made into a count. Vovk et al. 2021, Proposition 1 (from Vovk, Gammerman and Shafer 2005, Theorem 8.2).

### M-2 — the betting product, Ville's inequality, and the e-value

**The link** (`lintegral_capital_snoc`, `lintegral_capital`).
- *The setting.* The p-values are i.i.d. with law μ, uniform on [0, 1] by M-1. The capital is Sₙ = Π f_i(p_i), with S₀ = 1 (`capital`).
- *The premises on f_i.* Each f_i is ≥ 0 and measurable, chosen from p₁, …, p_{i−1} only, and integrates to 1 under μ. "Only" is the type: `f n` takes the history `Fin.init h`.
- *The conclusion.* Given any history, the next capital integrates to the current one: the martingale property on the path space. So E Sₙ = 1 at every horizon.
- *The proof.* Tonelli, through `MeasurableEquiv.piFinSuccAbove` at the last coordinate.

**The e-value** (`lintegral_stoppedCapital`).
- Stopped by any measurable rule that reads only the past, and by a horizon at the latest, the capital has expectation 1.
- Stopping switches the bet off: the stopped bets `stoppedBet` are again betting functions. Until the stop the stopped capital is the capital (`capital_stoppedBet_of_notMem`), and after it the capital stays put (`capital_stoppedBet_of_mem`).
- An unbounded stopping time gives E[S_τ; τ < ∞] ≤ 1 by Fatou over τ ∧ n. That is not formalized.

**Ville's inequality** (`ville_horizon`, `ville`).
- *At a horizon,* P(S_k ≥ c for some k ≤ n) ≤ 1/c. Stop at the alarm: the stopped capital is at least c on the alarm event, and has expectation at most 1, so Markov's inequality gives the bound.
- *Over every horizon,* with the p-values on any space whose laws at every horizon are μⁿ (M-1's conclusion), P(∃ n, Sₙ ≥ c) ≤ 1/c. The alarm events increase with the horizon, and the bound passes to their union (`tendsto_measure_iUnion_atTop`).
- *Ville needs only ∫ f ≤ 1.* The product is then a supermartingale, and every step above holds with ≤. `ville` takes that weaker premise, so a betting function normalized by quadrature can be scaled down by its rounding and stay valid. The martingale link and the e-value's equality keep the exact premise.

**Not Mathlib's maximal inequality.** The brief asked for Doob's maximal inequality from Mathlib's martingale library. That needs the product as a martingale in conditional-expectation form, which needs a freezing lemma (independence, with the bet chosen from the past). Switching the bet off proves the finite-horizon bound with Tonelli and Markov alone. It states the code's premises directly, and needs no new Mathlib modules: `Probability.Martingale` was not cached, and nothing was fetched.

**Checks.**
- *Exact, on a 4-point grid standing in for [0, 1]:* every path of 4 p-values enumerated, with bets averaging 1 over the grid. Three strategies: Simple Jumper (Algorithm 1 of Vovk et al. 2021, J = 1/100), a power bet, and a power bet switched by the past. For each, E[Sₙ | history] = S_{n−1} on every history, E S₄ = 1, and c·P(max S ≥ c) ≤ 1.
- *Each premise broken:* a bet integrating to 11/10, and one that reads the current p, both break the bound.
- *Monte Carlo (seed 20261108):*
  - Going all in on p₁ ≤ 1/c attains 1/c: 0.4993, 0.1002 and 0.0098 for c = 2, 10 and 100, at 200,000 runs.
  - The bets 1 + λ(1 − 2p) run until absorbed approach it as λ shrinks: 0.0865 at λ = ½ and 0.0995 at λ = 1/5, against 0.1.
  - Over 500 uniform p-values, Simple Jumper (0.42, 0.083) and the power bet with ε = ½ (0.26, 0.059) stay under 1/2 and 1/10.

**Simple Jumper is in this form.** Its capital is a mixture over the Markov chain's paths, and its step ratio S_n/S_{n−1} = Σ_ε w_ε f_ε(pₙ), with w_ε from the past, is itself a betting function. The same holds for the power mixture, and for the plug-in martingale once its kernel density is renormalized on [0, 1].

**Posed, not proved: CUSUM and Shiryaev–Roberts** (Vovk, *Statistical Science* 2021, Propositions 4.2 and 4.4, Corollaries 4.3 and 4.6). For a positive conformal martingale and c > 1:
- the Shiryaev–Roberts alarm time satisfies E τ₁ ≥ c, from Rₙ − n being a martingale, Rₙ = Σ_{i<n} Sₙ/S_i;
- the number of alarms Aₙ satisfies lim sup Aₙ/n ≤ 1/c almost surely;
- both hold for CUSUM, which alarms no more often.

Each Sₙ/S_i is the betting product from i + 1 to n, so `lintegral_capital` gives E Rₙ = n. The rest needs optional stopping at unbounded times and a strong law. That is a measure-theoretic step that changes no premise, so it is left posed.

*Proved since (2026-11-15), with no optional stopping at unbounded times:* the mean times to a false alarm and, in expectation at every horizon, the alarm counts, under `ville`'s premises ("E35's follow-up" at the end of this page). The almost sure lim sup stays posed.

**Prior art, read.**
- Vovk, Nouretdinov and Gammerman 2003 (power martingales, the e-value reading).
- Fedorova, Gammerman, Nouretdinov and Vovk 2012 (plug-in).
- Volkhonskiy, Burnaev, Nouretdinov, Gammerman and Vovk 2017 (inductive conformal martingales against classical CUSUM).
- Vovk et al. 2021 (Simple Jumper).
- Vovk, *Statistical Science* 2021 (Ville's procedure, CUSUM, Shiryaev–Roberts).
- The text extractions were checked against the formulas used here: (6) of 2003, Algorithm 1 of 2021, and (12) to (15) of *Statistical Science*.

---

## Slice 9 — cross-fitted forms and asymmetric CQR

Posed and proved on 2026-11-09 for the [conformal layer's slice 9](../../../issues/regression/archive/bass-gbass-khaos-activespaces/briefs/conformal-layer-brief.md#e16s-completion-the-plan-slices-7-onward), before any C#, from the [chip brief's assignment](conformal-chip-brief.md#e16-slice-9--cross-fitted-intervals-and-asymmetric-cqr-2026-11-09). The point is the contract: the exact rank rules the layer implements, and the hypotheses it must meet.

**Status (2026-11-09).** The jackknife+, CV+ with both of its bounds, and asymmetric CQR are proved in `lean/Enthymemata/ConformalValidity.lean`, PROOF-CLOSED, and `#print axioms` shows the standard three only. CV+ was added the same day, once Barber et al. 2021 was on disk. Later the same day came the order-reading jackknife+, the one-fold CV+ bound, and the position-aware counterexample. [`conformal-checks.py`](conformal-checks.py), section "Slice 9", ran first, and every assertion holds. Not included: cross-conformal predictive distributions (below).

### J+ — the jackknife+ covers at least 1 − (2k − 1)/(n + 1), so at least 1 − 2α

**Conventions.** These extend CB-1's, with conformity scores throughout.
- *Units.* There are m = n + 1 units in `C`: n training units and the test unit, with no split. Units outside `C` stay in every fit's training set, and with `C` every unit the form is the standard one.
- *The leave-two-out score.* s(i, j) is unit i's conformity under the fit to every unit but i and j (`pairScore`). For the jackknife+ it is −|Yᵢ − μ̂₋ᵢⱼ(Xᵢ)| (`jackFit`). The fit leaving out i and the test unit is the jackknife+'s leave-one-out fit μ̂₋ᵢ.
- *Stranger.* Unit u is stranger than j when s(u, j) < s(j, u), strictly, so a tie makes neither stranger (`strangerThan`).
- **The cross-fitted rule at level k** refuses u when u is stranger than at least m − k units (`AcceptsPairwise`). With the leave-two-out scores of any conformity score, this is Vovk's n-fold cross-conformal set.
- **The jackknife+ at level k** (`JackRefuses`), with k = ⌊α(n + 1)⌋:
  - The interval is closed, from the k-th smallest lower end μ̂₋ᵢ(X) − Rᵢ to the k-th largest upper end μ̂₋ᵢ(X) + Rᵢ, where Rᵢ = |Yᵢ − μ̂₋ᵢ(Xᵢ)|.
  - Equivalently, y is refused when at least n + 1 − k upper ends lie strictly below it, or at least n + 1 − k lower ends strictly above it. A y at an end is covered.
  - At k = 0 the interval is the whole line (`acceptsPairwise_zero` for the rule). Past k = n it is empty.
  - These are Barber et al.'s ends: the largest α-quantile of the lower ends with an atom at −∞, and the smallest (1 − α)-quantile of the upper ends with an atom at +∞.

**Claims.**
- *The counting core* (`card_refused_pairwise_le`). Whatever the comparisons and their ties, at most 2k − 1 units are refused at level k, and none at k = 0.
  - Each refused unit is stranger than at least m − k units, at most m − |F| of them outside the refused set F.
  - Inside F each comparison goes at most one way, so the comparisons there sum to at most |F|(|F| − 1)/2.
  - So |F|(m − k) ≤ |F|(|F| − 1)/2 + |F|(m − |F|), which gives |F| ≤ 2k − 1. This is Landau's argument as Barber et al. use it.
- *Over orderings* (`coverage_pairwise_orderings`). For a matrix score equivariant on `C` (`PairEquivariant`), the test unit is accepted on at least 1 − (2k − 1)/m of the orderings of a fixed data set.
- *Under a law* (`coverage_pairwise`, `crossFit_coverage`). The bound holds with probability at least 1 − (2k − 1)/m, for any conformity score of a symmetric fit (`FitSymmetric`).
- *The jackknife+*. `not_acceptsPairwise_of_jackRefuses`: the interval refuses only a strange test unit, so it contains the cross-fitted set. Hence `jackknifePlus_orderings`, `jackknifePlus_coverage`, and `jackknifePlus_two_alpha`: 1 − 2α at any k ≤ α(n + 1).

**The hypotheses, as the layer must meet them.**
1. *Exchangeability of all n + 1 units,* the training units and the test unit together.
2. *A symmetric fit.* On the data reordered by σ, the fit to the units T is the fit to the units σ(T) of the original: it reads its training units as a set (`FitSymmetric`).
   - The law forms take this up to a relabeling `act σ` of the fits' randomness. The law of the randomness and the data together must be invariant under each reordering with its relabeling.
   - So a fit that reads the order of its units qualifies when it reads them in one uniformly random order drawn per run and shared by every leave-out fit, or when its seeds are attached to the units.
   - **A fit that reads its units' order also qualifies, with no randomization** (`OrderReading`, `jackknifePlus_orderReading`, below). A sampler sweeping the units in index order with a fixed seed is such a fit.
   - **A fit that knows its units' positions does not.** That includes seeds keyed by row index, weights by index, and a mask showing which row was left out. The parity construction below covers exactly 1/2.
3. *Nothing else.* The fit may change with its sample size, need not be stable, and need not mix. The theorem does not need the fit to read only its training units. That locality is what makes the interval computable from the training data and X.

**What changed from the task.**
- *The guarantee is 1 − (2k − 1)/(n + 1), above 1 − 2α.*
  - The count is Barber et al.'s own: their Step 1 gives |S| ≤ 2α(n + 1) − 1 before Theorem 1 rounds it to 2α. With k = ⌊α(n + 1)⌋ it is 2k − 1, and the bound is at least 1 − 2α + 1/(n + 1) for k ≥ 1.
  - At k = 1 it is 1 − 1/(n + 1), the split form's exact rate, so the factor of 2 first costs anything at k = 2.
  - At α = 0.1 and n = 19, the worst case is 0.85, not 0.80.
- *The bound is attained.*
  - A symmetric fit whose leave-two-out residuals form a regular tournament on 2k − 1 units covers exactly 1 − (2k − 1)/m. That holds over every ordering, for each (m, k) checked with 2k − 1 ≤ m up to m = 7.
  - Barber et al.'s Theorem 2 makes the same point asymptotically: an i.i.d. law and a fit covering at most 1 − 2α + 6√(log n/n). The case here is the exact finite one.
  - Ordinary fits (mean, median, least-squares line, nearest neighbour) meet the bound at k = 1.
  - So "near 1 − α", the engine brief's figure, is a measurement under stability, not a guarantee. Barber et al.'s section on stability bounds it there.
- *The law form generalizes the task's.* It is stated up to a relabeling of the randomness, which admits samplers symmetric only in law (hypothesis 2).

**Checks** (`conformal-checks.py`, slice 9).
- *The counting core:* all 92,577 comparison outcomes. That is every outcome, ties included, up to m = 5, plus every tournament at m = 6. At most 2k − 1 units are refused, and exactly 2k − 1 is attained, capped at m.
- *The ends:* 304,500 cases, n up to 7, α from 1/20 to 1/2, with tied ends and y at the ends. Barber et al.'s quantiles, the k-th order statistics and the counts agree.
- *Over orderings:* 4 data sets (m = 4 to 6, ties included) × 5 symmetric fits, 4,920 orderings.
  - The leave-two-out matrix is equivariant, the interval refuses only a strange test unit, and every bound holds.
  - The interval and the cross-fitted set covered alike in every case there.
- *Under a law, exactly:* i.i.d. units on four points and a Pólya urn, m = 4, k = 1, against the bound 3/4. The coverage is 0.817 to 0.843.
- *Each premise broken.*
  - The plain jackknife, centred on the fit to all n units, covers nothing at any k ≥ 1 when the fit shifts with its sample size. The jackknife+ on the same fit meets its bound.
  - A tie counted as strange refuses all 4 units of a tied set at k = 1, where the core allows 1.
- *Order-reading fits.* Every fit that picks the stranger of two missing units from the ordered training tuple was searched exhaustively at m = 4: 2¹² fits. None refuses more than (2k − 1)/m. Exact MILPs in [`conformal-searches.py`](conformal-searches.py) give the same up to m = 6.
- *Fits that know their units' positions,* exactly over every ordering for m = 3 to 7: the parity construction refuses exactly 1/2 at every level. That passes 1 − 2α from m = 5 at k = 1.

**What a fit may see** (2026-11-09; the open item on symmetric fits, resolved).
- *Its units' order: the bound holds* (`jackknifePlus_orderReading`, through `coverage_orderReading`).
  - A fit reads its units' order (`OrderReading`) when, on the data reordered by σ, its fit to the units T is the fit to σ(T) whenever σ keeps the units of T in order. It may depend on that order, not on where in the index its units sit.
  - The proof moves any unit p to the front, the others kept in order: `(Fin.cycleRange p).symm`, increasing off 0 by `Fin.strictMono_succAbove`.
  - The test unit's comparisons in the moved data are unit p's in the original (`pairScore_cycleRange`). So the test unit is strange exactly as often as any unit (`acceptsPairwise_cycleRange`), and the counting core gives (2k − 1)/(n + 1).
  - Only these moves are needed, not the full symmetry. The interval itself never reads where the test unit sits.
- *Its units' positions: no guarantee.*
  - Swapping the test unit with training unit i flips the ordering's parity.
  - The fit leaving out position i and the test knows which positions are missing and, over a fixed bag, which two points. So it can make stranger the point at the test position in the even completion.
  - Then in every even ordering the test unit is stranger than every training unit, and in every odd one than none. Coverage is exactly 1/2 at every level 1 ≤ k ≤ n, for every n.
  - The exact MILP worst cases match: 1/2 at k = 1 for m = 3 to 6, and 2/3 at m = 5, k = 2.
  - Barber et al. 2023 handle fits that are not symmetric by changing the method, swapping the test unit with a random index K. They make no claim for the plain jackknife+.
- *The engine's reading.* A leave-out fit receives its training rows in index order, row i dropped and the rest closed up, with no index, mask or row-keyed seed. A sweep in index order with a fixed seed is then covered with no randomization.

**Not included.**
- *CV+* has its own section below.
- *Cross-conformal predictive distributions* (Vovk et al. 2018) are valid only approximately, and are not posed.
- *The asymmetric jackknife+* (Barber et al.'s section on asymmetric forms: signed residuals, α₊ + α₋ = α, the proof run with two signed matrices) is read, not posed.

**Prior art, read.**
- Barber, Candès, Ramdas and Tibshirani 2021 (arXiv 1905.02928v3, on disk since 2026-11-09), read in TeX:
  - the quantiles q̂⁺ and q̂⁻, the ⌈(1 − α)(n + 1)⌉-th and ⌊α(n + 1)⌋-th smallest, with q̂⁺ infinite past n;
  - the symmetry condition, the fit invariant to reordering its data;
  - Theorem 1 and its proof by Landau's theorem, whose Step 1 has the 2α(n + 1) − 1;
  - Theorem 2, the worst cases, where the jackknife covers nothing and the jackknife+ near 1 − 2α;
  - the section on the asymmetric jackknife+ and CV+.
- Barber, Candès, Ramdas and Tibshirani 2023 (arXiv 2202.13415v5), the jackknife+ appendix:
  - the restatement, with its footnote on which quantile is taken;
  - the remark that the n-fold cross-conformal set with any score has the same guarantee and lies inside the jackknife+ interval.

### CV+ — K folds, by merging the folds and by a tournament across them

**Conventions.** These are the jackknife+'s, with folds.
- *Folds.* n = Km training units in K folds of m. Fold k's fit μ̂₋ₖ leaves out fold k and the test unit. Unit i's residual Rᵢ is under its fold's fit, and so is the test unit's comparison with i.
- *Level.* ℓ = ⌊α(n + 1)⌋. Ties count for neither unit, and the interval is closed.
- *The rule* (`CVRefuses`, `CVRefusesAcross`). y is refused when at least n + 1 − ℓ upper ends μ̂₋ₖ₍ᵢ₎(X) + Rᵢ lie strictly below it, or as many lower ends strictly above it. At K = n it is the jackknife+.

**CV-a, merging the folds' ranks** (Vovk; Barber et al., Theorem 4(a)).
- *The rule.* rₖ is the test unit's rank in fold k with itself, ties counted for it as in CB-1, so it lies in 1, …, m + 1. CV+ refuses only when Σₖ rₖ < ℓ + K (`sum_rankLe_lt_of_cvRefuses`).
- *The bound* (`crossFold_refused_le`, `cvPlus_merged`). P(refused) ≤ 2(ℓ + K − 1)/(n + K).
  - With ℓ ≤ α(n + 1) this is at most 2α + 2(1 − α)(K − 1)/(n + K) (`cvPlus_merged_alpha`). That is Vovk's term, as Barber et al. state it for cross-conformal, sharpened by the floor in ℓ.
- *The proof.*
  - Each rₖ is super-uniform, by CB-1 in its fold (`rank_superUniform`).
  - Σₖ rₖ ≤ L forces L ≤ Σₖ (2L/K − rₖ)₊ (`le_sum_max_of_sum_le`).
  - The layer cake bounds E(c − rₖ)₊ by its average over 1, …, m + 1 (`sum_antitone_le`), a right Riemann sum at most c²/(2(m + 1)) (`sum_max_sub_le`).
  - So P ≤ 2L/(K(m + 1)) (`merged_le`). This is the arithmetic mean of p-values, valid up to a factor 2, here in counts. Barber et al.'s route randomizes the p-values and adjusts the level; the counts need neither.
- *The hypotheses.*
  - The n + 1 units are exchangeable, each fold with the test unit given the rest.
  - Each fold's fit is unchanged when its fold and the test unit are reordered: it reads neither.
  - **No symmetric fit is needed.** A fit that reads its other units' order qualifies, and the checks include two. So does one that knows their positions: reordering fold k and the test unit moves nothing the fold's fit reads.

**CV+ by one fold** (`cvPlus_oneFold`, 2026-11-09).
- *The bound.* A refusal leaves the test unit's rank at most ℓ in every fold, so in any one fold k₀. That fold with the test unit is CB-2's split form, so P(refused) ≤ ℓ/(m + 1).
- *The hypotheses.* Fold k₀ is exchangeable with the test unit, and its fit does not read either. It holds on any bag.
- *With two folds* it beats both bounds above. ℓ/(m + 1) ≤ α(2m + 1)/(m + 1) < 2α, so two-fold CV+ covers above 1 − 2α with no extra term. For K ≥ 3 it is about Kα and the others serve.

**CV-b, the tournament across folds** (Barber et al., Theorem 4(b)).
- *The setting.* m − 1 further units from the same law fill the test unit's fold, and the interval never reads them. Units are (fold, slot). Unit i is compared with j under the fit leaving out both their folds (`foldScore`), and only across folds (`strangerAcross`, `AcceptsAcross`).
- *The counting core* (`card_refused_across_le`). At most 2ℓ + m − 2 units are refused, where Barber et al. count 2α(n + 1) + m − 2. The same-fold pairs play no game, so the refused units' wins can be packed into fewer games.
- *The bound* (`coverage_across`, `cvPlus_tournament`). The coverage is at least 1 − (2ℓ + m − 2)/(n + m), which is at least Barber et al.'s 1 − 2α − (1 − K/n)/(K + 1) (checked). At m = 1 it is the jackknife+'s bound.
- *The hypotheses.*
  - The n + m units are exchangeable under permuting the folds and the slots within them (`ratio_prodPerm` averages over those).
  - i.i.d. data give that. A fixed data set of n + 1 units, sampled without replacement, does not: its law need not extend to n + m units.
  - A symmetric fit up to a relabeling of its randomness, as for the jackknife+.
- *On a fixed bag of n + 1 units,* without the further test units (2026-11-09; the open item, partly resolved).
  - *The reduction.* For each left-out set L, a fold with the test unit, a symmetric fit ranks L's residuals. Let d_L(x) count the units of L ranked strictly above x; under a uniform ordering it is uniform on 0, …, m. CV+ refuses iff Σₖ d_{Fₖ ∪ x}(x) ≤ ℓ − 1.
  - *Settled at K = 2.* The one-fold bound ℓ/(m + 1) lies within CV-b's at every level for m ≥ 2, and m = 1 is the jackknife+.
  - *Settled at K = 3, m = 2, ℓ ≤ 2.* At ℓ = 2 a refusal needs the test unit on top in two of its three triples, against once on average, so P ≤ 1/2, CV-b's bound by Markov. At ℓ = 1 the exact worst case is 1/7, against CV-b's 1/4. To see this, label the seven points of a Fano plane uniformly and pick a uniform test point. Its three incident lines pair up the other six points uniformly into folds. Refusal requires the test point to win on every incident line. Two points share a line and cannot both win it, so at most one point per plane can be refused. Ranking a designated point highest in every triple containing it attains 1/7. The [exact checks](conformal-searches.py) enumerate the 30 labelled planes and verify that all 105 choices of test and folds occur twice.
  - *Searched, no counterexample.* Annealing over K = 3, 4, 5 with m = 2, and K = 3 with m = 3, at every level found nothing past CV-b's bound, and the worst found stays well inside it: 0.40 against 0.58 at K = 3, m = 3, ℓ = 3 ([`conformal-searches.py`](conformal-searches.py)).
  - *Still open for K ≥ 3,* with no decision waiting: on a fixed bag the layer has the merged and one-fold bounds, which need no extension.

**Which bound serves.** The miscoverage bounds at α = 0.1:

| n | K | merged | tournament |
|---|---|---|---|
| 20 | 2 | 0.273 | 0.400 |
| 20 | 5 | 0.480 | 0.250 |
| 100 | 5 | 0.267 | 0.317 |
| 100 | 10 | 0.345 | 0.255 |
| 1000 | 5 | 0.207 | 0.332 |
| 1000 | 10 | 0.216 | 0.271 |
| n | n | over 1 | the jackknife+'s |

The merged bound serves few folds and many units, and needs no symmetric fit. The tournament serves many folds. Barber et al.'s minimum of the two is at most 2α + √(2/n).

**Checks** (`conformal-checks.py`, slice 9, CV+).
- *The merging inequality:* 21,142 (K, M, rank vector, L) cases, and the layer cake on super-uniform laws.
- *With K = 2 and uniform ranks,* the worst coupling lies within 1/M of the merged bound: 1/2 against 2/3 at M = 6, L = 4.
- *CV-a over orderings:* 51,600 orderings, 4 data sets × 5 fits, two of them reading their units' order. Every level holds, and every refusal has Σₖ rₖ < ℓ + K.
- *CV-b over orderings of n + m units:* 4,392 orderings, 3 data sets × 3 symmetric fits. Every level holds.
- *CV-b's counting core:* 551,934 outcomes of the cross-fold comparisons, ties included. The bound is attained for 2 and 3 folds of 2 and for 4 folds of 1. For 2 folds of 3 the worst is 4 against 5 at ℓ = 2.
- *Both bounds inside Barber et al.'s α forms:* 8,400 (K, m, α) cases.

**Prior art, read.** Barber, Candès, Ramdas and Tibshirani 2021 (1905.02928v3), read in TeX:
- the CV+ section and Theorem 4;
- the subsection on cross-conformal predictors: the cross-conformal set inside CV+, and the remark that both proofs cover it;
- the proofs for CV+: the modified cross-conformal form, Vovk's averaging of p-values, and the CV+ tournament with its test fold.

Vovk 2015, Vovk et al. 2018 and Vovk and Wang 2012, on combining p-values, are cited through it and were not read.

### CQR-A — asymmetric CQR, a threshold per tail

**Setting.** This is the split form (CB-2).
- The training split w fits q_lo and q_hi. Given w, the calibration and test units are exchangeable.
- Each tail has its own conformity: y − q_lo(x) below (`cqrLoScore`) and q_hi(x) − y above (`cqrHiScore`), the negatives of Romano et al.'s E_lo and E_hi.
- Each tail is ranked by CB-1's accept rule at its own level, k_lo and k_hi, with ties counted for the test unit. At α/2 per tail, k_lo = k_hi = ⌊(α/2)(n₂ + 1)⌋.

**The ends.**
- The set is the closed interval [q_lo(x) + c_lo(k_lo), q_hi(x) − c_hi(k_hi)]. Here c(k) is the k-th smallest calibration conformity of that tail.
- At level 0 a tail's end is infinite (`accepts_zero`).
- Romano et al.'s Q_{1−α_lo}(E_lo) is the ⌈(1 − α_lo)(n₂ + 1)⌉-th smallest of q_lo(xᵢ) − yᵢ, infinite past n₂. It gives the same lower end, and the upper end likewise.

**Claims.**
- *Each tail and both* (`cqrAsym_coverage`). The test unit falls below the set with probability at most k_lo/m and above it with probability at most k_hi/m. The set covers with probability at least 1 − (k_lo + k_hi)/m, whatever the fits.
- *At α/2 per tail* (`cqrAsym_coverage_half`), coverage is at least 1 − α.
- *Each tail exact with distinct scores.* CB-1's `coverage_orderings_eq` with `splitScore_equivariant` gives exactly k/m per tail. Random tie-breaking (`coverage_tieBreak`) makes any data set so.
- *The two-sided coverage is at most 1 − max(k_lo, k_hi)/m.* Checked, not formalized: the set lies inside each tail's.
- *Each tail refuses on its own side* (`cqrLo_isUpperSet`, `cqrHi_isLowerSet`). The lower tail accepts an upper set of y and the upper tail a lower set, so the set is an interval (`cqrAsym_ordConnected`). A lower-tail refusal is y below the set, and an upper-tail refusal y above it.

**Small calibration sets and rounding.**
- *A tail's end is finite only once ⌊(α/2)m⌋ ≥ 1,* that is m ≥ 2/α: n₂ ≥ 19 at α = 0.1, against n₂ ≥ 9 for the symmetric form.
- *Rounding each tail.* Since 2⌊αm/2⌋ is ⌊αm⌋ or one less, the asymmetric lower bound equals the symmetric finite-rank bound or exceeds it by 1/m. Its excess over nominal 1 − α is (αm − 2⌊αm/2⌋)/m, in [0, 2/m). At α = 0.1 and m = 31 the asymmetric and symmetric lower bounds are 29/31 and 28/31, respectively, with nominal coverage 9/10.

**Checks** (`conformal-checks.py`, slice 9).
- *Under an i.i.d. law,* N = 4, every pair of levels, ties included: each tail and both hold. At k_lo = k_hi = 1, the set misses below 0.190 of the time and above 0.170, and covers 0.640.
- *Over orderings with distinct scores:* each tail refuses exactly k/m at every level, and two-sided coverage lies in [1 − (k_lo + k_hi)/m, 1 − max(k_lo, k_hi)/m].
- *The ends:* 16,384 accept sets on a grid match the closed form and Romano et al.'s quantile.
- *Symmetric CQR controls only the total.* With both quantile fits below the data, every refusal falls above the set: 0.363 above, against 0.127 for the asymmetric form at k = 1 per tail.

**Prior art, read.** Romano, Patterson and Candès 2019 (arXiv 1905.03222v1), read in TeX: Theorem 2 with its proof, and their quantile Q_{1−α}(E, I₂), the (1 − α)(1 + 1/|I₂|)-th empirical quantile. Linusson et al. 2014, which they cite, was not read.

## Slice 10 — the full step's error, completed

Posed and proved on 2026-11-10 for the [conformal layer's slice 10](../../../issues/regression/archive/bass-gbass-khaos-activespaces/briefs/conformal-layer-brief.md#e16s-completion-the-plan-slices-7-onward), before any C#, from the [chip brief's assignment](conformal-chip-brief.md#e16-slices-10-and-12--the-full-steps-error-hoffs-score-and-the-latent-posterior-2026-11-10). It closes CB-3's open item, the margin under bounds that hold only with high probability.

**Status (2026-11-10).** Proved in `lean/Enthymemata/ConformalValidity.lean`, section `FullStepError`, PROOF-CLOSED, and `#print axioms` shows the standard three only. [`conformal-checks.py`](conformal-checks.py), section "Slice 10", ran first, exhaustive in exact arithmetic, and every assertion holds. Checked, not formalized: the PIT between grid values.

**Setting.** As in CB-3:
- the exact score S is equivariant;
- the estimates ŝ are read per ordering;
- the bounds εⱼ are carried by the units.

The full step's CDF score F̃ is equivariant, since each unit's add-one-in predictive CDF is symmetric in the units, so everything here applies to it. The equi-tailed set at k per tail is CB-1's rule on F̃ (the lower tail) and on −F̃ (the upper tail), `AcceptsTwoSided`. The full step takes k = ⌊αm/2⌋.

### The equi-tailed set

- *Exact.* It refuses at most 2k units (`card_not_acceptsTwoSided_le`).
- *Margin form* (`coverage_close_twoSided`). Coverage is at least 1 − |R ∪ M_lo ∪ M_hi|/m, so at least 1 − 2k/m − |M_lo ∪ M_hi|/m.
  - R is the exact set's refusals.
  - M_lo and M_hi are CB-3's fragile units of F̃ and of −F̃ (`FragileTwoSided`).
- *Inflated form* (`coverage_inflated_twoSided`). With each tail inflated by the same bounds, coverage is at least 1 − 2k/m, with bounds by ordering and position.
- *Do the tails' bands add or share units?* They enter as their union, never more than their sum.
  - A unit fragile in both tails has at least m − 2k + 2 units within εⱼ + ε_u of it, itself included (`card_near_of_fragile_both`). At k = ⌊αm/2⌋ that is more than (1 − α)m of the set.
  - Short of that bunching, the bands add.
- *Dominance* (`card_fragileTwoSided_le_estimated`). For every ordering, |M_lo ∪ M_hi| is at most the union band read from that ordering's estimates (each tail's `FragileEstimated`, shifts doubled). As for the density score, any ordering's reading is a valid report.

### The distribution: P(PIT ≤ u) − u

- *The PIT.* Given the scores and the ordering, the test unit's PIT is (R − 1 + V)/m.
  - R = B + J, with J uniform on 1..T + 1 (`tieRankProb`), and V uniform on (0, 1].
  - So P(PIT ≤ k/m) = P(R ≤ k) (`pitRankLe`).
  - With exact scores it is k/m at every grid value (checked).
- *Bounds* (`pit_ge_lower`, `pit_le_upper`). With estimates within the bounds,

  #{u : U(u) ≤ k}/m ≤ P(PIT ≤ k/m) ≤ #{u : L(u) < k}/m,

  where L(u) = #{j : F̃ⱼ + εⱼ + ε_u < F̃_u} (`shiftedBelow`) and U(u) = #{j : F̃ⱼ ≤ F̃_u + εⱼ + ε_u} (`shiftedAtMost`).
- *The two counts bracket k* (`card_shiftedAtMost_le_le_card_shiftedBelow`): #{U ≤ k} ≤ k ≤ #{L < k} for k ≤ m. So

  |P(PIT ≤ k/m) − k/m| ≤ max(#{L < k} − k, k − #{U ≤ k})/m.

  - With distinct exact scores, #{L < k} − k counts the units above level k whose shifted count falls below it. Likewise k − #{U ≤ k} counts the units at or below level k whose shifted count rises above it. These are the PIT's fragile units at level k.
  - Under ties the exact counts already pass k, and the bound is loose by the tied mass, while the randomized PIT stays uniform.
- *Uniformly in u.* Given R, the PIT is uniform on ((R − 1)/m, R/m]. So P(PIT ≤ u) is linear between grid values, as u is, and the largest departure sits at a grid value: sup_u |P(PIT ≤ u) − u| is the largest of the grid bounds. Checked, not formalized.
- *Dominance* (`card_pit_counts_estimated`). For every ordering, the counts read from its estimates with the shifts doubled bound the exact ones: #{L < k} ≤ #{L̂ < k} and #{Û ≤ k} ≤ #{U ≤ k}.

### Leaving out the degenerate turns

**Setting.** The draws w, the Monte Carlo seeds, have weights p summing to 1. For each draw:
- D(w) ⊆ C is the set of turns whose weights degenerate. It is read once for the data set (by k̂ or the ESS of each turn's weights), so it is the same for every ordering.
- ε(w) is the set of bounds carried over the kept turns only.
- ν bounds the weighted share of kept orderings (test unit outside D) on which some unit's bound fails. It is joint over the m units of an ordering.

On D's turns, and on kept orderings whose bounds fail, the estimate may be anything.

**Claims.**
- *The accept rule* (`coverage_close_excluding`). Coverage is at least 1 − |D ∪ R ∪ M_D|/m − ν.
  - R is the exact refusals (at most k).
  - M_D is the set of fragile units under the kept bounds.
- *The share is |D \ (R ∪ M_D)|/m.* A degenerate turn costs only when its unit is neither refused nor fragile already. The brief's form, 1 − (k + |M_D|)/m − |D|/m − ν, follows, and is never tighter. A far-out unit sits in the exact rule's low tail, so its turn usually costs nothing.
- *The inflated rule* (`coverage_inflated_excluding`): at least 1 − |D ∪ R|/m − ν, so at least 1 − (k + |D|)/m − ν.
- *The equi-tailed set* (`coverage_close_twoSided_excluding`, `coverage_inflated_twoSided_excluding`): at least 1 − |D ∪ R ∪ M_lo ∪ M_hi|/m − ν, and the inflated set at least 1 − |D ∪ R|/m − ν, with R its refusals (at most 2k).
- *The PIT* (`pit_excluding`):

  #{u ∉ D : U(u) ≤ k}/m − ν ≤ P(PIT ≤ k/m) ≤ #{u : u ∈ D or L(u) < k}/m + ν.
- *The steps.*
  - `coverage_ge_of_good_fail` is for indicators.
  - `avg_le_of_fail` and `avg_ge_of_fail` are for any reading in [0, 1]. A reading controlled on the kept orderings whose bounds hold averages within the controlled share, give or take ν.

**What the engine reads.** The engine reads these from one ordering whose kept bounds hold: units at their positions, the observed ordering in the engine.
- *The margin* (`card_excluded_le_estimated`). |D ∪ R ∪ M_D| is at most k + |M̂_D| + |{u ∈ D : the inflated rule accepts u}|. A degenerate turn that the inflated rule refuses is refused by the exact rule too (the contrapositive of `acceptsInflated_of_accepts`), so it costs nothing.
- *The equi-tailed set* (`card_excluded_le_estimated_twoSided`). The same bound holds with 2k, the union band, and the units of D that the inflated set accepts.
- *The PIT* (`card_pit_excluding_estimated`, through `shiftedBelow_estimated_le` and `shiftedAtMost_le_estimated`). The counts with D, read with the shifts doubled, bound the exact ones from the safe side.

**What it rests on.**
- *D must be fixed for the data set.* A turn is labelled by its unit and does not depend on which unit is the test unit. If the test unit's own turn were flagged after its estimate was seen, the excluded orderings would no longer be a share of units, and the bound would not follow. The precise rule, and how the engine reads it, is under "Reading from the realized ordering" below.
- *ν must be joint per ordering.* The layer's z = Φ⁻¹(1 − ν/(2m)) per unit gives it by a union bound.

**Against Bhattacharyya, Zhang and Barber.** Their 1 − α − 3√ν (arXiv 2605.29200, appendix) is for the inflated rule, with ν the chance that one score's estimate misses its bound. Ours takes ν per ordering, jointly over its m units. It pays ν once, for the margin and the inflated rule alike.
- Theirs, for a coverage cost of 0.01: it needs ν ≈ 1.1·10⁻⁵ per point, so z ≈ 4.39 for a two-sided normal bound.
- Ours: ν = 0.01 per ordering gives z ≈ 3.89 per unit at m = 100, and 4.42 at m = 1000.
- At ν = 0.05 theirs costs 0.67.

**An illustration at m = 100** (checks). One unit far out sits both in D and in the exact low tail; k = 5 and ν = 1/20.
- With bounds carried over every turn (3, widened by the far unit's turn), the margin reads 0.00, with 95 units fragile.
- Leaving out the far unit's turn (bounds 0.002), it reads 0.88, that is 1 − (5 + 2)/100 − 0.05, and D costs nothing.
- These bounds are chosen, not measured. BASS's own reading waits for the C#.

**Checks** (`conformal-checks.py`, slice 10).
- *The equi-tailed set:* 768 cases (m = 3 to 5, three bound profiles, a symmetric squared score and an in-sample CDF score, every k), each over every ordering with an adversarial and a noisy estimate.
- *The PIT:* 276 cases. The exact PIT is uniform at every grid value. The estimate's P(R ≤ k) lies within the counts, and each ordering's estimated counts dominate. Between grid values the departure is at most the larger at the ends.
- *Leaving out D:* 261 one-sided cases, and 489 two-sided and PIT cases. These include the engine's reading from every ordering.

**Decisions served.**
- *Can the full step's distribution report its Monte Carlo error?* Yes:
  - for the PIT, max(#{L̂ < k} − k, k − #{Û ≤ k})/m at each grid level, read from any ordering, or the form with D plus ν;
  - for the equi-tailed set, the union band.

  Whether to show it depends on whether it reads anything on real fits, which the C# measures.
- *Can the margin read anything on fits like BASS's?* With D left out, yes, when the degenerate turns are the far units that the exact rule already refuses. Their turns cost nothing, and the band shrinks to what the kept turns' bounds give.

### Reading from the realized ordering (added 2026-11-10, after verification)

The main session's verification raised two points about the engine's reading. Both are settled here.

**1. D is each unit's own turn, the candidate's included.**
- *What the theorems need.* They take one set D per data set and draw, the same whichever unit is the test (`D w` in `coverage_close_excluding` is fixed for every ordering). The kept bounds must hold on every ordering whose test unit is outside D. A unit's turn is the ordering in which it is the test. So the set to use is D_own, the units whose own turns degenerate by the per-turn rule.
- *In the full step at grid value y,* the data set has m units, the candidate (x_new, y) among them.
  - The candidate's own turn is the held ordering at y, read by its held add-one-in weights' k̂ or ESS.
  - A data unit's own turn is the ordering in which that unit is the candidate and (x_new, y) is a data unit.
  - So D(y) is read per (data set, y) over all m units, by one rule.
- *A D read from the n data units' turns alone fails.* It never contains the test unit, so it changes with which unit is the test, and the hypothesis fails.
- *What the engine can read.* It reads the candidate's own flag at y directly. A data unit's own-turn flag it reads only through a proxy, since that turn's draws are not held.
  - A proxy set is sound if it contains D_own at every y. A reading made with a larger set bounds the count made with D_own, so extra flags only charge more.
  - A proxy that misses a unit of D_own is not sound.
  - At a far y the far candidate is a data unit of the other turns. So a proxy must read each data unit's own add-one-in weight under a fit that includes the candidate, not the candidate's held weight.

**2. The engine reads only the realized ordering.**
- *At a flagged y the held ordering gives no reading at all.* When the candidate's own turn at y is in D, that ordering's estimates carry no guarantee for any unit, since every unit's estimate there uses the candidate's add-one-in weights. So the candidate cannot be charged as one unit of D with the rest read from those estimates. Another kept ordering of the same data set would serve, but the engine holds none.
- *Per data set, proved.* Each ordering reads its own estimates, and a left-out ordering reads nothing.
  - *Coverage* (`coverage_read_excluding`, `coverage_read_twoSided_excluding`, through `le_read_excluding`): at least 1 − (the orderings' average reading)/m − |D|/m − 2ν. Here the reading is ĉ(σ) = k + |M̂_D| + |D's units the inflated rule accepts| (2k and the union band for the equi-tailed set), from `card_excluded_le_estimated`.
  - *The PIT* (`pit_read_excluding`, through `read_le_excluding` as well): P(PIT ≤ k/m) lies between (the average lower reading)/m − 2ν and (the average upper reading)/m + |D|/m + 2ν.
  - *D is charged twice in the upper forms:* inside the count, for D's own orderings, and again for the orderings that read nothing.
  - *The second ν* is the realized ordering's own bounds failing, which voids its reading.
  - *A lower reading needs no D charge,* since a left-out ordering reads 0.
- *Over the grid.* This is the layer's reading, and it is not formalized, as the maximum over the grid was not before. The true y is unknown. In each realization the reading at the true y is the reading at one kept grid value, or nothing when that value is flagged, and |D| at the true y is at most its largest value over the grid. So the sound reports are:
  - *the margin:* 1 − max over kept y of ĉ(y)/m − max over all grid y of |D(y)|/m − 2ν;
  - *the PIT's bounds:* [min over kept y of (lower reading)(y)/m − max_y |D(y)|/m − 2ν, max over kept y of (upper reading)(y)/m + max_y |D(y)|/m + 2ν].
- *Does the maximum collapse when the grid reaches far out? No.*
  - Flagged y drop out of the maximum, and they are paid for by D's share.
  - A far y costs its one flagged unit, the candidate, in max_y |D(y)|.
  - The report collapses only if the flagged count itself is large at some grid value, as it would be with proxies that inherit the far candidate's weights.
- *Other readings.* A reading per y is not a coverage statement, since the true y is unknown. A reading over the candidates outside D is sound only with D's share charged, as above.

**Checks** (`conformal-checks.py`, slice 10, "reading from the realized ordering").
- 400 random data sets with arbitrary readings on left-out and failed orderings: both reading bounds hold.
- On the 261 (data, D, k) cases of the exclusion check, the engine's own reading from every kept ordering gives coverage ≥ 1 − average reading/m − |D|/m − 2ν.

### The grid reading (added 2026-11-10)

Formalized for slice 10's C#. It covers what the realized-ordering section checked but did not prove. Sections `GridReading` in `ConformalValidity.lean`, and "Slice 10's grid reading" in the checks.

**The PIT between grid values.**
- *The CDF* (`pitCdf`): P(PIT ≤ x) = Σ_r P(R = r) · clamp(m x − r + 1, 0, 1), with PIT = (R − 1 + V)/m.
- *It is linear between grid values* (`pitCdf_interp`): at x = (k + θ)/m, with θ ∈ [0, 1], it equals (1 − θ) P(R ≤ k) + θ P(R ≤ k + 1).
  - The same holds averaged over the draws and the orderings (`avg_pitCdf_interp`).
- *So the largest departure sits at a grid value* (`interp_dev_le`): |P(PIT ≤ x) − x| ≤ max(|P(R ≤ k) − k/m|, |P(R ≤ k + 1) − (k + 1)/m|).
  - sup_x |P(PIT ≤ x) − x| is the largest departure over the grid values k/m.
  - Bounds at the grid values (`pit_excluding`, `pit_read_excluding`) interpolate.

**The margin over the grid** (`coverage_grid_excluding`).
- *The engine's inputs.* It reads `eng w x` from an arrangement x (the data, with the test unit's value at its position), and flags positions `flagged w x`. The realized ordering's flag on its test unit is the data set's D, and the flagged count covers |D|.
- *The report* (`gridReport`): the largest reading at a grid value whose turn is kept, plus the largest flagged count, both over the grid's values.
- *The claim:* coverage is at least 1 − (the orderings' average report)/m − 2ν.
- *The report is read from the data.* The reading at the true value is one of the grid's (`le_gridMax`). With the grid fixed, the report is free of the test unit's own value (`gridMax_update`).
- *What the maximum needs: the grid must hold the value the test unit takes.*
  - The grid may be built from the rest of the data, for example over y at the test unit's own x.
  - In the checks, a grid that missed one value missed its reading in 222 of 1,620 arrangements.
  - So on a response with finitely many values, the grid takes all of them.
  - On a continuous response, the maximum over a fixed grid is a reading at the grid's points, not a bound. Each count is constant between the y where some comparison between estimates flips, or a flag changes. So a sound reading needs a grid value in each such interval: the breakpoints of the readings and flags themselves. Slice 8's boundary refinement locates the crossings of π(y), not these.
- *The PIT's counts over the grid* follow by the same steps: `le_read_excluding` or `read_le_excluding`, with `le_gridMax` bounding the reading at the true value.

**Checks.**
- The interpolation, and the departure bound at the ends, on 4,339 (data, unit, level) cases, ties included, averaged over orderings with noisy estimates too.
- The grid report's coverage bound on 90 (data, k) cases.
- The report's freedom from the test unit's value on 1,620 arrangements.
- The grid that misses a value, as above.

## Slice 12 — Hoff's score and the latent posterior

Posed and proved on 2026-11-10 for the [conformal layer's slice 12](../../../issues/regression/archive/bass-gbass-khaos-activespaces/briefs/conformal-layer-brief.md#e16s-completion-the-plan-slices-7-onward), before any C#, from the [chip brief's assignment](conformal-chip-brief.md#e16-slices-10-and-12--the-full-steps-error-hoffs-score-and-the-latent-posterior-2026-11-10). It closes CB-4's Hoff finding.

**Status (2026-11-10).** Proved in `lean/Enthymemata/ConformalValidity.lean`, PROOF-CLOSED, and `#print axioms` shows the standard three only. The sections are `ShiftedBounds`, `Hoff`, `DrawnLatents`, `HoffOptimal` and `HoffLink`. [`conformal-checks.py`](conformal-checks.py), section "Slice 12", ran first, in CB-4's model and in exact arithmetic, and every assertion holds.

### Hoff's score

- *The score* (`hoffScore`, every latent integrated): cᵢ = p(zᵢ | z₋ᵢ) = Σ_θ p(θ | z₋ᵢ) g(zᵢ | θ), which equals p(z)/p(z₋ᵢ) (`hoffScore_eq`). The numerator is the same for every unit.
- *Validity.* The score is equivariant (`hoffScore_equivariant`), so CB-1's rule covers with it (`hoff_coverage`), ties counted for the candidate. Hoff's strict rule covers 0 when every unit ties (checked; the CB-4 finding).
- *From the n-unit draws:*
  - The candidate's score is the plain n-unit predictive, c₍ₙ₊₁₎ = A = E_n[g(z_new | θ)] (`hoffScore_none`).
  - Each data unit's score is cᵢ = A/Bᵢ, with Bᵢ = E_n[g(z_new | θ)/g(zᵢ | θ)] (`hoffScore_some`, `hoffMeanWeight`).
  - The posterior given every unit but i is the n-unit posterior times g(z_new | θ)/g(zᵢ | θ), which is add-one-in and leave-one-out at once (`postUExcept_some`).
- *The accept rule needs only the Bᵢ.* The candidate is accepted at level k iff at least k data units have Bᵢ ≥ 1 (`hoff_accepts_iff`).
- *The estimator's error, in a form slice 10 carries.* From S draws, ĉᵢ = Â/B̂ᵢ and ĉ₍ₙ₊₁₎ = Â, with Â common to every unit of the ordering.
  - Ranks see only differences (`rankLe_add_const`, `accepts_add_const_iff`, `acceptsInflated_add_const_iff`), so CB-3 holds with bounds up to one shift per ordering (`coverage_close_shift`, `coverage_inflated_shift`). Every slice 10 form transfers the same way.
  - On the log scale the bounds are |log B̂ᵢ − log Bᵢ| for data units and 0 for the candidate.
  - The weights g(z_new)/g(zᵢ) are where PSIS's k̂ applies. A unit whose weights degenerate is one of slice 10's D.

### Under latents

- *(b), every latent integrated:* `hoffScore` itself.
- *(c), read literally:* E[f(zᵢ | θ, vᵢ)] under the joint posterior of θ and every latent given the other units. It equals (b) (`hoffJointScore_eq`): given the others, unit i's latent is at its prior. The latent posterior adds nothing to a leave-one-out score.
- *(c), read as the draws read it:* θ given the others, and unit i's latent given θ and its own value, Σ_θ p(θ | z₋ᵢ) Σ_v q f(zᵢ)²/g(zᵢ) (`hoffCondScore`).
  - This is what the C# computes when the data units read their retained latents with Hoff's weights and the candidate draws its latent from p(v′ | θ, z_new), unweighted (`hoffCondScore_some_draws`, `hoffCondScore_none_draws`).
  - It is equivariant (`hoffCondScore_equivariant`).
  - It sits above (b) unit by unit (`hoffScore_le_hoffCondScore`, through `integrated_le_sq_div`).
- *Reading the units differently breaks it, as CB-4(d) does.* With the data units conditioned and the candidate integrated, the score accepts no candidate that the conditioned score refuses (`accepts_hoffCond_of_mixed`). CB-4(d)'s model gives coverage 0 against 1/2 (`hoffMixed_counterexample`).

In the checks' model at m = 4, under the model's own law and then an iid law outside it:

| score | k = 1 (bound 0.75) | k = 2 (bound 0.50) |
|---|---|---|
| (b) | 0.813, 0.825 | 0.595, 0.587 |
| (c), conditioned | 0.812, 0.816 | 0.594, 0.583 |
| mixed (data conditioned, candidate integrated) | 0.531, 0.485 | 0.260, 0.226 |
| reversed | 0.889, 0.909 | 0.714, 0.744 |

### Hoff's optimality

**Setting.**
- P is any exchangeable, positive law of the m units, and μ > 0 is a base measure on a unit's values.
- A procedure is φ, with φ(z) the chance that the candidate's value z_t is accepted given the others. It may be any function: randomized, and not necessarily symmetric in the data.
- Its expected size is Σ_z P₋ₜ(z) μ(z_t) φ(z) (`expectedSize`).
- Hoff's score for P is cᵢ = P(z)/(μ(zᵢ) P₋ᵢ(z)) (`hoffP`). It is equivariant under any exchangeable P (`hoffP_comp`), and it is the model's `hoffScore` when P is the model's law and g a μ-density (`hoffScore_eq_hoffP`). With covariates, μ(x, y) = h(x) gives sets over y, and c is p(yᵢ | xᵢ, the rest).

**Claims.**
- *Averaging over the orderings* (`sum_perm_expectedSize`): the size is Σ_z P(z) Σ_σ φ(z∘σ)/c_σ(t)(z), up to the number of orderings. For each data set this is a linear program in the orderings' readings with weights 1/c, and a threshold in c solves it (`sum_mul_le_of_threshold`).
- *The theorem* (`hoff_expectedSize_le`). Let ψ be a threshold in Hoff's score on every data set: CB-1's rule on it, or that rule with ties broken at random. Every φ that covers at least as much as ψ on every data set has an expected size at least ψ's.
  - For CB-1's rule this is `hoff_cb1_expectedSize_le`.
  - With distinct Hoff scores, CB-1's rule covers exactly 1 − k/m, so it is smallest among all procedures valid at 1 − k/m.
  - This is Hoff's Theorem 4.1 in finite form: there, constant coverage forces coverage given the multiset.

**The tie finding.** Under CB-1's rule, ties count for the candidate.
- So on a data set where Hoff's scores tie at the boundary, Hoff's set covers more than 1 − k/m, and another score's CB-1 set that covers less can be smaller.
- In the checks' model, where values repeat: at m = 3, k = 1, CB-1's Hoff set has expected size 3.719, against 3.678 for the in-sample score's. The in-sample set covers less than Hoff's on 2 of 56 data sets.
- With ties broken at random, every set is exactly 1 − k/m and Hoff's is the smallest.

**What it buys.** Expected sizes with ties broken at random, counting (x, y) values:

| m, k | Hoff | in-sample (Fong and Holmes) | (c) | (c), conditioned Hoff |
|---|---|---|---|---|
| 3, 1 | 3.2538 | 3.2594 (+0.2%) | 3.3368 (+2.6%) | 3.3151 (+1.9%) |
| 4, 1 | 3.7783 | 3.7840 (+0.2%) | 3.8224 (+1.2%) | 3.8199 (+1.1%) |
| 4, 2 | 2.2011 | 2.2123 (+0.5%) | 2.2953 (+4.3%) | 2.2699 (+3.1%) |

Hoff's set is also no larger than any of 112 random procedures repaired to cover as much on every data set.

The in-sample and leave-one-out predictives differ by one unit's weight, so their gap should shrink as n grows. That is not checked here.

### CB-4(c) read from draws, the candidate's latent drawn per draw

- *The estimator.*
  - It reads the n-unit joint draws (θ_s, v_s) with add-one-in weights w_s = g(z_new | θ_s).
  - The data units read their retained latents, f(zᵢ | θ_s, v_s,i).
  - The candidate draws v′_s from p(v′ | θ_s, z_new) and reads f(z_new | θ_s, v′_s).
- *Its expectation is (c) for every unit, the candidate included* (`latentScore_some_draws`, `latentScore_none_draws`). These are the weighted readings under the n-unit joint posterior.
- *What the drawn latent adds* (`drawnLatent_mean_var`; also `sum_pi_prod_mul_sum`, `sum_pi_prod_mul_sq`). For fixed draws, over the drawn latents:
  - the candidate's estimate has mean Σ_s w̄_s E[f | θ_s, z_new], the Rao–Blackwell reading;
  - its variance is Σ_s w̄_s² Var(f | θ_s, z_new), with w̄ = w/Σw.

  This variance is the term that the candidate's bound in each ordering carries on top of the draws' own error. So each unit's carried bound carries it from its own turn as the candidate. Reading the conditional mean Σ_v q f²/g in closed form removes it.
- *The normalizer drops out.* Σ_s w_s is common to every unit of an ordering, so on the log scale it drops out as a common shift (`coverage_close_shift`). The bounds are on the log numerators.

**Checks** (`conformal-checks.py`, slice 12).
- Hoff's ratio form, equivariance and CB-1 coverage on 342 data sets.
- The draw identity and the B-form accept rule on all 1,296 data sets of 3 units and a candidate.
- 300 Monte Carlo draw sets for the common shift.
- (c) read literally equals (b), and the conditioned reading's equivariance, draw identity and order above (b) hold, on all 216 data sets.
- The coverage table.
- The drawn latent's mean and variance on 300 draw sets.
- The optimality comparisons at m = 3 and 4.

**Decisions served.**
- *Does the full step offer Hoff's score?* The evidence:
  - it is valid, and optimal at equal coverage;
  - in the checks' model it is 0.2–0.5% smaller than the in-sample score at exact coverage;
  - under CB-1's tie rule on discrete data it can be larger;
  - it costs a weight per data unit, which adds PSIS instability and degenerate turns.
- *Does the full step read a conditioned observation latent?*
  - For Hoff's score the literal latent posterior is (b), so there is nothing to read.
  - The conditioned reading is valid only when every unit reads alike, with the candidate's latent drawn per draw.
  - For CB-4(c), the drawn latent costs Σ w̄² Var on the candidate, which a closed-form conditional mean removes.

## Slice 13 — the tournament correction

Bhattacharyya, Zhang and Barber, arXiv 2605.29200v5, read in TeX: §2–§3, Theorems 1 to 3 and their proofs. Proved on 2026-11-10 in `ConformalValidity.lean`, section `Tournament`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 13") ran first.

**The rule is slice 9's cross-fitted rule.**
- Each training unit i plays a game against the test unit. Both are scored by one approximate fit, D₋ᵢ (every training unit but i) plus the pair {zᵢ, z_test}. The test unit is refused when it wins at least (1 − α)(n + 1) games, a win being a strictly higher nonconformity.
- In conformity terms, a win is `strangerThan`, and ties win nothing. Their threshold is `AcceptsPairwise` at k = ⌊αm⌋, with m = n + 1 (`tournament_count_iff`).
- *The pair fit* (`tournFit`) is the fit to the rest, then the scored value and the pair as an unordered multiset of values. `pairScore` of it is their s_approx(zᵢ; D₋ᵢⱼ + {zᵢ, zⱼ}).

**Theorem 1 and its randomized form: coverage ≥ 1 − (2k − 1)/m, so ≥ 1 − 2α.**
- *With a symmetric fit* (`tournament_coverage`, through `crossFit_coverage`): the fit to the rest is symmetric in its training units given its randomness, up to a relabeling, so seeds that travel with the units qualify. This is their row and column exchangeability of the seed matrix.
- *With an order-reading fit* (`tournament_orderReading`, through slice 9's `coverage_orderReading`): it needs no symmetry and no mixing. A chain that sweeps its training units in index order with a fixed seed, with the left-out places closed up, qualifies (`tournFit_orderReading`).
- *What the pair needs.* The fit receives the pair as values, so it cannot know which of the two units is the test (`tournFit_comp`). A fit that hands the pair over as (training, test), or otherwise marks the test, breaks the guarantee: in the checks such a fit refuses the test unit on every data set.
- *Their symmetric-in-the-pair assumption* is one way to meet this condition. A fit with roles (scored unit, opponent) by position would also qualify; roles (training, test) would not.

**Theorem 2, the 2ε-inflated rule under stability** (`coverage_tournament_stable`, `card_refused_by_le`), per data set:
- The test unit is accepted with share at least 1 − (k + a + 2F/(a + 1))/m, over the orderings, for every a.
- F counts the ordered pairs whose score strays from a reference R by more than ε (`strayCount`).
- *The argument.* Where neither side of a game strays, an inflated win needs R(u) < R(j). So a refused unit is among the k + a lowest by R, or has more than a stray games, and a Markov count bounds those.
- *Finding: R can be any reference scores of the data set.* It need be neither symmetric nor a valid score. Their R is the one-out score s_approx(zᵢ; D₋ᵢ + {zᵢ}), and any choice that makes F small serves.
- *Under an exchangeable law,* with each ordered pair straying with chance at most ν (their Assumption 1), E F ≤ m(m − 1)ν. So the miscoverage is at most k/m + a/m + 2(m − 1)ν/(a + 1), for every a.
- *Finding: the constant.* Optimized over a, this is about α + 2√(2ν), below their α + 3√ν: 0.227 against 0.300 at m = 20, ν = 0.01, and 0.0273 against 0.0300 at m = 1000, ν = 10⁻⁴. Their proof fixes the Markov threshold at 2√ν(n + 1) rather than optimizing it.
- The step from the data set to the law is the averaging every CB statement uses, and it is not formalized here.

**What the C# may claim.**
- 1 − 2α, or 1 − (2k − 1)/m exactly, at k = ⌊α(n + 1)⌋, counting strictly higher nonconformity as a win. This needs:
  - each game's pair fit to receive the pair as values, without marking the test;
  - its fit to the rest to be symmetric given its seeds, or to read its training units in index order with the left-out places closed up. It need not mix.
- 1 − α − 3√ν, or the sharper form above, for the 2ε-inflated set, under their stability with the stated ν.

**Checks.**
- The index identity on 5,000 random (m, α, count) triples.
- Coverage over every ordering on 313 (data, k) cases, with an order-reading, non-symmetric pair fit. The least slack was 0: the bound is attained.
- The fit that marks the test unit, which refuses it on every data set.
- The count bound on 3,000 random matrices, with arbitrary references.

## Slice 14 — shift

**Sources, read in TeX:**
- Tibshirani, Barber, Candès and Ramdas 2019 (arXiv 1904.06019v3): §3, Appendix A, Corollary 1, Lemma 3, Theorem 2.
- Barber, Candès, Ramdas and Tibshirani 2023 (arXiv 2202.13415v5): §3–§4, Theorem 2 and its proof, Lemma 1.
- Podkopaev and Ramdas 2021 (arXiv 2103.03323v4): §2.2, Theorem 2, Corollary 1, label-conditional conformal.

Proved on 2026-11-10 in `ConformalValidity.lean`, section `Shift`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 14") ran first.

**The weighted rule, in conformity form.** A unit u is refused when the weight of the units scoring at most its score, itself included, is at most α (`weightLe`). This is their quantile with the candidate's mass at +∞, and ties count for the candidate.

**The weighted quantile lemma** (`sum_weight_refused_le`). The refused units carry weight at most α, for any nonnegative weights: they all score at most the highest refused unit, whose weight-below is at most α.
- *At equal weights it is CB-1* (`weighted_eq_cb1`): accepted iff ⌊αm⌋ < rankLe. This is the commission's check.

**Covariate and label shift** (`coverage_weighted`), per data set.
- *The setting.* The test unit comes from the shifted law, with likelihood ratio w against the training law: w(x) under covariate shift (their Corollary 1), w(y) under label shift (Podkopaev and Ramdas, Theorem 2). Over the orderings of a data set, the test unit lands on unit u with weight w(z_u).
- *The claim.* The rule with the units' normalized weights refuses it with weight at most α of the total.
- *Label-conditional conformal* is CB-1 within each label's stratum (`coverage_stratum`). It holds under any shift of the label proportions and needs no weights.

**The loss with an estimated ratio** (`sum_weight_refused_estimated_le`). Our statement, since the papers give none in finite samples; Podkopaev and Ramdas's Corollary 1 is asymptotic and needs no jump at 1 − α.
- *The claim.* Per data set, the rule run with estimated weights p̂ refuses true weight at most α + Σ_u (p_u − p̂_u)⁺, the total variation between the normalized true and estimated weights of the data set's units. So the coverage loss is at most the expected total variation over data sets.

**Beyond exchangeability** (`refused_le_swap_tv`, their Theorem 2, on a finite law). The law of the data is any law P, and the score is symmetric.
- *The weights* are fixed by position, nonnegative and summing to one, with the test position's weight the largest, which is their wᵢ ≤ 1.
- *The claim.* P(refused) ≤ α + Σ_K q_K · TV(P, P ∘ swap(K, test)).
- *The proof is theirs.* A refused test unit makes unit K strange once K and the test are swapped (`weightLe_swap_le`, using q_K ≤ q_test). The swap costs the total-variation distance (TV), and the weighted lemma bounds the strange units.
- *The form proved* uses TV between the swapped data laws. Their sharper form uses TV between the swapped score vectors, which needs pushforward laws, and is not formalized.

**Checks.**
- The weighted lemma on 4,000 random (scores, weights) with ties, and CB-1 at equal weights.
- Shift on 3,000 data sets, including the loss with estimated weights.
- The nonexchangeable bound on 300 random finite laws.

## Slice 15 — posterior conformal prediction

Zhang and Candès (arXiv 2409.19712v2), read in TeX: §1.2, §2, Theorems 1, 3 and 4, Proposition 2, and the proof of Theorem 3. Proved on 2026-11-10 in `ConformalValidity.lean`, section `PosteriorConformal`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 15") ran first.

**Marginal coverage: exact, and conditional on the draw** (their Theorem 1).
- *Their construction.* The test unit's value sends a draw L ~ κ(· | z) (their L* ~ Multinomial(m, π(X_{n+1}))). Given the data set and L, the test unit is unit u with chance ∝ κ(L | z_u), which is their posterior weight w*.
- *The claim.* The weighted rule with those weights refuses with at most α of that chance, given L (`coverage_pcp`). It is slice 14's `coverage_weighted` with the weight function κ(L | ·). Summed over the draws, the marginal statement follows (`coverage_pcp_marginal`).
- *Estimated π.* Their Proposition 2, with π̂ fitted symmetrically, is the same statement with the fitted kernel inside the equivariant score and weights.

**The per-cluster statement** (`pcp_cluster_gap`; their Theorem 3, on a finite space).
- *Hypotheses.*
  - Given the covariates, the residuals are independent, unit i's law being the mixture Σ_k π_k(i) f_k.
  - The weights are fixed given the covariates, nonnegative and summing to one, with the test's the largest. Their KL weights qualify.
  - The score is symmetric.
- *The claim.* The miscoverage given the covariates is at most α + Σ_K q_K Σ_k |π_k(K) − π_k(test)|.
- *Finding: half their constant.* They state α + 2 Σ_K q_K Σ_k |…|. Through slice 14's `refused_le_swap_tv`, a swap costs Σ_v |m_K(v) − m_t(v)|, the L1 distance of the two units' laws (`swap_tv_prodLaw_le`, Barber et al.'s Lemma 1 without its square term). That is at most the L1 distance of their mixing weights (`mixture_l1_le`), once rather than twice.
- *How loose it is.* In the checks, on 400 finite mixture laws, the largest miscoverage seen used 14% of even the halved gap.
- *Their Theorem 4,* the multiplicative form under positive memberships, is not formalized.

**Checks.**
- Conditional and marginal validity on 1,500 random (data, kernel) pairs.
- The halved cluster gap on 400 finite mixture laws, computed exactly.

## Slice 16 — batch predictive inference

Lee, Tchetgen Tchetgen and Dobriban (arXiv 2409.13990v5), read in TeX: §2, Condition 1, Algorithm 1, Theorem 1 and its proof sketch, the quantile corollary, the selection corollary, and the theorem and corollary on simultaneous bounds. Proved on 2026-11-10 in `ConformalValidity.lean`, section `Batch`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 16") ran first.

**Setting.**
- The n + m units are their positions in score order, S monotone between `bot` and `top`. Any bounds on all n + m scores do, as inf s and sup s do.
- Ties are broken uniformly at random. That leaves the test set uniform over the C(n + m, m) subsets and moves no calibration order statistic. The procedures read only those order statistics, so their sets are the paper's, ties included.
- Every statement holds given the n + m score values, with only the test set's position random.

**The ranks are uniform on H** (`card_testRank_mem`, `card_rankSeqs`).
- The test units' calibration ranks, r_j = a_j − j for the sorted test positions a_j (0-based), are a bijection from the test sets onto H = {0 ≤ r_1 ≤ … ≤ r_m ≤ n}. The inverse is r ↦ {r_j + j}.
- So the rank sequence is uniform on H, and |H| = C(n + m, m). Their sketch states this step as (R_(n+1), …, R_(n+m)) ~ Unif(H).

**The sandwich** (`orderEmbOfFin_lt_cal`, `cal_lt_orderEmbOfFin`, `batch_sandwich`).
- In their 1-based count, S_(r_j − 1) ≤ S^test_(j) ≤ S_(r_j) for every j at once.
- So for any set I of rank sequences, on at least |H ∩ I| of the test sets one r ∈ I sandwiches every sorted test score. Each result below is this sandwich with a particular I.

**Functions of a batch** (`batch_interval`, `batch_coverage`; their Theorem 1).
- *The claim.* For h monotone, one r ∈ I has h(S_(r − 1)) ≤ h(sorted test scores) ≤ h(S_(r)).
- *Their interval.* Take I = {q_L ≤ h̃ ≤ q_U}. By the definitions of Q′_β and Q_{1−γ}, at most β of H lies below q_L and at most γ above q_U, so |H ∩ I| ≥ (1 − α)|H| with α = β + γ. Their B_L is a min over {h̃ ≥ q_L} and B_U a max over {h̃ ≤ q_U}, both containing I. So [B_L, B_U] covers on at least 1 − α of the test sets (`frac_of_card`).
- *Their h.* Their h, defined on sorted vectors, extends to v ↦ h(v↑), which is monotone by their Condition 1.

**Simultaneous sets** (`batch_simultaneous`; their theorem on joint bounds).
- With I the box w_j ≤ r_j ≤ q_j, every sorted test score lies in [S_(w_j − 1), S_(q_j)] on at least |H ∩ I| of the test sets. A coordinate left free takes 0 and n.
- Their condition F_{n,m}(w; q) ≥ (1 − α) C(n + m, m) is |H ∩ I| ≥ (1 − α)|H|.
- *A typo in their count.* Their closed form for L(ρ) has its product over j = 1, …, n. It holds with the product over consecutive pairs, j = 1, …, l − 1, which the checks confirm.

**Selection with false claims controlled** (`batch_exceed`, `false_claim_iff`; their quantile and selection corollaries).
- *The claim.* With I = {r_ζ ≤ k}, at most m − 1 − ζ test scores exceed S_(k), with ζ 0-based.
- *Their rule.* Their selection score is s = μ̂ · 1{y ≤ c}. At a threshold T ≥ 0, a false claim (μ̂ > T while y ≤ c) is exactly a score above T. So with T̂ = S_(q_η) and 1-based ζ = m − η, the number of false claims is at most η on at least 1 − α of the test sets.
- *The quantile.* q_η is the (1 − α) quantile of the law of r_(ζ). The checks confirm their closed form for that law.

**Not formalized.**
- The closed forms for |{r ∈ H : r_ζ = k}| and L(ρ) as Lean statements; the checks confirm them.
- The upper bound 1 − α + ε in the quantile corollary.
- Their propositions on NP-hardness, on optimality for the quantile, and on k-FWER through p-values.

**Checks.**
- The bijection onto H and |H| = C(n + m, m), for n ≤ 5 and m ≤ 4.
- Their law of r_(ζ) and L(ρ) in closed form, on 1,626 (n, m, t, ρ).
- On 400 (n, m, scores, h̃) cases, half with tied scores, over every test set:
  - the sandwich;
  - [B_L, B_U], built as in their Algorithm 1, covering at least 1 − β − γ;
  - the interval, joint-box and exceedance counts each at least their rank counts.
- Their selection rule on 180 data sets with scores tied at zero: at most η false claims on at least 1 − α of the test sets.

## Slice 18 — the guarantee read out

Snell and Griffiths (arXiv 2502.13228v2), read in TeX: §3–§4, Propositions 3.1–3.2, Theorems 4.1 and 4.3, Lemma 4.2, Corollary 4.4, the appendix proofs, and §5.1. Gressmann et al. (arXiv 1801.00753v3), read in TeX: §2.5, Definitions 2.2–2.3, Lemma 2.4, Propositions 2.5–2.7. Proved on 2026-11-10 in `ConformalValidity.lean`, sections `LossPosterior` and `ProperScores`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 18") ran first. Checked on 2026-11-13 against the sources the first reading lacked:
- Vlassis and Thomas (arXiv 2607.08415v1), read in TeX: §1 and Theorem 1, with the two-point validity and its proof.
- Angelopoulos et al., *Conformal risk control* (arXiv 2208.02814v4), read in TeX: §2, Theorems 1 and 2, and Theorem 1's proof.
- Gneiting and Raftery (JASA 2007), read from the PDF: Examples 1 and 3, and the quadratic score (18).

**The posterior over the loss.**
- *Their construction.* Take n calibration losses ℓ_(1) ≤ … ≤ ℓ_(n), set ℓ_(n+1) = B, and let L⁺ = Σ U_i ℓ_(i) with U ~ Dir(1, …, 1).
- *Its exact law, for small n.* The checks compute it exactly: P(L⁺ > x) is the divided difference of (t − x)_+^n on the knots, confluent where they tie.
- *Its mean is conformal risk control* (`crc_orderings`; their Proposition 3.2, through Angelopoulos et al.'s Theorem 1).
  - E(L⁺) = (Σ ℓ_i + B)/(n + 1). The rule λ_crc takes the first threshold with E(L⁺) ≤ α.
  - Over the orderings of any data set, the test unit's loss at that threshold averages at most α. So the posterior mean is calibrated over replicated calibration sets.
  - *Hypotheses.* Losses are nonincreasing in λ on a grid, at most B, and at most α at the last threshold.
  - *Against Angelopoulos et al.* Their Theorem 1 has these hypotheses, and its proof runs the same route: an oracle threshold over all n + 1 losses sits below λ̂, and the test loss is uniform over the multiset. On a grid, their right-continuity is not needed.

**Finding: their Theorem 4.3 does not hold as stated.**
- *The gap.* Its proof replaces p(t_{1:n} | ℓ_{1:n}) with the Dirichlet law of the spacings. That is their law before ℓ is seen, not after.
- *A counterexample.*
  - Put a point-mass prior on the quantile function that is 0 on [0, 1/10] and B = 1 above, so J = 9/10. It is consistent with all-zero losses.
  - Its posterior gives Pr(L ≤ b | ℓ) = 0 for b < 9/10, while Pr(L⁺ ≤ b) = 1 − (1 − b)^n > 0.
  - So inf_π Pr(L ≤ b | ℓ) ≥ Pr(L⁺ ≤ b) fails.
- *What stands.* Their Theorem 4.1, the right-rectangle bound given the quantile levels.

**The β-bound's calibration is frequentist, and known.**
- *Identification.* b*_β is Gaffke's (2005) bound.
  - Vlassis and Thomas's statistic is K(x) = P{Σ x_i D_i ≤ 1}, with D ~ Dir(1, …, 1) over the n values and an added zero.
  - With X_i = (B − ℓ_i)/(B − J), P(L⁺ ≤ J) = 1 − K(X). The checks confirm this identity exactly on 300 cases.
  - Learned-Miller and Thomas rediscovered the bound as a confidence bound on [0, 1].
- *Status: proved.* The calibration over replicated calibration sets, P(J ≤ b*_β) ≥ β, is Vlassis and Thomas's Theorem 1. Formalized at every n since 2026-11-14 (`gaffke_calibration`, under "Slice 18's follow-up").
  - b*_β misses J exactly when K(X) < 1 − β, which has chance at most 1 − β.
  - Gaffke conjectured it, reduced it to mean-one two-point laws, proved n = 2, and checked it numerically to n = 15. Learned-Miller and Thomas proved Bernoulli and half-Bernoulli laws.
- *Wider than their footnote.* The theorem needs the losses independent and at most B, and needs neither identical laws nor continuity. Snell and Griffiths's Lemma 4.2 needs continuity; the calibration does not, and the checks' laws are discrete.
- *One unit* (`gaffke_one`, our proof).
  - With one loss, b*_β = ℓ + β (B − ℓ).
  - The miss set's mass P satisfies P (B − J) ≤ (1 − β)(B − J) by a Markov-type step, strictly when P > 0. So P ≤ 1 − β.
- *How tight it is.* The checks find exact calibration on 150 loss laws, n ≤ 5, β ∈ {1/2, 3/4, 9/10, 19/20}. The least slack is 0.0005, so the bound is nearly tight.

**Their HPD rule** (`fail_le_fixed`).
- λ^β_hpd takes the first threshold with b*_β(λ) ≤ α.
- Its risk exceeds α only when the fixed threshold λ₀, the last with risk above α, has b*_β(λ₀) < J(λ₀). So its failure chance is at most 1 − β, by the calibration at λ₀.
- Their §5.1 shows this in simulation and does not state the step from a fixed λ to a chosen one.
- End to end, with the calibration at every n: `highProbability_risk` (2026-11-14), under "Slice 18's follow-up".

**Proper scores** (Gressmann et al.).
- *The log-loss* (`logLoss_proper`, `logLoss_strict`; their Proposition 2.5 (i)).
  - It is Gibbs' inequality on finitely many outcomes, strict at q ≠ p for q a law.
  - q must be positive where p is: Lean's log 0 is 0, where their loss is +∞. Gneiting and Raftery's log score (their Example 3) is strictly proper relative to the measures dominated by the reference measure, which is this hypothesis on finitely many outcomes.
- *The integrated squared loss* (`gneitingLoss_gap`, `gneitingLoss_strict`; their Lemma 2.4 (ii) and Proposition 2.5 (ii)). The gap is Σ (p − q)².
  - Against Gneiting and Raftery, the loss is minus their quadratic score (18). Their Brier score (Example 1) is 2p_i − Σ p_j² − 1, and the gap is its Bregman divergence.
- *Conditional predictions* (`proper_conditional`; their Proposition 2.7).
- *Finding: two sign slips in their proof of Proposition 2.5 (i).* Jensen is written log E Z ≤ E log Z, and Gibbs as E[log p_Y(Y)] ≤ E[log p(Y)]. Both hold for the loss −log, which is the reverse of what is written. The statement itself stands.

**Not formalized.**
- The calibration for n ≥ 2: Vlassis and Thomas's Theorem 1, whose proof is a chain-measure induction. *Formalized since (2026-11-14):* `gaffke_calibration` and `highProbability_risk`, under "Slice 18's follow-up".
- Angelopoulos et al.'s Theorem 2, the lower bound α − 2B/(n + 1) for i.i.d. losses without jumps.
- Their Theorem 4.1 on [0, 1], and Lemma 4.2 (Dirichlet spacings).
- The density case of the proper scores, and their Proposition 2.6 (the log-loss is the only strictly local, strictly proper loss).

**Checks.**
- The exact law of L⁺ on 393 knot sets: its mean, its ends, and the distinct-knot formula.
- At n = 1, b*_β = ℓ + β(B − ℓ).
- The calibration over replicated calibration sets, as above. At n = 1 it is checked as the Lean states it.
- P(L⁺ ≤ J) = 1 − K(X) for Vlassis and Thomas's K, exactly, on 300 cases.
- Their Theorem 4.3 failing under the point-mass prior.
- Conformal risk control over orderings on 300 data sets.
- The HPD rule exactly, end to end, on 100 laws: it fails at most where the fixed threshold does, within 1 − β (worst failure 0.15).
- The squared loss's gap exactly, and Gibbs in floats, on 300 law pairs.
- The written directions in their proof failing at once.

## Slice 14's follow-ups

Commissioned on 2026-11-13 (the brief's top section), each to settle a C# decision named in slice 14's C# report.

### Follow-up 1 — the weighted predictive distribution's PIT

Sources: Tibshirani, Barber, Candès and Ramdas 2019 (arXiv 1904.06019v3) and Podkopaev and Ramdas 2021, as read for slice 14. Proved on 2026-11-13 in `ConformalValidity.lean`, section `WeightedPit` beside `Shift`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 14 follow-up 1") ran first.

- *The distribution.* At the test unit's value, the weighted predictive distribution puts weight A strictly below it and B tied with it, its own included (`weightLt`, `weightEq`).
  - With the tied mass randomized, the PIT is T = A + V B, with V uniform (`wPitCdf`).
  - With ties counted for the candidate (V = 1), it is `weightLe`.
- *Uniform given the data set* (`sum_mul_wPitCdf`). If the test unit is unit u with chance q_u and Σ q = 1, then P(T ≤ x) = x for x ∈ [0, 1]. The tie levels' intervals [A, A + B) tile [0, 1] (`levels_tile`).
- *Over orderings, split and full* (`pit_weighted`). With the true ratio w and an equivariant score, the PIT is exactly uniform over the orderings, each weighted by w(z_{σt}).
  - w is w(x) for covariate shift and w(y) for response shift.
  - The split form is the score S z = s ∘ z.
- *Unrandomized:* `weightLe`'s super-uniformity at every level is `coverage_weighted` (slice 14).
- *An estimated ratio* (`wPit_estimated`, our statement).
  - The test unit lands by the true normalized weights p, and the PIT is read with the estimates p̂.
  - Then P(T ≤ x) is within TV(p, p̂) = Σ (p − p̂)⁺ of x, from either side.
  - The checks find the bound attained: the largest share of TV used is 1.000.
- *For the C#:* the distribution may take shift weights. Its PIT is uniform with the true ratio, and moves by at most the weights' TV with an estimated one.

**Checks.**
- The tiling on 3,000 random cases (scores with ties, weights with zeros), with `weightLe`'s super-uniformity.
- The full form over every ordering, with an equivariant score, on 150 data sets.
- The estimated ratio's bound on 3,000 cases.

### Follow-up 2 — the cross-fitted forms under weights

Barber, Candès, Ramdas and Tibshirani 2023 (arXiv 2202.13415v5), read in TeX: Appendix B (the nonexchangeable jackknife+, their Theorem 5 and its proof), Appendix D (cross-conformal), and their Lemma 2 (Lei's) with its proof. Proved on 2026-11-13 in `ConformalValidity.lean`, section `WeightedCrossFit`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 14 follow-up 2") ran first.

- *The weights are fixed by position.* They are Barber et al.'s w̃, summing to one with the test position's the largest, not a likelihood ratio.
- *Weighted strange points* (`strangeWeight`, `sum_weight_strange_le`; their Lemma 2, Lei's).
  - Unit u is refused when the units it is stranger than carry at least 1 − α of the weight.
  - Whatever the comparisons, the refused units carry at most 2α of it, for α ≥ 0. Among the refused units each comparison goes at most one way, so they are stranger than at most half their own weight's square.
  - *A hypothesis the lemma needs.* With α < 0 nothing is refused, but the bound 2αW is negative, so α ≥ 0 is stated; their α lies in (0, 1).
  - *Against the counts.* Unweighted, `card_refused_pairwise_le` gives 2k − 1, sharper by one from integrality.
- *The swap step* (`strangeWeight_swap_ge`). A refusal at t in z is a refusal at K in z ∘ swap(K, t), because the test position's weight is the largest. This is the step w̃_{π_K(i)} ≤ w̃_i of their proof.
- *Beyond exchangeability* (`crossFit_weighted_swap_tv`; their Theorem 5 and Appendix D). On any finite law, with any pair-equivariant score, the weighted cross-fitted rule refuses with chance at most 2α + Σ_K q_K TV(P, P ∘ swap(K, t)).
- *Exchangeable data* (`crossFit_weighted_exchangeable`). The swaps cost nothing, giving 2α, whatever the fixed weights, so long as the test's is the largest.
- *The jackknife+* (`JackRefusesW`, `strangeWeight_of_jackRefusesW`, `jackknifePlus_weighted`).
  - Y outside their weighted interval, with the test unit's mass at ±∞, is a weighted cross-fitted refusal.
  - So with a symmetric fit, the jackknife+ misses with chance at most 2α plus the swap distances: their Theorem 5 with a symmetric algorithm.
- *The tournament correction under weights* (`tournament_weighted`). It holds: the tournament is the cross-fitted rule on its pair score, so it inherits 2α plus the swap distances. Seeds that travel with the units are part of each unit's data.
- *Likelihood-ratio weights* (Prinster, Liu and Saria 2022, JAWS, arXiv 2207.10716v2, read in TeX: §3.1, Theorem 1 and its appendix proof).
  - *The setting* (`strangeWeight_comp`, `crossFit_likelihoodRatio`). The weights w(z_u) travel with the units, the test unit's own included, and over the orderings of a data set the test unit lands on unit u with weight w(z_u).
  - *The claim.* Weighing each ordering so, the weighted cross-fitted rule refuses at most 2α of the weight, with no TV term. Given the bag, the test unit is a strange unit with chance equal to their weight, at most 2α by Lei's lemma.
  - *Which shifts.* This is covariate shift, w(x), their case, and label shift, w(y), in candidate form.
  - *The forms.* JAW, their jackknife+ (`jackknifePlus_likelihoodRatio`), and the tournament (`tournament_likelihoodRatio`) follow.
  - **Finding: their Step 2 argues from exchangeability, which covariate shift breaks.** It claims A^w =d ΠA^wΠᵀ for every permutation matrix Π, so that the test point is as likely strange as any training point. It then multiplies those probabilities by the random weights p_j outside the expectation. Their Theorem 1 stands by the weighted-exchangeability argument, given the bag, which is what `crossFit_likelihoodRatio` proves.
- *K-fold CV+ under weights* (our statements; Barber et al. 2023 and JAWS give the leave-one-out forms only). Folds of training units, each with a fit leaving out it and the test unit; the weighted rule's total includes the test unit's weight.
  - *A refusal is a refusal in each fold* (`fold_refused_of_refused`). If the units the test unit beats carry at least 1 − α of the weight W_tot, then in any fold k₀ the units scoring at most the test unit under that fold's fit, itself included, carry at most α W_tot.
  - *Likelihood-ratio weights, by one fold* (`cvPlus_oneFold_likelihoodRatio`).
    - Given the other folds and the bag of fold k₀ with the test unit, the miss chance is at most min(1, α W_tot/(W_{k₀} + w_t)). Report the fold of largest weight.
    - At equal weights this is slice 9's ℓ/(m + 1); with two balanced folds it is about 2α.
    - The checks find the bound attained.
  - *Merged folds with position weights* (checked; the Lean is the next item).
    - A refusal gives Σ_k λ_k π_k ≤ c. Here π_k is fold k's weighted p-value, λ_k = (W_k + w_t)/(W + K w_t) is fixed, and c = (α(W + w_t) + (K − 1) w_t)/(W + K w_t).
    - Averaging super-uniform p-values with fixed weights costs a factor 2, so P ≤ 2c.
    - At equal weights 2c is slice 9's merged bound.
  - *Not covered.*
    - A merged bound with likelihood-ratio weights: the λ_k then depend on the data, and the averaging argument needs them fixed. Proved later the same day by another route (follow-up 2e).
    - The tournament across folds under weights.
- *For the C#:*
  - the jackknife+, the n-fold cross-conformal set and the tournament can take weights fixed by position, at 2α plus the swap distances, with the test position's weight the largest;
  - they can take likelihood-ratio weights at 2α;
  - K-fold CV+ can take likelihood-ratio weights, reporting min(1, α W_tot/(W_k + w_t)) from the fold of largest weight.

**Checks.**
- Lei's lemma on 4,000 random comparison matrices with ties. The largest share of the bound used is 0.667.
- The swap step on 2,925 (data, K).
- The bound on 75 finite laws with the leave-two-out and tournament pair scores, exchangeable laws within 2α exactly.
- The weighted jackknife+ refusing only where the cross-fitted rule does, on 2,000 data sets.
- Likelihood-ratio weights over every ordering at most 2α, on 120 data sets with both pair scores. The largest share of the bound used is 0.476.
- K-fold CV+:
  - the one-fold likelihood-ratio bound on 68 (data, fold) cases over every reordering of the fold and the test unit, attained (share 1.000);
  - the merged position-weight bound 2c on 48 cases over every ordering, using at most 0.284 of it, with equal weights giving slice 9's merged bound.

### Follow-up 3 — CB-3's margin and inflated rule under weights

CB-3's statements with the weighted rule in place of counts. Proved on 2026-11-13 in `ConformalValidity.lean`, section `WeightedMargin`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 14 follow-up 3") ran first.

- *The weighted fragile units* (`weightShift`, `WFragile`).
  - A fragile unit is accepted by the weighted rule on the exact score (α < `weightLe`).
  - But at most α of the weight sits at least ε_j + ε_u below it, itself counted.
- *The deterministic part* (`waccepts_or_fragile`). Estimates within their bounds refuse an exactly accepted unit only when it is weighted-fragile.
- *The weighted inflated rule* (`WAcceptsInflated`, `wAcceptsInflated_of_accepts`).
  - It accepts u when more than α of the weight has an estimate, lowered by its bound, at most u's raised by its own.
  - It keeps every unit the exact score accepts.
- *The estimated band* (`WFragileEstimated`, `wFragileEstimated_of_fragile`). It holds every weighted-fragile unit, read from the estimates and bounds alone (2ε_j + 2ε_u).
- *Over orderings, with shift weights w(z_u) and the bounds travelling with the units:*
  - the estimate's weighted rule refuses at most α + M, where M is the weighted-fragile units' share of the weight (`coverage_close_weighted`);
  - the weighted inflated rule refuses at most α (`coverage_inflated_weighted`).
- *Its dominance from any ordering's estimates* (`sum_wfragile_le_estimated`). The exact score's weighted-fragile units weigh at most the band read from any ordering's estimates, with the weights travelling. So the engine can report α + M̂ with M̂ read from what it holds.
- *For the C#:* weighted full-step sets can report their Monte Carlo error. The margin is α + M̂, with M̂ the estimated band's weight, and the inflated rule keeps 1 − α.

**Checks.**
- The deterministic statements on 4,000 random cases (scores, estimates within their bounds, bounds, weights).
- Over every ordering of 116 data sets, with an equivariant score and per-ordering estimates:
  - the estimate refuses at most α + M, and the inflated rule at most α;
  - M is at most any ordering's estimated band.
  - The largest excess of the estimate's refusals over α is 0.133, within M.

### Follow-up 2, continued — K-fold CV+ under position weights, by merging the folds

Our statements, the weighted form of slice 9's merged bound, at the user's ask on 2026-11-13. Proved on 2026-11-13 in `ConformalValidity.lean`, section `WeightedMerge`, PROOF-CLOSED, with the standard three axioms only. The checks ran first: section "Slice 14 follow-up 2c", and part (b) of follow-up 2's.

- *The reduction* (`cv_merged_of_refused`). A weighted CV+ refusal makes Σ_k λ_k π_k ≤ c, where:
  - π_k = (W_k⁻ + w_t)/(W_k + w_t) is fold k's weighted p-value, W_k⁻ its unbeaten weight;
  - λ_k = (W_k + w_t)/(W + K w_t), fixed when the weights are fixed by position;
  - c = (α(W + w_t) + (K − 1) w_t)/(W + K w_t).
- *Averaging with fixed weights* (`merge_le_two`).
  - If each π_k ≥ 0 is super-uniform at the levels 0 ≤ x < 2c, the average is at most c with chance at most 2c: the arithmetic mean's factor 2, as in slice 9's merged bound. (As first committed, the levels included negative x, which made the three theorems vacuous; corrected in follow-up 2e.)
  - The proof goes through a finite layer cake: Σ_v m_v (b − v) ≤ b²/2 when the mass up to each level is at most the level (`sum_mass_mul_sub_le`). So E(b − π)⁺ ≤ b²/2 for a super-uniform π (`sum_mul_max_sub_le`).
- *The bound* (`cvPlus_merged_weighted`).
  - With weights fixed by position and each fold's p-value super-uniform, CV+ refuses with chance at most 2c.
  - At equal weights 2c is slice 9's 2(ℓ + K − 1)/(n + K).
- *A hypothesis taken, not derived here.* Each fold's super-uniformity is a hypothesis of the theorem. On exchangeable data with the test position's weight the largest in the fold, it is Barber et al.'s Theorem 2 within the fold. `refused_le_swap_tv` states that over all units (`univ`), so applying it needs the fold, with the test unit, as the index type.
- *Not covered.* A merged bound with likelihood-ratio weights: the λ_k then depend on the data, and the averaging argument needs them fixed. For likelihood ratios the one-fold bound in follow-up 2 serves. Proved later the same day by another route (follow-up 2e).
- *For the C#:* K-fold CV+ under position weights can report 2c.

**Checks.**
- The end-to-end bound 2c on 48 cases over every ordering (follow-up 2's block, part (b)).
- The finite layer cake on 3,000 random (levels, masses).
- The averaging lemma on 300 cases of super-uniform weighted-rank p-values. The largest share of 2c used is 0.500.

### Follow-up 2, continued — the open K-fold items: the tournament across folds under likelihood-ratio weights, and the merged bound

At the user's ask on 2026-11-13 ("lets work on the remaining still open"). Proved on 2026-11-13 in `ConformalValidity.lean`, section `WeightedAcross`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 14 follow-up 2d") ran first.

**The tournament across folds, weighted** (our statements; slice 9's CV-b with weights).
- *The setting.*
  - Units are (fold, slot). The test unit's fold holds it and further units from the training law, which the interval never reads.
  - Units play only across folds: u beats j when its conformity, under the fit leaving out both their folds, is below j's.
  - u is refused when the units it beats carry at least 1 − α of the weight outside its fold, plus its own. For the test unit that is weighted CV+, its total the training weight plus the test unit's own.
- *The counting core* (`strangeAcrossWeight`, `sum_weight_refused_across_le`).
  - Whatever the comparisons, the refused units weigh at most 2αW + W_max, for 0 ≤ α ≤ 1, with W_max bounding every fold's weight.
  - Each refused unit beats all but α W of the weight across folds. Among the refused units, each comparison across folds goes at most one way, and the same-fold pairs weigh at most W_max times the refused weight.
  - At m = 1 this is Lei's lemma's 2αW plus the heaviest unit. Unweighted, slice 9's count, 2ℓ + m − 2, is sharper by integrality.
- *The averaging* (`sum_prodPerm_apply`). `ratio_prodPerm` for any real function: averaging over the fold and slot permutations is averaging over the units.
- *The coverage* (`coverage_across_likelihoodRatio`).
  - Given the bag, the test unit is u with chance w(z_u)/W and the folds are uniform. These are the fold and slot permutations, each weighed by the test unit's likelihood ratio.
  - The pair score is equivariant under them: a symmetric fit to the other folds.
  - Then weighted CV+ refuses with at most 2α + W_max/W of the weight.
- *Where it serves.* With K + 1 folds of similar weight, W_max/W is about 1/(K + 1), so this is the many-folds bound under likelihood ratios. The one-fold bound (`cvPlus_oneFold_likelihoodRatio`) serves two folds.
- *What it reads.* W_max includes the test unit's imagined fold-mates. It is bounded by m times the largest weight the ratio can take; a ratio bounded by w̄ gives W_max ≤ m w̄.

**The merged bound under likelihood-ratio weights: open here; proved in follow-up 2e, with the test unit's weight replaced by a cap read from the bag.**
- *Why the proof fails.* Under position weights the folds' p-values are averaged with fixed weights λ_k (`cvPlus_merged_weighted`). Likelihood ratios make λ_k and c depend on which unit is the test, and the averaging argument needs them fixed.
- *What the search found.* `conformal-searches.py --weights` anneals over every symmetric fit in an abstract model, 24 configurations with K ≤ 4, n + 1 ≤ 7 and weights 1–20. It found nothing past slice 9's unweighted merged bound or past E[2c]: the largest shares were 0.326 and 0.292.
- *Status.* That is evidence, not a statement the layer may report. Under likelihood ratios the layer reports the one-fold bound for few folds and the tournament bound for many.

**Checks.**
- The counting core on 5,000 random fold-structured comparisons, with ties, weights with zeros and α up to 1. The largest share of the bound used is 0.560.
- The coverage over every fold and slot permutation, weighted by the test unit's likelihood ratio, on 59 data sets. The largest share of the bound used is 0.278.

### Follow-up 2e — the merged bound under likelihood-ratio weights, and the averaging lemmas' hypotheses corrected

At the user's ask on 2026-11-13 ("im asking for more work on it"). Proved on 2026-11-13 in `ConformalValidity.lean`, section `WeightedMergeLR`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 14 follow-up 2e") ran first.

**The merged bound under likelihood-ratio weights** (our statements). This closes the open item above, with the test unit's weight replaced by a cap read from the bag.
- *What the fixed weights were for.* Markov's thresholds can be θ W_k, proportional to each block's weight, where a block is a fold with the test unit. W_k is fixed given the block as a set. Given that, the test unit is u with chance w_u/W_k, so the block's weight V_k scoring at most the test unit is super-uniform there.
- *The weighted layer cake* (`sum_weight_mul_max_sub_le`). Σ_u q_u (b − weightLe(u))⁺ ≤ b²/2, for any scores and q ≥ 0.
- *One block* (`sum_permsOn_mul_max_sub_le`). Over the block's reorderings, each weighed by the test unit's w, E(θW_k − V_k)⁺ ≤ θ²/2 · E W_k.
- *The refusal's slack* (`sum_blocks_filter`, `sum_blocks`, `blocks_slack_of_refused`).
  - The blocks' weights sum to W + (K − 1)w_t.
  - A refusal gives Σ_k (θW_k − V_k) ≥ (θ − α)W − (1 − θ)(K − 1)γ for θ ≤ 1, with γ at least the test unit's weight.
  - A refusal needs w_t ≤ αW, so γ need bound only the weights at most αW.
- *Markov's form* (`cvPlus_merged_likelihoodRatio_markov`).
  - Over every ordering of the bag, each weighed by the test unit's likelihood ratio: P(refused) · ((θ − α)W − (1 − θ)(K − 1)γ) ≤ θ²/2 · (W + (K − 1)e).
  - Here e = Σw²/W is the test unit's mean weight.
  - Each block's average runs over the cosets of its reorderings.
- *The bound* (`cvPlus_merged_likelihoodRatio`).
  - Take θ = 2L/N, with L = αW + (K − 1)γ, N = W + (K − 1)γ and γ the bag's largest weight. Then P(refused) ≤ min(1, 2c(γ)), with c(γ) = L/N.
  - This is follow-up 2's position-weight 2c with the test unit's weight replaced by the bag's largest.
  - At equal weights it is slice 9's merged bound in its α form (`cvPlus_merged_alpha`).
- *The refined form* (`cvPlus_merged_likelihoodRatio_refined`).
  - Take γ the largest weight at most αW, with 2L ≤ N. Then P(refused) ≤ 2L(W + (K − 1)e)/N².
  - It helps when a unit is too heavy to be refused.
- *What it reads.* W and the bag's weights (the training units' and the test point's own) are observed, so the layer computes both bounds per test point. The bound holds given each bag, so the miss rate over bags is at most its average.
- *Still open.* Whether E[2c(w_t)] bounds it, with the test unit's own weight in c. The two are not ordered: the proven bound is often the smaller (see the checks).
- *Which bound serves.* Miscoverage bounds at α = 0.1, for weights of mean 1 with the ratio capped at w̄; the one-fold bound is at the fold's mean weight.

| w̄ | n | K | one fold | tournament | merged |
|---|---|---|---|---|---|
| 1 | 20 | 2 | 0.191 | 0.676 | 0.282 |
| 1 | 100 | 10 | 0.918 | 0.299 | 0.347 |
| 1 | 1000 | 5 | 0.498 | 0.400 | 0.207 |
| 3 | 20 | 2 | 0.191 | 1.629 | 0.425 |
| 3 | 100 | 5 | 0.481 | 0.794 | 0.391 |
| 3 | 100 | 10 | 0.918 | 0.497 | 0.580 |
| 3 | 1000 | 10 | 0.991 | 0.500 | 0.247 |

  Under likelihood ratios the merged bound serves many units, the tournament many folds of few units, and the one-fold bound two folds.
- *For the C#:* K-fold CV+ under likelihood-ratio weights can report the smaller of min(1, 2c(w_max)) and the refined form. Here w_max is the bag's largest weight, the test point's included. The layer reports the least of these, the one-fold bound and the tournament bound.

**The averaging lemmas' hypotheses, corrected.** This was a defect in follow-up 2c (2928b06).
- `sum_mul_max_sub_le`, `merge_le_two` and `cvPlus_merged_weighted` assumed P(π ≤ x) ≤ x for every x below the bound, negative x included.
  - No law meets that: at x = −1 the left side is at least 0.
  - So the three theorems held vacuously.
- They now take π ≥ 0, and super-uniformity for 0 ≤ x only. That is what a p-value gives and what the proofs use. `cvPlus_merged_weighted` derives π ≥ 0 from its weights.
- The conclusions are unchanged, and so are the numeric checks behind them.
- Found while reusing the layer cake here. A search of the module for super-uniformity hypotheses over real levels found no other.

**Checks.**
- The weighted layer cake on 4,000 random cases. The largest share used is 0.789.
- The refusal's slack on 1,288 refusals among 4,000 random fold-structured cases, at every θ in steps of 0.1 and three choices of γ.
- End to end over every ordering on 84 bags. The fits are the mean of the other folds and the abstract model of symmetric fits.
  - Every bag stays within both bounds, and within the Markov form at every θ in steps of 0.1.
  - The largest share used is 0.400.
  - Equal weights give slice 9's merged bound.
- The adversarial search (`conformal-searches.py --weights`), rerun against the proven bound.
  - The bound holds in all 24 configurations, and the worst case found uses at most 0.322 of it.
  - It is below E[2c] in 17 of them, above it in one, and both are trivial in 6. The refined γ makes the difference where a unit is too heavy to be refused. Neither bound dominates the other.

## Slice 18's follow-up — the calibration at n units

Asked on 2026-11-14: extend `gaffke_one` to n units by Vlassis and Thomas's chain-measure induction, so the Lean states what `HighProbability` reports. Source: Vlassis and Thomas (arXiv 2607.08415v1), read in TeX in full (§1–§3: Theorem 1, Lemmas 2 and 5–7, Proposition 3, Remark 4). Proved on 2026-11-14 in `ConformalValidity.lean`, section `VlassisThomas` after `LossPosterior`, PROOF-CLOSED, with the standard three axioms only. The checks (section "Slice 18 follow-up: the calibration at n units") ran first.

**The model.**
- *Exponentials as laws on the line.* `unitExp` is the unit exponential by its density, `expLaw r` the law of rE, and `sumLaw m` the law of Σ_{r∈m} r E_r, built by convolution. Independence is the convolution, so there is no product space. `stat m` is P(Σ r E_r ≥ 0).
- *Their statistic* (`gaffkeK`). K(x) = P(E₀ + Σ (1 − x_i) E_i ≥ 0), their (1): the Dirichlet weights are normalized exponentials.
- *The posterior's distribution function* (`lplusCdf`). Snell and Griffiths's P(L⁺ ≤ x) is P((x − B)E₀ + Σ (x − ℓ_i)E_i ≥ 0).
- *Gaffke's bound* (`gaffkeBound`). b*_β = inf{x : P(L⁺ ≤ x) ≥ β}, the least such x, which `LossBound.Gaffke` bisects for.

**The analytic core** (their Lemma 5, by another route).
- *The race identity* (`stat_race`, through `race_Ici`).
  - (a + b) P(Y + aE − bE′ ≥ 0) = a P(Y + aE ≥ 0) + b P(Y − bE ≥ 0), for any independent Y.
  - It is the one computation, from the unit exponential's memorylessness (`lintegral_unitExp_Ici`) and its Laplace transform (`lintegral_unitExp_exp`).
- *The same-sign identities* (`stat_same_pos`, `stat_same_neg`). (a − b) P(Y + aE + bE′ ≥ 0) = a P(Y + aE ≥ 0) − b P(Y + bE ≥ 0), and the same for −a and −b. They follow from the race identity with an auxiliary exponential, by algebra.
- *The windows* (`sumLaw_hi_add_winHi`, `sumLaw_lo_eq`). With X_k high, the statistic loses u = E[e^{−Z/d}; Z ≥ 0]. With X_k low, it gains v = E[e^{Z/c}; Z < 0]. These are their u and v.
- *The flux identities* (their (19)). They follow from the race and same-sign identities, with no Stein identity and no derivative.
- *The likelihood-ratio inequality* (their (22) and (15); `lrAbove_edge`, `window_cross`).
  - Y + aE lies above Y − bE in the likelihood-ratio order, tested on ordered pairs of nonnegative functions.
  - The proof adds Y's exponentials one at a time (`LRAbove.conv`). Memorylessness splits each smoothed integral where its two arguments meet (`smooth_ordered`).
  - So it needs no densities and no theorem that convolution keeps log-concavity: each step adds one exponential.
- *θ and η* (their (15), (14) and (31); `windowShare_edge`, `transfer_nonneg`). θ = u/(u + v) falls along an edge. When c ≤ a, η = (1 − θ₋)(A₊ − A₋) − p(F₊ − F₋) ≥ 0, through (c + d)η = (a − c) w₊ (θ₊ − θ₋).

**The chain-measure induction** (their §2 and Lemma 6).
- *The objects.* Chains are sequences of finite sets (`IsChain`), with their measures' means (`chainMean`). `cubeK` is their K_m(S), the two-point statistic, and `bern` the independent law on the high sets.
- *The insertion's algebra* (their (27)–(29); `insertion_algebra`, `chain_insert_sum`).
  - λ_J = θ_J − θ_{J+1} is a law.
  - Σ_J λ_J E_{C^J} g − E_ν̂ g = Σ_j η_j (g(H_j) − g(C_j)) ≥ 0, so the best insertion beats the independent reveal.
- *The local insertion* (`local_insertion`).
  - F, A and B run along the chain, closed by their sentinel G_k.
  - The new unit has the least γ.
  - A unit with γ = 0 is never high and goes last, so their deletion step is not needed.
- *The induction* (`chain_dominates`). Strong induction over the placed units, inserting the least γ last; the unplaced units stay independent.
- *Two-point validity* (`twoPoint_valid`, their Proposition 3). The rejected states are a terminal segment of the final chain, and their masses telescope (`chain_terminal`).

**From two-point systems to b*_β.**
- *The two-point mixture* (`twoPoint_mixture`; their Lemma 7, for a finite law). The pair x_a < 1 < x_b has weight (x_b − x_a) p_a p_b / M, and p(x = 1) goes on the point 1.
- *Theorem 1 for finite laws* (`gaffkeK_valid`, through `sum_prod_twoPoint`, the bridge to the cube). Rescaling to mean one is as theirs (`gaffkeK_anti`). For any real α, P(K(X) ≤ α) ≤ max(α, 0).
- *A miss puts K low* (`gaffkeK_le_of_miss`).
  - If b*_β < J < B, then K(X) = 1 − P(L⁺ ≤ J) ≤ 1 − β, with X_i = (B − ℓ_i)/(B − J).
  - P(L⁺ ≤ J) + K(X) = 1 because a sum with a nonzero coefficient has no atom at 0 (`stat_add_stat_neg`), and `stat` is blind to a positive scale (`stat_map_mul`).
- *The calibration* (`gaffke_miss_le`, `gaffke_calibration`).
  - P(J ≤ b*_β) ≥ β, for independent units with their own finite loss laws at most B, each with mean at least J ≤ B.
  - J = B is covered: every loss with chance is B, and b*_β = B (`le_gaffkeBound_of_eq`).
  - J ≤ B matters only with no units.
  - `gaffke_one` is the case n = 1.
- *What HighProbability reports* (`highProbability_risk`, with `hpdRule` and `gaffkeBound_mono`).
  - The losses lie on a grid, nonincreasing in the threshold.
  - The rule meets α at a threshold whose risk exceeds α with chance at most 1 − β.
  - A rule that meets no threshold returns the last threshold, unmet (`RiskThreshold.Met` false). The guarantee is about the met case, as the C# reads it.

**Against the source.**
- *Their scope.* Their Theorem 1 is for any nonnegative laws. The Lean's is for finitely many values per unit, which is what a calibration set's losses can show. Their Lemma 7 for general laws is not needed.
- *Their Lemma 5's monotone θ.* They derive it from Y's log-concave density, citing that convolution keeps log-concavity (Karlin). The Lean needs neither: it carries the pair's likelihood-ratio order, not Y's density, through each added exponential.
- *Their Remark 4,* one random ordering serving every α at once, is not formalized.

**For the C#:** `RiskControl.HighProbability`'s guarantee is now proved as it reads: for independent units, it meets α at a threshold whose risk exceeds α with chance at most 1 − β. The docstrings that cite Vlassis and Thomas's Theorem 1 can name `highProbability_risk` and `gaffke_calibration`. No code changes.

**Checks.**
- The race identity against the Dirichlet law's divided differences, and the same-sign identities, exactly, on 600 cases.
- The flux identities, the likelihood-ratio inequality and η's sign on 1,500 edges. η is negative on 647 unsorted edges, so the sorting is used.
- The local insertion's identity on 6,079 insertions.
- The chain induction on 140 two-point systems at 2,201 levels: the hybrid means never fall, and the final chain stays within α.
- The two-point mixture on 373 finite laws.
- b*_β's calibration at n ≤ 4 units with their own laws, on 336 cases. The least slack is 0.100.
- The high-probability rule's event as the Lean states it.
  - On 480 random cases, all within 1 − β.
  - On 72 targeted cases where the rule meets α at a risk above it, within 1 − β; the worst failure chance is 0.343.

## RoBAS-Full's intervals at n ≥ 4

Asked on 2026-11-14 (the C# chip's commission, at the user's word): prove that RoBAS-Full's conformal sets are intervals for n ≥ 4, or find a counterexample there. The user chose scope A: everything in Lean except four bounds on Dawson's function, which are hypotheses stated on the module's own function and certified in the checks. Sources: RoBAS (arXiv 2607.04236v2), in TeX:
- the RoBAS-Full derivation (appendix, "Additional Details on BWM");
- Bersson and Hoff's Lemmas 1–2, as their `lem:prediction_region_interval`;
- RoBAS-EB's interval theorem (main text §3.4) for the route.

Proved on 2026-11-14 in `ConformalValidity.lean`, section `RobasFull` (59 declarations), PROOF-CLOSED, with the standard three axioms only. The checks (section "RoBAS-Full's intervals at n >= 4") ran first.

**The statement** (`robasFull_interval`).
- *Hypotheses.* n ≥ 4, calibration residuals not all equal, and `RobasDawsonBounds`.
- *Conclusion.* For every α, the set {z : p(z) > α} is an interval (`Set.OrdConnected`).
- *The p-value.* `robasPValue` is (#{j : S_cand ≤ S_j} + 1)/(n + 1), the candidate's own tie counted. This is `ResidualShrinkage`'s `atOrAbove + 1 > level` with level ⌊α(n + 1)⌋.
- *The score.* `robasFullScore` is the paper's −n s²/(2σ̂²) + ln ₁F₁(1; 3/2; −n r̄²/(2σ̂²)), of the n others, at the shared plug-in σ̂² = T/n, where T = Σ_all (r − r̄)² (`robasSpread`).
- *Not covered.* Calibration residuals all equal: there T vanishes at z = r̄, and the C# ties every unit.

**Dawson's function in the module** (no hypotheses).
- *Definitions.*
  - `kummerMoment j x` = ∫₀¹ v^{2j} e^{xv²} dv.
  - `kummerHalf x` = e^{−x} K₀(x) = ₁F₁(1; 3/2; −x) = D(√x)/√x.
  - `robasG` = K₁/K₀.
- *Its series* (`kummerMoment_hasSum`). K_j(x) = Σ xᵏ/(k!(2k + 2j + 1)), the series `LogKummer` sums.
- *Its derivatives.*
  - K_j′ = K_{j+1}, by differentiating under the integral.
  - The ODE 2x K_{j+1} + (2j + 1) K_j = eˣ, by the fundamental theorem.
- *Its convexity.* K₁² ≤ K₀K₂ (Cauchy–Schwarz), so g is nondecreasing (`robasG_monotone`).
- *Ω and its derivatives.* With Ω(y) = −2 ln ₁F₁(1; 3/2; −y²/2):
  - Ω′ = ω = 2y(1 − g);
  - ω′ = 4g(1 − x(1 − g)) at x = y²/2, through the ODE;
  - κ = 2Ω − yω has κ′ = c = ω − yω′;
  - ω′ ≤ 4 (`robasOmegaDeriv_le`).

**The reduction** (exact).
- *The score identity* (`robas_score_diff`). 2T(S_j − S_cand) = Δ_j(z) = (r_j − z)L_j(z) − T(Ω(p) − Ω(q)), where:
  - L_j(z) = (n − 1)z + (n + 1)r_j − 2S;
  - p and q are the leave-one-out sums over √T.
- *Single crossing* (`robasDelta_rises`). At each zero z ≠ r_j, Δ_j′(z)(r_j − z) = T·Ξ. Here d = p − q, and the slope at the zero is (`robas_slope_algebra`):
  - Ξ = (n − 1)d² + (d ω(p) − (Ω(p) − Ω(q))) + t d (κ(p) − κ(q));
  - t is the candidate's standardized deviation, (d − Ω̄/(n + 1))/2 at a zero.
- *Admissibility.* The other n − 1 standardized deviations bound (t + u)² ≤ (n − 1)(1 − t² − u²), by Cauchy–Schwarz.
- *Each unit's set is an interval* (`robas_acceptance_ordConnected`).
  - Δ_j(r_j) = 0.
  - Left of r_j, Δ_j rises through its zeros; right of r_j, it falls (`nonneg_of_rising`, `nonneg_of_falling`).
  - So {Δ_j ≥ 0} is an interval. This needs no limits at ±∞.
- *A common point* (`robas_common_point`).
  - Every unit's score is F_z(r_i)/(2T) less a shared term, with F_z(ρ) = nρ² + (Σ − ρ)² − TΩ((Σ − ρ)/√T).
  - F_z″ = 2n + 2 − ω′ > 0.
  - At z₀ with n z₀ = S g(S²/(2T(z₀))), which exists by the IVT, the candidate sits at F_{z₀}'s minimum.
- *Bersson and Hoff's step.* Intervals through z₀ make every superlevel set of the count an interval.

**The core and the four bounds** (`robas_core`, `RobasDawsonBounds`).
- *Lower bound on Ξ.*
  - Ξ ≥ d²[(n − 1) + ω′_min/2 − |κ(p) − κ(q)|/2 − |Ω̄| c_max/(2(n + 1))].
  - This uses Taylor's bound from ω′ ≥ ω′_min, and the mean value bounds |Ω̄| ≤ ω_max and |κ(p) − κ(q)| ≤ c_max|d|.
- *The two lengths.* Admissibility gives d² ≤ 2 − Ω̄²/((n + 1)(n − 1)).
- *n ≥ 5.* The long window suffices, and the bracket is at least 1.078.
- *n = 4.* Split on (5s)² ≤ 15(2 − 1.351²):
  - with the long window and |Ω̄| ≤ 1.6193, the bracket is at least 0.062;
  - otherwise |d| ≤ 1.351 and the short window applies, again at least 0.062.
- *The constants,* defined once in the module and read by the checks:

  | Bound | Value |
  |---|---|
  | \|ω\| | ≤ 1.8841 |
  | ω′ | ≥ −0.5511 |
  | \|c\| | ≤ 3.2365 |
  | κ's increment over \|b − a\| ≤ 1.415 | ≤ 4.2766 |
  | κ's increment over \|b − a\| ≤ 1.351 | ≤ 4.1060 |

**Against the source.**
- *Their RoBAS-EB route carries over.* Bersson and Hoff's two conditions hold.
  - Each I_j is an interval, though only from n = 4 here.
  - The common point comes from the scores' common convex profile.
- *Below n = 4.* The slope at a zero can be negative: the exact minimum is −2.14 at n = 2 and −0.93 at n = 3. This matches the C#'s pinned counterexamples, where one unit's β_j has three zeros.

**For the C#.**
- *What the theorem gives.* `ResidualShrinkage.Intervals` can serve RoBAS-Full from n = 4, as it serves RoBAS-EB, when the calibration residuals are not all equal.
- *The all-equal case* is outside the theorem. There every unit ties.
- *A point always in the set.* Take z₀ with n z₀ = S·g(S²/(2T(z₀))). There every unit's score is at least the candidate's, so p(z₀) = 1. It is a safe start for the bracketed search, in place of RoBAS-EB's â·r̄.

**Checks** (section "RoBAS-Full's intervals at n >= 4").
- *The four bounds certified.*
  - The constants are parsed from the Lean file and certified in decimal interval arithmetic with directed rounding, on 5,325 cells of [0, 12].
  - On each cell, g lies between its ends' series values, by monotonicity.
  - Results: ω ≤ 1.884025, ω′ ≥ −0.551087, |c| ≤ 3.236498. The window increments are 4.2765 at 1.415 and 4.1060 at 1.351, in exact rationals.
- *Past y = 12.* Analytic bounds through E W and E W² for W = 1 − V², with V's density ∝ e^{xv²}: ω ≤ 0.343, ω′ ≥ −0.111, 0 ≤ c ≤ 1.389.
- *The score identity.* Against the C#'s score form, on 1,689 units, worst 4.8e-15.
- *The pinned n = 3 residuals.* They split unit 1's set.
- *Direct grids at n = 4 to 7.* One piece on 218 sets, over every k.

## E35's follow-up — CUSUM's and Shiryaev–Roberts's guarantees

Asked on 2026-11-15 by the user: "The CUSUM and Shiryaev–Roberts guarantees". It is the item the addendum lists for the conformal Lean thread, and it decides whether `MartingaleAlarm`'s CUSUM and Shiryaev–Roberts alarms report proved guarantees.

- *Source.* Vovk, *Statistical Science* 2021, Propositions 4.2 and 4.4 and Corollaries 4.3 and 4.6, as M-2's "Posed, not proved" records them. The paper's text was read for M-2 on 2026-11-08, but is not in the local library and was not re-read.
- *Status.* Proved on 2026-11-15 in `ConformalValidity.lean`, section `ChangeAlarms`, right after M-2's `Betting` (36 declarations). PROOF-CLOSED, with the standard three axioms only. The checks (section "E35 follow-up: CUSUM's and Shiryaev-Roberts's guarantees") ran first.

**The model**, as `MartingaleAlarm` runs it. Each step's bet `fₙ` is chosen from the earlier p-values and read at the newest.
- *Shiryaev–Roberts* (`srRun`): `Rₙ = fₙ (Rₙ₋₁ + 1)` from 0. An alarm fires at `c ≤ Rₙ`, then the statistic restarts at 0.
- *CUSUM* (`cusumRun`): Page's `γₙ = fₙ max(γₙ₋₁, 1)` from 0, with the same alarm and restart.
- *Up to the first alarm* (`srFirst`, `cusumFirst`): the same runs, frozen at the first alarm. They count `τ ∧ n`.
- *The premises* are `ville`'s: the p-values are independent with law `μ` (uniform under exchangeability, by M-1), and each bet is measurable, reads only the past, and integrates to at most 1.

**The statements.**
- *The false alarms at every horizon.*
  - `sr_alarms_le`, `cusum_alarms_le`, and `alarms_stream_le` on a stream: `c E[A_N] ≤ N`. The expected number of false alarms among the first `N` arrivals is at most `N/c`, for both procedures.
- *CUSUM alarms no more often*, on every path (`cusum_count_le`, through `cusum_sr_invariant`). Between two CUSUM alarms, Shiryaev–Roberts alarms at least once: either CUSUM has had fewer alarms, or as many with its statistic no larger.
- *The first alarm at a horizon* (`sr_firstAlarm_le`, `cusum_firstAlarm_le`): `c P(τ ≤ N) ≤ E[τ ∧ N]`. CUSUM's first alarm comes no sooner (`cusumFirst_sr_invariant`).
- *The mean time to a false alarm* (`sr_meanRunLength`, `cusum_meanRunLength`): on a stream, `E τ ≥ c`. Here `τ = ∑ₙ 1{no alarm by n}`, which is infinite if no alarm ever fires.

**The route.** It needs no optional stopping at unbounded times.
- *The alarm count.* Let `Q` be the statistic after a step (0 after an alarm) and `Σ` the sum of the statistic at the alarms.
  - Each step adds `fₙ (Q + 1)` to `Q + Σ` (`srStep_sum`), which integrates to at most `Q + 1`.
  - So `E[Qₙ + Σₙ] ≤ n` (`lintegral_srRun_le`), through Tonelli at the newest p-value (`lintegral_pi_succ`).
  - Each alarm adds at least `c` to `Σ` (`srRun_count`).
- *The first alarm.* Frozen there, a step adds at most the step it counts (`lintegral_srFirst_le`). At the alarm, the statistic is at least `c`.
- *The mean.* `le_tsum_of_alarm_bound` turns `c P(τ ≤ N) ≤ ∑_{n<N} P(τ > n)` at every horizon into `c ≤ ∑ₙ P(τ > n) = E τ`.
  - If the alarm fires almost surely, the left side tends to `c`.
  - Otherwise every term is at least `P(τ = ∞) > 0`, and the sum is infinite.

**Against the source.**
- *The premise is wider.* Bets need only integrate to at most 1, which makes the capital a supermartingale; Vovk states a martingale. `Rₙ − n` is then a supermartingale, and every inequality runs the right way.
- *The frequency is stated in expectation.* Vovk's Corollaries 4.3 and 4.6 state an almost sure long-run frequency, lim sup `Aₙ/n ≤ 1/c`. Proved here is `E A_N ≤ N/c` at every horizon `N`, a finite-sample bound. The almost sure form, which needs a strong law, is not formalized.
- *CUSUM's comparison holds with restarts.* "CUSUM alarms no more often" is proved for the restarted runs `MartingaleAlarm` uses, path by path. Each procedure restarts at its own alarms.

**For the C#.**
- *Both guarantees are proved.* `AlarmGuarantee` can read `Proved` for CUSUM and Shiryaev–Roberts.
- *The wording.* It should be the proved form: "The mean time to a false alarm is at least c, and the expected number of false alarms among the first N arrivals is at most N/c." This replaces "repeated, the long-run frequency of false alarms at most 1/c", which is the almost sure form that is not proved.

**Checks.**
- *Exact enumeration on M-2's 4-point grid.* N = 6, so 4,096 paths, under 3 strategies and 4 thresholds.
  - `E[Q_N + Σ_N] = N` exactly.
  - `c A_N ≤ Σ_N` and `A^CUSUM ≤ A^SR` on every path, so `E[A_N] ≤ N/c`.
  - Frozen at the first alarm, `E[R] = E[τ ∧ N]` exactly, so `E[τ ∧ N] ≥ c P(τ ≤ N)`. CUSUM's first alarm never comes sooner.
- *A broken premise.* A bet averaging 2 breaks the count.
- *A seeded run.* Simple Jumper over 60,000 uniform p's at c = 20 gives 0.0471 Shiryaev–Roberts alarms per arrival (bound 0.05) and 0.0001 CUSUM alarms.

## Posed: M-1 under covariate shift (2026-11-15)

Posed by the conformal C# chip at the user's choice ("Build, checks now", 2026-11-15). The C# is d4c97d4 (`CovariateShiftPValues`) and its report e885e6a, in the addendum. Not yet proved.

- *The setting.* The joint density is Π_{t>m} w(z_t) · g(z₁, …, z_N) with g symmetric: Tibshirani, Barber, Candès and Ramdas 2019, Definition 1.
  - Positions 1..m (the reference sample) have weight 1, and the stream positions m+1..N have weight w.
  - The main case: independent Z_t ~ P for t ≤ m and Q for t > m, with dQ/dP = w.
  - The score α is a fixed function.
- *The p-values.* For n > m, let k = n − m and Pⱼ⁽ⁿ⁾ = πⱼ/k, where πⱼ = wⱼ e_{k−1}(w₋ⱼ)/e_k(w). Equivalently, Pⱼ is the share of arrangements putting unit j at n, each arrangement weighed by Π_{t=m+1..n} w(z_σ(t)).
  - pₙ = Σⱼ Pⱼ[αⱼ < αₙ] + θₙ Σⱼ Pⱼ[αⱼ = αₙ], with θₙ i.i.d. uniform on (0, 1] and independent of the data.
- *The claim.* p_{m+1}, …, p_N are i.i.d. uniform.
- *The route.* This is `seqRank_uniform`'s reverse-filtration step, with each ordering weighed by Π w over the stream positions.
  - Given the bag at n and every later unit, the arrangement has chance ∝ Π_{m<t≤n} w(z_σ(t)) (their Lemma 3), so pₙ is uniform given them.
  - Later p-values read only that bag and the later units.
  - Finite core: over every ordering of a fixed bag, weighed and normalized, the stream's tie-randomized weighted ranks are jointly uniform.
  - Equal weights give back M-1.
- *Checks so far (the C#'s).* Enumeration over six units: P(every pₜ ≤ uₜ) = Π uₜ to within 2e-15 over a grid. Ties are included, and so is a zero weight at a reference position. A stream unit with w = 0 has probability zero under the null, and the C# refuses it.
