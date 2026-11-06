# Conformal validity — prediction sets valid whatever the model does

Pins the contract of the engine's conformal layer, one layer over retained draws and pointwise densities for every model with a pointwise density ([E16](../../../issues/regression/bass/extensions-brief.md#conformal-bayes-the-consumer)). The layer owns coverage and the score owns width. These are the finite statements of that split. They are the sibling of [GuardedEstimator](../BASS/guarded-estimator-lemmas.md): there, means valid whatever the emulator does; here, prediction sets valid whatever the model does. Written from the [chip brief](conformal-chip-brief.md).

**Status (2026-11-06).** Proved in `lean/Enthymemata/ConformalValidity.lean`, namespace `Conformal`, with no apologies (88 declarations, PROOF-CLOSED; `#print axioms` shows `propext`, `Classical.choice`, `Quot.sound` only):
- CB-1: `card_rankLe_le`, `coverage_orderings`, `coverage_orderings_eq`, `coverage_strict_le`, `coverage_orderings_le`, `rank_uniform`, `coverage_tieBreak`, `coverage_given`, `coverage_law`, `coverage_stratum`;
- CB-2: `split_coverage`, `exchangeable_given_training`, `cqr_ordConnected`;
- CB-4: `sum_jointU`, `post_addOneIn`, `integratedScore_equivariant`, `latentScore_equivariant`, `latentScore_eq`, `sum_jointU_addOneIn`, `integratedScore_le_latentScore`, `accepts_latent_of_mixed`, `mixed_counterexample`;
- CB-3: `accepts_or_fragile`, `accepts_shift_of_accepts_close`, `acceptsInflated_of_accepts`, `fragileEstimated_of_fragile`, `coverage_close`, `coverage_inflated`, `card_fragile_le_estimated`.

[`conformal-checks.py`](conformal-checks.py) holds the exact-arithmetic checks, run before the Lean. Every assertion holds. Changed from the brief: CB-1's upper bound is an equality at the grid levels, and Hoff's strict rule fails under ties; CB-3's band is one-sided and tight, and it gives a third option, the inflated rule. Nothing proved false. Open: the law forms on continuous spaces, CB-3 under an ε that holds only with high probability, and the Hoff finding under CB-4.

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

**Claim, the deterministic part.** Let ŝ be within ε of an exact score s on C. Then:
- ŝ refuses a unit that s accepts only when the unit is **fragile** (`accepts_or_fragile`). Fragile means accepted by s but with no more than k units of C at least 2ε below it, itself counted (`Fragile`).
- ŝ accepts a unit that s refuses only when the unit sits within 2ε of acceptance (`accepts_shift_of_accepts_close`).
- So a decision changes only inside a band around the threshold. With distinct scores, the loss side of that band is M = {u : s₍ₖ₊₁₎ ≤ s_u < s₍ₖ₎ + 2ε}: the accepted units within 2ε of the highest refused score.

**Claim, over orderings.** Let S be equivariant, and let the estimate ŝ_σ for each ordering σ be anything within ε of S (Monte Carlo error, an asymmetric computation). Then:
- coverage is at least 1 − (k + |M|)/m, M the exact score's fragile units (`coverage_close`). This is the brief's identity P(test ∈ band) = E|band|/(n + 1), with the one-sided band;
- the **inflated rule**, which accepts when more than k estimated scores sit at most 2ε above the candidate's, keeps 1 − k/m (`coverage_inflated`, from `acceptsInflated_of_accepts`).

**The engine's reading.** |M| is at most the count of the band read from any ordering's estimated scores alone (`card_fragile_le_estimated`, from `fragileEstimated_of_fragile`): k < #{ŝⱼ ≤ ŝ_u + 2ε} and at most k units 4ε below. With distinct estimates that is ŝ₍ₖ₊₁₎ − 2ε ≤ ŝ_u < ŝ₍ₖ₎ + 4ε, the brief's "within 4ε" on the high side. So the layer can report 1 − (k + m̂)/m from what it holds.

**Changed from the brief.** The brief's band is two-sided, |s_test − s₍ₖ₎| ≤ 2ε. The coverage loss needs only the side above the threshold, and the bound is tight: an adversary within ε, the test unit down and the rest up, refuses exactly the units the bound counts. The k-th order statistic moving at most ε is checked, not formalized; the count form replaces it.

**Checks.**
- The deterministic part holds on 57,240 cases, exhaustive on a grid and over a seeded random sweep. A loss lies only in M, the inflated rule keeps every exact acceptance, M lies inside the estimated band, and order statistics move at most ε.
- The adversary attains the bound for every unit.
- Over orderings, the bound holds for an adversarial and a noisy estimate, and is attained by the adversary with distinct scores.
- Conformal Bayes by add-one-in weights over T = 20 and 200 draws on a Beta–Bernoulli grid: ε is 0.007 to 0.073, and every bound holds. Here the draws' decisions match the exact ones.

**Decision: whether the full-conformal layer reports its margin, and when it falls back to the split form.** One lemma gives three policies.
- *Report the margin:* coverage ≥ 1 − (k + m̂)/m, with m̂ the estimated band's count, taken as the maximum over the candidates on the grid.
- *Inflate by 2ε:* keep 1 − k/m at the cost of width.
- *Fall back to the split form (CB-2):* when m̂/m is too large for the requested level, or when ε cannot be trusted.

ε is a sure bound in these statements. A Monte Carlo ε holds only with some probability, and its failure probability adds to the miscoverage. ε comes from the sampler's own error estimate for each unit's add-one-in score, so it is only as good as the mixing. A chain that fails its diagnostics sends the layer to the split form.

**Prior art** (abstracts only, found by web search on 2026-11-06; none of these is on disk).
- Barber, Candès, Ramdas and Tibshirani, *Conformal prediction beyond exchangeability* (AoS 2023, arXiv 2202.13415): coverage gaps bounded by weighted total-variation terms, and a randomization technique for algorithms that are not symmetric.
- Ndiaye, *Stable conformal prediction sets* (ICML 2022, arXiv 2112.10224): a single fit, with a stability bound added to keep coverage. This is the inflated rule with a stability bound in place of ε.
- Razafindrakoto, Celisse and Lacaille, *Approximate full conformal prediction in an RKHS* (arXiv 2601.13102): a slack τᵢ(y) added to the training scores, so the approximate region contains the full one. The same route, per point.
- Bhattacharyya, Zhang and Barber, *Approximating full conformal prediction: … the tournament correction* (arXiv 2605.29200): 1 − 2α guaranteed, about 1 − α under stability.
- Abad Martinez et al. (AAAI 2023, arXiv 2202.01315): influence-function approximations whose error vanishes as n grows.
- None of these abstracts states the one-sided fragile count or its attainment. Whether a paper does is open.
