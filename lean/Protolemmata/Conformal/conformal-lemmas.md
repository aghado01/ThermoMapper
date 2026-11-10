# Conformal validity — prediction sets valid whatever the model does

Pins the contract of the engine's conformal layer, one layer over retained draws and pointwise densities for every model with a pointwise density ([E16](../../../issues/regression/bass/extensions-brief.md#conformal-bayes-the-consumer)). The layer owns coverage and the score owns width. These are the finite statements of that split. They are the sibling of [GuardedEstimator](../BASS/guarded-estimator-lemmas.md): there, means valid whatever the emulator does; here, prediction sets valid whatever the model does. Written from the [chip brief](conformal-chip-brief.md).

**Status (2026-11-06).** Proved in `lean/Enthymemata/ConformalValidity.lean`, namespace `Conformal`, with no apologies (PROOF-CLOSED, 187 declarations counting E35's M-1 and M-2 and slice 9 at the end of this page; `#print axioms` shows `propext`, `Classical.choice`, `Quot.sound` only):
- CB-1: `card_rankLe_le`, `coverage_orderings`, `coverage_orderings_eq`, `coverage_strict_le`, `coverage_orderings_le`, `rank_uniform`, `coverage_tieBreak`, `coverage_given`, `coverage_law`, `coverage_stratum`;
- CB-2: `split_coverage`, `exchangeable_given_training`, `cqr_ordConnected`;
- CB-4: `sum_jointU`, `post_addOneIn`, `integratedScore_equivariant`, `latentScore_equivariant`, `latentScore_eq`, `sum_jointU_addOneIn`, `integratedScore_le_latentScore`, `accepts_latent_of_mixed`, `mixed_counterexample`;
- CB-3, with an error bound per unit: `accepts_or_fragile`, `accepts_shift_of_accepts_close`, `acceptsInflated_of_accepts`, `fragileEstimated_of_fragile`, `fragileEstimated_mono`, `coverage_close`, `coverage_inflated`, `card_fragile_le_estimated`.

[`conformal-checks.py`](conformal-checks.py) holds the exact-arithmetic checks, run before the Lean. Every assertion holds. Changed from the brief: CB-1's upper bound is an equality at the grid levels, and Hoff's strict rule fails under ties; CB-3's band is one-sided and tight, its bounds are per unit, and it gives a third option, the inflated rule. Nothing proved false. CB-3's prior art was read in full on 2026-11-06; the inflated rule is known, and the margin was not found. Open: the law forms on continuous spaces, CB-3's margin under an ε that holds only with high probability (for the inflated rule, Bhattacharyya et al. give 1 − α − 3√ν), and the Hoff finding under CB-4.

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

**Prior art, read.**
- Vovk, Nouretdinov and Gammerman 2003 (power martingales, the e-value reading).
- Fedorova, Gammerman, Nouretdinov and Vovk 2012 (plug-in).
- Volkhonskiy, Burnaev, Nouretdinov, Gammerman and Vovk 2017 (inductive conformal martingales against classical CUSUM).
- Vovk et al. 2021 (Simple Jumper).
- Vovk, *Statistical Science* 2021 (Ville's procedure, CUSUM, Shiryaev–Roberts).
- The text extractions were checked against the formulas used here: (6) of 2003, Algorithm 1 of 2021, and (12) to (15) of *Statistical Science*.

---

## Slice 9 — cross-fitted forms and asymmetric CQR

Posed and proved on 2026-11-09 for the [conformal layer's slice 9](../../../issues/regression/conformal-layer-brief.md#e16s-completion-the-plan-slices-7-onward), before any C#, from the [chip brief's assignment](conformal-chip-brief.md#e16-slice-9--cross-fitted-intervals-and-asymmetric-cqr-2026-11-09). The point is the contract: the exact rank rules the layer implements, and the hypotheses it must meet.

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
  - *Settled at K = 3, m = 2, ℓ ≤ 2.* At ℓ = 2 a refusal needs the test unit on top in two of its three triples, against once on average, so P ≤ 1/2, CV-b's bound by Markov. At ℓ = 1 the exact MILP worst case is 1/7, against 1/4.
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
