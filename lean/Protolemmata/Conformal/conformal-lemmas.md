# Conformal validity — prediction sets valid whatever the model does

Pins the contract of the engine's conformal layer, one layer over retained draws and pointwise densities for every model with a pointwise density ([E16](../../../issues/regression/bass/extensions-brief.md#conformal-bayes-the-consumer)). The layer owns coverage and the score owns width. These are the finite statements of that split. They are the sibling of [GuardedEstimator](../BASS/guarded-estimator-lemmas.md): there, means valid whatever the emulator does; here, prediction sets valid whatever the model does. Written from the [chip brief](conformal-chip-brief.md).

**Status (2026-11-06).** Proved in `lean/Enthymemata/ConformalValidity.lean`, namespace `Conformal`, with no apologies (PROOF-CLOSED, 119 declarations counting E35's M-1 and M-2 at the end of this page; `#print axioms` shows `propext`, `Classical.choice`, `Quot.sound` only):
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
