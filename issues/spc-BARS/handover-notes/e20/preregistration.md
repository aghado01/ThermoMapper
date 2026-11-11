# E20, the graph panel: what it reads and what I expect, written before any sweep ran

The question (stopping primer §13.17): which of round 2's readings survive across graphs, and which belong to one graph?

## The panel

- **Graphs, q = 20.** BlattHierarchy (1,200 × 2, 50 per fine cluster, leaf spread 0.15, seed 42), AnisotropicGaussian (1,200, 400 per class, separation 4, ratio 1, seed 55), the toy (3,200 × 2, seed 42), Landsat (6,435 × 36, raw). K = 5, 10, 20. Gaussian and Cauchy. 24 graphs.
- **q arm.** Gaussian at K = 10 on the same four datasets, q = 5 and 100. 8 more sweeps, added to E20 as written because §13.17's table rests the melt's q-dependence on one graph.
- **Control.** The uniform cloud (400 × 2, seed 20261101), k = 10, Gaussian, q = 20, must reproduce melt/T_ps = 0.87. The toy at k = 10 must reproduce S₂ = 0.858 at T = 0.00446 and S₂·N about 1.12–1.14 at T = 0.357. If either misses, the harness is wrong and the panel stops.
- **Fixed on every graph.** Euclidean on raw features, mutual KNN, MST repair, `PathNeighbor` (no effect), bandwidth = mean edge length, gauge 1/K̂ over edges after repair, no LMP, floor 10⁻¹⁰. Each graph's settings are read from `Resolved` into `graphs.tsv`.
- **Grid.** 257 points uniform in ln T over [10⁻⁶, 3] (h = 0.058 e-folds), both ends anchored, as the user decided, plus the bracket's two ends as flagged extras. Volleys are nested subsets: 9, 17, 33, 65, 129, 257.
- **Sampling.** 8 replicas × (500 burn-in, 1,000 cycles), each temperature and replica independent, random start. Base seed 20261102. One pass, one seed: 8 replicas give R̂ and a spread, not a second draw of the panel. That is half the draws of round 1's main budget, 8 × (1,000, 2,000), chosen for cost on a shared machine, so absolute quiet shares hold at this budget only.
- **Machine.** Below-normal priority, at most 8 workers. Rows append per temperature, so a rerun skips what's done.

## Fixed before running: the tests

- **Margin.** A floor of n = 5% of N: 2n²/N² = 0.005 in S₂, 2n·ln 2/N = 0.0693 in the information currency, n/N = 0.05 in F_n.
- **Quiet interval.** Δ + q_α·SE(Δ) < floor, α = 0.05 Bonferroni over the volley's intervals. S₂ is scored two ways: within-run SE with z (only where R̂ ≤ 1.1 at both ends), and Welch on replica means with t at df from Welch–Satterthwaite. The information currency and F_n have replica means only, so Welch only.
- **Split ratio.** ρ = Δ(left child)/Δ(parent) on each parent not certified quiet. Resolved if |ρ − ½| ≤ ¼.
- **Stop.** The first volley v whose parents at v − 1 are all quiet or resolved. If none, "no stop by v5".
- **Closure.** Energy: |Δ⟨E⟩ − trapezoid of Var(E)·Δβ| over |Δ⟨E⟩| under 10% on intervals not quiet. Var(m) likewise through (m_t, E_{t+1}).
- **Equilibration flag.** R̂ > 1.1 in S₂ or E. Replica outlyingness: max |x_r − median|/(1.4826·MAD) over replicas.
- **Endpoint decided.** Minimum pairwise ARI over replicas at θ = 0.5 of at least 0.9.
- **Melt.** T½ = where S₂ falls through half its value at the bracket's cold end, hottest crossing. Also the C_v peak and the Var(m) peak. Each over T_ps(q), and over s̄/(4 ln(1 + √q)), which removes the strength.
- **Structure's hot edge.** The hottest T with S₂·N ≥ 2.

## What I expect

These come from §13.17's arithmetic. Each can fail.

1. **Hot end, Gaussian.** S₂·N at T = 0.357 stays in 1.08–1.20 on all 12 Gaussian graphs at q = 20, at every K, since s̄ is pinned. The hottest T with S₂·N ≥ 2 is below 0.357 on all of them.
2. **Hot end, Cauchy.** Close to Gaussian in S₂·N at T = 0.357 (under 1.30), since s̄ is 0.52–0.56, even though the kernel's tail is heavier. Cauchy doesn't declare BWD units, so nothing here promises it, and its largest failure would be Landsat.
3. **Melt across K.** On each dataset, T½/T_ps moves by under 25% across K = 5, 10, 20 within a kernel, wherever one melt dominates. Across datasets it does not settle on 0.87: the toy and Landsat, with wide strength spreads, differ from the cloud by more than 20%.
4. **Melt across q.** On each dataset at K = 10, T½/T_ps moves by under 15% across q = 5, 20, 100. This is the one part of §13.14 I expect to hold on every graph.
5. **Cold end.** S₂ at T = 0.00446 moves by K as the envelopes said on the toy, and sits inside its envelopes everywhere. C43 is a theorem: S₂ outside its envelope by more than noise means a bug or a stuck sampler, not a surprise about the data.
6. **Gaussian against Cauchy, cold.** On BlattHierarchy and AnisotropicGaussian the Gaussian graphs sit at S₂ well under 1 from the cold end at every K, because their bridges can't be felt; Cauchy reaches S₂ near 1.
7. **Stop volley.** It moves with the graph: on at least two datasets the S₂ stop volley differs by at least one level across K or kernel. Cauchy on BlattHierarchy needs at least as many volleys as Gaussian there, since it has more events to resolve.
8. **Currency.** The information currency stops no later than S₂ on the toy and Landsat, where shedding dominates S₂'s active intervals (E11).
9. **Equilibration.** R̂ flags cluster at melts, not on plateaus, at every K. At q = 5, there are fewer flags than at q = 100.
10. **Spread and width.** Across graphs, the S₂ transition's width (the 10%–90% fall in e-folds, on single-melt graphs) rises with ln(s_p90/s_p10). Spearman above 0.5. This is the weakest expectation; §13.17 says only that spread is a candidate.

**Note added after the two control cells, before the panel ran, with nothing above changed.** Both controls pass: the cloud gives T½/T_ps = 0.871; the toy at k = 10 gives S₂ = 0.859 at T = 0.00446 and S₂·N = 1.115 at T = 0.357. The toy shows T½ isn't a melt on a graph with levels: half its cold S₂ goes at the stripe split, T½ = 0.16·T_ps, while the C_v peak sits at 1.21·T_ps, near the melts round 1 measured. So expectations 3 and 4, written for T½ "wherever one melt dominates", get read on both T½ and the C_v peak, and I'll say which each statement rests on.

What would change the picture: if expectations 1, 3, 4 and 7 all hold, §13.17's reading stands. If 1 fails on Gaussian graphs, the pinned-strength account is wrong or incomplete. If 7 fails, the stopping quantities are less graph-dependent than §13.17 says, which would be good news for a general rule.
