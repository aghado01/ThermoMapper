### 13.20 Addendum, 2026-11-02: E20, the graph panel (C61, C62, C63)

E20 asks which of round 2's readings survive across graphs and which belong to one graph. The tests and ten expectations were
written before any sweep ran. Three follow-ups weren't planned; each came from chasing an anomaly, and each is marked.

- **Survive.** The hot end. T_ps as the scale of the energy melt and of the hottest structure, once s̄ replaces e^(−1/2). The
  q-form. C43 at q ≤ 20. The quiet certificate, never contradicted. Closure's h² law on melts.
- **Belong to one graph.** Melt = 0.87·T_ps, which is the cloud's. T½ as "the melt". The cold multiplier 0.05, which is the
  toy's at k = 10. The stop volley. How hard the sampler has to work.
- **New.** R̂ and closure catch different failures, and at q = 100 they split cleanly (C62). The stop rule as written reverses
  itself, and a version with a noise band mostly doesn't (C63). The lower percolation envelope gets T½ from the graph alone
  to a factor of 1.3–1.7 at a given K, where T_ps is off by up to 13 (C61).

**What ran.** Harness, preregistration and every readout: `side-session-notes/.../e20/`.
- *Sampling.* 8 × (500, 1,000), each temperature and replica independent, random start, base seed 20261102. One pass and one
  seed, at half the draws of round 1's main budget, so absolute quiet shares hold at this budget only.
- *Grid.* 257 points uniform in ln T over [10⁻⁶, 3], both ends anchored, plus the bracket's ends as extras. Volleys nest: v0–v5
  are 9, 17, 33, 65, 129, 257 points at h = 1.86, 0.93, 0.47, 0.23, 0.12, 0.058 e-folds. The bracket holds 4.38 e-folds of it, so
  v2–v5 put about 10, 19, 38 and 76 intervals inside it.
- *Build.* Commit 04cfbcb, clean `src`. 3 h 54 min at below-normal priority with at most 8 workers on a shared machine.
- *Controls, both passed before the panel started.* The cloud gives T½/T_ps = 0.871. The toy at k = 10 gives S₂(0.00446) = 0.859 and
  S₂·N(0.357) = 1.115.
- *Scope of "survives".* K, kernel and q varied. Metric, filter, repair, bandwidth rule, gauge, `ScaleEdges`, LMP and floor did
  not, and §13.18 shows those move readings too. Everything here is about this recipe.

**The graphs as resolved.** Euclidean on raw features, KNN, mutual filter, MST repair, `PathNeighbor` (no effect), Gaussian or
Cauchy with a = the mean edge length, g = 1/K̂ over the edges after repair, no LMP, floor 10⁻¹⁰. BlattHierarchy: 50 per fine
cluster, leaf spread 0.15, seed 42. AnisotropicGaussian: 400 per class, separation 4, ratio 1, seed 55. The toy: seed 42.
Per graph, with in-degree skew and the weakest J: `e20/graphs.tsv`.

| Graph | k | K̂ | Largest before repair | Repair edges | a | s̄ Gaussian (p10, p90) | s̄ Cauchy (p10, p90) |
|---|---|---|---|---|---|---|---|
| Uniform cloud, 400 × 2 | 10 | 8.47 | 1.000 | 0 | 0.0584 | 0.603 (0.386, 0.816) | Not run |
| BlattHierarchy, 1,200 × 2 | 5 | 3.62 | 0.042 | 82 | 0.0831 | 0.710 (0.219, 1.164) | 0.634 (0.197, 1.048) |
| | 10 | 7.17 | 0.042 | 37 | 0.0956 | 0.650 (0.206, 1.051) | 0.574 (0.193, 0.922) |
| | 20 | 14.71 | 0.126 | 20 | 0.1269 | 0.620 (0.218, 0.969) | 0.546 (0.201, 0.852) |
| AnisotropicGaussian, 1,200 × 2 | 5 | 3.91 | 0.328 | 19 | 0.1717 | 0.631 (0.200, 1.030), 2 components | 0.564 (0.205, 0.910) |
| | 10 | 7.97 | 0.336 | 5 | 0.2273 | 0.621 (0.214, 0.955) | 0.552 (0.220, 0.844) |
| | 20 | 15.73 | 0.998 | 2 | 0.3063 | 0.616 (0.209, 0.921) | 0.545 (0.215, 0.808) |
| The toy, 3,200 × 2 | 5 | 3.86 | 0.192 | 83 | 0.1204 | 0.643 (0.150, 1.022), 2 components | 0.573 (0.169, 0.905) |
| | 10 | 7.91 | 0.993 | 12 | 0.1622 | 0.636 (0.159, 0.945) | 0.563 (0.180, 0.834) |
| | 20 | 15.98 | 0.999 | 1 | 0.2269 | 0.633 (0.158, 0.900) | 0.559 (0.178, 0.789) |
| Landsat, 6,435 × 36 | 5 | 2.81 | 0.779 | 983 | 22.71 | 0.614 (0.185, 1.170) | 0.524 (0.161, 0.989) |
| | 10 | 5.32 | 0.946 | 282 | 23.84 | 0.616 (0.129, 1.192) | 0.527 (0.113, 1.009) |
| | 20 | 10.64 | 0.991 | 52 | 25.27 | 0.620 (0.133, 1.217) | 0.531 (0.128, 1.030) |

**The readings, q = 20.** T_lo = 0.00446 and T_hi = 0.357 are the bracket's ends, T_ps = 0.0892. T½ is where S₂ falls through half its
value at T_lo, the hottest such crossing. The hot edge is the hottest T with mean S₂·N ≥ 2: §13.19's T_melt on means rather than
medians, and not C45's T_melt. "Stop" gives the rule as written with within-run SE, then with Welch, then with C63's band, and a
dash means none by v5; † marks a band stop that a later volley reverses. "Flags" counts grid temperatures with R̂ > 1.1 in S₂ or E. ARI is
the pooled partition at θ = 0.5 against truth, best over T, per level from coarse to fine.

| Graph | S₂·N(T_hi) | S₂(T_lo) | S₂(10⁻⁶) | T½/T_ps | C_v peak/T_ps | Hot edge/T_ps | Stop | Flags | Best ARI |
|---|---|---|---|---|---|---|---|---|---|
| Cloud k10 Gaussian | 1.104 | 1.000 | 1.000 | 0.871 | 1.15 | 1.47 | 5, 5, 4 | 2 | |
| Blatt k5 Gaussian | 1.161 | 0.041 | 0.075 | 0.689 | 1.45 | 2.07 | 2, 2, 2 | 0 | 0.145, 0.562, 0.997 |
| Blatt k5 Cauchy | 1.135 | 0.070 | 1.000 | 0.275 | 1.21 | 1.87 | –, –, 3 | 0 | 1.000, 1.000, 0.997 |
| Blatt k10 Gaussian | 1.120 | 0.041 | 0.100 | 1.00 | 1.45 | 1.85 | 2, 2, 2† | 2 | 0.193, 0.698, 0.997 |
| Blatt k10 Cauchy | 1.102 | 0.053 | 0.998 | 0.693 | 1.36 | 1.63 | –, –, 4 | 2 | 1.000, 1.000, 0.997 |
| Blatt k20 Gaussian | 1.105 | 0.042 | 0.167 | 1.10 | 1.45 | 1.67 | –, 4, 4 | 7 | 0.332, 1.000, 0.998 |
| Blatt k20 Cauchy | 1.090 | 0.051 | 0.980 | 0.836 | 1.21 | 1.47 | –, 3, 3 | 7 | 1.000, 1.000, 0.995 |
| Aniso k5 Gaussian | 1.133 | 0.304 | 0.332 | 0.409 | 1.29 | 1.83 | 3, 3, 3 | 0 | 0.998 |
| Aniso k5 Cauchy | 1.112 | 0.330 | 1.000 | 0.363 | 1.08 | 1.63 | –, –, 3 | 0 | 0.998 |
| Aniso k10 Gaussian | 1.111 | 0.323 | 0.335 | 0.945 | 1.36 | 1.67 | –, –, 5 | 4 | 0.998 |
| Aniso k10 Cauchy | 1.095 | 0.331 | 1.000 | 0.815 | 1.21 | 1.49 | –, –, 3 | 0 | 0.998 |
| Aniso k20 Gaussian | 1.103 | 0.328 | 0.999 | 1.01 | 1.36 | 1.58 | –, –, 2† | 8 | 0.998 |
| Aniso k20 Cauchy | 1.089 | 0.335 | 1.000 | 0.874 | 1.21 | 1.41 | –, –, – | 6 | 0.998 |
| Toy k5 Gaussian | 1.137 | 0.680 | 0.991 | 0.082 | 1.36 | 1.85 | 3, –, 3 | 0 | 0.696 |
| Toy k5 Cauchy | 1.115 | 0.948 | 1.000 | 0.129 | 1.15 | 1.64 | 3, 3, 2 | 0 | 0.694 |
| Toy k10 Gaussian | 1.115 | 0.859 | 0.999 | 0.160 | 1.21 | 1.69 | –, –, 5 | 4 | 0.708 |
| Toy k10 Cauchy | 1.098 | 0.984 | 1.000 | 0.218 | 1.15 | 1.50 | –, –, 4 | 3 | 0.707 |
| Toy k20 Gaussian | 1.106 | 0.982 | 1.000 | 0.292 | 1.29 | 1.53 | –, –, 5 | 5 | 0.710 |
| Toy k20 Cauchy | 1.091 | 0.994 | 1.000 | 0.322 | 1.15 | 1.36 | –, –, 5 | 6 | 0.710 |
| Landsat k5 Gaussian | 1.137 | 0.993 | 1.000 | 0.549 | 1.53 | 1.91 | 4, 4, 4 | 1 | 0.235 |
| Landsat k5 Cauchy | 1.107 | 1.000 | 1.000 | 0.464 | 1.02 | 1.62 | 3, 3, 3 | 1 | 0.208 |
| Landsat k10 Gaussian | 1.118 | 0.988 | 1.000 | 0.609 | 1.72 | 1.99 | –, 5, 5 | 5 | 0.236 |
| Landsat k10 Cauchy | 1.095 | 0.999 | 1.000 | 0.558 | 1.45 | 1.72 | –, 4, 4 | 8 | 0.230 |
| Landsat k20 Gaussian | 1.109 | 0.988 | 1.000 | 0.654 | 1.72 | 2.04 | –, –, 5 | 9 | 0.278 |
| Landsat k20 Cauchy | 1.090 | 0.997 | 1.000 | 0.561 | 1.45 | 1.75 | –, 4, 4 | 11 | 0.279 |

C_v's peak is quantized by the grid to about 6%.

**The ten expectations, scored.** Written before the panel, unchanged since, apart from a note before it started that 3 and 4 would be read on the C_v peak too.

| # | Expected | Measured | Verdict |
|---|---|---|---|
| 1 | Gaussian S₂·N(0.357) in 1.08–1.20, hot edge below 0.357 | 1.103–1.161, hot edge at most 0.185 | Held |
| 2 | Cauchy under 1.30 | 1.089–1.135. The largest is BlattHierarchy at k = 5, not Landsat | Held |
| 3 | T½/T_ps moves under 25% across K where one melt dominates, and doesn't settle on 0.87 | T½ max/min: Landsat 1.19, 1.21, the rest 1.59–3.57. AnisotropicGaussian has one melt and moves 2.4–2.5. C_v peak: at most 1.12 on 7 of 8, Landsat Cauchy 1.42. T½ runs 0.08–1.10 | Failed on T½, held on the C_v peak on 7 of 8. Second half held |
| 4 | T½/T_ps moves under 15% across q at k = 10 | 1.05–1.11 as run, but the q = 100 sweeps weren't stationary. At burn-in 5,000: 1.01–1.08 on three datasets, the toy's q = 100 T½ not re-measured. C_v peak: 1.04–1.10 on three, Landsat 1.24 | Held on T½ after the fix. C_v peak on 3 of 4 |
| 5 | The toy's cold end moves with K, S₂ inside its envelopes | Gaussian 0.68, 0.86, 0.98. Cauchy 0.95, 0.98, 0.99. Envelopes: see C62's section | Held at q ≤ 20. The q = 100 breach is the stuck sampler the expectation names |
| 6 | Gaussian BlattHierarchy and AnisotropicGaussian well under 1 at the cold end, Cauchy near 1 | Cauchy at least 0.98. Gaussian on 5 of 6: AnisotropicGaussian at k = 20 reaches 0.999, since mutual edges join its classes (largest before repair 0.998) | Held, one exception |
| 7 | The stop moves with the graph on at least two datasets, and Cauchy BlattHierarchy needs at least as many volleys | Moves on all four. At k = 5 and 10 Cauchy has no stop against Gaussian's v2. At k = 20 under Welch it's v3 against v4 | Held, one exception |
| 8 | The information currency stops no later than S₂ on the toy and Landsat | 11 of 12, Landsat k = 5 Cauchy the exception | Held, but see the stop rule below |
| 9 | R̂ flags at melts at every K, fewer at q = 5 than q = 100 | q = 20: 89 flags, 74 inside an event, sitting around the C_v peak. q = 5: 0. q = 100: 44–54, of which 64–80% lie below 0.3·T_ps | Held at q ≤ 20. At q = 100 they mark drift, not melts |
| 10 | The transition's width rises with ln(s_p90/s_p10), Spearman above 0.5 on single-melt graphs | 0.46 on AnisotropicGaussian and Landsat (n = 12), 0.17 on all 24 | Failed |

The preregistration's test: if 1, 3, 4 and 7 hold, §13.17's reading stands. 1, 4 and 7 held, and 3 split in a way that corrects it: the
T_ps scale holds for energy and the hot edge, and T½ doesn't follow it.

**The scale survives, on C_v and the hot edge (C61).**
- *Hot end.* 4·T_ps is safe on all 24 graphs, with the hot edge at most 2.07·T_ps. Across q it's 2.37–2.64 at q = 5, and 1.42–1.65 at q = 100 on
  the unequilibrated sweeps. S₂·N exceeds §13.17's first-order 1 + s̄/(qT) by 0.013–0.062 at q = 20 and 0.115–0.150 at q = 5. It stays one
  fact about the recipe, as C58 said, and is now measured at more K and on a second kernel.
- *What s̄ buys.* Over T_ps, the C_v peak is 1.02–1.72, with Cauchy running low. Over s̄/(4 ln(1 + √q)), it's 1.16–1.70, the same range for both
  kernels. The hot edge goes from 1.36–2.07 to 1.47–2.00 the same way. So the kernel's effect on the temperature scale is its s̄, which the graph gives
  with no SPC. This fits §13.19's candidate B. It isn't evidence off the recipe.
- *Three "melts", three things.* C_v's peak counts bonds, so it reads the loss of local order. T½ counts connected mass, so it reads necks and levels. They agree
  only where there's one homogeneous melt. On the toy T½ is the stripe split, and it moves 3.6 times across K while C_v's peak moves 1.12. Var(m)'s peak sits
  at 0.05–0.34·T_ps on BlattHierarchy, AnisotropicGaussian and the toy, and it's a melt only on the cloud and Landsat.
- *The order.* On BlattHierarchy and AnisotropicGaussian the finest level's ARI stops reading 0.9 at 0.09–0.57·T_ps, far below the C_v peak at 1.08–1.45, and the
  hot edge lies 1.15–1.59 times above the C_v peak. Clusters readable at θ = 0.5 are gone well before the energy melt finishes, so the C_v peak doesn't locate the
  dense pass. What does depends on open question 1.
- *K ≥ 10.* Most of T½'s failure across K is k = 5, K̂ near 4. From k = 10 to 20, T½ moves at most 1.21 except on the toy, where it isn't a melt.

**The cold multiplier is the toy's, and it can cost levels (C61).** S₂(T_lo) runs from 0.041 to 1.000 across the panel. On Landsat almost nothing lies below T_lo. On
Cauchy BlattHierarchy nearly all of S₂'s range does, along with two of its three levels:

| BlattHierarchy, Cauchy | L0, coarse: T with ARI ≥ 0.9 | L1 | L2, fine | T_lo |
|---|---|---|---|---|
| k = 5 | 4.6×10⁻⁶ – 2.6×10⁻⁵ | 3.5×10⁻⁵ – 0.0025 | 0.0056 – 0.028 | 0.00446 |
| k = 10 | 3.0×10⁻⁶ – 1.7×10⁻⁵ | 2.3×10⁻⁵ – 0.0017 | 0.0035 – 0.038 | |
| k = 20 | 2.7×10⁻⁶ – 1.5×10⁻⁵ | 2.0×10⁻⁵ – 0.0014 | 0.0028 – 0.045 | |

Every L0 and L1 temperature sits below T_lo, and each reads 1.000 at its best. The Gaussian graphs read L2 alone, plus L1 at k = 20, because their repair edges sit at the kernel's clamp, J about 10⁻³⁰¹, and can't be felt. This
bears on the user's grid decision. The uniform grid isn't in question: it resolved everything between its anchors. But it can only show what the anchors contain, and an anchor taken from T_ps would never have shown these levels. Where the
anchors come from is still open, and the next reading is the cheap way to set them.

**What the graph says first (C61).** No SPC, and it's C43 applied.
- *T½ from the lower envelope.* Find where Bernoulli(p̃) falls through half its own value at T_lo. T½ over that is 0.98–1.29 at k = 5, 1.48–2.51 at k = 10, 2.12–3.20 at k = 20. That is a factor of 1.3–1.7 at
  a given K and 3.3 overall, against 13 for T_ps. It rises with K for a reason: on a tree the FK bonds are exactly Bernoulli(p̃), so a near-tree graph sits on the lower envelope and loops lift it off.
- *The sandwich at T½.* 0.84–2.26 e-folds wide at q = 20, 0.91–1.17 at q = 5, 1.87–2.36 at q = 100, on a rigorous bound.
- *The cold end.* The upper envelope at T_lo bounds S₂ there: 0.12–0.15 on Cauchy BlattHierarchy, whose S₂ reaches 1. So the graph says, before a sweep, that most of S₂ lies below T_lo. The weakest edge says how far: T = J_min/ln(q − 1) is within a factor of 1.3–2.4 of
  where S₂ finishes rising on the Cauchy BlattHierarchy and AnisotropicGaussian graphs, where that edge carries a real split. Where it holds only a pendant point it predicts nothing: the toy at k = 20, Gaussian, gives 10⁻¹⁰ against 0.003. The envelope weighs that automatically; the edge alone doesn't.
- *Sampler-free isn't data-free.* E14's first question, whether every curve sits between the envelopes, is answered below. Its second, the width at each event, is answered for T½ only.

**R̂ and closure catch different failures (C62, unplanned).** Energy closure at q = 100 failed on every volley at every graph, by 4–15 times in the median over the melt window, and the factor didn't shrink with h. Discretization shrinks, so
this wasn't resolution. A trace of energy by quarter of the recorded run:

| Graph, q, T | Burn-in | Mean E by quarter | R̂_E | Trapezoid/ΔE |
|---|---|---|---|---|
| Landsat k10, 100, T_ps | 500 | −1,334, −1,419, −1,457, −1,471 | 1.003 | 39 |
| | 5,000 | −1,485, −1,482, −1,483, −1,485 | 1.11 | 0.73 |
| AnisotropicGaussian k10, 100, T_ps | 500 | −235, −261, −274, −283 | 1.02 | 18 |
| | 5,000 | −300, −298, −297, −299 | 1.49 | 0.82 |
| | 20,000 | −298, −299, −300, −302 | 1.54 | 0.69 |
| AnisotropicGaussian k10, 20, T_ps | 500 | −307, −307, −306, −307 | 1.03 | 0.84 |

- *The split.* At burn-in 500 every replica drifts together, so R̂ stays quiet while the curve is wrong, and closure sees it. At burn-in 5,000 closure's median over the window is 0.87–1.06, and R̂ flags *rise* from 10–16 to 26–29 of 36 temperatures, because the replicas have
  settled into different long-lived states. R̂ got worse as the sweep got better. This is C46's hole, measured: "R̂ can miss a shared metastable phase". Neither check alone would have told us.
- *The mechanism fits C46.* A random start opens an edge only if its ends share a colour, so its open degree is at most K̂/q, 0.05–0.08 at q = 100 and 0.27–0.40 at q = 20 on these k = 10 graphs. The disordered phase can't percolate from there, and order waits on local nuclei.
- *What moved.* T½/T_ps at q = 100: AnisotropicGaussian 0.864 to 0.939, BlattHierarchy 0.939 to 0.972, Landsat 0.576 to 0.593. The C_v peak moved 1.08–1.36 to 1.28–1.53. The burn-in-500 peaks at 0.05–0.07·T_ps were drift.
- *C43 caught it too.* After recomputing the lower envelope on 20,000 draws (4,000 above N = 2,000), the misses at q ≤ 20 are gone: they were the panel's 200-draw estimate on one stream, and the worst z left is 4 on points selected for being low, across 8,481 comparisons. At q = 100 Landsat sits at S₂ = 0.9985–0.9990 against a bound of 1.0000 at nine cold temperatures, z = 3.4–5.5, all independent runs. A theorem about the stationary measure, breached, means the chain isn't there.
- *Controls.* q = 20 at burn-in 5,000, same window: Landsat T½ 0.609 to 0.607, C_v peak unchanged. The toy's T½ lies below the window, its S₂(T_lo) is 0.8592 to 0.8596, and its C_v peak moves one grid step. Closure and R̂ flags barely move. So burn-in 500 holds at q = 20, tested at k = 10 on the two largest graphs and elsewhere only through closure's median staying near 1.
- *Still not mixing.* At burn-in 5,000 and 20,000, R̂ at q = 100 reaches 1.3–3.3 near and above the melt, and Landsat at 1.5·T_ps still drifts. Its numbers above are better, not certified. §13.18 and §13.19 ran q = 40 at burn-in 200, which nothing here traced: if a q arm there reads oddly, burn-in comes first.

R̂ also tracks K. At q = 20, flags total 2, 28, 59 at k = 5, 10, 20 and the largest R̂ runs 1.04–1.32, 1.08–1.61, 1.30–2.33. K sets how first-order the transition looks to the sampler, and the sampler's load is the graph's.

**The stop rule (C63).** As written, a parent is quiet or has ρ in [¼, ¾].
- *As written it fails.* S₂ stops on 7 of 24 graphs within-run and 11 under Welch, the information currency on 17, F_n on 19. It also reverses: of the stops with a later volley to check, 5 of 7 and 6 of 10 fail at a later one. ρ has no noise term, and at v5 one noisy parent in 128 sinks the volley: the toy at k = 10 stops failing only on 3 unresolved parents out of 128.
- *With a band (exploration, written after that failure).* A parent fails only if a jackknife band over replicas, t on 7 df, Bonferroni over the volley, excludes [¼, ¾]. S₂ then stops on 23 of 24: 4 at v2, 7 at v3, 6 at v4, 6 at v5. It reverses on 2 of the 17 with a later volley. S₂'s piecewise-linear error against the full grid at the stop has median 0.013, largest 0.053.
- *The density needed moves with the graph.* The stop spans v2 to v5, a factor of 8. It tends to rise with K, median v3 at k = 5 and v4 or v5 at k = 10 and 20. So the user's "a few volleys" holds on some graphs and not others, which is the abstract minimum density of §13.4 showing up.
- *An earlier stop isn't a better one.* The information currency held expectation 8, and with the band it stops on all 24, 10 of them at v2. But it reverses on 6 of 22, and at its stop S₂'s curve is off by up to 0.106. That's less sensitivity as much as less work. Open question 4 isn't settled by stop volleys.
- *Not tested.* Whether a BARS fit at the stop is one to trust. Interpolation error stands in for it here, and BARS smooths.

**Closure.**
- *C41's h² holds on energy.* Median trapezoid/ΔE on not-quiet intervals, across graphs at v2, v3, v4, v5: 1.259, 1.062, 1.007, 0.994 at q = 20, and 1.196, 1.056, 1.014, 0.998 at q = 5. The first halvings shrink the bias by 3.5–4.2 times.
- *The 10% gate belongs to melts.* Pass shares at v4: Landsat 123 of 144, the toy 152 of 207, AnisotropicGaussian 105 of 199, BlattHierarchy 55 of 176. On Cauchy BlattHierarchy the failures scatter both ways, 0.52–1.68, and refining makes them no better: 15 of 93 at v5 for k = 5. A bridge moves S₂ a lot and E barely, so Var(E) rests on a few flips. Closure needs a noise term there, and on bridges it's blind by design: it sees energy.
- *Var(m)'s closure doesn't work at this budget.* It passes on about one interval in ten at any volley, and its median trapezoid/Δ at v5 is 0.56–1.04, below 1 on 22 of 23 graphs at q = 20, nearer 1 at q = 5 (0.91–0.99). Cause not chased.

**Smaller readings.**
- *The outlyingness test is at chance.* z > 5 on 8 normal replicas fires 4.7% of the time by simulation, so 12 of 257. The panel reads 7–25 at q = 20 and 9–27 at q = 5, then 18–124 at q = 100. §13.6's threshold doesn't discriminate at this replica count.
- *The certificate.* 0 of 1,415 quiet certifications contradicted at v3 and 0 of 6,106 at v5, half the replicas against the other half, over all 33 cells. Four replicas make this a weak test, and it's the same seed.
- *Levels are the graph's, or the data's.* AnisotropicGaussian reads 0.998 on all 12, with its window 8.0–10.1 e-folds under Gaussian (open to the grid's cold end) and 2.3–3.1 under Cauchy. The toy reads 0.69–0.71 and Landsat 0.21–0.28 on every K and kernel, so neither moves their best ARI at θ = 0.5 on this recipe. That extends §12.29 to k = 20 and Cauchy. §13.18 says the filter is what moves Landsat.
- *Spread and the hot edge.* ln(s_p90/s_p10) doesn't predict the width, but it does track the hot edge over its anchor: Spearman 0.70, n = 24, exploratory. It fits the hottest structure living on the strongest nodes.
- *The compiler's floor (not changed).* On Gaussian k = 5 the 10⁻¹⁰ floor cuts one mutual-KNN edge that is a graph bridge, on the toy (345–671, d/a = 7.52, κ = 5.25×10⁻¹³, leaving 3,199 + 1) and AnisotropicGaussian (307–863, d/a = 7.15, κ = 8.05×10⁻¹², leaving 800 + 400). It logs this. "Repair edges are exempt" doesn't cover a natural bridge. J about 10⁻¹² can't be felt on this grid, so the sweeps don't care. Also stale: `RefinementKind.Euclidean`'s entry at `GraphCompilerConfig.cs:154` and the comment at `GraphCompiler.cs:159`, since it throws at `GraphCompiler.cs:436`.

**Corrections.**
- §13.17's "s̄ runs 0.603–0.643": on the Gaussian panel it's 0.603–0.710, at BlattHierarchy k = 5, whose long repair edges raise a. s̄ = ⟨κ⟩ stays exact, and ≈ e^(−1/2) is what loosens.
- §13.14's melt at 0.87·T_ps is the cloud's T½. Across q the form held, and tighter once q = 100 was equilibrated.
- §13.14's "no sign at round 1's sizes" for C46 still holds at q = 20 on this panel and doesn't at q = 100.
- My first analysis scored E10 on all 24 graphs. The preregistration scoped it to single-melt graphs, and that's the 0.46 above.

**What it means for the stack (§13.5). A recommendation, not a decision.**

| Stage | After E20 |
|---|---|
| Equilibration | Two checks, because they see different things: R̂ for stuck replicas, energy closure for shared drift. Drop outlyingness at z > 5 on 8 replicas |
| Endpoints | Untouched. Undecided bands run 7–72 grid points at q = 20, by graph |
| Quiet intervals | Stands, on a weak test |
| Active intervals | The rule as written is out. ρ with a band is the candidate. Closure is a resolution diagnostic on melts only |

Four of the things that failed were mine: T½ as the melt, Var(m)'s closure, the outlyingness threshold, and ρ without a noise term. The user's call on the uniform grid stands up, and its anchors are the part that isn't settled. The cheap read up front,
which the user asked about, is better than I expected: the envelopes place T½ and flag missing cold structure from the graph alone.

**E21. The band rule on a fresh seed, and BARS at the stop (not run).** Fix the band rule and its thresholds in advance. On a new seed at the panel's budget, on the 24 graphs, run volleys to the stop. Fit BARS to S₂ at the stop with its settings stated, since the settings pass hasn't happened. Against the panel's 257-point curve as the independent reference: does the fit find the same features at the same temperatures, and does the rule reverse? This is the experiment that answers the question round 2 set: when a fit on a grid can be trusted to say whether features exist and where. About 3–4 hours on this machine. Off-recipe graphs and E16 remain separate.
