# The conformal layer: addendum (chip brief)

Opened 2026-11-14 at the user's word ("that brief was getting way too long, addendum brief"). The [conformal-layer brief](archive/bass-gbass-khaos-activespaces/briefs/conformal-layer-brief.md), now archived, is the record of E16 through slice 15 and slice 14's follow-ups: its decisions 1 to 12, the plan for slices 7 to 18, the gates common to every slice, the working rules and every report. They all stand. This addendum carries what is open and the reports from here on.

**Chip:** the conformal chip, `local_69bae1e5`. **Verified by:** the main session, `local_b6732445`, from clean archives (`scripts/verify-csharp.ps1`, `scripts/verify-qualification.ps1` for gated tests).

## Where E16 stands

- Slices 0 to 15 are built and verified; slice 14's follow-ups 1 to 3 and its ratio estimator too, on their Lean (decision 12).
- Follow-up 2 is completed by `3f68330`, the merged K-fold CV+ bound under likelihood ratios, on the Lean at `47acb3c`. Its report is the first below.
- The defects chip owns `PredictiveChecks.cs` and `Pointwise.cs`'s weight check; a change there is agreed with it first (the archived brief's coordination note).

## Open, each with the decision it would change

- **Slices 16 to 18,** by the user's need. Slice 18's plan reads Snell and Griffiths's bound as Gaffke's confidence bound, never a posterior (corrected 2026-11-13).
- **Slice 15's precision rule** reads the training split's ESS: with 1,000 training and 399 calibration points it chose m of 294 to 500, and 7% of the sets were the whole space. It decides whether the fit takes the calibration split's size and scales the target itself.
- **Slice 15's alternation speed** along the factorization's flat directions. It decides the cost at large n.
- **τ through richer features,** the ratio estimator's open follow-up. It decides whether the layer offers a nonparametric τ.
- **Two follow-ups agreed with the GBASS chip** for E21's slice 5 (`127d64d`): new-group coverage by group-level units, and the layer's own check that strata lie inside cells.
- **Slice 12a's k̂ diagnostic,** opt-in or on by default (the handoff's waiting list).

## Working rules

As in the archived brief: each slice is reported here, then waits for the main session's verification before the next starts; anything that changes a default or removes an option goes to the user; chips leave `docs/` alone and each report ends with **Doc changes owed**.

## Reports

Newest first.

### 2026-11-14 — Slice 16, batch predictive inference

Built at the user's pick (2026-11-14, "Slice 16, batches"), under the Lean's section `Batch` (proved 2026-11-10) and Lee, Tchetgen Tchetgen and Dobriban (arXiv 2409.13990v5, read in TeX: §2, Condition 1, Algorithm 1, Theorem 1, the quantile and selection corollaries, the joint-bounds theorem and corollary).

**Landed** on `main`, nothing pushed:
- `60ba2cb`: `BatchConformal`, with
  - `RankLaw` and `Quantile`, giving a `BatchInterval`;
  - `BoxCoverage`, `Box` and `Simultaneous`, giving `BatchBounds`;
  - `Function` and `Select`, the latter giving a `BatchSelection`;
  - and `BatchConformalTests`.
- *Nothing existing changes.*

**What it computes.** Batch sets for functions of m test units' scores from n calibration scores, exchangeable together. The scores' floor and ceiling, S_(0) and S_(n+1), are ±∞ unless given.
- *The basis.* The test units' sorted ranks are uniform on H = {1 ≤ r_1 ≤ … ≤ r_m ≤ n + 1}, |H| = C(n + m, m), and every sorted test score lies between S_(r_j − 1) and S_(r_j) at once (`card_testRank_mem`, `batch_sandwich`). Each set reads that sandwich on rank sequences holding at least 1 − α of H.
- *`Quantile`,* the ζ-th smallest score of the batch: [S_(q_L − 1), S_(q_U)], from the exact law of r_(ζ) (their quantile corollary). q_U is the least rank with at most γ of the law above it; q_L the greatest with at most β below. At m = 1 it is split conformal's set.
- *`Box` and `BoxCoverage`,* several order statistics at once: [S_(w_j − 1), S_(q_j)] with the box's share of H (their joint-bounds theorem, `batch_simultaneous`). The count runs by dynamic programming on L(ρ), with the product over consecutive orders, the Lean's correction of their closed form.
- *`Simultaneous`:* their centred search, each order's box centred at round(t_j·n/m) and widened by the least a whose share reaches 1 − α, by bisection.
- *`Function`,* any h monotone in each sorted score: their Algorithm 1, with the rank order h̃ given (h on the ranks by default). It enumerates H twice, refused above 2,000,000 sequences, since computing it is NP-hard in general.
- *`Select`:* scores μ̂·1{y ≤ c}, the threshold T̂ = S_(q_η), and the units with μ̂ above it declared above c, at most η of them falsely with chance at least 1 − α (their selection corollary; `batch_exceed`, `false_claim_iff`).
- *Exact levels.* The rank counts are integers (`BigInteger`). A share of H at most γ is the count at most γ·|H| for the given γ read exactly, so whole levels such as (1 − γ)(n + 1) = 18 land on their rank and do not round away from it.

**Gates.**
- *The counts are enumerations of H:*
  - the law of r_(ζ) at every k, for n ≤ 6 and m ≤ 4;
  - 300 random boxes' shares, to 10⁻¹².
- *At one test unit:* the quantile set is split conformal's, q_U = ⌈(1 − γ)(n + 1)⌉, across five (n, γ), whole levels included.
- *The batch median* (n 499, m 20, β = γ = 0.05, 4,000 batches): it covers 0.900, within its bound 0.9 + max p = 0.907, with mean width 0.94. Bonferroni over single sets covers 1.000 at mean width 6.07.
- *The quartiles at once* (n 200, m 40, α 0.1, 3,000 batches): the box holds 0.911 of H at ranks [21, 79] and [121, 179] of 201, and covers 0.914.
- *The batch mean* (n 12, m 3, scores in [0, 1], 4,000 batches): it covers 0.983 at mean width 0.78 of the range. With unbounded scores the mean's set is the whole line, since a sequence reaching S_(0) = −∞ is always held.
- *`Function` with h the second smallest and h̃ = r_(2)* is the quantile set exactly.
- *Selection* (n 200, m 50, μ̂(x) = x, c 0.7, η 2, α 0.1, 3,000 batches): false claims stay within 2 in 0.916 of the batches, with 4.7 units selected on average.
- *The refusals:*
  - a batch of none, an order outside 1, …, m, levels summing to 1 or more;
  - a NaN score, scores outside the floor or ceiling;
  - orders that do not rise, a box out of order;
  - an H above the enumeration cap;
  - negative predictions, η = m.
- *Per commit,* the selected classes and the pinned tier: 49 of 49 (one class, the new one).

**Defects found: none.**

**Follow-ups, each with the decision it would change.**
- *The rank order h̃ for the mean.* The default, h on the ranks, is conservative here (0.983 at nominal 0.9). Their remark on choosing h̃ suggests orders closer to h's own. Decides whether `Function` offers a better default for sums and means.
- *Sparse functions* (their §5: h reading a few order statistics), as `Quantile` reads one. It decides whether functions of a few order statistics avoid the enumeration.
- *Their upper bound 1 − α + ε* for the quantile, not formalized. The Lean states the lower bound only. Decides whether the layer reports the set's tightness.

**Doc changes owed** (for the offline docs pass): `docs/Maths/Regression/ConformalPrediction.md`, a batch section:
- the sandwich and the four forms, with the Lean names;
- exact levels;
- the gates' numbers;
- API rows for `BatchConformal`, `BatchInterval`, `BatchBounds` and `BatchSelection`.

### 2026-11-14 — Slice 14's follow-up 2, completed: weighted CV+'s merged bound under likelihood ratios

Built at the main session's word once its Lean (`47acb3c`; `cvPlus_merged_likelihoodRatio`, `cvPlus_merged_likelihoodRatio_refined`) was verified. It closes the open item in follow-up 2's report.

**Landed** on `main`, nothing pushed:
- `3f68330`: weighted CV+ under likelihood ratios reports a merged bound (`CvPlusBounds.Merged`), where it reported none. `WeightedCrossFitTests` gains its checks.
- *Verified by the main session on 2026-11-14,* from a clean archive of `1ea5d6b` with zero skips: the Commit selection, 1,213 of 1,213 ([receipts](receipts/verify-conformal-fu2e-3f68330/)).
- *Nothing else changes.* Position weights' 2c stands: its Lean was restated validly at `47acb3c`, with the same conclusion.

**What it computes,** per test point, from the bag's weights (the training units' and the test point's own, W their total):
- *simple:* min(1, 2c(w_max)), with c(γ) = (αW + (K − 1)γ)/(W + (K − 1)γ), K the folds and w_max the bag's largest weight;
- *refined:* γ the largest weight at most αW, L = αW + (K − 1)γ, N = W + (K − 1)γ and e = Σw²/W. The miss is 2L(W + (K − 1)e)/N² where 2L ≤ N, and 1 otherwise;
- where no weight is at most αW, no unit can be refused: miss 0;
- the merged coverage is 1 minus the lesser miss. It applies under `MergedRanks`, as slice 9's merged bound does.
- `Coverage` reports the largest of the three bounds that apply: one fold's, across folds', and merged.

**Gates.**
- *By hand* (three folds of weights 3, 2 and 4, the test's 2, α 0.2): the simple 2·8.2/17, the refined at γ = 2 (L = 6.2, N = 15, e = 21/11), the merged 1 minus their lesser, and the coverage the largest of the three bounds. Every weight above αW gives a merged bound of 1.
- *Covariate shift* (follow-up 2's fixture, 10 folds of 10 units): the merged bound averages 0.553 against one fold's 0.282 over 2,000 replications, while weighted CV+ covers 0.914; the across-folds bound, given w̄ = 2, is about 0.6 there (1 − 2α − 20/W with W near 101, computed rather than read from the run).
- *Existing classes,* the CV+ tests: unchanged and passing.
- *Per commit,* the selected classes and the pinned tier: 1,213 of 1,213.

**Defects found: none.**

**Doc changes owed:** the cross-fitted section's weighted CV+ bounds gain the merged bound under likelihood ratios, its two forms and where it serves (many units), with the Lean's table of which bound serves where.
