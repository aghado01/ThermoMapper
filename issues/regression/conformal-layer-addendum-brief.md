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

### In progress — Slice 16, batch predictive inference (started 2026-11-14)

At the user's pick (2026-11-14, "Slice 16, batches"), under the Lean's section `Batch` (proved 2026-11-10) and Lee, Tchetgen Tchetgen and Dobriban (arXiv 2409.13990v5, read in TeX: §2, Condition 1, Algorithm 1, Theorem 1, the quantile corollary, the selection corollary, the simultaneous bounds).

**Plan.** A `BatchConformal` class on calibration scores, scalar units, the scores' bounds S_(0) = inf s and S_(n+1) = sup s given or ±∞:
- *The quantile of a batch* (their quantile corollary; `batch_interval` with r_(ζ)): [S_(q_L − 1), S_(q_U)] from the law of r_(ζ). Its counts are exact big integers, so m = 1 is split conformal bit for bit.
- *Simultaneous bounds* on several order statistics of the batch (their joint-bounds theorem; `batch_simultaneous`). F, the box's count, is exact by dynamic programming on L(ρ), with the product over consecutive pairs (the Lean's correction of their typo). Two readings: a given box with its coverage, and their centered search.
- *Functions of a batch* (Algorithm 1, `batch_coverage`), by enumerating H where C(n + m, m) is small, since it is NP-hard in general. The rank order h̃ is given, h on the ranks by default.
- *Selection with false claims controlled* (their selection corollary; `batch_exceed`, `false_claim_iff`): the threshold T̂ = S_(q_η) on the scores μ̂·1{y ≤ c}, with at most η false claims with chance at least 1 − α.
- *Gates* (the archived plan's): each guarantee on exchangeable batches, the counts against enumeration of H, and width against Bonferroni over single sets.

**State.** Being written. Nothing committed.

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
