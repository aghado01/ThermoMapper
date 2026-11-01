# BASS — estimators that stay valid whatever the emulator does

Grounds a design-based layer beside the emulators: it takes predictions from any of them and reports means whose validity comes from how the runs were drawn. It is general machinery for any expensive code with cheap surrogates (multifidelity Monte Carlo in UQ, the regression estimator in survey sampling, prediction-powered inference in ML). Linked from the [lineage-gaps register](../../../issues/regression/bass/lineage-gaps-register.md).

**Status (2026-11-01).** Posed. [`lineage-gap-checks.py`](lineage-gap-checks.py) checks each item by simulation. The Lean is uncompiled.

**Symbols.** A pool of N_p configurations and n_r runs. These stay apart from BASS's N (cases) on purpose. r = m − f̂ is a fit's error against the means m.

---

## GE-1 — cross-fitting on a finite pool has an exact bias

**Claim.** Let the runs be a simple random sample of a fixed pool. For a fold with training set A, fit f̂ on A, and evaluate on the remaining runs. Given A and its outcomes, the fold's remaining runs are a simple random sample of the pool outside A, not of the whole pool, so

E[mean_pool f̂ + mean_fold(y − f̂) | A] − mean_pool m = (|A|/N_p)·(mean_{pool∖A} r − mean_A r).

**Reading.** The bias is the training share times the gap between the fit's mean error off and on its training points. A fit whose mean error is about the same on and off (a ridge) shows nothing to Monte Carlo error at any n_r/N_p. One built to fit its training points and fall away elsewhere shows it plainly:

| spike fit: pool, runs | bias | formula |
|---|---|---|
| 60, 30 | +0.786 ± 0.004 | +0.783 |
| 60, 12 | +0.346 ± 0.005 | +0.347 |
| 600, 12 | +0.033 ± 0.005 | +0.034 |
| 6000, 12 | −0.002 ± 0.005 | +0.003 |

So "unbiased whatever f̂ is" holds exactly for an f̂ fitted on a separate pilot, and only up to this term under cross-fitting on the pool.

**Lean.** The conditional identity is algebra, once the simple-random-sampling expectation is taken.

```lean
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Real.Basic

open Finset

theorem crossfit_gap {P : Type*} [Fintype P] [DecidableEq P] (m f : P → ℝ) (A : Finset P)
    (hA : A.Nonempty) (hAc : Aᶜ.Nonempty) :
    (∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ Aᶜ, (m i - f i)) / (Aᶜ.card : ℝ)
        - (∑ i, m i) / (Fintype.card P : ℝ) =
      ((A.card : ℝ) / (Fintype.card P : ℝ)) *
        ((∑ i ∈ Aᶜ, (m i - f i)) / (Aᶜ.card : ℝ) - (∑ i ∈ A, (m i - f i)) / (A.card : ℝ)) := by
  sorry
```

---

## GE-2 — draws from the generator, with the integral, have none

**Claim.** If run configurations are independent draws from a declared measure μ, and the first term is ∫f̂ dμ (closed form from `BassMoments`, or an independent pool), then given a fold's training runs the remaining runs are draws from μ itself, so each fold is exactly unbiased for ∫m dμ, for any fit. The check reads +0.004 ± 0.005 and +0.004 ± 0.004 on the spike fit.

**Why it matters.** The closed-form integral is what makes the estimator exact, rather than a separate convenience. This is a concrete return on the product-measure closed forms beyond Sobol indices, and [SG-3](sensitivity-generator-lemmas.md)'s joint mixtures widen the measures it covers.

---

## GE-3 — telescoping levels, the tuned weight, known probabilities

**Telescoping.** For an emulator f̂, a cheap level Y₁ on n_r runs and a costly level Y₂ on a random m of those runs,

∫f̂ dμ + mean_{n_r}(Y₁ − f̂) + mean_m(Y₂ − Y₁)

is unbiased for E[Y₂], provided each subsample is random and each difference's model is fitted apart from the sample it corrects. A fixed Y₁ needs no cross-fitting; one calibrated on Y₂ does. The check reads +0.0007 ± 0.0008. Multilevel and multifidelity Monte Carlo (Giles 2008; Peherstorfer, Willcox & Gunzburger 2016, 2018) give the allocation of n_r and m from each level's cost and the correlation between neighbouring levels.

**The tuned weight.** With a fixed f̂ and a pool independent of the runs, the variance of λ·mean_pool f̂ + mean_runs(y − λf̂) is λ²Var(f̂)/N_p + Var(y − λf̂)/n_r, minimized at

λ* = Cov(y, f̂) / ((1 + n_r/N_p)·Var(f̂)).

λ = 1 is the difference estimator, λ = 0 the plain mean, and the tuned weight does no worse than either (asymptotically, once estimated); a poor emulator can make λ = 1 worse than the plain mean. This is survey sampling's regression estimator and PPI++ (Angelopoulos, Duchi & Zrnic 2023). The check's grid minimizer is 1.045 against λ* = 1.0446.

**Known probabilities.** With each configuration run independently with known π_t > 0 and f̂ fitted apart,

mean_pool f̂ + (1/N_p)·Σ_{t run}(y_t − f̂_t)/π_t

is unbiased (Horvitz–Thompson). The check reads +0.0002 ± 0.0004. This is what reconciles the estimator with an adaptive design, which otherwise breaks the random-sample premise. An audit stream (a fixed random share of runs) is the simple alternative. Keep both, in docs/Maths when built.
- *Audit stream:* simple, but those runs do less to locate boundaries.
- *Known probabilities (active statistical inference, Zrnic & Candès 2024):* every run serves both purposes, but the variance grows as the floor on π_t drops. The variance-optimal π_t is proportional to the root-mean-square residual expected at the configuration, which a heteroscedastic emulator's predictive sd supplies. Where π_t depends on earlier runs the argument is a martingale one, not the independent-draw one above.

**Scope.** These guarantees cover means and reweighted means, and through PPI regression coefficients and quantiles. Sobol shares, boundaries and subspaces rest on the emulator being right, and are checked by held-out prediction and posterior predictive checks instead.
