I'll measure how often the two rules disagree on BARS's real candidates by temporarily wrapping `SplineTarget` to record both statistics of G = ZᵀZ across every evaluated configuration, testing dense uniform x, x with a data gap, and small n—data gaps being where divergence is most likely.

The frequency data shows disagreements always run the same way—PivotRatio rejects while the relative pivot admits—suggesting PivotRatio may also be screening out thinly supported basis functions near the edge of their data range. I'll inspect what these rejected configurations do to the fit in the gap case to test that hypothesis.

The data confirms thin-support configurations have wildly inflated posterior-mean fits (hundreds to thousands) regardless of whether PivotRatio's 1e-10 threshold admits or rejects them, showing the cutoff splits statistically similar cases inconsistently. I'll now check BASS's `npart` argument in the pinned source and measure it on these configurations.

I've confirmed that BASS enforces `npart` (minimum non-zero points per basis function, defaulting to min(20, 0.1·n)) as an explicit restriction in every proposal function, unlike BARS which hides this inside PivotRatio. Next, I'll measure this same minimum-points-per-basis-function count on the probe's configurations to pinpoint where wild fits start occurring.

Fits get wild with fewer than ~9 points under a basis function but stay sane at 9+, and BASS's default min(20, 0.1·n) threshold would have filtered out every wild case here—so I have enough evidence to remove the probe.

I'd switch the numerical rule to `MinRelativePivot`, but not on its own. The measurements show `PivotRatio` has been doing two jobs, and dropping it would lose the second one quietly.

## What the probe found

I recorded both statistics on all 6,000 Grams a BARS chain evaluated, on each of four datasets:

| Data | Rejected by PivotRatio only | By the relative pivot only | By both |
|---|---|---|---|
| Dense x, n = 150 | 0 | 0 | 0 |
| x with a gap, n = 120 | 3 | 0 | 0 |
| Small, n = 40 | 2 | 0 | 1 |

Disagreements are rare, and they always run one way: PivotRatio rejects what the relative pivot admits. The one case both reject is a genuinely collinear column, with relative pivot 2.7e-12.

The configurations PivotRatio alone rejects are not ill-conditioned: their relative pivots are 0.1–0.2. They are **thinly supported basis functions**, with only one or two data points under them, at the edge of their support. They are numerically fine and statistically wild. Against data ranging 0–2, the conditional-mean fit inside the gap reached 290–1,310. The reason is that the g-prior's covariance mσ²(ZᵀZ)⁻¹ puts no constraint on a direction with no data.

PivotRatio's 1e-10 cut doesn't separate wild from sane either. It falls in the middle of a continuum: configurations at 2.8e-11 (rejected) and 1.4e-10 to 3.7e-9 (admitted) are equally wild.

## Why the numerical check should be the relative pivot

- **Accuracy follows the equilibrated condition number.** A Cholesky solve's accuracy is governed by κ of D^(−½)AD^(−½), not κ(A) (van der Sluis 1969; Demmel 1989; Higham, *Accuracy and Stability*, ch. 10). The relative pivot bounds exactly that quantity; PivotRatio bounds κ(A), which includes column scaling Cholesky handles exactly.
- **Consistency with the model.** The g-prior marginal depends only on Z's column space, so the model ignores column scaling. A rule that depends on it makes the posterior depend on how the basis happens to be normalized.
- **Weight-dependent truncation.** Robust weights rescale Pa = Zᵀ(W + I/m)Z, and Poisson's Hessian is weighted by μ̂. Under PivotRatio, admission therefore depends on the weights and on the count level. That's a prediction I haven't measured; it would matter most for low-rate regions of spike-train fits.

## The second job, made explicit

BASS already declares it. Its `npart` is the "minimum number of non-zero points in a basis function" (`bass.R:23`), with default min(20, 0.1·n) (`bass.R:259–261`), and every proposal checks it: `rjmcmc_des_*.R` reject a candidate when `sum(cand.des$basis != 0) < npart.des`.

On these runs, the fewest data points under any basis function tracked the wild fits:

| Fewest points under a basis function | Max fit in the gap |
|---|---|
| 1–2 | often 30–1,300 |
| 3–7 | sometimes wild (226 at both 3 and 6 points) |
| 9 or more | always 1.2–1.9 |

BASS's default would have removed every wild configuration here. A shared `npart` also gives BARS and BASS one admission vocabulary.

There is a third, smaller piece. The `FlooredPivots == 0` condition rests on `BandCholesky`'s absolute 1e-12 floor, which is also scale-dependent; it hit once in 24,000 evaluations. A strict band mode would remove it, consistent with step 2.2.

## The proposal

Split admission into two declared parts:
1. **Numerical:** a strict factorization (no floor) that comes out positive definite, with `MinRelativePivot` ≥ 1e-10. That keeps Q4's 6-digit rationale, now on the right condition number.
2. **Support:** at least `npart` data points under every basis function, as a declared prior restriction, configurable.

I'd do this after 2.3. Step 2.3 is a pure refactor checked against the 0.6 baseline, and changing admission first would move that baseline. Afterwards it becomes its own recorded change, measured before and after on the validation fits, the posterior mean in data gaps, and rejection counts by reason.

## Where I'm least sure

- **I measured proposals, not the posterior.** I don't know how much posterior weight the wild configurations get. Their Occam cost is about ½·log(m + 1) per basis function, so the weight is probably small but not zero. A spike of 1,000 at weight 10⁻³ would still show in the gap's posterior mean.
- **`npart` is pragmatic, not principled.** The root cause is the g-prior's lack of shrinkage in data-free directions. The principled fix is a coefficient prior with a fixed scale component. That would break the identity "BARS is BASS with τ = 1/m" unless BASS changes too, so it belongs with BASS's prior decisions rather than here.

## For you to decide

- **The `npart` default for BARS.** BASS's min(20, 0.1·n) means 15 points per basis function at n = 150. Each cubic basis function spans four knot intervals, so that limits how close knots can be, and it could cost the sharp-peak fits. The probe suggests something around 8–10 removes the wild fits. I'd keep BASS's default for consistency and check the peak validation fits, rather than tune a BARS-specific number from this one dataset.
- **Whether you accept the split**, rather than a straight swap of PivotRatio for the relative pivot.