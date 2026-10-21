
### Verdict on the Core Recommendation (Q5)

**Approve the recommendation:** adopt the **likelihood-tempered path** as the primary path for Step `3.2`, with the Step `3.3` round-trip efficiency check recorded as the revisit trigger.

Point 1 in the note is the clincher, and it is worth stating the exact measure-theoretic reason in your docs: **a Markov chain targets a measure $\Pi(dx)$, not a density function**. Raising a density $f(x) = \frac{d\Pi}{d\mu}$ to a power $\beta$ defines $d\Pi_\beta \propto f(x)^\beta d\mu$, which is only invariant under reparameterizations with unit Jacobian ($\vert{}J\vert{} = 1$). 
Because the map from labeled knots on $[0, 1]^k$ to ordered knots on the simplex $\Delta^k$ has Jacobian $k!$, naive powered tempering leaks a $(k!)^{\beta - 1}$ factor across dimensions—turning $\text{Poisson}(\lambda)$ into $\text{Poisson}(\lambda^\beta)$ ($1.08$ knots) in ordered coordinates and into an overdispersed $\text{CMP}(\lambda^\beta, \beta)$ ($14.17$ knots) in labeled coordinates. Likelihood tempering defines the hot measure via $\frac{d\Pi_b}{d\Pi_0} = L(y \mid x)^b$ directly against the **prior measure $\Pi_0$**, so it is strictly representation-invariant.

### Decisions on the Three Open Questions

#### 1. Default or Required?

* **In `ParallelTempering` (`Maths.Samplers`):** **Required (no default)**. The generic sampler shell has no business guessing a thermodynamic path.

* **In `BarsConfig` (`Maths.Regression.Spline.Bars`):** **Default to `Likelihood`, with strict fail-fast validation in the constructor** (never silent fallback).

* Why not make it required on `BarsConfig`? Because Q5's role in the brief is to establish BARS's canonical tempering semantics for Normal regression.

* **Watch out for the Step `0.5d` trap:** In Step `0.5d`, you noted that `TemperingTests` includes a **robust five-level tempered run** (where `WeightResampler` updates non-unit weights $w_i \ne 1$). Since the $O(1)$ likelihood-tempered block only applies to unit-weight Normal models, the ensemble constructor should reject `TemperPath.Likelihood` whenever a `WeightResampler` or a non-Normal `IObservationModel` (`Poisson`, `Whittle`, `GibbsLoss`) is supplied, requiring those callers (including `TemperingTests`' robust case) to pass `TemperPath.Powered` explicitly.


#### 2. Proper $\sigma^2 \sim \text{IG}(g_1, g_2)$ Prior: Flag Now or Park?

* **Parameterize it inside `CollapsedGPriorBlock` right now (Step `2.3`), but keep BARS's default parked at $(g_1 = 0, g_2 = 0)$**.

* Look at the algebra of the likelihood-tempered Normal marginal with prior $\sigma^2 \sim \text{IG}(g_1, g_2)$ (where $g_1 = g_2 = 0$ is BARS's improper $1/\sigma^2$ prior):

$$\log p_b(y \mid k, \xi) = \text{const} + \frac{\nu}{2}\log\!\left(\frac{\tau}{\tau + b}\right) + \log\Gamma\!\left(g_1 + \frac{b m}{2}\right) - \left(g_1 + \frac{b m}{2}\right)\log\!\left(g_2 + \frac{S_b}{2}\right)$$

* Adding `double priorShape = 0, double priorScale = 0` to `CollapsedGPriorBlock` in Step `2.3` costs two fields and one guard (`if (b == 0 && priorShape == 0) reject`):

* At `(0, 0)` and $b > 0$, it reproduces BARS's Step `0.6` baseline to the bit.

* When BASS and the SMC core arrive in Stages 3–4 (which require $b = 0$ and proper $\text{IG}(g_1, g_2)$ priors), `CollapsedGPriorBlock` is already complete and never has to be reopened.

* As for changing BARS's user-facing default away from $(0, 0)$: leave that parked with the BASS prior decisions, since PT never visits $b = 0$ (`temperBetaMin > 0`).

#### 3. Tempered Marginals for Poisson and Whittle?

* **Do NOT build $b$-scaled IRLS Laplace marginals for Poisson, Whittle, or GibbsLoss**.

* Running a separate IRLS loop at every temperature rung $b_\ell$ would multiply per-step cost by $L\times$ and has zero code path into BASS.

* **However, notice a free win from Step `3.1` (`IRjTarget` splitting `LogPrior` and `LogLikelihood`):**
* Look at *why* the powered path failed in Point 1: it multiplied `LogPrior` ($\log p(k) + \log k!$) by $\beta$.

* Once Step `3.1` separates `LogPrior` from `LogLikelihood` (the collapsed Laplace/Gibbs marginal $\log p(y \mid k, \xi)$) on `IRjTarget`, a replica that evaluates `LogPrior + b * LogLikelihood` (tempering the *collapsed marginal likelihood* while leaving the knot prior untempered) is **still $O(1)$ per swap**, requires **zero extra IRLS solves**, and completely fixes the Point 1 representation defect for Poisson, Whittle, GibbsLoss, and weighted Normal runs!

* Even if you keep `Powered` as `b * (LogPrior + LogLikelihood)` for backward compatibility in Step `3.2`, keeping `LogPrior` and `LogLikelihood` split on `IRjTarget` in Step `3.1` gives you that option for free.

### One Extra Note on Point 2 & Point 5 (`UniformPrior(50)` on Hot Rungs)

Point 2 worries that under `UniformPrior(50)`, likelihood-tempered hot rungs wander to $\bar{k} \approx 25$ knots. Notice how Point 5 (`temperBetaMin = 0.05`, where $b m = 6$ at $m = 120$) actually protects you from Point 2 in Parallel Tempering:

* On the likelihood path, the $g$-prior's Occam factor per knot at temperature $b$ is $(1 + b/\tau)^{-1/2} = (1 + bm)^{-1/2}$.
* At $b = 0$, that factor is $1$, so `UniformPrior(50)` indeed wanders to $25$ knots.


* At $b_{\min} = 0.05$ with $m = 120$, $(1 + bm)^{-1/2} = 7^{-1/2} \approx 0.378$ per knot—which still imposes a $\approx 1$ nat penalty per redundant knot, keeping hot-rung dimension way below $25$ as long as $b_{\min} \gtrsim 1/m$! That gives Step `3.3`'s ladder adaptation a very clear lower-bound regime ($b_{\min} \sim 1/m$) to watch.