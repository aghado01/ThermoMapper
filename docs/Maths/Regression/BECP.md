# BECP — Bayesian Exact Change Point

> Piecewise-constant regression computed by exact dynamic programming across Gaussian or Cauchy segment models. Part of `Maths.Regression`.

### Exact changepoint regression

**Model and prior.** The number of segments is uniform on $1..K$, and segmentations are uniform given the count. The segment model fills tables of per-segment log evidence $\log D(a,b)$ and level moments.

**Forward and backward recursions.**

$$F_s(b)=\operatorname*{LSE}_{a}\big[F_{s-1}(a)+\log D(a+1,b)\big],\qquad B_s(a)=\operatorname*{LSE}_{c}\big[\log D(a+1,c)+B_{s-1}(c)\big].$$

**Outputs.**

- The segment-count posterior and its MAP; ties go to the smaller count.
- The log evidence, which includes the $-\ln K$ prior term.
- Per-position break probabilities: $n-1$ values, where entry $i-1$ is "break after position i".
- A posterior-mean regression curve with its standard deviation.
- `Posterior` (`ChangepointPosterior`): the dynamic program's tables, kept. They are each segment's log evidence, level moments and posterior probability $P(\text{segment }[a..b]\mid y)$, and the forward and backward tables. They cost $4(n+1)^2+2(K+1)(n+1)$ numbers.

**Exact segmentation draws.** `Posterior.Sample` draws independent, exact segmentations by forward filtering and backward sampling:

1. k from $P(k\mid y)$.
2. The segments from the right end: the start of the s-th segment, given its end b, takes $a$ with probability $F_{s-1}(a)\,D(a+1,b)/F_s(b)$. Given k the prior is uniform over segmentations, so this is the posterior exactly.
3. Each segment's level from the segment model: the conjugate Normal draw for the Gaussian model, and a quadrature node in proportion to its weight for the Cauchy model, so expectations over the draws equal the quadrature's.

The draws give what the marginals cannot: which breaks occur together, per-draw features such as the largest jump, and per-draw curves (`ChangepointDraw.Curve`). A brute-force enumeration of every segmentation on a short series checks the tables and the draws.

**Segment models.**

- **Gaussian:** closed form for $\log N(y_{a..b};\nu\mathbf1,\sigma^2I+\rho^2\mathbf1\mathbf1^\top)$.
- **Cauchy (outlier-robust):** quadrature over the level on a shared uniform grid, with the per-point log-likelihoods prefix-summed at every node.
  - **Range** $[\min(\nu,y_{\min})-6\rho,\ \max(\nu,y_{\max})+6\rho]$. Beyond it the prior and every Cauchy factor both decrease, so a segment loses at most the prior's $6\rho$ tail ratio, about $2\cdot10^{-9}$ of its mass per side.
  - **Spacing** $h=\tfrac12/\sqrt C$ with $C=1/\rho^2+2n/\gamma^2$. A Cauchy factor's log bends by at most $2/\gamma^2$ either way, so $C$ bounds the log integrand's curvature, and the narrowest segment posterior, about $\gamma/\sqrt{2n}$ wide, spans at least two nodes. The count $G$ is never below 512.
  - **Windowing.** Per segment the log integrand is read at every 16th node. With its curvature bounded by $1/\rho^2+2L/\gamma^2$ for a segment of L points, the reads bound the maximum between them, and only intervals that can come within 50 nats of it are summed node by node.
  - `gridPoints` fixes the count on the same range, and `Resolution` ($h\sqrt C$, at most ½ on the automatic grid) reports a coarse one. A scale too small for a table of $2^{25}$ numbers is refused.
- **Why the Cauchy grid changed (2026-10-24).** The former grid had 512 nodes on $\nu\pm6\rho$. Against a fine reference on integer steps (n = 60):
  - It was exact to $10^{-8}$ for $\rho/\gamma\le6$.
  - At $\rho/\gamma\approx22$ it was off by 0.08 in log evidence and 22% in level sd.
  - At $\rho/\gamma\approx110$, light noise on integer levels as in $b_1(T)$, it was off by 6.1 in log evidence, 4 posterior sds in level means, 0.017 in $P(k)$ and 3.7 in the total evidence.
  - A user-set narrow prior put a level beyond $\nu+6\rho$ and lost it: 132 nats, a level of 2.95 against 5.02, and break probabilities off by 0.17.
  - The new grid agrees with the reference to $4\cdot10^{-11}$ in log evidence on every case probed.
- **Alternatives considered:**
  - *The fixed grid:* fast, but wrong in the regime the engine serves.
  - *Per-segment adaptive integration:* without the shared prefix sums each node costs O(L), so O(n³) work.
  - *Gauss–Hermite at a Laplace centre:* misses the second mode of a segment that straddles a jump.

**Default hyperparameters.**

- $\nu=\bar y$.
- $\sigma^2=\max\big(\operatorname{med}_i[(y_{i+1}-y_i)^2]/(2\,\operatorname{med}\chi^2_1),\ 10^{-6}\big)$, with $2\,\operatorname{med}\chi^2_1=2\,\Phi^{-1}(\tfrac34)^2\approx0.909873$. Within a segment a difference is $N(0,2\sigma^2)$, and the median ignores the few that straddle a jump. With an even number of differences it takes the upper of the two middle values; $\sigma^2=10^{-6}$ when $n=1$. The divisor was the truncated 0.9098 until 2026-10-24.
- $\rho^2=\max(\text{var}-\sigma^2,\ 0.1\,\text{var}+10^{-9})$.

**Complexity.** Filling the segment tables costs $O(n^2)$ for the Gaussian model. The Cauchy model costs $O(Gn)$ for its prefix table, then per segment $G/16$ reads and its window. The dynamic programs cost $O(Kn^2)$ and the regression curve $O(n^2K^2)$, where $K=\min(\text{maxSegments},n)$. Measured Cauchy fills on integer steps:

| n | $\rho/\gamma$ | G | Former grid (inaccurate above about 20) | Now |
|---|---|---|---|---|
| 60 | 5.5 | 1,869 | 12 ms | 20 ms |
| 60 | 22 | 7,153 | 10 ms | 55 ms |
| 60 | 110 | 35,332 | 10 ms | 196 ms |
| 200 | 5.5 | 3,441 | 126 ms | 194 ms |
| 200 | 110 | 64,522 | 721 ms | 4.4 s |

