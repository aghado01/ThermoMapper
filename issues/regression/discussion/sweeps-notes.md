### The Sequential Sweep "Win" Is a Single-Scan Cascade Artifact

The "Sequential Sweep" win on round trips is due to a fundamental mechanics flaw in how an in-order sequential loop interacts with local MCMC steps.

---

### 1. Why the Sequential Sweep’s $1.0$ Round-Trip Rate Is Fake Mixing

The report claims that Syed et al. (2022) only restricted their optimality theorem to even-odd schemes (`DEO`/`SEO`) because even-odd pairs can run on separate machines, whereas on a single thread "the sweep's in-order swaps cost nothing". **That is false.**

Syed et al. restrict Parallel Tempering schedules to **vertex-disjoint matchings** (like Even and Odd passes) because **replicas only take a local MCMC mutation step (`Step()`) once per scan, *between* swap phases**:

* **Under `DEO` (vertex-disjoint swaps):**
Every replica participates in **at most one swap per scan**. Therefore, a particle traveling from rung $0$ to rung $N$ and back ($0 \to N \to 0$) must take at least $2N$ scans and is **guaranteed to undergo at least one local MCMC mutation step (`Step()`) at every single intermediate rung $1, 2, \dots, N-1$**.
* **Under a Unidirectional `Sequential` Sweep (`for i = 0..N-1: Swap(i, i+1)`):**
The swap pairs $(0,1), (1,2), \dots, (N-1,N)$ overlap inside the same scan. Look at what happens when `Swap(0, 1)` succeeds: the state $x_0$ that started at rung $0$ is now sitting in slot $1$. Immediately—**in the very next iteration of the `for` loop, with ZERO local MCMC steps at rung $1$**—the loop executes `Swap(1, 2)` on **that exact same unmutated state $x_0$**!

Now trace the $r = 0$ (rejection-free, 100% swap acceptance) limit for $N = 3$ (4 rungs: $0, 1, 2, 3$) starting with particles $[p_0, p_1, p_2, p_3]$ to see why `Sequential` gets a rate of $1$ instead of `DEO`’s $1/2$:

1. **Inside a single `Sequential` scan (`i = 0, 1, 2`):**
* `i = 0`: `Swap(0, 1)` $\to [p_1, p_0, p_2, p_3]$
* `i = 1`: `Swap(1, 2)` $\to [p_1, p_2, p_0, p_3]$
* `i = 2`: `Swap(2, 3)` $\to [p_1, p_2, p_3, p_0]$


2. **Look at what particle $p_0$ just did:**
* In **one single scan**, $p_0$ teleported all the way from rung $0$ to rung $3$ while taking **ZERO local MCMC steps at rung $1$ and ZERO local MCMC steps at rung $2$**!
* Then, on the return leg ($3 \to 2 \to 1 \to 0$), $p_0$ moves against the loop direction at $1$ rung per scan ($N = 3$ scans).
* Total round-trip time per particle: $1 + N = N + 1$ scans (giving an ensemble round-trip rate across all $N+1$ particles of $\frac{N+1}{N+1} = 1$ per scan, vs. `DEO`’s $\frac{N+1}{2(N+1)} = \frac{1}{2}$ per scan).


* Worse still: at the turnaround endpoint (rung $N$), $p_0$ gets **only a single MCMC step** before the next scan’s `Swap(N-1, N)` immediately pulls it back down to $N-1$!



**Why this invalidates the `Sequential` round-trip comparison:**

* Even away from $r = 0$, as you add rungs $N$ for a fixed barrier $\Lambda$ (so $r_{i,i+1} = \Lambda / N$), the probability that a particle cascades multiple rungs in the loop direction during a single scan—**skipping local MCMC exploration completely at every intermediate rung it jumps over**—approaches $e^{-\Lambda}$ for a full-ladder teleport! (For $\Lambda = 0.24$, $e^{-0.24} \approx 78.6\%$ of upward trips jump the entire ladder in 0 intermediate MCMC steps!)


* In Syed et al. (2022), the **Local Exploration Assumption (LEA)** assumes a particle thermalizes via the local MCMC kernel *between every swap attempt*. When your `Sequential` local-exploration test resamples states only once per scan (outside the swap loop) or counts multi-rung single-scan cascades as round trips, it rewards `Sequential` for bypassing the intermediate temperatures without ever running the MCMC kernel on them.


* That is why you should **keep `DEO` as the default** (exactly as committed in `e7efd67`), with burn-in ladder adaptation enabled.



---

### 2. Why the Q5 Benchmark Actually Proves the Case for the `Likelihood` Path

Look at the numbers in `sweeps-results.md` comparing the `Powered` path against the `Likelihood` path:

* On `counts`, the `Powered` path gets **$5\times$ more round trips** than the `Likelihood` path.


* Yet on actual **cold-chain ESS per second**, the `Powered` path is **worse** (`0.85×` on knot-count ESS/s and `0.84×` on energy ESS/s)!



Why did $5\times$ more round trips produce *lower* cold-chain ESS?

1. **Round-trip rate only measures mixing when the two ladders share the same hot endpoint $\pi_{\text{hot}}$.**
Under Syed et al.’s Theorem 1, the asymptotic `DEO` round-trip rate is $\frac{1}{2 + 2\Lambda}$, where $\Lambda = \sum_{i} r_{i,i+1}$ is the thermodynamic barrier between the cold posterior and the hot distribution.
2. **On `counts`, the `Powered` path's hot rung collapsed to a degenerate distribution.**
Because of the $(k!)^{\beta - 1}$ simplex-Jacobian leak identified in `tempering-questions.md`, the `Powered` hot rung sits at $\bar{k} = 1.0$ knot with a tiny barrier ($\Lambda = 0.24$), whereas the `Likelihood` path bridges all the way to the true $\text{Poisson}$ prior at $\bar{k} = 4.4$ knots ($\Lambda = 1.4$).
Cycling a replica rapidly to a hot rung frozen at $k = 1$ knot and back (`5×` round trips) strips away its knots rather than helping the cold chain explore $k \in \{3, 4, 5, 6\}$—which is why its knot-count ESS/s drops to `0.85×`!


3. Why `Powered` edged ahead (`1.3×`) on `peak`:
On `peak`, the posterior only needs a small number of knots around a single bump, and dropping to $k = 1$ knot on the hot rung lets that single knot slide across the bump without bumping into neighboring ordered knots ($\xi_1 < \xi_2 < \dots < \xi_k$). A `1.3×` difference across 3 noisy chains (where single estimates swing $2\text{–}8\times$) is well inside Monte Carlo noise and nowhere near a reason to abandon measure-theoretic coherence.



---

### Summary of Decisions to Hand Back

1. **Swap Schedule Default:** Keep `DEO` (plus burn-in ladder adaptation) as the default in both `ParallelTempering` and `BARS` (keep `e7efd67` as-is). Do *not* switch to `Sequential`. Document why `Sequential` inflates round-trip counts: overlapping swaps in a single scan cascade states across multiple rungs without intervening MCMC `Step()` mutations, violating Syed et al.'s Local Exploration Assumption and skipping intermediate-temperature exploration.


2. **Q5 Tempering Path Default:** Keep `Likelihood` as BARS's default. Record in the docs that `counts` proved round-trip rate cannot be compared across different hot endpoints ($\Lambda = 0.24$ at $\bar{k}=1.0$ vs. $\Lambda = 1.4$ at $\bar{k}=4.4$, where $5\times$ more round trips yielded $0.85\times$ the ESS/s), and defer any `peak` re-measurement until Step `3.4`’s rank-normalized ESS and multi-chain accumulator land.