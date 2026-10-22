No local copy of Syed et al., and the tests will assert their round-trip formulas, so I'll verify them against the paper itself:

The PDF came down as a file; reading the relevant pages directly:

Theorem 1 holds: DEO gives E[T] = 2(N+1)(1+E) and SEO gives 2(N+1)(N+E), yielding round-trip rates 1/(2+2E) and 1/(2N+2E) respectively. Next I'll work through the schedule optimization and guidance on choosing N.

I've confirmed the DEO/SEO rates, barrier bound, schedule inversion method, chain count, and tuning scheme against the paper. Now I'll check how Algorithm 1 accumulates rejection statistics.

I've confirmed Algorithm 1 tracks a Rao–Blackwellized rejection estimate across all pairs each scan, even though it only swaps the scheduled ones—something our cached evaluations make costless. Now I'll look at how BARS moves chains through burn-in to find where ladder adaptation plugs in.

Design finalized: ParallelTempering will support DEO/SEO/sequential swap schedules, per-pair rejection stats, round-trip counting, spline-based ladder adaptation, and BARS-style burn-in adaptation. Starting now on the chain's retempering hook.

Now `ParallelTempering` with the schedules, per-pair statistics, round trips and ladder adaptation:

Engine builds. Wiring the schedule and burn-in adaptation into BARS:

BARS is wired. Now the samplers tests: Theorem 1's round-trip rates under exact local exploration, the rejection-free limits, rung scaling DEO vs SEO, ladder equalization, and retempering:

All 56 tests pass and the theoretical predictions line up closely across exact local exploration, rejection-free limits, rung scaling, and adaptation checks. Next I'll verify the swap-tie test still catches R16 now that DEO is the default.

Still caught. Now BARS: the probe with `Sequential` and no adaptation should reproduce the pre-3.3 tempered runs bit for bit, which checks the engine refactor itself:

Bit-identical: the Sequential schedule reproduces the pre-3.3 ladders exactly. Running the regression suite with the new defaults (DEO plus burn-in adaptation):

192/192 with the new defaults; both tempered BARS checks find their peaks. Committing the engine first, with BARS pinned to `Sequential` so that commit is zero-drift (verified above):

Engine committed (`e7efd67`). Next the benchmark the brief asks for, which also serves as Q5's revisit check: round trips per unit work on BARS fits. Checking how existing benchmarks are gated:

I'll commit this as a reproducible benchmark gated by a `TEMPER_BENCH` variable, following the BGP benchmark's pattern. Checking BARS's move palette so the benchmark matches it:

The benchmark is running in the background: 3 fits × 5 configurations × 2 ladder sizes × 3 seeds. Meanwhile, reading the doc anchors:

While the benchmark runs, I'm drafting the result-independent parts of the docs:

I'm waiting on the benchmark before writing up BARS's defaults; I'll pick up when it reports.

The benchmark overturns the brief's premise: the sequential sweep outperforms DEO on round trips per scan, and Q5's revisit trigger formally fires in two of three fits, though the powered ladder's "hot end" for count data turns out to be barely distinguishable from the posterior itself. Before concluding anything, I'll compute cold-chain ESS per unit work and verify the sweep's rejection-free rate is actually 1.

Local exploration confirms the sweep's rejection-free rate is exactly 1 per scan and outperforms DEO by ~1.7× on equalized ladders, with its only drawback being sequential swap ordering (irrelevant here since replicas run on one machine). Next, I'll extend the benchmark to include cold-chain ESS per unit work and adaptive sweeps.

The extended benchmark is running (about 8–9 minutes). When it reports I'll set out the schedule and path findings, so you can decide the defaults on cold-chain ESS rather than round trips alone.

The benchmark is done. Here is what holds up and what doesn't.

**Round trips: robust, about ±5% across seeds.** Swap rate is shown alongside to show why it misleads.

| Fit, 8 rungs, likelihood path | DEO, adapted | Sweep, adapted |
|---|---|---|
| sine | 0.047 | 0.086 |
| peak | 0.018 | 0.025 |
| counts | 0.081 | 0.148 |

- **The sweep wins everywhere**, by about 1.4–1.8×, at the same cost per scan (the replica steps dominate; swaps are free here). It also wins under exact local exploration, where theory is clean: about 1.7×, with a rejection-free rate of exactly 1 per scan against DEO's ½.
- **DEO's guarantee is narrower than the brief assumed.** Syed et al.'s optimality is over even-odd schemes, which suit replicas on separate machines. This engine steps a ladder's replicas on one thread, where the sweep's in-order swaps cost nothing.
- **Adaptation** helps clearly on sine (with 4 rungs, the sweep goes 0.045 → 0.073 and DEO 0.027 → 0.036) and is neutral on the other two fits. It never hurt.

**Cold-chain ESS per second: too noisy at this budget to rank configurations.** With 3 chains per configuration, single estimates swing 2–8×. Averaged like for like (the same schedule and adaptation settings on both paths), powered against likelihood:

| | knot-count ESS/s | energy ESS/s |
|---|---|---|
| sine | 1.5× | 0.97× |
| peak | 1.4× | 1.3× |
| counts | 0.85× | 0.84× |

**Q5's trigger.** It fires formally on round trips: peak about 2× and counts about 5× for the powered path. But those trips are to a closer destination. On counts the powered barrier is 0.24 against 1.4, and its hot end sits at 1.0 knots against the prior's 4.4. On ESS the powered path is consistently ahead only on peak, by about 1.3×. That isn't "clearly ahead per unit work". My recommendation is to keep Q5, record peak as the fit to watch, and re-measure once step 3.4's rank-normalized ESS lands, with more chains.

**My recommendation for the schedule:** make the sequential sweep the default, for both `ParallelTempering` and BARS, with ladder adaptation on in BARS. Keep DEO and SEO as options: DEO for future setups where replicas run in parallel, SEO as the theory baseline. This reverses the brief's DEO default on the step's own measure. Nothing BARS-side is committed yet; the engine commit `e7efd67` still defaults to DEO and would change with this.