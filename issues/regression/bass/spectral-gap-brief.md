# The structure chain's true spectral gap on enumerated fixtures (chip brief)

Date: 2026-11-01. A measurement, tests only: it changes no code under `src/` and no default. Gap G1 of the [lineage-gaps register](lineage-gaps-register.md), step 1; SC-0 and SC-2 of [`structure-chain-lemmas.md`](../../../lean/Protolemmata/BASS/structure-chain-lemmas.md).

- *Why.* All 144 of F-R1's untempered fits fail the predictive R̂ gate, and on Friedman the flat histogram finds a basin at 66–68% that MCMC and the references hold at no more than 1.3%. BART has mixing theory; none was found for Bayesian MARS. At fixed τ, with knots at distinct data values and the coefficients and σ² collapsed, the structure chain is a finite reversible chain, so its relaxation time can be computed outright on small fixtures. That is ground truth for the diagnostics, and for any later lower-bound proof.
- *The question.* What is 1 − λ₂ of each move kernel and of the default mix; which set is the bottleneck; how does the gap scale with N; and do R̂ and the flat histogram report what the exact kernel says?

This section and the next three were written before any gap, flow or barrier was computed.

## What is measured

- *The kernel.* K(S, S′) and π(S | τ, y) from `BassKernelEnumeration`, τ held at unit information, b = 1. K(S, S) = 1 − Σ K(S, S′), so refusals, null moves and rejections hold.
- *The gap.* γ = 1 − λ₂ of D^(1/2) K D^(−1/2); t_rel = 1/γ. λ_min beside it where the dense solver runs, since the absolute gap needs it.
- *Mixes.* `default` is what a run does: birth/death and redraw on equal-over-available, no options. `cancelling` is the same on the prior-cancelling schedule at C = 0.4. `steps` adds local steps (ρ = ½, ω = 0.02), `family` adds the fixed-M factor moves, `all` adds both.
- *Per move kernel.* One group's flows inside its mix's palette, everything else held. A fixed-M group cannot be irreducible, so it reads per communicating class: the class count, the worst class's gap, and the gap of the class with the most mass.
- *The sets (SC-2).* L(A) = π(A)(1 − π(A)) / flow(A) lower-bounds t_rel for any A. Four families:
  - each *key basin*, the Friedman band brief's key (`BassMixing.Fixture.Basin`, the multiset of terms' variable sets);
  - each M level;
  - the *merge tree*: every connected set formed while adding K's edges in descending conductance π(S)K(S, S′), so each local mode's basin and their unions;
  - *eigenvector cuts*: thresholds at the eight largest jumps of the slow eigenvector.

  The *bottleneck* is the set with the largest L(A). The *relevant bottleneck* is the same over sets with π(A) and 1 − π(A) both at least 1%, since the gap can be set by a basin too light to matter to an answer.
- *The barrier.* D(A) = ln π(A) − ln of the largest conductance leaving A: how deep A's way out is, in nats.
- *The solvers.* `DenseEigen.DecomposeSymmetric` on the symmetrized kernel, which fits while γ stays well above its absolute error near 10⁻¹⁴. Below that, and for speed on sweeps: subtraction-free elimination of I − K (Grassmann, Taksar & Heyman), eliminating ascending π, which keeps relative accuracy from K's off-diagonals, then Lanczos on its inverse orthogonal to √π. Rates below about e⁻⁷⁰⁰ underflow; a rung that hits this says t_rel ≥ e⁷⁰⁰, leaves the gap fits, and keeps its barrier, which is computed in logs.

## The fixtures

- *Base.* `BassSweepTests.Data` at M ≤ 3 (the 2,263 states, N = 18); `CategoricalData` under both subset measures and `FunctionalData` at M ≤ 2 (N = 24).
- *Ladder R: replicates, F held.* Each grid point repeated, so the knots and the support rule's outcome stay put and N moves alone. The continuous 3 × 3 grid at r = 2, 4, …, 256 (N = 18 … 2,304), with r = 2 being the base data and each rung extending the last. Two surfaces, since the problem class may decide the answer:
  - *inside F*, the base surface, 1 + 2a + 3(a − ½)₊(½ − b)₊ + 0.3ε;
  - *outside F*, smooth, 2 sin(π·a·b) + 4(b − ½)² + a + 0.3ε, which no M ≤ 3 structure can hold, closer to Friedman's kind.

  Caps M ≤ 2 and M ≤ 3; data seeds 5, 11 and 23. M ≤ 4 on seed 5 at r = 2, 8, 32, 128 if its F stays under about 12,000. The categorical and functional fixtures too at r = 1, 2, 4, …, 64, seed as built.
- *Ladder L: levels, F grows.* Grid side L = 3, 4, 5, 6 at two replicates, M ≤ 2 (N = 18, 32, 50, 72), both surfaces, the same seeds. Real data's knots grow with N; here F grows with them, so only four rungs are enumerable.

## How the numbers will be read

**The gaps.** A table per fixture and mix: |F|, transitions, classes, γ, t_rel, λ_min. *Within or between:* the default mix's slow mode is *between basins* if γ is under a tenth of the worst class gap of the redraw alone (whose classes are expected to be the key basins), *within* if it matches that gap to a factor of 2, otherwise *mixed*.

**The bottleneck.** The set, its mode's terms, π(A), flow, L(A), D(A), L(A)·γ (1 means SC-2 is tight), for the bottleneck and the relevant bottleneck, and which family found each.

**The scaling.** For each ladder cell (fixture, surface, cap, seed), y(N) = −ln γ(N), and the same for D of the bottleneck and for ln L of the relevant bottleneck.
- Fit y = a + bN and y = a + c·ln N, and read the doubling increments Δ_k = y(2N) − y(N) on R, where N doubles.
- *Linear:* over the last three increments each is at least 1.5 times the one before. *Logarithmic:* the last three stay within a factor of 1.25 of each other. *Flat:* the last three are all below 0.1. Otherwise *neither*, and said so. L's four rungs do not double, so L reads only the two fits' residuals and says less.
- *The register's expectation* (barriers about linear in N, SC-4) holds in a cell that reads linear and fails in one that reads logarithmic or flat. The count of cells by surface and cap is reported, contradictions included, with no cell dropped.
- *Its mechanism* is SC-4's: leaving a basin loses a load-bearing term. Predicted slope, at the top rung: b_pred = the smallest over the mode's terms of [ln m(mode) − ln m(mode without it)] / N. b is read against b_pred. A rival, that a birth goes first and only an Occam cost near ½ ln N per spare term is paid, predicts c near ½ per term. Which caps show which is what M ≤ 2, 3 and 4 are for. No outcome is predicted here beyond the register's.
- *Identity.* Whether the bottleneck is the same set up the ladder, and whether its π(A) is negligible.

**R̂.** `Bass.Fit` at held τ on the same fixtures, four chains, master seeds 7101–7104, at S10's lengths (5,000 + 2,000) and a tenth of them, and at ten times them on seed 5 only. The gate is F-R1's, predictive R̂ ≤ 1.05, with the level and structural groups beside it. What it is held against is exact and shares nothing with it:
- *The law at time t.* μ_t = μ₀K^t from the run's own start, the prior on F. TV(μ_burn, π) at burn-in's end and at the last sweep, and the bottleneck's mass under μ_burn.
- *The relevant clock.* T·γ and T / L(relevant bottleneck), T the kept sweeps.
- *Exact ESS.* 4T/τ_int(M), with τ_int from the Poisson equation, against the Count diagnostic's ESS.
- *Truth, per run:* **mixed** if TV(μ_burn, π) < 0.05 and T ≥ 10·L(relevant); **not mixed** if TV(μ_burn, π) > 0.25 or T < L(relevant); otherwise *between*.
- *Tally:* gate against truth. A *false pass* passes while not mixed; a *false alarm* fails while mixed. Each is listed by cell, per group.

**The flat histogram.** `Bass.Fit` with its pilot-started route (four chains; pilot burn 5,000, traced 2,000, the linear form's ladder with its reference rung; learning 30,000, production 30,000, `Thin` 10; master seeds 7201, 7202), on ladder R at r = 2, 8, 32, 128, both surfaces, M ≤ 3, seed 5, and on the categorical and functional base fixtures.
- ln Ẑ(1) against the enumeration's ln Z(1): *agrees* if within the larger of 3 of its errors and 0.1 nats.
- Its weighted key-basin masses against π: *agrees* if TV over key basins is at most 0.05 and the top basin is the same.
- Round trips beside them. A rung where it agrees while T·γ ≪ 1 for the chains' sweeps is one where it is not bound by K's gap, and that is said as a reading of that rung only.

**Kept apart.** What the runs show goes in one list, what is made of it in another.

## Checks

An assertion fails the part.
1. `BassKernelEnumeration.Check` on every kernel read: detailed balance and πK = π.
2. λ₁ = 1 within 10⁻¹⁰ and its eigenvector ∝ √π within 10⁻⁸ wherever the dense solver runs.
3. The two solvers agree to 10⁻⁶ relative wherever both run and γ ≥ 10⁻⁸.
4. SC-2: L(A) ≤ t_rel(1 + 10⁻⁶) for every set tried, every mix.
5. Dirichlet forms add, so γ_mix ≥ Σ γ_group within 10⁻⁸ relative.
6. Solver facts on chains with closed forms, always run: a two-state chain, a birth–death chain, and a two-well chain with γ near 10⁻⁴⁰.
7. On ladder R, |F| is the same on every rung of a cell.
8. One sweep at held τ is one step of K: on the base fixture the production chain's M autocovariance at lags 1, 2, 5, 20 matches K^k's within five batch errors. T·γ means nothing otherwise.

## Limits, stated in the report

- These fixtures are small and F grows fast with N. The readings show scaling, not Friedman's barrier itself. Ladder R's replicates differ from fresh covariates, whose knots would move with N.
- τ is held. A sampled τ adds a continuous coordinate outside this kernel. A sweep's block scan with τ, a composition, is outside it too, so R̂ here is read at held τ, not at F-R1's sampled τ.
- Unit-information τ grows with N on the ladders, which is the engine's own default and part of what is read.

## Working rules

- Work on `main`. Commit by explicit path (`git commit -- paths`), leaving nothing staged: the index is shared. Push nothing. Commits end with the `Co-Authored-By` line the task gave.
- New files only: `tests/regression/BassKernelSpectrum.cs` and `BassKernelSpectrumTests.cs`, heavy parts behind `BASS_QUALIFICATION` = 1, with an accessor elsewhere only if one is needed. Other sessions edit under `tests/regression/`; `lean/` and the protolemmata are `lean-leads`'s.
- Build and run into `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-spectral-gap\`, set on each `dotnet` command itself. Delete that folder when done. No subagents, no worktrees. Other sessions' test hosts are running, so heavy parts keep to a few threads.
- Leave the user's uncommitted notes under `issues/regression/` and the uncommitted `src/graphs` work alone.
- If gap computation earns a place in the shared mixing layer, say so and propose it; do not build it here.
- Ends with a fresh-eyes review, fixed before the report.

## Reports

Newest first.
