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

## Added after registration

Written on 2026-11-01 after the base kernels and the case and level ladders had been read, and before any rung of the grid's replicate ladder, whose first run stopped with the session. The registered readings above stand and are reported as written.

- *D is read from the cut's lighter side.* L(A) = L(Aᶜ), but D as written depends on which side is called A, and the first rows named the heavy side, whose depth says how far the top mode climbs and not what traps the chain. So a set is oriented to its lighter side, D is that side's, `bn_mode` is that side's mode, and the other side's D stays beside it as `*_D_other`. b_pred is read at the lighter side's mode and at the top mode.
- *Each local mode's depth, by its free slots.* After the brief was written, the walk of Ronen et al. (SC-4, SC-5) gave a reason for the rival named above: S → S ∪ T → T costs only the Occam term, unless the count cap closes the superset, and the default cap rarely binds. An enumerable fixture's cap closes it for every state at M_max, so a gap could fall exponentially for a reason the walk already excepts. To separate them: the merge tree's depth of each mode, ln π(cluster) − ln of the conductance joining a likelier one, deepest per M_max − M at the mode, 0, 1 and 2 free, each with the same scaling reading. The walk's expectation is met if depth with a free slot grows like ln N, and falsified if it grows linearly in N at a mode where k = 1 is open; depth at 0 free is the cap's. M ≤ 4 against M ≤ 3 reads a mode the wider cap opens.

- *The flat histogram's refusals are readings, and a second route is added.* The registered pilot-started route refused several rungs at once ("SAMC's limit does not exist", where the discrete ℓ₁ spectrum leaves a segment nearly empty). A refusal is logged as a row. So that those rungs still read something, a *declared* route runs beside it on every rung, its window from F's exact density of states (the prior's 1% quantile of ℓ₁ to the posterior's 99%) at K = 12. It uses the enumeration, so it says what the sampler does with a good window, not what a run would find unaided. Added after seeing the first refusals and rows, on the base fixtures and r = 2.

- *M ≤ 5, depths only.* Added after M ≤ 4 inside F came back with depth at one free slot reading linear, which is the falsifier named above unless two terms must be exchanged at once. M ≤ 5 gives that mode two free slots. Its F is out of reach of the elimination, but a mode's depth needs only the kernel's edges, so `Depths_AcrossCaps` reads each mode at caps 3, 4 and 5 on the same data, and follows the modes the narrower caps trap. *Read before the run:* if the followed mode's depth stops growing at two free slots, the linear growth at one free was the cap closing a two-term exchange, which the walk excepts, and the walk stands. If it keeps growing linearly at two free, there is a separator the walk says should not exist, and the walk is falsified on this fixture. The first run, to N = 1,152, left the followed mode with a single increment at two free slots, so the same reading was extended to r = 512 and 2,048 (N = 4,608 and 18,432) before anything was written from it. Those are the only rungs past N = 2,304, and only depths are read there.

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

### 2026-11-01: measured

**In one paragraph.** The gap falls exponentially in N on every replicate cell, so the register's first expectation holds as a statement about γ on these fixtures. But all of that growth sits at modes at the count cap. Modes with enough free slots stop deepening, which is what the walk (SC-4, SC-5) said would happen unless the cap closes the superset. A mode needing two terms exchanged at once deepens linearly at zero and at one free slot and nearly stops at two. So both readings stand in their own domains, and these fixtures' γ is mostly measuring the cap. The move mix matters more than expected: the fixed-M family raises γ 33- to 210-fold on the base fixtures. An R̂ failure is trustworthy here (3 of 364 on chains the exact kernel calls mixed), while a pass misses 9% of unmixed runs. The flat histogram with a good window reads the exact posterior where K's gap is about 10⁻¹³⁰, but its pilot-started route refuses most rungs, and on one fixture was wrong by more than its error bar.

Rows: [`kernels`](spectral-gap-kernels.csv), [`ladder`](spectral-gap-ladder.csv), [`scaling`](spectral-gap-scaling.csv), [`capdepths`](spectral-gap-capdepths.csv), [`rhat`](spectral-gap-rhat.csv), [`flat`](spectral-gap-flat.csv). Harness: `tests/regression/BassKernelSpectrum.cs` and `BassKernelSpectrumTests.cs`, behind `BASS_QUALIFICATION` = 1, with `BASS_GAP_WIDE` = 1 for M ≤ 4 and 5.

#### What the runs show

**The checks.** All eight held on every kernel read. Check 3 ran on the continuous default mix, functional's five mixes and move classes up to 400 states; the categorical fixtures (2,918 states) and the continuous options were past the dense solver's budget, about 19 minutes at 2,263 states. Check 8: M's autocovariance at lags 1, 2, 5, 20 reads 0.1626 ± 0.0019, 0.1365 ± 0.0018, 0.0815 ± 0.0020, 0.0099 ± 0.0020 against K^k's 0.1607, 0.1345, 0.0801, 0.0095, so one sweep at held τ is one step of K. A rerun reproduced the ladders' 141 rows exactly in the ten columns compared, and the flat-histogram rows it repeated. Where the dense solver ran, λ_min stayed positive (0.14–0.27), so the absolute gap is 1 − λ₂.

**The gaps,** base fixtures, τ at unit information:

| Fixture | \|F\| | default | cancelling | steps | family | all |
|---|---|---|---|---|---|---|
| continuous, N = 18 | 2,263 | 3.88e-4 | 3.35e-4 | 5.24e-4 | **2.74e-2** | 1.72e-2 |
| categorical, R's measure, N = 24 | 2,918 | 9.20e-5 | 6.29e-5 | 8.66e-5 | **3.11e-3** | 3.40e-3 |
| categorical, uniform, N = 24 | 2,918 | 9.09e-5 | 6.18e-5 | 8.47e-5 | **3.01e-3** | 3.25e-3 |
| functional, N = 24 | 593 | 5.00e-5 | 3.22e-5 | 4.98e-5 | **1.05e-2** | 1.08e-2 |

Per move kernel, inside its mix:

| Group | continuous | categorical (R's) | functional |
|---|---|---|---|
| birth/death, irreducible | γ 1.01e-5 | 7.39e-6 | 6.78e-6 |
| redraw: classes; worst; heaviest | 18; 1.23e-6; 0.111 | 20; 1.42e-4; 2.78e-2 | 20; 1.59e-4; 8.34e-2 |
| local step: classes; worst | 157; 1.74e-5 | 75; 3.37e-7 | 117; 2.37e-2 |
| prune/grow: classes; worst | 4; 1.30e-10 | 3; 5.40e-7 | 3; 1.64e-7 |
| variable change: classes; worst | 852; 6.94e-11 | 50; 8.89e-8 | 50; 9.05e-5 |

- The redraw's classes are exactly the key basins on all three fixtures, as expected.
- *Within or between,* by the registered rule: continuous and functional read *mixed*, categorical *within*. The rule compares with the worst class, which is a light one. Read after the fact against the heaviest class, the redraw's gap is 290, 300 and 1,670 times the default mix's γ.

**The bottleneck,** default mix, from the lighter side:

| Fixture | Set, family | π(A) | L·γ | D | Mode |
|---|---|---|---|---|---|
| continuous | 75 states, eigenvector cut | 0.16% | 0.42 | 8.8 | x0+0·x1+0 + x0+0·x1−2 |
| categorical, R's | 4 states, eigenvector cut | 13.2% | 0.77 | 11.0 | x0+1·z0{0} + z1{0,1} |
| categorical, uniform | 4 states, merge tree | 13.2% | 0.77 | 11.0 | the same |
| functional | 585 states, merge tree | 28% | 0.43 | 10.8 | x0+0·t1+0 + x0−2 |

On the continuous fixture the gap is set by a basin with 0.16% of the mass; the best relevant set (3.7%) gives only L·γ = 0.002. SC-2 was within a factor of 2.4 of t_rel for every bottleneck above, and within about 5 on the ladders (L·γ from 0.19 to 0.87).

**The scaling, ladder R,** F held. −ln γ against N:

| Cell | b per seed | −ln γ | Deepest at M_max, top rung | Deepest with a free slot, N ≥ 144 |
|---|---|---|---|---|
| grid M ≤ 2, inside F | 0.714–0.723 | linear | 1,666 (N = 2,304) | ≤ 7.6 |
| grid M ≤ 2, outside F | 0.319–0.337 | linear | 781 | ≤ 7.6 |
| grid M ≤ 3, inside F | 0.110–0.126 | linear | 293 | ≤ 7.6 |
| grid M ≤ 3, outside F | 0.248–0.256 | linear | 595 | ≤ 7.6 |
| grid M ≤ 4, inside F (N to 1,152) | 0.090 | linear | 109 | **61.2 at 1 free**, 6.3 at 2 |
| grid M ≤ 4, outside F (N to 1,152) | 0.253 | linear | 294 | ≤ 7.5 |
| categorical, either measure | 0.552–0.553 | linear | 831 (N = 1,536) | ≤ 7.6 |
| functional | 0.526 | linear | 835 | ≤ 7.6 |

- *The register's expectation:* −ln γ reads linear in 17 of 17 cells, and D in 17 of 17. None reads logarithmic or flat, so the rival as the brief worded it fails for γ. Rungs below the solver's range are left out of the fits and keep their depths.
- *Depth by free slots:* at 0 free, linear in 17 of 17. At 1 free, flat 4, neither 12, linear 1. At 2 free, flat 8, neither 9. "Neither" here is increments of mixed sign around zero, with fitted b between −0.003 and 0.001. Over all cells to M ≤ 3 no free-slot mode goes past 13.8 nats at any N.
- *Where the trap sits:* from N = 144 up, the bottleneck's mode is at M_max on every rung of every cell. Below that it can have a slot free: on the continuous grid at N = 18–72 it is an M = 2 mode, whose depth peaks at 11.4 and is gone by N = 144.
- *Identity:* inside F at M ≤ 3 the trap is x0+1·x1−1 + x0−2·x1+0 + x0−2·x1−2 from N = 144 to 1,152, then another. The trap's π(A) falls with N, to 10⁻⁷⁸ outside F, so the bottleneck stops being relevant: a set holding 1% each side exists up to about N = 144–576 outside F. Where one exists, ln L(relevant) reads linear in all five outside-F cells with rungs enough. Inside F and on the case fixtures it reads flat, neither, or too few.

*The followed mode.* The registered falsifier came up at M ≤ 4 inside F, so M ≤ 5 was added to read it, with the rule written first. The same state, x0+1·x1−1 + x0−2·x1+0 + x0−2·x1−2, on the same data:

| N | 0 free (cap 3) | 1 free (cap 4) | 2 free (cap 5) |
|---|---|---|---|
| 72 | 8.5 | 8.4 | 8.4 |
| 288 | 43.9 | 24.7 | 17.8 |
| 1,152 | 154.7 | 61.2 | 19.1 |
| 4,608 | 567.1 | 199.9 | 20.6 |
| 18,432 | 2,349.5 | 801.9 | 21.2 |
| slope per N, last pair | 0.129 | 0.044 | 0.00004 |

At two free slots the increments per quadrupling are 1.27, 1.52 and 0.58, around SC-5's k = 2 Occam term, ln 4 = 1.39. By the rule, the walk stands on this fixture. Outside F one slot is enough: its cap-3 trap reads 299.6 at N = 1,152 and 4.5 at cap 4. Each cap makes its own traps: the deepest at 0 free at N = 1,152 reads 154.7, 108.6, 62.6 inside F and 299.6, 293.9, 272.0 outside at caps 3, 4, 5.

*SC-4's slope,* b_pred at the trap against dD/dN over the last doubling: M ≤ 2 inside F 0.73–0.76 against 0.69–0.74; M ≤ 3 inside 0.10–0.19 against 0.12–0.14; M ≤ 3 outside 0.26–0.32 against 0.24–0.26; M ≤ 2 outside 0.12 against 0.34; categorical 0.37 against 0.54; functional 0.31 against 0.56. The right order of size everywhere: close in the first three, and about a third to two thirds of the slope in the last three.

**Ladder L,** F growing from 299 to 5,049 as N goes 18 to 72, M ≤ 2, four rungs. −ln γ goes from 10–13 to 33–43 inside F and 6–9 to 18–23 outside; the linear fit is lower in 5 of 6 cells. Depth at 1 free goes from 3–5 to 6–8, the logarithmic fit lower in 5 of 6, with c from 1.5 to 2.6, inside SC-4's range for k = 1 at K = 2 (½ to 2½). Four rungs say little.

**The exact law from the prior.** TV(μ_t, π) plateaus: at N ≥ 144 it is nearly the same at 7,000 and 70,000 sweeps, with medians over seeds at 7,000 of 0.19–0.45 (M ≤ 3 inside F), 0.94–0.96 (M ≤ 3 outside), 0.37–0.82 (M ≤ 2 inside), 0.83–0.91 (M ≤ 2 outside), and 0.84 on categorical and functional at N = 1,536. That is the prior's mass that falls into traps and stays.

**R̂,** 1,064 fits at held τ. Predictive gate against truth, fail / pass:

| Lengths | Not mixed | Mixed | Between |
|---|---|---|---|
| 500 + 200 | 361 / 19 | 0 / 12 | 45 / 31 |
| 5,000 + 2,000 | 334 / 34 | 3 / 21 | 27 / 49 |
| 50,000 + 20,000, seed 5 | 80 / 4 | 4 / 8 | 11 / 21 |

- *False passes* at S10's lengths: 34 of 368, 9%. Of these 15 have TV above 0.25, 28 have T < L(relevant), and 9 have both. One is M ≤ 3 outside F at N = 144, where TV is 0.93 and T·γ about 10⁻¹⁵, with predictive R̂ 1.001: four chains agreeing on the wrong surface.
- *False alarms:* all seven are M ≤ 3 outside F at N = 18, by σ², with R̂ from 1.06 to 1.22, where T/L is 10–12, at the rule's edge, and T·γ is 5.8–6.8. None away from that edge.
- *Other groups,* S10's lengths, not mixed, fail / pass: level 344 / 24; structural (M) **51 / 317**. M's R̂ passes 86% of unmixed runs, since trapped chains share M = M_max.
- *ESS:* on the four rungs at S10's lengths with mixed truth and τ_int(M) readable, the Count diagnostic's bulk ESS is 1.7, 1.8, 2.1 and 4.8 times the exact 4T/τ_int(M). τ_int is blank (44 of the ladders' 141 rungs) where Var M is under 10⁻⁸, so the `exact_count_ess` column is blank there too.

**The flat histogram,** 40 runs, against the enumeration.

| Route | Runs | Refused | ln Z within the rule | Basins within the rule |
|---|---|---|---|---|
| pilot-started, as registered | 20 | 13, every rung at r ≥ 8 and one at r = 2 | 6 of 7 | 5 of 7 |
| window declared from the exact density of states | 20 | 0 | 19 of 20 | 20 of 20 |

- The declared route holds basin TV ≤ 0.021 with the top basin right at every rung, including r = 128 outside F, where γ is 4.9e-131 and sweeps·γ is 3e-125.
- The pilot route misses on functional, both seeds: the top basin is right, but its mass reads 0.58 and 0.89 against 0.718 (basin TV 0.16 and 0.19), with 7 and 11 round trips across the chains. One of the two also reads −25.768 ± 0.031 against −25.530, 7.6 errors out.
- The error bar is too small on 6 of 27 runs (|z| > 3, worst 7.6). The rule's 0.1-nat floor passes most of them. The declared route's one miss is −69.507 ± 0.018 against −69.373.
- One pilot run crashed rather than refused: `OverflowException` from `Wham.Nodes` (`src/maths/samplers/flathistogram/Wham.cs:251`).

#### What is made of it

- **The two expectations.** The first draft's mechanism, load-bearing deaths, is right for a state with no free slot, and the walk's is right for one with enough. Which applies depends on k and the slots free, and M_max is what the enumerable fixtures can reach, so γ here reads the cap. The walk's own exception covers this, but it should not be read as confirmation for an uncapped fit: see the limits.
- **Polynomial is not fast.** The followed mode's k = 2 depth is about 20 nats at any N past 288, which is about e⁻²⁰ per escape, before anything larger.
- **A cap near the truth's size is a trap.** At M_max close to what a surface needs, most of the prior's mass ends up stuck (TV up to 0.96 outside F) with no slot to escape by. Real fits' default cap rarely binds, so this is a warning for anyone setting a small M_max, not a reading of Friedman.
- **Exchanging terms is what helps.** The family's 33–210-fold gain on γ, with its own kernels' gaps as low as 10⁻¹⁰, says the value is in the paths it opens, not its own speed. Together with k = 2's reading, this supports the walk's suggestion of a several-term exchange. That is a reading of small fixtures, and nothing here tests such a move. Whether the gain survives at Friedman's scale is not read here.
- **R̂.** Within these fixtures' limits, F-R1's kind of failure deserves belief: a failed gate here is almost always a chain the exact kernel does not call mixed. A pass is weaker evidence, and four chains can agree on one trap. M's R̂ agreeing beside a failing predictive R̂ can mean chains trapped at a shared M, which the diagnostic's note on non-identifiability does not mention.
- **The flat histogram.** With a good window it is not bound by K's gap, which lends support to the method behind Friedman's a1. It says nothing about a1 itself. The functional runs caution that few crossings can give a confident wrong answer. Friedman's runs made 47–62 round trips, 6–21 per chain, against 7–11 in all there; those are different fixtures, so this is a question to keep open, not a verdict. The pilot's refusals may be largely a small-F effect, since the discrete ℓ₁ spectrum leaves segments empty and Friedman's is dense.

#### Limits

- These fixtures are small and F grows fast with N. The readings show scaling, not Friedman's barrier itself. Each cell's truth has few terms, and M_max sits close to it.
- **Ladder R nearly forces the free-slot result.** On a fixed grid, any term with real partial correlation gains in proportion to N, so a mode with a slot free and a useful single term cannot stay a mode for long. Friedman's weak directions, where the gain is comparable to ½ ln N, are the regime this ladder leaves quickly. Its low rungs, and ladder L, are the nearest reading, and they are small. So this supports the walk without testing it where it matters most.
- Replicates differ from fresh covariates, whose knots would move with N. Ladder L has only four rungs.
- k = 2 rests on one mode on one data seed, on one surface. k ≥ 3 is unread.
- τ is held at unit information, which grows with N, and Var M is too small to read τ_int past some rungs. A sampled τ, and a sweep's block scan with it, lie outside this kernel, so R̂ here is at held τ, not at F-R1's sampled τ.
- The truth classes for R̂ are the brief's thresholds, and "relevant" means 1%. The registered within-or-between rule read a light class and said little.
- Reading the sets is not an exhaustive search. L(A) is a lower bound over the sets tried, and L·γ shows how close it came.

#### Changes outside the new files

`BassKernelEnumeration.Space.Key` packs ids at 12 bits instead of 20, so five terms and their count fit a `long`, with a guard on the term count. It is only an identity key; the enumeration's existing tests pass unchanged.

#### Proposed, not built

- **The shared mixing layer.** `ReversibleChain`, `EliminatedChain` and the mode depths would serve any engine with an enumerable kernel, as an exact reference for its diagnostics. Propose moving them once a second consumer exists; the mean-field Potts toy (the register's G1, item 3) is the natural one. Until then they stay in tests.
- **`Wham.Nodes`.** A b = 1 rung sitting on one level has a rounding-size sd, which passes the `sd > 0` guard and overflows `checked((int)…)`. It should give the existing "no spread" refusal.
- **For the main session's record:** the register's G1 and SC-4 can take these readings, and F-M2's record the family's gains. Nothing here changes a default.
