# The route waiting times against exact hitting times and first passages (chip brief)

Date: 2026-11-04. A measurement, tests only: it changes no code under `src/` and no default. Gap G1 of the [lineage-gaps register](lineage-gaps-register.md), candidate reading (d). It follows the [path profile](path-profile-brief.md), whose report owes this check.

- *Why.* The path profile reads Σ_e exp(π(S) − ln c_e) sweeps from a start to a1 along one route: 29–45 nats by class flow from the chains' states on S10, against a 2,000-sweep phase. The claim that follows, "consistent with the chains staying put", rests on that series. Nothing has compared it with a real crossing time.
- *The reading has two errors of opposite sign.* A reversible chain's mean hitting time is exactly

  E_S τ_T = R_eff(S, T) · M, with M = Σ_x π(x) v(x),

  where v is the voltage: v(x) is the chance of reaching S before T from x, so v(S) = 1 and v(T) = 0. The route reading W = π(S) · R_route uses one route's series resistance and the start's own π. So

  ln W − ln E_S τ_T = ln(R_route / R_eff) − ln(M / π(S)),

  and both logs are non-negative: the first by Rayleigh monotonicity, the second because v(S) = 1 and M includes S. Parallel routes make the reading too slow. The basin's mass makes it too fast. The path profile's class credit addresses only the first.
- *The mass term is unlikely to be small.* a1's exemplar is 1 of 23,288 basin draws. Between a2 and a1, π differs by 14.8 nats and the basin shares by 5.8, and two-basin balance makes the ratio of the two waits follow the masses. The route reading's ratio follows π, so it is about 9 nats off on these draws.
- *The question.* Does the route reading predict crossing times, and which error wins? Past what the phase length already shows, how long do chains from these states really stay put?
- *Not in scope.*
  - Any move, kernel, schedule, option or default.
  - Choosing a lever, which is the user's, on this report.
  - Routes between modes found at N = 2000 (G1's reading (e)), the fixed-M family and the flat histogram's error.
  - Editing the register or the protolemmas. Put what they should say under "For the main session" in the report.

This section and the next three were written before any wait was computed. The harness's choices, where this leaves them open, go in the report and its docstrings.

**What each test can decide.** Only some of these can confirm the reading, and the report keeps them apart.

| Test | Barrier and size | What it can say |
|---|---|---|
| Part 1: exact, enumerated | any depth, small N | two-sided, on spaces with few parallel routes |
| Part 2: lattice resistance | S10, no chains | part of the parallel-route credit, as a bound |
| Part 3: chains at b = 1, shallow starts | S10, little barrier | two-sided, mostly about proposal cost |
| Part 3: chains at b = 1, deep starts | S10, full barrier | one-sided: refutes the reading if a chain crosses, and cannot confirm it, since the readings sit far past the budget |
| Part 4: a ladder in b | S10, barrier scaled down | two-sided at the real N, on a changed target |
| Part 5: chains at b = 1, N = 125 | carried modes, real barrier | two-sided, lower confidence |

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md), the path profile's brief and report, and [spectral-gap-brief.md](spectral-gap-brief.md), whose report found the small fixtures' depth sitting at the count cap.

In the code:
- *The kernel:* `tests/regression/BassKernelEnumeration.cs` (`Enumerate`, `Kernel`, with its b).
- *The solver:* `tests/regression/BassKernelSpectrum.cs` (`ReversibleChain`, and `EliminatedChain`'s state reduction, which holds digits where a plain LU loses them all once a wait passes about e³⁰).
- *The readings:* `tests/regression/BassPathProfile.cs` and its CSVs, `path-profile-{routes,steps,modes,carried,groundtruth}.csv`.
- *The chains and basins:* `BassMixing.cs` (`Basin`), `BassFlatHistogramQualificationTests.cs` (`ArmConfig`, seed 8017), `BassSweep.cs`.

## What is measured

**Part 1: exact hitting times.** On the path profile's enumerated pairs (127 at M ≤ 3, 40 at M ≤ 4), both directions, on the whole structure kernel at held τ with the production move mix. Per direction:
- ln h, the exact mean hitting time of T from S, by state reduction;
- ln R_eff and ln M, and so the two corrections, ln(R_route / R_eff) and ln(M / π(S));
- ln W by the exact route, the aligned route and class flow, as the path profile reads them;
- the same for hitting T's `Basin` as a set, against the route's series cut at its first state in that basin;
- whether the route holds M at the cap. Every summary is read at-cap and off-cap separately.

**Part 2: the lattice's resistance, no chains.** On each S10 pair's aligned lattice (d ≤ 11), R_lattice by elimination, with slides as the route takes them. By Rayleigh R_eff ≤ R_lattice ≤ R_route, so ln(R_route / R_lattice) is a part of the parallel-route credit at the real N that is measured rather than approximated. The exact lattices (d up to 19) are out of scope.

**Part 3: long chains at b = 1, N = 500.** S10's fixture and data seed, seed 8017's modes and chain states, exemplar 1. The target is a1's `Basin`.
- *Shallow starts:* a2 and a1+x6x8, with aligned readings of 15.3 and 14.5. Eight chains each, run until they cross or reach 10⁸ sweeps.
- *Deep starts:* A, B and chain0–chain3, with class readings of 29–45. Four chains each, in rungs of 10⁶, 10⁷ and 10⁸ sweeps, finishing each rung across all starts before the next.
- *Arm H:* no tempering, no flat histogram, τ held at the reading's τ. This is the chain the reading describes, and it is the primary arm.
- *Arm S:* τ drawn as in production, on the deep starts, to 10⁷ sweeps. This is the chain the sampler runs.
- *Per chain:* the sweep of first entry to the target, and whether it stays for 10⁴ sweeps; the sweep and destination of first exit from the start's basin; the first hit of a1's exact structure; the highest-π structure held and the fraction of sweeps spent on it; and a term-level trace of the last 10⁴ sweeps before first entry.
- *A start already in a1's `Basin`* has no basin crossing. Read it to the exact structure and say so.
- *The mass, read from the chains.* While a chain stays in the start's basin, the fraction of sweeps on the top structure estimates ln(mass / π(top)). So ln W_mass = ln mass + ln R_route is a corrected reading. Report it beside W wherever it can be read, with the flat histogram's estimate where the exemplar has draws, stating how noisy each is. M ≈ the basin's mass only while v ≈ 1 across the basin, and Part 1 checks that on small spaces.

**Part 4: a ladder in b, N = 500.** The deep starts, target a1's `Basin`, arm H, with π_b as the kernel's b defines it.
- At each b, re-search each pair's scored aligned lattice, so no new fits, and read W(b) by the route and by class flow, cut at the target's basin.
- Choose the rungs from those readings alone, before any chain at b runs: for each pair, the b at which the class reading is nearest ln = 10, 12, 14 and 16.
- Eight chains per rung, each stopped at crossing or at 10× the reading. Read the mean wait with its interval, and first exit as in Part 3.

**Part 5: long chains at b = 1, N = 125.** Two starts where a real barrier falls within reach and both ends hold by the carried reading: a1 → A on `stream` (route reading 19.2, depth 4.1), and B → a1 on `fresh` (`partial+redraw`, 20.8, depth 9.6). Arm H, eight chains, to 10⁸ sweeps. These modes are carried rather than found. Say so wherever this part is quoted.

**The budget.** At most 36 core-hours of chains and 4 threads, below normal priority, one test host.
- Measure sweeps per second first and state the plan's cost before starting the chains.
- The order is Parts 1 and 2; Part 3's shallow starts; deep starts to 10⁷, both arms; Part 4; Part 5; then deep starts to 10⁸. If the budget or the machine runs out, stop there and say what was not run.
- Everything appends and skips on rerun. Save each chain's structure every 10⁷ sweeps. A stopped chain resumes from it with a fresh derived seed, which leaves first passage unbiased because the chain is Markov.

## What must hold before any reading counts

- **The solver.** On the enumerated pairs, h by state reduction agrees with R_eff · M from the voltage solve, and with a dense solve wherever h ≤ e²⁰, to six digits.
- **The kernel against itself.** On at least 3 enumerated pairs with h between e⁶ and e¹², simulated passages from the enumerated kernel (at least 1,000 per pair) agree with h within their interval.
- **The production chain against the kernel.** On those pairs, at the same cap, the harness's production chain's first passages agree with h within their interval. A sweep must equal one attempt, as the path profile found. This check is what lets Parts 3–5 speak about the reading, so a failure stops them.
- **Rayleigh.** R_eff ≤ R_lattice ≤ R_route on every row it applies to. A violation stops the run and is reported.
- **The starts.** Each chain's starting π equals `path-profile-modes.csv`'s to rounding, on a fresh design. Recheck when the report is written, since other sessions are editing `src/`.
- **Exponential waits.** Where a start has at least 8 crossings, report the waits' coefficient of variation. A value far from 1 means the wait is not one escape rate, and Z/n is then a summary rather than a rate.

## How the outcomes read

Fixed now. The error is err = ln W − ln h, so a positive err means the reading is too slow. The tolerance is 3 nats, a factor of 20, which is below the 4.6 nats separating the path profile's two rule lines.

For chains: Z is the pooled sweeps to crossing or censoring and n the crossings. With n ≥ 3, the reading is ln(Z/n) with an exact Poisson interval. With n < 3 it is a bound, and with n = 0 the bound is mean wait ≥ Z/3 at 95%.

**Part 1, off-cap pairs, per reading:**
- **90th percentile |err| ≤ 3:** the reading predicts waits on enumerated fixtures. This does not carry to N = 500 alone, since these spaces hold few parallel routes.
- **Wider, but rank correlation with h ≥ 0.9:** it orders pairs without predicting waits. Report the offset and how it moves with depth.
- **Otherwise:** it does not predict waits even here.
- In every case, say which correction is larger, by how much, and whether each grows with depth.

**Chains that can reach their reading** (shallow starts, Part 5, ladder rungs):
- **|err| ≤ 3, widened by the interval:** agrees.
- **err > 3:** too slow, as parallel routes would make it. **err < −3:** too fast, as the basin's mass would make it.

**Deep starts at b = 1:**
- **A crossing whose interval lies wholly below W_class − 3:** the reading is refuted for that start, and a faster way exists. This is the informative outcome. Compare the traced path with the route before reporting anything else.
- **No crossing:** the reading survives and is not confirmed. Class flow is refuted on the slow side only where W_class < ln(Z/3) − 3.
- **Arm S crosses where arm H does not:** τ's movement opens the way, and the held reading is not the sampler's.

**The ladder:**
- **|err| ≤ 3 on at least two-thirds of the rungs with n ≥ 3:** the reading tracks crossing times at the real N with a real barrier, at b < 1.
- **Within a pair, err's spread ≤ 3 with its mean outside ±3:** a constant offset. Report it.
- **err growing toward b = 1:** extrapolating is not supported, and the report says so.
- A rung where chains leave the start's basin by 1% of the first passage has no held start, so the reading does not apply. It drops out and is counted.
- A fit of ln h against b, carried to b = 1, is reported only with at least 3 rungs and residuals ≤ 1 nat. It is labelled an extrapolation and not counted as a measurement.

## Working rules (a shared checkout)

- Work on `main`. Commit by explicit path (`git commit -- paths`), leaving nothing staged: the index is shared. Push nothing. Commits end with the `Co-Authored-By` line the task gave.
- New files only: `tests/regression/BassWaitingTime.cs` and `BassWaitingTimeTests.cs`, heavy parts behind `BASS_QUALIFICATION` = 1, with an accessor elsewhere only if one is needed. Data goes to `issues/regression/bass/waiting-time-*.csv`.
- Build and run into `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-waiting-time\`, set on each `dotnet` command itself. Delete that folder when done. No subagents, no worktrees.
- Leave other sessions' uncommitted notes under `issues/` and the uncommitted `lean/` work alone.
- If a hitting-time solver earns a place in the shared mixing layer, say so and propose it; do not build it here.

## Review before reporting

Read the harness and the report fresh, for correctness and for coverage gaps, and fix what turns up before reporting. Check in particular that W and h are read to the same target, in the same units and at the same τ on every comparison, and that each number quoted comes from the CSVs.

## Report

Append under "Reports": the checks; each part's table; which outcome each rule returned; what was not run and why; and the limits, stating plainly what the deep starts could not confirm. End with "For the main session": what G1's waiting-time bullet, SC-4 and the path profile's "consistent with the chains staying put" should now say.

## Added after registration

Written 2026-11-04, after Parts 1 and 2 and Part 3's deep starts to 10⁷, and before any ladder chain or Part 5 chain ran. The rules above stand and are read as written. These add to them.

- **A second "held start" beside the literal one.** The ladder's drop-out rule, read literally, would drop every rung, since first exit from the start's `Basin` comes within 34–317 sweeps at b = 1 for chains that cross near 10⁶: one spurious term is an exit, and it is gone a few dozen sweeps later. So the chains also record a *departure*, the start of the first stretch of at least 1,000 sweeps outside the start's basin. The report counts drop-outs both ways: by first exit, as written, and by departure, the same 1% against it.
- **More chains on the deep starts (`3x`).** Twelve more per deep start, arm H, to 10⁸ or the budget. Four chains can never reach the 8 crossings the exponential check asks for, and they cost about 1.4 core-hours, since the deep starts turned out to cross. These rows are kept apart from the registered four.
- **Each crossing's gateway.** Each chain from here on records the basin it entered the target from and π at entry, since the deep starts' traces show two doors into a1's basin that stay very differently.
- **The basins' shares.** The flat-histogram rerun also records its largest basins' shares, since those traces lean on `3 × x1·x2 + 2 × x3 + x4 + x5`, which the path profile did not read.

## Reports

### Measured — 2026-11-04

**The route reading does not predict waits, and the chains cross.** From all six deep starts at N = 500 and b = 1, chains enter a1's basin in a mean of about 4.5 million sweeps, e^{15.0}–e^{15.7} by start, against class readings of e^{29}–e^{45}. By the brief's rule that refutes the reading on all six, by 13–29 nats. "Too slow to have been taken in a 2,000-sweep phase" stands, but on the measured waits, about 2,270 times the phase. The tests that could have confirmed the reading did not: the ladder loses 19–20 of its 22 rungs to the held-start rule, and both N = 125 starts read too slow by about 5 nats. Of the two errors, parallel routes win by at least 18–32 nats at N = 500, and the aligned lattice's own are a small part of that on five of the six starts.

| Test | What it can say | What the rule returned |
|---|---|---|
| Part 1: exact, enumerated, off-cap | two-sided, on small spaces | **It does not predict waits even here,** on all eight readings: 90th percentile \|err\| 6.45–8.32, rank correlation 0.51–0.67 |
| Part 2: lattice resistance | part of the credit, as a bound | At most 1.02 nats on 10 of 12 rows; 3.63 on chain3, 12.99 on chain1 |
| Part 3: shallow starts | two-sided, mostly proposal cost | Aligned and partial readings **too slow,** by 6.9–11.5 nats. Exact cuts **agree,** at −0.53 and +1.44 |
| Part 3: deep starts | one-sided | **Refuted on all six** in arm H. Arm S never crosses where arm H does not |
| Part 4: ladder in b | two-sided at the real N, at b < 1 | **Nothing to read on tracking:** only chain0's 2–3 rungs survive the held-start rule. err grows toward b = 1, so **extrapolating is not supported** |
| Part 5: N = 125, carried modes | two-sided, lower confidence | **Too slow** on both, +4.87 and +4.85, with the interval's edge still 4.0 away |

#### Commits and state

`2a704c8` (harness and exact half), `825c185` (additions, before any ladder chain), `bbaf7c9` (all chains, the tally's pins). Tests only, nothing under `src/` or the register touched. Beside the two new files, `BassPathProfile.cs` gained accessors: `Slide`'s optional `visited`, `Neighbours` made internal, and `PathCsv`'s file prefix. `artifacts\chip-waiting-time\` is deleted. The CSVs are `waiting-time-{exact,lattice,readings,validation,chains,traces,ladder,rungs,flathistogram}.csv`, beside this brief. In `chains`, a chain's last row is current.

#### Checks

- **The solver.** On all 616 read rows, h by state reduction agrees with R_eff · M from the voltage, worst 1e-6 at the CSV's six decimals. A dense solve agrees on the 541 rows it ran on. It skips the 72 M ≤ 4 basin rows, by choice at about a minute a target, and the 3 rows past e²⁰. Closed forms hold on two states, a path, and two wells joined by a rate of 10⁻⁴⁰.
- **The kernel, and the production chain against it.** Each row is 1,000 passages, all inside the 99% interval. b = ½ is there for the ladder. This gate let Parts 3–5 run.

| Case | ln h | Kernel walk | Production chain |
|---|---|---|---|
| seed 11, 2033 → 273, structure | 6.994 | 6.997 | 6.978 |
| seed 15, 2033 → 2201, structure | 9.002 | 9.009 | 9.003 |
| seed 13, 2033 → 2201, structure | 11.163 | 11.120 | 11.126 |
| seed 5, 2201 → 273, basin | 7.938 | 7.880 | 7.925 |
| seed 15, 2033 → 2201, structure, b = ½ | 9.424 | 9.441 | 9.436 |
| seed 5, 2201 → 273, basin, b = ½ | 4.877 | 4.915 | 4.891 |

- **A sweep is one attempt.** Every segment of every chain checks it and throws otherwise.
- **Rayleigh.** R_eff ≤ R_lattice ≤ R_route held on every Part 1 row and all 24 of Part 2's.
- **The starts.** 13 structures on a fresh design, worst 4.8e-7, rechecked as this was written. Nothing under `src/` moved since the brief.
- **Exponential waits.** The added chains' cv by start is 0.72–1.32, and 1.05 over all 72, so the deep starts read as one escape rate. Shallow starts 1.01 and 0.91. N = 125: 0.70 and 1.41. Rungs with 8 crossings read 0.46–1.08, except B at ln 16 (1.64) and chain0 at ln 14 and ln 16 (2.72, 2.54), where Z/n is a summary and not a rate.

#### Choices the brief left open

The harness's are in its docstrings. These belong to the analysis:
- A rung drops when the median chain has left by 1% of Z/n. A chain with no departure before crossing counts at its last sweep. Where n = 0, Z/3 stands in for the mean.
- A's class reading never falls below 15.59 on the grid, so its ln 10, 12 and 14 are one rung at b = 0.005. That leaves 22 rungs.
- No start lies in a1's basin, so no reading runs to the exact structure.
- Part `3d`'s departure is not read: those rows predate the column, and a chain resumed outside its basin restarts the stretch. Departures here are from `3x`, the ladder and Part 5.

#### Part 1: exact hitting times

167 pairs both ways, on 12 fixtures of 2,263 states and 4 of 12,051, with ln h from 3.53 to 20.32. err = ln W − ln h:

| Target | Reading | Off-cap: n | p90 \|err\| | Median err | Rank corr. | At-cap: n | p90 \|err\| | Median err | Rank corr. |
|---|---|---|---|---|---|---|---|---|---|
| Structure | exact | 116 | 6.73 | 3.32 | 0.546 | 218 | 4.79 | −0.94 | 0.669 |
| | exact, class | 116 | 6.45 | 2.73 | 0.507 | 218 | 5.18 | −1.18 | 0.642 |
| | aligned | 114 | 6.94 | 3.39 | 0.546 | 220 | 5.82 | −1.06 | 0.593 |
| | aligned, class | 114 | 6.46 | 2.97 | 0.521 | 220 | 6.22 | −1.37 | 0.586 |
| Basin | exact | 93 | 7.90 | 5.20 | 0.618 | 189 | 6.95 | 2.31 | 0.421 |
| | exact, class | 93 | 7.48 | 4.97 | 0.607 | 189 | 6.13 | 1.74 | 0.385 |
| | aligned | 91 | 8.32 | 5.33 | 0.670 | 191 | 6.51 | 1.57 | 0.425 |
| | aligned, class | 91 | 7.73 | 4.99 | 0.659 | 191 | 5.90 | 1.33 | 0.401 |

Every reading returns the rule's third line. Class flow narrows the spread by under a nat, and on the structure target at-cap it does not narrow it at all.

The two corrections, parallel = ln(R_route / R_eff) and mass = ln(M / π(S)):

| Target, route | Cap | Parallel: median (p90) | Mass: median (p90) | Parallel larger on | Lattice's credit, median | Parallel against depth: rank, slope | Mass against depth: rank, slope |
|---|---|---|---|---|---|---|---|
| Structure, exact | off | 6.38 (8.75) | 1.86 (7.06) | 77% | 0.53 | 0.67, 0.30 | −0.48, −0.34 |
| | at | 4.55 (6.38) | 5.77 (9.12) | 41% | 0.54 | 0.39, 0.17 | −0.24, −0.23 |
| Structure, aligned | off | 6.47 (8.82) | 1.84 (6.94) | 75% | 0.40 | 0.70, 0.37 | −0.46, −0.30 |
| | at | 4.32 (7.13) | 5.86 (9.13) | 39% | 0.22 | 0.47, 0.27 | −0.26, −0.22 |
| Basin, exact | off | 6.73 (9.29) | 1.53 (6.00) | 84% | 0.55 | 0.81, 0.39 | −0.53, −0.24 |
| | at | 4.98 (9.92) | 2.43 (6.58) | 77% | 0.72 | 0.77, 0.64 | −0.03, −0.05 |
| Basin, aligned | off | 6.70 (9.44) | 1.53 (5.62) | 85% | 0.32 | 0.86, 0.46 | −0.47, −0.21 |
| | at | 4.59 (9.90) | 2.46 (6.63) | 70% | 0.21 | 0.76, 0.66 | −0.04, −0.06 |

Off-cap, parallel routes are larger by about 4.5–5.2 nats at the median, and they grow with depth. The mass falls with depth. At-cap on the structure target, the mass is larger by about 1.2–1.5. The sign of err turns on depth (aligned, structure, off-cap):

| Route's depth | n | Median err | Parallel | Mass |
|---|---|---|---|---|
| 0–4 | 24 | −2.78 | 3.11 | 6.49 |
| 4–8 | 41 | 3.43 | 6.29 | 2.06 |
| 8–12 | 30 | 3.94 | 7.28 | 1.68 |
| 12 and up | 19 | 5.93 | 7.45 | 0.94 |

So shallow routes read too fast, and deep ones too slow. M ≈ the start basin's mass fails often on small spaces: the ln ratio has median 2.39 (p10 −0.81, p90 5.51), and v over the basin has median 0.71, p10 0.02. W_mass does not rescue the reading: off-cap p90 is 7.65–7.81 on the structure and 8.47–9.02 on the basin, with rank correlation 0.51–0.70.

#### Part 2: the lattice's resistance

From S, on each S10 pair's aligned lattice:

| Pair | d | ln W, route | ln W, lattice | Credit |
|---|---|---|---|---|
| a2 | 5 | 15.30 | 15.30 | 0.00 |
| a1+x6x8 | 5 | 14.47 | 14.47 | 0.00 |
| A | 8 | 44.62 | 43.83 | 0.80 |
| B | 7 | 47.17 | 46.37 | 0.80 |
| chain0 | 10 | 42.89 | 42.23 | 0.66 |
| chain1 | 9 | 36.09 | 23.10 | 12.99 |
| chain2 | 11 | 33.68 | 33.22 | 0.45 |
| chain3 | 9 | 45.51 | 41.87 | 3.63 |
| A, 2/2 | 6 | 48.63 | 48.47 | 0.15 |
| A, d/1 | 10 | 44.44 | 43.42 | 1.02 |
| A, 1/d | 6 | 44.01 | 43.81 | 0.20 |
| A, d/d | 8 | 45.56 | 44.73 | 0.83 |

The lattice's parallel routes are a small part of the credit. Part 3 shows how small.

#### Part 3: shallow starts

Target a1's basin, arm H, 8 chains each, all crossed.

| Start | Reading | W, cut | W, to the exact structure | Hops, slides in the cut | n / Z | ln(Z/n) | cv | err | Outcome |
|---|---|---|---|---|---|---|---|---|---|
| a2 | aligned | 15.30 | 15.30 | 2, 1 | 8 / 367 | 3.83 [3.15, 4.67] | 1.01 | 11.48 | too slow |
| | exact | 3.30 | 11.60 | 1, 0 | | | | −0.53 | agrees |
| a1+x6x8 | aligned | 14.47 | 14.47 | 2, 1 | 8 / 14,503 | 7.50 [6.82, 8.34] | 0.91 | 6.97 | too slow |
| | partial | 14.44 | 14.44 | 2, 1 | | | | 6.94 | too slow |
| | exact | 8.94 | 27.94 | 1, 0 | | | | 1.44 | agrees |

a2's partial reading is its aligned one. The aligned series keeps a slide of 240 and 45 knot steps that no crossing needs. The exact cuts are one hop with no depth, so this is proposal cost and not a barrier: it says the reading's arithmetic holds where nothing else competes, and that the cut matters. Reading to the exact structure, the same starts read 8–19 nats higher.

#### Part 3: deep starts

The reading's side and the registered four chains a start, target a1's basin:

| Start | π(S) | W, route | W, class | Lattice | Arm H: n / Z | ln(Z/n) | W_class − 3 − interval's top | Arm S: n / Z | ln(Z/n) |
|---|---|---|---|---|---|---|---|---|---|
| A | −884.89 | 44.62 | 44.60 | 43.83 | 4 / 18,205,573 | 15.33 [14.39, 16.63] | 24.97 | 4 / 24,801,045 | 15.64 [14.70, 16.94] |
| B | −889.29 | 47.17 | 40.12 | 46.37 | 4 / 7,535,450 | 14.45 [13.51, 15.75] | 21.37 | 4 / 18,985,496 | 15.37 [14.43, 16.67] |
| chain0 | −933.98 | 42.89 | 36.55 | 42.23 | 4 / 24,881,312 | 15.64 [14.70, 16.94] | 16.61 | 3 / 28,547,204 | 16.07 [15.00, 17.65] |
| chain1 | −912.06 | 36.09 | 28.89 | 23.10 | 4 / 21,552,211 | 15.50 [14.56, 16.80] | 9.09 | 2 / 25,144,883 | bound [15.06, 18.46] |
| chain2 | −925.84 | 33.68 | 33.65 | 33.22 | 4 / 9,397,761 | 14.67 [13.73, 15.97] | 14.68 | 4 / 13,765,359 | 15.05 [14.11, 16.35] |
| chain3 | −917.22 | 45.51 | 44.81 | 41.87 | 4 / 34,974,950 | 15.98 [15.04, 17.28] | 24.53 | 3 / 16,593,025 | 15.53 [14.45, 17.10] |

**The rule's outcome.** Every interval lies wholly below W_class − 3, by 9–25 nats, so the reading is refuted for every start and a faster way exists. The lowest class reading over all three alignments changes nothing (B 36.64, chain0 36.22, chain1 27.60, chain2 31.73, chain3 40.77). Arm S's waits sit with arm H's. It runs to 10⁷ only and has 4 chains censored there, so it never crosses where arm H does not, and τ's movement is not what opens the way.

**The added chains (`3x`), kept apart from the four,** and pooled with them. All 96 crossed. The last did so at 25,488,101 sweeps.

| Start | `3x`: n / Z | ln(Z/n) | cv | Pooled: n | ln(Z/n) | cv | err, route | err, class | Shortest wait | First exit, median | Departure, median | Pre-crossing share in the start's basin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 12 / 41,743,963 | 15.06 [14.50, 15.72] | 1.26 | 16 | 15.14 [14.65, 15.70] | 1.15 | 29.49 | 29.47 | 225,866 | 109 | 4,803 | 0.085 |
| B | 12 / 45,638,688 | 15.15 [14.59, 15.81] | 1.32 | 16 | 15.02 [14.53, 15.58] | 1.33 | 32.16 | 25.10 | 269,129 | 80 | 23,741 | 0.646 |
| chain0 | 12 / 44,833,372 | 15.13 [14.58, 15.79] | 0.72 | 16 | 15.29 [14.80, 15.85] | 0.74 | 27.60 | 21.27 | 61,768 | 72 | 7,908 | 0.005 |
| chain1 | 12 / 41,633,155 | 15.06 [14.50, 15.72] | 1.08 | 16 | 15.19 [14.70, 15.75] | 1.11 | 20.90 | 13.70 | 539,121 | 130 | 2,046 | 0.625 |
| chain2 | 12 / 78,920,321 | 15.70 [15.14, 16.36] | 0.82 | 16 | 15.52 [15.04, 16.08] | 0.93 | 18.15 | 18.13 | 266,480 | 144 | 2,063 | 0.001 |
| chain3 | 12 / 65,935,574 | 15.52 [14.96, 16.18] | 1.08 | 16 | 15.66 [15.17, 16.22] | 1.18 | 29.85 | 29.16 | 561,527 | 106 | 5,796 | 0.498 |

Over all six: Z = 435,252,330, n = 96, mean 4,533,878 sweeps, ln 15.33 [15.13, 15.54], cv 1.09. That is 2,267 phases, and about 52 seconds of one core at the harness's 87,324 sweeps a second. On that rate a 2,000-sweep phase crosses with chance 4.4e-4 a chain. The shortest of 96 took 61,768 sweeps, 31 phases.

**Which error wins.** By the identity, ln(R_route / R_eff) = err + ln(M / π(S)) ≥ err, so on these starts parallel routes are worth at least 18–32 nats, whatever the mass is. The mass can only add to that. Since h ≥ π(S)R_eff, the interval's top also bounds the whole network, which sits under the lattice's resistance by at least 27, 31, 25, 6, 17 and 25 nats on A to chain3. Against that, the lattice's own credit is 0.45–13 and class flow's 0.02–7.2. So on every start the way past lies outside L: the lattice's whole resistance cannot reach the wait.

**The traced paths against the routes,** as the rule asks first. On the 44 crossings with a trace (both arms), 2 of 3,256 additions and changes in the last 10⁴ sweeps before entry hold a term of the start's aligned route, so the chains do not walk it. Nor do they find a shallower pass, on the evidence: on arm H's crossings, the widest path through the basins traced after the last state holding both main effects has a median low of −919 (21 chains, range −1,000 to −902; the other 3 had lost both before the window), against the routes' lows of −922 on A, −920 on B and −929 to −959 on the rest. The rows hold basin-changing moves only, so that low is approximate. My reading, not measured: what differs is the number of ways. An unpaid term costs one route its whole measure, 12.9 or 21.3 nats, and the class credit returns at most ln 51 or ln 2,601 of it, while a crossing to a basin can use any knot of that term.

**All six waits agree, whatever the depth.** Routes run 16.5–37.3 nats deep, and the waits differ by under 0.7. Chains leave the start's basin early: median first exit 72–144 sweeps, departure 2,046–23,741. The three starts in B's basin (B, chain1, chain3) spend 50–65% of their pre-crossing sweeps there. So all starts seem to fall in a few thousand sweeps into one common region and leave it at one rate. That is an inference from the waits and shares; no chain records the region.

**The gateways,** over the 72 added chains:

| Entered from | Entries | π at entry, median | Stay over 10⁴ sweeps, median | Wait, median |
|---|---|---|---|---|
| `3 × x1·x2 + 2 × x3 + x4 + x5` | 18 | −899.3 | 0.010 | 1,942,766 |
| a2's basin, `4 × x1·x2 + x2 + …` | 18 | −892.4 | 0.671 | 4,079,445 |
| `x1 + 4 × x1·x2 + …` | 7 | −888.9 | 0.573 | 692,995 |
| `5 × x1·x2 + …` | 2 | −904.5 | 0.543 | 335,162 |
| a1's basin and one term more (16 basins) | 27 | −894.4 | 0.654 | 3,505,832 |

First entry is not always holding: 17 of the 72 stay under 10% of the window, and the door through `3 × x1·x2 + …` stays 1% at the median. Yet none has gone back: 10⁴ sweeps after entry, the 96 arm-H chains sit at (x1, x2, x1·x2) = (0,0,4) on 68, (0,0,3) on 24, (0,1,4) on 3 and (1,0,4) on 1, and none holds both main effects. No chain of any part hit a1's exact structure (0 of 328). `3 × x1·x2 + 2 × x3 + x4 + x5` holds the highest π any arm-H deep chain reached, −877.64, above a1's exemplar at −883.02. The flat histogram gives that basin 0.594% of the mass and 7,906 draws, against a1's 65.914% and 23,288, B's 1.294% and 965, and A's 0.007% and 7,261. The path profile did not read it.

**The mass, read from the chains and the flat histogram.** ln W_mass = W + ln(mass / π(S)):

| Start | W | Chains: ln mass, median (range) | ln(mass / π(S)) | W_mass | Flat histogram: ln(mass / π(S)) | W_mass |
|---|---|---|---|---|---|---|
| a2 | 15.30 | −897.63 (−897.83 to −894.87) | 0.20 | 15.51 | 4.20 | 19.51 |
| a1+x6x8 | 14.47 | −890.89 (−893.90 to −886.06) | 3.35 | 17.82 | 8.17 | 22.65 |
| A | 44.62 | −875.96 (−879.04 to −873.24) | 8.94 | 53.56 | 3.61 | 48.23 |
| B | 47.17 | −878.98 (−881.33 to −875.55) | 10.30 | 57.48 | 6.07 | 53.25 |
| chain0 | 42.89 | −922.03 (−933.77 to −908.84) | 11.95 | 54.84 | no draw | |
| chain1 | 36.09 | −878.80 (−881.64 to −875.80) | 33.26 | 69.35 | no draw | |
| chain2 | 33.68 | −909.67 (−918.91 to −904.12) | 16.17 | 49.85 | no draw | |
| chain3 | 45.51 | −880.36 (−911.11 to −875.54) | 36.86 | 82.37 | no draw | |

The correction moves every reading further from the wait, so it cannot rescue W on these starts. The noise is large:
- A single chain's figure spans about 6 nats on A and B. Shallow chains cross before the mass can be read at all.
- Three starts in one basin read −878.98, −878.80 and −880.36, so the median repeats.
- The flat histogram rests on one draw of each exemplar, a structure selected by having been drawn, so its ratio is biased low. It sits 4.2–5.3 nats under the chains on A and B.
- Where S sits far under its basin's top, as on chain0–3, v at the top is well under 1 and the basin's mass is not M. Those rows do not read M.

#### Part 4: the ladder in b

Rungs were chosen from the readings alone (`waiting-time-rungs.csv`). Each ran 8 chains, stopped at crossing or 10× the class reading. W(b) is not monotone, since the maximin route changes with b: A's route reads 28.2 at b = 0.1, 49.1 at 0.5, 56.2 at 0.6, 43.4 at 0.7 and 44.6 at 1. A \* marks a departure median that half or more of the chains reached only by crossing first.

| Pair | ln | b | Route | Class | n / Z | ln(Z/n) | err, route | err, class | cv | First exit | Departure | Z/n | Drops: exit, departure | Stay |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 10 | 0.005 | 20.00 | 15.59 | 0 / 473,678,128 | ≥ 18.88 | | | | 4 | 4 | ≥ 157,892,709 | yes, yes | |
| A | 16 | 0.030 | 19.79 | 16.01 | 0 / 714,738,128 | ≥ 19.29 | | | | 9 | 9 | ≥ 238,246,043 | yes, yes | |
| B | 10 | 0.070 | 9.74 | 9.74 | 0 / 1,360,952 | ≥ 13.03 | | | | 10 | 10 | ≥ 453,651 | yes, yes | |
| B | 12 | 0.090 | 11.20 | 11.20 | 0 / 5,844,960 | ≥ 14.48 | | | | 10 | 10 | ≥ 1,948,320 | yes, yes | |
| B | 14 | 0.095 | 13.44 | 13.44 | 4 / 42,931,459 | 16.19 [15.25, 17.49] | −2.75 | −2.75 | | 6 | 11 | 10,732,865 | yes, yes | 0.00 |
| B | 16 | 0.100 | 15.69 | 15.69 | 8 / 62,039,299 | 15.86 [15.19, 16.70] | −0.17 | −0.17 | 1.64 | 6 | 6 | 7,754,912 | yes, yes | 0.00 |
| chain0 | 10 | 0.245 | 10.09 | 9.88 | 8 / 11,777 | 7.29 [6.62, 8.13] | 2.80 | 2.59 | 0.77 | 184 | 10,521\* | 1,472 | no, no | 0.82 |
| chain0 | 12 | 0.290 | 11.93 | 11.93 | 8 / 27,866 | 8.16 [7.48, 9.00] | 3.77 | 3.77 | 0.85 | 156 | 6,360\* | 3,483 | no, no | 0.84 |
| chain0 | 14 | 0.330 | 13.99 | 13.99 | 8 / 2,670,385 | 12.72 [12.04, 13.56] | 1.28 | 1.28 | 2.72 | 62 | 1,750 | 333,798 | yes, yes | 0.84 |
| chain0 | 16 | 0.370 | 16.09 | 16.09 | 8 / 225,986 | 10.25 [9.57, 11.09] | 5.84 | 5.84 | 2.54 | 138 | 2,418 | 28,248 | yes, no | 0.85 |
| chain1 | 10 | 0.230 | 16.23 | 9.75 | 5 / 1,116,686 | 12.32 [11.47, 13.44] | 3.91 | −2.57 | | 22 | 22 | 223,337 | yes, yes | 0.00 |
| chain1 | 12 | 0.465 | 12.03 | 12.03 | 4 / 11,108,648 | 14.84 [13.90, 16.14] | −2.81 | −2.81 | | 33 | 33 | 2,777,162 | yes, yes | 0.01 |
| chain1 | 14 | 0.525 | 14.07 | 14.07 | 8 / 5,621,103 | 13.46 [12.78, 14.30] | 0.61 | 0.61 | 0.77 | 76 | 84 | 702,638 | yes, yes | 0.01 |
| chain1 | 16 | 0.575 | 16.07 | 16.07 | 8 / 6,907,124 | 13.67 [12.99, 14.51] | 2.40 | 2.40 | 0.46 | 30 | 54 | 863,390 | yes, yes | 0.00 |
| chain2 | 10 | 0.480 | 14.63 | 9.93 | 1 / 1,577,070 | bound [12.55, 17.95] | | | | 140 | 424 | 1,577,070 | yes, yes | 0.00 |
| chain2 | 12 | 0.530 | 17.74 | 12.05 | 7 / 6,095,893 | 13.68 [12.95, 14.59] | 4.07 | −1.63 | | 83 | 352 | 870,842 | yes, yes | 0.00 |
| chain2 | 14 | 0.375 | 20.28 | 13.81 | 8 / 14,721,684 | 14.43 [13.75, 15.27] | 5.86 | −0.62 | 0.80 | 68 | 151 | 1,840,210 | yes, yes | 0.00 |
| chain2 | 16 | 0.380 | 20.94 | 14.46 | 8 / 6,691,239 | 13.64 [12.96, 14.48] | 7.30 | 0.82 | 1.08 | 84 | 84 | 836,405 | yes, yes | 0.01 |
| chain3 | 10 | 0.200 | 10.42 | 9.80 | 5 / 840,017 | 12.03 [11.18, 13.16] | −1.62 | −2.23 | | 140 | 612 | 168,003 | yes, yes | 0.00 |
| chain3 | 12 | 0.215 | 12.91 | 12.23 | 8 / 2,616,718 | 12.70 [12.02, 13.54] | 0.21 | −0.47 | 0.98 | 78 | 1,664 | 327,090 | yes, yes | 0.00 |
| chain3 | 14 | 0.225 | 14.62 | 13.93 | 8 / 2,179,069 | 12.51 [11.84, 13.35] | 2.10 | 1.41 | 0.98 | 138 | 813 | 272,384 | yes, yes | 0.00 |
| chain3 | 16 | 0.330 | 16.18 | 15.99 | 8 / 12,405,649 | 14.25 [13.58, 15.09] | 1.92 | 1.74 | 0.89 | 240 | 11,117 | 1,550,706 | yes, yes | 0.00 |

**By the rules:**
- **Drop-outs.** By first exit, as the brief wrote it, 20 of 22 rungs drop. By departure, 19. Only chain0's rungs survive, 2 and 3 of them. Within ±3 widened by the interval, the route and class readings each agree on 2 of 2 and 2 of 3 (strictly, 1 of 2 and 1 of 3). That meets the letter of two-thirds on one pair, whose start is two deletions from a1's basin. I do not read it as tracking at the real N with a real barrier: the ladder cannot say, because at b < 1 these starts are not held.
- **Rungs with n < 3.** Five, not counted. A's two and B's lowest two have n = 0, each bounding err ≤ −3.28, which is what no crossing in 8 chains at 10× the reading gives. They lean too fast, as the mass would make them, and they hold no start, with first exit in 4–10 sweeps.
- **A constant offset.** None. chain2's route comes closest, with spread 3.23 and mean 5.74.
- **err toward b = 1.** It grows. Within the rungs err by class rises on B, chain1 and chain3 (rank with b 1.00, 0.80, 1.00), is mixed on chain0 (0.40), and falls on chain2 (−0.50). From each pair's highest rung to b = 1 it jumps: −0.17 to 25.10, 5.84 to 21.27, 2.40 to 13.70, −1.63 to 18.13, 1.74 to 29.16. So extrapolating is not supported.
- **The fit, an extrapolation and not a measurement.** Two pairs qualify, on rungs the held-start rule drops: chain2 (3 rungs, worst residual 0.39) carries to ln h = 12.51 at b = 1, and chain3 (4 rungs, 0.26) to 24.94. Measured there: 15.52 [15.04, 16.08] and 15.66 [15.17, 16.22]. Both miss, by −3.0 and +9.3. chain0 and chain1 fail the residual (2.42, 1.20), and B has 2 rungs.

**If nothing were dropped,** which no rule asks: of the 17 rungs with n ≥ 3, the route is within ±3 on 11 (13 widened) and class flow on 15 (16 widened). All four route misses are too slow (chain0 ln 16, chain2's three). This says less than it looks, because the waits barely follow the readings. Rank correlation with class over the 17 is 0.34, and the slope of ln h on the class reading within a pair is −0.14, 0.65, 0.14, 0.08 and 0.32. The rungs were set at 10–16, and apart from chain0 the waits sit at e^{12}–e^{16} wherever the reading is.

The wait hardly depends on b, and the reading does. On chain1, from b = 0.575 to 1, ln(Z/n) goes from 13.67 to 15.19 while class goes from 16.07 to 28.89. This fits a wait set by finding the terms more than by the likelihood's barrier. That is an inference. At the lowest b the waits are longer again, with A at ≥ 18.9 near b = 0, so they are not monotone in b.

At b < 1, first entry is rarely holding. The median stay is at most 0.013 on 13 of the 17 rungs, and 84 of 122 entries come through `3 × x1·x2 + …`. chain0's four rungs, entering through `x1 + 4 × x1·x2 + …` and `5 × x1·x2 + …`, stay 0.82–0.85.

#### Part 5: N = 125, carried modes

Arm H, 8 chains a start, all crossed. **These modes are carried, not found.**

| Start → target | Data, reading | W, cut (whole route) | Class, cut | Depth, cut | n / Z | ln(Z/n) | cv | err | Interval's edge | Outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| a1 → A | `stream`, aligned from T | 15.18 (19.21) | 15.18 | 4.14 | 8 / 241,804 | 10.32 [9.64, 11.16] | 0.70 | 4.87 | 4.03 | too slow |
| B → a1 | `fresh`, partial+redraw | 20.85 (20.85) | 20.84 | 9.64 | 8 / 71,274,033 | 16.00 [15.32, 16.84] | 1.41 | 4.85 | 4.01 | too slow |

Both err the way parallel routes make them, by much less than at N = 500. But neither start is much of a held mode. Chains climb 23 nats above each carried start (a1 −314.03 to −291.34, B −307.16 to −283.67), so the exemplars are not this surface's modes. The share of pre-crossing sweeps in the start's basin has median 0.011 and 0.001. By the ladder's standard, B → a1 would drop (departure median 1,907 against 1% of 8.9 million) and a1 → A sits at the line (318 against 302). The stays are low, with median 0.07 and 0.003 and 7 of B's 8 entries through `3 × x1·x2 + …`. The mass adds to both: ln W_mass is 17.55 and 26.40. This is the one part that could read a barrier held at b = 1 two-sidedly, and these caveats limit it to the direction.

#### What was not run, and why

The whole plan ran, in the brief's order, on 5.69 core-hours of the ledger's 36: ladder 4.01, deep 0.84, added chains 0.77, N = 125 0.07. The flat-histogram rerun took 11 minutes outside the ledger. The worst case of about 46 that I stated beforehand assumed chains that run to their caps, and these stop at crossing. The deep starts' rung to 10⁸ was not needed past 2.5 × 10⁷.

Not run:
- arm S past 10⁷, by the brief;
- the exact lattices (d up to 19), by the brief;
- the mass on ladder rungs, since no start is held there;
- the wait to *hold* a1's basin, since chains stop 10⁴ sweeps after first entry.

#### Limits

- **What the deep starts cannot do.** They refute, and they did. They confirm nothing: no reading here is shown right for any deep start.
- **What could confirm did not.** The ladder lost its rungs. N = 125 is carried and unheld, and reads too slow. The two agreements at b = 1 are the exact cuts on shallow starts, one hop with no depth. **No test here confirms a route reading across a real barrier.**
- **A crossing is first entry to a1's basin.** It is not holding it, and it says nothing of a1's share. Nothing here tests the 66%. The chains reach the basin and never the exact structure, which is why the cut matters.
- **The scope.** One fixture, one data seed, one flat-histogram run, one exemplar a label, six starts, N = 500, with 2,000 taken as the phase. τ is held in the reading and in arm H.
- **The rates.** Intervals assume exponential waits, which the cv's support at b = 1. Pooling with `3x` is not a registered reading, and the rule's outcomes are read on the four. They agree.
- **The traces** hold only the last 10⁴ sweeps before entry and only basin-changing moves. "No shallower pass" and "the number of ways" are readings of them, not measurements of R_eff.

#### For the main session

**G1, the waiting-time bullet.** Suggested in place of "The waiting time on those routes reads e^{28}–e^{45} …":

> *The waiting time, measured* ([report](waiting-time-brief.md#reports); Friedman, N = 500, one flat-histogram run, one data seed). Untempered chains from the chains' states and from A and B enter a1's basin in a mean of about 4.5 million sweeps, e^{15.0}–e^{15.7} by start (16 chains a start, none censored, waits' cv 0.74–1.33), and τ drawn as in production reads the same. That is about 2,270 times a 2,000-sweep phase and about a minute of one core. The route series reads e^{29}–e^{45} by class flow on these pairs, so it is not a waiting time: it is 13–29 nats too slow, refuted on every start. On enumerated fixtures it does not predict waits either (90th percentile |err| 6.5–8.3 nats off-cap, rank correlation at most 0.67). There its two errors move opposite ways with depth: parallel routes grow, the basin's mass falls, so shallow routes read too fast and deep ones too slow. At N = 500 parallel routes win by at least 18–32 nats. The aligned lattice's own are worth about a nat on most pairs, and the chains cross on other knots, about as low in π as the route. First entry is not always holding. No chain hit a1's exact structure, and no test confirmed a route reading across a real barrier.

Under "Candidate readings", (d) is done. The bullet's "This is a weak reading …" goes.

**The path profile's "consistent with the chains staying put"** (its rule at line 70, its report at lines 470 and 512, and the register). The rule's label survives: every measured wait is past its 9.9 line, with the lowest interval starting at 13.51. "Staying put" no longer describes what the chains do, and the route is no longer the evidence. Suggested:

> Measured first passages are about 2,270 phases, so chains do not reach a1's basin within a phase: about 4 in 10,000 a chain, and the shortest of 96 took 31 phases. They do not stay in the start's basin: median first exit is 72–144 sweeps and departure 2,046–23,741. The route's e^{28}–e^{45} is not the reason. It is refuted by 13–29 nats.

**SC-4.**
- *"… and the best routes found on Friedman avoid them by passing a likelihood valley instead."* Add: those are the best routes through L, read as one series. Chains cross at e^{15} on other knots, and about as low in π, so one route's depth or series does not bound the time to a basin.
- *"Polynomial does not mean fast. … k = 4 costs e^{−20} before any proposal entropy."* That is arithmetic for one path. A's route adds four terms and net k = 2, and the measured rate to a1's basin is e^{−15.1} a sweep. So it should not be read as a rate.
- *The test* (the flow over a basin against N) stands, and this does not read it. N = 125's B → a1 at e^{16.0} beside N = 500's at e^{15.0} is one carried pair and not a slope. The polynomial expectation is neither broken nor confirmed.

SC-5 is untouched.

**The solver.** `HittingReduction` earns a place beside `EliminatedChain` and `ReversibleChain` in `BassKernelSpectrum.cs`. It reads h, R_eff, M and v on any reversible chain far past where a plain solve loses its digits, with stages shared across targets and starts, and it is checked against closed forms, dense solves and both walkers. I propose moving it when a second reading needs it. It is not built there.

**Next readings, for you to weigh:**
- *A series credited with a term's whole measure,* not its ±25-knot class. This tests "the number of ways" directly against Part 1's exact h, with no chains.
- *The wait to hold a1's basin,* with chains run past first entry, since a quarter of entries are touches. This is the figure a lever would have to beat.
- *Waits against N from modes found at each N,* for SC-4's slope.

Which lever, if any, stays yours. The fact it would rest on is above: a crossing costs about 4.5 million sweeps, with a held or a drawn τ.

#### The main session's check, 2026-11-04

Recomputed apart from the harness, from the CSVs: every chain's pooled wait and interval, arm S, the 96's mean and shortest, the shallow starts and N = 125; all of Part 1's table and its depth bins, the identity and the dense solves; the ladder's 20 and 19 drops; the gateways; the 96's final (x1, x2, x1·x2) counts; the 21 traced lows' median and range. All agree. Not recomputed: 2 of 3,256 and the flat histogram's masses against the run itself.

Two limits on the wording, so the register differs from the suggested text:
- **"About as low in π as the route" holds on A and B only.** The traced low by start, against the aligned route's: A −925 to −902 against −922, B −959 to −913 against −920, chain0 −949 to −904 against −959, chain1 −958 to −904 against −929, chain2 −1,000 to −912 against −942, chain3 −938 and −920 against −943. Those on chain0–3 sit mostly above the route's low, and the comparison is not like for like: a trace holds the last 10⁴ sweeps, after the chain has climbed from its start. So "no shallower pass" is not read.
- **2 of 3,256 matches by exact knot.** A chain walking the route's shape on other knots would score the same, so this says the chains are not on these knots, not that they take another path.
