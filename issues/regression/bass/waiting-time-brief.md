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

## Reports
