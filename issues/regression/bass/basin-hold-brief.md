# The route series credited with a term's whole measure, and the wait to hold a1's basin (chip brief)

Date: 2026-11-04. A measurement, tests only: it changes no code under `src/` and no default. Gap G1 of the [lineage-gaps register](lineage-gaps-register.md), candidate readings (f) and (g). It follows the [waiting-time report](waiting-time-brief.md#reports).

- *Why (f).* Parallel routes are worth at least 18–32 nats at N = 500, and the aligned lattice's own are about a nat on most pairs. The report's reading is that the difference is the number of ways: a crossing to a basin can use any knot of a term. Nothing has measured that. Crediting every addition's ±25-knot class to all later edges reaches the waits on chains 1 and 2 and stays 7–20 nats above them on chains 0 and 3.
- *Why (g).* A crossing is first entry, and 17 of 72 entries stay under a tenth of the next 10⁴ sweeps. A lever would have to beat the wait to hold, and nothing has read it. Chains now reach a1's basin, so they can also read a1's share, which rests on one flat-histogram method with 6–21 round trips a chain.
- *The question.* Does the number of ways account for the parallel routes? How long until a chain holds a1's basin, and do chains then agree with the flat histogram?
- *Not in scope.*
  - Any move, kernel, schedule, option or default.
  - Choosing a lever, which is the user's, on this report.
  - Modes found at each N (G1's reading (e)), the fixed-M family, a second data seed.
  - Editing the register or the protolemmas. Put what they should say under "For the main session".
  - Moving `HittingReduction` out of the tests. It stays where it is for this brief.

This section and the next three were written before any of this was computed. The harness's choices, where this leaves them open, go in the report and its docstrings.

**What each test can decide.**

| Test | Size | What it can say |
|---|---|---|
| Part 1: exact resistances by basin | enumerated | two-sided against R_eff, on spaces with small dictionaries |
| Part 2: whole-measure series | S10, no chains | one-sided against the waits, since h ≥ π(S)·R_eff and M is not known; two-sided on A and B only, with a noisy mass |
| Part 3: the route's shape | existing traces | whether chains walk the route's basins on other knots |
| Part 4: the wait to hold | S10, chains | two-sided: chains reach it or give a bound |
| Part 5: a1's share and the balance | S10, chains | two-sided against the flat histogram, if chains hold and return |

Read first: [AGENTS.md](../../../AGENTS.md), [CONTRIBUTING.md](../../../CONTRIBUTING.md), the waiting-time brief with its report and "The main session's check". The harness is yours: `tests/regression/BassWaitingTime.cs`, its CSVs, and `BassPathProfile.cs`.

## What is measured

Every reading is cut at the route's first state in the target's `Basin` and read at the held τ of the earlier readings, unless it says otherwise. Basins are the signature `Basin` writes, so a term's *class under whole measure* is its variable set, with every knot and sign.

**Part 1: exact resistances by basin, no chains.** On Part 1's enumerated pairs, both directions, with target the structure and the basin, off-cap and at-cap apart. Lump the kernel by basin, summing conductances.
- ln R_short, the effective resistance with each basin shorted. By Rayleigh R_short ≤ R_eff.
- ln R_lump, the series over the route's basin sequence with exact lumped flows. A hop inside a basin is not a step of this sequence, so slides drop out. It is neither bound: it is the best one route can read at the level of basins, so it is the ceiling for Part 2's method.
- ln W_own and ln W_held, as Part 2 defines them, on the route's own states.
- The dictionary's size on each fixture.
- The decomposition, exact: ln(R_eff / R_short) is resistance inside basins, finding the terms; ln(R_lump / R_short) is parallel routes between basins. Report each one's median and its rank with depth. The report's inference that the wait is set by finding the terms more than by the likelihood's barrier is read here.

**Part 2: the whole-measure series at N = 500, no chains.** S10, seed 8017's modes and chain states, exemplar 1. The six deep starts and the two shallow ones, on the aligned, partial and exact routes, target a1's `Basin`.
- W_held,25: the path profile's "every addition's class credited to all later edges", recomputed to the same cut. Its earlier values were read to a1's structure.
- W_own: on the route's basin sequence, at each hop, c_e sums the edge's own term over its whole class, with the route's other terms fixed.
- W_held: W_own, with c_e also multiplied, for each term the route added earlier and still holds, by Σ_{t'} π(x with t → t') / π(x) over that term's whole class, taken one term at a time. This evaluates the π ratios, which W_held,25 does not.
- The assumption these make: each basin is at local equilibrium between hops. Part 1 is where that can visibly fail. State how each kind of edge is credited.
- A class past about 10⁵ members may be estimated from a sample, with an interval. A sum of exponentials is carried by its largest terms, so a plain sample reads low: refine around the top, and gate it below.

**Part 3: the route's shape, no chains.** The report's 2 of 3,256 matches by exact knot, which cannot tell another path from the same path on other knots. Read it by basin, per start, on the aligned, partial and exact routes:
- the share of the 72 added chains' gateways (`entry_from`) on the start's route's basin sequence;
- the share of the 44 traced crossings' basin-changing moves with both basins on it.

**Part 4: the wait to hold, N = 500, b = 1.** Target a1's `Basin`. Chains run past first entry to the end of their rung.
- *Arm S, primary:* τ drawn as in production, six deep starts, 8 chains each. This is the chain the sampler runs, and a lever would act on it.
- *Arm H:* τ held, on A, B and chain0, 8 chains each, because Parts 1–3 read that kernel.
- *Definitions.* s_W(t) is the fraction of sweeps in [t, t + W) in the target. A chain *holds at (W, c)* at the first t with s_W(t) ≥ c. Primary: W = 10⁵, 50 phases, and c = 0.5, below a1's share and far above what a touch delivers. Also W = 2,000 and 10⁴, and c = 0.25. T_hold is that t.
- *Per chain:* first entry and its gateway; every entry's stay until hold, so the count of touches; T_hold at each (W, c); occupancy by (x1, x2, x1·x2) family and by the flat histogram's largest basins, before first entry and after hold; the sweeps per second.
- *The common region.* Before first entry and past departure, each chain's sweeps by family and by its 5 largest basins. The report inferred one region; this records it.

**Part 5: a1's share and the balance, from Part 4's chains.**
- *The share.* Past T_hold + W, leaving out the window that selected the chain, each chain's share of a1's basin, with chains as replicates. Read it conditional on holding neither both main effects: the fraction of sweeps in a1's basin among sweeps not holding both an x1 term and an x2 term. Compare with the flat histogram's same conditional, from the rerun's draws. Beside it, the same for its 12 largest basins, and each chain's first half against its second.
- *Which law.* Read the flat histogram's τ first. Arm S's share is the one to compare if its τ is drawn, and arm H's stationary law is a different one, so arm H's share is reported and not compared.
- *The balance.* Past hold, the sweep of the first return to a state holding both main effects, pooled as Z / n. Reversibility and two wells give ln(return) ≈ ln(first entry) + ln(mass of a1's side / mass of the both-main-effects side). Compute the mass ratio from the rerun's draws **before any Part 4 chain runs**, and write it in the CSV.

**The budget.** Parts 1–3 at most 4 core-hours. Parts 4 and 5 at most 24, on 4 threads, below normal priority, one test host. The machine is shared.
- Measure sweeps per second for both arms and state the plan's cost before starting the chains.
- The order is Parts 1–3; the gates; arm S to 2.5 × 10⁷; arm H to 2.5 × 10⁷; both to 10⁸ while the budget lasts, arm S first. If the budget or the machine runs out, stop there and say what was not run.
- Everything appends and skips on rerun. Save each chain's structure and τ every 10⁷ sweeps. A stopped chain resumes from it with a fresh derived seed.
- Keep the data under a few MB: online tallies and coarse blocks, not per-sweep rows.

## What must hold before any reading counts

- **The credit reduces to what exists.** At a window of 0 the credited series returns the route's reading, and at ±25 the path profile's class reading, to rounding.
- **Rayleigh.** R_short ≤ R_eff ≤ R_route on every Part 1 row. A violation stops the run and is reported.
- **The sampled sum.** On at least one one-factor and one two-factor class at N = 500 where the whole sum runs, the estimate's interval covers it.
- **The hold tally.** On an enumerated pair, hold and return tallies on the production chain, arm H, agree with a tally on simulated kernel walks. The long-run share of a basin agrees with exact π(basin) within its interval on at least 3 cases.
- **Arm S's share.** On one enumerated fixture, integrate τ by quadrature over π(S, τ) and check the production chain's share of a basin against it. If the harness cannot reach that, say so: arm S's share then has no check of its own.
- **The starts.** As before, each start's π to rounding on a fresh design, and rechecked when the report is written.
- **Exponential waits.** Report the cv wherever a start has at least 8 holds or returns.

## How the outcomes read

Fixed now. The tolerance is 3 nats, as before. Z, n and the intervals follow the waiting-time brief.

**The waits Part 2 reads against, fixed from the report's pooled chains:**

| Start | ln(Z/n) | Interval | ln(mass / π(S)), chains |
|---|---|---|---|
| A | 15.14 | [14.65, 15.70] | 8.94 |
| B | 15.02 | [14.53, 15.58] | 10.30 |
| chain0 | 15.29 | [14.80, 15.85] | not read |
| chain1 | 15.19 | [14.70, 15.75] | not read |
| chain2 | 15.52 | [15.04, 16.08] | not read |
| chain3 | 15.66 | [15.17, 16.22] | not read |

**Part 1, per reading, against π(S)·R_eff.** err_R = ln W − ln(π(S)·R_eff), since a credit addresses only parallel routes. Read the basin target first, because N = 500 reads it. Report err against h beside it, and W_class as the baseline.
- **90th percentile |err_R| ≤ 3, off-cap:** the reading accounts for the parallel routes on these spaces. It carries to N = 500 only as a method, since the dictionaries are small.
- **Wider, with rank correlation ≥ 0.9:** it orders pairs and does not read the size. Report the offset and how it moves with depth.
- **Otherwise:** it does not account for them even here.
- **If R_lump itself fails the first line,** one basin route cannot read the resistance on these spaces, so say that before Part 2 is read: W_own and W_held have nothing to carry.

**Part 2, per reading and start, against the table.** A reading of π(S)·R_eff must sit at or under ln h.
- **Above the interval's top by more than 3:** whole measure on this route's own terms does not reach the wait, refuted on that start. The way past uses other terms or another shape.
- **At or under that line:** it reaches the wait, and this does not confirm it, since M is unread.
- **On A and B,** add ln(mass / π(S)). Within ±3 of the wait, widened by the interval and the mass's noise of about 6 nats, it agrees. Past that it is too slow or too fast. The other four starts have no mass to add, as the report's table says.
- **Across the starts.** The waits spread 0.64 nats and class flow 15.9. A spread ≤ 3 recovers one rate. A wider one leaves the start-dependence the waits do not have.
- **Overall.** Refuted on at least 5 of the 6: the number of ways on the route's terms is not the explanation. Reaching on all 6: it is sufficient by this reading, and not shown necessary.

**Part 3, per route.** Gateways on the route's basin sequence:
- **At least half:** the chains walk the route's shape on other knots.
- **Under a quarter:** they take another path, so the route's terms are not the way.
- **Between:** mixed, and report which doors account for it.

**Part 4, at (10⁵, 0.5), arm S, per start and pooled.** Hold's wait is ln(Z/n) over T_hold, against first entry in the same chains.
- **Within 1 nat of first entry:** entry is the cost, and touches are minor.
- **1–3 above:** a door's cost, and report it by gateway.
- **More than 3 above:** holding costs more than entering, and 4.5 million sweeps understates what a lever must beat.
- **n < 3:** a bound, by the same rule as before.
- **Across starts:** a spread ≤ 1 nat is one rate. **Arm S against arm H:** apart by more than 3 means τ's movement matters to holding, unlike entry.
- **The common region.** The same (x1, x2, x1·x2) family with a median share ≥ 0.5 past departure on all six starts: there is one region, so name it and give its flat-histogram mass. Otherwise report the table and say the inference is not supported.

**Part 5.**
- **The share,** with at least 8 chains read past hold and no drift between halves beyond their interval. Within 0.05 of the flat histogram's conditional, widened by the chains' interval: the two agree, and a1's share has a second reading. Past it: they disagree, so say which side and which basins carry it. 0.05 reflects how far the flat histogram is known to read on poor runs.
- **Drift, or fewer than 8:** unread, and say how many sweeps it would take.
- **The balance,** with n ≥ 3 returns. Within ±3, widened by the interval: rates and masses agree. Past it: they disagree, so the two-well picture or the masses are wrong. With n < 3: a bound on the return, and whether that already excludes the ratio.
- The share and the balance are read on one fixture, one data seed and one flat-histogram run. Say so wherever they are quoted.

## Working rules (a shared checkout)

- Work on `main`. Commit by explicit path (`git commit -- paths`), leaving nothing staged: the index is shared. Push nothing. Commits end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- New code goes in `tests/regression/BassBasinHold.cs` and `BassBasinHoldTests.cs`, heavy parts behind `BASS_QUALIFICATION` = 1, with accessors in your existing harness where needed. Data goes to `issues/regression/bass/basin-hold-*.csv`.
- Build and run into `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-basin-hold\`, set on each `dotnet` command itself. Delete that folder when done. No subagents, no worktrees.
- Leave other sessions' uncommitted notes under `issues/` and the uncommitted `lean/` work alone.

## Review before reporting

Read the harness and the report fresh, for correctness and for coverage gaps, and fix what turns up before reporting. Check in particular that every comparison is read at one τ, to one cut and in one unit; that the share leaves out the selecting window; that the mass ratio's commit precedes the first Part 4 chain; and that each number quoted comes from the CSVs.

## Report

Append under "Reports": the checks; each part's table; which outcome each rule returned; what was not run and why; and the limits. End with "For the main session": what G1's waiting-time bullets and candidate readings (f) and (g) should now say, whether SC-4 changes, and the fact a lever would rest on. Which lever, if any, stays the user's.

## Reports
