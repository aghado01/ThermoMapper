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

### Measured — 2026-11-04

**Entry is the cost, a1's share has a second reading, and the number of ways is not shown to be the explanation.** With τ drawn as in production, all 48 chains from the six deep starts hold a1's basin, in a mean of about 5.6 million sweeps, e^{15.54} [e^{15.25}, e^{15.84}]. That is first entry's own wait, to 0.01 nats. Past the hold they put a1's share, among sweeps not holding both main effects, at 0.6723 ± 0.0018 against the flat histogram's 0.6712, with no drift. The two-well balance does not close as the brief wrote it: first returns come at e^{13.24}, 6.29 nats under the e^{19.53} it predicts, because a first return is a touch. Crediting an edge's own term with its whole measure is refuted on all six starts on every route. Also crediting the terms held earlier reaches the waits on all six on the exact route only, with a spread of 12–23 nats by route where the waits spread 0.64. And R_lump itself fails Part 1's first line, so by the brief's rule W_own and W_held have nothing to carry. The share and the balance are read on one fixture, one data seed and one flat-histogram run.

| Test | What it can say | What the rule returned |
|---|---|---|
| Part 1: exact resistances, off-cap, basin target | two-sided, 48-term dictionaries | **R_lump fails the first line:** 90th percentile \|err_R\| 5.21 and 5.24, rank 0.92 and 0.95. W_route, W_class, W_own and W_held are on the second line, W_held,25 on the third. At cap, all on the third |
| Part 1: decomposition | exact, same spaces | Inside basins, median 1.97 nats, falling with depth. Between basins, 4.37–4.64, growing with depth. **The inference that finding the terms sets the wait is not supported here,** and these dictionaries cannot read it at N = 500 |
| Part 2: whole-measure series | one-sided | **W_own refuted on 6 of 6** on every route. W_held refuted on 2, 1 and 0 (aligned, partial, exact): **reaches on all six only on the exact route,** spread 12.10. W_held,25 refuted on 6, 5 and 4 |
| Part 3: the route's shape | existing traces | No start reaches half. chain0's aligned and partial routes are **under a quarter** (2 of 12); the other 16 rows are **between** |
| Part 4: the wait to hold, arm S | two-sided | **Within 1 nat on every start and pooled, +0.01: entry is the cost.** Spread across starts 1.15, so the rule's one rate is not returned, on intervals that all overlap. Arm S against arm H within 0.87. **Common region not supported** |
| Part 5: the share | two-sided | **The two agree,** 0.6723 ± 0.0018 against 0.6712, 48 chains, drift −0.0016 ± 0.0036 |
| Part 5: the balance | two-sided | **They disagree:** 13.24 [12.96, 13.55], n = 48, against 19.53, 6.29 under |

#### Commits and state

`f4918af` holds the harness, Parts 1–3, the gates and the mass ratio. The later commit holds this report, the chains' CSVs, the budget's ceiling and the rungs. Tests only: nothing under `src/`, the register or the protolemmas. `StepRow` in `BassWaitingTime.cs` gained `ClassSize`; the rest is in `BassBasinHold.cs` and `BassBasinHoldTests.cs`. `artifacts\chip-basin-hold\` is deleted. The CSVs are `basin-hold-{classes,exact,hops,series,shape,flathistogram,validation,chains,visits,blocks}.csv`, beside this brief. In `chains`, a chain's last row is current.

#### Checks

- **The credit reduces.** 48 series. At a window of 0 it returns the route, and at ±25 the path profile's class reading, to 10⁻⁶. Credited live against from the stored rows agrees to 10⁻⁶. Against the path profile's printed column it is within 0.046, its one-decimal rounding.
- **Rayleigh.** Held on all 668 rows: least margins 0.0169 nats for R_eff over R_short and 0.406 for R_route over R_eff. R_lump over R_short reaches −0.000000, equal to rounding.
- **The sampled sum.** Nothing to read: no class was sampled. All 352 class sums ran whole, the largest over 10⁶ members, and no route holds a three-factor term.
- **The hold tally.** Two kernel pairs (seeds 5 and 6, 2033 → 273), 1,000 passages each, kernel walks against the production chain at held τ: first entry, hold, hold at (200, ¼), return, entries and touches agree, the largest z 2.11, on entries. These spaces have almost no touches (0, and 0.001 against 0.005 a passage), so the touch count rests on `Tally_AgreesWithABruteForceRead`, which also covers a stop and resume. The gate's first run failed with z = NaN on touches, all zero at windows of 20, 100 and 400 with a gap of 50. That was the check and not the tally. It now uses windows of 50, 200 and 1,000 with a gap of 20, and z = 0 for equal means, and the failed run's 5 rows are dropped.
- **Shares, arm H.** Exact π(basin) against 16 chains of 3 × 10⁶: 0.8153 and 0.8155, 0.8319 and 0.8319, 0.3718 and 0.3714, 0.4509 and 0.4507, all within 99%.
- **Arm S's share has its own check.** On seed 5, τ integrated over 1,201 points of ln τ gives 0.8542, and the production chain 0.8542 ± 0.0004.
- **The starts.** 9 structures against their stored π on a fresh design, to six decimals, before the chains and again when this was written, on the later build.
- **The mass ratio came first.** `f4918af` is stamped 17:04:17 and the first Part 4 chain started at 17:04:26. Two 10⁶-sweep rate chains per arm ran after the ratio was computed and before that commit, and nothing is read from them but the rate.
- **The chains' books.** On all 72 final rows, block sums equal the post-hold tallies, family counts equal their sweeps, pre-entry sweeps equal first entry less one, and post-hold sweeps start exactly at T_hold + W. No file holds a duplicate key.
- **Which τ.** The flat histogram's τ is drawn (its CSV's first row), so arm S's share is the one compared, and arm H's is reported only.
- **A host this chip did not start.** At 17:50 a test host with `BASS_QUALIFICATION` on reran Part 1, Part 2, the credit gate and three share gates against these files. The chain files show no trace of it. Its 7 validation rows matched or read nothing new (a Rayleigh row of `inf`, since every row was already there), and they are dropped.

#### Choices the brief left open

- *Classes* ran whole rather than sampled, since the cost allowed it: 124,449 members a second on four threads, and Part 2 took 0.55 core-hours.
- *How each edge is credited.* W_own: addition and deletion alike, c_e = Σ_{t′} min(π(lower) + ln q_b, π(lower + t′) + ln q_d) over the term's class, the lower state fixed, slides dropped. W_held: that, times Σ_{t′} π(x with t → t′) / π(x) at the lower state for each term added earlier and still held, one term at a time. W_held,25: the path profile's arithmetic, class flow on an addition and ln(class size) on every later step, slides included.
- *Off-cap* means the route's highest M is under the cap, as before.
- *The return's side.* The primary ratio is ln(not both / both) = 4.003, the same side the share conditions on. ln(neither / both) = 3.995 and ln(a1's basin / both) = 3.604 sit beside it in the CSV.
- *A visit* ends at 1,000 sweeps outside, and a touch is a visit closed before the hold's window.
- *T_hold* is the window's start, as the brief defines it, so it can precede first entry, by at most 39,831 sweeps here.
- *Halves* are by 10⁶-sweep blocks past T_hold + W. Waits take the exact Poisson 95% interval on Z / n, and shares Student's t over chains.
- *Past 2.5 × 10⁷* chains advance in rungs of 5 × 10⁶, not straight to 5 × 10⁷. A straight rung would have run four starts' chains and left the rest behind. The budget became a ceiling at the same time, counting the segments a check leaves in flight.

#### Part 1: exact resistances by basin

668 rows on Part 1's pairs, both directions, 564 read. 52 start in the target's basin, and 52 structure targets lie in the start's basin, so they have no R_short. The dictionary is 48 terms on both caps, with 2,263 states and 18 basins at M ≤ 3, 12,051 and 27 at M ≤ 4.

**R_lump fails the first line, so one basin route cannot read the resistance on these spaces, and W_own and W_held have nothing to carry to Part 2.** Part 2 below is read with that said.

Basin target, off-cap. err_R = ln W − ln(π(S)·R_eff). Each cell is the 90th percentile of |err|, then the median of err, then the rank correlation; the last column is the rank of err_R with depth on the exact route, then the aligned.

| Reading | err_R, exact route, n = 93 | err_R, aligned route, n = 91 | err_h, exact route | Line | err_R with depth |
|---|---|---|---|---|---|
| W_route | 9.29, 6.73, 0.90 | 9.44, 6.70, 0.91 | 7.90, 5.20, 0.62 | second | +0.81, +0.86 |
| W_class | 8.72, 6.52, 0.90 | 8.76, 6.48, 0.90 | 7.48, 4.97, 0.61 | second | +0.77, +0.81 |
| W_held,25 | 7.53, 5.29, 0.86 | 7.86, 5.27, 0.84 | 6.45, 3.57, 0.54 | third | +0.43, +0.51 |
| R_lump | 5.21, 2.79, 0.92 | 5.24, 2.49, 0.95 | 5.72, 1.38, 0.64 | second | +0.82, +0.83 |
| W_own | 7.39, 4.83, 0.90 | 7.39, 4.82, 0.91 | 6.46, 3.31, 0.60 | second | +0.81, +0.86 |
| W_held | 5.73, 3.58, 0.91 | 6.12, 3.57, 0.91 | 5.31, 2.18, 0.59 | second | +0.58, +0.62 |

Every reading orders the pairs at best, and its offset grows with depth over 0–17.9 nats, so none reads the size. W_held comes closest to R_lump: 0.78 nats above it at the median on the exact route (90th percentile 3.80), against W_own's 2.26 (4.80). So the assumption of local equilibrium costs about a nat there. At cap every reading is on the third line, with ranks 0.42–0.54 and 90th percentiles 4.91–9.92, and R_lump reads under R_eff at the median (−0.15 and −0.76), which it may, being neither bound. With the structure as target, off-cap, R_lump still fails the first line (4.44 and 5.01).

**The decomposition.** ln(R_eff / R_short) is resistance inside basins, finding the terms. ln(R_lump / R_short) is parallel routes between basins. Median [p10, p90], then rank with depth:

| Target, cap | Route | Inside | Between |
|---|---|---|---|
| Basin, off-cap | exact | 1.97 [0.68, 4.11], −0.44 | 4.64 [0.46, 7.23], +0.57 |
| Basin, off-cap | aligned | 1.97 [0.67, 3.88], −0.43 | 4.37 [0.25, 7.17], +0.69 |
| Basin, at cap | exact | 2.17 [0.27, 6.95], −0.66 | 3.10 [0.01, 5.52], +0.17 |
| Basin, at cap | aligned | 2.17 [0.27, 6.95], −0.63 | 0.92 [0.00, 5.10], +0.22 |
| Structure, off-cap | exact | 2.82 [1.76, 5.12], −0.26 | 4.66 [1.49, 7.23], +0.52 |
| Structure, off-cap | aligned | 2.82 [1.75, 5.10], −0.26 | 4.39 [0.44, 7.17], +0.64 |
| Structure, at cap | exact | 4.81 [2.46, 7.63], +0.08 | 3.79 [0.05, 6.16], +0.13 |
| Structure, at cap | aligned | 4.90 [2.47, 7.61], +0.05 | 2.89 [0.00, 5.37], +0.25 |

Off-cap, the deeper the pair, the less the resistance inside basins and the more the parallel routes between them. So the waiting-time report's inference, that the wait is set by finding the terms more than by the likelihood's barrier, is not supported on these spaces. It is not refuted at N = 500 either: a class here holds at most 48 terms, and a two-factor class at N = 500 holds 10⁶. This also reads G1's candidate (h), a whole-term redraw: inside basins there are about 2 nats to save on these spaces, and fewer at depth.

#### Part 2: the whole-measure series at N = 500

Read with Part 1's warning. Cut at the route's first state in a1's `Basin`, at the held τ, against the fixed table. The line is the interval's top plus 3. **Bold** is refuted.

| Start | Line | Aligned: W_held,25, W_own, W_held | Partial: W_held,25, W_own, W_held | Exact: W_held,25, W_own, W_held |
|---|---|---|---|---|
| A | 18.70 | **30.10**, **33.06**, 15.63 | **22.78**, **33.15**, 9.02 | **22.78**, **33.15**, 9.02 |
| B | 18.58 | **33.46**, **35.54**, **26.21** | **22.34**, **32.17**, 14.09 | **22.34**, **32.17**, 14.09 |
| chain0 | 18.85 | **35.44**, **30.96**, **29.73** | **35.02**, **30.61**, **29.32** | **25.55**, **30.86**, 9.09 |
| chain1 | 18.75 | **18.76**, **23.13**, 10.58 | 13.74, **21.84**, 5.83 | 10.99, **21.84**, 1.99 |
| chain2 | 19.08 | **23.52**, **21.13**, 9.80 | **21.60**, **19.22**, 8.38 | 13.86, **19.22**, 8.27 |
| chain3 | 19.22 | **34.69**, **33.47**, 9.23 | **33.16**, **29.40**, 6.26 | **22.78**, **29.40**, 9.35 |
| a2 | | 15.30, −5.23, −5.23 | 15.30, −5.23, −5.23 | 3.30, −3.28, −3.28 |
| a1x | | 14.47, −3.40, −3.40 | 14.44, −3.41, −3.41 | 8.94, −3.15, −3.15 |

| Reading | Aligned: refuted, spread | Partial: refuted, spread | Exact: refuted, spread | Overall |
|---|---|---|---|---|
| W_class, the baseline | 6, 15.92 | 6, 17.13 | 6, 17.13 | refuted |
| W_held,25 | 6, 16.68 | 5, 21.27 | 4, 14.56 | refuted on two routes, 2 reach on exact |
| W_own | 6, 14.41 | 6, 13.94 | 6, 13.94 | **refuted on every route** |
| W_held | 2, 20.49 | 1, 23.49 | 0, 12.10 | **reaches on all six on exact only** |

- **W_own.** Refuted on 6 of 6: the number of ways of the edge's own term, with the route's other terms fixed, is not the explanation. Two margins are thin: chain2's partial and exact are over by 0.14, and W_held,25 on chain1's aligned route by 0.01.
- **W_held.** Sufficient by this reading on the exact route, not shown necessary, and it does not confirm anything since M is unread. It moves by up to 20.64 nats with the route on one start (chain0, 29.73 aligned and 9.09 exact), so it is not a property of the pair.
- **One rate.** No reading recovers one: the smallest spread is 12.10, where the waits spread 0.64. The start-dependence the waits do not have stays.
- **A and B with the mass.** The window is the wait ± 3, widened by the interval and about 6 nats of the mass's noise: up to 24.70 on A and 24.58 on B. W_held agrees on partial and exact (A 17.96, B 24.39) and on A's aligned route (24.57); B's 24.39 and A's 24.57 lie within 0.2 of the window's top. B on aligned is too slow at 36.51. W_own and W_held,25 are too slow throughout. A window this wide tests little.
- **The shallow starts.** W_own and W_held read below 0 on a2 and a1x. Those are not waits: a class's credit can run past π(S), and h ≥ 1. W_held,25 there equals W_route: the cut leaves no later edge to credit.

#### Part 3: the route's shape

The 72 added chains' gateways, 12 a start, and the 44 traced crossings' basin-changing moves.

| Start | Route | Basins on the sequence | Gateways on it | At the last door | Moves on it / all | Rule |
|---|---|---|---|---|---|---|
| A | aligned | 7 | 4 | 4 | 78 / 1,216 | between |
| A | partial | 9 | 4 | 4 | 78 / 1,216 | between |
| A | exact | 8 | 4 | 4 | 76 / 1,216 | between |
| B | aligned | 5 | 3 | 3 | 12 / 1,094 | between |
| B | partial | 5 | 3 | 3 | 12 / 1,094 | between |
| B | exact | 5 | 3 | 3 | 12 / 1,094 | between |
| chain0 | aligned | 5 | 2 | 2 | 2 / 1,223 | **under a quarter** |
| chain0 | partial | 5 | 2 | 2 | 2 / 1,223 | **under a quarter** |
| chain0 | exact | 11 | 3 | 1 | 23 / 1,223 | between |
| chain1 | aligned | 7 | 5 | 3 | 7 / 846 | between |
| chain1 | partial | 8 | 5 | 3 | 7 / 846 | between |
| chain1 | exact | 9 | 5 | 3 | 7 / 846 | between |
| chain2 | aligned | 9 | 3 | 3 | 4 / 1,132 | between |
| chain2 | partial | 9 | 3 | 3 | 4 / 1,132 | between |
| chain2 | exact | 10 | 3 | 3 | 4 / 1,132 | between |
| chain3 | aligned | 6 | 3 | 3 | 2 / 1,053 | between |
| chain3 | partial | 6 | 3 | 3 | 2 / 1,053 | between |
| chain3 | exact | 8 | 3 | 3 | 2 / 1,053 | between |

**Which doors account for it.** Three doors, one hop from a1's basin, carry the mix: 3 × x1·x2 + 2 × x3 + x4 + x5, x1 + 4 × x1·x2 + 2 × x3 + x4 + x5, and 4 × x1·x2 + x2 + 2 × x3 + x4 + x5. No route holds all three, and a route's count is mostly whichever it holds as its last door: 18 of the 20 gateways on the aligned routes. The doors off the routes are mostly those same three and a1 + x6·x8. Earlier on the route the chains are not there: at most 6.4% of basin-changing moves have both basins on the route, and 1.1% or fewer on five of the six starts' aligned routes. So chains share the route's last step and not its shape.

#### Part 4: the wait to hold

Arm S, primary, τ drawn, 8 chains a start, none censored. Hold is at (10⁵, ½). Touches are the chains with any, then their total.

| Start | First entry, ln(Z/n) | Hold, ln(Z/n) | Hold − entry | cv of hold | cv of first entry | Touches |
|---|---|---|---|---|---|---|
| A | 16.00 [15.33, 16.84] | 16.00 [15.32, 16.84] | −0.00 | 0.97 | 0.97 | 1, 1 |
| B | 15.37 [14.69, 16.21] | 15.38 [14.70, 16.22] | +0.01 | 0.67 | 0.68 | 3, 53 |
| chain0 | 15.69 [15.01, 16.53] | 15.73 [15.05, 16.57] | +0.04 | 0.48 | 0.57 | 4, 28 |
| chain1 | 15.03 [14.35, 15.87] | 15.03 [14.35, 15.87] | +0.00 | 1.53 | 1.54 | 3, 17 |
| chain2 | 14.85 [14.18, 15.69] | 14.85 [14.17, 15.69] | −0.00 | 0.87 | 0.87 | 1, 4 |
| chain3 | 15.75 [15.07, 16.59] | 15.75 [15.07, 16.59] | −0.00 | 1.15 | 1.15 | 2, 11 |
| Pooled, n = 48 | 15.53 [15.25, 15.83] | 15.54 [15.25, 15.84] | +0.01 | 1.03 | 1.04 | 14, 114 |

Arm H, held τ, on the three starts Parts 1–3 read, 8 chains each, none censored:

| Start | First entry | Hold | Hold − entry | cv of hold | Arm S's hold − arm H's |
|---|---|---|---|---|---|
| A | 15.13 [14.46, 15.97] | 15.13 [14.46, 15.97] | +0.00 | 1.23 | +0.87 |
| B | 15.25 [14.57, 16.09] | 15.26 [14.58, 16.10] | +0.01 | 0.76 | +0.12 |
| chain0 | 15.53 [14.85, 16.37] | 15.75 [15.07, 16.59] | +0.22 | 0.66 | −0.02 |
| Pooled, n = 24 | 15.32 [14.92, 15.76] | 15.42 [15.02, 15.86] | +0.10 | 0.85 | |

- **By the rule: within 1 nat on every start and pooled, so entry is the cost and touches are minor.** The 4.5 million sweeps did not understate it. Arm S's mean hold is 5.58 × 10⁶ and its mean entry 5.54 × 10⁶; arm H's are 4.96 × 10⁶ and 4.50 × 10⁶.
- **Other windows and levels** do not move it. Arm S pooled at (2,000, ¼ and ½): 15.531, 15.538. At 10⁴: 15.539, 15.540. At 10⁵: 15.528, 15.535. Arm H: 15.373, 15.374; 15.373, 15.374; 15.406, 15.417.
- **Across starts** the spread is 1.15, over the rule's 1, so it does not return one rate. It does not show several either: all six intervals overlap, each about 1.5 nats wide at 8 chains.
- **Arm S against arm H:** at most 0.87 apart, so τ's movement does not matter to holding, as it does not to entry.
- **Touches.** 15 of 48 first entries stay under a tenth of the next 10⁴ sweeps (8 of 24 in arm H), as the earlier report found. They are not lost: a touch lasts a median of 31.5 in-target sweeps (39 in arm H; largest 5,475 and 26,902), and on 64 of 72 chains the hold starts within 10⁵ sweeps of first entry. The longest gaps are 1.95 × 10⁶ in arm S and 5.62 × 10⁶ in arm H.
- **The crossing itself is short.** From a chain's last sweep holding both main effects to its first entry, the median is 3,336 sweeps in arm S (largest 47,700) and 673 in arm H (largest 23,799). The other 5.5 million sweeps go before a crossing that commits starts.
- **Leaving the start.** Median first exit is 154.5 sweeps and departure 13,718 (largest 86,765) in arm S; 128.5 and 9,550.5 (72,745) in arm H.

**The gateways, arm S,** of 48: the three doors take 36, namely 3 × x1·x2 + 2 × x3 + x4 + x5 on 18, x1 + 4 × x1·x2 + 2 × x3 + x4 + x5 on 10, and 4 × x1·x2 + x2 + 2 × x3 + x4 + x5 on 8. Then a1 + x6·x8 on 4, and 8 others once each. Arm H: the same three on 19 of 24.

**The common region.** Sweeps before first entry and past departure, 2.65 × 10⁸ in arm S. Family is (x1, x2, x1·x2).

| Start | 1/1/3, median over chains | 1/1/3, pooled | Chains at ≥ 0.5 | Next family, median | Largest basin's pooled share |
|---|---|---|---|---|---|
| A | 0.874 | 0.899 | 7 | 1/1/2, 0.058 | 0.696 |
| B | 0.970 | 0.972 | 8 | 1/2/3, 0.011 | 0.719 |
| chain0 | **0.410** | 0.554 | 3 | 0/2/4, 0.038 | 0.422 |
| chain1 | 0.907 | 0.906 | 8 | 1/2/3, 0.013 | 0.693 |
| chain2 | **0.448** | 0.635 | 4 | 1/1/4, 0.204 | 0.489 |
| chain3 | 0.911 | 0.888 | 6 | 1/1/4, 0.033 | 0.684 |

By the rule, no family has a median share ≥ 0.5 on all six, so **the inference of one region is not supported.** Beside the rule: the largest basin on all six starts is the same one, x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5, which is the flat histogram's fifth largest at mass 0.0129; family 1/1/3 has mass 0.0166, most of the both-main-effects side's 0.0179. chain0's next basins are 4 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5 (0.111), x1 + 4 × x1·x2 + x2 + 2 × x3 + x4 + x5 (0.066) and 5 × x1·x2 + x2 + 2 × x3 + x4 + x5 (0.064), and chain2's is x1 + 4 × x1·x2 + x2 + 2 × x3 + x4 + x5 (0.121). The other four basins of each start's five are in `basin-hold-chains.csv`, as `region_basins`.

That basin is three hops from a1's: lose x1, lose x2, gain an x1·x2, in some order. The three doors are that exchange's last stops. Occupancy says chains spend the wait around it. No chain times that step, so that it is the slow one remains an inference.

#### Part 5: a1's share and the balance

**The share, arm S.** 48 chains read past T_hold + W, 2.11 × 10⁹ sweeps, at least 2.5 × 10⁷ a chain. Conditional on not holding both main effects: 0.6723 ± 0.0018, sd 0.0062 across chains, against the flat histogram's 0.671174. The first half reads 0.6731 and the second 0.6715, a drift of −0.0016 ± 0.0036, within its interval. **The difference is +0.0011, inside the chains' own interval and far inside 0.05: the two agree, and a1's share has a second reading.** By start: A 0.6691, B 0.6752, chain0 0.6731, chain1 0.6718, chain2 0.6698, chain3 0.6747, each ± 0.004–0.006.

The flat histogram's 12 largest beside it. Conditional on not both, but for the fifth, which holds both main effects and reads its share of all sweeps. Arm H is a different law: reported, not compared.

| # | Basin | Flat | Flat's four chains | Arm S ± 95% | Arm H |
|---|---|---|---|---|---|
| 1 | 4 × x1·x2 + 2 × x3 + x4 + x5 (a1) | 0.6712 | 0.636–0.685 | 0.6723 ± 0.0018 | 0.5996 ± 0.0030 |
| 2 | a1 + x6·x8 | 0.1079 | 0.092–0.137 | 0.1082 ± 0.0009 | 0.1181 |
| 3 | a1 + x10·x2·x9 | 0.0172 | 0.012–0.021 | 0.0156 ± 0.0005 | 0.0168 |
| 4 | a1 + x3 | 0.0139 | 0.009–0.023 | 0.0121 ± 0.0002 | 0.0158 |
| 5 | x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 0.0129 | 0.003–0.023 | 0.0019 ± 0.0024 | 0.0023 |
| 6 | a1 + x8 | 0.0104 | 0.010–0.011 | 0.0103 ± 0.0001 | 0.0130 |
| 7 | a1 + x5 | 0.0072 | 0.006–0.008 | 0.0075 ± 0.0000 | 0.0098 |
| 8 | a1 + x10·x9 | 0.0071 | 0.006–0.008 | 0.0075 ± 0.0001 | 0.0090 |
| 9 | a1 + x10·x6·x8 | 0.0067 | 0.004–0.009 | 0.0096 ± 0.0003 | 0.0108 |
| 10 | a1 + x6·x8·x9 | 0.0066 | 0.003–0.014 | 0.0047 ± 0.0001 | 0.0051 |
| 11 | 3 × x1·x2 + 2 × x3 + x4 + x5 | 0.0060 | 0.005–0.008 | 0.0066 ± 0.0007 | 0.0043 |
| 12 | 3 × x1·x2 + 2 × x3 + x4 + x5 + x6·x8 | 0.0058 | 0.001–0.010 | 0.0052 ± 0.0007 | 0.0045 |

"a1 + t" is a1's basin with one more term of t's variables. On the not-both side the chains sit inside the flat histogram's own spread on 10 of 11, the eleventh (#9) just above it. The fifth, on the other side, sits below all four of its chains.

**Arm H's share** is 0.5996 ± 0.0030, 24 chains, drift −0.0037 ± 0.0074, by start 0.5967, 0.6018 and 0.6004.

**The balance, as the brief wrote it.** The first return past T_hold + W to a state holding an x1 term and an x2 term, Z / n, against ln(first entry) + 4.003:

| | Return, ln(Z/n) | n | cv | First entry | Predicted | Return − predicted |
|---|---|---|---|---|---|---|
| A | 13.32 [12.64, 14.16] | 8 | 0.92 | 16.00 | 20.01 | −6.69 |
| B | 12.83 [12.15, 13.67] | 8 | 0.75 | 15.37 | 19.37 | −6.55 |
| chain0 | 13.29 [12.61, 14.13] | 8 | 0.65 | 15.69 | 19.69 | −6.41 |
| chain1 | 13.38 [12.70, 14.22] | 8 | 0.37 | 15.03 | 19.03 | −5.66 |
| chain2 | 13.54 [12.87, 14.38] | 8 | 1.36 | 14.85 | 18.86 | −5.31 |
| chain3 | 12.90 [12.23, 13.74] | 8 | 1.17 | 15.75 | 19.75 | −6.84 |
| Pooled | 13.24 [12.96, 13.55] | 48 | 0.97 | 15.53 | 19.53 | −6.29 |

**By the rule, past ±3 widened by the interval: they disagree, so the two-well picture or the masses are wrong.** Arm H's first return reads 12.95 [12.56, 13.40], n = 24, cv 1.18, with no ratio to set it against. I leave the rule as written and read it as the two-well picture failing for this set, for these reasons:

- A first return is a touch. Outside long stays, both main effects take 290,434 of 2.11 × 10⁹ sweeps, and with a first return every 5.6 × 10⁵ sweeps that is about 80 sweeps a touch, if later touches are like the first. On chains with no long stay, 98% of those sweeps are in families 1/1/3 and 1/1/4. Holding both main effects is a thin shell around a1's side as well as a well, so its first-hit time is not the reverse of a crossing from deep in that well.

**Beside the rule, and after seeing it fail,** so this does not replace it. Long stays read from the blocks, a 10⁶-sweep block at least a tenth in both main effects, in runs:

| Arm | Stays, by chain and the first block's end: sweeps in both | Count | ln(sweeps past hold / count) |
|---|---|---|---|
| S | start chain1's chain 5, 44 × 10⁶: 3,205,480; start chain2's chain 2, 7 × 10⁶: 954,516; start chain3's chain 1, 49 × 10⁶: 1,790,842, open at the end | 3 | 20.37 [19.30, 21.95] |
| H | start A's chain 5, 6 × 10⁶: 784,171; start B's chain 7, 18 × 10⁶: 439,607; start A's chain 1, 25 × 10⁶: 639,734, open at the end | 3 | 18.89 [17.82, 20.47] |

Arm S's interval covers the predicted 19.53. At a threshold of a half it holds the same 3, and at a hundredth 4, 20.08 [19.14, 21.38]. Arm H's count moves from 2 to 5 over the same thresholds. So in arm S stays that last are about as rare as the balance says, on 3 events.

**The both-main-effects side's own mass is open.** Past the hold it takes 0.00296 of arm S's sweeps, 6.24 × 10⁶ of 2.11 × 10⁹, where the flat histogram's 0.017931 would make 3.78 × 10⁷. That ratio, 5.82 against 4.00, is 1.8 nats apart, and 0.00296 sits below all four of the flat histogram's chains (0.005–0.032). Under the balance's own rate, with stays exponential at first entry's mean and the open stay taken as closed, a total this small has a probability of about 0.02. Against it, the chains were selected on a1's side and stays are few, and a ratio of 5.82 would predict a return of 21.35, also inside the interval. This does not touch a1's conditional share, which leaves the side out by design, and a1's own unconditional share moves little: 0.659 on the flat histogram's mass, 0.670 on the chains'. What it leaves open is the balance's ratio, which rests on a both-side mass only the flat histogram has read.

#### What was not run, and why

- **Arm S past 5 × 10⁷, and arm H past 2.5 × 10⁷: no chain reached 10⁸.** The budget stopped there at 23.83 of 24 core-hours, with 44 arm S chains at 5 × 10⁷ and chain3's chains 4–7 at 4.5 × 10⁷. Every rule had its n by 2.5 × 10⁷ except A's eighth return, which came later, and what the sweeps past that added was the share's precision and two more long stays.
- **The plan's cost was stated at 45.6 core-hours before the chains,** and the budget order applied. It read low: arm S to 2.5 × 10⁷ cost 13.55 where 7.95 was planned, at 24,601 sweeps a second a thread against 41,900, and arm H 4.22 where 3.45 was planned, at 39,517 against 48,300, with other sessions' timed passes on the machine.
- **Past 2.5 × 10⁷, arm S ran at 54,036 sweeps a second a thread.** The machine was quieter, and the rebuild took in the efficiency pass from `main` at `d3fb924`. The two are not separated. That report has its levers bit for bit or, for lever 1, to a tolerance, so the chain's law is unchanged. Most of the second half's sweeps are on the later build, and the halves show no drift.
- **Parts 1–3 and the gates** took 1.77 core-hours on the ledger, plus about 0.4 for the share gate, of the 4 allowed.
- **The sampled sum's gate,** since nothing was sampled.
- Modes at each N, the fixed-M family, a second data seed, and moving `HittingReduction`, as the brief says.

#### Limits

- One fixture, one data seed, one flat-histogram run, wherever the share or the balance is quoted.
- 8 chains a start: a start's interval is about 1.5 nats wide, and the pooled one 0.6.
- Part 1's dictionaries hold 48 terms, so the decomposition and every class reading carry to N = 500 only as a method.
- Part 2 is one-sided. M is unread on four of the six starts, and the mass's noise on A and B makes that window about 19 nats wide.
- Part 3's traces are the earlier report's: the last 10⁴ sweeps of 44 crossings, at held and drawn τ.
- The long stays' threshold was chosen after the balance's rule failed, and they number 3.
- Touches have almost no check on the kernel gate. They rest on the brute-force test.
- No tally reads a1's exact structure, the flow out of the largest pre-entry basin, or how many touches follow the first return.

#### For the main session

**G1, the waiting-time bullet.** Its last two sentences, from "First entry is not always holding" on, are replaced by this bullet:

> *The wait to hold, measured* ([report](basin-hold-brief.md#reports); Friedman, N = 500, one flat-histogram run, one data seed). With τ drawn as in production, 8 chains from each of the six deep starts hold a1's basin, at half of a 10⁵-sweep window, in a mean of about 5.6 million sweeps, e^{15.54} [e^{15.25}, e^{15.84}], none censored, cv 1.03. That is first entry to within 0.04 nats on every start, at every window and level read, and held τ reads within 0.87 nats. So entry is the cost. Entries still touch (15 of 48 stay under a tenth of the next 10⁴ sweeps), but a touch lasts a median of about 30 sweeps and the chain is back. The crossing itself is short: a median of 3,336 sweeps from last holding both main effects to first entry. Past the hold, a1's share among sweeps not holding both main effects reads 0.6723 ± 0.0018, with no drift, against the flat histogram's 0.6712, so that share has a second reading. Both main effects' own mass does not: chains give it 0.3%, on 3 long stays, against the flat histogram's 1.8%. The two-well balance, read as the first return to any state holding both, is 6.3 nats under its e^{19.53}, because those returns are touches of about 80 sweeps; stays that last, read after the fact, come at e^{20.4} [e^{19.3}, e^{22.0}]. No chain hit a1's exact structure in the earlier traces, and no test confirmed a route reading across a real barrier.

**G1, "What a phase can do".** "so the starts seem to fall into one common region and leave it at one rate. That is an inference; no chain records the region" becomes:

> Chains now record it. Before first entry, family 1/1/3 takes a median share of 0.87–0.97 on four starts and 0.41 and 0.45 on chain0 and chain2, so one region by family is not supported. One basin, x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5, is the largest on all six, at 0.42–0.72.

**G1, "How the chains cross is unread"** becomes:

> *How the chains cross, read by basin* ([report](basin-hold-brief.md#reports)). Three doors, one hop from a1's basin, take 36 of 48 entries: 3 × x1·x2 + 2 × x3 + x4 + x5, x1 + 4 × x1·x2 + 2 × x3 + x4 + x5, and 4 × x1·x2 + x2 + 2 × x3 + x4 + x5. They are the last stops of an exchange of x1 and x2 for one x1·x2, from the basin above. No route holds all three. 2–5 of 12 gateways lie on a start's route, nearly all at its last door, and at most 6.4% of basin-changing moves do, so chains share a route's last step and not its shape. *The number of ways.* A term's whole measure on an edge's own term is refuted on all six starts on every route, by 0.14–17 nats past the rule's line. Also crediting each term held earlier reaches the waits on all six on the exact route, five on partial and four on aligned, with a spread of 12–23 nats where the waits spread 0.64. So it can be sufficient and does not predict. On enumerated fixtures a basin series with exact lumped flows itself misses R_eff by 5.2 nats at the 90th percentile off-cap, so no one basin route carries this. There, resistance inside basins has a median of 2 nats and falls with depth, and parallel routes between basins have 4.4–4.6 and grow with depth. That does not support "finding the terms sets the wait", and those 48-term dictionaries cannot read it at N = 500, where a class holds 10⁶.

**Candidate readings.** (f) and (g) are done. (c) gains a target: the both side's mass, 0.3% on chains against 1.8%. (h) gains Part 1's reading: on the enumerated fixtures there are about 2 nats inside basins to save, fewer at depth, and N = 500 is unread. (e) stands.

**SC-4 does not change in substance.** The polynomial expectation stays neither broken nor confirmed, since nothing here reads N. Its line on k = 4 costing e^{−20}, which the waiting-time report read against the rate to a1's basin, can now say that rate is a rate of holding and not only of touching: e^{−15.5} a sweep. If anything is added, the exchange above is a concrete case for its k on this fixture: two terms out and one in, with the order untimed. SC-5 is untouched.

**Next readings, for you to weigh:**
- *The both side's mass.* Chains would need about 7 × 10⁹ sweeps past hold for 10 long stays, about 36 core-hours at 54,000 sweeps a second a thread. Reading (c) or a second flat-histogram run may be cheaper.
- *A flow tally out of the largest pre-entry basin.* It would time the exchange, and say whether it is the slow step, on chains already this cheap. Nothing has measured that.

**Which lever, if any, stays the user's.** The fact it would rest on: a chain from any of these six deep starts needs about 5.6 million sweeps to hold a1's basin, drawn τ or held, and holding costs nothing past entry. That is about 2,790 times a 2,000-sweep phase, and about 2–4 minutes of one thread at the rates measured here, which were read under shifting load and are not a benchmark. The crossing that succeeds takes about 3,300 sweeps. Once held, a chain leaves for a long stay about once in 7 × 10⁸ sweeps, on 3 events, and agrees with the flat histogram on a1's side. So a lever has one thing to beat, the rate at which a crossing that commits starts, and three quarters of those crossings end in one exchange. All of this is on one fixture and one data seed.
