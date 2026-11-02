# The path profile between Friedman's real modes (chip brief)

Date: 2026-11-02. A measurement, tests only: it changes no code under `src/` and no default. Gap G1 of the [lineage-gaps register](lineage-gaps-register.md), candidate reading (a); SC-4 and SC-5 of [`structure-chain-lemmas.md`](../../../lean/Protolemmata/BASS/structure-chain-lemmas.md).

- *Why.* On S10's Friedman fixture (N = 500) the flat histogram puts 66–68% of the posterior in basin a1, while untempered chains and the reference arms hold basins at no more than 1.3%, about 45 nats lower. SC-5 proves that the route S → S ∪ T → T falls at most k·c, with c the Occam cost per column, about ½ ln N. So the crossing cost should be polynomial in N and exponential in k, where k counts the terms that do not cover their Occam cost when added. The [exact-gap sweep](spectral-gap-brief.md) cannot read this: its growth sits at the count cap, its fixed grid makes free-slot modes shallow almost by construction, and it reaches only k = 2 on one mode.
- *The question.* Between real modes on Friedman, what are k and the depth of the best route found? How much of the crossing cost is the barrier, and how much is the chance of proposing the right terms? And does k grow with N? The register's polynomial expectation fails only if it does.
- *Not in scope.*
  - Any move, kernel, schedule, option or default.
  - Choosing a lever, which is the user's, on this report.
  - G1's other candidate readings: the flat histogram's error against round trips, and the fixed-M family at wider caps.
  - Editing the register or the protolemmas. Put what they should say under "For the main session" in the report.

This section and the next three were written before any route was computed. The harness's choices, where this leaves them open, go in the report and its docstrings.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the register's G1, and SC-4 and SC-5, including SC-4's note on Yang, Wainwright & Jordan: near-duplicate knots are the expected hard part;
- [mode-exit-pricing-brief.md](mode-exit-pricing-brief.md), whose reports already read the neighbouring facts: no pair death below 57 nats (P1), and drift endpoints within 20 nats for the x1·x2 block in every chain (P3);
- [local-move-measurements-brief.md](local-move-measurements-brief.md) for F-M2's readings: every term load-bearing, births losing to about ½ ln N, a local knot step accepted about 0.95;
- [flat-histogram-warmstart-brief.md](flat-histogram-warmstart-brief.md) and [friedman-band-brief.md](friedman-band-brief.md) for the basins and their shares.

In the code:
- *Modes:* `tests/regression/BassFlatHistogramQualificationTests.cs` (`ArmConfig`, `RunArm`, seed 8017) and `BassMixing.cs` (`Basin`, `Describe`, `FixedStates`).
- *Scoring:* `BassTarget.Evaluate(fit, logTau)`, `BassTarget.LogStructurePrior`, and `BassMixing.cs`'s `Removals`, `With`, `Without`, `Stepped` and `Score`.
- *Proposals:* `BassTermProposal.cs` (`LogBirthRatio`, `LogDeathRatio`) and `BassSweep.cs` (what one sweep attempts).
- *Enumeration, for validation:* `BassKernelEnumeration.cs`, and `BassKernelSpectrum.cs` (`ReversibleChain`, `ModeDepths`).

## What is measured

**The modes.** Basins are keyed by `Basin`, the multiset of terms' variable sets.
- From one flat-histogram run of arm (i) at seed 8017, with retained draws, take each basin's two highest-posterior structures at b = 1, for: a1 (4 × x1·x2 + 2 × x3 + x4 + x5), a1 + x6·x8, a2 (4 × x1·x2 + x2 + 2 × x3 + x4 + x5), A (x1 + 2 × x1·x2 + x2 + 2 × x3 + x4 + x5) and B (x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5).
- From `FixedStates`, take the four untempered chains' structures (the states P1–P3 read).
- Store every exemplar in `path-profile-modes.csv` beside this brief, term by term with exact variable, sign and knot index, so that later readings need no rerun. If a basin has no draw in this run, say so and go on.

**The pairs,** each read in both directions: A–a1, B–a1, a2–a1, (a1 + x6·x8)–a1 and each chain's structure–a1. A–a1 is also read on the second exemplars, to show how far k depends on the exemplar. Signature differences only bound how many terms differ, since knots may differ too.

**π.** The log posterior of a structure at b = 1, prior included, from `BassTarget`. Two readings of τ:
- *held,* at unit information, SC-0's setting and the spectral sweep's. This is the primary reading;
- *profiled,* at each state's own best τ, closer to a chain whose τ moves every sweep.

**The route.** For a pair (S, T), let D be the symmetric difference and L the 2^|D| structures between them, each holding S ∩ T plus a subset of D. One step adds or deletes one term.
- *Maximin route.* The path through L whose lowest π is highest, by a widest-path search. This is exact within L. If |D| is too large to enumerate, use a beam search and say so.
- *Admissibility.* A state counts only if production would admit it: under the cap `MaxBasis`, and past the support, repeated-term and angle checks. Report how many states each check refused, and how far they change the route. Also read each pair with the cap lifted, if that needs no change under `src/`.
- *Alignment.* Exact-knot comparison overstates D, because a local knot step costs almost nothing. So read each pair twice. *Exact:* terms match only if identical. *Aligned:* match terms with the same variables and signs by least knot distance, slide them together one knot step at a time, recording π along the slide, and run the route on what remains. The slide goes first or last, whichever gives the higher minimum. Aligned is the primary reading, and it is a route, so it bounds rather than decides.

**The readings per route.**
- *Depth.* d_S = π(S) − min π and d_T = π(T) − min π along the route.
- *k.* The number of additions with Δπ < 0 between S and the route's minimum. Also report the number of deletions with Δπ < 0, and the highest M reached.
- *The profile.* At each step: the term, Δ log marginal, Δ log prior, Δπ, and c from SC-5's form.
- *Flow.* For each edge, ln c_e = min(π(x) + ln q(x → y), π(y) + ln q(y → x)), with q from the production proposal, move selection included.
- *Class flow.* The exact term's q understates the flow, since a neighbouring knot often serves as well. For each addition, also sum the flow over terms with the same variables and signs within ±25 knot positions per factor, enumerated. This is an approximation, and it is labelled as one.
- *The split.* How much of the route's total −ln c_e is the fall in π and how much is proposal cost, exact and by class.
- *Waiting time.* Σ_e exp(π(S) − ln c_e) attempts, converted to sweeps by what a sweep attempts, on the assumption that the structures in S's basin have like exits. Compare it with the sampling phase's length in `BassMixing`.

**Part 2: k against N.** On A–a1, B–a1 and the best chain's structure–a1, at N = 125, 250, 1000 and 2000, with fresh covariates and noise from the fixture's generator.
- Carry each exemplar over by knot quantile, polish it with P3's greedy step, and record whether its `Basin` survives. Check with a short chain at b = 1 whether it holds as a mode, and report its rank.
- Then read the same route quantities, aligned, held and profiled.
- This part is lower confidence than Part 1, because these modes are carried over rather than found. Say so wherever it is quoted.

## What must hold before any reading counts

- **The search.** Against brute force over all orderings on |D| ≤ 7: the same maximin value.
- **SC-5's inequality.** On every superset edge read, L(M + 1, R²′) ≥ L(M, R²) − c, using the form proved in `collapsedLogMarginal_superset_ge`. A violation stops the run and is reported.
- **Ground truth on enumerated fixtures.** On fixtures `BassKernelEnumeration` can hold, compute the exact maximin between pairs of local modes over the whole state graph, independently of L. Report how far L's route falls short of it, with alignment and without. That says how loose "one route" is, on small fixtures only.
- **The scorer.** Every exemplar's π from a fresh design equals its value in the retained draw, to rounding. Recheck the stored exemplars when the report is written, since other sessions are editing `src/`.
- **The flow.** On a small fixture, ln c_e from the harness equals the enumerated kernel's π(x)·K(x, y) for the same edge. Where the two differ by move selection, name the factor.

## How the outcomes read

Fixed now, per chain-held pair against a1, using the class flow and the better τ reading.
- **Waiting time at least 10× the sampling phase:** this route is too slow to have been taken, which is consistent with the chains staying put. It does not show that no faster route exists.
- **Waiting time at most 1/10 of the sampling phase:** the chains should have crossed this way, and they did not. That is the informative outcome. Check whether a chain's states reach a1, and what the flow calculation missed, before reporting anything else.
- **In between:** inconclusive, and reported as that.

On k, reported without choosing a lever:
- **k ≤ 2 on every pair, with proposal cost the larger share:** the barrier is not the main obstacle, and finding the terms is.
- **k ≥ 3 on some pair, with the fall in π the larger share:** one-term moves meet SC-5's exponential in k.
- **k that moves between exemplars by more than 1:** k is not yet a stable reading of the basin, and the report says so before anything else.
- **Across N (Part 2):** k falling, flat or growing. Growing is the case that would break the register's polynomial expectation, so it needs the most care: first check that the carried-over structures are still modes. Depth at free slots cannot grow past SC-5's bound, so a depth that seems to is a bug or a closed cap.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Other sessions are editing here,** including `BassTarget.cs`, `BassFlatHistogramRun.cs` and files under `src/maths/samplers/`. Leave their files alone. If a build breaks in a file you did not touch, wait and retry, and report it if it persists. Record `git rev-parse HEAD` and `git status --short -- src/` in the report.
- **Build and test into your own artifacts root,** set inline on each command: `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-pathprofile\`. Do not run `scripts/fact-harness.ps1`. Delete that folder when you finish.
- **The machine is shared and has jammed under load.** Run long or heavy work at below-normal priority so children inherit it, with `-m:2` on builds and `-- xUnit.MaxParallelThreads=4` on tests.
- **Runs can die.** The harness appends each row as it is computed and skips rows already present on a rerun. CSVs go beside this brief as `path-profile-*.csv`. Heavy tests sit behind `BASS_QUALIFICATION` = 1.
- **Commit by pathspec** (`git commit -- paths`, after `git add` for new files), with targeted commits and nothing left staged. End commit messages with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Tooling here.** The Grep tool is denied by a repo hook, so use `rg` through Bash. The default `python` has no numpy; use `C:/Users/azrie/PDenv/PyenvPython/pyenv-win/versions/3.12.6/python.exe`, and put patterns in script files, since the `python` shim breaks on `|` in arguments.
- **New writing takes the module's symbols** (b, τ, ℓ, M, N), with an author's symbol named once where a source needs it.

## Review before reporting

Each step ends with a fresh-eyes review of its own code and readings, for correctness and for coverage gaps, fixed before the next step. In particular:
- recompute the headline readings from the CSVs, apart from the harness;
- read the profile of at least one route term by term, and check that its k agrees with what its Δπ column shows;
- name each reading that one mode, one data seed or one surface carries.

## Report

Append it under "Reports" below, newest first: commits, HEAD and the dirty `src/` files, run times, the validation results, then in order:
- the modes found and any basin missing;
- a table per pair (k, depths, highest M, refused states, flow, split, waiting time; aligned and exact; held and profiled);
- one route's full profile;
- Part 2's readings;
- the reading of each outcome above;
- the limits;
- "For the main session": what the register's G1 and SC-4 should now say.

## Reports

### Measured — 2026-11-02

**k is not yet a stable reading of the basin, by the rule above.** On A–a1, the aligned route at held τ gives k from A of 2 on the brief's first and second exemplars, and 4 on an added sign-distinct exemplar of A. The exact route gives k from a1 of 5, 5, 4, 7 and 5 across the five exemplar pairs. On one pair of exemplars, k from A is 2 aligned and 4 exact. The cause shows in the profiles: k counts swaps, an addition followed by the deletion of its near twin, which come and go with the knots. Net k, the terms held at the route's lowest state beyond the start's, is 2 from either end on every A–a1 row. Depth moves by about 1 nat across them. So read k below with net k beside it.

#### Commits and state

- `4fb7dc5`: the harness, `tests/regression/BassPathProfile.cs` and `BassPathProfileTests.cs`, with a few private test members opened to it and a `seed` argument on `Friedman`. `dc16dd8`: the rest of the harness, the CSVs and this report.
- HEAD was `bbcf2ee` at the start and `2858598` when this was written. Other sessions landed `a772fd9` (the weighted core) and the truncated-draws commits meanwhile. The dirty `src/` files now: `src/maths/distributions/HalfCauchy.cs`, `src/maths/distributions/Mixture.cs`, `src/maths/distributions/TailInversion.cs`, `src/maths/distributions/Weibull.cs`. Nothing under `src/` was changed here.
- The modes ran on a build made at the start, and Parts 1 and 2 on a later one. So the scorer was rechecked on a build of `2858598` with those dirty files: it holds, below.
- Run times, at below-normal priority on a machine shared with two other test hosts: arm (i) 4,882 s, against about 30 minutes alone; the ground truth under a minute; Part 1 9 m 33 s; Part 2 10 m 38 s, and 41 s for the redraw rows.

#### Validation

All of it held, and no reading stands on a failed check.

- **The search.** Against brute force over every ordering, on the 162 routes read: L itself on the 50 with |D| ≤ 7, and a subcube of up to seven of L's terms on the rest. The widest path equals a separate fixed-point relaxation exactly. It exceeds the best ordering on 0 of them. On 400 random cubes it matches the relaxation, and 32 are beaten only by a non-monotone route, as random π allows. None of the 324 route rows here is non-monotone.
- **SC-5's inequality.** 32,577,747 superset edges read at the held τ across every N, and each priced route's edges at the lower state's best τ. No violation; the least slack is 0.000000, so the bound is tight where a term adds nothing. A violation would have stopped the run.
- **Ground truth,** on the enumerated grid (N = 18), 167 pairs of local modes over 16 data seeds and caps. The third fixture, the five-value grid with main effects alone, has one local mode per seed, so no pairs. It says "one route" is loose by a few nats on these small fixtures, and nothing about Friedman.

| fixture | pairs | whole graph above birth/death (median, max) | exact route short of the whole graph | aligned route short of it | exact route short of the birth/death graph | reaches birth/death |
|---|---|---|---|---|---|---|
| M ≤ 3 | 127 | 1.97, 3.95 | 2.86, 9.64 | 2.69, 9.64 | 0.00, 7.41 | 65 |
| M ≤ 4 | 40 | 0.80, 2.80 | 2.06, 6.18 | 1.80, 9.64 | 1.05, 5.15 | 18 |

  With the cap lifted, the M ≤ 3 exact route's shortfall falls to 2.53 median, and M ≤ 4's does not move. Most of what the routes miss against the whole graph is the redraw edge, which L does not use.
- **The scorer.** All 19 exemplars from a fresh design: the log marginal at the retained draw's τ against that draw's ℓ₁, worst 1.2·10⁻⁶; against the stored π, worst 5·10⁻⁷. Both are the rows' six decimals. Rechecked when this was written, on `2858598`.
- **The flow.** On the 2,263-state fixture, 6,464 birth edges and 29,600 redraw edges: ln c_e equals the enumerated kernel's ln π(x)K(x, y) to 0.0 in doubles. Both hold the move's selection, so no factor is left over.
- **Recomputed apart from the harness.** k, the depths, the lowest π, the highest M and both log waiting times, from `path-profile-steps.csv` by a separate script: 324 direction rows, no differences. One route is read term by term below.

#### Choices the brief left open

Those that shape the harness are in the test class's docstring too.

- Exemplars rank by π at held τ. I added each basin's highest structure whose signs differ from the top one's (rank d), since the brief's rule on k needs exemplars that differ. It equals rank 2 except on A.
- The maximin value does not decide the route past its lowest state. The route takes the fewest steps, then the least Σ_e exp(−min π). `k_least` is the least k over all maximin routes.
- A third route, **partial**, beside aligned and exact: a matched pair goes to the route as an exchange when that raises the lowest π, greedily, farthest pair first. A slide takes its pairs together however far apart, and that can cost more than exchanging the far ones. Partial is at least as high as both on all eight pairs in the tables below, so it is the best route found. Aligned stays the primary reading.
- The split is read at the costliest edge, whose cost is the log waiting time to within ln(edges). The Σ over all edges is in the CSV. Slide edges swamp it on aligned routes, and it goes negative where a route rises above its start, so it is not read here.
- A slide's edges are D01's redraw at one knot's distance, 1/(2n_v) each, since S10's configuration has the local steps off. Their share of the waiting time is its own column, and on the chain-held pairs from the chain the route's edges dominate.
- The class flow applies to additions in the direction read. A deletion has one term to delete.
- Part 2 reads fresh cases from another seed of the generator, and also the generator's own stream, which truncates or extends S10's cases. At N = 125 and 250 the stream's cases are ones the exemplars were found on.
- "Holds" in Part 2 means at least half of a 2,000-sweep chain stays in the basin and the start ranks at or above that chain's median.
- *The one the brief's rule forced.* A hinge at its variable's lowest knot and its mirror at the highest are one linear column. D06 refuses a state holding both, which closes every superset above it. On those lattices D06 refuses exactly a quarter of the states, the ones holding both, 31,360 in all, and nothing else refuses anything. The aligned routes then delete a load-bearing term first, and read 46 to 1,029 nats deep. The brief reads such a depth as a bug or a closed cap, and it is the closed gate. `aligned+redraw` and `partial+redraw` take such a pair by one redraw, the move production makes it by. The plain rows stay in the CSV. No Part 1 pair holds one.

#### The modes

Every basin has draws in this run. a1's share is 65.9%, beside the register's 66–68%. Each exemplar was drawn once: 477,632 of 500,000 draws are distinct structures, so the top of a basin is an extreme order statistic of its draws.

basin | exemplar | M | π held | π profiled | ln τ best | share at b = 1 | draws in the basin | structure, each hinge's sign and knot index
--- | --- | --- | --- | --- | --- | --- | --- | ---
a1 | 1 | 8 | -883.02 | -882.19 | -7.05 | 65.914% | 23288 | x1+3·x2+1 + x1+110·x2+93 + x1+185·x2+321 + x1+342·x2+170 + x3−383 + x3+205 + x4+6 + x5+9
a1 | 2 | 8 | -883.04 | -882.22 | -7.05 | 65.914% | 23288 | x1+0·x2+2 + x1+99·x2+97 + x1+198·x2+317 + x1+335·x2+161 + x3−382 + x3+202 + x4−486 + x5+9
a1 + x6·x8 | 1 | 9 | -894.24 | -893.55 | -6.99 | 10.593% | 3559 | x1+4·x2+7 + x1+102·x2+92 + x1+197·x2+321 + x1+339·x2+181 + x3−385 + x3+206 + x4−499 + x5−494 + x6+208·x8−135
a1 + x6·x8 | 2 | 9 | -894.38 | -893.69 | -6.99 | 10.593% | 3559 | x1+1·x2+5 + x1+99·x2+83 + x1+199·x2+310 + x1+342·x2+146 + x3−377 + x3+196 + x4−498 + x5+20 + x6+172·x8−164
a2 | 1 | 9 | -897.83 | -897.44 | -6.94 | 0.199% | 112 | x1+3·x2+5 + x1+132·x2+341 + x1+200·x2+137 + x1+335·x2+181 + x2−86 + x3−393 + x3+206 + x4−491 + x5−483
a2 | 2 | 9 | -898.00 | -897.61 | -6.94 | 0.199% | 112 | x1+3·x2+1 + x1+136·x2+126 + x1+201·x2+321 + x1+338·x2+160 + x2+177 + x3−356 + x3+200 + x4−494 + x5+9
A | 1 | 8 | -884.89 | -884.52 | -6.97 | 0.007% | 7261 | x1−370 + x1−360·x2−408 + x1+203·x2+228 + x2−441 + x3−392 + x3+203 + x4−496 + x5+11
A | 2 | 8 | -885.18 | -884.81 | -6.97 | 0.007% | 7261 | x1−370 + x1−360·x2−408 + x1+203·x2+216 + x2−435 + x3−392 + x3+201 + x4−496 + x5+11
A | d | 8 | -885.18 | -884.81 | -6.97 | 0.007% | 7261 | x1−370 + x1−360·x2−413 + x1+203·x2+228 + x2−441 + x3−390 + x3+197 + x4−475 + x5−489
B | 1 | 9 | -889.29 | -888.89 | -6.94 | 1.294% | 965 | x1−475 + x1+6·x2−336 + x1+193·x2+129 + x1+340·x2+186 + x2+326 + x3−378 + x3+182 + x4+5 + x5−491
B | 2 | 9 | -889.51 | -889.12 | -6.94 | 1.294% | 965 | x1−475 + x1+7·x2−336 + x1+188·x2+123 + x1+340·x2+197 + x2+326 + x3−399 + x3+199 + x4+9 + x5+19

chain | M | π held | basin | rank in its phase | phase median | gap below the best median | read at
--- | --- | --- | --- | --- | --- | --- | ---
0 | 10 | -933.98 | x1 + 5 × x1·x2 + 2 × x3 + x4 + x5 | 73.2% | -933.91 | 41.52 | the burned state
1 (best) | 9 | -912.06 | x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 8.0% | -892.39 | 0.00 | the burned state
2 | 11 | -925.84 | 2 × x1 + 3 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5 | 56.7% | -926.19 | 33.80 | the burned state
3 | 9 | -917.22 | x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 0.0% | -916.01 | 23.63 | the phase's most frequent structure

- Chain 0: x1−291 + x1−303·x2−242 + x1−450·x2−456 + x1−154·x2+247 + x1+24·x2−191 + x1+333·x2+296 + x3−385 + x3+200 + x4−493 + x5+20
- Chain 1: x1−447 + x1−377·x2−450 + x1+60·x2−304 + x1+265·x2+192 + x2−479 + x3−393 + x3+178 + x4−476 + x5+4
- Chain 2: x1−430 + x1−319·x2−226 + x1−93·x2+245 + x1+278 + x1+279·x2+272 + x2−158 + x2−434 + x3−354 + x3+196 + x4−476 + x5+12
- Chain 3: x1−421 + x1−151·x2+125 + x1+104·x2−139 + x1+305·x2+259 + x2−455 + x3−159 + x3+341 + x4−476 + x5+10

- *a1's share is not a higher peak.* Its top structure stands 1.9 nats above A's and 6.3 above B's at held τ, while its share is 9,400 times A's and 51 times B's. So on these draws the mass is breadth: a1 has many high structures, and A's top is a rare one.
- *The chains' states are not their basins' tops.* They sit 29–51 nats below a1's top and 23–45 below B's. Chain 1's burned state ranks 8% in its own phase, 20 nats under that phase's median, so its route starts below the mode the chain then climbs to. Chain 3 is read at its phase's most frequent structure, as P1–P3 read it.

#### The pairs

Term names here use the fixture's x1…x10. `path-profile-steps.csv` holds them from x0, as `BassMixing.Key` writes them. Each row is one route read from one end. *k (net, least)* is k, net k and the least k over maximin routes. *Fall + proposal* is the costliest edge's π(from) − π(side) and its −ln q. *Route's, slide's* are the class waiting time's two parts. *Superset d_from* is the depth of S → S ∪ T → T on the same lattice, each half at its best order.

Every Part 1 route row has 0 refused states, by the cap, the support, repeated terms, the angle or the exact fit. The highest M on any row is 12 against a cap of 499. That cap is N − 1, the design's own bound, so lifting it changes nothing here. Partial's profiled rows are in the CSV.

**A – a1**

route | τ | from | \|D\| | slide | d_from | d_to | k (net, least) | deletions | top M | fall + proposal, exact | fall + proposal, class | ln wait exact | ln wait class | route's, slide's | superset d_from
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
aligned | held | A | 8 | last 124 | 37.3 | 39.2 | 2 (2, 2) | 0 | 10 | 25.6 + 19.0 | 25.6 + 19.0 | 44.6 | 44.6 | 44.6, 32.2 | 66.8
aligned | held | a1 | 8 | last 124 | 39.2 | 37.3 | 4 (2, 4) | 0 | 10 | 27.5 + 19.0 | 39.2 + 3.4 | 46.5 | 42.6 | 42.6, 34.0 | 68.7
aligned | profiled | A | 8 | last 124 | 37.9 | 40.2 | 2 (2, 2) | 0 | 10 | 26.0 + 19.0 | 26.0 + 19.0 | 45.1 | 45.0 | 45.0, 32.2 | 67.8
aligned | profiled | a1 | 8 | last 124 | 40.2 | 37.9 | 4 (2, 4) | 0 | 10 | 28.3 + 19.0 | 40.2 + 3.4 | 47.4 | 43.6 | 43.6, 34.6 | 70.2
partial | held | A | 12 | last 120 | 37.2 | 39.1 | 4 (2, 4) | 0 | 10 | 25.7 + 19.0 | 25.7 + 19.0 | 44.7 | 44.7 | 44.7, 32.3 | 92.1
partial | held | a1 | 12 | last 120 | 39.1 | 37.2 | 4 (2, 4) | 0 | 10 | 27.6 + 19.0 | 39.1 + 3.4 | 46.6 | 42.5 | 42.5, 34.2 | 93.9
exact | held | A | 16 | — | 37.2 | 39.1 | 4 (2, 4) | 0 | 10 | 25.7 + 19.0 | 25.7 + 19.0 | 44.7 | 44.7 | — | 111.4
exact | held | a1 | 16 | — | 39.1 | 37.2 | 5 (2, 5) | 1 | 10 | 27.6 + 19.0 | 39.1 + 3.4 | 46.6 | 42.5 | — | 113.2
exact | profiled | A | 16 | — | 37.8 | 40.1 | 4 (2, 4) | 0 | 10 | 26.1 + 19.0 | 26.1 + 19.0 | 45.2 | 45.2 | — | 112.8
exact | profiled | a1 | 16 | — | 40.1 | 37.8 | 5 (2, 5) | 1 | 10 | 28.5 + 19.0 | 40.1 + 3.4 | 47.5 | 43.5 | — | 115.2

**B – a1**

route | τ | from | \|D\| | slide | d_from | d_to | k (net, least) | deletions | top M | fall + proposal, exact | fall + proposal, class | ln wait exact | ln wait class | route's, slide's | superset d_from
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
aligned | held | B | 7 | first 166 | 30.3 | 36.6 | 1 (1, 1) | 0 | 10 | 28.3 + 18.8 | 28.3 + 11.2 | 47.2 | 40.1 | 40.1, 23.7 | 53.3
aligned | held | a1 | 7 | first 166 | 36.6 | 30.3 | 4 (2, 4) | 1 | 10 | 34.6 + 18.8 | 34.6 + 18.8 | 53.4 | 53.4 | 53.4, 29.9 | 59.6
aligned | profiled | B | 7 | first 166 | 30.9 | 37.6 | 1 (1, 1) | 0 | 10 | 28.8 + 18.8 | 28.8 + 11.2 | 47.6 | 40.6 | 40.6, 23.9 | 54.2
aligned | profiled | a1 | 7 | first 166 | 37.6 | 30.9 | 4 (2, 4) | 1 | 10 | 35.5 + 18.8 | 35.5 + 18.8 | 54.3 | 54.3 | 54.3, 30.6 | 60.9
partial | held | B | 9 | last 47 | 27.9 | 34.2 | 3 (1, 3) | 0 | 10 | 24.9 + 18.8 | 24.9 + 11.2 | 43.7 | 36.6 | 36.6, 7.7 | 71.9
partial | held | a1 | 9 | last 47 | 34.2 | 27.9 | 3 (2, 3) | 0 | 10 | 31.1 + 18.8 | 31.1 + 18.8 | 50.0 | 50.0 | 50.0, 14.0 | 78.1
exact | held | B | 17 | — | 27.9 | 34.2 | 3 (1, 3) | 0 | 10 | 24.9 + 18.8 | 24.9 + 11.2 | 43.7 | 36.6 | — | 129.6
exact | held | a1 | 17 | — | 34.2 | 27.9 | 7 (2, 7) | 0 | 10 | 31.1 + 18.8 | 31.1 + 18.8 | 50.0 | 50.0 | — | 135.9
exact | profiled | B | 17 | — | 28.3 | 35.0 | 3 (1, 3) | 0 | 10 | 25.2 + 18.8 | 25.2 + 11.2 | 44.0 | 37.0 | — | 131.2
exact | profiled | a1 | 17 | — | 35.0 | 28.3 | 7 (2, 7) | 0 | 10 | 31.9 + 18.8 | 31.9 + 18.8 | 50.7 | 50.7 | — | 137.9

**a2 – a1**

route | τ | from | \|D\| | slide | d_from | d_to | k (net, least) | deletions | top M | fall + proposal, exact | fall + proposal, class | ln wait exact | ln wait class | route's, slide's | superset d_from
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
aligned | held | a2 | 5 | first 240 | 0.0 | 14.8 | 0 (0, 0) | 0 | 9 | 0.0 + 10.9 | 0.0 + 10.9 | 15.3 | 15.3 | 2.6, 15.3 | 23.5
aligned | held | a1 | 5 | first 240 | 14.8 | 0.0 | 3 (1, 3) | 0 | 9 | 14.8 + 10.9 | 14.8 + 10.9 | 30.1 | 30.1 | 17.0, 30.1 | 38.3
aligned | profiled | a2 | 5 | first 240 | 0.0 | 15.3 | 0 (0, 0) | 0 | 9 | 0.0 + 10.9 | 0.0 + 10.9 | 15.3 | 15.3 | 2.5, 15.3 | 24.1
aligned | profiled | a1 | 5 | first 240 | 15.3 | 0.0 | 3 (1, 3) | 0 | 9 | 15.3 + 10.9 | 15.3 + 10.9 | 30.5 | 30.5 | 17.4, 30.5 | 39.4
partial | held | a2 | 5 | first 240 | 0.0 | 14.8 | 0 (0, 0) | 0 | 9 | 0.0 + 10.9 | 0.0 + 10.9 | 15.3 | 15.3 | 2.6, 15.3 | 23.5
partial | held | a1 | 5 | first 240 | 14.8 | 0.0 | 3 (1, 3) | 0 | 9 | 14.8 + 10.9 | 14.8 + 10.9 | 30.1 | 30.1 | 17.0, 30.1 | 38.3
exact | held | a2 | 17 | — | 6.5 | 21.3 | 6 (0, 6) | 0 | 9 | 6.5 + 3.3 | 6.5 + 3.3 | 11.6 | 10.9 | — | 130.8
exact | held | a1 | 17 | — | 21.3 | 6.5 | 3 (1, 3) | 0 | 9 | 21.3 + 3.3 | 21.3 + 3.3 | 26.4 | 25.7 | — | 145.6
exact | profiled | a2 | 17 | — | 6.5 | 21.7 | 6 (0, 6) | 0 | 9 | 6.5 + 3.3 | 6.5 + 3.3 | 11.6 | 10.9 | — | 132.5
exact | profiled | a1 | 17 | — | 21.7 | 6.5 | 3 (1, 3) | 0 | 9 | 21.7 + 3.3 | 21.7 + 3.3 | 26.8 | 26.1 | — | 147.7

**a1 + x6·x8 – a1**

route | τ | from | \|D\| | slide | d_from | d_to | k (net, least) | deletions | top M | fall + proposal, exact | fall + proposal, class | ln wait exact | ln wait class | route's, slide's | superset d_from
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
aligned | held | a1+x | 5 | first 45 | 1.5 | 12.7 | 1 (0, 1) | 0 | 9 | 0.0 + 10.9 | 0.0 + 10.9 | 14.5 | 14.5 | 8.7, 14.5 | 25.4
aligned | held | a1 | 5 | first 45 | 12.7 | 1.5 | 2 (1, 2) | 0 | 9 | 11.3 + 10.9 | 11.3 + 10.9 | 25.7 | 25.7 | 16.5, 25.7 | 36.6
aligned | profiled | a1+x | 5 | first 45 | 1.7 | 13.1 | 1 (0, 1) | 0 | 9 | 0.0 + 10.9 | 0.0 + 10.9 | 14.5 | 14.5 | 8.6, 14.5 | 26.1
aligned | profiled | a1 | 5 | first 45 | 13.1 | 1.7 | 2 (1, 2) | 0 | 9 | 11.4 + 10.9 | 11.4 + 10.9 | 25.8 | 25.8 | 16.9, 25.8 | 37.5
partial | held | a1+x | 7 | first 43 | 1.5 | 12.7 | 1 (0, 1) | 0 | 9 | 0.0 + 10.9 | 0.0 + 10.9 | 14.4 | 14.4 | 8.7, 14.4 | 38.1
partial | held | a1 | 7 | first 43 | 12.7 | 1.5 | 3 (1, 3) | 0 | 9 | 11.3 + 10.9 | 11.3 + 10.9 | 25.7 | 25.7 | 17.0, 25.7 | 49.4
exact | held | a1+x | 17 | — | 9.7 | 20.9 | 3 (0, 3) | 0 | 9 | 8.1 + 19.1 | 8.1 + 19.1 | 27.9 | 27.2 | — | 134.7
exact | held | a1 | 17 | — | 20.9 | 9.7 | 5 (1, 5) | 1 | 9 | 19.3 + 19.1 | 19.3 + 19.1 | 39.2 | 38.5 | — | 145.9
exact | profiled | a1+x | 17 | — | 10.0 | 21.3 | 3 (0, 3) | 0 | 9 | 8.5 + 19.1 | 8.5 + 19.1 | 28.3 | 27.7 | — | 136.6
exact | profiled | a1 | 17 | — | 21.3 | 10.0 | 5 (1, 5) | 1 | 9 | 19.9 + 19.1 | 19.9 + 19.1 | 39.7 | 39.0 | — | 147.9

**chain 0 – a1**

route | τ | from | \|D\| | slide | d_from | d_to | k (net, least) | deletions | top M | fall + proposal, exact | fall + proposal, class | ln wait exact | ln wait class | route's, slide's | superset d_from
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
aligned | held | chain | 10 | first 153 | 25.4 | 76.4 | 2 (2, 2) | 0 | 12 | 24.3 + 18.4 | 24.3 + 11.8 | 42.9 | 36.6 | 36.2, 35.3 | 54.7
aligned | held | a1 | 10 | first 153 | 76.4 | 25.4 | 6 (4, 5) | 1 | 12 | 75.2 + 18.4 | 75.2 + 18.4 | 93.9 | 93.9 | 93.9, 86.3 | 105.7
aligned | profiled | chain | 10 | first 153 | 25.6 | 77.7 | 2 (2, 2) | 0 | 12 | 24.6 + 18.4 | 24.6 + 11.8 | 43.2 | 36.9 | 36.5, 35.7 | 55.3
aligned | profiled | a1 | 10 | first 153 | 77.7 | 25.6 | 6 (4, 5) | 1 | 12 | 76.6 + 18.4 | 76.6 + 18.4 | 95.3 | 95.3 | 95.3, 87.7 | 107.3
partial | held | chain | 12 | first 148 | 25.2 | 76.2 | 2 (2, 2) | 0 | 12 | 23.9 + 18.4 | 23.9 + 11.8 | 42.5 | 36.2 | 35.9, 34.9 | 67.2
partial | held | a1 | 12 | first 148 | 76.2 | 25.2 | 7 (4, 6) | 1 | 12 | 74.9 + 18.4 | 74.9 + 18.4 | 93.5 | 93.5 | 93.5, 85.8 | 118.2
exact | held | chain | 18 | — | 25.2 | 76.2 | 5 (2, 5) | 1 | 12 | 23.9 + 18.4 | 22.2 + 18.6 | 42.7 | 40.9 | — | 112.7
exact | held | a1 | 18 | — | 76.2 | 25.2 | 7 (4, 6) | 1 | 12 | 74.9 + 18.4 | 74.9 + 18.4 | 93.7 | 93.5 | — | 163.7
exact | profiled | chain | 18 | — | 25.4 | 77.5 | 5 (2, 5) | 1 | 12 | 24.2 + 18.4 | 22.6 + 18.6 | 43.0 | 41.2 | — | 113.8
exact | profiled | a1 | 18 | — | 77.5 | 25.4 | 7 (4, 6) | 1 | 12 | 76.3 + 18.4 | 76.3 + 18.4 | 95.1 | 94.9 | — | 165.8

**chain 1 – a1**

route | τ | from | \|D\| | slide | d_from | d_to | k (net, least) | deletions | top M | fall + proposal, exact | fall + proposal, class | ln wait exact | ln wait class | route's, slide's | superset d_from
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
aligned | held | chain | 9 | last 141 | 17.3 | 46.4 | 3 (0, 3) | 1 | 9 | 17.0 + 19.1 | 17.0 + 11.2 | 36.1 | 28.9 | 28.9, 6.1 | 67.5
aligned | held | a1 | 9 | last 141 | 46.4 | 17.3 | 2 (1, 2) | 0 | 9 | 46.1 + 19.1 | 46.1 + 19.1 | 65.1 | 65.1 | 65.1, 35.2 | 96.5
aligned | profiled | chain | 9 | last 141 | 17.5 | 47.4 | 3 (0, 3) | 1 | 9 | 17.2 + 19.1 | 17.2 + 11.2 | 36.3 | 29.1 | 29.1, 5.7 | 68.3
aligned | profiled | a1 | 9 | last 141 | 47.4 | 17.5 | 2 (1, 2) | 0 | 9 | 47.2 + 19.1 | 47.2 + 19.1 | 66.3 | 66.3 | 66.3, 35.7 | 98.3
partial | held | chain | 15 | first 5 | 15.8 | 44.9 | 4 (0, 4) | 1 | 9 | 15.7 + 19.1 | 15.7 + 11.2 | 34.8 | 27.6 | 27.6, 11.8 | 91.7
partial | held | a1 | 15 | first 5 | 44.9 | 15.8 | 4 (1, 3) | 0 | 9 | 44.8 + 19.1 | 44.8 + 19.1 | 63.9 | 63.9 | 63.9, 40.8 | 120.8
exact | held | chain | 17 | — | 15.8 | 44.9 | 5 (0, 5) | 1 | 9 | 15.7 + 19.1 | 15.7 + 11.2 | 34.8 | 27.6 | — | 104.3
exact | held | a1 | 17 | — | 44.9 | 15.8 | 4 (1, 3) | 0 | 9 | 44.8 + 19.1 | 44.8 + 19.1 | 63.9 | 63.9 | — | 133.3
exact | profiled | chain | 17 | — | 16.0 | 45.9 | 5 (0, 5) | 1 | 9 | 15.9 + 19.1 | 15.9 + 11.2 | 35.0 | 27.8 | — | 105.4
exact | profiled | a1 | 17 | — | 45.9 | 16.0 | 4 (1, 3) | 0 | 9 | 45.9 + 19.1 | 45.9 + 19.1 | 65.0 | 65.0 | — | 135.4

**chain 2 – a1**

route | τ | from | \|D\| | slide | d_from | d_to | k (net, least) | deletions | top M | fall + proposal, exact | fall + proposal, class | ln wait exact | ln wait class | route's, slide's | superset d_from
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
aligned | held | chain | 11 | last 184 | 16.5 | 59.3 | 2 (0, 2) | 1 | 12 | 14.5 + 19.2 | 14.5 + 19.2 | 33.7 | 33.7 | 33.7, -8.0 | 69.2
aligned | held | a1 | 11 | last 184 | 59.3 | 16.5 | 5 (3, 4) | 1 | 12 | 57.3 + 19.2 | 53.5 + 19.2 | 76.5 | 72.8 | 72.8, 34.8 | 112.0
aligned | profiled | chain | 11 | last 184 | 16.8 | 60.9 | 2 (0, 2) | 1 | 12 | 14.8 + 19.2 | 14.8 + 19.2 | 34.0 | 34.0 | 34.0, -8.7 | 69.8
aligned | profiled | a1 | 11 | last 184 | 60.9 | 16.8 | 5 (3, 4) | 1 | 12 | 58.9 + 19.2 | 55.0 + 19.2 | 78.1 | 74.3 | 74.3, 35.4 | 113.9
partial | held | chain | 15 | first 38 | 14.5 | 57.3 | 2 (0, 2) | 1 | 12 | 12.5 + 19.2 | 12.5 + 19.2 | 31.7 | 31.7 | 31.7, 13.3 | 89.0
partial | held | a1 | 15 | first 38 | 57.3 | 14.5 | 7 (3, 6) | 2 | 12 | 55.4 + 19.2 | 51.5 + 19.2 | 74.6 | 70.7 | 70.7, 56.1 | 131.9
exact | held | chain | 19 | — | 14.5 | 57.3 | 4 (0, 4) | 1 | 12 | 12.5 + 19.2 | 12.5 + 19.2 | 31.7 | 31.7 | — | 113.7
exact | held | a1 | 19 | — | 57.3 | 14.5 | 7 (3, 6) | 2 | 12 | 55.4 + 19.2 | 51.5 + 19.2 | 74.6 | 70.7 | — | 156.5
exact | profiled | chain | 19 | — | 14.7 | 58.8 | 4 (0, 4) | 1 | 12 | 12.8 + 19.2 | 12.8 + 19.2 | 32.0 | 32.0 | — | 114.6
exact | profiled | a1 | 19 | — | 58.8 | 14.7 | 7 (3, 6) | 2 | 12 | 56.9 + 19.2 | 52.9 + 19.2 | 76.1 | 72.2 | — | 158.7

**chain 3 – a1**

route | τ | from | \|D\| | slide | d_from | d_to | k (net, least) | deletions | top M | fall + proposal, exact | fall + proposal, class | ln wait exact | ln wait class | route's, slide's | superset d_from
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
aligned | held | chain | 9 | last 487 | 25.8 | 60.0 | 2 (0, 2) | 1 | 10 | 25.8 + 19.0 | 25.8 + 19.0 | 45.5 | 44.8 | 44.8, 25.8 | 63.1
aligned | held | a1 | 9 | last 487 | 60.0 | 25.8 | 3 (1, 3) | 1 | 10 | 60.0 + 19.0 | 60.0 + 19.0 | 79.7 | 79.0 | 79.0, 60.0 | 97.3
aligned | profiled | chain | 9 | last 487 | 26.2 | 61.4 | 2 (0, 2) | 1 | 10 | 26.2 + 19.0 | 26.2 + 19.0 | 45.9 | 45.2 | 45.2, 25.8 | 63.9
aligned | profiled | a1 | 9 | last 487 | 61.4 | 26.2 | 3 (1, 3) | 1 | 10 | 61.4 + 19.0 | 61.4 + 19.0 | 81.2 | 80.5 | 80.5, 61.0 | 99.1
partial | held | chain | 13 | first 360 | 21.7 | 56.0 | 2 (0, 2) | 1 | 10 | 21.7 + 19.0 | 21.7 + 19.0 | 41.5 | 40.8 | 40.8, 33.1 | 77.4
partial | held | a1 | 13 | first 360 | 56.0 | 21.7 | 5 (1, 4) | 2 | 10 | 56.0 + 19.0 | 56.0 + 19.0 | 75.7 | 75.0 | 75.0, 67.3 | 111.6
exact | held | chain | 17 | — | 21.7 | 56.0 | 4 (0, 4) | 1 | 11 | 21.7 + 19.0 | 21.7 + 19.0 | 41.5 | 40.8 | — | 101.3
exact | held | a1 | 17 | — | 56.0 | 21.7 | 5 (1, 4) | 2 | 11 | 56.0 + 19.0 | 56.0 + 19.0 | 75.7 | 75.0 | — | 135.5
exact | profiled | chain | 17 | — | 22.1 | 57.3 | 4 (0, 4) | 1 | 11 | 22.1 + 19.0 | 22.1 + 19.0 | 41.8 | 41.1 | — | 102.3
exact | profiled | a1 | 17 | — | 57.3 | 22.1 | 5 (1, 4) | 2 | 11 | 57.3 + 19.0 | 57.3 + 19.0 | 77.0 | 76.4 | — | 137.5

A–a1 across exemplars, held τ:

exemplars (A/a1) | route | k from A (net) | k from a1 (net) | d from A | d from a1 | ln wait from A | ln wait from a1
--- | --- | --- | --- | --- | --- | --- | ---
1/1 | aligned | 2 (2) | 4 (2) | 37.3 | 39.2 | 44.6 | 42.6
1/1 | partial | 4 (2) | 4 (2) | 37.2 | 39.1 | 44.7 | 42.5
1/1 | exact | 4 (2) | 5 (2) | 37.2 | 39.1 | 44.7 | 42.5
2/2 | aligned | 2 (2) | 3 (2) | 37.8 | 39.9 | 48.6 | 43.8
2/2 | partial | 4 (2) | 4 (2) | 37.7 | 39.9 | 48.6 | 43.8
2/2 | exact | 4 (2) | 5 (2) | 37.7 | 39.9 | 48.6 | 43.8
d/1 | aligned | 4 (2) | 3 (2) | 37.0 | 39.1 | 44.4 | 42.6
d/1 | partial | 5 (2) | 3 (2) | 36.8 | 39.0 | 44.6 | 42.4
d/1 | exact | 5 (2) | 4 (2) | 36.8 | 39.0 | 44.6 | 42.4
1/d | aligned | 2 (2) | 3 (2) | 37.0 | 38.9 | 44.0 | 42.3
1/d | partial | 3 (2) | 3 (2) | 37.0 | 38.9 | 44.0 | 42.3
1/d | exact | 3 (2) | 7 (2) | 37.0 | 38.9 | 44.0 | 42.3
d/d | aligned | 3 (2) | 3 (2) | 37.8 | 39.9 | 45.5 | 43.4
d/d | partial | 3 (2) | 5 (2) | 36.7 | 38.8 | 44.0 | 42.2
d/d | exact | 5 (2) | 5 (2) | 36.7 | 38.8 | 44.0 | 42.2

The costliest edge three ways, on the primary reading. The brief's fall in π is the first two columns together, and the prior's part of it is the entropy of the terms held. The likelihood's part is the barrier in the sense SC-5 bounds.

pair | from | k (net) | Δ log marginal | Δ log prior | proposal | sum | the edge | fall's share (brief's split, class)
--- | --- | --- | --- | --- | --- | --- | --- | ---
A – a1 | the first | 2 (2) | 7.4 | 18.2 | 19.0 | 44.6 | delete x1-360·x2-408 | 57%
A – a1 | a1 | 4 (2) | 19.6 | 19.6 | 3.4 | 42.6 | delete x1+342·x2+170 | 92%
B – a1 | the first | 1 (1) | 19.9 | 8.4 | 11.2 | 39.5 | delete x2+326 | 72%
B – a1 | a1 | 4 (2) | 24.8 | 9.8 | 18.8 | 53.4 | delete x1+185·x2+321 | 65%
a2 – a1 | the first | 0 (0) | -0.0 | -0.0 | 15.3 | 15.3 | slide of 240 | 0%
a2 – a1 | a1 | 3 (1) | 3.2 | 9.8 | 17.1 | 30.1 | slide of 240 | 58%
a1 + x6·x8 – a1 | the first | 1 (0) | -0.0 | -0.0 | 14.5 | 14.5 | slide of 45 | 0%
a1 + x6·x8 – a1 | a1 | 2 (1) | -7.0 | 18.2 | 14.4 | 25.7 | slide of 45 | 51%
chain 0 – a1 | the first | 2 (2) | 6.1 | 18.1 | 11.8 | 36.0 | add x1+3·x2+1 | 67%
chain 0 – a1 | a1 | 6 (4) | 29.1 | 46.1 | 18.4 | 93.6 | delete x1+3·x2+1 | 80%
chain 1 – a1 | the first | 3 (0) | 26.8 | -9.8 | 11.2 | 28.2 | add x1+185·x2+321 | 60%
chain 1 – a1 | a1 | 2 (1) | 54.5 | -8.4 | 19.1 | 65.1 | delete x1+185·x2+321 | 71%
chain 2 – a1 | the first | 2 (0) | 32.6 | -18.1 | 19.2 | 33.7 | delete x1-319·x2-226 | 43%
chain 2 – a1 | a1 | 5 (3) | 60.5 | -7.0 | 19.2 | 72.8 | delete x1+110·x2+93 | 74%
chain 3 – a1 | the first | 2 (0) | 25.8 | -0.0 | 19.0 | 44.8 | delete x1-151·x2+125 | 58%
chain 3 – a1 | a1 | 3 (1) | 58.6 | 1.4 | 19.0 | 79.0 | delete x1+110·x2+93 | 76%

- *Barrier against proposal.* From the chains' states on the primary reading, the costliest edge costs 28–45 nats. Of that, 11–19 is the proposal and the rest the fall in π. By the brief's split the fall is 43–67% by class flow.
- *The fall is not all likelihood.* It is 6–33 nats of log marginal from the chains, with the prior between −18 and +18, because the route deletes before it adds.
- *SC-5 bounds the log marginal, not π.* In π an unpaying addition costs c plus its term's measure: 12.9 for one factor and 21.3 for two at N = 500, of which c is 3.11. So the superset route is not the best route. It is 54–69 nats deep from the chains on the aligned lattice and 101–114 on the exact one, against 15–26 for the maximin route. The maximin route swaps a smaller likelihood valley for not holding unpaid terms.

#### One route's profile

A–a1, aligned, held τ, from A. Read term by term: the additions with Δπ < 0 before the lowest state (π = −922.20, after step 1) are steps 0 and 1, so k = 2, as the table has it. From a1 the lowest state is reached through steps 6, 4, 3 and 2 read backward, each an addition with Δπ < 0, so k = 4. Steps 5 and 6 are a swap: x4+6 for x4−496, nearly one column under two signs. Each addition's Δ log marginal respects SC-5's −3.11.

step | move | term | M | Δ log marginal | Δ log prior | Δπ | c | π after | ln q forward | ln q reverse | ln c_e | ln c_e, class | class size
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
0 | add | x1+3·x2+1 | 8 → 9 | -2.32 | -18.23 | -20.55 | 3.11 | -905.44 | -19.30 | -3.30 | -908.74 | -902.26 | 783
1 | add | x1+342·x2+170 | 9 → 10 | 1.42 | -18.17 | -16.76 | 3.11 | -922.20 | -19.02 | -3.40 | -925.60 | -917.93 | 2601
2 | delete | x1−360·x2−408 | 10 → 9 | -6.45 | 18.17 | 11.72 | 3.11 | -910.47 | -3.40 | -19.02 | -929.50 | -921.64 | 2601
3 | delete | x1−370 | 9 → 8 | -0.40 | 9.82 | 9.42 | 3.11 | -901.06 | -3.30 | -11.14 | -913.77 | -909.85 | 51
4 | delete | x2−441 | 8 → 7 | 0.66 | 9.87 | 10.54 | 3.11 | -890.52 | -3.18 | -11.20 | -904.23 | -900.31 | 51
5 | add | x4+6 | 7 → 8 | -2.90 | -9.87 | -12.78 | 3.11 | -903.29 | -11.20 | -3.18 | -906.47 | -902.80 | 32
6 | delete | x4−496 | 8 → 7 | 2.74 | 9.87 | 12.61 | 3.11 | -890.68 | -3.18 | -11.20 | -906.47 | -903.37 | 29
7 | add | x1+110·x2+93 | 7 → 8 | 3.58 | -18.29 | -14.71 | 3.11 | -905.38 | -19.12 | -3.18 | -909.80 | -902.06 | 2601
8 | slide, 124 knot steps | lowest π -905.38 | 8 | 22.37 | 0.00 | 22.37 | | -883.02 | -10.78 | -10.78 | | |

Chain 1–a1, aligned, held τ. This one has a real likelihood valley: step 4 deletes x2−479 at −24.84 before step 5's replacement pays most of it back. Net k is 0 at the lowest state.

step | move | term | M | Δ log marginal | Δ log prior | Δπ | c | π after | ln q forward | ln q reverse | ln c_e | ln c_e, class | class size
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
0 | delete | x1−377·x2−450 | 9 → 8 | -9.44 | 18.23 | 8.79 | 3.11 | -903.27 | -3.30 | -19.30 | -922.57 | -914.70 | 2601
1 | add | x4+6 | 8 → 9 | -1.77 | -9.82 | -11.59 | 3.11 | -914.86 | -11.05 | -3.30 | -918.15 | -915.21 | 32
2 | delete | x4−476 | 9 → 8 | 3.11 | 9.82 | 12.92 | 3.11 | -901.93 | -3.30 | -11.05 | -918.15 | -914.11 | 49
3 | add | x1+3·x2+1 | 8 → 9 | 6.10 | -18.23 | -12.12 | 3.11 | -914.06 | -19.30 | -3.30 | -921.23 | -914.57 | 783
4 | delete | x2−479 | 9 → 8 | -24.84 | 9.82 | -15.02 | 3.11 | -929.08 | -3.30 | -11.14 | -940.22 | -936.39 | 46
5 | add | x1+185·x2+321 | 8 → 9 | 17.94 | -18.23 | -0.29 | 3.11 | -929.37 | -19.07 | -3.30 | -948.15 | -940.29 | 2601
6 | delete | x1+60·x2−304 | 9 → 8 | 3.09 | 18.23 | 21.32 | 3.11 | -908.05 | -3.30 | -19.07 | -932.66 | -924.78 | 2601
7 | delete | x1−447 | 8 → 7 | 3.11 | 9.87 | 12.98 | 3.11 | -895.07 | -3.18 | -11.20 | -911.23 | -907.29 | 51
8 | add | x1+110·x2+93 | 7 → 8 | 6.65 | -18.29 | -11.64 | 3.11 | -906.71 | -19.12 | -3.18 | -914.19 | -906.32 | 2601
9 | slide, 141 knot steps | lowest π -906.71 | 8 | 23.69 | 0.00 | 23.69 | | -883.02 | -10.78 | -10.78 | | |

#### Part 2: k against N

Lower confidence than Part 1 throughout: these modes are carried over, not found. Each row is aligned, held τ, and `aligned+redraw` where a closed pair exists. *First* is A, B or the best chain's structure.

pair | data | N | c | holds (first, a1) | closed pair? | d_first | d_a1 | k from first (net) | k from a1 (net) | ln wait from first | ln wait from a1 | superset d_first | d_first, partial
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
A-a1 | fresh | 125 | 2.42 | yes, yes | yes | 29.2 | 23.8 | 1 (0) | 2 (0) | 45.7 | 33.3 | 49.2 | 27.9
A-a1 | fresh | 250 | 2.76 | yes, yes |  | 26.6 | 21.9 | 2 (0) | 2 (0) | 44.5 | 32.0 | 63.1 | 26.6
A-a1 | fresh | 500 | 3.11 | yes, yes |  | 38.0 | 20.0 | 4 (0) | 0 (0) | 56.2 | 31.8 | 70.2 | 36.5
A-a1 | fresh | 1000 | 3.45 | yes, no | yes | 51.1 | 43.7 | 3 (0) | 0 (0) | 70.9 | 56.4 | 73.4 | 45.1
A-a1 | fresh | 2000 | 3.80 | no, yes |  | 58.9 | 114.4 | 1 (1) | 3 (1) | 80.7 | 128.3 | 75.8 | 13.2
A-a1 | stream | 125 | 2.42 | yes, yes |  | 15.4 | 4.1 | 3 (0) | 0 (0) | 29.8 | 18.5 | 59.3 | —
A-a1 | stream | 250 | 2.76 | yes, no | yes | 20.0 | 15.7 | 1 (0) | 2 (0) | 37.9 | 25.8 | 50.9 | —
A-a1 | stream | 1000 | 3.45 | yes, no | yes | 34.9 | 51.9 | 2 (2) | 3 (2) | 55.2 | 64.4 | 56.2 | —
A-a1 | stream | 2000 | 3.80 | no, yes |  | 31.1 | 90.1 | 3 (1) | 2 (1) | 52.9 | 104.0 | 52.4 | —
B-a1 | fresh | 125 | 2.42 | yes, yes | yes | 16.8 | 24.1 | 1 (-1) | 1 (0) | 33.1 | 34.2 | 35.6 | 9.6
B-a1 | fresh | 250 | 2.76 | yes, yes |  | 9.0 | 16.0 | 2 (-1) | 1 (0) | 25.0 | 27.7 | 49.1 | 8.2
B-a1 | fresh | 500 | 3.11 | no, yes |  | 38.7 | 35.7 | 3 (-1) | 0 (0) | 53.0 | 47.8 | 53.3 | 29.9
B-a1 | fresh | 1000 | 3.45 | yes, no |  | 32.4 | 32.2 | 3 (1) | 2 (2) | 38.8 | 39.3 | 56.5 | 32.4
B-a1 | fresh | 2000 | 3.80 | no, yes |  | 44.7 | 73.8 | 2 (0) | 1 (1) | 66.3 | 93.5 | 57.9 | 44.5
B-a1 | stream | 125 | 2.42 | no, yes |  | 8.9 | 15.7 | 0 (-1) | 3 (0) | 20.6 | 32.0 | 42.3 | —
B-a1 | stream | 250 | 2.76 | yes, no |  | 17.9 | 25.9 | 0 (-1) | 2 (0) | 35.5 | 43.6 | 48.1 | —
B-a1 | stream | 1000 | 3.45 | yes, no |  | 34.8 | 41.7 | 2 (1) | 3 (2) | 43.2 | 57.4 | 56.8 | —
B-a1 | stream | 2000 | 3.80 | no, yes |  | 42.9 | 68.4 | 1 (0) | 3 (1) | 58.6 | 90.0 | 58.9 | —
best chain-a1 | fresh | 125 | 2.42 | no, yes | yes | 10.0 | 18.4 | 3 (-1) | 0 (0) | 20.1 | 28.5 | 49.2 | 0.0
best chain-a1 | fresh | 250 | 2.76 | no, yes |  | 1.9 | 15.5 | 2 (0) | 3 (1) | 13.4 | 23.7 | 67.6 | 1.9
best chain-a1 | fresh | 500 | 3.11 | no, yes |  | 31.1 | 36.3 | 4 (-1) | 0 (0) | 44.0 | 50.4 | 70.7 | 21.0
best chain-a1 | fresh | 1000 | 3.45 | yes, no | yes | 16.5 | 23.8 | 2 (0) | 2 (1) | 24.0 | 34.2 | 62.4 | 16.5
best chain-a1 | fresh | 2000 | 3.80 | yes, yes |  | 28.6 | 82.9 | 2 (1) | 3 (2) | 50.1 | 104.1 | 54.5 | 25.4
best chain-a1 | stream | 125 | 2.42 | no, yes |  | 0.4 | 3.9 | 1 (0) | 4 (1) | 13.3 | 12.9 | 61.7 | —
best chain-a1 | stream | 250 | 2.76 | no, no |  | 2.9 | 12.9 | 1 (-1) | 3 (0) | 20.8 | 23.5 | 65.8 | —
best chain-a1 | stream | 1000 | 3.45 | yes, no | yes | 30.5 | 60.8 | 2 (1) | 3 (2) | 44.8 | 67.5 | 48.9 | —
best chain-a1 | stream | 2000 | 3.80 | yes, yes |  | 30.5 | 104.4 | 2 (0) | 3 (1) | 52.3 | 126.2 | 46.1 | —

The carried structures, polished at held τ, then 2,000 sweeps at b = 1:

N | data | label | M | π carried | π polished | passes | share in basin | top basin (share) | rank | start − median
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
125 | fresh | a1 | 8 | -309.9 | -299.9 | 6 | 92% | a1 (92%) | 0.996 | 3.1
125 | fresh | A | 8 | -304.3 | -294.5 | 9 | 92% | A (92%) | 1.000 | 4.6
125 | fresh | B | 9 | -313.5 | -307.2 | 4 | 53% | B (53%) | 0.670 | 1.7
125 | fresh | chain1 | 9 | -321.2 | -308.2 | 7 | 1% | A (90%) | 0.086 | -11.5
125 | stream | a1 | 8 | -322.6 | -314.0 | 6 | 91% | a1 (91%) | 0.987 | 2.8
125 | stream | A | 8 | -314.4 | -302.7 | 7 | 91% | A (91%) | 0.990 | 6.0
125 | stream | B | 9 | -328.9 | -320.8 | 5 | 2% | A (84%) | 0.100 | -15.6
125 | stream | chain1 | 9 | -329.1 | -317.5 | 9 | 0% | A (79%) | 0.129 | -12.7
250 | fresh | a1 | 8 | -518.8 | -504.3 | 8 | 89% | a1 (89%) | 0.996 | 3.9
250 | fresh | A | 8 | -511.4 | -499.7 | 6 | 90% | A (90%) | 0.996 | 5.0
250 | fresh | B | 9 | -523.0 | -511.3 | 11 | 78% | B (78%) | 0.993 | 2.8
250 | fresh | chain1 | 9 | -525.6 | -517.9 | 7 | 11% | A (70%) | 0.139 | -11.3
250 | stream | a1 | 8 | -521.5 | -515.7 | 5 | 27% | 3 × x1·x2 + 2 × x3 + x4 + x5 + x5·x6 (33%) | 1.000 | 6.2
250 | stream | A | 8 | -517.2 | -511.4 | 4 | 71% | A (71%) | 0.768 | 2.7
250 | stream | B | 9 | -527.6 | -523.7 | 5 | 58% | B (58%) | 0.998 | 6.2
250 | stream | chain1 | 9 | -533.5 | -525.6 | 5 | 1% | A (77%) | 0.156 | -13.0
500 | fresh | a1 | 8 | -927.8 | -901.1 | 17 | 91% | a1 (91%) | 0.999 | 5.0
500 | fresh | A | 8 | -931.3 | -883.1 | 10 | 87% | A (87%) | 0.997 | 2.1
500 | fresh | B | 9 | -927.5 | -898.1 | 16 | 13% | x1 + 3 × x1·x2 + 2 × x3 + x4 + x5 (40%) | 0.676 | 3.3
500 | fresh | chain1 | 9 | -946.9 | -906.3 | 13 | 20% | A (69%) | 0.260 | -15.3
1000 | fresh | a1 | 8 | -1629.9 | -1598.5 | 17 | 12% | 3 × x1·x2 + 2 × x3 + x4 + x5 (77%) | 0.191 | -14.4
1000 | fresh | A | 8 | -1677.9 | -1591.1 | 24 | 94% | A (94%) | 0.995 | 1.6
1000 | fresh | B | 9 | -1644.5 | -1598.2 | 15 | 92% | B (92%) | 0.914 | 0.8
1000 | fresh | chain1 | 9 | -1651.2 | -1605.8 | 18 | 92% | B (92%) | 0.997 | 3.2
1000 | stream | a1 | 8 | -1642.6 | -1636.0 | 10 | 40% | 4 × x1·x2 + x2·x6·x8 + 2 × x3 + x4 + x5 (52%) | 0.999 | 23.1
1000 | stream | A | 8 | -1672.7 | -1653.0 | 9 | 87% | A (87%) | 0.974 | 5.3
1000 | stream | B | 9 | -1653.2 | -1642.9 | 13 | 97% | B (97%) | 1.000 | 3.6
1000 | stream | chain1 | 9 | -1704.0 | -1666.4 | 19 | 65% | B (65%) | 0.998 | 8.6
2000 | fresh | a1 | 8 | -3132.1 | -3079.0 | 50 | 57% | a1 (57%) | 0.998 | 3.8
2000 | fresh | A | 8 | -3213.7 | -3134.4 | 25 | 2% | 3 × x1 + 2 × x1·x2 + x2 + x2·x5 + 3 × x3 + x4 + x5 (31%) | 0.586 | 1.9
2000 | fresh | B | 9 | -3145.4 | -3108.1 | 19 | 30% | x1 + 3 × x1·x2 + x2 + 4 × x3 + x4 + x5 (41%) | 0.989 | 6.7
2000 | fresh | chain1 | 9 | -3203.0 | -3133.3 | 39 | 87% | B (87%) | 0.991 | 2.6
2000 | stream | a1 | 8 | -3131.5 | -3077.7 | 50 | 97% | a1 (97%) | 0.930 | 0.9
2000 | stream | A | 8 | -3217.3 | -3136.7 | 31 | 1% | 3 × x1 + 2 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5 (77%) | 0.999 | 12.2
2000 | stream | B | 9 | -3142.2 | -3103.2 | 25 | 37% | x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5 + x6·x8 (48%) | 0.985 | 15.1
2000 | stream | chain1 | 9 | -3225.3 | -3151.5 | 32 | 88% | B (88%) | 0.991 | 1.0

- *Only rows where both ends hold count.* That leaves A–a1 at N = 125, 250 and 500 on fresh cases and N = 125 on the stream, B–a1 at N = 125 and 250 on fresh cases, and the best chain's structure against a1 at N = 2000 on both. Beyond N = 500 only that last pair holds at both ends. a1 fails to hold at N = 1000 on both data sets, and A and B at N = 2000.
- *Depth at the closed gate.* The literal aligned rows read 97, 46 and 83 nats at N = 125, 209 at N = 250 and 1,009–1,029 at N = 1000. Taken by one redraw they read 10–51. No depth exceeds the superset route's on its own lattice, and the superset route's depth from the first structure stays in 35–76 at every N.
- *Waiting* on the rows where both ends hold runs from ln 25 to ln 56 sweeps from the first structure, with no trend so few rows can carry.

#### How the outcomes read

*Waiting time, per chain-held pair against a1, class flow, the better τ:* **at least 10× the sampling phase on all four.** So each route is too slow to have been taken, which is consistent with the chains staying put. It does not show that no faster route exists. The better τ is held on every row. A sweep attempts one structure move (`MoveStep`), and τ's slice proposes none, so attempts are sweeps. The sampling phase is `BassMixing.Samples`, 2,000, so the rule's lines are ln 20,000 = 9.9 and ln 200 = 5.3.

chain | route | ln wait held | ln wait profiled | less ln 2,000 (better τ) | reading | best in L, exact flow | every addition's class credited
--- | --- | --- | --- | --- | --- | --- | ---
0 | aligned | 36.6 | 36.9 | 29.0 | ≥ 10× | 42.9 | 35.4
0 | partial | 36.2 | 36.5 | 28.6 | ≥ 10× | 42.5 | 35.0
0 | exact | 40.9 | 41.2 | 33.3 | ≥ 10× | — | 25.5
1 | aligned | 28.9 | 29.1 | 21.3 | ≥ 10× | 23.4 | 18.8
1 | partial | 27.6 | 27.8 | 20.0 | ≥ 10× | — | 13.7
1 | exact | 27.6 | 27.8 | 20.0 | ≥ 10× | — | 11.0
2 | aligned | 33.7 | 34.0 | 26.0 | ≥ 10× | 33.7 | 23.5
2 | partial | 31.7 | 32.0 | 24.1 | ≥ 10× | — | 21.6
2 | exact | 31.7 | 32.0 | 24.1 | ≥ 10× | — | 13.9
3 | aligned | 44.8 | 45.2 | 37.2 | ≥ 10× | 42.4 | 34.7
3 | partial | 40.8 | 41.1 | 33.2 | ≥ 10× | 34.1 | 33.2
3 | exact | 40.8 | 41.1 | 33.2 | ≥ 10× | — | 22.8

Two readings beside the rule, neither of which moves it:
- *The least waiting time over every route in L* on exact flows, where |D| ≤ 14. It lands well under the maximin route's class reading on two rows, 23.4 for chain 1 aligned and 34.1 for chain 3 partial, and both are still far past the line.
- *Every addition's class credited to all later edges,* which the brief's class flow does not do. This is a generous bound: it credits swaps too, and treats a class's members as equal in π. Its lowest is 11.0 for chain 1 on the exact route, about 60,000 sweeps, against the rule's 20,000.

*On k:* **k ≥ 3 on some pair with the fall in π the larger share,** so one-term moves meet SC-5's exponential in k. On the primary reading, from a1: A 4 (92% by class), B 4 (65%), chain 0 6 (80%), chain 2 5 (74%), chain 3 3 (76%), a2 3 (58%). From a chain's state: chain 1, 3 (60%). "k ≤ 2 on every pair" fails. Three things qualify it:
- k is not a stable reading, as the first paragraph says.
- Net k from the first structure is 0–2 on every pair, and from a1 at most 4. So most of the falls from a chain's state are one term lost before its replacement pays, not several unpaid terms held together.
- From the chains' states the split is close to even, 43–67%. Chain 2's costliest edge is proposal-led.

*Across N (Part 2):* **flat on what can be read.** On A–a1 from A, k is 1, 2 and 4 at N = 125, 250 and 500 on fresh cases, 3 at N = 125 on the stream and 2 on S10's own N = 500, with net k 0, 0, 0, 0 and 2. Its spread within one N is as wide as the rise. At N = 2000 the best chain's structure reads k = 2 from itself and 3 from a1. No rising k stands where both ends hold, and the one depth that seemed to grow is the closed gate. So nothing here breaks the register's polynomial expectation. It does not confirm it either: see the first limit below.

#### Limits

- **k against N is flat partly by construction.** Carrying an exemplar fixes its term count, so |D| stays 5–9 at every N and k cannot grow far. The chains say the fixed count stops being the mode: at N = 2000, A and B leave for basins with more terms (up to 3 × x1 + 2 × x1·x2 + x2 + x2·x5 + 3 × x3 + x4 + x5). That is SC-4's "M_n grows with N", the case that would matter, and it is unread. Reading it needs modes found at N = 2000, not carried.
- **One flat-histogram run, seed 8017, on one data seed, one surface.** All of Part 1 rests on these. a2's and a1 + x6·x8's exemplars come from 112 and 3,559 draws.
- **One mode or state each.** Each chain pair rests on one state, and chain 1's is not its chain's mode. A's second exemplar shares 5 of 8 terms with its first, a1's only 1 of 8. Only A–a1 was read across exemplars.
- **L is one family of routes.** It uses only the terms of S and T, plus knot steps on the slide. No other knot or term, and no redraw except a closed pair's, though production's redraw can replace a factor in one step. On small fixtures that gap is 2–3 nats at the median and at most 9.6. On Friedman it is unread. The routes bound, they do not decide. The partial search is greedy and stops at the first set of dropped pairs no single drop improves.
- **The waiting time is a series reading of one route,** on the assumption of like exits, not a mixing bound. The class flow is an approximation, one edge at a time, within ±25 knots.
- **τ.** Moves are read at fixed τ. The profile sets each state at its own best τ, which no chain does. It moves the depths by under 2 nats on Part 1's eight main pairs.
- **Part 2's carry is crude.** It maps knots by quantile and polishes greedily at held τ. a1 on fresh N = 500 ends 18 nats under A there, so whether a1 is that surface's mode at all is unread.

#### For the main session

**G1,** under "How to close it", candidate reading (a) is now measured. Suggested wording:
- Between Friedman's modes at N = 500 (one run, one data seed), the best routes found through L hold at most 2 terms above the start at their lowest state, and at most 4 from a1. They reach 15–26 nats below the chains' states and 28–37 below B and A, with 45–76 from a1 to the chains.
- The waiting time on them is e^{28}–e^{45} sweeps from the chains' states against a 2,000-sweep phase, consistent with the chains staying put. At the costliest edge the proposal is 11–19 nats of that, and the fall the rest.
- k counted by additions is not stable across exemplars; net k is.
- "The real k between a1 and the reference basins is unread" should go. So should "at N = 500 each such term costs about 3.1 nats, so eight cost about e^{−25}", since 3.1 is the Occam part alone: in π each unpaid term costs 12.9 or 21.3 nats by its order, before the ±25-knot class's credit of at most ln 51 or ln 2,601.
- k against N reads flat on carried modes, partly by construction. What remains open is k at modes found at large N.

**SC-4:**
- *"Polynomial does not mean fast. At N = 500 the births alone cost about e^{−5} each, so k = 4 costs e^{−20} before any proposal entropy."* The route reading adds that a term's proposal entropy does not stay out of it: the term's measure is in π for as long as the term is held unpaid, so it accumulates along a route. Four unpaid two-factor terms are up to 85 nats by the superset route, and the maximin route avoids them by passing a likelihood valley instead.
- *"D06 refuses only near-exact collinearity, where a removal costs little anyway."* This is wrong for a term and its mirror: a hinge at the lowest knot and its reflection at the highest are one column, D06 refuses them together, and removing the one held is load-bearing, −991 nats on x4 at N = 1000. Production crosses it with one redraw, so the sampler is not stuck, but the superset path is closed there. Any canonical-path argument built on SC-5, including the positive theorem G1 proposes, needs a redraw edge for such pairs or a rule that excludes them.
- *On Yang, Wainwright & Jordan.* Near-duplicate knots are not what sets the depth here: at N = 500 no state is refused, and swaps between near twins add steps and k without lowering the route's minimum. Their cost shows in the waiting time, as 1/(2n_v) per redraw with the local steps off.

**SC-5** stands as proved, with 32,577,747 edges read and no violation. Its statement is about the log marginal, and the register's and SC-4's sentences that read it as π's depth are the ones above.

**A next reading, for you to weigh:** routes between modes found at N = 2000, which is the case Part 2 could not reach.

`artifacts\chip-pathprofile\` is deleted. The CSVs are `path-profile-modes.csv`, `-routes.csv`, `-steps.csv`, `-validation.csv`, `-groundtruth.csv` and `-carried.csv`, beside this brief.
