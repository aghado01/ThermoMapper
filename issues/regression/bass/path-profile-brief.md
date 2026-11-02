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
