# BASS F-M2: pricing the mode exits (chip brief)

Date: 2026-10-29. Follows Phase 0 ([its brief and report](local-move-measurements-brief.md#reports)) and the design draft's "After Phase 0: proposals" section in [local-move-mixing-brief.md](local-move-mixing-brief.md). The user decided on 2026-10-29 to price these hypotheses before choosing Phase 1's levers; this brief is that work. Like Phase 0, it measures and changes nothing.
- *What Phase 0 found.*
  - The chains hold local modes 12–42 nats apart (M4). Each worse mode carries co-adapted terms that are each load-bearing: on Friedman, 38 of 41 terms have death α below 10⁻⁶ (M2).
  - A local knot step accepts about 0.95 (M3).
  - With the steep band opened, the ladder's ideal round-trip rate reached 0.03 per scan and no round trip completed (M5).
- *The hypotheses this chip prices,* P1–P4 below:
  - co-adapted pairs of terms that are jointly redundant;
  - structural moves at fixed M that re-represent a term;
  - a drift path along which sliding knots make a term redundant;
  - mid-ladder rungs frozen by b·Δℓ.
- *Not in scope.*
  - Building any move into production code.
  - Any default.
  - The design draft's questions, which are the user's.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the design draft's sections "Phase 0's readings" and "After Phase 0: proposals", which define the moves and their ratios;
- Phase 0's brief and report in full: its harness, its no-commit check, and its fixed states;
- the fitting brief's D01 (the term prior, the birth proposal, the move schedule).

In the code:
- *Phase 0's harness,* `tests/regression/BassMixing.cs`: fixtures, burned chains and ladders, `Score`, and the no-commit check. Its tests `BassMixing*Tests.cs` show how each was used.
- *The moves and the design:* `BassTermProposal.cs` (`LogBirth` through `LogBirthRatio` and `LogDeathRatio`), `BassSweep.cs` (`BassStructureMoves`), `BassDesign.cs` (`ProposeChange`, `ProposeDeath`, `Load`) and `BassTarget.cs` (`Evaluate`, `LogStructurePrior`; the tempered forms' likelihood at b).

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-exits\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-exits\` when you finish.
- **The main session works beside you,** caching the composition's column means in `src/maths/regression/spline/bass/BassSweep.cs`. That change is bit-identical and does not touch the structure moves. Leave the file alone. If a build breaks while its edit is in progress, wait and retry.
- **Your files.**
  - New test files `tests/regression/BassMixingExit*Tests.cs`.
  - Additions to `tests/regression/BassMixing.cs`, without changing what Phase 0's tests read: rerun M2–M4 and check that their readings are unchanged.
  - This brief's Reports section.
  - Leave everything else alone, the design draft included; name any change you need in your report.
- **Production code stays as it is.** The candidate moves are test-side scorers, not moves. Where a measurement cannot be taken without a change, the only one allowed is an internal accessor that changes no behaviour; name it, and rerun the Regression suite.
- **The harness changes nothing it measures.** Apply Phase 0's no-commit check after every scorer.
- **Commit your own files only, by explicit path,** one commit per measurement or a coherent pair. Use the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF.
- **Read, do not judge.** Each harness prints its readings and asserts only its own soundness. Recommend no lever; the readings feed the user's design pass.

## The measurements

- *Fixed states.* Phase 0's: each chain's state at the end of a 5,000-sweep burn-in, on S10's Friedman, mixed categorical and direct functional fits, four chains each.
- *Better and worse modes,* by M4's medians: Friedman's best is chain 1, and chains 0, 2 and 3 sit 24–42 nats below; on the mixed fit chain 2 sits 22 below the rest; on the functional fit chains 1 and 3 sit about 29 below chains 0 and 2.
- *First, check those states are in their chains' modes.* Compare each state's log posterior with its chain's M4 median. Report it, and if a state lies outside its chain's sampling range, use the chain's most frequent structure instead.
- *Report everything per chain,* marking which are the worse modes.

**P1. Pair deaths** (`BassMixingExitPairTests`).
- *For every unordered pair of terms in each state,* record:
  - Δℓ for removing both;
  - the pair death's full log ratio. That is the target's ratio (log marginal and log structure prior, through the target) plus the proposal part ln q_pair(t₁, t₂ | S′) + ln C(M, 2), where S′ is the reduced structure and q_pair(t₁, t₂ | S′) = q(t₁ | S′)·q(t₂ | S′ + t₁) + q(t₂ | S′)·q(t₁ | S′ + t₂), each q the D01 birth probability (`LogBirth`). Take the move-selection ratio as 1, and say so.
- *Against the single deaths:* for each pair, the two terms' single-death Δℓ. A pair whose joint removal costs far less than either alone is co-adapted.
- *Where M ≤ 12, also every triple's Δℓ* for the states whose best pair costs under 20 nats.
- *Report:* each state's best pairs (Δℓ, full ratio, the terms by name) and the count of pairs within 10 and 20 nats. Say whether the worse modes hold jointly redundant pairs the best mode does not.

**P2. Structural moves at fixed M** (`BassMixingExitFactorTests`).
- *Prune:* every factor of every term with K ≥ 2 in its part. Record Δℓ and the full log ratio, which the design draft derives as Δℓ alone.
  - Assert that the target's term-prior ratio and the proposal ratio cancel, to 1e-12 relative. This checks the derivation, so report it if it fails.
  - A prune onto a duplicate term, or below the support count, is a refusal: report its share.
- *Grow and variable change:* 5,000 proposals each per state, drawn as the draft defines them (grow: the variable uniform among the term's unused ones, the sign and knot from the prior; change: v → v′ uniform among the unused ones, the knot from v′'s prior). Record Δℓ, the cancellation check, and refusals.
- *Group:*
  - by part (continuous, categorical, functional);
  - on Friedman, by whether the factor's variable is active or inert;
  - by the term's order before and after.
- *Two-step profile:* after each prune or variable change costing under 20 nats, record the best single death and the best pair death in the new state. The question is whether a cheap re-representation opens a cheap removal.
- *Report:* acceptance (min(1, e^{ratio})) and Δℓ quantiles per group; each worse-mode state's ten cheapest exits by name; the two-step profile.

**P3. Profile deaths** (`BassMixingExitProfileTests`).
- *For each term of each state:* remove it, then re-optimize the remaining hinges' knots greedily.
  - In passes over the factors, try j ± δ for δ ∈ {1, …, 8}. Keep any step that raises the log marginal plus the log structure prior; the knot prior is uniform, so this is the marginal.
  - Stop after 5 passes or a pass without gain.
  - Categorical subsets stay fixed.
- *Record:* the plain death's Δℓ; the profile's Δℓ (after re-optimization, against the original state); the passes taken; and the knots' total displacement in distinct-value steps.
- *Report:* per state, plain against profile Δℓ for each term, and how many terms' profile Δℓ lies within 10 and 20 nats.
  - This is an optimization diagnostic: it bounds where a drift path could end, not that one exists or how long it takes.
  - Say so in the reading.

**P4. Structure moves along the dense ladder** (`BassMixingExitLadderTests`).
- *The ladder:* M5's Part 2 dense fixed ladder on Friedman, each chain, at the end of its burn-in.
- *At each rung's state, scored under that rung's tempered target:*
  - births (2,000 per rung);
  - every term's death;
  - change redraws (2,000);
  - local knot steps, δ ∈ {±1, ±2, ±3}, for every hinge factor;
  - and P2's prunes.
  - Use the ladder's own tempered form (the integrated form by default), through the target's evaluation at b, with the no-commit check per rung.
- *Report:* per rung (b, M), the mean α of each move kind, and the deaths' median Δℓ scaled by b. Say at which b deaths stop being negligible (α above 10⁻³), and whether the band's rungs are frozen.

## Review before reporting

Before you write the report, take a fresh-eyes pass over the harnesses and the readings, and fix what it finds:
- *Correctness:*
  - the no-commit check holds;
  - P1's q_pair sums both draw orders;
  - P2's cancellation check passes, or its failure is reported;
  - P4 scores at the rung's b and not at 1;
  - the counts add up.
- *Coverage:*
  - a state not in its chain's mode;
  - a reading that rests on one chain;
  - P3's greedy search stopping at its pass limit rather than at a local optimum.

Name the passes and their fixes in the report.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- each measurement's tables and its reading, kept apart: the numbers, then what they show, and how strongly;
- the run times;
- any accessor added, and the suites rerun;
- the review passes;
- anything the main session must know.

## Reports

Newest first.

### Mode exits priced: P1–P4 — 2026-10-29

*Landed:* tests only. No production change, no internal accessor, and `BassSweep.cs` untouched.
- `17f8b5f`: the additions to `BassMixing.cs`, and P1, `BassMixingExitPairTests`.
- `108cefe`: P2, `BassMixingExitFactorTests`. `b735551`: P3, `BassMixingExitProfileTests`.
- `bc1e6da`: P4, `BassMixingExitLadderTests`, with `811360e`, its header fix.
- Each runs behind `BASS_QUALIFICATION` = 1, reads, and asserts only its own soundness.
- *The additions to the helper.* None changes what Phase 0's tests read.
  - `FixedStates` burns each chain as before. A twin restored from the chain's snapshot continues it bit for bit through the 2,000-sweep sampling phase.
    - The run's retained draws must hold the twins' M, draw for draw, so the states are S10's and the phases M4's.
    - A state outside its phase's range is replaced by the phase's most frequent structure.
  - `Removals` removes k terms at once on a fresh design. Its reverse is k D01 births from the reduced structure, summed over their k! orders, over 1/C(M, k).
  - The design draft's fixed-M moves are test-side scorers: `Without`, `With`, `Replaced` and `Stepped`, with their proposal ratios (`LogPruneRatio`, `LogGrowRatio`). Each factor's prior is computed from the inputs, not through `BassStructurePrior`.
  - `Score` takes a rung's temperature, the posterior by default (the same bits at b = 1). `Untouched` takes a ladder replica's moves.
  - M2–M4 reran on the final helper, and their readings are Phase 0's, number for number: M2's losses and deaths, M3's buckets, M4's medians and gaps.
- *Choices the brief left open,* stated here and in the harnesses' docstrings:
  - *Move selection.* It is taken as matched (ratio 1) for the pair death, and for the factor moves' choice of move and of part.
  - *Several parts.* The prune, grow and change act within one part.
    - Beside the brief's prunes (K ≥ 2 in a part), P2 also scores a part's only factor where the term keeps another part (K = 1 → 0), tagged "prune, part emptied". It is the reverse of a grow into an empty part.
    - On the functional fit no state holds a part with K ≥ 2, so this is its only prune.
  - *P3's greedy step:* per factor, the best of its 16 knot steps, kept if it raises the objective by more than 10⁻⁹.
  - *P4's change α* leaves out redraws equal to the current factor, which are null moves.

*The fixed states,* each against its own chain's sampling phase. The medians are M4's, exactly.

| fixture | chain | M | log post | rank in its phase | phase median | phase range | gap below the best median | read at |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman | 0 (worse) | 10 | -932.8 | 73.2% | -933.9 | [-981.5, -930.7] | 41.5 | the burned state |
| Friedman | 1 | 9 | -910.8 | 8.0% | -892.4 | [-932.0, -890.3] | 0.0 | the burned state |
| Friedman | 2 (worse) | 11 | -925.2 | 56.7% | -926.2 | [-953.9, -913.4] | 33.8 | the burned state |
| Friedman | 3 (worse) | 9 | -950.4 | 0.0% | -916.0 | [-950.4, -912.4] | 23.6 | its most frequent structure (-916.1) |
| mixed | 0 | 6 | -310.2 | 46.2% | -309.9 | [-326.1, -307.8] | 0.0 | the burned state |
| mixed | 1 | 7 | -332.0 | 5.2% | -312.5 | [-334.7, -308.0] | 2.5 | the burned state |
| mixed | 2 (worse) | 6 | -332.2 | 33.0% | -331.9 | [-358.9, -329.8] | 22.0 | the burned state |
| mixed | 3 | 6 | -317.6 | 10.1% | -309.9 | [-340.2, -306.9] | 0.0 | the burned state |
| functional | 0 | 5 | 180.2 | 27.4% | 181.1 | [163.1, 183.7] | 0.0 | the burned state |
| functional | 1 (worse) | 7 | 152.0 | 31.5% | 152.7 | [136.4, 156.1] | 28.4 | the burned state |
| functional | 2 | 5 | 168.3 | 45.3% | 168.6 | [153.4, 179.3] | 12.5 | the burned state |
| functional | 3 (worse) | 8 | 143.6 | 12.4% | 151.7 | [124.8, 154.4] | 29.4 | the burned state |

- *Friedman chain 3* is the one state outside its phase: its burned state lies below all 2,000 of the phase's states, carrying a just-born inert x6 · x8. P1–P3 read that chain at its most frequent structure, M4's chain 3 top (−916.1, M = 9).
- *The other eleven* lie inside their phases, at ranks 5–73%.
- *"Worse"* marks a chain whose median sits 20 nats or more below the best chain's: Friedman 0, 2 and 3; mixed 2; functional 1 and 3, as the brief names them.

#### P1. Pair deaths

*Per state,* at fixed τ:
- The ratio is the target's ratio, plus ln q_pair, plus ln C(M, 2), with the selection taken as 1.
- The relief is the pair's Δℓ less its cheaper single's.

| fixture | chain | M | best single Δℓ | best pair Δℓ | its ratio | pairs within 10, 20 nats | of C(M, 2) | pairs cheaper than their cheaper single | largest relief | largest Δℓ_pair − Δℓ₁ − Δℓ₂ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman | 0 (worse) | 10 | -49.46 | -81.58 | -76.38 | 0, 0 | 45 | 0 | -15.60 | 251.89 |
| Friedman | 1 | 9 | -9.41 | -138.78 | -134.65 | 0, 0 | 36 | 0 | -11.35 | 238.48 |
| Friedman | 2 (worse) | 11 | -17.71 | -56.68 | -52.14 | 0, 0 | 55 | 0 | -14.18 | 242.55 |
| Friedman | 3 (worse) | 9 | -49.96 | -169.48 | -165.41 | 0, 0 | 36 | 0 | -61.33 | 256.24 |
| mixed | 0 | 6 | -68.75 | -98.40 | -94.58 | 0, 0 | 15 | 0 | -29.65 | 165.18 |
| mixed | 1 | 7 | 2.15 | -62.95 | -59.90 | 0, 0 | 21 | 0 | -35.11 | 174.41 |
| mixed | 2 (worse) | 6 | -17.60 | -69.87 | -66.70 | 0, 0 | 15 | 0 | -52.27 | 154.27 |
| mixed | 3 | 6 | -56.29 | -81.96 | -78.32 | 0, 0 | 15 | 0 | -25.66 | 162.91 |
| functional | 0 | 5 | -329.27 | -353.62 | -350.48 | 0, 0 | 10 | 0 | -24.36 | 346.97 |
| functional | 1 (worse) | 7 | -99.39 | -269.24 | -265.35 | 0, 0 | 21 | 0 | -74.16 | 464.70 |
| functional | 2 | 5 | -310.26 | -346.28 | -343.13 | 0, 0 | 10 | 0 | -36.02 | 333.27 |
| functional | 3 (worse) | 8 | -2.56 | -110.14 | -106.58 | 0, 0 | 28 | 0 | -37.63 | 641.07 |

*Each state's cheapest pair:*

| fixture | chain | cheapest pair | Δℓ | ratio | single Δℓ, each |
| --- | --- | --- | --- | --- | --- |
| Friedman | 0 (worse) | x1−0.31 · x2+0.49 + x1−0.61 · x2−0.48 | -81.58 | -76.38 | -49.46, -56.48 |
| Friedman | 1 | x1−0.76 · x2−0.90 + x2−0.96 | -138.78 | -134.65 | -9.41, -46.95 |
| Friedman | 2 (worse) | x1+0.56 + x1−0.19 · x2+0.49 | -56.68 | -52.14 | -17.71, -20.08 |
| Friedman | 3 (worse) | x1−0.30 · x2+0.25 + x1+0.21 · x2−0.28 | -169.48 | -165.41 | -49.96, -90.14 |
| mixed | 0 | x3+0.71 · z{2} + x3−0.30 · z{2} | -98.40 | -94.58 | -68.75, -92.34 |
| mixed | 1 | x3+0.70 · z{2} + x1−0.62 · x2−0.38 · x3−0.50 | -62.95 | -59.90 | -65.27, 2.15 |
| mixed | 2 (worse) | x2+0.01 · x3−0.28 · z{2} + x2−0.84 · x3−0.45 · z{2} | -69.87 | -66.70 | -54.55, -17.60 |
| mixed | 3 | x3−0.29 + x3−0.27 · z{0,1} | -81.96 | -78.32 | -84.38, -56.29 |
| functional | 0 | x0+0.68 · t+0.37 + x0+0.64 · t−0.74 | -353.62 | -350.48 | -347.84, -329.27 |
| functional | 1 (worse) | x0+0.71 · t−0.53 + x0+0.69 · t+0.37 | -269.24 | -265.35 | -262.71, -195.08 |
| functional | 2 | x0+0.59 · t+0.42 + x0+0.58 · t−0.68 | -346.28 | -343.13 | -331.63, -310.26 |
| functional | 3 (worse) | x0+0.00 · t−0.84 + t+0.79 | -110.14 | -106.58 | -110.36, -2.56 |

- No state's best pair costs under 20 nats, so no triple was read.
- The largest Δℓ_pair − Δℓ₁ − Δℓ₂ belongs to terms sharing one effect.
  - Friedman's two x3 hinges (chain 0: x3+0.40 and x3−0.77): singles −256.4 and −238.7, the pair −254.3, so +240.8.
  - On the functional fit, two x0 main effects (chain 3: x0+0.66 and x0−0.61): +641.

*What P1 shows.*
- *No jointly redundant pair exists in any state.*
  - The cheapest pair costs 57 nats (Friedman chain 2) to 354 (functional chain 0). Its full ratio stands 3–5 nats above its Δℓ.
  - Each of the 307 pairs costs more than its cheaper single; the largest relief is −11.4.
- *The co-adaptation is substitution, not redundancy.* Pairs carry shared structure, up to 641 nats against their singles' sum. Yet removing both always costs more than removing the cheaper one alone: each term does what the other cannot do without it.
- *The worse modes hold no jointly redundant pair the best mode lacks.* Best and worse states read alike on every fixture.
- *Strength: strong.*
  - Every pair of every state is enumerated, and the twelve states agree.
  - The single deaths match the chain's design and the death move's ratio.
  - One fixed state per chain, at fixed τ.

#### P2. Structural moves at fixed M

*The cancellation holds.*
- Over 89,209 proposals, the term prior's ratio and the proposal ratio cancel to at most 1.0 · 10⁻¹⁶, relative to the log priors they difference.
- Through the target, over 78,192 admitted candidates, they cancel to at most 1.4 · 10⁻¹⁶.
- None exceeds 10⁻¹², so each move's log ratio is Δℓ alone, as the draft derived.

*Per state.*
- *Refused:* for support, as a repeated term, by the angle, or as a grow at its part's cap.
- *Unavailable:* a change with no unused variable in its part.

| fixture | chain | M | kind | proposals | unavailable | refused | mean α | largest α | α > 0.01 | largest Δℓ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman | 0 (worse) | 10 | prune | 10 | 0.0% | 0.0% | 1.4E-17 | 8.38E-17 | 0 | -37.02 |
| Friedman | 0 (worse) | 10 | grow | 5000 | 0.0% | 11.3% | 3.14E-08 | 5.85E-05 | 0 | -9.75 |
| Friedman | 0 (worse) | 10 | change | 5000 | 0.0% | 6.5% | 1.4E-22 | 2.16E-19 | 0 | -42.98 |
| Friedman | 1 | 9 | prune | 6 | 0.0% | 0.0% | 3.08E-06 | 1.41E-05 | 0 | -11.17 |
| Friedman | 1 | 9 | grow | 5000 | 0.0% | 5.9% | 4.03E-06 | 0.00135 | 0 | -6.61 |
| Friedman | 1 | 9 | change | 5000 | 0.0% | 4.7% | 9.33E-07 | 0.000523 | 0 | -7.56 |
| Friedman | 2 (worse) | 11 | prune | 6 | 0.0% | 0.0% | 1.96E-07 | 1.18E-06 | 0 | -13.65 |
| Friedman | 2 (worse) | 11 | grow | 5000 | 0.0% | 12.1% | 0.000138 | 0.132 | 12 | -2.02 |
| Friedman | 2 (worse) | 11 | change | 5000 | 0.0% | 5.2% | 7.08E-10 | 3.87E-07 | 0 | -14.76 |
| Friedman | 3 (worse) | 9 | prune | 6 | 0.0% | 0.0% | 8.5E-19 | 5.1E-18 | 0 | -39.82 |
| Friedman | 3 (worse) | 9 | grow | 5000 | 0.0% | 10.0% | 7.72E-12 | 4.56E-09 | 0 | -19.21 |
| Friedman | 3 (worse) | 9 | change | 5000 | 0.0% | 5.7% | 1.19E-20 | 6.62E-18 | 0 | -39.56 |
| mixed | 0 | 6 | prune, part emptied | 8 | 0.0% | 50.0% | 1.4E-28 | 1.12E-27 | 0 | -62.05 |
| mixed | 0 | 6 | grow | 5000 | 0.0% | 55.1% | 0.000517 | 0.793 | 4 | -0.23 |
| mixed | 0 | 6 | change | 5000 | 50.2% | 6.3% | 1.08E-32 | 8.16E-31 | 0 | -69.28 |
| mixed | 1 | 7 | prune, part emptied | 8 | 0.0% | 50.0% | 5.25E-27 | 4.2E-26 | 0 | -58.43 |
| mixed | 1 | 7 | prune | 3 | 0.0% | 0.0% | 0.507 | 0.517 | 3 | -0.66 |
| mixed | 1 | 7 | grow | 5000 | 0.0% | 58.0% | 0.0286 | 1 | 178 | 1.11 |
| mixed | 1 | 7 | change | 5000 | 56.3% | 5.6% | 2.43E-31 | 1.81E-29 | 0 | -66.18 |
| mixed | 2 (worse) | 6 | prune, part emptied | 8 | 0.0% | 0.0% | 1.17E-05 | 9.38E-05 | 0 | -9.27 |
| mixed | 2 (worse) | 6 | prune | 4 | 0.0% | 0.0% | 0.0144 | 0.0321 | 2 | -3.44 |
| mixed | 2 (worse) | 6 | grow | 5000 | 0.0% | 58.6% | 6.61E-05 | 0.0089 | 0 | -4.72 |
| mixed | 2 (worse) | 6 | change | 5000 | 35.7% | 12.4% | 6.91E-08 | 1.89E-05 | 0 | -10.88 |
| mixed | 3 | 6 | prune, part emptied | 8 | 0.0% | 0.0% | 1.05E-19 | 8.44E-19 | 0 | -41.62 |
| mixed | 3 | 6 | grow | 5000 | 0.0% | 47.7% | 0.00128 | 0.117 | 55 | -2.15 |
| mixed | 3 | 6 | change | 5000 | 33.1% | 7.0% | 2.35E-21 | 8.57E-19 | 0 | -41.60 |
| functional | 0 | 5 | prune, part emptied | 6 | 0.0% | 0.0% | 1.39E-140 | 8.36E-140 | 0 | -320.24 |
| functional | 0 | 5 | grow | 5000 | 0.0% | 42.3% | 1.03E-23 | 4.29E-21 | 0 | -46.90 |
| functional | 0 | 5 | change | 5000 | 31.5% | 6.8% | 4.39E-143 | 5.34E-141 | 0 | -322.99 |
| functional | 1 (worse) | 7 | prune, part emptied | 8 | 0.0% | 12.5% | 0.00596 | 0.0477 | 1 | -3.04 |
| functional | 1 (worse) | 7 | grow | 5000 | 0.0% | 48.1% | 1.05E-14 | 2.14E-11 | 0 | -24.57 |
| functional | 1 (worse) | 7 | change | 5000 | 41.6% | 5.5% | 3.27E-15 | 5.81E-12 | 0 | -25.87 |
| functional | 2 | 5 | prune, part emptied | 6 | 0.0% | 0.0% | 9.15E-129 | 5.49E-128 | 0 | -293.03 |
| functional | 2 | 5 | grow | 5000 | 0.0% | 39.5% | 1.59E-21 | 7.94E-19 | 0 | -41.68 |
| functional | 2 | 5 | change | 5000 | 30.5% | 6.1% | 1.95E-132 | 4.59E-130 | 0 | -297.81 |
| functional | 3 (worse) | 8 | prune, part emptied | 10 | 0.0% | 0.0% | 0.00126 | 0.0105 | 1 | -4.56 |
| functional | 3 (worse) | 8 | grow | 5000 | 0.0% | 48.3% | 0.00802 | 1 | 166 | 0.19 |
| functional | 3 (worse) | 8 | change | 5000 | 44.6% | 5.4% | 3.58E-05 | 0.0377 | 5 | -3.28 |

*Pooled, by kind and part:*

| fixture | kind, part | proposals | unavailable | refused: support | repeated | angle | cap | mean α | Δℓ 5% | 25% | median | 75% | 95% | max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman | change, continuous | 20000 | 0.0% | 5.5% | 0.0% | 0.0% | 0.0% | 2.33E-07 | -517.33 | -259.80 | -239.40 | -71.31 | -21.12 | -7.56 |
| Friedman | grow, continuous | 20000 | 0.0% | 9.8% | 0.0% | 0.0% | 0.0% | 3.56E-05 | -476.40 | -246.96 | -199.85 | -69.19 | -16.37 | -2.02 |
| Friedman | prune, continuous | 28 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 7.03E-07 | -232.86 | -174.43 | -78.68 | -44.23 | -12.79 | -11.17 |
| mixed | change, categorical | 8039 | 100.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0 | — | — | — | — | — | — |
| mixed | change, continuous | 11961 | 6.1% | 13.1% | 0.0% | 0.0% | 0.0% | 2.89E-08 | -248.51 | -143.26 | -129.52 | -87.56 | -37.93 | -10.88 |
| mixed | grow, categorical | 10054 | 0.0% | 1.7% | 0.0% | 0.0% | 76.5% | 0.0149 | -125.20 | -107.85 | -87.61 | -66.98 | -0.68 | 1.11 |
| mixed | grow, continuous | 9946 | 0.0% | 27.5% | 0.0% | 0.0% | 3.8% | 0.000297 | -259.36 | -161.71 | -114.16 | -74.02 | -26.84 | -0.23 |
| mixed | prune, continuous | 7 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.225 | -32.28 | -6.71 | -3.44 | -0.69 | -0.67 | -0.66 |
| mixed | prune, part emptied, categorical | 17 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 5.52E-06 | -268.00 | -102.11 | -80.15 | -61.61 | -24.79 | -9.27 |
| mixed | prune, part emptied, continuous | 15 | 0.0% | 0.0% | 26.7% | 26.7% | 0.0% | 5.63E-20 | -176.44 | -154.96 | -144.58 | -92.41 | -56.85 | -41.62 |
| functional | change, continuous | 12591 | 0.0% | 9.5% | 0.0% | 0.0% | 0.0% | 1.42E-05 | -882.08 | -705.84 | -507.15 | -326.65 | -113.39 | -3.28 |
| functional | change, functional | 7409 | 100.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0 | — | — | — | — | — | — |
| functional | grow, continuous | 10028 | 0.0% | 19.2% | 0.0% | 0.0% | 0.0% | 0.004 | -817.67 | -553.78 | -376.54 | -229.46 | -60.75 | 0.19 |
| functional | grow, functional | 9972 | 0.0% | 3.3% | 0.1% | 0.0% | 66.7% | 3.27E-19 | -1113.90 | -585.73 | -485.11 | -295.76 | -91.78 | -35.56 |
| functional | prune, part emptied, continuous | 15 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.00402 | -580.29 | -451.98 | -293.03 | -130.77 | -4.11 | -3.04 |
| functional | prune, part emptied, functional | 15 | 0.0% | 0.0% | 6.7% | 0.0% | 0.0% | 7.52E-47 | -524.56 | -333.43 | -285.54 | -163.02 | -110.05 | -103.50 |

*Friedman, by the factor's variable:*

| fixture | kind, the factor's variable | proposals | unavailable | refused: support | repeated | angle | cap | mean α | Δℓ 5% | 25% | median | 75% | 95% | max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman | change, active → active | 8539 | 0.0% | 5.2% | 0.0% | 0.0% | 0.0% | 1.17E-07 | -521.94 | -260.24 | -239.67 | -93.07 | -21.26 | -10.24 |
| Friedman | change, active → inert | 11461 | 0.0% | 5.8% | 0.0% | 0.0% | 0.0% | 3.2E-07 | -517.32 | -259.31 | -238.42 | -71.02 | -20.89 | -7.56 |
| Friedman | grow, active | 8423 | 0.0% | 9.1% | 0.0% | 0.0% | 0.0% | 7.95E-05 | -456.28 | -244.72 | -191.66 | -70.42 | -16.32 | -2.02 |
| Friedman | grow, inert | 11577 | 0.0% | 10.3% | 0.0% | 0.0% | 0.0% | 3.64E-06 | -482.17 | -248.44 | -206.75 | -68.21 | -16.49 | -4.83 |
| Friedman | prune, active | 28 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 7.03E-07 | -232.86 | -174.43 | -78.68 | -44.23 | -12.79 | -11.17 |

*Pooled, by the term's orders before and after.* Rows wholly unavailable or at the cap are left out.
- The mixed fit's categorical part and the functional fit's output part each have one variable, so a change there never has an unused variable.
- A grow into a part already at its cap of 1 is refused.

| fixture | kind, orders before → after | proposals | unavailable | refused: support | repeated | angle | cap | mean α | Δℓ 5% | 25% | median | 75% | 95% | max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman | change, K1 → K1 | 12752 | 0.0% | 3.8% | 0.0% | 0.0% | 0.0% | 4.65E-11 | -523.96 | -277.45 | -249.72 | -226.08 | -49.43 | -18.27 |
| Friedman | change, K2 → K2 | 7248 | 0.0% | 8.6% | 0.0% | 0.0% | 0.0% | 6.44E-07 | -281.78 | -221.63 | -93.71 | -52.96 | -12.61 | -7.56 |
| Friedman | grow, K1 → K2 | 12662 | 0.0% | 6.1% | 0.0% | 0.0% | 0.0% | 4.87E-05 | -495.51 | -266.08 | -233.42 | -166.29 | -47.25 | -2.02 |
| Friedman | grow, K2 → K3 | 7338 | 0.0% | 16.2% | 0.0% | 0.0% | 0.0% | 1.29E-05 | -234.54 | -179.55 | -74.99 | -43.84 | -11.91 | -4.16 |
| Friedman | prune, K2 → K1 | 28 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 7.03E-07 | -232.86 | -174.43 | -78.68 | -44.23 | -12.79 | -11.17 |
| mixed | change, K1 → K1 | 4095 | 0.0% | 6.4% | 0.0% | 0.0% | 0.0% | 2.14E-36 | -133.06 | -131.48 | -129.94 | -127.92 | -81.45 | -78.12 |
| mixed | change, K1Z1 → K1Z1 | 6010 | 0.0% | 15.7% | 0.0% | 0.0% | 0.0% | 1.96E-21 | -262.32 | -172.41 | -143.21 | -94.89 | -56.76 | -41.60 |
| mixed | change, K2Z1 → K2Z1 | 1128 | 0.0% | 31.8% | 0.0% | 0.0% | 0.0% | 3.06E-07 | -54.57 | -49.58 | -26.33 | -15.60 | -12.90 | -10.88 |
| mixed | grow, K1 → K1Z1 | 2015 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.00319 | -125.20 | -107.85 | -92.64 | -73.94 | -26.37 | -2.15 |
| mixed | grow, K1 → K2 | 2029 | 0.0% | 10.3% | 0.0% | 0.0% | 0.0% | 4.7E-22 | -128.28 | -120.89 | -104.78 | -80.21 | -60.55 | -42.81 |
| mixed | grow, K1Z1 → K2Z1 | 5903 | 0.0% | 34.0% | 0.0% | 0.0% | 0.0% | 1.26E-14 | -255.41 | -168.92 | -124.31 | -77.92 | -39.58 | -25.19 |
| mixed | grow, K2Z1 → K3Z1 | 850 | 0.0% | 48.2% | 0.0% | 0.0% | 0.0% | 0.000389 | -35.61 | -29.83 | -19.77 | -8.56 | -5.23 | -4.72 |
| mixed | grow, K3 → K3Z1 | 350 | 0.0% | 49.7% | 0.0% | 0.0% | 0.0% | 0.409 | -0.68 | -0.68 | 0.31 | 1.11 | 1.11 | 1.11 |
| mixed | grow, Z1 → K1Z1 | 791 | 0.0% | 13.1% | 0.0% | 0.0% | 0.0% | 0.00332 | -265.89 | -260.06 | -251.39 | -229.45 | -27.65 | -0.23 |
| mixed | prune, K2Z1 → K1Z1 | 4 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0144 | -37.11 | -17.79 | -6.71 | -3.61 | -3.47 | -3.44 |
| mixed | prune, K3 → K2 | 3 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.507 | -0.69 | -0.69 | -0.68 | -0.67 | -0.66 | -0.66 |
| mixed | prune, part emptied, K1Z1 → K1 | 15 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 3.89E-27 | -269.10 | -131.14 | -80.86 | -68.18 | -59.23 | -58.43 |
| mixed | prune, part emptied, K1Z1 → Z1 | 15 | 0.0% | 0.0% | 26.7% | 26.7% | 0.0% | 5.63E-20 | -176.44 | -154.96 | -144.58 | -92.41 | -56.85 | -41.62 |
| mixed | prune, part emptied, K2Z1 → K2 | 2 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 4.69E-05 | -27.70 | -23.82 | -18.97 | -14.12 | -10.24 | -9.27 |
| functional | change, K1 → K1 | 6631 | 0.0% | 9.4% | 0.0% | 0.0% | 0.0% | 3.6E-216 | -1357.06 | -800.87 | -593.96 | -506.10 | -495.28 | -491.71 |
| functional | change, K1T1 → K1T1 | 5960 | 0.0% | 9.5% | 0.0% | 0.0% | 0.0% | 3E-05 | -740.83 | -566.22 | -323.98 | -188.39 | -87.44 | -3.28 |
| functional | grow, K1 → K1T1 | 3322 | 0.0% | 10.0% | 0.3% | 0.0% | 0.0% | 9.81E-19 | -1113.90 | -585.73 | -485.11 | -295.76 | -91.78 | -35.56 |
| functional | grow, K1 → K2 | 3312 | 0.0% | 19.2% | 0.0% | 0.0% | 0.0% | 3.53E-142 | -1211.94 | -678.44 | -557.56 | -498.05 | -418.51 | -317.60 |
| functional | grow, K1T1 → K2T1 | 6032 | 0.0% | 19.8% | 0.0% | 0.0% | 0.0% | 8.67E-15 | -678.81 | -350.72 | -285.41 | -157.37 | -77.45 | -24.57 |
| functional | grow, T1 → K1T1 | 684 | 0.0% | 13.2% | 0.0% | 0.0% | 0.0% | 0.0586 | -564.25 | -513.34 | -238.74 | -4.46 | -0.31 | 0.19 |
| functional | prune, part emptied, K1T1 → K1 | 15 | 0.0% | 0.0% | 6.7% | 0.0% | 0.0% | 7.52E-47 | -524.56 | -333.43 | -285.54 | -163.02 | -110.05 | -103.50 |
| functional | prune, part emptied, K1T1 → T1 | 15 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.00402 | -580.29 | -451.98 | -293.03 | -130.77 | -4.11 | -3.04 |

*The worse states' ten cheapest exits, by family:* the kind, the old term and the new term's variables, each family by its best member.

| fixture, chain | exit | kind | Δℓ | ratio | α | the family's proposals | their mean α |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman 0 (worse) | x1−0.58 → x1−0.58 · x2+0.03 | grow | -9.75 | -9.75 | 5.85E-05 | 56 | 2.8E-06 |
| Friedman 0 (worse) | x1−0.31 · x2+0.49 → x1−0.31 · x2+0.49 · x9+0.06 | grow | -20.17 | -20.17 | 1.74E-09 | 62 | 5.11E-11 |
| Friedman 0 (worse) | x1−0.31 · x2+0.49 → x1−0.31 · x2+0.49 · x5−0.94 | grow | -21.31 | -21.31 | 5.57E-10 | 50 | 2.2E-11 |
| Friedman 0 (worse) | x1−0.31 · x2+0.49 → x1−0.31 · x2+0.49 · x4+0.02 | grow | -23.83 | -23.83 | 4.49E-11 | 60 | 1.03E-12 |
| Friedman 0 (worse) | x1−0.31 · x2+0.49 → x1−0.31 · x2+0.49 · x10+0.00 | grow | -24.51 | -24.51 | 2.27E-11 | 65 | 1.06E-12 |
| Friedman 0 (worse) | x1−0.31 · x2+0.49 → x1−0.31 · x2+0.49 · x6+0.01 | grow | -24.57 | -24.57 | 2.13E-11 | 69 | 7.73E-13 |
| Friedman 0 (worse) | x1−0.31 · x2+0.49 → x1−0.31 · x2+0.49 · x3−0.98 | grow | -26.21 | -26.21 | 4.14E-12 | 62 | 9.7E-14 |
| Friedman 0 (worse) | x1−0.31 · x2+0.49 → x1−0.31 · x2+0.49 · x7+0.01 | grow | -26.31 | -26.31 | 3.75E-12 | 59 | 1.79E-13 |
| Friedman 0 (worse) | x1−0.31 · x2+0.49 → x1−0.31 · x2+0.49 · x8+0.06 | grow | -26.68 | -26.68 | 2.59E-12 | 55 | 8E-14 |
| Friedman 0 (worse) | x1−0.58 → x1−0.58 · x5+0.03 | grow | -27.52 | -27.52 | 1.12E-12 | 55 | 2.03E-14 |
| Friedman 2 (worse) | x1+0.56 → x1+0.56 · x2−0.98 | grow | -2.02 | -2.02 | 0.132 | 47 | 0.00696 |
| Friedman 2 (worse) | x2−0.32 → x1+0.01 · x2−0.32 | grow | -2.30 | -2.30 | 0.1 | 42 | 0.00686 |
| Friedman 2 (worse) | x1−0.19 · x2+0.49 → x1−0.19 · x2+0.49 · x5−0.99 | grow | -4.16 | -4.16 | 0.0156 | 58 | 0.000554 |
| Friedman 2 (worse) | x1−0.19 · x2+0.49 → x1−0.19 · x2+0.49 · x9+0.01 | grow | -4.83 | -4.83 | 0.00796 | 58 | 0.000586 |
| Friedman 2 (worse) | x1−0.19 · x2+0.49 → x1−0.19 · x2+0.49 · x3−1.00 | grow | -7.06 | -7.06 | 0.000857 | 47 | 2.97E-05 |
| Friedman 2 (worse) | x1−0.19 · x2+0.49 → x1−0.19 · x2+0.49 · x10+0.01 | grow | -7.21 | -7.21 | 0.000739 | 55 | 2.75E-05 |
| Friedman 2 (worse) | x1−0.19 · x2+0.49 → x1−0.19 · x2+0.49 · x4−0.99 | grow | -7.29 | -7.29 | 0.000682 | 55 | 2.1E-05 |
| Friedman 2 (worse) | x1+0.56 → x1+0.56 · x3+0.10 | grow | -7.56 | -7.56 | 0.000521 | 51 | 2.96E-05 |
| Friedman 2 (worse) | x1−0.19 · x2+0.49 → x1−0.19 · x2+0.49 · x7−0.94 | grow | -7.73 | -7.73 | 0.000439 | 65 | 2.61E-05 |
| Friedman 2 (worse) | x1−0.19 · x2+0.49 → x1−0.19 · x2+0.49 · x6+0.04 | grow | -8.11 | -8.11 | 0.0003 | 61 | 2.68E-05 |
| Friedman 3 (worse) | x1−0.30 · x2+0.25 → x1−0.30 · x2+0.25 · x6+0.04 | grow | -19.21 | -19.21 | 4.56E-09 | 64 | 1.9E-10 |
| Friedman 3 (worse) | x1−0.30 · x2+0.25 → x1−0.30 · x2+0.25 · x9+0.05 | grow | -19.32 | -19.32 | 4.07E-09 | 73 | 1.86E-10 |
| Friedman 3 (worse) | x1−0.30 · x2+0.25 → x1−0.30 · x2+0.25 · x10+0.01 | grow | -19.68 | -19.68 | 2.83E-09 | 79 | 7.39E-11 |
| Friedman 3 (worse) | x1−0.30 · x2+0.25 → x1−0.30 · x2+0.25 · x7−1.00 | grow | -19.84 | -19.84 | 2.42E-09 | 60 | 6.66E-11 |
| Friedman 3 (worse) | x1−0.30 · x2+0.25 → x1−0.30 · x2+0.25 · x5−0.98 | grow | -20.76 | -20.76 | 9.68E-10 | 70 | 1.94E-11 |
| Friedman 3 (worse) | x1−0.30 · x2+0.25 → x1−0.30 · x2+0.25 · x4+0.05 | grow | -21.12 | -21.12 | 6.76E-10 | 76 | 1.15E-11 |
| Friedman 3 (worse) | x1+0.21 · x2−0.28 → x1+0.21 · x2−0.28 · x5+0.02 | grow | -21.72 | -21.72 | 3.69E-10 | 73 | 5.51E-12 |
| Friedman 3 (worse) | x1+0.21 · x2−0.28 → x1+0.21 · x2−0.28 · x4+0.02 | grow | -23.01 | -23.01 | 1.02E-10 | 73 | 3.98E-12 |
| Friedman 3 (worse) | x1−0.30 · x2+0.25 → x1−0.30 · x2+0.25 · x8−0.96 | grow | -24.38 | -24.38 | 2.58E-11 | 70 | 8.38E-13 |
| Friedman 3 (worse) | x1−0.30 · x2+0.25 → x1−0.30 · x2+0.25 · x3−0.91 | grow | -25.30 | -25.30 | 1.03E-11 | 57 | 2.92E-13 |
| mixed 2 (worse) | x2−0.84 · x3−0.45 · z{2} → x3−0.45 · z{2} | prune | -3.44 | -3.44 | 0.0321 | 1 | 0.0321 |
| mixed 2 (worse) | x2+0.01 · x3−0.28 · z{2} → x3−0.28 · z{2} | prune | -3.67 | -3.67 | 0.0256 | 1 | 0.0256 |
| mixed 2 (worse) | x2−0.84 · x3−0.45 · z{2} → x1−1.00 · x2−0.84 · x3−0.45 · z{2} | grow | -4.72 | -4.72 | 0.0089 | 397 | 0.000833 |
| mixed 2 (worse) | x2−0.84 · x3−0.45 · z{2} → x2−0.84 · x3−0.45 | prune, part emptied | -9.27 | -9.27 | 9.38E-05 | 1 | 9.38E-05 |
| mixed 2 (worse) | x2−0.84 · x3−0.45 · z{2} → x2−0.84 · z{2} | prune | -9.75 | -9.75 | 5.85E-05 | 1 | 5.85E-05 |
| mixed 2 (worse) | x2−0.84 · x3−0.45 · z{2} → x1+0.00 · x3−0.45 · z{2} | change | -10.88 | -10.88 | 1.89E-05 | 277 | 7.21E-07 |
| mixed 2 (worse) | x2−0.84 · x3−0.45 · z{2} → x1+0.00 · x2−0.84 · z{2} | change | -11.95 | -11.95 | 6.49E-06 | 274 | 5.32E-07 |
| mixed 2 (worse) | x2+0.01 · x3−0.28 · z{2} → x1+0.00 · x2+0.01 · x3−0.28 · z{2} | grow | -17.37 | -17.37 | 2.86E-08 | 453 | 3.4E-10 |
| mixed 2 (worse) | x2+0.01 · x3−0.28 · z{2} → x1+0.00 · x3−0.28 · z{2} | change | -26.16 | -26.16 | 4.34E-12 | 298 | 9.08E-14 |
| mixed 2 (worse) | x2+0.01 · x3−0.28 · z{2} → x2+0.01 · x3−0.28 | prune, part emptied | -28.67 | -28.67 | 3.53E-13 | 1 | 3.53E-13 |
| functional 1 (worse) | x0−0.95 · t+0.74 → t+0.74 | prune, part emptied | -3.04 | -3.04 | 0.0477 | 1 | 0.0477 |
| functional 1 (worse) | x0−0.95 · t+0.74 → x0−0.95 · x1+0.00 · t+0.74 | grow | -24.57 | -24.57 | 2.14E-11 | 354 | 1.48E-13 |
| functional 1 (worse) | x0−0.95 · t+0.74 → x1+0.00 · t+0.74 | change | -25.87 | -25.87 | 5.81E-12 | 370 | 4.41E-14 |
| functional 1 (worse) | x0+0.69 → x0+0.69 · t+0.00 | grow | -35.56 | -35.56 | 3.61E-16 | 339 | 9.59E-18 |
| functional 1 (worse) | x0+0.69 · t+0.37 → x0+0.69 · x1−1.00 · t+0.37 | grow | -94.37 | -94.37 | 1.04E-41 | 368 | 6.27E-44 |
| functional 1 (worse) | x0−0.95 · t+0.74 → x0−0.95 | prune, part emptied | -103.50 | -103.50 | 1.13E-45 | 1 | 1.13E-45 |
| functional 1 (worse) | x0+0.71 · t−0.53 → x0+0.71 · x1+0.00 · t−0.53 | grow | -125.64 | -125.64 | 2.71E-55 | 370 | 1.48E-57 |
| functional 1 (worse) | t−0.74 → x0−0.98 · t−0.74 | grow | -135.24 | -135.24 | 1.85E-59 | 178 | 1.29E-61 |
| functional 1 (worse) | x0+0.69 · t+0.37 → t+0.37 | prune, part emptied | -135.42 | -135.42 | 1.54E-59 | 1 | 1.54E-59 |
| functional 1 (worse) | x0+0.69 · t+0.37 → x1+0.00 · t+0.37 | change | -153.25 | -153.25 | 2.78E-67 | 376 | 3.18E-69 |
| functional 3 (worse) | t+0.79 → x1−0.98 · t+0.79 | grow | 0.19 | 0.19 | 1 | 171 | 0.149 |
| functional 3 (worse) | t+0.79 → x0−0.90 · t+0.79 | grow | 0.00 | 0.00 | 1 | 156 | 0.0937 |
| functional 3 (worse) | x0−1.00 · t+0.05 → x1+0.00 · t+0.05 | change | -3.28 | -3.28 | 0.0377 | 333 | 0.000537 |
| functional 3 (worse) | x0−1.00 · t+0.05 → t+0.05 | prune, part emptied | -4.56 | -4.56 | 0.0105 | 1 | 0.0105 |
| functional 3 (worse) | x0+0.00 · t−0.84 → t−0.84 | prune, part emptied | -6.15 | -6.15 | 0.00214 | 1 | 0.00214 |
| functional 3 (worse) | x0+0.76 · t+0.63 → x0+0.76 · x1+0.00 · t+0.63 | grow | -53.26 | -53.26 | 7.42E-24 | 296 | 1.56E-25 |
| functional 3 (worse) | x0+0.64 · t−0.37 → x0+0.64 · x1+0.00 · t−0.37 | grow | -66.56 | -66.56 | 1.23E-29 | 313 | 1.22E-31 |
| functional 3 (worse) | x0+0.00 · t−0.84 → x0+0.00 · x1−1.00 · t−0.84 | grow | -82.68 | -82.68 | 1.23E-36 | 325 | 1.62E-38 |
| functional 3 (worse) | x0−0.61 → x0−0.61 · t−1.00 | grow | -87.12 | -87.12 | 1.46E-38 | 337 | 2.6E-40 |
| functional 3 (worse) | x0−1.00 · t+0.05 → x0−1.00 · x1+0.00 · t+0.05 | grow | -96.78 | -96.78 | 9.35E-43 | 290 | 3.22E-45 |

*The two-step profile.* For each distinct state reached by a prune or change costing under 20 nats, its best removal of the other terms. Removing the re-represented term itself reaches the fixed state's own single death of the old term, so it is left out.

| fixture | chain | states reached under 20 nats (prunes, changes) | best other single there | at the fixed state | best total, step + single | fixed state's best single | best other pair there | at the fixed state | best total, step + pair | fixed state's best pair | totals over −10 (single, pair) | over −20 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman | 0 (worse) | 0 | — | — | — | -49.46 | — | — | — | -81.58 | — | — |
| Friedman | 1 | 492 (2, 490) | 1.06 | -46.95 | -11.27 | -9.41 | -120.42 | -246.90 | -131.67 | -138.78 | 0, 0 | 1, 0 |
| Friedman | 2 (worse) | 128 (1, 127) | -2.98 | -243.71 | -21.56 | -17.71 | -39.51 | -245.72 | -58.09 | -56.68 | 0, 0 | 0, 0 |
| Friedman | 3 (worse) | 0 | — | — | — | -49.96 | — | — | — | -169.48 | — | — |
| mixed | 0 | 0 | — | — | — | -68.75 | — | — | — | -98.40 | — | — |
| mixed | 1 | 3 (3, 0) | -64.82 | -65.27 | -65.48 | 2.15 | -99.02 | -100.38 | -99.70 | -62.95 | 0, 0 | 0, 0 |
| mixed | 2 (worse) | 303 (4, 299) | 2.49 | -17.60 | -1.18 | -17.60 | -103.71 | -112.04 | -107.38 | -69.87 | 2, 0 | 2, 0 |
| mixed | 3 | 0 | — | — | — | -56.29 | — | — | — | -81.96 | — | — |
| functional | 0 | 0 | — | — | — | -329.27 | — | — | — | -353.62 | — | — |
| functional | 1 (worse) | 1 (1, 0) | -265.29 | -195.08 | -268.33 | -99.39 | -306.48 | -269.24 | -309.52 | -269.24 | 0, 0 | 0, 0 |
| functional | 2 | 0 | — | — | — | -310.26 | — | — | — | -346.28 | — | — |
| functional | 3 (worse) | 12 (2, 10) | 3.86 | -110.36 | 0.55 | -2.56 | 7.67 | -573.85 | 4.27 | -110.14 | 10, 10 | 12, 11 |

*The best paths,* three per state that has any:

| fixture, chain | step | step Δℓ | then the death | its Δℓ | the same at the fixed state | total Δℓ | the death's ratio | end state's log post less the fixed state's |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman 1 | prune: x1−0.76 · x2−0.90 → x2−0.90 | -12.32 | x2−0.96 | 1.06 | -46.95 | -11.27 | 3.12 | 6.96 |
| Friedman 1 | prune: x1−0.76 · x2−0.90 → x1−0.76 | -11.17 | x1−0.90 | -32.61 | -264.64 | -43.79 | -30.55 | -25.56 |
| Friedman 1 | change: x1−0.76 · x2−0.90 → x2−0.90 · x5+0.03 | -11.64 | x2−0.96 | -52.22 | -46.95 | -63.85 | -50.25 | -54.04 |
| Friedman 2 (worse) | change: x1+0.56 → x3−0.84 | -18.58 | x3−0.71 | -2.98 | -243.71 | -21.56 | -0.82 | -11.85 |
| Friedman 2 (worse) | change: x1+0.56 → x3−0.85 | -18.79 | x3−0.71 | -5.69 | -243.71 | -24.48 | -3.53 | -14.77 |
| Friedman 2 (worse) | change: x1+0.56 → x3−0.86 | -18.81 | x3−0.71 | -5.89 | -243.71 | -24.69 | -3.72 | -14.98 |
| mixed 1 | prune: x1−0.62 · x2−0.38 · x3−0.50 → x2−0.38 · x3−0.50 | -0.66 | x3+0.70 · z{2} | -64.82 | -65.27 | -65.48 | -62.80 | -49.19 |
| mixed 1 | prune: x1−0.62 · x2−0.38 · x3−0.50 → x1−0.62 · x3−0.50 | -0.68 | x3+0.70 · z{2} | -64.95 | -65.27 | -65.63 | -62.94 | -49.34 |
| mixed 1 | prune: x1−0.62 · x2−0.38 · x3−0.50 → x1−0.62 · x2−0.38 | -0.70 | x3+0.70 · z{2} | -65.05 | -65.27 | -65.75 | -63.04 | -49.46 |
| mixed 2 (worse) | prune: x2+0.01 · x3−0.28 · z{2} → x3−0.28 · z{2} | -3.67 | x2−0.84 · x3−0.45 · z{2} | 2.49 | -17.60 | -1.18 | 4.04 | 22.67 |
| mixed 2 (worse) | prune: x2−0.84 · x3−0.45 · z{2} → x3−0.45 · z{2} | -3.44 | x2+0.01 · x3−0.28 · z{2} | -4.10 | -54.55 | -7.54 | -2.55 | 16.31 |
| mixed 2 (worse) | change: x2−0.84 · x3−0.45 · z{2} → x1+0.00 · x3−0.45 · z{2} | -10.88 | x2+0.01 · x3−0.28 · z{2} | -18.59 | -54.55 | -29.47 | -16.98 | -12.01 |
| functional 1 (worse) | prune, part emptied: x0−0.95 · t+0.74 → t+0.74 | -3.04 | x0+0.69 · t+0.37 | -265.29 | -195.08 | -268.33 | -263.35 | -252.26 |
| functional 1 (worse) | prune, part emptied: x0−0.95 · t+0.74 → t+0.74 | -3.04 | x0+0.71 · t−0.53 + x0+0.69 · t+0.37 | -306.48 | -269.24 | -309.52 | -302.81 | -282.80 |
| functional 3 (worse) | change: x0−1.00 · t+0.05 → x1+0.02 · t+0.05 | -3.40 | x1−0.92 · t+0.11 + t+0.79 | 7.67 | -573.85 | 4.27 | 11.23 | 19.90 |
| functional 3 (worse) | change: x0−1.00 · t+0.05 → x1+0.00 · t+0.05 | -3.28 | x0+0.00 · t−0.84 + t+0.79 | 7.49 | -110.14 | 4.21 | 11.06 | 19.85 |
| functional 3 (worse) | change: x0−1.00 · t+0.05 → x1+0.03 · t+0.05 | -4.12 | x1−0.92 · t+0.11 + x0+0.00 · t−0.84 | 7.59 | -591.27 | 3.47 | 11.85 | 24.58 |

*What P2 shows.*
- *On Friedman's worse modes 0 and 3, no fixed-M move is a cheap exit.*
  - Their cheapest prune, grow or change costs 9.8 and 19.2 nats (α ≤ 6 · 10⁻⁵).
  - No prune or change costs under 20 nats there, so their two-step profiles are empty.
  - Chain 2's cheapest exits are grows onto edge hinges: a knot at its variable's end, pointing inward, so the new factor is nearly linear. x1+0.56 → x1+0.56 · x2−0.98 accepts 0.13.
  - Friedman's 28 prunes (all K2 → K1) cost 11.2 nats or more.
- *On the mixed and functional fits, each fixture has one worse state with a near-free two-step exit.*
  - *Mixed chain 2.*
    - Pruning x2 from x2+0.01 · x3−0.28 · z{2} costs 3.7 nats (α 0.026). Its twin x2−0.84 · x3−0.45 · z{2} then dies with Δℓ +2.5.
    - The path costs 1.2 nats in all and ends 22.7 above the fixed state in log posterior: the other chains' level, since chain 2 sits 22.0 below them.
  - *Functional chain 3.*
    - Changing x0−1.00 · t+0.05 to x1 · t+0.05, its knot near x1's low end, costs 3.3–3.4 nats (α up to 0.038). On the path through x1+0.02 · t+0.05, x1−0.92 · t+0.11 and t+0.79 then die together with Δℓ +7.7.
    - The path gains 4.3 and ends 19.9 above the fixed state: 5 below chain 2's median and 18 below chain 0's.
    - Ten of its twelve reached states have totals above −10.
  - In both, the second step is a death only the first opens. At the fixed state the same death costs 17.6 (mixed), and the same pair 573.9 (functional).
- *On the better states, the cheap re-representations open nothing.* Friedman chain 1 reaches 492 states under 20 nats, 490 of them by variable changes. Their best two-step total is −11.3, against that state's own best single death at −9.4.
- *Grows of a near-surplus term accept freely.*
  - Mixed chain 1's x1 · x2 · x3 (death Δℓ +2.2) grows a z factor at a mean α of 0.41.
  - Functional chain 3's t+0.79 (death −2.6) grows an x factor at up to α = 1.
- *Active and inert inputs, Friedman:* grows onto an active variable accept 8.0 · 10⁻⁵ on average, against 3.6 · 10⁻⁶ for an inert one. Changes accept alike (1.2 · 10⁻⁷ and 3.2 · 10⁻⁷).
- *Strength.*
  - The cancellation is exact, and the two-step paths are exact computations.
  - Each exit rests on one state. Mixed chain 2 is that fixture's only worse chain. Functional chain 1, the other worse state there, has none: its one cheap prune opens only deaths costing 265 nats or more.
  - The acceptance readings pool 5,000 draws per state and agree across chains.

#### P3. Profile deaths

*Per state.*
- The profile Δℓ is measured against the fixed state, and "within" counts Δℓ above −10 and −20.
- The control is the same search on the full state, with no term removed.

| fixture | chain | M | best plain Δℓ | best profile Δℓ (5 passes) | plain within 10, 20 | profile within 10, 20 | searches stopped at 5 passes | run on to 50: best Δℓ, within 10, 20, still at the limit | full state's own profile, 5 passes | run on (passes) | best death less the full profile, run on | deaths within 10, 20 of the full profile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman | 0 (worse) | 10 | -49.46 | 4.91 | 0, 0 | 1, 2 | 10 | 6.44, 4, 4, 0 of 10 | 10.87 | 11.02 (8) | -4.58 | 2, 4 |
| Friedman | 1 | 9 | -9.41 | 6.18 | 1, 1 | 3, 3 | 9 | 7.96, 3, 3, 0 of 9 | 10.42 | 11.11 (17) | -3.15 | 3, 3 |
| Friedman | 2 (worse) | 11 | -17.71 | 6.69 | 0, 1 | 2, 2 | 11 | 9.17, 5, 6, 1 of 11 | 5.54 | 5.74 (7) | 3.43 | 4, 6 |
| Friedman | 3 (worse) | 9 | -49.96 | -40.50 | 0, 0 | 0, 0 | 9 | 9.27, 4, 4, 0 of 9 | 5.98 | 7.51 (9) | 1.76 | 4, 4 |
| mixed | 0 | 6 | -68.75 | -20.51 | 0, 0 | 0, 0 | 5 | -18.83, 0, 1, 0 of 5 | 1.43 | 1.43 (5) | -20.26 | 0, 0 |
| mixed | 1 | 7 | 2.15 | 3.28 | 1, 1 | 1, 2 | 6 | 3.28, 1, 2, 0 of 6 | 1.92 | 2.60 (13) | 0.68 | 1, 1 |
| mixed | 2 (worse) | 6 | -17.60 | -13.08 | 0, 1 | 0, 1 | 6 | -13.06, 0, 1, 0 of 6 | 4.23 | 4.27 (8) | -17.32 | 0, 1 |
| mixed | 3 | 6 | -56.29 | -43.99 | 0, 0 | 0, 0 | 5 | -40.22, 0, 0, 0 of 5 | 1.88 | 1.88 (3) | -42.10 | 0, 0 |
| functional | 0 | 5 | -329.27 | -237.27 | 0, 0 | 0, 0 | 2 | -222.47, 0, 0, 0 of 2 | 5.88 | 5.88 (3) | -228.35 | 0, 0 |
| functional | 1 (worse) | 7 | -99.39 | 6.72 | 0, 0 | 1, 1 | 2 | 6.72, 1, 1, 0 of 2 | 5.10 | 5.10 (3) | 1.62 | 1, 1 |
| functional | 2 | 5 | -310.26 | -209.14 | 0, 0 | 0, 0 | 2 | -209.14, 0, 0, 0 of 2 | 17.18 | 17.18 (6) | -226.33 | 0, 0 |
| functional | 3 (worse) | 8 | -2.56 | 6.85 | 1, 1 | 3, 4 | 4 | 7.76, 3, 4, 0 of 4 | 11.66 | 11.66 (5) | -3.90 | 3, 3 |

*Each term.* The displacement is the knots' total distinct-value steps at the search's end, run on where it was.

| fixture, chain | term removed | plain Δℓ | profile Δℓ, 5 passes | run on (passes) | less the full profile | displacement |
| --- | --- | --- | --- | --- | --- | --- |
| Friedman 0 (worse) | x1−0.58 | -61.49 | 4.91 | 6.44 (12) | -4.58 | 352 |
| Friedman 0 (worse) | x1−0.31 · x2+0.49 | -49.46 | -18.25 | 3.79 (32) | -7.23 | 532 |
| Friedman 0 (worse) | x1−0.90 · x2−0.91 | -183.83 | -138.57 | -0.69 (34) | -11.72 | 1237 |
| Friedman 0 (worse) | x1−0.61 · x2−0.48 | -56.48 | -23.52 | -8.02 (27) | -19.04 | 531 |
| Friedman 0 (worse) | x1+0.05 · x2−0.38 | -292.53 | -234.07 | -42.02 (39) | -53.04 | 890 |
| Friedman 0 (worse) | x1+0.67 · x2+0.59 | -153.44 | -133.05 | -130.65 (15) | -141.67 | 371 |
| Friedman 0 (worse) | x3+0.40 | -256.38 | -247.65 | -184.55 (35) | -195.57 | 443 |
| Friedman 0 (worse) | x3−0.77 | -238.66 | -222.81 | -185.59 (26) | -196.61 | 349 |
| Friedman 0 (worse) | x5+0.04 | -246.93 | -238.35 | -238.31 (9) | -249.33 | 236 |
| Friedman 0 (worse) | x4−0.99 | -524.23 | -514.95 | -508.62 (20) | -519.64 | 624 |
| Friedman 1 | x1+0.12 · x2−0.61 | -68.02 | 3.46 | 7.96 (11) | -3.15 | 244 |
| Friedman 1 | x1−0.76 · x2−0.90 | -9.41 | 6.18 | 6.25 (7) | -4.86 | 137 |
| Friedman 1 | x2−0.96 | -46.95 | -3.51 | 2.28 (17) | -8.83 | 316 |
| Friedman 1 | x1−0.90 | -264.64 | -188.76 | -119.14 (25) | -130.25 | 464 |
| Friedman 1 | x1+0.53 · x2+0.38 | -237.73 | -207.45 | -155.38 (28) | -166.49 | 502 |
| Friedman 1 | x3+0.36 | -249.82 | -242.80 | -176.48 (39) | -187.59 | 487 |
| Friedman 1 | x3−0.79 | -236.40 | -226.10 | -184.35 (31) | -195.46 | 381 |
| Friedman 1 | x5+0.01 | -236.48 | -228.45 | -227.73 (14) | -238.84 | 298 |
| Friedman 1 | x4−0.95 | -514.03 | -507.98 | -506.36 (15) | -517.47 | 404 |
| Friedman 2 (worse) | x2−0.87 | -224.21 | -120.36 | 9.17 (50 (at the limit)) | 3.43 | 1221 |
| Friedman 2 (worse) | x1+0.56 | -17.71 | 6.69 | 7.57 (10) | 1.83 | 225 |
| Friedman 2 (worse) | x2−0.32 | -67.43 | -29.81 | 6.37 (24) | 0.63 | 545 |
| Friedman 2 (worse) | x1−0.64 · x2−0.45 | -107.25 | -80.56 | 1.55 (39) | -4.19 | 851 |
| Friedman 2 (worse) | x1−0.19 · x2+0.49 | -20.08 | -9.46 | -6.84 (11) | -12.59 | 279 |
| Friedman 2 (worse) | x1−0.86 | -274.44 | -169.78 | -10.84 (43) | -16.58 | 1033 |
| Friedman 2 (worse) | x1+0.56 · x2+0.55 | -218.59 | -195.95 | -185.11 (21) | -190.85 | 612 |
| Friedman 2 (worse) | x3+0.39 | -256.73 | -244.72 | -186.13 (30) | -191.87 | 639 |
| Friedman 2 (worse) | x3−0.71 | -243.71 | -229.09 | -189.39 (26) | -195.13 | 403 |
| Friedman 2 (worse) | x5+0.02 | -247.54 | -239.49 | -238.70 (9) | -244.44 | 241 |
| Friedman 2 (worse) | x4−0.95 | -520.43 | -516.29 | -512.80 (20) | -518.54 | 481 |
| Friedman 3 (worse) | x1−0.84 | -322.03 | -250.21 | 9.27 (33) | 1.76 | 934 |
| Friedman 3 (worse) | x2−0.91 | -274.21 | -201.72 | 9.01 (46) | 1.51 | 1097 |
| Friedman 3 (worse) | x1+0.21 · x2−0.28 | -90.14 | -63.34 | 8.86 (28) | 1.36 | 511 |
| Friedman 3 (worse) | x1−0.30 · x2+0.25 | -49.96 | -40.50 | 6.43 (30) | -1.08 | 451 |
| Friedman 3 (worse) | x3+0.68 | -184.17 | -173.18 | -171.70 (13) | -179.21 | 332 |
| Friedman 3 (worse) | x3−0.32 | -192.73 | -184.04 | -180.33 (17) | -187.84 | 297 |
| Friedman 3 (worse) | x1+0.61 · x2+0.52 | -241.94 | -215.67 | -185.21 (17) | -192.71 | 610 |
| Friedman 3 (worse) | x5+0.02 | -238.18 | -227.90 | -227.74 (10) | -235.25 | 257 |
| Friedman 3 (worse) | x4−0.95 | -509.42 | -501.71 | -498.99 (19) | -506.50 | 496 |
| mixed 0 | x3+0.71 · z{2} | -68.75 | -20.51 | -18.83 (8) | -20.26 | 93 |
| mixed 0 | z{0,1} | -259.98 | -152.76 | -42.84 (20) | -44.27 | 295 |
| mixed 0 | x3−0.30 · z{2} | -92.34 | -85.93 | -85.93 (6) | -87.36 | 66 |
| mixed 0 | x3−0.84 · z{2} | -170.84 | -127.62 | -89.09 (16) | -90.52 | 264 |
| mixed 0 | x1+0.11 | -127.70 | -126.69 | — | -128.12 | 28 |
| mixed 0 | x2−0.96 · z{0,1} | -143.68 | -142.22 | -142.22 (7) | -143.65 | 37 |
| mixed 1 | x1−0.62 · x2−0.38 · x3−0.50 | 2.15 | 3.28 | — | 0.68 | 38 |
| mixed 1 | x3+0.70 · z{2} | -65.27 | -19.92 | -18.32 (9) | -20.91 | 205 |
| mixed 1 | z{2} | -265.51 | -167.29 | -44.28 (21) | -46.88 | 405 |
| mixed 1 | x3−0.28 · z{2} | -94.59 | -87.55 | -82.02 (26) | -84.62 | 387 |
| mixed 1 | x3−0.83 · z{2} | -177.24 | -135.84 | -88.43 (35) | -91.03 | 592 |
| mixed 1 | x1+0.14 | -126.32 | -120.12 | -91.20 (25) | -93.80 | 378 |
| mixed 1 | x2−0.99 · z{0,1} | -140.40 | -131.84 | -113.07 (23) | -115.67 | 384 |
| mixed 2 (worse) | x2−0.84 · x3−0.45 · z{2} | -17.60 | -13.08 | -13.06 (7) | -17.32 | 84 |
| mixed 2 (worse) | x2+0.01 · x3−0.28 · z{2} | -54.55 | -40.28 | -37.87 (15) | -42.14 | 221 |
| mixed 2 (worse) | x3+0.72 · z{2} | -109.34 | -88.01 | -81.13 (14) | -85.40 | 313 |
| mixed 2 (worse) | x3+0.27 · z{2} | -266.04 | -220.08 | -82.85 (29) | -87.11 | 476 |
| mixed 2 (worse) | x1−0.97 | -130.31 | -125.52 | -125.50 (7) | -129.77 | 114 |
| mixed 2 (worse) | x2−0.87 · z{0,1} | -155.27 | -143.93 | -143.24 (11) | -147.50 | 252 |
| mixed 3 | x3−0.27 · z{0,1} | -56.29 | -43.99 | -40.22 (12) | -42.10 | 127 |
| mixed 3 | x3−0.29 | -84.38 | -65.65 | -60.22 (13) | -62.10 | 176 |
| mixed 3 | x3+0.75 · z{2} | -105.11 | -85.74 | -83.18 (13) | -85.06 | 126 |
| mixed 3 | x3+0.21 · z{2} | -276.45 | -238.49 | -98.99 (29) | -100.87 | 484 |
| mixed 3 | x1+0.03 | -128.30 | -124.02 | — | -125.90 | 46 |
| mixed 3 | x2−0.95 · z{0,1} | -161.52 | -147.06 | -146.43 (8) | -148.30 | 116 |
| functional 0 | x0−0.97 | -593.91 | -248.28 | -222.47 (8) | -228.35 | 58 |
| functional 0 | x0+0.68 · t+0.37 | -347.84 | -237.27 | — | -243.15 | 22 |
| functional 0 | x0+0.64 · t−0.74 | -329.27 | -242.55 | — | -248.43 | 14 |
| functional 0 | x0+0.66 | -502.79 | -357.43 | — | -363.31 | 15 |
| functional 0 | x1+0.00 · t+0.00 | -765.51 | -610.65 | -610.65 (6) | -616.53 | 51 |
| functional 1 (worse) | x0−0.95 · t+0.74 | -99.39 | 6.72 | — | 1.62 | 14 |
| functional 1 (worse) | t−0.74 | -571.82 | -40.01 | -40.01 (6) | -45.11 | 50 |
| functional 1 (worse) | x0+0.69 · t+0.37 | -195.08 | -74.25 | — | -79.35 | 32 |
| functional 1 (worse) | x0−0.76 | -713.47 | -146.13 | — | -151.23 | 77 |
| functional 1 (worse) | x0+0.71 · t−0.53 | -262.71 | -202.42 | — | -207.52 | 29 |
| functional 1 (worse) | x0+0.69 | -877.94 | -211.22 | -211.22 (6) | -216.32 | 46 |
| functional 1 (worse) | x1−1.00 · t+0.00 | -564.97 | -543.20 | — | -548.30 | 35 |
| functional 2 | x0+0.61 | -494.77 | -209.14 | — | -226.33 | 53 |
| functional 2 | x0+0.02 | -572.74 | -257.47 | -210.74 (9) | -227.92 | 44 |
| functional 2 | x0+0.59 · t+0.42 | -331.63 | -225.84 | — | -243.02 | 20 |
| functional 2 | x0+0.58 · t−0.68 | -310.26 | -231.15 | — | -248.33 | 13 |
| functional 2 | x1+0.00 · t+0.00 | -741.49 | -599.51 | -594.13 (8) | -611.31 | 57 |
| functional 3 (worse) | x0−1.00 · t+0.05 | -300.94 | 6.35 | 7.76 (7) | -3.90 | 53 |
| functional 3 (worse) | x0+0.00 · t−0.84 | -110.36 | 2.33 | 7.21 (9) | -4.45 | 39 |
| functional 3 (worse) | t+0.79 | -2.56 | 6.85 | — | -4.81 | 20 |
| functional 3 (worse) | x0+0.64 · t−0.37 | -137.08 | -19.76 | — | -31.42 | 45 |
| functional 3 (worse) | x0−0.61 | -867.29 | -74.93 | -56.68 (11) | -68.34 | 61 |
| functional 3 (worse) | x0+0.76 · t+0.63 | -130.63 | -60.28 | -60.23 (7) | -71.89 | 35 |
| functional 3 (worse) | x0+0.66 | -1394.03 | -90.52 | — | -102.18 | 56 |
| functional 3 (worse) | x1−0.92 · t+0.11 | -577.09 | -543.87 | — | -555.53 | 41 |

*What P3 shows.*
- *The 5-pass limit binds.*
  - It stopped 71 of the 89 searches: all 39 on Friedman, 22 of 25 on the mixed fit and 10 of 25 on the functional fit.
  - Run on, all but one converge within 50 passes. The exception, Friedman chain 2's x2−0.87, is still gaining at 50 and already 3.4 above its full profile.
  - The readings below use the converged searches.
- *A sampled state is not its knots' optimum.*
  - The same search on the full state gains 5.7–11.1 nats on Friedman, 1.4–4.3 on the mixed fit and 5.1–17.2 on the functional fit.
  - So a profile death above the fixed state (Δℓ > 0) does not by itself make a term redundant. The full profile is the fair baseline.
- *Friedman: the x1 · x2 block is re-representable, in every chain.*
  - Against the full profile, 4, 3, 6 and 4 terms of chains 0–3 die within 20 nats, all of them x1 or x2 main effects or x1 · x2 interactions. The best are −4.6, −3.2, +3.4 and +1.8.
  - The x3, x4 and x5 terms stay load-bearing, 179–520 nats below it.
  - The best mode, chain 1, reads as the worse ones do.
  - The ends are far: 137–1,237 distinct-value steps in total over a state's hinges, with 500 values per input.
- *Mixed: little is re-representable.* Only chain 1's surplus x1 · x2 · x3 (+0.7 against its full profile) and chain 2's x2−0.84 · x3−0.45 · z{2} (−17.3) come within 20 nats. The closest of the rest, x3+0.7 · z{2}, sits 20.3 and 20.9 below on chains 0 and 1.
- *Functional: the worse modes are re-representable, the better ones not.*
  - Chain 1 holds one term within 10 nats of its full profile: x0−0.95 · t+0.74, at +1.6.
  - Chain 3 holds three: x0−1.00 · t+0.05, x0+0.00 · t−0.84 and t+0.79, at −3.9 to −4.8.
  - Chains 0 and 2 hold none above −226.
  - The displacements there are short, 14–53 steps over x's 60 values and t's 20.
- *Strength: an optimization diagnostic only.*
  - It bounds where a drift could end, with the term removed first and the knots re-optimized after.
  - It does not show a path along which the term is still present while the knots slide, nor how long one would take.
  - The Friedman reading holds on all four chains. The functional one holds on both worse chains against both better ones.

#### P4. Structure moves along the dense ladder

The dense ladders hold 45, 44, 43 and 45 rungs, as in M5's Part 2. Each rung's state is scored under its replica's own temperature:
- The replica's cached log density equals its temperature's density of the state, bit for bit.
- The tempered Δℓ_b differs from Δℓ₁ on every rung below b = 0.99 that holds a term.

*Chain 0, every rung.*
- α is the mean over the kind's proposals, and a rung's deaths are every term's.
- Δℓ_b is the change in the integrated form's tempered log marginal.

| rung | b | M | births α | deaths α: mean (largest) | median | terms over 10⁻³ | changes α | steps α at δ = ±1, ±2, ±3 | prunes α (n) | deaths' median Δℓ_b | median b·Δℓ₁ | median Δℓ₁ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1 | 10 | 0.00669 | 5.29E-08 (5.29E-07) | 2.89E-75 | 0 | 0.0469 | 0.930, 0.879, 0.833 | 2.24E-27 (6) | -207.35 | -207.35 | -207.35 |
| 1 | 0.8943 | 9 | 0.0202 | 2.94E-20 (2.65E-19) | 5.88E-77 | 0 | 0.0425 | 0.978, 0.948, 0.897 | 5.42E-74 (4) | -177.59 | -178.35 | -199.42 |
| 2 | 0.7949 | 11 | 0.0118 | 8.54E-05 (0.000939) | 8.65E-38 | 0 | 0.0813 | 0.982, 0.962, 0.940 | 1E-07 (8) | -87.42 | -88.61 | -111.48 |
| 3 | 0.7245 | 9 | 0.0291 | 9.29E-11 (8.36E-10) | 2.36E-73 | 0 | 0.0658 | 0.935, 0.905, 0.876 | 3.04E-56 (4) | -169.29 | -170.74 | -235.65 |
| 4 | 0.6456 | 9 | 0.0115 | 2.91E-08 (2.62E-07) | 5.22E-65 | 0 | 0.0611 | 0.967, 0.930, 0.911 | 3.14E-07 (6) | -149.99 | -151.99 | -235.42 |
| 5 | 0.5717 | 10 | 0.0112 | 1.18E-10 (1.18E-09) | 1.34E-57 | 0 | 0.0803 | 0.988, 0.978, 0.945 | 4.16E-06 (6) | -133.16 | -136.18 | -238.20 |
| 6 | 0.4971 | 8 | 0.0269 | 9.46E-36 (7.57E-35) | 8.87E-52 | 0 | 0.0632 | 0.953, 0.919, 0.896 | 3.41E-36 (4) | -128.16 | -130.72 | -262.98 |
| 7 | 0.4376 | 8 | 0.0127 | 9.82E-36 (4.45E-35) | 3.29E-43 | 0 | 0.0504 | 0.975, 0.938, 0.904 | 8.1E-36 (4) | -112.84 | -116.58 | -266.39 |
| 8 | 0.4 | 8 | 0.00864 | 9.38E-09 (7.51E-08) | 9.81E-40 | 0 | 0.0815 | 0.978, 0.954, 0.927 | 3.21E-37 (4) | -93.43 | -96.46 | -241.14 |
| 9 | 0.39 | 7 | 0.0113 | 1.84E-08 (1.29E-07) | 1.23E-37 | 0 | 0.0645 | 0.948, 0.917, 0.878 | 2.96E-10 (6) | -86.69 | -89.60 | -229.75 |
| 10 | 0.38 | 9 | 0.01 | 0.111 (1) | 3.98E-39 | 1 | 0.115 | 0.976, 0.950, 0.910 | 9.81E-33 (4) | -90.49 | -93.80 | -246.83 |
| 11 | 0.37 | 10 | 0.0124 | 1.02E-06 (1.02E-05) | 8.7E-30 | 0 | 0.0884 | 0.986, 0.970, 0.956 | 2.96E-19 (4) | -69.12 | -73.61 | -198.95 |
| 12 | 0.36 | 9 | 0.0296 | 1.47E-08 (1.32E-07) | 2.2E-30 | 0 | 0.0984 | 0.967, 0.948, 0.928 | 3.2E-17 (4) | -70.36 | -74.12 | -205.89 |
| 13 | 0.35 | 7 | 0.0482 | 1.99E-32 (1.3E-31) | 4.88E-48 | 0 | 0.0614 | 0.952, 0.938, 0.913 | 5.11E-26 (4) | -110.63 | -114.83 | -328.08 |
| 14 | 0.34 | 8 | 0.0144 | 3.92E-25 (3.05E-24) | 8.38E-31 | 0 | 0.0942 | 0.979, 0.961, 0.933 | 2.66E-23 (4) | -71.34 | -74.97 | -220.49 |
| 15 | 0.33 | 8 | 0.00922 | 5.54E-25 (4.43E-24) | 3.55E-34 | 0 | 0.0553 | 0.958, 0.925, 0.883 | 7.69E-26 (4) | -84.80 | -88.12 | -267.03 |
| 16 | 0.32 | 8 | 0.0184 | 1.81E-20 (1.45E-19) | 3.98E-31 | 0 | 0.0744 | 0.957, 0.926, 0.903 | 3.11E-21 (4) | -76.50 | -80.29 | -250.90 |
| 17 | 0.31 | 8 | 0.0325 | 9.39E-17 (7.5E-16) | 1.05E-21 | 0 | 0.125 | 0.982, 0.939, 0.920 | 9.53E-10 (4) | -51.69 | -56.24 | -181.44 |
| 18 | 0.3 | 9 | 0.0202 | 0.111 (1) | 5.35E-25 | 1 | 0.189 | 0.982, 0.965, 0.949 | 0.429 (7) | -57.86 | -65.20 | -217.34 |
| 19 | 0.29 | 8 | 0.0481 | 0.0362 (0.29) | 2.14E-18 | 1 | 0.132 | 0.975, 0.934, 0.906 | 3.71E-19 (4) | -48.52 | -52.51 | -181.08 |
| 20 | 0.28 | 9 | 0.0545 | 0.111 (1) | 5.61E-18 | 1 | 0.191 | 0.980, 0.967, 0.953 | 0.333 (6) | -41.69 | -46.60 | -166.44 |
| 21 | 0.27 | 6 | 0.0724 | 2.74E-17 (1.57E-16) | 1.37E-21 | 0 | 0.0811 | 0.922, 0.908, 0.868 | 7.81E-38 (4) | -65.64 | -69.14 | -256.09 |
| 22 | 0.26 | 7 | 0.0237 | 4.37E-12 (1.33E-11) | 2.94E-15 | 0 | 0.099 | 0.992, 0.980, 0.975 | 3.76E-12 (2) | -35.39 | -38.11 | -146.58 |
| 23 | 0.25 | 7 | 0.0195 | 5.08E-08 (3.53E-07) | 1.68E-13 | 0 | 0.124 | 0.979, 0.970, 0.960 | 1.58E-08 (2) | -31.34 | -34.06 | -136.25 |
| 24 | 0.24 | 8 | 0.0233 | 0.125 (1) | 1.52E-13 | 1 | 0.167 | 0.996, 0.989, 0.971 | 1.43E-10 (2) | -31.86 | -35.25 | -146.90 |
| 25 | 0.23 | 6 | 0.0257 | 5.79E-09 (2.85E-08) | 1.01E-10 | 0 | 0.184 | 0.992, 0.984, 0.979 | — (0) | -34.48 | -37.65 | -163.70 |
| 26 | 0.22 | 6 | 0.0184 | 7.53E-10 (3.65E-09) | 2.16E-10 | 0 | 0.117 | 0.994, 0.991, 0.982 | — (0) | -34.31 | -37.15 | -168.88 |
| 27 | 0.21 | 6 | 0.0147 | 7.22E-07 (3.76E-06) | 9.52E-10 | 0 | 0.153 | 0.995, 0.987, 0.979 | — (0) | -31.74 | -34.34 | -163.53 |
| 28 | 0.2 | 7 | 0.0177 | 0.286 (1) | 2.32E-05 | 2 | 0.308 | 0.993, 0.988, 0.977 | — (0) | -12.69 | -15.05 | -75.27 |
| 29 | 0.19 | 6 | 0.0167 | 2.34E-07 (1.32E-06) | 6.32E-09 | 0 | 0.178 | 0.995, 0.992, 0.966 | — (0) | -28.71 | -31.41 | -165.33 |
| 30 | 0.18 | 6 | 0.0363 | 3.49E-06 (1.44E-05) | 6.03E-07 | 0 | 0.23 | 0.998, 0.993, 0.987 | — (0) | -23.78 | -27.43 | -152.38 |
| 31 | 0.17 | 7 | 0.0412 | 0.143 (1) | 5.53E-06 | 1 | 0.332 | 0.986, 0.980, 0.966 | 1 (2) | -14.03 | -16.98 | -99.88 |
| 32 | 0.16 | 7 | 0.0305 | 0.156 (1) | 6.36E-06 | 2 | 0.332 | 0.995, 0.991, 0.987 | 0.976 (2) | -13.89 | -17.22 | -107.62 |
| 33 | 0.15 | 6 | 0.0269 | 0.000481 (0.00272) | 7.45E-06 | 1 | 0.245 | 0.996, 0.985, 0.963 | — (0) | -19.38 | -22.16 | -147.75 |
| 34 | 0.14 | 6 | 0.0227 | 8.96E-06 (2.43E-05) | 5.01E-06 | 0 | 0.184 | 0.993, 0.988, 0.985 | — (0) | -19.92 | -23.11 | -165.05 |
| 35 | 0.13 | 8 | 0.027 | 0.251 (1) | 0.000861 | 4 | 0.381 | 0.997, 0.994, 0.988 | 0.955 (2) | -9.19 | -12.25 | -94.23 |
| 36 | 0.12 | 6 | 0.0764 | 0.334 (1) | 0.00278 | 3 | 0.369 | 0.996, 0.988, 0.987 | 0.992 (3) | -10.17 | -12.58 | -104.85 |
| 37 | 0.11 | 4 | 0.0387 | 0.00144 (0.00575) | 1.05E-05 | 1 | 0.224 | 0.995, 0.985, 0.978 | — (0) | -14.15 | -16.62 | -151.09 |
| 38 | 0.1 | 4 | 0.027 | 0.00677 (0.0261) | 0.000489 | 1 | 0.291 | 0.995, 0.981, 0.977 | — (0) | -10.48 | -12.90 | -129.01 |
| 39 | 0.09877 | 4 | 0.0595 | 0.0113 (0.0452) | 9.81E-05 | 1 | 0.296 | 0.992, 0.982, 0.968 | — (0) | -11.36 | -14.37 | -145.51 |
| 40 | 0.08322 | 3 | 0.148 | 0.0475 (0.135) | 0.00736 | 2 | 0.421 | 0.994, 0.992, 0.987 | — (0) | -6.14 | -8.60 | -103.38 |
| 41 | 0.07214 | 5 | 0.0622 | 0.428 (1) | 0.129 | 4 | 0.537 | 0.996, 0.991, 0.987 | 0.507 (2) | -3.37 | -5.34 | -74.02 |
| 42 | 0.06131 | 2 | 0.161 | 0.0699 (0.137) | 0.0699 | 2 | 0.425 | 0.998, 0.998, 0.996 | — (0) | -4.68 | -6.68 | -108.98 |
| 43 | 0.04533 | 1 | 0.197 | 0.14 (0.14) | 0.14 | 1 | 0.179 | 0.995, 0.993, 0.991 | 0.532 (2) | -3.16 | -4.76 | -104.99 |
| 44 | 0.002 | 0 | 0.174 | — (—) | — | 0 | — | —, —, — | — (0) | — | — | — |

*All four chains at the window's shared rungs:*

| b | chain 0: M; median death α; median Δℓ_b; changes α | chain 1: M; median death α; median Δℓ_b; changes α | chain 2: M; median death α; median Δℓ_b; changes α | chain 3: M; median death α; median Δℓ_b; changes α |
| --- | --- | --- | --- | --- |
| 1 | 10; 2.89E-75; -207.35; 0.0469 | 13; 3.37E-41; -95.35; 0.0517 | 11; 2.31E-87; -201.56; 0.0429 | 9; 8E-111; -255.57; 0.0454 |
| 0.4 | 8; 9.81E-40; -93.43; 0.0815 | 10; 2.11E-28; -66.83; 0.115 | 8; 6.23E-31; -79.85; 0.0593 | 8; 8.3E-41; -101.85; 0.0553 |
| 0.35 | 7; 4.88E-48; -110.63; 0.0614 | 10; 6.28E-28; -65.31; 0.138 | 9; 9.99E-32; -73.45; 0.108 | 9; 4.35E-27; -62.77; 0.1 |
| 0.3 | 9; 5.35E-25; -57.86; 0.189 | 8; 1.24E-22; -53.99; 0.123 | 8; 5.41E-24; -60.60; 0.0937 | 7; 6.98E-41; -94.15; 0.0875 |
| 0.25 | 7; 1.68E-13; -31.34; 0.124 | 7; 2.25E-12; -28.75; 0.116 | 8; 1.77E-15; -36.45; 0.12 | 9; 3.82E-15; -35.26; 0.239 |
| 0.2 | 7; 2.32E-05; -12.69; 0.308 | 6; 5.01E-09; -29.70; 0.131 | 6; 5.06E-09; -29.64; 0.133 | 6; 5.11E-09; -29.71; 0.147 |
| 0.15 | 6; 7.45E-06; -19.38; 0.245 | 5; 9.61E-08; -17.76; 0.358 | 4; 4.85E-08; -20.32; 0.194 | 6; 1.26E-06; -21.65; 0.157 |
| 0.12 | 6; 0.00278; -10.17; 0.369 | 4; 5.76E-06; -14.46; 0.244 | 5; 8.19E-08; -17.92; 0.313 | 6; 0.000292; -15.23; 0.256 |
| 0.1 | 4; 0.000489; -10.48; 0.291 | 4; 0.000282; -10.43; 0.34 | 4; 2.25E-05; -12.99; 0.249 | 4; 0.000223; -10.27; 0.356 |

*Per chain:* the coldest rung at which deaths stop being negligible, by three readings, and the band's rungs pooled.

| chain | rungs | coldest b with deaths' mean α > 0.001 | with a death's α > 0.001 | with the median term's α > 0.001 | band [0.1, 0.4], pooled: M range | births α | deaths α | term deaths over 0.001 | changes α | steps α at δ = ±1, ±2, ±3 | prunes α (n) | rungs whose median term's death α ≤ 0.001 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 45 | 0.38 | 0.38 | 0.12 | 4–10 | 0.0283 | 0.0562 | 20 of 221 | 0.173 | 0.980, 0.963, 0.945 | 0.177 (78) | 30 of 31 |
| 1 | 44 | 0.4 | 0.4 | 0.14 | 4–10 | 0.0272 | 0.0243 | 14 of 217 | 0.163 | 0.981, 0.965, 0.945 | 0.154 (87) | 30 of 31 |
| 2 | 43 | 0.778 | 0.778 | 0.0947 | 4–10 | 0.0247 | 0.0232 | 10 of 225 | 0.155 | 0.981, 0.965, 0.950 | 0.0857 (77) | 31 of 31 |
| 3 | 45 | 0.925 | 0.925 | 0.0783 | 4–9 | 0.0266 | 0.0153 | 9 of 218 | 0.146 | 0.983, 0.969, 0.952 | 0.0536 (76) | 31 of 31 |

The hottest rungs, below b = 0.09, hold M = 0–5. There the median term's death α is 0.0008–1, and births accept 0.06–0.34.

*What P4 shows.*
- *The band's rungs are frozen for deaths.*
  - In [0.10, 0.40], the median term's death α stays at or below 10⁻³ on 30, 30, 31 and 31 of the 31 rungs.
  - The typical term's death becomes non-negligible only at b = 0.08–0.14, the band's hot edge and below.
- *A rung's mean death α misleads.*
  - One surplus term at a time, often just born, dies with α = 1 and sets the mean.
  - By the mean, deaths stop being negligible at b = 0.38–0.93; by the median term, at 0.08–0.14.
  - Pooled over the band, 4–9% of term deaths accept above 10⁻³.
- *The draft's b · Δℓ arithmetic holds.*
  - The median death's Δℓ_b runs −63 to −111 at b = 0.40 and 0.35, −29 to −36 at 0.25, and −10 to −18 at 0.10–0.12.
  - It sits 0.5–7.3 nats above b · Δℓ₁ on every rung above b = 0.05, so the integrated form is slightly gentler than the linear one.
  - On those rungs the median death's untempered Δℓ₁ stays between −74 and −345: the hot rungs' terms are as load-bearing, untempered, as the cold ones'.
- *The other moves loosen as b falls, and the local steps are free everywhere.*
  - *Births* accept 0.025–0.028, pooled over the band, against 0.007 at chain 0's b = 1.
  - *Change redraws* accept 0.15–0.17, against 0.04–0.05 at b = 1.
  - *Local steps* accept 0.92–1.00 at ±1 and 0.83–1.00 at ±3, on every rung.
  - *Prunes* accept 0.05–0.18, pooled over the band, but only where a prune lands on a surplus term; elsewhere they are negligible.
- *Strength:* strong for the deaths, on four chains and 177 rung states, though with one state per rung, at the end of burn-in.

#### Run times

Release, four chains in parallel within each harness. Each harness ran alone, on a machine shared with the main session.
- *P1:* the states and their phases 0.3–0.7 s per fixture, with each fixture's S10 fit for the twins' check; the removals under 0.05 s.
- *P2:* the states 0.2–1.0 s; the moves and two-step profiles 0.2–1.4 s per fixture.
- *P3:* the states 0.3–0.9 s; the searches 0.2–1.2 s.
- *P4:* 25 s. Each chain, in parallel, burns its coarse ladder and its dense ladder, then scores its rungs.
- *The Regression suite* without the switch: 1 min 32 s.

#### Accessors and suites

- No internal accessor was added, no production file changed, and `BassSweep.cs` was not touched.
- M2–M4 reran with identical readings.
- The Regression suite, rerun without the switch: 838/838 (828 before, plus the 10 new cases).

#### The review passes and their fixes

1. *P1's first read* labelled three rows per state "co-adapted" by their relief, which is negative for every pair. Those rows are now the most redundant pairs by Δℓ_pair − Δℓ₁ − Δℓ₂, and the summary counts the pairs cheaper than their cheaper single (none).
2. *P2's first two-step profile* counted the re-represented term's own death, which only reaches the fixed state's single death of the old term. It now reads the other terms' removals, beside:
   - the same removals at the fixed state;
   - the fixed state's own best;
   - the end state's log posterior.
3. *P2's cheapest exits* listed ten knots of one family. They are now by family, with the family's proposals and mean α.
4. *P2's unavailable variable changes* carried a default score that read as admitted. The two-step profile now takes only available proposals; this was found before the first run.
5. *P3's 5-pass limit* stopped 71 of its 89 searches. They now run on to 50 passes, and all but one converge.
6. *P3's profile deaths above the fixed state* looked like redundancy. A sampled state is not its knots' optimum, so the full state's own profile is now the control.
7. *P4's first read* had three faults:
   - Its hottest rungs can hold M = 0, with no change or death to propose. Those are now counted as unavailable.
   - One surplus term set a rung's mean death α. The median term's α and the count above 10⁻³ now stand beside it.
   - The band's per-rung means went blank wherever a rung held no prune. The band's readings are now pooled over its proposals.
8. *P4's header* printed |δ|, which split its markdown tables (`811360e`). The readings are the same.
9. *Correctness held.*
   - *No commit:* checked after every scorer: P1 after its removals, P2 after its moves and after its two-step profile, P3 after its searches, and P4 per rung.
   - *P1's q_pair* was recomputed from both orders' counts, to 10⁻¹⁰.
   - *The single deaths* matched the chain's design (RSS and Q to 10⁻⁹) and the death move's ratio (to 10⁻¹²).
   - *P2's cancellation* holds, as reported above.
   - *P4* scores each rung at its own b.
   - *The counts add up:*
     - the removals number C(M, k);
     - each prune is enumerated once;
     - each state draws 5,000 grows and 5,000 changes;
     - P4's changes are the scored ones, the null redraws and the unavailable ones.
   - *The twins* reproduce the run's draws.
10. *Coverage.*
    - *A state outside its mode:* Friedman chain 3, replaced.
    - *Readings resting on one chain:* P2's two-step exits rest on mixed chain 2 and functional chain 3, and P3's functional reading on chains 1 and 3.
    - *Readings that hold on every chain:* P1, P3 on Friedman, and P4.
    - *Triples were never read,* since no best pair came under 20 nats.

#### For the main session

- *Phase 0's Friedman chain 3.* Its state at the end of burn-in lies below its whole sampling phase, holding a just-born inert x6 · x8. M2's exception there (that term's death at α 0.20) and its births came from that transient edge state.
- *The helper's signatures.* `Score` takes an optional temperature, and `Untouched` has a replica constructor. Phase 0's tests read the same.
- *Design points the harnesses had to fix, for the user's pass* (not decided here):
  - the move and part selection for factor moves under several parts;
  - whether a part's only factor may be pruned while the term keeps another part, the reverse of a grow into an empty part, and the only prune on the functional fit;
  - a variable change in a part with one variable, which has no unused variable to move to.
- *Readings that bear on the design draft's "After Phase 0: proposals"* (not mine to edit):
  - No pair death in any state costs less than 57 nats (P1).
  - The fixed-M moves offer no cheap exit from Friedman's worse modes 0 and 3. In the mixed and functional worse states, a prune or variable change followed by a death is a near-free path to the better modes' level (P2).
  - Drift endpoints within 20 nats exist for Friedman's x1 · x2 block in every chain, and for the functional fit's worse chains (P3).
  - The band's rungs are frozen for deaths, as the draft's b · Δℓ arithmetic predicted (P4).
