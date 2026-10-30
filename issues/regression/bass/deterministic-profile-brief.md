# F-R1: the deterministic-simulator profile's comparison (design draft and chip brief)

Date: 2026-10-30. Status: a draft for the user's pass. "The decisions" below lays out the choices with recommendations. Once they are taken, "The comparison, fixed before any run" is the chip's brief, and nothing in it changes after a run.

## Why

The fitting segment's regime note (2026-10-26, [port-fitting-brief.md](port-fitting-brief.md), "Regime note: deterministic simulators") proposed a named profile, `DeterministicSimulator`, never a default:
- D03's scaled, proper σ² prior (already the default);
- a wider λ prior, with h₁ and h₂ below the usual 10 (GBMARS's advice for interpolating surrogates, which need more terms);
- the ridge as the alternative coefficient prior (KHAOS found its fixed ridge best on noise-free simulators and its g-prior on most noisy ones);
- in direct functional fits, the residual-dependence check read before any interval, since a simulator's residual is all misfit.

Which coefficient prior the profile takes was left to a pre-registered comparison of our own. It was waiting on the ridge, which is built (2026-10-29, `docs/Maths/BayesianRegression.md`, "BASS: the coefficient prior's two axes, and the ridge"). The Opus chat's column-scale concern, that a common τ on raw hinge columns shrinks terms by where their knots sit, is handled by building the ridge on unit-RMS columns. So the comparison measures the prior, not a column scaling. The Bissiri reading (the extensions brief) adds that σ² on deterministic data is a learned loss weight, so the profile keeps σ² sampled.

## The decisions (for the user's pass, with recommendations)

1. **The coefficient priors compared.** *Recommended:* the 2 × 2 of form and scale, so the ridge-against-g question and the sampled-against-held question are read apart:
   - the g-form under Zellner–Siow, today's default;
   - the g-form with τ held at τ₀ = 1/N, D03's second mode;
   - the ridge with τ sampled under the unit-information Gamma at a = ½ on τ₀ = 1;
   - the ridge with τ held at 1. The held scale is part of the model, so it is declared here; τ = 1 is one observation's worth per term, and on orthogonal columns it is the g-form at τ₀.
   - *Not compared:* hyper-g and hyper-g/n, which would double the arms for a question the 2 × 2 already frames.
2. **The λ prior.** *Recommended:* the default (h₁, h₂) = (10, 10) and one wider value, (1, 1).
   - Both have prior mean λ = 1; (1, 1) has the heavier tail. In the count law's ratio P(M)/P(M − 1) = (h₁ + M − 1)/(M·(h₂ + 1)):
     - under (10, 10), a term costs about 2 nats near M = 20, against about 0.1 for the first;
     - under (1, 1), λ ~ Exp(1), the law is geometric with ratio ½, about 0.69 nats per term at every M.
   - The alternative, (1, 0.1), costs about 0.1 nats per term, close to no count penalty at all. It is not recommended: it hands model size to the term measure and the Occam factor alone.
3. **The benchmarks** (noise-free, each with a truth computed in the harness):
   - Friedman's f₁ with no noise, p = 10, five of them inert (S10's generator);
   - Ishigami, a = 7 and b = 0.1, on [−π, π]³;
   - the borehole function, p = 8, on its standard ranges (Morris, Mitchell and Ylvisaker 1993);
   - S10's direct functional fixture's mean, f = 1 + 2x₀ + 1.5·t·x₁ + 3·(x₀ − 0.5)₊·sin(πt), with no noise (n = 60, m = 20).
   - *Recommended:* these four. Two are smooth (Ishigami, borehole), one has inert inputs (Friedman), and one is functional.
4. **The designs and the holdout.** *Recommended:*
   - N = 500 for the scalar benchmarks, from three independent random Latin hypercubes each;
   - for the functional benchmark, n = 60 cases from three independent draws, as its generator makes them;
   - a holdout of 10,000 uniform points for each scalar benchmark (200 cases × m for the functional one), fixed across arms and designs;
   - one training size only, so the reading does not cover how the answer moves with N.
5. **The sampler and the mixing gate.**
   - *Recommended:* S10's MCMC configuration: four chains, `BurnIn` 5,000, `Samples` 2,000, untempered. That is what a user would run, and it matches F-M2's pass.
   - *The trap.* Noise-free posteriors are sharp, σ² runs small, and BASS's structure mixing is already slow (F-M2). An arm that wins only because its chains stall in one good structure is no win.
   - *The gate, recommended:* a (benchmark, design) pair where any arm's fit has a predictive-group rank-normalized R̂ above 1.05 is reported and left out of the score for every arm alike, so the arms stay paired. If two or more of a benchmark's three designs are left out, the benchmark is unscored, and the rule runs on the rest.
   - *Why 1.05 and not D10's verdict:* D10's reference target is 1.01, and S10's untempered Friedman chains read structural R̂ of 1.5–2.5 (F-M2), so its verdict would likely leave Friedman out for every arm. The predictive group agrees better than the structure. D10's verdict is reported for every fit beside the gate.
   - *The alternative:* the self-sized tempered ladder for every arm, which mixes better but is not what a default fit runs and costs several times more.
6. **The metrics.** *Recommended:*
   - *primary:* the holdout error of the posterior mean surface, as NRMSE, the root mean square error over the truth's sd on the holdout;
   - *the second gate:* the coverage of the truth by the pointwise 95% latent band, over the holdout points;
   - *context:* the posterior's M, the time, and D10's verdict. For the functional benchmark, also the residual-dependence check on the training residuals, read before its intervals.
7. **The rule.** *Recommended*, in F-M2's form:
   - per benchmark f, an arm's ratio is r_f = exp(mean over designs of ln NRMSE(arm) − the same for A0), with A0 today's default (the g-form under Zellner–Siow, λ at (10, 10));
   - the score is the geometric mean of r_f over the scored benchmarks; lower is better;
   - an arm clears when its score is at most 0.9, no r_f exceeds 1.05, and its coverage, averaged over designs, is no more than 0.05 below A0's on any scored benchmark;
   - among arms that clear, the lowest score wins, and its coefficient prior and λ become the profile's;
   - if none clears, no profile is built, every option stays an option, and F-R1 closes as read.
   - The thresholds ask for a 10% gain with no benchmark more than 5% worse. The coverage condition keeps an arm from winning by a fit more overconfident than the default's. It is relative because a noise-free band's absolute coverage depends on how σ² settles, for every arm alike.
8. **The noisy controls** (read, not ruled). *Recommended:* the same arms on two noisy twins, on the same designs: Friedman's f₁ with N(0, 1) noise, S10's level, and Ishigami with Gaussian noise of sd 0.37, about a tenth of its output's sd. The report says whether the winner loses to A0 there, as KHAOS's regime split predicts. That is what keeps the profile a profile and not a default.
9. **What a verdict does.** *Recommended:*
   - a winner becomes the named profile, a `BassConfig` factory in production built after the verdict, with its λ, its coefficient prior and D03's σ² default. Its docs say that direct functional fits read the residual check before intervals.
   - No default changes either way.

## The comparison, fixed before any run (as recommended)

**Arms.** Each is S10's configuration (four chains, `BurnIn` 5,000, `Samples` 2,000, `Thin` 1, untempered) with:

| arm | `Coefficients` | `Count` |
| --- | --- | --- |
| A0, the default | `BassCoefficientPrior.G()` | `BassCountPrior(10, 10)` |
| A1 | `BassCoefficientPrior.G(TauPrior.Held())` | `BassCountPrior(10, 10)` |
| A2 | `BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma())` | `BassCountPrior(10, 10)` |
| A3 | `BassCoefficientPrior.Ridge(TauPrior.Held(1.0))` | `BassCountPrior(10, 10)` |
| A4 | `BassCoefficientPrior.G()` | `BassCountPrior(1, 1)` |
| A5 | `BassCoefficientPrior.G(TauPrior.Held())` | `BassCountPrior(1, 1)` |
| A6 | `BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma())` | `BassCountPrior(1, 1)` |
| A7 | `BassCoefficientPrior.Ridge(TauPrior.Held(1.0))` | `BassCountPrior(1, 1)` |

**Benchmarks,** as decision 3, noise-free. The scalar ones have N = 500 from three random Latin hypercubes each; the functional one has three draws of its 60 cases. The holdouts are as decision 4.

**Seeds.** Each (benchmark, design) has one data seed and one master seed, shared by all eight arms, so the arms are paired. Starts are the fit's own exact prior draws (S7). The seeds are declared in the harness before the first run.

**Running.**
- Release, one fit at a time, `MaxDegreeOfParallelism` = 1, the harness run alone (`--filter` on its class).
- Within each (benchmark, design), the arms run in an order shuffled by a fixed seed.
- The time is the whole `Bass.Fit`.
- *A budget:* if the 1/10 run (below) projects the full comparison past about 6 hours, stop there and report the projection, before any timed run.

**Readings, per fit.**
- NRMSE of the posterior mean surface on the holdout;
- the 95% latent band's coverage of the truth on the holdout;
- the predictive group's rank-normalized R̂ and ESS (the gate), and D10's verdict;
- the posterior median and 95% interval of M, and whether any draw reached the basis cap;
- the time;
- for the functional benchmark, the residual-dependence check (lag one along t) on the training residuals.

**The rule,** decision 7 verbatim, with decision 5's gate applied first and reported. No arm, benchmark or design is dropped after the runs except by that gate.

**Context,** reported beside the verdict and not part of it:
- each arm's per-design log ratios and their spread;
- each arm's score with any one benchmark left out;
- the noisy controls' r_f, by the same arithmetic.

## Checks before the runs

- *The benchmarks are what they say:*
  - Ishigami's mean a/2 = 3.5 and variance a²/8 + bπ⁴/5 + b²π⁸/18 + ½ ≈ 13.84, by Monte Carlo over 10⁶ points within 5 SE;
  - Friedman's f₁ against S10's generator at its own inputs, to rounding;
  - the borehole function's monotone directions (increasing in r_w and H_u, decreasing in H_l and L) at random points;
  - the functional benchmark against `BassQualificationRecoveryTests.Functional`'s truth, to rounding.
- *The arms are what they say:* each fit's resolved coefficient prior, τ₀ and count prior, and the fingerprints differ across arms exactly as the table does.
- *The pairing holds:* the eight arms of one (benchmark, design) see the same data and start from the same exact prior draws.
- *The metrics are right:* NRMSE and coverage on a fit whose posterior mean is set by hand.
- Run the whole comparison once at 1/10 of the lengths on one design, and check that every reading is produced and the rule's arithmetic runs, before the real run.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Tests only.** The harness is `tests/regression/BassDeterministicProfileTests.cs`, behind `BASS_QUALIFICATION` = 1. It changes no production file; the profile itself is built after the verdict.
- **The flat-histogram Phase 2 chip works beside you** on BASS's production files (`BassConfig.cs`, `Bass.cs`, `BassFit.cs` and others; [its brief](flat-histogram-phase2-brief.md)). Build against the checkout as it is. If its commits change an API your harness uses, adapt your harness; never edit its files. If a build fails mid-change in a file you do not own, wait and retry.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-fr1\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-fr1\` when you finish.
- **Your files:** the harness, any test-side helpers it needs under `tests/regression/`, and this brief's Reports section.
- **Commit your own files only, by explicit path,** in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. The user's untracked notes and other sessions' edits stay unstaged. Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - the benchmarks' formulas and domains;
  - the holdout fixed across arms;
  - NRMSE and coverage computed on the latent surface, not the predictive;
  - the gate applied before the rule, and reported;
  - the rule's arithmetic, against a hand computation on the per-fit table.
- *Coverage:*
  - a reading resting on one design;
  - an arm whose fits reach the basis cap;
  - a benchmark the gate leaves unscored.

## Report

Under "Reports" below, newest first:
- the verdict first: each arm's score, its r_f per benchmark, its coverage, and which arm clears, if any;
- the gate's exclusions;
- the per-fit table;
- the noisy controls;
- the times;
- the checks with their tolerances;
- the review passes;
- anything the main session must know.

Recommend no default; the rule decides the profile.

## Reports

Newest first.
