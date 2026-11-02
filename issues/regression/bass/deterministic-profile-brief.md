# F-R1: the deterministic-simulator profile's comparison (design draft and chip brief)

Date: 2026-10-30. Status: decided. The user took the nine recommendations below as drafted on 2026-10-30, so "The comparison, fixed before any run" is the chip's brief, and nothing in it changes after a run.

## Why

The fitting segment's regime note (2026-10-26, [port-fitting-brief.md](port-fitting-brief.md), "Regime note: deterministic simulators") proposed a named profile, `DeterministicSimulator`, never a default:
- D03's scaled, proper σ² prior (already the default);
- a wider λ prior, with h₁ and h₂ below the usual 10 (GBMARS's advice for interpolating surrogates, which need more terms);
- the ridge as the alternative coefficient prior (KHAOS found its fixed ridge best on noise-free simulators and its g-prior on most noisy ones);
- in direct functional fits, the residual-dependence check read before any interval, since a simulator's residual is all misfit.

Which coefficient prior the profile takes was left to a pre-registered comparison of our own. It was waiting on the ridge, which is built (2026-10-29, `docs/Maths/BayesianRegression.md`, "BASS: the coefficient prior's two axes, and the ridge"). The Opus chat's column-scale concern, that a common τ on raw hinge columns shrinks terms by where their knots sit, is handled by building the ridge on unit-RMS columns. So the comparison measures the prior, not a column scaling. The Bissiri reading (the extensions brief) adds that σ² on deterministic data is a learned loss weight, so the profile keeps σ² sampled.

## The decisions (decided 2026-10-30, as recommended)

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

## The comparison, fixed before any run

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
- **Other chips may work beside you,** the flat-histogram warm-start chip among them ([its brief](flat-histogram-warmstart-brief.md)), in `src/maths/samplers/flathistogram/`, `BassFlatHistogramRun.cs` and the flat-histogram tests. Build against the checkout as it is. If a commit changes an API your harness uses, adapt your harness; never edit another chip's files. If a build fails mid-change in a file you do not own, wait and retry. Its long runs share the machine, so your times are context, as the brief says.
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

## For a rerun

Whether F-R1 reruns is the user's call, and nothing here changes the comparison that ran. *The user's ruling, 2026-11-01:* when it does, a truncated σ² prior, R's `s2.lower`, joins as a pre-registered arm.
- *Why:* the contract shapes left `s2.lower` out on grounds of form, before this comparison ran. This brief names σ² running small as the trap on noise-free posteriors, and its arms varied the coefficient prior and the count alone. So the one control R offers for that was never read.
- *It changes the target.* With σ² ≥ s under a prior IG(a₀, b₀), the collapsed marginal gains the factor P(a, b/s)/P(a₀, b₀/s), P the regularized lower incomplete gamma and (a, b) the σ² posterior's shape and scale, and the composition draws σ² from the truncated inverse gamma. So it needs its own derivation, and checks against enumeration, before a comparison reads it.
- *To declare before the run:* the floor s or its rule, since it sets how much of the response the fit may leave as noise, and that is the arm's meaning.

## Reports

Newest first.

### 2026-10-31: the run (harness `8349663`, review fixes `fb35f51`)

**The verdict: no arm clears.** The gate leaves out all twelve noise-free (benchmark, design) pairs, so no benchmark is scored, and no arm has a score, r_f or coverage under the rule. By decision 7 no profile is built, every option stays an option, and F-R1 closes as read.

| arm | score | r_f per benchmark | coverage | clears |
| --- | --- | --- | --- | --- |
| A1, g, τ held at 1/N, λ (10, 10) | none | none scored | none scored | no |
| A2, ridge, τ sampled, λ (10, 10) | none | none scored | none scored | no |
| A3, ridge, τ held at 1, λ (10, 10) | none | none scored | none scored | no |
| A4, g, Zellner–Siow, λ (1, 1) | none | none scored | none scored | no |
| A5, g, τ held at 1/N, λ (1, 1) | none | none scored | none scored | no |
| A6, ridge, τ sampled, λ (1, 1) | none | none scored | none scored | no |
| A7, ridge, τ held at 1, λ (1, 1) | none | none scored | none scored | no |

**The gate's exclusions.** Every one of the 144 fits reads a predictive R̂ above 1.05. That excludes all 12 noise-free pairs and all 6 noisy ones for every arm, leaving Friedman, Ishigami, borehole and the functional benchmark unscored, each with none of its three designs kept.
- Noise-free: median 2.14, range 1.34–3.44. σ² sets the maximum in 55 of 96 fits.
- Noisy controls: median 1.74, range 1.23–2.56. The fitted values set it in 45 of 48.
- The least predictive ESS is 4.4–12.2 (median 5.4) of 8,000 draws.
- D10's verdict reads N in all three groups (predictive, level, structural) for all 144 fits.

| arm | predictive R̂, median [range] over 18 fits | least ESS, median | set by σ² |
| --- | --- | --- | --- |
| A0 | 2.16 [1.68, 3.35] | 5.1 | 6 |
| A1 | 1.92 [1.34, 2.78] | 5.6 | 3 |
| A2 | 2.20 [1.23, 3.44] | 5.1 | 10 |
| A3 | 1.85 [1.55, 2.66] | 5.7 | 9 |
| A4 | 2.43 [1.55, 3.38] | 4.9 | 7 |
| A5 | 1.76 [1.36, 2.43] | 6.0 | 3 |
| A6 | 2.29 [1.46, 3.21] | 5.0 | 11 |
| A7 | 1.72 [1.29, 2.17] | 6.1 | 9 |

Each fit's R̂, and the quantity setting it, are in the per-fit CSV.

**Context outside the rule: the ratios without the gate.** These are not a verdict, and not part of the comparison as fixed. Every fit behind them is a chain stalled in the way decision 5's trap describes. They are here because the gate decided F-R1, not the arms.

| arm | score | r_f Friedman | r_f Ishigami | r_f borehole | r_f functional | coverage − A0's | would clear |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A1 | 1.653 | 1.473 | 1.242 | 2.042 | 2.001 | −0.002 to +0.017 | no |
| A2 | 0.870 | 0.766 | 1.113 | 0.871 | 0.771 | −0.011 to +0.008 | no |
| A3 | 1.276 | 1.109 | 1.343 | 1.085 | 1.639 | −0.040 to +0.014 | no |
| A4 | 0.825 | 0.740 | 0.978 | 0.920 | 0.694 | −0.001 to +0.007 | yes |
| A5 | 1.545 | 1.271 | 1.162 | 1.951 | 1.979 | −0.009 to +0.024 | no |
| A6 | 0.598 | 0.619 | 0.943 | 0.770 | 0.285 | −0.003 to +0.022 | yes |
| A7 | 1.093 | 0.905 | 1.183 | 0.888 | 1.500 | −0.030 to +0.024 | no |

- *A0's own readings.* Its NRMSE, as a geometric mean over designs, is 0.0118 on Friedman, 0.0584 on Ishigami, 0.0084 on borehole and 0.0027 on the functional benchmark. Its coverage is 0.95–0.98.
- *The functional ratios rest on one design.* Per design, ln r runs from −3.03 to −0.37 for A6 and from −2.01 to +1.00 for A4 (sd 1.53 for both). A0's functional NRMSE is 0.0085 on design 0, against 0.0016 and 0.0015 on the other two.
- *One design left out.* A6 would still clear (0.524–0.730). A4 would not without design 0 (1.025).
- *The noisy controls, ungated.* The scores run from 0.943 (A6) to 1.111 (A2). None is at or below 0.9, and A4 reads 0.992.

**The per-fit table** is [`deterministic-profile-results.csv`](deterministic-profile-results.csv), beside this brief: 144 rows, one per fit.
- *What it holds:* every reading this brief lists, plus the run order, the holdout hash, the draws hash and the fingerprint.
- *One edit:* the run wrote D10's verdict unquoted. The committed file re-quotes that one field and changes nothing else, and the harness quotes it from `fb35f51` on.
- *The basis cap:* no fit's draw reaches it. The largest M is 128, against a cap of 499, or 1,000 for the functional benchmark.
- *M's median by arm, noise-free:*

  | arm | A0 | A1 | A2 | A3 | A4 | A5 | A6 | A7 |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | M | 24–42 | 9–34 | 28–41 | 35–76 | 27–49 | 11–42 | 31–58 | 52–108 |
- *The functional residual check:* the mean lag-one r₁ is 0.07–0.55, against −0.05 expected for independent residuals, in every arm and design. The noise-free functional fits' training residuals are dependent along t, so their intervals read after that check, as decision 9 says.

**The noisy controls.** The gate leaves out all six pairs, so both controls are unscored. With no winner, nothing is compared against A0. Their ungated scores are in the context above.

**Times.** The run took 60.4 minutes in all: the fits 20.0 minutes (1,202 s, timed), and the readings, the prediction on the holdouts, 39.5 minutes.
- *The machine:* Release, 18 logical processors, workstation GC, CPU over wall 0.770. The flat-histogram warm-start chip and other sessions shared it.
- *Against the trial:* the 1/10 trial projected about 0.45 hours of fits and readings, and the run took 0.99. Its fits took 3.3 times the linear projection (20.0 minutes against 6), and its readings 1.9 times (39.5 against 21). The run stayed well inside the 6-hour budget.
- *Seconds per arm, over its 12 noise-free fits:*

  | arm | A0 | A1 | A2 | A3 | A4 | A5 | A6 | A7 |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | seconds | 15.2 | 10.5 | 48.7 | 138.4 | 17.7 | 12.3 | 121.1 | 676.0 |

  A7 costs 44 times A0, and its functional fits take 114 s each.

**The checks, with their tolerances.** All passed before the run, and again after the review fixes.
- *Ishigami,* by Monte Carlo over 10⁶ points, within 5 SE: the mean reads 3.49477 against 3.5 (1.40 SE), and the variance 13.8558 against 13.8446 (0.51 SE).
- *Friedman's f₁ and the functional truth* against S10's generators: they differ by 0, within 1e-12.
- *The borehole function:* increasing in r_w and H_u, decreasing in H_l and L, at 1,000 random points.
- *The designs:* each Latin hypercube has one point per cell per variable.
- *The holdouts* are reproducible:

  | benchmark | hash | sd |
  | --- | --- | --- |
  | Friedman | a66dc51962b0142e | 4.884 |
  | Ishigami | 3a1cd6180768a91c | 3.703 |
  | borehole | 1a32494249af75d2 | 45.89 |
  | functional | c802ebdb47e355c5 | 0.8927 |
- *The noisy twins:* the noise sd is within 0.85–1.15 of nominal, reading 1.018 for Friedman and 0.365 against 0.37 for Ishigami.
- *The arms,* as the table says:
  - each fit's resolved form, τ₀ and held τ;
  - the Gamma(0.5, 0.5·n₀) τ prior;
  - the count prior and S10's sampling;
  - eight distinct fingerprints per (benchmark, design);
  - an arm with its priors reset gives A0's fingerprint.
- *The metrics,* on a posterior set by hand: NRMSE 1/√1.25 = 0.894427…, coverage ½, and the bands' type-7 quantiles.
- *The rule,* against a hand computation: a clearing arm, an r_f failure at 1.1696, a coverage failure, a gated design, an unscored benchmark, and the leave-one-out scores.
- *The holdout path reads the latent surface:* at the training cases it matches `BassQualification.Fitted` within 1.35e-15 relative, against a tolerance of 1e-9.
- *The 1/10 run* on one design produced every reading, and the rule's arithmetic ran. Every pair was gated there too, at R̂ 1.9–4.1.

**The review passes.**
- *Correctness.*
  - The benchmarks' formulas and domains are covered by the checks above.
  - The holdout is fixed: one hash per surface over all 48 of its fits, the noisy twins sharing their noise-free twin's.
  - The readings are latent, not predictive: they take `Predictive = false`, and the holdout-path check pins it.
  - The gate runs before the rule: the harness excludes first, and the report lists every exclusion.
  - The rule's arithmetic: an independent recomputation from the CSV, in Python, finds the same 12 + 6 exclusions, no benchmark scored and no winner. It also computes the ungated context above.
- *Coverage.*
  - A reading resting on one design: none is scored. Ungated, the functional ratios rest on design 0, as above.
  - An arm at the basis cap: no fit reaches it.
  - A benchmark the gate leaves unscored: all four, and both controls.
- *Found and fixed* in `fb35f51`, checked on a 1/50 run. None changes a reading.
  - The CSV's verdict was unquoted.
  - For an unscored arm, the noisy controls' "loses to A0" read "no"; it now reads "unscored".
  - When nothing is scored, the context tables printed empty headers; they now give a line.

**What the main session must know.**
- *The gate decided, not the arms.* Decision 5 put the 1.05 gate on the predictive group, expecting it to agree better than the structure. At S10's lengths it does not, on noise-free or noisy data: R̂ medians of 2.14 and 1.74, and a least ESS of about 5 in 8,000 draws.
  - F-M2 read the same configuration the same way: a least ESS median of 4.96, and M's R̂ of 1.96–2.45. So the cause is BASS's untempered mixing at these lengths, not this harness.
  - Every option F-R1 compared stays an option.
  - A re-run under decision 5's alternative, the self-sized tempered ladder for every arm, or with longer chains, would be a new pre-registration, and the user's decision.
- *The ungated context* shows what lies behind the gate: A6 and A4 would clear without it, and A6's lead rests most on the functional benchmark's design 0. It is reported as context, not as a result.
- *Pairing, as checked.* The arms cannot share literal starts, since each fit starts from its own prior's exact draw (S7). These were checked instead:
  - one data fingerprint and one master seed across all eight arms of a (benchmark, design);
  - the same start structures within each count prior;
  - equal τ/τ₀ starts for the two sampled arms of a count prior, to 1e-12;
  - the held arms starting at their held τ.
- *The rule's edges.*
  - Comparisons against 0.9, 1.05 and 0.05 take a 1e-12 slack, so rounding does not decide a reading at a threshold.
  - The posterior mean is the mean of the draws.
  - A gate quantity that is constant is skipped. A non-estimate reads +∞ and gates.
- *The checkout.* The run started at 13:35 on `8349663`. Its build included other sessions' uncommitted SMC and flat-histogram edits, none on the untempered MCMC path:
  - `BassResolved.cs`, `BassSmcRun.cs`, `CutSmc.cs`, `SmcRun.cs` and `SmcWeights.cs`.

  The edits to `BassBasis.cs`, `BassDesign.cs` and `BassModel.cs`, and the new `DesignKernel*.cs` and `BassObservationModel.cs`, came after 14:07 and are not in the run. The post-review checks and the 1/50 run did build with them.
- *A stale message.* A message saying none of this work was in the checkout reached this chip after `8349663` had been committed.
