# The BASS completion audit (chip brief)

Date: 2026-11-01. Status: chipped on the user's word. An audit of the port against its contract, ending in the completion report the hub asks for ([bass-full-port-brief.md](../bass-full-port-brief.md), the last paragraph of "References and documentation updates").

## Why

Every stage is built, and each was verified when it landed, but each against its own brief. Nobody has read the whole against the contract. The first note on what remained ([opus-bass-completion-remainder.md](../opus-bass-completion-remainder.md)) named edge cases nobody could vouch for, the categorical representation above 64 levels among them. This audit turns "close to done" into a concrete list: what is fulfilled and on what evidence, what is built but untested, what is missing, and what was left out by decision.

## The standard of evidence

**Start from R and the contract.** The contract is the hub's "Scope and completion contract". The R source is BASS 1.3.1 at `D:/aghado01/packages/BASS`, commit `26710146c623a258854162367a939041efedd091`, and [the public inventory](../bass-public-inventory-and-backend-review.md) maps it.

**Our records say where to look; the evidence is code and tests you have read.** The checklists' reports and the "Checked" bullets in `docs/Maths` were written by the work they describe. A verdict that rests on a record alone is not yet a verdict.

Each row takes one verdict:
- **Fulfilled:** the C# location (`file:line`); a named test you have read, which would fail if the capability broke; and its reference: R parity where R is right, the corrected value where R is wrong (the hub's parity rule), exact enumeration, closed form or controlled recovery.
- **Fulfilled, limited:** as above, with the measured limitation stated.
- **Built, untested:** the code path exists and no test exercises it.
- **Gap:** missing, or it differs from the contract.
- **Left out by decision:** cite the decision (a D-number, the parity rule, the hub's list of extensions) and where it is recorded.

## The audit, in order

Each part lands in `issues/regression/bass/completion-audit-matrix.md`, committed as it finishes, so a restart loses little. Record the audited revision (`git rev-parse HEAD`) first. Other sessions commit beside you, so note any later commit that touches BASS.

1. **Capabilities P01–P15.**
   - For each exported function and S3 method, read R's signature and return. Map every argument to its C# field, and every returned element to its C# member, or give the other verdicts. The inventory's "Fitting controls that must remain representable", "Three response workflows" and "Sensitivity, calibration, and operational coverage" are the checklist.
   - Cross the three response representations (scalar, direct functional, supplied basis or PCA) with continuous, categorical and mixed inputs, and with each output: prediction at new inputs and output coordinates, sensitivity, effects, reconstruction, restart. A cell that cannot apply says why.
   - The plot, print and summary methods (P09–P15): the contract asks for their data in regression, with rendering above maths. Say which data exists; rendering is listed as remaining, not as a gap.
2. **The baseline.**
   - The qualified existing-engine revision and the closure evidence for B01–B06 ([the prerequisites review](../bass-adaptation-prerequisites-review.md)), keeping baseline repairs distinguishable from port changes.
   - The inventory's "Observed contract defects and qualification gaps": each R defect has its fixture, the C# test asserts the corrected value, and the R oracle keeps R's value beside it, as the parity rule says.
3. **The decisions.**
   - D01–D10, from the hub's register and the segment briefs: each decision as recorded, its code, and its test.
   - **The defaults, mechanically:** every row of the BASS configuration tables in `docs/Maths/BayesianRegression.md` and the BASS-facing rows in `docs/Maths/Samplers.md`, against the value in code. List each difference.
   - The register's **Corrected** entries.
4. **The validation policy** (the hub's "Validation and performance qualification"), each item with its test, across the representations where it applies:
   - the independent checks: boundary proposal ratios, enumerable model spaces, Gaussian conditionals, precision sampling, empirical covariance, functional contractions, sensitivity identities, calibration targets;
   - the SMC tests: a known initial distribution, incremental weights, normalizer ratios, resampling, state ownership, populations whose modes are known;
   - the essential edge cases: intercept-only models; one predictor, one response component, one case, one draw; constant columns and constant responses; nearly dependent terms; more than 64 category levels; category probabilities of zero; bounded and truncated input measures; full PCA rank; rejected proposals leaving the current state untouched.
5. **The benchmark list.** The hub's axes and measures against [the cost half's report](performance-cost-brief.md#reports). Known open: statistical error per unit time ([its brief](performance-statistical-brief.md), not run) and the sensitivity moment work alone. Run no benchmark.
6. **The records.**
   - The open checklist boxes (N15, N17, F-M2, F-X4 at the time of writing): why each is open, and whether its box is stale.
   - The API and source-map rows in `docs/Maths/BayesianRegression.md` against the public surface and the files: rows naming what no longer exists, and public types no row names.
   - Statements you found false along the way. The extension seams' write-up is already known to be owed, so list it once.
7. **What remains, and the limits.** The optional backends and extensions ([the extensions brief](extensions-brief.md); the hub tracks the SMC backends apart from capability completion), and the measured limitations, each with where it was measured: Friedman's posterior, the deterministic profile (F-R1), N15's cost and any others you meet.

## Probes

Where the contract names something and no test exists (edge cases above all), a small probe separates "built, untested" from "gap".
- They go in one file, `tests/regression/BassCompletionAuditProbes.cs`, each a few seconds at most, each named for what it pins.
- A probe that passes stays as coverage.
- A probe that exposes a defect goes in as `[Fact(Skip = "audit gap G-n")]`, with its reproduction in the report.
- Nothing under `src/` changes. A gap is reported for the user to schedule, not fixed here.

## The machine

The user's machine is shared with other sessions and applications, and it jammed under load on 2026-11-01. So this chip stays light:
- Use targeted filters while you work, and run the suites once at the end, one at a time: `Maths.Samplers.Tests`, `Maths.Regression.Tests`, and the BASS classes of `Maths.Oracle.Tests` (`FullyQualifiedName~Maths.Oracle.Tests.Bass`).
- Run every build and test at below-normal priority, so its children inherit it, and cap the work: `-m:2` on builds, `-- xUnit.MaxParallelThreads=4` on tests.
- Leave `BASS_QUALIFICATION` and `SEED_SWEEP` unset; those runs are the user's to schedule.
- Check that a run has really stalled before repeating it. These suites took 5–6 minutes each on an idle machine, and take far longer on a busy one.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing, and spawn no subagents.
- **Your files:** the matrix, the probes, the report and this brief's Reports section. `docs/Maths`, the hub and the checklists are the main session's, which verifies your report and records it. Leave other sessions' uncommitted edits and the user's notes under `issues/regression/` alone.
- **Commit your own files only, by explicit path** (`git commit -- paths`), never leaving anything staged, since the index is shared, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Build and test into your own artifacts root:** `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-audit\`, on each command itself rather than through a variable set earlier in a chain. Delete `artifacts\chip-audit\` when you finish.
- **Search with `rg` through Bash;** the Grep tool is denied here. `docs/Maths/BayesianRegression.md` runs past 2,000 lines, so read it by section; the `mdnav` skill helps.
- **Line endings:** briefs and C# are LF.
- **Recommend no default.** Which gaps to fix, and in what order, is the user's call; give each gap its size and leave the choice.

## Review before reporting

A fresh-eyes pass, its fixes made before you report:
- **Completeness of the audit:** every P-ID, B-ID, D-number and edge case above has a row, and so does every R argument and return.
- **A re-trace:** take every tenth Fulfilled row and trace it again from R to the test, without reading your note first. Ask whether that test would fail if the capability were removed. If a row fails this, widen the check to its neighbours and report how many moved.
- **Demotion:** any verdict that rests on a record alone drops to what you can show.

## Report

`issues/regression/bass-completion-report.md`, beside the hub and linking back to it, with what the hub lists:
- the audited revision, and the suites' counts there;
- the qualified existing-engine revision and the B01–B06 closure evidence, baseline repairs apart from port changes;
- the fulfilled capability IDs;
- the resolved statistical decisions;
- the remaining optional backends and extensions;
- the measured limitations;
- **the gap list, ranked:** contract gaps first, then built-but-untested, then record drift, each with a one-line reading of its size.

Open it by saying plainly whether the contract is fulfilled, and if not, what stands between. Keep what you found apart from what you make of it. The matrix holds the row-by-row detail. Put a short note under "Reports" below, newest first, with what the main session should check.

## Not in this brief

- **Any fix under `src/`,** and any change to a default.
- **Performance qualification's statistical half** and the efficiency pass.
- **The extensions** (E1–E18), the notation pass, and rendering above maths.
- **The Friedman band reading,** a separate chip ([its brief](friedman-smc-band-brief.md)).

## Reports

Newest first.

### The closing, checked by the main session — 2026-11-02

- *Read:* `src/` across `0ab1e4d`, `4c34948`, `1364cbb`, `b844f53`, `3a83d0a` and `26befff`, and the tests that carry them.
  - **G-2 against R:** `vf[vf==func.var]<-NA` (`sobol_normTmix.R:684`) removes only the variable conditioned on, so R keeps the others as inputs, as `ConditionOn` does. The normalizer is one constant over the part, carried once: `ConditionOnOne_HandBuiltTerms_CarryTheNormalizerOnce` builds the terms where a split would show, a coordinate past a hinge's end among them. Terms with no output factor get a normalizer of 1, so skipping `Reads` there is harmless.
  - **G-3 on a weighted fit:** the selected draws keep raw weights that no longer sum to 1, which holds because the means divide by their mass and the summaries read a weighted sample.
  - **G-4:** `ListFault` refuses NaN, a list not starting at 1, one not descending, and the wrong end for the reference rung.
  - **G-6:** the default draw is bit-identical, a part left out draws nothing, and `Reconstruct` takes an empty residual.
- *The matrix, by my script again:* F 192, FL 5, D 18, G 0, total 215. Every R argument and return still has a row.
- *The suites, rerun by the main session at `971c84b` (Release, below-normal priority, four threads): Samplers 284 of 284, Regression 1,544 of 1,544 with none skipped, the BASS oracle classes 90 of 90, the chip's counts.*
- *The chip's correction stands:* a list moves only a tempered config's fingerprint, since `Tempering` is null by default. So only a tempered checkpoint from before `3a83d0a` stops resuming, a narrower break than I told the user. No hash re-pinned.
- *One choice put to the user:* a declared list on the default `Adaptive = true` re-spaces its interior rungs. `Levels` and `HotEnd` declare a count and an end, but a list declares the spacing itself, and on this default it does not run as declared.
- *Written by me from the chip's list:* the API, configuration and prose rows for `InverseTemperatures`, `ConditionOn`, `Draws`, `BasisParts`, `OtherVariance` and `OtherIndex`, in `docs/Maths/BayesianRegression.md`.
- *Open:* G-7 and G-8, waiting on a quiet window; the extension seams' write-up, mine, owed.

### The gaps closed (chip, 2026-11-02)

G-1 `0ab1e4d`, G-2 `4c34948`, G-3 `1364cbb` with `b844f53`, G-4 `3a83d0a`, G-6 `26befff`; G-5 left out, as ruled. [The matrix](completion-audit-matrix.md) and [the report](../bass-completion-report.md) are flipped: F 192, FL 5, D 18, G 0, with P04, P05 and P08 Fulfilled, and the report's open questions turned into the rulings. G-7 and G-8 wait on a quiet window. At `b844f53`: Samplers 284 of 284, Regression 1,544 of 1,544 with none skipped, BASS oracle classes 90 of 90. `artifacts\chip-sens\` is deleted.

**For the main session to check:**
- *`BassMovePassTests`' hash did not move, against your note.* `new BassConfig()`'s JSON has `"Tempering":null`, so the new field sits inside a null. The pin and the bit-identity tests pass unchanged, so I re-pinned nothing. Only a tempered config's fingerprint moves, so only a tempered checkpoint from before `3a83d0a` stops resuming, a flat-histogram pilot's ladder included.
- *Choices the rulings did not spell out:*
  - **The list and `Adaptive`.** They stay independent, so on the default `true` a list is a start whose interior rungs re-space. R's as-given list is `Adaptive = false`.
  - **`BasisParts` is one enum of three,** not R's two booleans. R's "neither" is `Latent`, and R's default pair is `Noise`. D07's full draw stays the default.
  - **`ConditionOn` takes several names,** wider than R's one.
  - **G-2's default coordinates** are the held variables' distinct training values. This also applies when every output variable is held, where it was one coordinate per location, and differs only if locations repeat.
  - **`other` snaps to 0** within the rounding tolerance and on a zero-variance draw.
- *What the review changed:* basis effects evaluated every component draw and then picked the selected ones. They now evaluate only those (`b844f53`), so a selection fits a budget the whole run would exceed.
- *Not compared with R:* the new options have no R fixtures. Their references are quadrature, the whole analysis's draws, bit identity with the built ladder, and hand-assembled covariances.
- *Docs owed, yours:* rows in `docs/Maths/BayesianRegression.md` for `TemperingSpec.InverseTemperatures`; `SensitivitySpec.Draws` and `ConditionOn`; `BassPredictOptions.BasisParts` (its row lists the record's fields) and `BassBasisPrediction.Parts`; `CoordinateVariables` on `SobolAnalysis` and `EffectAnalysis`; `OtherVariance`, `OtherIndex` and `MeanOtherIndex` on `SobolDecomposition`, and `OtherIndex` and `OtherVariance` on `SobolSummary`.

### Checked by the main session — 2026-11-01

- *The gaps against R and the code:* each of G-1 to G-6 holds at the lines cited: `V.other` and its `names.ind` (`bassPCA.R:590-603`), `func.var` (`sobol_normTmix.R:215-225`) against the refusal at `sensitivity/BassSensitivity.cs:391`, no draw field on `SensitivitySpec`, the ladder's list on an internal constructor (`BassLadder.cs:56`), `start.temper` (`bass.R:39,78`), and `trunc.error` and `nugget` (`bassPCA.R:183`) against the one `Predictive`.
- *The audit's completeness, by a script of my own:* `NAMESPACE` holds the six exports and nine S3 methods the hub counts; every argument of the eleven R functions and every name in `bass`'s returned lists appears in the matrix. Its verdicts recount to the report's 215: F 183, FL 12, D 17, G 3.
- *Rows re-traced:* the defaults (a dozen rows' values against `BassConfig.cs`, no difference), the 64-level representation, and `DeclaredControls_ReachTheModelAndTheRetainedTerms`, which reads the proposals, the count ratio and the retained terms, not the config's values alone.
- *The suites, rerun by the main session at `e1eea5d`* (Release, below-normal priority, four threads, before the chip began on G-1): `Maths.LinAlg.Tests` 139 of 139 and `Maths.Distributions.Tests` 148 of 148, which the chip had read but not run; Samplers 284 of 284; Regression 1,506 of 1,507, the one skip being G-2's probe, as the chip read. The BASS oracle classes I take from its run (90 of 90); no `src/maths` change separates the two runs.
- *Two additions to the report:*
  - **G-1's definition.** R's `other` is the variance less the subsets it keeps, after `nind` keeps the largest (`bassPCA.R:582-590`). The port has no `nind`, by decision (`port-outputs-brief.md:88`), so here it is the variance less every subset up to `MaxOrder`. They agree where `nind` is null. `0ab1e4d` defines it that way, the variance less the subsets reported.
  - **`s2.lower` is the one Left-out row worth reopening** (A14, D). It was left out at the contract shapes on grounds of form ("a truncated σ² prior, which would have to be stated as one"), not of merit, and before F-R1. The inventory lists it among the controls that must remain representable. F-R1 then read every noise-free fit failing its predictive gate, its brief naming σ² running small as the trap, and its arms varied the coefficient prior and the count alone. So a truncated σ² prior is untested where it would matter. It changes the target (the collapsed marginal gains a regularized incomplete gamma factor), so it is a pre-registered arm, not one more field. The user's ruling.
- *My readings on what the chip asked:* a refusal with its reason meets "one case": one case makes every input constant, and two fit. The per-iteration records stay limited, not gaps: they are not a draw's generative state.
- *Since the report:* the user set this chip to closing the surface gaps in its own session (G-1 `0ab1e4d`, G-2 next), on its plan: build G-3, G-4 and G-6, record G-5 as left out, since the self-sizing ladder reads its swap rejections during burn-in and a late start would starve it. I agree with those leans. They are checked when it reports.
- *The user's rulings, the same night* ("Yes, take them", on the leans in my message): G-3, G-4 and G-6 are built; G-5 is left out; the per-iteration records stay limited; `s2.lower` is a pre-registered arm at F-R1's rerun. Recorded in `port-fitting-brief.md` ("Not carried"), `deterministic-profile-brief.md` ("For a rerun") and `docs/Maths/BayesianRegression.md`. The extension seams' write-up in `docs/Maths` is still the main session's, owed.
- *Fixed from "Record drift":* the source map's three files, `BassComponents.Assemble` in the API row, the five false statements (annotated in place, or fixed in `BassRunOptions.cs` and this brief). G-8 joined the statistical half's cost blocks. The 71 public types no API row names stay as listed: the rows name entry points, and most of those are readouts.

### The audit done (chip, 2026-11-01)

[The report](../bass-completion-report.md) and [the matrix](completion-audit-matrix.md); probes in `tests/regression/BassCompletionAuditProbes.cs` (`348ad80`). Audited at `8f4d2da`. Samplers 284/284, Regression 1,507 (1 skipped: G-2), BASS oracle classes 90/90. Nothing under `src/` changed. `artifacts\chip-audit\` is deleted.

**Found:** the contract is not yet fulfilled in full. Every capability exists and verifies, the 215 rows found no wrong number, and the defaults match the doc at every row. Six small items of R's surface differ (`other`, `func.var` on several output variables, `mcmc.use` on Sobol, `temp.ladder` as a list, `start.temper`, the two prediction switches), and two benchmark measures stay open, as known.

**For the main session to check:**
- *Three rulings the brief leaves open* (the report's "What the brief does not rule on"): they decide whether G-3 to G-6 are gaps or recorded choices. I decided none.
- *Whether "one case" is met by a refusal* (Ed04, Fulfilled, limited). D01 and D07 imply it; the hub's list names the case without saying.
- *The re-trace's finding:* declared controls read only at their defaults (seven rows, now covered by a probe). A similar pattern may exist beyond BASS's config.
- *The run's tree:* the main checkout, with another session's `cb0400b` in it. Its one test run was not a clean snapshot, because the oracle's R library is ignored by git.
- *Records to fix,* yours by the brief: the report's "Record drift", including the source map's three files and five false statements, one of them in this brief (E1–E15).
- *`Maths.LinAlg.Tests` and `Maths.Distributions.Tests`* were read, not run here, and four rows cite LinAlg tests.
