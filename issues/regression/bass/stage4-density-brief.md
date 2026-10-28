# BASS stage 4: the pointwise density and residual dependence, S9 (chip B brief)

Date: 2026-10-28. A chip of the [fitting segment's stage-4 plan](port-fitting-brief.md#stage-4-plan-decided-2026-10-28), which the user decided on 2026-10-28. S1–S6 are built and reviewed (the [checklist's reports](port-fitting-checklist.md#reports)). This chip builds S9 ahead of the run: the pointwise log density per draw (OQ17's seam) and the residual-dependence check (F11). Both read only retained draws and the model, whose shapes are fixed (`BassDraw`, `BassModel`, `BassTermTable`). So you build them as static entry points over draws, and the main session wires them into `BassFit` after chip A, which builds the run (S7), lands. Chip C writes the docs section meanwhile.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then, in the [fitting segment](port-fitting-brief.md):
- "Per-observation log density per draw (OQ17)": the seam's design, units and latents' scope;
- D04's "Pointwise density default (OQ17)";
- the residual-dependence proposal at the segment's head and F11 in the checklist;
- the contract shapes, which give `BassFit`'s `ResidualDependence?` field.

In the code:
- the shared seam: `src/maths/regression/Pointwise.cs` (`UnitPartition`, `PointwiseContent`, `LatentReport`, `PointwiseLogDensity`) and `PredictiveChecks.cs` (WAIC, PSIS-LOO, r_eff from the chains);
- BARS's producer, the template: `BarsEnsemble.PointwiseLogDensity`;
- the collapsed content's block: `GPriorConditional.LeaveOutLogDensity` in `src/maths/regression/gprior/CollapsedGPrior.cs`;
- the draws' basis: `BassDraw`, whose coefficients are uncentred, a₀ + Σ aⱼBⱼ (`BassState.cs`).

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-density\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-density\` when you finish.
- **Run oracle tests filtered** to BASS if you run any; the TDA oracle tests time out when several builds run at once.
- **Your files:** new `src/maths/regression/spline/bass/BassPointwise*.cs` and `BassResidualDependence*.cs`, and their new tests (`tests/regression/BassPointwise*Tests.cs`, `tests/regression/BassResidualDependence*Tests.cs`). Generate test draws with `BassChain` and `BassStructureMoves.Compose`, which exist. Do not edit other BASS files: chip A edits `Bass.cs`, `BassSweep.cs`, `BassState.cs` and `BassStructurePrior.cs`. If you need a hook there, say so in your report. A change to shared code stays generic, is tested in its own host, and is named in your report.
- **Commit your own files only, by explicit path,** one commit per task, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** C# and briefs are LF.
- **Decisions are not yours to settle.** Where the decided docs leave something open (below), implement the recommendation only if it is reversible, and put the question in your report.
- **Report** under "Reports" at the end of this brief, newest first: what landed, the commits, the tests and suite counts, the review passes (below) and anything the main session must know.

## Task 1: the pointwise log density

**The entry point.** A static entry point over the model, the term table and the draws by chain, giving a `PointwiseLogDensity` with its chains for PSIS's r_eff. Follow BARS's producer for the shape.

**Units.** The caller's `UnitPartition`. The default is one observation for a scalar response, and one case (its whole curve of m values) for a direct functional response: the hub's unit for functional data. Validate the partition against N.

**Contents.**
- *Retained parameters, the default (D04).* The Normal density of each unit's observations at the draw's coefficients and σ². The noise is independent across output locations, so a case's density is the product over its curve. The structure, τ, a and σ² are global latents, conditioned on.
- *The collapsed predictive, selectable.* For each draw's structure and τ, the leave-unit-out predictive with a and σ² integrated, through the shared block's `LeaveOutLogDensity` on the structure's explicit design. Under this content a and σ² are integrated.

  Its cost is a design per distinct structure. Note it, and cache by structure if the draws repeat.

**The check D04 asks for before a comparison relies on the collapsed content.** BARS's brute-force leave-one-out found its collapsed k̂ pessimistic under clumping. So, on a small fixture, check the collapsed content against explicit refits that leave each unit out, at the same structure and τ.

**Checks.**
- Both contents against dense computation: the explicit design and Normal or Student-t densities.
- Units of observations and of cases on direct functional data.
- WAIC and PSIS-LOO (`PredictiveChecks`) read the result, with r_eff from the chains.
- The latent reports say what was conditioned on and what was integrated.

## Task 2: the residual-dependence check (F11)

**What it reports.** For a direct functional fit with one output variable: each case's residual under the posterior mean, the fitted values averaged over draws, and its lag-one autocorrelation along the output index ordered by t. Report a summary across cases (mean and quantiles) beside the per-case values: the `ResidualDependence` type the contract shapes name.

**Where it is null.** A scalar response gives null.

**Several output variables are open.** No order along the outputs is decided. Return null there, and put the question in your report with a recommendation.

**Checks.**
- A deterministic smooth misfit, a direct fit whose structure cannot reach a smooth curve in t, reads a high lag-one value.
- Independent noise reads about zero, within its sampling error.
- The values match a direct computation.

## Review before reporting

A task is not done when its tests pass. Before the report, review the work as a separate pass, as the main session reviewed S5 and S6 (each of those reviews found real defects or untested paths):

1. **Reread the whole change** as a reviewer seeing it for the first time: `git diff <the commit before your first>..HEAD -- <your paths>`. For each change, ask:
   - Is it correct? Check the densities and their normalizing constants, the unit bookkeeping, the edge cases (a unit with one observation, a structure with no terms, the flat intercept and Jeffreys under the collapsed content), the exception paths, the parameter names in exceptions, thread safety and the seam's contract.
   - Does it read like the surrounding code?
   - Do its docstrings say what it does?
2. **Hunt the coverage gaps.** Look for each option and convention the code accepts that no test runs:
   - the intercept conventions (level centring, flat, zero);
   - the σ² priors (scaled, Jeffreys, a declared IG);
   - sampled or fixed g;
   - scalar, categorical and direct functional data;
   - both contents, and each unit partition.

   Look also for each refusal path. Prefer independent references (dense computation, explicit refits) to self-consistency.
3. **Fix what you find,** add the tests and rerun. Commit the fixes separately (`fix(bass): …` or `test(bass): …`).
4. **Review again** if the pass found anything substantive. Stop when a pass finds nothing.
5. **Record it** in your report: what each pass found, what was fixed, and anything left, with the reason.

## Reports

### Chip B: the pointwise density and residual dependence (S9) — 2026-10-28

Both tasks are built, tested and reviewed in three passes. Commits: `60dba78` (task 1), `b819126` (task 2), `fd107f4` (the first review's fixes), `2db92dd` (the second review's tests). Only the four new files changed, and no shared code.

**Task 1: `BassPointwise`** (`src/maths/regression/spline/bass/BassPointwise.cs`).
- *The entry point.* `BassPointwise.LogDensity(model, terms, draws, units?, content?)` returns a `PointwiseLogDensity` with the chains' lengths, which PSIS reads for r_eff. `draws` is per chain, and a fit's `IReadOnlyList<DrawStore<BassDraw>>` passes as it is (covariance, tested).
- *Units.* `BassPointwise.DefaultUnits(model)` gives each observation of a scalar response, and each case of a direct one. Observations run case by case, so a case is `ConsecutiveBlocks(n, m)`. The partition is checked against N.
- *Retained parameters* (the default): the Normal density at the draw's uncentred coefficients and σ². The latents structure, τ, a and σ² are all conditioned on. Under fixed g τ is held, so it is not reported as a latent.
- *The collapsed predictive.* Draws are grouped by structure. Each distinct structure gets one `BassDesign.Load`, one block (`FromCentredLeastSquares`) and one explicit [1, Z_c], and only one structure is held in memory at a time. Each draw then costs N solves through `LeaveOutLogDensity` at its τ. a and σ² are reported `IntegratedGivenOthers`.
- *Refusals.* A partition that does not cover N; content out of range. Draws: no chain, an empty or null chain, a null draw, ids outside the table or not ascending and distinct, a coefficient count other than M + 1 or a non-finite coefficient, a σ² that is not positive and finite. τ is checked the same way, for the collapsed content only. A term that does not fit the model's inputs (named `terms`).
- *Collapsed-content refusals.* A structure outside F: the design's refusal, or an exact fit under Jeffreys. A unit that leaves no observation outside it, beyond the one the flat intercept spends. Under Jeffreys that is σ²'s shape given the rest reaching 0. Under a proper prior the shared block would otherwise throw on S₋U = 0, with an exception that names `qr`.

**Task 2: `ResidualDependence`** (`BassResidualDependence.cs`).
- *The record.* `ResidualDependence(LagOne, Mean, Probabilities, Quantiles)`: r₁ per case, and the mean and type-7 quantiles over the cases where r₁ is defined. The default probabilities are 5, 25, 50, 75 and 95%.
- *The estimator.* `ResidualDependence.Of(residuals, index, probabilities?)` is generic. It orders the locations by the index, keeping ties in location order, and computes r₁ = Σ(eₖ − ē)(eₖ₊₁ − ē)/Σ(eₖ − ē)². A case whose residual does not vary gets NaN, and the summary leaves it out.
- *The BASS check.* `BassResidualDependence.Of(model, terms, draws, weights?, probabilities?)` takes the residual under the posterior mean, with SMC's weights when they are given, ordered by the output variable. It returns null for a scalar response and for several output variables. The arguments are checked even where the result is null.

**Checks** (`BassPointwiseTests` 62, `BassResidualDependenceTests` 6, all passing).
- *Retained content.* Against `Normal.LogPdf` at the explicit uncentred design, to 10⁻¹⁰ relative.
- *Collapsed content.* Against explicit refits that leave each unit out, at each draw's structure and τ, to 10⁻⁸ relative. Each refit is the normal–inverse-gamma posterior from the rest's rows under the full data's prior, then the unit's multivariate t. This is D04's check.
- *Coverage.* Every combination of the three intercepts, three σ² priors (scaled, Jeffreys, IG(2, ½)), two g modes and three kinds of data (scalar, categorical, direct functional): 54 cases, each with 16 `BassChain` draws over 2–12 structures. Also:
  - partitions: scalar observations, and pairs with a one-observation unit; functional cases, observations, and locations across cases (not contiguous);
  - the empty structure under each convention;
  - exact fits: refused under Jeffreys, matched under the scaled default (the deterministic-simulator regime);
  - `DrawStore`s with unequal chain lengths;
  - a `Bass.Fit`'s own merged table and stores.
- *WAIC and PSIS-LOO* read both contents, with r_eff estimated from the chains. WAIC refuses the collapsed content. As observations, not assertions (3 chains × 150 draws):
  - scalar fixture: retained elpd −9.89 with 3 units above k̂ 0.62, collapsed −10.14 with 9;
  - functional fixture: −0.12 with 3, and −0.71 with 2.
- *The residual check.*
  - A hand-worked case, with ties and an unordered index.
  - A direct computation with t stored out of order, under equal and given weights (one of them zero), to 10⁻¹².
  - A deterministic sine fitted by a case-only structure reads exactly the sine's own r₁, 0.9458.
  - Chains capped at two terms read a mean of 0.78, with a 5% quantile of 0.71.
  - Independent noise reads −0.089 ± 0.044, against −1/m = −0.062, with a per-case sd of 0.251 against 1/√m = 0.250.

**Suite.** `Maths.Regression.Tests` in `artifacts\chip-density\`: 651 of 651 at `b819126`, and 654 of 654 after the review passes (HEAD `a8cfa36`, which includes `2db92dd`). Both runs saw other sessions' uncommitted edits in the shared checkout.

**Review passes.**
1. *Pass 1, after both tasks.*
   - Defects fixed: a non-finite coefficient produced NaN densities, and is now refused; the generic check read a non-finite residual as "does not vary", and now refuses it; the null check on `draws` ran after the partition's; the record mixed array and list fields.
   - Gaps closed: the collapsed content's exact-fit refusal under Jeffreys; exact fits under the scaled default; `DrawStore` input with unequal chains.
   - Building the tests found the unit refusal above: a unit with nothing outside it reached the shared block's `qr` exception under proper priors. It is now one rule for every prior.
   - Committed as `fd107f4`.
2. *Pass 2.* No source defects. Gaps closed: the refusals of a null chain and of a structure without ids, and the end-to-end test on `Bass.Fit`'s output; one vacuous assertion fixed. Committed as `2db92dd`.
3. *Pass 3.* Nothing found.

**For the main session (wiring).**
- *The fit holds no data or model.* `BassFit.PointwiseLogDensity(units, content)` needs one: `BassModel.Create(fit.Config, data)` rebuilds the inputs, basis and target that the table's knot and level indices refer to. The end-to-end test does exactly that.
- *The `ResidualDependence` field* is `BassResidualDependence.Of(model, fit.Terms, fit.Draws, fit.DrawWeights)`, at the end of the run.
- *Ids must be ascending.* Both entry points require it; `BassRun`'s remap keeps it by re-sorting each draw's ids with its coefficients.
- *Weights.* `PointwiseLogDensity` and its consumers treat draws as equally weighted, so SMC's weighted draws (stage 8) will need weights there. Not built.
- *Cost.* The shared block re-solves the τ-free rows for every draw at a structure. Caching them per structure would need a hook on `GPriorConditional`; not needed at this chip's sizes.
- *The shared domain.* `LeaveOutPredictive.StudentT` refuses S₋U ≤ 0 under every σ² prior, though a proper prior keeps the predictive finite there. BASS's unit rule keeps every reachable case inside its domain, so the shared code is unchanged.

**Questions** (each implemented reversibly, or left null):
1. *Several output variables (F11).* The check returns null, as the brief says. The output variables mix coordinates and covariates, and the data contract does not say which is which. Recommendation:
   - declare the index on the call, an output variable or an explicit location order, and read r₁ along it;
   - for d ≥ 2 coordinates, read Moran's I per case over each location's k nearest neighbours in the scaled coordinates, reported like r₁.
2. *r₁'s centring.* The estimator is the textbook one, about each case's own mean residual, so a constant offset in a case does not read as dependence along t. An uncentred r₁ would read it. The choice sits in one line.
3. *D04's elpd-level check.* This chip checked each collapsed term against explicit refits. BARS's brute-force leave-one-out, the check that found the collapsed k̂ pessimistic, has not been run for BASS. That reference samples structures without the unit, on a target that subtracts the leave-out term. `BassChain` takes a custom target, so it can run now. Recommended before a comparison relies on the collapsed content's k̂.
