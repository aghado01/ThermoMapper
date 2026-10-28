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
