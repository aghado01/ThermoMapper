# BASS port: reference fixtures that do not wait on the design pass (chip brief)

Date: 2026-10-27. Issued from the BASS planning session, which continues the D01–D06 design discussion while this runs. A second chip, [shared-infrastructure-brief.md](shared-infrastructure-brief.md), runs at the same time in the same checkout; the file ownership below keeps the two apart.

## Goal

Build the stage-1 reference fixtures whose expected values do not depend on how D01–D09 settle. They are R oracles in `r/oracles/`, run by `r/scripts/bass_stage1_checks.R`, in the style of the four already built (F-O1, F-O2, F-O7, O-O4). Each one either computes an independent mathematical reference or reproduces an R BASS behaviour from the pinned source next to its corrected expectation, so it can serve as a defect fixture. No C# yet: parity tests come with the build items that consume each fixture (F3, F4, F6, O1–O7, C1–C3).

Done means every item below has its oracle, its checks pass under `Rscript scripts/bass_stage1_checks.R` (cwd = `r/`), its checklist item is ticked with a report, and this brief's Report section says what landed.

## Items

| Item | Checklist | What the fixture establishes | Where the spec is |
|---|---|---|---|
| F-O5 | [fitting](port-fitting-checklist.md) | R's `dmwnchBass` (`mcmc_funcs.R`) against brute-force sums over draw orders: each subset's probability as the sum, over every order of its variables, of the product of sequential weighted selections without replacement. Interaction orders 1–3 (R's range) and 4–5 (exactness beyond it); equal weights (the central case) and strongly unequal ones; several predictor counts. D01 states that R's formula is exact, not an approximation; this is the check. | D01 draft in [port-fitting-brief.md](port-fitting-brief.md), the variable-proposal bullets |
| F-O6 | fitting | Three defect fixtures. (a) R's tempering: births and deaths multiply prior and proposal factors by the inverse temperature (`rjmcmc_des.R`), and the swap density `lp` (`mcmc_funcs.R`) omits the term-prior blocks; on an enumerable structure space, show that R's swap fails detailed balance for every single tempered target, beside D05's corrected ratio. (b) R's Jeffreys mode: reproduce its per-term factors and its σ² shape from the source and set them against the flat and Jeffreys coefficient priors' coherent factors (the derivation is in D03). (c) KHAOS's double-counted (M+1) in the birth ratio: under a constant likelihood, the count law each birth ratio implies, beside the stated prior's. | D03 ("R's Jeffreys mode has no single target"), D05 ("R's tempering is a defect fixture"), D01 (term measure; KHAOS's form) |
| O-O1 | [outputs](port-outputs-checklist.md) | R output-path defects: the truncation sampler's indexing, recycled noise, the residual's sign and units, 100% retention, constant columns, draw order and dropped axes, interval quantiles. Each beside its corrected expectation. | D07 draft in [port-outputs-brief.md](port-outputs-brief.md); the inventory review's [defect section](../bass-public-inventory-and-backend-review.md#observed-contract-defects-and-qualification-gaps) |
| O-O2 | outputs | R sensitivity defects: ignored category weights, unclipped truncation bounds, decomposition units, hard-coded effects. | D08 draft; inventory review, same section |
| O-O3 | outputs | Hinge-product moments against quadrature and category enumeration: hinge degrees above 1, the truncated-normal mixtures R's `sobol` supports (`sobol_normTmix.R`), bounded measures, and category probabilities including zeros. The reference is numerical integration or enumeration, never the closed form being checked. | D08 draft |
| C-O1 | [calibration](port-calibration-checklist.md) | R calibration defects (`calib.R`): the stale likelihood after a discrepancy draw, the branch chosen by shape, normalized returned coordinates, the unit mismatch, the all-groups `s2.df` switch. | D09 draft in [port-calibration-brief.md](port-calibration-brief.md); inventory review |
| C-O2 | calibration | Reference values for both covariance formulations: dense solves and log-determinants against the low-rank route, and the physical posterior unchanged between scaled and unscaled y (the corrected statement; R's unit mismatch is C-O1's). | D09 draft |

Out of scope: F-O3 and F-O4 (count and order laws, boundary ratios), which test D01's and D02's drafted measure and wait for the design pass.

## How to work

- **The R source is pinned.** BASS 1.3.1 at `D:/aghado01/packages/BASS`, commit `26710146c623a258854162367a939041efedd091`; verify the commit before relying on it. R BASS is not installed in the project's renv. `source()` the pinned files you need (they are plain R) instead of installing the package.
- **Packages.** Base R and `jsonlite` only. If a fixture would need another package, stop and record it in the report as a blocker: installing one is a download that needs the user's approval, and the other chip owns `r/DESCRIPTION` and `r/renv.lock`.
- **Oracle conventions.** Follow the existing BASS oracles: read a fixture file, emit JSON through `_common.R`'s helpers, keep every check in `bass_stage1_checks.R`, and say in each oracle's header what it is the reference for. Judge numerical ranks against one tolerance for the whole problem (O-O4's lesson). Compute residuals directly, not as yᵀy − Q (F-O1's lesson). Sum cross-products pairwise where cancellation matters (F-O7's lesson).
- **Independence.** A reference must come from outside the formula under test: enumeration, brute force, quadrature, or an explicit refit. Where R's behaviour is the defect, the fixture records both R's number and the corrected one.
- **A draft that a fixture contradicts.** The decision drafts are under discussion in the planning session. Do not edit the hub or the segment briefs; put the finding in your report with the evidence, and the planning session folds it in.
- **Scope of claims.** These are fixtures, not studies. Keep each to the cases that pin the behaviour, and say in the report what each one does not cover.

## File ownership

This chip writes only:
- new files `r/oracles/bass_*_oracle.R`;
- `r/scripts/bass_stage1_checks.R`;
- `r/oracles/_common.R` (this chip owns it; the other chip does not edit it);
- the fitting, outputs and calibration checklists: tick items and prepend reports under "## Reports";
- this brief's Report section.

It does not touch C#, `docs/`, the hub, the segment briefs, `r/DESCRIPTION`, `r/renv.lock`, or the numerics checklist. It does not run `dotnet` (the other chip builds in this checkout).

## Commits

Work on `main` with targeted commits of your own paths only: `git add <your new files>`, then `git commit -m "..." -- <your paths>`, so the other chip's staged work is never swept in. Never `git add -A` or `git add .`. If `.git/index.lock` exists, the other session is committing: wait and retry. Leave the user's untracked notes alone (`issues/regression/clumped-problems.md`, `issues/regression/error-capture-proposal.md`, `issues/spc/opus-SMC-notes.md`, `opus-gvdv-ideas.md`, `issues/spred/pca-degeneracy.md`). End each commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

Briefs and checklists use LF line endings; R files keep theirs. For multi-line splices into markdown, the repo's habit is a Python script with raw strings, checked with `git diff --stat`.

## Report

Append the report here when the work lands, newest first: date, commits, each fixture's evidence and limits, anything that contradicts a draft, and blockers.

### Report — 2026-10-27

`a6960f0` (oracles, `_common.R`, checks). Every item has its oracle and its checks pass under `Rscript scripts/bass_stage1_checks.R`:
- F-O5: `bass_wallenius_oracle.R`.
- F-O6: `bass_enumeration_oracle.R` (modes `swap` and `count`) and `bass_jeffreys_oracle.R`.
- O-O1 to O-O3: `bass_output_defects_oracle.R`, `bass_sensitivity_defects_oracle.R` and `bass_moments_oracle.R`.
- C-O1 and C-O2: `bass_calibration_defects_oracle.R` and `bass_covariance_oracle.R`.

The full evidence and each fixture's limits are in the [fitting](port-fitting-checklist.md), [outputs](port-outputs-checklist.md) and [calibration](port-calibration-checklist.md) checklists' reports of today.

- **How the oracles reach R.** `_common.R` gains `bass_source`, which checks that the pinned checkout is at `2671014` and sources the named files into an environment. It also gains `source_calls` and `source_assignment`, which pull an expression out of a function body so that R's code is evaluated as written, and `gauss_legendre`. Base R and `jsonlite` only.
- **No draft is contradicted.** Each R behaviour the drafts call a defect reproduces, and each corrected expectation holds. Refinements for the planning session:
  - *D01.* The inclusion–exclusion cut at K ≤ 3 is safe at R's default w₂ but not with skewed weights, which lose digits from K = 2 (1·10⁻⁷ at K = 2, 3·10⁻⁵ at K = 3 with weights spanning five decades). An error estimate would choose better than K alone. `dmwnchBass` doubles the probability at K = 1, a path R never takes.
  - *D05.* R's `lp` is the log joint at a's conditional mode, so beyond the missing term priors it carries ½ log|G| − (k/2) log(2πσ²/(1 + τ)). R's tempered kernel alone, before any swap, targets none of the candidate laws.
  - *D09.* R's Woodbury helpers match the low-rank route's accuracy, so the case against m × m inverses is cost alone. At κ = 9·10¹³ no dense reference certifies either route (the two dense algorithms disagree by 10⁻³); both stay backward stable.
  - *F-O4,* by the way: R's fixed thirds put the b = 1 chain 0.020 in total variation off the posterior on F-O6's enumerable space.
- **Blockers:** none. R's Student moments could not be compared (they need `hypergeo`); the quadrature references cover them.
