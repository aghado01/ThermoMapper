# Regression coordination: handoff (2026-11-10)

Written by the coordinating session (`local_b6732445`) at the user's quota pause. It supersedes the [2026-11-09 handoff](coordination-handoff-2026-11-09.md) for the current state. That file keeps the earlier record and the verification routine. The briefs are the record; this file points into them.

## November 11 review and finish order

The active Sol worker chats are [E21 grouped candidate densities](codex://threads/01a1cc3f-c2ad-7b82-a47a-177f777f4518) and [Conformal slice 10](codex://threads/01a1cc40-0194-7411-bf64-b2150765c9eb). Both run in the shared `main` checkout. Each maintains its brief's **In progress** entry and sends verification events to the [coordinating chat](codex://threads/01a1cc28-f4d3-7951-882a-65aeb6945f8c).

The reviewed code is `4254f5d` on `main`. Its eight changed files cover both E21 slice 5 and conformal slice 10a. The E21 collapse opt-in had already landed in `aa0d0ee`. Form B is the default on every rung. The saved Commit selection at `7b0a85d` passed **2,280/2,280**. [Opt-in receipt](../../artifacts/verify-e21-optin-7b0a85d/commit.trx)

The regression test assembly built in Release from archived sources at `4254f5d`. Of 11 selected existing cases, 10 passed. The failing test still expects a blanket refusal of grouped candidates. Six additional checks exercised the unfinished public paths against the same product code. The box-vertex check passed and five checks failed. These were focused runs. To repeat the six checks, copy the probe source into an archive's `tests/regression/` and filter the regression test project on `WipReadinessChecks`. [Existing tests](../../artifacts/verify-wip-4254f5d/baseline.trx) · [Probe source](../../artifacts/verify-wip-4254f5d/WipReadinessChecks.cs) · [Probe results](../../artifacts/verify-wip-4254f5d/probes.trx)

1. **E21 slice 5 — grouped candidate densities.** Start at `BassCandidates.LogDensity`. Its call to `PredictionCases.Prepare` omits the fit's groupings. Even with valid cell strata, a known school is rejected as an unknown predictor. The mixed-group check stops at the same point. Pass the groupings as prediction does, then check bitwise agreement with the training column and the mixed latent report. Repair cell identity as well. The crossed pairs `("a|new:b", "c")` and `("a", "b|new:c")` both received cell 11. Replace the old blanket-refusal assertion with the intended success and refusal cases. The remaining gates include the Gaussian and independent quadrature references for new groups and simulated coverage within known cells. [E21 agreement and gates](bass/ridge-blocks-brief.md#slice-5-the-candidate-seam-under-random-effects-codex-continuation-2026-11-11)
2. **Conformal 10a — the density-score margin.** With 128 draws, nine training units and an ESS floor of 10,000, every turn must be flagged. `FullStep` returned a null flagged count. It omits the configured floor when constructing the scorer and omits the flag arrays when constructing the returned margin. The grid refinement callback is declared and unused. Finish those paths, compare the proxy to exact conjugate refits, and record coverage and the BASS comparison across ESS floors. The approved sequence continues with 10b (CDF), 11 and 12b. The BASS reading informs the order of 13–18. [Conformal slice 10 plan](conformal-layer-brief.md#in-progress--slice-10-the-full-steps-error-completed-started-2026-11-10)
3. **E8 slice 7a — the projection.** The four vertices of a two-dimensional box reached every projected face within padding. Pairing the analysis with another transform that kept the predictor names but changed the first domain from `[0,1]` to `[-2,2]` succeeded. Enforce the analysis's input-scale contract, then check pushforward moments, dimension refusals and the high-dimensional domain scaling before proceeding to the projected fit (7b) and comparison (7c). [Approved projection gates](bass/active-subspace-brief.md#slice-7s-plan-rotate-and-refit-2026-11-10-for-the-users-call)

**Independent verification still due.** Conformal Lean 18 landed in `4383916` and is ready for its clean-archive checks, module build, axiom inspection and meta-ci. E34 slice 5 has the worker's **3,015/3,015** at `4135eb9`. The coordinator's verification of the slice and hint remains open. E21 slice 4 has cleared E34 slice 4's prerequisite. The u–s funnel remedy and BPPR commissioning remain user decisions. [Lean 18 report](../../lean/Protolemmata/Conformal/conformal-chip-brief.md#2026-11-10--slice-18-the-guarantee-read-out-proved) · [E34 report](bass/multi-output-brief.md#slice-5-the-pca-modes-joint-option-chip-2026-11-10)

### Independent checks on the two active slices

The four failures reproduced for these slices at `4254f5d` passed on the workers’ first corrected sources. Two additional checks passed in the same build. One exercised grouped candidate densities through the full conformal step. The other compared the returned flag counts with independent Gaussian conjugate calculations. At ESS floors 10, 100 and 400, the returned count fell below the analytic count at 0, 3 and 31 of 240 grid points respectively. The largest shortfall was three units. The Kept margin is an empirical diagnostic whose coverage reading depends on the flags containing the truly degenerate turns and the grid covering the test response. [Six checks](../../artifacts/verify-e21-conformal-nov11/review.trx) · [Source manifest](../../artifacts/verify-e21-conformal-nov11/review-inputs.json)

On a clean archive of conformal 10a commit `4395b18`, the coordinator reran both Kept-margin test classes with `BASS_QUALIFICATION=1`. All 15 passed, including the exact-refit comparison, BASS readings and the finite-endpoint refinement regression. [Coordinator receipt](../../artifacts/verify-wip-nov11-10a-focused-4395b18/qualification.trx) · [Worker report](conformal-layer-brief.md#in-progress--slice-10-the-full-steps-error-completed-started-2026-11-10)

**10a verified for its scope, 10b released.** The same archive’s Commit selection passed 2,558 of 2,559 with zero skips. Its only failure was the existing random-intercept test expecting a blanket refusal of grouped predictions. It failed in `BassFit` before reaching conformal code, and the same failure was reproduced at `4254f5d`. E21 replaced that assertion in `40158dc`. [Coordinator Commit selection](../../artifacts/verify-wip-nov11-10a-focused-4395b18/commit-selection.trx)

E21 code and tests landed in `40158dc`, with the report in `f7cee05`. All six files in the worker’s corrected qualification snapshot match their committed Git blobs. Its 25 gated cases passed. The far-tail tests exercise the Student-t peak at log mixing variance 40 with and without extra Gaussian variance, where the original scan could loop indefinitely. At nominal 0.90, marginal split coverage across the three known cells was 0.90528 for the Gaussian fit and 0.90611 for the t fit. New-group coverage remains open. [Source comparison](../../artifacts/verify-wip-nov11-e21-10a-40158dc/e21-source-check.json) · [Qualification](../../artifacts/worker-e21-s5-boundary/results/qualification.trx) · [Worker report](bass/ridge-blocks-brief.md#slice-5-the-candidate-seam-under-random-effects-codex-continuation-2026-11-11)

**E21 slice 5 and conformal 10a verified together at `40158dc`.** The clean archive's Commit selection covered the 12 changed production and test paths since the original WIP. All 2,584 tests passed with zero skips, including the replacement random-intercept assertion. The coordinator then ran all three new test classes on that archive's DLL with `BASS_QUALIFICATION=1`. All 40 passed with zero skips. The two pre-correction worker Slice runs were stopped so the final whole-assembly run can cover 10b as well. [Changed paths](../../artifacts/verify-wip-nov11-e21-10a-40158dc/selected-inputs.txt) · [Commit selection](../../artifacts/verify-wip-nov11-e21-10a-40158dc/selection.trx) · [Gated qualification](../../artifacts/verify-wip-nov11-e21-10a-40158dc/qualification.trx)

## Standing rules (the user, 2026-11-10; in [AGENTS.md](../../AGENTS.md#docops--conventions-for-issues))

- **Docs.** Chips leave `docs/` alone. Each report lists the doc changes its work calls for under "Doc changes owed". `docs/` is updated offline in batches. The user is reorganizing `docs/Maths/Regression`; the [tracking audit](../docs/regression-docs-tracking-audit.md) lists what moves out.
- **Messages.** Messages between sessions carry events and coordination, tersely, and results go in reports. Chips message each other to coordinate.
- **No transcript reading.** No session reads another's transcript. An event points the orchestrator to the report it reads.
- **Notes stay in the repository.** A chip's hand-over state goes in its brief's report as an "In progress" entry; probe code and data go to `artifacts/`. The app's side-session notes were imported into each workstream's `issues/<workstream>/handover-notes/`.
- **Carried over:**
  - Work on `main`, commit by pathspec, leave nothing staged, push nothing, and add no attribution lines.
  - Leave `basin-hold-validation.csv`, `path-profile-validation.csv` and the user's notes out of every commit.
  - Shared files (`BassConfig`, `BassModel`, `BassResolved`, `BassSweep`, `BassFit`, the moment service) are announced to the other chips before editing and committed by themselves. Where another session holds a hunk in the same file, commit your own hunk through the index.
- **Verification.** Each slice is verified before the next starts.
  - *C#:* `git archive HEAD` into `artifacts/verify-main/snap`, then the Commit selection by `scripts/test-select.ps1` called in-process with the slice's changed files.
  - *Lean:* `git archive <commit> lean` with a junction to the real `lean/.lake/packages`; then the checks script, the module build, `#print axioms` and `meta-ci -NoBuild`. **Remove the junction with `cmd /c rmdir <junction>` before deleting the snapshot**, or the real packages go with it.
  - The scripts used are recorded in the earlier handoff.

## Paused at the quota: resume these first

| Session | Work | State in the tree | Next |
|---|---|---|---|
| [E8](claude://claude.ai/epitaxy/local_66d45fa0-e3f5-4429-b588-1786ae46c0ef) ([brief](bass/active-subspace-brief.md)) | Slice 7a, the frame and the projected data. The plan was approved as recommended (`72ec20f`) | Committed as WIP at the user's word (`bd5dc9c`): `sensitivity/BassActiveSubspace.cs` and the new `spline/bass/projection/`. It compiles; it is untested | 7a, then 7b and 7c, each verified |
| [Conformal slices](claude://claude.ai/epitaxy/local_69bae1e5-39ad-45c0-8df2-e8e705b21e42) ([brief](conformal-layer-brief.md)) | Slice 10, the full step's Monte Carlo error, which had just started | Nothing in the tree | Slice 10 under the rules beside decision 9, the grid rule (`1d76a22`) included; then slices 11 and 12b (decision 10) |
| [GBASS](claude://claude.ai/epitaxy/local_add4be38-7ddb-4981-8b14-6f8ac64f7dab) ([brief](bass/ridge-blocks-brief.md)) | E21 slice 4 is done and verified (`4f5bf3f`, `94133f8`, `1940cc0`, `397ca88`, report `2cd9711`; 2,278 of 2,278 on a clean snapshot of `2cd9711`). The b = 1 collapse becomes an opt-in, at the user's word (`2d97551`) | The opt-in change, if the chip has started it | The opt-in as its own commit, verified; then slice 5, the candidate seam, agreed with the conformal chip |
| [Conformal Lean](claude://claude.ai/epitaxy/local_21649695-833f-41ac-8d7c-326b4f691dda) ([brief](../../lean/Protolemmata/Conformal/conformal-chip-brief.md)) | **Done.** Slice 18 landed in `4383916` after these rows were written, completing all six items of the commission | Nothing | Verify slice 18 (unverified). The sources it wants fetched are in "Waiting on the user" |

**E21 slice 4's finding, settled.** The first, successive-conditional Geweke design failed on ln s and is on record as unqualified.
- *The cause* is the u–s funnel: ln s's integrated autocorrelation time is about 19,000 sweeps, and z falls as the chain grows.
- *The gate* is now exact-start rows: 9 of 9 pass, and they catch a bias of the failing size at z −5.00.
- *The fix* for the funnel is the user's call (waiting list, item 1).

**E34 can start its slice 4** (noise by row and random intercepts), now that E21 slice 4 is verified. Its chip is idle.

## Landed on 2026-11-10

| Work | Commits | Verified here |
|---|---|---|
| Receipts: reproduction inputs tracked | `c4ea4fd` | — |
| Conformal Lean, slice 10 and its follow-ups (the excluded set; reading from the realized ordering) | `dfe5c99`, `a708eb2`, `f82ae97` | Yes, on clean archives of `dfe5c99` and `3b61763` |
| Conformal Lean, slice 12 (Hoff's score; the latent posterior from draws) | `f692e37` | Yes (`3b61763`) |
| Conformal Lean, slice 10's grid reading | `bb5df90` | Yes |
| Conformal Lean, slice 13, the tournament correction | `f998131` | Yes |
| Conformal Lean, slice 14, shift | `27a6662` | Yes |
| Conformal Lean, slices 15 (PCP) and 16 (batches) | `60d1c8c`, `77773c8` | Yes: checks pass, the build is clean, all 32 results on the standard axioms (one uses only two of them), and meta-ci passes |
| Conformal slice 12a, Hoff's score as an opt-in | `ab21104`, `cbd9b7d` | Yes: 2,313 of 2,313 with weighted fits, at `cbd9b7d` |
| Weighted fits through `Bass.Fit`, and their candidates | `7343bc8`, `c1bbf51`, `97b2fc4` | Yes (the same run) |
| E8 slice 6, the readout follow-ups | `276af70`, `e9d8827`, `30df59e`, `9a669b6`, `3c26660` | Yes: 2,278 of 2,278 at `3c26660` |
| E8 slice 7's plan; the BPPR draft brief | `5816337`, `72ec20f`; `57e1ef6` | Plan read and approved by the user. The brief awaits the user |
| E34 slice 5, the PCA mode's joint option; decision 2's hint | `4135eb9`, `cb4849b`, `6f2d3ea` | **Not yet.** The chip ran the whole assembly, 3,015 of 3,015 |
| E21 slice 4, random effects under the noises | `4f5bf3f`, `94133f8`, `1940cc0`, `397ca88`, `2cd9711` | Yes: 2,278 of 2,278 at `2cd9711` |
| KHAOS: SMC's step cap, the residual birth as raw-ridge default | `0ffc0e0`, `5a6cb1d`, `0738183`, `233c0e0`, `a980d1f` | Not rerun (the user: "Not now") |

## Decided on 2026-11-10

- **Conformal:**
  - decision 9, Hoff's score opt-in;
  - decision 10, the queue: 12a, then 10, 11 and 12b, and slices 13 to 18 ordered after slice 10's BASS reading.
- **E21:** decision 7 amended. Under a noise, u is a latent on every rung and is integrated out in the structure step at b = 1 only.
- **E34:** decision 2. `ScaledToResponse` stays the default; Huang–Wand is opt-in and named below N/K = 10.
- **E8:** slice 7's five decisions, as recommended. The matrix Poincaré check stays off.
- **The standing rules above.**

## Waiting on the user

1. **E21's u–s funnel fix:** commission it or not. ln s's integrated autocorrelation time is about 19,000 sweeps under every random-effects fit with a noise. The remedies are a step on s with u integrated, exact at b = 1 only, or interweaving u's centred and non-centred forms, exact at every rung (the coordinating session's recommendation).
2. **The BPPR brief** ([bppr-brief.md](bass/bppr-brief.md)): six open decisions, each with a recommendation the coordinating session agrees with.
   - The coordinating session also recommends that the moves' balance go Lean first.
   - The power spherical paper (arXiv 2006.04437) is not on disk.
   - Approval gives the go-ahead to install `BayesPPR` 0.2.0 from CRAN.
3. **Sources for the conformal Lean's slice 18.**
   - *Wanted:* Vlassis and Thomas (arXiv 2607.08415), who prove Gaffke's conjecture, and Angelopoulos et al., Conformal Risk Control (arXiv 2208.02814).
   - *Optional:* Gneiting and Raftery 2007.
   - *Its findings:* Snell and Griffiths's Theorem 4.3 is false as stated (a point-mass prior on all-zero losses breaks it), and Gressmann's proof of Proposition 2.5(i) has two sign slips while the statement stands.
4. **Slice 12a's k̂ diagnostic:** opt-in (recommended; it costs about 200× the score) or on by default.
5. **The original side-session notes** in `%APPDATA%\Claude\side-session-notes\`: delete them or keep them. Copies are in the repository; three PyTeXML notes were left there.
6. **Conformal slices 13–18's C# order,** after slice 10's BASS reading. The tournament correction (13) is the candidate if the margin still reads nothing on BASS.
7. **Lean for the active-subspace geometry** (proposed 2026-11-10):
   - E36's direction moves;
   - the matrix Poincaré check's reduction and the product law's constant;
   - P̄'s readings under mixed ranks;
   - rotate-and-refit's non-invariance.
8. **Carried from 2026-11-09:**
   - GBASS's open decisions;
   - `RungWorkers`' default and the GC pass;
   - the Lean promotion review;
   - `DenseEigen`'s relative threshold;
   - the next wave (E25 and E15, then E3);
   - test tooling's follow-ups;
   - the offline docs pass with each report's "Doc changes owed", and the ~200 links into `docs/.archive/BayesianRegression.md` once the docs layout settles.

## Local, untracked

- `artifacts/side-session-imports/` (62 MB): the probes' code, logs and data from the side-session notes.
- `artifacts/verify-*`: verification receipts, small except the TRX logs.
