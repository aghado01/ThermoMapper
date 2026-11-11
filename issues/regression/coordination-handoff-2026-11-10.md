# Regression coordination: handoff (2026-11-10)

Written by the coordinating session (`local_b6732445`) at the user's quota pause. It supersedes the [2026-11-09 handoff](coordination-handoff-2026-11-09.md) for the current state. That file keeps the earlier record and the verification routine. The briefs are the record; this file points into them.

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
| [E8](claude://claude.ai/epitaxy/local_66d45fa0-e3f5-4429-b588-1786ae46c0ef) ([brief](bass/active-subspace-brief.md)) | Slice 7a, the frame and the projected data. The plan was approved as recommended (`72ec20f`) | Uncommitted: `sensitivity/BassActiveSubspace.cs` and the new `spline/bass/projection/` | 7a, then 7b and 7c, each verified |
| [Conformal slices](claude://claude.ai/epitaxy/local_69bae1e5-39ad-45c0-8df2-e8e705b21e42) ([brief](conformal-layer-brief.md)) | Slice 10, the full step's Monte Carlo error, which had just started | Nothing in the tree | Slice 10 under the rules beside decision 9, the grid rule (`1d76a22`) included; then slices 11 and 12b (decision 10) |
| [GBASS](claude://claude.ai/epitaxy/local_add4be38-7ddb-4981-8b14-6f8ac64f7dab) ([brief](bass/ridge-blocks-brief.md)) | E21 slice 4, random effects under the noises. The code is committed (`4f5bf3f`, ledger `94133f8`) | Uncommitted: the brief's "In progress" entry and `RandomEffectsNoiseQualificationTests.cs`, whose Geweke rows were rewritten as exact-start tests | Finish the qualification on the four conditions below, report, be verified; then slice 5, the candidate seam, agreed with the conformal chip |
| [Conformal Lean](claude://claude.ai/epitaxy/local_21649695-833f-41ac-8d7c-326b4f691dda) ([brief](../../lean/Protolemmata/Conformal/conformal-chip-brief.md)) | **Done.** Slice 18 landed in `4383916` after these rows were written, completing all six items of the commission | Nothing | Verify slice 18 (unverified). The sources it wants fetched are in "Waiting on the user" |

**E21 slice 4's open finding.** The first qualification run failed Geweke on ln s in nearly every row, with z up to 5.1. The chip reads it as slow mixing of the half-Cauchy scale, not a kernel defect. Verification needs four things:
1. the failing successive-chain rows kept on record as unqualified;
2. exact-start rows covering the same functions, ln s included, with enough replicates to see a bias of that size;
3. evidence of slow mixing: ln s's autocorrelation, and z shrinking with the chain's length;
4. if the cause is the u–s funnel under form B, the remedy goes to the user as a design change (a step on s with u integrated at b = 1, or interweaving).

**E34 waits.** Its slice 4, noise by row and random intercepts, starts after E21 slice 4 is verified.

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
| E21 slice 4's code | `4f5bf3f`, `94133f8` | **Not yet:** its qualification is unfinished |
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

1. **The BPPR brief** ([bppr-brief.md](bass/bppr-brief.md)): six open decisions, each with a recommendation the coordinating session agrees with.
   - The coordinating session also recommends that the moves' balance go Lean first.
   - The power spherical paper (arXiv 2006.04437) is not on disk.
   - Approval gives the go-ahead to install `BayesPPR` 0.2.0 from CRAN.
2. **Sources for the conformal Lean's slice 18.**
   - *Wanted:* Vlassis and Thomas (arXiv 2607.08415), who prove Gaffke's conjecture, and Angelopoulos et al., Conformal Risk Control (arXiv 2208.02814).
   - *Optional:* Gneiting and Raftery 2007.
   - *Its findings:* Snell and Griffiths's Theorem 4.3 is false as stated (a point-mass prior on all-zero losses breaks it), and Gressmann's proof of Proposition 2.5(i) has two sign slips while the statement stands.
3. **Slice 12a's k̂ diagnostic:** opt-in (recommended; it costs about 200× the score) or on by default.
4. **The original side-session notes** in `%APPDATA%\Claude\side-session-notes\`: delete them or keep them. Copies are in the repository; three PyTeXML notes were left there.
5. **Conformal slices 13–18's C# order,** after slice 10's BASS reading. The tournament correction (13) is the candidate if the margin still reads nothing on BASS.
6. **Lean for the active-subspace geometry** (proposed 2026-11-10):
   - E36's direction moves;
   - the matrix Poincaré check's reduction and the product law's constant;
   - P̄'s readings under mixed ranks;
   - rotate-and-refit's non-invariance.
7. **Carried from 2026-11-09:**
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
