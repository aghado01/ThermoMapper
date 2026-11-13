# Regression coordination: handoff (2026-11-11)

Written by the coordinating session (`local_b6732445`) at close-out, for a fresh coordinating session. It supersedes the [2026-11-10 handoff](coordination-handoff-2026-11-10.md), which keeps the November 11 Sol review and the table of what landed on 2026-11-10. The [2026-11-09 handoff](coordination-handoff-2026-11-09.md) keeps the older sessions and resources. The briefs are the record; this file points into them.

## Start here

1. Read [AGENTS.md](../../AGENTS.md) and [CONTRIBUTING.md](../../CONTRIBUTING.md), then this file. Read a workstream's brief before acting on it.
2. **No chip is building.** E8's slice 7d landed and was verified on 2026-11-13 (below).
3. Then take the [waiting list](#waiting-on-the-user) to the user. No other chip has work commissioned.

## Standing rules

They are in [AGENTS.md](../../AGENTS.md#docops--conventions-for-issues) (the user, 2026-11-10):

- **Docs.** Chips leave `docs/` alone. Each report ends with "Doc changes owed". The user updates `docs/` offline in batches and is reorganizing `docs/Maths/Regression`. The [tracking audit](../docs/regression-docs-tracking-audit.md) lists what moves out.
- **Messages** carry events and coordination, tersely. Results go in reports. Chips message each other to coordinate.
- **No transcript reading.** No session reads another's transcript. An event points the coordinator to the report it reads.
- **Notes stay in the repository.** A chip's hand-over state is an "In progress" entry in its brief. Probe code and data go to `artifacts/`. Nothing is written to the app's side-session notes.
- **Git.**
  - Work on `main`, commit by pathspec, leave nothing staged, push nothing, add no attribution lines.
  - Keep these out of every commit: `basin-hold-validation.csv`, `path-profile-validation.csv`, the user's notes (`issues/regression/discussion/opus-open-items-11082026.md`, `issues/geometry/gemini-abliteralus-thermomapper.md`), `issues/.obsidian/workspace.json` and `issues/Untitled.canvas`.
- **Shared files** (`BassConfig`, `BassModel`, `BassResolved`, `BassSweep`, `BassFit`, the moment service) are announced to the other chips before editing and committed on their own.
- **Sub-agents** need the user's permission. **Recorded decisions are settled:** execute them; don't re-ask.

## Sessions

| Chip | Brief | State at handoff |
|---|---|---|
| [E8: active subspaces](claude://claude.ai/epitaxy/local_66d45fa0-e3f5-4429-b588-1786ae46c0ef) `local_66d45fa0` | [active-subspace-brief.md](bass/active-subspace-brief.md) | **Running:** the projected fit under the training-range rule (below) |
| [Conformal slices](claude://claude.ai/epitaxy/local_69bae1e5-39ad-45c0-8df2-e8e705b21e42) `local_69bae1e5` | [conformal-layer-brief.md](conformal-layer-brief.md) | Idle; slices through 13 done |
| [Conformal Lean](claude://claude.ai/epitaxy/local_21649695-833f-41ac-8d7c-326b4f691dda) `local_21649695` | [conformal-chip-brief.md](../../lean/Protolemmata/Conformal/conformal-chip-brief.md) | Idle; commission complete |
| [E34: several outputs](claude://claude.ai/epitaxy/local_26b4e7ee-1be4-482b-a19d-995477d32291) `local_26b4e7ee` | [multi-output-brief.md](bass/multi-output-brief.md) | Idle; E34 complete |
| [GBASS slices](claude://claude.ai/epitaxy/local_add4be38-7ddb-4981-8b14-6f8ac64f7dab) `local_add4be38` | [ridge-blocks-brief.md](bass/ridge-blocks-brief.md) (E21) | Idle; E21 slices 0–5 done |
| [KHAOS slices](claude://claude.ai/epitaxy/local_d665e4ce-cd64-427f-9c2a-8c94d96497dd), [weighted fits](claude://claude.ai/epitaxy/local_57bad71b-a245-4e83-8cb0-52c46cd48b60) | — | Idle; their work landed |

Older chips are listed in the 2026-11-09 handoff. Message a chip by its `local_…` id.

## Where each workstream stands

| Workstream | Done and verified | Next |
|---|---|---|
| **E8** | Through slice 7. Slice 7 verified: 7a `7b2ca7a` (2,349 of 2,349), 7b `70a2aa2` (48 of 48), 7c `015d977` (46 of 46); 7c's report is `5064e50`. Slice 7d, the training-range rule (`566c1f1` to `6aff134`, report `f72039b`), verified 2026-11-13: 2,502 of 2,502 from a clean archive of `788d249` ([inputs](receipts/verify-e8-s7d-788d249/)) | The BPPR brief, if the user approves it |
| **Conformal (E16)** | 12a, 10a/10b, 11 (149 of 149 at `9248e33`), 12b (2,659 of 2,659 at `6ef6c98`), 13 (1,120 of 1,120 at `60694f0`) | The queue is empty. Slices 14 to 18 come by the user's need (waiting list) |
| **Conformal Lean** | Slices 10, 12, 13–16, 18 and the grid reading, all verified | Slice 18's sources are now on disk (waiting list) |
| **E21** | Slices 0–5. The b = 1 collapse is opt-in (`BassSampling.CollapseEffectsInMoves`), with form B on every rung by default | The u–s funnel fix (waiting list). New-group coverage stays open as the layer's follow-up; see the [two-layer addendum](bass/ridge-blocks-brief.md) (`8bc0d55`) |
| **E34** | Complete: slices 0–5. Slice 5 at 2,396 of 2,396 (`e5f142f`); slice 4 at 2,651 of 2,651 (report `fc9c343`) | Its slice 4 report's follow-ups wait for the user. Doc changes owed: `MultiOutput.md`, `BASS.md` |
| **E36, BPPR** | Draft brief `57e1ef6` | The user's approval (waiting list) |
| **KHAOS** | SMC's step cap; the residual birth as the raw-ridge default (`0ffc0e0` … `a980d1f`) | Not rerun (the user: "Not now") |

**E8's training-range rule (the user, 2026-11-11, in E8's thread).** Projected fits normalize each hinge over the training range. This settles the first item of the 2026-11-10 waiting list: column-mass support and the raw ridge under a projected fit.
- *The problem it settles.* The projected bounding box is 4.06× the data's range at p = 40. Under it, those two rules' means differ by 33% and 2.4% between the two domains.
- *Where it is recorded.* The E8 chip records the decision and its gates in its brief with the report.
- *The working tree at handoff.* E8 is the only chip running. The tree holds uncommitted edits to 14 BASS files (the design and kernel files, `BassProjectedFit`, and the sensitivity files) plus the user's two Obsidian files.

## Decided since the 2026-11-10 handoff

- **Conformal, decision 11** ("Yes, 13 first"). The tournament correction came first among 13–18, after 11 and 12b. It landed and is verified.
- **E8, the projected fit:** normalize each hinge over the training range (above).
- **The requested papers are on disk** (the user, 2026-11-11):

| Paper | Location |
|---|---|
| The power spherical law (De Cao and Aziz, arXiv 2006.04437), for BPPR | `D:\aghado01\codex-scientiae\supellex-private\gauntlet\bayesian\2006.04437v2\` (PDF, TeX, HTML) |
| Vlassis and Thomas (arXiv 2607.08415), Gaffke's conjecture | `…\gauntlet\conformal-bayes\2607.08415v1\` |
| Angelopoulos et al., Conformal Risk Control (arXiv 2208.02814) | `…\gauntlet\conformal-bayes\2208.02814v4\` |
| Gneiting and Raftery 2007, JASA | `…\supellex-private\PDFs\Bayesian\Gneiting2007jasa.pdf` |

## Waiting on the user

1. **E21's u–s funnel fix.** Commission it or not. ln s's integrated autocorrelation time is about 19,000 sweeps under every random-effects fit with a noise. There are two remedies:
   - a step on s with u integrated out, exact at b = 1 only;
   - interweaving u's centred and non-centred forms, exact at every rung. This is the coordinating session's recommendation.
2. **The BPPR brief** ([bppr-brief.md](bass/bppr-brief.md)) has six open decisions. The coordinating session agrees with each of their recommendations, and also recommends that the moves' balance go Lean first. The power spherical paper is now on disk. Approval also gives the go-ahead to install `BayesPPR` 0.2.0 from CRAN.
3. **Conformal Lean, slice 18's sources.** Should the Lean chip read them against its slice 18 report?
   - *Gaffke.* It proved one unit of Gaffke's bound (`gaffke_one`); Vlassis and Thomas prove the whole conjecture.
   - *Conformal risk control.* Its `crc_orderings` reads Snell and Griffiths's Proposition 3.2 as conformal risk control.
   - *Proper scores.* Gneiting and Raftery underlie that part.
   - *Its findings to confirm.* Snell and Griffiths's Theorem 4.3 is false as stated. Gressmann's proof of Proposition 2.5(i) has two sign slips, while the statement stands.
4. **The next conformal C# slice among 14–18.** Slice 15 (posterior conformal prediction, for conditional coverage) is recommended. Its Lean is verified, as is the Lean for 14, 16 and 18.
5. **Slice 12a's k̂ diagnostic:** opt-in (recommended; it costs about 200× the score) or on by default.
6. **E34 slice 4's follow-ups**, listed in its report (`fc9c343`).
7. **The original side-session notes** in `%APPDATA%\Claude\side-session-notes\`: delete them or keep them. Copies are in each workstream's `issues/<workstream>/handover-notes/`, except three PyTeXML notes, which were left there.
8. **Lean for the active-subspace geometry** (proposed 2026-11-10):
   - E36's direction moves;
   - the matrix Poincaré check's reduction and the product law's constant;
   - P̄'s readings under mixed ranks;
   - rotate-and-refit's non-invariance.
9. **Carried from 2026-11-09:**
   - GBASS's open decisions;
   - `RungWorkers`' default and the GC pass;
   - the Lean promotion review;
   - `DenseEigen`'s relative threshold;
   - the next wave (E25 and E15, then E3);
   - test tooling's follow-ups;
   - the offline docs pass with each report's "Doc changes owed", and the ~200 links into `docs/.archive/BayesianRegression.md` once the docs layout settles.

## Open defects in shared code

- **`CholeskyDecomposition` on a near-singular covariance.** Found by conformal slice 11 (`9248e33`).
  - The factor's entries go infinite while `LogDet` stays finite, and only `Status.Finite` reports it.
  - A caller that does not read the status gets wrong numbers silently. BGP's drawn function uses a pivoted factor instead.
  - This is for LinAlg's owner.
- **The shared weight reader at extreme weights** (`PredictiveChecks.Rows`). Weights like `[1e300, 1e-300]` and `[1e308, 1e308]` read as −∞ log weights. The acceptance criteria are in conformal slice 10's report.

## Verification

Each slice is verified from a clean archive before the next one starts on the same code. The two scripts used all week are tracked in [receipts/verification/](receipts/verification/). Both were dry-run tested at close-out (`-NoRun`): the snapshot matched the tracked files, and the junction came and went with the real packages intact.

- **C#:** [verify-csharp.ps1](receipts/verification/verify-csharp.ps1).
  - It archives HEAD into `artifacts/verify-main/snap`, leaving out `r/`, then runs `scripts/test-select.ps1 -Scope Commit` over the changed files listed in `-ChangedFile`. That list is `git diff --name-only <last verified>..HEAD -- src tests projects`.
  - It copies the TRX, its counters and the log to `artifacts/verify-<Name>/`, then removes `verify-main`.
  - It runs below normal priority on two processors.
- **Lean:** [verify-lean.ps1](receipts/verification/verify-lean.ps1).
  - It archives `<Commit>`'s `lean/` into `artifacts/verify-lean-<Commit>/` and links `.lake/packages` to the real packages with a junction.
  - It then runs the checks script, `lake build <Module>`, `#print axioms` through `-AxiomsFile` (`import <Module>`, then one `#print axioms` line per new result) and `meta-ci -NoBuild`.
  - **It removes the junction with `cmd /c rmdir` before deleting anything**, and checks that the real packages (9 entries) survived.
- A slice counts as verified when every test in its selection passes with zero skips. Gated qualification classes run with `BASS_QUALIFICATION=1` where a report names them.

## Resources

- **Oracle pins, read-only clones, and the BPPR and E34 papers:** listed under Resources in the [2026-11-09 handoff](coordination-handoff-2026-11-09.md#resources).
- **The paper library:**
  - `D:\aghado01\codex-scientiae\supellex-private\gauntlet\{bayesian,conformal-bayes,stats}\`, each with an `inventory.jsonl`;
  - `…\supellex-private\PDFs\Bayesian\`, also with an `inventory.jsonl`.
- **Tracked reproduction inputs:** [receipts/](receipts/README.md).

## Local, untracked

- `artifacts/side-session-imports/` (62 MB): the probes' code, logs and data from the side-session notes.
- `artifacts/verify-*`: verification receipts. They are small except for the TRX logs.
