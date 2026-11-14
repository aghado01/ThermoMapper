# Regression coordination: handoff (2026-11-11)

Written by the coordinating session (`local_b6732445`) at close-out, for a fresh coordinating session. It supersedes the [2026-11-10 handoff](coordination-handoff-2026-11-10.md), which keeps the November 11 Sol review and the table of what landed on 2026-11-10. The [2026-11-09 handoff](coordination-handoff-2026-11-09.md) keeps the older sessions and resources. The briefs are the record; this file points into them.

Updated 2026-11-13 with that day's round-up: E8's slice 7d verified, the LinAlg consolidation complete, and three new items for the user (the first three on the waiting list).

## Start here

1. Read [AGENTS.md](../../AGENTS.md) and [CONTRIBUTING.md](../../CONTRIBUTING.md), then this file. Read a workstream's brief before acting on it.
2. **This session's chips:** the conformal chip is building the C# for slice 14's follow-ups 1 and 2, with slices 14 and 15 verified on 2026-11-13. The conformal Lean chip's three follow-ups and its K-fold bounds are proved and verified through `7e2b905`. E21's slice 6 landed (`5cde502`) and was verified on 2026-11-13; the GBASS chip is building slice 7's item 2, item 1 verified. E8's slice 7e, the frame's column signs, landed and was verified in seven test projects. This session (`local_b6732445`) verifies those two; the second coordinating session (`local_b84fa8ad`) verifies its own two chips. E8's slice 7d landed and was verified on 2026-11-13 (below).
3. Two Opus chips were commissioned the same day by a second coordinating session, at the user's request: the defects chip and the active-subspace Lean chip ([Sessions](#sessions)). Both have landed their work. All three pieces are verified from clean archives: the Lean, the weight reader and the Cholesky refusal. Each chip has since been given more by the user (the defects chip a third slice, the Lean chip an extension), and each addition needs that session's clean-archive run when it lands.
4. Then take the [waiting list](#waiting-on-the-user) to the user.

## Standing rules

They are in [AGENTS.md](../../AGENTS.md#docops--conventions-for-issues) (the user, 2026-11-10):

- **Docs.** Chips leave `docs/` alone. Each report ends with "Doc changes owed". The user updates `docs/` offline in batches and is reorganizing `docs/Maths/Regression`. The [tracking audit](../docs/regression-docs-tracking-audit.md) lists what moves out.
- **Messages** carry events and coordination, tersely. Results go in reports. Chips message each other to coordinate.
- **No transcript reading.** No session reads another's transcript. An event points the coordinator to the report it reads.
- **Notes stay in the repository.** A chip's hand-over state is an "In progress" entry in its brief. Probe code and data go to `artifacts/`. Nothing is written to the app's side-session notes.
- **Git.**
  - Work on `main`, commit by pathspec, leave nothing staged, push nothing, add no attribution lines.
  - Keep these out of every commit:
    - `basin-hold-validation.csv`, `path-profile-validation.csv`;
    - the user's notes: `issues/regression/discussion/opus-open-items-11082026.md`, `issues/geometry/gemini-abliteralus-thermomapper.md`, `issues/regression/discussion/gemini-linalg-regression-notes.md`;
    - `issues/.obsidian/workspace.json` and `issues/Untitled.canvas`;
    - the user's reorganization in progress: `AGENTS-draft.md`, `CONTRIBUTING-draft.md`, `DOCS-OWED-draft.md`, the `* copy.md` files, the edits to `AGENTS.md` and `CONTRIBUTING.md`, and the move of the BASS reports and `opus-bayesian-evals.md` out of `issues/regression/`.
- **Shared files** (`BassConfig`, `BassModel`, `BassResolved`, `BassSweep`, `BassFit`, the moment service) are announced to the other chips before editing and committed on their own.
- **Sub-agents** need the user's permission. **Recorded decisions are settled:** execute them; don't re-ask.

## Sessions

| Chip | Brief | State at handoff |
|---|---|---|
| [E8: active subspaces](claude://claude.ai/epitaxy/local_66d45fa0-e3f5-4429-b588-1786ae46c0ef) `local_66d45fa0` | [active-subspace-brief.md](bass/active-subspace-brief.md) | Idle; slice 7e landed and verified |
| [Conformal slices](claude://claude.ai/epitaxy/local_69bae1e5-39ad-45c0-8df2-e8e705b21e42) `local_69bae1e5` | [conformal-layer-brief.md](conformal-layer-brief.md) | **Building:** slice 14's follow-ups 1 and 2 in C#; slices 14 and 15 verified |
| [Conformal Lean](claude://claude.ai/epitaxy/local_21649695-833f-41ac-8d7c-326b4f691dda) `local_21649695` | [conformal-chip-brief.md](../../lean/Protolemmata/Conformal/conformal-chip-brief.md) | Idle; slice 14's follow-ups proved and verified through `7e2b905` |
| [E34: several outputs](claude://claude.ai/epitaxy/local_26b4e7ee-1be4-482b-a19d-995477d32291) `local_26b4e7ee` | [multi-output-brief.md](bass/multi-output-brief.md) | Idle; E34 complete |
| [GBASS slices](claude://claude.ai/epitaxy/local_add4be38-7ddb-4981-8b14-6f8ac64f7dab) `local_add4be38` | [ridge-blocks-brief.md](bass/ridge-blocks-brief.md) (E21) | **Building:** E21 slice 7's item 2, the cases-per-level grid; item 1 verified |
| [KHAOS slices](claude://claude.ai/epitaxy/local_d665e4ce-cd64-427f-9c2a-8c94d96497dd), [weighted fits](claude://claude.ai/epitaxy/local_57bad71b-a245-4e83-8cb0-52c46cd48b60) | — | Idle; their work landed |
| LinAlg consolidation: the user's Sol sessions, listed in its brief | [spectral-solver-consolidation-brief.md](../linalg/spectral-solver-consolidation-brief.md) | Complete |
| [Shared-code defects](claude://claude.ai/epitaxy/local_0e31b37f-31de-42a6-a527-9c1daa3819bc) `local_0e31b37f` (Opus) | [shared-code-defects-brief.md](shared-code-defects-brief.md) | **Both slices landed:** the weight reader (`f50753f`) and the factor's refusal (`37eb63c`), each with its report. Both are verified from clean archives (the brief's two verification entries). **Slice 3,** BASS's draw weights over their whole range (`c9f9a84`), was added by the user and is verified too (2,659 of 2,659 at `badc0d9`); its R-oracle selection stands on the chip's own run. Three decisions wait on the user ([item 13](#waiting-on-the-user)) |
| [Active-subspace Lean](claude://claude.ai/epitaxy/local_eaa7e8ed-844f-440c-ae5f-b72b72d6d68a) `local_eaa7e8ed` (Opus) | [active-subspace-lean-chip-brief.md](../../lean/Protolemmata/BASS/active-subspace-lean-chip-brief.md) | **Done and verified:** all four proved (`cbfd0ec` to `85d9d53`), one hypothesis corrected after review (`1682d4a`); verified from a clean archive at `42c980d` ([item 10](#waiting-on-the-user)). **Extending,** by its own account at the user's request: items 1, 2, 3 and 5 of its report's open list (the continuous AS-2 pieces, AS-3's continuous form, AS-4(viii) in general, and trace equals rank). Those commits need a new `verify-lean` run |

Older chips are listed in the 2026-11-09 handoff. Message a chip by its `local_…` id.

## Where each workstream stands

| Workstream | Done and verified | Next |
|---|---|---|
| **E8** | Through slice 7. Slice 7 verified: 7a `7b2ca7a` (2,349 of 2,349), 7b `70a2aa2` (48 of 48), 7c `015d977` (46 of 46); 7c's report is `5064e50`. Slice 7d, the training-range rule (`566c1f1` to `6aff134`, report `f72039b`), verified 2026-11-13: 2,502 of 2,502 from a clean archive of `788d249` ([inputs](receipts/verify-e8-s7d-788d249/)). Slice 7e, the frame's column signs by `Pca.Orient` (`9b9e4e1` to `199abd7`, report `9283a91`), verified 2026-11-13 in seven test projects: regression 2,575 of 2,575 at `9283a91`, Oracle 72 of 72, and the rest ([receipts](receipts/verify-e8-s7e-9283a91/)) | Slice 8's plan (waiting list, item 15); the BPPR brief, if the user approves it. `StiefelManifold` was repaired ahead of a consumer at the user's request (`6a55e2f`, report `8d64d5b`); its verification is queued |
| **Conformal (E16)** | 12a, 10a/10b, 11 (149 of 149 at `9248e33`), 12b (2,659 of 2,659 at `6ef6c98`), 13 (1,120 of 1,120 at `60694f0`). Slice 14, shift (`c19ce98`): Commit 1,160 of 1,160 at `6200d6d` and the gated covariate-shift measurement 1 of 1 at `59dd8b4` ([receipts](receipts/verify-conformal-s14-6200d6d/)). Its ratio estimator (`4aeebfc`, `DensityRatio.Logistic`): 1,174 of 1,174 at `59dd8b4` ([receipts](receipts/verify-ratio-estimator-59dd8b4/)). Slice 15, posterior conformal prediction (`55f4b93`): Commit 1,184 of 1,184 at `3980720` ([receipts](receipts/verify-conformal-s15-21dc9f8/)) | Slice 14's follow-ups 1 and 2 in C#, building, on their Lean verified at `38f2e3f` ([receipts](receipts/verify-lean-conformal-38f2e3f/)); follow-up 3 after its Lean's run (decision 12). Slice 15's report leaves a follow-up for the user: the precision rule reads the training split's ESS, so with 1,000 training and 399 calibration points it chose m of 294 to 500 and 7% of the sets were the whole space. The brief's [coordination note](conformal-layer-brief.md#coordination-note-2026-11-13) (`4094ec7`) gives the defects chip `PredictiveChecks.cs` and `Pointwise.cs`'s weight check; a change there is agreed with it first. Slices 16 to 18 come by the user's need. Two follow-ups agreed with the GBASS chip for E21's slice 5 are recorded (`127d64d`): new-group coverage by group-level units, and the layer's check that strata lie inside cells |
| **Conformal Lean** | Slices 10, 12, 13–16, 18 and the grid reading, all verified. Slice 18 checked against its sources (`f9bcbb1`, 2026-11-13): every result stands. The β-bound's calibration is Vlassis and Thomas's Theorem 1 (P(L⁺ ≤ J) = 1 − K(X), checked exactly on 300 cases), `crc_orderings` is Angelopoulos et al.'s Theorem 1 on a grid, and the scores match Gneiting and Raftery. Three docstrings changed and no proof did. The chip verified it on a clean archive (checks 1,577 s, build, axioms, meta-ci); the coordinating session confirmed the Lean diff is docstrings only and did not rerun it. Slice 14's follow-ups 1 to 3 and the K-fold bounds (`d81f293` to `7e2b905`), verified 2026-11-13 from clean archives of `38f2e3f` and `7e2b905`, the barrel built ([receipts](receipts/verify-lean-conformal-7e2b905/)). JAWS's Step 2 argued from exchangeability and was repaired; the merged K-fold bound under likelihood ratios is recorded open, its proof failing and its search finding no counterexample | Idle |
| **E21** | Slices 0–5. The b = 1 collapse is opt-in (`BassSampling.CollapseEffectsInMoves`), with form B on every rung by default | Slice 6 landed (`5cde502`, report `6125ed3`): s's interweaving, an opt-in off by default (`BassSampling.InterweaveGroupScales`). On weak data it gives ln s 3.6 times the ESS per sweep; on the n = 300 fit it gains nothing and costs the other quantities. Verified 2026-11-13: Commit 2,536 of 2,536 at `4cba94d`, and the Geweke rows with qualification on, 13 of 13 at `1ef4846` (worst \|z\| 3.17 against the limit 4) ([receipts](receipts/verify-e21-s6-4cba94d/)). Slice 7, commissioned 2026-11-13 (`2f13572`): slice 6's two follow-ups. Item 1, verified 2026-11-13 (68 of 68 at `113532b`, the gated rows 8 of 8 at `b651bc6`; [receipts](receipts/verify-e21-s7-item1-5154bc0/)): M's loss with the step on runs through σ²'s block, and the joint (s, σ) move recovered nothing over twelve seeds, so the step is unchanged. Item 2, building: a grid of cases per level for an automatic default, whose switch goes to the user. New-group coverage stays open as the layer's follow-up; see the [two-layer addendum](bass/ridge-blocks-brief.md) (`8bc0d55`) |
| **E34** | Complete: slices 0–5. Slice 5 at 2,396 of 2,396 (`e5f142f`); slice 4 at 2,651 of 2,651 (report `fc9c343`) | Its slice 4 report's follow-ups wait for the user. Doc changes owed: `MultiOutput.md`, `BASS.md` |
| **E36, BPPR** | Draft brief `57e1ef6` | The user's approval (waiting list) |
| **KHAOS** | SMC's step cap; the residual birth as the raw-ridge default (`0ffc0e0` … `a980d1f`) | Not rerun (the user: "Not now") |
| **E3, BAMRS** | Preparation brief `f5f8d3f`; its Lean section rewritten against the revised BARS draft (`1fe6d8c`) | The user's seven decisions (waiting list). E38, a replicate layer in BAPS (`419f547`), is parked after it |
| **LinAlg consolidation** | Complete per its [brief](../linalg/spectral-solver-consolidation-brief.md), reviewed by its own coordinating chat. The element-wise relative Jacobi rule is the dense solvers' only rule (`75523c4`). The `Fma` switch and `DenseEigenOptions` are gone. `ThinSvd` carries column norms and skips the left vectors where callers read none (`9649ccf`). A leading-patterns readout by LOBPCG is opt-in (`5ddea64`). BGP factors with `DenseCholesky` (`5db2e40`). The eigensolver performance brief is closed (`9a54710`). This session verified none of it separately; 7d's run at `788d249` ran on top of it | Docs owed only (waiting list) |

**E8's training-range rule (the user, 2026-11-11; built in slice 7d).** Projected fits normalize each hinge over the training range. This settles the first item of the 2026-11-10 waiting list: column-mass support and the raw ridge under a projected fit.
- *The problem it settles.* The projected bounding box is 4.06× the data's range at p = 40. Under it, those two rules' means differ by 33% and 2.4% between the two domains. Under the training range, 0 of 60 draws differ.
- *What was built.* `BassStructureSpec.HingeNormalization` takes `Box` or `TrainingRange`. Null means the box and is not serialized, so existing configs keep their fingerprints. `BassProjectedFit.Of` resolves an unset setting to the training range and keeps a stated one. κ still moves with the box; the evaluated column does not.
- *Its follow-up,* the training range as BASS's default wherever the declared domain is wider than the data, is on the waiting list.

## Decided since the 2026-11-10 handoff

- **Conformal, decision 11** ("Yes, 13 first"). The tournament correction came first among 13–18, after 11 and 12b. It landed and is verified.
- **E8, the projected fit:** normalize each hinge over the training range (above).
- **LinAlg, the owner's five decisions** (2026-11-12, in the consolidation brief):
  - the `EIGEN_REFERENCE` build stays;
  - the `Fma` switch is removed;
  - B1 migrated its callers in stages;
  - one element-wise relative rule serves positive definite and indefinite inputs alike;
  - the matrix Poincaré check reports one ratio per draw and stays opt-in.
- **E21, the u–s funnel's fix: interweaving** (the user, 2026-11-13, in the GBASS thread; decision 7 in the [brief](bass/ridge-blocks-brief.md), `3d46d04`). It is exact at every rung. The 19,000-sweep autocorrelation time came from Geweke's toy chain (12 cases in 3 levels, y redrawn each step); on the n = 300 fit ln s's is about 5, so the funnel bites where the data pin each group's effect weakly.
- **Conformal: slice 14, then 15** (the user, 2026-11-13, in the conformal chip's thread: "ok 14 then 15"). One chip builds both, and 14 is verified before 15 starts. A commission of 15 to a new chip (`0a1c648`) was withdrawn in `4094ec7`; that chip was stopped before it wrote anything.
- **The Lean chip's two findings against slice 18's papers** (the user, 2026-11-13). Neither touches built code or docs.
  - *Snell and Griffiths's Theorem 4.3 fails as stated,* and Corollary 4.4 with it. A point-mass prior consistent with the losses defeats the claim that L⁺ bounds the posterior under every prior. What holds is frequentist: b*_β is Gaffke's bound, calibrated over replicated calibration sets (Vlassis and Thomas's Theorem 1). Conformal C# slice 18's plan row now reads it out as that confidence bound, never as a posterior.
  - *Gressmann et al.'s proof of Proposition 2.5(i) has two sign slips;* the statement stands and the Lean proves it independently. Recorded in the Lean ledger; nothing to do.
  - Both belong in the corrections of a conformal page under `docs/Maths` when one is written.
- **E8's hinge normalization: no change** (the user, 2026-11-13, in E8's thread; recorded in `fd68f68`). `Box` stays the default for ordinary fits and `TrainingRange` for projected hinge fits, with an explicit setting winning in either case. There is no width threshold, since width cannot tell a prediction allowance from a choice of where the raw ridge's prior is calibrated. For an ordinary fit whose domain is widened only to predict beyond the data, state `TrainingRange` explicitly. Doc owed: `BASS.md`'s basis section.
- **E38 recorded** (`419f547`): a replicate layer in BAPS, subsuming the hierarchical P-spline, parked after E3.
- **The requested papers are on disk** (the user, 2026-11-11):

| Paper | Location |
|---|---|
| The power spherical law (De Cao and Aziz, arXiv 2006.04437), for BPPR | `D:\aghado01\codex-scientiae\supellex-private\gauntlet\bayesian\2006.04437v2\` (PDF, TeX, HTML) |
| Vlassis and Thomas (arXiv 2607.08415), Gaffke's conjecture | `…\gauntlet\conformal-bayes\2607.08415v1\` |
| Angelopoulos et al., Conformal Risk Control (arXiv 2208.02814) | `…\gauntlet\conformal-bayes\2208.02814v4\` |
| Gneiting and Raftery 2007, JASA | `…\supellex-private\PDFs\Bayesian\Gneiting2007jasa.pdf` |

## Waiting on the user

New on 2026-11-13:

1. **BAMRS (E3), free-knot tensor-product surfaces.** Its working directory is [abmrs-bamrs/](abmrs-bamrs/) (the user, 2026-11-13): `briefs-slices/`, `papers/`, `probes/`, `receipts/`, `verification/` and the authors' supplement. The [preparation brief](abmrs-bamrs/briefs-slices/bamrs-brief.md) has seven decisions, each with a recommendation:
   - one engine over the product carrier, with BARS as its one-axis case and BARS's pinned tier as the acceptance (the alternatives are a sibling ensemble or a BASS term type);
   - the g-prior with the thin-column gate by default; a ridge prior only if slice 2 measures the gate refusing configurations the posterior would hold;
   - the count law's default left open until the benchmark;
   - fixed per-axis counts as a mode;
   - readouts in three tiers;
   - the full-grid shortcut as its own slice;
   - non-normal responses after the Normal model qualifies.

   Approval commissions slice 0, which only reads and measures: `BarsEnsemble`'s body, and a bicubic design through `WeightedNormalModel` against dense algebra.
2. **E8's follow-up, the training range as BASS's default:** decided on 2026-11-13, no change (under Decided).
    - *Lean (2026-11-13):* the invariance this rests on is proved. Under the rule the evaluated hinges do not read the declared box (`trainingRange_hinge`); κ alone does (`kappa_moves`); the g-form reads no column scale under either rule (`hat_invariant`).
3. **PBX 9501, the application in 2311.18146:** deferred by the user on 2026-11-13 to a broader pass over the applications in the papers behind the regression engine (E8's brief, `f978a86`). The data stay in `artifacts/datasets/PBX9501`; the authors' script's mixed draws (numerator from draws 1 to 10, denominator from 100 to 1000) are recorded there.

Carried:

4. **The BPPR brief** ([bppr-brief.md](bass/bppr-brief.md)) has six open decisions. The coordinating session agrees with each of their recommendations, and also recommends that the moves' balance go Lean first. The power spherical paper is now on disk. Approval also gives the go-ahead to install `BayesPPR` 0.2.0 from CRAN.
5. **Conformal Lean, slice 18's sources:** done on 2026-11-13 at the user's request (`f9bcbb1`; see the Conformal Lean row). The chip's two earlier findings are settled (under Decided).
6. **The next conformal C# slice:** decided on 2026-11-13, 14 then 15 (under Decided). The number is kept so other notes' item numbers still hold.
7. **Slice 12a's k̂ diagnostic:** opt-in (recommended; it costs about 200× the score) or on by default.
8. **E34 slice 4's follow-ups**, listed in its report (`fc9c343`).
9. **The original side-session notes** in `%APPDATA%\Claude\side-session-notes\`: delete them or keep them. Copies are in each workstream's `issues/<workstream>/handover-notes/`, except three PyTeXML notes, which were left there.
10. **Lean for the active-subspace geometry** (proposed 2026-11-10; commissioned by the user on 2026-11-13). All four targets are proved in `lean/Enthymemata/ActiveSubspace/` and verified from a clean archive at `42c980d`. The chip's [brief](../../lean/Protolemmata/BASS/active-subspace-lean-chip-brief.md) holds the reports. What comes back:
    - *for item 2:* the training-range rule's invariance is proved;
    - *for E8's brief and BPPR's brief:* the report gives each a line to take (slice 5's readouts and slice 7; decision 5 and slice 2's ratios);
    - *BPPR's ratios:* the paper's printed birth ratio slips at one active coordinate; the change move never flips a sign there, so the sampler treats that case apart; a further 1/(M + 1) in a birth targets λ^M/(M!)², not the Poisson prior;
    - *two attributions in E8's brief change, the readouts standing:* the cross-rank distance is not Ye and Lim's metric, and Barbieri and Berger prove the ½ rule for another loss;
    - *the frame:* the refit depends on the frame at every frame, since only signed permutations keep the hinge class. The eigengap report warns of ties only;
    - *still open, as posed:* the one-variable Poincaré inequality is a hypothesis (Mathlib has none), and the direction moves are proved in finite algebraic form.
11. **Carried from 2026-11-09:**
    - GBASS's open decisions;
    - `RungWorkers`' default and the GC pass;
    - the Lean promotion review;
    - the next wave: E25 and E15 (E3 now has its brief, item 1);
    - test tooling's follow-ups.
12. **The offline docs pass.** Each report's "Doc changes owed" lists its sections. The owed pages are piling up:
    - `LinearAlgebra.md` and `ActiveSubspaces.md`, from the LinAlg consolidation;
    - `BASS.md` and `ActiveSubspaces.md`, from 7d;
    - `MultiOutput.md` and `BASS.md`, from E34.

    The ~200 links into `docs/.archive/BayesianRegression.md` wait until the docs layout settles.
13. **The defects chip's four decisions** (2026-11-13), with measurements in its [reports](shared-code-defects-brief.md#reports):
    - *the floor's rule:* the absolute 10⁻¹² floor fires on scale alone, and once it fires `LogDet` and the solves are wrong while every entry stays finite. The chip recommends no floor, with callers regularizing; a floor relative to the largest diagonal is the alternative. GMM leans on the floor where its ridge is 0 or a warm start is singular;
    - *a pivoted factor in `Maths.LinAlg`:* promote `BgpDrawnFunction`'s private factor, with a low-rank solve and sample;
    - *retiring `CholeskyDecomposition` into `DenseCholesky`:* a removal, so its own question; it follows the floor's decision;
    - *BASS's draw-weight readers* (`Draws3.Select`, `PosteriorSummary.Over`): **decided by the user,** as the chip's brief records it: the full-range rule, built as the chip's slice 3.
    - *A review of the other three* is in `archive/bass-gbass-khaos-activespaces/discussion/astra-shared-code-cleanup.md`. It recommends dropping the floor, promoting the pivoted factor with a narrow first contract, and retiring the class after its callers migrate. The chip corrected slice 2's report from it (`77f80ba`).
14. **E8: a sign convention for the projection frame** (raised by the active-subspace Lean chip, 2026-11-13). **Decided by the user the same day,** as the first coordinating session reports it ("Yes, commission it"): E8 slice 7e, commissioned to the E8 chip in `08ed07b` ([brief](bass/active-subspace-brief.md)) and verified by that session. The record of the finding:
    - *The finding, checked in source:* `BassProjection.Of` copies `SubspacePosterior.MeanCVectors` as `DenseEigen` returns them (`BassProjection.cs:136`), and neither fixes a column's sign. A projected refit is invariant in law under a sign flip of a frame column, but not draw for draw. The invariance in law rests on an informal argument: AS-4(iii)'s class-level converse is not formalized. A solver change that flips a column (`75523c4` is the kind) moves the refit's draws and any pin on a projected fit.
    - *The proposal:* make each frame column's largest-magnitude entry positive, ties to the lowest index. Then "W₁ fixes the refit" holds bit for bit across solver versions.
    - *The rule already exists* (found by the first coordinating session, checked in source): `Pca.Orient` (`src/maths/linalg/PCA.cs:52`) makes an axis's largest-magnitude entry positive, the first on ties, and its doc names it as shared with estimators that return principal axes of their own. `MoMPCA` calls it, and `BassDecomposition.cs:418` holds an inline copy.
    - *Placement:* `BassProjectionComparison.cs:97` compares the frame with `MeanCVectors` by exact inequality. Applying `Pca.Orient` where `SubspacePosterior` forms `MeanCVectors` gives E8 the repo's convention, and the frame and that comparison agree with no further change.
    - *Caveat:* where two entries of opposite sign nearly tie in magnitude, the sign stays sensitive to rounding under any rule of this kind.
15. **E8 slice 8's plan** (drafted at the user's request, 2026-11-13, `a6da53b`; [brief](bass/active-subspace-brief.md)): the five unbuilt readout and comparison items as slices 8a (co-active functionals at new locations; Poincaré constants for mixture and Student t laws), 8b (a held-out log predictive density; cross-fitting with groupings and case weights) and 8c (basis fits through a projection). Nine decisions, A1 to E3, each with a recommendation.
    - *What moves:* pins on projected fits, once, wherever the solver's sign differs from the convention.

Closed since the 2026-11-11 handoff: `DenseEigen`'s relative threshold (`75523c4`), the projected fit's prior rules (7d), and regression's use of the legacy `CholeskyDecomposition` (BGP moved, `5db2e40`).

## Open defects in shared code

Both were commissioned by the user on 2026-11-13, to the defects chip ([brief](shared-code-defects-brief.md)): the weight reader first, then the factor's refusal of a non-finite result. Both repairs have landed: the reader keeps every finite weight's normalized log (`f50753f`), and the factor refuses a result that is not finite, with its three callers reading the refusal (`37eb63c`). Both repairs are verified from clean archives. A floored factor that stays finite is still read, and still wrong once a pivot floors: that is the first decision of [item 13](#waiting-on-the-user).

- **`CholeskyDecomposition` on a near-singular covariance.** Found by conformal slice 11 (`9248e33`).
  - The factor's entries go infinite while `LogDet` stays finite, and only `Status.Finite` reports it.
  - A caller that does not read the status gets wrong numbers silently.
  - Regression no longer calls it: BGP moved to `DenseCholesky` (`5db2e40`). Its remaining callers are GMM (`GaussianComponent`, `BhattacharyyaCoefficient`) and viz-core's serializer.
  - This is for LinAlg's owner.
- **The shared weight reader at extreme weights** (`PredictiveChecks.Rows`). Weights like `[1e300, 1e-300]` and `[1e308, 1e308]` read as −∞ log weights. The acceptance criteria are in conformal slice 10's report.

## Verification

Each slice is verified from a clean archive before the next one starts on the same code. The two scripts used all week are in [scripts/](../../scripts/) (moved there by the user on 2026-11-13). Both were dry-run tested at close-out (`-NoRun`): the snapshot matched the tracked files, and the junction came and went with the real packages intact. The C# script's first real run was 7d's (2,502 tests; the selection took 890 s). Since `959ceae` it checks its arguments before writing anything, and reads a relative `-ChangedFile` from the caller's directory.

- **C#:** [verify-csharp.ps1](../../scripts/verify-csharp.ps1).
  - It archives HEAD into `artifacts/verify-main/snap`, leaving out `r/`, then runs `scripts/test-select.ps1 -Scope Commit` over the changed files listed in `-ChangedFile`. That list is `git diff --name-only <last verified>..HEAD -- src tests projects`.
  - It copies the TRX, its counters and the log to `artifacts/verify-<Name>/`, then removes `verify-main`.
  - It runs below normal priority on two processors.
  - **The selection runs in one test project,** the regression project unless `-Project` names another (the second coordinating session found this on 2026-11-13). Run it once for each test project the changed sources reach: a change to `src/maths/linalg/PCA.cs`, for one, selects `MoMPCATests` and `SubspaceAnnealerTests` in `Maths.Geometry.Tests`.
  - Gated tests behind `BASS_QUALIFICATION` are run separately, from a fresh snapshot with a `dotnet test --filter` over the named tests, since `test-select` takes no filter.
- **Lean:** [verify-lean.ps1](../../scripts/verify-lean.ps1).
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
