# Regression docs: development tracking to move out

Date: 2026-11-10. A read-only audit, at the user's request ("Yes, flag them"), of `docs/Maths/Regression/*.md` during the user's reorganization of those pages. It flags passages that track development rather than state the method. For each passage it notes whether the content already sits in an `issues/` brief or report. Nothing in `docs/` was edited. The pages were read as they stood in the working tree, the user's uncommitted edits included. Line numbers are as read and will drift, so each flag is anchored by its heading.

## The rule applied

A page states the settled method: the model, the math, the contract, the options, usage, and the decided trade-offs with their alternatives. Development tracking belongs in the brief's report under `issues/`:
- verification results and test classes;
- measured costs;
- build dates and slice numbers;
- commit hashes;
- what is not built yet;
- open decisions;
- exploratory readings.

## How it was checked

- *Markers, per section:* slice references, dates, commit hashes, test counts, timings, words like checked, measured, built and waits, and links into `issues/`.
- *Overlap:* each section's distinctive numbers (five characters or more, round numbers excluded) were searched in `issues/**/*.md` and `lean/Protolemmata/**/*.md`. The file holding the most of them was taken as the likely home.
- *What the overlap shows:* that the numbers exist in that file, not that the prose matches. Where the best match is `issues/spc-BARS/stopping-theory-primer.md`, a long unrelated file, read it as chance.
- *Reading:* I read the flagged lines and sampled the sections, not every line of 5,212.

## The patterns

| | Pattern | Example | Disposition |
|---|---|---|---|
| A | **Verification sections**: "What was checked", "Checks", "Checks and cost", "Checked", qualification tables | ActiveSubspaces' five "What was checked (slice N)" tables | Delete where the report holds them, about 400 lines across the pages. A one-line pointer to the report can stay if wanted |
| B | **Measurements inside method sections**: timings, readings over seeds, benchmark tables | VarianceMeanMixtures' GIG draw "took 0.08 to 0.36 µs a draw"; BARS's `MixingBenchmark` and `BirthBenchmark` paragraphs | Keep the conclusion as guidance ("judge a run by ESS"); move the numbers to the report |
| C | **Build provenance and status**: section leads such as "Built 2026-11-05 ([the build brief]…), slice 1"; headings with "(slice N)"; "Not built yet"; "Refused until their slices"; status columns in tables | Every VarianceMeanMixtures family section; BASS's "Stage N (date)" leads; the extension seams' "Today" and "Filled by" columns | Strip the leads, rename headings by topic, and keep "refuses X" only as a contract statement. Build status goes to the checklists |
| D | **Readings and open design calls**: the investigation behind a decision, or a question still open | ConformalPrediction's "The predictive distribution (slice 4): the reading"; BASS's "this is a design call for the user's review" in "a design whose working response moves" | Move to the brief as a design note; the page keeps the decided outcome and its trade-offs |
| E | **Change history and provenance tags**: "Until 2026-10-27 it did X", "(found by the SPC × BARS work, 2026-10-24)", "the user's decision of 2026-11-07", inline commit hashes | `5ce7453` in BASS's sweep; Regression's "joined the shell on 2026-10-23" | Drop the history (git has it), or move a real change to Sources&Changes. Keep credits (see the gray zones) |

**Stale already.** These show why status goes stale on a method page:
- PolynomialChaos says KHAOS's SMC trap remedy "is the user's open decision". It was decided on 2026-11-09: the step cap, opt-in, and the residual birth as the raw-ridge default.
- ConformalPrediction's "Not built yet" lists CB-4's latent-posterior reading, Hoff's score and BGP's drawn f. Slice 12's Lean is done and slice 11 is queued.
- Status&Caveats opens with "This section describes BASS through stage 4", written for the old single page.

## By page

The share given is the section's distinctive numbers found in its best-matching `issues/` file. "Any" is the share found anywhere in `issues/`.

**[ActiveSubspaces.md](../../docs/Maths/Regression/ActiveSubspaces.md)**
- A: "What was checked" for slices 1 to 5, about 105 lines. All in [active-subspace-brief.md](../regression/bass/active-subspace-brief.md): 3/3, 5/5, 35/35, 16/16 and 11/11.
- C: section headings carry "(slice N)". "The decisions (the brief's, taken in commissioning)" is the brief's record. Keep any choice not already stated in its method section, as a trade-off.

**[BASS.md](../../docs/Maths/Regression/BASS.md)**
- A and B: "qualification" holds the recovery table and the k̂ trust readings (17/24 in [port-fitting-checklist.md](../regression/bass/port-fitting-checklist.md)). Keep its conclusion: read the retained content's k̂ under either g.
- B:
  - "the residual birth (E4)": "The measurements", 22/22 in [proposal-kernels-brief.md](../regression/bass/proposal-kernels-brief.md);
  - "the local-move options (F-M2)": "The pass's verdict" ([local-move-pass-brief.md](../regression/bass/local-move-pass-brief.md));
  - "SMC" and "the flat-histogram sampler": readings, 14/34 and 42/69, with any 34/34 and 69/69.
- C: the "Built …" or "Stage N (date)" lead of these sections:
  - the terms and the prior, F-M2, the birth mixture, the residual birth, SMC, the flat histogram (its heading too);
  - decomposition, sensitivity, calibration, presentation;
  - the extension seams (and its "Today" and "Filled by" columns), the weighted core;
  - priors that stop reading a moving response, a design whose working response moves, random-effect intercepts.
- D: "a design whose working response moves" holds an open design call for the user's review.
- E: "The design pass's correction (2026-10-27)" in design and admission; "Until 2026-11-06 it was reported `IntegratedGivenOthers`" in the pointwise density; commit hashes in the sweep and presentation.

**[BARS.md](../../docs/Maths/Regression/BARS.md)**
- B: in "target and moves", the relocation kernel's and births' measurements (`MixingBenchmark`, `BirthBenchmark`); 7/8 in [bass-preparation-brief.md](../regression/bass-preparation-brief.md).
- E: "Until 2026-10-21 the Normal model …" in the collapsed g-prior block; the dated "Centring the Normal prior (2026-10-22, after step 3.10)"; "found by the SPC × BARS work, 2026-10-24" in orchestration.
- *Check before deleting:* orchestration and readouts matches only 14/33 (any 26/33), so some of its numbers may live only here.

**[BAPS.md](../../docs/Maths/Regression/BAPS.md), [BGP.md](../../docs/Maths/Regression/BGP.md), [BECP.md](../../docs/Maths/Regression/BECP.md), [BAMRS.md](../../docs/Maths/Regression/BAMRS.md)**
- B:
  - BAPS's "Anisotropic smoothing" carries measurement runs within its 139 lines (18/20 in the preparation brief);
  - BGP's single section has readings (13/13 there);
  - BECP's timings match only 3/7, so move them rather than delete.
- C: BAMRS's "Overview & Development Arc" is framed as history.
- E: BAMRS's "Source Analysis & Corrections" belongs with Sources&Changes.

**[ConformalPrediction.md](../../docs/Maths/Regression/ConformalPrediction.md)**
- A, about 145 lines:
  - "Checks (slice 0)": 14/16 in [conformal-layer-brief.md](../regression/conformal-layer-brief.md);
  - "Cost (slice 0)";
  - "Checks and cost (slice 2)": 27/48, any 39;
  - "Checks and cost (slice 3)": any 21/24.
  - *Move, don't delete,* slices 2 and 3: 9 and 3 of their numbers are found nowhere in `issues/`.
- B:
  - "The predictive distribution, built (slice 4)": what a case reports is method; its readings and costs match only partly (any 41/58), so move them;
  - "The per-unit bounds (between slices 4 and 5)": the carried bound's definition stays and the readings go (23/37, any 31);
  - E35's three sections: µs costs and false-alarm runs;
  - "CB-3": its readings (43/44 in the conformal brief).
- C:
  - "(slice N)" in headings;
  - the overview's "The build is the layer's brief, slice by slice";
  - "Not built yet" (stale, as above), whose plan is the brief's slices 10 to 18.
- D: "The predictive distribution (slice 4): the reading" is the investigation behind the 2026-11-08 decision (any 11/15).
- E: "the user's decision of 2026-11-07" in the split form.

**[MultiOutput.md](../../docs/Maths/Regression/MultiOutput.md)**
- A: "Checked", 6/6 in [multi-output-brief.md](../regression/bass/multi-output-brief.md).
- C:
  - headings carry "(decision N; slice N)";
  - "The remaining boundaries" is build status unless it states refusals.
  - The method content is otherwise clean.

**[PolynomialChaos.md](../../docs/Maths/Regression/PolynomialChaos.md)**
- A: "What was checked (slices 1, 2, 4 and 5)", 62/71 in [khaos-brief.md](../regression/bass/khaos-brief.md). Seven of its numbers are found nowhere in `issues/`.
- B and C: "The engine's reach (slice 5)" carries the SMC trap's measurement and its stale open decision.
- C:
  - the first line and the overview: "Slice 1 was built on 2026-11-07 …" and the build's slice list;
  - "not built" remarks in the term space, the proposal and the coefficient priors (contract boundary or status, the user's call).
- B, a gray item: "The basis (E19)": "No column cache … a degree-15 column still costs about 4.5× to 6× …" is a trade-off with its measurement.

**[RandomEffects.md](../../docs/Maths/Regression/RandomEffects.md)**
- C: headings carry "(slice N)", and the title carries "(E21's ridge blocks)".
- Otherwise method. "Migrating a saved one-group configuration" is user-facing and stays.

**[Regression.md](../../docs/Maths/Regression/Regression.md)**
- E:
  - "Configuration" carries 8 dates and 14 decision IDs in its table;
  - "Design": "joined the shell on 2026-10-23";
  - pointwise densities: "came from NaN arithmetic until 2026-10-27".
- Otherwise method.

**[VarianceMeanMixtures.md](../../docs/Maths/Regression/VarianceMeanMixtures.md)**
- C:
  - each family section opens "Built <date> ([the build brief]…), slice N", nine leads in all;
  - headings carry "(E1's slice N)";
  - the overview table has a "Built" column;
  - "Status and caveats" opens "Built: slices 0 to 7 of the brief".
- B:
  - "The GIG's draw" holds timing readings;
  - the family sections carry readings, matching only partly (13/24, 17/26, 19/27 in [gbass-noise-brief.md](../regression/bass/gbass-noise-brief.md); any up to 29/36). Move these rather than delete.
- Keep the credits ("the user's derivation of 2026-10-26", "the user's u and z parameterization").

**[PosteriorContraction.md](../../docs/Maths/Regression/PosteriorContraction.md)**: method throughout. Nothing to flag beyond a few scope remarks.

**Your two new pages**
- [Status&Caveats.md](../../docs/Maths/Regression/Status&Caveats.md) splits in two:
  - *Status, to `issues/`:* the completion contract, "through stage 4", "Built / Built since / Not yet built / Meanwhile", and the efficiency defaults waiting on measurement. Its homes are [bass-completion-report.md](../regression/bass-completion-report.md) and the [extensions checklist](../regression/bass/extensions-checklist.md).
  - *Caveats, which stay:*
    - `EbicPrior`'s count law against its location measure;
    - BAPS's heavy-tailed λ;
    - local births stranding knots;
    - units (OQ19);
    - reporting limits;
    - the knot count mixing slowly. Keep its guidance ("judge a run by ESS"); its τ table is a measurement (9/9 in the preparation brief).
- [Sources&Changes.md](../../docs/Maths/Regression/Sources&Changes.md) is legitimate docs: the register of corrections to sources, with evidence, that you decided on 2026-10-27. Strip the commit hashes and the "(E8's slice 0)"-style references from its cells; the evidence stays.

## Your calls

1. **A validation line.** Delete the verification sections outright, or keep one line per capability, such as "checked against the R oracle to 10⁻¹⁴; the report", linking the report.
2. **Decision IDs inline** (D03, OQ17, F-M2, E-numbers). Keep them as stable cross-references into `issues/`, or strip them from the method text. They are densest in Regression's configuration table and BASS.
3. **Credits with dates.** Keep the credit; the date is optional.
4. **"Not built" remarks.** One that tells a user what the API refuses is contract and stays. One that tracks the build goes to a checklist.
5. **Unique numbers.** Where a section's numbers are found nowhere in `issues/`, moving them into the brief's report before deleting keeps the record whole. This covers ConformalPrediction's checks and costs for slices 2 and 3, PolynomialChaos's checks, the VarianceMeanMixtures family readings, BARS's orchestration and BECP. I can do that move once you've settled the layout, if wanted.
6. **The structure.** Most pages are ordered by build slice. Renaming headings by topic is the larger fix, and the reorganization you're doing already covers it.
