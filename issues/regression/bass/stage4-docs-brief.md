# BASS stage 4: the docs section (chip C brief)

Date: 2026-10-28. A chip of the [fitting segment's stage-4 plan](port-fitting-brief.md#stage-4-plan-decided-2026-10-28), which the user decided on 2026-10-28. S10's docs section, written now: D01–D06 are decided and S1–S6 are built and reviewed, so their decisions, alternatives and trade-offs are settled. The run (S7, chip A), tempering (S8) and the pointwise density (S9, chip B) are still being built, and their paragraphs are added when they land. Decided statistical choices land in `docs/Maths` with their alternatives, trade-offs and implications, so a reader can see what was chosen, what was not, and what follows.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md) (DocOps, the vocabulary, positive declarative prescriptions). Then:
- the [hub](../bass-full-port-brief.md): its scope, the lineage register, the parity rule and the representation table;
- the [fitting segment](port-fitting-brief.md): D01–D06 in full, with D03's proper-defaults subsection, the contract shapes, the stage-4 plan and the regime note;
- the [checklist's reports](port-fitting-checklist.md#reports): S1–S6 with their reviews, F-C2, the fixtures F-O1–F-O9 and the design pass;
- `docs/Maths/BayesianRegression.md` itself, whose BARS sections are the model for structure, depth and voice.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Switch no branch, push nothing.
- **Your file:** `docs/Maths/BayesianRegression.md` only, besides this brief. Chips A and B edit code; the main session owns the fitting segment and its checklist.
- **Line endings: `docs/Maths/*.md` are CRLF.** Edit with tools that keep them, and check `git diff --stat` after every scripted edit. A whole-file line count in the diff means the endings flipped.
- **Commit by explicit path,** in the repo's message style (`docs(bass): …`), ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Write what is decided and built.** Invent nothing. Where something is not built yet, say so in "Status and caveats" rather than describing it. Where the sources disagree with what was decided, follow the decision and say what changed.
- **Report** under "Reports" at the end of this brief, newest first: what landed, the commit, the review passes (below), and anything in the briefs you found inconsistent, which the main session fixes.

## The task: BASS in `BayesianRegression.md`

BASS's sections sit beside BARS's, in the document's own structure.
- **Overview and key concepts.** BASS's place:
  - Bayesian adaptive spline surfaces, reversible jump over hinge-product terms;
  - its relation to BARS through the shared collapsed g-prior block;
  - its response kinds: scalar, categorical predictors, direct functional responses.
- **Algorithms, as built,** one subsection per area. Each decision gets its alternatives, its trade-offs and its implications.
  - *Terms and the prior* (D01, D02):
    - the term representation;
    - the term measure, with the mixed order prior decided 2026-10-28, the categorical subset measures and the complement pairs;
    - the count prior with λ integrated out, under global truncation, and its tilt;
    - the proposals: Wallenius' variable sets, by the draw-order recursion, not inclusion–exclusion; the pseudo-counts w₁ and w₂;
    - the move schedules.
  - *The design and admission* (D06):
    - the centred basis and why the marginal reads only the column space of [1, Z];
    - the order-free angle gate with its strict confirmation;
    - support per part;
    - the refusal reasons;
    - the residual's explicit recomputation on every accept.
  - *The collapsed target* (D03, D05):
    - the proper defaults and the trilemma table;
    - τ's unit-information family, Zellner–Siow, and the correction of Francom et al.'s shape and R's rate;
    - the Liang convention and the zero centre as options;
    - fixed g;
    - propriety and the exact-fit rule;
    - the integrated and linear forms, with their endpoints and where each is defined.
  - *The sweep* (D04):
    - the collapsed order rule;
    - τ's slice step on ln τ;
    - composition at retention;
    - the draw's uncentred coefficients;
    - how the design follows the chain.
- **Configuration.** `BassConfig`'s fields with their defaults, and what `BassResolved` picks from the data.
- **Corrections to the sources.** R BASS's defects that the port corrects rather than copies (the parity rule):
  - the fixed-thirds move ratios;
  - the τ prior's rate;
  - the Jeffreys mode;
  - the mixed order scoring;
  - `dmwnchBass`'s lost digits;
  - `getQf`'s pivot rule;
  - the swap density and the tempered prior factors (cite D05).

  Then the papers' statements read differently: Francom et al. 2018's τ shape, and the per-λ truncation.
- **Status and caveats.** What is built (S1–S6) and what is not yet (the run, tempering, the density, prediction and the later stages). Say what to read meanwhile.
- **Source map.** The files under `src/maths/regression/spline/bass/` and the shared pieces they use: `CollapsedGPrior`, `Wallenius`, `BorderedCholesky`, `KhatriRao`, `ReversibleJumpChain`, the per-replica blocks.

## Review before reporting

The section is not done when it is written. Before the report, review it as a separate pass, as the main session reviewed each build step:

1. **Reread the whole change** as a reader seeing it for the first time: `git diff <the commit before yours>..HEAD -- docs/Maths/BayesianRegression.md`.
   - Check every claim against its source: the brief's decision, the report, or the code at `file:line`.
   - Check every number against the report it came from.
   - Make sure the intent of every decision survives your wording, and that nothing is invented.
2. **Hunt the gaps:**
   - a decided choice with no alternatives or trade-offs beside it;
   - a built behaviour the section does not mention;
   - a statement about something not built yet;
   - a link that does not resolve;
   - a term used before it is defined.
3. **Fix what you find,** and recheck the line endings.
4. **Review again** if the pass found anything substantive. Stop when a pass finds nothing.
5. **Record it** in your report: what each pass found, what was fixed, and any inconsistency in the briefs themselves, for the main session.

## Reports
