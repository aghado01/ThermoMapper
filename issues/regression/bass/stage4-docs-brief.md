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

### The docs section — 2026-10-28

`a8cfa36`: BASS's section in `docs/Maths/BayesianRegression.md`, 270 lines added and 5 changed, CRLF kept throughout.
- **Four algorithm sections beside BARS's:** *terms and the prior* (D01, D02), *the design and admission* (D06), *the collapsed target* (D03, D05) and *the sweep* (D04). Each decision carries its alternatives, trade-offs and implications, the code that holds it and the checks behind it, with the numbers taken from the S1–S6 reports and the fixture reports. A "references" paragraph defines F-O1–F-O9 and the posterior mode's fixtures A, B, C and E before their first use.
- **Around them:**
  - BASS in the header, the overview (with a mermaid subgraph), the key concepts, the design section (values throughout, against `BarsConfig`'s instances) and the API table;
  - `BassConfig`'s rows in the configuration table, and a `BassResolved` paragraph;
  - nine rows in "Corrections to the sources": R's fixed thirds, τ rate, Jeffreys mode, mixed order scoring, `dmwnchBass`, `getQf` and tempering; Francom et al. 2018's τ shape; KHAOS's double count. Three "changed, not corrected" items: the per-λ truncation, R's intercept and σ² defaults, the subset measure;
  - two status bullets (built through the sweep, not the run; what waits on measurement) and two source-map rows.
- **Described as not built:** the run, tempering and evidence, the pointwise density, qualification, and stages 5–8. The status bullet lists them, and each section names the step (S7–S10) where a decision's consumer lands. When those steps land, their paragraphs replace these pointers: the status bullet, the config rows marked "read by the run (S7)" and the tempering row's "the ladder lands with S8".

**Review.** Four passes over the diff, each claim checked against its decision, its report or the code.
1. *Accuracy.* Fixed:
   - N was used before it was defined;
   - the three-hinge example was misquoted (1.6·10⁻³ is h₂ after h₁, not the third column);
   - the pairwise-sum errors were attributed to one column (D06's wording restored);
   - R's statistics were quoted at 2·10⁻¹⁶, where S/2 is 4·10⁻¹⁶;
   - the ABMRS ridge was tied to the wrong correction;
   - the equal schedule was called provisional, where it stands until F-M1;
   - "a normalizer it cannot compute" became D02's "intractable";
   - the support threshold's symbol b clashed with the inverse temperature.

   Added: the mixed order prior's implication (the empty combination cut by the truncation), the cost of uncentring the coefficients at composition, the Liang caveat on the g-prior display, and S7–S10 named in the intro.
2. *Gaps and definitions.* κ(E), M_max, K_max, p, n_v, Z_U and G(U) were used before they were defined; interning's "storage choice" was missing; "what it buys" read as built (now "neither is built yet"); the stage-4 plan was unlinked. Fixed.
3. *A fresh read of the final text.* Fixed: the Wallenius set S clashed with the structure S (now V, with Z̄ for the weight outside); G_c, b_c, α, |G|, ‖Rβ‖²'s R and Q were undefined in the BASS sections; F-O1's quadrature was phrased loosely; "He, Yang and Kang's code" was not D01's attribution.
4. *Nothing further.* A script checked that every added link resolves, file and heading slug, and the file is CRLF on every line.

**For the main session: inconsistencies in the briefs.**
1. *D06's "Shared state" bullet* still says the residual is updated on accepted births and recomputed after deaths and on a cadence. S4 recomputes it on every acceptance, as its report says, but D06 was not amended, as D01 was for Wallenius. The doc describes the built rule and names the change.
2. *Level centring's extra credit at M = 2:* D03's proper-defaults subsection says 0.67–0.69 nats, the F-O9 report 0.68–0.69. The doc quotes M = 1, M = 5 and structures one term apart only.
3. *D06's column-mass fixture* ("a fixture compares the predictive variance of edge terms under the two") has no checklist item and no test; S4 tested column mass as an option only. The doc does not claim it.
4. *Fixture A's posteriors differ between reports:* for M = 1–3 at M_max = 3, the S6 report's oracle values are 0.176, 0.751 and 0.073, the oracle-modes report's 0.7719, 0.1757 and 0.0525. The C# bridge writes its own data, so these are two draws of A's recipe rather than a conflict, but neither report says so. The doc quotes the C# test's numbers.
5. *D04's "τ by slice sampling on log τ, with its Jacobian"* against S6's "no Jacobian is added": both hold, since the Jacobian sits inside the target's density on u, but D04 reads as though the block adds it.

Not mine and untouched: chip A's uncommitted change to `src/maths/samplers/ensemble/ChainEnsemble.cs`, in the working tree during this work.
