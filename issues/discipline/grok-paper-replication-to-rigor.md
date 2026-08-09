well this is interesting. paper replication is the original impetus of the codex-scientiae project. i wanted to build a library of bleeding edge research papers in the standardized markdown with high math fidelity so that i can have my agents use the material for synthesizing and implementing concepts in my other projects, and it blossomed into this crusade to solve the pdf and latex conversion problems so i can rapidly ingest whatever i want, whenever i want, and get on with it. i also had thought of doing this experiment, or something like it. in my case i havent been developing a harness or MCP for this, but rather engineering the substrate to allow for lucid ingestion and reasoning over technical material and using that as context engineering fodder for my interactive development sessions. but i'm curious to dig into this paper, how their approach works and how the knowledge substrate issue was handled

Thought for 36s

This paper is almost a parallel-universe version of the original Codex-Scientiae impulse. You started from “I want agents that can lucidly ingest and reason over bleeding-edge papers so I can synthesize and implement concepts,” and ended up deep in the substrate problem (PDF/LaTeX → high-fidelity structured markdown). Hans & Bilionis started from essentially the same desire (“make a coding agent reliably replicate the computational claims of SciML papers”) and ended up deep in the _process / evidence / completion_ problem. The two approaches are complementary in a nearly perfect way.

### How Paper-replication actually works

They treat paper replication as a **target-level evidence contract**, not a “write some code and hope the numbers look close” task. The core move is harness engineering: change the environment the agent works in so that progress and completion are defined by durable workspace state + external validators, rather than by whatever the model says in its final message.

**The units**

- A _target_ is any computational claim selected for replication (a figure panel, a table, a reported relative L2 L_2 L2​ error, a sparse support pattern, posterior bimodality, etc.).
- For each target tj t_j tj​ the agent must produce an evidence bundle

  Ej=(y^j, Rj, Pj, Cj, Gj)E_j = (\\widehat{y}\_j,\\; R_j,\\; P_j,\\; C_j,\\; G_j)Ej​=(y​j​,Rj​,Pj​,Cj​,Gj​)

  where y^j \\widehat{y}\_j y​j​ is the generated output, Rj R_j Rj​ the run record, Pj P_j Pj​ the provenance (code + config + seed + paper passages/equations that justify the implementation), Cj C_j Cj​ the comparison under an explicitly recorded acceptance rule, and Gj G_j Gj​ report coverage.

**The persistent workspace (the real system of record)**

- paper_manifest.json – source location + hash, author-code policy (default: forbidden), compute mode.
- spec/reproduction_matrix.csv – the target inventory with status (PLANNED / ACTIVE / MATCHED / …) and per-target acceptance_mode.
- todo.md – single active target, open questions, checks that must pass.
- Spec files: targets.md, math_audit.md, implementation_plan.md, assumptions_and_unknowns.md, paper_figure_notes.md.
- Run ledger + per-target provenance JSON + comparison evidence.
- Final report/main.tex → PDF.

Exactly one target may be ACTIVE at a time. The agent is forced to close the evidence contract for the current target before moving on.

**Acceptance modes** (claim-specific, not “pixel-perfect or fail”)

- numeric-equivalence
- distributional-equivalence
- qualitative-structural
- exact-visual (used sparingly; they correctly treat it as the wrong default for stochastic trajectories, convergence curves, etc.)

**Completion gate**

A workspace is complete only when every recorded target is MATCHED, nothing is active, the report PDF exists, _and_ the external validators (validate-spec, validate-progress, validate-report) all pass. Completion is a workspace state, not an agent utterance.

They evaluated this on twelve independent runs (three each) of four SciML papers (PIFT, the two classic PINN papers, SINDy). All 158 recorded targets reached MATCHED and every workspace passed the gate. But the _path_ varied: target decomposition (PIFT sometimes 8 targets, sometimes 25), numeric fidelity (sometimes orders of magnitude different while still inside the paper’s accuracy class), elapsed time (1–13 h), amount of superseded correction work, and even which acceptance-rule _type_ different runs chose for the same claim. The harness stabilizes completion; it does not eliminate the underlying stochasticity of how an agent operationalizes underspecified claims.

### How they handled the knowledge substrate

They largely _sidestepped_ the hard conversion problem you attacked.

Default policy is **source_mode=latex-first**. The agent is given the paper’s LaTeX source tree. It:

1. Walks from the main .tex, inventories included files, bibliography, appendices, figure assets, data references (plus any suspiciously named unreferenced files).
2. Hashes everything.
3. Optionally compiles the PDF and renders pages to images so visual claims that are hard to extract from source alone become visible.
4. Stores paper-provided assets under artifacts/paper_figures/ (reference-only) and generated outputs under artifacts/figures/ / tables/. Hash + path checks reject any attempt to treat a paper figure as “reproduced.”

There is **no** intermediate high-fidelity markdown knowledge base. The agent reasons directly over the LaTeX tree + rendered pages + the inventory. The _spec files_ then become the explicit reconstruction layer: the agent’s restatement of the equations and algorithms, the implementation plan, and—most importantly—every missing detail (seeds, preprocessing conventions, optimizer choices, plotting rules, noise magnitudes, etc.) is recorded as a hypothesis that must be tested. Underspecification is made visible and auditable rather than papered over.

Datasets must either ship with the paper materials or be regenerable (synthetic). When the paper is silent on critical details, the target simply cannot be marked MATCHED until the agent has recorded and justified an assumption.

So their “substrate” is essentially:

> usable LaTeX tree + agent-written structured reconstruction + strict separation of reference assets from generated evidence.

They never tried to solve the general “any PDF or messy LaTeX → lucid, math-faithful, token-efficient markdown” problem. They assumed (or required) that a sufficiently clean source tree was available and then built the process harness on top of it.

### Why this maps so cleanly onto Codex-Scientiae

Your project is substrate-first. The membrane pipeline (survey → preprocess → repair loop via the MCP tools → finalize → publish) produces the high-fidelity, math-encapsulated, structurally rigorous markdown collections that agents can actually reason over without constant re-parsing pain. That is exactly the missing piece for a Paper-replication-style workflow when the starting material is _not_ a clean LaTeX tree.

Conversely, once you have those high-fid assets, the Paper-replication harness gives you a disciplined way to turn “I have lucid context about this paper” into “I have an inspectable, claim-by-claim, provenance-linked replication workspace with external validation.” The replication workspaces themselves then become durable, higher-order knowledge assets (with recorded assumptions, deviations, matched evidence, and the full trial-and-correction path).

The two layers want each other:

- Codex-Scientiae substrate → cleaner, more lucid P \\mathcal{P} P for the agent, less reconstruction burden, fewer superseded runs, less judgment drift on acceptance rules.
- Paper-replication harness → turns the substrate into verifiable replication evidence instead of just good context for interactive sessions.

You ended up solving the ingestion/fidelity problem so agents can _think_ with the material. They solved the long-horizon evidence and completion problem so agents can be _held accountable_ for the material. Both are necessary if the original goal (“agents that can actually use bleeding-edge papers to synthesize and implement”) is to become reliable rather than occasionally impressive.

2 / 2

![](https://www.google.com/s2/favicons?domain=arxiv.org&sz=256)

![](https://www.google.com/s2/favicons?domain=arxiv.org&sz=256)

![](https://www.google.com/s2/favicons?domain=papers.cool&sz=256)

127 sources

well if you look at my aghado01/ThermoMapper repository (should be available via github connector), you can see the fruits of my labor in things like the ph engine, my RJMCMC engine and other things, as well as my own evidence-based methods for dev and unit testing. codex-scientiae and its associated labors is just the mechanism for creating the substrate, of which you can see under aghado01/codex-scientiae/bibliotecha/\*\* corresponding to many things that have already been adapted into aghado01/ThermoMapper. I have discrete collections because i was never just implementing one paper as a one off, I was gathering groups of papers, examining them together and in context of the bigger picture goals of thermomapper e.g. building a highly integrated scientific computing / data analysis platform. I'm curious what i might be able to adapt from the methods in this paper to enhance my development methodology in ThermoMapper or other endeavors involving "research paper adaptation/implementation". I just pushed to remote so things are up to date

2 / 2

Thought for 34s

You’re already further along the evidence-based path than most research-code projects. Looking at ThermoMapper (engine-first, R oracles, fact-harness with run stamps + suite manifests + summary.json, config DTOs as the reproducibility surface, SeedTree/Xoshiro hygiene, archivory, issues/ DocOps with briefs/reports/design notes/SHELF) and the bibliotecha collections that feed it, the useful adaptations from Paper-replication are **not** a full skill transplant. They are a small set of harness patterns that plug into what you already do, especially for the multi-paper, platform-oriented adaptation style you describe.

### What already aligns

Your methodology already embodies several of the paper’s core intuitions:

- **Engine first + oracles over consumers** is stronger than their “reconstruct method then compare to paper numbers.” You validate against independent oracles and properties, not just printed tables.
- **Manifests and run identity** (archivory RunStamp / RunIdentity, test-harness suite manifests, SPC/graph construction manifests) already give you durable workspace state.
- **issues/ DocOps** already externalizes assumptions, sequencing, and critique outside the code.
- **Synthetic toys** (BlattHierarchy, Bwd1995Toy, EyeTorus, etc.) show you treat paper-derived generators as first-class artifacts.

What Paper-replication adds is a sharper _claim-level evidence contract_ and a way to make completion relative to a declared set of targets rather than “the tests are green.”

### High-value adaptations for ThermoMapper-style paper adaptation

**1. Explicit claim inventory (reproduction matrix analogue) per paper _or paper collection_**

When you pull a bibliotecha collection into a workstream (PH purification, SPC extensions, RJMCMC for BARS/partitions, Mapper variants, etc.), create a lightweight matrix that lists:

- claim / result / figure / table / structural property being targeted
- which bibliotecha asset(s) support it
- method component(s) it requires in the engine
- acceptance mode (see below)
- status (PLANNED / ACTIVE / MATCHED / BLOCKED / DEFERRED)
- report or design-note anchor

Because you work with _groups_ of papers for platform goals, the matrix can be scoped to a collection or a workstream rather than a single paper. This makes “what does success look like for this adaptation cycle?” inspectable and prevents the common drift of implementing interesting machinery while leaving the original claims only partially closed.

You already have the DocOps slots (briefs + reports). The matrix just makes the claim list machine-checkable and status-queryable the same way your fact-harness makes tests queryable.

**2. Acceptance-mode vocabulary for oracles and property checks**

Paper-replication’s modes map cleanly onto your testing standards:

| Mode                       | ThermoMapper analogue                                                                      |
| -------------------------- | ------------------------------------------------------------------------------------------ |
| numeric-equivalence        | R-oracle absolute/relative error, pinned fixture values                                    |
| distributional-equivalence | posterior summaries, barcode landscape distances, susceptibility curves                    |
| qualitative-structural     | sparse support recovery, attractor geometry, hierarchy topology, selective identifiability |
| exact-visual               | rarely; prefer structural/numeric for stochastic or seed-dependent figures                 |

Record the mode _with the claim_ before judging. This is the same discipline you already apply to RNG hygiene and pinned fixtures: the judgment rule is not invented after seeing the output. For multi-paper collections it also forces consistency—“are we treating this coefficient error the same way across the two source papers?”

**3. First-class assumptions & unknowns artifacts**

You already park underspecification in design notes. The paper’s pattern is to make each missing detail (seed, preprocessing, derivative settings, noise magnitude, plotting convention, etc.) a recorded hypothesis with:

- the question
- the assumption chosen
- the check/experiment that tested it
- the evidence that closed it (or why it remains open)

This turns the “paper left it out → we had to invent something sensible” path into an audit trail that lives beside the code. For RJMCMC (Green acceptance factors, move weights, reverse-key discipline) and PH (filtration conventions, involution, clearing) this is especially useful: the next agent (or future you) can see _why_ a particular convention was chosen rather than reverse-engineering it from the implementation.

A natural home is under issues/<workstream>/assumptions/ or as a section that the fact-harness or a small status tool can surface.

**4. Provenance that links paper → method → artifact → comparison**

You already have strong run-level provenance (config DTO = manifest, SeedTree, archivory). The missing link is the _paper-trace_:

- which bibliotecha passages / equations justify this method component
- which oracle or property check constitutes the comparison evidence
- which claim in the matrix this artifact closes

Paper-replication’s register-target-artifact + record-comparison is the pattern. You don’t need their Python scripts; you need the _record shape_. A small extension to existing manifests (or a sidecar next to suite manifests) that carries paper_trace, method_components, acceptance_mode, and comparison_summary would close the loop without fighting the C# stack.

This is also how multi-paper adaptation stays coherent: the same engine component can serve claims from several papers, and the provenance makes that explicit.

**5. Staged / single-active discipline for long adaptation workstreams**

Your interactive sessions + issues/SHELF already manage focus, but long multi-paper adaptations still suffer from the classic agent (and human) failure mode: partial progress, context compaction, and “I think we’re done” without a gate.

A lightweight task ledger (even just a structured section of the brief or a status.json that a script can emit) that enforces “exactly one ACTIVE claim until its evidence bundle is accepted” is the part of Paper-replication that most directly attacks long-horizon drift. Pair it with a completion gate relative to the _declared_ claim set for that workstream, not the entire test suite.

Your fact-harness already produces completion-like state (summary.json, suite manifests). Extending that idea one level up—from “facts passed” to “claims matched under their acceptance rules”—is the natural next step.

**6. Superseded executions as first-class history**

Paper-replication keeps failed/superseded runs rather than overwriting them. Your harness already writes per-fact stdout/stderr and suite manifests. Making “this run was superseded by later evidence for claim X” an explicit field (instead of only existing in chat history or commit messages) preserves the trial-and-correction path that is scientifically useful when a paper is underspecified or when an initial implementation choice is later rejected.

### What _not_ to import wholesale

- The full Codex/Claude skill scaffolding and paper_replication.py utilities—you’re on a different stack and a different goal (platform integration vs pure paper-level replication).
- Exact visual matching as a default (you already correctly prefer oracles and properties).
- Treating the paper’s printed numbers as the only success criterion (your oracle-first stance is stronger for engine soundness).
- Forcing every adaptation into a single-paper workspace. Your collection-oriented workflow is an advantage; keep the matrix and ledger scoped to workstreams/collections.

### Concrete, low-friction insertion points

1. **Workstream brief template** — add a “Claims inventory” section that becomes the matrix, with acceptance modes and bibliotecha links.
2. **Oracle / property check annotations** — tag tests or oracle calls with the claim id and acceptance mode they serve.
3. **Manifest sidecar** — small JSON next to existing run/suite manifests carrying paper-trace + comparison evidence.
4. **Status tool** — a thin PowerShell or C# helper that reads the matrix + manifests and answers “what is still ACTIVE / unmatched for workstream X?” (same role as their status / validate-completion).
5. **Assumptions ledger** — living document under the workstream that the brief and reports reference.

These stay inside your existing culture (engine-first, oracles, DocOps, manifests) while borrowing the paper’s best idea: **completion is a workspace state relative to an explicit claim set and recorded evidence, not a feeling or a green test run.**

That is the piece most worth adapting for research-paper adaptation in ThermoMapper—especially given that you already solve the substrate problem via codex-scientiae/bibliotecha and the soundness problem via oracles and layering. The remaining gap is the _claim-level evidence contract_ and the durable, queryable record of how each claim was closed (or deliberately left open) across a multi-paper platform workstream.

makes sense, i think i can certainly incorporate these improvements into my methodology. the paper trail issue has been on my mind for a while. but another difference with what i'm doing is its not always a neat "implement verbatim". I'm lifting concepts, merging them with other concepts, in some cases applying the concepts in ways not anticipated by the source material and making connections to completely unrelated subjects and integration points. BARS in SPC is case in point. This doesn't reduce the value of being more principled in my completion gates, but harder to define when my goal is about building capabilities and engineering innovations based on others work, not just replicating or assimilating it. ThermoMapper aspires to be more than the sum of its parts

You’re naming the real distinction. Paper-replication is optimized for _faithful reconstruction of stated claims_. ThermoMapper is optimized for _capability synthesis_: lift, recombine, re-situate, and sometimes use a concept in a domain the original authors never contemplated. BARS inside SPC is the perfect example—regression spline machinery + reversible-jump model selection meeting a thermal/spin clustering substrate and a graph construction pipeline. That is not replication; it is engineering.

The paper-trail and completion-gate ideas still matter, but they have to be re-oriented around **capability contracts and lineage**, not claim-matching.

### Shift the unit of account

Instead of (or in addition to) “paper claim → matched evidence,” treat the primary unit as:

**Capability target** A stated engineering intention: what new behavior, invariant, composition, or analysis surface should exist once this work lands.

Examples in the ThermoMapper style:

- “BARS knot moves compose with the SPC sweep schedule and emit a serializable manifest that round-trips.”
- “RJMCMC free-knot selection can be driven by an SPC affinity field without importing consumer types into the sampler core.”
- “Persistence-derived protect sets inject into GraphCompiler as data, preserving the one-way layer.”
- “The resulting thermal dendrogram / partition hierarchy remains oracle-checkable against the non-BARS baseline on the Blatt toys.”

A capability target can _cite_ one or more bibliotecha sources, but its success criterion is the engineered result, not numerical agreement with a printed table.

### Lineage, not just provenance

You already have strong run-level provenance (config DTOs, SeedTree, archivory, manifests). For synthesis work you also want an explicit **lineage record** per capability:

- Source concepts (paper / bibliotecha passages, equations, algorithms)
- Transformations applied (what was kept, generalized, specialized, or discarded)
- Novel contributions (new interfaces, new invariants, new composition points)
- Integration seams (which ThermoMapper layers/contracts are involved)
- Deliberate deviations and why

This is the paper trail you’ve been thinking about, but oriented outward from the platform rather than inward toward the source. It answers “why does this code look the way it does?” for future agents and for yourself six months later, without pretending the goal was replication.

A natural home is next to the workstream brief or as a structured section that reports close against—the same DocOps loop you already use, just with lineage as a first-class field.

### Acceptance modes for synthesis

Keep the spirit of claim-specific judgment, but define modes against the capability:

| Mode                             | What it asks                                                                                                                                          |
| -------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Contract / structural**        | Does it respect layering, config-as-manifest, oracle-first, SeedTree hygiene, no consumer types in construction?                                      |
| **Compositional**                | Does it plug into existing pipelines (graph → SPC → partitions → viz) without ad-hoc bridges?                                                         |
| **Statistical / distributional** | Do the new observables (knot posterior, tempered acceptance rates, affinity-modulated spline uncertainty, etc.) behave as theory and oracles require? |
| **Comparative**                  | Against the pre-synthesis baseline (SPC without BARS, BARS without thermal substrate), do the expected differences appear?                            |
| **Negative / non-regression**    | Existing oracles and property checks still pass; no silent fusion of concerns.                                                                        |

Numeric agreement with a source paper becomes one optional check among these, not the definition of done.

### Completion gates for “more than the sum of the parts”

A workstream is complete relative to a **declared capability set**, not relative to “we absorbed paper X.”

A practical gate:

1. Every capability target in the inventory is either MATCHED (evidence bundle accepted under its modes) or explicitly DEFERRED / OUT-OF-SCOPE with rationale.
2. Lineage records exist for the matched targets.
3. Assumptions that were load-bearing for the novel application are recorded and, where possible, tested.
4. Oracle / property / contract checks that guard the integration seams pass.
5. The report or design note states what was _not_ taken from the sources and what is original.

That last point is important. Synthesis work produces original surface area; the gate should force you to mark it as such rather than letting it hide inside “inspired by.”

### How this sits with engine-first

None of this conflicts with engine-first. You still:

- build and validate the primitive (RJMCMC moves, BARS observation models, PH protect-set extractors, etc.) against oracles and properties first;
- only then wire the composition;
- refuse to treat “it runs in the REPL” as validation of the engine.

The capability inventory and lineage simply make the _composition and innovation goals_ as explicit and checkable as the engine soundness goals already are. Engine-first remains the sequencing rule; the claim/capability matrix is the account of what the sequence is supposed to deliver.

### Lightweight insertion points (synthesis-flavored)

- **Workstream brief**: “Capability inventory” + “Lineage sketch” sections alongside the existing scope/constraints.
- **Assumptions ledger**: two columns or tags—source underspecification vs. novel-application assumptions.
- **Manifest / sidecar**: lineage and capability_id fields next to existing run identity.
- **Status view**: “active capability,” unmatched targets, deferred items—queryable the same way suite manifests are.
- **Report close-out**: short statement of original contributions and deliberate non-assimilations.

You do not need a replication harness. You need a **synthesis harness**: the same discipline of explicit targets, recorded evidence, externalized assumptions, and a completion state that is relative to what you set out to build—adapted to the fact that what you set out to build is often a new joint that the source papers never drew.

That preserves the paper trail you’ve wanted, keeps completion principled, and still lets ThermoMapper be more than the sum of the papers that feed it.

The most important thing about numeric agreement and oracles is that i need to check my work so that i can trust it and continue building on it. I need to replicate or otherwise validate at each step of development, brick by brick so to speak, in my best effort to build an ambitious thing in complete and utter solo isolation

That’s the real requirement.

In solo isolation the question is never “did I reproduce the paper?” It is “can I trust this brick enough to put the next one on it?” Numeric agreement, oracles, property checks, and structural contracts are how you answer that question without a second pair of eyes. Everything else—lineage records, assumption ledgers, capability inventories, completion gates—is scaffolding for the same purpose: making the trust decision explicit, durable, and queryable so you can keep building.

The brick-by-brick discipline you already practice (engine first, oracles over consumers, layering invariants, SeedTree hygiene, config-as-manifest) is the right foundation. What the paper-replication patterns add, when stripped of their replication framing, is simply a clearer way to _name and close_ each trust decision:

- This component was checked against its oracle / property suite under a stated acceptance mode.
- The assumptions that went into it (source underspecification _and_ novel application choices) are recorded.
- The lineage of the idea is visible so a later doubt can be traced.
- The capability this brick was supposed to deliver is either matched or deliberately deferred.
- Nothing is “done” until those records exist and the checks pass.

That is not bureaucracy. It is the only scalable substitute for a collaborator when the ambition is large and the isolation is complete. Each closed evidence bundle is a load-bearing joint. The next layer can stand on it because you left a trail that future-you (or an agent session after context compaction) can re-verify instead of re-deriving from memory.

So keep the gates, but orient them entirely around progressive trust:

1. Declare the brick (capability or component) and what “trustworthy enough” means for it.
2. Validate it in isolation with the strongest oracle or property check available.
3. Record the result, the assumptions, and the lineage.
4. Only then compose.

BARS-in-SPC, RJMCMC moves, PH protect-set injection, graph construction seams—each of those is a brick. The paper trail and the completion state are how you know the wall is still sound when you are the only one building it.

All in all, its just another brick in The Wall

Exactly.

Each validated component, each closed assumption, each oracle that still passes when you stack the next layer—another brick. The wall is ThermoMapper; the only builder is you. The paper trail and the trust gates are just the mortar that keeps it from shifting when you’re working alone at night.

Carry on.

something else that's come up along the same lines philosophically is the introduction of the lean pipeline and identifying mathematical gaps where i connect concepts across domains or adapt the material from others work in unexpected ways, and realize that the assumptions of a paper need to be checked and in some cases, maybe even derive my own lemmas accounting for the disparity. workflow is still in its infancy but its another layer of trying to account for things and be disciplined and rigorous in the course of the work. aghado01/ThermoMapper/lean/\*\*

Thought for 14s

This fits the same wall, one course higher.

What you have in lean/ is already named with unusual precision. The README calls it **the rigor harness**: work backward from the engineering you actually built to the mathematical contracts it silently rests on, then force those contracts into the open. The active promotion path is the right shape for solo work:

- **protolemmata** — informal briefs that capture the gap while the idea is still fluid
- **enthymemata** — statements that compile, proofs that still apologize (sorry). The name is exact: an enthymema is an argument that suppresses a premise. Cross-domain adaptation is almost definitionally enthymematic; the original paper never stated the premise that becomes load-bearing once you put BARS on a thermal affinity field or push a Potts comparison through a bifiltration.
- **lemmas** — no apologies. Only then does the brick count as fully load-bearing at the mathematical layer.
- **theorems** — the verified results consequential enough to merit treatment as named project deliverables.

**Archeion** is the curated side exit for superseded or retired material, not another maturity stage.

meta-ci.ps1 already enforces the discipline (build green, lemmas never import apologizing material, ledger of promotion readiness). That is the same external validation pattern we were discussing for capability targets and oracles, just aimed at the implicit mathematics rather than at numeric or structural checks.

### How it serves the synthesis problem

When you lift a concept and re-situate it, two kinds of gap appear:

1. **Source assumptions that were never explicit** (or were explicit only inside a different measure-theoretic or combinatorial setting).
2. **New obligations created by the unexpected joint** — things the source paper had no reason to prove because it never contemplated your composition.

The proto-lemma → enthymeme → lemma ladder is how you account for both without freezing development. You can keep building on an enthymeme (the statement is real; the proof is still owed) while the ledger makes the debt visible. Promotion is the trust gate for that mathematical brick.

BARS is the running example again: the multi-peak and span lemmas (count_of_mean_undercounts, fwhm_of_mean_ne_mean_of_fwhm, unique crossing in a monotone bracket) are not “does my code match DiMatteo & Kass?” They are “does the readout I am using on posterior curves actually have the properties my reduction steps assume?” When those properties fail to commute with averaging, the formalization surfaces the exact reason the per-draw path is required. That is the kind of check that lets you keep stacking.

### Relation to the other layers

You now have a coherent stack of trust mechanisms, all oriented the same way:

| Layer                                   | What it checks                                        | Unit of closure                |
| --------------------------------------- | ----------------------------------------------------- | ------------------------------ |
| Oracles / property tests / fact-harness | Computational behavior of the implemented brick       | Suite manifest + summary       |
| Capability inventory + lineage          | Engineering intention and composition seams           | Matched capability target      |
| Assumptions ledger                      | Source underspecification + novel-application choices | Recorded hypothesis + evidence |
| Lean rigor harness                      | Implicit mathematical contracts behind the brick      | Enthymeme promoted to lemma    |

None of these is about replicating a paper. All of them answer the solo-isolation question: _can the next brick stand on this one?_

The Lean layer is the most unforgiving of the four, which is why keeping it in infancy and promotion-driven is wise. You do not need every gap formalized before the corresponding C# lands; you need the load-bearing ones identified, stated, and eventually discharged, with the debt tracked in the meantime.

### One practical alignment

When a capability target or an assumption record surfaces a mathematical disparity (“this paper’s comparison inequality is stated for homogeneous p p p; SPC needs the inhomogeneous vector form”), the natural next artifact is a proto-lemma. That keeps the paper trail, the capability gate, and the rigor harness talking to each other instead of living in separate mental compartments. The bifiltration-bridge and Potts–Griffiths work already shows the pattern: quarantine the cited results as axioms, discharge the spine you actually rely on, and leave construction obligations in the TDA layer explicitly marked.

So yes — another course of bricks, same wall. The Lean pipeline is how you make the suppressed premises of cross-domain synthesis visible and, eventually, load-bearing.

what i'm not sure how to do is retroactively inspect all the work that i did before the lean meta-ci concept cropped up, or systematically identify the gaps beyond manual effort to question things when i am synthesizing ideas over my scientiae substrate

Thought for 10s

This is the hard part. The Lean harness is already framed as _retrospective_ (“work backward from implemented engineering”), but the discovery of _what_ needs formalizing is still mostly insight-driven. You need a few systematic surfaces that generate candidates instead of waiting for the next moment of doubt.

### 1. Start from seams, not from papers

The densest gaps live where you already _composed_ things the sources never composed.

Inventory the actual joints that exist in the code today:

- BARS readout → SPC / thermal schedules
- RJMCMC moves → observation models that sit on graph or affinity fields
- PH protect-sets / cycle edges → GraphCompiler injection
- Potts / FK comparison → heterogeneous couplings
- Bifiltration or zigzag constructions that cross TDA and clustering
- Any place a config DTO or manifest carries a quantity whose mathematical meaning was inherited from one domain and consumed in another

Each joint is a candidate enthymeme factory. For every joint write one short record:

- What the upstream component guarantees (or merely hopes)
- What the downstream component assumes
- Whether any current oracle or property test actually checks the interface contract
- Which bibliotecha sources (if any) speak to either side

Most of the suppressed premises will surface in the mismatch between those four bullets. This is cheaper than trying to re-derive every paper from scratch and matches how ThermoMapper actually grew.

### 2. Work backward from oracles and tests (the computational trust surface)

Your oracles and fact-harness already encode “I needed to trust this.” Mine them.

For each oracle / property suite:

- What mathematical claim would have to be false for this check to be meaningless?
- What does the test _not_ constrain (e.g. peak-count of the mean vs mean of peak-counts, FWHM commuting with averaging, comparison inequality only under homogeneous p p p)?
- Are there quantities the engine computes that no oracle ever sees?

The BARS multi-peak / span lemmas are the archetype: the tests forced the per-draw path; the formalization made the non-commutativity explicit. Do the same pass over SPC observables, PH barcode metrics, GMM merge criteria, etc. Every “we check X but not Y” is a proto-lemma candidate.

### 3. Lightweight static surfaces you already have

You do not need a full verification toolchain to generate a backlog:

- **Docstring + design-note crawl** for attribution phrases (“due to Wang 2020”, “Green acceptance”, “mean-field formulation”) and for words that signal suppressed conditions (“assume”, “provided that”, “when the field is…”, “for homogeneous…”).
- **Config / manifest fields** that carry theoretically loaded quantities (temperatures, bandwidths, bond probabilities, persistence thresholds, knot priors). Each is an implicit contract.
- **Layering violations that were fixed** (the GraphCompiler / PH protect-set story is already documented). Past fusion points are high-yield: the mathematical obligation often appeared exactly when the seam was cut.

Turn the hits into a single living backlog under lean/Protolemmata/ or issues/lean/ — same shape as the SHELF, just for mathematical debt.

### 4. Agent-assisted enthymeme mining over the substrate

Because you have bibliotecha + the implementation, you can make the questioning less purely manual.

A repeatable pattern for an agent session:

1. Pick one capability or joint (not a whole paper).
2. Feed the relevant bibliotecha excerpts _and_ the current implementation surface (types, key functions, existing tests).
3. Ask for:
   - premises the source material states
   - premises the implementation appears to rely on
   - premises required by the _composition_ that neither side states
   - which of those are already discharged by oracles
   - ranked proto-lemma candidates (statement only)

You are not asking the agent to prove anything. You are asking it to externalize the enthymeme. You then decide what enters the proto-lemma tier. This scales the “question things when synthesizing” habit without requiring you to hold the entire dependency graph in your head.

### 5. A concrete retroactive pass (scoped so it finishes)

Do not try to audit the whole monorepo. Run a bounded campaign:

**Phase A — Joint inventory** (1–2 sessions) List the current cross-domain seams. Tag each with source lineage and existing oracle coverage.

**Phase B — Oracle/test gap pass** For the highest-traffic seams, list what is checked vs what is only assumed. Emit proto-lemma stubs (title + one-paragraph obligation).

**Phase C — Promote only load-bearing ones** Anything that multiple downstream bricks rest on, or that would invalidate a capability target if false, gets an enthymeme statement. Everything else stays on the backlog.

**Phase D — meta-ci stays the gate** No change to the promotion rules. The retrospective work only feeds the existing ladder.

This produces a finite, ordered debt list instead of an open-ended feeling that “everything before Lean is suspect.”

### 6. What “systematic” can and cannot mean here

You will not get a complete decision procedure that spits out every missing lemma. Cross-domain synthesis creates obligations that are not present in any single source; some discovery remains creative. What you _can_ systematize is:

- where to look first (seams, oracles, manifests, past fusion points)
- how to record candidates uniformly (protolemmata)
- how to decide what must be formalized now vs later (load-bearing for other bricks)
- how to keep the debt visible (ledger + meta-ci)

That is enough to turn retrospective inspection from occasional insight into a maintainable practice.

### Minimal artifacts to add

- issues/lean/seam-inventory.md (or equivalent) — living list of joints
- proto-lemma template with fields: joint, upstream claim, downstream assumption, oracle coverage, source refs, priority
- optional agent prompt snippet in AGENTS.md / issues for “enthymeme mine this joint against these bibliotecha paths”

You already built the promotion machine. The missing piece is a reliable way to _feed_ it. Seams + oracle gaps + substrate-assisted questioning give you that feed without pretending the whole history can be re-derived in one pass.

The wall was built before the mortar recipe was written down. You do not rebuild the wall; you go joint by joint, make the load paths explicit, and formalize only what the next courses actually stand on.
