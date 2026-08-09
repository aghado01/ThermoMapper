# Literature follow-up — the PH confidence program

The citation agenda for [confidence-pushforward-lemmas.md](confidence-pushforward-lemmas.md) and its candidate
families (SL / ZZ / NV / OP / TR), produced by walking the **codex-scientiae** corpus inventory against the claims
the chain leans on. Corpus root: `codex-scientiae/compendia/{ph,mapper,statistics,bars,...}`.

> **Headline — the corpus is far more complete than the candidate cites implied.** Reading the four relevant
> compendia turned what looked like ~8 speculative acquisitions into **two real ones, one in-corpus paper to
> elevate, and one self-correction**. The discussion-arc threads are mostly *grounded already*; the gaps are
> narrow and named.

Status: **HELD** (in corpus) · **GAP→ACQUIRE** (named, absent) · **ELEVATE** (held, under-exploited) · **PARK**
(held-enough, or not yet load-bearing).

---

## The three actions that matter

1. **GAP→ACQUIRE — frequentist confidence sets for diagrams.** Fasy, Lecci, Rinaldo, Wasserman, Balakrishnan,
   Singh, *Confidence sets for persistence diagrams* (Ann. Stat. 2014). This is the **direct prior art for L3**:
   a confidence band on the function pushed through stability (L2) to a bottleneck confidence *set* on the diagram
   — the chain's L3 is the Bayesian-posterior twin of exactly this. The `statistics` compendium holds only
   TKH2022 (topological ABC); this canonical paper is absent. Cite it at L3/L4 the way BCKL2010 is cited at L2.
   *(Subsampling variant: Chazal et al, *Subsampling methods for persistent homology*.)* **Priority: HIGH.**

2. **GAP→ACQUIRE — random-Mapper stability (grounds NV-2).** The `mapper` compendium holds only PNV20XX (NV-1,
   the deterministic-cover faithfulness) and GLL2026 (applied). The **stochastic-cover** leg — the SW-built cover
   is random, each realization ε-good yet the ensemble nerve varies — has no source. Acquire **Carrière–Michel–
   Oudot**, *Statistical analysis and parameter selection for Mapper* (JMLR 2018) and **Brown–Bobrowski–Munch–
   Wang**, *Probabilistic convergence and stability of random Mapper graphs* (J. Appl. Comput. Topol. 2021).
   ThermoMapper layer-B/C correctness rides on this. **Priority: MEDIUM-HIGH.**

3. **ELEVATE — DS2026 is already paid for and under-exploited.** `ph/DS2026.md`, *Quasi Zigzag Persistence: A
   Topological Framework for Analyzing Time-Varying Data* (ZZ-GRIL). Two payloads the arc needs and didn't cite:
   - its **quasi-zigzag *bifiltration*** is the *principled formal object* for "a warped curve through two ordered
     axes" — the thing you flagged as **not** a proper bifiltration. DS2026 is the honest middle path between a
     single slice and the fraught full bipersistence (REF-MPH §9.2.1).
   - **§3.1 carries the zigzag stability** (ZZ-GRIL), so together with **FH2024 §4** (block-extension/interleaving
     for zigzag modules — already the L5 cite) it **closes ZZ-1 from inside the corpus**. This is the
     self-correction below. **Action: read it; weave into ZZ-1 (done), the vocabulary warped-curve note, and the
     SL family.**

### Self-correction (already applied)
ZZ-1 in the proto-lemma was flagged as needing **Botnan–Lesnick** (zigzag algebraic stability) acquired. That was
written before the inventory: **FH2024 §4 + DS2026 §3.1 are in the corpus and supply the leg.** Botnan–Lesnick is
demoted to *optional* (cleanest general constant only). The ZZ-1 line is fixed accordingly.

---

## Corpus status by thread

| Thread | Claim | Status | In-corpus sources |
|---|---|---|---|
| **L1–L4 inference spine** | BARS curve posterior; PH posterior consistency | **HELD** | `bars/*` (BD2005, DMGK2001, **GRE1995** RJ-MCMC, WLK2008, HYK2024, MRA2015); WRD2025, MNO2019, MMO2019 |
| L2 stability (function) | `d_b ≤ ‖·‖∞` + minimax rate | **HELD** | REF-PH §5.2 + §6; BCKL2010 |
| L3 confidence pushforward | band → bottleneck set | **GAP→ACQUIRE** | *(Fasy 2014 — frequentist twin, absent)* |
| L5 summary | landscape mean; KDE density | **HELD** | FH2024, MMO2019; REF-PH §7.3 |
| **SL** slice / multiparam | essential cover; which-slice | **HELD** | **AL2026 §4** (essential cover); REF-MPH; SGW2025; FH2024 §2.2.2 |
| SL fibered middle-ground | vary γ, watch barcode | **PARK** | *(vineyards / RIVET — absent, optional)* |
| **ZZ** zigzag structural | decomposition + stability | **HELD** | CDSM2009, DLST2026, AL2026 (decomp); **FH2024 §4 + DS2026 §3.1** (stability) |
| **NV-1** nerve faithfulness | nerve ≈ space | **HELD** | PNV20XX (`mapper`) |
| **NV-2** random cover | stochastic-cover stability | **GAP→ACQUIRE** | *(Carrière–Michel–Oudot; Brown et al — absent)* |
| **OP-2** harmonic = barcode | `dim ker Δ_k = b_k` | **HELD** | QW2024, STGW2024 §4.1, Eckmann (classical) |
| **OP-1** non-harmonic unstable | spectrum jumps | **PARK** | QW2024 (likely covers); *(Mémoli–Wan–Wang 2022 if a crisp statement is wanted)* |
| ⊕ phase-transition target | feature = transition / bifurcation | **HELD** | MR2026 (entropy); DLST2026 (Conley-Morse); SGL2022 (XY witness) |
| directed far rung | path / Mayer homology | **HELD** (theory) | KGW2026, WGW2023, QW2024 ch4 |
| **TR** triangulation | independent-signal corroboration | **OPEN** | *(no TDA-specific source; general statistics)* |

---

## Park list (named, but not yet)

- **OP-1** — read QW2024's persistent-Laplacian chapters for a spectral-instability statement before acquiring
  **Mémoli–Wan–Wang**, *Persistent Laplacians: properties, algorithms and implications* (SIAM JMDS 2022). The
  dissertation may already carry the non-harmonic behavior OP-1 needs.
- **SL fibered barcode** — **vineyards** (Cohen-Steiner–Edelsbrunner–Morozov, SoCG 2006) and **RIVET**
  (Lesnick–Wright) only when one slice can't be justified and you watch the barcode vary with γ.
- **Directed backbone** — **Chowdhury–Mémoli**, *Persistent path homology of directed networks* (SODA 2018) +
  **Grigor'yan–Lin–Muranov–Yau** (path-homology foundations) supply the *persistence + stability* the KGW2026
  far rung would need. Park until a causal/temporal backbone is load-bearing.
- **TR / triangulation** — **OPEN.** No TDA-specific "corroboration of independent topological summaries" source
  is known; the grounding is general statistics (independence / conditional-independence / evidence combination).
  Worth a scout, but this is a place the program may be **ahead of the literature** — record it as such rather
  than forcing a cite. (See the corrected TR-1/TR-2: coincidence of *independent* signals is evidence; the failure
  is *dependence*, not coincidence.)

## Read-don't-acquire (held, worth the read)

- **DS2026** (above — the highest-ROI read; already in `ph`).
- **TKH2022** (`statistics`) — *Topological approximate Bayesian computation*: an **alternative inference route**
  to the barcode posterior (ABC on topological summaries instead of the BARS-curve pushforward). Read as a
  cousin/contrast to the L1–L4 path — it reaches a posterior over topology without the function-stability spine,
  which is a useful adversarial check on whether the BARS route's structure is load-bearing or merely convenient.

---

## Sequencing

Only two "soon" acquisitions — **Fasy 2014** (the frequentist-L3, nearest to already-needed) and the
**random-Mapper pair** (NV-2, as ThermoMapper comes online). Everything else is **read** (DS2026, TKH2022) or
**park**. Before acquiring any named paper, grep the corpus keys — the inventory was read at title level and a
paper may already be held under an unguessed key.

**Cross-refs:** [confidence-pushforward-lemmas.md](confidence-pushforward-lemmas.md) (the candidate families this
services) · [vocabulary.md](../../.discussion/issues/architecture-overhaul/vocabulary.md) (the warped-curve /
ordered-backbone note — a DS2026 cite belongs there) · [[project_thermomapper_architecture]].
