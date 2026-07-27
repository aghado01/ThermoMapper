# The Arc Seen Whole — filtrations, readouts, quotients

**Status:** capstone synthesis; discussion-grade, nothing scheduled. Sixth and closing
doc of the arc — companions:
[claude-HOPE-paper-analysis.md](claude-HOPE-paper-analysis.md),
[claude-hilbert-synthesis.md](claude-hilbert-synthesis.md),
[claude-heat-semigroup-engines.md](claude-heat-semigroup-engines.md),
[claude-repeated-units-bundles.md](claude-repeated-units-bundles.md),
[claude-spc-repetition-nerve.md](claude-spc-repetition-nerve.md).

**Seeds (Azriel):** circle back on how the whole line connects — repeated structure in
graphs, HOPE, NLD, and expanded-SPC + BARS related to Mapper theoretically; the
Mapper / PH / SPC concepts are still being untangled, with a sense of greenfield
pursued implicitly for a while.

Three levels: **generate** a filtration, **read** its observable curves, **quotient**
by what repeats. §II is the untangling.

---

## I. Generate — three filtration axes, one dataset

Every engine in the arc produces a one-parameter operator or measure family. Choosing
the axis is choosing the physics:

| axis | family | physics | native to |
|---|---|---|---|
| **heat** `t` | `{e^{−tL}}` | free / Gaussian field | NLD, diffusion distances |
| **thermal** `T` | `{Gibbs_T}` (Potts) | interacting field | SPC |
| **metric** `ε` | `{Rips_ε}` | none — combinatorial | PH |
| **capacity** | progressive compression | (component domain) | HOPE |

The heat↔thermal relation is doc 3 §II's bridge: `e^{−tL}` is the Gibbs propagator of
the *free* field, Potts correlations the interacting cousin — diffusion is the
Gaussian sector of the theory the stack already runs. And the H₀ content is
**literally shared**: the single-linkage dendrogram is the metric H₀ barcode; the SPC
dendrogram is the *thermal* H₀ barcode. SPC has been a persistence theory in
disguise. NLD's `∫dt` is the integrate-the-whole-filtration instinct on the heat axis.

## II. Untangling Mapper / Reeb / merge tree / PH — and where zigzag enters

**The distinction that organizes everything: sublevel sets vs level sets.**

| construction | filters by | tracks | output | sees loops? |
|---|---|---|---|---|
| **merge tree** | sublevel `f⁻¹(−∞,t]` | components of sublevel sets | tree | **no** |
| **H₀ barcode** | sublevel | births/deaths (elder rule) | multiset of bars | no |
| **Reeb graph** | level `f⁻¹(t)` | components of level sets | graph | **yes** |
| **Mapper** | interval cover `f⁻¹(U_i)` | clusters per cell → nerve | graph/complex | yes |
| **levelset zigzag** | same interval cover | H_k across cells/intersections | barcode | yes |

**Established relations** (not my synthesis):

- **Dendrogram = merge tree = ultrametric** (Carlsson–Mémoli; single linkage is the
  stable one). So **the SPC dendrogram is a thermal merge tree**, and the whole
  hierarchy machinery is merge-tree machinery.
- **Merge tree → H₀ barcode is lossy**: the tree strictly refines the barcode.
  Likewise **Reeb graph → merge tree** (take the merge tree of the Reeb graph under
  its induced function). A chain of quotients: Reeb ⟶ merge ⟶ barcode.
- **Mapper is the statistical Reeb graph**: it replaces level sets with interval-cover
  preimages and connected components with a clusterer, then takes the nerve.
  Convergence to the Reeb graph and stability under cover refinement are the
  Munch–Wang / Dey–Mémoli–Wang (multiscale Mapper) results.
- **Levelset zigzag** (Carlsson–de Silva–Morozov) runs the alternating diagram
  `f⁻¹(I₁) ← f⁻¹(I₁∩I₂) → f⁻¹(I₂) ← …` over an interval cover of the range and
  decomposes it into bars.

**The consequence worth pinning — Mapper and levelset zigzag are two readouts of the
same interval cover.** Mapper takes the **nerve** (an object); zigzag takes the
**persistence module** (its decomposition into bars). Same cover, same function, two
faces. Since the repo already has a zigzag engine, **the machinery for the
Reeb-theoretic content of Mapper is already built** — the implicit pursuit made
explicit. (Caveat: the barcode is the decomposition; the Reeb graph as an object —
in the cosheaf view of de Silva–Munch–Patel — is strictly richer than its bars.
Mapper approximates the object; zigzag decomposes it.)

**The second consequence — SPC alone is structurally blind to loops.** A merge tree
is a tree by construction; no amount of temperature sweeping makes a dendrogram show
genus. Loops enter only via the level-set side (Reeb/Mapper) or via H₁ (PH). This is
not a defect of SPC — it is the merge-tree/Reeb-graph distinction — and it explains
*why the stack needs Mapper as a peer and not a wrapper*.

**Third — the loops are exactly the holonomy sites.** Doc 4's bundle framing measures
twist along cycles of the base graph. The base is a nerve; its cycles are the Reeb
graph's `b₁` generators. So **Reeb/Mapper supplies the loops; the connection Laplacian
supplies the twist along them** (pinwheel-grade structure). PH's H₁ certifies which
loops are real rather than cover artifacts.

### What ThermoMapper's compositions actually are

The three SPC × Mapper compositions are mathematically *different objects* — worth
naming separately rather than treating as one feature:

1. **SPC observable field as the lens** (`T_melt(i)`, local susceptibility-peak
   location): Mapper then approximates the **Reeb graph of a thermodynamic
   function**. A genuinely unusual lens family — the standard menu is
   eccentricity/density/PCA.
2. **SPC as the clusterer inside cells**: this changes what is being approximated.
   Mapper's convergence theory assumes the clusterer stands in for *connected
   components*; a thermal/density clusterer instead yields the Reeb graph of the lens
   **restricted to the density-supported region** (the same caveat that applies to
   DBSCAN-in-Mapper, standard practice but rarely stated). Declare which object you
   are estimating — this is a faithfulness question, not a tuning question.
3. **Global SPC over the nerve**: clustering the nerve graph itself — a second-order
   operation on the quotient, not a Mapper variant.

## III. Read — SPCX + BARS is an axis-agnostic readout, and it is persistence one level up

Nothing in the BARS layer cares whether the x-axis is `T`, `t`, or `ε`: the Lean slate
theorems are about spline curves. So "expanded SPC" is really a **general theory of
inferring features of filtration-indexed observable curves from noisy samples** —
applies verbatim to `χ(T)`, `tr(e^{−tL})`, `−dS/d log t`, `Q_S(t)`.

**The identity worth pinning: peak prominence is 0-dimensional persistence of the
curve.** Under the superlevel-set filtration of `f(T)`, maxima are births, cols are
deaths, and prominence *is* the persistence. The θ-thresholded `SignificantPeakCount`
is therefore a persistence cutoff. The multi-peak readout is not persistence-*like*;
it **is** PH applied one level up — persistence of the *observable of* the filtration
rather than of the data.

Consequences for the Lean slate (docs 5 §VI):

- **MP-3**'s "mass split between k and k+1 near a bifurcation is correct ambiguity"
  is barcode instability read correctly — the mature statement is the PH stability
  theorem.
- **MP-4**'s intensity `λ(T)` is kin to **expected persistence diagrams / persistence
  intensity functions** (Chazal–Divol): theirs over (birth, death), yours the
  location marginal with span marks. The slate has been formalizing baby persistence
  theory without the vocabulary — which says both that it is well-founded and where
  the mature theorems live.
- Schedule-neutral interleaved-uniform sampling (doc 5 §VI) is what makes
  **multi-feature** curves first-class; hunting schedules are single-`T_c`-shaped and
  would distort exactly the multi-transition curves repetition produces.

## IV. Quotient — repetition is symmetry along the two directions of one field

Everything in the arc is a field `f(x, s)` over **(data × scale)** — the function
algebra's domain grammar built for this sentence:

- **Fractillitude = invariance along the scale direction** (dilation on `s`):
  power-law windows, flat `−dS/d log t`.
- **Repeated units = invariance along the data direction** (translation/deck symmetry
  on `x`): spectral near-degeneracies, `F ≈ m·f_unit(T)`, co-melting.
- **Hierarchies of units** (the mixed case): a *ladder* of peaks — discrete scale
  invariance, whose known signature is log-periodicity.

Detection happens at level III: symmetries of the data manifest as **degeneracies and
self-similarities of the curves**. Construction happens here: **Mapper is the quotient
functor of the whole story** — a nerve is precisely the combinatorial quotient by
"same cell", and when the cells are discovered units, the nerve is the bundle's base.

**The pipeline, stated once:** SPC proposes covers (thermal, multiscale) → BARS
certifies which scales are stable (plateaus, peak posteriors) → Mapper quotients
(nerve = base) → couplings glue fibers (connection) → connection-Laplacian harmonics
measure the twist (holonomy). Each stage carries a declared fidelity per the
(gauge, measure, metric) audit; doc 4's manifold-null warning gates the entry.

## V. Closure — the original question

The thread opened on recalling Hilbert spaces and operators for attention and
activations. The frame answers it: **model internals are analyzed exactly like data**
— operator families under declared measures (HOPE's contribution), curve readouts,
symmetry quotients. Repeated computational units in transformers (circuit motifs
recurring across layers and positions) are the **translation-symmetric case in the
component domain**: doc 4's program applied to HOPE's objects. HOPE's own merging is
quotienting by approximate neuron-equality — the parent neuron is a fiber template,
DEFT's core the quotient's stable part — and superposition has a bundle reading
(features as sections not aligned with the neuron basis). The interpretability motif
hunt, which currently runs on bespoke probes, is a special case of the screens →
couplings → holonomy toolkit this arc assembled.

## VI. Greenfield register (labelled by maturity)

| # | item | status |
|---|---|---|
| 1 | **Zigzag ↔ Mapper on one cover** — use the existing zigzag engine for the Reeb content of Mapper | established theory (CdSM), unexploited here |
| 2 | **Thermal lenses** — Reeb graph of `T_melt(i)` / susceptibility-peak fields | synthesis; buildable on existing parts |
| 3 | **Declare which object SPC-in-cells estimates** (density-supported Reeb graph) | correctness note, cheap, do first |
| 4 | **(T, cover-resolution) bifiltration** — ThermoMapper is natively 2-parameter | established frontier: no complete discrete invariant; rank invariant / fibered barcodes / RIVET-style tooling |
| 5 | **Soft nerve from SPC soft memberships** — weighted/fuzzy Mapper from `Groups` | thin literature — real greenfield |
| 6 | **Loops as holonomy sites** (Reeb `b₁` = base cycles = defect carriers) | my synthesis; testable via doc 4 §II |
| 7 | **Merge-tree comparison for repeated units** — interleaving distance / ultrametric GW between subtree dendrograms | established but computationally hard; approximations exist; feeds doc 5 §III hashing |
| 8 | **Reeb *space*** for multivariate thermal lenses | established (Edelsbrunner–Harer–Patel); pulls back to item 4 |

## Open edges

- Item 3 first: does swapping SPC into Mapper cells break the convergence guarantees
  in a way that matters at realistic cover resolutions, or only in principle?
- Does the zigzag engine's existing interval/cover machinery accept a Mapper cover
  unchanged, or is there an impedance mismatch in how covers are specified?
- Two-parameter persistence: is the honest near-term move a *fibered* readout (slices
  along fixed `T`, barcode per slice — which the BARS layer then reads as curves over
  the slice parameter) rather than a genuine 2-parameter invariant?
- Merge-tree interleaving is NP-hard exactly; is the subtree-hash screen (doc 5 §III)
  plus a cheap approximation enough for the repetition program?
- Prominence-as-persistence: does adopting the PH stability theorem wholesale
  supersede any of MP-2/MP-3, or are the slate's elementary statements still the
  right granularity for the Lean harness?
- Soft nerve: what is the right notion of "intersection" when cell membership is
  fuzzy — and does the nerve lemma survive it in any usable form?
