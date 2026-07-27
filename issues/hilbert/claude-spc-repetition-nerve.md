# SPC, PH, Mapper and Repetition — multiplicity vs correspondence

**Status:** analysis and forward sketch; discussion-grade, nothing scheduled. Fifth doc
of the arc — companions:
[claude-HOPE-paper-analysis.md](claude-HOPE-paper-analysis.md),
[claude-hilbert-synthesis.md](claude-hilbert-synthesis.md),
[claude-heat-semigroup-engines.md](claude-heat-semigroup-engines.md),
[claude-repeated-units-bundles.md](claude-repeated-units-bundles.md).

**Seed (Azriel):** could discovering intrinsic repeating structure in large graphs be
approached thermodynamically, via SPC's systematic exploration over temperature and
the hierarchy/dendrogram machinery already built — sensing a connection between PH,
SPC and Mapper.

---

## I. The division of labor (the load-bearing limit)

**SPC segments; it does not recognize.** A partition assigns nodes to groups. Two
copies of a template in different locations land in *different* clusters, and nothing
in the dendrogram states they are the same shape. The same holds for PH and Mapper.

| gives you | machinery |
|---|---|
| **multiplicity** — how many things, what size/shape, at what scale, how stable | SPC sweep, PH barcodes, Mapper nerve |
| **correspondence** — which node of copy A maps to which node of copy B | functional maps, GW couplings (doc 3 §III) |

So the answer to the seed is a split, not a yes/no: the thermodynamic/topological side
is a strong **detector and proposer**; the coupling side remains the only **matcher**.

**Upgrade to doc 4's ladder:** SPC replaces the seed-and-grow heuristic at rung 3.
The dendrogram over `T` is a principled multi-scale candidate pool with stability
already attached — strictly better than growing from a rarest-role seed.

## II. Thermodynamic signatures of repetition

**Degeneracy is the fingerprint.** `m` near-copies melt at the same temperature. The
signal is therefore not the location of the χ(T) peak but the **cluster-size
distribution at the transition**: a spike of `m` near-equal-sized clusters breaking
together. This is the thermal twin of doc 4 §III.1's spectral near-degeneracy screen —
two independent screens for one phenomenon, which is worth having.

**Extensivity of free energy — the stronger, fittable claim.** If the graph is `m`
weakly-coupled copies, the partition function approximately factorizes:

```
F(T) ≈ m · f_unit(T) + coupling
```

so free energy is **extensive in the number of units**, and `m` can be *fit* from the
free-energy curve rather than counted from clusters. The sweep machinery already
computes the ingredients. Caveat to test: the coupling term is exactly what
horizontal/overlap connections contribute, so the fit degrades as units stop being
weakly coupled — the residual is itself a coupling-strength estimate.

**Thermal autonomy** (twin of doc 3 §II's dynamical autonomy): does a unit melt at the
same `T` in isolation as in situ? Units that do are genuinely modular. Same
boundary-menu question, thermal currency.

## III. The dendrogram trick — repetition without the hard matching problem

If the data has repeated units, the merge tree contains `m` isomorphic subtrees. **Tree
isomorphism is linear-time** (AHU canonical form) — unlike graph isomorphism.

> Run SPC → canonically hash every subtree → look for repeated hashes.

Candidate repeated units fall out with no matching solver. Approximate variants: tree
edit distance, or hashing merge-height profiles rather than exact shape.

**Honest limit:** the merge tree discards cycle information, so identical subtrees are
*necessary-ish, not sufficient* — a screen, not a proof. That is the correct role for
it: cheap enough to run on everything, and its output feeds the correspondence ladder
directly.

## IV. The PH / SPC / Mapper connection — one construction, three faces

The sensed connection is structural, not vibes:

- The single-linkage dendrogram **is** the H₀ barcode of the Rips filtration. SPC is
  thermal single-linkage, so the SPC dendrogram is an H₀ barcode of a **thermal**
  rather than metric filtration.
- Mapper is the **nerve of a cover** with a clusterer in each cell.
- Therefore: **SPC over `T` generates the cover; the nerve of that cover is Mapper;
  and that nerve is exactly doc 4's base graph.** The bundle framing and the existing
  engine are the same construction approached from two sides.
- **PH re-enters as the stability layer**: multiscale Mapper's interleaving theory is
  what certifies that the tower of nerves over `T` is trustworthy rather than an
  artifact of the sweep.
- **Barcode multiplicity** — `m` bars with near-identical birth/death — is itself a
  repetition signature, consistent with §II's degeneracy reading.

## V. Combined pipeline (screens → proposals → correspondence)

1. Screens (cheap, run always): entropy-susceptibility peak vs plateau (doc 3),
   spectral near-degeneracy (doc 4), cluster-size degeneracy at transition + free-energy
   extensivity (§II), barcode multiplicity (§IV).
2. Proposals: SPC dendrogram subtree hashing (§III) — candidate units with scale and
   stability attached.
3. Cover/nerve: candidates as an overlapping cover → Mapper nerve = base graph.
4. Correspondence: rung 0 → 1 → 2 (doc 3 §III) on the proposed units.
5. Bundle: connection Laplacian on base + correspondences → defects/holonomy (doc 4 §II).

Doc 4's manifold-null warning applies unchanged and hardest at step 1: on point-cloud
kNN graphs, local repetition is the null, so screens must fire on **mesoscale**
structure against a **matched random-geometric** baseline.

## VI. The SPCX layer — thermal curves as inference objects

**Seed (Azriel):** the ideas above dovetail with the spirit of **SPCX** — not "run SPC
across T and find T_c," but advanced analytical treatment of thermal observable
curves: BARS estimation of the curve as a **joint posterior**, features extracted
*analytically* from the fitted object; groundwork already laid (adaptive BARS
scheduling — uniform first pass, iterative RJMCMC curve fitting; Lean lemmas and
protolemmas attached to this line). Plus the SPC × Mapper applications (ThermoMapper
proper: SPC-derived Mapper lenses; global SPC over the Mapper nerve) — discussed, not
yet implemented.

**The reframe this buys the whole arc:** every screen in docs 3–5 is a *curve
feature* — χ(T) peak location/width, `−dS/d log t` peak-vs-plateau, cluster-size
degeneracy at transition, plateau extents, extensivity slope. SPCX-as-curve-inference
upgrades each from a grid-read to a posterior quantity:

1. **Analytic derivatives with uncertainty.** The screens need derivatives and
   curvature of noisy MC-estimated curves; finite differences amplify exactly the
   noise MC produces. Splines differentiate exactly — the joint posterior hands every
   derivative-based screen its credible interval for free.
2. **The knot posterior is a transition detector.** RJMCMC knot placement concentrates
   where the curve has structure — knot density over T is a "where is the action"
   observable at zero extra cost, and it is what the adaptive scheduling refines on:
   uniform pass → knots cluster → resample there → refit. An active-learning loop
   over T.
3. **Extensivity becomes Bayesian model comparison.** `F(T) ≈ m·f_unit(T)` vs a free
   curve is a shared-shape-times-multiplier model; `m` gets a posterior, and the
   Bayes-factor machinery already in the K.You track applies. The §II fit stops being
   a heuristic.
4. **The thermal curve is a behavioral signature** (rung-0 descriptor). Repeated units
   have the *same* `f_unit(T)`; clustering per-subgraph thermal-response curves is the
   thermal twin of HKS-profile clustering. Node-level version: per-node melting
   profile / `T_melt(i)` as a field over nodes.
5. **Thermal autonomy becomes a two-curve comparison with a likelihood** (in-situ vs
   island melting curves), not an eyeball.

**Engine note — MCMC feeding MCMC.** SW sampling yields heteroscedastic noisy
observables per T; the robust-by-augmentation BARS design is built for exactly that
noise model, and SPC + BARS already sit in the planned shared mixing/diagnostics
family. The Lean harness has natural targets here of the boring-load-bearing kind
(spline-derivative exactness; validity conditions for the extensivity decomposition).

**ThermoMapper proper** (destination markers, not scheduled): (a) curve-feature
fields as Mapper lenses — `T_melt(i)`, local susceptibility-peak location — now
analytic-with-uncertainty rather than grid-read; (b) global SPC over the Mapper
nerve; (c) §IV's SPC-cover → nerve. The repetition program and ThermoMapper proper
consume the **same nerve**.

## Open edges

- Does the free-energy extensivity fit survive realistic coupling (overlap +
  horizontal edges), and can the residual be read as coupling strength?
- Subtree-hash sensitivity: how much noise/overlap before hashes stop colliding —
  and is a height-profile hash materially more robust?
- Cluster-size-degeneracy screen vs spectral-degeneracy screen: do they fire on the
  same cases, and is either strictly stronger?
- Thermal autonomy: which boundary-menu entry corresponds to "in isolation" for a
  Potts subsystem?
- Does the SPC-cover → Mapper-nerve construction need the interleaving guarantees
  before it is usable, or is that a later rigor pass?
- Heteroscedastic MC-noise propagation into the BARS likelihood: is the augmentation
  scheme sufficient as-is, or does per-T error need explicit modeling?
- Does the knot posterior actually concentrate at transitions in practice — testable
  cheaply on synthetic m-copy fixtures before trusting it as a detector.
- `m`-posterior identifiability as coupling strengthens: where does the extensivity
  model comparison stop being able to distinguish m from m±1?
- Which curve-feature extraction claims are Lean-lemma-ready now vs protolemma-stage.
