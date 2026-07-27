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
