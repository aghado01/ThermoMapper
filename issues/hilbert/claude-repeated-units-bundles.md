# Repeated Units — translation symmetry, bundles, and the discovery problem

**Status:** analysis and forward sketch; discussion-grade, nothing scheduled. Fourth
doc of the arc — companions:
[claude-HOPE-paper-analysis.md](claude-HOPE-paper-analysis.md) (first read),
[claude-hilbert-synthesis.md](claude-hilbert-synthesis.md) (HOPE × NLD × stack),
[claude-heat-semigroup-engines.md](claude-heat-semigroup-engines.md) (affordability,
subgraph reframing, correspondence ladder).

**Seed (Azriel):** not only scale-invariant self-similarity ("fractillitude"), but
*flatly repeated subgraph structure at a roughly fixed scale* — the receptive-field
case: each unit has similar internal connectivity, processes a different part of the
visual field, with **partial overlap** and **horizontal connections**. Two asks: (1)
**test** a graph for this kind of structure, (2) **discover** the repeated unit
patterns embedded in the larger graph, so correspondence-based methods can work with
them.

---

## I. A different symmetry — and the screen already separates them

| | fractillitude | repeated units |
|---|---|---|
| symmetry | **dilation** — object into itself under rescaling | **translation / permutation (deck)** — one template, many copies |
| index | scale | position |
| machinery | RG, coarsening ladders, multiscale GW | bundles, quotients/equitable partitions, motif alignment |
| screen signature | `−dS/d log t` **flat plateau** (no preferred scale) | `−dS/d log t` **isolated peak** (a characteristic scale) |

The same entropy-susceptibility curve from the engines doc distinguishes them with no
new machinery: repeated units of roughly fixed size **have** a preferred scale, so
they announce as a bump where the fractal case announces as a plateau. One Chebyshev
basis buys the first test. The two symmetries also compose — a hierarchy of repeated
units (retina → columns → areas) shows a *ladder of peaks*, which is a third,
distinguishable reading of the same curve.

## II. The right object: a bundle over a base graph

"Repeated units with partial overlap and horizontal connections" is not a bag of
isomorphic subgraphs. It is:

- **base** — the arrangement of units (retinotopic position);
- **fiber** — the canonical unit's internal structure;
- **connection** — the rule identifying one unit's internal coordinates with its
  neighbor's, carried by the overlap and the horizontal edges.

Because the units **overlap rather than partition**, the combinatorial home of the
base is a **nerve of a cover**: Mapper is not an analogy here, it is the data
structure. And the operator for a graph with fiber structure plus an identification
rule is the **connection Laplacian** (Singer–Wu vector diffusion maps) or its
**cellular-sheaf** generalization — both already on the field-ladder aspiration list
and adjacent to the magnetic/sheaf persistent-Laplacian track. This is the natural
convergence of the two: repeated-unit discovery *produces* the sheaf that the
higher-degree Laplacian machinery *consumes*.

**Pinwheels are holonomy — the payoff that makes the framing more than elegant.**
Traverse a loop of overlapping orientation columns and preference rotates by 180°:
the orientation map is a section of a bundle over the cortical surface with
half-integer topological defects (established neuroscience, not extrapolation). In
bundle language, *"does the correspondence compose back to the identity around a
cycle?"* is exactly what the connection Laplacian's harmonic space measures. So
discovery and consistency are not two projects — the **defects are the interesting
structure**, and they are an eigen-computation once correspondences exist.

**A sharp testable question about horizontal connections.** Cortex wires like-to-like
across columns; in bundle terms, the lateral connection is *compatible* with the
fiber identification. That is directly checkable: compare the correspondence induced
by **overlap** against the one induced by **horizontal edges** and ask whether they
agree. Disagreement is a measurable, falsifiable claim about the circuit — not a
metaphor.

## III. Discovery ladder (cheap → committal)

Each rung reuses the prior rung's computed descriptors; nothing is recomputed.

1. **Spectral-density screen** (KPM, no eigensolve). `m` near-copies of a unit produce
   ~`m`-fold **near-degeneracies**: anomalous multiplicity clusters and spikes in the
   density of states are the fingerprint. Cheap, global, checked against a null (§IV).
   Pairs with the `−dS/d log t` peak from §I.
2. **Role discovery, not unit discovery.** Cluster nodes by **HKS profile over a
   t-window**. The framing worth stating cleanly: **HKS is a soft, noise-tolerant
   Weisfeiler–Leman coloring** — WL refinement computes the coarsest **equitable
   partition**, whose **quotient is the template**. This is the Hilbert-side statement
   of "repeated units": an equitable partition ⇔ an invariant subspace of the
   Laplacian. Established theory and a validation regime come via **cluster
   synchronization** (Pecora–Sorrentino, Schaub) — the same partitions govern which
   node groups can synchronize, which is an independent handle on correctness.
3. **Unit localization.** Seed on the **rarest role** (highest information per seed),
   grow until the role complement closes — one candidate unit per seed, with the
   boundary menu (island / Dirichlet / Schur, engines doc §II) declaring how the
   environment is treated.
4. **Alignment.** Feed candidates to the correspondence ladder (engines doc §III):
   rung 0 invariant descriptors → rung 1 functional maps → rung 2 entropic GW. Cheap
   at this point precisely because the screen already computed the descriptors.
5. **Bundle assembly.** Nerve of the discovered (overlapping) units = base; pairwise
   correspondences = connection; then connection-Laplacian harmonics = defects and
   global consistency (§II).

## IV. Two honest warnings

**1. On a manifold point-cloud kNN graph, local repetition is the NULL, not the
signal.** Locally-Euclidean *means* every neighborhood looks alike; a naive
repeated-unit detector fires everywhere and reports nothing. Consequences:

- the signal must be **mesoscale structure not explained by local homogeneity**;
- the null must be a **matched random geometric graph**, *not* a degree-preserving
  rewiring — rewiring nulls make trivially-homogeneous geometry look anomalous and
  will manufacture false positives wholesale;
- the connectome case is genuinely different: repetition in wiring is non-trivial
  because the null there is not "smooth manifold."

This is the translation-symmetry analogue of the false-fractal hazard, and it bites
harder because point clouds are the target data.

**2. Approximate symmetry detection is less mature than the fractal-screen side.**
Exact automorphism (nauty) is fast but exact-only; frequent-subgraph mining is brittle
and combinatorial; the descriptor-clustering route recommended above is a **heuristic
that presumes roles are distinguishable**. Build for **soft, partial, approximate**
units — which is also what the biology argues: the columnar story itself is contested
(Horton & Adams; salt-and-pepper organization in rodents), so crisp units are not a
safe presumption even in the motivating case.

## V. Nearer application than visual cortex

In a **trajectory graph**, repeated subgraph structure is **recurrence** — the same
state-space region revisited. Motif discovery wearing SIFTS clothing, on data already
in hand ([[project_spcx_telos]]'s NHP warm-up), where the null question is better
posed and the ground truth is closer.

## Open edges

- Sharpen the mesoscale null: what exactly does "not explained by local homogeneity"
  mean operationally for a matched random-geometric baseline?
- Which discovery rung is the honest stopping point for noisy weighted graphs — is
  role-clustering + rung-1 alignment enough without GW?
- Degeneracy screen sensitivity: how much noise/overlap before near-degeneracy
  clusters wash out?
- Overlap-vs-horizontal correspondence agreement (§II) — formalize as a statistic with
  a null.
- Relationship to the equitable-partition literature's exact algorithms: is there a
  soft relaxation with guarantees, or is descriptor clustering the practical ceiling?
