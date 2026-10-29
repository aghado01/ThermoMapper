# Spin-glass threads: frustration, overlaps, ultrametricity, and the regression bridge

Discussion digest, 2026-10-28. Seeded by Domany's later papers (local copies in
`codex-scientiae/supellex-private/PDFs/Domany`). Status: ideation. Nothing here is
decided; the "Candidate work" list is what could be briefed.

## Sources read

| File | Paper | Used here for |
|---|---|---|
| `PhysRevLett.92.157201` | Hed, Young, Domany 2004, *Lack of ultrametricity in the low-T phase of 3D Ising spin glasses* | ultrametricity test conditioned on a clustering |
| `PhysRevB.76.132408` | Hed, Domany 2007, *Nontrivial link overlap distribution in 3D Ising spin glasses* | link overlap, cross-cluster conditioning |
| `j.physa.2013.07.036` | Farhi, Hed, …, Domany 2013, *Temperature Integration* | ln Z along a tempering ladder; ΔF between systems with equal Z(0); the energy cutoff |
| `bioinformatics_19_9_1153` | Getz, Domany 2003, *Coupled two-way clustering server* | SPC stability index in biclustering (not pursued) |
| `PhysRevD.98.063511` | Leclercq 2018, *Bayesian optimization for likelihood-free cosmological inference* (not Domany) | GP-surrogate acquisition under a simulation budget |
| `gauntlet/bayesian/1306.6430v2` | Bissiri, Holmes, Walker 2016, *A general framework for updating belief distributions* (*JRSS B* 78, not Domany; read 2026-10-29) | SPC as a general-Bayes update over partitions; calibrating the temperature (Thread 7) |

**The shared move in both spin-glass papers.** A global statistic over an
equilibrium ensemble is dominated by its least informative mode (pairs inside the
dominant state). They cluster the sampled *states* first (average linkage on
d = (1 − q)/2), then compute the statistic only across the tree's top split. The
clustering conditions the diagnostic.

## Thread 1: ultrametricity as a faithfulness test for any dendrogram

- Single linkage outputs the subdominant ultrametric of its dissimilarity, so the
  cumulative field's hierarchy (and any single-linkage tree) is a tree by
  construction, whether or not the data are hierarchical at that scale.
- Hed–Young–Domany test the *underlying* distances against the tree:
  - **Triplet index.** Take α ∈ C1a, β ∈ C1b, γ ∈ C2 (across the top splits,
    never at random: random triplets mostly fall in one state and are
    equilateral, i.e. uninformative). Sort the three distances and set
    K = (largest − middle) / smallest ∈ [0, 1]. Ultrametric ⇒ K → 0. They discard
    collinear triangles (K = 1) as a finite-size effect.
  - **Block uniformity.** Under ultrametricity the C1 × C2 distance block depends
    only on the common ancestor, so it is uniform. Statistic: std of the block,
    normalized by the Bernoulli bound √(m(1 − m)).
  - **Result.** SK (infinite range) becomes ultrametric as N grows; the 3D
    short-range model does not (K peaks near 0.8 and stays).
- **For ThermoMapper.** A per-split faithfulness score for any `Dendrogram`
  against an original or graph-geodesic distance: sharper and more local than a
  cophenetic correlation. It is the hierarchy-level counterpart of the FKG
  nesting gate in `HierarchyEom`. Fits the faithfulness axis.
- **Open placement question.** `Clustering.Evaluation` (needs a `Dendrogram`
  reference) or `Clustering.Dendrograms` (where the trees live).

## Thread 2: link overlap = edge-level partition comparison

- Site overlap q counts flipped *volume*; link overlap
  Q = N_b⁻¹ Σ_⟨ij⟩ s_i s_j t_i t_j counts agreement on *bonds*, i.e. domain-wall
  area. In SK, Q is a function of q; in short-range systems they differ (surface
  vs volume).
- **Exact identity.** For two partitions A, B, the Rand index is the fraction of
  *all* node pairs on which they agree (same/same or different/different). That is
  the Potts (δ-form) link overlap on the complete graph, the mean-field case.
- **Restricted to the coupling graph's edges** (optionally J-weighted), the same
  count is the short-range link overlap. It measures whether two partitions place
  their boundaries in the same places, which ARI averages away.
- **Uses.** Comparing cumulative-field vs Swendsen–Wang partitions; comparing the
  Mutual/Inclusive symmetrization members; and, under Thread 4, replica overlaps.
- **Inventory.** `Clustering.Evaluation` has `AdjustedRandIndex` but no plain Rand
  and no graph-restricted variant. The graph variant compares two partitions on a
  graph, a shape neither `IExternalClusterEvaluator` (two labelings) nor
  `IGraphPartitionEvaluator` (one labeling plus a field) covers exactly.

## Thread 3: clustering the states, not the data

- They cluster spin configurations to recover the hierarchy of Gibbs (pure)
  states. The ThermoMapper analogue is clustering retained Swendsen–Wang draws by
  a partition distance (VI, or Thread 2's link overlap) to expose the structure of
  the sampler's own posterior over partitions near a transition. This touches the
  parked consensus and soft-membership work. It needs per-sweep partitions: SW
  retention S5 (`issues/spc/sw-retention-brief.md`).

## Thread 4: thermodynamic clustering with frustration

Where it stops being an analogy. Sketch only.

- **Signed couplings.** J_ij > 0 attract, J_ij < 0 repel. Sources:
  - null-model subtraction J = A − γP (every dataset becomes frustrated, dense, SK-like);
  - explicit repulsion of non-neighbors;
  - signed similarity data (correlation matrices, the f-SPC/Yelibi domain);
  - cannot-link constraints.
- **Ground state.** With H = −Σ J_ij δ(σ_i, σ_j) and q at least the number of
  clusters, the T → 0 problem is correlation clustering: minimize positive edges
  cut plus negative edges kept. The number of clusters emerges and is not an input.
- **What frustration means for clustering.** For Potts with many states it is
  cycles with *exactly one* negative edge (weak balance), not cycles with an odd
  number of negative edges (Ising/strong balance): a triangle of three repulsions is
  satisfiable with three clusters. The minimum number of such violations (a
  frustration index) is a data-level clusterability diagnostic, available before
  any sampling.
- **Sampling.** The FK/Swendsen–Wang representation needs J ≥ 0 (a negative bond
  probability is meaningless). Options:
  - antiferromagnetic-Potts cluster moves (WSK);
  - two-replica cluster moves (Houdayer);
  - local updates inside parallel tempering, the spin-glass workhorse. This is
    where `Maths.Samplers` already helps (Thread 5).
- **Observables.** Magnetization and FK susceptibility lose their meaning. Use
  replica overlaps at each T from two or more independent chains. Thread 2's
  link overlap on co-membership is label-free, so P(Q) per temperature is
  well-defined without aligning Potts labels.
  - A single sharp peak means one clustering.
  - A broad or multimodal P(Q) means the data support several incompatible
    clusterings: an ambiguity readout ferromagnetic SPC cannot express.
- **Output object.** Without FKG the temperature stack need not nest, so the
  existing non-nested resolvers apply directly: zigzag persistence (`TDA.Ph`) and
  `LineagePersistence`. Thread 3 plus Thread 1 decide whether the competing
  clusterings form a tree at all.
- **Range matters (speculative but grounded in the 2004 result).** Dense or
  null-model (SK-like) frustrated clustering may yield an ultrametric hierarchy of
  alternative clusterings; sparse kNN frustration (short-range) may not. If so,
  graph topology decides whether a dendrogram is the right output.
- **Out of scope for the cumulative field.** Its ladder needs J ≥ 0 to be
  monotone; signed couplings break the filtration.

## Thread 5: the regression bridge (exact, not analogy)

- The Potts Boltzmann family π_β(σ) ∝ exp(−βH(σ)) over the qᴺ configurations is a
  power-posterior path: uniform prior at β = 0 (ln Z = N ln q) and "likelihood"
  exp(−H). Thermodynamic integration (Farhi et al.; Friel and Pettitt) and
  stepping-stone apply unchanged.
- `Maths.Samplers` already has the generic machinery, built for BARS/BASS:
  - `ParallelTempering<TState>` with non-reversible DEO swaps (Syed et al. 2022);
  - `TemperingLadder` with rung adaptation from swap-rejection rates;
  - `LadderEvidence` (stepping-stone, and path sampling with the Friel–Hurn–Wyse
    correction);
  - SMC.
- Three planks:
  1. **BARS reads f(T)** (exists).
  2. **The SPC temperature grid is a tempering ladder without swaps.** Adding
     replica exchange makes it parallel tempering, which frustrated SPC needs and
     which also speeds ferromagnetic SPC near transitions.
  3. **Free energy by temperature integration** (Farhi et al. 2013).
     ln Z(β) − N ln q = −∫₀^β ⟨H⟩ dβ′, with the SPC grid's temperatures as the
     intermediates. Z(0) = qᴺ for every coupling graph on the same nodes and the
     Potts energy is bounded, so no cutoff is needed and the result is the
     absolute free energy. The SW cycle records H = −Σ J δ once per undirected
     edge, the Hamiltonian it samples.
     - A trapezoid over the running ⟨E⟩(T) sums works from existing output.
       `LadderEvidence` (stepping-stone, path sampling with a Monte Carlo error)
       needs per-sweep energies: SW retention S1. The Potts slope is bounded
       (⟨H⟩ = −Σ J/q at β = 0), so path sampling is safe here, unlike BASS's
       integrated form, whose slope at b = 0 is infinite under Zellner–Siow.
     - **Grid closure: a reference rung**, as BASS has (`BassEvidence`: an exact
       prior draw every scan at b = 0). The regression side already runs this
       shape of evidence, so the design flows from there to SPC. The Potts β = 0
       member is iid uniform spins, drawn exactly, so stepping-stone from it up
       through a few cheap hot rungs (SW forms few bonds there and mixes in a
       sweep) closes the integral with no approximation.
       - Under uniform spins the edge agreements are pairwise independent, so
         Var(H) = Σₑ Jₑ²(1/q)(1 − 1/q). A step from β = 0 to β stays well
         conditioned while β²·Var(H) is of order 1: the analogue of the hot end
         BASS auto-picks at 1/N.
     - **Forest term, as a check.** ln Z ≈ N ln q + Σₑ ln(1 + (e^{βJₑ} − 1)/q) is
       exact on a forest; cycles add corrections at order (e^{βJ} − 1)^girth. It
       checks the reference-rung closure on trees and approximates it elsewhere.
     - **Equilibrium.** The method assumes converged parallel tempering. SPC's
       fixed-T chains have no swaps, so ⟨E⟩ can be biased near transitions. R-hat
       across replicas (S1) marks the rungs to distrust; plank 2 removes the cause.
- **Temperature placement.** The rejection-rate ladder concentrates rungs where
  the energy distribution changes fastest (the communication barrier tracks energy
  dispersion ∝ √C_v). That is SPC's transitions, which ties tempering theory to the
  SPC acquisition/bracket question.
- **Caution.** ΔF between two coupling graphs is not model evidence for the graph
  construction: there is no likelihood of the data given the couplings. It compares
  thermodynamic potentials, not fits.
- **Two systems.** Two Potts systems on the same nodes and q share Z(0), so ΔF(T)
  between them is two temperature integrals and needs no λ-path. Between coupling
  graphs it is dominated by the total coupling scale unless the couplings are
  normalized (the parked coupling-normalization family). Two pairings may mean
  more; both are speculative:
  - **Signed vs unsigned couplings** on the same nodes (Thread 4): the free-energy
    cost of the repulsions at each T.
  - **Constrained vs free.** A must-link or cannot-link constraint is an infinite
    energy term; the paper's cutoff makes it finite so both systems share Z at
    high T. ΔF = −T ln P(the constraint holds at T): the cost of forcing a split or
    merge. It reaches probabilities too small to count from samples, as at low T,
    where direct counting fails.

## Thread 6: acquisition under a simulation budget (Leclercq)

- BOLFI-style: GP surrogate of the expensive signal, acquire where expected
  posterior uncertainty drops most. Maps onto choosing Swendsen–Wang temperatures
  (and BASS stopping/readout) under a fixed sweep budget. Complements Thread 5's
  ladder rule: tempering decides rung spacing for mixing; acquisition decides where
  the *readout* needs resolution.

## Thread 7: SPC as a general-Bayes update over partitions (Bissiri et al.)

- **The reading.** Bissiri, Holmes and Walker's model-free clustering (2016, §5.2)
  is P(S | x) ∝ π(S)·exp{−w·l(S, x)} over partitions S, with a loss l and a
  weight w. Its readout is the co-clustering matrix, P(i and j in one cluster).
  - The Potts family at β is that update, with the loss H and w = β. SPC's pair
    correlation (`EdgeCorrelation`) is their co-clustering probability.
  - Thread 5's power-posterior path is their annealing (§3.1).
  - Under their framework each SPC temperature is a valid belief update about the
    partition, not only a physical ensemble, and that needs no likelihood of the
    data given the couplings (Thread 5's caution).
- **SPC already has a prior on partitions, implicitly.** In the Fortuin–Kasteleyn
  random-cluster form, Z ∝ Σ_ω Π_{e∈ω} p_e Π_{e∉ω} (1 − p_e) · q^{k(ω)}, with
  p_e = 1 − e^{−βJ_e} and k(ω) the bond configuration's cluster count.
  - The q^{k(ω)} factor is a prior weight on the number of SW clusters, so q is a
    prior choice as well as a resolution knob.
  - Bissiri et al. state theirs explicitly: a Poisson prior on the number of
    groups.
- **What it adds:**
  1. **A temperature rule for discrete structures** (§5.2). Choose w so that a
     partition whose loss reduction from the single cluster sits at the α tail of
     random allocation's null gets posterior odds R against the single cluster.
     - It fixes T by what reduction should count as a real clustering: an
       evidence-based alternative to the susceptibility peak and to the
       fuzzy-entropy plateau the SW→f-SPC handoff reads.
     - *The printed rule is inverted.* The last bullet of §5.2 reads
       w = −log(R)·l(S, x)·f*(k − 1)/(n − k). From the odds equation it is
       w = −ln R·(1 + c)/(c·l(S₀, x)), with c = f*_α(k − 1)/(n − k) and S₀ the
       single cluster. Only this form reproduces the paper's own numbers: w falls
       from 0.0036 to 0.0012 as α rises from 0.99 to 0.999.
     - *The null does not carry over as it stands.* The paper's F reference holds
       for a within-cluster sum of squares. SPC's loss is H on a coupling graph, so
       the null of H's reduction under random allocation would come from
       permutations, empirically. Its cost is unmeasured.
  2. **An explicit prior on partitions,** beside or instead of q's implicit one.
  3. **The joint allocation's uncertainty as the primary readout:** the
     co-clustering matrix at a calibrated temperature, read as a posterior
     probability rather than a thermodynamic observable. It is the same number
     `EdgeCorrelation` already gives; what changes is how the temperature is
     chosen and what the readout is taken to mean.
- **Caution.** A weight calibrated on one coupling graph does not transfer to
  another, since the loss's scale moves with the couplings (the parked
  coupling-normalization family). Comparisons between graphs stay thermodynamic
  (Thread 5).

## Candidate work (unbriefed)

1. **Rand index + graph-restricted link overlap** in `Clustering.Evaluation`
   (Thread 2). **Done 2026-10-28:** `RandIndex` is an external quality index in
   `Clustering.Evaluation`. `LinkOverlap` is an observable, not an evaluator: it
   compares any two node labelings, so it lives in `Graphs.Observables`
   (`LinkOverlap.Compute(graph, a, b, weighted)`, symmetric, negative labels =
   unplaced). Replica overlaps (Thread 4) reuse it directly.
2. **Ultrametricity test** (triplet K distribution across top splits, block
   uniformity) as a dendrogram faithfulness diagnostic (Thread 1). Moderate;
   placement decision.
3. **Free-energy curve by temperature integration** (Thread 5). A trapezoid over
   the existing ⟨E⟩(T) sums needs no new sampling; the `LadderEvidence` route
   (Monte Carlo error, stepping-stone) depends on SW retention S1. Either needs the
   reference-rung closure at β = 0. Small once S1 lands.
4. **SW as a replica inside `ParallelTempering`** (Thread 5, prerequisite for
   Thread 4). Moderate; needs a Potts state type fitting the ladder's replica
   contract.
5. **Frustrated SPC prototype** (Thread 4). Research-scale: signed couplings from
   one source (null-model subtraction is the simplest), local updates in PT,
   replica link-overlap observables, zigzag/lineage output.
6. **Temperature calibration by posterior odds** (Thread 7). Fix T by the odds a
   significant clustering should get against the single cluster, with a
   permutation null of H's reduction. Compare it with the susceptibility peak and
   the entropy plateau on the Domany battery. Research-scale; first measure the
   null's cost.

## Reading leads (verify before citing)

- Bansal, Blum, Chawla, *Correlation Clustering* (Machine Learning 2004).
- Reichardt, Bornholdt, *Statistical mechanics of community detection* (PRE 2006):
  modularity as a Potts spin-glass ground state with resolution γ.
- Ronhovde, Nussinov: replica-correlation multiresolution community detection
  (PRE 2009); local resolution-limit-free "Absolute Potts Model" (PRE 2010).
- Wang, Swendsen, Kotecký, antiferromagnetic Potts cluster algorithm (PRL 1989).
- Swendsen, Wang, *Replica Monte Carlo simulation of spin glasses* (PRL 1986);
  Houdayer, two-replica cluster algorithm (EPJ B 2001).
- Davis 1967 (weak balance / clusterability); Harary 1953 (structural balance).
- Domany, Hed, Palassini, Young, PRB 64, 224406 (2001): the earlier
  state-clustering paper the 2004 and 2007 work builds on (not in the folder).
