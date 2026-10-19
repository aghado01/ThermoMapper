### Conceptual Overview

**Conditioned backbone persistence** reframes persistent homology from asking *"what is alike?"* (unconditioned Rips clustering) to asking **"what returns, once a known reference scaffold or prior is already in place?"**

1. **The Core Mechanism:**
   - A **backbone** $K_0$ is anchored at filtration value $\varepsilon_0 = 0$ (absorbing the base 0-skeleton / path connectivity).
   - Content edges (similarity chords or predictive relations) enter at their distance or prediction residual.
   - Any 1-cycle formed against the backbone is an $H_1$ generator representing a **return relative to the prior**. Chords born early are expected by the prior; chords born late represent persistent anomalies or surprises.

2. **The Lineage & Generalization:**
   - **SIFTS** (Sequential Information-Filtered Topological Summaries) was the degenerate precursor: reading-order backbone + static TF-IDF cosine similarity.
   - **Li (2025, [arXiv:2508.11646](https://arxiv.org/abs/2508.11646)):** The backbone was generalized from a static path into a **generative edge prior** $t_j \approx t_i + \tau_{ij}$. Edges are admitted when their residual $r_{ij} = |t_j - (t_i + \tau_{ij})| < \delta$. SIFTS is recovered as the degenerate prior $\tau \equiv 0$.
   - **The $\Phi \perp \Psi$ Duality:**
     - **$\Phi$ (Homology / Content):** Evaluated over filtration axes: slack $\delta$ and prediction horizon/reach $\Delta$ ($|\tau| \le \Delta$). Answers *which topological returns exist*.
     - **$\Psi$ (Sheaf / Context):** Evaluated over a connection / gauge field $q$ (e.g., $U(1)$ magnetic phase). Answers *whether directional returns cohere and what flux they enclose*. Per Wei–Wei (2021), a magnetic Laplacian on a path is a cellular sheaf.

---

### Work Completed So Far

| Phase | Component | Location | Validation & Status |
|---|---|---|---|
| **P0** | **Undirected conditioned core** | [`ConditionedFiltration.cs`](file:///d:/aghado01/ThermoMapper/src/tda/ph/ConditionedFiltration.cs) | ✅ **Landed & green** in [`ConditionedFiltrationTests.cs`](file:///d:/aghado01/ThermoMapper/tests/TDA.Ph.Tests/ConditionedFiltrationTests.cs). `BuildGraph` unions $K_0$ at $\varepsilon_0$ with similarity edges; `ComputeBarcode` routes through the existing Rips $\to$ involuted $Z_2$ reducer without altering reduction code. Confirms SIFTS ($\tau \equiv 0$) as degenerate case. |
| **P1a** | **Residual-prior admission (monotone $\delta$)** | [`ResidualPrior.cs`](file:///d:/aghado01/ThermoMapper/src/tda/ph/ResidualPrior.cs#L46-L77) | ✅ **Landed & green** in [`ResidualPriorTests.cs`](file:///d:/aghado01/ThermoMapper/tests/TDA.Ph.Tests/ResidualPriorTests.cs). Admits candidate edges by symmetrized residual $r_{ij} = \|t_j - (t_i + \tau_{ij})\|$. Confirms that predicted returns are born early ($\approx 0$), while unpredicted returns are born late (persistent surprise). |
| **P1b** | **$\Delta$ reach horizon** | [`ResidualPrior.ReachSlices`](file:///d:/aghado01/ThermoMapper/src/tda/ph/ResidualPrior.cs#L86-L99) | ✅ **Landed & green** in [`ResidualPriorTests.cs`](file:///d:/aghado01/ThermoMapper/tests/TDA.Ph.Tests/ResidualPriorTests.cs). Adds a prior reach cutoff $|\tau| \le \Delta$, turning $(\delta, \Delta)$ into a monotone bifiltration read via nested 1D slices in ordinary PH. |
| **Spectral Substrate** | **Magnetic Laplacian ($U(1)$ operator)** | [`MagneticLaplacianOperator.cs`](file:///d:/aghado01/ThermoMapper/src/graphs/spectral/MagneticLaplacianOperator.cs) | ✅ **Landed & green** in [`MagneticFluxTests.cs`](file:///d:/aghado01/ThermoMapper/tests/VizCore.Tests/MagneticFluxTests.cs). Matrix-free $L^{(q)} = D - A \circ \exp(i 2\pi q \Theta)$ embedded into a real $2n \times 2n$ operator solved via [`LOBPCG`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/Lobpcg.cs). Proves flux grades returns by backbone span. |

---

### Direction & Next Phases

The roadmap follows the synthesis laid out in [`dev-sequencing.md`](file:///d:/aghado01/ThermoMapper/issues/ph/dev-sequencing.md) and [`sol-ph-dev-discussion.md`](file:///d:/aghado01/ThermoMapper/issues/ph/sol-ph-dev-discussion.md) (*"Ship $\Phi$, earn $\Psi$"*):

```
P0 (Core) ──> P1a/P1b (Residual & Reach) ──> P1c (Zigzag H0 Driver) ──> P2 (Directedness & Flux) ──> P3/P4 (Sheaf / PSL)
                                                   │
                                                   └──> P1d (Filled H1 Zigzag) / Z6 (Multiparameter)
```

1. **P1c · Non-Monotone Reader ($H_0$ Zigzag) — [Immediate Next Step]:**
   - **Brief:** [`p1c-nonmonotone-zigzag-brief.md`](file:///d:/aghado01/ThermoMapper/issues/ph/p1c-nonmonotone-zigzag-brief.md).
   - **Goal:** When reading dynamic sweeps (trajectories, sliding windows, or temperature $T$-sweeps), edges enter *and* leave. Monotone PH cannot represent component splits.
   - **Action:** Build `ConditionedZigzag.cs` as a sweep-to-churn diffing driver that feeds consecutive frame diffs directly into the parity-complete [`GraphZigzag.Compute`](file:///d:/aghado01/ThermoMapper/src/tda/ph/GraphZigzag.cs) engine.

2. **P1d · Filled-$H_1$ over Zigzag (Deferred):**
   - `GraphZigzag` computes 1-complex (unfilled graph cycle) persistence. Reconciling triangle-filled Rips $H_1$ over non-monotone sweeps requires cell-level [`ZigzagFiltration`](file:///d:/aghado01/ThermoMapper/src/tda/ph/ZigzagFiltration.cs) and resolving whether cycles die by triangle fill or edge departure.

3. **P2 · Directedness v1.5:**
   - **Combinatorial:** Implement directed flag complexes (Lütgehetmann et al. 2020) over ordered simplices, riding the existing $Z_2$ column reduction.
   - **Spectral:** Execute experiments M1–M3 using [`MagneticLaplacianOperator`](file:///d:/aghado01/ThermoMapper/src/graphs/spectral/MagneticLaplacianOperator.cs). The key validation gate is **M2**: two cycles with identical persistence differentiated by directional magnetic flux.

4. **P3 & P4 · Sheaf Cohomology & Persistent Sheaf Laplacians:**
   - **P3:** Formalize cellular sheaves over cell complexes, treating $H^1(\mathcal{F}) \ne 0$ as a cohomological obstruction / coherence failure.
   - **P4:** Compute continuous $\lambda_q$ spectral surfaces and persistent sheaf Laplacians ($\ker \Delta_q^{K, L} \cong \text{persistent sheaf Betti}$) using the Wei–Wei framework.
   - **Z6:** Combine multiparameter persistence with non-monotone zigzag sweeps (*quasi-zigzag*).

   Searched for "relative"
Searched for "VSMJ2011"
Searched for "VSMJ"
Viewed sol-ph-dev-discussion.md:355-370
Viewed PersistentInvolutedHomology.cs:1-60
Viewed PersistentCohomology.cs:1-55
Searched for "relative"

**No, you do not currently have relative homology implemented.** 

What you currently have in code is **ordinary (absolute) persistent homology of the union graph**, $H_*(K)$, with the backbone edges injected at birth $\varepsilon_0 = 0$.

---

### 1. What the code actually does

In [`ConditionedFiltration.cs`](file:///d:/aghado01/ThermoMapper/src/tda/ph/ConditionedFiltration.cs#L86-L97):
```csharp
CsrGraph graph = BuildGraph(n, backbone, similarity, baseWeight);
SimplicialFiltration filtration = RipsFiltration.GraphRips(
    graph, FiltrationWeights.RawDistance, maxDimension, label: "conditioned");
return PersistentInvolutedHomology.Compute(filtration, representatives: true);
```
[`BuildGraph`](file:///d:/aghado01/ThermoMapper/src/tda/ph/ConditionedFiltration.cs#L35-L68) simply assigns weight $\varepsilon_0$ (default `0.0`) to backbone edges and distance $d$ to similarity chords. The resulting complex is fed straight into standard involuted persistent homology.

There is no subcomplex quotienting, no relative chain complex $C_*(K, K_0) = C_*(K)/C_*(K_0)$, and no relative reduction in [`FiltrationCohomologyReducer`](file:///d:/aghado01/ThermoMapper/src/tda/ph/FiltrationCohomologyReducer.cs).

---

### 2. Why your formula is mathematically exact

Your algebraic derivation from the long exact sequence of the pair $(K, K_0)$ is completely correct:

$$\dots \to H_1(K_0) \xrightarrow{i_*} H_1(K) \xrightarrow{j_*} H_1(K, K_0) \xrightarrow{\partial} H_0(K_0) \xrightarrow{i_*} H_0(K) \to \dots$$

When the backbone $K_0$ spans all vertices ($V(K_0) = V(K)$) and is connected:
- $H_0(K_0) \cong H_0(K) \cong \mathbb{Z}_2$, making $i_*$ an isomorphism on $H_0$.
- Hence $\ker(i_*) = 0$, which forces $\operatorname{im}(\partial) = 0$, so $\ker(\partial) = H_1(K, K_0)$.
- By exactness, $\operatorname{im}(j_*) = \ker(\partial) = H_1(K, K_0)$.
- By the first isomorphism theorem:
  $$H_1(K, K_0) \cong H_1(K) / \ker(j_*) = H_1(K) / \operatorname{im}(H_1(K_0))$$

---

### 3. Why the current shortcut works (and where it breaks)

* **When $K_0$ is a simple path or tree (all current tests):**
  $K_0$ is acyclic, so $H_1(K_0) = 0 \implies \operatorname{im}(H_1(K_0)) = 0$.
  Therefore, $H_1(K, K_0) \cong H_1(K)$.
  Because all test fixtures (reading-order paths, token trajectories) are 1D paths without loops, **absolute persistence of the union happens to coincide with relative persistence**.

* **If $K_0$ contains cycles of its own:**
  Because all backbone edges enter at $\varepsilon_0 = 0$, any intrinsic cycle in $K_0$ will be born at filtration value $0.0$.
  In the current engine, **those cycles will appear as bars born at 0** (persisting as $[0, \infty)$ or until filled by 2-simplices). They are **not** quotiented out.

---

### 4. Where this sits in the design

In the synthesis discussion ([`sol-ph-dev-discussion.md`](file:///d:/aghado01/ThermoMapper/issues/ph/sol-ph-dev-discussion.md#L362)), this distinction was noted:
> *"The $\Phi$ / homology half and the absolute-vs-relative pairing is your existing pCoh + the **VSMJ2011 relative-(co)homology** feature already scoped in the port — not new."*

When P0 was built, it took the reuse-maximal path: since linear sequences and trajectories are acyclic paths ($H_1(K_0) = 0$), absolute persistence was used as a drop-in stand-in to avoid modifying the reduction core.

If you ever pass a backbone with internal loops (e.g. a prior containing cycles, ring networks, or state-machine loops), you would need to implement true relative persistence—either by clearing $K_0$ columns in the boundary matrix or via VSMJ2011 relative persistent cohomology—to quotient out $\operatorname{im} H_1(K_0)$.