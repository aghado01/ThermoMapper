Searched for "backbone persistence"
Searched for "backbone"
Searched for "ConditionedFiltration"
Listed directory ph
Viewed dev-sequencing.md:1-102
Viewed p1a-residual-prior-brief.md:1-100
Viewed p1b-delta-reach-brief.md:1-91
Viewed p1c-nonmonotone-zigzag-brief.md:1-85
Viewed sol-ph-analysis.md:1-100
Searched for "backbone"
Viewed sol-ph-dev-discussion.md:200-250
Viewed sol-ph-dev-discussion.md:330-400
Searched for "conditioned-persistence-synthesis"
Viewed sol-ph-dev-discussion.md:401-500
Viewed sol-ph-dev-discussion.md:501-560
Used MCP tool: git/git_log
Used MCP tool: git/git_status
Listed directory ph
Viewed ResidualPrior.cs:1-100
Viewed ConditionedFiltration.cs:1-99
Viewed MagneticLaplacianOperator.cs:1-60

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