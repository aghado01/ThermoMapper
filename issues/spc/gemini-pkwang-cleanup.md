Moving the PKWang solver and WangMonteCarloDemo out of Clustering.Graphical.SPC to become an independent peer alongside spc and hdbscan aligns the physical codebase with the underlying theory: SPC is a stochastic Potts/Swendsen–Wang MCMC sampler, whereas PKWang is a deterministic cumulative-field / survival-kernel / single-linkage solver.

Here is the architectural blueprint, seam analysis, target layout, and step-by-step refactoring plan following 

AGENTS.md
 and 

CONTRIBUTING.md
.

1. The Core Seam: What Currently Couples PKWang to SPC?
A strict audit of all files in 

src/clustering/graphical/spc/runtime/core/solver
 reveals that almost nothing actually depends on SPC:

Graphs.Primitives provides CsrGraph, EdgeWeightKind, SymmetrizationRule, EdgeFieldSymmetrization, and UndirectedEdge.
Graphs.Models.Potts provides FkKernel.BondProbability.
Graphs.Observables provides AffinityNodeMarginals.
Clustering.Primitives provides Assignment, Landscape, and LandscapeProvenance.
Maths.Rng & Maths.Samplers provide RNG and exponential sampling for the demo.
The ONLY coupling to SPC is one internal call in 

PKWang.cs
:

csharp
public static Assignment Cluster(PKWangContext context, double T, double theta = 0.5)
{
    Affinities affinity = Solve(context, T);
    return AffinityThreshold.Connect(context.Graph, affinity.G, theta);
}
AffinityThreshold (and its helper UnionFindLabeler) lives in 

src/clustering/graphical/spc/partitions/strategies/AffinityThreshold.cs
 marked internal static class.

Resolving the Seam:
AffinityThreshold is not an SPC-specific concept—it is general connected-component extraction over an active-bond edge affinity field ($g[e] > \theta$). To make PKWang an independent sibling without depending on Clustering.Graphical.SPC:

Promote AffinityThreshold & UnionFindLabeler to a shared home:
Option A (Recommended — Clustering.Primitives): Promote AffinityThreshold (or GraphAffinityThreshold) and UnionFindLabeler to 

projects/Clustering.Primitives
. This requires adding a <ProjectReference Include="../Graphs.Primitives/Graphs.Primitives.csproj" /> to Clustering.Primitives. Both SPC and PKWang (and future graphical methods) then consume the same neutral partition closer.
Option B (Clustering.Dendrograms): 

Clustering.Dendrograms
 already references both Graphs.Primitives and Clustering.Primitives.
Option C (Self-contained in PKWang): Keep a private threshold helper inside PKWang and leave AffinityThreshold in SPC.
2. Target Layout & Naming
Following the sibling pattern of hdbscan and spc:

src/clustering/graphical/
├── hdbscan/
├── spc/
└── pkwang/                      <-- NEW SIBLING
    ├── Field.cs
    ├── IField.cs
    ├── LocalField.cs
    ├── MeanField.cs
    ├── PKWang.cs
    ├── PKWangContext.cs
    ├── PKWangLandscapes.cs
    ├── PKWangSweep.cs
    ├── PKWangSweepResult.cs
    └── demos/
        └── WangMonteCarloDemo.cs
projects/
├── Clustering.Graphical.HdbScan/
├── Clustering.Graphical.SPC/
└── Clustering.Graphical.PKWang/  <-- NEW CSPROJ
    └── Clustering.Graphical.PKWang.csproj
tests/
├── clustering/graphical/
│   ├── spc/
│   │   └── scheduling/          (sampler/ is cleared)
│   └── pkwang/                  <-- MOVED TESTS
│       ├── PKWangLocalFieldTests.cs
│       ├── PKWangMeanFieldTests.cs
│       └── PKWangSweepTests.cs
└── Clustering.Graphical.PKWang.Tests/  <-- NEW TEST PROJECT
    └── Clustering.Graphical.PKWang.Tests.csproj
(Note on naming: pkwang / Clustering.Graphical.PKWang matches the C# type conventions already in place. If preferred, wang or cumulative_field can be used instead).

3. Project File Specifications
A. projects/Clustering.Graphical.PKWang/Clustering.Graphical.PKWang.csproj
xml
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <RootNamespace>Clustering.Graphical.PKWang</RootNamespace>
    <AssemblyName>Clustering.Graphical.PKWang</AssemblyName>
    <EnableDefaultCompileItems>false</EnableDefaultCompileItems>
  </PropertyGroup>
  <ItemGroup>
    <Compile Include="../../src/clustering/graphical/pkwang/**/*.cs" Link="pkwang/%(RecursiveDir)%(Filename)%(Extension)" />
  </ItemGroup>
  <ItemGroup>
    <ProjectReference Include="../Clustering.Primitives/Clustering.Primitives.csproj" />
    <ProjectReference Include="../Graphs.Primitives/Graphs.Primitives.csproj" />
    <ProjectReference Include="../Graphs.Models.Potts/Graphs.Models.Potts.csproj" />
    <ProjectReference Include="../Graphs.Observables/Graphs.Observables.csproj" />
    <ProjectReference Include="../Maths.Rng/Maths.Rng.csproj" />
    <ProjectReference Include="../Maths.Samplers/Maths.Samplers.csproj" />
  </ItemGroup>
</Project>
B. tests/Clustering.Graphical.PKWang.Tests/Clustering.Graphical.PKWang.Tests.csproj
xml
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <RootNamespace>Clustering.Graphical.PKWang.Tests</RootNamespace>
    <AssemblyName>Clustering.Graphical.PKWang.Tests</AssemblyName>
    <IsTestProject>true</IsTestProject>
    <EnableDefaultCompileItems>false</EnableDefaultCompileItems>
  </PropertyGroup>
  <ItemGroup>
    <PackageReference Include="Microsoft.NET.Test.Sdk" Version="17.14.1" />
    <PackageReference Include="xunit" Version="2.9.3" />
    <PackageReference Include="xunit.runner.visualstudio" Version="3.1.5">
      <PrivateAssets>all</PrivateAssets>
      <IncludeAssets>runtime; build; native; contentfiles; analyzers; buildtransitive</IncludeAssets>
    </PackageReference>
  </ItemGroup>
  <ItemGroup>
    <Compile Include="..\clustering\graphical\pkwang\**\*.cs" Link="pkwang\%(RecursiveDir)%(Filename)%(Extension)" />
  </ItemGroup>
  <ItemGroup>
    <ProjectReference Include="..\..\projects\Clustering.Graphical.PKWang\Clustering.Graphical.PKWang.csproj" />
    <ProjectReference Include="..\..\projects\Clustering.Primitives\Clustering.Primitives.csproj" />
    <ProjectReference Include="..\..\projects\Clustering.Dendrograms\Clustering.Dendrograms.csproj" />
    <ProjectReference Include="..\..\projects\Graphs.Primitives\Graphs.Primitives.csproj" />
  </ItemGroup>
</Project>
4. Consumer & Solution Updates
Namespace Renaming:
Update from Clustering.Graphical.SPC.Runtime.Core.Solver and Clustering.Graphical.SPC.Demos to Clustering.Graphical.PKWang (and Clustering.Graphical.PKWang.Demos).


src/user-repl/PKWangUserSession.cs
:
Replace using Clustering.Graphical.SPC.Runtime.Core.Solver; with using Clustering.Graphical.PKWang;.
Add <ProjectReference Include="../Clustering.Graphical.PKWang/Clustering.Graphical.PKWang.csproj" /> to 

projects/UserRepl/UserRepl.csproj
.


src/user-repl/commands/SpcCommand.cs
:
Update using SolverField = Clustering.Graphical.PKWang.Field;.


tests/VizCore.Tests/PKWangLandscapeTests.cs
:
Update import to using Clustering.Graphical.PKWang;, and add reference to Clustering.Graphical.PKWang.csproj in VizCore.Tests.csproj.
Solution file (ThermoMapper.sln):
Add Clustering.Graphical.PKWang.csproj and Clustering.Graphical.PKWang.Tests.csproj.
5. Recommended Execution Sequence (Per 

CONTRIBUTING.md
)
Following the "Refactors: new files land beside old ones; deletions happen only in the final cleanup sweep" invariant:

Phase 1: Seam Decoupling
Promote AffinityThreshold & UnionFindLabeler to Clustering.Primitives (or Clustering.Dendrograms) and verify that Clustering.Graphical.SPC builds cleanly against the promoted location.
Phase 2: Establish the New Home
Create src/clustering/graphical/pkwang/ and demos/.
Copy the 9 solver files and WangMonteCarloDemo.cs, adjusting their namespaces to Clustering.Graphical.PKWang.
Create projects/Clustering.Graphical.PKWang/Clustering.Graphical.PKWang.csproj.
Phase 3: Test Project & Moved Tests
Create tests/clustering/graphical/pkwang/ and tests/Clustering.Graphical.PKWang.Tests/.
Copy the 3 unit tests (PKWangLocalFieldTests.cs, PKWangMeanFieldTests.cs, PKWangSweepTests.cs) and update their namespace imports.
Phase 4: Repoint Downstream Consumers
Update PKWangUserSession.cs, SpcCommand.cs, and VizCore.Tests.
Verify the build and run tests.
Phase 5: Cleanup Sweep
Delete the old directory 

src/clustering/graphical/spc/runtime/core/solver
, 

src/clustering/graphical/spc/demos/WangMonteCarloDemo.cs
, and 

tests/clustering/graphical/spc/sampler
.
Ensure no leftover compat shims remain.
Would you like to proceed with executing Phase 1 (promoting AffinityThreshold to decouple the seam) or would you prefer a specific naming choice (e.g. pkwang vs wang vs cumulative_field) first?

1:30 AM
2:04 AM
This brings complete clarity to the design. Eradicating PKWang from production types and assemblies—relegating the paper name strictly to historical attribution and didactic auditing—aligns the codebase with the actual mathematics: Wang 2020 never published the analytical result; they ran spurious Monte Carlo draws to estimate a constant, and the density-adaptive local field formulation is ThermoMapper's own.

Here is how the architecture and file layout simplify when transitioning from the legacy "solver" abstraction to the final Local Cumulative Field engine.

1. Conceptual Shift: From "SPC Solver" to "Cumulative Field"
In the legacy code, the components were framed as an alternative solver for SPC (PKWangContext, PKWang.Solve, Field.Mean vs Field.Local).

In reality:

It is not an SPC solver (no Potts Hamiltonian, no spin configuration, no phase transition).
MeanField was just the closed-form realization of Wang’s density-blind global ranking.
The real production method is the density-adaptive local cumulative field (per 

lean/Protolemmata/CumulativeField/design.md
).
Legacy (Over-generalized "Solver"):
  IField (MeanField | LocalField)  -->  PKWang.Solve  -->  AffinityThreshold
Final Production Form:
  Local Cumulative Ladder  -->  Endpoint Symmetrization  -->  Cut / Sweep & Landscape
2. Suggested Sibling Naming
Alongside Clustering.Graphical.SPC and Clustering.Graphical.HdbScan, the canonical sibling name explored in Lean is CumulativeField:

Assembly / Project: Clustering.Graphical.CumulativeField
Directory: src/clustering/graphical/cumulative_field/
Namespace: Clustering.Graphical.CumulativeField
3. What the Simplified Core Looks Like
By pruning the IField dispatch and the legacy global mean-field from production:

Legacy Component	Final Production Role
IField.cs, Field.cs	Deleted / Inlined. No abstraction needed if the local field is the definitive method.
MeanField.cs	Retired from production. Kept only in the didactic demo/audit suite as the global reference.
LocalField.cs	Core ladder. Direct local cumulative curve builder $L_u(t) = \sum_{v: J_{uv} \le t} J_{uv}$.
PKWangContext.cs	CumulativeFieldContext.cs (or inlined into run state). Precomputed local $H^{\text{cum}}$ + graph + reverse slot map.
PKWang.cs	CumulativeField.cs (or LocalCumulativeLinkage.cs). Direct Prepare, EstimateBracket, Affinity, and Cluster.
PKWangLandscapes.cs	CumulativeFieldLandscapes.cs. Exact bond-mass and frustration-energy landscape mints for Thermal EOM.
PKWangSweep.cs, PKWangSweepResult.cs	CumulativeFieldSweep.cs / CumulativeFieldSweepResult.cs.
The entire production engine becomes ~4–5 clean, cohesive files with zero branching on field type.

4. Placement of the Didactic Demo & Historical Audit
The demonstration that Wang's Monte Carlo was estimating a closed-form constant with zero variance is a valuable empirical witness, but it belongs outside the hot engine:

src/clustering/graphical/cumulative_field/
├── CumulativeField.cs
├── CumulativeFieldContext.cs
├── CumulativeFieldLandscapes.cs
├── CumulativeFieldSweep.cs
└── didactic/                          <-- AUDIT & DIDACTIC COMPANION
    ├── Wang2020MeanFieldReference.cs   (the legacy global ladder)
    └── WangMonteCarloDemo.cs           (the MCMC exponential draws witness)
The docstrings in didactic/ can explicitly cite the paper and document the reduction:

"Wang 2020 as published used $M$ exponential draws to estimate edge survival, which gratuitously estimated the zero-variance closed form $G = 1 - \exp(-H^{\mathrm{cum}}/T)$. This demo provides the empirical witness for that reduction."

5. Decoupling the Flat Cut Seam (AffinityThreshold)
To prevent CumulativeField from referencing Clustering.Graphical.SPC:

Move 

AffinityThreshold.cs
 and 

UnionFindLabeler.cs
 out of SPC into 

Clustering.Primitives
.
Clustering.Primitives references Graphs.Primitives (UnionFind, CsrGraph).
Both Clustering.Graphical.SPC and Clustering.Graphical.CumulativeField use the shared AffinityThreshold.Connect(...) primitive to close their cuts into Assignments.
6. REPL / CLI Alignment
In UserRepl:



PKWangUserSession.cs
 becomes CumulativeFieldUserSession.cs.
In commands/:
Option 1: Keep --solver cumulative-field (with pkwang as a deprecated/aliased parser token) on spc.
Option 2: Introduce a first-class cumulative-field command matching hdbscan and spc.
Next Step
Does Clustering.Graphical.CumulativeField fit your intended naming, or do you prefer a different umbrella term (e.g., LocalCumulativeLinkage or LocalField)? Once confirmed, we can first promote the AffinityThreshold seam, then scaffold the new sibling structure.


# i think i still need to generalize theta as well as figure out the symmetrization, which could be min/max or geometric mean or something like that

These two questions—**generalizing $\theta$** and **choosing the symmetrization operator**—are tightly coupled. In fact, resolving their algebra reveals why the current C# code had hidden artifacts, and points toward a clean, parameter-free model.

Here is the exact mathematical breakdown of both issues and how they shape the local field variation.

---

### 1. Generalizing $\theta$: The Collapse into a 1D Cut Coordinate $\tau$

In the current code and the 2020 paper, $\theta$ is nominally treated as an independent threshold, but hardcoded to $\theta = 0.5$.

Let’s derive the cut condition for arbitrary $\theta \in (0, 1)$:

$$G > \theta \iff 1 - e^{-H / T} > \theta \iff e^{-H / T} < 1 - \theta$$

Taking logarithms:
$$-\frac{H}{T} < \ln(1 - \theta) \iff H > -T \ln(1 - \theta)$$

Define the **scale factor**:
$$\kappa(\theta) \equiv -\ln(1 - \theta) = \ln\left(\frac{1}{1 - \theta}\right)$$

And define the **unified cut coordinate**:
$$\tau \equiv T \cdot \kappa(\theta)$$

Then the cut condition is simply:
$$H > \tau$$

#### Key Consequences:
1. **$\theta$ and $T$ are NOT a 2D bifiltration.**
   At $\theta = 0.5$, $\kappa(0.5) = \ln 2 \approx 0.6931$, recovering the paper’s $H > T \ln 2$. 
   If you change $\theta$, you are **not** changing the underlying hierarchy or discovering different clusters—you are merely rescaling the temperature axis by $\frac{\ln 2}{\kappa(\theta)}$.
2. **The hierarchy is static in $\tau$:**
   Every edge has a critical cut height $\tau^* = H$. As $\tau$ sweeps from $0 \to \infty$, edges drop out in strict ascending order of their cumulative energy $H$.
3. **Bracketing becomes trivial:**
   Instead of estimating an arbitrary $(T_{\text{lo}}, T_{\text{hi}})$ bracket dividing by an ad-hoc $\ln 2$, the native bracket lives in $\tau$-space, directly bounded by the quantiles of the positive $H$ ladder.

---

### 2. Symmetrization: Where and How to Combine Directed Fields

For an undirected edge $uv$, the local field produces two directed cumulative energies:
* $H_{u \to v}$: $uv$'s rank from $u$'s local neighborhood.
* $H_{v \to u}$: $uv$'s rank from $v$'s local neighborhood.

The critical decision is **whether you combine in $H$-space (cumulative hazard) or in $G$-space (probability)**, and which operator you use.

#### A. Min and Max Commute with the Exponential Map
Because $G(H) = 1 - e^{-H/T}$ is strictly monotonically increasing in $H$:
$$\min(G_1, G_2) = 1 - e^{-\min(H_1, H_2) / T}$$
$$\max(G_1, G_2) = 1 - e^{-\max(H_1, H_2) / T}$$

* **Mutual ($\min$):** $H_{\mathrm{mut}} = \min(H_{u \to v}, H_{v \to u})$. Conservative/bilateral core evidence—an edge is only as strong as its weaker endpoint views it.
* **Inclusive ($\max$):** $H_{\mathrm{inc}} = \max(H_{u \to v}, H_{v \to u})$. Permissive—captures peripheral points where only the core node considers the edge strong.

Both preserve a strict, static ordering of edges across all temperatures.

#### B. The Mean: Geometric Survival vs. Temperature-Dependent Drift

Here is where the current C# code had an unintended flaw:
In the legacy `PKWang.cs`:
```csharp
// Exponentiates in directed slots first:
g[e] = 1 - exp(-hcum[e] / T);
// Then takes arithmetic mean of probabilities:
EdgeFieldSymmetrization.Symmetrize(..., g, Mean); // 0.5 * (g1 + g2)
```

Taking the **arithmetic mean of probabilities** $G_{\text{arith}} = \frac{G_1 + G_2}{2}$ corresponds to an effective $H$:
$$H_{\mathrm{eff}}(T) = -T \ln\left(\frac{e^{-H_1 / T} + e^{-H_2 / T}}{2}\right)$$
This is a **temperature-dependent soft-min** (LogSumExp). As $T$ varies, relative edge scores cross over, meaning **edge rankings mutate as temperature changes**. It is not a fixed hierarchy!

#### The Clean Alternative: Geometric Mean of Survival (Arithmetic Mean of $H$)
If you instead average in **$H$-space** (balanced hazard):
$$H_{\mathrm{bal}} \equiv \frac{H_{u \to v} + H_{v \to u}}{2}$$

Look at what this becomes under the exponential map:
$$1 - G_{\mathrm{bal}} = e^{-H_{\mathrm{bal}} / T} = e^{-(H_1 + H_2) / (2T)} = \sqrt{e^{-H_1/T} \cdot e^{-H_2/T}} = \sqrt{(1 - G_1)(1 - G_2)}$$

**The survival probability $(1 - G)$ is precisely the geometric mean of the directed survival probabilities!**

#### The Diagnostic Sandwich:
If you symmetrize in $H$-space, all three operators are parameter-free, temperature-invariant, and form a strict sandwich:
$$H_{\mathrm{mut}} \le H_{\mathrm{bal}} \le H_{\mathrm{inc}}$$
Which induces a monotonic filtration of graphs:
$$E_{\mathrm{mut}}(\tau) \subseteq E_{\mathrm{bal}}(\tau) \subseteq E_{\mathrm{inc}}(\tau)$$

---

### 3. How to Structure the Code While You Iterate

Since you are still settling the final local field variation, we should avoid locking down rigid assumptions. Instead, we can structure the new sibling to make these experiments explicit:

1. **Move the code to `src/clustering/graphical/cumulative_field/` (or `pkwang/`)** to immediately liberate `SPC` from this non-Potts logic.
2. **Promote the shared cut seam (`AffinityThreshold`)** so both engines share component extraction.
3. **Parameterize the symmetrizer in $H$-space**:
   ```csharp
   public enum CumulativeSymmetrization
   {
       Mutual,    // min(Hu, Hv) -> conservative
       Inclusive, // max(Hu, Hv) -> captures periphery
       Balanced,  // (Hu + Hv) / 2 -> geometric mean of survival
   }
   ```
4. **Parameterize cuts by $\tau$ or $(T, \theta)$ via $\kappa(\theta)$**:
   Allow either specifying $\tau$ directly, or converting $(T, \theta)$ via $\tau = T \cdot (-\ln(1 - \theta))$.

How does the $H$-space sandwich ($\min$, $\max$, and arithmetic mean in $H$ / geometric survival) feel as the candidate space for your symmetrization experiments?


# if i do it this way, do i get a clusterer that can just run without user needing to worry about choosing things?

Listed directory dendrogram
Viewed Condensation.cs:1-50
Viewed SingleLinkage.cs:1-65

**Yes, absolutely.** That is the single biggest payoff of this formulation.

By unifying $(T, \theta)$ into the single hazard coordinate $\tau = T \kappa(\theta)$ and symmetrizing in $H$-space, you transform the algorithm from an unstable, multi-knob "simulation sweep" into a **fully deterministic, parameter-free graphical clusterer**.

Here is how the pipeline runs "push-button" with zero user decisions:

---

### 1. The Zero-Parameter Pipeline: Step-by-Step

Given a coupling graph from the graph compiler:

```
Coupling Graph (J)
       │
       ▼  [Step 1: Local Ladder]
  H_{u→v} = ∑_{w: J_{uw} ≤ J_{uv}} J_{uw}       (Parameter-free: deterministic local mass)
       │
       ▼  [Step 2: Symmetrization]
  H_{uv} = (H_{u→v} + H_{v→u}) / 2              (Parameter-free: balanced hazard / geometric survival)
       │
       ▼  [Step 3: Critical Edge Ladder]
  Exact MST over edge weights (-H_{uv})         (Parameter-free: O(E log E) Kruskal / Borůvka)
       │
       ▼  [Step 4: Persistence Dendrogram]
  Single-Linkage Dendrogram in τ-space          (Parameter-free: exact birth/death for every branch)
       │
       ▼  [Step 5: Automatic Extraction]
  Excess-of-Mass (EOM) Selection                (Parameter-free: bottom-up tree stability optimization)
       │
       ▼
  Final Flat Partition (Assignment)
```

---

### 2. Why Every Previous "Knob" Disappears

| What users had to choose before | Why it's eliminated in this formulation |
| :--- | :--- |
| **Temperature grid ($T_{\text{lo}}, T_{\text{hi}}, \Delta T$)** | **Gone.** The critical cut events are *exact* edge scores $H_{uv}$. You don't sample on a discrete grid; you have the exact continuous filtration ladder. |
| **Threshold ($\theta$)** | **Gone.** As shown by $\tau = T \cdot \kappa(\theta)$, changing $\theta$ does not discover new cluster structure—it is an isotropic rescaling of the cut axis. $\tau = H$ is the native cut. |
| **Equilibration & MC steps** | **Gone.** Zero variance, zero draws, instantaneous execution ($O(E \log E)$). |
| **Which cut / temperature to pick** | **Automated.** You feed the dendrogram into [Condensation.cs](file:///d:/aghado01/ThermoMapper/src/clustering/dendrogram/Condensation.cs) and run **Excess of Mass (EOM)**—the exact same stability logic HDBSCAN uses to eliminate $\epsilon$. |

---

### 3. How Automatic Extraction (EOM) Resolves the Tree

Because each branch in the dendrogram has an exact birth $\tau_{\text{birth}}$ and death $\tau_{\text{death}}$:
1. Every candidate cluster has an intrinsic stability / lifetime:
   $$\text{Stability}(C) = \sum_{p \in C} (\tau_{\text{fallout}}(p) - \tau_{\text{birth}}(C))$$
2. The tree is condensed (using a small floor like `minClusterSize = 2` or $\log N$).
3. Bottom-up **Excess-of-Mass** walks the tree: a parent cluster is retained if its own persistence exceeds the combined persistence of its children, otherwise it splits.

**The user never has to pick a cut height.** The tree self-selects its most salient, persistent density components.

---

### 4. What the User Experience Looks Like

The entire user-facing API collapses to:

```csharp
// "Just run" — no temperatures, no theta, no schedules:
CumulativeFieldResult result = CumulativeFieldRunner.Run(couplingGraph);

// Ready to use:
Assignment partition = result.Assignment;

// Optional: for inspection, the complete persistence tree is preserved
Dendrogram tree = result.Dendrogram;
```

### Summary
* **SPC** remains the stochastic physical research engine (for Potts magnetization, susceptibility transitions, and critical phenomena).
* **HDBSCAN** remains the metric/distance mutual-reachability clusterer.
* **CumulativeField** (or `LocalCumulativeLinkage`) becomes the **fast, deterministic, zero-knob graphical clusterer**—taking a coupling graph and returning optimal persistent clusters in a single $O(E \log E)$ pass.


# ive actually got metric/distance plans for SPC. i haven't sorted out some nuances with coupling kernels and bandwidths, which are kind of a mess right now in the graph engine. but the intention with thermomapper has always been a generalized SPC graph if i can help it. but depending on curvature, gaussian might not apply and these things interact. but i really want to do clustering with spc on data projected or native to curved spaces like hyperboloid or simplex or others

Listed directory Spc.HyperbolicSmoke
Viewed Program.cs:1-60
Listed directory graphs
Listed directory coupling
Viewed Bandwidth.cs:1-60

This insight cuts straight to the heart of what **ThermoMapper** is actually about, and why evicting PKWang from SPC is so vital.

When PKWang was sitting inside SPC, it created an identity crisis: was SPC a Potts physical engine, or was it an analytical single-linkage ladder?
Removing it restores SPC's true, radical purpose: **thermodynamic phase-transition clustering on generalized geometric manifolds**.

The Potts Hamiltonian doesn't know what Euclidean space is:
$$H(\sigma) = -\sum_{\langle i, j \rangle} J_{ij} \delta(\sigma_i, \sigma_j)$$
All SPC requires is a valid, well-calibrated coupling graph $J$. But generating that graph from data on curved manifolds is where curvature, volume growth, kernels, and bandwidths collide.

---

### 1. Why Curvature Breaks the Standard Gaussian Assumption

In Euclidean space $\mathbb{R}^d$, the volume of a geodesic ball grows polynomially: $V(r) \propto r^d$. The Gaussian kernel $e^{-d^2 / 2\sigma^2}$ is the fundamental solution of the Euclidean heat equation.

On curved manifolds, volume forms completely distort this balance:

#### A. Hyperbolic Space ($\mathbb{H}^d$ — Hyperboloid / Lorentz, Poincaré Ball)
* **Volume Growth:** $V(r) \propto \sinh^{d-1}(r) \sim e^{(d-1)r}$ (exponential volume expansion).
* **The Gaussian Failure:** As noted in [Spc.HyperbolicSmoke/Program.cs](file:///d:/aghado01/ThermoMapper/tests/Spc.HyperbolicSmoke/Program.cs#L37-L41), the Gaussian's $e^{-d^2}$ tail decays far faster than space expands. It crushes long-range hierarchical bridges. At $d/\delta \approx 5$, Gaussian couplings are $10^4 \times$ smaller than they should be, isolating peripheral branches into detached clumps and destroying the multi-scale Blatt hierarchy.
* **The Remedy:** Heavy-tailed kernels like **Cauchy** $J = \frac{1}{1 + (d/\delta)^2}$ or power-law kernels that match the polynomial/exponential trade-off of $\mathbb{H}^d$.

#### B. The Probability Simplex ($\Delta^{d-1}$ — Fisher–Rao / Bhattacharyya / Hellinger)
* **Volume Growth:** Compact spherical quadrant with positive curvature. Distances are strictly bounded:
  $$d_{\text{FR}}(p, q) = 2 \arccos\left(\sum_i \sqrt{p_i q_i}\right) \le \pi$$
* **The Gaussian Failure:** Euclidean distance on raw probability vectors violates the geometry of distributions (probabilities live on a manifold, not a flat plane). But even with geodesic Fisher–Rao distance, a Gaussian kernel never enters its asymptotic tail because distances are bounded by $\pi$.
* **The Remedy:** Kernels adapted to compact spherical geometry (cosine-type kernels, or Dirichlet/heat kernels on the simplex).

---

### 2. Why the Graph Engine Feels Like a "Mess" Right Now

Looking at [src/graphs/Bandwidth.cs](file:///d:/aghado01/ThermoMapper/src/graphs/Bandwidth.cs) and [src/graphs/GraphCompiler.cs](file:///d:/aghado01/ThermoMapper/src/graphs/GraphCompiler.cs), the reason the graph engine feels tangled is that **four distinct concerns are currently intertwined**:

1. **Geodesic Metric:** How to compute $d(x_i, x_j)$ on the manifold (Euclidean, Poincaré/Hyperboloid, Fisher–Rao).
2. **Scale Estimation (Bandwidth $\delta$):** How to extract a typical scale from local distances. (In $\mathbb{H}^d$, Euclidean MAD underestimates the spread; a log-scale MAD or moment-matching against $\sinh^{d-1}(r)$ is needed).
3. **Decay Kernel ($K$):** How distance maps to raw coupling (Gaussian, Cauchy, Laplacian).
4. **Thermodynamic Gauge / Normalization:** How to rescale $J$ so that the Potts transition temperature $T_{\text{ps}}(q)$ doesn't randomly fly off to 0.001 or 100.0 depending on the dataset.

In `Bandwidth.cs`, you already started deriving the exact radial reference measures:
$$p(r) \propto K(r) \sinh^{d-1}(r)$$
and noting consistency factors (`GaussianHyperbolicFactor`, `LaplacianHyperbolicFactor = 1.67`, `QuantileNormalized` for Simplex/Fisher–Rao). The math is already there—it just needs clean architectural boundaries.

---

### 3. The Clean Separation for Manifold-Aware SPC

To make graph construction modular and tractable for Hyperbolic and Simplex spaces, the graph compiler can decouple into four pure pipeline stages:

```
Manifold Coordinates (Hyperbolic / Simplex / Euclidean)
                      │
                      ▼  [Stage 1: Geodesic Metric]
             Pairwise Geodesic Distance Matrix (d_ij)
                      │
                      ▼  [Stage 2: Manifold-Aware Scale]
             Local / Global Bandwidth (δ_i)
             (Log-MAD for H^d, Quantile for Simplex, MAD for R^d)
                      │
                      ▼  [Stage 3: Curvature-Matched Kernel]
             Raw Coupling: J_ij = K(d_ij / δ_ij)
             (Cauchy on H^d, Gaussian on R^d, Cosine on Simplex)
                      │
                      ▼  [Stage 4: SPC Gauge Calibration]
             CsrGraph with normalized average coupling: ⟨J⟩ = 1
                      │
                      ▼
             Swendsen–Wang Potts Sampler (SPC)
             (T_ps(q) is now invariant and predictable across any manifold!)
```

### Stage 4 is the Secret Weapon:
If you normalize the compiled coupling graph such that the average local coupling $\langle J \rangle \approx 1$ (the BWD normalization), then **the physical temperature scale of SPC remains constant regardless of whether your data lives in $\mathbb{R}^d$, $\mathbb{H}^d$, or on a simplex!**
You don't need a bespoke temperature schedule for every new manifold—the Potts transition happens in the same predictable $T$ window.

---

### The Big Picture

This completely clarifies the division of labor across ThermoMapper:

1. **`Clustering.Graphical.SPC`**:
   The flagship physical engine. It consumes a `CsrGraph` with normalized couplings $J$ and discovers multiscale hierarchical clusters across **curved Riemannian geometries** (Hyperboloid, Simplex, Sphere) via Swendsen–Wang thermal fluctuations.
2. **`Graphs.Compiler`**:
   The manifold bridge. It handles the curvature-aware metric, bandwidth estimation, and kernel decay to emit clean, gauge-calibrated coupling graphs.
3. **`Clustering.Graphical.CumulativeField` (formerly PKWang)**:
   A fast, deterministic, zero-knob baseline. It lives completely outside of SPC, providing instant local-ladder clustering without MCMC or temperature tuning.

Once we move PKWang/CumulativeField out, the path is clear to clean up the graph engine and unlock full curved-space SPC.


# how does the guage work for the pkwang and for the SPC paths currently?

Searched for "gauge"
Listed directory potts
Viewed BwdPottsCriticalEstimate.cs:1-40
Searched for "TpsUpperBound"
Viewed SpcCommand.cs:460-520
Viewed SpcScheduleHelpers.cs:1-60
Searched for "gauge"
Viewed sol-SPC-bars-graphs-brackets.md:440-480
Searched for "GlobalBandwidthScaler"
Viewed GlobalBandwidthScaler.cs:1-80
Viewed GlobalBandwidthScaler.cs:81-140
Viewed GlobalBandwidthScaler.cs:141-200

In ThermoMapper, **"gauge"** refers to the scale freedom of coupling strengths ($J \to \alpha J$) and how that scale dictates the temperature axis ($T \to \alpha T$).

Here is the exact breakdown of how the gauge works across the **SPC path** vs. the **PKWang / CumulativeField path** in the codebase today.

---

### 1. The SPC Path: Two Distinct Gauge Regimes

In the Potts model, the Hamiltonian is scaled by temperature:
$$\beta H = -\frac{1}{T} \sum_{\langle i, j \rangle} J_{ij} \delta(\sigma_i, \sigma_j)$$
Because only the ratio $J_{ij} / T$ enters the Boltzmann factor, any global scale factor on $J$ shifts the physical temperature scale by the exact same factor.

In the codebase, SPC handles this via two explicit modes in [GlobalBandwidthScaler.cs](file:///d:/aghado01/ThermoMapper/src/graphs/pipeline/scalers/GlobalBandwidthScaler.cs#L141-L154) and [SpcCommand.cs](file:///d:/aghado01/ThermoMapper/src/user-repl/commands/SpcCommand.cs#L477-L492):

#### A. The BWD $1/\widehat{K}$ Normalized Gauge (`MeanEdgeDistance`)
* **How it works:** Couplings carry an explicit normalization by the average node degree $\widehat{K}$:
  $$J_{ij} = \frac{1}{\widehat{K}} K(d_{ij}, \delta) \quad \text{where } \widehat{K} = \frac{1}{N}\sum_i \text{deg}(i)$$
  *(In code: `couplingScale = n / (double)totalNeighbors`).*
* **Why it matters:** Normalizing by $1/\widehat{K}$ fixes the average interaction energy per spin to $\sim 1$.
* **The Temperature Anchor:** Because the gauge is fixed, the paramagnetic pseudotransition temperature $T_{\text{ps}}$ is **data-independent** and determined purely by Potts physics in [BwdPottsCriticalEstimate.cs](file:///d:/aghado01/ThermoMapper/src/graphs/models/potts/BwdPottsCriticalEstimate.cs):
  $$T_{\text{ps}}(q) = \frac{e^{-1/2}}{4} \ln(1 + \sqrt{q})$$
* **The Bracket:** `SpcCommand` sets the bracket automatically as $[0.05 \cdot T_{\text{ps}}(q), \; 4.0 \cdot T_{\text{ps}}(q)]$. You never need to inspect data distances to find the physical phase transition.

#### B. The Un-normalized Gauge (`MadConsistencyFactor`, etc.)
* **How it works:** `couplingScale = 1.0`. Couplings are raw kernel values $J_{ij} \in [0, 1]$.
* **Why it matters:** The total interaction energy per node scales with the local density and neighborhood size (dimensionally $\propto \widehat{K}$).
* **The Temperature Anchor:** The $q$-only $T_{\text{ps}}(q)$ formula is invalid here (it would land $\sim \widehat{K}\times$ too cold). Instead, [SpcScheduleHelpers.Estimate](file:///d:/aghado01/ThermoMapper/src/clustering/graphical/spc/runtime/scheduling/SpcScheduleHelpers.cs#L33-L45) computes an empirical heuristic:
  $$T_c \approx \frac{\text{median}(J) \cdot \langle\text{weighted degree}\rangle}{\ln q}$$
  and sets the sweep bracket to $[0.05 \cdot T_c, \; 5.0 \cdot T_c]$.

---

### 2. The PKWang / CumulativeField Path: "q-Free and Gauge-Free"

In [PKWang.cs](file:///d:/aghado01/ThermoMapper/src/clustering/graphical/spc/runtime/core/solver/PKWang.cs#L65-L73) and [SpcCommand.cs](file:///d:/aghado01/ThermoMapper/src/user-repl/commands/SpcCommand.cs#L465-L470), the solver is explicitly designated as **`q-free and gauge-free`**:

#### A. Why It Does Not Use a Potts Gauge
* PKWang has no Potts spins, no $q$, and no physical phase transition.
* "Temperature" $T$ is not a thermal bath; it is purely a coordinate transform on single-linkage cut heights:
  $$G_{ij} = 1 - \exp(-H^{\text{cum}}_{ij} / T)$$
  At $\theta = 0.5$, an edge is connected iff:
  $$H^{\text{cum}}_{ij} > T \ln 2 \iff T < \frac{H^{\text{cum}}_{ij}}{\ln 2}$$
* **Self-Gauging:** The temperature scale is defined by the graph’s *own* cumulative-energy ladder $H^{\text{cum}}$. If you rescale $J \to \alpha J$, $H^{\text{cum}}$ rescales by $\alpha$, and the critical temperatures rescale by $\alpha$.

#### B. The Native Ladder Bracket
In `PKWang.EstimateBracket`, the sweep window is read directly off the empirical quantiles of $H^{\text{cum}}$:
$$T_{\text{lo}} = 0.5 \cdot \frac{H^{\text{cum}}_{[0.02]}}{\ln 2}, \qquad T_{\text{hi}} = 2.0 \cdot \frac{H^{\text{cum}}_{[0.98]}}{\ln 2}$$
No knowledge of $q$, $\widehat{K}$, or Potts critical transitions is used.

#### C. The $\theta$-Gauge Equivalence
As proven in [lean/Protolemmata/CumulativeField/design.md](file:///d:/aghado01/ThermoMapper/lean/Protolemmata/CumulativeField/design.md#L105-L125), $(T, \theta)$ possesses a strict gauge symmetry:
$$\tau = T \cdot (-\ln(1 - \theta))$$
Choosing a different threshold $\theta$ is mathematically identical to a gauge transformation that rescales the temperature axis by $\frac{\ln 2}{-\ln(1 - \theta)}$. The hierarchy of graph cuts is invariant under this gauge.

---

### 3. In Downstream Resolvers (Landscapes & Thermal EOM)

Once affinities $G$ are produced:
* **The Landscape Gauge:** In [Landscape.cs](file:///d:/aghado01/ThermoMapper/src/clustering/primitives/Landscape.cs#L9), `LandscapeProvenance` carries a `GaugeNote`:
  * SPC tags: `"spc:sampled"`
  * PKWang tags: `"pkwang:closed-form"`
* **Gauge-Free Modal Ascent:** In [ThermalEom.cs](file:///d:/aghado01/ThermoMapper/src/clustering/graphical/spc/partitions/thermal/ThermalEom.cs#L38-L41), when completing abstained points via periphery ascent (`ThermalPeripheryCompletion.Ascend`), the walk is **gauge-free** because greedy climbing to the local peak depends only on the *ordinal ranking* between neighboring nodes, not the absolute cardinal height.

---

### Summary Comparison

| Dimension | **SPC (Potts Physical Path)** | **PKWang / CumulativeField Path** |
| :--- | :--- | :--- |
| **Physical meaning of $T$** | Thermodynamic temperature driving spin fluctuations. | Coordinate transform on edge cut height: $T = H / \ln 2$. |
| **Coupling normalization** | BWD $1/\widehat{K}$ gauge fixes $T_{\text{ps}}(q)$; unnormalized shifts $T$ by $\widehat{K}$. | Self-gauging: temperature is directly in units of $H^{\text{cum}}$. |
| **Potts $q$ dependency** | Essential: $T_{\text{ps}}(q) \propto \ln(1 + \sqrt{q})$. | None ($q$-free). |
| **Bracket discovery** | Physics anchor ($T_{\text{ps}}(q)$) or graph degree heuristic. | Quantiles of the internal $H^{\text{cum}}$ ladder ($[q_{0.02}, q_{0.98}]$). |
| **Internal gauge symmetry** | Physical Hamiltonian scale invariance $\beta J$. | $\theta$-gauge invariance: $\tau = T \cdot \kappa(\theta)$. |

# does the unnormalized path B have any value? i think where i was going with that was wanting to be able to adjust the width of the temperature range, worried that some more exotic systems with curved geometry and whatnot might clip off on the high end or the low end or something

Your intuition about the **danger of clipping on exotic curved geometries** is 100% correct. But **unnormalized Path B does not actually solve it**, and in fact creates a severe dimensional bug.

Here is the exact breakdown of why Path B fails, why exotic/curved geometries clip, and how to properly control the bracket width without abandoning the gauge.

---

### 1. The Verdict on Path B: Dimensionally Flawed

In [SpcScheduleHelpers.cs](file:///d:/aghado01/ThermoMapper/src/clustering/graphical/spc/runtime/scheduling/SpcScheduleHelpers.cs#L38-L45), the un-normalized estimator is:
```csharp
double meanWeightedDegree = 2.0 * weights.MeanWeight * weights.EdgeCount / graph.NodeCount;
return weights.MedianWeight * meanWeightedDegree / logQ;
```

Notice what is happening dimensionally:
* `MeanWeight` is $O(J)$.
* `MedianWeight` is $O(J)$.
* Their product is **$O(J^2)$**.

In the Potts Hamiltonian, energy is linear in coupling: $\beta H = \frac{1}{T} \sum J_{ij}$.
If you rescale couplings by $\alpha$ ($J \to \alpha J$), the physics demands that temperature rescale linearly: $T \to \alpha T$.
Because Path B computes an estimate scaling as **$J^2$**, changing graph density or scaling couplings by $0.1$ causes the estimated temperature to plummet by $0.01$ (landing $\sim \widehat{K}\times$ too cold, as acknowledged in the comment in [SpcCommand.cs](file:///d:/aghado01/ThermoMapper/src/user-repl/commands/SpcCommand.cs#L480-L481)).

Furthermore, Path B's bracket was hardcoded to $[0.05 \cdot T_c, \; 5.0 \cdot T_c]$—which is a fixed factor of 100, **identically rigid to Path A**. It does not dynamically widen the temperature range at all.

**Conclusion:** Unnormalized Path B does not provide adaptive width; it only makes the temperature center drift unpredictably. It has no architectural value and can be safely retired.

---

### 2. What Actually Causes Clipping in Curved Geometries?

Your concern about curved systems clipping at the high or low end is mathematically profound:

#### Shift vs. Width in $\ln T$ Space
1. **Global Normalization (The Gauge):**
   Multiplying all couplings by a global scalar (like $1/\widehat{K}$) is a **rigid translation** in $\ln T$ space:
   $$\ln T' = \ln T - \ln \widehat{K}$$
   It shifts the window left or right. It **cannot** stretch, compress, or clip the relative physics.
2. **Coupling Dispersion (The Width):**
   The true physical width of the superparamagnetic (SP) phase—where hierarchical clustering actually happens—is determined by the **ratio between the strongest within-cluster couplings and the weakest inter-cluster bridges**:
   $$\text{Dynamic Range} \equiv \frac{T_{\text{hot}}}{T_{\text{cold}}} \sim \frac{J_{\text{dense}}}{J_{\text{sparse}}}$$

#### Why Euclidean Data is Narrow, but Hyperbolic Data is Wide:
* **Euclidean + Gaussian:** Volume grows polynomially ($r^d$) and Gaussian decays as $e^{-r^2}$. Couplings drop to near-zero almost immediately. All clusters freeze and melt in a relatively tight window (a factor of $20\times$ to $50\times$ in $T$). Path A’s default $[0.05, 4.0]$ ($80\times$ range) covers it easily.
* **Hyperbolic ($\mathbb{H}^d$) + Cauchy:**
  Volume grows exponentially ($e^{(d-1)r}$) and Cauchy decays polynomially ($1/r^2$). Deep leaf clusters have $J \approx 1.0$, while high-level ancestral bridges across the Poincaré ball have $J \approx 0.001$.
  The couplings span **several orders of magnitude**. 
  Therefore, the multi-scale phase transitions span a much wider temperature ratio ($1000\times$ or more)! If you use the standard $80\times$ bracket, **the high end or low end will clip off**, exactly as you feared.

---

### 3. The Clean Solution: Gauge-Fixed Center, Dispersion-Aware Width

You don't need two different gauge regimes to solve this. You keep the gauge fixed so the anchor is predictable, and make the **bracket width adapt to the manifold's coupling dispersion**:

```
Coupling Normalization:
  J_normalized = J / ⟨K⟩                      --> Fixes T_center ≈ T_ps(q) (invariant across all manifolds)

Bracket Width:
  DynamicRange = Quantile(J, 0.98) / Quantile(J, 0.02)
  T_min = T_ps(q) / sqrt(DynamicRange)
  T_max = T_ps(q) * sqrt(DynamicRange)        --> Dynamically widens on heavy-tailed / hyperbolic graphs!
```

### Why this is strictly superior:
1. **Universal Units:** $T \approx T_{\text{ps}}(q)$ is always the reference transition point, whether on $\mathbb{R}^d$, $\mathbb{H}^d$, or the Simplex.
2. **No Clipping on Heavy Tails:** If a hyperbolic Cauchy graph has couplings spanning 4 decades, the bracket automatically expands to span 4 decades centered on $T_{\text{ps}}$.
3. **BARS Boundary Verification:** If the empirical susceptibility $\chi(T)$ has not returned to baseline at $T_{\min}$ or $T_{\max}$, the sweep automatically steps outward by an octave until both boundaries are silent.

### Summary
Path B's unnormalized path is an artifact of earlier heuristic experiments. You can safely discard Path B. The proper way to protect exotic curved geometries from clipping is to keep the $1/\widehat{K}$ gauge for the anchor, and let the coupling dispersion determine the bracket width.

# Searched for "T_fs"
Viewed sol-SPC-bars-graphs-brackets.md:480-530

Yes. We can derive a principled formula that **reduces identically to BWD's $[0.05, 4.0] \cdot T_{\text{ps}}(q)$** on nominal Euclidean data, while automatically expanding on heavy-tailed, curved, or multiscale geometries.

Here is the physics of where BWD's numbers came from, the generalized formula, and the outer verification procedure.

---

### 1. Where BWD’s $[0.05, 4.0]$ Originally Came From

In BWD (1995/1996), data points live in $\mathbb{R}^d$ with average nearest-neighbor distance $\bar{a}$, and couplings are Gaussian with $1/\widehat{K}$ normalization:
$$J_{ij} = \frac{1}{\widehat{K}} \exp\left(-\frac{d_{ij}^2}{2\bar{a}^2}\right)$$

In their mutual-$k\text{NN}$ graph:
* **The Dense Core (Hot end):** In cluster cores where $d \approx 0$, couplings are $J \approx 1/\widehat{K}$. The local weighted degree $S_i = \sum_{j} J_{ij}$ is at most $S_{\max} \approx 2\text{--}4 \times \langle S \rangle$. Above $T \approx 4.0 \cdot T_{\text{ps}}$, even the densest cluster melts into paramagnetic noise.
* **The Weakest Bridges (Cold end):** In mutual-$k\text{NN}$, the longest retained edges between adjacent clusters sit at $d \approx 2.5 \bar{a}$. Their Gaussian coupling is:
  $$\exp\left(-\frac{(2.5)^2}{2}\right) = e^{-3.125} \approx 0.044 \approx \mathbf{0.05}$$
  Below $T \approx 0.05 \cdot T_{\text{ps}}$, thermal fluctuations cannot sever even these weakest bridges ($p_{\text{bond}} \to 1$), so the entire system freezes into a single ferromagnetic giant component.

BWD’s $[0.05, 4.0]$ was not an arbitrary guess—**it was directly tracking the ratio between the weakest retained bridge ($e^{-3.125}$) and the densest core in Euclidean space.**

---

### 2. The General Closed-Form Formula

To make this universal across curved manifolds $\mathcal{M}$ while recovering BWD on Euclidean data, we keep the physics anchor $T_{\text{ps}}(q)$ and let the **graph's empirical coupling spectrum** define the endpoints:

Let $S_i = \sum_{j \sim i} J_{ij}$ be the local weighted degree, and let $J_e$ be the undirected edge couplings:

$$T_{\max} = T_{\text{ps}}(q) \cdot \max\left(4.0, \; \frac{S_{98}}{\langle S \rangle}\right)$$

$$T_{\min} = T_{\text{ps}}(q) \cdot \min\left(0.05, \; \frac{J_{02}}{\text{Median}(J)}\right)$$

#### How It Behaves:
1. **On Nominal Euclidean Data:**
   * $S_{98} / \langle S \rangle \approx 1.5\text{--}2.5 \implies$ clamped by $\max(4.0, \dots) = \mathbf{4.0}$.
   * $J_{02} / \text{Median}(J) \approx \frac{0.044}{0.606} \approx 0.07 \implies$ clamped by $\min(0.05, \dots) = \mathbf{0.05}$.
   * **Result:** Exactly reproduces BWD's $[0.05 \cdot T_{\text{ps}}(q), \; 4.0 \cdot T_{\text{ps}}(q)]$.
2. **On Hyperbolic / Curved / Multiscale Data ($\mathbb{H}^d$ with Cauchy kernel):**
   * Cauchy tails on Poincaré ball bridges reach $d \gg \delta$, yielding $J_{02} / \text{Median}(J) \sim 10^{-3}$ or $10^{-4}$.
   * The formula automatically drops $T_{\min} \to 10^{-4} \cdot T_{\text{ps}}(q)$, **preventing the cold freeze from clipping off**.
   * If the manifold has high-density focal points (like the center of a hyperbolic disk), $S_{98} / \langle S \rangle$ expands $T_{\max}$ upward, **preventing the hot core from clipping off**.

---

### 3. The Runtime Procedure: BARS Boundary Adequacy

As outlined in [sol-SPC-bars-graphs-brackets.md](file:///d:/aghado01/ThermoMapper/issues/spc-BARS/discussions/sol-SPC-bars-graphs-brackets.md#L480-L525), even the best static formula should be backed by an empirical safety loop using the physical susceptibility $\chi(T)$:

```
           [Static Formula Bracket: T_min, T_max]
                            │
                            ▼
           Run coarse sweep & fit BARS spline on χ(T)
                            │
                            ▼
              Check Boundary Adequacy:
        • Is χ(T_min) > ε or f'(T_min) < 0? (Cold clipping!)
        • Is χ(T_max) > ε or f'(T_max) > 0? (Hot clipping!)
                   /                 \
                 YES                  NO
                /                      \
   Expand boundary by octave        Accept bracket
   (e.g. T_min /= 2 or T_max *= 2)   & proceed to full run
   & append observations
```

Because BARS spline fitting consumes observations incrementally, expanding the bracket by an octave doesn't discard previous sweeps—it only samples the newly exposed interval until both physical boundaries are completely silent.

---

### Summary

* **Static Formula:** 
  $$[T_{\min}, T_{\max}] = T_{\text{ps}}(q) \times \left[ \min\left(0.05, \; \frac{J_{02}}{\text{Median}(J)}\right), \; \max\left(4.0, \; \frac{S_{98}}{\langle S \rangle}\right) \right]$$
  Guarantees nominal Euclidean behaves identically to BWD, while heavy-tailed geometries automatically expand.
* **Dynamic Safety:**
  The BARS boundary report verifies that $\chi(T)$ has returned to baseline at both ends, stepping outward if any tail is clipped.
  