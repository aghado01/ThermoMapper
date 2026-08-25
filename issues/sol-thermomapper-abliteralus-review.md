ThermoMapper resolves into eight distinct cross-examination territories. The strongest first wave is the neutral graph substrate plus graph observables/spectral diagnostics, paired with metric geometry. It avoids duplicating `feat/grassmann`, does not reopen SPC, and sits beneath many consumers.

## Architectural map

The checkout contains roughly 518 C# source files across 40 production assemblies and 21 test projects. Its intended dependency order is explicit in [AGENTS.md](/D:/aghado01/ThermoMapper/AGENTS.md:24):

```text
Maths.LinAlg / Rng / Distance / Information / Geometry
                         ↓
       Graphs.Primitives → Distance / Proximity / Compiler
                         → Diagnostics / Observables / Spectral
                         ↓
          TDA.Ph     clustering     VizApi     UserRepl
             ↘          ↓             ↙
                    TDA.Mapper
               high-level integration
```

Direct project references reinforce that structure:

- `Maths.LinAlg` is consumed by 12 production assemblies.
- `Graphs.Primitives` is consumed by 11.
- `TDA.Mapper` has 13 dependencies and is therefore an integration assembly, not a primitive.
- SPC and `UserRepl.Commands` each have 15 dependencies and sit firmly in the application/control plane.

ThermoMapper’s production algorithms are also largely locally implemented rather than wrappers around a large scientific framework. That makes the later C#→PyTorch work an algorithm-and-contract translation, not merely an API substitution.

## Cross-examination territories

| Territory | Implemented substance | OBLITERATUS questions |
|---|---|---|
| 1. Numerical, metric, and Riemannian foundations | Eigensolvers, LOBPCG, PCA/ICA/FFT, Mahalanobis; Wasserstein, Fisher–Rao, hyperbolic and information distances; robust manifold means and medians | How should weights, activations, subspaces, and surgery deltas be compared? Which operations are tensor-native, batched, differentiable, and device-safe? |
| 2. Neutral graph construction | Typed metrics, kNN/epsilon graphs, mutual filtering, connectivity repair, path refinement, graph projections, serializable `CsrGraph` | What should nodes and edges represent—layers, neurons, directions, prompts, or interventions? How do graph recipe and provenance affect every downstream conclusion? |
| 3. Graph diagnostics and spectral operators | Hubness, neighborhood scale, graph health, typed graph signals, entropy, cycles, Laplacians, effective resistance, coherent fields, magnetic Laplacians | Can these reveal intervention propagation, harmful communities, unstable layers, directional coherence, or graph pathologies before and after surgery? |
| 4. PH kernels and diagram currency | Explicit/lazy/full Rips, homology/cohomology, cycle representatives, Wasserstein diagram metrics, static and dynamic/zigzag PH | Can topology provide stable before/after diagnostics? Is zigzag appropriate for non-monotone surgery passes or layer-depth trajectories? |
| 5. Mapper and persistent Mapper | Scalar/multidimensional lenses, covers, local clustering, graph nerve skeletons, node lineage and merge/split tracking | Could Mapper produce useful atlases of prompts, modules, experts, or weight updates? The current semantics need review before application. |
| 6. Sweep and uncertainty inference | BARS/BAPS, splines, GP regression, changepoints, RJMCMC, parallel tempering, chain diagnostics | Could these characterize response curves across surgery strength, identify regime changes, or quantify uncertainty in intervention selection? |
| 7. Synthetic and rigor infrastructure | Manifold/hierarchical fixtures with retained truth and tangent geometry; R oracles; xUnit parity paths; Lean theorem hierarchy | How can PyTorch ports be falsified on controlled cases, rather than merely demonstrated on a model checkpoint? |
| 8. Visualization, provenance, and supporting tools | Runnable legacy graph/field visualization; emerging evidence-aware Viz v2; Archivory manifests; Hashish fingerprints; repo/test harnesses | Can OBLITERATUS gain reproducible surgery studies, before/after scene packages, direction-field views, artifact identity, or checkpoint/layer fingerprints? |

Representative seams include the declarative [GraphCompiler](/D:/aghado01/ThermoMapper/src/graphs/GraphCompiler.cs:49), durable [CsrGraph](/D:/aghado01/ThermoMapper/src/graphs/primitives/CsrGraph.cs:33), typed [graph-signal contract](/D:/aghado01/ThermoMapper/src/graphs/observables/IGraphSignal.cs:59), [graph Laplacian machinery](/D:/aghado01/ThermoMapper/src/graphs/spectral/GraphLaplacian.cs:36), PH’s [IFiltration](/D:/aghado01/ThermoMapper/src/tda/ph/IFiltration.cs:7), and the evidence-bearing [VizStudy](/D:/aghado01/ThermoMapper/src/viz/contracts/VizStudy.cs:10).

Grassmann/SPRED should remain a calibration case because it is already in flight. SPC’s algorithm should remain deferred, while its control-plane, provenance, and graph-health patterns remain reviewable.

## Important cautions

- PyTorch should not be forced onto every algorithm. Dense linear algebra, distances, and many spectral operations are natural tensor ports; union-find, dynamic trees, filtration enumeration, and exact PH reduction may remain CPU/Python code or require compiled extensions.
- The `issues/hilbert` corpus contains promising heat-semigroup, wavelet, sheaf, and correspondence ideas, but these are research proposals built around a smaller implemented spectral substrate—not already-implemented engines.
- Legacy Viz is runnable; Viz v2 has stronger epistemic/provenance contracts but is still foundational.
- Mapper currently emits a graph 1-skeleton. Its reported loop count is graph cycle rank, not automatically homology of the full cover nerve.
- Persistent Mapper’s correspondence currently selects one best successor despite overlapping Mapper covers. Its output is better treated as event lineage pending a correctness review, not assumed to be established persistence.
- Graph-restricted Rips measures the graph-construction recipe as well as the underlying points. Metric, \(k\), thresholds, filtering, and symmetrization must therefore be first-class experimental provenance.
- R parity tests can return successfully when R is unavailable. A green test exit alone does not establish that an external oracle ran.
- Lean currently formalizes selected cut, symmetrization, stability, and audit claims—not the correctness of the C# PH reducers or Mapper.

## Proposed deep-dive plan

Use a repeating four-role review cell for each bounded wave:

1. Root synthesis: maintain the OBLITERATUS object model, culture, constraints, and cross-territory comparison.
2. ThermoMapper reviewer: inspect implementation, tests, invariants, approximations, complexity, and maturity.
3. Mathematical research agent: verify cited algorithms and novelty against primary sources.
4. PyTorch/evidence reviewer: classify tensor/device/autograd/sparse/quantization implications and design falsifiable parity tests.

Each candidate should produce one compact concept dossier:

- mathematical object and invariants;
- exact ThermoMapper implementation boundary;
- implemented, provisional, and proposed portions;
- possible OBLITERATUS attachment point;
- PyTorch-native, CPU-combinatorial, or extension-backed translation;
- numerical and computational risks;
- existing evidence and missing oracle;
- smallest useful experiment and negative control.

Recommended waves:

1. Graph construction, graph health/observables, and spectral operators; pair with metric semantics.
2. PH kernels, diagram metrics, and zigzag/dynamic persistence.
3. Mapper correctness and Mapper-as-diagnostic research—only after Wave 2 establishes the topology currency.
4. Synthetic fixtures, legacy visualization diagnostics, Viz v2 provenance, and artifact schemas.
5. Regression/sampling, Hashish, non-SPC clustering, and research tooling as secondary opportunity areas.
6. Score candidates by OBLITERATUS relevance, leverage, PyTorch fit, evidence maturity, computational cost, isolation, and conceptual novelty; then authorize only small spikes.

The best next bounded assignment is therefore:

> Cross-examine ThermoMapper’s graph compiler, graph currencies, graph-health observables, and spectral operators against OBLITERATUS’s layer/weight/intervention objects, while a research agent verifies the mathematical lineage and a PyTorch reviewer designs the minimal diagnostic prototype.

This survey was read-only. It reflects `main` at `56c9ccab4f637c9be0bf3760a5c3165cbbd5e602`, two commits ahead of `origin/main`, with an existing dirty working tree containing an active Lean/design reorganization. I did not run full build or test gates, and none of those working-tree changes were modified.