# Visualization v2 program

## North star

ThermoMapper Viz is an experimental composition and inspection environment over
the real ThermoMapper modules. It supports:

- canonical, provenance-complete method runs;
- partial execution and inspection of intermediate artifacts;
- comparative sweeps and linked A/B views;
- frozen-stage branches and controlled oracle substitutions;
- diagnostic investigation of why a result occurred;
- exploratory development of methodological intuition;
- curated demonstrations and didactic sequences using the same artifacts.

These are three epistemic postures of one workbench:

- **Exploratory:** what changes when an assumption or intervention changes?
- **Diagnostic:** why did this result occur, and at which stage?
- **Demonstrative:** what does the pinned method do, and how can it be explained?

Didactic value is a presentation purpose, not a fourth execution mode.

## System boundary

The complete system has two halves which must remain separable.

### Research orchestration

The orchestrator invokes public APIs from Synthetic, Graphs, Clustering, TDA,
and future producers. It owns typed recipes, jobs, cancellation, provenance,
branching, sweeps, and artifact registration. It may execute stages partially or
in a different valid order when a recorded experiment calls for it.

Its intervention vocabulary is deliberately constrained:

- `ReplaceInput`
- `BypassStage`
- `FreezeArtifact`
- `Fork`
- `Sweep`
- `Compare`

This should become a typed operation graph, not an unrestricted visual node
editor. Operation ports and artifacts enforce scientifically meaningful
composition.

### Visualization workbench

The workbench compiles registered artifacts into visual scenes, renders linked
panels, accepts direct manipulation, and issues typed actions. It never hides a
scientific recomputation inside a color, camera, or DOM callback.

```text
domain APIs -> recipes/jobs -> study artifacts -> scene compiler -> renderers
                     ^                |                 ^             |
                     |                v                 |             v
                 interventions   linked selection   view actions   pixels/UI
```

## The study model

`VizStudy` is the durable research document. It can contain:

- **entity sets:** samples, graph nodes, Mapper nodes, Gaussian components,
  persistence bars, simplices, hierarchy nodes;
- **coordinate spaces:** feature, display, Euclidean, Poincaré-ball, spherical,
  simplex, projected, and manifold-local coordinates;
- **artifacts:** graphs, fields, assignments, mixtures, complexes, barcodes,
  hierarchies, series, reports, and domain-specific extensions;
- **relations:** cluster membership, Mapper-node membership, bar-to-generator,
  sample correspondence, node matching, and lineage;
- **panels:** spatial, chart, barcode, hierarchy, report, and table views;
- **links:** shared selection, axis cursor, camera, filtering, and hover policy;
- **presets and films:** saved research views and curated sequences;
- **provenance:** producer identity, configuration, inputs, fingerprints,
  warnings, and evidence role. Configuration is recorded as a nested config
  trace with requested/resolved/source per auto-filled knob (the archivory
  provenance discipline), never flattened into a parameter list.

The core uses mathematical carrier families rather than algorithm names:

- anchored categorical, scalar, vector, line, and edge fields;
- graphs, complexes, mixtures, hierarchies, barcodes, series, and reports;
- indexed artifacts over a scientific axis;
- relations between entity sets.

Algorithm-specific meaning is expressed through extensible semantic IDs and
producer adapters, not an ever-growing core enum.

## Evidence and interventions

Every artifact has an evidence role:

- `Observed` — produced by the canonical method from admissible inputs;
- `Oracle` — ground truth or privileged knowledge supplied by a fixture or
  external reference;
- `Counterfactual` — produced after a recorded intervention;
- `DiagnosticDerived` — analysis of other artifacts used to explain behavior;
- `PresentationDerived` — projection, thinning, interpolation, or other view
  preparation that makes no scientific claim.

Ground truth is not a global mode. It may be used as:

- a reference overlay;
- a classifier of an observed result, such as oracle cross-group graph edges;
- an upstream substitution in a counterfactual branch;
- calibration for an empirical field;
- one side of an oracle-versus-observed explanation.

An oracle artifact remains available to the study while being invisible and
inadmissible to a canonical recipe stage. This separation must be enforced above
the renderer.

## Client architecture

The mature client is divided into the following systems.

### Render kernel

Owns the Three.js/WebGL lifecycle, GPU resources, render scheduling, cameras,
picking buffers, resizing, instancing, and disposal. It renders only when dirty
except during animation.

### Visual scene

A declarative retained model of coordinate bindings, marks, groups, styles,
visibility, and stable entity references. Geometry and appearance are separate
so selection or recoloring does not rebuild topology.

### Interaction services

Picking, hover, selection, lasso, brushing, measurement, probing, clipping, and
camera manipulation. Interactions emit typed actions rather than mutating
Three.js objects throughout UI callbacks.

### View controllers

Spatial, chart, barcode, hierarchy, report, and table panels. Each controller
translates study artifacts and view state into a renderer-specific presentation.

### Workbench shell

Artifact outline, contextual inspector, recipes/jobs, linked panel layout,
provenance trail, history, saved views, and film playback.

## State separation

- **Scientific/document state:** artifacts, recipes, branches, provenance.
- **View state:** cameras, layer visibility, encodings, clipping, axis cursor,
  linked-selection policy.
- **Transient state:** hover, open tooltip, drag, provisional lasso, animation
  interpolation.

Scientific and durable view state can be saved. Transient state cannot leak into
the study package.

## Renderer design rules

- Build a narrow visual engine around ThermoMapper's required marks rather than
  wrapping all of Three.js.
- Use stable scientific entity references for picking and linked selection.
- Preserve explicit coordinate-space and projection provenance.
- Use packed buffers and instancing for large point, edge, arrow, and glyph sets.
- Update scene resources incrementally and centralize ownership/disposal.
- Render native geometry where it changes interpretation, including Poincaré
  boundaries and geodesic arcs. Geometry-model-faithful drawing math (geodesic
  interpolation, boundary construction) is renderer competence — the analogue
  of drawing a straight Euclidean line — not scientific computation.
- Keep legends, inspectors, controls, and ordinary annotations in the UI layer;
  use WebGL text only when spatial anchoring requires it.
- Advertise renderer capabilities instead of scattering geometry-specific
  conditionals through every layer.

## Transport

Transport composes on `Archivory`, the repository's centralized persistence
machinery — chunk encoding, atomic writes, artifact scopes, JSON artifact
conventions, and config-trace provenance. The study package is one Archivory
container profile: a JSON manifest plus independently compressed chunks, chosen
so a browser can read it directly and load lazily. Other consumers (REPL swap,
checkpoint archives) use different container and compression profiles of the
same chunk machinery; the profiles share standards, not a single layout.

The durable package will be a versioned manifest plus typed chunks:

- dtype, shape, byte order, compression, checksum, and semantic role are explicit;
- CSR row offsets, columns, weights, and slots remain structurally identifiable;
- temporal frames and large fields may be loaded independently;
- embedded, directory/static, and live sources implement the same artifact-source
  contract;
- package writing is asynchronous and chunk-oriented, never one giant JSON
  string copied into an HTML template;
- chunks are independently decompressible, so a partial load seeks without
  whole-stream decompression;
- chunk locations are package-relative; checksums identify chunks, and package
  fingerprints exclude location;
- the client rejects unsupported schema versions rather than guessing.

Offline means offline: a self-contained export cannot depend on a CDN.

## Assembly direction

The planned assemblies are introduced only as needed:

```text
Viz.Contracts       BCL-only durable research contracts
      |
Viz.Scene           renderer-facing scene and interaction model
      |
Viz.Transport       manifests, chunks, readers, writers, version validation
      |             (composes on Archivory, the repository persistence layer)
Viz.Web             browser renderer and workbench client assets

Viz.Adapters.*      producer-specific adapters into Viz.Contracts
Viz.Orchestration   typed recipes, jobs, branches, and interventions
Viz.Host            live/static application composition
```

`Viz.Adapters.Graphs` may reference Graphs; `Viz.Adapters.Tda` may reference TDA;
and so on. Those references never point back into the domain modules.
`Viz.Transport` references `Archivory` — itself BCL-only — rather than defining
a bespoke format. The contracts and scene assemblies stay neutral.

## First vertical slice

The EyeTorus graph laboratory is the architectural proof:

1. Generate a seeded EyeTorus case with explicit sample identities, feature and
   display coordinates, oracle hierarchy, and generative provenance.
2. Run two real graph recipes from the Graphs public API.
3. Preserve each complete `GraphBuildResult`, including CSR slots, metric and
   weight semantics, manifest, repair information, health report, and logs.
4. Derive an oracle cross-group edge classification without modifying either
   graph.
5. Present linked A/B spatial panels plus graph diagnostics and an inspector.
6. Support linked selection from points and edges into provenance and reports.
7. Render the Poincaré variant with a boundary and native geodesic arcs.
8. Open the exact same study from a static package and a live host.

No scientific computation may be added to JavaScript to make this slice work.

## Non-goals for the foundation

- Feature parity with `viewer.html`.
- A generic visual programming environment.
- A comprehensive enum of every future ThermoMapper artifact.
- A framework decision for the final workbench UI before the artifact and
  interaction contracts are exercised.
- Migrating `LocalTangent` into the new Viz tree; scientific producers must move
  to their owning maths or graph-analysis modules.
- A compatibility adapter that makes `VizDataset` appear to be `VizStudy`.

## Quality gates

- Contract round-trip tests and schema-version rejection.
- A study-level admissibility audit: no `Observed` artifact transitively
  consumes an `Oracle` input through its provenance chain.
- Browser contract and screenshot tests for representative scenes.
- Stable entity-picking tests across incremental updates.
- Resource lifecycle tests for repeated scene replacement.
- Static/live parity tests against the same package fingerprint.
- Scientific fixture tests remain in producer modules; visualization tests assert
  faithful adaptation and representation, not the producer's mathematics.
- Every film declares whether a transition represents computation, a scientific
  axis, or presentation only.
