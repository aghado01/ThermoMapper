# Viz v2 migration roadmap

## Strategy

Use a strangler migration with independent v2 contracts. `viz-core` remains
runnable while v2 grows through complete vertical slices, but no compatibility
surface is placed between `VizDataset` and `VizStudy`. Useful rendering ideas may
be reimplemented or extracted deliberately; the old schema and monolithic viewer
do not dictate the new design.

The migration has three rules:

1. Do not break existing demos merely to create an empty v2 seam.
2. Do not add new scientific capabilities to `viz-core` when they belong in v2.
3. Delete the old projects wholesale only after objective replacement gates pass.

## Phase 0 — foundation

Deliver:

- `Viz.Contracts`, with stable IDs, evidence roles, provenance, entities,
  coordinate spaces, artifacts, relations, runs, panels, and links; the typed
  presentation-sequence contract is deferred until two scripts stabilize the
  grammar (see the screenplay format's vocabulary boundary);
- `Viz.Scene`, with declarative marks, stable entity bindings, camera/view state,
  typed view actions, and renderer/compiler interfaces;
- the planning and screenplay library under `issues/viz`;
- focused builds proving both assemblies are neutral.

Exit gate: both projects compile without references to ThermoMapper domain
assemblies or browser technology.

## Phase 1 — EyeTorus graph laboratory

Deliver:

- a richer synthetic-case seam for EyeTorus oracle anatomy, latent sample IDs,
  realization provenance, and cross-realization correspondence;
- a Synthetic adapter and Graphs adapter;
- two complete graph artifacts and graph-health reports;
- a static study manifest with inline development payloads or initial binary
  chunks;
- a minimal spatial renderer with point/edge layers, picking, visibility, and
  linked A/B selection;
- the first executable version of `scripts/eye-torus-pilot.md`.

Exit gate: the user can inspect the same EyeTorus samples under two graph
constructions, select a scientific edge by preserved CSR slot, and toggle an
oracle cross-group classification without changing the graphs.

## Phase 2 — package and host parity

Deliver:

- versioned manifest and typed chunk format, composed on Archivory's chunk
  machinery — the chunk codec lands in Archivory first, and the study package
  is its browser-readable container profile (spec: MarkBrain vault,
  `ThermoMapper/issues/archivory/format-spec.md`);
- checksums, schema validation, lazy chunk loading, and package fingerprints;
- directory/static, embedded, and live artifact sources;
- asynchronous C# package writer;
- live job protocol with progress, cancellation, request identity, and stale
  result rejection;
- fully local browser dependencies for offline export, with chunk codecs
  restricted to what the browser decodes without a server (raw or gzip;
  brotli stays available to backend-only profiles).

Exit gate: a live-generated study can be saved, reopened without the server, and
produce the same artifact inventory and visible initial state.

## Phase 3 — mature spatial engine

Deliver:

- packed point, segment, polyline, triangle, vector, line-field, ellipsoid, and
  annotation marks;
- centralized GPU resource ownership and incremental scene updates;
- GPU or equivalently scalable stable-ID picking;
- lasso, probing, measurement, clipping, and selection framing;
- renderer capability negotiation;
- native Poincaré boundary and geodesic rendering;
- explicit 4D-to-3D projection artifacts and controls.

Exit gate: EyeTorus Euclidean/Poincaré/4D takes run without special-purpose code
inside the render kernel, and repeated study replacement shows no resource leak.

## Phase 4 — coordinated research workbench

Deliver:

- spatial, chart, barcode, hierarchy, report, and table panels;
- shared selections and scientific-axis cursors;
- artifact outline and contextual inspector;
- view history, undo/redo, saved presets, and comparison layouts;
- provenance trail and visible evidence-role disclosure;
- film playback that the user may pause, inspect, and leave without entering a
  separate application mode.

Exit gate: an SPC temperature cursor can coordinate spatial assignments, edge
currencies, curves, and lineage while keeping cinematic time independent.

## Phase 5 — orchestration and producer breadth

Deliver typed recipes and adapters in evidence-driven order:

1. Synthetic and Graphs.
2. Graph fields and spectral diagnostics.
3. SPC sweeps, landscapes, hierarchies, and lineage.
4. Mapper, PH, filtrations, barcodes, and nerve diffs.
5. GMM, HDBSCAN, mixtures, responsibilities, dendrograms, and condensed trees.
6. Additional producers as concrete films require them.

Recipe execution adds `FreezeArtifact`, `Fork`, `Sweep`, `Compare`,
`ReplaceInput`, and `BypassStage` only with typed validity rules and complete
provenance.

Exit gate: at least one canonical, partial, comparative, and oracle-intervention
recipe can be reproduced from a saved study.

## Phase 6 — `viz-core` sunset

Delete `src/viz-core`, `projects/VizCore`, the old VizApi composition path, and
their dedicated smoke/test projects only when all of the following are true:

- every still-valued existing demo has a v2 film, preset, or documented reason
  for retirement;
- v2 covers the existing spatial mark vocabulary;
- serializer/viewer contract checks and browser tests are green;
- temporal data is synchronized through scientific axes rather than the old
  label-only temporal sequence;
- static export is genuinely offline and packages large studies without a giant
  JSON copy;
- live regeneration has cancellation, progress, and stale-response protection;
- graph identity, including CSR slots and edge currencies, survives round-trip;
- feature/display coordinate separation and projection provenance are enforced;
- native hyperbolic rendering is available;
- `LocalTangent` and other scientific computation no longer live in viz;
- operational documentation identifies v2 as the only supported path.

No `[Obsolete]` aliases or long-lived forwarding shims remain after the cut.

## Risk controls

- **Architecture without feedback:** each phase ends in a visible vertical slice,
  not a collection of unused abstractions.
- **Browser framework churn:** keep `Viz.Contracts` and `Viz.Scene` independent of
  the workbench framework.
- **Over-general contracts:** add semantic carriers when a second real producer
  demonstrates the common shape.
- **Scientific leakage:** domain adapters are one-way and JavaScript never grows
  substitute implementations.
- **Beautiful but false films:** every cue declares its timeline and evidence
  semantics; presentation-derived thinning and interpolation remain visible in
  provenance.
- **Performance surprises:** preserve chunks and incremental updates from the
  start, then test with canonical fixture sizes rather than only miniature scenes.

## Current checkpoint

Phase 0 is active. The initial `Viz.Contracts` and `Viz.Scene` projects establish
the dependency direction but are expected to evolve when the EyeTorus vertical
slice exercises them. Phase 1 is the next implementation target; later phases are
planning, not promises that their current decomposition is final.
