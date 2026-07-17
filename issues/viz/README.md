# ThermoMapper visualization program

Status: inception / v2 foundation

This directory is the planning room and screenplay library for the visualization
system that will replace `src/viz-core`. The replacement is a research workbench,
not merely a demo viewer: it runs real ThermoMapper producers, preserves their
artifacts and provenance, supports controlled interventions, and presents those
artifacts through reusable linked views.

`src/viz-core` and `projects/VizCore` remain the reference implementation during
the migration. They are not the architectural base of v2 and should receive only
necessary maintenance while the replacement proves its own vertical slices.

## Documents

- [program.md](program.md) — north star, boundaries, subsystem design, and
  non-negotiable invariants.
- [migration-roadmap.md](migration-roadmap.md) — staged delivery plan and the
  objective gates for sunsetting `viz-core`.
- [cinematic-universe.md](cinematic-universe.md) — the broader research and
  explanatory film slate, organized by scientific question.
- [screenplay-format.md](screenplay-format.md) — the authoring grammar for takes,
  shots, cues, interactions, and epistemic disclosure.
- [scripts/eye-torus-pilot.md](scripts/eye-torus-pilot.md) — the first concrete
  film and v2 vertical-slice specification.
- [sol-viz-v2-filmography.md](sol-viz-v2-filmography.md) — source conversation and
  early design material. It is retained as input, not treated as an authoritative
  specification.

## Code foundation

New implementation lives under `src/viz` and is compiled by projects under
`projects/Viz.*`:

```text
src/viz/
  contracts/       neutral study, artifact, provenance, and panel contracts
  scene/           renderer-facing visual scene and interaction state
```

The first assemblies intentionally contain no references to Synthetic, Graphs,
Clustering, TDA, or the browser host. Producer-specific adapters and transport
will be added only when a vertical slice requires them.

## Working rules

1. Scientific algorithms execute in their owning C# modules. Visualization code
   invokes public APIs and adapts returned artifacts; it does not reimplement the
   science.
2. The renderer consumes a visual scene. It does not own experiment, algorithm,
   or oracle state.
3. Feature coordinates, display coordinates, and their projection provenance are
   distinct artifacts.
4. Ground truth is evidence with an explicit role. It never enters a canonical
   method run unless a recorded oracle intervention requests it.
5. Scientific time, pipeline order, and cinematic time are independent axes.
6. Static replay and live execution consume the same versioned study package.
7. Domain artifacts preserve identity across the seam: sample IDs, CSR slots,
   Mapper membership, persistence back-references, and cross-frame matches must
   not be flattened away.
8. The v2 cut introduces no compatibility façade over the old `VizDataset`
   schema. Necessary content is migrated deliberately and obsolete surfaces are
   deleted at sunset.

## Immediate checkpoint

The current milestone is the contract-and-scene foundation followed by the
EyeTorus graph laboratory. It is complete only when the same study can be opened
as a static package and driven by a live host without scientific recomputation in
the browser.
