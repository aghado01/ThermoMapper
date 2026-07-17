# Viz v2 source

This tree is the replacement visualization system. It is intentionally separate
from `src/viz-core`, which remains a reference implementation until the migration
gates in `issues/viz/migration-roadmap.md` pass.

Current source areas:

- `contracts/` — neutral durable study, artifact, evidence, coordinate, panel,
  and relation contracts (`Viz.Contracts`).
- `scene/` — renderer-facing visual scene, durable panel view state, typed view
  actions, and renderer lifecycle contract (`Viz.Scene`).

Dependency direction:

```text
Viz.Contracts <- Viz.Scene
```

Neither assembly references Synthetic, Graphs, Clustering, TDA, Three.js, or a
web framework. Producer adapters, orchestration, transport, host, and browser
implementation will be introduced as separate projects when the EyeTorus
vertical slice requires them.

Scientific algorithms do not belong in this tree. Viz invokes producers in their
own modules and preserves the returned artifacts and identities.
