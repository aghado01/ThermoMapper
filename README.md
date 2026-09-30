# ThermoMapper

Research code for clustering, graph geometry, and topological data analysis in C#.

![EyeTorusToy rendered by viz-core: gold upper lid, cyan iris, slate lower lid, and violet pupil](docs/images/eye-human.png)

**EyeTorusToy.** A 3D toroidal interpretation of the toy data in
[Blatt, Wiseman & Domany (1996)](https://doi.org/10.1103/PhysRevLett.76.3251).
The live editor's eye preset, with full tube depth and a solid pupil, seed 42.
Colors identify the generator's components.

[More fixture screenshots and how to reproduce them](docs/gallery.md).

## Code

- [Clustering](src/clustering/) — superparamagnetic clustering (SPC), HDBSCAN,
  Gaussian mixtures, k-means, and dendrograms.
- [Graphs](src/graphs/) — construction, distances, couplings, and observables.
- [Topology](src/tda/) — persistent homology, Mapper, and dimensionality reduction.
- [Synthetic data](src/synthetic/) — Euclidean and manifold fixtures.
- [Visualization](src/viz-core/) — the Three.js viewer used above.
- [R oracles](r/) and [Lean](lean/) — reference calculations and formal work.

Source lives in `src/`; assembly definitions live in `projects/`; tests live in
`tests/`. The visualization replacement has its own [status note](src/viz/README.md).

## Open the eye editor

Requires the .NET 10 SDK. From the repository root:

```sh
dotnet run --project projects/VizApi/VizApi.csproj -c Release -- --urls http://127.0.0.1:8766
```

Open [the editor](http://127.0.0.1:8766/eye). Set each component's radii,
position, point count, and color. The lids also have angle and arc-span controls.
Give each object its own spin, brightness, glow, and shine. **Photo view** pauses
motion and hides the panels for screenshots.

## Render the fixtures

```sh
dotnet run --project tests/VizCoreSmoke/VizCoreSmoke.csproj -c Release -- --gallery artifacts/gallery
```

Open `artifacts/gallery/eye-torus.html` in a browser. Drag to orbit; scroll to
zoom. The HTML contains the sampled data and loads Three.js 0.160.0 from a CDN,
so it needs an internet connection.

[Development practices](CONTRIBUTING.md) · [Repository guide](AGENTS.md)
