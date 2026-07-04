# AGENTS.md — layering invariants

## The one-way layer order

```
maths  →  graphs  →  { tda, clustering, viz, user-repl }
maths/topology  →  tda/ph  →  tda/mapper  →  tda/pipelines
```

Dependencies point left (downward). `graphs` is a neutral substrate: it emits a
weighted `CsrGraph` + provenance manifest, and SW/SPC/PH/Mapper are *consumers*
of that artifact. Placement authority for the topology stack lives in
`tda-placement.md` (MarkBrain: `ThermoMapper/issues/tda-purification/persistent-homology/`).
There is no `tda/primitives` — everything PH lives in `tda/ph`.

## The construction invariant

Every graph-construction stage (`src/graphs/**`) references only `graphs` and
`maths` namespaces. Consumer and engine semantics enter construction as **data
through typed contracts** — config values, injected delegates, precomputed edge
scores — never as an import of a consumer or PH type.

Worked example (violation cleared 2026-07-03): `GraphCompiler` used to import
the PH engine and run involuted persistence mid-build
(`H1CycleEdges.FromDistanceGraph`) to compute its LMP protect-set. The fix is
the pattern to repeat: `GraphCompiler.Build(..., ProtectedEdgeSource? protectedEdges)`
— the caller (UserRepl commands) computes the H1 protect-set above the seam and
injects it as a delegate over the distance graph. Future topology-scored passes
(the graph sculptor's persistence/spectral criteria) flow in the same way:
scores in, neutral substrate out.

## The litmus test

`CsrGraph` round-trips to disk (`WriteTo`/`FromBinary`). The compiler is
decoupled exactly when a non-SPC reader — a pure PH run, a Mapper nerve, a
spectral embedding — consumes the output unchanged. Needing anything
SPC-specific from construction means fusion has crept back in: re-cut the seam
instead of threading the consumer concern into a build stage.

## Keep the boundary artifacts boring

- The compile boundary emits an immutable, serializable artifact + manifest;
  runs snapshot their full parameter set (config-artifact provenance).
- Configs are declarative, JSON-serializable DTOs. Delegates and live objects
  ride the call, not the config.
- Fluent/chained ergonomics live at the CLI/REPL boundary; the backend consumes
  the pure DTO.
- Pre-release discipline: no back-compat shims, no `[Obsolete]` aliases —
  superseded surfaces are deleted wholesale.
