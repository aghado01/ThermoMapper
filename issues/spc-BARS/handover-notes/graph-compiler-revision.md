# Graph compiler revision: hand-over note (2026-11-01)

Full report: `D:\aghado01\ThermoMapper\issues\spc-BARS\graph-compiler-brief.md`, §4.2 (decisions) and "Report". Design and knob register: `docs/Graphs/GraphConstruction.md`. Commits `7a35ffd`..`dfcb817` on main.

What other sessions sharing this checkout should know:
- **Repair edges are never cut by the coupling floor** and LMP leaves them their kernel coupling. BlattHierarchy under a Gaussian kernel is now one component (it was six). A Gaussian bridge still activates only near T ≈ its coupling, so SW at a reachable T sees what it saw before.
- **The gauge is its own field** (`CouplingProjection.Gauge`, applied last). `MeanEdgeDistance` no longer implies the 1/K̂ factor; the three BWD presets and `BwdParityTests` now set `gauge` and `scale_edges: after-repair`. Their graphs and expected values are unchanged.
- **Mean-edge bandwidth reads the pre-repair edges by default** (`ScaleEdges`); a graph built with `mean-edge-distance` and no preset moves.
- **`auto` temperature grids** come from `SpcScheduleHelpers.AutoBracket`, reading the build's `Resolved` declaration. Cauchy or LMP graphs with the mean-edge bandwidth no longer take the q-only T_ps bracket. Mapper patches keep the heuristic.
- **Every cached graph rebuilds** (the config fingerprint covers the new projection fields). Older run manifests with `mean-edge-distance` re-materialize without the gauge in `extract`.
- **LMP default is still the strict rank**, which writes 0 on each node's weakest edge and fragments sparse graphs; `LmpRank = Inclusive` is the option. Flipping the default is the user's call.
- `GraphBuildResult.Resolved` and `GraphConstructionManifest.Resolved` record what each `Auto` picked and the units of the returned graph.
- `docs/Clustering/SpcEngine.md` line 163 still says the REPL's `auto` reads `MeanEdgeDistance`; it was left because another session has uncommitted edits in that file.
- Nothing in `src/maths`, `issues/regression` or `stopping-theory-primer.md` was touched.
