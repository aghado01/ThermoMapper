# Graph settings behind round 2's SPC results

From the stopping round-2 session, 2026-11-02. The record is `issues/spc-BARS/stopping-theory-primer.md` §13.17 and C58. Harness and raw output: `graph-settings-round2/` beside this file (`Program.cs`, `graphline.csproj`, `graphline.out`). No SPC ran; graphs only, on the current compiler.

The user's rule: "these experiments that involve running SPC aren't taking into account what the settings on the graph are, these things needs to be made explicit so at least we know and can consider if there settings used might explain the results. default settings != correct settings".

## What any session running SPC needs

- **State the graph as resolved**, from `GraphBuildResult.Resolved`, in one line. Example:
  `iris 150 × 4, Euclidean, raw; KNN k = 5, mutual (skew 0.18), MST repair (15 edges, largest 0.527 before), PathNeighbor; Gaussian a = 0.307 (mean edge), g = 0.286 (1/K̂, after repair), no LMP, floor 1e-10; K̂ = 3.49, s̄ = 0.616 (0.151, 1.126), one component, weakest bridge J = 1.7e-7; q = 20, budget, grid, seeds, commit.`
- **Refinement has no effect under Euclidean.** `Auto` → `PathNeighbor` takes min(edge, shortest path). Measured difference from the true distance: at most 5.6e-17 on all five graphs, repair edges included.
- **The recipe pins mean strength.** With g = 1/K̂ and a = the mean edge length, s̄ = ⟨κ⟩ ≈ e^(−1/2) = 0.607 on any data at any K (0.603–0.643 measured). So agreement of T_ps or the hot end across datasets on this recipe counts as one check.
- **Compiler defaults are a different graph:** OR filter, MAD bandwidth (2.5–3.5× smaller), no gauge, repair and LMP chosen per dataset, s̄ from 0.09 to 2.58, 57 coupled components on the toy and 177 on Landsat.

## What the graph alone fixes on iris

Exactly one edge joins setosa to the rest at k = 5, 10, 20, 40, a repair edge about 1.64 long. J = 1.7e-7, 5.7e-6, 1.2e-4, 9.5e-4. As a bridge it's open with probability p̃ = p/(p + q(1 − p)), at most 2.0e-6, 6.4e-5, 1.3e-3, 1.2e-2, at T = 0.00446. So on the Gaussian recipe the setosa split is present from the cold end, and a resolver can't mix it.

The versicolor–virginica cut: 13, 34, 76, 285 edges, ΣJ 1.68, 2.35, 2.65, 4.38.

## K axis (parity recipe, graph only)

S₂ at T = 0.00446 from the percolation envelopes:

| k | Toy | Iris | Landsat |
|---|---|---|---|
| 5 | [0.631, 0.831] | [0.538, 0.552] | [0.992, 0.998] |
| 10 | [0.838, 0.936] | [0.540, 0.555] | [0.986, 0.997] |
| 20 | [0.949, 0.995] | [0.554, 0.571] | [0.982, 0.998] |
| 40 | [0.987, 1.000] | [0.557, 0.640] | [0.967, 0.997] |
