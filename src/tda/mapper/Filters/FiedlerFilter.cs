// ============================================================================
// TDA.Mapper — FiedlerFilter.cs
// ============================================================================
// Fiedler-vector filter implementation. Separated from GraphFilters.cs because
// it carries non-trivial dependencies — Graphs.Spectral (Laplacian materialiser),
// Maths.LinAlg (dense eigendecomposition), and Graphs.Diagnostics (connectivity
// validation) — mirroring how Pca1DFilter is split out from DefaultFilters.cs.
// ============================================================================

#nullable enable
using System;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.LinAlg;
using TDA.Mapper;

namespace TDA.Mapper.Filters;

/// <summary>
/// Computes the Fiedler vector (eigenvector of the second-smallest eigenvalue)
/// of the weighted graph Laplacian L = D − W, where D is the diagonal of
/// weighted degrees and W is the symmetric edge-weight matrix.
///
/// Mathematical preliminaries: for a connected weighted graph, L is positive
/// semi-definite. Its smallest eigenvalue is 0 (eigenvector = constant). The
/// Fiedler value (second-smallest eigenvalue) measures algebraic connectivity;
/// its eigenvector (the Fiedler vector) partitions the graph along the
/// dominant connectivity axis.
///
/// Implementation: delegates to <c>Spectral.ComputeBottomK(k: 2, SolverKind.Auto)</c>
/// on the combinatorial Laplacian (non-positive edge weights clamped to 1.0) and
/// takes the second (Fiedler) mode. <see cref="SolverKind.Auto"/> keeps small
/// graphs on the dense path and routes large ones to the matrix-free LOBPCG
/// solver, so this is no longer a fixed O(N³)/O(N²) full decomposition computed
/// just to keep one vector.
/// </summary>
internal sealed class FiedlerVectorFilter : IGraphFilter
{
    public string Name => "Fiedler vector (Laplacian PC2, dense)";

    public double[] Apply(CsrGraph graph, double[][]? features = null)
    {
        int n = graph.NodeCount;
        if (n == 0) return Array.Empty<double>();
        if (n == 1) return new double[] { 0.0 };

        // Validate connectivity. Disconnected graphs have multi-dimensional
        // null space (one constant eigenvector per component); "the Fiedler
        // vector" is ill-defined.
        var diag = Connectivity.Validate(graph);
        if (diag.ComponentCount > 1)
            throw new InvalidOperationException(
                $"FiedlerVectorFilter requires a connected graph; found {diag.ComponentCount} components. " +
                "Use ConnectivityRepair.EnsureConnected to repair connectivity via MST bridges before computing Fiedler.");

        // Bottom-2 modes of the combinatorial Laplacian: pair[0] is the trivial
        // constant mode (λ≈0), pair[1] is the Fiedler vector. Routed through the
        // orchestrator so SolverKind.Auto can hand large graphs to matrix-free
        // LOBPCG instead of a full dense decomposition for a single vector.
        var pairs = Spectral.ComputeBottomK(
            graph,
            seed: 0,
            k: 2,
            lapType: LaplacianType.Combinatorial,
            solverKind: SolverKind.Auto);

        if (pairs.Count < 2)
            throw new InvalidOperationException(
                $"Fiedler extraction requires 2 eigenpairs; got {pairs.Count}.");

        // EigenPair.Vector is a freshly allocated copy — safe to hand to the caller.
        return pairs[1].Vector;
    }
}
