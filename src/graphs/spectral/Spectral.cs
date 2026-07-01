// src/graphs/spectral/Spectral.cs
#nullable enable
using System;
using System.Buffers;
using System.Collections.Generic;
using Graphs.Primitives;
using Maths.LinAlg;

namespace Graphs.Spectral;

public enum SolverKind { Dense, Iterative, Auto }

/// <summary>
/// Graph-spectral analysis orchestrator. Materialises the graph Laplacian via
/// <see cref="GraphLaplacian"/> and delegates dense eigendecomposition +
/// bottom-K selection to <see cref="SpectralMath"/> in <c>Maths.LinAlg</c>.
/// This class owns the CSR-to-dense seam; <see cref="SpectralMath"/> owns the
/// graph-free spectral math.
///
/// "Bottom" is by eigenvalue ascending; eigenvalue 0 corresponds to the
/// constant mode (one per connected component) and is included in the
/// returned list. Callers that want only non-trivial modes can drop entries
/// where <c>Lambda &lt; ε</c>.
///
/// <see cref="SolverKind.Dense"/> materialises the full Laplacian and runs a
/// dense decomposition; <see cref="SolverKind.Iterative"/> runs the matrix-free
/// <see cref="Maths.LinAlg.LOBPCG"/> primitive (see
/// <see cref="GraphSpectral.ComputeBottomK"/>) and honours the <c>seed</c>.
/// <see cref="SolverKind.Auto"/> resolves to one of those via
/// <see cref="SpectralSolverPolicy"/> from the request shape <c>(n, k)</c>.
/// All kinds return the same trivial-inclusive bottom-K set.
/// </summary>
public static class Spectral
{
    public static IReadOnlyList<EigenPair> ComputeBottomK(
        CsrGraph graph,
        int seed = 0,
        int k = 8,
        LaplacianType lapType = LaplacianType.Combinatorial,
        SolverKind solverKind = SolverKind.Dense,
        DenseEigenOptions denseOptions = default,
        DenseLaplacianMaterialization denseMaterialization = DenseLaplacianMaterialization.Rectangular)
    {
        int n = graph.NodeCount;

        SolverKind resolved = solverKind == SolverKind.Auto
            ? SpectralSolverPolicy.Resolve(n, k)
            : solverKind;

        if (resolved == SolverKind.Iterative)
            return GraphSpectral.ComputeBottomK(
                graph, k: k, lapType: lapType, seed: seed, deflateNullSpace: false).Eigenpairs;

        if (n == 0 || k <= 0)
            return Array.Empty<EigenPair>();

        if (denseMaterialization == DenseLaplacianMaterialization.FlatColumnMajor)
        {
            int matrixSize = n * n;
            double[] flatLaplacian = ArrayPool<double>.Shared.Rent(matrixSize);
            try
            {
                GraphLaplacian.BuildDenseColumnMajor(graph, lapType, flatLaplacian.AsSpan(0, matrixSize));
                return SpectralMath.BottomK(flatLaplacian.AsSpan(0, matrixSize), n, k, denseOptions);
            }
            finally
            {
                ArrayPool<double>.Shared.Return(flatLaplacian);
            }
        }

        double[,] laplacian = GraphLaplacian.BuildDense(graph, lapType);
        return SpectralMath.BottomK(laplacian, k, denseOptions);
    }
}
