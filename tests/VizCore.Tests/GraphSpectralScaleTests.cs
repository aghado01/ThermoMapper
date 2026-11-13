using System;
using System.Collections.Generic;
using System.Linq;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The graph callers of the dense eigensolver under rescaled weights (spectral-solver brief, B1 stage three). A Laplacian's spectrum
/// scales with its weights and its eigenvectors do not move, so every result at weight scale c is c times, or the same as, the result
/// at unit scale. An absolute stopping bound breaks this: below it no pivot rotates, and at 10⁻²⁰ the solve returns the degrees.
/// </summary>
public sealed class GraphSpectralScaleTests
{
    public static IEnumerable<object[]> Scales => new[] { 1e-40, 1e-20, 1e-8, 1.0, 1e8, 1e20, 1e40 }.Select(c => new object[] { c });

    [Theory]
    [MemberData(nameof(Scales))]
    public void AlgebraicConnectivity_ScalesWithTheWeights(double scale)
    {
        double unit = AlgebraicConnectivity.Compute(WeightedPath(1.0)).Lambda2;
        AlgebraicConnectivityReport report = AlgebraicConnectivity.Compute(WeightedPath(scale));

        Assert.True(report.Computed);
        Assert.InRange(Math.Abs(report.Lambda2 / scale - unit), 0.0, 1e-12 * unit);
    }

    [Theory]
    [MemberData(nameof(Scales))]
    public void ComputeBottomK_Dense_ScalesWithTheWeights(double scale)
    {
        foreach (DenseLaplacianMaterialization layout in new[] { DenseLaplacianMaterialization.Rectangular, DenseLaplacianMaterialization.FlatColumnMajor })
        {
            IReadOnlyList<EigenPair> unit = Spectral.ComputeBottomK(WeightedPath(1.0), k: 5, denseMaterialization: layout);
            IReadOnlyList<EigenPair> scaled = Spectral.ComputeBottomK(WeightedPath(scale), k: 5, denseMaterialization: layout);
            double top = unit[^1].Lambda;

            Assert.Equal(5, scaled.Count);
            for (int j = 0; j < 5; j++)
            {
                // Eigenvalues to 10⁻¹² of the largest, eigenvectors (distinct eigenvalues) to the same up to sign.
                Assert.InRange(Math.Abs(scaled[j].Lambda / scale - unit[j].Lambda), 0.0, 1e-12 * top);
                double overlap = Math.Abs(scaled[j].Vector.Zip(unit[j].Vector, (a, b) => a * b).Sum());
                Assert.InRange(1.0 - overlap, -1e-12, 1e-12);
            }
        }
    }

    /// <summary>
    /// A star of five leaves around node 0, every alignment a: G = I + a·A, whose top eigenvector is (2, 1, 1, 1, 1, 1)/√20 for
    /// eigenvalue 1 + a√5, the hub weighted by √(leaves).
    /// </summary>
    [Fact]
    public void AlignmentEigenCentrality_StarHasTheKnownPerronVector()
    {
        const int leaves = 5;
        CsrGraph star = CsrGraph.FromEdges([.. Enumerable.Range(1, leaves).Select(i => new Edge(0, i, 1.0))], leaves + 1);
        var alignments = new Alignments { Temperature = 1.0, G = [.. Enumerable.Repeat(0.3, star.Targets.Length)] };

        double[] centrality = new AlignmentEigenCentrality().Compute(alignments, star);

        double norm = Math.Sqrt(2.0 * leaves);
        Assert.InRange(Math.Abs(centrality[0] - Math.Sqrt(leaves) / norm), 0.0, 1e-12);
        for (int i = 1; i <= leaves; i++) Assert.InRange(Math.Abs(centrality[i] - 1.0 / norm), 0.0, 1e-12);
    }

    // A path of five nodes with weights 1, 2, 3, 4, all scaled: its Laplacian's five eigenvalues are distinct.
    private static CsrGraph WeightedPath(double scale) =>
        CsrGraph.FromEdges([.. Enumerable.Range(0, 4).Select(i => new Edge(i, i + 1, (i + 1) * scale))], 5);
}
