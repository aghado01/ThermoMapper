#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Graphs.Primitives;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

/// <summary>
/// Random graph-restricted Rips filtrations — parity of pCoh vs pHcol.
/// </summary>
public sealed class PersistentCohomologyRipsParityTests
{
    static IEnumerable<(double Birth, double Death, int Dim)> Signatures(Barcode bc) =>
        bc.Bars
            .Select(b => (b.Birth, b.Death, b.Dimension))
            .OrderBy(x => x.Dimension)
            .ThenBy(x => x.Birth)
            .ThenBy(x => x.Death);

    static void AssertBarcodeParity(SimplicialFiltration filtration, int maxDimension = int.MaxValue)
    {
        Barcode ph = PersistentHomology.Compute(filtration, maxDimension);
        Barcode pcoh = PersistentCohomology.Compute(filtration, maxDimension);
        Assert.Equal(Signatures(ph), Signatures(pcoh));
    }

    [Theory]
    [InlineData(42, 0)]
    [InlineData(42, 1)]
    [InlineData(42, 2)]
    [InlineData(99, 0)]
    [InlineData(99, 1)]
    [InlineData(99, 2)]
    public void RandomRipsGraph_ParityWithHomology(int seed, int trial)
    {
        var rng = new Random(seed * 31 + trial);
        int n = rng.Next(4, 10);
        var edges = new List<Edge>();
        double nextWeight = 1.0;
        for (int u = 0; u < n; u++)
        {
            for (int v = u + 1; v < n; v++)
            {
                if (rng.NextDouble() < 0.45)
                {
                    edges.Add(new Edge(u, v, nextWeight));
                    nextWeight += 1.0;
                }
            }
        }

        if (edges.Count == 0)
            edges.Add(new Edge(0, 1, 1.0));

        var g = CsrGraph.FromEdges(edges.ToArray(), nodeCount: n);
        var filtration = RipsFiltration.GraphRips(g, FiltrationWeights.RawDistance, maxDimension: 2);
        AssertBarcodeParity(filtration);
    }
}
