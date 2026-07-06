#nullable enable
using System.Collections.Generic;
using System.Linq;
using Graphs.Primitives;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

public sealed class LazyRipsFiltrationTests
{
    static IEnumerable<(double Birth, double Death, int Dimension)> Signatures(Barcode bc) =>
        bc.Bars
            .Select(b => (b.Birth, b.Death, b.Dimension))
            .OrderBy(x => x.Dimension)
            .ThenBy(x => x.Birth)
            .ThenBy(x => x.Death);

    [Fact]
    public void MatchesMaterializedRipsFiltration()
    {
        Edge[] edges = new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(0, 2, 1.0),
            new Edge(1, 2, 1.0),
            new Edge(0, 3, 2.0),
            new Edge(1, 3, 2.0),
            new Edge(2, 3, 2.0),
        };
        var g = CsrGraph.FromEdges(edges, 4);

        var materialized = RipsFiltration.GraphRips(g, FiltrationWeights.RawDistance, maxDimension: 2);
        var lazy = new LazyRipsFiltration(g, FiltrationWeights.RawDistance);

        Barcode phMat = PersistentHomology.Compute(materialized);
        Barcode phLazy = PersistentHomology.Compute(lazy);
        Barcode pcohMat = PersistentCohomology.Compute(materialized);
        Barcode pcohLazy = PersistentCohomology.Compute(lazy);

        Assert.Equal(Signatures(phMat), Signatures(phLazy));
        Assert.Equal(Signatures(pcohMat), Signatures(pcohLazy));
    }
}
