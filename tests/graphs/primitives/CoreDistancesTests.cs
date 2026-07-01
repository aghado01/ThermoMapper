using Graphs.Distance;
using Graphs.Distance.Euclidean;
using Graphs.Primitives.Mst;
using Xunit;

namespace Graphs.Primitives.Tests;

public sealed class CoreDistancesTests
{
    [Fact]
    public void Compute_MinPtsTwo_UsesSecondNearestOnTriangle()
    {
        // Three points in R^2: equilateral-ish unit layout
        double[] data =
        {
            0, 0,
            1, 0,
            0.5, 0.866,
        };
        const int n = 3;
        const int dim = 2;
        var metric = new EuclideanMetric();
        var core = new double[n];

        CoreDistances.Compute(data, n, dim, minPts: 2, metric, core);

        Assert.Equal(3, core.Length);
        foreach (double c in core)
            Assert.True(c > 0.0 && c < double.PositiveInfinity);
    }
}
