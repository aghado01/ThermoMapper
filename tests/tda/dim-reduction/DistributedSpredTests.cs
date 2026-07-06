using System;
using Graphs;
using Maths.Geometry;
using Xunit;

namespace TDA.DimReduction.Tests;

public sealed class DistributedSpredTests
{
    [Fact]
    public void AggregateProjections_DuplicateSubspaceWinsMedian()
    {
        double[][] xy =
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        double[][] xyRotated =
        [
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ];
        double[][] xz =
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ];

        double[][] aggregated = DistributedSpred.AggregateProjections(
            new[] { xy, xyRotated, xz },
            ambientDim: 3,
            targetDim: 2);

        var grass = new GrassmannManifold(ambientN: 3, subspaceR: 2);
        double distance = grass.Distance(PackFrame(aggregated), PackFrame(xy));
        Assert.InRange(distance, 0.0, 1e-8);
    }

    [Fact]
    public void Compute_SingleBlock_ReturnsOrthonormalProjection()
    {
        double[][] projection = DistributedSpred.Compute(
            Circle3D(32),
            targetDim: 2,
            blockCount: 1,
            SmallConfig(),
            maxIters: 8,
            seed: 11);

        AssertOrthonormalRows(projection);
    }

    private static PersistenceObjectiveConfig SmallConfig() => new()
    {
        Graph = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 6 },
            Filter = new FilterConfig { Kind = FilterKind.OrRule },
            Repair = new RepairConfig { Kind = RepairKind.NoRepair },
            Projection = new DistanceProjection(),
        },
        Dimensions = [(1, 1.0)],
        MaxDimension = 2,
    };

    private static double[][] Circle3D(int n)
    {
        var pts = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double t = 2.0 * Math.PI * i / n;
            pts[i] = [Math.Cos(t), Math.Sin(t), 0.0];
        }
        return pts;
    }

    private static void AssertOrthonormalRows(double[][] projection)
    {
        Assert.Equal(2, projection.Length);
        for (int i = 0; i < projection.Length; i++)
        {
            double self = Dot(projection[i], projection[i]);
            Assert.InRange(self, 1.0 - 1e-9, 1.0 + 1e-9);

            for (int j = i + 1; j < projection.Length; j++)
                Assert.InRange(Dot(projection[i], projection[j]), -1e-9, 1e-9);
        }
    }

    private static double Dot(double[] a, double[] b)
    {
        double sum = 0.0;
        for (int i = 0; i < a.Length; i++) sum += a[i] * b[i];
        return sum;
    }

    private static double[] PackFrame(double[][] projection)
    {
        int targetDim = projection.Length;
        int ambientDim = projection[0].Length;
        var frame = new double[ambientDim * targetDim];
        for (int col = 0; col < targetDim; col++)
            Array.Copy(projection[col], 0, frame, col * ambientDim, ambientDim);
        return frame;
    }
}
