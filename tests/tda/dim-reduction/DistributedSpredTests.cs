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
    public void AggregateProjections_CleanMajorityResistsCorruptedBlocks()
    {
        double[][] xy =
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        double[][] xyTiltX =
        [
            [Math.Cos(0.04), 0.0, Math.Sin(0.04)],
            [0.0, 1.0, 0.0],
        ];
        double[][] xyTiltY =
        [
            [1.0, 0.0, 0.0],
            [0.0, Math.Cos(0.05), Math.Sin(0.05)],
        ];
        double[][] xz =
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        double[][] yz =
        [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];

        double[][] aggregated = DistributedSpred.AggregateProjections(
            new[] { xy, xyTiltX, xyTiltY, xz, yz },
            ambientDim: 3,
            targetDim: 2);

        var grass = new GrassmannManifold(ambientN: 3, subspaceR: 2);
        double cleanDistance = grass.Distance(PackFrame(aggregated), PackFrame(xy));
        double corruptedXzDistance = grass.Distance(PackFrame(aggregated), PackFrame(xz));
        double corruptedYzDistance = grass.Distance(PackFrame(aggregated), PackFrame(yz));

        Assert.InRange(cleanDistance, 0.0, 0.1);
        Assert.True(cleanDistance < corruptedXzDistance);
        Assert.True(cleanDistance < corruptedYzDistance);
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

    [Fact]
    public void Compute_MultipleBlocks_RunsSplitAndAggregatesProjection()
    {
        double[][] projection = DistributedSpred.Compute(
            RepeatedCircleBlocks(blockCount: 2, pointsPerBlock: 24),
            targetDim: 2,
            blockCount: 2,
            SmallConfig(),
            maxIters: 0,
            seed: 17);

        AssertOrthonormalRows(projection);

        double[][] xy =
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        var grass = new GrassmannManifold(ambientN: 3, subspaceR: 2);
        double distance = grass.Distance(PackFrame(projection), PackFrame(xy));
        Assert.InRange(distance, 0.0, 1e-8);
    }

    [Fact]
    public void ComputeWithDiagnostics_MultipleBlocks_ReportsBlockMetadataAndProjections()
    {
        const int seed = 23;
        const int blockCount = 3;
        const int pointsPerBlock = 18;

        DistributedSpredResult result = DistributedSpred.ComputeWithDiagnostics(
            RepeatedCircleBlocks(blockCount, pointsPerBlock),
            targetDim: 2,
            blockCount,
            SmallConfig(),
            maxIters: 0,
            seed);

        Assert.Equal(3, result.AmbientDimension);
        Assert.Equal(2, result.TargetDimension);
        Assert.Equal(blockCount, result.BlockCount);
        Assert.Equal(blockCount, result.Blocks.Count);
        AssertOrthonormalRows(result.Projection);
        AssertFiniteFullObjective(result);

        double[][] xy =
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        var grass = new GrassmannManifold(ambientN: 3, subspaceR: 2);
        Assert.InRange(grass.Distance(PackFrame(result.Projection), PackFrame(xy)), 0.0, 1e-8);

        for (int block = 0; block < blockCount; block++)
        {
            DistributedSpredBlockResult info = result.Blocks[block];
            Assert.Equal(block, info.Index);
            Assert.Equal(block * pointsPerBlock, info.Start);
            Assert.Equal(pointsPerBlock, info.Count);
            Assert.Equal(seed + 1009 * block, info.Seed);
            AssertOrthonormalRows(info.Projection);
            Assert.InRange(grass.Distance(PackFrame(info.Projection), PackFrame(xy)), 0.0, 1e-8);
            AssertFiniteObjectives(info);
            Assert.Equal(info.LocalObjective, info.AggregateObjective, precision: 10);
        }
    }

    [Fact]
    public void ComputeWithDiagnostics_UnevenBlocks_CoversEveryRowWithoutTruncation()
    {
        const int rowCount = 50;
        const int blockCount = 3;

        DistributedSpredResult result = DistributedSpred.ComputeWithDiagnostics(
            Circle3D(rowCount),
            targetDim: 2,
            blockCount,
            SmallConfig(),
            maxIters: 0,
            seed: 29);

        Assert.Equal(blockCount, result.BlockCount);
        Assert.Equal(0, result.Blocks[0].Start);

        int covered = 0;
        for (int i = 0; i < result.Blocks.Count; i++)
        {
            DistributedSpredBlockResult block = result.Blocks[i];
            int expectedStart = i * rowCount / blockCount;
            int expectedEnd = (i + 1) * rowCount / blockCount;

            Assert.Equal(i, block.Index);
            Assert.Equal(expectedStart, block.Start);
            Assert.Equal(expectedEnd - expectedStart, block.Count);
            Assert.Equal(covered, block.Start);
            AssertOrthonormalRows(block.Projection);
            AssertFiniteObjectives(block);

            covered += block.Count;
        }

        Assert.Equal(rowCount, covered);
        Assert.Equal(rowCount, result.Blocks[^1].Start + result.Blocks[^1].Count);
        AssertOrthonormalRows(result.Projection);
        AssertFiniteFullObjective(result);
    }

    [Fact]
    public void ComputeWithDiagnostics_CorruptedBlocks_AggregatesCleanMajority()
    {
        const int pointsPerBlock = 24;

        DistributedSpredResult result = DistributedSpred.ComputeWithDiagnostics(
            CorruptedCircleBlocks(pointsPerBlock),
            targetDim: 2,
            blockCount: 5,
            SmallConfig(),
            maxIters: 0,
            seed: 31);

        double[][] xy =
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        double[][] xz =
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        double[][] yz =
        [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];

        var grass = new GrassmannManifold(ambientN: 3, subspaceR: 2);
        double aggregateToClean = grass.Distance(PackFrame(result.Projection), PackFrame(xy));
        double aggregateToXz = grass.Distance(PackFrame(result.Projection), PackFrame(xz));
        double aggregateToYz = grass.Distance(PackFrame(result.Projection), PackFrame(yz));

        Assert.Equal(5, result.BlockCount);
        AssertFiniteFullObjective(result);
        Assert.InRange(aggregateToClean, 0.0, 1e-8);
        Assert.True(aggregateToClean < aggregateToXz);
        Assert.True(aggregateToClean < aggregateToYz);

        for (int block = 0; block < 3; block++)
            Assert.InRange(grass.Distance(PackFrame(result.Blocks[block].Projection), PackFrame(xy)), 0.0, 1e-8);

        Assert.InRange(grass.Distance(PackFrame(result.Blocks[3].Projection), PackFrame(xz)), 0.0, 1e-8);
        Assert.InRange(grass.Distance(PackFrame(result.Blocks[4].Projection), PackFrame(yz)), 0.0, 1e-8);
        Assert.True(grass.Distance(PackFrame(result.Blocks[3].Projection), PackFrame(xy)) > 0.5);
        Assert.True(grass.Distance(PackFrame(result.Blocks[4].Projection), PackFrame(xy)) > 0.5);

        for (int block = 0; block < 3; block++)
        {
            AssertFiniteObjectives(result.Blocks[block]);
            Assert.Equal(result.Blocks[block].LocalObjective, result.Blocks[block].AggregateObjective, precision: 10);
        }

        AssertFiniteObjectives(result.Blocks[3]);
        AssertFiniteObjectives(result.Blocks[4]);
        Assert.True(result.Blocks[3].AggregateObjective > result.Blocks[3].LocalObjective);
        Assert.True(result.Blocks[4].AggregateObjective > result.Blocks[4].LocalObjective);
    }

    [Fact]
    public void ComputeWithDiagnostics_CleanFixture_MatchesGlobalSpredBaseline()
    {
        double[][] data = RepeatedCircleBlocks(blockCount: 3, pointsPerBlock: 24);
        PersistenceObjectiveConfig config = SmallConfig();

        double[][] global = Spred.Compute(data, targetDim: 2, config, maxIters: 0, seed: 41);
        DistributedSpredResult distributed = DistributedSpred.ComputeWithDiagnostics(
            data,
            targetDim: 2,
            blockCount: 3,
            config,
            maxIters: 0,
            seed: 41);

        double globalObjective = EvaluateFullDataObjective(data, config, global);

        AssertOrthonormalRows(global);
        AssertOrthonormalRows(distributed.Projection);
        AssertFiniteFullObjective(distributed);
        Assert.Equal(globalObjective, distributed.FullDataObjective, precision: 10);

        var grass = new GrassmannManifold(ambientN: 3, subspaceR: 2);
        Assert.InRange(grass.Distance(PackFrame(global), PackFrame(distributed.Projection)), 0.0, 1e-8);
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

    private static double[][] RepeatedCircleBlocks(int blockCount, int pointsPerBlock)
    {
        var data = new double[blockCount * pointsPerBlock][];
        for (int block = 0; block < blockCount; block++)
        {
            double[][] circle = Circle3D(pointsPerBlock);
            Array.Copy(circle, 0, data, block * pointsPerBlock, pointsPerBlock);
        }
        return data;
    }

    private static double[][] CorruptedCircleBlocks(int pointsPerBlock)
    {
        var data = new double[5 * pointsPerBlock][];
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 1), data, 0, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 1), data, 1, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 1), data, 2, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 2), data, 3, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 1, axisB: 2), data, 4, pointsPerBlock);
        return data;
    }

    private static void CopyBlock(double[][] block, double[][] data, int blockIndex, int pointsPerBlock)
    {
        Array.Copy(block, 0, data, blockIndex * pointsPerBlock, pointsPerBlock);
    }

    private static double[][] CircleInPlane(int n, int axisA, int axisB)
    {
        var pts = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double t = 2.0 * Math.PI * i / n;
            pts[i] = new double[3];
            pts[i][axisA] = Math.Cos(t);
            pts[i][axisB] = Math.Sin(t);
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

    private static void AssertFiniteObjectives(DistributedSpredBlockResult block)
    {
        Assert.True(double.IsFinite(block.LocalObjective));
        Assert.True(double.IsFinite(block.AggregateObjective));
        Assert.True(block.LocalObjective >= 0.0);
        Assert.True(block.AggregateObjective >= 0.0);
    }

    private static void AssertFiniteFullObjective(DistributedSpredResult result)
    {
        Assert.True(double.IsFinite(result.FullDataObjective));
        Assert.True(result.FullDataObjective >= 0.0);
    }

    private static double EvaluateFullDataObjective(
        double[][] data,
        PersistenceObjectiveConfig config,
        double[][] projection)
    {
        return new PersistenceObjective(data, config).Evaluate(projection);
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
