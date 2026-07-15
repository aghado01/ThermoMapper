using System;
using System.Diagnostics;
using Graphs;
using Maths.Geometry;
using Xunit;
using Xunit.Abstractions;

namespace TDA.DimReduction.Tests;

/// <summary>Opt-in scale profile for global, serial-block, and parallel-block SPRED execution.</summary>
public sealed class DistributedSpredProfileTests
{
    private const string ProfileEnvironmentVariable = "THERMOMAPPER_SPRED_PROFILE";
    private readonly ITestOutputHelper _output;

    public DistributedSpredProfileTests(ITestOutputHelper output) => _output = output;

    [Fact]
    [Trait("Category", "Benchmark")]
    public void Profile_GlobalVsSerialAndParallelBlocks()
    {
        if (!string.Equals(
                Environment.GetEnvironmentVariable(ProfileEnvironmentVariable),
                "1",
                StringComparison.Ordinal))
        {
            _output.WriteLine($"Set {ProfileEnvironmentVariable}=1 to run the distributed SPRED profile.");
            return;
        }

        const int blockCount = 6;
        const int pointsPerBlock = 32;
        const int maxIters = 4;
        const int seed = 71;
        int parallelWorkers = Math.Min(blockCount, Environment.ProcessorCount);
        double[][] data = NoisyCircleBlocks(blockCount, pointsPerBlock);
        PersistenceObjectiveConfig config = ProfileConfig();

        // Warm the block path before collecting timings.
        DistributedSpred.Compute(
            NoisyCircleBlocks(blockCount: 2, pointsPerBlock: 12),
            targetDim: 2,
            blockCount: 2,
            config,
            maxIters: 0,
            seed,
            maxDegreeOfParallelism: 1);

        TimedProjection global = Measure(() =>
            Spred.Compute(data, targetDim: 2, config, maxIters, seed));
        TimedProjection serial = Measure(() =>
            DistributedSpred.Compute(
                data,
                targetDim: 2,
                blockCount,
                config,
                maxIters,
                seed,
                maxDegreeOfParallelism: 1));
        TimedProjection parallel = Measure(() =>
            DistributedSpred.Compute(
                data,
                targetDim: 2,
                blockCount,
                config,
                maxIters,
                seed,
                maxDegreeOfParallelism: parallelWorkers));

        var fullObjective = new PersistenceObjective(data, config);
        TimedObjective globalObjective = Evaluate(fullObjective, global.Projection);
        TimedObjective serialObjective = Evaluate(fullObjective, serial.Projection);
        TimedObjective parallelObjective = Evaluate(fullObjective, parallel.Projection);
        double serialParallelDistance = DistanceBetweenProjections(serial.Projection, parallel.Projection);

        _output.WriteLine("mode                 workers  fit(ms)  full-eval(ms)  full-objective");
        WriteResult("global", 1, global, globalObjective);
        WriteResult("distributed-serial", 1, serial, serialObjective);
        WriteResult("distributed-parallel", parallelWorkers, parallel, parallelObjective);
        _output.WriteLine($"serial/parallel Grassmann distance: {serialParallelDistance:G6}");

        Assert.True(double.IsFinite(globalObjective.Value));
        Assert.True(double.IsFinite(serialObjective.Value));
        Assert.True(double.IsFinite(parallelObjective.Value));
        Assert.InRange(serialParallelDistance, 0.0, 1e-10);
        Assert.Equal(serialObjective.Value, parallelObjective.Value, precision: 10);
    }

    private void WriteResult(string mode, int workers, TimedProjection fit, TimedObjective objective)
    {
        _output.WriteLine(
            $"{mode,-20} {workers,7} {fit.ElapsedMilliseconds,8:F1} " +
            $"{objective.ElapsedMilliseconds,14:F1} {objective.Value,15:G8}");
    }

    private static TimedProjection Measure(Func<double[][]> fit)
    {
        var stopwatch = Stopwatch.StartNew();
        double[][] projection = fit();
        stopwatch.Stop();
        return new TimedProjection(projection, stopwatch.Elapsed.TotalMilliseconds);
    }

    private static TimedObjective Evaluate(PersistenceObjective objective, double[][] projection)
    {
        var stopwatch = Stopwatch.StartNew();
        double value = objective.Evaluate(projection);
        stopwatch.Stop();
        return new TimedObjective(value, stopwatch.Elapsed.TotalMilliseconds);
    }

    private static PersistenceObjectiveConfig ProfileConfig() => new()
    {
        Graph = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 8 },
            Filter = new FilterConfig { Kind = FilterKind.OrRule },
            Repair = new RepairConfig { Kind = RepairKind.NoRepair },
            Projection = new DistanceProjection(),
        },
        Dimensions = [(1, 1.0)],
        MaxDimension = 2,
        MinPersistence = 0.01,
    };

    private static double[][] NoisyCircleBlocks(int blockCount, int pointsPerBlock)
    {
        var data = new double[blockCount * pointsPerBlock][];
        for (int block = 0; block < blockCount; block++)
        {
            for (int point = 0; point < pointsPerBlock; point++)
            {
                double angle = 2.0 * Math.PI * point / pointsPerBlock;
                int row = block * pointsPerBlock + point;
                data[row] =
                [
                    Math.Cos(angle),
                    Math.Sin(angle),
                    0.05 * Math.Sin(3.0 * angle + 0.2 * block),
                ];
            }
        }
        return data;
    }

    private static double DistanceBetweenProjections(double[][] a, double[][] b)
    {
        var grass = new GrassmannManifold(ambientN: a[0].Length, subspaceR: a.Length);
        return grass.Distance(PackFrame(a), PackFrame(b));
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

    private sealed record TimedProjection(double[][] Projection, double ElapsedMilliseconds);
    private sealed record TimedObjective(double Value, double ElapsedMilliseconds);
}
