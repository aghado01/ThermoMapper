using System;
using System.Collections.Generic;
using Maths.Geometry.DimReduction;
using Maths.Geometry;
using Maths.Geometry.Estimators.Intrinsic;

namespace TDA.DimReduction;

public sealed record DistributedSpredBlockResult(
    int Index,
    int Start,
    int Count,
    int? Seed,
    double[][] Projection,
    double LocalObjective,
    double AggregateObjective);

public sealed record DistributedSpredResult(
    int AmbientDimension,
    int TargetDimension,
    double[][] Projection,
    IReadOnlyList<DistributedSpredBlockResult> Blocks)
{
    public int BlockCount => Blocks.Count;
}

/// <summary>
/// Distributed SPRED (§3.2): run SPRED independently on contiguous data blocks, then aggregate the
/// resulting projection subspaces by geometric median on the Grassmann manifold.
/// </summary>
public static class DistributedSpred
{
    public static double[][] Compute(
        double[][] data,
        int targetDim,
        int blockCount,
        PersistenceObjectiveConfig objective,
        int maxIters = 1000,
        int? seed = null)
    {
        return ComputeWithDiagnostics(data, targetDim, blockCount, objective, maxIters, seed).Projection;
    }

    public static DistributedSpredResult ComputeWithDiagnostics(
        double[][] data,
        int targetDim,
        int blockCount,
        PersistenceObjectiveConfig objective,
        int maxIters = 1000,
        int? seed = null)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(objective);
        if (data.Length == 0) throw new ArgumentException("Empty data", nameof(data));
        if (blockCount < 1 || blockCount > data.Length)
            throw new ArgumentOutOfRangeException(nameof(blockCount), "Block count must be between 1 and data.Length.");

        int ambientDim = data[0].Length;
        if (targetDim < 1 || targetDim > ambientDim)
            throw new ArgumentOutOfRangeException(nameof(targetDim), "Target dimension must satisfy 1 <= targetDim <= ambient dimension.");
        for (int i = 1; i < data.Length; i++)
            if (data[i].Length != ambientDim)
                throw new ArgumentException("All data rows must have the same dimension.", nameof(data));

        if (blockCount == 1)
        {
            BlockRun run = RunBlock(0, 0, data, targetDim, objective, maxIters, seed);
            var blocks = new[]
            {
                run.ToResult(run.Projection),
            };
            return new DistributedSpredResult(ambientDim, targetDim, run.Projection, blocks);
        }

        var projections = new List<double[][]>(blockCount);
        var blockRuns = new List<BlockRun>(blockCount);
        for (int block = 0; block < blockCount; block++)
        {
            int start = block * data.Length / blockCount;
            int end = (block + 1) * data.Length / blockCount;
            double[][] slice = SliceRows(data, start, end);
            int? blockSeed = BlockSeed(seed, block);
            BlockRun run = RunBlock(block, start, slice, targetDim, objective, maxIters, blockSeed);
            projections.Add(run.Projection);
            blockRuns.Add(run);
        }

        double[][] aggregate = AggregateProjections(projections, ambientDim, targetDim);
        var blockResults = new List<DistributedSpredBlockResult>(blockCount);
        foreach (BlockRun run in blockRuns)
            blockResults.Add(run.ToResult(aggregate));

        return new DistributedSpredResult(ambientDim, targetDim, aggregate, blockResults);
    }

    private sealed class BlockRun
    {
        private readonly PersistenceObjective _objective;

        public BlockRun(
            int index,
            int start,
            int count,
            int? seed,
            double[][] projection,
            double localObjective,
            PersistenceObjective objective)
        {
            Index = index;
            Start = start;
            Count = count;
            Seed = seed;
            Projection = projection;
            LocalObjective = localObjective;
            _objective = objective;
        }

        public int Index { get; }
        public int Start { get; }
        public int Count { get; }
        public int? Seed { get; }
        public double[][] Projection { get; }
        public double LocalObjective { get; }

        public DistributedSpredBlockResult ToResult(double[][] aggregate)
        {
            return new DistributedSpredBlockResult(
                Index,
                Start,
                Count,
                Seed,
                Projection,
                LocalObjective,
                _objective.Evaluate(aggregate));
        }
    }

    private static BlockRun RunBlock(
        int index,
        int start,
        double[][] data,
        int targetDim,
        PersistenceObjectiveConfig objective,
        int maxIters,
        int? seed)
    {
        var ph = new PersistenceObjective(data, objective);
        double[][] projection = SubspaceAnnealer.Compute(data, targetDim, ph.Evaluate, maxIters, seed);
        double localObjective = ph.Evaluate(projection);
        return new BlockRun(index, start, data.Length, seed, projection, localObjective, ph);
    }

    internal static double[][] AggregateProjections(
        IReadOnlyList<double[][]> projections,
        int ambientDim,
        int targetDim)
    {
        if (projections.Count == 0)
            throw new ArgumentException("At least one projection is required.", nameof(projections));

        var frames = new double[projections.Count][];
        for (int i = 0; i < projections.Count; i++)
            frames[i] = ProjectionToFrame(projections[i], ambientDim, targetDim);

        double[] median = (double[])frames[0].Clone();
        double[] weights = new double[frames.Length];
        Array.Fill(weights, 1.0);

        var grass = new GrassmannManifold(ambientDim, targetDim);
        GeometricMedian.Compute(grass, frames, weights, median);

        return FrameToProjection(median, ambientDim, targetDim);
    }

    private static double[][] SliceRows(double[][] data, int start, int end)
    {
        var slice = new double[end - start][];
        Array.Copy(data, start, slice, 0, slice.Length);
        return slice;
    }

    private static int? BlockSeed(int? seed, int block)
    {
        return seed is null ? null : unchecked(seed.Value + 1009 * block);
    }

    private static double[] ProjectionToFrame(double[][] projection, int ambientDim, int targetDim)
    {
        if (projection.Length != targetDim)
            throw new ArgumentException("Projection row count must match target dimension.", nameof(projection));

        var frame = new double[ambientDim * targetDim];
        for (int col = 0; col < targetDim; col++)
        {
            if (projection[col].Length != ambientDim)
                throw new ArgumentException("Projection row length must match ambient dimension.", nameof(projection));
            Array.Copy(projection[col], 0, frame, col * ambientDim, ambientDim);
        }
        return frame;
    }

    private static double[][] FrameToProjection(double[] frame, int ambientDim, int targetDim)
    {
        var projection = new double[targetDim][];
        for (int row = 0; row < targetDim; row++)
        {
            projection[row] = new double[ambientDim];
            Array.Copy(frame, row * ambientDim, projection[row], 0, ambientDim);
        }
        return projection;
    }
}
