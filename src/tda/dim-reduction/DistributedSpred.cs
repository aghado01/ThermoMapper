using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using Maths.Geometry.DimReduction;
using Maths.Geometry;
using Maths.Geometry.Estimators.Intrinsic;

namespace TDA.DimReduction;

/// <summary>Diagnostics for one local SPRED block.</summary>
/// <param name="Index">Zero-based block index.</param>
/// <param name="Start">Zero-based index of the block's first row in the original input.</param>
/// <param name="Count">Number of input rows assigned to the block.</param>
/// <param name="Seed">Seed passed to the block annealer; <c>null</c> draws OS entropy.</param>
/// <param name="Projection">The block's locally optimized k x d orthonormal projection.</param>
/// <param name="LocalObjective">The block objective evaluated at <paramref name="Projection"/>.</param>
/// <param name="AggregateObjective">The block objective evaluated at the final aggregate projection.</param>
public sealed record DistributedSpredBlockResult(
    int Index,
    int Start,
    int Count,
    int? Seed,
    double[][] Projection,
    double LocalObjective,
    double AggregateObjective);

/// <summary>The aggregate projection and diagnostics from a distributed SPRED run.</summary>
/// <param name="AmbientDimension">Input dimension d.</param>
/// <param name="TargetDimension">Projection dimension k.</param>
/// <param name="Projection">The k x d Grassmann-median projection aggregated from all blocks.</param>
/// <param name="FullDataObjective">The full input objective evaluated at <paramref name="Projection"/>.</param>
/// <param name="Blocks">Per-block results in ascending input-partition order.</param>
public sealed record DistributedSpredResult(
    int AmbientDimension,
    int TargetDimension,
    double[][] Projection,
    double FullDataObjective,
    IReadOnlyList<DistributedSpredBlockResult> Blocks)
{
    /// <summary>Number of input blocks represented by <see cref="Blocks"/>.</summary>
    public int BlockCount => Blocks.Count;
}

/// <summary>
/// Distributed SPRED (§3.2): run SPRED independently on contiguous data blocks, then aggregate the
/// resulting projection subspaces by geometric median on the Grassmann manifold.
/// </summary>
public static class DistributedSpred
{
    /// <summary>
    /// Run SPRED on contiguous, non-overlapping input blocks and return their Grassmann-median projection.
    /// </summary>
    /// <param name="data">Row-major ambient samples.</param>
    /// <param name="targetDim">Projection dimension k.</param>
    /// <param name="blockCount">Number of contiguous blocks; must be between one and the row count.</param>
    /// <param name="objective">Persistent-homology objective shared by all local runs.</param>
    /// <param name="maxIters">Simulated-annealing steps per block.</param>
    /// <param name="seed">Base RNG seed; block i receives <c>seed + 1009 * i</c>. Null draws OS entropy.</param>
    /// <param name="maxDegreeOfParallelism">Maximum concurrent block runs. One preserves serial execution.</param>
    /// <param name="cancellationToken">Cancellation observed between annealing iterations and pipeline phases.</param>
    /// <returns>The aggregate k x d orthonormal projection.</returns>
    public static double[][] Compute(
        double[][] data,
        int targetDim,
        int blockCount,
        PersistenceObjectiveConfig objective,
        int maxIters = 1000,
        int? seed = null,
        int maxDegreeOfParallelism = 1,
        CancellationToken cancellationToken = default)
    {
        int ambientDim = ValidateInputs(data, targetDim, blockCount, objective, maxDegreeOfParallelism);
        cancellationToken.ThrowIfCancellationRequested();
        if (blockCount == 1)
            return Spred.Compute(data, targetDim, objective, maxIters, seed, cancellationToken);

        var projections = new double[blockCount][][];
        RunBlocks(blockCount, maxDegreeOfParallelism, cancellationToken, block =>
        {
            int start = block * data.Length / blockCount;
            int end = (block + 1) * data.Length / blockCount;
            double[][] slice = SliceRows(data, start, end);
            projections[block] = Spred.Compute(
                slice,
                targetDim,
                objective,
                maxIters,
                BlockSeed(seed, block),
                cancellationToken);
        });

        cancellationToken.ThrowIfCancellationRequested();
        return AggregateProjections(projections, ambientDim, targetDim);
    }

    /// <summary>
    /// Run distributed SPRED and also evaluate each local projection, the aggregate on every block,
    /// and the aggregate on the full input.
    /// </summary>
    /// <param name="data">Row-major ambient samples.</param>
    /// <param name="targetDim">Projection dimension k.</param>
    /// <param name="blockCount">Number of contiguous blocks; must be between one and the row count.</param>
    /// <param name="objective">Persistent-homology objective shared by all local runs.</param>
    /// <param name="maxIters">Simulated-annealing steps per block.</param>
    /// <param name="seed">Base RNG seed; block i receives <c>seed + 1009 * i</c>. Null draws OS entropy.</param>
    /// <param name="maxDegreeOfParallelism">Maximum concurrent block runs. One preserves serial execution.</param>
    /// <param name="cancellationToken">Cancellation observed between annealing iterations and pipeline phases.</param>
    /// <returns>The aggregate projection and ordered block diagnostics.</returns>
    /// <remarks>This path performs additional objective evaluations and is more expensive than <see cref="Compute"/>.</remarks>
    public static DistributedSpredResult ComputeWithDiagnostics(
        double[][] data,
        int targetDim,
        int blockCount,
        PersistenceObjectiveConfig objective,
        int maxIters = 1000,
        int? seed = null,
        int maxDegreeOfParallelism = 1,
        CancellationToken cancellationToken = default)
    {
        int ambientDim = ValidateInputs(data, targetDim, blockCount, objective, maxDegreeOfParallelism);
        cancellationToken.ThrowIfCancellationRequested();

        if (blockCount == 1)
        {
            BlockRun run = RunBlock(0, 0, data, targetDim, objective, maxIters, seed, cancellationToken);
            double singleBlockObjective = run.AggregateObjective(run.Projection, cancellationToken);
            var blocks = new[]
            {
                run.ToResult(run.Projection, cancellationToken),
            };
            return new DistributedSpredResult(ambientDim, targetDim, run.Projection, singleBlockObjective, blocks);
        }

        var blockRuns = new BlockRun[blockCount];
        RunBlocks(blockCount, maxDegreeOfParallelism, cancellationToken, block =>
        {
            int start = block * data.Length / blockCount;
            int end = (block + 1) * data.Length / blockCount;
            double[][] slice = SliceRows(data, start, end);
            int? blockSeed = BlockSeed(seed, block);
            blockRuns[block] = RunBlock(
                block,
                start,
                slice,
                targetDim,
                objective,
                maxIters,
                blockSeed,
                cancellationToken);
        });

        var projections = new double[blockCount][][];
        for (int block = 0; block < blockCount; block++)
            projections[block] = blockRuns[block].Projection;

        cancellationToken.ThrowIfCancellationRequested();
        double[][] aggregate = AggregateProjections(projections, ambientDim, targetDim);
        var blockResults = new DistributedSpredBlockResult[blockCount];
        RunBlocks(blockCount, maxDegreeOfParallelism, cancellationToken, block =>
            blockResults[block] = blockRuns[block].ToResult(aggregate, cancellationToken));

        cancellationToken.ThrowIfCancellationRequested();
        var fullObjective = new PersistenceObjective(data, objective);
        cancellationToken.ThrowIfCancellationRequested();
        double fullDataObjective = fullObjective.Evaluate(aggregate);
        cancellationToken.ThrowIfCancellationRequested();
        return new DistributedSpredResult(ambientDim, targetDim, aggregate, fullDataObjective, blockResults);
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

        public double AggregateObjective(double[][] aggregate, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            double value = _objective.Evaluate(aggregate);
            cancellationToken.ThrowIfCancellationRequested();
            return value;
        }

        public DistributedSpredBlockResult ToResult(double[][] aggregate, CancellationToken cancellationToken)
        {
            return new DistributedSpredBlockResult(
                Index,
                Start,
                Count,
                Seed,
                Projection,
                LocalObjective,
                AggregateObjective(aggregate, cancellationToken));
        }
    }

    private static int ValidateInputs(
        double[][] data,
        int targetDim,
        int blockCount,
        PersistenceObjectiveConfig objective,
        int maxDegreeOfParallelism)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(objective);
        if (data.Length == 0) throw new ArgumentException("Empty data", nameof(data));
        if (blockCount < 1 || blockCount > data.Length)
            throw new ArgumentOutOfRangeException(nameof(blockCount), "Block count must be between 1 and data.Length.");
        if (maxDegreeOfParallelism < 1)
            throw new ArgumentOutOfRangeException(
                nameof(maxDegreeOfParallelism),
                "Maximum degree of parallelism must be at least one.");

        int ambientDim = data[0].Length;
        if (targetDim < 1 || targetDim > ambientDim)
            throw new ArgumentOutOfRangeException(nameof(targetDim), "Target dimension must satisfy 1 <= targetDim <= ambient dimension.");
        for (int i = 1; i < data.Length; i++)
            if (data[i].Length != ambientDim)
                throw new ArgumentException("All data rows must have the same dimension.", nameof(data));

        return ambientDim;
    }

    private static void RunBlocks(
        int blockCount,
        int maxDegreeOfParallelism,
        CancellationToken cancellationToken,
        Action<int> body)
    {
        cancellationToken.ThrowIfCancellationRequested();
        if (maxDegreeOfParallelism == 1 || blockCount == 1)
        {
            for (int block = 0; block < blockCount; block++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                body(block);
            }
            return;
        }

        Parallel.For(
            0,
            blockCount,
            new ParallelOptions
            {
                MaxDegreeOfParallelism = maxDegreeOfParallelism,
                CancellationToken = cancellationToken,
            },
            body);
    }

    private static BlockRun RunBlock(
        int index,
        int start,
        double[][] data,
        int targetDim,
        PersistenceObjectiveConfig objective,
        int maxIters,
        int? seed,
        CancellationToken cancellationToken)
    {
        cancellationToken.ThrowIfCancellationRequested();
        var ph = new PersistenceObjective(data, objective);
        cancellationToken.ThrowIfCancellationRequested();
        double[][] projection = SubspaceAnnealer.Compute(
            data,
            targetDim,
            ph.Evaluate,
            maxIters,
            seed,
            cancellationToken);
        cancellationToken.ThrowIfCancellationRequested();
        double localObjective = ph.Evaluate(projection);
        cancellationToken.ThrowIfCancellationRequested();
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
