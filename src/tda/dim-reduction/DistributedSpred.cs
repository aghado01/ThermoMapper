using System;
using System.Collections.Generic;
using System.Runtime.ExceptionServices;
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
    /// <param name="annealerOptions">Proposal mixture, step adaptation, and cooling shared by every
    /// block's <see cref="SubspaceAnnealer"/>; null takes the engine defaults. Validated before any
    /// block work, so a bad record fails fast rather than inside a block wrap.</param>
    /// <param name="cancellationToken">Cancellation observed between annealing iterations and pipeline phases.</param>
    /// <returns>The aggregate k x d orthonormal projection.</returns>
    /// <remarks>
    /// With more than one block, a block whose objective construction or annealing fails surfaces as an
    /// <see cref="InvalidOperationException"/> naming the block index and row count, with the original
    /// failure as <see cref="Exception.InnerException"/> — validation cannot see the objective's graph
    /// recipe (e.g. its kNN K), so blocks too small for the recipe fail here rather than up front.
    /// Serial and parallel runs throw identically; cancellation surfaces as
    /// <see cref="OperationCanceledException"/>, never <see cref="AggregateException"/>.
    /// </remarks>
    public static double[][] Compute(
        double[][] data,
        int targetDim,
        int blockCount,
        PersistenceObjectiveConfig objective,
        int maxIters = 1000,
        int? seed = null,
        int maxDegreeOfParallelism = 1,
        SubspaceAnnealerOptions? annealerOptions = null,
        CancellationToken cancellationToken = default)
    {
        int ambientDim = ValidateInputs(data, targetDim, blockCount, objective, maxDegreeOfParallelism);
        annealerOptions?.Validate();
        cancellationToken.ThrowIfCancellationRequested();
        if (blockCount == 1)
            return Spred.Compute(data, targetDim, objective, maxIters, seed, annealerOptions, cancellationToken);

        var projections = new double[blockCount][][];
        RunBlocks(blockCount, maxDegreeOfParallelism, cancellationToken, block =>
        {
            int start = block * data.Length / blockCount;
            int end = (block + 1) * data.Length / blockCount;
            double[][] slice = SliceRows(data, start, end);
            try
            {
                projections[block] = Spred.Compute(
                    slice,
                    targetDim,
                    objective,
                    maxIters,
                    BlockSeed(seed, block),
                    annealerOptions,
                    cancellationToken);
            }
            catch (Exception failure) when (failure is not OperationCanceledException)
            {
                throw BlockSetupFailure(block, slice.Length, failure);
            }
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
    /// <param name="annealerOptions">Proposal mixture, step adaptation, and cooling shared by every
    /// block's <see cref="SubspaceAnnealer"/>; null takes the engine defaults. Validated before any
    /// block work, so a bad record fails fast rather than inside a block wrap.</param>
    /// <param name="cancellationToken">Cancellation observed between annealing iterations and pipeline phases.</param>
    /// <returns>The aggregate projection and ordered block diagnostics.</returns>
    /// <remarks>
    /// This path performs additional objective evaluations and is more expensive than <see cref="Compute"/>.
    /// The failure surface matches <see cref="Compute"/>: with more than one block, a block whose objective
    /// construction or annealing fails surfaces as an <see cref="InvalidOperationException"/> naming the
    /// block index and row count, with the original failure as <see cref="Exception.InnerException"/>,
    /// identically for serial and parallel runs; cancellation surfaces as
    /// <see cref="OperationCanceledException"/>, never <see cref="AggregateException"/>.
    /// </remarks>
    public static DistributedSpredResult ComputeWithDiagnostics(
        double[][] data,
        int targetDim,
        int blockCount,
        PersistenceObjectiveConfig objective,
        int maxIters = 1000,
        int? seed = null,
        int maxDegreeOfParallelism = 1,
        SubspaceAnnealerOptions? annealerOptions = null,
        CancellationToken cancellationToken = default)
    {
        int ambientDim = ValidateInputs(data, targetDim, blockCount, objective, maxDegreeOfParallelism);
        annealerOptions?.Validate();
        cancellationToken.ThrowIfCancellationRequested();

        if (blockCount == 1)
        {
            // The aggregate is the block projection and the block data is the full input, so the
            // deterministic local objective serves all three reported values without re-evaluation.
            BlockRun run = RunBlock(0, 0, data, targetDim, objective, maxIters, seed, annealerOptions, cancellationToken);
            var blocks = new[]
            {
                new DistributedSpredBlockResult(
                    run.Index, run.Start, run.Count, run.Seed, run.Projection,
                    run.LocalObjective, run.LocalObjective),
            };
            return new DistributedSpredResult(ambientDim, targetDim, run.Projection, run.LocalObjective, blocks);
        }

        var blockRuns = new BlockRun[blockCount];
        RunBlocks(blockCount, maxDegreeOfParallelism, cancellationToken, block =>
        {
            int start = block * data.Length / blockCount;
            int end = (block + 1) * data.Length / blockCount;
            double[][] slice = SliceRows(data, start, end);
            int? blockSeed = BlockSeed(seed, block);
            try
            {
                blockRuns[block] = RunBlock(
                    block,
                    start,
                    slice,
                    targetDim,
                    objective,
                    maxIters,
                    blockSeed,
                    annealerOptions,
                    cancellationToken);
            }
            catch (Exception failure) when (failure is not OperationCanceledException)
            {
                throw BlockSetupFailure(block, slice.Length, failure);
            }
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

    /// <summary>
    /// Run <paramref name="body"/> once per block index — serially when
    /// <paramref name="maxDegreeOfParallelism"/> is one, otherwise via
    /// <see cref="Parallel.For(int, int, ParallelOptions, Action{int})"/>. Both paths share one failure
    /// surface: the first body exception is rethrown with its original stack (never
    /// <see cref="AggregateException"/>), and cancellation surfaces as <see cref="OperationCanceledException"/>.
    /// </summary>
    internal static void RunBlocks(
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

        try
        {
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
        catch (AggregateException aggregate)
        {
            ExceptionDispatchInfo.Capture(aggregate.InnerExceptions[0]).Throw();
        }
    }

    private static InvalidOperationException BlockSetupFailure(int index, int rowCount, Exception failure)
    {
        return new InvalidOperationException(
            $"Distributed SPRED block {index} ({rowCount} rows) failed objective construction or annealing; " +
            "the block may be too small for the objective's graph recipe.",
            failure);
    }

    private static BlockRun RunBlock(
        int index,
        int start,
        double[][] data,
        int targetDim,
        PersistenceObjectiveConfig objective,
        int maxIters,
        int? seed,
        SubspaceAnnealerOptions? annealerOptions,
        CancellationToken cancellationToken)
    {
        cancellationToken.ThrowIfCancellationRequested();
        var ph = new PersistenceObjective(data, objective);
        cancellationToken.ThrowIfCancellationRequested();
        SubspaceAnnealerResult annealed = SubspaceAnnealer.Compute(
            data,
            targetDim,
            ph.Evaluate,
            maxIters,
            seed,
            annealerOptions,
            cancellationToken);
        cancellationToken.ThrowIfCancellationRequested();
        return new BlockRun(index, start, data.Length, seed, annealed.Projection, annealed.Objective, ph);
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

        var grass = new GrassmannManifold(ambientDim, targetDim);

        // Warm-start at the medoid frame rather than frames[0]: an arbitrary leading block can sit
        // at principal angle exactly pi/2 from the clean majority (Y^T Z singular — the Grassmann
        // cut locus, where LogMap degenerates and Weiszfeld cannot move off the initialization;
        // see GeometricMedian.MedoidIndex). Also makes the aggregate invariant to block order.
        double[] median = (double[])frames[GeometricMedian.MedoidIndex(grass, frames)].Clone();
        double[] weights = new double[frames.Length];
        Array.Fill(weights, 1.0);

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
