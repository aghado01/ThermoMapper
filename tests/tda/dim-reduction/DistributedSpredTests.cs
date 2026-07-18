using System;
using System.Threading;
using Graphs;
using Maths.Geometry;
using Maths.Geometry.DimReduction;
using Maths.Rng;
using Xunit;

namespace TDA.DimReduction.Tests;

public sealed class DistributedSpredTests
{
    [Theory]
    [InlineData(0)]
    [InlineData(9)]
    public void Compute_BlockCountOutsideRowRange_Throws(int blockCount)
    {
        ArgumentOutOfRangeException error = Assert.Throws<ArgumentOutOfRangeException>(() =>
            DistributedSpred.Compute(
                Circle3D(8),
                targetDim: 2,
                blockCount,
                SmallConfig(),
                maxIters: 0,
                seed: 1));

        Assert.Equal("blockCount", error.ParamName);
    }

    [Theory]
    [InlineData(0)]
    [InlineData(4)]
    public void Compute_TargetDimensionOutsideAmbientRange_Throws(int targetDim)
    {
        ArgumentOutOfRangeException error = Assert.Throws<ArgumentOutOfRangeException>(() =>
            DistributedSpred.Compute(
                Circle3D(8),
                targetDim,
                blockCount: 2,
                SmallConfig(),
                maxIters: 0,
                seed: 1));

        Assert.Equal("targetDim", error.ParamName);
    }

    [Fact]
    public void Compute_RaggedInput_Throws()
    {
        double[][] data =
        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0],
        ];

        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            DistributedSpred.Compute(
                data,
                targetDim: 2,
                blockCount: 1,
                SmallConfig(),
                maxIters: 0,
                seed: 1));

        Assert.Equal("data", error.ParamName);
    }

    [Theory]
    [InlineData(-1)]
    [InlineData(0)]
    public void Compute_NonPositiveParallelism_Throws(int maxDegreeOfParallelism)
    {
        ArgumentOutOfRangeException error = Assert.Throws<ArgumentOutOfRangeException>(() =>
            DistributedSpred.Compute(
                Circle3D(8),
                targetDim: 2,
                blockCount: 2,
                SmallConfig(),
                maxIters: 0,
                seed: 1,
                maxDegreeOfParallelism));

        Assert.Equal("maxDegreeOfParallelism", error.ParamName);
    }

    [Fact]
    public void Compute_PreCanceledToken_Throws()
    {
        using var cancellation = new CancellationTokenSource();
        cancellation.Cancel();

        Assert.Throws<OperationCanceledException>(() =>
            DistributedSpred.Compute(
                Circle3D(16),
                targetDim: 2,
                blockCount: 2,
                SmallConfig(),
                maxIters: 8,
                seed: 1,
                maxDegreeOfParallelism: 2,
                cancellationToken: cancellation.Token));
    }

    // A bad options record is a caller error, not a block failure: it must surface raw with the
    // offending property's name, before any block work — never wrapped in the block-context
    // InvalidOperationException.
    [Fact]
    public void Compute_NonFiniteAnnealerOption_FailsFastUnwrapped()
    {
        ArgumentOutOfRangeException error = Assert.Throws<ArgumentOutOfRangeException>(() =>
            DistributedSpred.Compute(
                Circle3D(16),
                targetDim: 2,
                blockCount: 2,
                SmallConfig(),
                maxIters: 1,
                seed: 1,
                maxDegreeOfParallelism: 1,
                new SubspaceAnnealerOptions { CoolingRate = double.NaN }));

        Assert.Equal("CoolingRate", error.ParamName);
    }

    [Fact]
    public void ComputeWithDiagnostics_NonFiniteAnnealerOption_FailsFastUnwrapped()
    {
        ArgumentOutOfRangeException error = Assert.Throws<ArgumentOutOfRangeException>(() =>
            DistributedSpred.ComputeWithDiagnostics(
                Circle3D(16),
                targetDim: 2,
                blockCount: 2,
                SmallConfig(),
                maxIters: 1,
                seed: 1,
                maxDegreeOfParallelism: 1,
                new SubspaceAnnealerOptions { InitialTemperature = double.PositiveInfinity }));

        Assert.Equal("InitialTemperature", error.ParamName);
    }

    // blockCount == data.Length gives 1-row blocks — far below what SmallConfig's kNN (K=6) recipe and
    // the annealer's PCA warm start can support. Validation cannot see the recipe's K, so the failure
    // surfaces at block setup, wrapped with the block index and row count — identically for serial and
    // parallel runs.
    [Theory]
    [InlineData(1)]
    [InlineData(4)]
    public void Compute_BlocksTooSmallForGraphRecipe_ThrowsWithBlockContext(int maxDegreeOfParallelism)
    {
        double[][] data = Circle3D(8);

        InvalidOperationException error = Assert.Throws<InvalidOperationException>(() =>
            DistributedSpred.Compute(
                data,
                targetDim: 2,
                blockCount: data.Length,
                SmallConfig(),
                maxIters: 0,
                seed: 1,
                maxDegreeOfParallelism));

        Assert.Matches(@"block \d+ \(1 rows\)", error.Message);
        Assert.NotNull(error.InnerException);
        if (maxDegreeOfParallelism == 1)
            Assert.Contains("block 0 (1 rows)", error.Message);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(4)]
    public void ComputeWithDiagnostics_BlocksTooSmallForGraphRecipe_ThrowsWithBlockContext(int maxDegreeOfParallelism)
    {
        double[][] data = Circle3D(8);

        InvalidOperationException error = Assert.Throws<InvalidOperationException>(() =>
            DistributedSpred.ComputeWithDiagnostics(
                data,
                targetDim: 2,
                blockCount: data.Length,
                SmallConfig(),
                maxIters: 0,
                seed: 1,
                maxDegreeOfParallelism));

        Assert.Matches(@"block \d+ \(1 rows\)", error.Message);
        Assert.NotNull(error.InnerException);
    }

    // Deterministic mid-run cancellation: the body cancels the shared source and throws with that
    // token, so the parallel path must collapse to OperationCanceledException rather than surfacing
    // Parallel.For's AggregateException.
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void RunBlocks_BodyCancellation_SurfacesOperationCanceledException(int maxDegreeOfParallelism)
    {
        using var cancellation = new CancellationTokenSource();

        OperationCanceledException error = Assert.Throws<OperationCanceledException>(() =>
            DistributedSpred.RunBlocks(blockCount: 4, maxDegreeOfParallelism, cancellation.Token, _ =>
            {
                cancellation.Cancel();
                throw new OperationCanceledException(cancellation.Token);
            }));

        Assert.Equal(cancellation.Token, error.CancellationToken);
    }

    [Fact]
    public void ComputeWithDiagnostics_ParallelMatchesSeededSerialRun()
    {
        const int blockCount = 3;
        double[][] data = RepeatedCircleBlocks(blockCount, pointsPerBlock: 18);
        PersistenceObjectiveConfig config = SmallConfig();

        DistributedSpredResult serial = DistributedSpred.ComputeWithDiagnostics(
            data,
            targetDim: 2,
            blockCount,
            config,
            maxIters: 3,
            seed: 37,
            maxDegreeOfParallelism: 1);
        DistributedSpredResult parallel = DistributedSpred.ComputeWithDiagnostics(
            data,
            targetDim: 2,
            blockCount,
            config,
            maxIters: 3,
            seed: 37,
            maxDegreeOfParallelism: blockCount);

        Assert.InRange(DistanceBetweenProjections(serial.Projection, parallel.Projection), 0.0, 1e-10);
        Assert.Equal(serial.FullDataObjective, parallel.FullDataObjective, precision: 10);
        Assert.Equal(serial.BlockCount, parallel.BlockCount);

        for (int block = 0; block < blockCount; block++)
        {
            DistributedSpredBlockResult expected = serial.Blocks[block];
            DistributedSpredBlockResult actual = parallel.Blocks[block];

            Assert.Equal(expected.Index, actual.Index);
            Assert.Equal(expected.Start, actual.Start);
            Assert.Equal(expected.Count, actual.Count);
            Assert.Equal(expected.Seed, actual.Seed);
            Assert.InRange(DistanceBetweenProjections(expected.Projection, actual.Projection), 0.0, 1e-10);
            Assert.Equal(expected.LocalObjective, actual.LocalObjective, precision: 10);
            Assert.Equal(expected.AggregateObjective, actual.AggregateObjective, precision: 10);
        }
    }

    // One serial-vs-parallel comparison can miss a race that only fires under a particular
    // thread interleaving; five seeded repetitions on a small fixture raise the odds of
    // surfacing one while keeping the fact cheap.
    [Fact]
    public void ComputeWithDiagnostics_RepeatedParallelRuns_StayDeterministic()
    {
        const int blockCount = 3;
        double[][] data = RepeatedCircleBlocks(blockCount, pointsPerBlock: 14);
        PersistenceObjectiveConfig config = SmallConfig();

        for (int rep = 0; rep < 5; rep++)
        {
            int seed = 61 + 17 * rep;

            DistributedSpredResult serial = DistributedSpred.ComputeWithDiagnostics(
                data,
                targetDim: 2,
                blockCount,
                config,
                maxIters: 2,
                seed,
                maxDegreeOfParallelism: 1);
            DistributedSpredResult parallel = DistributedSpred.ComputeWithDiagnostics(
                data,
                targetDim: 2,
                blockCount,
                config,
                maxIters: 2,
                seed,
                maxDegreeOfParallelism: blockCount);

            Assert.InRange(DistanceBetweenProjections(serial.Projection, parallel.Projection), 0.0, 1e-10);
            Assert.Equal(serial.FullDataObjective, parallel.FullDataObjective, precision: 10);

            for (int block = 0; block < blockCount; block++)
            {
                DistributedSpredBlockResult expected = serial.Blocks[block];
                DistributedSpredBlockResult actual = parallel.Blocks[block];

                Assert.Equal(expected.Seed, actual.Seed);
                Assert.InRange(DistanceBetweenProjections(expected.Projection, actual.Projection), 0.0, 1e-10);
                Assert.Equal(expected.LocalObjective, actual.LocalObjective, precision: 10);
                Assert.Equal(expected.AggregateObjective, actual.AggregateObjective, precision: 10);
            }
        }
    }

    [Fact]
    public void AggregateProjections_DuplicateSubspaceWinsMedian()
    {
        double[][] xy = XyPlane();
        double[][] xyRotated =
        [
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ];
        double[][] xz = XzPlane();

        double[][] aggregated = DistributedSpred.AggregateProjections(
            new[] { xy, xyRotated, xz },
            ambientDim: 3,
            targetDim: 2);

        AssertNearPlane(aggregated, xy, 1e-8);
    }

    [Fact]
    public void AggregateProjections_CleanMajorityResistsCorruptedBlocks()
    {
        double[][] xy = XyPlane();
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
        double[][] xz = XzPlane();
        double[][] yz = YzPlane();

        double[][] aggregated = DistributedSpred.AggregateProjections(
            new[] { xy, xyTiltX, xyTiltY, xz, yz },
            ambientDim: 3,
            targetDim: 2);

        double cleanDistance = DistanceToPlane(aggregated, xy);
        double corruptedXzDistance = DistanceToPlane(aggregated, xz);
        double corruptedYzDistance = DistanceToPlane(aggregated, yz);

        Assert.InRange(cleanDistance, 0.0, 0.1);
        Assert.True(cleanDistance < corruptedXzDistance);
        Assert.True(cleanDistance < corruptedYzDistance);
    }

    [Fact]
    public void AggregateProjections_CorruptedFirstBlock_StillFindsCleanMajority()
    {
        double[][] xz = XzPlane();
        double[][] xy = XyPlane();
        double[][] xyRotated =
        [
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ];
        double[][] xyTilted =
        [
            [Math.Cos(0.04), 0.0, Math.Sin(0.04)],
            [0.0, 1.0, 0.0],
        ];

        // Regression guard for the old frames[0] warm start: every clean frame contains
        // the y axis while xz is normal to it, so from xz each clean frame sits at
        // principal angle exactly pi/2 (Y^T Z singular — the Grassmann cut locus, where
        // LogMap degenerates and Weiszfeld stalls on the corrupted initialization).
        // The medoid warm start must land inside the clean majority instead.
        double[][] aggregated = DistributedSpred.AggregateProjections(
            new[] { xz, xy, xyRotated, xyTilted },
            ambientDim: 3,
            targetDim: 2);

        double cleanDistance = DistanceToPlane(aggregated, xy);

        Assert.InRange(cleanDistance, 0.0, 0.1);
        Assert.True(cleanDistance < DistanceToPlane(aggregated, xz));
    }

    [Fact]
    public void AggregateProjections_TiltedCorruptedFirstBlock_ConvergesToCleanMajority()
    {
        // Corrupted first frame tilted 0.3 rad off xz so every principal angle to the
        // clean frames stays strictly below pi/2 — all pairwise logs inside
        // GrassmannManifold.LogMap's documented domain. Weiszfeld converges to the clean
        // majority from any warm start here; contrast with
        // AggregateProjections_CorruptedFirstBlock_StillFindsCleanMajority, whose exactly
        // orthogonal frames sit on the cut locus and require the medoid warm start.
        double[][] xzTilted =
        [
            [1.0, 0.0, 0.0],
            [0.0, Math.Sin(0.3), Math.Cos(0.3)],
        ];
        double[][] xy = XyPlane();
        double[][] xyRotated =
        [
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ];
        double[][] xyTilted =
        [
            [Math.Cos(0.04), 0.0, Math.Sin(0.04)],
            [0.0, 1.0, 0.0],
        ];

        double[][] aggregated = DistributedSpred.AggregateProjections(
            new[] { xzTilted, xy, xyRotated, xyTilted },
            ambientDim: 3,
            targetDim: 2);

        double cleanDistance = DistanceToPlane(aggregated, xy);

        Assert.InRange(cleanDistance, 0.0, 0.1);
        Assert.True(cleanDistance < DistanceToPlane(aggregated, XzPlane()));
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

    // With maxIters: 0 the annealer returns the PCA warm start, so this validates the
    // PCA-per-block + Grassmann-median aggregation path, not the annealed SPRED estimator
    // (the annealed path is covered by ComputeWithDiagnostics_ParallelMatchesSeededSerialRun
    // and ComputeWithDiagnostics_AdversarialLowIteration_RunStaysInterpretable).
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
        AssertNearPlane(projection, XyPlane(), 1e-8);
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

        double[][] xy = XyPlane();
        AssertNearPlane(result.Projection, xy, 1e-8);

        for (int block = 0; block < blockCount; block++)
        {
            DistributedSpredBlockResult info = result.Blocks[block];
            Assert.Equal(block, info.Index);
            Assert.Equal(block * pointsPerBlock, info.Start);
            Assert.Equal(pointsPerBlock, info.Count);
            Assert.Equal(SeedTree.Derive(seed, blockCount)[block], info.Seed);
            AssertOrthonormalRows(info.Projection);
            AssertNearPlane(info.Projection, xy, 1e-8);
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

    // With maxIters: 0 the annealer returns the PCA warm start, so this validates the
    // PCA-per-block + Grassmann-median aggregation path, not the annealed SPRED estimator
    // (the annealed path is covered by ComputeWithDiagnostics_ParallelMatchesSeededSerialRun
    // and ComputeWithDiagnostics_AdversarialLowIteration_RunStaysInterpretable).
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

        double[][] xy = XyPlane();
        double[][] xz = XzPlane();
        double[][] yz = YzPlane();

        double aggregateToClean = DistanceToPlane(result.Projection, xy);
        double aggregateToXz = DistanceToPlane(result.Projection, xz);
        double aggregateToYz = DistanceToPlane(result.Projection, yz);

        Assert.Equal(5, result.BlockCount);
        AssertFiniteFullObjective(result);
        Assert.InRange(aggregateToClean, 0.0, 1e-8);
        Assert.True(aggregateToClean < aggregateToXz);
        Assert.True(aggregateToClean < aggregateToYz);

        for (int block = 0; block < 3; block++)
            AssertNearPlane(result.Blocks[block].Projection, xy, 1e-8);

        AssertNearPlane(result.Blocks[3].Projection, xz, 1e-8);
        AssertNearPlane(result.Blocks[4].Projection, yz, 1e-8);
        Assert.True(DistanceToPlane(result.Blocks[3].Projection, xy) > 0.5);
        Assert.True(DistanceToPlane(result.Blocks[4].Projection, xy) > 0.5);

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

    // With maxIters: 0 the annealer returns the PCA warm start, so this validates the
    // PCA-per-block + Grassmann-median aggregation path, not the annealed SPRED estimator
    // (the annealed path is covered by ComputeWithDiagnostics_ParallelMatchesSeededSerialRun
    // and ComputeWithDiagnostics_AdversarialLowIteration_RunStaysInterpretable).
    [Fact]
    public void ComputeWithDiagnostics_CorruptedFirstBlock_AggregatesCleanMajority()
    {
        const int pointsPerBlock = 24;

        // Regression guard for the old frames[0] warm start: the block PCA frames span
        // exact coordinate planes, so the clean xy frames sit at principal angle exactly
        // pi/2 from the leading xz frame (the cut-locus configuration where Weiszfeld
        // stalls). The medoid warm start must recover the clean majority end to end.
        DistributedSpredResult result = DistributedSpred.ComputeWithDiagnostics(
            CorruptedFirstCircleBlocks(pointsPerBlock),
            targetDim: 2,
            blockCount: 5,
            SmallConfig(),
            maxIters: 0,
            seed: 59);

        double[][] xy = XyPlane();
        double[][] xz = XzPlane();
        double[][] yz = YzPlane();

        double aggregateToClean = DistanceToPlane(result.Projection, xy);

        Assert.Equal(5, result.BlockCount);
        AssertFiniteFullObjective(result);
        Assert.InRange(aggregateToClean, 0.0, 0.1);
        Assert.True(aggregateToClean < DistanceToPlane(result.Projection, xz));
        Assert.True(aggregateToClean < DistanceToPlane(result.Projection, yz));

        AssertNearPlane(result.Blocks[0].Projection, xz, 1e-8);
        for (int block = 1; block < 4; block++)
            AssertNearPlane(result.Blocks[block].Projection, xy, 1e-8);
        AssertNearPlane(result.Blocks[4].Projection, yz, 1e-8);
    }

    // With maxIters: 0 the annealer returns the PCA warm start, so this validates the
    // PCA-per-block + Grassmann-median aggregation path, not the annealed SPRED estimator
    // (the annealed path is covered by ComputeWithDiagnostics_ParallelMatchesSeededSerialRun
    // and ComputeWithDiagnostics_AdversarialLowIteration_RunStaysInterpretable).
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

        Assert.InRange(DistanceBetweenProjections(global, distributed.Projection), 0.0, 1e-8);
    }

    // With maxIters: 0 the annealer returns the PCA warm start, so this validates the
    // PCA-per-block + Grassmann-median aggregation path, not the annealed SPRED estimator
    // (the annealed path is covered by ComputeWithDiagnostics_ParallelMatchesSeededSerialRun
    // and ComputeWithDiagnostics_AdversarialLowIteration_RunStaysInterpretable).
    [Fact]
    public void ComputeWithDiagnostics_CorruptedFixture_DocumentsGlobalVsDistributedTradeoff()
    {
        double[][] data = ScaledCorruptedCircleBlocks(pointsPerBlock: 24, corruptedRadius: 4.0);
        PersistenceObjectiveConfig config = SmallConfig();

        double[][] global = Spred.Compute(data, targetDim: 2, config, maxIters: 0, seed: 47);
        DistributedSpredResult distributed = DistributedSpred.ComputeWithDiagnostics(
            data,
            targetDim: 2,
            blockCount: 5,
            config,
            maxIters: 0,
            seed: 47);

        double globalObjective = EvaluateFullDataObjective(data, config, global);

        AssertOrthonormalRows(global);
        AssertOrthonormalRows(distributed.Projection);
        Assert.True(double.IsFinite(globalObjective));
        AssertFiniteFullObjective(distributed);

        double[][] xy = XyPlane();
        double globalToClean = DistanceToPlane(global, xy);
        double distributedToClean = DistanceToPlane(distributed.Projection, xy);

        Assert.InRange(distributedToClean, 0.0, 1e-8);
        Assert.True(globalToClean > 0.5);
        Assert.True(distributedToClean < globalToClean);
        Assert.True(globalObjective >= 0.0);
    }

    [Fact]
    [Trait("Category", "SeededAnnealer")]
    public void ComputeWithDiagnostics_AdversarialLowIteration_RunStaysInterpretable()
    {
        double[][] data = CorruptedCircleBlocks(pointsPerBlock: 20);
        PersistenceObjectiveConfig config = SmallConfig();

        DistributedSpredResult result = DistributedSpred.ComputeWithDiagnostics(
            data,
            targetDim: 2,
            blockCount: 5,
            config,
            maxIters: 6,
            seed: 53);

        double[][] xy = XyPlane();
        double[][] xz = XzPlane();
        double[][] yz = YzPlane();

        Assert.Equal(5, result.BlockCount);
        AssertOrthonormalRows(result.Projection);
        AssertFiniteFullObjective(result);

        double aggregateToClean = DistanceToPlane(result.Projection, xy);
        Assert.True(aggregateToClean < DistanceToPlane(result.Projection, xz));
        Assert.True(aggregateToClean < DistanceToPlane(result.Projection, yz));

        foreach (DistributedSpredBlockResult block in result.Blocks)
        {
            AssertOrthonormalRows(block.Projection);
            AssertFiniteObjectives(block);
        }

        Assert.True(result.Blocks[3].AggregateObjective > result.Blocks[3].LocalObjective);
        Assert.True(result.Blocks[4].AggregateObjective > result.Blocks[4].LocalObjective);
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

    private static double[][] CorruptedFirstCircleBlocks(int pointsPerBlock)
    {
        var data = new double[5 * pointsPerBlock][];
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 2), data, 0, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 1), data, 1, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 1), data, 2, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 1), data, 3, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 1, axisB: 2), data, 4, pointsPerBlock);
        return data;
    }

    private static double[][] ScaledCorruptedCircleBlocks(int pointsPerBlock, double corruptedRadius)
    {
        var data = new double[5 * pointsPerBlock][];
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 1), data, 0, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 1), data, 1, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 1), data, 2, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 0, axisB: 2, corruptedRadius), data, 3, pointsPerBlock);
        CopyBlock(CircleInPlane(pointsPerBlock, axisA: 1, axisB: 2, corruptedRadius), data, 4, pointsPerBlock);
        return data;
    }

    private static void CopyBlock(double[][] block, double[][] data, int blockIndex, int pointsPerBlock)
    {
        Array.Copy(block, 0, data, blockIndex * pointsPerBlock, pointsPerBlock);
    }

    private static double[][] CircleInPlane(int n, int axisA, int axisB, double radius = 1.0)
    {
        var pts = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double t = 2.0 * Math.PI * i / n;
            pts[i] = new double[3];
            pts[i][axisA] = radius * Math.Cos(t);
            pts[i][axisB] = radius * Math.Sin(t);
        }
        return pts;
    }

    private static double[][] XyPlane() =>
    [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ];

    private static double[][] XzPlane() =>
    [
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ];

    private static double[][] YzPlane() =>
    [
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];

    private static void AssertNearPlane(double[][] projection, double[][] plane, double maxDistance)
    {
        Assert.InRange(DistanceToPlane(projection, plane), 0.0, maxDistance);
    }

    private static double DistanceToPlane(double[][] projection, double[][] plane)
    {
        return DistanceBetweenProjections(projection, plane);
    }

    private static double DistanceBetweenProjections(double[][] a, double[][] b)
    {
        var grass = new GrassmannManifold(ambientN: a[0].Length, subspaceR: a.Length);
        return grass.Distance(PackFrame(a), PackFrame(b));
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
