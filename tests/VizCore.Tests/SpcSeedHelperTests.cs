using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Execution.Sinks;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Models.Potts;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// <see cref="SpcSeedHelper.Derive"/> folds the temperature's own bits, so a grid's temperatures run on distinct
/// streams at every scale. The ⌊1000·T⌋ it folded before gave three pairs among the q-only bracket's 48 points a shared
/// seed, and all 158 points below 10⁻³ on the stopping schedule one seed per replica (stopping primer §13.14).
/// </summary>
public sealed class SpcSeedHelperTests
{
    [Fact]
    public void TheQOnlyBracket_GivesEachTaskItsOwnSeed()
    {
        (double lo, double hi) = BwdPottsCriticalEstimate.Bracket(20);
        AssertEachTaskHasItsOwnSeed(SpcScheduleHelpers.LogSpaceGrid(lo, hi, 48));
    }

    [Theory]
    [InlineData(1e-6, 0.3, 289)]   // the stopping schedule
    [InlineData(1e-6, 1e-3, 97)]   // all of it below the old resolution
    public void AGridReachingAMillionth_GivesEachTaskItsOwnSeed(double lo, double hi, int steps)
        => AssertEachTaskHasItsOwnSeed(SpcScheduleHelpers.LogSpaceGrid(lo, hi, steps));

    private static void AssertEachTaskHasItsOwnSeed(double[] grid)
    {
        foreach (int baseSeed in new[] { 42, 4242, 20261029 })
            foreach (int round in new[] { 0, -1 })
            {
                var seeds = (from t in grid
                             from r in Enumerable.Range(0, 8)
                             select SpcSeedHelper.Derive(baseSeed, t, r, round)).ToList();

                Assert.True(seeds.Distinct().Count() == seeds.Count, $"base seed {baseSeed}, round {round}");
            }
    }

    /// <summary>The literals are an independent evaluation of the fold, outside this code base.</summary>
    [Theory]
    [InlineData(42, 0.1, 0, 0, -1968973061)]
    [InlineData(42, 1e-6, 3, -1, 873334969)]
    [InlineData(-7, 0.00446, 7, 0, 1801028353)]
    [InlineData(20261029, 2.269, 1, 2, -1146591670)]
    public void TheFold_IsPinned(int baseSeed, double temperature, int replica, int round, int expected)
        => Assert.Equal(expected, SpcSeedHelper.Derive(baseSeed, temperature, replica, round));

    [Fact]
    public void AdjacentDoubles_AreDifferentTemperatures()
        => Assert.NotEqual(SpcSeedHelper.Derive(42, 1e-6, 0), SpcSeedHelper.Derive(42, Math.BitIncrement(1e-6), 0));

    [Fact]
    public void NoBaseSeed_MeansNoDerivedSeed()
        => Assert.Null(SpcSeedHelper.Derive(null, 0.1, 0));

    // Far below the triangle's couplings every bond holds, so a shared seed repeats one trajectory exactly, and the runs
    // end on one RNG state. Two entry points feed the helper, and both are read here.
    private static readonly double[] ColdGrid = SpcScheduleHelpers.LogSpaceGrid(1e-6, 1e-4, 9);
    private static readonly RunBudget ColdBudget = new(10, 20);
    private const int ColdReplicas = 2;

    [Fact]
    public void FixedGridSweep_RunsEachColdTaskOnItsOwnStream()
    {
        SweepResult result = new FixedGridSweepStrategy(new FixedGridSweepConfig
        {
            Temperatures = ColdGrid, Replicas = ColdReplicas, SweepBudget = ColdBudget, EquilibriumBudget = ColdBudget,
            Sampler = new PottsModelConfig { Q = 4 }, BaseSeed = 42,
        }).Run(SwAccumulatorProbeTests.Triangle());

        AssertDistinctStreams(result.SweepRuns);
    }

    [Fact]
    public void Executor_RunsEachColdTaskOnItsOwnStream()
    {
        var sink = new InMemoryFrameSink();
        var executor = new SpcExecutor();
        List<SpcTaskSpec> tasks = executor.BuildTaskList(
            ColdGrid, ColdReplicas, q: 4, AccumulationSpec.None, ColdBudget, checkpointDirectory: "unused");

        executor.Run(SwAccumulatorProbeTests.Triangle(), tasks,
            executionOptions: new SpcExecutionOptions { BaseSeed = 42, FrameSink = sink });

        AssertDistinctStreams(sink.Results);
    }

    private static void AssertDistinctStreams(IReadOnlyList<SpcRunResult> runs)
    {
        Assert.Equal(ColdGrid.Length * ColdReplicas, runs.Count);
        var finalStates = runs.Select(r => r.Accumulator)
            .Select(a => (a.RngState0, a.RngState1, a.RngState2, a.RngState3))
            .ToHashSet();
        Assert.Equal(runs.Count, finalStates.Count);
    }
}
