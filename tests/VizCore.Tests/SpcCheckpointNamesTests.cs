using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Execution.Sinks;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Models.Potts;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// A checkpoint's name carries T in its round-trip form, so distinct temperatures get distinct files at every scale. The
/// 5-decimal form gave 117 of the 289 points on the stopping schedule a shared name: a task then loaded another
/// temperature's frame as its own, or threw at the second writer.
/// </summary>
public sealed class SpcCheckpointNamesTests
{
    private static readonly double[] Schedule = SpcScheduleHelpers.LogSpaceGrid(1e-6, 0.3, 289);
    private static readonly RunBudget Budget = new(2, 4);
    private const int Replicas = 2;

    [Fact]
    public void TheName_CarriesTExactly()
    {
        foreach (double t in Schedule)
        {
            string name = SpcExecutor.GetCheckpointFileName(t, replicaIndex: 3);
            string carried = name["T_".Length..name.IndexOf("_rep_", StringComparison.Ordinal)];
            Assert.Equal(t, double.Parse(carried, CultureInfo.InvariantCulture));
        }
    }

    [Fact]
    public void AdjacentDoubles_GetDistinctNames()
        => Assert.NotEqual(SpcExecutor.GetCheckpointFileName(1e-6, 0), SpcExecutor.GetCheckpointFileName(Math.BitIncrement(1e-6), 0));

    [Fact]
    public void Executor_WritesAndResumesEachTasksOwnFrame()
    {
        string directory = NewDirectory();
        try
        {
            var executor = new SpcExecutor();
            List<SpcTaskSpec> tasks = executor.BuildTaskList(Schedule, Replicas, q: 4, AccumulationSpec.None, Budget, directory);
            var options = new SpcExecutionOptions { BaseSeed = 42, FrameSink = SpcxDiskFrameSink.Instance };

            Assert.Equal(tasks.Count, tasks.Select(t => t.CheckpointPath).Distinct().Count());

            SpcBatchResult first = executor.RunBatch(SwAccumulatorProbeTests.Triangle(), tasks, default, options);
            Assert.Equal(tasks.Count, first.ExecutedTaskCount);
            Assert.Equal(tasks.Count, Directory.GetFiles(directory, "*.spcx").Length);
            AssertEachFrameIsItsTasks(tasks, first.AlignedRuns);

            SpcBatchResult resumed = executor.RunBatch(SwAccumulatorProbeTests.Triangle(), tasks, default, options);
            Assert.Equal(tasks.Count, resumed.SkippedTaskCount);
            AssertEachFrameIsItsTasks(tasks, resumed.AlignedRuns);
        }
        finally
        {
            Directory.Delete(directory, recursive: true);
        }
    }

    [Fact]
    public void FixedGridSweep_WritesAndResumesEachTasksOwnFrame()
    {
        string directory = NewDirectory();
        try
        {
            var config = new FixedGridSweepConfig
            {
                Temperatures = Schedule, Replicas = Replicas, SweepBudget = Budget, EquilibriumBudget = Budget,
                Sampler = new PottsModelConfig { Q = 4 }, BaseSeed = 42, CheckpointDirectory = directory,
            };
            // The strategy lists its tasks replica by replica.
            var expected = (from r in Enumerable.Range(0, Replicas)
                            from t in Schedule
                            select (Temperature: t, Replica: r)).ToList();

            SweepResult first = new FixedGridSweepStrategy(config).Run(SwAccumulatorProbeTests.Triangle());
            SweepResult resumed = new FixedGridSweepStrategy(config).Run(SwAccumulatorProbeTests.Triangle());

            // One file per sweep task, and one for the equilibrium pass.
            Assert.Equal(expected.Count + 1, Directory.GetFiles(directory, "*.spcx").Length);
            Assert.Equal(expected, first.SweepRuns.Select(r => (r.Accumulator.Temperature, r.Accumulator.ReplicaIndex)));
            Assert.Equal(expected, resumed.SweepRuns.Select(r => (r.Accumulator.Temperature, r.Accumulator.ReplicaIndex)));
            Assert.Equal(first.SweepRuns.Select(r => r.Accumulator.RngState0), resumed.SweepRuns.Select(r => r.Accumulator.RngState0));
        }
        finally
        {
            Directory.Delete(directory, recursive: true);
        }
    }

    private static void AssertEachFrameIsItsTasks(IReadOnlyList<SpcTaskSpec> tasks, IReadOnlyList<SpcRunResult?> runs)
    {
        Assert.Equal(tasks.Count, runs.Count);
        for (int i = 0; i < tasks.Count; i++)
        {
            Assert.Equal(tasks[i].Temperature, runs[i]!.Accumulator.Temperature);
            Assert.Equal(tasks[i].ReplicaIndex, runs[i]!.Accumulator.ReplicaIndex);
        }
    }

    private static string NewDirectory()
    {
        string directory = Path.Combine(Path.GetTempPath(), "spc-checkpoint-names-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(directory);
        return directory;
    }
}
