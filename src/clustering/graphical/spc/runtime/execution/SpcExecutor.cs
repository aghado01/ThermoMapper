using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using Clustering.Graphical.SPC.Runtime.Execution.Sinks;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Execution;

/// <summary>
/// SPC executor — schedule builder, runspace pool manager, and
/// orchestration of per-task <see cref="IFrameSink"/> dispatch. Owns
/// the parallel loop and seed derivation; defers persistence semantics
/// to the configured sink.
/// </summary>
/// <remarks>
/// <para>This class subsumes the old <c>SpcScheduler</c> and
/// <c>SpcRunner</c> execution semantics. Persistence used to be
/// controlled by a boolean (<c>PersistArtifacts</c>); it is now
/// controlled by the sink on
/// <see cref="SpcExecutionOptions.FrameSink"/>. Default sink is
/// <see cref="SpcxDiskFrameSink"/>, which preserves the prior
/// write-checkpoints-to-disk behavior.</para>
/// </remarks>
public sealed class SpcExecutor : ISpcExecutor
{
    /// <summary>
    /// Builds a flat schedule of independent SPC tasks from the requested
    /// temperatures, replicas, cycles, and output directory.
    /// </summary>
    public List<SpcTaskSpec> BuildTaskList(
        IReadOnlyList<double> temperatures,
        int numReplicas,
        int q,
        AccumulationSpec accumulation,
        RunBudget budget,
        string checkpointDirectory,
        int? baseSeed = null)
    {
        if (temperatures is null || temperatures.Count == 0)
            throw new ArgumentException("Temperatures must be non-empty.", nameof(temperatures));
        if (numReplicas <= 0)
            throw new ArgumentOutOfRangeException(nameof(numReplicas), "Replica count must be positive.");
        if (budget.Cycles <= 0)
            throw new ArgumentOutOfRangeException(nameof(budget), "Budget.Cycles must be positive.");
        if (budget.BurnIn < 0)
            throw new ArgumentOutOfRangeException(nameof(budget), "Budget.BurnIn cannot be negative.");
        if (string.IsNullOrWhiteSpace(checkpointDirectory))
            throw new ArgumentException("Checkpoint directory must be provided.", nameof(checkpointDirectory));

        var tasks = new List<SpcTaskSpec>(temperatures.Count * numReplicas);
        foreach (double temperature in temperatures)
        {
            for (int replica = 0; replica < numReplicas; replica++)
            {
                tasks.Add(new SpcTaskSpec
                {
                    Temperature    = temperature,
                    ReplicaIndex   = replica,
                    Budget         = budget,
                    Q              = q,
                    Accumulation   = accumulation,
                    CheckpointPath = Path.Combine(checkpointDirectory, FileNameFor(temperature, replica)),
                    BaseSeed       = baseSeed,
                });
            }
        }

        return tasks;
    }

    /// <summary>
    /// Execute the schedule and dispatch each completed result to the
    /// configured <see cref="IFrameSink"/>. Tasks the sink reports as
    /// already-completed are skipped. The default sink
    /// (<see cref="SpcxDiskFrameSink"/>) writes canonical SPCX/SPCE
    /// files; pass <see cref="NullFrameSink.Instance"/> to skip
    /// persistence.
    /// </summary>
    public void Run(
        CsrGraph graph,
        IReadOnlyList<SpcTaskSpec> tasks,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null)
    {
        if (tasks is null)
            throw new ArgumentNullException(nameof(tasks));
        if (tasks.Count == 0)
            return;

        IFrameSink sink = executionOptions?.FrameSink ?? SpcxDiskFrameSink.Instance;

        WorkerBudget budget = ResolveBudget(tasks.Count, executionOptions);
        DispatchParallel(tasks.Count, CreateParallelOptions(budget.Workers, ct), i =>
        {
            var task = tasks[i];
            if (sink.TaskAlreadyCompleted(task))
                return;

            var result = ExecuteTask(graph, task, executionOptions);
            sink.Accept(task, result);
        });
    }

    /// <summary>
    /// Execute the schedule and return a batch envelope containing
    /// results and execution counts. Sink behavior is controlled by
    /// <see cref="SpcExecutionOptions.FrameSink"/>; the returned
    /// envelope is independent of the sink (task-list-aligned results
    /// for tasks that were executed).
    /// </summary>
    public SpcBatchResult RunBatch(
        CsrGraph graph,
        IReadOnlyList<SpcTaskSpec> tasks,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null)
    {
        if (tasks is null)
            throw new ArgumentNullException(nameof(tasks));

        IFrameSink sink = executionOptions?.FrameSink ?? SpcxDiskFrameSink.Instance;

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var aligned = new SpcRunResult?[tasks.Count];
        // Tracks whether the slot was freshly executed this call (true)
        // vs hydrated from sink cache (false). Slots that stay null
        // remain null in both views.
        var executedThisCall = new bool[tasks.Count];
        var executedCount = 0;
        var skippedCount = 0;

        WorkerBudget budget = ResolveBudget(tasks.Count, executionOptions);
        DispatchParallel(tasks.Count, CreateParallelOptions(budget.Workers, ct), i =>
        {
            var task = tasks[i];

            // Resume path: ask the sink whether the task is cached, and if
            // so try to hydrate the prior result so AlignedRuns is dense.
            // If the hydration fails (e.g. tier-2 sidecar missing despite
            // the spcx existing), fall through to re-execute.
            if (sink.TaskAlreadyCompleted(task))
            {
                SpcRunResult? cached = sink.TryLoad(task, graph);
                if (cached is not null)
                {
                    aligned[i] = cached;
                    Interlocked.Increment(ref skippedCount);
                    return;
                }
            }

            var result = ExecuteTask(graph, task, executionOptions);
            sink.Accept(task, result);
            aligned[i] = result;
            executedThisCall[i] = true;
            Interlocked.Increment(ref executedCount);
        });

        sw.Stop();

        // Dense Runs view: freshly-executed results in task-list order.
        // Preserves the historical "Runs = work done this call" contract;
        // resumed-from-cache slots show up only in AlignedRuns.
        var executedRuns = new List<SpcRunResult>(executedCount);
        for (int i = 0; i < aligned.Length; i++)
        {
            if (executedThisCall[i] && aligned[i] is not null)
                executedRuns.Add(aligned[i]!);
        }

        return new SpcBatchResult
        {
            Runs = executedRuns,
            AlignedRuns = aligned,
            RequestedTaskCount = tasks.Count,
            ExecutedTaskCount = executedCount,
            SkippedTaskCount = skippedCount,
            Elapsed = sw.Elapsed,
            WorkerBudget = budget,
        };
    }

    /// <summary>
    /// Execute the full schedule in memory and return the results.
    /// Pure execution path — bypasses the sink entirely; useful when
    /// checkpoint persistence is not required and the caller wants
    /// task-list-aligned in-memory results.
    /// </summary>
    public IReadOnlyList<SpcRunResult> ExecuteAll(
        CsrGraph graph,
        IReadOnlyList<SpcTaskSpec> tasks,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null)
    {
        if (tasks is null)
            throw new ArgumentNullException(nameof(tasks));

        var results = new SpcRunResult[tasks.Count];
        WorkerBudget budget = ResolveBudget(tasks.Count, executionOptions);
        DispatchParallel(tasks.Count, CreateParallelOptions(budget.Workers, ct), i =>
        {
            results[i] = ExecuteTask(graph, tasks[i], executionOptions);
        });

        return results;
    }

    /// <summary>
    /// Execute the full schedule in memory and return the results for the
    /// provided Potts run specs. Pure execution path — does not invoke
    /// the sink (SwRunSpec carries no checkpoint metadata).
    /// </summary>
    public IReadOnlyList<SpcRunResult> ExecuteAll(
        IReadOnlyList<SwRunSpec> specs,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null)
    {
        if (specs is null)
            throw new ArgumentNullException(nameof(specs));

        var results = new SpcRunResult[specs.Count];
        WorkerBudget budget = ResolveBudget(specs.Count, executionOptions);
        DispatchParallel(specs.Count, CreateParallelOptions(budget.Workers, ct), i =>
        {
            results[i] = ExecuteTask(specs[i], executionOptions);
        });

        return results;
    }

    /// <summary>
    /// Execute a batch of Potts run specs and return an execution
    /// envelope. Pure execution path — does not invoke the sink.
    /// </summary>
    public SpcBatchResult ExecuteBatch(
        IReadOnlyList<SwRunSpec> specs,
        CancellationToken ct = default,
        SpcExecutionOptions? executionOptions = null)
    {
        if (specs is null)
            throw new ArgumentNullException(nameof(specs));

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var results = new SpcRunResult[specs.Count];

        WorkerBudget budget = ResolveBudget(specs.Count, executionOptions);
        DispatchParallel(specs.Count, CreateParallelOptions(budget.Workers, ct), i =>
        {
            results[i] = ExecuteTask(specs[i], executionOptions);
        });

        sw.Stop();
        return new SpcBatchResult
        {
            Runs = results,
            // Every SwRunSpec gets executed (no sink consultation in
            // this path), so AlignedRuns == Runs.
            AlignedRuns = results,
            RequestedTaskCount = specs.Count,
            ExecutedTaskCount = specs.Count,
            SkippedTaskCount = 0,
            Elapsed = sw.Elapsed,
            WorkerBudget = budget,
        };
    }

    /// <summary>
    /// Execute a single schedule task by mapping it to a Potts run spec.
    /// Seed derivation is controlled by <paramref name="executionOptions"/>.
    /// </summary>
    public SpcRunResult ExecuteTask(CsrGraph graph, SpcTaskSpec task, SpcExecutionOptions? executionOptions = null)
    {
        if (task is null)
            throw new ArgumentNullException(nameof(task));

        return ExecuteTask(new SwRunSpec
        {
            Graph        = graph,
            Temperature  = task.Temperature,
            Q            = task.Q,
            Accumulation = task.Accumulation,
            Seed         = ResolveSeed(task.BaseSeed, task.ReplicaIndex, task.Temperature, executionOptions),
            Budget       = task.Budget,
            ReplicaIndex = task.ReplicaIndex,
        }, executionOptions);
    }

    /// <summary>
    /// Execute a single Potts run spec and return the unified SPC result.
    /// Seed derivation is controlled by <paramref name="executionOptions"/>.
    /// </summary>
    public SpcRunResult ExecuteTask(SwRunSpec spec, SpcExecutionOptions? executionOptions = null)
    {
        if (spec is null)
            throw new ArgumentNullException(nameof(spec));

        var resolvedSpec = spec with
        {
            Seed = ResolveSeed(spec.Seed, spec.ReplicaIndex, spec.Temperature, executionOptions),
        };

        var pottsResult = SwRunner.Run(resolvedSpec);
        return new SpcRunResult
        {
            Graph = spec.Graph,
            Accumulator = pottsResult.Accumulator,
        };
    }

    private static WorkerBudget ResolveBudget(int taskCount, SpcExecutionOptions? executionOptions)
        => WorkerBudgetResolver.Resolve(taskCount, executionOptions?.Parallelism ?? new WorkerBudgetPolicy());

    private static ParallelOptions CreateParallelOptions(int workers, CancellationToken ct)
        => new()
        {
            CancellationToken = ct,
            MaxDegreeOfParallelism = workers,
        };

    // Single parallel-dispatch seam: a future cache-aware partitioner or
    // affinity-pinned pool swaps in here without touching the call sites.
    private static void DispatchParallel(int count, ParallelOptions options, Action<int> body)
        => Parallel.For(0, count, options, body);

    private static int? ResolveSeed(
        int? explicitSeed,
        int replicaIndex,
        double temperature,
        SpcExecutionOptions? executionOptions)
    {
        if (explicitSeed.HasValue)
            return explicitSeed;
        if (executionOptions?.BaseSeed is not int baseSeed)
            return null;

        // round=0: SpcExecutor's flat task list has no refinement-round
        // dimension. Shares the canonical formula with the sweep
        // strategies so identical (baseSeed, T, replica) tuples hash to
        // the same seed regardless of which entry point built the task.
        return SpcSeedHelper.Derive(
            baseSeed,
            temperature,
            replica: replicaIndex,
            round: 0,
            quantizationFactor: executionOptions.TemperatureQuantizationFactor);
    }

    public static string GetCheckpointFileName(double temperature, int replicaIndex)
        => $"T_{temperature.ToString("F5", CultureInfo.InvariantCulture)}_rep_{replicaIndex}.spcx";

    private static string FileNameFor(double temperature, int replicaIndex)
        => GetCheckpointFileName(temperature, replicaIndex);
}
