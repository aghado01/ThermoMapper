using System;

namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>
/// Pure resolver from a <see cref="WorkerBudgetPolicy"/> + task count to a concrete
/// <see cref="WorkerBudget"/>. The ambient overload's only environmental read is
/// <see cref="Environment.ProcessorCount"/>; the explicit-core overload is fully deterministic
/// in its inputs, so the truth table is host-independent and unit-testable.
/// </summary>
/// <remarks>
/// This is the planner half of the planner/dumb-runner split: it decides <em>how many</em>
/// workers; the executor's dispatch decides nothing about the count. Auto mode reserves cores
/// (do not hog the machine); a batch smaller than the ceiling takes one worker per task.
/// </remarks>
public static class WorkerBudgetResolver
{
    /// <summary>
    /// Resolve the worker count for a batch of <paramref name="taskCount"/> independent tasks
    /// against the ambient <see cref="Environment.ProcessorCount"/>.
    /// </summary>
    public static WorkerBudget Resolve(int taskCount, WorkerBudgetPolicy policy)
        => Resolve(taskCount, policy, Environment.ProcessorCount);

    /// <summary>
    /// Resolve against an explicit <paramref name="logicalCores"/> count (deterministic, no
    /// environment read). A <paramref name="taskCount"/> of 0 yields 1 (a quiescent pool),
    /// never 0.
    /// </summary>
    public static WorkerBudget Resolve(int taskCount, WorkerBudgetPolicy policy, int logicalCores)
    {
        ArgumentNullException.ThrowIfNull(policy);
        if (taskCount < 0)
            throw new ArgumentOutOfRangeException(nameof(taskCount), taskCount, "Task count cannot be negative.");

        int logical = Math.Max(1, logicalCores);

        WorkerBudgetPolicyKind kind;
        int ceiling;
        if (policy.MaxWorkers is int explicitMax)
        {
            kind = WorkerBudgetPolicyKind.Explicit;
            ceiling = Math.Clamp(explicitMax, 1, logical);
        }
        else
        {
            kind = WorkerBudgetPolicyKind.Auto;
            int reserved = Math.Min(Math.Max(0, policy.ReservedCores), logical - 1); // always leave ≥ 1 core
            ceiling = Math.Max(1, logical - reserved);
        }

        // One worker per task, bounded by the ceiling. No spin-up amortization (the colonel's
        // grade-down): Parallel.For workers are warm pooled threads, so capping a small batch
        // below its task count would only idle cores. An empty batch yields a quiescent 1.
        int workers = taskCount == 0 ? 1 : Math.Min(ceiling, taskCount);

        return new WorkerBudget(workers, kind, taskCount, logical, ceiling);
    }
}
