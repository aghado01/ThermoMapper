using System;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Xunit;

namespace Clustering.Graphical.SPC.Tests.Scheduling;

public sealed class WorkerBudgetResolverTests
{
    // Auto mode: one worker per task up to the ceiling (16 logical cores, reserved 2 → 14).
    // Small batches take a worker each (warm pooled threads cost ~nothing to spin); the ceiling
    // bounds oversubscription on large ones. No spin-up grade-down — that was the colonel's
    // PowerShell-runspace artifact and does not transfer to C# TPL.
    [Theory]
    [InlineData(1, 1)]
    [InlineData(4, 4)]
    [InlineData(8, 8)]
    [InlineData(20, 14)]
    [InlineData(60, 14)]
    [InlineData(200, 14)]
    public void Auto_UsesOneWorkerPerTask_UpToCeiling(int taskCount, int expectedWorkers)
    {
        var policy = new WorkerBudgetPolicy { ReservedCores = 2 };

        WorkerBudget budget = WorkerBudgetResolver.Resolve(taskCount, policy, logicalCores: 16);

        Assert.Equal(expectedWorkers, budget.Workers);
        Assert.Equal(WorkerBudgetPolicyKind.Auto, budget.Policy);
    }

    [Fact]
    public void EmptyBatch_YieldsOneWorker_NeverZero()
    {
        WorkerBudget budget = WorkerBudgetResolver.Resolve(0, new WorkerBudgetPolicy(), logicalCores: 16);

        Assert.Equal(1, budget.Workers);
    }

    [Fact]
    public void ExplicitMax_ClampsToLogicalCores()
    {
        var policy = new WorkerBudgetPolicy { MaxWorkers = 100 };

        WorkerBudget budget = WorkerBudgetResolver.Resolve(1000, policy, logicalCores: 8);

        Assert.Equal(WorkerBudgetPolicyKind.Explicit, budget.Policy);
        Assert.Equal(8, budget.Workers);
    }

    [Fact]
    public void ExplicitMax_BypassesReservedCores()
    {
        var policy = new WorkerBudgetPolicy { MaxWorkers = 6, ReservedCores = 2 };

        WorkerBudget budget = WorkerBudgetResolver.Resolve(1000, policy, logicalCores: 16);

        Assert.Equal(6, budget.Workers); // explicit ceiling 6, not the auto ceiling 14
    }

    [Fact]
    public void ReservedCores_NeverStarvesBelowOneCore()
    {
        var policy = new WorkerBudgetPolicy { ReservedCores = 99 };

        WorkerBudget budget = WorkerBudgetResolver.Resolve(1000, policy, logicalCores: 4);

        Assert.Equal(1, budget.Workers); // reserved clamped so ≥ 1 core always remains
    }

    [Fact]
    public void SmallBatch_NeverExceedsTaskCount()
    {
        var policy = new WorkerBudgetPolicy { ReservedCores = 2 };

        WorkerBudget budget = WorkerBudgetResolver.Resolve(3, policy, logicalCores: 16);

        Assert.Equal(3, budget.Workers); // never more workers than tasks, even with ceiling 14
    }

    [Fact]
    public void NegativeTaskCount_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => WorkerBudgetResolver.Resolve(-1, new WorkerBudgetPolicy(), logicalCores: 8));
    }

    [Fact]
    public void Provenance_RecordsInputsAndCeiling()
    {
        var policy = new WorkerBudgetPolicy { ReservedCores = 2 };

        WorkerBudget budget = WorkerBudgetResolver.Resolve(60, policy, logicalCores: 16);

        Assert.Equal(60, budget.TaskCount);
        Assert.Equal(16, budget.LogicalCores);
        Assert.Equal(14, budget.Ceiling);  // logical 16 − reserved 2
        Assert.Equal(14, budget.Workers);  // 60 tasks capped at the ceiling
    }
}
