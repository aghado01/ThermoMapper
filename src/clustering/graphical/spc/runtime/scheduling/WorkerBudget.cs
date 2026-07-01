namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>How the worker ceiling was derived.</summary>
public enum WorkerBudgetPolicyKind
{
    /// <summary>Ceiling derived from logical cores minus reserved cores.</summary>
    Auto,

    /// <summary>Ceiling taken from an explicit <see cref="WorkerBudgetPolicy.MaxWorkers"/>.</summary>
    Explicit,
}

/// <summary>
/// Resolved worker budget plus the provenance of the decision — the requested inputs and the
/// intermediate ceiling — so a run can record why it used <c>Workers</c> threads
/// (requested-vs-resolved, per the config-artifact-provenance convention).
/// </summary>
/// <param name="Workers">Final worker count handed to the parallel dispatch. Always ≥ 1.</param>
/// <param name="Policy">Auto vs Explicit ceiling derivation.</param>
/// <param name="TaskCount">Number of tasks the budget was resolved against.</param>
/// <param name="LogicalCores">Logical processor count observed at resolve time.</param>
/// <param name="Ceiling">Upper bound on workers (explicit max, or logical − reserved).</param>
public readonly record struct WorkerBudget(
    int Workers,
    WorkerBudgetPolicyKind Policy,
    int TaskCount,
    int LogicalCores,
    int Ceiling);
