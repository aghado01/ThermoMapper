namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>
/// Declarative worker-budget policy for the SPC executor's flat-task pool. Pairs with the
/// pure <see cref="WorkerBudgetResolver"/>; holds no resolved state — only the knobs a caller
/// sets. Supersedes the bare <c>int MaxDegreeOfParallelism</c> the executor used to read,
/// which saturated every logical core with no headroom.
/// </summary>
public sealed record WorkerBudgetPolicy
{
    /// <summary>
    /// Explicit worker ceiling. <see langword="null"/> selects auto mode (with
    /// <see cref="ReservedCores"/> withheld from the logical core count). Clamped to
    /// <c>[1, logical cores]</c> at resolve time.
    /// </summary>
    public int? MaxWorkers { get; init; }

    /// <summary>
    /// Cores withheld from the ceiling in auto mode so a sweep does not saturate the machine.
    /// Default 2 leaves headroom for the caller / OS; set 0 for maximum throughput. Ignored
    /// when <see cref="MaxWorkers"/> is set.
    /// </summary>
    public int ReservedCores { get; init; } = 2;
}
