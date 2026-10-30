namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Result of one Swendsen–Wang sampler execution — the run's final
/// <see cref="Accumulator"/> (scalar moments, the per-draw record + resume state,
/// plus the optional arrays the run tracked).
/// </summary>
public sealed record SwRunResult
{
    /// <summary>The sufficient-statistics the run folded its draws into.</summary>
    public required Accumulator Accumulator { get; init; }

    /// <summary>
    /// The FK susceptibility measured by this run result.
    /// </summary>
    public double FkSusceptibility
        => Accumulator.RunningSumSqClusterSizes / (Accumulator.DrawCount * (double)Accumulator.Spins.Length);
}
