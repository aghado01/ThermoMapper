namespace Clustering.Graphical.CumulativeField;

/// <summary>
/// The cumulative-field sub-axis: how the cumulative-energy ladder that drives the
/// survival kernel is pooled.
/// </summary>
/// <remarks>
/// <see cref="Mean"/> is one global ladder — density-blind: a dense region's
/// couplings crowd the ladder and lift the effective cut everywhere.
/// <see cref="Local"/> is a per-site ladder — locally adaptive, robust to
/// density variation (not frustration). This is a fidelity ordering, not a type
/// name.
/// </remarks>
public enum Field
{
    /// <summary>One global ascending sort + cumulative sum (the mean field; global ladder after Wang et al. 2020).</summary>
    Mean,

    /// <summary>Per-site ascending sort + cumulative sum (density-adaptive).</summary>
    Local,
}
