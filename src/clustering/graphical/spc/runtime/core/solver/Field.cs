namespace Clustering.Graphical.SPC.Runtime.Core.Solver;

/// <summary>
/// The PKWang field sub-axis: how the cumulative-energy ladder that drives the
/// survival kernel is pooled.
/// </summary>
/// <remarks>
/// <see cref="Mean"/> is one global ladder (Wang 2020) — density-blind: a dense
/// region's couplings crowd the ladder and lift the effective cut everywhere.
/// <see cref="Local"/> is a per-site ladder — locally adaptive, robust to
/// density variation (not frustration). This is a fidelity ordering, not a type
/// name; see <c>memory/project_wang2020_spc.md</c>.
/// </remarks>
public enum Field
{
    /// <summary>One global ascending sort + cumulative sum (Wang 2020 mean field).</summary>
    Mean,

    /// <summary>Per-site ascending sort + cumulative sum (density-adaptive).</summary>
    Local,
}
