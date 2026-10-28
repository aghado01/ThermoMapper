using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Solver;

/// <summary>
/// Result of a focused PKWang temperature sweep: the per-temperature partitions
/// produced by one prepared context, plus the configuration that produced them
/// (declared/traced provenance). PKWang-native — it carries no Swendsen–Wang
/// thermodynamic observables (χ, specific heat, magnetization), which are moot by
/// construction for PKWang: there is no global spin configuration and no Monte
/// Carlo to average.
/// </summary>
public sealed record PKWangSweepResult
{
    /// <summary>Field the sweep was run with.</summary>
    public required Field Field { get; init; }

    /// <summary>Symmetrization rule (inert for the symmetric MeanField).</summary>
    public required SymmetrizationRule Symmetrization { get; init; }

    /// <summary>Bond-activity threshold applied at each temperature.</summary>
    public required double Theta { get; init; }

    /// <summary>Temperatures in caller order, parallel to <see cref="Partitions"/>.</summary>
    public required double[] Temperatures { get; init; }

    /// <summary>Partition at each temperature, parallel to <see cref="Temperatures"/>.</summary>
    public required Assignment[] Partitions { get; init; }

    /// <summary>Cluster count at each temperature (fresh array per call).</summary>
    public int[] ClusterCounts()
    {
        var counts = new int[Partitions.Length];
        for (int i = 0; i < Partitions.Length; i++)
            counts[i] = Partitions[i].Count;
        return counts;
    }
}
