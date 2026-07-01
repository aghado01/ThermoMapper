using System;
using System.Collections.Generic;

namespace Clustering.Primitives;

/// <summary>
/// The clustering-layer currency: a resolved assignment of points to clusters — the universal output
/// any clustering algorithm (statistical / graphical / geometric) reduces to. <see cref="Labels"/>[i]
/// is the cluster id of point <c>i</c>, dense in <c>[0, Count)</c>, or <see cref="Unassigned"/> when the
/// resolution strategy declined to place the point. A crisp, disjoint, fully-covering assignment is a
/// <i>partition</i> — the special case where no label is <see cref="Unassigned"/>.
/// </summary>
/// <remarks>
/// <para><b>Unassigned is first-class.</b> Any resolution strategy may abstain on a point it cannot
/// confidently place (a frustrated/boundary point) — HDBSCAN's noise generalized across the framework.
/// That outcome is <see cref="Unassigned"/>.</para>
///
/// <para><b>Storage vs semantics.</b> The store is a flat <c>int[]</c> — HPC- and interop-friendly, with
/// <see cref="Unassigned"/> = <c>-1</c> following the sklearn/MATLAB convention (a negative sentinel
/// fails loud if indexed naively, rather than silently mis-attributing). The "unassigned is explicit"
/// semantics live in this type's API (<see cref="IsAssigned"/>, <see cref="Coverage"/>,
/// <see cref="Assigned"/>) — not in a per-element <c>int?</c> (which doubles memory, kills vectorization,
/// and doesn't even close the <c>null == null</c> comparison trap). Consumers iterate the assigned subset
/// rather than hand-rolling the sentinel, mirroring MATLAB's <c>~Missing</c> / <c>PrivX</c> pattern.</para>
/// </remarks>
public sealed record Assignment
{
    /// <summary>Sentinel label for a point the resolution strategy declined to assign.</summary>
    public const int Unassigned = -1;

    /// <summary>Per-point cluster id; dense in <c>[0, Count)</c>, or <see cref="Unassigned"/>.</summary>
    public required int[] Labels { get; init; }

    /// <summary>Number of clusters; valid labels live in <c>[0, Count)</c>.</summary>
    public required int Count { get; init; }

    /// <summary>Total number of points (assigned + unassigned).</summary>
    public int PointCount => Labels.Length;

    /// <summary>True when point <paramref name="i"/> carries a real cluster (not <see cref="Unassigned"/>).</summary>
    public bool IsAssigned(int i) => Labels[i] != Unassigned;

    /// <summary>Number of points that received a real cluster. O(n) scan.</summary>
    public int AssignedCount
    {
        get
        {
            int n = 0;
            for (int i = 0; i < Labels.Length; i++)
                if (Labels[i] != Unassigned) n++;
            return n;
        }
    }

    /// <summary>
    /// Fraction of points assigned, in <c>[0, 1]</c> (0 for an empty assignment) — DeWolfe's coverage,
    /// reported separately so "refused 40%, nailed the 60%" and "assigned everything at 60%" don't collapse.
    /// </summary>
    public double Coverage => Labels.Length == 0 ? 0.0 : (double)AssignedCount / Labels.Length;

    /// <summary>
    /// The assigned <c>(index, label)</c> pairs — the unassigned-safe way to iterate (MATLAB's
    /// <c>~Missing</c> subset). Skips every <see cref="Unassigned"/> point.
    /// </summary>
    public IEnumerable<(int Index, int Label)> Assigned
    {
        get
        {
            for (int i = 0; i < Labels.Length; i++)
                if (Labels[i] != Unassigned)
                    yield return (i, Labels[i]);
        }
    }

    /// <summary>
    /// Build an assignment from dense labels, deriving <see cref="Count"/> as <c>max(label) + 1</c> over the
    /// assigned labels (0 when all unassigned). The caller asserts labels are dense in
    /// <c>[0, Count)</c> ∪ {<see cref="Unassigned"/>}.
    /// </summary>
    public static Assignment FromLabels(int[] labels)
    {
        ArgumentNullException.ThrowIfNull(labels);
        int max = -1;
        for (int i = 0; i < labels.Length; i++)
            if (labels[i] > max) max = labels[i];
        return new Assignment { Labels = labels, Count = max + 1 };
    }
}
