#nullable enable
using System.Collections.Generic;

namespace TDA.Ph;

/// <summary>
/// The <c>V</c> of a reduction <c>R = ∂V</c>, kept as the earlier columns each column took, and read a column at a time. A
/// run keeps a few integers per column for this and pays for a chain only when it is read, in time set by the cells in it.
/// </summary>
internal sealed class ReductionChains
{
    readonly int[]?[] _took;

    public ReductionChains(int count) => _took = new int[]?[count];

    public void Record(int column, List<int> took)
    {
        if (took.Count > 0)
            _took[column] = took.ToArray();
    }

    /// <summary>The column's chain, ascending: the cells whose boundaries sum to its reduced column. Its last cell is the column.</summary>
    public int[] ChainOf(int column)
    {
        // V_j = e_j + the sum of V_i over the columns j took, all earlier. So once the later columns are expanded, the latest
        // one pending has its final coefficient, and it is pending exactly when its own cell is in the chain.
        var pending = new SortedSet<int> { column };
        var chain = new List<int>();
        while (pending.Count > 0)
        {
            int j = pending.Max;
            pending.Remove(j);
            chain.Add(j);

            if (_took[j] is { } took)
                foreach (int i in took)
                {
                    if (!pending.Remove(i)) pending.Add(i);
                }
        }

        chain.Reverse();
        return chain.ToArray();
    }
}
