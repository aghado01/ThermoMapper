using System;
using System.Collections.Generic;
using System.Linq;

namespace Graphs.Pipeline;

/// <summary>
/// The undirected edges the repair stage added: present after repair, in neither direction
/// before it. The scalers keep them past the coupling floor and Local Mutual Proximity
/// leaves them their kernel coupling, so a repaired graph stays connected.
/// </summary>
public sealed class RepairEdgeSet
{
    private readonly HashSet<long> _keys;
    private readonly int _n;
    private IReadOnlySet<(int Lo, int Hi)>? _pairs;

    private RepairEdgeSet(HashSet<long> keys, int n)
    {
        _keys = keys;
        _n = n;
    }

    /// <summary>The edges in <paramref name="after"/> that <paramref name="before"/> holds in neither direction.</summary>
    public static RepairEdgeSet Between(NeighborSelection before, NeighborSelection after, int n)
    {
        var existing = new HashSet<long>();
        for (int i = 0; i < n; i++)
            foreach (var nb in before.AllNeighbors[i])
                existing.Add(KeyOf(i, nb.Index, n));

        var added = new HashSet<long>();
        for (int i = 0; i < n; i++)
            foreach (var nb in after.AllNeighbors[i])
            {
                long key = KeyOf(i, nb.Index, n);
                if (!existing.Contains(key)) added.Add(key);
            }

        return new RepairEdgeSet(added, n);
    }

    public int Count => _keys.Count;

    public bool Contains(int a, int b) => _keys.Count > 0 && _keys.Contains(KeyOf(a, b, _n));

    /// <summary>The edges as normalized <c>(Lo, Hi)</c> vertex pairs, the form the LMP protect set takes.</summary>
    public IReadOnlySet<(int Lo, int Hi)> Pairs =>
        _pairs ??= new HashSet<(int Lo, int Hi)>(_keys.Select(key => ((int)(key / Math.Max(_n, 1)), (int)(key % Math.Max(_n, 1)))));

    /// <summary>An edge's key, <c>lo · n + hi</c>; the same for either orientation.</summary>
    public static long KeyOf(int a, int b, int n) => a < b ? (long)a * n + b : (long)b * n + a;

    public IReadOnlyCollection<long> Keys => _keys;
}
