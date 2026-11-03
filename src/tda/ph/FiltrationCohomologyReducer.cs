#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>Dim-by-dim persistent cohomology reduction with H0 clearing (Ripser-style).</summary>
internal sealed class FiltrationCohomologyReducer
{
    private readonly IFiltration _filtration;
    private readonly int _n;
    private readonly Dictionary<int, int> _pivotOwner = new();
    private SortedSet<int>[] _reduced;
    private readonly HashSet<int> _reducedKeys = new();
    private readonly HashSet<int> _deathSimplices = new();
    private List<int> _columnsToReduce;

    public FiltrationCohomologyReducer(
        IFiltration filtration,
        HashSet<int> toReduce)
    {
        _filtration = filtration;
        _n = filtration.Count;
        _reduced = new SortedSet<int>[_n];
        _columnsToReduce = toReduce.ToList();
    }

    public void RunDimension(
        int dim,
        int maxDimension,
        double cutoff,
        bool representatives,
        List<Bar>? bars,
        List<int>? deaths,
        List<int>? infiniteBirths)
    {
        var columns = new List<int>(_columnsToReduce);
        columns.Sort((a, b) => b.CompareTo(a));

        EnsureReducedCapacity(_filtration.Count);

        foreach (int j in columns)
        {
            int dimJ = _filtration.GetDimension(j);

            if (TryEmergentPair(j, dimJ, maxDimension, cutoff, bars, representatives))
                continue;

            // Cohomology is the anti-transpose of the boundary reduction: the pivot is the LOWEST
            // index in the reduced coboundary column (forward-filtration indexing), not the highest.
            // Since simplices sort by (filtration, dim, combinatorial index), col.Min is the
            // earliest-filtration cofacet = the correct death simplex.
            var col = new SortedSet<int>(_filtration.GetCoboundaryIndices(j));
            while (col.Count > 0 && _pivotOwner.TryGetValue(col.Min, out int owner))
            {
                SortedSet<int>? chain = _reduced[owner];
                if (chain is null)
                    break;

                foreach (int r in chain)
                {
                    if (!col.Remove(r)) col.Add(r);
                }
            }

            _reduced[j] = col;

            if (col.Count == 0)
            {
                if (!_deathSimplices.Contains(j)
                    && dimJ <= maxDimension
                    && PersistenceClearing.PassesCutoff(
                        _filtration.GetBirth(j), double.PositiveInfinity, cutoff))
                {
                    infiniteBirths?.Add(j);
                    if (bars != null)
                        bars.Add(new Bar(
                            _filtration.GetBirth(j),
                            double.PositiveInfinity,
                            dimJ, j,
                            CocycleRep(representatives, col, j, infinite: true)));
                }
                continue;
            }

            int pivot = col.Min;
            _pivotOwner[pivot] = j;
            _reducedKeys.Add(pivot);
            _deathSimplices.Add(pivot);

            double birth = _filtration.GetBirth(j);
            double death = _filtration.GetBirth(pivot);
            if (dimJ <= maxDimension)
            {
                // Every death, whatever the cutoff: a boundary reduction needs all of a degree's columns to reach the right
                // pivots, and it applies the cutoff to the bars it finds.
                deaths?.Add(pivot);
                if (bars != null && PersistenceClearing.PassesCutoff(birth, death, cutoff))
                    bars.Add(new Bar(birth, death, dimJ, j, CocycleRep(representatives, col, j, infinite: false)));
            }
        }
    }

    /// <summary>The column a death simplex from <see cref="RunDimension"/> is the pivot of: the bar's birth cell.</summary>
    public int BirthOf(int deathSimplex) => _pivotOwner[deathSimplex];

    void EnsureReducedCapacity(int size)
    {
        if (size <= _reduced.Length)
            return;

        var grown = new SortedSet<int>[size];
        Array.Copy(_reduced, grown, _reduced.Length);
        _reduced = grown;
    }

    bool TryEmergentPair(
        int j,
        int dimJ,
        int maxDimension,
        double cutoff,
        List<Bar>? bars,
        bool representatives)
    {
        if (!_filtration.EmergentPairs)
            return false;

        int[] cob = _filtration.GetCoboundaryIndices(j);
        if (cob.Length == 0)
            return false;

        // The emergent pivot is the min-index cofacet (same min-convention as the main reduction),
        // computed order-independently so this holds for any IFiltration's cofacet ordering.
        int tau = cob[0];
        for (int k = 1; k < cob.Length; k++)
            if (cob[k] < tau) tau = cob[k];

        if (_filtration.GetBirth(tau) != _filtration.GetBirth(j))
            return false;

        if (_reducedKeys.Contains(tau) || _pivotOwner.ContainsKey(tau))
            return false;

        EnsureReducedCapacity(Math.Max(_n, Math.Max(j, tau) + 1));
        _deathSimplices.Add(tau);
        _reducedKeys.Add(tau);
        _pivotOwner[tau] = j;
        // Reduced columns are keyed by COLUMN index (like _reduced[j] = col in the main loop), with
        // _pivotOwner keyed by pivot. Column j's reduced coboundary is exactly {tau}.
        _reduced[j] = new SortedSet<int> { tau };

        if (dimJ <= maxDimension
            && PersistenceClearing.PassesCutoff(_filtration.GetBirth(j), _filtration.GetBirth(tau), cutoff)
            && bars != null)
        {
            bars.Add(new Bar(
                _filtration.GetBirth(j),
                _filtration.GetBirth(tau),
                dimJ,
                j,
                CocycleRep(representatives, _reduced[j], j, infinite: false)));
        }

        return true;
    }

    /// <summary>
    /// Moves to the next degree: its cells that were not cleared as a death. They are found by dimension, not as cofaces of
    /// this degree, since a cell of a general complex can have no faces (a torus's 2-cell over Z/2, a relative cell whose
    /// faces all lie in the subcomplex), and a degree can be empty with cells above it.
    /// </summary>
    public bool AdvanceDimension(int dim, int dimCap)
    {
        if (dim >= dimCap)
            return false;

        var next = new List<int>();
        for (int tau = 0; tau < _n; tau++)
            if (_filtration.GetDimension(tau) == dim + 1 && !_reducedKeys.Contains(tau))
                next.Add(tau);

        _columnsToReduce = next;
        return true;
    }

    static int[]? CocycleRep(bool representatives, SortedSet<int> col, int birthIndex, bool infinite)
    {
        if (!representatives)
            return null;

        if (infinite)
            return new[] { birthIndex };

        var indices = new int[col.Count];
        col.CopyTo(indices);
        return indices;
    }
}
