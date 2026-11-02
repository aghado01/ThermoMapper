#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace TDA.Ph.Tests;

/// <summary>A subspace of Z/2 chains over global cell ids, kept reduced by each vector's highest bit.</summary>
sealed class ChainSpan
{
    readonly Dictionary<int, BigInteger> _byTop = new();

    public int Dimension => _byTop.Count;

    public IEnumerable<BigInteger> Basis => _byTop.Values;

    BigInteger Reduce(BigInteger v)
    {
        while (!v.IsZero)
        {
            int top = (int)v.GetBitLength() - 1;
            if (!_byTop.TryGetValue(top, out BigInteger b)) return v;
            v ^= b;
        }
        return v;
    }

    public bool Contains(BigInteger v) => Reduce(v).IsZero;

    public void Add(BigInteger v)
    {
        v = Reduce(v);
        if (!v.IsZero) _byTop[(int)v.GetBitLength() - 1] = v;
    }

    public ChainSpan Plus(ChainSpan other)
    {
        var sum = new ChainSpan();
        foreach (BigInteger v in Basis) sum.Add(v);
        foreach (BigInteger v in other.Basis) sum.Add(v);
        return sum;
    }

    /// <summary>This space intersected with the chains that avoid one cell.</summary>
    public ChainSpan Avoiding(int cell)
    {
        var kept = new ChainSpan();
        BigInteger? holder = null;
        foreach (BigInteger v in Basis)
        {
            if (!ZigzagModuleOracle.HasCell(v, cell)) kept.Add(v);
            else if (holder is null) holder = v;
            else kept.Add(v ^ holder.Value);
        }
        return kept;
    }
}

/// <summary>
/// The zigzag module H_p of a filtration, by exact Z/2 linear algebra over the complex at each level, sharing no code
/// with the engines. A subspace of H_p(K_ℓ) is kept as the cycles it holds, a space between B(K_ℓ) and Z(K_ℓ).
/// <para>The right filtration at a level (Carlsson and de Silva 2010, mirrored) says where each class there ends: walking
/// left from the last level, a forward arrow puts its kernel first and pulls the rest back, and a backward arrow pushes
/// every step forward and puts what is outside the image last. A class's label is the first step that holds it. A class
/// with label d has a compatible sequence up to level d that ends there as Dey and Hou's Definition 5 asks, and none that
/// runs further or ends sooner, so this is the "set of representatives" condition without building the sequence.</para>
/// </summary>
sealed class ZigzagModuleOracle
{
    readonly ZigzagFiltration _f;
    readonly int _p;
    readonly Dictionary<int, int[]>[] _cells;
    readonly Dictionary<int, int>[] _dims;
    readonly ChainSpan[] _z;
    readonly ChainSpan[] _b;
    readonly List<(int Label, ChainSpan Space)>[] _right;

    public static bool HasCell(BigInteger chain, int cell) => !((chain >> cell) & BigInteger.One).IsZero;

    public static BigInteger ChainOf(IEnumerable<int> cells)
    {
        BigInteger chain = BigInteger.Zero;
        foreach (int c in cells)
        {
            if (c < 0) throw new ArgumentException($"{c} is not a cell id.", nameof(cells));
            chain ^= BigInteger.One << c;
        }
        return chain;
    }

    public ZigzagModuleOracle(ZigzagFiltration f, int p)
    {
        _f = f;
        _p = p;
        int m = f.Count;
        _cells = new Dictionary<int, int[]>[m + 1];
        _dims = new Dictionary<int, int>[m + 1];
        _z = new ChainSpan[m + 1];
        _b = new ChainSpan[m + 1];
        _right = new List<(int, ChainSpan)>[m + 1];

        _cells[0] = new Dictionary<int, int[]>();
        _dims[0] = new Dictionary<int, int>();
        for (int s = 0; s < m; s++)
        {
            ZigzagStep step = f[s];
            _cells[s + 1] = new Dictionary<int, int[]>(_cells[s]);
            _dims[s + 1] = new Dictionary<int, int>(_dims[s]);
            if (step.Direction == ZigzagDirection.Add)
            {
                int[] boundary = step.BoundaryAtAdd!;
                _cells[s + 1][step.GlobalCellId] = boundary;
                _dims[s + 1][step.GlobalCellId] = boundary.Length == 0 ? 0 : _dims[s][boundary[0]] + 1;
            }
            else
            {
                _cells[s + 1].Remove(step.GlobalCellId);
                _dims[s + 1].Remove(step.GlobalCellId);
            }
        }

        for (int l = 0; l <= m; l++)
        {
            _z[l] = Cycles(l);
            _b[l] = Boundaries(l);
        }

        _right[m] = new List<(int, ChainSpan)> { (-1, _b[m]) };
        Push(_right[m], m, _z[m]);
        for (int l = m - 1; l >= 0; l--)
            _right[l] = WalkLeft(l);
    }

    public int LastLevel => _f.Count;

    public ChainSpan Z(int level) => _z[level];

    public ChainSpan B(int level) => _b[level];

    public bool Holds(int level, int cell) => _cells[level].ContainsKey(cell) && _dims[level][cell] == _p;

    public BigInteger BoundaryOf(int level, BigInteger chain)
    {
        BigInteger boundary = BigInteger.Zero;
        foreach (int cell in _cells[level].Keys)
            if (HasCell(chain, cell))
                boundary ^= ChainOf(_cells[level][cell]);
        return boundary;
    }

    /// <summary>Where a class at this level ends: −1 for the zero class, otherwise the last level of its bar.</summary>
    public int EndOf(int level, BigInteger cycle)
    {
        foreach (var (label, space) in _right[level])
            if (space.Contains(cycle))
                return label;

        throw new ArgumentException($"Not a cycle of level {level}.");
    }

    /// <summary>The barcode by definition: each step that creates a class, and where its interval ends.</summary>
    public (double Birth, double Death)[] Intervals()
    {
        var bars = new List<(double, double)>();
        for (int s = 0; s < _f.Count; s++)
        {
            ZigzagStep step = _f[s];
            int level = s + 1;
            int cell = step.GlobalCellId;
            if (step.Direction == ZigzagDirection.Add)
            {
                if (_dims[level][cell] != _p || _z[level].Dimension == _z[s].Dimension) continue;

                // The first step of the right filtration holding a cycle through the new cell, the one class the image misses.
                foreach (var (label, space) in _right[level])
                    if (space.Basis.Any(v => HasCell(v, cell)))
                    {
                        bars.Add((s, label));
                        break;
                    }
            }
            else
            {
                if (_dims[s][cell] != _p + 1) continue;

                BigInteger boundary = ChainOf(_cells[s][cell]);
                if (!_b[level].Contains(boundary))
                    bars.Add((s, EndOf(level, boundary)));
            }
        }

        return bars.OrderBy(x => x).ToArray();
    }

    /// <summary>
    /// Whether a bar's cycle is a representative of it: null if so, otherwise the first condition it fails. The necessary
    /// conditions come first; the last is the strong one, that the class runs the bar's whole length and ends at its end.
    /// </summary>
    public string? Check(Bar bar, BigInteger cycle)
    {
        string? necessary = CheckAtBirth(bar, cycle);
        if (necessary is not null) return necessary;

        int level = (int)bar.Birth + 1;
        int end = EndOf(level, cycle);
        return end == (int)bar.Death ? null : $"its class ends at level {end}, not {bar.Death}";
    }

    /// <summary>The conditions that read only the birth step, which alone cannot tell a class from one that outlives it.</summary>
    public string? CheckAtBirth(Bar bar, BigInteger cycle)
    {
        int s = (int)bar.Birth;
        int level = s + 1;

        for (int cell = 0; cell < (int)cycle.GetBitLength(); cell++)
            if (HasCell(cycle, cell) && !Holds(level, cell))
                return $"cell {cell} is not a {_p}-cell of level {level}";

        if (!BoundaryOf(level, cycle).IsZero) return "it does not close";
        if (_b[level].Contains(cycle)) return "it bounds at its birth level";

        if (bar.BirthEnd == IntervalEnd.Closed)
        {
            if (_f[s].Direction != ZigzagDirection.Add) return "a closed birth on a deletion";
            if (_z[s].Plus(_b[level]).Contains(cycle)) return "its class comes from the previous level";
        }
        else
        {
            if (_f[s].Direction != ZigzagDirection.Delete) return "an open birth on an addition";
            if (!_b[s].Contains(cycle)) return "its class survives into the previous level";
        }

        return null;
    }

    List<(int Label, ChainSpan Space)> WalkLeft(int l)
    {
        var filtration = new List<(int, ChainSpan)> { (-1, _b[l]) };
        ZigzagStep step = _f[l];

        if (step.Direction == ZigzagDirection.Add)
        {
            int cell = step.GlobalCellId;
            bool leavesCycles = _dims[l + 1][cell] == _p;
            foreach (var (label, space) in _right[l + 1])
                Push(filtration, label < 0 ? l : label, leavesCycles ? space.Avoiding(cell) : space);
        }
        else
        {
            foreach (var (label, space) in _right[l + 1])
                if (label >= 0)
                    Push(filtration, label, space.Plus(_b[l]));
            Push(filtration, l, _z[l]);
        }

        return filtration;
    }

    static void Push(List<(int Label, ChainSpan Space)> filtration, int label, ChainSpan space)
    {
        if (space.Dimension > filtration[^1].Space.Dimension)
            filtration.Add((label, space));
    }

    ChainSpan Cycles(int level)
    {
        var cycles = new ChainSpan();
        var pivots = new Dictionary<int, (BigInteger Boundary, BigInteger Chain)>();
        foreach (var (cell, boundary) in _cells[level])
        {
            if (_dims[level][cell] != _p) continue;

            BigInteger b = ChainOf(boundary);
            BigInteger chain = BigInteger.One << cell;
            while (!b.IsZero)
            {
                int top = (int)b.GetBitLength() - 1;
                if (!pivots.TryGetValue(top, out var pivot)) break;
                b ^= pivot.Boundary;
                chain ^= pivot.Chain;
            }

            if (b.IsZero) cycles.Add(chain);
            else pivots[(int)b.GetBitLength() - 1] = (b, chain);
        }

        return cycles;
    }

    ChainSpan Boundaries(int level)
    {
        var boundaries = new ChainSpan();
        foreach (var (cell, boundary) in _cells[level])
            if (_dims[level][cell] == _p + 1)
                boundaries.Add(ChainOf(boundary));
        return boundaries;
    }
}
