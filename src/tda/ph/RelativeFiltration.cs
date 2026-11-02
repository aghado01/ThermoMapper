#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// The relative chain complex of a filtration against a fixed subcomplex <c>K₀</c>: the cells of <c>K</c>
/// outside <c>K₀</c>, in <c>K</c>'s order, each boundary keeping only its faces outside <c>K₀</c>. Reducing it
/// is the persistence of the pair <c>H(K_δ, K₀)</c>, the features <c>K₀</c> does not explain, where a read of
/// <c>K</c> alone also reports <c>K₀</c>'s own. A chord whose endpoints both lie in <c>K₀</c> has an empty
/// relative boundary, so it is a cycle by itself; a <c>K₀</c> cycle that <c>K</c> fills reappears one degree up,
/// born where it is filled.
/// <para>For a connected spanning <c>K₀</c>, <c>H₁(K, K₀) ≅ H₁(K) / im H₁(K₀)</c> (Azriel's derivation from the
/// pair's long exact sequence), so the relative H₁ bars are the absolute ones born after <c>K₀</c>. This is fixed
/// <c>K₀</c>, not the moving <c>H(K, K_i)</c> of de Silva, Morozov and Vejdemo-Johansson's dualities: an ordinary
/// filtration of the quotient, read by <see cref="PersistentHomology"/> with no new reduction code. The clearing
/// and involuted paths assume absolute structure, so they do not apply.</para>
/// </summary>
public sealed class RelativeFiltration : IFiltration
{
    readonly bool[] _inK0;
    readonly int[] _ambientIndex;
    readonly int[][] _boundaries;
    readonly int[][] _coboundaries;

    /// <summary>The filtration of <c>K</c> this complex is relative to <c>K₀</c> within.</summary>
    public SimplicialFiltration Ambient { get; }

    /// <summary>For each relative cell, its index in <see cref="Ambient"/>.</summary>
    public IReadOnlyList<int> AmbientIndices => _ambientIndex;

    public int Count => _ambientIndex.Length;
    public string Label { get; }
    public bool EmergentPairs => false;

    /// <param name="ambient">The filtration of <c>K</c>.</param>
    /// <param name="inSubcomplex">Membership in <c>K₀</c>. Throws when it is not closed under taking faces.</param>
    /// <param name="label">Axis label of a reduction over this complex.</param>
    public RelativeFiltration(SimplicialFiltration ambient, Func<Simplex, bool> inSubcomplex, string label = "relative")
    {
        ArgumentNullException.ThrowIfNull(ambient);
        ArgumentNullException.ThrowIfNull(inSubcomplex);

        Ambient = ambient;
        Label = label;

        int m = ambient.Count;
        var inK0 = new bool[m];
        for (int a = 0; a < m; a++)
            inK0[a] = inSubcomplex(ambient.Simplices[a]);
        _inK0 = inK0;

        var relativeOf = new int[m];
        var ambientIndex = new List<int>(m);
        for (int a = 0; a < m; a++)
        {
            if (inK0[a])
            {
                relativeOf[a] = -1;
                foreach (int face in ambient.GetBoundaryIndices(a))
                    if (!inK0[face])
                        throw new ArgumentException(
                            $"K₀ is not a subcomplex: a simplex of dimension {ambient.GetDimension(a)} at index {a} "
                            + $"is in K₀ and its face at index {face} is not.", nameof(inSubcomplex));
            }
            else
            {
                relativeOf[a] = ambientIndex.Count;
                ambientIndex.Add(a);
            }
        }

        _ambientIndex = ambientIndex.ToArray();
        _boundaries = new int[_ambientIndex.Length][];
        var coboundaries = new List<int>[_ambientIndex.Length];
        for (int r = 0; r < coboundaries.Length; r++)
            coboundaries[r] = new List<int>();

        for (int r = 0; r < _ambientIndex.Length; r++)
        {
            var kept = new List<int>();
            foreach (int face in ambient.GetBoundaryIndices(_ambientIndex[r]))
            {
                int f = relativeOf[face];
                if (f < 0) continue;
                kept.Add(f);
                coboundaries[f].Add(r);
            }
            _boundaries[r] = kept.ToArray();
        }

        _coboundaries = new int[coboundaries.Length][];
        for (int r = 0; r < coboundaries.Length; r++)
            _coboundaries[r] = coboundaries[r].ToArray();
    }

    public int GetDimension(int simplexIndex) => Ambient.GetDimension(_ambientIndex[simplexIndex]);

    public double GetBirth(int simplexIndex) => Ambient.GetBirth(_ambientIndex[simplexIndex]);

    public ReadOnlySpan<int> GetVertices(int simplexIndex) => Ambient.GetVertices(_ambientIndex[simplexIndex]);

    public int[] GetBoundaryIndices(int simplexIndex) => _boundaries[simplexIndex];

    public int[] GetCoboundaryIndices(int simplexIndex) => _coboundaries[simplexIndex];

    /// <summary>
    /// Reduces this complex through <see cref="PersistentHomology"/>, with each bar's
    /// <see cref="Bar.Generator"/> and <see cref="Bar.Cycle"/> indexing <see cref="Ambient"/>. A
    /// <see cref="Bar.Cycle"/> is a <em>relative</em> cycle, a chain whose boundary lies in <c>K₀</c>, so it need
    /// not close; for a connected spanning <c>K₀</c>, an H₁ bar's generator is the chord that closed the return.
    /// Zero-length bars stay, as the reducer emits them.
    /// </summary>
    public Barcode ComputeBarcode(int maxDimension = int.MaxValue, bool representatives = false)
    {
        Barcode relative = PersistentHomology.Compute(this, maxDimension, representatives);

        var bars = new Bar[relative.Bars.Count];
        for (int b = 0; b < bars.Length; b++)
        {
            Bar bar = relative.Bars[b];
            bars[b] = bar with
            {
                Generator = bar.Generator is int g ? _ambientIndex[g] : null,
                Cycle = bar.Cycle is { } cycle ? Array.ConvertAll(cycle, c => _ambientIndex[c]) : null,
            };
        }

        return new Barcode(bars, Label);
    }

    /// <summary>
    /// Both readings from one absolute run of the clearing engine with <c>K₀</c> strictly first, which
    /// <c>relative_barcode</c> in <c>lean/Lemmas/RelativePersistence.lean</c> proves determines the relative barcode:
    /// <list type="table">
    /// <listheader><term>Absolute bar</term><description>Relative bar</description></listheader>
    /// <item><term>Born outside <c>K₀</c></term><description>The same bar</description></item>
    /// <item><term>Born in <c>K₀</c>, killed outside it</term><description>Essential, born at that death, one degree up</description></item>
    /// <item><term>Born in <c>K₀</c>, essential or killed within it</term><description>None</description></item>
    /// </list>
    /// <c>Relative</c> holds the same pairs as <see cref="ComputeBarcode"/>, each with the same
    /// <see cref="Bar.Generator"/>. This path gives no <see cref="Bar.Cycle"/>: a moved-up class's relative cycle is
    /// the chain <c>V_j</c>, with <c>∂V_j</c> in <c>K₀</c>, and no engine keeps <c>V</c>. Ask
    /// <see cref="ComputeBarcode"/> for representatives.
    /// <para><c>Absolute</c> is <c>K</c> with <c>K₀</c> at the base: its cells read at the ambient's smallest value.
    /// When they already sit there, as a backbone at ε₀ = 0 does, it is the ambient's own barcode.</para>
    /// <para>Internal, so <see cref="ComputeBarcode"/> is the one public path. This is for wanting both: for the relative
    /// barcode alone <see cref="ComputeBarcode"/> reduces a smaller problem, and cost less at the largest size measured.
    /// Only a busy machine has timed this; figures in <c>issues/ph/ph-loose-ends-brief.md</c>.</para>
    /// </summary>
    internal (Barcode Absolute, Barcode Relative) ComputeBothReadings(int maxDimension = int.MaxValue)
    {
        var absolute = new List<Bar>();
        var relative = new List<Bar>();
        double baseValue = Ambient.Count == 0 ? 0.0 : Ambient.GetBirth(0);

        double Value(int cell) =>
            cell < 0 ? double.PositiveInfinity : _inK0[cell] ? baseValue : Ambient.GetBirth(cell);

        foreach (AbsolutePair pair in AbsolutePairsWithK0First(maxDimension))
        {
            var bar = new Bar(Value(pair.Birth), Value(pair.Death), pair.Dimension, pair.Birth);
            absolute.Add(bar);

            if (pair.Row == RelativeRow.Kept)
                relative.Add(bar);
            else if (pair.Row == RelativeRow.MovedUp && pair.Dimension < maxDimension)
                relative.Add(new Bar(Value(pair.Death), double.PositiveInfinity, pair.Dimension + 1, pair.Death));
        }

        return (new Barcode(absolute, Ambient.Label), new Barcode(relative, Label));
    }

    /// <summary>The absolute run with <c>K₀</c> strictly first, as pairs of ambient cells, each with its table row.</summary>
    internal List<AbsolutePair> AbsolutePairsWithK0First(int maxDimension)
    {
        var first = new K0FirstOrder(Ambient, _inK0);
        List<(int Birth, int Death, int Dimension)> pairs = PersistentCohomology.ComputePairs(first, maxDimension);

        var read = new List<AbsolutePair>(pairs.Count);
        foreach ((int birth, int death, int dimension) in pairs)
        {
            int ambientBirth = first.AmbientOf(birth);
            int ambientDeath = death < 0 ? -1 : first.AmbientOf(death);
            RelativeRow row = !_inK0[ambientBirth] ? RelativeRow.Kept
                : ambientDeath >= 0 && !_inK0[ambientDeath] ? RelativeRow.MovedUp
                : RelativeRow.Dropped;
            read.Add(new AbsolutePair(ambientBirth, ambientDeath, dimension, row));
        }

        return read;
    }

    /// <summary>
    /// The ambient's cells in a new order: <c>K₀</c> first, then the rest, each in the ambient's order. Faces stay before
    /// cofaces since <c>K₀</c> is closed under faces. <c>K₀</c>'s cells read just below every value in the ambient, since
    /// the engine's values must not decrease along the order and a bar's cutoff reads them.
    /// </summary>
    sealed class K0FirstOrder : IFiltration
    {
        readonly SimplicialFiltration _ambient;
        readonly int[] _ambientOf;
        readonly int[] _positionOf;
        readonly int _k0Count;
        readonly double _k0Value;

        public K0FirstOrder(SimplicialFiltration ambient, bool[] inK0)
        {
            _ambient = ambient;
            int m = ambient.Count;
            _ambientOf = new int[m];
            _positionOf = new int[m];
            _k0Value = m == 0 ? 0.0 : Math.BitDecrement(ambient.GetBirth(0));

            for (int a = 0; a < m; a++)
                if (inK0[a])
                    _ambientOf[_k0Count++] = a;

            int next = _k0Count;
            for (int a = 0; a < m; a++)
                if (!inK0[a])
                    _ambientOf[next++] = a;

            for (int p = 0; p < m; p++)
                _positionOf[_ambientOf[p]] = p;
        }

        public int AmbientOf(int position) => _ambientOf[position];

        public int Count => _ambientOf.Length;
        public string Label => _ambient.Label;
        public bool EmergentPairs => false;

        public int GetDimension(int position) => _ambient.GetDimension(_ambientOf[position]);

        public double GetBirth(int position) => position < _k0Count ? _k0Value : _ambient.GetBirth(_ambientOf[position]);

        public ReadOnlySpan<int> GetVertices(int position) => _ambient.GetVertices(_ambientOf[position]);

        public int[] GetBoundaryIndices(int position) =>
            Array.ConvertAll(_ambient.GetBoundaryIndices(_ambientOf[position]), a => _positionOf[a]);

        public int[] GetCoboundaryIndices(int position) =>
            Array.ConvertAll(_ambient.GetCoboundaryIndices(_ambientOf[position]), a => _positionOf[a]);
    }
}

/// <summary>Which row of the relative table an absolute pair falls in.</summary>
internal enum RelativeRow { Kept, MovedUp, Dropped }

/// <summary>A pair of the absolute run, as ambient cell indices; <c>Death</c> is −1 for an essential class.</summary>
internal readonly record struct AbsolutePair(int Birth, int Death, int Dimension, RelativeRow Row);
