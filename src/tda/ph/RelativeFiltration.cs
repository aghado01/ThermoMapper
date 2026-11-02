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
}
