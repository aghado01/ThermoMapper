#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Ordered cells, each carrying its explicit Z/2 boundary as indices into earlier cells.
/// This serves as the currency for Zigzag Z0 and acts as a standard-reduction target for Z2.
/// </summary>
public sealed class CellFiltration : IFiltration
{
    private readonly int[] _dimensions;
    private readonly double[] _births;
    private readonly int[][] _boundaries;
    private readonly int[][] _coboundaries;

    public int Count => _dimensions.Length;
    public string Label { get; }
    public bool EmergentPairs => false; // CellFiltration doesn't support emergent pairs optimization by default

    public CellFiltration(IReadOnlyList<(int Dimension, double Birth, int[] Boundary)> cells, string label = "CellFiltration")
    {
        Label = label;
        int n = cells.Count;
        _dimensions = new int[n];
        _births = new double[n];
        _boundaries = new int[n][];

        // First pass: copy cells
        for (int i = 0; i < n; i++)
        {
            _dimensions[i] = cells[i].Dimension;
            _births[i] = cells[i].Birth;
            _boundaries[i] = cells[i].Boundary ?? Array.Empty<int>();
        }

        // Second pass: compute coboundaries
        var coboundaryLists = new List<int>[n];
        for (int i = 0; i < n; i++)
        {
            coboundaryLists[i] = new List<int>();
        }

        for (int i = 0; i < n; i++)
        {
            foreach (int face in _boundaries[i])
            {
                if (face >= 0 && face < n)
                {
                    coboundaryLists[face].Add(i);
                }
            }
        }

        _coboundaries = new int[n][];
        for (int i = 0; i < n; i++)
        {
            _coboundaries[i] = coboundaryLists[i].ToArray();
        }
    }

    public int GetDimension(int simplexIndex) => _dimensions[simplexIndex];

    public double GetBirth(int simplexIndex) => _births[simplexIndex];

    public ReadOnlySpan<int> GetVertices(int simplexIndex)
    {
        // CellFiltration relies on pure Z/2 boundaries and doesn't store explicit vertex string-keys.
        // Consumers requesting explicit vertices (e.g. for LMP cycle reconstruction) are generally unsupported
        // for pure abstract cells. We return empty here.
        return ReadOnlySpan<int>.Empty;
    }

    public int[] GetBoundaryIndices(int simplexIndex) => _boundaries[simplexIndex];

    public int[] GetCoboundaryIndices(int simplexIndex) => _coboundaries[simplexIndex];
}
