#nullable enable
using System;

namespace Maths.Topology;

/// <summary>
/// A single simplex in a filtered simplicial complex.
/// </summary>
public readonly struct Simplex
{
    public int[] Vertices { get; }
    public double FiltrationValue { get; }
    public int Dimension => Vertices.Length - 1;

    public Simplex(double filtrationValue, params int[] vertices)
    {
        var v = (int[])vertices.Clone();
        Array.Sort(v);
        Vertices = v;
        FiltrationValue = filtrationValue;
    }
}
