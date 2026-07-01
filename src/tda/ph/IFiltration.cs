#nullable enable
using System;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Pure filtration seam — coboundary/boundary enumeration and births for PH engines.
/// Graph-coupled implementors (lazy Rips) live in <c>tda</c>; explicit complexes here.
/// </summary>
public interface IFiltration
{
    int Count { get; }
    string Label { get; }
    bool EmergentPairs { get; }

    int GetDimension(int simplexIndex);
    double GetBirth(int simplexIndex);
    ReadOnlySpan<int> GetVertices(int simplexIndex);
    int[] GetBoundaryIndices(int simplexIndex);
    int[] GetCoboundaryIndices(int simplexIndex);
}
