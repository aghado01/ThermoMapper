using System;

namespace Graphs.Primitives.Mst;

public readonly record struct MstEdge(int U, int V, double Weight) : IComparable<MstEdge>
{
    public int CompareTo(MstEdge other) => Weight.CompareTo(other.Weight);
}
