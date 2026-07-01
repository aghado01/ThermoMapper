using System;

namespace Graphs.Primitives;

public struct Neighbor : IComparable<Neighbor>
{
    public int Index;
    public double Distance;

    public int CompareTo(Neighbor other) => other.Distance.CompareTo(Distance);
}
