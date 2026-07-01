using System;
using System.Runtime.InteropServices;

namespace Maths.Distance;

/// <summary>
/// Core distance contract for hot-path graph algorithms (MST, core distances).
/// Graph bandwidth metadata lives on <c>Graphs.Distance.IDistanceMetric</c>.
/// </summary>
public interface IDistanceMetric
{
    double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b);

    double Distance(ref double a, ref double b, int dim)
        => Distance(
            MemoryMarshal.CreateReadOnlySpan(ref a, dim),
            MemoryMarshal.CreateReadOnlySpan(ref b, dim));
}
