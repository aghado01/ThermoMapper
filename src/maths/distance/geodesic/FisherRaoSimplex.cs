using System;
using System.Buffers;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Threading;

namespace Maths.Distance.Geodesic;

public static class FisherRaoSimplex
{
    private static readonly ThreadLocal<double[]> _frScratch =
        new(() => new double[512]);

    public static double Distance(double[] p, double[] q)
        => Distance((ReadOnlySpan<double>)p, q);

    public static double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
    {
        int len = p.Length;
        if (len == 0) return 0.0;
        if (q.Length != len)
            throw new ArgumentException(
                $"FisherRaoSimplex requires equal-length PMF vectors (p={len}, q={q.Length}).");

        if (len <= 64)
        {
            Span<double> buf = stackalloc double[len];
            return SimplexCore(p, q, buf, len);
        }

        if (len <= 512)
            return SimplexCore(p, q, _frScratch.Value!.AsSpan(0, len), len);

        var pool = ArrayPool<double>.Shared;
        double[] rented = pool.Rent(len);
        try
        {
            return SimplexCore(p, q, rented.AsSpan(0, len), len);
        }
        finally
        {
            pool.Return(rented);
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double SimplexCore(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> buf, int len)
    {
        var pSpan = p.Slice(0, len);
        var qSpan = q.Slice(0, len);

        TensorPrimitives.Multiply<double>(pSpan, qSpan, buf);
        TensorPrimitives.Max<double>(buf, 0.0, buf);
        TensorPrimitives.Sqrt<double>(buf, buf);
        double bhat = TensorPrimitives.Sum<double>(buf);

        if (bhat > 1.0) bhat = 1.0;
        else if (bhat < -1.0) bhat = -1.0;
        return 2.0 * Math.Acos(bhat);
    }
}
