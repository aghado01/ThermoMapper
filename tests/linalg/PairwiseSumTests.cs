using System;
using System.Numerics;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// Pairwise summation against exact sums: each term is an integer times a power of two, so the sum is exact as a
/// BigInteger. The error stays within the pairwise bound, (Block/4 + ⌈log₂(N/Block)⌉ + 3)·u·Σ|terms|, including on a
/// centred near-constant column, where left-to-right summation's bound, N·u, is a thousand times looser.
/// </summary>
public sealed class PairwiseSumTests
{
    private const double UnitRoundoff = 1.1102230246251565e-16;   // 2⁻⁵³

    [Theory]
    [InlineData(0)]
    [InlineData(1)]
    [InlineData(3)]
    [InlineData(127)]
    [InlineData(128)]
    [InlineData(129)]
    [InlineData(1000)]
    [InlineData(100_003)]
    public void SumAndDot_StayWithinThePairwiseBound(int n)
    {
        var rng = new Xoshiro256PlusPlus(seed: 900 + n);
        var a = new double[n];
        var b = new double[n];
        for (int i = 0; i < n; i++)
        {
            a[i] = (rng.NextDouble() - 0.5) * Math.Pow(10.0, 4.0 * rng.NextDouble() - 2.0);
            b[i] = (rng.NextDouble() - 0.5) * Math.Pow(10.0, 4.0 * rng.NextDouble() - 2.0);
        }

        double sumScale = 0.0, dotScale = 0.0;
        for (int i = 0; i < n; i++)
        {
            sumScale += Math.Abs(a[i]);
            dotScale += Math.Abs(a[i] * b[i]);
        }
        double sumError = Error(PairwiseSum.Of(a), ExactDot(a, null));
        double dotError = Error(PairwiseSum.Dot(a, b), ExactDot(a, b));
        Assert.True(sumError <= Bound(n) * sumScale, $"sum error {sumError:E3} over bound {Bound(n) * sumScale:E3}");
        Assert.True(dotError <= Bound(n) * dotScale, $"dot error {dotError:E3} over bound {Bound(n) * dotScale:E3}");
    }

    [Fact]
    public void NearConstantColumn_IsSummedWithinItsBound()
    {
        // A categorical indicator true for all but its first case, centred: one entry near −1 first, then 99,999
        // equal entries near 10⁻⁵, whose squares round the same way each time they are added to a running sum near 1.
        const int n = 100_000;
        double mean = (n - 1.0) / n;
        var x = new double[n];
        x[0] = -mean;
        for (int i = 1; i < n; i++) x[i] = 1.0 - mean;

        double scale = 0.0;
        for (int i = 0; i < n; i++) scale += x[i] * x[i];
        double error = Error(PairwiseSum.Dot(x, x), ExactDot(x, x));
        Assert.True(error <= Bound(n) * scale, $"error {error:E3} over bound {Bound(n) * scale:E3}");
        Assert.True(Bound(n) <= n * UnitRoundoff / 1000.0);
    }

    private static double Bound(int n) =>
        (PairwiseSum.Block / 4 + Math.Ceiling(Math.Log2(Math.Max(1.0, (double)n / PairwiseSum.Block))) + 3) * UnitRoundoff;

    // Σ a_i·b_i exactly, as Mantissa·2^Exponent; b = null sums a.
    private static (BigInteger Mantissa, int Exponent) ExactDot(double[] a, double[]? b)
    {
        var terms = new (BigInteger M, int E)[a.Length];
        int min = int.MaxValue;
        for (int i = 0; i < a.Length; i++)
        {
            var (ma, ea) = Exact(a[i]);
            var (mb, eb) = b is null ? (BigInteger.One, 0) : Exact(b[i]);
            terms[i] = (ma * mb, ea + eb);
            if (!terms[i].M.IsZero) min = Math.Min(min, terms[i].E);
        }
        if (min == int.MaxValue) return (BigInteger.Zero, 0);
        BigInteger sum = BigInteger.Zero;
        foreach (var (m, e) in terms) if (!m.IsZero) sum += m << (e - min);
        return (sum, min);
    }

    private static (BigInteger Mantissa, int Exponent) Exact(double x)
    {
        long bits = BitConverter.DoubleToInt64Bits(x);
        int biased = (int)((bits >> 52) & 0x7FF);
        long mantissa = bits & 0xF_FFFF_FFFF_FFFF;
        if (biased == 0) biased = 1;
        else mantissa |= 1L << 52;
        var m = new BigInteger(mantissa);
        return (bits < 0 ? -m : m, biased - 1075);
    }

    // |computed − exact|, exactly, then rounded once.
    private static double Error(double computed, (BigInteger M, int E) exact)
    {
        var (mc, ec) = Exact(computed);
        if (mc.IsZero) ec = exact.E;
        if (exact.M.IsZero) exact = (BigInteger.Zero, ec);
        int e = Math.Min(ec, exact.E);
        BigInteger difference = (mc << (ec - e)) - (exact.M << (exact.E - e));
        return Math.ScaleB((double)BigInteger.Abs(difference), e);
    }
}
