#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Diagram distances on <see cref="Barcode"/> / <see cref="Bar"/>.
/// Unrelated to <c>Maths.Distance.Geodesic.Wasserstein1</c> (1-D PMF transport).
/// </summary>
public static class DiagramMetrics
{
    internal enum EssentialKind { InfiniteOnMismatch, FinitePenalty }

    /// <summary>
    /// Policy for essential (infinite-death) bars, handled before finite matching.
    /// <see cref="InfiniteOnMismatch"/> is the default via <c>default(EssentialPolicy)</c>.
    /// </summary>
    public readonly record struct EssentialPolicy
    {
        internal EssentialKind Kind { get; init; }
        internal double PerBar { get; init; }

        public static EssentialPolicy InfiniteOnMismatch => default;

        public static EssentialPolicy FinitePenalty(double perBar) =>
            new() { Kind = EssentialKind.FinitePenalty, PerBar = perBar };
    }

    /// <summary>
    /// Wasserstein distance W_p between two barcodes in the given homological dimension.
    /// L∞ ground metric; balanced (n+m) assignment (Kerber–Morozov–Nigmetov).
    /// </summary>
    public static double Wasserstein(
        Barcode a,
        Barcode b,
        int dimension,
        double p = 2.0,
        EssentialPolicy essential = default)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        ValidateP(p);

        int essA = CountEssential(a, dimension);
        int essB = CountEssential(b, dimension);

        if (essential.Kind == EssentialKind.InfiniteOnMismatch && essA != essB)
            return double.PositiveInfinity;

        var finiteA = CollectFinite(a, dimension);
        var finiteB = CollectFinite(b, dimension);

        double[,] cost = BuildCost(finiteA, finiteB, p, out double big);
        int s = finiteA.Count + finiteB.Count;
        double assignmentSum = s == 0 ? 0.0 : MinAssignment(cost, s);

        if (essential.Kind == EssentialKind.FinitePenalty)
        {
            int surplus = Math.Abs(essA - essB);
            if (surplus > 0)
                assignmentSum += surplus * Math.Pow(essential.PerBar, p);
        }

        return Math.Pow(assignmentSum, 1.0 / p);
    }

    /// <summary>Bottleneck distance — deferred; gate uses <see cref="Wasserstein"/>.</summary>
    public static double Bottleneck(
        Barcode a,
        Barcode b,
        int dimension,
        EssentialPolicy essential = default) =>
        throw new NotImplementedException("Bottleneck (Hopcroft–Karp threshold matching) is P1.");

    static void ValidateP(double p)
    {
        if (p < 1.0)
            throw new ArgumentOutOfRangeException(nameof(p), "Wasserstein requires p >= 1.");
        if (double.IsPositiveInfinity(p))
            throw new ArgumentOutOfRangeException(nameof(p), "Use Bottleneck for p = infinity.");
    }

    static int CountEssential(Barcode barcode, int dimension)
    {
        int count = 0;
        foreach (Bar bar in barcode.Bars)
            if (bar.Dimension == dimension && bar.IsInfinite)
                count++;
        return count;
    }

    static List<Bar> CollectFinite(Barcode barcode, int dimension)
    {
        var list = new List<Bar>();
        foreach (Bar bar in barcode.Bars)
            if (bar.Dimension == dimension && !bar.IsInfinite)
                list.Add(bar);
        return list;
    }

    static double Dinf(in Bar p, in Bar q) =>
        Math.Max(Math.Abs(p.Birth - q.Birth), Math.Abs(p.Death - q.Death));

    static double Diag(in Bar p) => 0.5 * (p.Death - p.Birth);

    static double[,] BuildCost(IReadOnlyList<Bar> a, IReadOnlyList<Bar> b, double p, out double big)
    {
        int n = a.Count, m = b.Count, s = n + m;
        var c = new double[s, s];
        double maxReal = 0.0;

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                maxReal = Math.Max(maxReal, Dinf(a[i], b[j]));
        for (int i = 0; i < n; i++)
            maxReal = Math.Max(maxReal, Diag(a[i]));
        for (int j = 0; j < m; j++)
            maxReal = Math.Max(maxReal, Diag(b[j]));

        big = Math.Pow(maxReal + 1.0, p) * (s + 1) + 1.0;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
                c[i, j] = Math.Pow(Dinf(a[i], b[j]), p);
            for (int k = 0; k < n; k++)
                c[i, m + k] = k == i ? Math.Pow(Diag(a[i]), p) : big;
        }

        for (int j = 0; j < m; j++)
        {
            for (int jj = 0; jj < m; jj++)
                c[n + j, jj] = jj == j ? Math.Pow(Diag(b[j]), p) : big;
            for (int k = 0; k < n; k++)
                c[n + j, m + k] = 0.0;
        }

        return c;
    }

    /// <summary>Potential-based Hungarian (Kuhn–Munkres), O(n³). Square n×n, min total cost.</summary>
    static double MinAssignment(double[,] a, int n)
    {
        const double inf = double.PositiveInfinity;
        var u = new double[n + 1];
        var v = new double[n + 1];
        var p = new int[n + 1];
        var way = new int[n + 1];

        for (int i = 1; i <= n; i++)
        {
            p[0] = i;
            int j0 = 0;
            var minv = new double[n + 1];
            var used = new bool[n + 1];
            Array.Fill(minv, inf);

            do
            {
                used[j0] = true;
                int i0 = p[j0];
                int j1 = -1;
                double delta = inf;

                for (int j = 1; j <= n; j++)
                {
                    if (used[j]) continue;
                    double cur = a[i0 - 1, j - 1] - u[i0] - v[j];
                    if (cur < minv[j])
                    {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta)
                    {
                        delta = minv[j];
                        j1 = j;
                    }
                }

                for (int j = 0; j <= n; j++)
                {
                    if (used[j])
                    {
                        u[p[j]] += delta;
                        v[j] -= delta;
                    }
                    else
                    {
                        minv[j] -= delta;
                    }
                }

                j0 = j1;
            } while (p[j0] != 0);

            do
            {
                int j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
            } while (j0 != 0);
        }

        double res = 0.0;
        for (int j = 1; j <= n; j++)
            res += a[p[j] - 1, j - 1];
        return res;
    }
}
