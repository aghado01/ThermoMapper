using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Runtime.Intrinsics;
using System.Threading;
using Maths.LinAlg;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;

namespace Maths.Regression.Tests;

internal interface IBandwidth { static abstract int P { get; } }
internal readonly struct Bw1 : IBandwidth { public static int P => 1; }
internal readonly struct Bw2 : IBandwidth { public static int P => 2; }
internal readonly struct Bw3 : IBandwidth { public static int P => 3; }

/// <summary>
/// The shared factorizations and the reference rates (the brief's <c>machine</c> and <c>kernels:*</c> blocks). Replicas are
/// test-side and compiled fully optimized: the production loop's arithmetic in a contiguous row-major band (bits identical,
/// asserted), a compile-time bandwidth by <c>static abstract</c> (bits identical), and a fixed-order four-lane dot on the
/// same layout through <see cref="PairwiseSum.Dot"/> (a reduction: the bits differ, and the difference is read).
/// </summary>
internal static class CostKernels
{
    internal static double Sink;
    internal static object? ObjectSink;

    // ---- Matrices ----

    /// <summary>A diagonally dominant SPD band in LAPACK lower-band storage <c>band[d, j] = A(j+d, j)</c>.</summary>
    internal static double[,] BandSpd(int n, int p, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var band = new double[p + 1, n];
        var rowSum = new double[n];
        for (int j = 0; j < n; j++)
            for (int d = 1; d <= Math.Min(p, n - 1 - j); d++)
            {
                double v = (2.0 * rng.NextDouble() - 1.0) * Math.Pow(0.7, d);
                band[d, j] = v;
                rowSum[j] += Math.Abs(v);
                rowSum[j + d] += Math.Abs(v);
            }
        for (int j = 0; j < n; j++) band[0, j] = 1.0 + 1.1 * rowSum[j];
        return band;
    }

    internal static double[,] BandToDense(double[,] band, int n, int p)
    {
        var a = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            a[j, j] = band[0, j];
            for (int d = 1; d <= Math.Min(p, n - 1 - j); d++) a[j + d, j] = a[j, j + d] = band[d, j];
        }
        return a;
    }

    internal static double[,] DenseSpd(int m, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var a = new double[m, m];
        var rowSum = new double[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < i; j++)
            {
                double v = (2.0 * rng.NextDouble() - 1.0) / (1.0 + 0.1 * (i - j));
                a[i, j] = a[j, i] = v;
                rowSum[i] += Math.Abs(v);
                rowSum[j] += Math.Abs(v);
            }
        for (int i = 0; i < m; i++) a[i, i] = 1.0 + 1.1 * rowSum[i];
        return a;
    }

    internal static double[] Random(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var v = new double[n];
        for (int i = 0; i < n; i++) v[i] = 2.0 * rng.NextDouble() - 1.0;
        return v;
    }

    // ---- Operation counts (the loops' own bounds) ----

    /// <summary>The 3-operand multiply-subtracts (L·L·D) of the banded LDLᵀ/Cholesky factorization, counted by the loops' bounds.</summary>
    internal static long BandFactorTerms(int n, int p)
    {
        long c = 0;
        for (int j = 0; j < n; j++)
        {
            c += j - Math.Max(0, j - p);
            int dmax = Math.Min(p, n - 1 - j);
            for (int d = 1; d <= dmax; d++) c += Math.Max(0, j - Math.Max(0, j + d - p));
        }
        return c;
    }

    /// <summary>The closed form for an interior column: n·p(p+1)/2. The count differs from it by the first p columns' shortfall only.</summary>
    internal static long BandFactorTermsClosed(int n, int p) => (long)n * p * (p + 1) / 2;

    internal static long BandSolveTerms(int n, int p)
    {
        long c = 0;
        for (int i = 0; i < n; i++) c += i - Math.Max(0, i - p);
        return 2 * c;
    }

    internal static long DenseFactorTerms(int m) { long c = 0; for (int i = 0; i < m; i++) c += (long)i * (i + 1) / 2 + i; return c; }

    // ---- Band replicas: row-major layout, L(i, k) at l[i·p + k − i + p], k ∈ [i − p, i − 1] ----

    internal static void ToRows(double[,] band, int n, int p, double[] l, double[] d)
    {
        Array.Clear(l);
        for (int i = 0; i < n; i++)
        {
            d[i] = band[0, i];
            for (int k = Math.Max(0, i - p); k < i; k++) l[i * p + k - i + p] = band[i - k, k];
        }
    }

    private static bool Row(double[] l, double[] d, int i, int p, ref double logDet)
    {
        int ri = i * p, klo = Math.Max(0, i - p);
        for (int j = klo; j < i; j++)
        {
            int rj = j * p;
            double s = l[ri + j - i + p];
            for (int k = klo; k < j; k++) s -= l[ri + k - i + p] * l[rj + k - j + p] * d[k];
            l[ri + j - i + p] = s / d[j];
        }
        double dj = d[i];
        for (int k = klo; k < i; k++) { double lik = l[ri + k - i + p]; dj -= lik * lik * d[k]; }
        if (!(dj > 0.0 && double.IsFinite(dj))) return false;
        d[i] = dj;
        logDet += Math.Log(dj);
        return true;
    }

    /// <summary>The production loop's arithmetic, in the same order, on a contiguous layout: the factor and ln|A|, strict.</summary>
    internal static bool FactorRows(double[] l, double[] d, int n, int p, out double logDet)
    {
        logDet = 0.0;
        for (int i = 0; i < n; i++)
            if (!Row(l, d, i, p, ref logDet)) return false;
        return true;
    }

    /// <summary>The same factor with the bandwidth a compile-time constant: the interior rows are straight-line code.</summary>
    internal static bool FactorFixed<TP>(double[] l, double[] d, int n, out double logDet) where TP : struct, IBandwidth
    {
        int p = TP.P;
        logDet = 0.0;
        int i = 0;
        for (; i < Math.Min(n, p); i++)
            if (!Row(l, d, i, p, ref logDet)) return false;
        for (; i < n; i++)
        {
            int ri = i * p;
            double dj;
            if (TP.P == 1)
            {
                double l0 = l[ri] / d[i - 1];
                l[ri] = l0;
                dj = d[i] - l0 * l0 * d[i - 1];
            }
            else if (TP.P == 2)
            {
                double d2 = d[i - 2], d1 = d[i - 1];
                double l0 = l[ri] / d2;
                double s1 = l[ri + 1];
                s1 -= l0 * l[(i - 1) * 2 + 1] * d2;
                double l1 = s1 / d1;
                l[ri] = l0; l[ri + 1] = l1;
                dj = d[i];
                dj -= l0 * l0 * d2;
                dj -= l1 * l1 * d1;
            }
            else
            {
                double d3 = d[i - 3], d2 = d[i - 2], d1 = d[i - 1];
                double l0 = l[ri] / d3;
                double s1 = l[ri + 1];
                s1 -= l0 * l[(i - 2) * 3 + 2] * d3;
                double l1 = s1 / d2;
                double s2 = l[ri + 2];
                s2 -= l0 * l[(i - 1) * 3 + 1] * d3;
                s2 -= l1 * l[(i - 1) * 3 + 2] * d2;
                double l2 = s2 / d1;
                l[ri] = l0; l[ri + 1] = l1; l[ri + 2] = l2;
                dj = d[i];
                dj -= l0 * l0 * d3;
                dj -= l1 * l1 * d2;
                dj -= l2 * l2 * d1;
            }
            if (!(dj > 0.0 && double.IsFinite(dj))) return false;
            d[i] = dj;
            logDet += Math.Log(dj);
        }
        return true;
    }

    /// <summary>The same factor with every inner sum a fixed-order four-lane dot (<see cref="PairwiseSum.Dot"/>) against W = L·D.</summary>
    internal static bool FactorVector(double[] l, double[] w, double[] d, int n, int p, out double logDet)
    {
        logDet = 0.0;
        for (int i = 0; i < n; i++)
        {
            int ri = i * p, klo = Math.Max(0, i - p);
            for (int j = klo; j < i; j++)
            {
                int rj = j * p, len = j - klo;
                double s = l[ri + j - i + p];
                if (len > 0) s -= PairwiseSum.Dot(l.AsSpan(ri + klo - i + p, len), w.AsSpan(rj + klo - j + p, len));
                l[ri + j - i + p] = s / d[j];
            }
            int width = i - klo;
            for (int k = klo; k < i; k++) w[ri + k - i + p] = l[ri + k - i + p] * d[k];
            double dj = d[i];
            if (width > 0) dj -= PairwiseSum.Dot(l.AsSpan(ri + klo - i + p, width), w.AsSpan(ri + klo - i + p, width));
            if (!(dj > 0.0 && double.IsFinite(dj))) return false;
            d[i] = dj;
            logDet += Math.Log(dj);
        }
        return true;
    }

    /// <summary>The largest relative difference of two factors' entries, and whether every entry is identical.</summary>
    internal static (double MaxRelative, bool Identical) Compare(double[] x, double[] y)
    {
        double max = 0.0;
        bool same = true;
        for (int i = 0; i < x.Length; i++)
        {
            if (BitConverter.DoubleToInt64Bits(x[i]) != BitConverter.DoubleToInt64Bits(y[i])) same = false;
            double scale = Math.Max(Math.Abs(x[i]), Math.Abs(y[i]));
            if (scale > 0) max = Math.Max(max, Math.Abs(x[i] - y[i]) / scale);
        }
        return (max, same);
    }

    // ---- Dense replicas: packed lower triangle, L(i, j) at i(i+1)/2 + j ----

    private static int RowOf(int i) => i * (i + 1) / 2;

    /// <summary>DenseCholesky's loops, without its condition norm: the packed factor and ln|A|.</summary>
    internal static bool FactorPacked(double[,] a, double[] l, int m, out double logDet)
    {
        logDet = 0.0;
        for (int i = 0; i < m; i++)
        {
            int ri = RowOf(i);
            for (int j = 0; j < i; j++)
            {
                int rj = RowOf(j);
                double sum = a[i, j];
                for (int k = 0; k < j; k++) sum -= l[ri + k] * l[rj + k];
                l[ri + j] = sum / l[rj + j];
            }
            double d = a[i, i];
            for (int k = 0; k < i; k++) d -= l[ri + k] * l[ri + k];
            if (!(d > 0.0) || !double.IsFinite(d)) return false;
            double lii = Math.Sqrt(d);
            l[ri + i] = lii;
            logDet += Math.Log(lii);
        }
        logDet *= 2.0;
        return true;
    }

    /// <summary>The same factor with every row sum a fixed-order four-lane dot (<see cref="PairwiseSum.Dot"/>).</summary>
    internal static bool FactorPackedVector(double[,] a, double[] l, int m, out double logDet, bool machineDependent = false)
    {
        logDet = 0.0;
        for (int i = 0; i < m; i++)
        {
            int ri = RowOf(i);
            for (int j = 0; j < i; j++)
            {
                int rj = RowOf(j);
                double dot = j == 0 ? 0.0 : machineDependent
                    ? MatrixOps.Dot(l.AsSpan(ri, j), l.AsSpan(rj, j))
                    : PairwiseSum.Dot(l.AsSpan(ri, j), l.AsSpan(rj, j));
                l[ri + j] = (a[i, j] - dot) / l[rj + j];
            }
            double sq = i == 0 ? 0.0 : machineDependent ? MatrixOps.Dot(l.AsSpan(ri, i), l.AsSpan(ri, i)) : PairwiseSum.Dot(l.AsSpan(ri, i), l.AsSpan(ri, i));
            double d = a[i, i] - sq;
            if (!(d > 0.0) || !double.IsFinite(d)) return false;
            double lii = Math.Sqrt(d);
            l[ri + i] = lii;
            logDet += Math.Log(lii);
        }
        logDet *= 2.0;
        return true;
    }

    /// <summary><see cref="CholeskyDecomposition"/>'s first pass, in its own arrays: the factor, floored.</summary>
    internal static void CholeskyStep1(double[,] sigma, double[,] l, int dim)
    {
        for (int i = 0; i < dim; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = sigma[i, j];
                for (int k = 0; k < j; k++) sum -= l[i, k] * l[j, k];
                if (i == j)
                {
                    if (sum <= 1e-12) sum = 1e-12;
                    l[i, i] = Math.Sqrt(sum);
                }
                else
                {
                    l[i, j] = sum / l[j, j];
                    l[j, i] = 0.0;
                }
            }
    }

    /// <summary><see cref="CholeskyDecomposition"/>'s second pass: L⁻¹ by forward substitution, which BGP's fit never reads.</summary>
    internal static void CholeskyStep2(double[,] l, double[,] lInv, int dim)
    {
        for (int i = 0; i < dim; i++)
        {
            lInv[i, i] = 1.0 / l[i, i];
            for (int j = 0; j < i; j++)
            {
                double sum = 0.0;
                for (int k = j; k < i; k++) sum -= l[i, k] * lInv[k, j];
                lInv[i, j] = sum / l[i, i];
                lInv[j, i] = 0.0;
            }
        }
    }

    internal static long CholeskyStep1Terms(int n) { long c = 0; for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) c += j; return c; }

    internal static long CholeskyStep2Terms(int n) { long c = 0; for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) c += i - j; return c; }

    // ---- The reference rates ----

    private static double ScalarDot(double[] a, double[] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.Length; i++) s += a[i] * b[i];
        return s;
    }

    private static double VectorDot(double[] a, double[] b)
    {
        int w = Vector<double>.Count, i = 0;
        var acc = Vector<double>.Zero;
        for (; i <= a.Length - w; i += w) acc += new Vector<double>(a, i) * new Vector<double>(b, i);
        double s = Vector.Sum(acc);
        for (; i < a.Length; i++) s += a[i] * b[i];
        return s;
    }

    // The in-box vector kernels (System.Runtime.Intrinsics): the test project does not reference System.Numerics.Tensors, whose
    // TensorPrimitives is built on them.
    private static double Vector256Dot4(double[] a, double[] b)
    {
        var s0 = Vector256<double>.Zero; var s1 = s0; var s2 = s0; var s3 = s0;
        int i = 0;
        for (; i + 15 < a.Length; i += 16)
        {
            s0 += Vector256.LoadUnsafe(ref a[0], (nuint)i) * Vector256.LoadUnsafe(ref b[0], (nuint)i);
            s1 += Vector256.LoadUnsafe(ref a[0], (nuint)(i + 4)) * Vector256.LoadUnsafe(ref b[0], (nuint)(i + 4));
            s2 += Vector256.LoadUnsafe(ref a[0], (nuint)(i + 8)) * Vector256.LoadUnsafe(ref b[0], (nuint)(i + 8));
            s3 += Vector256.LoadUnsafe(ref a[0], (nuint)(i + 12)) * Vector256.LoadUnsafe(ref b[0], (nuint)(i + 12));
        }
        double s = Vector256.Sum((s0 + s1) + (s2 + s3));
        for (; i < a.Length; i++) s += a[i] * b[i];
        return s;
    }

    private static void Vector256Map(double[] x, double[] y, bool exp)
    {
        int i = 0;
        for (; i + 3 < x.Length; i += 4)
        {
            Vector256<double> v = Vector256.LoadUnsafe(ref x[0], (nuint)i);
            (exp ? Vector256.Exp(v) : Vector256.Log(v)).StoreUnsafe(ref y[0], (nuint)i);
        }
        for (; i < x.Length; i++) y[i] = exp ? Math.Exp(x[i]) : Math.Log(x[i]);
    }

    private interface IShape { double Eval(double x); }
    private sealed class Square : IShape { public double Eval(double x) => x * x; }
    private sealed class Cube : IShape { public double Eval(double x) => x * x * x; }

    internal static void MachineBlock()
    {
        CostProfile.MachineRows();
        foreach (int n in new[] { 1000, 10_000, 100_000, 1_000_000 })
        {
            double[] a = Random(n, 1), b = Random(n, 2);
            var parts = new List<CostProfile.Part>
            {
                new("scalar_madd", c => { for (int i = 0; i < c; i++) Sink = ScalarDot(a, b); }, n, "madd"),
                new("vector_madd", c => { for (int i = 0; i < c; i++) Sink = VectorDot(a, b); }, n, "madd"),
                new("pairwise_dot", c => { for (int i = 0; i < c; i++) Sink = PairwiseSum.Dot(a, b); }, n, "madd"),
                new("matrixops_dot", c => { for (int i = 0; i < c; i++) Sink = MatrixOps.Dot(a, b); }, n, "madd"),
                new("vector256_dot_4_accumulators", c => { for (int i = 0; i < c; i++) Sink = Vector256Dot4(a, b); }, n, "madd"),
            };
            CostProfile.Point("machine", $"dot n={n}", parts);
        }

        {
            const int n = 1024;
            double[] x = Random(n, 3).Select(v => 0.5 + 0.5 * Math.Abs(v)).ToArray(), y = new double[n];
            var shapes = new IShape[] { new Square(), new Cube() };
            Func<double, double> square = v => v * v;
            Func<double, double> cube = v => v * v * v;
            Func<double, double>[] delegates = [square, cube];
            var table = new Dictionary<string, int> { ["birth"] = 0, ["death"] = 1, ["relocate"] = 2 };
            string[] keys = ["birth", "death", "relocate", "death"];
            var gate = new Lock();
            long counter = 0;
            var parts = new List<CostProfile.Part>
            {
                new("exp", c => { for (int r = 0; r < c; r++) { double s = 0; for (int i = 0; i < n; i++) s += Math.Exp(-x[i]); Sink = s; } }, n, "call"),
                new("log", c => { for (int r = 0; r < c; r++) { double s = 0; for (int i = 0; i < n; i++) s += Math.Log(x[i]); Sink = s; } }, n, "call"),
                new("sqrt", c => { for (int r = 0; r < c; r++) { double s = 0; for (int i = 0; i < n; i++) s += Math.Sqrt(x[i]); Sink = s; } }, n, "call"),
                new("divide", c => { for (int r = 0; r < c; r++) { double s = 0; for (int i = 0; i < n; i++) s += 1.0 / x[i]; Sink = s; } }, n, "call"),
                new("vector256_exp", c => { for (int r = 0; r < c; r++) { Vector256Map(x, y, exp: true); Sink = y[0]; } }, n, "element"),
                new("vector256_log", c => { for (int r = 0; r < c; r++) { Vector256Map(x, y, exp: false); Sink = y[0]; } }, n, "element"),
                new("delegate_call", c => { for (int r = 0; r < c; r++) { double s = 0; var f = delegates[r & 1]; for (int i = 0; i < n; i++) s += f(x[i]); Sink = s; } }, n, "call"),
                new("interface_call_monomorphic", c => { IShape sh = shapes[0]; for (int r = 0; r < c; r++) { double s = 0; for (int i = 0; i < n; i++) s += sh.Eval(x[i]); Sink = s; } }, n, "call"),
                new("interface_call_bimorphic", c => { for (int r = 0; r < c; r++) { double s = 0; for (int i = 0; i < n; i++) s += shapes[i & 1].Eval(x[i]); Sink = s; } }, n, "call"),
                new("dictionary_string_lookup", c => { for (int r = 0; r < c; r++) { int s = 0; for (int i = 0; i < n; i++) s += table[keys[i & 3]]; Sink = s; } }, n, "lookup"),
                new("lock_enter_exit", c => { for (int r = 0; r < c; r++) { for (int i = 0; i < n; i++) lock (gate) counter++; } Sink = counter; }, n, "lock"),
            };
            CostProfile.Point("machine", "scalar functions n=1024", parts);
        }

        foreach (int len in new[] { 8, 64, 512, 4096 })
            CostProfile.Point("machine", $"allocation double[{len}]", new List<CostProfile.Part>
            {
                new("new_array", c => { for (int i = 0; i < c; i++) ObjectSink = new double[len]; }, 1, "array"),
            });

        StreamRows();
    }

    /// <summary>Memory streams (a copy and a triad over 64 MB arrays, past the L3) at 1, 2, 4 and 8 threads, in GB/s.</summary>
    private static void StreamRows()
    {
        if (CostProfile.FreeMemoryGb() < 1.5) { CostProfile.Row("machine", "stream", "", "skipped", "free memory below 1.5 GB"); return; }
        const int n = 8 << 20;
        double[] a = new double[n], b = Random(n, 5), c = Random(n, 6);
        foreach (int threads in new[] { 1, 2, 4, 8 })
            foreach (string kind in new[] { "copy", "triad" })
            {
                var rates = new List<double>();
                for (int rep = 0; rep < 6; rep++)
                {
                    var barrier = new Barrier(threads);
                    var workers = new Thread[threads];
                    long start = 0, stop = 0;
                    for (int t = 0; t < threads; t++)
                    {
                        int k = t;
                        workers[t] = new Thread(() =>
                        {
                            int lo = (int)((long)k * n / threads), hi = (int)((long)(k + 1) * n / threads);
                            barrier.SignalAndWait();
                            if (k == 0) start = Stopwatch.GetTimestamp();
                            if (kind == "copy") Array.Copy(b, lo, a, lo, hi - lo);
                            else for (int i = lo; i < hi; i++) a[i] = b[i] + 3.0 * c[i];
                            barrier.SignalAndWait();
                            if (k == 0) stop = Stopwatch.GetTimestamp();
                        });
                        workers[t].Start();
                    }
                    foreach (Thread w in workers) w.Join();
                    double seconds = (stop - start) / (double)Stopwatch.Frequency;
                    rates.Add((kind == "copy" ? 16.0 : 24.0) * n / seconds / 1e9);
                }
                if (rates.Count > 1) rates.RemoveAt(0);   // the first touch
                CostProfile.Row("machine", $"stream {kind}", $"threads={threads}", "GB_per_s", CostProfile.Median(rates));
            }
        Sink = a[17];
    }

    // ---- kernels:band and kernels:engines ----

    private static void BandPoint(string block, string name, int n, int p)
    {
        double[,] band = BandSpd(n, p, 100 + n + p);
        var strict = new BandCholesky(n, p, BandFactorization.Ldlt, strict: true);
        var floored = new BandCholesky(n, p, BandFactorization.Ldlt);
        var chol = new BandCholesky(n, p, BandFactorization.Cholesky);
        double[] lt = new double[n * p], dt = new double[n], l = new double[n * p], d = new double[n], w = new double[n * p];
        ToRows(band, n, p, lt, dt);
        double[] rhs = Random(n, 9), x = new double[n];
        long terms = BandFactorTerms(n, p), solveTerms = BandSolveTerms(n, p);
        long closed = BandFactorTermsClosed(n, p);
        if (n >= 4 * p && Math.Abs(terms - closed) > (long)p * p * p + p * p)
            throw new InvalidOperationException($"The count {terms} differs from n·p(p+1)/2 = {closed} by more than the first columns' shortfall (n = {n}, p = {p}).");

        // Replicas against production: bits.
        strict.DecomposeBanded(band);
        Array.Copy(lt, l, l.Length); Array.Copy(dt, d, n);
        if (!FactorRows(l, d, n, p, out double ldRows)) throw new InvalidOperationException("The row-major replica refused an SPD band.");
        if (BitConverter.DoubleToInt64Bits(ldRows) != BitConverter.DoubleToInt64Bits(strict.LogDet))
            throw new InvalidOperationException($"The row-major replica's ln|A| differs from production's ({ldRows:R} against {strict.LogDet:R}).");
        double[] lRows = (double[])l.Clone(), dRows = (double[])d.Clone();
        string fixedBits = "n/a";
        if (p <= 3)
        {
            Array.Copy(lt, l, l.Length); Array.Copy(dt, d, n);
            bool ok = p == 1 ? FactorFixed<Bw1>(l, d, n, out double ld) : p == 2 ? FactorFixed<Bw2>(l, d, n, out ld) : FactorFixed<Bw3>(l, d, n, out ld);
            (double rel, bool same) = Compare(l, lRows);
            fixedBits = ok && same && Compare(d, dRows).Identical ? "identical" : $"differs {rel:E2}";
            if (!(ok && same)) throw new InvalidOperationException($"The compile-time bandwidth replica differs from the row replica (p = {p}, {rel:E2}).");
        }
        Array.Copy(lt, l, l.Length); Array.Copy(dt, d, n); Array.Clear(w);
        FactorVector(l, w, d, n, p, out double ldVec);
        (double relVec, bool sameVec) = Compare(l, lRows);
        double relD = Compare(d, dRows).MaxRelative;
        CostProfile.Row(block, name, "vector_dot_factor", "bits_identical", sameVec ? 1 : 0);
        CostProfile.Row(block, name, "vector_dot_factor", "max_relative_difference_L", relVec);
        CostProfile.Row(block, name, "vector_dot_factor", "max_relative_difference_D", relD);
        CostProfile.Row(block, name, "vector_dot_factor", "relative_difference_logdet", Math.Abs(ldVec - ldRows) / Math.Max(1e-300, Math.Abs(ldRows)));
        CostProfile.Row(block, name, "row_major_factor", "bits_identical_to_production", 1);
        CostProfile.Row(block, name, "compile_time_bandwidth_factor", "bits_vs_row_major", fixedBits);
        CostProfile.Row(block, name, "(count)", "factor_terms", terms);
        CostProfile.Row(block, name, "(count)", "factor_terms_closed_form", closed);
        CostProfile.Row(block, name, "(count)", "solve_terms", solveTerms);

        var parts = new List<CostProfile.Part>
        {
            new("production_ldlt_strict_factor", c => { for (int i = 0; i < c; i++) strict.DecomposeBanded(band); }, terms, "term"),
            new("production_ldlt_floored_factor", c => { for (int i = 0; i < c; i++) floored.DecomposeBanded(band); }, terms, "term"),
            new("production_cholesky_factor", c => { for (int i = 0; i < c; i++) chol.DecomposeBanded(band); }, terms, "term"),
            new("row_major_factor", c => { for (int i = 0; i < c; i++) { Array.Copy(lt, l, l.Length); Array.Copy(dt, d, n); FactorRows(l, d, n, p, out double ld); Sink = ld; } }, terms, "term"),
            new("vector_dot_factor", c => { for (int i = 0; i < c; i++) { Array.Copy(lt, l, l.Length); Array.Copy(dt, d, n); FactorVector(l, w, d, n, p, out double ld); Sink = ld; } }, terms, "term"),
            new("production_solve", c => { for (int i = 0; i < c; i++) { Array.Copy(rhs, x, n); strict.SolveInPlace(x); } }, solveTerms, "term"),
            new("copy_only", c => { for (int i = 0; i < c; i++) { Array.Copy(lt, l, l.Length); Array.Copy(dt, d, n); } }, 0),
        };
        if (p <= 3)
            parts.Add(new("compile_time_bandwidth_factor", c =>
            {
                for (int i = 0; i < c; i++)
                {
                    Array.Copy(lt, l, l.Length); Array.Copy(dt, d, n);
                    double ld;
                    if (p == 1) FactorFixed<Bw1>(l, d, n, out ld); else if (p == 2) FactorFixed<Bw2>(l, d, n, out ld); else FactorFixed<Bw3>(l, d, n, out ld);
                    Sink = ld;
                }
            }, terms, "term"));
        CostProfile.Point(block, name, parts);
    }

    internal static void BandBlock()
    {
        foreach (int n in new[] { 8, 16, 32, 64, 128, 256, 512, 1024 })
            foreach (int p in new[] { 1, 2, 3, 5, 8, 16, 32, 64 })
                if (p < n) BandPoint("kernels:band", $"nu={n} p={p}", n, p);
    }

    internal static void EngineBandsBlock()
    {
        const string block = "kernels:engines";
        foreach ((string name, int n, int p) in new[]
        {
            ("BARS k=4", 8, 3), ("BARS k=16", 20, 3), ("BARS k=32", 36, 3),
            ("BAPS 1-D 25 knots", 29, 3), ("BAPS 1-D 100 knots", 104, 3), ("BAPS 1-D 500 knots", 504, 3),
        })
            BandPoint(block, name, n, p);

        foreach (int knots in new[] { 6, 10 })
        {
            (int nu, int bandwidth) = TensorBand(knots);
            CostProfile.Row(block, $"BAPS tensor {knots} knots", "(count)", "nu", nu);
            CostProfile.Row(block, $"BAPS tensor {knots} knots", "(count)", "bandwidth_read_from_design_and_penalty", bandwidth);
            BandPoint(block, $"BAPS tensor {knots} knots", nu, bandwidth);
        }
    }

    /// <summary>The size and the half-bandwidth BAPS's tensor fit factors: read from <see cref="TensorProductDesign"/> and
    /// <see cref="TensorPenalty"/> as <c>Baps</c> does.</summary>
    internal static (int Nu, int Bandwidth) TensorBand(int knots)
    {
        int g = 3 * (knots + 4) / 2;
        var xs = new double[g * g];
        var ys = new double[g * g];
        for (int a = 0; a < g; a++)
            for (int b = 0; b < g; b++) { xs[a * g + b] = (a + 0.5) / g; ys[a * g + b] = (b + 0.5) / g; }
        var k = new KnotConfig(Enumerable.Range(1, knots).Select(i => i / (knots + 1.0)).ToArray());
        TensorDesign td = TensorProductDesign.Build(new SplineBasis(3), k, xs, new SplineBasis(3), k, ys);
        var penalty = new TensorPenalty(td.NuX, td.NuY);
        return (td.NuX * td.NuY, Math.Max(new BandedDesign(td.Design).Bandwidth, penalty.Bandwidth));
    }

    // ---- kernels:dense, kernels:bgp, kernels:wide ----

    internal static void DenseBlock()
    {
        foreach (int m in new[] { 5, 10, 20, 50, 100, 200, 500, 1000 })
        {
            double[,] a = DenseSpd(m, 200 + m);
            var factor = new DenseCholesky(m);
            var packed = new double[m * (m + 1) / 2];
            var packedV = new double[packed.Length];
            double[] rhs = Random(m, 11), x = new double[m];
            long terms = DenseFactorTerms(m);
            long closedForm = (long)m * (m - 1) * (m + 1) / 6 + (long)m * (m - 1) / 2;   // Σ_i (i(i+1)/2 + i)
            if (terms != closedForm) throw new InvalidOperationException($"Dense count {terms} differs from the closed form {closedForm} at m = {m}.");

            factor.Decompose(a);
            FactorPacked(a, packed, m, out double ld);
            FactorPackedVector(a, packedV, m, out double ldV);
            bool identical = true;
            double maxRel = 0.0;
            for (int i = 0; i < m; i++)
                for (int j = 0; j <= i; j++)
                {
                    double ref0 = factor.Factor(i, j), mine = packed[RowOf(i) + j], vec = packedV[RowOf(i) + j];
                    if (BitConverter.DoubleToInt64Bits(ref0) != BitConverter.DoubleToInt64Bits(mine)) identical = false;
                    maxRel = Math.Max(maxRel, Math.Abs(vec - ref0) / Math.Max(1e-300, Math.Abs(ref0)));
                }
            if (!identical) throw new InvalidOperationException($"The packed replica differs from DenseCholesky at m = {m}.");
            string name = $"m={m}";
            CostProfile.Row("kernels:dense", name, "vector_dot_factor", "max_relative_difference_L", maxRel);
            CostProfile.Row("kernels:dense", name, "vector_dot_factor", "relative_difference_logdet", Math.Abs(ldV - ld) / Math.Abs(ld));
            CostProfile.Row("kernels:dense", name, "(count)", "factor_terms", terms);
            long solveTerms = (long)m * (m - 1) / 2;
            var parts = new List<CostProfile.Part>
            {
                new("production_factor", c => { for (int i = 0; i < c; i++) factor.Decompose(a); }, terms, "term"),
                new("packed_scalar_factor", c => { for (int i = 0; i < c; i++) { FactorPacked(a, packed, m, out double l0); Sink = l0; } }, terms, "term"),
                new("packed_vector_dot_factor", c => { for (int i = 0; i < c; i++) { FactorPackedVector(a, packedV, m, out double l0); Sink = l0; } }, terms, "term"),
                new("packed_matrixops_dot_factor", c => { for (int i = 0; i < c; i++) { FactorPackedVector(a, packedV, m, out double l0, machineDependent: true); Sink = l0; } }, terms, "term"),
                new("production_forward_solve", c => { for (int i = 0; i < c; i++) { Array.Copy(rhs, x, m); factor.ForwardSolveInPlace(x); } }, solveTerms, "term"),
                new("production_back_solve", c => { for (int i = 0; i < c; i++) { Array.Copy(rhs, x, m); factor.BackSolveInPlace(x); } }, solveTerms, "term"),
            };
            CostProfile.Point("kernels:dense", name, parts);
        }
    }

    internal static void BgpCholeskyBlock()
    {
        foreach (int n in new[] { 50, 100, 200, 400, 800, 1600 })
        {
            double[,] a = DenseSpd(n, 300 + n);
            var chol = new CholeskyDecomposition(n);
            var strict = new DenseCholesky(n);
            var l = new double[n, n];
            var lInv = new double[n, n];
            double[] rhs = Random(n, 13);
            long t1 = CholeskyStep1Terms(n), t2 = CholeskyStep2Terms(n);
            chol.Decompose(a);
            CholeskyStep1(a, l, n);
            CholeskyStep2(l, lInv, n);
            CostProfile.Row("kernels:bgp", $"n={n}", "(count)", "step1_terms", t1);
            CostProfile.Row("kernels:bgp", $"n={n}", "(count)", "step2_terms", t2);
            var parts = new List<CostProfile.Part>
            {
                new("production_decompose", c => { for (int i = 0; i < c; i++) chol.Decompose(a); }, t1 + t2, "term"),
                new("replica_step1_factor", c => { for (int i = 0; i < c; i++) CholeskyStep1(a, l, n); }, t1, "term"),
                new("replica_step2_inverse_of_L", c => { for (int i = 0; i < c; i++) CholeskyStep2(l, lInv, n); }, t2, "term"),
                new("dense_cholesky_strict_factor", c => { for (int i = 0; i < c; i++) strict.Decompose(a); }, DenseFactorTerms(n), "term"),
                new("production_solve", c => { for (int i = 0; i < c; i++) Sink = chol.Solve(rhs)[0]; }, 2L * n * (n + 1) / 2, "term"),
                new("new_n_by_n_array", c => { for (int i = 0; i < c; i++) ObjectSink = new double[n, n]; }, 1, "array"),
            };
            CostProfile.Point("kernels:bgp", $"n={n}", parts);
        }
    }

    internal static void WideBlock()
    {
        foreach ((int n, int p) in new[] { (200, 30), (500, 60), (1000, 100), (500, 250) })
        {
            double[,] band = BandSpd(n, p, 400 + n + p);
            double[,] dense = BandToDense(band, n, p);
            var wide = new WideBandCholesky(n, p);
            var narrow = new BandCholesky(n, p, BandFactorization.Cholesky, strict: true);
            var full = new DenseCholesky(n);
            double[] rhs = Random(n, 15), x = new double[n];
            long terms = BandFactorTerms(n, p), solve = BandSolveTerms(n, p);
            wide.Decompose(dense);
            if (!wide.Status.PositiveDefinite) throw new InvalidOperationException($"WideBandCholesky refused an SPD band (n = {n}, p = {p}).");
            string name = $"n={n} p={p}";
            var parts = new List<CostProfile.Part>
            {
                new("wide_band_factor", c => { for (int i = 0; i < c; i++) wide.Decompose(dense); }, terms, "term"),
                new("band_cholesky_factor_from_dense", c => { for (int i = 0; i < c; i++) narrow.Decompose(dense); }, terms, "term"),
                new("band_cholesky_factor_from_band", c => { for (int i = 0; i < c; i++) narrow.DecomposeBanded(band); }, terms, "term"),
                new("dense_cholesky_factor", c => { for (int i = 0; i < c; i++) full.Decompose(dense); }, DenseFactorTerms(n), "term"),
                new("wide_band_solve", c => { for (int i = 0; i < c; i++) { Array.Copy(rhs, x, n); wide.SolveInPlace(x); } }, solve, "term"),
                new("band_cholesky_solve", c => { for (int i = 0; i < c; i++) { Array.Copy(rhs, x, n); narrow.SolveInPlace(x); } }, solve, "term"),
            };
            CostProfile.Point("kernels:wide", name, parts);
        }
    }

    // ---- kernels:bordered, kernels:lowrank, kernels:khatri ----

    internal static void BorderedBlock()
    {
        foreach (int m in new[] { 5, 10, 20, 50, 100, 200 })
        {
            double[,] big = DenseSpd(m + 1, 500 + m);
            var a = new double[m, m];
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) a[i, j] = big[i, j];
            var border = new double[m];
            for (int i = 0; i < m; i++) border[i] = big[m, i];
            double diagonal = big[m, m];
            var shortBorder = new double[m - 1];
            for (int i = 0; i < m - 1; i++) shortBorder[i] = big[m, i < m / 2 ? i : i + 1];
            double[] rhs = Random(m + 1, 17);
            var bc = new BorderedCholesky(m + 2, rebuildInterval: 0);
            bc.Decompose(a);
            int middle = m / 2;
            double mm = (double)m * m;
            var parts = new List<CostProfile.Part>
            {
                new("decompose_strict", c => { for (int i = 0; i < c; i++) bc.Decompose(a); }, DenseFactorTerms(m), "term"),
                new("propose_birth", c => { for (int i = 0; i < c; i++) Sink = bc.ProposeBirth(border, diagonal).LogDetChange; }, 1.5 * mm, "1.5M^2"),
                new("propose_birth_and_quadratic_form", c => { for (int i = 0; i < c; i++) { bc.ProposeBirth(border, diagonal); Sink = bc.ProposedInverseQuadraticForm(rhs); } }, 2.5 * mm, "2.5M^2"),
                new("propose_death", c => { for (int i = 0; i < c; i++) Sink = bc.ProposeDeath(middle).LogDetChange; }, 2.0 * mm, "2M^2"),
                new("propose_change", c => { for (int i = 0; i < c; i++) Sink = bc.ProposeChange(middle, shortBorder, diagonal).LogDetChange; }, 3.5 * mm, "3.5M^2"),
                new("birth_accept_death_accept", c =>
                {
                    for (int i = 0; i < c; i++) { bc.ProposeBirth(border, diagonal); bc.Accept(); bc.ProposeDeath(m); bc.Accept(); }
                }, 3.5 * mm, "3.5M^2"),
                new("rebuild", c => { for (int i = 0; i < c; i++) Sink = bc.Rebuild().GetHashCode(); }, mm * m / 3.0, "M^3/3"),
                new("solve", c => { for (int i = 0; i < c; i++) Sink = bc.Solve(rhs.AsSpan(0, m))[0]; }, mm, "M^2"),
            };
            CostProfile.Point("kernels:bordered", $"M={m}", parts);
        }
    }

    internal static void LowRankBlock()
    {
        foreach (int m in new[] { 20, 50, 200 })
            foreach (int k in new[] { 10, 39, 100 })
            {
                var rng = new Xoshiro256PlusPlus(600 + m + k);
                double[] diagonal = Enumerable.Range(0, m).Select(_ => 0.5 + rng.NextDouble()).ToArray();
                double[] weights = Enumerable.Range(0, k).Select(_ => 0.1 + rng.NextDouble()).ToArray();
                var factor = new double[m, k];
                for (int i = 0; i < m; i++) for (int j = 0; j < k; j++) factor[i, j] = 2.0 * rng.NextDouble() - 1.0;
                var w = new double[m, k];
                for (int i = 0; i < m; i++) for (int j = 0; j < k; j++) w[i, j] = factor[i, j] * Math.Sqrt(weights[j]) / Math.Sqrt(diagonal[i]);
                var sigma = new DiagonalPlusLowRank(diagonal, factor, weights);
                double[] rhs = Random(m, 19), x = new double[m];
                double svdSweepUnits = 0.5 * k * (k - 1) * (7.0 * m + 4.0 * k);
                var parts = new List<CostProfile.Part>
                {
                    new("construct", c => { for (int i = 0; i < c; i++) ObjectSink = new DiagonalPlusLowRank(diagonal, factor, weights); }, svdSweepUnits, "one-sided Jacobi sweep"),
                    new("thin_svd_alone", c => { for (int i = 0; i < c; i++) ObjectSink = ThinSvd.Decompose(w); }, svdSweepUnits, "one-sided Jacobi sweep"),
                    new("solve_in_place", c => { for (int i = 0; i < c; i++) { Array.Copy(rhs, x, m); sigma.SolveInPlace(x); } }, 2.0 * m * Math.Min(m, k), "madd"),
                    new("inverse_quadratic_form", c => { for (int i = 0; i < c; i++) Sink = sigma.InverseQuadraticForm(rhs); }, 4.0 * m * Math.Min(m, k), "madd"),
                };
                CostProfile.Point("kernels:lowrank", $"m={m} k={k}", parts);
            }
    }

    internal static void KhatriBlock()
    {
        const int count = 20;
        foreach (int nx in new[] { 100, 500 })
            foreach (int m in new[] { 10, 50, 200, 1000 })
            {
                var a = new double[count][];
                var h = new double[count][];
                for (int j = 0; j < count; j++) { a[j] = Random(nx, 700 + j); h[j] = Random(m, 800 + j); }
                double[] aNew = Random(nx, 901), hNew = Random(m, 902);
                var y = new double[nx, m];
                double[] flat = Random(nx * m, 903);
                for (int i = 0; i < nx; i++) for (int l = 0; l < m; l++) y[i, l] = flat[i * m + l];
                double borderUnits = (double)count * (nx + m) + nx + m;
                var parts = new List<CostProfile.Part>
                {
                    new("gram", c => { for (int i = 0; i < c; i++) ObjectSink = KhatriRao.Gram(a, h); }, 0.5 * count * count * (nx + m), "madd"),
                    new("border", c => { for (int i = 0; i < c; i++) ObjectSink = KhatriRao.Border(a, h, aNew, hNew, out double _); }, borderUnits, "madd"),
                    new("centred_border", c => { for (int i = 0; i < c; i++) ObjectSink = KhatriRao.CentredBorder(a, h, aNew, hNew, out double _); }, 3.0 * borderUnits, "pass element"),
                    new("cross_products", c => { for (int i = 0; i < c; i++) ObjectSink = KhatriRao.CrossProducts(a, h, y); }, (double)count * nx * m + count * nx, "madd"),
                    new("centred_cross_products", c => { for (int i = 0; i < c; i++) ObjectSink = KhatriRao.CentredCrossProducts(a, h, y); }, (double)count * nx * m + 6.0 * nx * m, "madd"),
                };
                CostProfile.Point("kernels:khatri", $"n_x={nx} m={m} M={count}", parts);
            }
    }
}
