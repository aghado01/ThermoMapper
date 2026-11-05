using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The efficiency pass's primitives, old form against new in one process (<c>efficiency-pass-brief.md</c>, levers 2 to 4): the
/// old scalar forms are the verbatim copies <see cref="BassColumnBitsTests"/> holds, and each pair runs round-robin with the other
/// (a warm-up round, then nine batches of calls, the median), so load falls on both alike. The pairs: a term's column
/// (evaluation, support count, centring), a subset's membership over the cases, and the commit's explicit residual.
/// Behind <c>BASS_QUALIFICATION</c> = 1; rows go to <c>BASS_COST_OUT</c>; times wait for a quiet machine.
/// </summary>
public sealed class BassPrimitiveCostTests
{
    private readonly ITestOutputHelper _out;
    public BassPrimitiveCostTests(ITestOutputHelper output) => _out = output;

    private static Func<int, double> Timer(int calls, Func<double> call)
    {
        double sink = 0.0;
        return [MethodImpl(MethodImplOptions.AggressiveOptimization)] (int _) =>
        {
            long t0 = Stopwatch.GetTimestamp();
            for (int c = 0; c < calls; c++) sink += call();
            return (Stopwatch.GetTimestamp() - t0) * 1e9 / Stopwatch.Frequency / calls + sink * 0.0;
        };
    }

    private static int Calls(int unitsPerCall, long target = 4_000_000) => (int)Math.Clamp(target / Math.Max(1, unitsPerCall), 2, 200_000);

    [Fact]
    public void Column_OldAgainstNew()
    {
        if (!BassQualification.Enabled) return;
        var rows = new BassCost.Rows("primitives", "column");
        foreach ((string q, string v) in BassCost.Machine()) rows.Add("machine", q, v);
        var load = new BassCost.Load();
        var rng = new Xoshiro256PlusPlus(seed: 811);
        foreach (int n in new[] { 200, 1000, 10000, 50000 })
            foreach (int degree in new[] { 1, 2 })
                foreach (int order in new[] { 1, 2 })
                {
                    var columns = new List<BassPredictor>();
                    for (int v = 0; v < 4; v++) columns.Add(new ContinuousPredictor($"x{v}", [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())]));
                    var data = new BassData([.. columns], new ScalarResponse([.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())]));
                    BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
                    var basis = new BassBasis(inputs, degree);
                    HingeFactor[] factors = [.. Enumerable.Range(0, order).Select(v => new HingeFactor(v, v % 2 == 0 ? 1 : -1, inputs.ContinuousKnots(v).Length / (2 + v)))];
                    var term = new BassTerm(factors);
                    double[][] x = [.. Enumerable.Range(0, 4).Select(v => inputs.Cases.Continuous(v).ToArray())];
                    double[][] knots = [.. Enumerable.Range(0, 4).Select(v => inputs.ContinuousKnots(v).ToArray())];
                    var column = new double[n];
                    double minimum = 5.0;

                    // As BassDesign.Evaluate formed a scalar column: fill, kernels, support count, mean, centring.
                    double Old()
                    {
                        Array.Fill(column, 1.0);
                        OldHinge(degree, term.Continuous, knots, x, column);
                        int nonzero = 0;
                        foreach (double v in column) if (v != 0.0) nonzero++;
                        if (nonzero < minimum) return 0.0;
                        double mean = PairwiseSum.Of(column) / n;
                        for (int i = 0; i < column.Length; i++) column[i] -= mean;
                        return column[n / 2];
                    }
                    double New()
                    {
                        basis.CaseColumn(term, inputs.Cases, column);
                        if (ColumnOps.NonzeroCount(column) < minimum) return 0.0;
                        ColumnOps.Subtract(column, PairwiseSum.Of(column) / n);
                        return column[n / 2];
                    }
                    double[] before = (double[])ColumnOf(Old, column).Clone();
                    double[] after = (double[])ColumnOf(New, column).Clone();
                    for (int i = 0; i < n; i++) Assert.Equal(BitConverter.DoubleToInt64Bits(before[i]), BitConverter.DoubleToInt64Bits(after[i]));

                    int calls = Calls(n * (order + 4));
                    Dictionary<string, double[]> times = BassCost.Interleave([("old", Timer(calls, Old)), ("new", Timer(calls, New))]);
                    double old = BassCost.Median(times["old"]), now = BassCost.Median(times["new"]);
                    string key = $"n={n} q={degree} order={order}";
                    rows.Add(key, "old ns per case", old / n);
                    rows.Add(key, "new ns per case", now / n);
                    rows.Add(key, "old ns per element pass (n(order + 4))", old / (n * (order + 4)));
                    rows.Add(key, "new ns per element pass (n(order + 4))", now / (n * (order + 4)));
                    rows.Add(key, "speedup", old / now);
                    _out.WriteLine($"column {key}: old {old / n:F2} ns/case, new {now / n:F2}, ×{old / now:F1}");
                }
        Finish(rows, load);
    }

    private static double[] ColumnOf(Func<double> call, double[] column)
    {
        call();
        return column;
    }

    // HingeProductKernel.Apply's values route as it was.
    private static void OldHinge(int degree, IReadOnlyList<HingeFactor> factors, double[][] knots, double[][] x, double[] column)
    {
        double peak = 1.0;
        foreach (HingeFactor f in factors)
        {
            double t = knots[f.Variable][f.Knot], s = f.Sign;
            peak *= (s + 1.0) / 2.0 - s * t;
            double[] xs = x[f.Variable];
            for (int i = 0; i < column.Length; i++) column[i] *= Value(degree, f.Sign, t, xs[i]);
        }
        double c = degree switch { 1 => peak, 2 => peak * peak, _ => Math.Pow(peak, degree) };
        if (c == 0.0) return;
        for (int i = 0; i < column.Length; i++) column[i] /= c;
    }

    private static double Value(int degree, int sign, double knot, double x)
    {
        double u = Math.Max(sign * (x - knot), 0.0);
        return degree switch { 1 => u, 2 => u * u, _ => Math.Pow(u, degree) };
    }

    [Fact]
    public void Membership_OldAgainstNew()
    {
        if (!BassQualification.Enabled) return;
        var rows = new BassCost.Rows("primitives", "membership");
        foreach ((string q, string v) in BassCost.Machine()) rows.Add("machine", q, v);
        var load = new BassCost.Load();
        var rng = new Xoshiro256PlusPlus(seed: 821);
        const int n = 4096;
        foreach (int levels in new[] { 2, 4, 16, 64, 256, 1024 })
        {
            var factor = new CategoricalFactor(0, Enumerable.Range(0, levels).Where(l => l % 2 == 0));   // a subset of L/2 levels, as the cost half read it
            int[] sorted = [.. factor.Levels];
            int[] codes = [.. Enumerable.Range(0, n).Select(i => i % levels)];
            var column = new double[n];
            double Old()
            {
                Array.Fill(column, 1.0);
                foreach (int c in codes) { }
                for (int i = 0; i < n; i++)
                    if (!(codes[i] >= 0 && Array.BinarySearch(sorted, codes[i]) >= 0)) column[i] = 0.0;
                return column[n / 3];
            }
            double New()
            {
                Array.Fill(column, 1.0);
                factor.ZeroOutside(codes, column);
                return column[n / 3];
            }
            int calls = Calls(n);
            Dictionary<string, double[]> times = BassCost.Interleave([("old", Timer(calls, Old)), ("new", Timer(calls, New))]);
            double old = BassCost.Median(times["old"]), now = BassCost.Median(times["new"]);
            string key = $"L={levels} n={n}";
            rows.Add(key, "old ns per case", old / n);
            rows.Add(key, "new ns per case", now / n);
            rows.Add(key, "speedup", old / now);
            _out.WriteLine($"membership {key}: old {old / n:F2} ns/case, new {now / n:F2}, ×{old / now:F1}");
        }
        Finish(rows, load);
    }

    [Fact]
    public void Refit_OldAgainstNew()
    {
        if (!BassQualification.Enabled) return;
        var rows = new BassCost.Rows("primitives", "refit");
        foreach ((string q, string v) in BassCost.Machine()) rows.Add("machine", q, v);
        var load = new BassCost.Load();
        var rng = new Xoshiro256PlusPlus(seed: 831);
        foreach ((int n, int m) in new[] { (1000, 20), (10000, 5), (10000, 200), (50000, 20) })
        {
            double[][] columns = [.. Enumerable.Range(0, m).Select(_ => Enumerable.Range(0, n).Select(_ => rng.NextDouble() - 0.5).ToArray())];
            double[] beta = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble())];
            double[] y = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            var r = new double[n];
            double Old()
            {
                Array.Copy(y, r, n);
                for (int j = 0; j < m; j++)
                {
                    double bj = beta[j];
                    double[] z = columns[j];
                    for (int i = 0; i < n; i++) r[i] -= bj * z[i];
                }
                return r[n / 2];
            }
            double New()
            {
                Array.Copy(y, r, n);
                ColumnOps.SubtractScaledColumns(r, columns, beta);
                return r[n / 2];
            }
            Old();
            double[] before = (double[])r.Clone();
            New();
            for (int i = 0; i < n; i++) Assert.Equal(BitConverter.DoubleToInt64Bits(before[i]), BitConverter.DoubleToInt64Bits(r[i]));
            int calls = Calls(n * m, 8_000_000);
            Dictionary<string, double[]> times = BassCost.Interleave([("old", Timer(calls, Old)), ("new", Timer(calls, New))]);
            double old = BassCost.Median(times["old"]), now = BassCost.Median(times["new"]);
            string key = $"scalar N={n} M={m}";
            rows.Add(key, "old ns per N·M", old / ((double)n * m));
            rows.Add(key, "new ns per N·M", now / ((double)n * m));
            rows.Add(key, "speedup", old / now);
            _out.WriteLine($"refit {key}: old {old / ((double)n * m):F2} ns/NM, new {now / ((double)n * m):F2}, ×{old / now:F1}");
        }
        foreach ((int cases, int locations, int m) in new[] { (100, 10, 20), (100, 200, 20), (500, 1000, 20) })
        {
            double[][] a = [.. Enumerable.Range(0, m).Select(_ => Enumerable.Range(0, cases).Select(_ => rng.NextDouble()).ToArray())];
            double[][] h = [.. Enumerable.Range(0, m).Select(_ => Enumerable.Range(0, locations).Select(_ => rng.NextDouble()).ToArray())];
            double[] beta = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble())];
            double[] y = [.. Enumerable.Range(0, cases * locations).Select(_ => rng.NextDouble())];
            var r = new double[cases * locations];
            double[] level = [.. a.Zip(h, (aj, hj) => PairwiseSum.Of(aj) / cases * (PairwiseSum.Of(hj) / locations))];
            double Old()
            {
                Array.Copy(y, r, r.Length);
                for (int j = 0; j < m; j++)
                {
                    double bj = beta[j];
                    double[] aj = a[j], hj = h[j];
                    for (int i = 0; i < cases; i++)
                    {
                        double ai = aj[i];
                        int row = i * locations;
                        for (int l = 0; l < locations; l++) r[row + l] -= bj * (ai * hj[l] - level[j]);
                    }
                }
                return r[r.Length / 2];
            }
            double New()
            {
                Array.Copy(y, r, r.Length);
                for (int j = 0; j < m; j++)
                    for (int i = 0; i < cases; i++) ColumnOps.SubtractScaledProduct(r.AsSpan(i * locations, locations), beta[j], a[j][i], h[j], level[j]);
                return r[r.Length / 2];
            }
            Old();
            double[] before = (double[])r.Clone();
            New();
            for (int i = 0; i < r.Length; i++) Assert.Equal(BitConverter.DoubleToInt64Bits(before[i]), BitConverter.DoubleToInt64Bits(r[i]));
            int calls = Calls(cases * locations * m, 8_000_000);
            Dictionary<string, double[]> times = BassCost.Interleave([("old", Timer(calls, Old)), ("new", Timer(calls, New))]);
            double old = BassCost.Median(times["old"]), now = BassCost.Median(times["new"]);
            double units = (double)cases * locations * m;
            string key = $"direct n_x={cases} m={locations} M={m}";
            rows.Add(key, "old ns per N·M", old / units);
            rows.Add(key, "new ns per N·M", now / units);
            rows.Add(key, "speedup", old / now);
            _out.WriteLine($"refit {key}: old {old / units:F2} ns/NM, new {now / units:F2}, ×{old / now:F1}");
        }
        Finish(rows, load);
    }

    private void Finish(BassCost.Rows rows, BassCost.Load load)
    {
        (string[] end, double cpu, double seconds, double others) = load.End();
        rows.Add("load", "other test hosts at the start", string.Join(" ", load.AtStart));
        rows.Add("load", "other test hosts at the end", string.Join(" ", end));
        rows.Add("load", "own CPU over wall", cpu);
        rows.Add("load", "other processes' cores busy", others);
        rows.Add("load", "block seconds", seconds);
        rows.Add("load", "flagged", load.AtStart.Length > 0 || end.Length > 0 ? "yes" : "no");
        _out.WriteLine($"load: others {others:F2} cores, {(load.AtStart.Length > 0 || end.Length > 0 ? "FLAGGED" : "quiet")}; rows {rows.Flush()}");
    }
}
