using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The efficiency pass's levers 2, 3 and 4 (<c>efficiency-pass-brief.md</c>): the vectorized column kernels, the categorical
/// subsets' bitset and the commit's residual passes are element-wise and fused-multiply-add free, so each result is the old
/// scalar form's bit for bit. Each test holds a verbatim copy of the old form and compares bit patterns (every NaN taken as
/// equal to every NaN) over chosen cases and a seeded random parameter sweep, far ranges, ±0, ±∞ and NaN included. Run it
/// again with <c>DOTNET_EnableHWIntrinsic=0</c> for the unaccelerated path.
/// </summary>
public sealed class BassColumnBitsTests
{
    private readonly ITestOutputHelper _out;
    public BassColumnBitsTests(ITestOutputHelper output) => _out = output;

    // ---- The old scalar forms, verbatim ----

    private static double OldPower(int degree, double x) => degree switch { 1 => x, 2 => x * x, _ => Math.Pow(x, degree) };

    private static double OldValue(int degree, int sign, double knot, double x) => OldPower(degree, Math.Max(sign * (x - knot), 0.0));

    // HingeProductKernel.Apply as it was, the values route.
    private static void OldMultiply(int degree, IReadOnlyList<HingeFactor> factors, Func<int, double> knots, Func<int, double[]> values, double[] column)
    {
        if (factors.Count == 0) return;
        double peak = 1.0;
        foreach (HingeFactor f in factors)
        {
            double t = knots(f.Variable * 1000 + f.Knot), s = f.Sign;
            peak *= (s + 1.0) / 2.0 - s * t;
            double[] x = values(f.Variable);
            for (int i = 0; i < column.Length; i++) column[i] *= OldValue(degree, f.Sign, t, x[i]);
        }
        double c = OldPower(degree, peak);
        if (c == 0.0) return;
        for (int i = 0; i < column.Length; i++) column[i] /= c;
    }

    private static bool Same(double a, double b) => (double.IsNaN(a) && double.IsNaN(b)) || BitConverter.DoubleToInt64Bits(a) == BitConverter.DoubleToInt64Bits(b);

    private static void AssertSame(double[] expected, double[] actual, string what)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++)
            if (!Same(expected[i], actual[i]))
                Assert.Fail($"{what}: entry {i}: expected {expected[i]:R} ({BitConverter.DoubleToInt64Bits(expected[i]):X16}), got {actual[i]:R} ({BitConverter.DoubleToInt64Bits(actual[i]):X16})");
    }

    // A value drawn from the far ranges and the specials, or an ordinary one in [0, 1].
    private static double Draw(Xoshiro256PlusPlus rng, bool far)
    {
        if (!far || rng.NextDouble() < 0.4) return rng.NextDouble();
        int kind = (int)(rng.NextUInt64() % 12);
        double sign = rng.NextDouble() < 0.5 ? -1.0 : 1.0;
        return kind switch
        {
            0 => 0.0,
            1 => -0.0,
            2 => double.PositiveInfinity,
            3 => double.NegativeInfinity,
            4 => double.NaN,
            5 => sign * double.Epsilon * (1 + rng.NextUInt64() % 1000),
            6 => sign * double.MaxValue,
            _ => sign * Math.Pow(10.0, -300.0 + 600.0 * rng.NextDouble()),
        };
    }

    // ---- Lever 2: the hinge kernel ----

    private static (BassBasis Basis, double[][] Knots) Basis(int variables, int cases, int degree, Xoshiro256PlusPlus rng, bool ties)
    {
        var columns = new List<BassPredictor>();
        for (int v = 0; v < variables; v++)
        {
            double[] x = new double[cases];
            for (int i = 0; i < cases; i++) x[i] = ties ? Math.Floor(rng.NextDouble() * 7.0) : rng.NextDouble() * (1 + v);
            columns.Add(new ContinuousPredictor($"x{v}", x));
        }
        double[] y = [.. Enumerable.Range(0, cases).Select(_ => rng.NextDouble())];
        var data = new BassData([.. columns], new ScalarResponse(y));
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        return (new BassBasis(inputs, degree), [.. Enumerable.Range(0, variables).Select(v => inputs.ContinuousKnots(v).ToArray())]);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void HingeKernel_IsTheScalarForm_BitForBit(int degree)
    {
        var rng = new Xoshiro256PlusPlus(seed: 7000 + degree);
        int checkedTerms = 0;
        foreach (bool ties in new[] { false, true })
        {
            (BassBasis basis, double[][] knots) = Basis(5, 60, degree, rng, ties);
            foreach (int length in new[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 33, 63, 64, 67, 100, 1001 })
                foreach (bool far in new[] { false, true })
                {
                    double[][] x = [.. Enumerable.Range(0, 5).Select(_ => Enumerable.Range(0, length).Select(_ => Draw(rng, far)).ToArray())];
                    var at = new CaseInputs(x, []);
                    for (int trial = 0; trial < 40; trial++)
                    {
                        int order = 1 + (int)(rng.NextUInt64() % 4);
                        int[] variables = [.. Enumerable.Range(0, 5).OrderBy(_ => rng.NextUInt64()).Take(order)];
                        HingeFactor[] factors = [.. variables.Select(v => new HingeFactor(v, rng.NextDouble() < 0.5 ? 1 : -1, (int)(rng.NextUInt64() % (ulong)knots[v].Length)))];
                        var term = new BassTerm(factors);
                        var got = new double[length];
                        basis.CaseColumn(term, at, got);
                        var want = new double[length];
                        Array.Fill(want, 1.0);
                        OldMultiply(degree, term.Continuous, k => knots[k / 1000][k % 1000], v => x[v], want);
                        AssertSame(want, got, $"degree {degree}, length {length}, far {far}, {term}");
                        checkedTerms++;
                    }
                }
        }
        _out.WriteLine($"degree {degree}: {checkedTerms} terms bit-identical");
    }

    [Fact]
    public void HingeKernel_KnotsAtTheCasesValues_AndTheEnds()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7100);
        foreach (int degree in new[] { 1, 2 })
        {
            (BassBasis basis, double[][] knots) = Basis(3, 80, degree, rng, ties: false);
            // The training inputs themselves are the knots: x = t exactly, where the hinge is zero and max(−0, 0) matters.
            double[][] x = [.. knots.Select(k => k.Concat(k).ToArray())];
            var at = new CaseInputs(x, []);
            foreach (int v in new[] { 0, 1, 2 })
                foreach (int knot in new[] { 0, 1, knots[v].Length / 2, knots[v].Length - 2, knots[v].Length - 1 })
                    foreach (int sign in new[] { 1, -1 })
                    {
                        var term = new BassTerm([new HingeFactor(v, sign, knot)]);
                        var got = new double[at.Cases];
                        basis.CaseColumn(term, at, got);
                        var want = new double[at.Cases];
                        Array.Fill(want, 1.0);
                        OldMultiply(degree, term.Continuous, k => knots[k / 1000][k % 1000], vv => x[vv], want);
                        AssertSame(want, got, $"degree {degree}, {term}");
                    }
        }
    }

    // ---- Lever 3: the subset ----

    [Fact]
    public void Subset_Contains_IsTheBinarySearch()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7200);
        foreach (int levels in new[] { 2, 3, 63, 64, 65, 127, 128, 129, 1000, 1024, 5000 })
            for (int trial = 0; trial < 60; trial++)
            {
                int size = 1 + (int)(rng.NextUInt64() % (ulong)(levels - 1));
                int[] members = [.. Enumerable.Range(0, levels).OrderBy(_ => rng.NextUInt64()).Take(size)];
                var factor = new CategoricalFactor(trial % 3, members);
                int[] sorted = [.. members.OrderBy(l => l)];
                foreach (int level in Enumerable.Range(-3, levels + 130).Concat([int.MinValue, int.MaxValue, 1 << 20]))
                    Assert.Equal(level >= 0 && Array.BinarySearch(sorted, level) >= 0, factor.Contains(level));
            }
    }

    [Fact]
    public void Subset_ZeroOutside_IsTheCasePerCaseLoop()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7210);
        foreach (int levels in new[] { 2, 7, 64, 65, 300, 1024 })
            foreach (int length in new[] { 0, 1, 5, 64, 1000 })
            {
                int size = 1 + (int)(rng.NextUInt64() % (ulong)(levels - 1));
                var factor = new CategoricalFactor(0, Enumerable.Range(0, levels).OrderBy(_ => rng.NextUInt64()).Take(size));
                int[] codes = [.. Enumerable.Range(0, length).Select(_ => (int)(rng.NextUInt64() % (ulong)(levels + 3)) - 1)];   // −1 and past the table included
                double[] start = [.. Enumerable.Range(0, length).Select(_ => rng.NextDouble() - 0.5)];
                double[] want = (double[])start.Clone(), got = (double[])start.Clone();
                for (int i = 0; i < length; i++)
                    if (!factor.Contains(codes[i])) want[i] = 0.0;
                factor.ZeroOutside(codes, got);
                AssertSame(want, got, $"levels {levels}, length {length}");
            }
    }

    // ---- Lever 2 and 4: the column passes ----

    [Fact]
    public void ColumnOps_AreTheScalarLoops_BitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7300);
        foreach (int n in new[] { 0, 1, 3, 4, 5, 7, 8, 9, 16, 17, 31, 63, 64, 65, 129, 1000, 4099 })
            foreach (bool far in new[] { false, true })
            {
                double[] z = [.. Enumerable.Range(0, n).Select(_ => Draw(rng, far))];
                double[] w = [.. Enumerable.Range(0, n).Select(_ => Draw(rng, far))];
                double amount = Draw(rng, far), factor = Draw(rng, far);

                // NonzeroCount
                int want = 0;
                foreach (double v in z) if (v != 0.0) want++;
                Assert.Equal(want, ColumnOps.NonzeroCount(z));

                // Subtract
                double[] a = (double[])z.Clone(), b = (double[])z.Clone();
                for (int i = 0; i < n; i++) a[i] -= amount;
                ColumnOps.Subtract(b, amount);
                AssertSame(a, b, $"Subtract n={n}");

                // CentreAndScale
                a = (double[])z.Clone(); b = (double[])z.Clone();
                for (int i = 0; i < n; i++) a[i] = w[i] * (a[i] - amount);
                ColumnOps.CentreAndScale(b, amount, w);
                AssertSame(a, b, $"CentreAndScale n={n}");

                // CentreInto, with and without roots
                double[] into = new double[n], wantInto = new double[n];
                for (int i = 0; i < n; i++) wantInto[i] = z[i] - amount;
                ColumnOps.CentreInto(into, z, amount, default);
                AssertSame(wantInto, into, $"CentreInto n={n}");
                for (int i = 0; i < n; i++) wantInto[i] *= w[i];
                ColumnOps.CentreInto(into, z, amount, w);
                AssertSame(wantInto, into, $"CentreInto with roots n={n}");

                // Scale
                for (int i = 0; i < n; i++) wantInto[i] = factor * z[i];
                ColumnOps.Scale(into, factor, z);
                AssertSame(wantInto, into, $"Scale n={n}");

                // SubtractScaledProduct
                double level = Draw(rng, far), beta = Draw(rng, far), ai = Draw(rng, far);
                a = (double[])z.Clone(); b = (double[])z.Clone();
                for (int l = 0; l < n; l++) a[l] -= beta * (ai * w[l] - level);
                ColumnOps.SubtractScaledProduct(b, beta, ai, w, level);
                AssertSame(a, b, $"SubtractScaledProduct n={n}");
            }
    }

    [Fact]
    public void ColumnOps_SubtractScaledColumns_IsColumnByColumn()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7400);
        foreach (int n in new[] { 0, 1, 3, 4, 5, 8, 17, 100, 1000 })
            foreach (int count in new[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 63 })
                foreach (bool far in new[] { false, true })
                {
                    double[][] columns = [.. Enumerable.Range(0, count).Select(_ => Enumerable.Range(0, n).Select(_ => Draw(rng, far)).ToArray())];
                    double[] betas = [.. Enumerable.Range(0, count).Select(_ => Draw(rng, far))];
                    double[] start = [.. Enumerable.Range(0, n).Select(_ => Draw(rng, far))];
                    double[] want = (double[])start.Clone(), got = (double[])start.Clone();
                    for (int j = 0; j < count; j++)
                        for (int i = 0; i < n; i++) want[i] -= betas[j] * columns[j][i];
                    ColumnOps.SubtractScaledColumns(got, columns, betas);
                    AssertSame(want, got, $"n={n}, {count} columns, far {far}");
                }
    }
}
