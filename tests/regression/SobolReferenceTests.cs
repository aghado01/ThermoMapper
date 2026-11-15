using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The sampled-sensitivity brief's slice 0 gates on its references (<see cref="SobolReferences"/>): each closed form against tensor
/// Gauss–Legendre quadrature (<see cref="SobolQuadrature"/>), a method apart from it and from the estimators under test. The mean and
/// every closed variance τ_u to 10⁻¹² of the variance; conditional means at seeded points to 10⁻¹² of the function's scale; the main and
/// total variances a family states directly against their Möbius sums; and the series behind Friedman's against quadrature.
/// </summary>
public sealed class SobolReferenceTests(ITestOutputHelper output)
{
    private const double Tolerance = 1e-12;

    public static TheoryData<string> Names()
    {
        var names = new TheoryData<string>();
        foreach (SobolReference r in SobolReferences.All()) names.Add(r.Name);
        return names;
    }

    private static SobolReference Named(string name) => SobolReferences.All().Single(r => r.Name == name);

    // Every nonempty subset of the inputs f varies with.
    private static IEnumerable<int[]> Subsets(SobolReference f)
    {
        int[] varying = [.. Enumerable.Range(0, f.Dimension).Where(f.Varies)];
        for (int mask = 1; mask < 1 << varying.Length; mask++)
            yield return [.. Enumerable.Range(0, varying.Length).Where(k => (mask >> k & 1) != 0).Select(k => varying[k])];
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void ClosedVariances_MatchQuadrature(string name)
    {
        SobolReference f = Named(name);
        IReadOnlyList<QuadratureRule> rules = SobolQuadrature.Rules(f);
        double v = f.Variance, worst = Math.Abs(f.Mean - SobolQuadrature.Mean(f, rules)) / Math.Max(1.0, Math.Abs(f.Mean));
        Assert.True(v > 0.0);
        foreach (int[] u in Subsets(f))
            worst = Math.Max(worst, Math.Abs(f.ClosedVariance(u) - SobolQuadrature.ClosedVariance(f, rules, u)) / v);
        output.WriteLine($"{name}: V = {v:R}, worst |τ − quadrature| / V = {worst:E2}");
        Assert.True(worst <= Tolerance, $"{name}: {worst:E2}");
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void ConditionalMeans_MatchQuadrature(string name)
    {
        SobolReference f = Named(name);
        IReadOnlyList<QuadratureRule> rules = SobolQuadrature.Rules(f);
        double scale = Math.Abs(f.Mean) + Math.Sqrt(f.Variance), worst = 0.0;
        var rng = new Xoshiro256PlusPlus(20261115);
        foreach (int[] u in Subsets(f).Where(u => u.Length <= 2).Append([.. Enumerable.Range(0, f.Dimension)]))
            for (int point = 0; point < 6; point++)
            {
                double[] at = [.. u.Select(v => f.Law.Lower[v] + (f.Law.Upper[v] - f.Law.Lower[v]) * rng.NextDouble())];
                worst = Math.Max(worst, Math.Abs(f.ConditionalMean(u, at) - SobolQuadrature.ConditionalMean(f, rules, u, at)) / scale);
            }
        output.WriteLine($"{name}: worst |E[f | x_u] − quadrature| / scale = {worst:E2}");
        Assert.True(worst <= Tolerance, $"{name}: {worst:E2}");
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void StatedPieces_MatchTheirMobiusSums(string name)
    {
        SobolReference f = Named(name);
        double v = f.Variance, worst = 0.0, mains = 0.0;
        int[] all = [.. Enumerable.Range(0, f.Dimension)];
        foreach (int[] u in Subsets(f))
        {
            double mobius = 0.0;
            for (int mask = 0; mask < 1 << u.Length; mask++)
            {
                int[] w = [.. Enumerable.Range(0, u.Length).Where(k => (mask >> k & 1) != 0).Select(k => u[k])];
                if (w.Length > 0) mobius += ((u.Length - w.Length) % 2 == 0 ? 1.0 : -1.0) * f.ClosedVariance(w);
            }
            double main = f.MainVariance(u);
            mains += main;
            worst = Math.Max(worst, Math.Abs(main - mobius) / v);
            Assert.True(main >= -Tolerance * v, $"{name}: V_[{string.Join(",", u)}] = {main:R}");
        }
        for (int i = 0; i < f.Dimension; i++)
        {
            double rest = f.Dimension == 1 ? 0.0 : f.ClosedVariance([.. all.Where(j => j != i)]);
            worst = Math.Max(worst, Math.Abs(f.TotalVariance(i) - (v - rest)) / v);
        }
        output.WriteLine($"{name}: worst |stated − Möbius| / V = {worst:E2}; Σ mains / V − 1 = {mains / v - 1.0:E2}");
        Assert.True(worst <= Tolerance, $"{name}: {worst:E2}");
        Assert.True(Math.Abs(mains / v - 1.0) <= Tolerance, $"{name}: the mains sum to {mains / v:R} of V");
    }

    [Fact]
    public void Effects_AreTheAnovaComponents()
    {
        // f₀ + Σ_u f_u(x_u) over every subset rebuilds f at seeded points, and each main effect averages to 0 over its input.
        var rng = new Xoshiro256PlusPlus(15);
        foreach (SobolReference f in SobolReferences.All())
        {
            IReadOnlyList<QuadratureRule> rules = SobolQuadrature.Rules(f);
            double scale = Math.Abs(f.Mean) + Math.Sqrt(f.Variance), worst = 0.0;
            for (int point = 0; point < 8; point++)
            {
                double[] x = [.. Enumerable.Range(0, f.Dimension).Select(v => f.Law.Lower[v] + (f.Law.Upper[v] - f.Law.Lower[v]) * rng.NextDouble())];
                double sum = f.Mean;
                foreach (int[] u in Subsets(f)) sum += f.Effect(u, [.. u.Select(v => x[v])]);
                worst = Math.Max(worst, Math.Abs(sum - f.Value(x)) / scale);
            }
            for (int v = 0; v < f.Dimension; v++)
            {
                if (!f.Varies(v)) continue;
                QuadratureRule rule = rules[v];
                double mean = 0.0;
                for (int q = 0; q < rule.Nodes.Length; q++) mean += rule.Weights[q] * f.Effect([v], [rule.Nodes[q]]);
                worst = Math.Max(worst, Math.Abs(mean) / scale);
            }
            output.WriteLine($"{f.Name}: worst ANOVA identity error / scale = {worst:E2}");
            Assert.True(worst <= Tolerance, $"{f.Name}: {worst:E2}");
        }
    }

    [Fact]
    public void FriedmansSeries_MatchQuadrature()
    {
        QuadratureRule Rule(double upper) => QuadratureRule.Uniform(0.0, upper, 24, 8);
        double Integral(double upper, Func<double, double> g)
        {
            QuadratureRule rule = Rule(upper);
            double s = 0.0;
            for (int q = 0; q < rule.Nodes.Length; q++) s += rule.Weights[q] * g(rule.Nodes[q]);
            return s * upper;
        }
        double si = Integral(Math.PI, t => Math.Sin(t) / t), si2 = Integral(2.0 * Math.PI, t => Math.Sin(t) / t);
        double cin = Integral(Math.PI, t => 2.0 * Math.Pow(Math.Sin(0.5 * t), 2) / t);
        const decimal Pi = 3.1415926535897932384626433833m;
        double[] errors =
        [
            Math.Abs(FriedmanSobolReference.SineIntegral(Pi) - si),
            Math.Abs(FriedmanSobolReference.SineIntegral(2m * Pi) - si2),
            Math.Abs(FriedmanSobolReference.CosineIntegral(Pi) - cin),
            Math.Abs(FriedmanSobolReference.CosineIntegral(2m * Pi) - FriedmanReference.CinTwoPi()),
        ];
        output.WriteLine($"Si(π) = {si:R}, Si(2π) = {si2:R}, Cin(π) = {cin:R}; errors {string.Join(", ", errors.Select(e => e.ToString("E2")))}");
        Assert.All(errors, e => Assert.True(e <= 4e-15, $"{e:E2}"));
    }

    [Fact]
    public void Ishigami_VarianceIsItsKnownSum()
    {
        // V = a²/8 + bπ⁴/5 + b²π⁸/18 + ½ (cited from memory) for the chosen and a few other parameters.
        foreach ((double a, double b) in new[] { (7.0, 0.1), (0.0, 0.0), (5.0, 0.0), (13.0, 0.4) })
        {
            var f = new IshigamiReference(a, b);
            double known = a * a / 8.0 + b * Math.Pow(Math.PI, 4) / 5.0 + b * b * Math.Pow(Math.PI, 8) / 18.0 + 0.5;
            Assert.True(Math.Abs(f.Variance - known) <= 1e-14 * known, $"a = {a}, b = {b}: {f.Variance:R} against {known:R}");
        }
    }
}
