using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// E8's oracle and its defects (<c>active_subspace_oracle.R</c>; the active-subspace brief, slice 0): coactivity's <c>C_bass</c>
/// (knrumsey/coactivity@35e4573, the closed form of Rumsey et al., arXiv 2307.11241) on hinge draws given here, beside the mode's
/// reference, which integrates each draw's own gradient by tensor Gauss–Legendre split at its knots under the stated law. Each defect
/// is confirmed or struck, with its size computed on this side where it can be. These are corrections, not parity variants: the engine
/// takes the stated law.
/// <list type="bullet">
/// <item><b>Parity:</b> on degree-one draws under uniform laws, <c>C_bass</c> is the reference.</item>
/// <item><b>Confirmed:</b> <c>C_bass</c> reads no hinge degree; the truncated normal's second moment flips the lower end's term, in the
/// code and in the paper's ξ_TNorm; a normal mixture's components are normalized twice; <c>build_prior</c>'s gamma prior reaches
/// <c>C_bass</c> as NA, and a hand-set one is not renormalized on the box; additive fits (maxInt = 1) stop, as the paper's supplement
/// says.</item>
/// <item><b>Readings:</b> <c>C_bass</c> without <c>mcmc.use</c> reads the last draw only; <c>use_native_scale</c> defaults to false where
/// the documentation says true.</item>
/// </list>
/// </summary>
public sealed class ActiveSubspaceDefectParityTests(ITestOutputHelper output)
{
    // f = 0.5 + 1.2·h(x₁; +, 0.3) − 0.8·h(x₂; −, 0.6) + 2·h(x₁; +, 0.2)·h(x₂; +, 0.4), BASS's normalized hinges.
    private static readonly (int Input, int Sign, double Knot)[][] TwoInputTerms =
    [
        [(1, 1, 0.3)], [(2, -1, 0.6)], [(1, 1, 0.2), (2, 1, 0.4)],
    ];

    private static readonly double[] TwoInputCoefficients = [0.5, 1.2, -0.8, 2.0];

    private static readonly (int Input, int Sign, double Knot)[][] ThreeInputTerms =
    [
        [(1, 1, 0.3)], [(2, -1, 0.7), (3, 1, 0.4)], [(1, -1, 0.6), (2, 1, 0.25), (3, -1, 0.8)], [(3, 1, 0.55)],
    ];

    private static readonly double[] ThreeInputCoefficients = [-0.2, 1.5, 0.9, -2.5, 0.7];

    // A law per input: uniform (lower, upper), or normal components (mean, sd), equal weights, truncated to [0, 1].
    private abstract record Law;
    private sealed record UniformLaw(double Lower, double Upper) : Law;
    private sealed record NormalLaw(params (double Mean, double Sd)[] Components) : Law;

    private sealed record Draw(int Id, int Degree, (int Input, int Sign, double Knot)[][] Terms, double[] Coefficients, Law[] Laws);

    private static readonly Draw[] Draws =
    [
        new(1, 1, TwoInputTerms, TwoInputCoefficients, [new UniformLaw(0, 1), new UniformLaw(0, 1)]),
        new(2, 2, TwoInputTerms, TwoInputCoefficients, [new UniformLaw(0, 1), new UniformLaw(0, 1)]),
        new(3, 1, TwoInputTerms, TwoInputCoefficients, [new NormalLaw((0.5, 0.5)), new NormalLaw((0.2, 0.3))]),
        new(4, 1, TwoInputTerms, TwoInputCoefficients, [new NormalLaw((0.2, 0.1), (0.9, 0.4)), new UniformLaw(0.1, 0.8)]),
        new(5, 1, ThreeInputTerms, ThreeInputCoefficients, [new UniformLaw(0, 1), new UniformLaw(0, 1), new UniformLaw(0.25, 0.75)]),
    ];

    private static IEnumerable<double[]> Rows(Draw d)
    {
        yield return [3, d.Id, 0, d.Laws.Length, d.Degree, 0];
        for (int m = 0; m < d.Coefficients.Length; m++) yield return [1, d.Id, m, d.Coefficients[m], 0, 0];
        for (int m = 0; m < d.Terms.Length; m++)
            foreach ((int input, int sign, double knot) in d.Terms[m]) yield return [2, d.Id, m + 1, input, sign, knot];
        for (int v = 0; v < d.Laws.Length; v++)
            switch (d.Laws[v])
            {
                case UniformLaw u:
                    yield return [4, d.Id, v + 1, 1, u.Lower, u.Upper];
                    break;
                case NormalLaw n:
                    foreach ((double mean, double sd) in n.Components) yield return [4, d.Id, v + 1, 2, mean, sd];
                    break;
            }
    }

    private static readonly Lazy<JsonElement?> Fixture = new(() => GbassOracle.Run("active_subspace_oracle.R", Draws.SelectMany(Rows)));

    private static JsonElement DrawResult(JsonElement r, int id) => r.GetProperty("draws").EnumerateArray().Single(d => d.GetProperty("draw").GetInt32() == id);

    private static double[,] Matrix(JsonElement e)
    {
        double[][] rows = [.. e.EnumerateArray().Select(GbassOracle.Numbers)];
        var m = new double[rows.Length, rows[0].Length];
        for (int i = 0; i < rows.Length; i++)
            for (int j = 0; j < rows[i].Length; j++) m[i, j] = rows[i][j];
        return m;
    }

    private static double Gap(double[,] a, double[,] b)
    {
        double worst = 0.0, scale = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
            {
                worst = Math.Max(worst, Math.Abs(a[i, j] - b[i, j]));
                scale = Math.Max(scale, Math.Abs(b[i, j]));
            }
        return worst / scale;
    }

    // A normal component's mass on [0, 1].
    private static double Mass(double mean, double sd) => Normal.MeanSd(mean, sd).Cdf(1.0) - Normal.MeanSd(mean, sd).Cdf(0.0);

    [Fact]
    public void CBass_OnDegreeOneDrawsUnderUniformLaws_IsTheReference()
    {
        if (Fixture.Value is not { } r) return;
        foreach (int id in new[] { 1, 5 })
        {
            JsonElement d = DrawResult(r, id);
            double gap = Gap(Matrix(d.GetProperty("C")), Matrix(d.GetProperty("reference")));
            output.WriteLine($"draw {id} ({d.GetProperty("p").GetInt32()} inputs): C_bass against the reference {gap:E2}");
            Assert.True(gap <= 1e-12, $"draw {id}: {gap:E2}");
        }
    }

    /// <summary>C_bass's integrals and hinge scale are degree one's, and it reads no degree: a degree-two draw gets the same C as
    /// the same knots at degree one.</summary>
    [Fact]
    public void CBass_ReadsNoDegree()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement one = DrawResult(r, 1), two = DrawResult(r, 2);
        double same = Gap(Matrix(two.GetProperty("C")), Matrix(one.GetProperty("C")));
        double missed = Gap(Matrix(two.GetProperty("C")), Matrix(two.GetProperty("reference")));
        double patched = Gap(Matrix(two.GetProperty("corrected")), Matrix(two.GetProperty("reference")));
        output.WriteLine($"degree two against degree one in C_bass {same:E2}; against degree two's reference {missed:E3} (patched normal terms change nothing: {patched:E3})");
        Assert.Equal(0.0, same);
        Assert.True(missed > 0.05, $"{missed:E2}");
    }

    /// <summary>
    /// ∫_a^b z²φ(z)dz = Δ + Aφ(A) − Bφ(B) for the standardized ends A, B and Δ = Φ(B) − Φ(A). XI_FUNC_TNORM_SVW and the paper's
    /// ξ_TNorm write Δ − Aφ(A) − Bφ(B), wrong wherever the lower end is finite; patched, the code meets the integral.
    /// </summary>
    [Fact]
    public void NormalSecondMoment_FlipsTheLowerEndsTerm()
    {
        if (Fixture.Value is not { } r) return;
        foreach (JsonElement c in r.GetProperty("normal_second_moment").EnumerateArray())
        {
            double mu = c.GetProperty("mean").GetDouble(), sigma = c.GetProperty("sd").GetDouble();
            double lo = c.GetProperty("lower").GetDouble(), hi = c.GetProperty("upper").GetDouble();   // clipped to [0, 1] by the mode
            Normal standard = Normal.Standard;
            double a = (lo - mu) / sigma, b = (hi - mu) / sigma, delta = standard.Cdf(b) - standard.Cdf(a);
            double pa = Math.Exp(standard.LogPdf(a)), pb = Math.Exp(standard.LogPdf(b));
            double here = (mu * mu * delta + 2 * sigma * mu * (pa - pb) + sigma * sigma * (delta + a * pa - b * pb)) / Mass(mu, sigma);
            double code = c.GetProperty("code").GetDouble(), corrected = c.GetProperty("corrected").GetDouble(), integral = c.GetProperty("integral").GetDouble();
            output.WriteLine($"N({mu}, {sigma}²) on [{lo}, {hi}]: E[x²; a, b] here {here:F12}, integrate {integral:F12}, patched {corrected:F12}, as written {code:F12}");
            Assert.Equal(1.0, integral / here, 12);
            Assert.Equal(1.0, corrected / here, 12);
            Assert.True(Math.Abs(code / here - 1.0) > 1e-3, $"{code} against {here}");
        }
    }

    /// <summary>
    /// C_bass weights each normal component w/cc, cc = Σ w_l z_l, then takes the component's truncated moment, which divides by z_l
    /// again: the law it reads is Σ (w_l/cc)φ_l/z_l on [0, 1], of mass Σ w_l/(cc z_l). With the second moment patched, C_bass meets
    /// that law (as_read); with the weights w z/cc as well, the stated law. One component per input scales C by Π 1/z_k.
    /// </summary>
    [Fact]
    public void NormalMixture_ComponentsAreNormalizedTwice()
    {
        if (Fixture.Value is not { } r) return;
        foreach (int id in new[] { 3, 4 })
        {
            JsonElement d = DrawResult(r, id);
            double[,] reference = Matrix(d.GetProperty("reference"));
            double asWritten = Gap(Matrix(d.GetProperty("C")), reference);
            double secondFixed = Gap(Matrix(d.GetProperty("second_moment_fixed")), Matrix(d.GetProperty("as_read")));
            double bothFixed = Gap(Matrix(d.GetProperty("corrected")), reference);
            output.WriteLine($"draw {id}: as written against the stated law {asWritten:E3}; second moment patched against the law it reads {secondFixed:E2}; both patched against the stated law {bothFixed:E2}");
            Assert.True(asWritten > 0.1, $"{asWritten:E2}");
            Assert.True(secondFixed <= 1e-12, $"{secondFixed:E2}");
            Assert.True(bothFixed <= 1e-12, $"{bothFixed:E2}");
        }

        // One component per input: the law read is the stated one times 1/z per input, so C is the reference times Π 1/z_k.
        JsonElement single = DrawResult(r, 3);
        double factor = 1.0 / (Mass(0.5, 0.5) * Mass(0.2, 0.3));
        double[,] scaled = Matrix(single.GetProperty("reference"));
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++) scaled[i, j] *= factor;
        double gap = Gap(Matrix(single.GetProperty("second_moment_fixed")), scaled);
        output.WriteLine($"single components: Π 1/z = {factor:F6}; second moment patched against the reference times it {gap:E2}");
        Assert.True(gap <= 1e-12, $"{gap:E2}");
    }

    [Fact]
    public void Readings_GammaAdditiveDefaultDrawsAndScale()
    {
        if (Fixture.Value is not { } r) return;
        JsonElement gamma = r.GetProperty("gamma_prior");
        JsonElement built = gamma.GetProperty("built"), hand = gamma.GetProperty("hand");
        bool builtIsNa = !built.GetProperty("ran").GetBoolean()
                         || built.GetProperty("C").EnumerateArray().SelectMany(row => row.EnumerateArray()).All(v => v.ValueKind == JsonValueKind.Null);
        double handC = Matrix(hand.GetProperty("C"))[0, 0], mass = hand.GetProperty("mass").GetDouble();
        JsonElement additive = r.GetProperty("additive");
        double[,] last = Matrix(r.GetProperty("default_draws").GetProperty("last")), first = Matrix(r.GetProperty("default_draws").GetProperty("first"));
        double[,] byDefault = Matrix(r.GetProperty("default_draws").GetProperty("default"));
        output.WriteLine($"build_prior's gamma through C_bass: {(builtIsNa ? "NA or an error" : "a number")}");
        output.WriteLine($"a hand-set Gamma(2, 3) on f = x: C_bass {handC:F12}, the law's mass on [0, 1] {mass:F12} (renormalized, 1)");
        output.WriteLine($"maxInt = 1: {(additive.GetProperty("ran").GetBoolean() ? "ran" : $"stops: {additive.GetProperty("error").GetString()}")}");
        output.WriteLine($"without mcmc.use: the last draw's C {Gap(byDefault, last):E1} away, the first's {Gap(byDefault, first):E2}");
        output.WriteLine($"use_native_scale defaults to {r.GetProperty("native_scale_default")}");
        Assert.True(builtIsNa);
        Assert.Equal(1.0, handC / mass, 12);
        Assert.False(additive.GetProperty("ran").GetBoolean());
        Assert.Equal(0.0, Gap(byDefault, last));
        Assert.True(Gap(byDefault, first) > 0.5);
        Assert.False(r.GetProperty("native_scale_default").GetBoolean());
    }
}
