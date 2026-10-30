using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The shared τ family (<see cref="TauPrior"/>; decided 2026-10-29), each member on its form's unit-information precision
/// τ₀ = 1/n₀: the density on u = ln τ normalizes by quadrature on u and carries its Jacobian, its moments meet their
/// closed forms, the exact draws meet the law within 5 SE, and the improper and held members refuse what they cannot
/// give.
/// </summary>
public sealed class TauPriorTests
{
    private readonly ITestOutputHelper _out;
    public TauPriorTests(ITestOutputHelper output) => _out = output;

    public static TheoryData<string, double, double> Members => new()
    {
        { "gamma", 0.1, 1 }, { "gamma", 0.5, 40 }, { "gamma", 3, 1 }, { "gamma", 0.5, 7500 },
        { "hyper-g", 2.5, 40 }, { "hyper-g", 4, 1 }, { "hyper-g/n", 3, 40 }, { "hyper-g/n", 2.2, 1 },
    };

    internal static TauPrior Member(string kind, double a) => kind switch
    {
        "gamma" => TauPrior.UnitInformationGamma(a),
        "hyper-g" => TauPrior.HyperG(a),
        _ => TauPrior.HyperGOverN(a),
    };

    // A range of u holding all but about e⁻⁴⁰ of the member's mass: the Gamma's left tail falls as e^(au), hyper-g's as
    // e^((a/2 − 1)u), both right tails at least as e^(−u).
    internal static (double Lo, double Hi) Range(ResolvedTauPrior prior)
    {
        double centre = Math.Log(prior.UnitInformation), a = prior.Declared.Shape!.Value;
        double rate = prior.Gamma is null ? 0.5 * a - 1.0 : a;
        double offset = prior.Declared.Kind == TauPriorKind.HyperG ? 0.0 : centre;
        return (offset - 45 / rate, offset + 45);
    }

    // E over the member's law of f(u), by composite Gauss–Legendre on u.
    internal static double Expect(ResolvedTauPrior prior, Func<double, double> f)
    {
        (double lo, double hi) = Range(prior);
        double sum = 0.0;
        foreach ((double u, double w) in BassTargetTests.Nodes(lo, hi, 0.25)) sum += w * Math.Exp(prior.LogDensity(u)) * f(u);
        return sum;
    }

    [Theory]
    [MemberData(nameof(Members))]
    public void Density_NormalizesOnU_AndMeetsItsClosedForms(string kind, double a, double n0)
    {
        ResolvedTauPrior prior = Member(kind, a).For(n0);
        double mass = Expect(prior, _ => 1.0);
        _out.WriteLine($"{prior}: mass {mass:R}");
        Assert.True(Math.Abs(mass - 1.0) <= 1e-12, $"mass {mass:R}");
        if (kind == "gamma")
        {
            // E[τ] = τ₀ at every shape, and E[ln τ] = ψ(a) − ln(a·n₀) (F-O8's identities).
            Assert.Equal(1.0, Expect(prior, u => Math.Exp(u)) * n0, 12);
            double log = SpecialFunctions.Digamma(a) - Math.Log(a * n0);
            Assert.True(Math.Abs(Expect(prior, u => u) - log) <= 1e-11 * Math.Max(1.0, Math.Abs(log)));
            return;
        }
        // g′/(1 + g′) ~ Beta(1, a/2 − 1), with g′ = s/τ, s = 1 (hyper-g) or τ₀ (hyper-g/n): its mean is 2/a and its
        // second moment E[B²] = 2/((β + 1)(β + 2)), β = a/2 − 1.
        double s = kind == "hyper-g" ? 1.0 : prior.UnitInformation, beta = 0.5 * a - 1.0;
        Assert.Equal(1.0 / (1.0 + beta), Expect(prior, u => s / (Math.Exp(u) + s)), 11);
        Assert.Equal(2.0 / ((beta + 1.0) * (beta + 2.0)), Expect(prior, u => Math.Pow(s / (Math.Exp(u) + s), 2)), 11);
    }

    [Theory]
    [MemberData(nameof(Members))]
    public void Density_CarriesItsJacobian(string kind, double a, double n0)
    {
        ResolvedTauPrior prior = Member(kind, a).For(n0);
        foreach (double u in new[] { -30.0, -8.0, Math.Log(prior.UnitInformation), -0.3, 0.0, 2.0, 12.0 })
        {
            double tau = Math.Exp(u), expected;
            if (kind == "gamma")
            {
                // The Gamma's density on τ, times |dτ/du| = τ.
                double rate = a * n0;
                expected = a * Math.Log(rate) - SpecialFunctions.LogGamma(a) + (a - 1) * Math.Log(tau) - rate * tau + u;
            }
            else
            {
                // Hyper-g's density on g′ = s/τ, times |dg′/du| = g′.
                double g = (kind == "hyper-g" ? 1.0 : 1.0 / n0) / tau;
                expected = Math.Log(0.5 * (a - 2)) - 0.5 * a * Math.Log(1 + g) + Math.Log(g);
            }
            Assert.True(Math.Abs(prior.LogDensity(u) - expected) <= 1e-12 * Math.Max(1.0, Math.Abs(expected)),
                        $"{prior} at u = {u}: {prior.LogDensity(u):R} against {expected:R}");
        }
        Assert.Equal(double.NegativeInfinity, prior.LogDensity(-800));
        Assert.Equal(double.NegativeInfinity, prior.LogDensity(800));
    }

    [Theory]
    [MemberData(nameof(Members))]
    public void Draws_MeetTheLaw_Within5SE(string kind, double a, double n0)
    {
        ResolvedTauPrior prior = Member(kind, a).For(n0);
        var rng = new Xoshiro256PlusPlus(17);
        const int draws = 40_000;
        double[] u = [.. Enumerable.Range(0, draws).Select(_ => prior.SampleLogTau(rng))];
        Assert.All(u, v => Assert.True(double.IsFinite(prior.LogDensity(v))));

        // Every bounded statistic's mean against its quadrature, 5 SE: the CDF at three of the law's quantiles, and a
        // smooth one, s/(τ + s), on the scale the member is stated on.
        double s = kind == "hyper-g" ? 1.0 : prior.UnitInformation;
        var statistics = new List<(string Name, Func<double, double> F, double Expected)>
        {
            ("s/(τ + s)", v => s / (Math.Exp(v) + s), Expect(prior, v => s / (Math.Exp(v) + s))),
        };
        foreach (double level in new[] { 0.1, 0.5, 0.9 })
        {
            double q = Quantile(prior, level);
            statistics.Add(($"P(u ≤ q{level})", v => v <= q ? 1.0 : 0.0, level));
        }
        foreach ((string name, Func<double, double> f, double expected) in statistics)
        {
            double[] values = [.. u.Select(f)];
            double mean = values.Average(), sd = Math.Sqrt(values.Select(x => (x - mean) * (x - mean)).Sum() / (draws - 1));
            double error = sd / Math.Sqrt(draws);
            _out.WriteLine($"{prior} {name}: draws {mean:F5}, law {expected:F5} (± {error:G2})");
            Assert.True(Math.Abs(mean - expected) <= 5 * error, $"{prior} {name}: {mean:F5} against {expected:F5} (± {error:G2})");
        }
    }

    // The level-quantile of u, by bisection on the CDF, each by quadrature on a grid that ends at its argument.
    internal static double Quantile(ResolvedTauPrior prior, double level)
    {
        (double lo, double hi) = Range(prior);
        double start = lo;
        for (int i = 0; i < 60; i++)
        {
            double mid = 0.5 * (lo + hi);
            double cdf = BassTargetTests.Nodes(start, mid, 0.25).Sum(p => p.W * Math.Exp(prior.LogDensity(p.U)));
            if (cdf < level) lo = mid;
            else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    [Theory]
    [InlineData(1.5)]
    [InlineData(2.0)]
    public void HyperG_AtOrBelowTwo_IsImproper_StatedUpToItsConstant(double a)
    {
        foreach (TauPrior declared in new[] { TauPrior.HyperG(a), TauPrior.HyperGOverN(a) })
        {
            Assert.False(declared.IsProper);
            ResolvedTauPrior prior = declared.For(40);
            Assert.False(prior.IsProper);
            Assert.Throws<InvalidOperationException>(() => prior.SampleLogTau(new Xoshiro256PlusPlus(1)));
            double s = declared.Kind == TauPriorKind.HyperG ? 1.0 : 1.0 / 40;
            foreach (double u in new[] { -20.0, 0.0, 5.0 })
            {
                double g = s / Math.Exp(u);
                Assert.Equal(-0.5 * a * Math.Log(1 + g) + Math.Log(g), prior.LogDensity(u), 12);
            }
            // Its left tail falls as e^((a/2 − 1)u), no faster: not integrable at τ = 0.
            Assert.Equal(0.5 * a - 1.0, prior.LogDensity(-200) - prior.LogDensity(-201), 9);
        }
    }

    [Fact]
    public void Held_HoldsItsValue_OrTauZero()
    {
        ResolvedTauPrior unit = TauPrior.Held().For(40), declared = TauPrior.Held(0.3).For(40), ridge = TauPrior.Held(2.5).For(1);
        Assert.Equal((1.0 / 40, 0.3, 2.5), (unit.HeldTau, declared.HeldTau, ridge.HeldTau));
        Assert.True(unit.IsHeld && unit.IsProper && TauPrior.Held().IsHeld);
        Assert.Throws<InvalidOperationException>(() => unit.LogDensity(0.0));
        Assert.Throws<InvalidOperationException>(() => unit.SampleLogTau(new Xoshiro256PlusPlus(1)));
        Assert.False(TauPrior.UnitInformationGamma().IsHeld);
        Assert.Throws<ArgumentOutOfRangeException>(() => TauPrior.Held().For(0));
    }

    [Theory]
    [InlineData(40, 1, 0.2, 3.0)]
    [InlineData(40, 2, 0.5, 4.0)]
    [InlineData(1000, 1, 0.2, 3.0)]
    [InlineData(1000, 3, 0.05, 2.5)]
    public void HyperG_BayesFactor_MeetsLiangsClosedForm(int n, int k, double r2, double a)
    {
        // Liang et al.'s eq. (17), under their flat intercept: BF[k columns : null] = (a − 2)/(k + a − 2)·₂F₁((n − 1)/2, 1;
        // (k + a)/2; R²), against the quadrature on u of hyper-g's density times (1 + g)^((n−1−k)/2)(1 + g(1 − R²))^(−(n−1)/2).
        ResolvedTauPrior prior = TauPrior.HyperG(a).For(n);
        double LogIntegrand(double u)
        {
            double g = Math.Exp(-u);
            return prior.LogDensity(u) + 0.5 * (n - 1 - k) * Math.Log(1 + g) - 0.5 * (n - 1) * Math.Log(1 + g * (1 - r2));
        }
        var nodes = BassTargetTests.Nodes(-90, 45, 0.1).ToList();
        double[] logw = [.. nodes.Select(p => Math.Log(p.W) + LogIntegrand(p.U))];
        double top = logw.Max(), logBf = top + Math.Log(logw.Sum(l => Math.Exp(l - top)));

        double alpha = 0.5 * (n - 1), gamma = 0.5 * (k + a), logZ = Math.Log(r2);
        var terms = new List<double>();
        for (int j = 0; ; j++)
        {
            terms.Add(SpecialFunctions.LogGamma(alpha + j) - SpecialFunctions.LogGamma(alpha) - SpecialFunctions.LogGamma(gamma + j)
                      + SpecialFunctions.LogGamma(gamma) + j * logZ);
            if (j > 10 && terms[^1] < terms.Max() - 50) break;
        }
        double peak = terms.Max();
        double closed = Math.Log((a - 2) / (k + a - 2)) + peak + Math.Log(terms.Sum(t => Math.Exp(t - peak)));
        _out.WriteLine($"N {n}, k {k}, R² {r2}, a {a}: ln BF quadrature {logBf:R}, closed {closed:R}");
        Assert.True(Math.Abs(logBf - closed) <= 1e-9 * Math.Max(1.0, Math.Abs(closed)), $"{logBf:R} against {closed:R}");
    }

    [Fact]
    public void HyperGOverN_OnTheRidge_IsHyperG()
    {
        // τ₀ = 1 under the normalized ridge, so the two members coincide there.
        ResolvedTauPrior g = TauPrior.HyperG(3).For(1), gn = TauPrior.HyperGOverN(3).For(1);
        foreach (double u in new[] { -9.0, -1.0, 0.0, 3.5 }) Assert.Equal(g.LogDensity(u), gn.LogDensity(u));
        var r1 = new Xoshiro256PlusPlus(5);
        var r2 = new Xoshiro256PlusPlus(5);
        for (int i = 0; i < 100; i++) Assert.Equal(g.SampleLogTau(r1), gn.SampleLogTau(r2));
    }
}
