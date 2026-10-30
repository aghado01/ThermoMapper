using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The coefficient prior's two axes in BASS (decided 2026-10-29): every member of τ's family under each covariance form.
/// τ's slice step on a fixed structure against one-dimensional quadrature of its prior on u times the collapsed marginal
/// (D03's quadrature check).
/// </summary>
public sealed class BassCoefficientPriorTests
{
    private readonly ITestOutputHelper _out;
    public BassCoefficientPriorTests(ITestOutputHelper output) => _out = output;

    internal static readonly BassTerm[] TwoTerms =
        [new BassTerm([new HingeFactor(0, 1, 2)]), new BassTerm([new HingeFactor(0, 1, 2), new HingeFactor(1, -1, 2)])];

    internal static BassCoefficientPrior Prior(string form, string kind, double a)
    {
        TauPrior scale = TauPriorTests.Member(kind, a);
        return form == "ridge" ? BassCoefficientPrior.Ridge(scale) : BassCoefficientPrior.G(scale);
    }

    internal sealed class Stay : IRjMove<BassState>
    {
        public string Key => "stay";
        public string ReverseKey => "stay";
        public double Weight(BassState state) => 1.0;
        public Proposal<BassState>? Propose(BassState current, Xoshiro256PlusPlus rng) => null;
    }

    // The batch-means error of a trace's mean.
    internal static double BatchError(double[] trace, int batches = 25)
    {
        int size = trace.Length / batches;
        double[] means = [.. Enumerable.Range(0, batches).Select(b => trace.Skip(b * size).Take(size).Average())];
        double mean = means.Average();
        return Math.Sqrt(means.Sum(m => (m - mean) * (m - mean)) / (batches - 1) / batches);
    }

    public static TheoryData<string, string, double> TauCases => new()
    {
        { "g", "gamma", 0.5 }, { "g", "gamma", 2 }, { "g", "hyper-g", 3 }, { "g", "hyper-g/n", 3 }, { "g", "hyper-g", 1.5 },
    };

    [Theory]
    [MemberData(nameof(TauCases))]
    public void TauPosterior_OnAFixedStructure_MatchesQuadrature(string form, string kind, double a)
    {
        BassConfig config = BassSweepTests.Config(maxBasis: 10) with { Coefficients = Prior(form, kind, a) };
        var model = BassModel.Create(config, BassSweepTests.Data(large: true));
        BassTarget target = model.Target;
        var moves = new BassStructureMoves(model, TwoTerms, Math.Log(target.Scale.UnitInformation));
        var chain = new ReversibleJumpChain<BassState>([new Stay()], target, moves.Start, new Xoshiro256PlusPlus(7));
        var scan = new BlockScan<ReversibleJumpChain<BassState>>([new BassTauBlock(target)], chain, new Xoshiro256PlusPlus(8), c => c);
        const int sweeps = 60_000;
        for (int i = 0; i < 1_000; i++) scan.Step();
        double[] u = [.. Enumerable.Range(0, sweeps).Select(_ => scan.Step().Current.LogTau)];

        // τ | structure, y on u: its prior's density (Jacobian included) times the collapsed marginal.
        BassStructureFit fit = moves.Start.Fit;
        // The likelihood vanishes as τ^(k/2) at τ = 0 and levels off as τ grows, so 60 below τ₀ and 45 above the prior's
        // own centre hold the posterior, improper hyper-g's included; the end weights confirm it.
        double centre = target.Scale.Declared.Kind == TauPriorKind.HyperG ? 0.0 : Math.Log(target.Scale.UnitInformation);
        double lo = Math.Min(centre, Math.Log(target.Scale.UnitInformation)) - 60, hi = Math.Max(centre, Math.Log(target.Scale.UnitInformation)) + 45;
        var nodes = BassTargetTests.Nodes(lo, hi, 0.25).ToList();
        double LogPosterior(double v) => target.LogTauDensity(v) + target.Evaluate(fit, v).LogMarginal;
        double[] logw = [.. nodes.Select(p => Math.Log(p.W) + LogPosterior(p.U))];
        double top = logw.Max();
        double[] w = [.. logw.Select(l => Math.Exp(l - top))];
        Assert.True(w[0] < 1e-20 && w[^1] < 1e-15, $"end weights {w[0]:G3}, {w[^1]:G3}");
        double norm = w.Sum();
        double Expect(Func<double, double> f) => nodes.Select((p, i) => w[i] * f(p.U)).Sum() / norm;
        // The CDF on a grid that ends at its argument, so it is continuous there.
        double Cdf(double q) => BassTargetTests.Nodes(lo, q, 0.25).Sum(p => p.W * Math.Exp(LogPosterior(p.U) - top)) / norm;

        double median = Bisect(v => Cdf(v) - 0.5, lo, hi);
        var checks = new List<(string Name, Func<double, double> F, double Exact)>
        {
            ("E[u]", x => x, Expect(x => x)), ("P(u ≤ median)", x => x <= median ? 1.0 : 0.0, 0.5),
            ("E[τ/(τ + τ₀)]", x => Math.Exp(x) / (Math.Exp(x) + target.Scale.UnitInformation),
             Expect(x => Math.Exp(x) / (Math.Exp(x) + target.Scale.UnitInformation))),
        };
        foreach ((string name, Func<double, double> f, double exact) in checks)
        {
            double[] trace = [.. u.Select(f)];
            double mean = trace.Average(), error = BatchError(trace);
            _out.WriteLine($"{form}, {target.Scale.Declared}: {name} {mean:F5}, quadrature {exact:F5} (± {error:G2})");
            Assert.True(Math.Abs(mean - exact) <= 5 * error, $"{name}: {mean:F5} against {exact:F5} (± {error:G2})");
        }
    }

    private static double Bisect(Func<double, double> f, double lo, double hi)
    {
        for (int i = 0; i < 60; i++)
        {
            double mid = 0.5 * (lo + hi);
            if (f(mid) < 0) lo = mid;
            else hi = mid;
        }
        return 0.5 * (lo + hi);
    }
}
