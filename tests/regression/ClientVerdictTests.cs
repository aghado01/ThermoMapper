using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Bgp;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// D10's verdict on the ensemble's other clients (N14): BARS's two groups, predictive (the fitted grid values) and
/// structural (the knot count), beside its R̂ readouts; and BGP, BAPS and the hierarchical P-spline, whose results carry the
/// verdict on their one group with their functionals named. Each verdict is the shared rule read on the result's own
/// diagnostics.
/// </summary>
public sealed class ClientVerdictTests
{
    private readonly ITestOutputHelper _out;
    public ClientVerdictTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Bars_ReadsItsFitAndItsKnotCount_InTwoGroups()
    {
        (double[] x, double[] y) = Curve(80, seed: 4);
        var config = new BarsConfig { Chains = 3, MasterSeed = 9, BurnIn = 300, MaxSamples = 300 };
        BarsResult r = Bars.Run(config, x, y, x);
        ConvergenceVerdict v = r.Verdict;
        Assert.Equal(new[] { BarsEnsemble.PredictiveGroup, BarsEnsemble.StructuralGroup }, v.Groups.Select(g => g.Name));
        Assert.Equal(ChainEnsemble.ReferenceRHatTarget, v.Target);

        var readings = new Dictionary<string, (RHatDiagnostic, EssDiagnostic?)> { ["knots"] = (r.Knots.RHat, r.Knots.Ess) };
        for (int g = 0; g < x.Length; g++) readings[$"fit[{g}]"] = (r.RHatFit[g], null);
        VerdictAssert.Consistent(v, readings, config.Chains);
        Assert.All(v.Group(BarsEnsemble.PredictiveGroup).Failing, f => Assert.Null(f.Ess));   // the grid is unmonitored
        _out.WriteLine($"BARS: {v.Outcome}; " + string.Join("; ", v.Groups.Select(g => $"{g.Name} {g.Outcome}, {g.Failing.Count} failing")));

        // An R̂ stop ends Converged; the stop and the verdict read the same groups.
        BarsResult stopped = Bars.Run(config with { RHatTarget = 1.5, BatchSize = 50, MinEssPerChain = 10, MaxSamples = 3000 }, x, y, x);
        Assert.True(stopped.SamplesUsed < 3000, $"ran {stopped.SamplesUsed}");
        Assert.Equal(ConvergenceOutcome.Converged, stopped.Verdict.Outcome);
        Assert.Equal(1.5, stopped.Verdict.Target);
    }

    [Fact]
    public void Bgp_CarriesTheVerdict_OnItsNamedFunctionals()
    {
        (double[] x1, double[] y) = Curve(30, seed: 2);
        var x = new double[x1.Length, 1];
        for (int i = 0; i < x1.Length; i++) x[i, 0] = x1[i];
        var known = new GpRegression(x, y, 0.01, new SquaredExponentialKernel());
        BgpResult r = new BgpSampler(known, new EmpiricalBayesBandwidthPrior(known)).Run(x, draws: 300, burn: 100, proposalSd: 0.5, seed: 5, chains: 3);
        VerdictAssert.Consistent(r.Verdict, new() { ["bandwidth"] = (r.BandwidthConvergence.RHat, r.BandwidthConvergence.Ess) }, 3);
        Assert.Equal(FunctionalGroup.AllName, Assert.Single(r.Verdict.Groups).Name);

        GpRegression unknown = GpRegression.WithUnknownScale(x, y, new SquaredExponentialKernel());
        BgpResult u = new BgpSampler(unknown, new EmpiricalBayesBandwidthPrior(unknown)).Run(x, draws: 300, burn: 100, proposalSd: 0.5, seed: 5, chains: 3);
        ConvergenceDiagnostic eta = u.NoiseRatioConvergence!.Value;
        VerdictAssert.Consistent(u.Verdict, new()
        {
            ["bandwidth"] = (u.BandwidthConvergence.RHat, u.BandwidthConvergence.Ess), ["noise_ratio"] = (eta.RHat, eta.Ess),
        }, 3);
    }

    [Theory]
    [InlineData(BapsScheme.Gibbs)]
    [InlineData(BapsScheme.Collapsed)]
    public void Baps_CarriesTheVerdict_OnItsNamedFunctionals(BapsScheme scheme)
    {
        (double[] x, double[] y) = Curve(60, seed: 3);
        var knots = new double[12];
        for (int k = 0; k < knots.Length; k++) knots[k] = (k + 1.0) / (knots.Length + 1);
        double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
        BapsResult r = new Baps(z, y, new DifferencePenalty(2)).Run(scheme, chains: 3, masterSeed: 5, burn: 200, samples: 300);
        VerdictAssert.Consistent(r.Verdict, new()
        {
            ["log_lambda"] = (r.LogLambdaConvergence.RHat, r.LogLambdaConvergence.Ess),
            ["log_sigma2"] = (r.LogSigma2Convergence.RHat, r.LogSigma2Convergence.Ess),
        }, 3);
    }

    [Fact]
    public void Hierarchical_CarriesTheVerdict_OnItsNamedFunctionals()
    {
        (double[] x, double[] y0) = Curve(40, seed: 6);
        var rng = new Xoshiro256PlusPlus(12);
        double[][] y = [.. Enumerable.Range(0, 4).Select(_ => y0.Select(v => v + 0.2 * (rng.NextDouble() - 0.5)).ToArray())];
        var knots = new double[10];
        for (int k = 0; k < knots.Length; k++) knots[k] = (k + 1.0) / (knots.Length + 1);
        double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
        HierarchicalResult r = new HierarchicalPSpline(z, y, new DifferencePenalty(2)).Run(chains: 3, masterSeed: 3, burn: 200, samples: 300);
        var readings = new Dictionary<string, (RHatDiagnostic, EssDiagnostic?)>
        {
            ["log_sigma2"] = (r.NoiseVarianceConvergence.RHat, r.NoiseVarianceConvergence.Ess),
            ["log_tau_u2"] = (r.DeviationVarianceConvergence.RHat, r.DeviationVarianceConvergence.Ess),
            ["log_tau_02"] = (r.PopulationSmoothingVarianceConvergence.RHat, r.PopulationSmoothingVarianceConvergence.Ess),
        };
        for (int j = 0; j < r.RHatPopulation.Length; j++) readings[$"population[{j}]"] = (r.RHatPopulation[j], null);
        VerdictAssert.Consistent(r.Verdict, readings, 3);
    }

    // y = sin 2πx plus noise of sd 0.1 on n points of (0, 1).
    private static (double[] X, double[] Y) Curve(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.1 * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        return (x, y);
    }
}
