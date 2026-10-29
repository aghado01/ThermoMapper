using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// S10's controlled recovery: BASS's defaults on simulated surfaces of each data kind, the truth known. R's posteriors are
/// exploration only (the parity rule); the reference here is the generating surface.
/// <list type="bullet">
///   <item><b>Scalar continuous:</b> Friedman's surface, 10 sin(πx₁x₂) + 20(x₃ − ½)² + 10x₄ + 5x₅, over p = 10 uniform
///   inputs of which five are inert, n = 500, σ = 1.</item>
///   <item><b>Mixed categorical:</b> 2x₁ + 3x₂·[z ≠ c] + 2 sin(2πx₃)·[z = c] over three uniform inputs and a factor z on
///   {a, b, c}, n = 300, σ = 0.5.</item>
///   <item><b>Direct functional:</b> 1 + 2x₀ + 1.5·t·x₁ + 3(x₀ − ½)₊·sin(πt) over n = 60 cases and t on 20 locations
///   (N = 1,200), σ = 0.2.</item>
///   <item><b>Fits:</b> the defaults (interactions to 3, the count prior, the proper σ² and level centring, sampled g),
///   four chains of 5,000 burn-in and 2,000 draws.</item>
///   <item><b>Read:</b> the posterior mean's RMSE against the truth at the training points, the 95% pointwise
///   intervals' coverage of it, σ's 95% interval, the mean M, the largest R̂ the diagnostics report, and for Friedman's
///   surface the share of draws using an inert input.</item>
/// </list>
/// Pre-registered criteria: σ's interval covers the truth, and the pointwise intervals cover at least 80% of the
/// training points; the rest is read. A qualification run (<see cref="BassQualification.Enabled"/>). Untempered
/// Friedman fails its σ criterion, a failure the user kept as registered on 2026-10-28: its test stays red until the
/// local moves mix better (the fitting checklist's F-M2).
/// </summary>
public sealed class BassQualificationRecoveryTests
{
    private readonly ITestOutputHelper _out;
    public BassQualificationRecoveryTests(ITestOutputHelper output) => _out = output;

    private static BassConfig Config(int seed, TemperingSpec? tempering = null) => new()
    {
        Sampling = new BassSampling { Chains = 4, MasterSeed = seed, BurnIn = 5000, Samples = 2000, Tempering = tempering },
    };

    [Fact]
    public void Friedman_ScalarContinuous()
    {
        if (!BassQualification.Enabled) return;   // heavy
        (BassData data, double[] truth) = Friedman();
        BassFit fit = Recover("Friedman", data, truth, sigma: 1.0, Config(53));
        double inert = fit.Draws.SelectMany(c => c)
                          .Average(d => d.Structure.TermIds.Any(id => fit.Terms[id].Continuous.Any(h => h.Variable >= 5)) ? 1.0 : 0.0);
        _out.WriteLine($"Friedman: draws using an inert input {inert:P1}");
    }

    // Friedman's surface again, tempered (S8 with chip E's self-sizing ladder), the defaults otherwise: whether the
    // ladder repairs the chains' disagreement in M that the untempered run shows. Added after that run's read.
    [Fact]
    public void Friedman_Tempered()
    {
        if (!BassQualification.Enabled) return;   // heavy
        (BassData data, double[] truth) = Friedman();
        Recover("Friedman, tempered", data, truth, sigma: 1.0, Config(53, new TemperingSpec()));
    }

    private static (BassData Data, double[] Truth) Friedman()
    {
        var rng = new Xoshiro256PlusPlus(seed: 51);
        const int n = 500, p = 10;
        var x = new double[p][];
        for (int v = 0; v < p; v++) x[v] = new double[n];
        var truth = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            for (int v = 0; v < p; v++) x[v][i] = rng.NextDouble();
            truth[i] = 10 * Math.Sin(Math.PI * x[0][i] * x[1][i]) + 20 * Math.Pow(x[2][i] - 0.5, 2) + 10 * x[3][i] + 5 * x[4][i];
            y[i] = truth[i] + BassQualification.Gaussian(rng);
        }
        return (new BassData([.. Enumerable.Range(0, p).Select(v => (BassPredictor)new ContinuousPredictor($"x{v + 1}", x[v]))], new ScalarResponse(y)), truth);
    }

    [Fact]
    public void MixedCategorical()
    {
        if (!BassQualification.Enabled) return;   // heavy
        var rng = new Xoshiro256PlusPlus(seed: 55);
        const int n = 300;
        var x = new double[3][];
        for (int v = 0; v < 3; v++) x[v] = new double[n];
        var z = new int[n];
        var truth = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            for (int v = 0; v < 3; v++) x[v][i] = rng.NextDouble();
            z[i] = (int)(3 * rng.NextDouble());
            truth[i] = 2 * x[0][i] + (z[i] != 2 ? 3 * x[1][i] : 2 * Math.Sin(2 * Math.PI * x[2][i]));
            y[i] = truth[i] + 0.5 * BassQualification.Gaussian(rng);
        }
        var data = new BassData(
            [new ContinuousPredictor("x1", x[0]), new ContinuousPredictor("x2", x[1]), new ContinuousPredictor("x3", x[2]),
             new CategoricalPredictor("z", z, ["a", "b", "c"])],
            new ScalarResponse(y));
        Recover("mixed categorical", data, truth, sigma: 0.5, Config(57));
    }

    [Fact]
    public void DirectFunctional()
    {
        if (!BassQualification.Enabled) return;   // heavy
        var rng = new Xoshiro256PlusPlus(seed: 59);
        const int n = 60, m = 20;
        double[] t = [.. Enumerable.Range(0, m).Select(r => r / (m - 1.0))];
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n, m];
        var truth = new double[n * m];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            for (int r = 0; r < m; r++)
            {
                double f = 1 + 2 * x0[i] + 1.5 * t[r] * x1[i] + 3 * Math.Max(x0[i] - 0.5, 0) * Math.Sin(Math.PI * t[r]);
                truth[i * m + r] = f;
                y[i, r] = f + 0.2 * BassQualification.Gaussian(rng);
            }
        }
        var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                                new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
        BassFit fit = Recover("direct functional", data, truth, sigma: 0.2, Config(61));
        ResidualDependence? check = fit.ResidualDependence;
        _out.WriteLine($"direct functional: residual r₁ mean {check?.Mean:F3} against {check?.Expected:F3}");
    }

    // Fits the defaults, reads the recovery and checks the pre-registered criteria.
    private BassFit Recover(string name, BassData data, double[] truth, double sigma, BassConfig config)
    {
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassFit fit = Bass.Fit(config, data);
        double seconds = clock.Elapsed.TotalSeconds;
        double[][] fitted = BassQualification.Fitted(BassModel.Create(config, data), fit);
        int points = truth.Length, covered = 0;
        double squares = 0.0;
        for (int k = 0; k < points; k++)
        {
            (double mean, _, double lower, double upper) = BassQualification.At(fitted, k);
            squares += (mean - truth[k]) * (mean - truth[k]);
            if (lower <= truth[k] && truth[k] <= upper) covered++;
        }
        double[] sd = [.. fit.Draws.SelectMany(c => c).Select(d => Math.Sqrt(d.Sigma2)).Order()];
        double sdLower = PosteriorInterval.Quantile(sd, 0.025), sdUpper = PosteriorInterval.Quantile(sd, 0.975);
        BassDiagnostics diagnostics = fit.Diagnostics;
        double rHat = new[] { diagnostics.Count, diagnostics.Sigma2 }.Concat(diagnostics.Fitted.Select(f => f.Diagnostic))
                                                                    .Max(d => d.RHat.Value);
        double meanM = fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length);
        double coverage = covered / (double)points;
        _out.WriteLine($"{name}: RMSE {Math.Sqrt(squares / points):F4}, 95% coverage {coverage:P1} of {points}, σ in " +
                       $"[{sdLower:F4}, {sdUpper:F4}] against {sigma}, mean M {meanM:F2}, largest R̂ {rHat:F3}, {seconds:F1} s");
        _out.WriteLine($"{name}: R̂ of M {diagnostics.Count.RHat.Value:F3}, σ² {diagnostics.Sigma2.RHat.Value:F3}, τ {diagnostics.Tau?.RHat.Value:F3}, " +
                       $"fitted {string.Join(", ", diagnostics.Fitted.Select(f => f.Diagnostic.RHat.Value.ToString("F3")))}");
        _out.WriteLine($"{name}: moves {string.Join(", ", diagnostics.Moves.Select(mv => $"{mv.Key} {mv.Accepted}/{mv.Attempts}"))}; " +
                       $"refusals {string.Join(", ", Enum.GetValues<BassRefusal>().Select(r => $"{r} {diagnostics.Refusals.Count(r)}"))}");
        if (diagnostics.Tempering is { } t)
            _out.WriteLine($"{name}: levels {string.Join(", ", t.Levels)}, largest pair rejection {t.MaxPairRejection:F3}, barrier {t.Barrier:F2}, " +
                           $"round trips per scan {t.RoundTripsPerScan:G3} against the ideal {t.IdealRoundTripRate:G3}");
        // Where a ladder's barrier sits: each chain's pairs above 0.9 rejection, by their rungs b.
        for (int c = 0; c < (diagnostics.Tempering?.PairRejection.Length ?? 0); c++)
        {
            double[] r = diagnostics.Tempering!.PairRejection[c], b = diagnostics.Tempering.Ladders[c];
            string steep = string.Join(", ", Enumerable.Range(0, r.Length).Where(i => r[i] > 0.9).Select(i => $"[{b[i]:G3}, {b[i + 1]:G3}] {r[i]:F3}"));
            _out.WriteLine($"{name}: chain {c}: pair rejection median {r.Order().ElementAt(r.Length / 2):F3}, hot end b {b[^1]:G3}; " +
                           $"pairs above 0.9 {(steep.Length == 0 ? "none" : steep)}");
        }
        foreach ((DrawStore<BassDraw> chain, int c) in fit.Draws.Select((chain, c) => (chain, c)))
            _out.WriteLine($"{name}: chain {c}: mean M {chain.Average(d => d.Structure.TermIds.Length):F2}, mean σ {chain.Average(d => Math.Sqrt(d.Sigma2)):F4}, " +
                           $"mean τ {chain.Average(d => d.Tau):G3}, M range {chain.Min(d => d.Structure.TermIds.Length)}–{chain.Max(d => d.Structure.TermIds.Length)}");
        Assert.True(sdLower <= sigma && sigma <= sdUpper, $"{name}: σ's interval [{sdLower:F4}, {sdUpper:F4}] misses {sigma}");
        Assert.True(coverage >= 0.80, $"{name}: the pointwise intervals cover {coverage:P1}");
        return fit;
    }
}
