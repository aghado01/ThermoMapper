using System;
using System.Diagnostics;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// Rank-normalized split R̂ and bulk and tail ESS (Vehtari, Gelman, Simpson, Carpenter and Bürkner 2021), against cases
/// with known answers:
/// <list type="bullet">
/// <item>independent draws: ESS ≈ S, R̂ ≈ 1;</item>
/// <item>AR(1) with coefficient φ: ESS ≈ S(1 − φ)/(1 + φ), above S when φ &lt; 0;</item>
/// <item>chains with one mean but different scales, which classic R̂ misses and the folded R̂ catches;</item>
/// <item>Cauchy draws, which have no variance for classic R̂ to use, and which rank normalization handles;</item>
/// <item>the failure policy: constant, frozen, too short, non-finite.</item>
/// </list>
/// </summary>
public sealed class ConvergenceDiagnosticsTests
{
    private readonly ITestOutputHelper _out;
    public ConvergenceDiagnosticsTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // C chains of n draws of a stationary Gaussian AR(1) with coefficient phi (phi = 0: independent).
    private static double[][] Ar1(int c, int n, double phi, int seed, double scale = 1.0, double shift = 0.0)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var chains = new double[c][];
        double sd = Math.Sqrt(1.0 - phi * phi);
        for (int j = 0; j < c; j++)
        {
            chains[j] = new double[n];
            double x = Gaussian(rng);
            for (int i = 0; i < n; i++)
            {
                x = phi * x + sd * Gaussian(rng);
                chains[j][i] = shift + scale * x;
            }
        }
        return chains;
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(0.9)]
    [InlineData(-0.5)]
    public void Ess_MatchesTheAr1Theory(double phi)
    {
        double[][] chains = Ar1(4, 20_000, phi, seed: 3);
        double total = 4 * 20_000.0;
        double expected = total * (1.0 - phi) / (1.0 + phi);
        ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(chains);
        _out.WriteLine($"[ess] φ = {phi}: bulk {d.BulkEss:F0}, tail {d.TailEss:F0}, theory {expected:F0}; R̂ {d.RHat:F4}");

        Assert.Equal(DiagnosticState.Estimated, d.BulkEss.State);
        Assert.InRange(d.BulkEss.Value / expected, 0.88, 1.12);
        Assert.InRange(d.RHat.Value, 0.99, 1.01);
        // Tail ESS is an exceedance indicator's: independent draws give S, and a Gaussian process's exceedances
        // decorrelate faster than the process, so for φ > 0 it lies between the bulk ESS and S.
        if (phi == 0.0) Assert.InRange(d.TailEss.Value / total, 0.8, 1.2);
        if (phi > 0.0) Assert.InRange(d.TailEss.Value, d.BulkEss.Value, total);
    }

    [Fact]
    public void FoldedRHat_CatchesAScaleMismatch_ThatClassicRHatMisses()
    {
        double[][] chains = Ar1(4, 2_000, 0.0, seed: 8);
        double[][] wide = Ar1(1, 2_000, 0.0, seed: 9, scale: 3.0);
        chains[3] = wide[0];   // same mean, three times the spread

        ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(chains);
        var means = new double[4];
        var variances = new double[4];
        for (int j = 0; j < 4; j++) { RunningMoments m = RunningMoments.Of(chains[j]); means[j] = m.Mean; variances[j] = m.Variance; }
        RHatDiagnostic classic = ChainDiagnostics.RHat(means.AsSpan(), variances.AsSpan(), 2_000);
        _out.WriteLine($"[rhat] classic {classic:F4}, bulk {d.BulkRHat:F4}, tail {d.TailRHat:F4}");

        Assert.True(classic.Value < 1.01, $"classic R̂ {classic.Value:F4}");
        Assert.True(d.BulkRHat.Value < 1.01, $"bulk R̂ {d.BulkRHat.Value:F4}");
        Assert.True(d.TailRHat.Value > 1.05, $"tail R̂ {d.TailRHat.Value:F4}");
        Assert.Equal(d.TailRHat, d.RHat);
    }

    [Fact]
    public void BulkRHat_CatchesAShiftedChain()
    {
        double[][] chains = Ar1(4, 2_000, 0.5, seed: 12);
        double[][] shifted = Ar1(1, 2_000, 0.5, seed: 13, shift: 1.0);
        chains[1] = shifted[0];
        ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(chains);
        Assert.True(d.BulkRHat.Value > 1.1, $"bulk R̂ {d.BulkRHat.Value:F4}");
    }

    [Fact]
    public void RankNormalization_HandlesCauchyDraws()
    {
        var rng = new Xoshiro256PlusPlus(seed: 4);
        var chains = new double[4][];
        for (int j = 0; j < 4; j++)
        {
            chains[j] = new double[5_000];
            for (int i = 0; i < chains[j].Length; i++) chains[j][i] = Math.Tan(Math.PI * (rng.NextDouble() - 0.5));
        }
        ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(chains);
        _out.WriteLine($"[cauchy] bulk {d.BulkEss:F0}, tail {d.TailEss:F0}, R̂ {d.RHat:F4}");
        Assert.InRange(d.BulkEss.Value / 20_000.0, 0.88, 1.12);
        Assert.InRange(d.TailEss.Value / 20_000.0, 0.8, 1.2);
        Assert.InRange(d.RHat.Value, 0.99, 1.01);
    }

    [Fact]
    public void DiscreteDraws_WithTies_AreDiagnosed()
    {
        var rng = new Xoshiro256PlusPlus(seed: 6);
        var chains = new double[4][];
        for (int j = 0; j < 4; j++)
        {
            chains[j] = new double[4_000];
            for (int i = 0; i < chains[j].Length; i++) chains[j][i] = Math.Floor(4.0 * rng.NextDouble());   // knot-count-like
        }
        ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(chains);
        Assert.Equal(DiagnosticState.Estimated, d.BulkEss.State);
        Assert.InRange(d.BulkEss.Value / 16_000.0, 0.88, 1.12);
        Assert.InRange(d.RHat.Value, 0.99, 1.01);
    }

    [Fact]
    public void FailurePolicy_MatchesRHats()
    {
        double[][] Fill(int c, int n, Func<int, int, double> f)
        {
            var x = new double[c][];
            for (int j = 0; j < c; j++) { x[j] = new double[n]; for (int i = 0; i < n; i++) x[j][i] = f(j, i); }
            return x;
        }

        ConvergenceDiagnostic constant = ChainDiagnostics.Diagnose(Fill(4, 100, (_, _) => 0.1));
        Assert.Equal(DiagnosticState.Constant, constant.RHat.State);
        Assert.Equal(DiagnosticState.Constant, constant.BulkEss.State);
        Assert.True(constant.Ess.MeetsMinimum(1e9));

        ConvergenceDiagnostic frozen = ChainDiagnostics.Diagnose(Fill(4, 100, (j, _) => j));
        Assert.Equal(DiagnosticState.Frozen, frozen.RHat.State);
        Assert.False(frozen.Ess.MeetsMinimum(0.0));

        ConvergenceDiagnostic jumps = ChainDiagnostics.Diagnose(Fill(1, 100, (_, i) => i < 50 ? 0.0 : 1.0));   // one chain, halves differ
        Assert.Equal(DiagnosticState.Frozen, jumps.RHat.State);

        Assert.Equal(DiagnosticState.Unavailable, ChainDiagnostics.Diagnose(Fill(4, 3, (_, i) => i)).RHat.State);
        Assert.Equal(DiagnosticState.Unavailable,
                     ChainDiagnostics.Diagnose(Fill(4, 100, (j, i) => j == 2 && i == 7 ? double.NaN : i)).BulkEss.State);
        Assert.Throws<ArgumentException>(() => ChainDiagnostics.Diagnose([new double[10], new double[11]]));
    }

    [Fact]
    public void OneChain_IsDiagnosedThroughItsHalves()
    {
        ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(Ar1(1, 10_000, 0.0, seed: 2));
        Assert.Equal(DiagnosticState.Estimated, d.RHat.State);
        Assert.InRange(d.BulkEss.Value / 10_000.0, 0.85, 1.15);
    }

    [Fact]
    public void Ess_IsLinearithmic()
    {
        double[][] chains = Ar1(4, 200_000, 0.9, seed: 5);
        var clock = Stopwatch.StartNew();
        ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(chains);
        double seconds = clock.Elapsed.TotalSeconds;
        _out.WriteLine($"[ess-cost] 800,000 draws in {seconds:F2} s; bulk {d.BulkEss:F0}");
        Assert.True(seconds < 10.0, $"diagnostics took {seconds:F1} s");   // the lag-by-lag variogram took minutes here
    }
}
