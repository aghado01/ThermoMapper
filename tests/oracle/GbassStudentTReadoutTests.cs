using System;
using System.Linq;
using System.Text.Json;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// Student t BASS against GBASS's <c>tbass</c> (E1's slice 1 gate), posterior readouts on <see cref="GbassReadoutParityTests"/>' data and
/// grid, ν = 3 and c = (ν − 2)/ν, so σ² is the noise variance on both sides, which is <c>tbass</c>'s w. The chains differ by design and
/// so do the priors, each difference stated: GBASS draws knots from U(0, 1) where the engine takes D01's data grid; mutates a whole term
/// where <c>Change</c> moves one factor; weights orders and variables with pseudo-counts 1 where the engine takes 5; holds a ridge
/// a ~ N(0, wτI) on the raw columns with the intercept inside and τ ~ IG(½, N/2), where the engine's default is D03's g-prior on centred
/// columns with Zellner–Siow τ (the ridge option is read beside it); floors w at var(y)/N, which the engine does not; and samples λ
/// with a Gamma(1, 1) prior where D02 integrates it under Gamma(10, 10). So the surfaces meet within their bands, w within a fifth, and
/// M is read, not held.
/// </summary>
public sealed class GbassStudentTReadoutTests(ITestOutputHelper output)
{
    private static readonly Lazy<JsonElement?> Tbass = new(() => GbassReadoutParityTests.Readout("tbass", 6000, 3000, "df=3"));

    private static BassData Data()
    {
        double[][] rows = GbassReadoutParityTests.Data();
        return new BassData([new ContinuousPredictor("x1", [.. rows.Select(r => r[0])]), new ContinuousPredictor("x2", [.. rows.Select(r => r[1])])],
                            new ScalarResponse([.. rows.Select(r => r[2])]));
    }

    private static (double[] Mean, double[] Lo, double[] Hi, double W, double WLo, double WHi, double M) Engine(BassCoefficientPrior form, int seed)
    {
        BassData data = Data();
        double[][] grid = GbassReadoutParityTests.Grid();
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2 },
            Coefficients = form,
            Noise = BassNoise.StudentT(3, 1.0 / 3),
            Sampling = new BassSampling { Chains = 4, MasterSeed = seed, BurnIn = 3000, Samples = 1500, Thin = 2 },
        };
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        var cases = new PredictionCases([new ContinuousColumn("x1", [.. grid.Select(g => g[0])]), new ContinuousColumn("x2", [.. grid.Select(g => g[1])])]);
        BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Predictive = false });
        int points = grid.Length;
        var mean = new double[points];
        var lo = new double[points];
        var hi = new double[points];
        double[,] m = p.Mean();
        for (int i = 0; i < points; i++)
        {
            double[] f = [.. Enumerable.Range(0, p.DrawCount).Select(d => p.Latent[d, i, 0]).Order()];
            mean[i] = m[i, 0];
            lo[i] = f[(int)(0.05 * f.Length)];
            hi[i] = f[(int)(0.95 * f.Length) - 1];
        }
        double[] w = [.. fit.Draws.SelectMany(c => c).Select(d => d.Sigma2).Order()];
        return (mean, lo, hi, w.Average(), w[(int)(0.05 * w.Length)], w[(int)(0.95 * w.Length) - 1],
                fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length));
    }

    [Fact]
    public void TBass_MeetsTbass_OnTheSurfaceAndW()
    {
        if (Tbass.Value is not { } r) return;
        JsonElement surface = r.GetProperty("surface");
        double[] rMean = GbassOracle.Numbers(surface.GetProperty("mean")), rLo = GbassOracle.Numbers(surface.GetProperty("q05")),
                 rHi = GbassOracle.Numbers(surface.GetProperty("q95"));
        JsonElement w = r.GetProperty("w");
        double rW = w.GetProperty("mean").GetDouble();
        output.WriteLine($"tbass: {r.GetProperty("seconds").GetDouble():F1} s; w {rW:G4} [{w.GetProperty("q05").GetDouble():G4}, {w.GetProperty("q95").GetDouble():G4}], " +
                         $"M {r.GetProperty("M").GetProperty("mean").GetDouble():F2}");
        foreach ((string name, BassCoefficientPrior form) in new[] { ("g-prior", BassCoefficientPrior.G()), ("ridge", BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma())) })
        {
            var e = Engine(form, 41);
            double worst = 0, sumSq = 0;
            for (int i = 0; i < rMean.Length; i++)
            {
                double half = Math.Max(0.5 * (rHi[i] - rLo[i]), 0.5 * (e.Hi[i] - e.Lo[i])), gap = Math.Abs(e.Mean[i] - rMean[i]) / half;
                worst = Math.Max(worst, gap);
                sumSq += gap * gap;
            }
            output.WriteLine($"engine, {name}: w {e.W:G4} [{e.WLo:G4}, {e.WHi:G4}], M {e.M:F2}; surface gap over the wider 90% half-band: " +
                             $"rms {Math.Sqrt(sumSq / rMean.Length):F3}, worst {worst:F3}, at {rMean.Length} grid points");
            Assert.True(worst <= 1.0, $"{name}: a surface mean lies outside the other's 90% band ({worst:F3})");
            Assert.InRange(e.W / rW, 0.8, 1.25);
        }
    }
}
