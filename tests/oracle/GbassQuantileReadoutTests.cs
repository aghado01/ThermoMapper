using System;
using System.Linq;
using System.Text.Json;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// Quantile BASS against GBASS's <c>qbass</c> (E1's slice 2 gate), posterior readouts on <see cref="GbassReadoutParityTests"/>' data and
/// grid at q = 0.9. Both hold β at 1/q − 1/(1 − q) and c at 2/(q(1 − q)), so σ² is <c>qbass</c>'s w on both sides, with w's prior
/// IG(0.1, 0.05), <c>qbass</c>'s. The differences stated: GBASS draws knots from U(0, 1) where the engine takes D01's grid, mutates a whole
/// term, weights orders and variables with pseudo-counts 1 against 5, holds a raw ridge with the intercept inside and τ ~ IG(½, N/2) where
/// the engine's default is D03's g-prior (the ridge read beside it), and samples λ; and its w step at β ≠ 0 keeps no floor, which the
/// engine has none of either. The priors read the sample 0.9-quantile on the engine's side, level and scale, where GBASS reads none.
/// <para>The ridge arm, GBASS's own form, holds the parity bar: each surface mean within the wider 90% half-band. The g-prior arm carries a
/// stated model difference, so its gap is not expected to be zero; its bar, 1.5 half-bands, catches gross defects and is no parity bar.
/// Over engine seeds 43 to 48 its worst gap read 0.49 to 1.23 (mean 0.85) with the GIG's tangent hat drawing v and 0.61 to 1.26 (mean
/// 0.92) with its ratio-of-uniforms draw, two of six above 1.0 under either: the 1.0 bar it held at first measured seed 43's stream more
/// than the sampler. The ridge arm read 0.51 to 0.81 on all six (set 2026-11-06, with the faster GIG draw).</para>
/// </summary>
public sealed class GbassQuantileReadoutTests(ITestOutputHelper output)
{
    private const double Level = 0.9;

    private static readonly Lazy<JsonElement?> Qbass = new(() => GbassReadoutParityTests.Readout("qbass", 3000, 1500, "q=0.9"));

    private static BassData Data()
    {
        double[][] rows = GbassReadoutParityTests.Data();
        return new BassData([new ContinuousPredictor("x1", [.. rows.Select(r => r[0])]), new ContinuousPredictor("x2", [.. rows.Select(r => r[1])])],
                            new ScalarResponse([.. rows.Select(r => r[2])]));
    }

    [Fact]
    public void QuantileBass_MeetsQbass_OnTheSurfaceAndW()
    {
        if (Qbass.Value is not { } r) return;
        JsonElement surface = r.GetProperty("surface");
        double[] rMean = GbassOracle.Numbers(surface.GetProperty("mean")), rLo = GbassOracle.Numbers(surface.GetProperty("q05")),
                 rHi = GbassOracle.Numbers(surface.GetProperty("q95"));
        JsonElement w = r.GetProperty("w");
        double rW = w.GetProperty("mean").GetDouble();
        output.WriteLine($"qbass: {r.GetProperty("seconds").GetDouble():F1} s; w {rW:G4} [{w.GetProperty("q05").GetDouble():G4}, {w.GetProperty("q95").GetDouble():G4}], " +
                         $"M {r.GetProperty("M").GetProperty("mean").GetDouble():F2}, share below its fit {r.GetProperty("below_fit").GetDouble():F3}");
        BassData data = Data();
        double[][] grid = GbassReadoutParityTests.Grid();
        var cases = new PredictionCases([new ContinuousColumn("x1", [.. grid.Select(g => g[0])]), new ContinuousColumn("x2", [.. grid.Select(g => g[1])])]);
        foreach ((string name, BassCoefficientPrior form) in new[] { ("g-prior", BassCoefficientPrior.G()), ("ridge", BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma())) })
        {
            var config = new BassConfig
            {
                Structure = new BassStructureSpec { MaxInteraction = 2 },
                Coefficients = form,
                Sigma2 = Sigma2Prior.InverseGamma(0.1, 0.05),
                Noise = BassNoise.AsymmetricLaplace(Level),
                Sampling = new BassSampling { Chains = 4, MasterSeed = 43, BurnIn = 3000, Samples = 1500, Thin = 2 },
            };
            BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
            BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Predictive = false });
            double[,] mean = p.Mean();
            double worst = 0, sumSq = 0;
            for (int i = 0; i < grid.Length; i++)
            {
                double[] f = [.. Enumerable.Range(0, p.DrawCount).Select(d => p.Latent[d, i, 0]).Order()];
                double half = Math.Max(0.5 * (rHi[i] - rLo[i]), 0.5 * (f[(int)(0.95 * f.Length) - 1] - f[(int)(0.05 * f.Length)]));
                double gap = Math.Abs(mean[i, 0] - rMean[i]) / half;
                worst = Math.Max(worst, gap);
                sumSq += gap * gap;
            }
            double[] sigma2 = [.. fit.Draws.SelectMany(c => c).Select(d => d.Sigma2).Order()];
            double m = fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length);
            output.WriteLine($"engine, {name}: w {sigma2.Average():G4} [{sigma2[(int)(0.05 * sigma2.Length)]:G4}, {sigma2[(int)(0.95 * sigma2.Length) - 1]:G4}], M {m:F2}; " +
                             $"surface gap over the wider 90% half-band: rms {Math.Sqrt(sumSq / grid.Length):F3}, worst {worst:F3}");
            double bar = name == "ridge" ? 1.0 : 1.5;
            Assert.True(worst <= bar, $"{name}: a surface mean lies {worst:F3} of the wider 90% half-band from the other's, past {bar}");
            Assert.InRange(sigma2.Average() / rW, 0.8, 1.25);
        }
    }
}
