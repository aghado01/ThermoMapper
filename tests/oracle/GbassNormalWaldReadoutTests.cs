using System;
using System.Linq;
using System.Text.Json;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// Normal-Wald BASS against GBASS's <c>nwbass</c> (E1's slice 3 gate), posterior readouts on <see cref="GbassReadoutParityTests"/>' data and
/// grid. Both draw the skew β ~ N(0, 1) and γ under N(2, 1), with w's prior IG(0.1, 0.05), <c>nwbass</c>'s, and c = 1, and R runs with
/// <c>lag_beta</c> = 0, so neither holds the skew for a burn-in. The differences stated: GBASS draws γ from its conditional's normal on the
/// whole line, where the engine takes its part on γ &gt; 0 (E1's slice 0); and as for <c>qbass</c>, GBASS's knots from U(0, 1) against D01's
/// grid, its whole-term mutation, pseudo-counts 1 against 5, its raw ridge with the intercept inside, λ sampled, and the engine's priors
/// reading the response's moments where GBASS reads none.
/// </summary>
public sealed class GbassNormalWaldReadoutTests(ITestOutputHelper output)
{
    private static readonly Lazy<JsonElement?> Nwbass =
        new(() => GbassReadoutParityTests.Readout("nwbass", 3000, 1500, "m_beta=0", "s_beta=1", "m_gamma=2", "s_gamma=1", "lag_beta=0"));

    private static BassData Data()
    {
        double[][] rows = GbassReadoutParityTests.Data();
        return new BassData([new ContinuousPredictor("x1", [.. rows.Select(r => r[0])]), new ContinuousPredictor("x2", [.. rows.Select(r => r[1])])],
                            new ScalarResponse([.. rows.Select(r => r[2])]));
    }

    [Fact]
    public void NormalWaldBass_MeetsNwbass_OnTheSurfaceWBetaAndGamma()
    {
        if (Nwbass.Value is not { } r) return;
        JsonElement surface = r.GetProperty("surface");
        double[] rMean = GbassOracle.Numbers(surface.GetProperty("mean")), rLo = GbassOracle.Numbers(surface.GetProperty("q05")),
                 rHi = GbassOracle.Numbers(surface.GetProperty("q95"));
        (double Mean, double Lo, double Hi) Summary(string key)
        {
            JsonElement s = r.GetProperty(key);
            return (s.GetProperty("mean").GetDouble(), s.GetProperty("q05").GetDouble(), s.GetProperty("q95").GetDouble());
        }
        var (rW, rWLo, rWHi) = Summary("w");
        var (rNoise, _, _) = Summary("s2");
        var (rBeta, rBetaLo, rBetaHi) = Summary("beta");
        var (rGamma, rGammaLo, rGammaHi) = Summary("gamma");
        output.WriteLine($"nwbass: {r.GetProperty("seconds").GetDouble():F1} s; w {rW:G4} [{rWLo:G4}, {rWHi:G4}], β {rBeta:G4} [{rBetaLo:G4}, {rBetaHi:G4}], " +
                         $"γ {rGamma:G4} [{rGammaLo:G4}, {rGammaHi:G4}], M {r.GetProperty("M").GetProperty("mean").GetDouble():F2}");
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
                Noise = BassNoise.NormalWald(new BassNoisePrior(2.0, 1.0), new BassNoisePrior(0.0, 1.0)),
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
            var draws = fit.Draws.SelectMany(c => c).ToList();
            (double Mean, double Lo, double Hi) Band(Func<BassDraw, double> read)
            {
                double[] s = [.. draws.Select(read).Order()];
                return (s.Average(), s[(int)(0.05 * s.Length)], s[(int)(0.95 * s.Length) - 1]);
            }
            var w = Band(d => d.Sigma2);
            var beta = Band(d => d.Latents!.Skew);
            var gamma = Band(d => d.Latents!.Gamma!.Value);
            var law = new VarianceMeanMixture(config.Noise!);
            var noise = Band(d => law.Variance(d.Sigma2, d.Latents));
            double m = draws.Average(d => d.Structure.TermIds.Length);
            output.WriteLine($"engine, {name}: w {w.Mean:G4} [{w.Lo:G4}, {w.Hi:G4}], β {beta.Mean:G4} [{beta.Lo:G4}, {beta.Hi:G4}], γ {gamma.Mean:G4} " +
                             $"[{gamma.Lo:G4}, {gamma.Hi:G4}], noise variance {noise.Mean:G4} against nwbass's {rNoise:G4}, M {m:F2}; surface gap over the wider " +
                             $"90% half-band: rms {Math.Sqrt(sumSq / grid.Length):F3}, worst {worst:F3}");
            Assert.True(worst <= 1.0, $"{name}: a surface mean lies outside the other's 90% band ({worst:F3})");
            // The noise's variance w(E[v] + β²Var[v]) is what the data fix; w and γ trade along it, so each is read for overlap: each side's
            // mean inside the other's 90% interval, as β's.
            Assert.InRange(noise.Mean / rNoise, 0.8, 1.25);
            Assert.InRange(w.Mean, rWLo, rWHi);
            Assert.InRange(rW, w.Lo, w.Hi);
            Assert.InRange(beta.Mean, rBetaLo, rBetaHi);
            Assert.InRange(rBeta, beta.Lo, beta.Hi);
            Assert.InRange(gamma.Mean, rGammaLo, rGammaHi);
            Assert.InRange(rGamma, gamma.Lo, gamma.Hi);
        }
    }
}
