using System;
using System.Linq;
using System.Text.Json;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// Horseshoe BASS against GBASS's <c>hbass</c> (E1's slice 4 gate), posterior readouts on <see cref="GbassReadoutParityTests"/>' data and grid.
/// <c>hbass</c> takes v ~ GBP(1, ½, ½), c = var(y) and √w ~ C⁺(0, 1); the engine takes the same v, c = 1 and σ ~ C⁺(0, s) at the
/// reference's mean square s² by the hierarchical route (<see cref="Sigma2Prior.HalfCauchy"/>), so its σ² is <c>hbass</c>'s w·var(y), the
/// noise's scale per unit v on both sides, up to the reference's divisor N where var(y) takes N − 1. <c>hbass</c> as written adapts v's
/// random-walk step from mean(ln v)'s spread (slice 0's defect), which here shrinks it to 0.12 at an acceptance of 0.98: v mixes slowly, and
/// against the same sampler with v's step held at 2.5 and its adaptation off, the correction and the reference, it reads M 4.16 against 5.36
/// and w·var(y) 0.082 against 0.068. It is read for the record. The engine moves v by a slice step on ln v. The chains differ by design otherwise as
/// under <c>tbass</c> (<see cref="GbassStudentTReadoutTests"/>): knots, mutation, pseudo-counts, the ridge's intercept and λ. The bars are
/// <c>qbass</c>'s (<see cref="GbassQuantileReadoutTests"/>): the ridge, GBASS's own form, holds parity, the surface within the wider 90% band
/// (0.52 to 0.82 over engine seeds 41 to 46); the g-prior, a stated model difference, a bar for gross defects at 1.5 (0.55 to 1.20). σ²
/// meets w·var(y) within a fifth (1.10 to 1.24 of it over those seeds), and v's posterior meets the reference's: the mean of ln v within 0.1
/// (−0.21 to −0.16 against −0.16) and the share of v below 0.1 within 0.03 (0.159 to 0.171 against 0.172).
/// </summary>
public sealed class GbassHorseshoeReadoutTests(ITestOutputHelper output)
{
    private static readonly Lazy<JsonElement?> AsWritten = new(() => GbassReadoutParityTests.Readout("hbass", 6000, 3000));
    private static readonly Lazy<JsonElement?> Corrected = new(() => GbassReadoutParityTests.Readout("hbass", 6000, 3000, "prop_sigma=2.5"));

    private static BassData Data()
    {
        double[][] rows = GbassReadoutParityTests.Data();
        return new BassData([new ContinuousPredictor("x1", [.. rows.Select(r => r[0])]), new ContinuousPredictor("x2", [.. rows.Select(r => r[1])])],
                            new ScalarResponse([.. rows.Select(r => r[2])]));
    }

    internal sealed record EngineReadout(double[] Mean, double[] Lo, double[] Hi, double S2, double S2Lo, double S2Hi, double M, double LnV, double VSmall);

    internal static EngineReadout Engine(BassCoefficientPrior form, int seed)
    {
        BassData data = Data();
        double[][] grid = GbassReadoutParityTests.Grid();
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2 },
            Coefficients = form,
            Noise = BassNoise.Horseshoe(),
            Sigma2 = Sigma2Prior.HalfCauchy(),
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
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        double[] s2 = [.. draws.Select(d => d.Sigma2).Order()];
        double[][] v = [.. draws.Select(d => d.Latents!.Values.ToArray())];
        double lnv = v.Average(x => x.Average(Math.Log)), small = v.Average(x => x.Average(vi => vi < 0.1 ? 1.0 : 0.0));
        return new EngineReadout(mean, lo, hi, s2.Average(), s2[(int)(0.05 * s2.Length)], s2[(int)(0.95 * s2.Length) - 1],
                                 draws.Average(d => d.Structure.TermIds.Length), lnv, small);
    }

    // The worst gap of the engine's surface means from the reference's over the wider 90% half-band, and its rms.
    internal static (double Worst, double Rms) Gap(JsonElement r, EngineReadout e)
    {
        JsonElement surface = r.GetProperty("surface");
        double[] rMean = GbassOracle.Numbers(surface.GetProperty("mean")), rLo = GbassOracle.Numbers(surface.GetProperty("q05")),
                 rHi = GbassOracle.Numbers(surface.GetProperty("q95"));
        double worst = 0, sumSq = 0;
        for (int i = 0; i < rMean.Length; i++)
        {
            double half = Math.Max(0.5 * (rHi[i] - rLo[i]), 0.5 * (e.Hi[i] - e.Lo[i])), gap = Math.Abs(e.Mean[i] - rMean[i]) / half;
            worst = Math.Max(worst, gap);
            sumSq += gap * gap;
        }
        return (worst, Math.Sqrt(sumSq / rMean.Length));
    }

    internal static double VarY()
    {
        double[] y = [.. GbassReadoutParityTests.Data().Select(row => row[2])];
        double ybar = y.Average();
        return y.Sum(v => (v - ybar) * (v - ybar)) / (y.Length - 1);
    }

    private string Describe(string name, JsonElement r)
    {
        double varY = VarY();
        JsonElement w = r.GetProperty("w");
        return $"{name}: {r.GetProperty("seconds").GetDouble():F1} s; v's step {r.GetProperty("v_step").GetDouble():G3}, accepted {r.GetProperty("v_accept").GetDouble():F3}; " +
               $"w·var(y) {w.GetProperty("mean").GetDouble() * varY:G4} [{w.GetProperty("q05").GetDouble() * varY:G4}, {w.GetProperty("q95").GetDouble() * varY:G4}]; " +
               $"mean ln v {r.GetProperty("lnv").GetProperty("mean").GetDouble():F3}, v below 0.1 {r.GetProperty("v_small").GetDouble():F3}; " +
               $"M {r.GetProperty("M").GetProperty("mean").GetDouble():F2}";
    }

    [Fact]
    public void HBass_MeetsCorrectedHbass_OnTheSurfaceTheNoiseScaleAndV()
    {
        if (Corrected.Value is not { } r || AsWritten.Value is not { } plain) return;
        output.WriteLine(Describe("hbass as written", plain));
        output.WriteLine(Describe("hbass, v's step held at 2.5", r));
        double rS2 = r.GetProperty("w").GetProperty("mean").GetDouble() * VarY(), rLnV = r.GetProperty("lnv").GetProperty("mean").GetDouble(),
               rSmall = r.GetProperty("v_small").GetDouble();
        foreach ((string name, BassCoefficientPrior form, double bar) in new[]
                 { ("g-prior", BassCoefficientPrior.G(), 1.5), ("ridge", BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()), 1.0) })
        {
            EngineReadout e = Engine(form, 41);
            (double worst, double rms) = Gap(r, e);
            output.WriteLine($"engine, {name}: σ² {e.S2:G4} [{e.S2Lo:G4}, {e.S2Hi:G4}], mean ln v {e.LnV:F3}, v below 0.1 {e.VSmall:F3}, M {e.M:F2}; " +
                             $"surface gap over the wider 90% half-band: rms {rms:F3}, worst {worst:F3}");
            Assert.True(worst <= bar, $"{name}: a surface mean lies {worst:F3} of the wider half-band from the other's, past {bar}");
            Assert.InRange(e.S2 / rS2, 0.8, 1.25);
            Assert.True(Math.Abs(e.LnV - rLnV) <= 0.1, $"{name}: mean ln v {e.LnV:F3} against {rLnV:F3}");
            Assert.True(Math.Abs(e.VSmall - rSmall) <= 0.03, $"{name}: v below 0.1 {e.VSmall:F3} against {rSmall:F3}");
        }
    }
}
