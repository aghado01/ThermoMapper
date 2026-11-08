using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// KHAOS slice 4's coefficient route (<see cref="BassSensitivity.FromCoefficients"/>) against <c>sobol_khaos</c> on the same draws:
/// an engine fit to Ishigami's function under KHAOS's defaults, its draws' multi-indices, coefficients and σ² packed as khaos 2.1.1's
/// <c>adaptive_khaos</c> object holds them (<c>khaos_sobol_oracle.R</c>). <c>sobol_khaos</c> divides by the explained variance plus σ²,
/// so each of its indices is the engine's main or total variance over Var f + σ², draw for draw. Opt-in: it runs where R is present.
/// </summary>
public sealed class KhaosSobolParityTests(ITestOutputHelper output)
{
    [Fact]
    public void SobolKhaos_OnTheEnginesDraws_IsTheCoefficientRoute()
    {
        if (!ROracle.IsAvailable) return;
        double[][] rows = KhaosReadoutParityTests.Data(120, 20261108);
        int p = rows[0].Length - 1;
        var domain = new InputDomain(0.0, 1.0);
        var data = new BassData([.. Enumerable.Range(0, p).Select(j => new ContinuousPredictor($"x{j}", [.. rows.Select(r => r[j])], domain))],
                                new ScalarResponse([.. rows.Select(r => r[p])]));
        BassConfig config = new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 3, Legendre = new LegendreBasis { MaxDegree = 8 } },
            Count = new BassCountPrior(4, 4.0 / rows.Length),
            Coefficients = BassCoefficientPrior.ModifiedG(TauPrior.Gamma(0.001, 1000.0)),
            Intercept = GPriorCentre.Zero,
            Sigma2 = Sigma2Prior.Jeffreys,
            Sampling = new BassSampling { Chains = 2, MasterSeed = 9, BurnIn = 2_000, Samples = 40, Thin = 25 },
        };
        BassFit fit = Bass.Fit(config, data);
        SobolDecomposition s = BassSensitivity.FromCoefficients(fit, data).Sobol.Decompositions[0];
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        Assert.Equal(pooled.Length, s.Variance.Length);

        // The draws as khaos holds them: coefficients (the intercept as term 0), each term's factors, and σ².
        var packed = new List<double[]>();
        for (int d = 0; d < pooled.Length; d++)
        {
            BassDraw draw = pooled[d];
            for (int j = 0; j < draw.Coefficients.Length; j++) packed.Add([1, d + 1, j, draw.Coefficients[j], 0]);
            for (int j = 0; j < draw.Structure.TermIds.Length; j++)
                foreach (LegendreFactor f in fit.Terms[draw.Structure.TermIds[j]].Legendre) packed.Add([2, d + 1, j + 1, f.Variable + 1, f.Degree]);
            packed.Add([3, d + 1, 0, draw.Sigma2, 0]);
        }
        JsonElement? result = GbassOracle.Run("khaos_sobol_oracle.R", packed, p.ToString());
        if (result is not { } r) return;
        string[] labels = [.. r.GetProperty("labels").EnumerateArray().Select(e => e.GetString()!)];
        JsonElement[] sRows = [.. r.GetProperty("S").EnumerateArray()], tRows = [.. r.GetProperty("T").EnumerateArray()];
        double[] leftover = GbassOracle.Numbers(r.GetProperty("leftover"));
        string Label(int[] u) => string.Join(":", u.Select(v => $"x{v + 1}"));
        Assert.All(labels, label => Assert.Contains(s.Subsets, u => Label(u) == label));

        double worst = 0.0;
        for (int d = 0; d < pooled.Length; d++)
        {
            double denominator = s.Variance[d] + pooled[d].Sigma2;
            double[] total = GbassOracle.Numbers(tRows[d]), main = GbassOracle.Numbers(sRows[d]);
            for (int v = 0; v < p; v++) worst = Math.Max(worst, Math.Abs(total[v] - s.TotalVariance[d, v] / denominator));
            for (int u = 0; u < s.Subsets.Count; u++)
            {
                int at = Array.IndexOf(labels, Label(s.Subsets[u]));
                worst = Math.Max(worst, Math.Abs((at < 0 ? 0.0 : main[at]) - s.MainVariance[d, u] / denominator));
            }
            worst = Math.Max(worst, Math.Abs(leftover[d] - pooled[d].Sigma2 / denominator));
        }
        output.WriteLine($"{pooled.Length} draws, {labels.Length} of sobol_khaos's labels, {s.Subsets.Count} engine subsets: largest index gap {worst:E2}");
        Assert.True(worst <= 1e-12, $"{worst:E2}");
    }
}
