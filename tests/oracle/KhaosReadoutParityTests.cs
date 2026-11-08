using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// KHAOS's posterior readouts (<c>khaos_readouts_oracle.R</c>; E19's build brief, slice 0): <c>adaptive_khaos</c> under the ridge and
/// the modified g-prior, and <c>sobol_khaos</c>, on Ishigami's function, a smooth simulator whose Sobol indices are known in closed
/// form, with one inert input. khaos's chain is not the engine's, so the engine meets these readouts from slice 1 on, never their
/// traces, each difference stated beside the comparison. Slice 0 checks that each mode runs from <c>r/</c> and reads what its quantities
/// mean. Under the ridge (degree 8, order 3, 6,000 iterations) the surface follows the function and the total indices are Ishigami's
/// to a few hundredths, the x₁:x₃ interaction included, on every seed tried. Under the g-prior as written the chain carries the defects
/// <see cref="KhaosDefectParityTests"/> confirms: it keeps about 70 terms on these 120 cases, a quarter of them repeating a variable, and
/// its indices move with the seed (x₁'s total 0.69 on one, 0.19 on another, at 4,000 iterations), so its readouts are checked for their
/// form only, never for their values. The corrected samplers (<c>corrected=1</c>), which reproduce the exact structure law
/// (<see cref="KhaosExactLawTests"/>), are the reference the engine meets: under both priors they read Ishigami's indices to a few
/// hundredths, with no term repeating a variable.
/// </summary>
public sealed class KhaosReadoutParityTests(ITestOutputHelper output)
{
    internal const int Inputs = 4, Observations = 120;
    internal const double NoiseSd = 0.5;

    /// <summary>Ishigami's function (a = 7, b = 0.1) on [0, 1]³ mapped to [−π, π]³.</summary>
    internal static double Ishigami(IReadOnlyList<double> x)
    {
        double z1 = Math.PI * (2 * x[0] - 1), z2 = Math.PI * (2 * x[1] - 1), z3 = Math.PI * (2 * x[2] - 1);
        return Math.Sin(z1) + 7 * Math.Pow(Math.Sin(z2), 2) + 0.1 * Math.Pow(z3, 4) * Math.Sin(z1);
    }

    /// <summary>Ishigami's total Sobol indices under independent uniform inputs: V₁ = ½(1 + bπ⁴/5)², V₂ = a²/8, V₁₃ = b²π⁸(1/18 − 1/50).</summary>
    internal static double[] IshigamiTotals()
    {
        double v1 = 0.5 * Math.Pow(1 + 0.1 * Math.Pow(Math.PI, 4) / 5, 2), v2 = 49.0 / 8, v13 = 0.01 * Math.Pow(Math.PI, 8) * (1.0 / 18 - 1.0 / 50);
        double v = v1 + v2 + v13;
        return [(v1 + v13) / v, v2 / v, v13 / v, 0.0];
    }

    /// <summary>n cases on [0, 1]^<see cref="Inputs"/>, Ishigami's function of the first three and noise 0.5·N(0, 1).</summary>
    internal static double[][] Data(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        return [.. Enumerable.Range(0, n).Select(_ =>
        {
            double[] x = [.. Enumerable.Range(0, Inputs).Select(_ => rng.NextDouble())];
            return x.Append(Ishigami(x) + NoiseSd * Normal.SampleStandard(rng)).ToArray();
        })];
    }

    /// <summary>The surface's grid: 16 points on [0, 1]^<see cref="Inputs"/>.</summary>
    internal static double[][] Grid()
    {
        var rng = new Xoshiro256PlusPlus(20261108);
        return [.. Enumerable.Range(0, 16).Select(_ => Enumerable.Range(0, Inputs).Select(_ => rng.NextDouble()).ToArray())];
    }

    internal static JsonElement? Readout(string prior, int nmcmc, int nburn, params string[] args)
    {
        if (!ROracle.IsAvailable) return null;
        string grid = GbassOracle.WriteCsv(Grid());
        try { return GbassOracle.Run("khaos_readouts_oracle.R", Data(Observations, 20261107), [prior, grid, nmcmc.ToString(), nburn.ToString(), "11", .. args]); }
        finally { System.IO.File.Delete(grid); }
    }

    private static readonly Dictionary<string, Lazy<JsonElement?>> Fits = new()
    {
        ["ridge"] = new(() => Readout("ridge", 6000, 3000, "degree=8", "order=3")),
        // The g-prior as written keeps about 70 terms here, so its chain, predict and sobol_khaos cost several times the ridge's: a
        // shorter chain, thinned, keeps the mode inside the bridge's 120 s on a loaded machine.
        ["gprior"] = new(() => Readout("gprior", 2500, 1250, "degree=8", "order=3", "thin=2")),
        // The corrected reference, khaos's samplers with slice 0's corrections (khaos_corrected), which reproduce the exact structure
        // law (KhaosExactLawTests): what the engine meets.
        ["ridge corrected"] = new(() => Readout("ridge", 6000, 3000, "degree=8", "order=3", "corrected=1")),
        ["gprior corrected"] = new(() => Readout("gprior", 6000, 3000, "degree=8", "order=3", "corrected=1")),
    };

    private static string PriorOf(string fit) => fit.Split(' ')[0];

    [Theory]
    [InlineData("ridge")]
    [InlineData("gprior")]
    [InlineData("ridge corrected")]
    [InlineData("gprior corrected")]
    public void EachPrior_RunsFromR_AndEmitsItsReadouts(string fit)
    {
        if (Fits[fit].Value is not { } r) return;
        string prior = PriorOf(fit);
        bool corrected = fit.EndsWith("corrected");
        Assert.Equal(prior, r.GetProperty("prior").GetString());
        Assert.Equal(corrected, r.GetProperty("corrected").GetBoolean());
        Assert.Equal(Observations, r.GetProperty("N").GetInt32());
        Assert.Equal(Inputs, r.GetProperty("p").GetInt32());
        JsonElement surface = r.GetProperty("surface");
        double[][] grid = Grid();
        foreach (string k in new[] { "mean", "q05", "q50", "q95" })
            Assert.Equal(grid.Length, GbassOracle.Numbers(surface.GetProperty(k)).Count(double.IsFinite));
        double[] mean = GbassOracle.Numbers(surface.GetProperty("mean")), lo = GbassOracle.Numbers(surface.GetProperty("q05")), hi = GbassOracle.Numbers(surface.GetProperty("q95"));
        Assert.All(Enumerable.Range(0, grid.Length), i => Assert.True(lo[i] <= hi[i]));
        double[] truth = [.. grid.Select(Ishigami)];
        double rmse = Math.Sqrt(truth.Zip(mean, (t, m) => (t - m) * (t - m)).Average());
        int covered = Enumerable.Range(0, grid.Length).Count(i => lo[i] <= truth[i] && truth[i] <= hi[i]);
        double[] proposed = GbassOracle.Numbers(r.GetProperty("count_propose")), accepted = GbassOracle.Numbers(r.GetProperty("count_accept"));
        Assert.True(accepted.Sum() > 0 && proposed.Zip(accepted).All(t => t.First >= t.Second));
        JsonElement terms = r.GetProperty("terms");
        output.WriteLine($"{fit}: {r.GetProperty("seconds").GetDouble():F1} s; M {r.GetProperty("M").GetProperty("mean").GetDouble():F2}, s2 {r.GetProperty("s2").GetProperty("mean").GetDouble():G4} " +
                         $"(noise {NoiseSd * NoiseSd}); the surface's RMSE against the function {rmse:F3} (its sd 3.72), 90% bands cover {covered} of {grid.Length}; " +
                         $"acceptance B/D/degree/swap {string.Join("/", accepted.Zip(proposed, (a, p) => p > 0 ? (a / p).ToString("F3") : "-"))}; " +
                         $"{terms.GetProperty("repeated").GetInt32()} of {terms.GetProperty("count").GetInt32()} kept terms repeat a variable");
        if (prior == "gprior")
            output.WriteLine($"g₀²: mean {r.GetProperty("g0sq").GetProperty("mean").GetDouble():G4}, 90% [{r.GetProperty("g0sq").GetProperty("q05").GetDouble():G4}, {r.GetProperty("g0sq").GetProperty("q95").GetDouble():G4}]");
        if (corrected) Assert.Equal(0, terms.GetProperty("repeated").GetInt32());
        if (fit != "gprior") Assert.True(rmse < 1.0, $"{fit}'s surface is {rmse:F2} from the function");
    }

    [Theory]
    [InlineData("ridge")]
    [InlineData("gprior")]
    [InlineData("ridge corrected")]
    [InlineData("gprior corrected")]
    public void SobolKhaos_ReadsTheIndices_AndIshigamisWhereTheChainHoldsTheLaw(string fit)
    {
        if (Fits[fit].Value is not { } r) return;
        JsonElement s = r.GetProperty("sobol");
        double[] total = GbassOracle.Numbers(s.GetProperty("T_f_mean")), truth = IshigamiTotals();
        double[] lo = GbassOracle.Numbers(s.GetProperty("T_q05")), hi = GbassOracle.Numbers(s.GetProperty("T_q95"));
        string[] labels = [.. s.GetProperty("labels").EnumerateArray().Select(e => e.GetString()!)];
        double[] first = GbassOracle.Numbers(s.GetProperty("S_mean"));
        output.WriteLine($"{fit}: total indices over Var f {string.Join(", ", total.Select((t, j) => $"x{j + 1} {t:F3} (true {truth[j]:F3})"))}; " +
                         $"sobol_khaos's own, over Var f + s2: {string.Join(", ", GbassOracle.Numbers(s.GetProperty("T_mean")).Select(t => t.ToString("F3")))}, " +
                         $"the noise's share {s.GetProperty("leftover_mean").GetDouble():F3}");
        output.WriteLine($"  partial indices: {string.Join(", ", labels.Zip(first).Where(t => t.Second > 0.005).Select(t => $"{t.First} {t.Second:F3}"))}");
        Assert.Equal(Inputs, total.Length);
        Assert.All(Enumerable.Range(0, Inputs), j => Assert.True(lo[j] <= hi[j] && total[j] >= 0));
        Assert.InRange(s.GetProperty("leftover_mean").GetDouble(), 0.0, 1.0);
        if (fit == "gprior") return;
        // Under the ridge, and under either corrected sampler: each total index within 0.05 of Ishigami's, the inert input's below 0.02.
        for (int j = 0; j < 3; j++) Assert.InRange(total[j], truth[j] - 0.05, truth[j] + 0.05);
        Assert.True(total[3] < 0.02, $"the inert input's total index {total[3]:F3}");
    }
}
