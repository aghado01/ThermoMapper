using System;
using System.Globalization;
using System.IO;
using System.Text;
using System.Text.Json;
using Maths.Distributions;
using Maths.Regression;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// The rank-normalized diagnostics and PSIS-LOO against R (OQ9, N13): <see cref="ChainDiagnostics"/> against the
/// posterior package's <c>rhat</c>, <c>ess_bulk</c>, <c>ess_tail</c> and <c>ess_mean</c>, and
/// <see cref="PredictiveChecks.PsisLoo"/> against the loo package's k̂, <c>relative_eff</c> and Monte Carlo errors, on the
/// same fixtures (<c>oracles/diagnostics_oracle.R</c>). Both sides implement the same algorithms, so they agree to
/// rounding, at odd chain lengths too: the tail diagnostics take their median and quantiles over all draws, as posterior
/// does. One convention differs on purpose, and is asserted as such: a PSIS tail whose largest ratios tie is k̂ = −∞ here
/// (<see cref="ParetoTail.Tied"/>) and +∞ in loo, which then withholds the total's Monte Carlo error.
/// </summary>
public sealed class DiagnosticsParityTests
{
    private const double Relative = 1e-9;
    private readonly ITestOutputHelper _out;
    public DiagnosticsParityTests(ITestOutputHelper output) => _out = output;

    public static TheoryData<string> ChainFixtures => new() { "ar1", "antithetic", "scales", "shift", "cauchy", "odd" };

    [Theory]
    [MemberData(nameof(ChainFixtures))]
    public void RankNormalizedDiagnostics_MatchPosterior(string name)
    {
        if (!ROracle.IsAvailable) return;
        double[][] chains = Chains(name);
        string csv = WriteCsv(chains[0].Length, chains.Length, (i, c) => chains[c][i]);
        try
        {
            JsonElement r = ROracle.Run("oracles/diagnostics_oracle.R", csv, "chains");
            ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(chains);
            _out.WriteLine($"{name}: R̂ {d.RHat.Value:R} (bulk {d.BulkRHat.Value:R}, tail {d.TailRHat.Value:R}), bulk ESS {d.BulkEss.Value:R}, tail ESS {d.TailEss.Value:R}");
            _out.WriteLine($"  posterior: rhat_tail {Number(r, "rhat_tail"):R}, ess_tail {Number(r, "ess_tail"):R}");

            Close(Number(r, "rhat"), d.RHat.Value, "R̂");
            Close(Number(r, "rhat_bulk"), d.BulkRHat.Value, "bulk R̂");
            Close(Number(r, "rhat_tail"), d.TailRHat.Value, "tail R̂");
            Close(Number(r, "ess_bulk"), d.BulkEss.Value, "bulk ESS");
            Close(Number(r, "ess_tail"), d.TailEss.Value, "tail ESS");
            Close(Number(r, "ess_mean"), ChainDiagnostics.MeanEss(chains).Value, "ESS of the mean");
        }
        finally
        {
            File.Delete(csv);
        }
    }

    public static TheoryData<string> LooFixtures => new() { "independent", "chains", "influential", "tied", "nofit" };

    [Theory]
    [MemberData(nameof(LooFixtures))]
    public void PsisLoo_MatchesLoo(string name)
    {
        if (!ROracle.IsAvailable) return;
        (double[,] ll, int chains) = LogLikelihood(name);
        int s = ll.GetLength(0), units = ll.GetLength(1);
        string csv = WriteCsv(s, units, (d, u) => ll[d, u]);
        try
        {
            JsonElement r = ROracle.Run("oracles/diagnostics_oracle.R", csv, "loo", chains.ToString(CultureInfo.InvariantCulture));
            int[]? lengths = chains > 0 ? Repeat(s / chains, chains) : null;
            var density = new PointwiseLogDensity(ll, UnitPartition.EachObservation(units), Array.Empty<LatentReport>(), lengths);
            PredictiveCheck c = PredictiveChecks.PsisLoo(density);
            double[] k = Numbers(r, "pareto_k");
            _out.WriteLine($"{name}: elpd {c.Elpd:R} vs {Number(r, "elpd"):R}; mcse {c.MonteCarloError:R} vs {Number(r, "mcse"):R}; largest k̂ {Max(c.ParetoK!):R} vs {Max(k):R} (threshold {Number(r, "threshold"):R})");

            Close(Number(r, "threshold"), PredictiveChecks.ParetoKThreshold(s), "k̂ threshold");
            Close(Number(r, "elpd"), c.Elpd, "elpd");
            Close(Number(r, "se"), c.StandardError, "elpd's standard error");
            Close(Number(r, "p_loo"), c.EffectiveParameters, "p_loo");
            double[] elpd = Numbers(r, "elpd_i"), mcse = Numbers(r, "mcse_i"), rEff = Numbers(r, "r_eff");
            for (int u = 0; u < units; u++)
            {
                Close(elpd[u], c.PerUnit[u], $"unit {u}'s elpd");
                Close(mcse[u], c.MonteCarloErrors![u], $"unit {u}'s Monte Carlo error");
                Close(rEff[u], c.RelativeEfficiency![u], $"unit {u}'s r_eff");
                if (c.ParetoTails![u] == ParetoTail.Tied)
                {
                    // The convention that differs: loo cannot fit a tied tail and calls it +∞.
                    Assert.Equal(double.PositiveInfinity, k[u]);
                    Assert.Equal(double.NegativeInfinity, c.ParetoK![u]);
                }
                else Close(k[u], c.ParetoK![u], $"unit {u}'s k̂", absolute: 1e-9);
            }
            ParetoTail[] expected = name switch
            {
                "tied" => [ParetoTail.Tied, ParetoTail.Fitted],
                "nofit" => [ParetoTail.NoFit, ParetoTail.Fitted],
                _ => Array.ConvertAll(new int[units], _ => ParetoTail.Fitted),
            };
            Assert.Equal(expected, c.ParetoTails);
            double total = Number(r, "mcse");
            if (name == "tied")
            {
                Assert.True(double.IsNaN(total), "loo withholds the total's error, its k̂ being +∞");
                Assert.True(double.IsFinite(c.MonteCarloError), "a tied tail counts as reliable here");
            }
            else if (double.IsNaN(total)) Assert.True(double.IsNaN(c.MonteCarloError), "loo withholds the total's error above the threshold");
            else Close(total, c.MonteCarloError, "the total's Monte Carlo error");
            if (name is "influential" or "nofit") Assert.True(Max(k) > Number(r, "threshold"), "a unit's k̂ should exceed the threshold");
            else if (name != "tied") Assert.True(double.IsFinite(total), "the total's Monte Carlo error should be compared");
            if (name == "chains") Assert.True(Array.Exists(rEff, v => v < 1.0), "some unit's r_eff should be below 1");
        }
        finally
        {
            File.Delete(csv);
        }
    }

    // Four chains: AR(1) at φ = 0.9 (ESS ≈ S/19) and φ = −0.3 (antithetic, ESS above S), one chain at three times the
    // others' scale, one shifted by half a standard deviation, Cauchy draws, and an odd length.
    private static double[][] Chains(string name)
    {
        var rng = new Xoshiro256PlusPlus(Seed(name));
        int length = name == "odd" ? 999 : 1000;
        var chains = new double[4][];
        for (int c = 0; c < 4; c++)
        {
            chains[c] = new double[length];
            double phi = name switch { "ar1" or "odd" => 0.9, "antithetic" => -0.3, _ => 0.0 };
            double x = Normal.SampleStandard(rng);
            for (int i = 0; i < length; i++)
            {
                double z = Normal.SampleStandard(rng);
                x = name == "cauchy" ? z / Normal.SampleStandard(rng) : phi * x + Math.Sqrt(1.0 - phi * phi) * z;
                chains[c][i] = name switch
                {
                    "scales" when c == 3 => 3.0 * x,
                    "shift" when c == 3 => x + 0.5,
                    _ => x,
                };
            }
        }
        return chains;
    }

    // y_u ~ N(0, 1) for 29 units and a 30th at 3.5 (2.5 for the chains, 9 for the influential case), with draws of (μ, σ) from their posterior under a flat
    // prior on (μ, log σ): independent, or four AR(1) chains in μ at ρ = 0.8, laid out chain by chain.
    private static (double[,] LogLikelihood, int Chains) LogLikelihood(string name)
    {
        if (name is "tied" or "nofit") return (TwoUnits(name), 0);
        var rng = new Xoshiro256PlusPlus(Seed(name));
        const int units = 30, draws = 4000;
        var y = new double[units];
        for (int u = 0; u < units - 1; u++) y[u] = Normal.SampleStandard(rng);
        y[units - 1] = name switch { "influential" => 9.0, "chains" => 2.5, _ => 3.5 };
        double mean = 0.0, ss = 0.0;
        foreach (double v in y) mean += v / units;
        foreach (double v in y) ss += (v - mean) * (v - mean);

        int chains = name == "chains" ? 4 : 0;
        var ll = new double[draws, units];
        double z = Normal.SampleStandard(rng);
        for (int d = 0; d < draws; d++)
        {
            double chi2 = 0.0;
            for (int j = 0; j < units - 1; j++) { double g = Normal.SampleStandard(rng); chi2 += g * g; }
            double sigma = Math.Sqrt(ss / chi2);
            if (chains > 0 && d % (draws / chains) == 0) z = Normal.SampleStandard(rng);
            z = chains > 0 ? 0.8 * z + 0.6 * Normal.SampleStandard(rng) : Normal.SampleStandard(rng);
            double mu = mean + sigma / Math.Sqrt(units) * z;
            for (int u = 0; u < units; u++)
                ll[d, u] = -0.5 * Math.Log(2.0 * Math.PI * sigma * sigma) - (y[u] - mu) * (y[u] - mu) / (2.0 * sigma * sigma);
        }
        return (ll, chains);
    }

    // Two units over 400 draws, whose tail of M = 60 log ratios the Pareto fit cannot read; the second unit varies smoothly.
    // "tied": three values, as a collapsed predictive over three structures takes; the largest ratio, 1/p = e², is held by
    // 200 draws, so the whole tail ties. "nofit": 355 draws share ratio e and 45 exceed it, so the tail's lowest 15 sit at
    // the cutoff and its first quartile is its minimum (posterior's gpdfit guard).
    private static double[,] TwoUnits(string name)
    {
        var rng = new Xoshiro256PlusPlus(Seed(name));
        var ll = new double[400, 2];
        for (int d = 0; d < 400; d++)
        {
            ll[d, 0] = name == "tied"
                ? d < 200 ? -2.0 : d < 300 ? -1.0 : -0.5
                : d < 355 ? -1.0 : -1.0 - 0.02 * (d - 354);
            double mu = 0.1 * Normal.SampleStandard(rng);
            ll[d, 1] = -0.5 * Math.Log(2.0 * Math.PI) - 0.5 * (0.3 - mu) * (0.3 - mu);
        }
        return ll;
    }

    private static int Seed(string name)
    {
        int h = 17;
        foreach (char ch in name) h = unchecked(h * 31 + ch);
        return h & 0x7fffffff;
    }

    private static int[] Repeat(int value, int count)
    {
        var a = new int[count];
        Array.Fill(a, value);
        return a;
    }

    private static double Max(double[] values)
    {
        double m = double.NegativeInfinity;
        foreach (double v in values) m = Math.Max(m, v);
        return m;
    }

    private static void Close(double expected, double actual, string what, double absolute = 0.0)
        => Assert.True(double.IsInfinity(expected)
                           ? actual == expected
                           : Math.Abs(actual - expected) <= Math.Max(absolute, Relative * Math.Abs(expected)),
                       $"{what}: C# {actual:R} vs R {expected:R}");

    // A JSON number; the oracle writes non-finite values as "Inf", "-Inf", "NaN" and "NA" (NA is NaN here).
    private static double Number(JsonElement root, string name) => Value(root.GetProperty(name));

    private static double[] Numbers(JsonElement root, string name)
    {
        JsonElement e = root.GetProperty(name);
        var a = new double[e.GetArrayLength()];
        int i = 0;
        foreach (JsonElement v in e.EnumerateArray()) a[i++] = Value(v);
        return a;
    }

    private static double Value(JsonElement e) => e.ValueKind switch
    {
        JsonValueKind.Number => e.GetDouble(),
        JsonValueKind.String => e.GetString() switch
        {
            "Inf" => double.PositiveInfinity,
            "-Inf" => double.NegativeInfinity,
            "NaN" or "NA" => double.NaN,
            var other => throw new FormatException($"Unexpected oracle value '{other}'."),
        },
        _ => throw new FormatException($"Unexpected oracle value kind {e.ValueKind}."),
    };

    private static string WriteCsv(int rows, int columns, Func<int, int, double> value)
    {
        string path = Path.Combine(Path.GetTempPath(), $"diagnostics_{Guid.NewGuid():N}.csv");
        var sb = new StringBuilder();
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                if (j > 0) sb.Append(',');
                sb.Append(value(i, j).ToString("R", CultureInfo.InvariantCulture));
            }
            sb.Append('\n');
        }
        File.WriteAllText(path, sb.ToString());
        return path;
    }
}
