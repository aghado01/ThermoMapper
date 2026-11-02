using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// The distributions' quantiles against base R (<c>r/oracles/distribution_quantile_oracle.R</c>), and truncated draws
/// against exact truncated moments in R's log tails. Each quantile round-trips through R's <c>p*</c>, which shares nothing
/// with it. R's own <c>q*</c> is the target only where it passes the same round trip; where it does not, the count is
/// printed for the record, as R's qt stops short in the far tail. Opt-in: it runs where the R toolchain is present.
/// </summary>
public sealed class DistributionQuantileParityTests
{
    private const double Ln2 = 0.69314718055994530942;
    private const int Draws = 400_000;
    private static readonly double[] LogTails = [-Ln2, -1.0, Math.Log(0.1), Math.Log(1e-3), Math.Log(1e-10), Math.Log(1e-100), -700.0];
    private readonly ITestOutputHelper _out;

    public DistributionQuantileParityTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Quantiles_MatchBaseR_WhereRSolvesItsCase()
    {
        if (!ROracle.IsAvailable) return;
        var cases = new List<double[]>();
        foreach (double shape in new[] { 1e-3, 0.01, 0.1, 0.5, 1.0, 2.5, 10.0, 100.0, 1000.0, 1e5 }) Add(cases, 1, shape, 0.0);
        foreach ((double a, double b) in new[] { (0.01, 0.01), (0.1, 5.0), (5.0, 0.1), (0.5, 0.5), (1.0, 3.0), (3.0, 1.0), (2.0, 5.0), (5.0, 2.0), (20.0, 30.0), (100.0, 1000.0), (1000.0, 100.0), (1000.0, 1000.0), (1e4, 100.0) })
            Add(cases, 2, a, b);
        foreach (double nu in new[] { 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 1000.0, 1e5 }) Add(cases, 3, nu, 0.0);
        Add(cases, 4, 0.0, 0.0);
        foreach (double tail in new[] { -1000.0, -1e4, -1e6 })
        {
            cases.Add([4.0, 0.0, 0.0, tail, 1.0, 0.0]);
            cases.Add([4.0, 0.0, 0.0, tail, 0.0, 0.0]);
        }

        foreach (double[] row in cases) row[5] = Ours(row[0], row[1], row[2], row[3], row[4] == 1.0);

        JsonElement r = Run(cases, "quantile");
        double[] rq = Doubles(r, "q"), back = Doubles(r, "r_log_tail"), atOurs = Doubles(r, "x_log_tail");
        double worstGap = 0.0, worstTrip = 0.0;
        int gaps = 0, trips = 0, unsolved = 0;
        string worstCase = "", worstTripCase = "";
        for (int i = 0; i < cases.Count; i++)
        {
            (double family, double a, double b, double tail, bool lower, double ours) = (cases[i][0], cases[i][1], cases[i][2], cases[i][3], cases[i][4] == 1.0, cases[i][5]);
            string what = $"family {family} ({a}, {b}) tail {tail} lower {lower}: {ours:R}, R's {rq[i]:R}, R's tail there {atOurs[i]:R}, R's at its own {back[i]:R}";
            double tolerance = 1e-12 * Math.Max(1.0, Math.Abs(tail)) + Prefactor(family, a, b, ours);

            // The round trip: an x past the double's reach of an end, or rounded onto 1, no longer holds the tail.
            bool holdsTail = double.IsFinite(ours) && Math.Abs(ours) > 1e-290 && !(family == 2.0 && 1.0 - ours < 1e-3);
            if (holdsTail)
            {
                double trip = Math.Abs(atOurs[i] - tail);
                Assert.True(trip <= tolerance, $"{what}: tolerance {tolerance:E1}");
                trips++;
                if (trip / tolerance > worstTrip)
                {
                    worstTrip = trip / tolerance;
                    worstTripCase = what;
                }
            }

            // R's root rounded onto 1 is no fault of R's; a NaN, an overflow, a floor at 2⁻¹⁰²³ or a far-tail t that stopped short is.
            bool rHoldsTail = !(family == 2.0 && 1.0 - rq[i] < 1e-3);
            if (double.IsNaN(rq[i]) || (rq[i] != ours && rHoldsTail && !(Math.Abs(back[i] - tail) <= tolerance)))
            {
                unsolved++;
                _out.WriteLine($"R short of its own round trip: {what}");
                continue;
            }
            // Where both saturate or underflow, the value is the end; at the centre of a symmetric law R's root is a rounding
            // error about 0; otherwise relative, with the flat conditioning a small shape gives x.
            if (ours == rq[i] || (Math.Abs(rq[i]) < 1e-290 && Math.Abs(ours) < 1e-290)) continue;
            if (family >= 3.0 && tail == -Ln2 && Math.Abs(ours) < 1e-14 && Math.Abs(rq[i]) < 1e-14) continue;
            double relative = Math.Abs(ours - rq[i]) / Math.Abs(rq[i]);
            double flat = family == 2.0 ? Math.Min(1.0, Math.Min(a, b)) : family == 4.0 ? 1.0 : Math.Min(1.0, a);
            Assert.True(relative <= 2e-10 / flat, what);
            gaps++;
            if (relative > worstGap)
            {
                worstGap = relative;
                worstCase = what;
            }
        }
        _out.WriteLine($"{cases.Count} cases: {trips} round trips through R's tails, worst at {worstTrip:P0} of tolerance ({worstTripCase})");
        _out.WriteLine($"{gaps} gaps from R's own quantile, worst {worstGap:E1} ({worstCase})");
        _out.WriteLine($"{unsolved} cases R's quantile fails its own round trip");
        Assert.True(trips >= 0.8 * cases.Count, $"only {trips} of {cases.Count} cases held their tail");
    }

    // A log tail's prefactor carries a relative 2·10⁻¹⁴ from ln Γ on its own size, and a quantile inherits it.
    private static double Prefactor(double family, double a, double b, double x)
    {
        double logX = Math.Log(Math.Abs(x) + 1e-300);
        return 2e-14 * family switch
        {
            1.0 => Math.Abs(SpecialFunctions.LogGamma(a)) + Math.Abs(a * logX),
            2.0 => Math.Abs(SpecialFunctions.LogBeta(a, b)) + Math.Abs(a * logX),
            3.0 => Math.Abs(SpecialFunctions.LogBeta(0.5 * a, 0.5)) + Math.Abs(a * logX),
            _ => 0.0,
        };
    }

    private static double Ours(double family, double a, double b, double tail, bool lower)
    {
        IContinuousDistribution d = family switch
        {
            1.0 => Gamma.ShapeRate(a, 1.0),
            2.0 => Beta.Shapes(a, b),
            3.0 => StudentT.Standard(a),
            _ => Normal.Standard,
        };
        return lower ? d.InverseLogCdf(tail) : d.InverseLogSurvival(tail);
    }

    private static void Add(List<double[]> cases, double family, double a, double b)
    {
        foreach (double tail in LogTails)
        {
            cases.Add([family, a, b, tail, 1.0, 0.0]);
            cases.Add([family, a, b, tail, 0.0, 0.0]);
        }
    }

    [Theory]
    [InlineData(1, 0.0, 1.0, 10.0, double.PositiveInfinity)]
    [InlineData(1, 0.0, 1.0, 40.0, double.PositiveInfinity)]
    [InlineData(1, 0.0, 1.0, double.NegativeInfinity, -10.0)]
    [InlineData(1, 0.0, 1.0, 10.0, 10.05)]
    [InlineData(1, 5.0, 0.1, 0.0, 1.0)]
    [InlineData(1, 0.0, 1.0, -1.0, 2.0)]
    [InlineData(1, 0.0, 1.0, 0.0, double.PositiveInfinity)]
    [InlineData(2, 2.0, 1.0, 40.0, double.PositiveInfinity)]
    [InlineData(2, 2.0, 1.0, 1000.0, double.PositiveInfinity)]
    [InlineData(2, 2.0, 1.0, 40.0, 41.0)]
    [InlineData(2, 1000.0, 1.0, 0.0, 800.0)]
    [InlineData(2, 1000.0, 1.0, 1000.0, 1001.0)]
    [InlineData(2, 0.5, 1.0, 0.0, 1e-3)]
    [InlineData(2, 0.01, 1.0, 5.0, 6.0)]
    [InlineData(2, 3.5, 100.0, 0.5, double.PositiveInfinity)]
    [InlineData(3, 50.0, 1e-3, 1.0, double.PositiveInfinity)]
    [InlineData(3, 500.0, 1.0, 1.0, double.PositiveInfinity)]
    [InlineData(3, 3.0, 2.0, 0.0, 0.05)]
    [InlineData(3, 3.0, 2.0, 0.05, 0.051)]
    [InlineData(3, 2.5, 1e-4, 10.0, double.PositiveInfinity)]
    [InlineData(3, 2.5, 1.0, 0.1, 10.0)]
    public void Draws_MatchExactTruncatedMoments(int family, double a, double b, double lower, double upper)
    {
        if (!ROracle.IsAvailable) return;
        JsonElement r = Run([[family, a, b, lower, upper]], "moments");
        double mean = Doubles(r, "mean")[0], variance = Doubles(r, "variance")[0], logMass = Doubles(r, "log_mass")[0];

        IContinuousDistribution parent = family switch
        {
            1 => Normal.MeanSd(a, b),
            2 => Gamma.ShapeRate(a, b),
            _ => InverseGamma.ShapeScale(a, b),
        };
        // The normalizer, in R's tails against this side's.
        Assert.True(Math.Abs(parent.LogIntervalProbability(lower, upper) - logMass) <= 1e-11 * Math.Max(1.0, Math.Abs(logMass)),
                    $"log mass {parent.LogIntervalProbability(lower, upper):R} against R's {logMass:R}");

        var rng = new Xoshiro256PlusPlus(20261102);
        var draws = new double[Draws];
        for (int i = 0; i < Draws; i++) draws[i] = parent.SampleOn(rng, lower, upper);
        // Taken about R's mean, so the sample's own rounding does not swamp a far-tail spread.
        double sum = 0.0, squares = 0.0, fourth = 0.0;
        foreach (double x in draws)
        {
            double d = x - mean;
            sum += d;
            squares += d * d;
            fourth += d * d * d * d;
        }
        double gap = sum / Draws, spread = squares / Draws;
        double kurtosisTerm = Math.Sqrt(Math.Max(fourth / Draws - spread * spread, 0.0));
        _out.WriteLine($"mean gap {gap:E2} (se {Math.Sqrt(variance / Draws):E2}), variance {spread:R} against R's {variance:R}");
        Assert.True(Math.Abs(gap) <= 5.0 * Math.Sqrt(variance / Draws) + 1e-13 * Math.Abs(mean), $"mean gap {gap:R}, variance {variance:R}");
        Assert.True(Math.Abs(spread - variance) <= 5.0 * kurtosisTerm / Math.Sqrt(Draws) + 1e-10 * mean * mean * 1e-3, $"variance {spread:R} against {variance:R}");
    }

    private static JsonElement Run(IEnumerable<double[]> rows, string mode)
    {
        string fixture = Path.Combine(Path.GetTempPath(), $"distribution_quantile_{Guid.NewGuid():N}.csv");
        File.WriteAllLines(fixture, rows.Select(row => string.Join(",", row.Select(Format))));
        try
        {
            return ROracle.Run("oracles/distribution_quantile_oracle.R", fixture, mode);
        }
        finally
        {
            File.Delete(fixture);
        }
    }

    private static string Format(double v)
        => double.IsPositiveInfinity(v) ? "Inf" : double.IsNegativeInfinity(v) ? "-Inf" : v.ToString("R", CultureInfo.InvariantCulture);

    private static double[] Doubles(JsonElement r, string name)
    {
        JsonElement e = r.GetProperty(name);
        if (e.ValueKind == JsonValueKind.Number) return [e.GetDouble()];
        return [.. e.EnumerateArray().Select(v => v.ValueKind == JsonValueKind.Number ? v.GetDouble()
                                                 : v.ValueKind == JsonValueKind.String && v.GetString() == "Inf" ? double.PositiveInfinity
                                                 : v.ValueKind == JsonValueKind.String && v.GetString() == "-Inf" ? double.NegativeInfinity
                                                 : double.NaN)];
    }
}
