using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's truncation residuals (stage 5, O1; D07) against the O-O4 oracle (<c>bass_truncation_oracle.R</c>), which builds
/// the leave-one-case-out residuals by n explicit refits in output space and the in-sample ones by projection, each rank
/// judged at max(n, m)·ε·d₁ for the whole data: the centred rank, r, both residual sets, the leave-one-out ranks, their RMS
/// and the basis-stability ratio, at n = 2, 5, 10 and 40 under the 99% rule and at full retention. The curves follow the
/// stage-1 checks' recipe with this side's noise: the oracle reads the data it is given. Opt-in: it runs where the R
/// toolchain is present.
/// </summary>
public sealed class BassOutputTruncationParityTests
{
    private readonly ITestOutputHelper _out;
    public BassOutputTruncationParityTests(ITestOutputHelper output) => _out = output;

    // n curves at m locations whose phase and amplitude move with two inputs, with noise of sd 0.01.
    internal static double[,] Curves(int n, int m, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            double z1 = rng.NextDouble(), z2 = rng.NextDouble();
            for (int l = 0; l < m; l++)
            {
                double t = l / (m - 1.0);
                double noise = Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());
                y[i, l] = Math.Sin(2 * Math.PI * (t - 0.3 * z1)) * (1 + z2) + 0.5 * z2 * t * t + 0.01 * noise;
            }
        }
        return y;
    }

    internal static JsonElement Run(string oracle, double[,] y, params string[] options)
    {
        string csv = Path.Combine(Path.GetTempPath(), $"bass_output_{Guid.NewGuid():N}.csv");
        var text = new StringBuilder();
        for (int i = 0; i < y.GetLength(0); i++)
            text.AppendLine(string.Join(",", Enumerable.Range(0, y.GetLength(1)).Select(l => y[i, l].ToString("R", CultureInfo.InvariantCulture))));
        File.WriteAllText(csv, text.ToString());
        try
        {
            return ROracle.Run(oracle, csv, options);
        }
        finally
        {
            File.Delete(csv);
        }
    }

    // jsonlite writes a matrix by rows.
    internal static double[,] Matrix(JsonElement e)
    {
        JsonElement[] rows = [.. e.EnumerateArray()];
        int columns = rows[0].GetArrayLength();
        var a = new double[rows.Length, columns];
        for (int i = 0; i < rows.Length; i++)
        {
            int l = 0;
            foreach (JsonElement v in rows[i].EnumerateArray()) a[i, l++] = v.GetDouble();
        }
        return a;
    }

    private static double MaxGap(double[,] a, double[,] b)
    {
        Assert.Equal(a.GetLength(0), b.GetLength(0));
        Assert.Equal(a.GetLength(1), b.GetLength(1));
        double gap = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int l = 0; l < a.GetLength(1); l++) gap = Math.Max(gap, Math.Abs(a[i, l] - b[i, l]));
        return gap;
    }

    [Theory]
    [InlineData(2, 99.0)]
    [InlineData(5, 99.0)]
    [InlineData(10, 99.0)]
    [InlineData(40, 99.0)]
    [InlineData(5, 100.0)]
    public void Residuals_MatchTheOracle(int n, double percent)
    {
        if (!ROracle.IsAvailable) return;
        const int m = 200;
        double[,] y = Curves(n, m, seed: 100 + n);
        JsonElement r = Run("oracles/bass_truncation_oracle.R", y, percent.ToString("R", CultureInfo.InvariantCulture));
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = percent / 100 });

        double scale = y.Cast<double>().Max(Math.Abs);
        double inSampleGap = MaxGap(Matrix(r.GetProperty("in_sample")), d.Truncation.InSample!);
        double looGap = MaxGap(Matrix(r.GetProperty("loo")), d.Truncation.LeaveOneOut!);
        int[] looRanks = [.. r.GetProperty("loo_rank").EnumerateArray().Select(v => v.GetInt32())];
        double? ratio = r.GetProperty("ratio").ValueKind == JsonValueKind.Null ? null : r.GetProperty("ratio").GetDouble();
        _out.WriteLine($"n = {n}, {percent}%: r = {d.Rank} of {d.Pca!.NumericalRank} (R {r.GetProperty("r").GetInt32()} of " +
                       $"{r.GetProperty("centred_rank").GetInt32()}); in-sample gap {inSampleGap:E1}, leave-one-out gap {looGap:E1}; " +
                       $"ratio {d.Truncation.StabilityRatio?.ToString("F3") ?? "null"} (R {ratio?.ToString("F3") ?? "null"})");

        Assert.Equal(r.GetProperty("centred_rank").GetInt32(), d.Pca.NumericalRank);
        Assert.Equal(r.GetProperty("r").GetInt32(), d.Rank);
        Assert.Equal(r.GetProperty("full_retention").GetBoolean(), d.Rank == d.Pca.NumericalRank);
        Assert.Equal(r.GetProperty("zero_tolerance").GetDouble(), d.Pca.RankTolerance, 1e-12 * d.Pca.RankTolerance);
        Assert.Equal(looRanks, d.Truncation.LeaveOneOutRanks!);
        Assert.InRange(inSampleGap, 0, 1e-12 * scale);
        Assert.InRange(looGap, 0, 1e-10 * scale);
        Assert.Equal(r.GetProperty("loo_rms").GetDouble(), d.Truncation.LeaveOneOutRms!.Value, 1e-10 * scale);
        Assert.Equal(ratio is null, d.Truncation.StabilityRatio is null);
        if (ratio is { } expected) Assert.Equal(expected, d.Truncation.StabilityRatio!.Value, 1e-8 * expected);
        if (n == 2)
            for (int l = 0; l < m; l++) Assert.Equal(y[0, l] - y[1, l], d.Truncation.LeaveOneOut![0, l], 1e-12 * scale);
    }
}
