using System;
using System.IO;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// Wallenius' set probabilities, the variable-set law of BASS's births (D01), against F-O5
/// (<c>r/oracles/bass_wallenius_oracle.R</c>; stage 4's S3): brute-force sums over every draw order, each term and
/// denominator a sum of positives, for every subset of each order. The weights are BASS's z = w₂ + appearances, at the
/// default w₂ = 5 with skewed counts and at a small w₂ with counts in the hundreds, where R's own inclusion–exclusion
/// (<c>dmwnchBass</c>) loses digits; its worst error is printed for the record. Opt-in: it runs where the R toolchain
/// and the BASS checkout are present.
/// </summary>
public sealed class WalleniusParityTests
{
    private static readonly string BassSource = Environment.GetEnvironmentVariable("BASS_R_SOURCE") ?? "D:/aghado01/packages/BASS";
    private readonly ITestOutputHelper _out;
    public WalleniusParityTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData("5,45,5,8,125,5,12,5,6", "1,2,3,4,5,6")]
    [InlineData("0.01,200.01,0.01,50.01,0.01,3.01,0.01,900.01", "1,2,3,4,5")]
    public void SetProbabilities_MatchBruteForce(string weights, string orders)
    {
        if (!ROracle.IsAvailable || !Directory.Exists(BassSource)) return;
        double[] z = [.. weights.Split(',').Select(double.Parse)];
        string fixture = Path.Combine(Path.GetTempPath(), $"bass_wallenius_{Guid.NewGuid():N}.csv");
        File.WriteAllText(fixture, weights + "\n");
        try
        {
            JsonElement r = ROracle.Run("oracles/bass_wallenius_oracle.R", fixture, orders, "2000");
            foreach (JsonElement order in r.GetProperty("orders").EnumerateArray())
            {
                int k = order.GetProperty("K").GetInt32();
                Assert.True(order.GetProperty("enumerated").GetBoolean());
                JsonElement[] subsets = [.. order.GetProperty("subsets").EnumerateArray()];
                double[] brute = [.. order.GetProperty("brute").EnumerateArray().Select(v => v.GetDouble())];
                double worst = 0.0;
                for (int s = 0; s < subsets.Length; s++)
                {
                    int[] set = [.. subsets[s].EnumerateArray().Select(v => v.GetInt32() - 1)];
                    Assert.Equal(k, set.Length);
                    double ours = Math.Exp(Wallenius.LogProbability(z, set));
                    double relative = Math.Abs(ours / brute[s] - 1.0);
                    worst = Math.Max(worst, relative);
                    Assert.True(relative <= 1e-12, $"K = {k}, set {{{string.Join(",", set)}}}: {ours:R} against {brute[s]:R}");
                }
                _out.WriteLine($"K = {k}: {subsets.Length} sets, worst {worst:E1}; R's dmwnchBass worst {order.GetProperty("max_rel_r").GetDouble():E1}");
            }
        }
        finally
        {
            File.Delete(fixture);
        }
    }
}
