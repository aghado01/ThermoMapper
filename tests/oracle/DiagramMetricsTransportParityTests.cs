using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using TDA.Ph;
using Xunit;

namespace Maths.Oracle.Tests;

/// <summary>
/// T4transport / lpSolve parity for the DiagramMetrics distance backends (the "T4transport
/// oracle" clause of the ISOLET H0 matching-cost gate). The R side reconstructs the
/// diagonal-augmented balanced geometry independently from the raw bars; T4transport masses are
/// probability (1/s per point) vs the C# unit masses, so its sliced/Sinkhorn distances carry the
/// factor s^(1/p). The fixture is sized (s = 100) so T4transport's interpolated-quantile
/// smoothing of the per-slice 1-D transport is a few-percent residual (measured 3-4%), absorbed
/// by the sliced tolerance. Skips silently when the R toolchain is absent.
/// </summary>
public sealed class DiagramMetricsTransportParityTests
{
    private const int NumProjections = 20000;   // Monte Carlo on the R side; ~1% mean error
    private const int Seed = 73;

    // Deterministic mixed-scale finite bars, unequal counts (n = 60 vs m = 40), s = 100 — large
    // enough that T4transport's quantile interpolation bias is small (see class doc).
    private static (double Birth, double Death)[] BarsA()
    {
        var bars = new (double, double)[60];
        for (int i = 0; i < bars.Length; i++)
        {
            double birth = 0.05 * i + 0.35 * Math.Sin(1.1 * i);
            bars[i] = (birth, birth + 0.15 + 0.8 * Math.Abs(Math.Cos(0.9 * i)));
        }
        return bars;
    }

    private static (double Birth, double Death)[] BarsB()
    {
        var bars = new (double, double)[40];
        for (int j = 0; j < bars.Length; j++)
        {
            double birth = 0.07 * j + 0.30 * Math.Cos(0.8 * j);
            bars[j] = (birth, birth + 0.20 + 0.7 * Math.Abs(Math.Sin(1.3 * j)));
        }
        return bars;
    }

    [Theory]
    [InlineData(1.0)]
    [InlineData(2.0)]
    public void DiagramDistances_MatchT4TransportAndLpAssign(double p)
    {
        if (!ROracle.IsAvailable) return;

        const double epsilon = 0.05;   // matched smoothing: lambda = epsilon * cMax on the R side
        string csv = WriteBarsCsv();
        try
        {
            JsonElement r = ROracle.Run(
                "oracles/transport_oracle.R",
                csv,
                p.ToString(CultureInfo.InvariantCulture),
                NumProjections.ToString(CultureInfo.InvariantCulture),
                Seed.ToString(CultureInfo.InvariantCulture),
                epsilon.ToString(CultureInfo.InvariantCulture));

            double massScale = Math.Pow(r.GetProperty("s").GetDouble(), 1.0 / p);
            Barcode a = Diagram(BarsA());
            Barcode b = Diagram(BarsB());

            // Exact Hungarian vs the external assignment LP (both exact; the LP is integral).
            double exact = DiagramMetrics.Wasserstein(a, b, dimension: 0, p);
            double lp = Math.Pow(r.GetProperty("lp_cost").GetDouble(), 1.0 / p);
            Assert.InRange(Math.Abs(exact - lp) / lp, 0.0, 1e-8);

            // Sliced: deterministic slices vs T4transport's Monte Carlo of the same integral —
            // slicing sees only 1-D projected distances, so this is full-semantics parity. The
            // tolerance covers T4transport's interpolated-quantile smoothing (measured 3-4% at
            // s = 100) plus ~1% Monte Carlo error; structural errors (mass scale, combination
            // convention, metric slips) are all far larger.
            double sliced = DiagramMetrics.SlicedWasserstein(a, b, dimension: 0, p, directions: 2000);
            double swdist = massScale * r.GetProperty("swdist").GetDouble();
            Assert.InRange(Math.Abs(sliced - swdist) / swdist, 0.0, 0.08);

            // Sinkhorn at matched smoothing: same kernel exp(-C/(epsilon*cMax)), plans agree up
            // to the mass scale. This is the external Sinkhorn check; the eps -> 0 limit is
            // pinned by an externally closed chain instead — lp_cost above ties our Hungarian to
            // lpSolve at 1e-8, and the unit suite ties Sinkhorn's small-eps limit to that
            // Hungarian at <= 1% (driving eps -> 0 at s = 100 would need far more than 2e4
            // iterations — the slow-convergence regime already recorded in the ISOLET brief).
            double sink = DiagramMetrics.SinkhornWasserstein(
                a, b, dimension: 0, p, epsilon: epsilon, maxIters: 20000);
            double sinkhornD = massScale * r.GetProperty("sinkhorn").GetDouble();
            Assert.InRange(Math.Abs(sink - sinkhornD) / sinkhornD, 0.0, 1e-2);
        }
        finally
        {
            File.Delete(csv);
        }
    }

    private static Barcode Diagram((double Birth, double Death)[] bars) =>
        new(bars.Select(bar => new Bar(bar.Birth, bar.Death, 0)).ToArray());

    private static string WriteBarsCsv()
    {
        string path = Path.Combine(Path.GetTempPath(), $"diagrams_{Guid.NewGuid():N}.csv");
        var sb = new StringBuilder();
        foreach ((double birth, double death) in BarsA())
            sb.AppendLine(string.Join(",", "0",
                birth.ToString("R", CultureInfo.InvariantCulture),
                death.ToString("R", CultureInfo.InvariantCulture)));
        foreach ((double birth, double death) in BarsB())
            sb.AppendLine(string.Join(",", "1",
                birth.ToString("R", CultureInfo.InvariantCulture),
                death.ToString("R", CultureInfo.InvariantCulture)));
        File.WriteAllText(path, sb.ToString());
        return path;
    }
}
