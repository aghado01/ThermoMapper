using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Profiling;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// P1: SpPlateau — BWD1996 §C SP-plateau midpoint detection.
/// Temperatures in SweepProfile are ascending (SweepCurves.ByTemperature uses SortedDictionary).
/// </summary>
public sealed class SpPlateauTests
{
    private static SweepProfile MakeProfile(double[] temperatures, double[] chi)
        => new(
            temperatures,
            chi,
            Array.Empty<double>(),
            Array.Empty<double>(),
            null,
            new Dictionary<string, IReadOnlyList<double>>());

    [Fact]
    public void SpPlateau_NominalCase_ReturnsMidpoint()
    {
        // χ: peak at T=0.10, plateau, then cliff descending toward T=0.25.
        // Paper semantics (BWD1996 §C): T_ps = departure from the plateau —
        // the first smoothed value below 80% of the peak (0.83 < 1.147 at
        // index 3, T=0.20).
        double[] temps = { 0.05, 0.10, 0.15, 0.20, 0.25 };
        double[] chi   = { 0.50, 2.00, 1.80, 0.40, 0.30 };

        var result = SpcProfileAnalysis.SpPlateau(MakeProfile(temps, chi));

        Assert.True(result.CliffFound);
        Assert.Equal(0.10, result.TFs,   precision: 6);
        Assert.Equal(0.20, result.TPs,   precision: 6);  // first smoothed value below 0.8 × peak
        Assert.Equal(0.15, result.TClus, precision: 6);  // (0.10 + 0.20) / 2
        Assert.True(result.TFs < result.TClus && result.TClus < result.TPs);
    }

    [Fact]
    public void SpPlateau_FlatAbovePeak_FallsBackToTFs()
    {
        // χ is flat above the peak — no ratio strictly > 1.0 → no cliff.
        double[] temps = { 0.05, 0.10, 0.15, 0.20, 0.25 };
        double[] chi   = { 0.50, 2.00, 2.00, 2.00, 2.00 };

        var result = SpcProfileAnalysis.SpPlateau(MakeProfile(temps, chi));

        Assert.False(result.CliffFound);
        Assert.Equal(result.TFs, result.TClus, precision: 10);
    }

    [Fact]
    public void SpPlateau_PeakAtSecondToLastPoint_FallsBackToTFs()
    {
        // Peak at index n-2: only one point above it, too few to find a cliff.
        double[] temps = { 0.05, 0.10, 0.20 };
        double[] chi   = { 0.30, 0.50, 2.00 };

        var result = SpcProfileAnalysis.SpPlateau(MakeProfile(temps, chi));

        Assert.False(result.CliffFound);
        Assert.Equal(0.20, result.TFs, precision: 6);
        Assert.Equal(result.TFs, result.TClus, precision: 10);
    }

    [Fact]
    public void SpPlateau_SinglePoint_FallsBackToTFs()
    {
        double[] temps = { 0.10 };
        double[] chi   = { 1.50 };

        var result = SpcProfileAnalysis.SpPlateau(MakeProfile(temps, chi));

        Assert.False(result.CliffFound);
        Assert.Equal(0.10, result.TFs, precision: 6);
        Assert.Equal(0.10, result.TClus, precision: 6);
    }

    [Fact]
    public void SpPlateau_MultipleDrops_PicksLargestAbsoluteDrop()
    {
        // Anti-ratio-explosion: the plateau departure (smoothed 1.40 < 0.8 ×
        // 1.90 at index 4) must win over the larger drop-RATIO sitting in
        // the near-zero tail (0.80/0.21 ≈ 3.9 at index 5) — BWD's "χ
        // abruptly diminishes" is the END of the near-constant plateau.
        // Ratio semantics would pick T=0.30 here.
        double[] temps = { 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40 };
        double[] chi   = { 0.40, 1.60, 2.00, 1.90, 1.80, 0.50, 0.10, 0.02 };

        var result = SpcProfileAnalysis.SpPlateau(MakeProfile(temps, chi));

        Assert.True(result.CliffFound);
        Assert.Equal(0.20,  result.TFs, precision: 6);   // smoothed peak (3-pt smoothing shifts the raw 0.15 peak)
        Assert.Equal(0.25,  result.TPs, precision: 6);   // largest absolute smoothed drop
        Assert.Equal(0.225, result.TClus, precision: 6); // (0.20 + 0.25) / 2
        Assert.True(result.TFs < result.TClus && result.TClus < result.TPs);
    }
}
