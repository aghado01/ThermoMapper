// ============================================================================
// TDA.Mapper — ICover.cs
// ============================================================================
// 1-D cover interfaces and standard implementations.
//
// A 1-D cover partitions the range of a scalar filter into overlapping intervals;
// each interval defines a preimage that MAPPER then locally clusters.
//
// Two cover families:
//   UniformCover            — equal-width intervals. Standard default.
//   BalancedHistogramCover  — equal-population bins via filter-value quantiles.
//                             Adaptive to skewed filter distributions; recommended
//                             for hyperbolic data where filter values can be
//                             heavy-tailed (e.g., Poincaré radial near boundary).
//
// Multi-D covers (product cubes) live in Multi.cs.
// ============================================================================

#nullable enable
using System;
using System.Collections.Generic;

namespace TDA.Mapper.Cover;

// ── 1-D cover interface ─────────────────────────────────────────────────────

public interface ICover
{
    /// <summary>Generate cover bins from a vector of N scalar filter values.</summary>
    CoverResult Generate(double[] filterValues);

    string Name { get; }
}

public sealed record CoverResult(
    IReadOnlyList<CoverBin> Bins,
    double FilterMin,
    double FilterMax);

public sealed record CoverBin(
    int BinId,
    double Lower,
    double Upper,
    List<int> PointIndices);

// ── Multi-D cover interface (sibling) ───────────────────────────────────────

/// <summary>
/// Multi-D cover for product-cover MAPPER. Implementations live in Multi.cs.
/// </summary>
public interface IMultiCover
{
    int Dimension { get; }
    CoverResult Generate(double[][] filterValues);   // N × Dimension input
    string Name { get; }
}

// ── Factory ─────────────────────────────────────────────────────────────────

public static class Cover
{
    /// <summary>Equal-width intervals over [min, max] of filter values, with
    /// each adjacent pair overlapping by <paramref name="overlapPercent"/> of
    /// the interval width. Standard MAPPER default.</summary>
    public static UniformCover Uniform(int numIntervals = 12, double overlapPercent = 0.4)
        => new(numIntervals, overlapPercent);

    /// <summary>Equal-population bins via filter-value quantiles. Each bin holds
    /// roughly N/numIntervals points, with neighboring bins overlapping by
    /// <paramref name="overlapPercent"/> of their populations. Recommended for
    /// skewed filter distributions (e.g., Poincaré radial near boundary).</summary>
    public static BalancedHistogramCover Balanced(int numIntervals = 12, double overlapPercent = 0.4)
        => new(numIntervals, overlapPercent);
}

// ── UniformCover ────────────────────────────────────────────────────────────

public sealed class UniformCover : ICover
{
    public int NumIntervals { get; }
    public double OverlapPercent { get; }

    public string Name => $"Uniform (n={NumIntervals}, overlap={OverlapPercent:P0})";

    public UniformCover(int numIntervals, double overlapPercent)
    {
        if (numIntervals < 2)
            throw new ArgumentOutOfRangeException(nameof(numIntervals), "numIntervals must be >= 2");
        if (overlapPercent < 0 || overlapPercent >= 1)
            throw new ArgumentOutOfRangeException(nameof(overlapPercent), "overlapPercent must be in [0, 1)");

        NumIntervals = numIntervals;
        OverlapPercent = overlapPercent;
    }

    public CoverResult Generate(double[] filterValues)
    {
        ArgumentNullException.ThrowIfNull(filterValues);
        if (filterValues.Length == 0)
            return new CoverResult(Array.Empty<CoverBin>(), 0.0, 0.0);

        // Compute finite range; ignore +/-Infinity (e.g., Poincaré radial at boundary).
        double minVal = double.PositiveInfinity, maxVal = double.NegativeInfinity;
        for (int i = 0; i < filterValues.Length; i++)
        {
            double v = filterValues[i];
            if (!double.IsFinite(v)) continue;
            if (v < minVal) minVal = v;
            if (v > maxVal) maxVal = v;
        }

        if (!double.IsFinite(minVal) || !double.IsFinite(maxVal))
            return new CoverResult(Array.Empty<CoverBin>(), 0.0, 0.0);

        // Degenerate case: all finite values identical.
        if (maxVal - minVal < 1e-12)
        {
            var idx = new List<int>(filterValues.Length);
            for (int i = 0; i < filterValues.Length; i++)
                if (double.IsFinite(filterValues[i])) idx.Add(i);
            return new CoverResult(new[] { new CoverBin(0, minVal, maxVal, idx) }, minVal, maxVal);
        }

        double range = maxVal - minVal;
        double intervalWidth = range / NumIntervals;
        double overlapWidth = intervalWidth * OverlapPercent;
        double step = intervalWidth - overlapWidth;

        var bins = new List<CoverBin>(NumIntervals);

        for (int i = 0; i < NumIntervals; i++)
        {
            double lower = minVal + i * step;
            double upper = lower + intervalWidth;
            if (i == NumIntervals - 1) upper = maxVal + 1e-9;

            var indices = new List<int>(filterValues.Length / NumIntervals + 8);
            for (int p = 0; p < filterValues.Length; p++)
            {
                double v = filterValues[p];
                if (!double.IsFinite(v)) continue;
                if (v >= lower && v <= upper)
                    indices.Add(p);
            }
            if (indices.Count > 0)
                bins.Add(new CoverBin(i, lower, upper, indices));
        }

        return new CoverResult(bins, minVal, maxVal);
    }
}

// ── BalancedHistogramCover (adaptive, quantile-based) ───────────────────────

/// <summary>
/// Equal-population bins via filter-value quantiles. Each bin contains
/// approximately N/NumIntervals points, with neighboring bins overlapping by
/// <see cref="OverlapPercent"/> of their populations (not their value-width).
///
/// Adaptive to skewed filter distributions — useful when filter values are
/// non-uniform (e.g., Poincaré radial concentrates points near origin if data
/// has dense root cluster; hyperbolic-tail filters in general).
///
/// Compare to UniformCover, which uses equal-width intervals on the value axis
/// and can produce wildly variable bin populations on skewed data.
/// </summary>
public sealed class BalancedHistogramCover : ICover
{
    public int NumIntervals { get; }
    public double OverlapPercent { get; }

    public string Name => $"BalancedHistogram (n={NumIntervals}, overlap={OverlapPercent:P0})";

    public BalancedHistogramCover(int numIntervals, double overlapPercent)
    {
        if (numIntervals < 2)
            throw new ArgumentOutOfRangeException(nameof(numIntervals), "numIntervals must be >= 2");
        if (overlapPercent < 0 || overlapPercent >= 1)
            throw new ArgumentOutOfRangeException(nameof(overlapPercent), "overlapPercent must be in [0, 1)");

        NumIntervals = numIntervals;
        OverlapPercent = overlapPercent;
    }

    public CoverResult Generate(double[] filterValues)
    {
        ArgumentNullException.ThrowIfNull(filterValues);
        if (filterValues.Length == 0)
            return new CoverResult(Array.Empty<CoverBin>(), 0.0, 0.0);

        // Collect finite values with their original indices, sorted by value.
        var pairs = new List<(double Value, int Index)>(filterValues.Length);
        for (int i = 0; i < filterValues.Length; i++)
        {
            double v = filterValues[i];
            if (double.IsFinite(v)) pairs.Add((v, i));
        }
        if (pairs.Count == 0)
            return new CoverResult(Array.Empty<CoverBin>(), 0.0, 0.0);

        pairs.Sort((a, b) => a.Value.CompareTo(b.Value));

        double minVal = pairs[0].Value;
        double maxVal = pairs[^1].Value;

        // Each bin holds roughly nPerBin points; neighbors overlap by overlapCount
        // points on each side (effectively widening the bin by 2·overlapCount).
        int n = pairs.Count;
        int nPerBin = Math.Max(1, n / NumIntervals);
        int overlapCount = (int)Math.Round(nPerBin * OverlapPercent);

        var bins = new List<CoverBin>(NumIntervals);

        for (int i = 0; i < NumIntervals; i++)
        {
            int centerStart = i * nPerBin;
            int centerEnd = (i == NumIntervals - 1) ? n : (i + 1) * nPerBin;

            int lo = Math.Max(0, centerStart - overlapCount);
            int hi = Math.Min(n, centerEnd + overlapCount);

            if (lo >= hi) continue;

            var indices = new List<int>(hi - lo);
            for (int k = lo; k < hi; k++) indices.Add(pairs[k].Index);

            bins.Add(new CoverBin(
                BinId: i,
                Lower: pairs[lo].Value,
                Upper: pairs[hi - 1].Value,
                PointIndices: indices));
        }

        return new CoverResult(bins, minVal, maxVal);
    }
}
