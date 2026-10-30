using System;
using System.Collections.Generic;

namespace Maths.Samplers.FlatHistogram;

/// <summary>
/// The flat-histogram diagnostics of one phase's visits, learning or production: the hats' occupancy over the window
/// against its flat target, the visits outside the window, and tunnelling, the round trips each chain makes between the
/// window's lowest and highest bins (a round trip leaves the lowest bin, reaches the highest and returns, as a tempering
/// ladder's leaves the hot end; a bin is a node's half-spacing cell, the edge bins taking the extension beyond).
/// </summary>
/// <param name="Occupancy">Per node, the visits' summed hat values φ_k(ℓ).</param>
/// <param name="Desired">Per node, the flat occupancy π_k (<see cref="MulticanonicalWeights.FlatOccupancy"/>).</param>
/// <param name="InWindow">Visits inside the window.</param>
/// <param name="BelowWindow">Visits below it.</param>
/// <param name="AboveWindow">Visits above it.</param>
/// <param name="RoundTrips">Per chain, the round trips completed.</param>
/// <param name="LastEnd">Per chain, the end it was last at: −1 the lowest bin, +1 the highest, 0 neither yet; restart state.</param>
public sealed record FlatHistogramDiagnostics(
    double[] Occupancy, double[] Desired, long InWindow, long BelowWindow, long AboveWindow, long[] RoundTrips, sbyte[] LastEnd)
{
    /// <summary>All visits.</summary>
    public long Visits => InWindow + BelowWindow + AboveWindow;

    /// <summary>Each node's occupancy share inside the window over its flat share: 1 everywhere on a flat histogram;
    /// NaN with no visit inside.</summary>
    public double[] Flatness
    {
        get
        {
            var ratio = new double[Occupancy.Length];
            for (int k = 0; k < ratio.Length; k++) ratio[k] = InWindow > 0 ? Occupancy[k] / (InWindow * Desired[k]) : double.NaN;
            return ratio;
        }
    }

    /// <summary>The histogram's largest departure from flat, max_k |occupancy share / π_k − 1|, Wang–Landau's flatness
    /// criterion read as a diagnostic: about 0 when flat, 1 when some node is never visited; NaN with no visit inside.</summary>
    public double MaxFlatnessDeviation
    {
        get
        {
            double worst = InWindow > 0 ? 0.0 : double.NaN;
            foreach (double r in Flatness) worst = Math.Max(worst, Math.Abs(r - 1.0));
            return worst;
        }
    }

    /// <summary>The chains that completed no round trip: they did not cross the window within the phase.</summary>
    public int[] ChainsWithoutRoundTrip
    {
        get
        {
            var none = new List<int>();
            for (int c = 0; c < RoundTrips.Length; c++) if (RoundTrips[c] == 0) none.Add(c);
            return [.. none];
        }
    }

    /// <summary>The round trips summed over chains.</summary>
    public long TotalRoundTrips
    {
        get
        {
            long total = 0;
            foreach (long r in RoundTrips) total += r;
            return total;
        }
    }

    /// <summary>The chains' diagnostics merged in chain order, each holding one chain.</summary>
    public static FlatHistogramDiagnostics Merge(IReadOnlyList<FlatHistogramDiagnostics> chains)
    {
        ArgumentNullException.ThrowIfNull(chains);
        if (chains.Count == 0) throw new ArgumentException("Nothing to merge.", nameof(chains));
        int k = chains[0].Occupancy.Length;
        var occupancy = new double[k];
        long inside = 0, below = 0, above = 0;
        var trips = new List<long>();
        var last = new List<sbyte>();
        foreach (FlatHistogramDiagnostics d in chains)
        {
            for (int i = 0; i < k; i++) occupancy[i] += d.Occupancy[i];
            inside += d.InWindow;
            below += d.BelowWindow;
            above += d.AboveWindow;
            trips.AddRange(d.RoundTrips);
            last.AddRange(d.LastEnd);
        }
        return new FlatHistogramDiagnostics(occupancy, (double[])chains[0].Desired.Clone(), inside, below, above, [.. trips], [.. last]);
    }
}

// The running form of the diagnostics over a fixed set of chains: a visit at a time, in chain order.
internal sealed class HistogramAccumulator
{
    private readonly MulticanonicalWeights _weights;   // the window and nodes; only the shape is read
    private readonly double[] _occupancy;
    private readonly long[] _roundTrips;
    private readonly sbyte[] _lastEnd;
    private long _inside, _below, _above;

    public HistogramAccumulator(MulticanonicalWeights weights, int chains)
    {
        _weights = weights;
        _occupancy = new double[weights.Nodes];
        _roundTrips = new long[chains];
        _lastEnd = new sbyte[chains];
    }

    public HistogramAccumulator(MulticanonicalWeights weights, FlatHistogramDiagnostics saved)
    {
        _weights = weights;
        if (saved.Occupancy.Length != weights.Nodes) throw new ArgumentException("The saved diagnostics have another node count.", nameof(saved));
        _occupancy = (double[])saved.Occupancy.Clone();
        _roundTrips = (long[])saved.RoundTrips.Clone();
        _lastEnd = (sbyte[])saved.LastEnd.Clone();
        _inside = saved.InWindow;
        _below = saved.BelowWindow;
        _above = saved.AboveWindow;
    }

    public void Add(int chain, double level)
    {
        if (_weights.TryLocate(level, out int segment, out double fraction))
        {
            _occupancy[segment] += 1.0 - fraction;
            _occupancy[segment + 1] += fraction;
            _inside++;
        }
        else if (level < _weights.Low) _below++;
        else if (level > _weights.High) _above++;

        double half = 0.5 * _weights.Spacing;
        if (level < _weights.Low + half)
        {
            if (_lastEnd[chain] == +1) _roundTrips[chain]++;
            _lastEnd[chain] = -1;
        }
        else if (level > _weights.High - half && _lastEnd[chain] == -1) _lastEnd[chain] = +1;
    }

    public FlatHistogramDiagnostics Snapshot()
        => new((double[])_occupancy.Clone(), _weights.FlatOccupancy(), _inside, _below, _above,
               (long[])_roundTrips.Clone(), (sbyte[])_lastEnd.Clone());
}
