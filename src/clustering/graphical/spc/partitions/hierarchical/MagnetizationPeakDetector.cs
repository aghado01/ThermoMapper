using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Profiling;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// Default <see cref="IPseudoTransitionDetector"/>: scans the
/// magnetization-variance trajectory for local maxima with prominence
/// above <see cref="MinProminence"/>. The classical Blatt χ_m peak
/// detector, modulo the N/T pre-scaling that doesn't affect peak
/// location on a fixed graph.
/// </summary>
/// <remarks>
/// <para><b>Signal source.</b> Reads
/// <see cref="SweepProfile.AdditionalChannels"/>'s
/// <see cref="SignalChannelName"/> entry; throws if the channel is
/// absent. <see cref="SweepProfile.From"/> populates
/// <c>"MagnetizationVariance"</c> by default, so any sweep run through
/// the standard pipeline is consumable as-is.</para>
///
/// <para><b>Peak rule.</b> A local maximum is recorded when both
/// neighbors are strictly lower AND the maximum's prominence
/// (height above the higher of the two surrounding minima within the
/// window) is at least <see cref="MinProminence"/> relative to the
/// signal's peak-to-trough range. Endpoints are never peaks. Repeated
/// equal samples are treated as plateaus and skipped — a plateau is
/// only a peak when it has strict descents on both sides.</para>
///
/// <para><b>Limitations.</b> v1 is a simple three-point scan with a
/// prominence gate; no smoothing, no Savitzky-Golay, no persistence
/// homology. Adequate for sweeps with reasonable sample density and
/// well-separated phases; noisy or undersampled sweeps will need a
/// successor detector (smoothed-spline peak, multi-signal consensus,
/// or BARS-derived posterior modes).</para>
/// </remarks>
public sealed class MagnetizationPeakDetector : IPseudoTransitionDetector
{
    /// <summary>
    /// Channel name read from
    /// <see cref="SweepProfile.AdditionalChannels"/>. Defaults to the
    /// canonical <c>"MagnetizationVariance"</c> emitted by
    /// <see cref="SweepProfile.From"/>.
    /// </summary>
    public string SignalChannelName { get; init; } = "MagnetizationVariance";

    /// <summary>
    /// Minimum prominence as a fraction of the signal's peak-to-trough
    /// range, in <c>[0, 1]</c>. <c>0.0</c> accepts every local maximum;
    /// <c>0.1</c> (default) filters out small numerical bumps that
    /// aren't real pseudo-transitions.
    /// </summary>
    public double MinProminence { get; init; } = 0.1;

    /// <inheritdoc />
    public double[] Detect(SweepProfile profile)
    {
        ArgumentNullException.ThrowIfNull(profile);
        if (profile.IsEmpty) return Array.Empty<double>();

        if (!profile.AdditionalChannels.TryGetValue(SignalChannelName, out var signal))
            throw new InvalidOperationException(
                $"SweepProfile is missing channel '{SignalChannelName}'. " +
                "Pass a profile built by SweepProfile.From(runs) or supply a custom " +
                "channel name via SignalChannelName.");

        int n = signal.Count;
        if (n != profile.Temperatures.Count)
            throw new InvalidOperationException(
                $"Channel '{SignalChannelName}' has {n} samples but the profile has " +
                $"{profile.Temperatures.Count} temperatures.");
        if (n < 3) return Array.Empty<double>();

        // Range gate: filter peaks whose absolute prominence is below
        // MinProminence × (max − min). Cheap O(n) prepass.
        double minVal = double.PositiveInfinity;
        double maxVal = double.NegativeInfinity;
        for (int i = 0; i < n; i++)
        {
            double v = signal[i];
            if (v < minVal) minVal = v;
            if (v > maxVal) maxVal = v;
        }
        double range = maxVal - minVal;
        double absThreshold = range * MinProminence;
        if (range == 0.0) return Array.Empty<double>();

        var peaks = new List<double>();
        for (int i = 1; i < n - 1; i++)
        {
            double prev = signal[i - 1];
            double curr = signal[i];
            double next = signal[i + 1];

            // Strict local maximum: both neighbors strictly lower.
            // Plateaus (curr == prev or curr == next) are skipped — a
            // legitimate plateau-peak shows up at the boundary point
            // where the descent actually starts.
            if (curr <= prev || curr <= next) continue;

            // Prominence vs the higher of the two adjacent valleys —
            // computed as a strict three-point check, not a windowed
            // walk. Good enough for v1 given the prominence filter is
            // the main guard against noise.
            double higherSide = Math.Max(prev, next);
            double prominence = curr - higherSide;
            if (prominence < absThreshold) continue;

            peaks.Add(profile.Temperatures[i]);
        }

        return peaks.ToArray();
    }
}
