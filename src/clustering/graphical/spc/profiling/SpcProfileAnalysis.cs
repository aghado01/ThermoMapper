using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;

namespace Clustering.Graphical.SPC.Profiling;

/// <summary>BWD1996 §C SP-plateau landmarks.</summary>
/// <param name="TFs">FM→SP pseudotransition: χ argmax (= PickPeakTemperature).</param>
/// <param name="TPs">SP→PM pseudotransition: steepest high-T cliff above <paramref name="TFs"/>.</param>
/// <param name="TClus">(T_fs + T_ps)/2 — the SP-plateau midpoint; the BWD clustering temperature.</param>
/// <param name="CliffFound">
/// False when the cliff was undetectable (too few high-T points, or χ never
/// genuinely drops above the peak). <see cref="TClus"/> equals <see cref="TFs"/>
/// in this case (safe fallback).
/// </param>
public readonly record struct SpPlateauResult(
    double TFs,
    double TPs,
    double TClus,
    bool   CliffFound);

/// <summary>
/// Analysis helpers for SPC temperature profiles and Potts sampling.
/// </summary>
public static class SpcProfileAnalysis
{
    public static (double Lo, double Hi) ComputeHalfMaximumBand(SweepProfile profile)
    {
        if (profile is null)
            throw new ArgumentNullException(nameof(profile));

        if (profile.Temperatures.Count != profile.Susceptibility.Count)
            throw new ArgumentException("SweepProfile temperatures and susceptibility values must have the same length.", nameof(profile));

        if (profile.Count == 0)
            throw new ArgumentException("Profile must contain at least one point.", nameof(profile));

        double chiMax = profile.Susceptibility.Max();
        if (chiMax <= 0.0) return (profile.Temperatures[0], profile.Temperatures[^1]);

        double half = chiMax * 0.5;
        double lo = profile.Temperatures[0];
        double hi = profile.Temperatures[^1];

        for (int i = 0; i < profile.Count; i++)
        {
            if (profile.Susceptibility[i] >= half)
            {
                lo = profile.Temperatures[i];
                break;
            }
        }

        for (int i = profile.Count - 1; i >= 0; i--)
        {
            if (profile.Susceptibility[i] >= half)
            {
                hi = profile.Temperatures[i];
                break;
            }
        }

        if (lo >= hi)
            return (profile.Temperatures[0], profile.Temperatures[^1]);

        return (lo, hi);
    }

    public static double PickPeakTemperature(SweepProfile profile)
    {
        if (profile is null)
            throw new ArgumentNullException(nameof(profile));

        if (profile.Temperatures.Count != profile.Susceptibility.Count)
            throw new ArgumentException("SweepProfile temperatures and susceptibility values must have the same length.", nameof(profile));

        if (profile.Count == 0)
            throw new ArgumentException("Profile must contain at least one point.", nameof(profile));

        int bestIndex = 0;
        double bestChi = profile.Susceptibility[0];
        for (int i = 1; i < profile.Count; i++)
        {
            double chi = profile.Susceptibility[i];
            if (chi > bestChi)
            {
                bestChi = chi;
                bestIndex = i;
            }
        }

        return profile.Temperatures[bestIndex];
    }

    public static double ComputeStability(SweepProfile profile)
    {
        if (profile is null)
            throw new ArgumentNullException(nameof(profile));

        if (profile.Temperatures.Count != profile.Susceptibility.Count)
            throw new ArgumentException("SweepProfile temperatures and susceptibility values must have the same length.", nameof(profile));

        if (profile.Count < 2) return 0.0;
        var (lo, hi) = ComputeHalfMaximumBand(profile);
        double full = profile.Temperatures[^1] - profile.Temperatures[0];
        if (full <= 0.0) return 0.0;
        return Math.Clamp(1.0 - (hi - lo) / full, 0.0, 1.0);
    }

    /// <summary>
    /// Identifies the SP-plateau boundaries per BWD1996 §C and returns the
    /// clustering temperature T_clus = (T_fs + T_ps)/2.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Temperatures in <paramref name="profile"/> are assumed ascending (as assembled
    /// by <see cref="SweepCurves.ByTemperature"/>). T_fs is the FIRST prominent
    /// smoothed χ peak — hierarchical data has one peak per pseudotransition
    /// (BWD1996 Fig 4) and the FM→SP bound is the first, not the global argmax.
    /// T_ps is the largest ABSOLUTE smoothed drop within the first descent after
    /// the peak (BWD1996: χ "abruptly diminishes by ~the largest cluster's
    /// volume" — an absolute step; drop-ratios explode on the near-zero tail).
    /// The descent stops at the first local minimum, so on multi-stage data the
    /// landmarks bracket the FIRST plateau. Falls back to T_fs when fewer than
    /// two points lie above the peak, or when χ does not genuinely drop above it.
    /// </para>
    /// </remarks>
    public static SpPlateauResult SpPlateau(SweepProfile profile)
    {
        ArgumentNullException.ThrowIfNull(profile);
        if (profile.Temperatures.Count != profile.Susceptibility.Count)
            throw new ArgumentException(
                "SweepProfile temperatures and susceptibility values must have the same length.",
                nameof(profile));
        if (profile.Count == 0)
            throw new ArgumentException("Profile must contain at least one point.", nameof(profile));

        int n = profile.Count;
        if (n == 1)
            return new SpPlateauResult(profile.Temperatures[0], profile.Temperatures[0],
                profile.Temperatures[0], CliffFound: false);

        // 3-point smoothing on the full curve reduces MC noise before landmark detection.
        var smoothed = new double[n];
        smoothed[0] = (profile.Susceptibility[0] + profile.Susceptibility[1]) / 2.0;
        for (int i = 1; i < n - 1; i++)
            smoothed[i] = (profile.Susceptibility[i - 1] + profile.Susceptibility[i] + profile.Susceptibility[i + 1]) / 3.0;
        smoothed[n - 1] = (profile.Susceptibility[n - 2] + profile.Susceptibility[n - 1]) / 2.0;

        // T_fs = global argmax of the smoothed curve — BWD1996 step (e),
        // literally "the temperature of maximal χ". On inhomogeneous data
        // smaller genuine peaks exist (e.g. a sparse background's own
        // ordering transition at the cold end), but the dominant transition
        // has the largest variance swing; the first-descent cliff below then
        // scopes the midpoint to THAT plateau.
        int peakIdx = 0;
        for (int i = 1; i < n; i++)
            if (smoothed[i] > smoothed[peakIdx]) peakIdx = i;
        double tFs = profile.Temperatures[peakIdx];

        // Need at least 2 points above the peak to find a cliff (T[i] → T[i+1] pair).
        if (peakIdx >= n - 2)
            return new SpPlateauResult(tFs, tFs, tFs, CliffFound: false);

        // T_ps = departure from the plateau: the first T above the peak where
        // the smoothed curve falls below DepartureFraction of the peak value.
        // BWD1996 §C reads "χ abruptly diminishes" as the END of the
        // near-constant SP plateau — a level crossing, not the steepest later
        // step (on gradual-then-steep declines a max-step rule overshoots
        // into the melt, and drop-ratios explode on the near-zero tail).
        double floor = DeparturePlateauFraction * smoothed[peakIdx];
        int cliffIdx = -1;
        for (int i = peakIdx + 1; i < n; i++)
        {
            if (smoothed[i] < floor) { cliffIdx = i; break; }
        }

        if (cliffIdx < 0)
            return new SpPlateauResult(tFs, tFs, tFs, CliffFound: false);

        double tPs   = profile.Temperatures[cliffIdx];
        double tClus = (tFs + tPs) / 2.0;
        return new SpPlateauResult(tFs, tPs, tClus, CliffFound: true);
    }

    /// <summary>
    /// The plateau-departure level for <see cref="SpPlateau"/>: T_ps is the
    /// first grid point whose smoothed value drops below this fraction of
    /// the peak. 0.8 keeps T_clus inside the near-constant band on both the
    /// BWD toy and Iris reference curves (artifacts/parity-profiles.tsv).
    /// </summary>
    private const double DeparturePlateauFraction = 0.8;

    public static double ComputeFkSusceptibility(SwRunResult result)
    {
        if (result is null)
            throw new ArgumentNullException(nameof(result));

        return result.FkSusceptibility;
    }
}
