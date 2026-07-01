using System;
using System.Collections.Generic;
using System.Globalization;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// Configuration for <see cref="BlattPartitionStrategy"/>. Carries the
/// per-phase sampler knobs plus the pluggable detector and cut policy.
/// </summary>
public sealed class BlattPartitionConfig
{
    /// <summary>Number of Potts colors (q) for per-phase equilibria.
    /// Should match the q used to build the sweep — different q's would
    /// place transitions at different temperatures.</summary>
    public int Q { get; init; } = 20;

    /// <summary>Burn-in and measurement cycles per phase equilibrium.</summary>
    public RunBudget EquilibriumBudget { get; init; } = new(1000, 5000);

    /// <summary>
    /// What sufficient-statistics the per-phase equilibrium accumulates.
    /// Defaults to <see cref="AccumulationSpec.Currencies"/> — the canonical
    /// Blatt cut (<see cref="ThresholdSpinAgreement"/>) requires the
    /// <c>Alignments</c> currency, so both per-edge precursors are collected.
    /// </summary>
    public AccumulationSpec Accumulation { get; init; } = AccumulationSpec.Currencies;

    /// <summary>
    /// Root seed for reproducible per-phase RNG derivation. Null draws
    /// from OS entropy. Per-phase seeds are derived via
    /// <see cref="SpcSeedHelper.Derive"/> with <c>round=-1</c> and
    /// <c>replica=phaseIndex</c> to keep phase seed streams disjoint
    /// from any companion sweep's probe seeds.
    /// </summary>
    public int? BaseSeed { get; init; }

    /// <summary>
    /// Detector that scans <see cref="Profiling.SweepProfile"/> for
    /// pseudo-transitions. Defaults to
    /// <see cref="MagnetizationPeakDetector"/> with its default
    /// prominence threshold (0.1 of the χ_m range).
    /// </summary>
    public IPseudoTransitionDetector? Detector { get; init; }

    /// <summary>
    /// Cut policy applied to each per-phase equilibrium. Defaults to
    /// <see cref="ThresholdSpinAgreement"/> at θ=0.5 — the canonical
    /// Blatt 1996 friends-of-friends cut.
    /// </summary>
    public IPartitionStrategy? CutPolicy { get; init; }
}

/// <summary>
/// Hierarchical SPC partitioner in the Blatt 1996 / Blatt-Wiseman-Domany
/// 1997 picture: detect pseudo-transitions on the χ_m trajectory, treat
/// the intervals between consecutive peaks as stable super-paramagnetic
/// phases, run a fresh Tier-1 equilibrium at each phase's representative
/// temperature, and apply the friends-of-friends cut to read off the
/// per-phase partition. The ordered sequence of partitions is the
/// hierarchy.
/// </summary>
/// <remarks>
/// <para><b>Inputs.</b> A <see cref="Runtime.Scheduling.SweepResult"/>
/// from any sweep strategy — typically <see cref="FixedGridSweepStrategy"/>
/// with a dense user-supplied grid that brackets the expected transition
/// range. The strategy reads only <c>sweep.Summary.Profile</c> (for peak
/// detection) and re-samples the equilibrium at each phase representative
/// on its own terms.</para>
///
/// <para><b>Phase representative T.</b> v1 uses the geometric midpoint
/// of each interval between consecutive peaks, with the cold extremum
/// = <c>min(sweep.Temperatures)</c> below the first peak and the hot
/// extremum = <c>max(sweep.Temperatures)</c> above the last peak. The
/// geometric mean is the natural choice on log-spaced T grids and the
/// arithmetic mean degenerates to it on linear grids that don't span a
/// decade. Future variants may pick the local χ_m minimum inside each
/// interval (sharper but noisier).</para>
///
/// <para><b>Nesting check.</b> The Blatt picture predicts strict
/// nesting between consecutive levels — every cluster at a hotter level
/// is a subset of some cluster at the colder level. The strategy
/// validates this against the actual per-phase partitions and reports
/// the verdict in <see cref="PartitionHierarchy.NestingHolds"/>; the
/// boolean is informational, not a failure (an undersampled or noisy
/// sweep may produce non-nesting levels that are still useful as
/// individual partitions).</para>
/// </remarks>
public sealed class BlattPartitionStrategy : IHierarchicalPartitionStrategy
{
    private readonly BlattPartitionConfig _config;

    public BlattPartitionStrategy(BlattPartitionConfig? config = null)
    {
        _config = config ?? new BlattPartitionConfig();
    }

    /// <inheritdoc />
    public PartitionHierarchy Apply(SweepResult sweep, CsrGraph graph)
    {
        var detector = _config.Detector ?? new MagnetizationPeakDetector();
        var cut = _config.CutPolicy ?? new ThresholdSpinAgreement { Theta = 0.5 };

        var profile = sweep.Summary.Profile;
        if (profile.IsEmpty || profile.Temperatures.Count < 2)
            return new PartitionHierarchy(Array.Empty<HierarchyLevel>(), NestingHolds: true);

        double[] peakTemperatures = detector.Detect(profile);
        double tMin = profile.Temperatures[0];
        double tMax = profile.Temperatures[profile.Temperatures.Count - 1];

        // Phase representatives: cold extremum | midpoints between peaks | hot extremum.
        // No peaks → single phase at the sweep midpoint.
        double[] phaseTs = ComputePhaseTemperatures(peakTemperatures, tMin, tMax);

        var levels = new List<HierarchyLevel>(phaseTs.Length);
        for (int phaseIdx = 0; phaseIdx < phaseTs.Length; phaseIdx++)
        {
            double T = phaseTs[phaseIdx];
            int? seed = SpcSeedHelper.Derive(_config.BaseSeed, T, replica: phaseIdx, round: -1);
            SpcRunResult phaseRun = SweepKernel.RunEquilibrium(
                graph, T, _config.Q,
                _config.EquilibriumBudget,
                _config.Accumulation, seed);

            Affinities  phaseAffinities = SwCurrencies.ToAffinities(phaseRun.Accumulator);
            Alignments? phaseAlignments = phaseRun.Accumulator.SpinAgreementCount is null
                ? null
                : SwCurrencies.ToAlignments(phaseRun.Accumulator);
            Assignment partition = cut.Apply(graph, phaseAffinities, phaseAlignments);
            levels.Add(new HierarchyLevel(
                Temperature: T,
                Partition:   partition,
                Provenance:  BuildProvenance(phaseIdx, phaseTs.Length, peakTemperatures)));
        }

        bool nestingHolds = PartitionNesting.Holds(levels);

        return new PartitionHierarchy(levels, NestingHolds: nestingHolds);
    }

    /// <summary>
    /// Computes phase-representative temperatures from the detected
    /// peak list. Geometric midpoints are used so log-spaced and
    /// linear sweep grids both yield natural in-phase samples.
    /// </summary>
    private static double[] ComputePhaseTemperatures(
        double[] peakTemperatures, double tMin, double tMax)
    {
        if (peakTemperatures.Length == 0)
        {
            // No peaks detected — fall back to a single representative
            // at the sweep's geometric midpoint.
            return new[] { GeometricMidpoint(tMin, tMax) };
        }

        var sorted = (double[])peakTemperatures.Clone();
        Array.Sort(sorted);

        var phases = new double[sorted.Length + 1];
        phases[0] = GeometricMidpoint(tMin, sorted[0]);
        for (int i = 1; i < sorted.Length; i++)
            phases[i] = GeometricMidpoint(sorted[i - 1], sorted[i]);
        phases[sorted.Length] = GeometricMidpoint(sorted[sorted.Length - 1], tMax);
        return phases;
    }

    /// <summary>
    /// Geometric mean of two positive temperatures, with arithmetic
    /// fallback when either endpoint is non-positive (defensive — the
    /// SPC temperature axis should be strictly positive but the
    /// fallback keeps us out of NaN territory in degenerate sweeps).
    /// </summary>
    private static double GeometricMidpoint(double a, double b)
    {
        if (a <= 0.0 || b <= 0.0) return 0.5 * (a + b);
        return Math.Sqrt(a * b);
    }

    private static string BuildProvenance(int phaseIndex, int phaseCount, double[] peaks)
    {
        // Levels[0] is the cold-end phase; Levels[^1] is the hot-end
        // phase. Provenance documents which detected peak (if any) flanks
        // the phase on each side.
        var ci = CultureInfo.InvariantCulture;
        if (peaks.Length == 0) return "no-peaks: single phase across full sweep range";
        if (phaseIndex == 0) return $"cold of T_peak={peaks[0].ToString("G6", ci)}";
        if (phaseIndex == phaseCount - 1) return $"hot of T_peak={peaks[peaks.Length - 1].ToString("G6", ci)}";
        return $"between T_peak={peaks[phaseIndex - 1].ToString("G6", ci)} and " +
               $"T_peak={peaks[phaseIndex].ToString("G6", ci)}";
    }
}
