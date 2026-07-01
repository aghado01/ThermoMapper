using System;

namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>
/// Canonical per-task seed derivation for SPC schedules and ad-hoc
/// task lists. Centralizes what used to be two near-identical formulas
/// (one in <see cref="Execution.SpcExecutor"/>, one in the
/// adaptive scheduler) so any task with the same logical
/// <c>(baseSeed, T, replica, round)</c> tuple resolves to the same seed
/// regardless of which strategy or executor builds it.
/// </summary>
/// <remarks>
/// <para><b>Quantization.</b> Temperatures are quantized to roughly 3
/// significant digits (factor 1000 by default) so callers that compute
/// <c>T</c> via slightly different float paths still hash to the same
/// bucket. Bump <see cref="DefaultTemperatureQuantizationFactor"/> on
/// the per-call <c>quantizationFactor</c> argument when more precision
/// is needed.</para>
///
/// <para><b>Determinism.</b> The mixer is a fixed SplitMix64 fold —
/// <em>not</em> <see cref="System.HashCode"/>, whose per-process random
/// seed silently breaks cross-process reproducibility (same
/// <c>(baseSeed, T, replica, round)</c> gave a different engine seed
/// every run; fixed 2026-06-10). The same tuple now resolves to the
/// same seed on any process, machine, and .NET version. This change is
/// itself a one-time seed-stream identity break relative to the
/// HashCode-era streams (which were never stable to begin with).</para>
/// </remarks>
public static class SpcSeedHelper
{
    /// <summary>
    /// Default temperature quantization factor (3 sig-fig granularity).
    /// </summary>
    public const int DefaultTemperatureQuantizationFactor = 1000;

    /// <summary>
    /// Derive a per-task seed from a base seed and the discriminating
    /// task coordinates. Returns <see langword="null"/> when
    /// <paramref name="baseSeed"/> is null — that path means "no
    /// reproducibility; draw from OS entropy at run time."
    /// </summary>
    /// <param name="baseSeed">Root seed shared across the whole
    /// schedule, or null for non-reproducible runs.</param>
    /// <param name="temperature">Task temperature.</param>
    /// <param name="replica">Replica index inside the (T, round)
    /// bucket. Use 0 when the workload has no replica dimension.</param>
    /// <param name="round">Refinement-round index. Use 0 when the
    /// workload has no round dimension (most fixed-grid sweeps and
    /// raw task lists).</param>
    /// <param name="quantizationFactor">Temperature quantization
    /// factor; defaults to
    /// <see cref="DefaultTemperatureQuantizationFactor"/>.</param>
    public static int? Derive(
        int? baseSeed,
        double temperature,
        int replica,
        int round = 0,
        int quantizationFactor = DefaultTemperatureQuantizationFactor)
    {
        if (baseSeed is null) return null;
        int tQuantized = (int)(temperature * quantizationFactor);

        ulong h = Mix((uint)baseSeed.Value);
        h = Mix(h ^ (uint)tQuantized);
        h = Mix(h ^ (uint)replica);
        h = Mix(h ^ (uint)round);
        return (int)(h ^ (h >> 32));
    }

    /// <summary>SplitMix64 finalizer (Steele, Lea &amp; Flood 2014) — a fixed,
    /// platform-stable 64-bit avalanche.</summary>
    private static ulong Mix(ulong z)
    {
        z += 0x9E3779B97F4A7C15UL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
        return z ^ (z >> 31);
    }
}
