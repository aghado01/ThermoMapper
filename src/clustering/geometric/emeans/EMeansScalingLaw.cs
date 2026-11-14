using System;
using System.Collections.Generic;

namespace Clustering.Geometric.EMeans;

/// <summary>
/// The two power laws of Bessa et al. 2026 (Table 2) that set the batch size and the peak step from the problem size:
/// least-squares fits in log space to the validation error of one corpus (DEEP), reported to transfer to four others.
/// Outside the grid the fits were made on, the laws are extrapolations; <see cref="Extrapolations"/> names each bound
/// a problem crosses.
/// </summary>
public static class EMeansScalingLaw
{
    /// <summary>Largest <c>n</c> at which the batch law applies undivided; above it the paper divides by 8.</summary>
    public const int UndividedBatchLimit = 2_000_000;

    /// <summary>
    /// Batch law: <c>B = 0.035 · k^0.13 · n^0.88 · T^0.81</c>, divided by 8 when <c>n</c> exceeds
    /// <see cref="UndividedBatchLimit"/>, rounded to the nearest power of two in log space and capped at <c>n</c>.
    /// </summary>
    public static int BatchSize(int n, int k, int epochs)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n));
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));
        if (epochs < 1) throw new ArgumentOutOfRangeException(nameof(epochs));

        double coefficient = n > UndividedBatchLimit ? 0.035 / 8.0 : 0.035;
        double raw = coefficient * Math.Pow(k, 0.13) * Math.Pow(n, 0.88) * Math.Pow(epochs, 0.81);
        double exponent = Math.Round(Math.Log2(raw), MidpointRounding.AwayFromZero);
        exponent = Math.Clamp(exponent, 0.0, 30.0);
        return Math.Min(n, 1 << (int)exponent);
    }

    /// <summary>
    /// Step law: <c>α₀ = min(0.824 · k^−0.46 · B^0.90 · n^−0.45 · T^−0.11, 1)</c>, the peak step at batch size
    /// <paramref name="batch"/>.
    /// </summary>
    public static double PeakStep(int n, int k, int epochs, int batch)
        => Math.Min(UnclippedPeakStep(n, k, epochs, batch), 1.0);

    /// <summary>The step law before its cap at 1.</summary>
    public static double UnclippedPeakStep(int n, int k, int epochs, int batch)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n));
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));
        if (epochs < 1) throw new ArgumentOutOfRangeException(nameof(epochs));
        if (batch < 1) throw new ArgumentOutOfRangeException(nameof(batch));

        return 0.824 * Math.Pow(k, -0.46) * Math.Pow(batch, 0.90) * Math.Pow(n, -0.45) * Math.Pow(epochs, -0.11);
    }

    /// <summary>
    /// The bounds of the fitting grid (the paper's Appendix B) that a problem crosses, one entry per bound; empty when
    /// the problem lies inside it. The grid: <c>n</c> from 10⁵ to 10⁸ (10⁷ for the batch law), <c>k ≥ 256</c> with
    /// <c>n/k ≥ 8</c>, step-law batches from 256 to <c>n/4</c>, and epoch budgets up to 25.
    /// </summary>
    public static IReadOnlyList<string> Extrapolations(int n, int k, int epochs, int batch)
    {
        var crossed = new List<string>();
        if (n < 100_000) crossed.Add($"n={n} below the fitted 1e5");
        if (n > 100_000_000) crossed.Add($"n={n} above the fitted 1e8");
        if (k < 256) crossed.Add($"k={k} below the fitted 256");
        if ((double)n / k < 8.0) crossed.Add($"n/k={(double)n / k:G4} below the fitted 8");
        if (batch < 256) crossed.Add($"B={batch} below the step law's fitted 256");
        if (batch > n / 4) crossed.Add($"B={batch} above the step law's fitted n/4");
        if (epochs > 25) crossed.Add($"T={epochs} above the fitted 25");
        return crossed;
    }
}
