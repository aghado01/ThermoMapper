using System;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// Convergence diagnostics on the ensemble's reduced functionals (R̂ and ESS are undefined or meaningless on a
/// kernel's raw state — varying-dimension (k,ξ), rotation-degenerate Z, label-symmetric partitions — so they are
/// computed on the fixed-dimensional invariant reductions the run reports). Pure functional reductions with no
/// coupling to any kernel: the most reusable asset in the engine family, shared across every member.
/// </summary>
public static partial class ChainDiagnostics
{
    /// <summary>
    /// Gelman–Rubin R̂ for one scalar functional from each chain's streaming moments (all chains must hold the same
    /// number of draws). See <see cref="RHat(ReadOnlySpan{double}, ReadOnlySpan{double}, double)"/> for the states.
    /// </summary>
    public static RHatDiagnostic RHat(ReadOnlySpan<RunningMoments> chains)
    {
        int c = chains.Length;
        if (c < 2) return RHatDiagnostic.Unavailable;
        long n = chains[0].Count;
        var means = new double[c];
        var variances = new double[c];
        for (int i = 0; i < c; i++)
        {
            if (chains[i].Count != n)
                throw new ArgumentException("Every chain must hold the same number of draws.", nameof(chains));
            means[i] = chains[i].Mean;
            variances[i] = chains[i].Variance;
        }
        return RHat(means, variances, n);
    }

    /// <summary>
    /// Gelman–Rubin potential scale reduction R̂ for one scalar functional from per-chain means and sample variances
    /// over <paramref name="drawsPerChain"/> draws each; R̂ → 1 at convergence. Zero within-chain variance in every
    /// chain gives <see cref="DiagnosticState.Constant"/> when the chain means coincide and <see cref="DiagnosticState.Frozen"/>
    /// when they differ. Fewer than two chains or draws, or a non-finite input, gives
    /// <see cref="DiagnosticState.Unavailable"/>.
    /// </summary>
    public static RHatDiagnostic RHat(ReadOnlySpan<double> chainMeans, ReadOnlySpan<double> chainVariances,
                                      double drawsPerChain)
    {
        int c = chainMeans.Length;
        if (chainVariances.Length != c)
            throw new ArgumentException("One variance per chain mean is required.", nameof(chainVariances));
        if (c < 2 || !(drawsPerChain >= 2.0)) return RHatDiagnostic.Unavailable;
        for (int i = 0; i < c; i++)
            if (!double.IsFinite(chainMeans[i]) || !double.IsFinite(chainVariances[i]) || chainVariances[i] < 0.0)
                return RHatDiagnostic.Unavailable;

        double grand = 0.0, w = 0.0;
        for (int i = 0; i < c; i++) { grand += chainMeans[i]; w += chainVariances[i]; }
        grand /= c;
        w /= c;

        if (w == 0.0)
        {
            for (int i = 1; i < c; i++)
                if (chainMeans[i] != chainMeans[0]) return RHatDiagnostic.Frozen;
            return RHatDiagnostic.Constant;
        }

        double b = 0.0;
        for (int i = 0; i < c; i++) { double d = chainMeans[i] - grand; b += d * d; }
        b *= drawsPerChain / (c - 1);
        double vHat = (drawsPerChain - 1.0) / drawsPerChain * w + b / drawsPerChain;
        return RHatDiagnostic.Estimate(Math.Sqrt(vHat / w));
    }

}
