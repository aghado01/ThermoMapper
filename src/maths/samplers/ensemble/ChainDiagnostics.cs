using System;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// Convergence diagnostics on the ensemble's reduced functionals (R̂ and ESS are undefined or meaningless on a
/// kernel's raw state — varying-dimension (k,ξ), rotation-degenerate Z, label-symmetric partitions — so they are
/// computed on the fixed-dimensional invariant reductions the run reports). Pure functional reductions with no
/// coupling to any kernel: the most reusable asset in the engine family, shared across every member.
/// </summary>
public static class ChainDiagnostics
{
    /// <summary>
    /// Gelman–Rubin potential scale reduction R̂ for one scalar functional, from per-chain moments: <c>C</c>
    /// chains of <c>n</c> draws each, given each chain's Σ and Σ² of the functional. R̂ → 1 at convergence.
    /// </summary>
    public static double RHat(double[] chainSums, double[] chainSumSquares, int n)
    {
        ArgumentNullException.ThrowIfNull(chainSums);
        ArgumentNullException.ThrowIfNull(chainSumSquares);
        int c = chainSums.Length;
        if (c < 2 || n < 2) return double.NaN;

        double grand = 0.0;
        for (int i = 0; i < c; i++) grand += chainSums[i];
        grand /= c * (double)n;

        double w = 0.0, b = 0.0;
        for (int i = 0; i < c; i++)
        {
            double mean = chainSums[i] / n;
            double within = (chainSumSquares[i] - chainSums[i] * chainSums[i] / n) / (n - 1);
            w += within;
            double dm = mean - grand;
            b += dm * dm;
        }
        w /= c;
        b *= (double)n / (c - 1);

        if (w <= 1e-300) return b <= 1e-300 ? 1.0 : double.PositiveInfinity;
        double vHat = (n - 1.0) / n * w + b / n;
        return Math.Sqrt(vHat / w);
    }

    /// <summary>
    /// Effective sample size of a scalar functional across chains (Gelman BDA3 variogram estimator with
    /// Geyer's initial-positive-sequence truncation): <c>ESS = C·N / (1 + 2 Σ ρ̂_t)</c>. Each row of
    /// <paramref name="chains"/> is one chain's draw sequence. Reports the information content behind the
    /// pooled posterior, discounting autocorrelation.
    /// </summary>
    public static double Ess(double[][] chains)
    {
        ArgumentNullException.ThrowIfNull(chains);
        int c = chains.Length;
        if (c < 1) return double.NaN;
        int n = chains[0].Length;
        if (n < 4) return double.NaN;

        var means = new double[c];
        double grand = 0.0;
        for (int j = 0; j < c; j++)
        {
            double s = 0.0;
            for (int i = 0; i < n; i++) s += chains[j][i];
            means[j] = s / n;
            grand += s;
        }
        grand /= c * (double)n;

        double w = 0.0, b = 0.0;
        for (int j = 0; j < c; j++)
        {
            double s = 0.0;
            for (int i = 0; i < n; i++) { double d = chains[j][i] - means[j]; s += d * d; }
            w += s / (n - 1);
            double dm = means[j] - grand;
            b += dm * dm;
        }
        w /= c;
        b *= (double)n / (c - 1);

        double vHat = (n - 1.0) / n * w + b / n;
        double total = c * (double)n;
        if (vHat <= 0.0) return total;

        // ρ̂_t from the across-chain variogram; truncate on the first non-positive consecutive pair.
        var rho = new double[n];
        for (int t = 1; t < n; t++)
        {
            double vt = 0.0;
            long cnt = 0;
            for (int j = 0; j < c; j++)
                for (int i = t; i < n; i++) { double d = chains[j][i] - chains[j][i - t]; vt += d * d; cnt++; }
            rho[t] = 1.0 - vt / cnt / (2.0 * vHat);
        }

        double rhoSum = 0.0;
        for (int t = 1; t + 1 < n; t += 2)
        {
            double pair = rho[t] + rho[t + 1];
            if (pair < 0.0) break;
            rhoSum += pair;
        }

        double tau = 1.0 + 2.0 * rhoSum;
        double ess = total / Math.Max(tau, 1e-6);
        return Math.Min(ess, total);
    }
}
