using System;
using Maths.Distributions;

namespace Maths.Regression.Changepoint;

/// <summary>The exact piecewise-constant posterior summaries (Hutter 2005): the segment-number posterior, its
/// MAP, the model evidence, the per-position break probabilities, and the regression curve.</summary>
/// <param name="SegmentNumberPosterior"><c>P(k | y)</c> for <c>k = 1..maxSegments</c> (index <c>k−1</c>).</param>
/// <param name="MapSegments">k̂ = argmax_k P(k | y) — the number of segments (k̂−1 changepoints).</param>
/// <param name="LogEvidence">log P(y), the model evidence (up to the k-uniform prior constant).</param>
/// <param name="BreakProbability"><c>P(break after position i | y)</c> for <c>i = 1..n−1</c> (index <c>i−1</c>).</param>
/// <param name="RegressionCurve">Posterior-mean level <c>E[f_i | y]</c> at each point (length n) — Bayes-averaged
/// over all segmentations: flat in clear segments, wiggling in uncertain ones, jumping at boundaries.</param>
/// <param name="RegressionStd">Posterior standard deviation of <c>f_i</c> (length n) — the curve's error band.</param>
public sealed record ChangepointResult(
    double[] SegmentNumberPosterior,
    int MapSegments,
    double LogEvidence,
    double[] BreakProbability,
    double[] RegressionCurve,
    double[] RegressionStd);

/// <summary>
/// Exact Bayesian regression of piecewise-constant functions (Marcus Hutter 2005) — change-point detection by
/// dynamic programming, the exact (no-MCMC) counterpart to the RJ-MCMC <c>StepBasis</c> and the "how many
/// transitions, where" answer for b₁(T) over a thermal grid. An <see cref="ISegmentModel"/> supplies the
/// single-segment evidence and posterior level (Gaussian closed-form, or Cauchy/robust by quadrature); a
/// log-space forward DP sums the per-segment evidences over every k-segmentation, a backward pass gives the break
/// probabilities, and the segment posteriors give the regression curve. Everything is exact and polynomial —
/// the DP is O(n²·K) — so for a coarse thermal grid (small n) it is cheap and MCMC-error-free. The
/// uniform-over-segmentations boundary prior supplies the Occam penalty (`1/C(n−1,k−1)`); k is uniform over
/// 1..maxSegments.
/// </summary>
public static class ExactChangepoint
{
    /// <summary>Gaussian-model convenience overload (σ² noise, N(ν, ρ²) level prior).</summary>
    public static ChangepointResult Infer(double[] y, double sigma2, double nu, double rho2, int maxSegments)
        => Infer(y, new GaussianSegmentModel(y, sigma2, nu, rho2), maxSegments);

    /// <summary>Infer with an explicit segment model (Gaussian, Cauchy/robust, …).</summary>
    public static ChangepointResult Infer(double[] y, ISegmentModel model, int maxSegments)
    {
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(model);
        int n = y.Length;
        if (n < 1) throw new ArgumentException("Need at least one data point.", nameof(y));
        if (model.Count != n) throw new ArgumentException("Model length must match y.", nameof(model));
        int k = Math.Clamp(maxSegments, 1, n);

        // Per-segment evidence + posterior level tables (upper triangle 1 ≤ a ≤ b ≤ n).
        var logD = new double[n + 1, n + 1];
        var mean = new double[n + 1, n + 1];
        var var = new double[n + 1, n + 1];
        model.Fill(logD, mean, var);

        // Forward DP: fwd[s, b] = log Σ over s-segmentations of y_{1..b} of Π D(segment).
        var fwd = Fill(k + 1, n + 1, double.NegativeInfinity);
        fwd[0, 0] = 0.0;                                                    // empty prefix: 0 segments, 0 points
        for (int b = 1; b <= n; b++) fwd[1, b] = logD[1, b];
        for (int s = 2; s <= k; s++)
            for (int b = s; b <= n; b++)
            {
                double acc = double.NegativeInfinity;
                for (int a = s - 1; a < b; a++)                              // first s−1 segments cover 1..a
                    acc = LogSumExp(acc, fwd[s - 1, a] + logD[a + 1, b]);
                fwd[s, b] = acc;
            }

        // Backward DP: bwd[s, a] = log evidence of y_{a+1..n} with s segments (a = 0..n−1).
        var bwd = Fill(k + 1, n + 1, double.NegativeInfinity);
        bwd[0, n] = 0.0;                                                    // empty suffix: 0 segments after n
        for (int a = 0; a < n; a++) bwd[1, a] = logD[a + 1, n];
        for (int s = 2; s <= k; s++)
            for (int a = 0; a <= n - s; a++)
            {
                double acc = double.NegativeInfinity;
                for (int c = a + 1; c <= n - s + 1; c++)                     // first segment a+1..c
                    acc = LogSumExp(acc, logD[a + 1, c] + bwd[s - 1, c]);
                bwd[s, a] = acc;
            }

        // P(k|y) ∝ P(k) · A[k][n] / C(n−1, k−1); P(k) uniform on 1..k ⇒ logP(k) = −log k (cancels in ratios).
        double logPk = -Math.Log(k);
        var logJoint = new double[k + 1];
        double logEvidence = double.NegativeInfinity;
        for (int s = 1; s <= k; s++)
        {
            logJoint[s] = logPk - LogBinom(n - 1, s - 1) + fwd[s, n];
            logEvidence = LogSumExp(logEvidence, logJoint[s]);
        }

        var posterior = new double[k];
        for (int s = 1; s <= k; s++) posterior[s - 1] = Math.Exp(logJoint[s] - logEvidence);
        int map = 1;
        for (int s = 2; s <= k; s++) if (posterior[s - 1] > posterior[map - 1]) map = s;

        // Break after position i: any segmentation with a boundary at i — j segments on 1..i, (s−j) on i+1..n.
        var breaks = new double[Math.Max(0, n - 1)];
        for (int i = 1; i <= n - 1; i++)
        {
            double acc = double.NegativeInfinity;
            for (int s = 2; s <= k; s++)
            {
                double pre = logPk - LogBinom(n - 1, s - 1);
                int jMax = Math.Min(i, s - 1);
                for (int j = 1; j <= jMax; j++)
                {
                    int after = s - j;
                    if (n - i < after) continue;
                    acc = LogSumExp(acc, pre + fwd[j, i] + bwd[after, i]);
                }
            }
            breaks[i - 1] = Math.Exp(acc - logEvidence);
        }

        // Regression curve E[f_i|y] and its variance: each segment [a..b] contributes its posterior level
        // (mean[a,b], var[a,b]) weighted by P(seg|y) = exp(logD + logΦ − logEvidence), range-added over [a..b].
        var diffMean = new double[n + 2];
        var diffM2 = new double[n + 2];
        for (int a = 1; a <= n; a++)
            for (int b = a; b <= n; b++)
            {
                double logPhi = double.NegativeInfinity;                    // before(j)/after(m) split context
                for (int j = 0; j < k; j++)
                {
                    double fj = fwd[j, a - 1];
                    if (double.IsNegativeInfinity(fj)) continue;
                    for (int m = 0; j + 1 + m <= k; m++)
                    {
                        double bm = bwd[m, b];
                        if (double.IsNegativeInfinity(bm)) continue;
                        logPhi = LogSumExp(logPhi, logPk - LogBinom(n - 1, j + m) + fj + bm);   // k = j+1+m
                    }
                }
                if (double.IsNegativeInfinity(logPhi)) continue;
                double w = Math.Exp(logD[a, b] + logPhi - logEvidence);     // P(segment [a..b] | y)
                if (w <= 0.0) continue;

                double mab = mean[a, b];
                diffMean[a] += w * mab; diffMean[b + 1] -= w * mab;
                double m2 = w * (var[a, b] + mab * mab);
                diffM2[a] += m2; diffM2[b + 1] -= m2;
            }

        var curve = new double[n];
        var std = new double[n];
        double runMean = 0.0, runM2 = 0.0;
        for (int i = 1; i <= n; i++)
        {
            runMean += diffMean[i];
            runM2 += diffM2[i];
            curve[i - 1] = runMean;
            std[i - 1] = Math.Sqrt(Math.Max(0.0, runM2 - runMean * runMean));
        }

        return new ChangepointResult(posterior, map, logEvidence, breaks, curve, std);
    }

    /// <summary>
    /// Rough, semi-principled Gaussian hyperparameters (Hutter §7): ν = mean(y), ρ² = data variance (segment-level
    /// spread), σ² from the robust median of squared first differences (within-segment noise, immune to the few
    /// jumps). Override for a calibrated fit.
    /// </summary>
    public static (double Sigma2, double Nu, double Rho2) DefaultHyperparameters(double[] y)
    {
        ArgumentNullException.ThrowIfNull(y);
        int n = y.Length;
        double mean = 0.0;
        foreach (double v in y) mean += v;
        mean /= n;
        double variance = 0.0;
        foreach (double v in y) variance += (v - mean) * (v - mean);
        variance = n > 1 ? variance / n : 1.0;

        double sigma2 = 1e-6;
        if (n > 1)
        {
            var sqDiff = new double[n - 1];
            for (int i = 0; i < n - 1; i++) { double d = y[i + 1] - y[i]; sqDiff[i] = d * d; }
            Array.Sort(sqDiff);
            double median = sqDiff[(n - 1) / 2];
            sigma2 = Math.Max(median / 0.9098, 1e-6);   // E[(Δy)²]=2σ²·E[χ²₁], median(χ²₁)≈0.4549 ⇒ σ²≈median/0.91
        }
        double rho2 = Math.Max(variance - sigma2, 0.1 * variance + 1e-9);
        return (sigma2, mean, rho2);
    }

    private static double[,] Fill(int rows, int cols, double value)
    {
        var m = new double[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++) m[i, j] = value;
        return m;
    }

    private static double LogSumExp(double a, double b)
    {
        if (double.IsNegativeInfinity(a)) return b;
        if (double.IsNegativeInfinity(b)) return a;
        double max = Math.Max(a, b);
        return max + Math.Log(Math.Exp(a - max) + Math.Exp(b - max));
    }

    private static double LogBinom(int n, int kk)
        => SpecialFunctions.LogGamma(n + 1) - SpecialFunctions.LogGamma(kk + 1) - SpecialFunctions.LogGamma(n - kk + 1);
}
