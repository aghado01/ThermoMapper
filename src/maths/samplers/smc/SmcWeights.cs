using System;

namespace Maths.Samplers.Smc;

/// <summary>
/// Importance-weight arithmetic for a weighted particle population: log-weight normalization by log-sum-exp, the
/// effective sample size of the weights, and the conditional ESS of an incremental reweighting (Zhou, Johansen and
/// Aston 2016), which adaptive stage selection reads. Weight concentration is a property of this population's
/// importance weights; it is not the autocorrelation ESS of an independent chain, and equal weights after resampling
/// do not make the particles independent.
/// </summary>
public static class SmcWeights
{
    /// <summary>
    /// Writes the normalized weights W_i = exp(l_i)/Σ_j exp(l_j) into <paramref name="weights"/> and returns
    /// log Σ_j exp(l_j), computed about the largest l so nothing overflows. Entries of −∞ get weight 0; if every
    /// entry is −∞ the population has no support and this throws.
    /// </summary>
    public static double Normalize(ReadOnlySpan<double> logWeights, Span<double> weights)
    {
        if (weights.Length != logWeights.Length) throw new ArgumentException("The two spans must have the same length.", nameof(weights));
        double max = double.NegativeInfinity;
        foreach (double l in logWeights)
        {
            if (double.IsNaN(l) || double.IsPositiveInfinity(l)) throw new ArgumentException($"A log weight is {l}.", nameof(logWeights));
            if (l > max) max = l;
        }
        if (double.IsNegativeInfinity(max)) throw new InvalidOperationException("Every particle has zero weight: the stage density is zero on the whole population.");
        double sum = 0.0;
        for (int i = 0; i < logWeights.Length; i++)
        {
            weights[i] = Math.Exp(logWeights[i] - max);
            sum += weights[i];
        }
        for (int i = 0; i < weights.Length; i++) weights[i] /= sum;
        return max + Math.Log(sum);
    }

    /// <summary>1/Σ W_i² for normalized weights: N for equal weights, 1 when one particle holds all the mass.</summary>
    public static double EffectiveSampleSize(ReadOnlySpan<double> weights)
    {
        double s = 0.0;
        foreach (double w in weights) s += w * w;
        return 1.0 / s;
    }

    /// <summary>
    /// The conditional ESS of reweighting normalized weights W by incremental log weights g:
    /// N·(Σ W_i e^{g_i})² / Σ W_i e^{2g_i}, in [0, N]. Unlike the ESS of the updated weights it measures the
    /// reweighting step alone, whatever the weights' degeneracy before it, so it can choose the next stage whether or
    /// not the population was just resampled.
    /// </summary>
    public static double ConditionalEss(ReadOnlySpan<double> weights, ReadOnlySpan<double> logIncrements)
    {
        if (weights.Length != logIncrements.Length) throw new ArgumentException("The two spans must have the same length.", nameof(logIncrements));
        double max = double.NegativeInfinity;
        for (int i = 0; i < weights.Length; i++)
            if (weights[i] > 0.0 && logIncrements[i] > max) max = logIncrements[i];
        if (double.IsNegativeInfinity(max)) return 0.0;
        double first = 0.0, second = 0.0;
        for (int i = 0; i < weights.Length; i++)
        {
            if (weights[i] == 0.0) continue;
            double e = Math.Exp(logIncrements[i] - max);
            first += weights[i] * e;
            second += weights[i] * e * e;
        }
        return weights.Length * first * first / second;
    }

    /// <summary>log Σ W_i e^{g_i}: the log ratio of consecutive stages' normalizers the reweighting estimates.</summary>
    public static double LogNormalizerIncrement(ReadOnlySpan<double> weights, ReadOnlySpan<double> logIncrements)
    {
        double max = double.NegativeInfinity;
        for (int i = 0; i < weights.Length; i++)
            if (weights[i] > 0.0 && logIncrements[i] > max) max = logIncrements[i];
        if (double.IsNegativeInfinity(max)) return double.NegativeInfinity;
        double sum = 0.0;
        for (int i = 0; i < weights.Length; i++)
            if (weights[i] > 0.0) sum += weights[i] * Math.Exp(logIncrements[i] - max);
        return max + Math.Log(sum);
    }
}
