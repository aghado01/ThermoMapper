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

    /// <summary>
    /// Lee and Whiteley's (2018) genealogy estimate of the normalizer's relative variance Var(Ẑ)/Z², from normalized
    /// weights, each particle's initial ancestor (its Eve index, <paramref name="families"/>) and the number of
    /// resampling events behind them: <c>1 − (N/(N−1))^{R+1}·Σ_{i,j: E_i ≠ E_j} W_i W_j</c>. Particles that share an
    /// initial ancestor are correlated through it, which the weights alone do not see. Between resamplings the
    /// weights carry the accumulated potentials, so the stages between two resamplings are one step of Lee and
    /// Whiteley's particle system. Times Ẑ², it is unbiased for Var(Ẑ) under multinomial resampling at fixed times
    /// with kernels that do not read the population; Du and Guyader (2021) carry its consistency to adaptive
    /// schedules. Its N/(N−1) factors are multinomial resampling's: the lower-variance schemes (systematic, stratified,
    /// residual) merge lineages more slowly, and it then reads too small. At R = 0 it is importance sampling's
    /// (N·ΣW² − 1)/(N − 1). It can be negative, and when every lineage has collapsed into one it reads 1, whatever the
    /// run: the families count says when. NaN for one particle.
    /// </summary>
    public static double GenealogyRelativeVariance(ReadOnlySpan<double> weights, ReadOnlySpan<int> families, int resamplings)
    {
        if (weights.Length != families.Length) throw new ArgumentException("The two spans must have the same length.", nameof(families));
        if (resamplings < 0) throw new ArgumentOutOfRangeException(nameof(resamplings));
        int n = weights.Length;
        if (n < 2) return double.NaN;
        var family = new double[n];   // each family's weight; Eve indices are initial positions, 0 … N−1
        double total = 0.0;
        for (int i = 0; i < n; i++)
        {
            family[families[i]] += weights[i];
            total += weights[i];
        }
        double within = 0.0;
        foreach (double s in family) within += s * s;
        double across = (total * total - within) / (total * total);   // Σ_{E_i ≠ E_j} W_i W_j, renormalized
        return 1.0 - Math.Exp((resamplings + 1) * Math.Log((double)n / (n - 1))) * across;
    }

    /// <summary>The standard error of ln Ẑ from a relative variance v: √ln(1 + v), the log-normal reading (Bérard,
    /// Del Moral and Doucet 2014), which is the delta method's √v for small v. A negative v reads 0.</summary>
    public static double LogNormalizerError(double relativeVariance)
        => double.IsNaN(relativeVariance) ? double.NaN : Math.Sqrt(Math.Log(1.0 + Math.Max(relativeVariance, 0.0)));
}
