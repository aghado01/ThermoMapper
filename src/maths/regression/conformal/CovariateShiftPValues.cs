using System;
using System.Collections.Generic;
using Maths.Rng;

namespace Maths.Regression.Conformal;

/// <summary>One stream arrival's weighted smoothed p-value under covariate shift (<see cref="CovariateShiftPValues"/>): its conformity
/// <see cref="Score"/>, the weighted mass of the units so far strictly below it (<see cref="Below"/>) and tied with it, itself included
/// (<see cref="Tied"/>), the tie-breaking draw θ, and p = Below + θ·Tied. <see cref="Count"/> is the stream arrivals so far, itself
/// included; <see cref="Newest"/> is the chance the weights give the newest unit itself of holding the newest position.</summary>
public readonly record struct WeightedPValue(int Count, double Score, double Below, double Tied, double Theta, double P, double Newest);

/// <summary>
/// The smoothed conformal p-values of a stream under covariate shift with a known likelihood ratio (E35's streams under shift;
/// Tibshirani, Barber, Candès and Ramdas 2019, Lemma 3, carried into Vovk, Nouretdinov and Gammerman 2003's online p-values).
/// <para><b>The design.</b> A reference sample of m units from the source law P comes first, scored by a fixed fit. Then the stream
/// arrives from the target law Q, the inputs' law shifted and y given x unchanged. Every unit carries its likelihood ratio
/// w(x) = dQ/dP at its inputs: the reference units' as well as the stream's.</para>
/// <para><b>The p-value.</b> At the k-th stream arrival, the bag holds n = m + k units. Given the bag, and every unit after it, the newest
/// position holds unit j with chance Pⱼ = πⱼ/k. Here πⱼ = wⱼ e_(k−1)(w₋ⱼ)/e_k(w) is the chance that unit j holds one of the k stream
/// positions, the conditional Poisson inclusion probability, and e_r is the elementary symmetric polynomial of the bag's weights. The
/// newest score is ranked by those chances, p = Σ Pⱼ[αⱼ &lt; αₙ] + θ Σ Pⱼ[αⱼ = αₙ], with θ uniform on (0, 1].</para>
/// <para><b>The guarantee</b> (posed here, checked by enumeration, not yet in Lean). Suppose the units are independent, the reference
/// units from P and the stream from Q with dQ/dP = w, or more generally their joint law is weighted exchangeable with weights 1 then w
/// (Tibshirani et al.'s Definition 1), and the scores are one fixed function. Then the stream's p-values are i.i.d. uniform.
/// <list type="bullet">
/// <item>Given the bag at time n and every later unit, the units' arrangement has chance ∝ Π w over the stream positions (their
/// Lemma 3), so pₙ is uniform given them.</item>
/// <item>Every later p-value reads only that bag and the later units, so they are independent of pₙ, as in Vovk's reverse-filtration
/// proof of M-1.</item>
/// </list>
/// With equal weights it is M-1. Q may put no mass where P does: a reference unit may have w = 0. A stream unit may not, since the null
/// rules it out: such a unit lies outside the target law, itself proof of a shift beyond the stated one, and it is refused.</para>
/// <para><b>An estimated ratio</b> in place of the true one, fitted outside the stream (<see cref="DensityRatio"/>), costs validity by
/// how far the weights move, as the weighted split rule does; no estimate here certifies it.</para>
/// <para><b>Cost.</b> Each arrival reads the inclusion probabilities afresh, O(n·min(k, m)) by leave-one-out Poisson-binomial
/// probabilities. The weights are tilted so the expected count is k, every term is positive, and nothing cancels. Not thread-safe: one
/// per stream.</para>
/// </summary>
public sealed class CovariateShiftPValues
{
    private readonly List<double> _scores = [];
    private readonly List<double> _weights = [];
    private readonly Xoshiro256PlusPlus _theta;

    /// <param name="referenceScores">The reference units' conformity scores, higher conforming better, from the fit that scores the stream.</param>
    /// <param name="referenceWeights">The reference units' likelihood ratios w(x) = dQ/dP, each finite and nonnegative.</param>
    /// <param name="seed">The tie-breaking draws' seed; the draws are independent of the scores by construction.</param>
    public CovariateShiftPValues(IReadOnlyList<double> referenceScores, IReadOnlyList<double> referenceWeights, int seed)
    {
        ArgumentNullException.ThrowIfNull(referenceScores);
        ArgumentNullException.ThrowIfNull(referenceWeights);
        if (referenceScores.Count != referenceWeights.Count)
            throw new ArgumentException($"One weight per reference unit is required: {referenceWeights.Count} for {referenceScores.Count} units.", nameof(referenceWeights));
        for (int i = 0; i < referenceScores.Count; i++)
        {
            if (double.IsNaN(referenceScores[i])) throw new ArgumentException($"Reference unit {i}'s score is NaN; a NaN does not order.", nameof(referenceScores));
            if (!(referenceWeights[i] >= 0.0 && double.IsFinite(referenceWeights[i])))
                throw new ArgumentException($"Reference unit {i}'s likelihood ratio must be finite and nonnegative.", nameof(referenceWeights));
            _scores.Add(referenceScores[i]);
            _weights.Add(referenceWeights[i]);
        }
        ReferenceCount = referenceScores.Count;
        _theta = new Xoshiro256PlusPlus(seed);
    }

    /// <summary>The reference units m.</summary>
    public int ReferenceCount { get; }

    /// <summary>The stream arrivals so far.</summary>
    public int Count => _scores.Count - ReferenceCount;

    /// <summary>Ranks the newest stream unit's <paramref name="score"/> among every unit so far, each weighted by its chance of holding the
    /// newest position, and gives its weighted smoothed p-value.</summary>
    /// <param name="score">Its conformity score: finite or ±∞, never NaN.</param>
    /// <param name="weight">Its likelihood ratio w(x) = dQ/dP, positive and finite.</param>
    public WeightedPValue Add(double score, double weight)
    {
        if (double.IsNaN(score)) throw new ArgumentException("A conformity score is NaN; a NaN does not order.", nameof(score));
        if (!(weight > 0.0 && double.IsFinite(weight)))
            throw new ArgumentException("A stream unit's likelihood ratio must be positive and finite: a unit with w = 0 lies outside the target law, which the null rules out.", nameof(weight));
        _scores.Add(score);
        _weights.Add(weight);
        int n = _scores.Count, k = Count;
        double[] pi = StreamShares([.. _weights], k);
        double total = 0.0;
        foreach (double v in pi) total += v;
        double below = 0.0, tied = 0.0;
        for (int j = 0; j < n; j++)
        {
            if (_scores[j] < score) below += pi[j];
            else if (_scores[j] == score) tied += pi[j];
        }
        below /= total;
        tied /= total;
        double theta = 1.0 - _theta.NextDouble();
        // Rounded up to the least positive double where the newest unit's own share underflows: a larger p only weakens the evidence.
        return new WeightedPValue(k, score, below, tied, theta, Math.Clamp(below + theta * tied, double.Epsilon, 1.0), pi[n - 1] / total);
    }

    /// <summary>
    /// The chance πⱼ = wⱼ e_(k−1)(w₋ⱼ)/e_k(w) that each unit holds one of <paramref name="k"/> stream positions among the bag's, the
    /// conditional Poisson inclusion probabilities; they sum to k. Units of weight 0 hold none.
    /// </summary>
    /// <remarks>Tilting every weight by c leaves π unchanged, and with aᵢ = c wᵢ/(1 + c wᵢ) the count of units drawn, each with its own
    /// chance aᵢ, is Poisson-binomial with P(count = r) = Π(1 − a) c^r e_r(w). So πⱼ = aⱼ P₋ⱼ(count = k − 1)/P(count = k). c is solved
    /// so the expected count is k, which puts the probabilities read near the mode, where none underflows. Where k exceeds the
    /// reference positions, the count of units not drawn is tracked instead, so the tables run to min(k, n − k).</remarks>
    internal static double[] StreamShares(double[] weights, int k)
    {
        int n = weights.Length;
        var pi = new double[n];
        var positive = new List<int>();
        for (int j = 0; j < n; j++) if (weights[j] > 0.0) positive.Add(j);
        int np = positive.Count, reference = np - k;
        if (reference < 0) throw new InvalidOperationException("Fewer units of positive weight than stream positions.");
        if (reference == 0)
        {
            foreach (int j in positive) pi[j] = 1.0;
            return pi;
        }
        var logW = new double[np];
        for (int i = 0; i < np; i++) logW[i] = Math.Log(weights[positive[i]]);
        double t = Tilt(logW, k);
        var drawn = new double[np];    // aᵢ
        var spared = new double[np];   // 1 − aᵢ
        for (int i = 0; i < np; i++)
        {
            double x = t + logW[i];
            drawn[i] = x >= 0.0 ? 1.0 / (1.0 + Math.Exp(-x)) : Math.Exp(x) / (1.0 + Math.Exp(x));
            spared[i] = x >= 0.0 ? Math.Exp(-x) / (1.0 + Math.Exp(-x)) : 1.0 / (1.0 + Math.Exp(x));
        }
        // Count drawn units when k is the smaller side, else spared ones. Unit j drawn leaves the others' count at `need`.
        bool countDrawn = k <= reference;
        int top = countDrawn ? k : reference, need = countDrawn ? k - 1 : reference;
        double[] counted = countDrawn ? drawn : spared, other = countDrawn ? spared : drawn;
        // Suffix tables: the count over units i..np−1, to `top`.
        var suffix = new double[np + 1, top + 1];
        suffix[np, 0] = 1.0;
        for (int i = np - 1; i >= 0; i--)
            for (int r = 0; r <= top; r++)
                suffix[i, r] = other[i] * suffix[i + 1, r] + (r > 0 ? counted[i] * suffix[i + 1, r - 1] : 0.0);
        double whole = suffix[0, top];
        var prefix = new double[top + 1];
        prefix[0] = 1.0;
        for (int i = 0; i < np; i++)
        {
            double without = 0.0;
            for (int r = 0; r <= need; r++) without += prefix[r] * suffix[i + 1, need - r];
            pi[positive[i]] = drawn[i] * without / whole;
            for (int r = top; r >= 0; r--) prefix[r] = other[i] * prefix[r] + (r > 0 ? counted[i] * prefix[r - 1] : 0.0);
        }
        return pi;
    }

    // ln c so that Σ σ(ln c + ln wᵢ) = k, by Newton's method kept inside a bracket; the solve only conditions the tables, so a loose
    // root serves.
    private static double Tilt(double[] logW, int k)
    {
        double low = double.PositiveInfinity, high = double.NegativeInfinity, mean = 0.0;
        foreach (double l in logW) { low = Math.Min(low, l); high = Math.Max(high, l); mean += l / logW.Length; }
        double lo = -high - 40.0, hi = -low + 40.0;   // the count's mean is near 0 at lo and near the units' number at hi
        double t = Math.Clamp(Math.Log(k / (double)(logW.Length - k)) - mean, lo, hi);
        for (int iteration = 0; iteration < 100; iteration++)
        {
            double h = -k, slope = 0.0;
            foreach (double l in logW)
            {
                double x = t + l, a = x >= 0.0 ? 1.0 / (1.0 + Math.Exp(-x)) : Math.Exp(x) / (1.0 + Math.Exp(x));
                h += a;
                slope += a * (1.0 - a);
            }
            if (Math.Abs(h) < 1e-6) break;
            if (h < 0.0) lo = t; else hi = t;
            double next = slope > 0.0 ? t - h / slope : double.NaN;
            t = next > lo && next < hi ? next : 0.5 * (lo + hi);
        }
        return t;
    }
}
