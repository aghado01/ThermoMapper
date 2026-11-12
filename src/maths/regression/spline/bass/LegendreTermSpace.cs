using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// KHAOS's degree-split term space 𝒜_{p, d_max, q_max} (E19): every Legendre term over p continuous predictors whose order q,
/// its count of variables, is at most q_max and whose degree d, the sum of its factors' α ≥ 1, is at most d_max. A term is its
/// variable set with a composition of its degree over them, so
/// |𝒜| = Σ_{q ≤ q_max} Σ_{d ≤ d_max} C(p, q)·C(d − 1, q − 1) = Σ_{q ≤ q_max} C(p, q)·C(d_max, q), the inner sum by the hockey
/// stick. The structure prior's term measure is uniform over it, and a uniform composition of d into q positive parts, q − 1
/// distinct cuts from {1, …, d − 1}, has probability 1/C(d − 1, q − 1) (the brief's third decision, which corrects
/// <c>random_partition</c>). Immutable.
/// </summary>
public sealed class LegendreTermSpace
{
    private readonly double[] _orderWeights;   // C(p, q)·C(d_max, q) by q, the terms of order q
    private readonly long[]? _orderCounts;     // the same as integers when |𝒜| fits, for an exact draw
    private readonly long _total;              // |𝒜| as an integer when it fits; 0 otherwise

    /// <summary>The space over <paramref name="variables"/> predictors with degree at most <paramref name="maxDegree"/> and
    /// order at most <paramref name="maxOrder"/>, which may not exceed either.</summary>
    public LegendreTermSpace(int variables, int maxDegree, int maxOrder)
    {
        if (variables < 1) throw new ArgumentOutOfRangeException(nameof(variables), "A Legendre basis needs a continuous predictor.");
        if (maxOrder < 1 || maxOrder > variables) throw new ArgumentOutOfRangeException(nameof(maxOrder), $"q_max lies in [1, {variables}], the predictors.");
        if (maxDegree < maxOrder) throw new ArgumentOutOfRangeException(nameof(maxDegree), "d_max must reach q_max, or the top orders hold no term.");
        Variables = variables;
        MaxDegree = maxDegree;
        MaxOrder = maxOrder;
        _orderWeights = new double[maxOrder + 1];
        var counts = new long[maxOrder + 1];
        var logs = new double[maxOrder];
        bool exact = true;
        long total = 0;
        for (int q = 1; q <= maxOrder; q++)
        {
            logs[q - 1] = LogChoose(variables, q) + LogChoose(maxDegree, q);
            _orderWeights[q] = Math.Exp(logs[q - 1]);
            if (!exact) continue;
            bool fits = Choose(variables, q, out long a) & Choose(maxDegree, q, out long b);   // both, so both are assigned
            exact = fits && a <= (long.MaxValue >> 2) / b && a * b <= (long.MaxValue >> 2) - total;
            if (!exact) continue;
            counts[q] = a * b;
            total += counts[q];
        }
        _orderCounts = exact ? counts : null;
        _total = exact ? total : 0;
        double top = double.NegativeInfinity;
        foreach (double l in logs) top = Math.Max(top, l);
        double sum = 0.0;
        foreach (double l in logs) sum += Math.Exp(l - top);
        LogSize = exact ? Math.Log(total) : top + Math.Log(sum);
        Size = exact ? total : Math.Exp(LogSize);
    }

    /// <summary>p, the continuous predictors.</summary>
    public int Variables { get; }

    /// <summary>d_max, the largest degree.</summary>
    public int MaxDegree { get; }

    /// <summary>q_max, the largest order.</summary>
    public int MaxOrder { get; }

    /// <summary>|𝒜|, the terms: exact below 2⁵³, rounded above, where the draw still counts in integers while |𝒜| fits a long.</summary>
    public double Size { get; }

    /// <summary>ln |𝒜|: the uniform term measure's log probability is its negative.</summary>
    public double LogSize { get; }

    /// <summary>Refuses a term outside the space: another part, a variable beyond p, an order above q_max, a degree above d_max.</summary>
    public void Validate(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        if (term.Continuous.Count + term.Categorical.Count + term.Functional.Count > 0)
            throw new ArgumentException("A Legendre fit's terms hold Legendre factors alone.", nameof(term));
        int q = term.Legendre.Count;
        if (q == 0 || q > MaxOrder) throw new ArgumentException($"The term's order {q} lies outside [1, {MaxOrder}].", nameof(term));
        foreach (LegendreFactor f in term.Legendre)
            if (f.Variable >= Variables) throw new ArgumentException($"Continuous variable {f.Variable} does not exist.", nameof(term));
        if (term.LegendreDegree > MaxDegree) throw new ArgumentException($"The term's degree {term.LegendreDegree} exceeds d_max = {MaxDegree}.", nameof(term));
    }

    /// <summary>A term drawn uniformly from the space: its order q with probability C(p, q)·C(d_max, q)/|𝒜|, its variables a
    /// uniform q-subset, and its degrees a uniform q-subset s₁ &lt; … &lt; s_q of {1, …, d_max} read as α₁ = s₁, αₖ = sₖ − sₖ₋₁,
    /// which is uniform over the compositions of every d ≤ d_max into q parts.</summary>
    public BassTerm Draw(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        int q = DrawOrder(rng);
        int[] variables = SortedSubset(rng, Variables, q, 0);
        int[] sums = SortedSubset(rng, MaxDegree, q, 1);
        var factors = new LegendreFactor[q];
        for (int k = 0, previous = 0; k < q; k++)
        {
            factors[k] = new LegendreFactor(variables[k], sums[k] - previous);
            previous = sums[k];
        }
        return new BassTerm(legendre: factors);
    }

    /// <summary>Every term of the space, by order, then variable set, then degrees; for enumeration on small spaces.</summary>
    public IEnumerable<BassTerm> Terms()
    {
        for (int q = 1; q <= MaxOrder; q++)
            foreach (int[] variables in Subsets(Variables, q, 0))
                foreach (int[] sums in Subsets(MaxDegree, q, 1))
                {
                    var factors = new LegendreFactor[q];
                    for (int k = 0, previous = 0; k < q; k++)
                    {
                        factors[k] = new LegendreFactor(variables[k], sums[k] - previous);
                        previous = sums[k];
                    }
                    yield return new BassTerm(legendre: factors);
                }
    }

    /// <summary>A uniform composition of <paramref name="degree"/> into <paramref name="order"/> positive parts: q − 1 distinct cuts
    /// from {1, …, d − 1}, sorted, each part the gap between consecutive cuts.</summary>
    public static int[] DrawComposition(Xoshiro256PlusPlus rng, int degree, int order)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (order < 1 || degree < order) throw new ArgumentOutOfRangeException(nameof(degree), "A composition needs 1 ≤ q ≤ d.");
        int[] cuts = SortedSubset(rng, degree - 1, order - 1, 1);
        var parts = new int[order];
        for (int k = 0, previous = 0; k < order; k++)
        {
            int next = k < order - 1 ? cuts[k] : degree;
            parts[k] = next - previous;
            previous = next;
        }
        return parts;
    }

    /// <summary>ln 1/C(d − 1, q − 1), a uniform composition's log probability.</summary>
    public static double LogComposition(int degree, int order)
    {
        if (order < 1 || degree < order) throw new ArgumentOutOfRangeException(nameof(degree), "A composition needs 1 ≤ q ≤ d.");
        return -LogChoose(degree - 1, order - 1);
    }

    /// <summary>ln C(n, k), by the log gamma function.</summary>
    internal static double LogChoose(int n, int k)
        => k < 0 || k > n ? double.NegativeInfinity
           : k == 0 || k == n ? 0.0
           : SpecialFunctions.LogGamma(n + 1.0) - SpecialFunctions.LogGamma(k + 1.0) - SpecialFunctions.LogGamma(n - k + 1.0);

    private int DrawOrder(Xoshiro256PlusPlus rng)
    {
        if (MaxOrder == 1) return 1;
        if (_orderCounts is { } counts)
        {
            long draw = BassDraws.Below(rng, _total);
            for (int q = 1; q < MaxOrder; q++)
            {
                if (draw < counts[q]) return q;
                draw -= counts[q];
            }
            return MaxOrder;
        }
        double u = rng.NextDouble() * Size;   // past a long's range: doubles, whose rounding is far below any test's reach
        for (int q = 1; q < MaxOrder; q++)
        {
            if (u < _orderWeights[q]) return q;
            u -= _orderWeights[q];
        }
        return MaxOrder;
    }

    // k distinct values of {offset, …, offset + n − 1}, uniformly among the C(n, k), sorted: a partial Fisher–Yates shuffle.
    private static int[] SortedSubset(Xoshiro256PlusPlus rng, int n, int k, int offset)
    {
        var pool = new int[n];
        for (int v = 0; v < n; v++) pool[v] = v + offset;
        for (int i = 0; i < k; i++)
        {
            int j = i + BassDraws.Below(rng, n - i);
            (pool[i], pool[j]) = (pool[j], pool[i]);
        }
        int[] chosen = pool[..k];
        Array.Sort(chosen);
        return chosen;
    }

    // Every k-subset of {offset, …, offset + n − 1}, ascending, in lexicographic order.
    private static IEnumerable<int[]> Subsets(int n, int k, int offset)
    {
        var current = new int[k];
        for (int i = 0; i < k; i++) current[i] = i;
        while (true)
        {
            var subset = new int[k];
            for (int i = 0; i < k; i++) subset[i] = current[i] + offset;
            yield return subset;
            int t = k - 1;
            while (t >= 0 && current[t] == n - k + t) t--;
            if (t < 0) yield break;
            current[t]++;
            for (int i = t + 1; i < k; i++) current[i] = current[i - 1] + 1;
        }
    }

    // C(n, k) exactly when it fits a long.
    private static bool Choose(int n, int k, out long value)
    {
        value = 1;
        k = Math.Min(k, n - k);
        for (int i = 1; i <= k; i++)
        {
            long numerator = n - k + i;
            if (value > long.MaxValue / numerator) return false;
            value = value * numerator / i;
        }
        return true;
    }
}
