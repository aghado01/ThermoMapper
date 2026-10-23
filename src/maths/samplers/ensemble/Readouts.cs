using System;
using System.Collections.Generic;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// A readout statistic that pools by merging: folding another of the same shape into this one is associative and
/// commutative, so pooling across chains or particles does not depend on order. Readouts carry additive sufficient
/// statistics only; ordered draws are a <see cref="ChainTrace"/>'s, and restart state is a checkpoint's. Nonlinear
/// reductions (a standard deviation, a quantile, a mode) read pooled readouts once, where the estimand lives.
/// </summary>
public interface IMergeable<in TSelf>
{
    /// <summary>Fold <paramref name="other"/> into this readout.</summary>
    void MergeFrom(TSelf other);
}

/// <summary>
/// Weighted sums over a fixed index set (a field on a grid, or a few scalars) and the total weight folded in, so each
/// index's weighted mean is <c>Sums[i] / TotalWeight</c>. A fold adds <c>weight × value</c> at each index it touches
/// and adds its weight to the total once, whether it touches every index (<see cref="Add"/>) or a few
/// (<see cref="AddSparse"/>). Unit weights are the MCMC case; other weights pool weighted particles.
/// </summary>
public sealed class WeightedSums : IMergeable<WeightedSums>
{
    private readonly double[] _sums;

    public WeightedSums(int length)
    {
        if (length < 0) throw new ArgumentOutOfRangeException(nameof(length));
        _sums = new double[length];
    }

    /// <summary>The number of indices.</summary>
    public int Length => _sums.Length;

    /// <summary>The weighted sums, one per index.</summary>
    public ReadOnlySpan<double> Sums => _sums;

    /// <summary>The total weight folded in: the draw count under unit weights.</summary>
    public double TotalWeight { get; private set; }

    /// <summary>The number of folds.</summary>
    public long Folds { get; private set; }

    /// <summary>Fold one draw that has a value at every index.</summary>
    public void Add(ReadOnlySpan<double> values, double weight = 1.0)
    {
        if (values.Length != _sums.Length) throw new ArgumentException("One value per index is required.", nameof(values));
        CheckWeight(weight);
        for (int i = 0; i < _sums.Length; i++) _sums[i] += weight * values[i];
        TotalWeight += weight;
        Folds++;
    }

    /// <summary>Fold one draw that adds <c>weight × values[j]</c> at <c>indices[j]</c> (repeats accumulate) and zero
    /// elsewhere; an empty draw still adds its weight.</summary>
    public void AddSparse(ReadOnlySpan<int> indices, ReadOnlySpan<double> values, double weight = 1.0)
    {
        if (indices.Length != values.Length) throw new ArgumentException("One value per index is required.", nameof(values));
        CheckWeight(weight);
        for (int j = 0; j < indices.Length; j++) _sums[indices[j]] += weight * values[j];
        TotalWeight += weight;
        Folds++;
    }

    /// <summary>Index <paramref name="i"/>'s weighted mean; NaN when nothing has been folded in.</summary>
    public double Mean(int i) => TotalWeight > 0.0 ? _sums[i] / TotalWeight : double.NaN;

    /// <summary>An independent copy.</summary>
    public WeightedSums Copy()
    {
        var c = new WeightedSums(_sums.Length) { TotalWeight = TotalWeight, Folds = Folds };
        _sums.CopyTo(c._sums, 0);
        return c;
    }

    public void MergeFrom(WeightedSums other)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (other._sums.Length != _sums.Length) throw new ArgumentException("Readouts of different shapes do not merge.", nameof(other));
        for (int i = 0; i < _sums.Length; i++) _sums[i] += other._sums[i];
        TotalWeight += other.TotalWeight;
        Folds += other.Folds;
    }

    private static void CheckWeight(double weight)
    {
        if (!(weight >= 0.0) || double.IsPositiveInfinity(weight))
            throw new ArgumentOutOfRangeException(nameof(weight), "A fold weight must be finite and non-negative.");
    }
}

/// <summary>
/// A weighted multiset of scalar draws, for readouts no fixed-size statistic answers: quantiles and credible
/// intervals. Merging is multiset union. Its order of insertion means nothing; the ordered draws diagnostics need
/// are a <see cref="ChainTrace"/>.
/// </summary>
public sealed class DrawSample : IMergeable<DrawSample>
{
    private readonly List<double> _values = new();
    private List<double>? _weights;   // null while every weight is 1

    /// <summary>The draws held.</summary>
    public int Count => _values.Count;

    /// <summary>The total weight: the draw count under unit weights.</summary>
    public double TotalWeight { get; private set; }

    /// <summary>An independent copy.</summary>
    public DrawSample Copy()
    {
        var c = new DrawSample { TotalWeight = TotalWeight };
        c._values.AddRange(_values);
        if (_weights is not null) c._weights = new List<double>(_weights);
        return c;
    }

    /// <summary>Fold one draw.</summary>
    public void Add(double value, double weight = 1.0)
    {
        if (!(weight >= 0.0) || double.IsPositiveInfinity(weight))
            throw new ArgumentOutOfRangeException(nameof(weight), "A fold weight must be finite and non-negative.");
        if (weight != 1.0 && _weights is null) _weights = Ones(_values.Count);
        _values.Add(value);
        _weights?.Add(weight);
        TotalWeight += weight;
    }

    public void MergeFrom(DrawSample other)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (other._weights is not null && _weights is null) _weights = Ones(_values.Count);
        _values.AddRange(other._values);
        if (_weights is not null) _weights.AddRange(other._weights ?? Ones(other._values.Count));
        TotalWeight += other.TotalWeight;
    }

    /// <summary>
    /// The <paramref name="p"/> quantile. Under unit weights it is the type-7 estimator, linear interpolation between
    /// order statistics at position <c>p(n − 1)</c>. Under general weights it is linear interpolation of the weighted ECDF
    /// between the midpoints of each distinct value's weight: zero weights carry no mass and are dropped, and tied values
    /// pool their weight, so the answer depends on the weighted multiset alone, not on the order its draws arrived or
    /// merged in. NaN for an empty sample or one with no mass.
    /// </summary>
    public double Quantile(double p)
    {
        if (!(p >= 0.0 && p <= 1.0)) throw new ArgumentOutOfRangeException(nameof(p));
        int n = _values.Count;
        if (n == 0) return double.NaN;
        if (_weights is null)
        {
            var sorted = _values.ToArray();
            Array.Sort(sorted);
            if (n == 1) return sorted[0];
            double pos = p * (n - 1);
            int lo = (int)Math.Floor(pos);
            double frac = pos - lo;
            return lo + 1 < n ? sorted[lo] * (1.0 - frac) + sorted[lo + 1] * frac : sorted[lo];
        }

        // Sorted by value, then weight, so that tied weights are pooled in one order and the sums are reproducible.
        var draws = new (double Value, double Weight)[n];
        for (int i = 0; i < n; i++) draws[i] = (_values[i], _weights[i]);
        Array.Sort(draws, (a, b) => a.Value != b.Value ? a.Value.CompareTo(b.Value) : a.Weight.CompareTo(b.Weight));
        var distinct = new List<(double Value, double Weight)>(n);
        double total = 0.0;
        foreach ((double value, double weight) in draws)
        {
            if (weight == 0.0) continue;
            if (distinct.Count > 0 && distinct[^1].Value == value) distinct[^1] = (value, distinct[^1].Weight + weight);
            else distinct.Add((value, weight));
            total += weight;
        }
        if (distinct.Count == 0) return double.NaN;   // no mass

        double target = p * total, below = 0.0, previousMid = double.NaN;
        for (int i = 0; i < distinct.Count; i++)
        {
            double mid = below + 0.5 * distinct[i].Weight;
            if (target <= mid)
            {
                if (i == 0) return distinct[0].Value;
                double t = (target - previousMid) / (mid - previousMid);
                return distinct[i - 1].Value + t * (distinct[i].Value - distinct[i - 1].Value);
            }
            previousMid = mid;
            below += distinct[i].Weight;
        }
        return distinct[^1].Value;
    }

    private static List<double> Ones(int count)
    {
        var ones = new List<double>(count);
        for (int i = 0; i < count; i++) ones.Add(1.0);
        return ones;
    }
}

/// <summary>Counts of outcomes by category (for instance, why each proposal's design was or was not admitted).</summary>
public sealed class Tally<TCategory> : IMergeable<Tally<TCategory>> where TCategory : struct, Enum
{
    private static readonly TCategory[] Categories = Enum.GetValues<TCategory>();
    private readonly long[] _counts = new long[Categories.Length];

    /// <summary>Count one outcome.</summary>
    public void Add(TCategory category) => _counts[IndexOf(category)]++;

    /// <summary>How many outcomes fell in <paramref name="category"/>.</summary>
    public long Count(TCategory category) => _counts[IndexOf(category)];

    /// <summary>All outcomes counted.</summary>
    public long Total
    {
        get
        {
            long t = 0;
            foreach (long c in _counts) t += c;
            return t;
        }
    }

    /// <summary>Every category with its count, in declaration order.</summary>
    public IReadOnlyList<(TCategory Category, long Count)> Entries
    {
        get
        {
            var e = new (TCategory, long)[Categories.Length];
            for (int i = 0; i < e.Length; i++) e[i] = (Categories[i], _counts[i]);
            return e;
        }
    }

    public void MergeFrom(Tally<TCategory> other)
    {
        ArgumentNullException.ThrowIfNull(other);
        for (int i = 0; i < _counts.Length; i++) _counts[i] += other._counts[i];
    }

    /// <summary>An independent copy.</summary>
    public Tally<TCategory> Copy()
    {
        var c = new Tally<TCategory>();
        _counts.CopyTo(c._counts, 0);
        return c;
    }

    /// <summary>Overwrite these counts with <paramref name="other"/>'s, keeping this instance (a restore whose
    /// rebuilding evaluated states counts nothing).</summary>
    public void Assign(Tally<TCategory> other)
    {
        ArgumentNullException.ThrowIfNull(other);
        other._counts.CopyTo(_counts, 0);
    }

    private static int IndexOf(TCategory category)
    {
        int i = Array.IndexOf(Categories, category);
        if (i < 0) throw new ArgumentOutOfRangeException(nameof(category));
        return i;
    }
}
