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
/// intervals, by one weighted quantile that is type 7 at equal weights (<see cref="Quantile"/>). Merging is multiset
/// union. Its order of insertion means nothing; the ordered draws diagnostics need are a <see cref="ChainTrace"/>.
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
    /// The <paramref name="p"/> quantile, by Akinshin's (2023, §6) weighted type 7: with n* = (Σw)²/Σw², Kish's effective
    /// size, and h = (n* − 1)p + 1, the mean of the weighted quantile function over [(h − 1)/n*, h/n*]. Zero weights carry
    /// no mass and are dropped, and tied values pool their weight; n* counts the draws.
    /// <list type="bullet">
    /// <item>Under equal weights, unit or any other common weight, n* = n and it is exactly type 7, linear interpolation
    /// between order statistics at position <c>p(n − 1)</c>.</item>
    /// <item>It is continuous in the weights, a weight going to 0 included, so a negligible draw moves no quantile. Under
    /// unequal weights it is therefore not anchored: Q(0) and Q(1) lie inside the extremes once the end draws' weights fall
    /// below 1/n*.</item>
    /// <item>It depends on the weighted multiset alone, not on the order its draws arrived or merged in, nor on a common
    /// rescaling of the weights.</item>
    /// </list>
    /// NaN for an empty sample or one with no mass.
    /// </summary>
    public double Quantile(double p)
    {
        if (!(p >= 0.0 && p <= 1.0)) throw new ArgumentOutOfRangeException(nameof(p));
        int n = _values.Count;
        if (n == 0) return double.NaN;
        if (_weights is null) return Type7(_values.ToArray(), p);

        double common = 0.0;
        bool equal = true;
        int positive = 0;
        foreach (double w in _weights)
        {
            if (w == 0.0) continue;
            if (positive++ == 0) common = w;
            else equal &= w == common;
        }
        if (positive == 0) return double.NaN;   // no mass
        if (equal)
        {
            var kept = new double[positive];
            for (int i = 0, k = 0; i < n; i++)
                if (_weights[i] != 0.0) kept[k++] = _values[i];
            return Type7(kept, p);
        }

        // Sorted by value, then weight, so that tied weights are pooled in one order and the sums are reproducible; the
        // weights are taken relative to the largest, so their sums neither overflow nor depend on a common scale.
        var draws = new (double Value, double Weight)[n];
        double largest = 0.0;
        for (int i = 0; i < n; i++)
        {
            draws[i] = (_values[i], _weights[i]);
            largest = Math.Max(largest, _weights[i]);
        }
        Array.Sort(draws, (a, b) => a.Value != b.Value ? a.Value.CompareTo(b.Value) : a.Weight.CompareTo(b.Weight));
        var distinct = new List<(double Value, double Weight)>(n);
        double total = 0.0, squares = 0.0;
        foreach ((double value, double weight) in draws)
        {
            if (weight == 0.0) continue;
            double w = weight / largest;
            if (distinct.Count > 0 && distinct[^1].Value == value) distinct[^1] = (value, distinct[^1].Weight + w);
            else distinct.Add((value, w));
            total += w;
            squares += w * w;
        }

        // The window's uniform CDF at the cumulative weight t is F(t) = n*·t − (n* − 1)p, clamped to [0, 1]; each distinct
        // value takes the window's share of its step of the weighted quantile function.
        double effective = total * total / squares, offset = (effective - 1.0) * p;
        double quantile = 0.0, below = 0.0, previous = 0.0;
        for (int i = 0; i < distinct.Count && previous < 1.0; i++)
        {
            below += distinct[i].Weight;
            double t = i == distinct.Count - 1 ? 1.0 : below / total;
            double f = Math.Clamp(effective * t - offset, 0.0, 1.0);
            quantile += (f - previous) * distinct[i].Value;
            previous = f;
        }
        return Math.Clamp(quantile, distinct[0].Value, distinct[^1].Value);
    }

    // The type-7 quantile of `values`, sorted in place.
    private static double Type7(double[] values, double p)
    {
        Array.Sort(values);
        int n = values.Length;
        if (n == 1) return values[0];
        double pos = p * (n - 1);
        int lo = (int)Math.Floor(pos);
        double frac = pos - lo;
        return lo + 1 < n ? values[lo] * (1.0 - frac) + values[lo + 1] * frac : values[lo];
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
