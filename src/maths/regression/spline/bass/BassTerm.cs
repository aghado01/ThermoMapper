using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// One hinge factor of a term, (s·(x − t))₊ raised to the fit's degree: a continuous predictor, or an output variable
/// for a term's functional part. <see cref="Knot"/> indexes the variable's knots, its distinct scaled values (D01).
/// </summary>
public readonly record struct HingeFactor(int Variable, int Sign, int Knot);

/// <summary>
/// One categorical factor of a term, the indicator 1(z ∈ D): a categorical predictor and a proper, nonempty subset D
/// of its training levels, held as sorted level indices, so a factor of any number of levels is one value.
/// </summary>
public sealed class CategoricalFactor : IEquatable<CategoricalFactor>
{
    private readonly int[] _levels;

    /// <summary>The subset <paramref name="levels"/> of categorical predictor <paramref name="variable"/>'s training levels.</summary>
    public CategoricalFactor(int variable, IEnumerable<int> levels)
    {
        ArgumentNullException.ThrowIfNull(levels);
        if (variable < 0) throw new ArgumentOutOfRangeException(nameof(variable));
        _levels = levels.ToArray();
        Array.Sort(_levels);
        if (_levels.Length == 0) throw new ArgumentException("A categorical factor needs a level.", nameof(levels));
        if (_levels[0] < 0) throw new ArgumentOutOfRangeException(nameof(levels), "A level index must not be negative.");
        for (int k = 1; k < _levels.Length; k++)
            if (_levels[k] == _levels[k - 1]) throw new ArgumentException($"Level {_levels[k]} is listed twice.", nameof(levels));
        Variable = variable;
    }

    /// <summary>The categorical predictor, by its position among the categorical predictors.</summary>
    public int Variable { get; }

    /// <summary>The subset's training-level indices, ascending.</summary>
    public IReadOnlyList<int> Levels => _levels;

    /// <summary>Whether <paramref name="level"/> lies in the subset; a level outside the training levels (−1) never does.</summary>
    public bool Contains(int level) => level >= 0 && Array.BinarySearch(_levels, level) >= 0;

    /// <inheritdoc/>
    public bool Equals(CategoricalFactor? other)
        => other is not null && Variable == other.Variable && _levels.AsSpan().SequenceEqual(other._levels);

    /// <inheritdoc/>
    public override bool Equals(object? obj) => Equals(obj as CategoricalFactor);

    /// <inheritdoc/>
    public override int GetHashCode()
    {
        var hash = new HashCode();
        hash.Add(Variable);
        foreach (int l in _levels) hash.Add(l);
        return hash.ToHashCode();
    }

    /// <inheritdoc/>
    public override string ToString() => $"z{Variable}∈{{{string.Join(",", _levels)}}}";
}

/// <summary>
/// A BASS term (D01): a continuous part of hinge factors, a categorical part of subset indicators and, for a direct
/// functional response, a functional part of hinge factors over the output variables. Each part holds distinct
/// variables, and a term holds at least one factor; a part may be empty, which leaves it the constant 1.
/// <para>A term is held in canonical form, each part sorted by variable, so it does not depend on the order its
/// factors were drawn in, and terms compare by value. That makes a model's collection of terms a set of distinct values,
/// which D01's M! counting needs, and lets a <see cref="BassTermTable"/> intern them. The term knows indices only; a
/// <see cref="BassBasis"/> checks them against the inputs and evaluates the term.</para>
/// </summary>
public sealed class BassTerm : IEquatable<BassTerm>
{
    private readonly HingeFactor[] _continuous, _functional;
    private readonly CategoricalFactor[] _categorical;
    private readonly int _hash;

    /// <summary>The term with these factors, in any order.</summary>
    public BassTerm(IEnumerable<HingeFactor>? continuous = null, IEnumerable<CategoricalFactor>? categorical = null,
                    IEnumerable<HingeFactor>? functional = null)
    {
        _continuous = Hinges(continuous, "continuous");
        _functional = Hinges(functional, "functional");
        _categorical = categorical?.ToArray() ?? [];
        foreach (CategoricalFactor? f in _categorical)
            if (f is null) throw new ArgumentException("A categorical factor is null.", nameof(categorical));
        Array.Sort(_categorical, (a, b) => a.Variable.CompareTo(b.Variable));
        for (int k = 1; k < _categorical.Length; k++)
            if (_categorical[k].Variable == _categorical[k - 1].Variable)
                throw new ArgumentException($"Categorical variable {_categorical[k].Variable} appears twice in one term.", nameof(categorical));
        if (_continuous.Length + _categorical.Length + _functional.Length == 0)
            throw new ArgumentException("A term needs a factor; the intercept is not a term.");

        var hash = new HashCode();
        foreach (HingeFactor h in _continuous) hash.Add(h);
        hash.Add(-1);
        foreach (CategoricalFactor c in _categorical) hash.Add(c);
        hash.Add(-2);
        foreach (HingeFactor h in _functional) hash.Add(h);
        _hash = hash.ToHashCode();
    }

    private static HingeFactor[] Hinges(IEnumerable<HingeFactor>? factors, string part)
    {
        HingeFactor[] h = factors?.ToArray() ?? [];
        Array.Sort(h, (a, b) => a.Variable.CompareTo(b.Variable));
        for (int k = 0; k < h.Length; k++)
        {
            if (h[k].Variable < 0 || h[k].Knot < 0) throw new ArgumentOutOfRangeException(part, "Variables and knots are indices, not negative.");
            if (h[k].Sign is not (1 or -1)) throw new ArgumentOutOfRangeException(part, "A hinge's sign is +1 or −1.");
            if (k > 0 && h[k].Variable == h[k - 1].Variable)
                throw new ArgumentException($"The {part} variable {h[k].Variable} appears twice in one term.", part);
        }
        return h;
    }

    /// <summary>The continuous part, by variable; its count is the term's continuous order.</summary>
    public IReadOnlyList<HingeFactor> Continuous => _continuous;

    /// <summary>The categorical part, by variable.</summary>
    public IReadOnlyList<CategoricalFactor> Categorical => _categorical;

    /// <summary>The functional part, by output variable.</summary>
    public IReadOnlyList<HingeFactor> Functional => _functional;

    /// <summary>Whether the term has a factor on the case side, continuous or categorical.</summary>
    public bool HasCaseFactor => _continuous.Length + _categorical.Length > 0;

    /// <inheritdoc/>
    public bool Equals(BassTerm? other)
        => other is not null && (ReferenceEquals(this, other)
            || (_hash == other._hash && _continuous.AsSpan().SequenceEqual(other._continuous)
                && _functional.AsSpan().SequenceEqual(other._functional) && _categorical.SequenceEqual(other._categorical)));

    /// <inheritdoc/>
    public override bool Equals(object? obj) => Equals(obj as BassTerm);

    /// <inheritdoc/>
    public override int GetHashCode() => _hash;

    /// <inheritdoc/>
    public override string ToString()
    {
        var parts = new List<string>();
        parts.AddRange(_continuous.Select(h => $"x{h.Variable}{(h.Sign > 0 ? '+' : '−')}k{h.Knot}"));
        parts.AddRange(_categorical.Select(c => c.ToString()));
        parts.AddRange(_functional.Select(h => $"t{h.Variable}{(h.Sign > 0 ? '+' : '−')}k{h.Knot}"));
        return string.Join(" · ", parts);
    }
}

/// <summary>
/// Interned terms: each distinct term gets one id, in the order first seen, and a structure holds ids, so a retained
/// draw never aliases mutable state. One writer at a time: a run keeps one table per chain and merges them in chain
/// order at the end, so the ids do not depend on how the chains were scheduled.
/// </summary>
public sealed class BassTermTable
{
    private readonly List<BassTerm> _terms = [];
    private readonly Dictionary<BassTerm, int> _ids = [];

    /// <summary>The terms interned so far.</summary>
    public int Count => _terms.Count;

    /// <summary>The term with id <paramref name="id"/>.</summary>
    public BassTerm this[int id] => _terms[id];

    /// <summary>The id of <paramref name="term"/>, interning it if it is new.</summary>
    public int Intern(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        if (_ids.TryGetValue(term, out int id)) return id;
        id = _terms.Count;
        _terms.Add(term);
        _ids.Add(term, id);
        return id;
    }

    /// <summary>The id of <paramref name="term"/>, if it has been interned.</summary>
    public bool TryGetId(BassTerm term, out int id) => _ids.TryGetValue(term, out id);
}
