using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A sparse unit direction θ over the continuous predictors (E36): its active variables 𝒥, distinct and ascending, and their nonzero
/// weights θ_𝒥 with ‖θ_𝒥‖² = 1 to 10⁻¹⁰, a point of the coordinate subsphere the paper's direction prior is uniform on. The sign is the
/// direction's own, so θ and −θ are two directions. Directions compare by value, bit for bit, so a term holding one interns as any term.
/// </summary>
public sealed class ProjectionDirection : IEquatable<ProjectionDirection>
{
    private readonly int[] _variables;
    private readonly double[] _weights;
    private readonly int _hash;

    /// <summary>The direction with weight <paramref name="weights"/>[k] on continuous predictor <paramref name="variables"/>[k], given in
    /// any order; a variable the direction does not read is left out rather than given weight 0.</summary>
    public ProjectionDirection(IEnumerable<int> variables, IEnumerable<double> weights)
    {
        ArgumentNullException.ThrowIfNull(variables);
        ArgumentNullException.ThrowIfNull(weights);
        _variables = variables.ToArray();
        _weights = weights.ToArray();
        if (_variables.Length == 0) throw new ArgumentException("A direction needs an active variable.", nameof(variables));
        if (_weights.Length != _variables.Length) throw new ArgumentException("A direction takes one weight per active variable.", nameof(weights));
        Array.Sort(_variables, _weights);
        double squares = 0.0;
        var hash = new HashCode();
        for (int k = 0; k < _variables.Length; k++)
        {
            if (_variables[k] < 0) throw new ArgumentOutOfRangeException(nameof(variables), "Variables are indices, not negative.");
            if (k > 0 && _variables[k] == _variables[k - 1]) throw new ArgumentException($"Variable {_variables[k]} appears twice in one direction.", nameof(variables));
            if (!double.IsFinite(_weights[k]) || _weights[k] == 0.0)
                throw new ArgumentException($"The weight on variable {_variables[k]} must be finite and nonzero; leave out a variable the direction does not read.", nameof(weights));
            squares += _weights[k] * _weights[k];
            hash.Add(_variables[k]);
            hash.Add(_weights[k]);
        }
        if (!(Math.Abs(squares - 1.0) <= 1e-10)) throw new ArgumentException($"A direction is a unit vector; its squared norm is {squares:R}.", nameof(weights));
        _hash = hash.ToHashCode();
    }

    /// <summary>The active variables 𝒥, ascending, by position among the continuous predictors.</summary>
    public IReadOnlyList<int> Variables => _variables;

    /// <summary>The weights θ_𝒥, in <see cref="Variables"/>' order.</summary>
    public IReadOnlyList<double> Weights => _weights;

    /// <summary>The active count a = |𝒥|.</summary>
    public int Count => _variables.Length;

    /// <summary>θ's weight on continuous predictor <paramref name="variable"/>: 0 for one outside 𝒥.</summary>
    public double WeightOn(int variable)
    {
        int k = Array.BinarySearch(_variables, variable);
        return k >= 0 ? _weights[k] : 0.0;
    }

    /// <inheritdoc/>
    public bool Equals(ProjectionDirection? other)
        => other is not null && (ReferenceEquals(this, other)
            || (_hash == other._hash && _variables.AsSpan().SequenceEqual(other._variables) && _weights.AsSpan().SequenceEqual(other._weights)));

    /// <inheritdoc/>
    public override bool Equals(object? obj) => Equals(obj as ProjectionDirection);

    /// <inheritdoc/>
    public override int GetHashCode() => _hash;

    /// <inheritdoc/>
    public override string ToString()
    {
        var text = new StringBuilder();
        for (int k = 0; k < _variables.Length; k++)
        {
            double w = _weights[k];
            text.Append(k == 0 ? (w < 0.0 ? "−" : "") : (w < 0.0 ? " − " : " + "));
            text.Append(Math.Abs(w).ToString("R", CultureInfo.InvariantCulture)).Append("·z").Append(_variables[k]);
        }
        return text.ToString();
    }
}

/// <summary>
/// A directional hinge (θᵀz − t)₊ raised to the fit's hinge degree q (E36's term type; its brief's decision 1, the first stage): the
/// projection of the standardized continuous predictors z (<see cref="Standardization"/>, the training mean and sample deviation;
/// decision 3) onto <see cref="Direction"/>, and a knot <see cref="Knot"/> on that projection. The knot is continuous rather than an
/// index into the cases, since a knot indexed to a case would jump whenever θ moved the projections (decision 4). The reverse hinge
/// (t − θᵀz)₊ is the pair (−θ, −t). A term holds at most one, alone on its case side.
/// </summary>
public sealed record DirectionalHinge
{
    /// <summary>The hinge on <paramref name="direction"/> with knot <paramref name="knot"/>, in standardized projection units.</summary>
    public DirectionalHinge(ProjectionDirection direction, double knot)
    {
        Direction = direction ?? throw new ArgumentNullException(nameof(direction));
        if (!double.IsFinite(knot)) throw new ArgumentOutOfRangeException(nameof(knot), "The knot is a finite projection value.");
        Knot = knot;
    }

    /// <summary>The direction θ.</summary>
    public ProjectionDirection Direction { get; }

    /// <summary>The knot t, on the projection θᵀz.</summary>
    public double Knot { get; }

    /// <inheritdoc/>
    public override string ToString()
        => $"({Direction} {(Knot < 0.0 ? '+' : '−')} {Math.Abs(Knot).ToString("R", CultureInfo.InvariantCulture)})₊";
}
