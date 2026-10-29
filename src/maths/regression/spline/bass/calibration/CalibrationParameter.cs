using System;
using System.Collections.Generic;
using System.Text.Json.Serialization;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// A calibration parameter θ_j (D09), declared by the emulator input it sets: continuous on a bounded box, sampled in
/// [0, 1] and reported in physical coordinates, or categorical over declared levels, which R's calibration does not
/// have. Values only, so a calibration config serializes whole.
/// </summary>
[JsonPolymorphic(TypeDiscriminatorPropertyName = "Kind")]
[JsonDerivedType(typeof(ContinuousParameter), "Continuous")]
[JsonDerivedType(typeof(CategoricalParameter), "Categorical")]
public abstract record CalibrationParameter(string Name)
{
    /// <summary>Whether <paramref name="other"/> declares the same parameter: its kind, name, bounds or levels.</summary>
    public abstract bool SameAs(CalibrationParameter other);
}

/// <summary>
/// A continuous parameter on [<paramref name="Lower"/>, <paramref name="Upper"/>] in physical units, the prior uniform
/// there (R's). The sampler moves u in [0, 1]; the physical value is Lower + u·(Upper − Lower), the transform R applies
/// before calling its emulator and leaves off what it returns (C-O1).
/// </summary>
public sealed record ContinuousParameter(string Name, double Lower, double Upper) : CalibrationParameter(Name)
{
    /// <summary>The physical value at unit coordinate <paramref name="u"/>.</summary>
    public double ToPhysical(double u) => Lower + u * (Upper - Lower);

    /// <summary>The unit coordinate of physical value <paramref name="x"/>.</summary>
    public double ToUnit(double x) => (x - Lower) / (Upper - Lower);

    /// <inheritdoc/>
    public override bool SameAs(CalibrationParameter other)
        => other is ContinuousParameter c && c.Name == Name && c.Lower.Equals(Lower) && c.Upper.Equals(Upper);
}

/// <summary>A categorical parameter over <paramref name="Levels"/>, the prior uniform over them; a level is held as its
/// index into the declared levels.</summary>
public sealed record CategoricalParameter(string Name, string[] Levels) : CalibrationParameter(Name)
{
    /// <inheritdoc/>
    public override bool SameAs(CalibrationParameter other)
        => other is CategoricalParameter c && c.Name == Name && c.Levels.AsSpan().SequenceEqual(Levels);
}

/// <summary>
/// The parameters split by kind, in declaration order within each kind: θ's continuous part, the [0, 1] coordinates the
/// sampler moves, and its categorical part, a level index per parameter. Refuses a declaration that is empty, names a
/// parameter twice or holds an invalid one.
/// </summary>
public sealed class ParameterLayout
{
    private readonly ContinuousParameter[] _continuous;
    private readonly CategoricalParameter[] _categorical;

    public ParameterLayout(IReadOnlyList<CalibrationParameter> parameters)
    {
        ArgumentNullException.ThrowIfNull(parameters);
        if (parameters.Count == 0) throw new ArgumentException("Calibration needs at least one parameter.", nameof(parameters));
        var names = new HashSet<string>(StringComparer.Ordinal);
        var continuous = new List<ContinuousParameter>();
        var categorical = new List<CategoricalParameter>();
        long combinations = 1;
        foreach (CalibrationParameter? p in parameters)
        {
            if (p is null) throw new ArgumentException("A parameter is null.", nameof(parameters));
            if (string.IsNullOrEmpty(p.Name)) throw new ArgumentException("A parameter needs a name.", nameof(parameters));
            if (!names.Add(p.Name)) throw new ArgumentException($"The parameter '{p.Name}' is declared twice.", nameof(parameters));
            switch (p)
            {
                case ContinuousParameter c:
                    if (!(double.IsFinite(c.Lower) && double.IsFinite(c.Upper) && c.Lower < c.Upper))
                        throw new ArgumentException($"The parameter '{c.Name}' needs finite bounds with Lower below Upper.", nameof(parameters));
                    continuous.Add(c);
                    break;
                case CategoricalParameter f:
                    if (f.Levels is null || f.Levels.Length < 2)
                        throw new ArgumentException($"The categorical parameter '{f.Name}' needs at least two levels.", nameof(parameters));
                    var levels = new HashSet<string>(StringComparer.Ordinal);
                    foreach (string? level in f.Levels)
                        if (level is null || !levels.Add(level))
                            throw new ArgumentException($"The categorical parameter '{f.Name}' has a null or repeated level.", nameof(parameters));
                    combinations = combinations > long.MaxValue / f.Levels.Length ? long.MaxValue : combinations * f.Levels.Length;
                    categorical.Add(f);
                    break;
                default:
                    throw new ArgumentException($"The parameter '{p.Name}' is of no known kind.", nameof(parameters));
            }
        }
        _continuous = [.. continuous];
        _categorical = [.. categorical];
        Combinations = categorical.Count == 0 ? 1 : combinations;
    }

    /// <summary>The continuous parameters, in declaration order.</summary>
    public IReadOnlyList<ContinuousParameter> Continuous => _continuous;

    /// <summary>The categorical parameters, in declaration order.</summary>
    public IReadOnlyList<CategoricalParameter> Categorical => _categorical;

    /// <summary>The joint level combinations, Π_j L_j (saturating at <see cref="long.MaxValue"/>); 1 without categorical
    /// parameters.</summary>
    public long Combinations { get; }

    /// <summary>Combination <paramref name="index"/>'s level indices, the first parameter varying slowest, into
    /// <paramref name="levels"/>.</summary>
    public void Levels(long index, Span<int> levels)
    {
        if ((ulong)index >= (ulong)Combinations) throw new ArgumentOutOfRangeException(nameof(index));
        for (int j = _categorical.Length - 1; j >= 0; j--)
        {
            int count = _categorical[j].Levels.Length;
            levels[j] = (int)(index % count);
            index /= count;
        }
    }

    /// <summary>The physical values of unit coordinates <paramref name="unit"/>, into <paramref name="physical"/>.</summary>
    public void ToPhysical(ReadOnlySpan<double> unit, Span<double> physical)
    {
        for (int i = 0; i < _continuous.Length; i++) physical[i] = _continuous[i].ToPhysical(unit[i]);
    }

    /// <summary>Whether <paramref name="parameters"/> declares the same parameters in the same order.</summary>
    public static bool Same(IReadOnlyList<CalibrationParameter> parameters, IReadOnlyList<CalibrationParameter> other)
    {
        if (parameters.Count != other.Count) return false;
        for (int i = 0; i < parameters.Count; i++)
            if (!parameters[i].SameAs(other[i])) return false;
        return true;
    }
}
