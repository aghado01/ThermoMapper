using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// The data a BASS fit reads: named predictors in the caller's order, which the fit keeps, so results, sensitivity
/// indices and predictions refer to them by name and position, and one response. Mixed predictor types are
/// first-class. Every value is copied and checked at construction, so a fit, its fingerprint and its retained draws
/// never see one change under them (the contract shapes, decided 2026-10-28).
/// <para>Refused here, with a reason: no predictor; a name used twice, output variables included; a predictor whose
/// length is not the response's case count. Each predictor and response refuses its own defects. A zero or constant
/// response depends on the intercept convention, so <see cref="BassResolved"/> refuses it (D03).</para>
/// </summary>
public sealed record BassData
{
    /// <summary>The predictors, in the caller's order, and the response.</summary>
    public BassData(IReadOnlyList<BassPredictor> predictors, BassResponse response)
    {
        ArgumentNullException.ThrowIfNull(predictors);
        ArgumentNullException.ThrowIfNull(response);
        if (predictors.Count == 0) throw new ArgumentException("BASS needs at least one predictor.", nameof(predictors));
        var names = new HashSet<string>(StringComparer.Ordinal);
        var copy = new BassPredictor[predictors.Count];
        for (int j = 0; j < copy.Length; j++)
        {
            BassPredictor p = predictors[j] ?? throw new ArgumentException($"Predictor {j} is null.", nameof(predictors));
            if (p.Length != response.Cases)
                throw new ArgumentException($"Predictor '{p.Name}' has {p.Length} values for {response.Cases} cases.", nameof(predictors));
            if (!names.Add(p.Name)) throw new ArgumentException($"The name '{p.Name}' is used twice.", nameof(predictors));
            copy[j] = p;
        }
        if (response is FunctionalResponse functional)
            foreach (ContinuousPredictor v in functional.OutputVariables)
                if (!names.Add(v.Name)) throw new ArgumentException($"The name '{v.Name}' is used twice.", nameof(response));
        Predictors = Array.AsReadOnly(copy);
        Response = response;
    }

    /// <summary>The predictors, in the caller's order.</summary>
    public IReadOnlyList<BassPredictor> Predictors { get; }

    /// <summary>The response.</summary>
    public BassResponse Response { get; }
}

/// <summary>A named predictor column, continuous or categorical, one value per case.</summary>
public abstract record BassPredictor
{
    private protected BassPredictor(string name)
    {
        if (string.IsNullOrWhiteSpace(name)) throw new ArgumentException("A predictor needs a name.", nameof(name));
        Name = name;
    }

    /// <summary>The name results refer to it by.</summary>
    public string Name { get; }

    /// <summary>How many values it holds: the case count, or for an output variable the location count.</summary>
    public abstract int Length { get; }
}

/// <summary>
/// A continuous predictor, or an output variable of a <see cref="FunctionalResponse"/>. Refused: no values, a
/// non-finite value, and a constant column, which has no knot (D01's knots are its distinct values).
/// </summary>
public sealed record ContinuousPredictor : BassPredictor
{
    private readonly double[] _values;

    /// <summary>The column <paramref name="values"/>, copied.</summary>
    public ContinuousPredictor(string name, double[] values) : base(name)
    {
        ArgumentNullException.ThrowIfNull(values);
        if (values.Length == 0) throw new ArgumentException($"Predictor '{name}' has no values.", nameof(values));
        _values = (double[])values.Clone();
        double min = double.PositiveInfinity, max = double.NegativeInfinity;
        foreach (double v in _values)
        {
            if (!double.IsFinite(v)) throw new ArgumentException($"Predictor '{name}' has a non-finite value.", nameof(values));
            min = Math.Min(min, v);
            max = Math.Max(max, v);
        }
        if (!(max > min)) throw new ArgumentException($"Predictor '{name}' is constant, so it has no knot.", nameof(values));
    }

    /// <summary>The values, in case order.</summary>
    public ReadOnlyMemory<double> Values => _values;

    /// <inheritdoc/>
    public override int Length => _values.Length;
}

/// <summary>
/// A categorical predictor: each case's code indexes <see cref="Levels"/>. Refused: levels that are empty or repeated,
/// a code outside them, and a factor with fewer than two levels among its cases, which carries no information (D01).
/// A level no case takes is allowed here; the fit's training levels leave it out, and the resolved record says so.
/// </summary>
public sealed record CategoricalPredictor : BassPredictor
{
    private readonly int[] _codes;

    /// <summary>The codes <paramref name="codes"/> into <paramref name="levels"/>, both copied.</summary>
    public CategoricalPredictor(string name, int[] codes, IReadOnlyList<string> levels) : base(name)
    {
        ArgumentNullException.ThrowIfNull(codes);
        ArgumentNullException.ThrowIfNull(levels);
        if (codes.Length == 0) throw new ArgumentException($"Predictor '{name}' has no values.", nameof(codes));
        var levelCopy = new string[levels.Count];
        var seen = new HashSet<string>(StringComparer.Ordinal);
        for (int l = 0; l < levelCopy.Length; l++)
        {
            string level = levels[l];
            if (string.IsNullOrEmpty(level)) throw new ArgumentException($"Predictor '{name}' has an empty level.", nameof(levels));
            if (!seen.Add(level)) throw new ArgumentException($"Predictor '{name}' lists level '{level}' twice.", nameof(levels));
            levelCopy[l] = level;
        }
        _codes = (int[])codes.Clone();
        var observed = new bool[levelCopy.Length];
        foreach (int c in _codes)
        {
            if ((uint)c >= (uint)levelCopy.Length)
                throw new ArgumentException($"Predictor '{name}' has code {c}, outside its {levelCopy.Length} levels.", nameof(codes));
            observed[c] = true;
        }
        if (Array.FindAll(observed, o => o).Length < 2)
            throw new ArgumentException($"Predictor '{name}' takes fewer than two levels, so it carries no information.", nameof(codes));
        Levels = Array.AsReadOnly(levelCopy);
        Observed = observed;
    }

    /// <summary>The codes, in case order, each an index into <see cref="Levels"/>.</summary>
    public ReadOnlyMemory<int> Codes => _codes;

    /// <summary>The levels as declared, in the caller's order.</summary>
    public IReadOnlyList<string> Levels { get; }

    /// <summary>Whether some case takes each declared level.</summary>
    internal bool[] Observed { get; }

    /// <inheritdoc/>
    public override int Length => _codes.Length;
}

/// <summary>A BASS response: one value per case, or one curve over shared output locations per case.</summary>
public abstract record BassResponse
{
    private protected BassResponse() { }

    /// <summary>n_x: the cases.</summary>
    public abstract int Cases { get; }

    /// <summary>n_r: the output locations each case is observed at; 1 for a scalar response.</summary>
    public abstract int Locations { get; }

    /// <summary>N = n_x · n_r, the observations D03's unit information counts.</summary>
    public int Observations => Cases * Locations;

    /// <summary>Every value, case by case and, within a case, location by location.</summary>
    internal abstract ReadOnlySpan<double> Values { get; }
}

/// <summary>A scalar response, one value per case. Refused: no values, or a non-finite one.</summary>
public sealed record ScalarResponse : BassResponse
{
    private readonly double[] _y;

    /// <summary>The response <paramref name="y"/>, copied.</summary>
    public ScalarResponse(double[] y)
    {
        ArgumentNullException.ThrowIfNull(y);
        if (y.Length == 0) throw new ArgumentException("The response has no values.", nameof(y));
        foreach (double v in y)
            if (!double.IsFinite(v)) throw new ArgumentException("The response has a non-finite value.", nameof(y));
        _y = (double[])y.Clone();
    }

    /// <summary>The values, in case order.</summary>
    public ReadOnlyMemory<double> Y => _y;

    /// <inheritdoc/>
    public override int Cases => _y.Length;

    /// <inheritdoc/>
    public override int Locations => 1;

    internal override ReadOnlySpan<double> Values => _y;
}

/// <summary>
/// A direct functional response: Y is n_x × n_r, every case observed at the same n_r output locations, which need not
/// form a grid. The output variables describe the locations: their coordinates and any covariates known at each one,
/// each a further functional variable of length n_r. Refused: an empty Y, a non-finite value, no output variable, one
/// whose length is not n_r, and a name used twice among them.
/// </summary>
public sealed record FunctionalResponse : BassResponse
{
    private readonly double[] _y;   // case by case, location by location
    private readonly int _cases, _locations;

    /// <summary>The response <paramref name="y"/> (cases × locations) and its output variables, copied.</summary>
    public FunctionalResponse(double[,] y, IReadOnlyList<ContinuousPredictor> outputVariables)
    {
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(outputVariables);
        _cases = y.GetLength(0);
        _locations = y.GetLength(1);
        if (_cases == 0 || _locations == 0) throw new ArgumentException("The response has no values.", nameof(y));
        if ((long)_cases * _locations > int.MaxValue) throw new ArgumentException("The response has too many values.", nameof(y));
        _y = new double[_cases * _locations];
        for (int i = 0; i < _cases; i++)
            for (int r = 0; r < _locations; r++)
            {
                double v = y[i, r];
                if (!double.IsFinite(v)) throw new ArgumentException("The response has a non-finite value.", nameof(y));
                _y[i * _locations + r] = v;
            }
        if (outputVariables.Count == 0) throw new ArgumentException("A functional response needs an output variable.", nameof(outputVariables));
        var names = new HashSet<string>(StringComparer.Ordinal);
        var copy = new ContinuousPredictor[outputVariables.Count];
        for (int d = 0; d < copy.Length; d++)
        {
            ContinuousPredictor v = outputVariables[d] ?? throw new ArgumentException($"Output variable {d} is null.", nameof(outputVariables));
            if (v.Length != _locations)
                throw new ArgumentException($"Output variable '{v.Name}' has {v.Length} values for {_locations} locations.", nameof(outputVariables));
            if (!names.Add(v.Name)) throw new ArgumentException($"The name '{v.Name}' is used twice.", nameof(outputVariables));
            copy[d] = v;
        }
        OutputVariables = Array.AsReadOnly(copy);
    }

    /// <summary>The value of case <paramref name="case"/> at location <paramref name="location"/>.</summary>
    public double this[int @case, int location]
    {
        get
        {
            if ((uint)@case >= (uint)_cases) throw new ArgumentOutOfRangeException(nameof(@case));
            if ((uint)location >= (uint)_locations) throw new ArgumentOutOfRangeException(nameof(location));
            return _y[@case * _locations + location];
        }
    }

    /// <summary>The output variables, each of length <see cref="Locations"/>.</summary>
    public IReadOnlyList<ContinuousPredictor> OutputVariables { get; }

    /// <inheritdoc/>
    public override int Cases => _cases;

    /// <inheritdoc/>
    public override int Locations => _locations;

    internal override ReadOnlySpan<double> Values => _y;
}
