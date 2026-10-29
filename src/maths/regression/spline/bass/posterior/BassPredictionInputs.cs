using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass;

/// <summary>What prediction does with a categorical level no training case took (the contract shapes, choice 5).</summary>
public enum UnseenLevel
{
    /// <summary>Refuse the call, naming the predictor and the level: the default.</summary>
    Refuse,

    /// <summary>Place the level outside every subset, so every categorical factor on it is 0.</summary>
    OutsideEverySubset,
}

/// <summary>One named column of new values, for a predictor or an output variable.</summary>
public abstract record PredictionColumn
{
    private protected PredictionColumn(string name)
    {
        if (string.IsNullOrWhiteSpace(name)) throw new ArgumentException("A column needs a name.", nameof(name));
        Name = name;
    }

    /// <summary>The predictor or output variable it gives values for.</summary>
    public string Name { get; }

    /// <summary>The values it holds.</summary>
    public abstract int Length { get; }
}

/// <summary>New values of a continuous predictor or an output variable, in raw units; refused: none, or a non-finite
/// one. Unlike a training column it may be constant, as a single case is.</summary>
public sealed record ContinuousColumn : PredictionColumn
{
    private readonly double[] _values;

    /// <summary>The values <paramref name="values"/>, copied.</summary>
    public ContinuousColumn(string name, double[] values) : base(name)
    {
        ArgumentNullException.ThrowIfNull(values);
        if (values.Length == 0) throw new ArgumentException($"Column '{name}' has no values.", nameof(values));
        foreach (double v in values)
            if (!double.IsFinite(v)) throw new ArgumentException($"Column '{name}' has a non-finite value.", nameof(values));
        _values = (double[])values.Clone();
    }

    /// <summary>The values.</summary>
    public ReadOnlyMemory<double> Values => _values;

    /// <inheritdoc/>
    public override int Length => _values.Length;
}

/// <summary>New values of a categorical predictor, by level name (the training levels are kept by name).</summary>
public sealed record CategoricalColumn : PredictionColumn
{
    private readonly string[] _levels;

    /// <summary>Each case's level <paramref name="levels"/>, copied; refused: none, or a null one.</summary>
    public CategoricalColumn(string name, string[] levels) : base(name)
    {
        ArgumentNullException.ThrowIfNull(levels);
        if (levels.Length == 0) throw new ArgumentException($"Column '{name}' has no values.", nameof(levels));
        foreach (string level in levels)
            if (level is null) throw new ArgumentException($"Column '{name}' has a null level.", nameof(levels));
        _levels = (string[])levels.Clone();
    }

    /// <summary>Each case's level.</summary>
    public IReadOnlyList<string> Levels => _levels;

    /// <inheritdoc/>
    public override int Length => _levels.Length;
}

/// <summary>
/// The cases to predict at: one column per predictor of the fit, matched by name in any order, all of one length.
/// Refused: no column, a name given twice, columns of different lengths.
/// </summary>
public sealed record PredictionCases
{
    private readonly PredictionColumn[] _columns;

    /// <summary>The columns <paramref name="columns"/>.</summary>
    public PredictionCases(IReadOnlyList<PredictionColumn> columns)
    {
        ArgumentNullException.ThrowIfNull(columns);
        if (columns.Count == 0) throw new ArgumentException("Prediction needs a column per predictor.", nameof(columns));
        var names = new HashSet<string>(StringComparer.Ordinal);
        _columns = new PredictionColumn[columns.Count];
        for (int j = 0; j < _columns.Length; j++)
        {
            PredictionColumn c = columns[j] ?? throw new ArgumentException($"Column {j} is null.", nameof(columns));
            if (!names.Add(c.Name)) throw new ArgumentException($"The column '{c.Name}' is given twice.", nameof(columns));
            if (c.Length != columns[0].Length)
                throw new ArgumentException($"Column '{c.Name}' has {c.Length} values, column '{columns[0].Name}' {columns[0].Length}.", nameof(columns));
            _columns[j] = c;
        }
        Cases = _columns[0].Length;
    }

    /// <summary>The columns.</summary>
    public IReadOnlyList<PredictionColumn> Columns => _columns;

    /// <summary>The cases.</summary>
    public int Cases { get; }

    /// <summary>
    /// The cases' inputs through the fit's <paramref name="transform"/>, arranged as the training <paramref name="inputs"/>
    /// are, and whether each case lies outside the training range in some continuous predictor, where the hinges
    /// extrapolate linearly. An unseen level is refused or placed outside every subset, as <paramref name="unseen"/> says.
    /// Refused: a predictor with no column, a column of the wrong kind, and a column that names no predictor.
    /// </summary>
    internal (CaseInputs Inputs, bool[] OutsideRange) Prepare(BassInputTransform transform, BassInputs inputs, UnseenLevel unseen)
    {
        if (unseen is not (UnseenLevel.Refuse or UnseenLevel.OutsideEverySubset)) throw new ArgumentOutOfRangeException(nameof(unseen));
        var byName = new Dictionary<string, PredictionColumn>(StringComparer.Ordinal);
        foreach (PredictionColumn c in _columns) byName.Add(c.Name, c);
        var continuous = new double[inputs.Cases.ContinuousCount][];
        var categorical = new int[inputs.Cases.CategoricalCount][];
        var outside = new bool[Cases];
        for (int j = 0; j < transform.Predictors.Count; j++)
        {
            PredictorScale scale = transform.Predictors[j];
            if (!byName.Remove(scale.Name, out PredictionColumn? column))
                throw new ArgumentException($"No column gives the predictor '{scale.Name}'.", "cases");
            (bool isContinuous, int index) = inputs.PredictorPosition(j);
            switch ((scale, column))
            {
                case (RangeScale range, ContinuousColumn values) when isContinuous:
                {
                    ReadOnlySpan<double> raw = values.Values.Span;
                    var scaled = new double[Cases];
                    for (int i = 0; i < scaled.Length; i++)
                    {
                        scaled[i] = range.Scale(raw[i]);
                        if (!range.Contains(raw[i])) outside[i] = true;
                    }
                    continuous[index] = scaled;
                    break;
                }
                case (TrainingLevels training, CategoricalColumn levels) when !isContinuous:
                {
                    var indices = new int[Cases];
                    for (int i = 0; i < indices.Length; i++)
                    {
                        indices[i] = training.IndexOf(levels.Levels[i]);
                        if (indices[i] < 0 && unseen == UnseenLevel.Refuse)
                            throw new ArgumentException($"Case {i} has level '{levels.Levels[i]}' of '{scale.Name}', which no training case took; " +
                                                        "declare UnseenLevel.OutsideEverySubset to place it outside every subset.", "cases");
                    }
                    categorical[index] = indices;
                    break;
                }
                default:
                    throw new ArgumentException($"The column '{scale.Name}' is not of its predictor's kind.", "cases");
            }
        }
        foreach (string extra in byName.Keys) throw new ArgumentException($"The column '{extra}' names no predictor of the fit.", "cases");
        return (new CaseInputs(continuous, categorical), outside);
    }

    /// <summary>
    /// New output locations for a direct functional fit: one column per output variable, matched by name, all of one
    /// length, through the training ranges, with each location flagged when it lies outside them.
    /// </summary>
    internal static (OutputInputs Inputs, bool[] OutsideRange) PrepareLocations(BassInputTransform transform, IReadOnlyList<ContinuousColumn> locations)
    {
        if (locations.Count == 0) throw new ArgumentException("New output locations need a column per output variable.", "options");
        var byName = new Dictionary<string, ContinuousColumn>(StringComparer.Ordinal);
        foreach (ContinuousColumn c in locations)
        {
            if (c is null) throw new ArgumentException("An output-location column is null.", "options");
            if (!byName.TryAdd(c.Name, c)) throw new ArgumentException($"The output column '{c.Name}' is given twice.", "options");
            if (c.Length != locations[0].Length) throw new ArgumentException("The output-location columns must be of one length.", "options");
        }
        int m = locations[0].Length;
        var columns = new double[transform.OutputVariables.Count][];
        var outside = new bool[m];
        for (int d = 0; d < columns.Length; d++)
        {
            RangeScale range = transform.OutputVariables[d];
            if (!byName.Remove(range.Name, out ContinuousColumn? column))
                throw new ArgumentException($"No column gives the output variable '{range.Name}'.", "options");
            ReadOnlySpan<double> raw = column.Values.Span;
            columns[d] = new double[m];
            for (int r = 0; r < m; r++)
            {
                columns[d][r] = range.Scale(raw[r]);
                if (!range.Contains(raw[r])) outside[r] = true;
            }
        }
        foreach (string extra in byName.Keys) throw new ArgumentException($"The column '{extra}' names no output variable of the fit.", "options");
        return (new OutputInputs(columns), outside);
    }

    /// <summary>The training data's cases as prediction columns, by name: each continuous predictor's values and each
    /// categorical predictor's levels.</summary>
    public static PredictionCases Of(BassData data)
    {
        ArgumentNullException.ThrowIfNull(data);
        return Of(data.Predictors);
    }

    /// <summary>The cases of <paramref name="predictors"/>, as <see cref="Of(BassData)"/> reads a data set's.</summary>
    public static PredictionCases Of(IReadOnlyList<BassPredictor> predictors)
    {
        ArgumentNullException.ThrowIfNull(predictors);
        var columns = new PredictionColumn[predictors.Count];
        for (int j = 0; j < columns.Length; j++)
            columns[j] = predictors[j] switch
            {
                ContinuousPredictor c => new ContinuousColumn(c.Name, c.Values.ToArray()),
                CategoricalPredictor f => new CategoricalColumn(f.Name, Array.ConvertAll(f.Codes.ToArray(), code => f.Levels[code])),
                _ => throw new ArgumentException($"Predictor '{predictors[j]?.Name}' is of no known kind.", nameof(predictors)),
            };
        return new PredictionCases(columns);
    }
}
