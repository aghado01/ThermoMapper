using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// The input transform, fitted once from the training data and kept in the fit, so prediction reads inputs exactly as
/// fitting did. Continuous predictors and output variables take R's range scaling to [0, 1], stored as each column's
/// (min, max); knots are the scaled distinct values (D01). A categorical predictor keeps its training levels, by name:
/// the declared levels some case takes, in the declared order.
/// <para>At prediction a continuous value outside the training range scales past [0, 1], where the hinges extrapolate
/// linearly, and is flagged per case (<see cref="RangeScale.Contains"/>). A level outside the training levels has no
/// index (<see cref="TrainingLevels.IndexOf"/>); the prediction call says whether it is refused, the default, or placed
/// outside every subset.</para>
/// </summary>
public sealed record BassInputTransform
{
    private BassInputTransform(PredictorScale[] predictors, RangeScale[] outputVariables)
    {
        Predictors = Array.AsReadOnly(predictors);
        OutputVariables = Array.AsReadOnly(outputVariables);
    }

    /// <summary>One scale per predictor, in the data's order.</summary>
    public IReadOnlyList<PredictorScale> Predictors { get; }

    /// <summary>One range per output variable, in the response's order; none for a scalar response.</summary>
    public IReadOnlyList<RangeScale> OutputVariables { get; }

    /// <summary>The transform the training data <paramref name="data"/> defines.</summary>
    public static BassInputTransform Fit(BassData data)
    {
        ArgumentNullException.ThrowIfNull(data);
        var predictors = new PredictorScale[data.Predictors.Count];
        for (int j = 0; j < predictors.Length; j++)
            predictors[j] = data.Predictors[j] switch
            {
                ContinuousPredictor c => RangeScale.Of(c),
                CategoricalPredictor f => TrainingLevels.Of(f),
                _ => throw new ArgumentException($"Predictor '{data.Predictors[j].Name}' is of no known kind.", nameof(data)),
            };
        var outputs = data.Response is FunctionalResponse functional
            ? Array.ConvertAll([.. functional.OutputVariables], RangeScale.Of)
            : [];
        return new BassInputTransform(predictors, outputs);
    }
}

/// <summary>How one predictor's raw values become the fit's inputs.</summary>
public abstract record PredictorScale
{
    private protected PredictorScale(string name) => Name = name;

    /// <summary>The predictor's name.</summary>
    public string Name { get; }
}

/// <summary>R's range scaling, <c>(x − min)/(max − min)</c>: the training range onto [0, 1], ends exact.</summary>
public sealed record RangeScale : PredictorScale
{
    /// <summary>The range [<paramref name="min"/>, <paramref name="max"/>], with min below max.</summary>
    public RangeScale(string name, double min, double max) : base(name)
    {
        if (!double.IsFinite(min) || !double.IsFinite(max) || !(max > min))
            throw new ArgumentException($"The range of '{name}' must be finite with min below max.", nameof(max));
        Min = min;
        Max = max;
    }

    internal static RangeScale Of(ContinuousPredictor column)
    {
        double min = double.PositiveInfinity, max = double.NegativeInfinity;
        foreach (double v in column.Values.Span)
        {
            min = Math.Min(min, v);
            max = Math.Max(max, v);
        }
        return new RangeScale(column.Name, min, max);
    }

    /// <summary>The training minimum.</summary>
    public double Min { get; }

    /// <summary>The training maximum.</summary>
    public double Max { get; }

    /// <summary>The scaled value: 0 at the minimum and 1 at the maximum, exactly, and linear beyond them.</summary>
    public double Scale(double x) => (x - Min) / (Max - Min);

    /// <summary>Whether <paramref name="x"/> lies in the training range, where no extrapolation is needed.</summary>
    public bool Contains(double x) => x >= Min && x <= Max;
}

/// <summary>A categorical predictor's training levels: the declared levels some training case takes.</summary>
public sealed record TrainingLevels : PredictorScale
{
    private readonly Dictionary<string, int> _index;
    private readonly int[] _fromCode;

    private TrainingLevels(string name, string[] levels, int[] fromCode) : base(name)
    {
        Levels = Array.AsReadOnly(levels);
        _fromCode = fromCode;
        _index = new Dictionary<string, int>(levels.Length, StringComparer.Ordinal);
        for (int l = 0; l < levels.Length; l++) _index[levels[l]] = l;
    }

    internal static TrainingLevels Of(CategoricalPredictor factor)
    {
        var levels = new List<string>();
        var fromCode = new int[factor.Levels.Count];
        for (int c = 0; c < fromCode.Length; c++)
        {
            fromCode[c] = factor.Observed[c] ? levels.Count : -1;
            if (factor.Observed[c]) levels.Add(factor.Levels[c]);
        }
        return new TrainingLevels(factor.Name, [.. levels], fromCode);
    }

    /// <summary>The training levels, in the declared order.</summary>
    public IReadOnlyList<string> Levels { get; }

    /// <summary>The training index of <paramref name="level"/>, or −1 for a level no training case took.</summary>
    public int IndexOf(string level) => level is not null && _index.TryGetValue(level, out int l) ? l : -1;

    /// <summary>The training index of the training data's code <paramref name="code"/>; −1 for an unobserved level.</summary>
    internal int IndexOfCode(int code) => _fromCode[code];
}
