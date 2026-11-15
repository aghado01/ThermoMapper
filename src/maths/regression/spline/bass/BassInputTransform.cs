using System;
using System.Collections.Generic;
using Maths.LinAlg;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// The input transform, fitted once from the training data and kept in the fit, so prediction reads inputs exactly as
/// fitting did. Continuous predictors and output variables take R's range scaling to [0, 1], stored as each column's
/// (min, max): its declared domain's ends when it declares one (<see cref="InputDomain"/>), else its training range's; knots
/// are the scaled distinct values (D01). A categorical predictor keeps its training levels, by name: the declared levels some
/// case takes, in the declared order.
/// <para>At prediction a continuous value outside the training range lies past the last knot, where the hinges extrapolate
/// linearly, and is flagged per case (<see cref="RangeScale.Contains"/>); without a declared domain it scales past [0, 1].
/// A level outside the training levels has no index (<see cref="TrainingLevels.IndexOf"/>); the prediction call says
/// whether it is refused, the default, or placed outside every subset.</para>
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

/// <summary>
/// R's range scaling, <c>(x − min)/(max − min)</c>: the input's box onto [0, 1], ends exact. The box is the declared domain
/// when there is one (<see cref="InputDomain"/>), else the training range, which the scale keeps either way for its flags.
/// Beside the box it keeps the training data's <see cref="Standardization"/>, the metric E36's directions are drawn in.
/// </summary>
public sealed record RangeScale : PredictorScale
{
    /// <summary>The training range [<paramref name="min"/>, <paramref name="max"/>], with min below max, and the declared
    /// <paramref name="domain"/> that holds it, or null to scale by the range itself; and the training data's
    /// <paramref name="standardization"/>, or null for a scale built without them.</summary>
    public RangeScale(string name, double min, double max, InputDomain? domain = null, Standardization? standardization = null) : base(name)
    {
        if (!double.IsFinite(min) || !double.IsFinite(max) || !(max > min))
            throw new ArgumentException($"The range of '{name}' must be finite with min below max.", nameof(max));
        if (domain is not null && !(domain.Contains(min) && domain.Contains(max)))
            throw new ArgumentException($"The training range of '{name}' leaves its declared domain {domain}.", nameof(domain));
        TrainingMin = min;
        TrainingMax = max;
        Domain = domain;
        Min = domain?.Lower ?? min;
        Max = domain?.Upper ?? max;
        Standardization = standardization;
    }

    internal static RangeScale Of(ContinuousPredictor column)
    {
        ReadOnlySpan<double> values = column.Values.Span;
        double min = double.PositiveInfinity, max = double.NegativeInfinity;
        foreach (double v in values)
        {
            min = Math.Min(min, v);
            max = Math.Max(max, v);
        }
        // The data refuse a constant column first; were one to reach here, the range would refuse it before a deviation divided by zero.
        return new RangeScale(column.Name, min, max, column.Domain, max > min ? Standardization.Of(values) : null);
    }

    /// <summary>The box's lower end, which scales to 0: the declared domain's, or the training minimum.</summary>
    public double Min { get; }

    /// <summary>The box's upper end, which scales to 1: the declared domain's, or the training maximum.</summary>
    public double Max { get; }

    /// <summary>The training minimum.</summary>
    public double TrainingMin { get; }

    /// <summary>The training maximum.</summary>
    public double TrainingMax { get; }

    /// <summary>The declared domain the box is, or null when the box is the training range.</summary>
    public InputDomain? Domain { get; }

    /// <summary>The scaled value: 0 at the minimum and 1 at the maximum, exactly, and linear beyond them.</summary>
    public double Scale(double x) => (x - Min) / (Max - Min);

    /// <summary>Whether <paramref name="x"/> lies in the training range, where no extrapolation is needed.</summary>
    public bool Contains(double x) => x >= TrainingMin && x <= TrainingMax;

    /// <summary>The training data's mean and sample deviation, or null for a scale built without them.</summary>
    public Standardization? Standardization { get; }
}

/// <summary>
/// A continuous input's training mean and sample standard deviation, with n − 1: z = (x − <see cref="Mean"/>)/<see cref="Deviation"/>.
/// The metric E36's direction prior is uniform in (its decision 3): the paper's standardization, and <c>BayesPPR</c>'s, which uses the
/// sample deviation, so that a comparison with either is like for like. Under a uniform design it differs from the box by one common
/// factor, so the two give the same prior; under a skewed design they do not.
/// </summary>
public sealed record Standardization
{
    /// <summary>The mean <paramref name="mean"/> and deviation <paramref name="deviation"/>, finite, the deviation positive.</summary>
    public Standardization(double mean, double deviation)
    {
        if (!double.IsFinite(mean)) throw new ArgumentOutOfRangeException(nameof(mean), "The mean is finite.");
        if (!(deviation > 0.0) || double.IsPositiveInfinity(deviation)) throw new ArgumentOutOfRangeException(nameof(deviation), "The deviation is finite and positive.");
        Mean = mean;
        Deviation = deviation;
    }

    /// <summary>The training mean.</summary>
    public double Mean { get; }

    /// <summary>The training sample standard deviation, with n − 1.</summary>
    public double Deviation { get; }

    /// <summary>The standardized value of <paramref name="x"/>.</summary>
    public double Standardize(double x) => (x - Mean) / Deviation;

    // Two passes, each summed pairwise: the mean, then the squared deviations from it.
    internal static Standardization Of(ReadOnlySpan<double> values)
    {
        double mean = PairwiseSum.Of(values) / values.Length;
        var squares = new double[values.Length];
        for (int i = 0; i < squares.Length; i++)
        {
            double d = values[i] - mean;
            squares[i] = d * d;
        }
        return new Standardization(mean, Math.Sqrt(PairwiseSum.Of(squares) / (values.Length - 1)));
    }
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
