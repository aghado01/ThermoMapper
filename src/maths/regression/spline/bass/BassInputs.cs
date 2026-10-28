using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// Case-side inputs a basis is evaluated at: each continuous predictor's scaled values and each categorical
/// predictor's training-level indices, one per case. A level outside the training levels is −1, which lies outside
/// every subset. The training cases are one instance; prediction's new cases are another.
/// </summary>
public sealed class CaseInputs
{
    private readonly double[][] _continuous;
    private readonly int[][] _categorical;

    /// <summary>The columns, continuous then categorical, each by its position among its kind, all of one length.</summary>
    public CaseInputs(double[][] continuous, int[][] categorical)
    {
        ArgumentNullException.ThrowIfNull(continuous);
        ArgumentNullException.ThrowIfNull(categorical);
        int cases = continuous.Length > 0 ? continuous[0].Length : categorical.Length > 0 ? categorical[0].Length : 0;
        foreach (double[] c in continuous)
            if (c is null || c.Length != cases) throw new ArgumentException("Every column needs one value per case.", nameof(continuous));
        foreach (int[] c in categorical)
            if (c is null || c.Length != cases) throw new ArgumentException("Every column needs one value per case.", nameof(categorical));
        _continuous = continuous;
        _categorical = categorical;
        Cases = cases;
    }

    /// <summary>The cases.</summary>
    public int Cases { get; }

    /// <summary>Continuous predictor <paramref name="variable"/>'s scaled values.</summary>
    public ReadOnlySpan<double> Continuous(int variable) => _continuous[variable];

    /// <summary>Categorical predictor <paramref name="variable"/>'s training-level indices.</summary>
    public ReadOnlySpan<int> Categorical(int variable) => _categorical[variable];

    /// <summary>The continuous predictors.</summary>
    public int ContinuousCount => _continuous.Length;

    /// <summary>The categorical predictors.</summary>
    public int CategoricalCount => _categorical.Length;
}

/// <summary>Output-side inputs: each output variable's scaled values at the output locations.</summary>
public sealed class OutputInputs
{
    private readonly double[][] _variables;

    /// <summary>The output variables' scaled columns, all of one length.</summary>
    public OutputInputs(double[][] variables)
    {
        ArgumentNullException.ThrowIfNull(variables);
        if (variables.Length == 0) throw new ArgumentException("Output inputs need a variable.", nameof(variables));
        int locations = variables[0]?.Length ?? 0;
        foreach (double[] v in variables)
            if (v is null || v.Length != locations) throw new ArgumentException("Every output variable needs one value per location.", nameof(variables));
        _variables = variables;
        Locations = locations;
    }

    /// <summary>The output locations.</summary>
    public int Locations { get; }

    /// <summary>Output variable <paramref name="variable"/>'s scaled values.</summary>
    public ReadOnlySpan<double> Variable(int variable) => _variables[variable];

    /// <summary>The output variables.</summary>
    public int Count => _variables.Length;
}

/// <summary>
/// The training inputs, prepared once from the data through its transform: the scaled case-side and output-side
/// columns, each variable's knots, and each categorical predictor's level count. The data's predictors split by kind
/// here; <see cref="PredictorPosition"/> maps each back to its place in the data, so results can name it.
/// <para>Knots are each variable's distinct scaled values, ascending (D01; R's <c>unique</c> on the scaled column). A
/// hinge reads only the scaled knot, so two raw values that scale to the same double are one knot.</para>
/// </summary>
public sealed class BassInputs
{
    private readonly double[][] _continuousKnots, _functionalKnots;
    private readonly int[] _levelCounts;
    private readonly (bool Continuous, int Index)[] _positions;

    private BassInputs(CaseInputs cases, OutputInputs? outputs, double[][] continuousKnots, double[][] functionalKnots,
                       int[] levelCounts, (bool, int)[] positions)
    {
        Cases = cases;
        Outputs = outputs;
        _continuousKnots = continuousKnots;
        _functionalKnots = functionalKnots;
        _levelCounts = levelCounts;
        _positions = positions;
    }

    /// <summary>The training cases' inputs.</summary>
    public CaseInputs Cases { get; }

    /// <summary>The output locations' inputs; null for a scalar response.</summary>
    public OutputInputs? Outputs { get; }

    /// <summary>Continuous predictor <paramref name="variable"/>'s knots, ascending.</summary>
    public ReadOnlySpan<double> ContinuousKnots(int variable) => _continuousKnots[variable];

    /// <summary>Output variable <paramref name="variable"/>'s knots, ascending.</summary>
    public ReadOnlySpan<double> FunctionalKnots(int variable) => _functionalKnots[variable];

    /// <summary>Categorical predictor <paramref name="variable"/>'s training levels.</summary>
    public int LevelCount(int variable) => _levelCounts[variable];

    /// <summary>Data predictor <paramref name="predictor"/>'s kind and its position among its kind.</summary>
    public (bool Continuous, int Index) PredictorPosition(int predictor) => _positions[predictor];

    /// <summary>The inputs of <paramref name="data"/> through <paramref name="transform"/>, fitted on it.</summary>
    public static BassInputs Prepare(BassData data, BassInputTransform transform)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(transform);
        if (transform.Predictors.Count != data.Predictors.Count)
            throw new ArgumentException("The transform was not fitted on these data.", nameof(transform));
        var continuous = new List<double[]>();
        var knots = new List<double[]>();
        var categorical = new List<int[]>();
        var levels = new List<int>();
        var positions = new (bool, int)[data.Predictors.Count];
        for (int j = 0; j < positions.Length; j++)
        {
            BassPredictor p = data.Predictors[j];
            if (transform.Predictors[j].Name != p.Name)
                throw new ArgumentException("The transform was not fitted on these data.", nameof(transform));
            switch ((p, transform.Predictors[j]))
            {
                case (ContinuousPredictor c, RangeScale scale):
                    positions[j] = (true, continuous.Count);
                    double[] scaled = Scaled(c, scale);
                    continuous.Add(scaled);
                    knots.Add(Distinct(scaled));
                    break;
                case (CategoricalPredictor f, TrainingLevels training):
                    positions[j] = (false, categorical.Count);
                    ReadOnlySpan<int> codes = f.Codes.Span;
                    var indices = new int[codes.Length];
                    for (int i = 0; i < indices.Length; i++) indices[i] = training.IndexOfCode(codes[i]);
                    categorical.Add(indices);
                    levels.Add(training.Levels.Count);
                    break;
                default:
                    throw new ArgumentException($"The transform's scale for '{p.Name}' does not match its kind.", nameof(transform));
            }
        }

        OutputInputs? outputs = null;
        double[][] functionalKnots = [];
        if (data.Response is FunctionalResponse functional)
        {
            if (transform.OutputVariables.Count != functional.OutputVariables.Count)
                throw new ArgumentException("The transform was not fitted on these data.", nameof(transform));
            var columns = new double[functional.OutputVariables.Count][];
            functionalKnots = new double[columns.Length][];
            for (int d = 0; d < columns.Length; d++)
            {
                columns[d] = Scaled(functional.OutputVariables[d], transform.OutputVariables[d]);
                functionalKnots[d] = Distinct(columns[d]);
            }
            outputs = new OutputInputs(columns);
        }
        return new BassInputs(new CaseInputs([.. continuous], [.. categorical]), outputs, [.. knots], functionalKnots,
                              [.. levels], positions);
    }

    private static double[] Scaled(ContinuousPredictor column, RangeScale scale)
    {
        ReadOnlySpan<double> raw = column.Values.Span;
        var scaled = new double[raw.Length];
        for (int i = 0; i < scaled.Length; i++) scaled[i] = scale.Scale(raw[i]);
        return scaled;
    }

    private static double[] Distinct(double[] values)
    {
        var sorted = (double[])values.Clone();
        Array.Sort(sorted);
        int count = 0;
        for (int i = 0; i < sorted.Length; i++)
            if (i == 0 || sorted[i] != sorted[count - 1]) sorted[count++] = sorted[i];
        Array.Resize(ref sorted, count);
        return sorted;
    }
}
