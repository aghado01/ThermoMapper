using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json.Serialization;

namespace Maths.Regression.Spline.Bass;

/// <summary>The scale and named grouping shared by independent Gaussian effect blocks.</summary>
public interface IBassRandomEffect
{
    string Grouping { get; }
    double Scale { get; }
    BassGroupScalePrior? ScalePrior { get; }
    string? ParentGrouping => null;
}

/// <summary>The training input units in which a random slope's scale is declared.</summary>
public enum RandomSlopeScale
{
    /// <summary>The fitted input transform's unit interval.</summary>
    Unit,
    /// <summary>Training mean zero and mean square one, fitted after the unit transform.</summary>
    Standardized,
}

/// <summary>
/// Independent slopes on a named continuous input, with explicitly chosen units. The contribution zᵢb_g is orthogonal to
/// the intercept and population slope over the training cases. Effects and scales are independent of every intercept block.
/// </summary>
public sealed record BassRandomSlope(string Grouping, string Input, [property: JsonRequired] RandomSlopeScale InputScale, double Scale,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] BassGroupScalePrior? ScalePrior = null) : IBassRandomEffect;

internal sealed class RandomSlopeInput
{
    private RandomSlopeInput(string name, int index, double mean, double deviation, double[] values)
    { Name = name; Index = index; Mean = mean; Deviation = deviation; Values = values; }

    public string Name { get; }
    public int Index { get; }
    public double Mean { get; }
    public double Deviation { get; }
    public double[] Values { get; }

    public static RandomSlopeInput Of(BassRandomSlope slope, IReadOnlyList<ContinuousPredictor>? predictors)
    {
        int index = -1;
        if (predictors is not null)
            for (int j = 0; j < predictors.Count; j++) if (predictors[j].Name == slope.Input) index = j;
        if (index < 0) throw new ArgumentException($"Random slope '{slope.Grouping}' must name a continuous predictor; '{slope.Input}' is absent.", nameof(slope));
        ContinuousPredictor predictor = predictors![index];
        RangeScale transform = RangeScale.Of(predictor);
        double[] unit = [.. predictor.Values.ToArray().Select(transform.Scale)];
        double mean = 0, deviation = 1;
        if (slope.InputScale == RandomSlopeScale.Standardized)
        {
            mean = unit.Average();
            deviation = Math.Sqrt(unit.Select(z => (z-mean)*(z-mean)).Average());
            if (!(deviation > 0) || !double.IsFinite(deviation))
                throw new ArgumentException($"Random slope input '{slope.Input}' has zero training deviation; Standardized scaling is undefined.", nameof(slope));
        }
        return new(slope.Input, index, mean, deviation, [.. unit.Select(z => (z-mean)/deviation)]);
    }

    public double At(CaseInputs inputs, int row) => (inputs.Continuous(Index)[row]-Mean)/Deviation;
}
