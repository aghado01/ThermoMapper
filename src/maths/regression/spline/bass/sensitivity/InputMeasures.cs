using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// One input's measure for sensitivity analysis (D08), declared in the input's own units: a continuous predictor's or
/// output variable's <see cref="UniformMeasure"/>, <see cref="MixtureMeasure"/> or <see cref="EmpiricalMeasure"/>, or a
/// categorical predictor's <see cref="CategoricalMeasure"/>. The inputs' measures form a product: their dependence is not
/// modelled. <see cref="ProductMeasure.Resolve"/> carries each into the fit's transformed units.
/// </summary>
public abstract record InputMeasure
{
    private protected InputMeasure() { }
}

/// <summary>U(<paramref name="Lower"/>, <paramref name="Upper"/>), each end the training range's when omitted (R's default).</summary>
public sealed record UniformMeasure(double? Lower = null, double? Upper = null) : InputMeasure;

/// <summary>
/// A finite mixture of truncated normal and Student t components on [<see cref="Lower"/>, <see cref="Upper"/>], each end the
/// training range's when omitted. <see cref="Weighting"/> says what the weights are: by default the weights before
/// truncation, the whole mixture renormalized on the interval (R's convention), or the components' shares on it.
/// </summary>
public sealed record MixtureMeasure(IReadOnlyList<MixtureComponent> Components) : InputMeasure
{
    /// <summary>The truncation's lower end; the training minimum when null.</summary>
    public double? Lower { get; init; }

    /// <summary>The truncation's upper end; the training maximum when null.</summary>
    public double? Upper { get; init; }

    /// <summary>What the components' weights are.</summary>
    public MixtureWeighting Weighting { get; init; } = MixtureWeighting.BeforeTruncation;
}

/// <summary>What a truncated mixture's weights are.</summary>
public enum MixtureWeighting
{
    /// <summary>The weights before truncation: the whole mixture is renormalized on the interval, so a component's share
    /// follows its mass there (R's convention).</summary>
    BeforeTruncation,

    /// <summary>The components' shares on the interval: each component is renormalized on it.</summary>
    ConditionalOnTruncation,
}

/// <summary>One mixture component, with its weight.</summary>
public abstract record MixtureComponent(double Weight);

/// <summary>A normal component N(<paramref name="Mean"/>, <paramref name="Sd"/>²).</summary>
public sealed record NormalComponent(double Weight, double Mean, double Sd) : MixtureComponent(Weight);

/// <summary>
/// A Student t component, location + scale·T_ν: <paramref name="Scale"/> is the t's scale, not its standard deviation
/// (R's <c>sd</c> is this scale); <see cref="FromSd"/> converts a marginal standard deviation for ν &gt; 2.
/// </summary>
public sealed record StudentComponent(double Weight, double DegreesOfFreedom, double Location, double Scale) : MixtureComponent(Weight)
{
    /// <summary>The component whose marginal standard deviation is <paramref name="sd"/>: scale = sd·√((ν − 2)/ν), for ν &gt; 2.</summary>
    public static StudentComponent FromSd(double weight, double degreesOfFreedom, double location, double sd)
    {
        if (!(degreesOfFreedom > 2.0)) throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "A Student t has a standard deviation only for ν > 2.");
        if (!(sd > 0.0) || double.IsPositiveInfinity(sd)) throw new ArgumentOutOfRangeException(nameof(sd), "The standard deviation must be positive and finite.");
        return new StudentComponent(weight, degreesOfFreedom, location, sd * Math.Sqrt((degreesOfFreedom - 2.0) / degreesOfFreedom));
    }
}

/// <summary>
/// The empirical marginal: the input's training values, equally weighted, repeats kept. For an output variable in the
/// integrated functional form, the training output grid's values.
/// </summary>
public sealed record EmpiricalMeasure : InputMeasure;

/// <summary>
/// A categorical predictor's level probabilities by level name, over its training levels: a level left out has probability
/// 0, and the values are normalized. Uniform over the training levels when null.
/// </summary>
public sealed record CategoricalMeasure(IReadOnlyDictionary<string, double>? Probabilities = null) : InputMeasure;

/// <summary>A continuous input's measure in the fit's transformed units: its interval moments and its support.</summary>
/// <param name="Measure">The measure's interval moments.</param>
/// <param name="Lower">The support's lower end.</param>
/// <param name="Upper">The support's upper end.</param>
public sealed record ContinuousMeasure(IIntervalMoments Measure, double Lower, double Upper)
{
    /// <summary>Whether the support reaches past the training range [0, 1], where the hinges extrapolate.</summary>
    public bool Extrapolates => Lower < 0.0 || Upper > 1.0;
}

/// <summary>
/// The product measure over a fit's inputs, in its transformed units: one measure per predictor, in the data's order, then,
/// when <see cref="IncludesOutputs"/>, one per output variable, the integrated functional form's extra inputs. Its
/// variables are numbered in that order, which is how sensitivity results name inputs.
/// </summary>
public sealed class ProductMeasure
{
    private readonly ContinuousMeasure?[] _continuous;
    private readonly double[]?[] _probabilities;
    private readonly int[] _ofContinuous, _ofCategorical;
    private readonly int _outputOffset;

    private ProductMeasure(string[] names, InputMeasure[] declared, ContinuousMeasure?[] continuous, double[]?[] probabilities,
                           int[] ofContinuous, int[] ofCategorical, int outputOffset, bool includesOutputs)
    {
        Names = Array.AsReadOnly(names);
        Declared = Array.AsReadOnly(declared);
        _continuous = continuous;
        _probabilities = probabilities;
        _ofContinuous = ofContinuous;
        _ofCategorical = ofCategorical;
        _outputOffset = outputOffset;
        IncludesOutputs = includesOutputs;
    }

    /// <summary>The variables' names: the predictors', then the output variables' when included.</summary>
    public IReadOnlyList<string> Names { get; }

    /// <summary>Each variable's measure as declared, defaults filled in.</summary>
    public IReadOnlyList<InputMeasure> Declared { get; }

    /// <summary>The variables.</summary>
    public int Count => Names.Count;

    /// <summary>Whether the output variables are inputs (the integrated functional form).</summary>
    public bool IncludesOutputs { get; }

    /// <summary>Variable <paramref name="variable"/>'s measure if it is continuous; null for a categorical predictor.</summary>
    public ContinuousMeasure? Continuous(int variable) => _continuous[variable];

    /// <summary>Variable <paramref name="variable"/>'s level probabilities over its training levels if it is categorical.</summary>
    public IReadOnlyList<double>? Probabilities(int variable) => _probabilities[variable];

    /// <summary>The variable of continuous predictor <paramref name="index"/> (its position among the continuous ones).</summary>
    public int OfContinuous(int index) => _ofContinuous[index];

    /// <summary>The variable of categorical predictor <paramref name="index"/>.</summary>
    public int OfCategorical(int index) => _ofCategorical[index];

    /// <summary>The variable of output variable <paramref name="index"/>; −1 when the outputs are not inputs.</summary>
    public int OfOutput(int index) => IncludesOutputs ? _outputOffset + index : -1;

    /// <summary>
    /// The product of <paramref name="measures"/> over the predictors and, when <paramref name="includeOutputs"/>,
    /// <paramref name="outputMeasures"/> over the output variables, each keyed by name and declared in the variable's own
    /// units, carried into the units of <paramref name="transform"/>, which <paramref name="inputs"/> were prepared with. A
    /// variable left out takes the uniform measure on its training range, or the uniform level probabilities. A name that is
    /// not a variable's, or a measure of the wrong kind, is refused.
    /// </summary>
    public static ProductMeasure Resolve(BassInputTransform transform, BassInputs inputs,
                                         IReadOnlyDictionary<string, InputMeasure>? measures = null,
                                         IReadOnlyDictionary<string, InputMeasure>? outputMeasures = null,
                                         bool includeOutputs = false)
    {
        ArgumentNullException.ThrowIfNull(transform);
        ArgumentNullException.ThrowIfNull(inputs);
        measures ??= new Dictionary<string, InputMeasure>();
        outputMeasures ??= new Dictionary<string, InputMeasure>();
        if (includeOutputs && inputs.Outputs is null) throw new ArgumentException("A scalar response has no output variables to integrate.", nameof(includeOutputs));
        Unknown(measures.Keys, transform.Predictors.Select(p => p.Name), "predictor", nameof(measures));
        Unknown(outputMeasures.Keys, transform.OutputVariables.Select(p => p.Name), "output variable", nameof(outputMeasures));
        if (!includeOutputs && outputMeasures.Count > 0)
            throw new ArgumentException("Output-variable measures apply only when the output variables are integrated.", nameof(outputMeasures));

        int p = transform.Predictors.Count, outputs = includeOutputs ? transform.OutputVariables.Count : 0;
        var names = new string[p + outputs];
        var declared = new InputMeasure[names.Length];
        var continuous = new ContinuousMeasure?[names.Length];
        var probabilities = new double[]?[names.Length];
        var ofContinuous = new int[inputs.Cases.ContinuousCount];
        var ofCategorical = new int[inputs.Cases.CategoricalCount];
        for (int j = 0; j < p; j++)
        {
            PredictorScale scale = transform.Predictors[j];
            names[j] = scale.Name;
            measures.TryGetValue(scale.Name, out InputMeasure? given);
            (bool isContinuous, int index) = inputs.PredictorPosition(j);
            if (isContinuous)
            {
                ofContinuous[index] = j;
                declared[j] = given ?? new UniformMeasure();
                continuous[j] = ResolveContinuous(scale.Name, declared[j], (RangeScale)scale, inputs.Cases.Continuous(index));
            }
            else
            {
                ofCategorical[index] = j;
                declared[j] = given ?? new CategoricalMeasure();
                probabilities[j] = ResolveCategorical(scale.Name, declared[j], (TrainingLevels)scale);
            }
        }
        for (int d = 0; d < outputs; d++)
        {
            RangeScale scale = transform.OutputVariables[d];
            names[p + d] = scale.Name;
            outputMeasures.TryGetValue(scale.Name, out InputMeasure? given);
            declared[p + d] = given ?? new UniformMeasure();
            continuous[p + d] = ResolveContinuous(scale.Name, declared[p + d], scale, inputs.Outputs!.Variable(d));
        }
        return new ProductMeasure(names, declared, continuous, probabilities, ofContinuous, ofCategorical, p, includeOutputs);
    }

    private static void Unknown(IEnumerable<string> given, IEnumerable<string> known, string kind, string parameter)
    {
        string[] unknown = [.. given.Except(known, StringComparer.Ordinal)];
        if (unknown.Length > 0) throw new ArgumentException($"No {kind} is named {string.Join(", ", unknown.Select(n => $"'{n}'"))}.", parameter);
    }

    // The declared measure in transformed units: (v − min)/(max − min), scales divided by the range.
    private static ContinuousMeasure ResolveContinuous(string name, InputMeasure declared, RangeScale scale, ReadOnlySpan<double> scaledTraining)
    {
        double range = scale.Max - scale.Min;
        double At(double raw) => (raw - scale.Min) / range;
        switch (declared)
        {
            case UniformMeasure u:
            {
                double lo = u.Lower is { } l ? At(l) : 0.0, hi = u.Upper is { } h ? At(h) : 1.0;
                if (!double.IsFinite(lo) || !double.IsFinite(hi) || !(hi > lo))
                    throw new ArgumentException($"The uniform measure of '{name}' needs finite ends, lower below upper.");
                return new ContinuousMeasure(Uniform.On(lo, hi), lo, hi);
            }
            case MixtureMeasure m:
            {
                if (m.Components is null || m.Components.Count == 0) throw new ArgumentException($"The mixture measure of '{name}' needs a component.");
                double lo = m.Lower is { } l ? At(l) : 0.0, hi = m.Upper is { } h ? At(h) : 1.0;
                if (!double.IsFinite(lo) || !double.IsFinite(hi) || !(hi > lo))
                    throw new ArgumentException($"The mixture measure of '{name}' needs a finite truncation, lower below upper.");
                try
                {
                    var weights = m.Components.Select(c => c.Weight).ToArray();
                    if (m.Weighting == MixtureWeighting.BeforeTruncation)
                    {
                        IContinuousDistribution[] parts = [.. m.Components.Select(c => Component(c, scale.Min, range))];
                        return new ContinuousMeasure(Truncated.Of(Mixture.Of(parts, weights), lo, hi), lo, hi);
                    }
                    IContinuousDistribution[] truncated = [.. m.Components.Select(c => TruncatedComponent(c, scale.Min, range, lo, hi))];
                    return new ContinuousMeasure(Mixture.Of(truncated, weights), lo, hi);
                }
                catch (ArgumentException e)
                {
                    throw new ArgumentException($"The mixture measure of '{name}': {e.Message}", e);
                }
            }
            case EmpiricalMeasure:
            {
                Empirical empirical = Empirical.Of(scaledTraining);
                return new ContinuousMeasure(empirical, empirical.Lower, empirical.Upper);
            }
            default:
                throw new ArgumentException($"'{name}' is continuous; its measure is uniform, a mixture or empirical, not {declared?.GetType().Name ?? "null"}.");
        }
    }

    private static IContinuousDistribution Component(MixtureComponent c, double min, double range) => c switch
    {
        NormalComponent n => Normal.MeanSd((n.Mean - min) / range, n.Sd / range),
        StudentComponent s => StudentT.LocationScale(s.DegreesOfFreedom, (s.Location - min) / range, s.Scale / range),
        _ => throw new ArgumentException($"A mixture component is normal or Student t, not {c?.GetType().Name ?? "null"}."),
    };

    private static IContinuousDistribution TruncatedComponent(MixtureComponent c, double min, double range, double lo, double hi) => c switch
    {
        NormalComponent n => Truncated.Of(Normal.MeanSd((n.Mean - min) / range, n.Sd / range), lo, hi),
        StudentComponent s => Truncated.Of(StudentT.LocationScale(s.DegreesOfFreedom, (s.Location - min) / range, s.Scale / range), lo, hi),
        _ => throw new ArgumentException($"A mixture component is normal or Student t, not {c?.GetType().Name ?? "null"}."),
    };

    private static double[] ResolveCategorical(string name, InputMeasure declared, TrainingLevels levels)
    {
        if (declared is not CategoricalMeasure c)
            throw new ArgumentException($"'{name}' is categorical; its measure is level probabilities, not {declared?.GetType().Name ?? "null"}.");
        int count = levels.Levels.Count;
        var p = new double[count];
        if (c.Probabilities is null)
        {
            Array.Fill(p, 1.0 / count);
            return p;
        }
        double total = 0.0;
        foreach ((string level, double probability) in c.Probabilities)
        {
            int l = levels.IndexOf(level);
            if (l < 0) throw new ArgumentException($"'{level}' is not a training level of '{name}'.");
            if (!(probability >= 0.0) || double.IsPositiveInfinity(probability))
                throw new ArgumentException($"The probability of level '{level}' of '{name}' must be nonnegative and finite.");
            p[l] = probability;
            total += probability;
        }
        if (!(total > 0.0)) throw new ArgumentException($"The level probabilities of '{name}' need a positive sum.");
        for (int l = 0; l < count; l++) p[l] /= total;
        return p;
    }
}
