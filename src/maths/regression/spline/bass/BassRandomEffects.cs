using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json.Serialization;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>One retained grouping's scale and effects, in the observed levels' declaration order.</summary>
public sealed record BassGroupDraw(string Grouping, double Scale, double[] Effects)
{
    /// <summary>The slope's continuous input; absent for an intercept block.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public string? Input { get; init; }

    /// <summary>Several outputs' effects (E34's slice 4), one array per output over the observed levels, u_g ∈ ℝᴷ with one s for every output;
    /// <see cref="Effects"/> is then empty. Absent for one response.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public double[][]? OutputEffects { get; init; }
}

/// <summary>Named Gaussian intercept and slope blocks, jointly integrated through their row covariance.</summary>
public sealed class RandomEffectsBlock
{
    private readonly IReadOnlyList<IBassRandomEffect> _specifications;
    private readonly IReadOnlyList<CategoricalPredictor> _groupings;
    private readonly RandomEffectGrouping[] _layouts;
    private readonly RandomEffectGeometry? _geometry;
    private readonly RandomEffectWhitening? _general;
    private readonly RandomInterceptBlock? _single;
    private readonly double[] _scales;

    private RandomEffectsBlock(IBassRandomEffect[] specifications, CategoricalPredictor[] groupings, RandomEffectGrouping[] layouts,
                               RandomEffectGeometry? geometry, RandomInterceptBlock? single = null)
    {
        _specifications = Array.AsReadOnly(specifications);
        _groupings = Array.AsReadOnly(groupings);
        _layouts = layouts;
        _geometry = geometry;
        _scales = [.. specifications.Select(s => s.Scale)];
        if (geometry is null) _single = single ?? RandomInterceptBlock.Of(groupings[0], (BassRandomIntercept)specifications[0]);
        else _general = geometry.At(_scales);
        LogScalePrior = _single?.LogScalePrior ?? specifications.Sum(s => s.ScalePrior is { } prior ? prior.LogDensity(s.Scale) + Math.Log(s.Scale) : 0);
    }

    /// <summary>Each block's grouping layout, in the declaration order.</summary>
    internal IReadOnlyList<RandomEffectGrouping> Layouts => _layouts;

    public static RandomEffectsBlock Of(IReadOnlyList<CategoricalPredictor> groupings, IReadOnlyList<IBassRandomEffect> specifications,
                                        IReadOnlyList<ContinuousPredictor>? predictors = null)
    {
        ArgumentNullException.ThrowIfNull(groupings);
        ArgumentNullException.ThrowIfNull(specifications);
        if (specifications.Any(s => s is not (BassRandomIntercept or BassRandomSlope)))
            throw new ArgumentException("Each effect must be a declared intercept or slope.", nameof(specifications));
        BassConfig.CheckRandomIntercepts(specifications.OfType<BassRandomIntercept>().ToArray());
        BassConfig.CheckRandomSlopes(specifications.OfType<BassRandomSlope>().ToArray());
        if (specifications.Count == 0) throw new ArgumentException("At least one random-effect block is required.", nameof(specifications));
        var available = new Dictionary<string, CategoricalPredictor>(StringComparer.Ordinal);
        foreach (CategoricalPredictor grouping in groupings)
        {
            ArgumentNullException.ThrowIfNull(grouping);
            if (!available.TryAdd(grouping.Name, grouping)) throw new ArgumentException($"Grouping '{grouping.Name}' is declared twice.", nameof(groupings));
        }
        var columns = new CategoricalPredictor[specifications.Count];
        var layouts = new RandomEffectGrouping[specifications.Count];
        for (int k = 0; k < specifications.Count; k++)
        {
            IBassRandomEffect spec = specifications[k];
            if (!available.TryGetValue(spec.Grouping, out CategoricalPredictor? grouping))
                throw new ArgumentException($"Random-effect block '{spec.Grouping}' has no data grouping.", nameof(specifications));
            columns[k] = grouping;
            if (spec is BassRandomSlope slope)
                layouts[k] = RandomEffectGrouping.Slope(grouping, RandomSlopeInput.Of(slope, predictors));
            else if (((BassRandomIntercept)spec).Convention == RandomInterceptConvention.ObservationCentred)
                layouts[k] = RandomEffectGrouping.ObservationCentred(grouping);
            else
            {
                if (!available.TryGetValue(spec.ParentGrouping!, out CategoricalPredictor? parent))
                    throw new ArgumentException($"Parent grouping '{spec.ParentGrouping}' of '{spec.Grouping}' is missing.", nameof(specifications));
                layouts[k] = RandomEffectGrouping.WithinParent(grouping, parent);
            }
        }
        RandomEffectGeometry? geometry = specifications.Count == 1 && specifications[0] is BassRandomIntercept { Convention: RandomInterceptConvention.ObservationCentred }
            ? null : new RandomEffectGeometry(layouts);
        return new RandomEffectsBlock(specifications.ToArray(), columns, layouts, geometry);
    }

    public IReadOnlyList<IBassRandomEffect> Specifications => _specifications;
    public IReadOnlyList<CategoricalPredictor> Groupings => _groupings;
    public ReadOnlyMemory<double> Scales => _scales;
    public int Cases => _layouts[0].Cases;
    public int Rank => _geometry?.Rank ?? _single!.Levels - 1;
    public int Count => _layouts.Length;
    public double LogDeterminant => _single?.LogDeterminant ?? _general!.LogDeterminant;
    public double LogScalePrior { get; }
    internal IReadOnlyList<RandomEffectAlias> Aliases => _geometry?.Aliases ?? Array.Empty<RandomEffectAlias>();
    internal RandomEffectGrouping Layout(int group) => _layouts[group];

    public RandomEffectsBlock At(ReadOnlySpan<double> scales)
    {
        if (scales.Length != Count) throw new ArgumentException("Give one scale per random-effect block.", nameof(scales));
        if (scales.SequenceEqual(_scales)) return this;
        var specifications = new IBassRandomEffect[Count];
        for (int k = 0; k < Count; k++)
        {
            if (!(scales[k] > 0) || !double.IsFinite(scales[k])) throw new ArgumentOutOfRangeException(nameof(scales));
            specifications[k] = _specifications[k] switch
            {
                BassRandomIntercept intercept => intercept with { Scale = scales[k] },
                BassRandomSlope slope => slope with { Scale = scales[k] },
                _ => throw new InvalidOperationException("Unknown effect declaration."),
            };
        }
        return new RandomEffectsBlock(specifications, _groupings.ToArray(), _layouts, _geometry, _single?.At(scales[0]));
    }

    public RandomEffectsBlock At(int group, double scale)
    {
        if ((uint)group >= (uint)Count) throw new ArgumentOutOfRangeException(nameof(group));
        double[] scales = (double[])_scales.Clone();
        scales[group] = scale;
        return At(scales);
    }

    internal RandomEffectsBlock At(IReadOnlyList<BassGroupDraw> groups)
        => At(Ordered(groups).Select(g => g.Scale).ToArray());

    // Resolve retained groups by name; data, configuration and draw collections need not have the same order.
    internal BassGroupDraw[] Ordered(IReadOnlyList<BassGroupDraw> groups)
    {
        ArgumentNullException.ThrowIfNull(groups);
        if (groups.Count != Count) throw new ArgumentException("The draw must carry every configured grouping.", nameof(groups));
        var named = new Dictionary<(string Grouping, string? Input), BassGroupDraw>();
        foreach (BassGroupDraw group in groups)
        {
            ArgumentNullException.ThrowIfNull(group);
            if (!named.TryAdd((group.Grouping, group.Input), group)) throw new ArgumentException($"The draw repeats grouping '{group.Grouping}' on input '{group.Input}'.", nameof(groups));
        }
        var ordered = new BassGroupDraw[Count];
        for (int k = 0; k < Count; k++)
        {
            string name = _specifications[k].Grouping;
            if (!named.TryGetValue((name, _layouts[k].Input), out BassGroupDraw? group)) throw new ArgumentException($"The draw has no effects for grouping '{name}' on input '{_layouts[k].Input}'.", nameof(groups));
            if (group.OutputEffects is { } several ? Array.Exists(several, e => e is null || e.Length != _layouts[k].Levels)
                                                   : group.Effects is null || group.Effects.Length != _layouts[k].Levels)
                throw new ArgumentException($"The draw's effects do not match the observed levels of '{name}'.", nameof(groups));
            if (!(group.Scale > 0) || !double.IsFinite(group.Scale)) throw new ArgumentException($"The draw's scale of '{name}' must be finite and positive.", nameof(groups));
            ordered[k] = group;
        }
        return ordered;
    }

    public void Whiten(Span<double> x)
    {
        if (_single is not null) _single.Whiten(x);
        else _general!.Whiten(x);
    }

    /// <summary>The jointly whitened reading, under scalar Gaussian noise and unit row weights.</summary>
    public WorkingResponse Whitened(WorkingResponse working)
    {
        ArgumentNullException.ThrowIfNull(working);
        if (!working.IsUnweighted || working.Latents is not null || !working.Offset.IsEmpty || working.Response is not (ScalarResponse or MultiOutputResponse))
            throw new ArgumentException("Random effects whiten a scalar response or several outputs under unit weights, with no offset or latents.", nameof(working));
        return WhitenResponse(working);
    }

    // Replace the old whitening from raw observations, preserving the design's response and prior-reference identities.
    internal WorkingResponse WhitenResponse(WorkingResponse working)
    {
        ReadOnlySpan<double> y = working.Response.Values;
        double[] white = y.ToArray();
        if (working.Response is MultiOutputResponse outputs)
        {
            // Several outputs (E34's slice 4): V^(−½) on each output's column, the rows' covariance V ⊗ Σ under one s for every output.
            int k = outputs.Outputs, n = outputs.Cases;
            var column = new double[n];
            for (int l = 0; l < k; l++)
            {
                for (int i = 0; i < n; i++) column[i] = white[i * k + l];
                Whiten(column);
                for (int i = 0; i < n; i++) white[i * k + l] = column[i];
            }
        }
        else Whiten(white);
        var offset = new double[y.Length];
        for (int i = 0; i < y.Length; i++) offset[i] = y[i] - white[i];
        return new WorkingResponse(working.Response, offset, output: working.Output, meaning: working.Meaning, priorGram: working.PriorGram,
                                   reference: working.Reference, logJacobian: -LogDeterminant) { Whitening = this };
    }

    /// <summary>
    /// <paramref name="reading"/> whitened by this block at its scales. At unit weights V^(−½) acts on the raw response and the offset
    /// carries it (<see cref="WhitenResponse"/>). Under weights, a noise's, the reading keeps its own offset and weights and carries the
    /// weighted whitening beside them (<see cref="WeightedWhitening"/>; E21's slice 4), which the design applies to the response and each
    /// column.
    /// </summary>
    internal WorkingResponse Rewhiten(WorkingResponse reading)
    {
        ArgumentNullException.ThrowIfNull(reading);
        if (reading.IsUnweighted) return WhitenResponse(reading);
        if (reading.Response is not ScalarResponse) throw new ArgumentException("Random effects under weights take a scalar response.", nameof(reading));
        return reading.WithWhitening(this, new WeightedWhitening(this, reading.Weights));
    }

    internal double ScaleInformation(int group) => _single?.LogScaleInformation ?? _general!.ScaleInformation(group);

    internal double[][] MeanEffects(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients)
        => _single is { } one ? [one.EffectMean(one.ResidualSums(model, terms, coefficients))] : _general!.MeanEffects(Residual(model, terms, coefficients));

    internal double[][] DrawEffects(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients,
                                     double sigma2, Xoshiro256PlusPlus rng)
        => _single is { } one ? [one.DrawEffects(one.ResidualSums(model, terms, coefficients), sigma2, rng)]
            : _general!.DrawEffects(Residual(model, terms, coefficients), sigma2, rng);

    internal BassGroupDraw[] DrawGroups(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients,
                                       double sigma2, Xoshiro256PlusPlus rng, WorkingResponse? reading = null)
    {
        // Under a noise's weights, from u's conditional given the coefficients, σ² and the latents, on the reading's own y − σβv.
        double[][] effects = reading?.Weighted is { } weighted
            ? weighted.DrawEffects(Residual(model, terms, coefficients, reading.Values), sigma2, rng).Effects
            : DrawEffects(model, terms, coefficients, sigma2, rng);
        return Named(effects);
    }

    /// <summary>The effects' conditional mean given the coefficients: on <paramref name="reading"/>'s weights and y − σβv under a noise, where
    /// it carries u as a latent or integrates it; at unit weights the block's own.</summary>
    internal double[][] MeanEffects(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients, WorkingResponse? reading)
    {
        if (reading?.Weighted is { } weighted) return weighted.MeanEffects(Residual(model, terms, coefficients, reading.Values));
        if (reading is { IsUnweighted: false, Effects: { } latent })
        {
            double[] values = reading.Values.ToArray();   // y − σβv − Z̃u, with Z̃u put back
            latent.AddTo(values);
            return new WeightedWhitening(this, reading.Weights).MeanEffects(Residual(model, terms, coefficients, values));
        }
        return MeanEffects(model, terms, coefficients);
    }

    /// <summary>
    /// Several outputs' group effects at composition (E34's slice 4, decision 7): u (levels × K per grouping) from its conditional given the
    /// coefficients and Σ, MN(M, K_u, Σ), with u_g ∈ ℝᴷ ~ N(0, s²Σ) a priori. M's column l is the scalar conditional mean on output l's residual
    /// over the unwhitened columns, and the spread is the scalar block's standardized noise drawn once per output, in output order, mixed across
    /// the outputs by Σ's factor <paramref name="sigmaRoot"/>: at K = 1, mean + √Σ₀₀·noise on the scalar draw's stream.
    /// <paramref name="coefficients"/> is (M + 1) × K, the intercepts first and the terms in <paramref name="terms"/>' order.
    /// </summary>
    internal BassGroupDraw[] DrawOutputGroups(BassModel model, IReadOnlyList<BassTerm> terms, double[,] coefficients, double[,] sigmaRoot, Xoshiro256PlusPlus rng)
    {
        double[][] residuals = OutputResiduals(model, terms, coefficients);
        int k = residuals.Length;
        var means = new double[k][][];
        var noises = new double[k][][];
        for (int l = 0; l < k; l++)
        {
            means[l] = OutputMeanOf(residuals[l]);
            noises[l] = _single is { } one ? [one.DrawEffects(new double[one.Levels], 1.0, rng)] : _general!.DrawEffects(new double[residuals[l].Length], 1.0, rng);
        }
        var groups = new BassGroupDraw[Count];
        for (int g = 0; g < Count; g++)
        {
            int levels = means[0][g].Length;
            var effects = new double[k][];
            for (int l = 0; l < k; l++)
            {
                var u = new double[levels];
                for (int h = 0; h < levels; h++)
                {
                    double mixed = 0.0;
                    for (int q = 0; q <= l; q++) mixed += sigmaRoot[l, q] * noises[q][g][h];
                    u[h] = means[l][g][h] + mixed;
                }
                effects[l] = u;
            }
            groups[g] = new BassGroupDraw(_specifications[g].Grouping, _scales[g], []) { Input = _layouts[g].Input, OutputEffects = effects };
        }
        return groups;
    }

    /// <summary>Several outputs' effects' conditional mean given <paramref name="coefficients"/> ((M + 1) × K): per output the scalar block's, on
    /// its residual over the unwhitened columns; [output][grouping][level].</summary>
    internal double[][][] MeanOutputEffects(BassModel model, IReadOnlyList<BassTerm> terms, double[,] coefficients, WorkingResponse? reading = null)
    {
        // Under a noise U is the reading's latent (E34's slice 4): each output's mean at the reading's weights, on y − vγ (the reading's values
        // with Z̃U put back).
        if (reading is { IsUnweighted: false, Effects: { OutputEffects: not null } latent })
        {
            double[] values = reading.Values.ToArray();
            latent.AddTo(values);
            var whitening = new WeightedWhitening(this, reading.Weights);
            return [.. OutputResiduals(model, terms, coefficients, values).Select(r => whitening.MeanEffects(r))];
        }
        double[][] residuals = OutputResiduals(model, terms, coefficients);
        return [.. residuals.Select(OutputMeanOf)];
    }

    private double[][] OutputMeanOf(double[] residual)
        => _single is { } one ? [one.EffectMean(one.GroupSums(residual))] : _general!.MeanEffects(residual);

    // Each output's y_l − α_l − Σⱼ β_jl·Bⱼ on the unwhitened case columns; `of`, when given, in place of y (case by case, output by output).
    private static double[][] OutputResiduals(BassModel model, IReadOnlyList<BassTerm> terms, double[,] coefficients, double[]? of = null)
    {
        var response = (MultiOutputResponse)model.Data.Response;
        int k = response.Outputs, n = response.Cases;
        if (coefficients.GetLength(0) != terms.Count + 1 || coefficients.GetLength(1) != k)
            throw new ArgumentException("Give the intercepts and one row per term, one column per output.", nameof(coefficients));
        ReadOnlySpan<double> values = of ?? response.Values;
        var residuals = new double[k][];
        for (int l = 0; l < k; l++)
        {
            var r = new double[n];
            for (int i = 0; i < n; i++) r[i] = values[i * k + l] - coefficients[0, l];
            residuals[l] = r;
        }
        var column = new double[n];
        for (int j = 0; j < terms.Count; j++)
        {
            model.Basis.CaseColumn(terms[j], model.Inputs.Cases, column);
            for (int l = 0; l < k; l++)
            {
                double beta = coefficients[j + 1, l];
                for (int i = 0; i < n; i++) residuals[l][i] -= beta * column[i];
            }
        }
        return residuals;
    }

    // Several outputs' effects ([grouping][output][level]) as named retained records, at this block's scales.
    internal BassGroupDraw[] NamedOutputs(double[][][] effects)
    {
        var groups = new BassGroupDraw[Count];
        for (int k = 0; k < Count; k++) groups[k] = new BassGroupDraw(_specifications[k].Grouping, _scales[k], []) { Input = _layouts[k].Input, OutputEffects = effects[k] };
        return groups;
    }

    // The effects as named retained records, at this block's scales.
    internal BassGroupDraw[] Named(double[][] effects)
    {
        var groups = new BassGroupDraw[Count];
        for (int k = 0; k < Count; k++) groups[k] = new BassGroupDraw(_specifications[k].Grouping, _scales[k], effects[k]) { Input = _layouts[k].Input };
        return groups;
    }

    /// <summary>Adds Σ_k z_k·u_k at each training case, the effects' contribution, to <paramref name="fit"/>.</summary>
    internal void AddEffects(double[][] effects, Span<double> fit)
    {
        for (int k = 0; k < Count; k++)
            for (int i = 0; i < fit.Length; i++) fit[i] += _layouts[k].Multiplier(i) * effects[k][_layouts[k].GroupOf(i)];
    }

    private static double[] Residual(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients)
        => Residual(model, terms, coefficients, model.Data.Response.Values);

    // values − α − Σⱼ βⱼBⱼ on the unwhitened case columns.
    private static double[] Residual(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients, ReadOnlySpan<double> values)
    {
        if (coefficients.Length != terms.Count + 1) throw new ArgumentException("Give an intercept and one coefficient per term.", nameof(coefficients));
        double[] residual = values.ToArray();
        for (int i = 0; i < residual.Length; i++) residual[i] -= coefficients[0];
        var column = new double[residual.Length];
        for (int j = 0; j < terms.Count; j++)
        {
            model.Basis.CaseColumn(terms[j], model.Inputs.Cases, column);
            for (int i = 0; i < residual.Length; i++) residual[i] -= coefficients[j + 1] * column[i];
        }
        return residual;
    }
}

/// <summary>One grouping's observed levels and orthonormal level contrasts, under observation or within-parent centring.</summary>
internal sealed class RandomEffectGrouping
{
    private readonly int[] _group, _observed;
    private readonly double[] _counts;
    private readonly double[,] _contrast;
    private readonly string[] _levels;
    private readonly string[]? _parentLevels;
    private readonly HashSet<string>? _observedParents;
    private readonly RandomSlopeInput? _slope;

    private RandomEffectGrouping(CategoricalPredictor grouping, CategoricalPredictor? parent, RandomSlopeInput? slope = null)
    {
        _slope = slope;
        Name = grouping.Name;
        ParentName = parent?.Name;
        _observed = new int[grouping.Levels.Count];
        Array.Fill(_observed, -1);
        var levels = new List<string>();
        for (int g = 0; g < _observed.Length; g++)
            if (grouping.Observed[g])
            {
                _observed[g] = levels.Count;
                levels.Add(grouping.Levels[g]);
            }
        _levels = [.. levels];
        _group = [.. grouping.Codes.ToArray().Select(g => _observed[g])];
        _counts = new double[levels.Count];
        foreach (int g in _group) _counts[g]++;
        var parents = new int[levels.Count];
        Array.Fill(parents, -1);
        if (parent is not null)
        {
            if (parent.Length != grouping.Length) throw new ArgumentException("The grouping and its parent have different case counts.", nameof(parent));
            if (parent.Name == grouping.Name) throw new ArgumentException("A grouping cannot be its own parent.", nameof(parent));
            ReadOnlySpan<int> codes = parent.Codes.Span;
            for (int i = 0; i < _group.Length; i++)
            {
                int g = _group[i];
                if (parents[g] >= 0 && parents[g] != codes[i])
                    throw new ArgumentException($"Level '{_levels[g]}' of '{Name}' belongs to more than one level of parent '{parent.Name}'.", nameof(parent));
                parents[g] = codes[i];
            }
        }
        else Array.Fill(parents, 0);

        if (parent is not null)
        {
            _parentLevels = [.. parents.Select(p => parent.Levels[p])];
            _observedParents = new HashSet<string>(_parentLevels, StringComparer.Ordinal);
        }

        if (slope is not null)
        {
            if (slope.Values.Length != Cases) throw new ArgumentException("Slope input and grouping have different case counts.", nameof(slope));
            var first = new double[Levels];
            var second = new double[Levels];
            for (int i = 0; i < Cases; i++)
            {
                first[_group[i]] += slope.Values[i];
                second[_group[i]] += slope.Values[i]*slope.Values[i];
            }
            _contrast = Complement([first, second]);
            return;
        }

        // The constraint rows have disjoint supports. Build the Euclidean complement of each count vector separately.
        var basis = new List<double[]>();
        foreach (int p in parents.Distinct())
        {
            int[] children = [.. Enumerable.Range(0, levels.Count).Where(g => parents[g] == p)];
            double norm = Math.Sqrt(children.Sum(g => _counts[g] * _counts[g]));
            int start = basis.Count;
            foreach (int child in children)
            {
                if (basis.Count - start == children.Length - 1) break;
                var column = new double[levels.Count];
                column[child] = 1;
                for (int pass = 0; pass < 2; pass++)
                {
                    double along = children.Sum(g => column[g] * _counts[g]) / (norm * norm);
                    foreach (int g in children) column[g] -= along * _counts[g];
                    for (int k = start; k < basis.Count; k++)
                    {
                        along = PairwiseSum.Dot(column, basis[k]);
                        for (int g = 0; g < column.Length; g++) column[g] -= along * basis[k][g];
                    }
                }
                double length = Math.Sqrt(PairwiseSum.Dot(column, column));
                if (length <= 1e-10) continue;
                for (int g = 0; g < column.Length; g++) column[g] /= length;
                basis.Add(column);
            }
            if (basis.Count - start != children.Length - 1) throw new InvalidOperationException("The grouping's contrast basis lost rank.");
        }
        _contrast = new double[levels.Count, basis.Count];
        for (int k = 0; k < basis.Count; k++)
            for (int g = 0; g < levels.Count; g++) _contrast[g, k] = basis[k][g];
    }

    public static RandomEffectGrouping ObservationCentred(CategoricalPredictor grouping)
    {
        ArgumentNullException.ThrowIfNull(grouping);
        return new RandomEffectGrouping(grouping, null);
    }

    public static RandomEffectGrouping WithinParent(CategoricalPredictor grouping, CategoricalPredictor parent)
    {
        ArgumentNullException.ThrowIfNull(grouping);
        ArgumentNullException.ThrowIfNull(parent);
        return new RandomEffectGrouping(grouping, parent);
    }

    public static RandomEffectGrouping Slope(CategoricalPredictor grouping, RandomSlopeInput input) => new(grouping, null, input);

    private static double[,] Complement(double[][] constraints)
    {
        int levels = constraints[0].Length;
        var matrix = new double[constraints.Length,levels];
        for (int k = 0; k < constraints.Length; k++)
            for (int g = 0; g < levels; g++) matrix[k,g] = constraints[k][g];
        ThinSvdResult factor = ThinSvd.Decompose(matrix);
        // Reveal rank jointly: independently normalizing a nearly zero row would turn rounding noise into a constraint.
        List<double[]> rows = [.. factor.RightVectors.Take(factor.NumericalRank)];
        var basis = new List<double[]>();
        for (int g = 0; g < levels && basis.Count < levels-rows.Count; g++)
        {
            var column = new double[levels];
            column[g] = 1;
            for (int pass = 0; pass < 2; pass++)
            {
                RemoveRows(column, rows);
                RemoveRows(column, basis);
            }
            double norm = Math.Sqrt(PairwiseSum.Dot(column, column));
            if (norm <= 1e-12) continue;
            for (int h = 0; h < levels; h++) column[h] /= norm;
            basis.Add(column);
        }
        if (basis.Count != levels-rows.Count) throw new InvalidOperationException("The slope's constraint complement lost rank.");
        var result = new double[levels,basis.Count];
        for (int g = 0; g < levels; g++)
            for (int j = 0; j < basis.Count; j++) result[g,j] = basis[j][g];
        return result;
    }

    private static void RemoveRows(double[] column, List<double[]> rows)
    {
        for (int pass = 0; pass < 2; pass++)
            foreach (double[] row in rows)
            {
                double dot = PairwiseSum.Dot(column, row);
                for (int g = 0; g < column.Length; g++) column[g] -= dot*row[g];
            }
    }

    public string Name { get; }
    public string? ParentName { get; }
    public string? Input => _slope?.Name;
    public int Cases => _group.Length;
    public int Levels => _counts.Length;
    public int Dimension => _contrast.GetLength(1);
    public ReadOnlySpan<double> Counts => _counts;
    public int GroupOf(int i) => _group[i];
    public int ObservedOf(int declared) => (uint)declared < (uint)_observed.Length ? _observed[declared] : -1;
    public string LevelName(int g) => _levels[g];
    public string? ParentLevelOf(int g) => _parentLevels?[g];
    public bool IsObservedParentLevel(string level) => _observedParents?.Contains(level) ?? false;
    public double Contrast(int g, int k) => _contrast[g, k];
    public double Multiplier(int row) => _slope?.Values[row] ?? 1;
    public double Multiplier(CaseInputs inputs, int row) => _slope?.At(inputs, row) ?? 1;

    public double[] Effects(ReadOnlySpan<double> coordinates)
    {
        if (coordinates.Length != Dimension) throw new ArgumentException("The coordinates do not match the grouping's contrast dimension.", nameof(coordinates));
        var effects = new double[Levels];
        for (int g = 0; g < effects.Length; g++)
            for (int k = 0; k < coordinates.Length; k++) effects[g] += _contrast[g, k] * coordinates[k];
        return effects;
    }
}

internal sealed record RandomEffectAlias(string FirstGrouping, string FirstLevel, string SecondGrouping, string SecondLevel);

/// <summary>The scale-free thin factorization of all constrained grouping indicators, with aliased columns retained on the effects side.</summary>
internal sealed class RandomEffectGeometry
{
    private readonly double[][] _q;
    private readonly double[,] _r;
    private readonly int[] _coordinateGroup, _offsets;
    private readonly IReadOnlyList<RandomEffectGrouping> _groups;

    /// <summary>The factorization of the groupings' constrained indicators, each row scaled by <paramref name="roots"/> where given: W^½'s
    /// for a block under a noise's weights (<see cref="WeightedWhitening"/>), whose aliases the unweighted geometry has reported.</summary>
    public RandomEffectGeometry(IReadOnlyList<RandomEffectGrouping> groups, ReadOnlySpan<double> roots = default)
    {
        ArgumentNullException.ThrowIfNull(groups);
        if (groups.Count == 0) throw new ArgumentException("At least one grouping is required.", nameof(groups));
        ArgumentNullException.ThrowIfNull(groups[0]);
        _groups = Array.AsReadOnly(groups.ToArray());
        Cases = groups[0].Cases;
        var names = new HashSet<(string Grouping, string? Input)>();
        _offsets = new int[groups.Count + 1];
        for (int k = 0; k < groups.Count; k++)
        {
            ArgumentNullException.ThrowIfNull(groups[k]);
            if (groups[k].Cases != Cases) throw new ArgumentException("The groupings have different case counts.", nameof(groups));
            if (!names.Add((groups[k].Name, groups[k].Input))) throw new ArgumentException($"Grouping '{groups[k].Name}' on input '{groups[k].Input}' is declared twice.", nameof(groups));
            _offsets[k + 1] = _offsets[k] + groups[k].Dimension;
        }
        _coordinateGroup = new int[_offsets[^1]];
        var columns = new double[Cases, Dimension];
        for (int k = 0; k < groups.Count; k++)
            for (int j = 0; j < groups[k].Dimension; j++)
            {
                int column = _offsets[k] + j;
                _coordinateGroup[column] = k;
                for (int i = 0; i < Cases; i++) columns[i, column] = groups[k].Multiplier(i)*groups[k].Contrast(groups[k].GroupOf(i), j);
            }
        if (!roots.IsEmpty)
        {
            if (roots.Length != Cases) throw new ArgumentException("Give one root weight per case.", nameof(roots));
            for (int i = 0; i < Cases; i++)
                for (int column = 0; column < Dimension; column++) columns[i, column] *= roots[i];
        }
        ThinSvdResult factor = ThinSvd.Decompose(columns);
        Rank = factor.NumericalRank;
        _q = [.. factor.LeftVectors.Take(Rank)];
        _r = new double[Rank, Dimension];
        for (int i = 0; i < Rank; i++)
            for (int j = 0; j < Dimension; j++) _r[i, j] = factor.SingularValues[i] * factor.RightVectors[i][j];
        var aliases = new List<RandomEffectAlias>();
        for (int a = 0; a < groups.Count && roots.IsEmpty; a++)
            for (int b = a + 1; b < groups.Count; b++)
            {
                if (groups[a].Input is not null || groups[b].Input is not null) continue;
                if (groups[a].ParentName == groups[b].Name || groups[b].ParentName == groups[a].Name) continue;
                for (int g = 0; g < groups[a].Levels; g++)
                    for (int h = 0; h < groups[b].Levels; h++)
                    {
                        if (groups[a].Counts[g] != groups[b].Counts[h]) continue;
                        bool same = true;
                        for (int i = 0; i < Cases && same; i++) same = (groups[a].GroupOf(i) == g) == (groups[b].GroupOf(i) == h);
                        if (same) aliases.Add(new RandomEffectAlias(groups[a].Name, groups[a].LevelName(g), groups[b].Name, groups[b].LevelName(h)));
                    }
            }
        Aliases = aliases.AsReadOnly();
    }

    public int Cases { get; }
    public int Rank { get; }
    public int Dimension => _coordinateGroup.Length;
    public IReadOnlyList<RandomEffectGrouping> Groups => _groups;
    public IReadOnlyList<RandomEffectAlias> Aliases { get; }
    internal double[,] R => _r;
    internal int CoordinateGroup(int coordinate) => _coordinateGroup[coordinate];

    public RandomEffectWhitening At(ReadOnlySpan<double> scales) => new(this, scales);

    internal double[] Project(ReadOnlySpan<double> x)
    {
        if (x.Length != Cases) throw new ArgumentException("The vector does not match the groupings' case count.", nameof(x));
        var projected = new double[Rank];
        for (int j = 0; j < Rank; j++) projected[j] = PairwiseSum.Dot(_q[j], x);
        return projected;
    }

    internal void AddProjected(Span<double> x, ReadOnlySpan<double> projected)
    {
        for (int j = 0; j < Rank; j++)
            for (int i = 0; i < Cases; i++) x[i] += _q[j][i] * projected[j];
    }

    internal double[][] Effects(ReadOnlySpan<double> coordinates)
    {
        var effects = new double[_groups.Count][];
        for (int k = 0; k < effects.Length; k++)
            effects[k] = _groups[k].Effects(coordinates.Slice(_offsets[k], _offsets[k + 1] - _offsets[k]));
        return effects;
    }
}

/// <summary>The joint Gaussian row whitening and conditional group effects at one vector of scales.</summary>
internal sealed class RandomEffectWhitening
{
    private readonly RandomEffectGeometry _geometry;
    private readonly double[] _scales, _roots, _gain, _white;
    private readonly double[][] _e, _v;

    internal RandomEffectWhitening(RandomEffectGeometry geometry, ReadOnlySpan<double> scales)
    {
        if (scales.Length != geometry.Groups.Count) throw new ArgumentException("Give one scale per grouping.", nameof(scales));
        _geometry = geometry;
        _scales = scales.ToArray();
        foreach (double s in _scales)
            if (!(s > 0) || !double.IsFinite(s)) throw new ArgumentOutOfRangeException(nameof(scales), "Each scale is positive and finite.");
        int rank = geometry.Rank;
        if (rank == 0)
        {
            _e = [];
            _v = [];
            _roots = [];
            _gain = [];
            _white = [];
            return;
        }
        double normalizer = _scales.Max();
        // Diagonalize the covariance through its weighted factor. Forming B Bᵀ first would square
        // the condition number and erase a small variance component beside a much larger one.
        var weighted = new double[rank, geometry.Dimension];
        for (int a = 0; a < rank; a++)
            for (int j = 0; j < geometry.Dimension; j++)
                weighted[a, j] = geometry.R[a, j] * (_scales[geometry.CoordinateGroup(j)] / normalizer);
        ThinSvdResult factor = ThinSvd.Decompose(weighted);
        _e = factor.LeftVectors;
        _v = factor.RightVectors;
        _roots = new double[rank];
        _gain = new double[rank];
        _white = new double[rank];
        double logDeterminant = 0;
        for (int j = 0; j < rank; j++)
        {
            double value = factor.SingularValues[j];
            double logValue = value == 0 ? double.NegativeInfinity : Math.Log(normalizer) + Math.Log(value);
            double logVariance = 2 * logValue;
            double logRoot = 0.5 * (logVariance > 0 ? logVariance + SpecialFunctions.Log1p(Math.Exp(-logVariance))
                : SpecialFunctions.Log1p(Math.Exp(logVariance)));
            _roots[j] = Math.Exp(-logRoot);
            _white[j] = _roots[j] - 1;
            _gain[j] = Math.Exp(logValue - 2 * logRoot);
            logDeterminant += 2 * logRoot;
        }
        LogDeterminant = logDeterminant;
    }

    public RandomEffectGeometry Geometry => _geometry;
    public ReadOnlySpan<double> Scales => _scales;
    public double LogDeterminant { get; }

    internal double ScaleInformation(int group)
    {
        if ((uint)group >= (uint)_scales.Length) throw new ArgumentOutOfRangeException(nameof(group));
        int rank = _geometry.Rank;
        var gram = new double[rank, rank];
        for (int j = 0; j < _geometry.Dimension; j++)
        {
            if (_geometry.CoordinateGroup(j) != group) continue;
            var column = new double[rank];
            for (int a = 0; a < rank; a++) column[a] = _scales[group] * _geometry.R[a, j];
            double[] whitened = Transform(column, _roots);
            for (int a = 0; a < rank; a++)
                for (int b = 0; b < rank; b++) gram[a, b] += whitened[a] * whitened[b];
        }
        double information = 0;
        foreach (double entry in gram) information += entry * entry;
        return 2 * information;
    }

    public void Whiten(Span<double> x)
    {
        double[] projected = _geometry.Project(x);
        _geometry.AddProjected(x, Transform(projected, _white));
    }

    public double[][] MeanEffects(ReadOnlySpan<double> residual)
        => _geometry.Effects(ConditionalMean(_geometry.Project(residual)));

    // The posterior square root in standardized contrast coordinates is I + V(diag(roots)-I)Vᵀ.
    // Draw all q coordinates, retaining the prior variation in the q-r aliased directions.
    public double[][] DrawEffects(ReadOnlySpan<double> residual, double sigma2, Xoshiro256PlusPlus rng)
        => _geometry.Effects(DrawCoordinates(residual, sigma2, rng));

    // The draw's contrast coordinates, each a grouping's scale times its standardized coordinate.
    internal double[] DrawCoordinates(ReadOnlySpan<double> residual, double sigma2, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!(sigma2 > 0) || !double.IsFinite(sigma2)) throw new ArgumentOutOfRangeException(nameof(sigma2));
        double sd = Math.Sqrt(sigma2);
        double[] mean = ConditionalMean(_geometry.Project(residual));
        double[] z = [.. Enumerable.Range(0, _geometry.Dimension).Select(_ => Normal.SampleStandard(rng))];
        double[] noise = ConditionalNoise(z);
        for (int j = 0; j < noise.Length; j++) mean[j] += sd * noise[j];
        return mean;
    }

    // u's prior energy over σ², Σ_j (c_j/s_g(j))², the coefficients' Q's counterpart in σ²'s conditional: u | σ² ~ N(0, σ²S²) in contrasts.
    internal double Energy(ReadOnlySpan<double> coordinates)
    {
        double energy = 0;
        for (int j = 0; j < coordinates.Length; j++)
        {
            double standardized = coordinates[j] / _scales[_geometry.CoordinateGroup(j)];
            energy += standardized * standardized;
        }
        return energy;
    }

    // The physical contrast square root at unit σ, applied to supplied standard coordinates.
    internal double[] ConditionalNoise(ReadOnlySpan<double> z)
    {
        if (z.Length != _geometry.Dimension) throw new ArgumentException("Give one standard coordinate per contrast.", nameof(z));
        double[] noise = z.ToArray();
        for (int a = 0; a < _v.Length; a++)
        {
            double along = PairwiseSum.Dot(_v[a], z) * _white[a];
            for (int j = 0; j < noise.Length; j++) noise[j] += _v[a][j] * along;
        }
        for (int j = 0; j < noise.Length; j++) noise[j] *= _scales[_geometry.CoordinateGroup(j)];
        return noise;
    }

    private double[] ConditionalMean(ReadOnlySpan<double> projected)
    {
        var coordinates = new double[_geometry.Dimension];
        for (int a = 0; a < _e.Length; a++)
        {
            double along = PairwiseSum.Dot(_e[a], projected) * _gain[a];
            for (int j = 0; j < coordinates.Length; j++) coordinates[j] += _v[a][j] * along;
        }
        for (int j = 0; j < coordinates.Length; j++) coordinates[j] *= _scales[_geometry.CoordinateGroup(j)];
        return coordinates;
    }

    private double[] Transform(ReadOnlySpan<double> projected, ReadOnlySpan<double> diagonal)
    {
        var result = new double[projected.Length];
        for (int j = 0; j < _e.Length; j++)
        {
            double along = PairwiseSum.Dot(projected, _e[j]) * diagonal[j];
            for (int a = 0; a < result.Length; a++) result[a] += _e[j][a] * along;
        }
        return result;
    }
}

/// <summary>
/// A random-effect block's whitening under observation weights W (E21's slice 4, the block under GBASS's noises). With y's covariance
/// σ²(W⁻¹ + Z̃S²Z̃ᵀ), the map L = V_w^(−½)W^½, V_w = I + W^½Z̃S²Z̃ᵀW^½, takes it to σ²I, at the price −½ln|V_w| beside the weights'
/// ½Σ ln wᵢ. The blocks' conventions centre u with unweighted counts, so L·1 is not W^½·1: the intercept's direction is c = L·1, its mass
/// cᵀc, and the centring of a column or of the response is its GLS mean cᵀLz/cᵀc. At W = I this is the unweighted block's (V1 = 1). V_w is
/// factored by the thin SVD of W^½Z̃ (<see cref="RandomEffectGeometry"/>) at each reading, O(n·q²); u's conditional given the coefficients
/// and σ² is the general block's on the scaled rows.
/// </summary>
internal sealed class WeightedWhitening
{
    private readonly double[] _roots, _direction;
    private readonly RandomEffectWhitening _whitening;

    internal WeightedWhitening(RandomEffectsBlock block, ReadOnlySpan<double> weights)
    {
        ArgumentNullException.ThrowIfNull(block);
        int n = weights.Length;
        if (n != block.Cases) throw new ArgumentException("Give one weight per case.", nameof(weights));
        _roots = new double[n];
        for (int i = 0; i < n; i++) _roots[i] = Math.Sqrt(weights[i]);
        _whitening = new RandomEffectGeometry(block.Layouts, _roots).At(block.Scales.Span);
        var ones = new double[n];
        Array.Fill(ones, 1.0);
        _direction = Apply(ones);
        Mass = PairwiseSum.Dot(_direction, _direction);
    }

    /// <summary>cᵀc, the intercept's mass under L.</summary>
    internal double Mass { get; }

    /// <summary>ln|V_w|.</summary>
    internal double LogDeterminant => _whitening.LogDeterminant;

    /// <summary>c = L·1, the intercept's direction.</summary>
    internal ReadOnlySpan<double> Direction => _direction;

    /// <summary>L·z = V_w^(−½)(W^½∘z), a new array.</summary>
    internal double[] Apply(ReadOnlySpan<double> z)
    {
        var t = new double[z.Length];
        for (int i = 0; i < t.Length; i++) t[i] = _roots[i] * z[i];
        _whitening.Whiten(t);
        return t;
    }

    /// <summary>The GLS centring of a vector already taken through L: cᵀt/cᵀc.</summary>
    internal double CentringOf(ReadOnlySpan<double> applied) => PairwiseSum.Dot(_direction, applied) / Mass;

    /// <summary>The GLS mean of <paramref name="z"/>, 1ᵀC⁻¹z/1ᵀC⁻¹1.</summary>
    internal double Mean(ReadOnlySpan<double> z) => CentringOf(Apply(z));

    /// <summary>u's conditional mean given the residual e = y − σβv − f, through W^½∘e.</summary>
    internal double[][] MeanEffects(ReadOnlySpan<double> residual) => _whitening.MeanEffects(Scaled(residual));

    /// <summary>A draw of u given the residual and σ², with its prior energy over σ² (<see cref="RandomEffectWhitening.Energy"/>).</summary>
    internal (double[][] Effects, double Energy) DrawEffects(ReadOnlySpan<double> residual, double sigma2, Xoshiro256PlusPlus rng)
    {
        double[] coordinates = _whitening.DrawCoordinates(Scaled(residual), sigma2, rng);
        return (_whitening.Geometry.Effects(coordinates), _whitening.Energy(coordinates));
    }

    /// <summary>q, u's contrast coordinates, each drawn: the count u's prior adds to σ²'s conditional.</summary>
    internal int Dimension => _whitening.Geometry.Dimension;

    private double[] Scaled(ReadOnlySpan<double> x)
    {
        var scaled = new double[x.Length];
        for (int i = 0; i < x.Length; i++) scaled[i] = _roots[i] * x[i];
        return scaled;
    }
}

/// <summary>
/// Random effects carried as a latent under a noise (E21's decision 7, amended 2026-11-10): u at the block's scales, which the reading takes
/// as the offset Z̃u beside σβv. Every rung conditions on it, so ladders, SMC and the reference rung are exact; at b = 1 the structure step
/// may integrate it and draw it afresh (<see cref="BassCollapsedStructureBlock"/>). u | σ², s ~ N(0, σ²S²) in each grouping's orthonormal
/// contrasts, so its prior energy over σ² is Σ_k ‖u_k‖²/s_k² on Σ_k q_k coordinates.
/// </summary>
internal sealed record RandomEffectsLatent(RandomEffectsBlock Block, double[][] Effects)
{
    /// <summary>u = 0 at <paramref name="block"/>'s scales: a reading's start.</summary>
    internal static RandomEffectsLatent Zero(RandomEffectsBlock block)
        => new(block, [.. Enumerable.Range(0, block.Count).Select(k => new double[block.Layout(k).Levels])]);

    /// <summary>u from its prior given <paramref name="sigma2"/>, at <paramref name="block"/>'s scales: the reference rung's and SMC's draw.</summary>
    internal static RandomEffectsLatent DrawPrior(RandomEffectsBlock block, double sigma2, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double sigma = Math.Sqrt(sigma2);
        var effects = new double[block.Count][];
        for (int k = 0; k < block.Count; k++)
        {
            RandomEffectGrouping layout = block.Layout(k);
            double sd = sigma * block.Scales.Span[k];
            double[] coordinates = [.. Enumerable.Range(0, layout.Dimension).Select(_ => sd * Normal.SampleStandard(rng))];
            effects[k] = layout.Effects(coordinates);
        }
        return new RandomEffectsLatent(block, effects);
    }

    /// <summary>These effects at <paramref name="block"/>'s scales (a step of s, which moves no effect).</summary>
    internal RandomEffectsLatent At(RandomEffectsBlock block) => this with { Block = block };

    /// <summary>Several outputs' effects U (E34's slice 4, decision 7): [grouping][output][level], u_g ∈ ℝᴷ ~ N(0, s²Σ) in each grouping's
    /// orthonormal contrasts; null for one response, whose effects are <see cref="Effects"/>.</summary>
    internal double[][][]? OutputEffects { get; init; }

    /// <summary>Several outputs' effects <paramref name="effects"/> ([grouping][output][level]) at <paramref name="block"/>'s scales.</summary>
    internal static RandomEffectsLatent OfOutputs(RandomEffectsBlock block, double[][][] effects)
        => new(block, [.. Enumerable.Range(0, block.Count).Select(_ => Array.Empty<double>())]) { OutputEffects = effects };

    /// <summary>U = 0 for <paramref name="outputs"/> outputs at <paramref name="block"/>'s scales: a reading's start.</summary>
    internal static RandomEffectsLatent ZeroOutputs(RandomEffectsBlock block, int outputs)
        => OfOutputs(block, [.. Enumerable.Range(0, block.Count).Select(k => Enumerable.Range(0, outputs).Select(_ => new double[block.Layout(k).Levels]).ToArray())]);

    /// <summary>U from its prior given Σ (its factor <paramref name="sigmaRoot"/>), at <paramref name="block"/>'s scales: each grouping's
    /// contrasts N(0, s_k²·I ⊗ Σ), the reference rung's and SMC's draw.</summary>
    internal static RandomEffectsLatent DrawPriorOutputs(RandomEffectsBlock block, double[,] sigmaRoot, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        int outputs = sigmaRoot.GetLength(0);
        var effects = new double[block.Count][][];
        for (int k = 0; k < block.Count; k++)
        {
            RandomEffectGrouping layout = block.Layout(k);
            int q = layout.Dimension;
            double s = block.Scales.Span[k];
            var z = new double[outputs][];
            for (int l = 0; l < outputs; l++) z[l] = [.. Enumerable.Range(0, q).Select(_ => Normal.SampleStandard(rng))];
            effects[k] = new double[outputs][];
            for (int l = 0; l < outputs; l++)
            {
                var coordinates = new double[q];
                for (int c = 0; c < q; c++)
                {
                    double mixed = 0.0;
                    for (int m = 0; m <= l; m++) mixed += sigmaRoot[l, m] * z[m][c];
                    coordinates[c] = s * mixed;
                }
                effects[k][l] = layout.Effects(coordinates);
            }
        }
        return OfOutputs(block, effects);
    }

    /// <summary>Several outputs' pseudo-rows for Σ's integral: per grouping U_kᵀU_k (K × K, over its levels, its contrasts' cross product), its
    /// q_k contrasts and s_k².</summary>
    internal IEnumerable<(double[,] Cross, int Count, double Scale2)> OutputPseudoRows()
    {
        double[][][] effects = OutputEffects ?? throw new InvalidOperationException("One response's effects have no pseudo-rows.");
        for (int k = 0; k < effects.Length; k++)
        {
            int outputs = effects[k].Length;
            var cross = new double[outputs, outputs];
            for (int a = 0; a < outputs; a++)
                for (int d = 0; d <= a; d++)
                    cross[a, d] = cross[d, a] = PairwiseSum.Dot(effects[k][a], effects[k][d]);
            double s = Block.Scales.Span[k];
            yield return (cross, Block.Layout(k).Dimension, s * s);
        }
    }

    /// <summary>Σ_k ‖u_k‖²/s_k², the prior energy over σ², and Σ_k q_k, its coordinates.</summary>
    internal (double Energy, int Count) Prior
    {
        get
        {
            double energy = 0;
            int count = 0;
            for (int k = 0; k < Effects.Length; k++)
            {
                double s = Block.Scales.Span[k], squares = 0;
                foreach (double u in Effects[k]) squares += u * u;
                energy += squares / (s * s);
                count += Block.Layout(k).Dimension;
            }
            return (energy, count);
        }
    }

    /// <summary>Grouping k's ‖u_k‖².</summary>
    internal double Squares(int k) => Effects[k].Sum(u => u * u);

    /// <summary>Adds Z̃u at each training case to <paramref name="fit"/>; for several outputs Z̃U, case by case and output by output.</summary>
    internal void AddTo(Span<double> fit)
    {
        if (OutputEffects is not { } several)
        {
            Block.AddEffects(Effects, fit);
            return;
        }
        int outputs = several[0].Length, n = Block.Cases;
        if (fit.Length != n * outputs) throw new ArgumentException("Give one value per case and output.", nameof(fit));
        var column = new double[n];
        for (int l = 0; l < outputs; l++)
        {
            Array.Clear(column);
            Block.AddEffects([.. several.Select(g => g[l])], column);
            for (int i = 0; i < n; i++) fit[i * outputs + l] += column[i];
        }
    }

    /// <summary>The effects as named retained records; for several outputs each grouping's <see cref="BassGroupDraw.OutputEffects"/>.</summary>
    internal BassGroupDraw[] Groups => OutputEffects is { } several ? Block.NamedOutputs(several) : Block.Named(Effects);
}
