using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>One retained grouping's scale and effects, in the observed levels' declaration order.</summary>
public sealed record BassGroupDraw(string Grouping, double Scale, double[] Effects);

/// <summary>Named Gaussian intercept blocks on one set of cases, jointly integrated through their row covariance.</summary>
public sealed class RandomEffectsBlock
{
    private readonly IReadOnlyList<BassRandomIntercept> _specifications;
    private readonly IReadOnlyList<CategoricalPredictor> _groupings;
    private readonly RandomEffectGrouping[] _layouts;
    private readonly RandomEffectGeometry? _geometry;
    private readonly RandomEffectWhitening? _general;
    private readonly RandomInterceptBlock? _single;
    private readonly double[] _scales;

    private RandomEffectsBlock(BassRandomIntercept[] specifications, CategoricalPredictor[] groupings, RandomEffectGrouping[] layouts,
                               RandomEffectGeometry? geometry, RandomInterceptBlock? single = null)
    {
        _specifications = Array.AsReadOnly(specifications);
        _groupings = Array.AsReadOnly(groupings);
        _layouts = layouts;
        _geometry = geometry;
        _scales = [.. specifications.Select(s => s.Scale)];
        if (geometry is null) _single = single ?? RandomInterceptBlock.Of(groupings[0], specifications[0]);
        else _general = geometry.At(_scales);
        LogScalePrior = _single?.LogScalePrior ?? specifications.Sum(s => s.ScalePrior is { } prior ? prior.LogDensity(s.Scale) + Math.Log(s.Scale) : 0);
    }

    public static RandomEffectsBlock Of(IReadOnlyList<CategoricalPredictor> groupings, IReadOnlyList<BassRandomIntercept> specifications)
    {
        ArgumentNullException.ThrowIfNull(groupings);
        ArgumentNullException.ThrowIfNull(specifications);
        BassConfig.CheckRandomIntercepts(specifications);
        if (specifications.Count == 0) throw new ArgumentException("At least one random-intercept block is required.", nameof(specifications));
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
            BassRandomIntercept spec = specifications[k];
            if (!available.TryGetValue(spec.Grouping, out CategoricalPredictor? grouping))
                throw new ArgumentException($"Random-intercept block '{spec.Grouping}' has no data grouping.", nameof(specifications));
            columns[k] = grouping;
            if (spec.Convention == RandomInterceptConvention.ObservationCentred)
                layouts[k] = RandomEffectGrouping.ObservationCentred(grouping);
            else
            {
                if (!available.TryGetValue(spec.ParentGrouping!, out CategoricalPredictor? parent))
                    throw new ArgumentException($"Parent grouping '{spec.ParentGrouping}' of '{spec.Grouping}' is missing.", nameof(specifications));
                layouts[k] = RandomEffectGrouping.WithinParent(grouping, parent);
            }
        }
        RandomEffectGeometry? geometry = specifications.Count == 1 && specifications[0].Convention == RandomInterceptConvention.ObservationCentred
            ? null : new RandomEffectGeometry(layouts);
        return new RandomEffectsBlock(specifications.ToArray(), columns, layouts, geometry);
    }

    public IReadOnlyList<BassRandomIntercept> Specifications => _specifications;
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
        if (scales.Length != Count) throw new ArgumentException("Give one scale per random-intercept block.", nameof(scales));
        if (scales.SequenceEqual(_scales)) return this;
        var specifications = new BassRandomIntercept[Count];
        for (int k = 0; k < Count; k++)
        {
            if (!(scales[k] > 0) || !double.IsFinite(scales[k])) throw new ArgumentOutOfRangeException(nameof(scales));
            specifications[k] = _specifications[k] with { Scale = scales[k] };
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
        var named = new Dictionary<string, BassGroupDraw>(StringComparer.Ordinal);
        foreach (BassGroupDraw group in groups)
        {
            ArgumentNullException.ThrowIfNull(group);
            if (!named.TryAdd(group.Grouping, group)) throw new ArgumentException($"The draw repeats grouping '{group.Grouping}'.", nameof(groups));
        }
        var ordered = new BassGroupDraw[Count];
        for (int k = 0; k < Count; k++)
        {
            string name = _specifications[k].Grouping;
            if (!named.TryGetValue(name, out BassGroupDraw? group)) throw new ArgumentException($"The draw has no effects for grouping '{name}'.", nameof(groups));
            if (group.Effects is null || group.Effects.Length != _layouts[k].Levels)
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
        if (!working.IsUnweighted || working.Latents is not null || !working.Offset.IsEmpty || working.Response is not ScalarResponse)
            throw new ArgumentException("Random intercepts whiten a scalar response under unit weights, with no offset or latents.", nameof(working));
        return WhitenResponse(working);
    }

    // Replace the old whitening from raw observations, preserving the design's response and prior-reference identities.
    internal WorkingResponse WhitenResponse(WorkingResponse working)
    {
        ReadOnlySpan<double> y = working.Response.Values;
        double[] white = y.ToArray();
        Whiten(white);
        var offset = new double[y.Length];
        for (int i = 0; i < y.Length; i++) offset[i] = y[i] - white[i];
        return new WorkingResponse(working.Response, offset, output: working.Output, meaning: working.Meaning, priorGram: working.PriorGram,
                                   reference: working.Reference, logJacobian: -LogDeterminant) { Whitening = this };
    }

    internal double ScaleInformation(int group) => _single?.LogScaleInformation ?? _general!.ScaleInformation(group);

    internal double[][] MeanEffects(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients)
        => _single is { } one ? [one.EffectMean(one.ResidualSums(model, terms, coefficients))] : _general!.MeanEffects(Residual(model, terms, coefficients));

    internal double[][] DrawEffects(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients,
                                     double sigma2, Xoshiro256PlusPlus rng)
        => _single is { } one ? [one.DrawEffects(one.ResidualSums(model, terms, coefficients), sigma2, rng)]
            : _general!.DrawEffects(Residual(model, terms, coefficients), sigma2, rng);

    internal BassGroupDraw[] DrawGroups(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients,
                                       double sigma2, Xoshiro256PlusPlus rng)
    {
        double[][] effects = DrawEffects(model, terms, coefficients, sigma2, rng);
        var groups = new BassGroupDraw[Count];
        for (int k = 0; k < Count; k++) groups[k] = new BassGroupDraw(_specifications[k].Grouping, _scales[k], effects[k]);
        return groups;
    }

    private static double[] Residual(BassModel model, IReadOnlyList<BassTerm> terms, ReadOnlySpan<double> coefficients)
    {
        if (coefficients.Length != terms.Count + 1) throw new ArgumentException("Give an intercept and one coefficient per term.", nameof(coefficients));
        double[] residual = model.Data.Response.Values.ToArray();
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

    private RandomEffectGrouping(CategoricalPredictor grouping, CategoricalPredictor? parent)
    {
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

    public string Name { get; }
    public string? ParentName { get; }
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

    public RandomEffectGeometry(IReadOnlyList<RandomEffectGrouping> groups)
    {
        ArgumentNullException.ThrowIfNull(groups);
        if (groups.Count == 0) throw new ArgumentException("At least one grouping is required.", nameof(groups));
        ArgumentNullException.ThrowIfNull(groups[0]);
        _groups = Array.AsReadOnly(groups.ToArray());
        Cases = groups[0].Cases;
        var names = new HashSet<string>(StringComparer.Ordinal);
        _offsets = new int[groups.Count + 1];
        for (int k = 0; k < groups.Count; k++)
        {
            ArgumentNullException.ThrowIfNull(groups[k]);
            if (groups[k].Cases != Cases) throw new ArgumentException("The groupings have different case counts.", nameof(groups));
            if (!names.Add(groups[k].Name)) throw new ArgumentException($"Grouping '{groups[k].Name}' is declared twice.", nameof(groups));
            _offsets[k + 1] = _offsets[k] + groups[k].Dimension;
        }
        _coordinateGroup = new int[_offsets[^1]];
        var columns = new double[Cases, Dimension];
        for (int k = 0; k < groups.Count; k++)
            for (int j = 0; j < groups[k].Dimension; j++)
            {
                int column = _offsets[k] + j;
                _coordinateGroup[column] = k;
                for (int i = 0; i < Cases; i++) columns[i, column] = groups[k].Contrast(groups[k].GroupOf(i), j);
            }
        ThinSvdResult factor = ThinSvd.Decompose(columns);
        Rank = factor.NumericalRank;
        _q = [.. factor.LeftVectors.Take(Rank)];
        _r = new double[Rank, Dimension];
        for (int i = 0; i < Rank; i++)
            for (int j = 0; j < Dimension; j++) _r[i, j] = factor.SingularValues[i] * factor.RightVectors[i][j];
        var aliases = new List<RandomEffectAlias>();
        for (int a = 0; a < groups.Count; a++)
            for (int b = a + 1; b < groups.Count; b++)
            {
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
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!(sigma2 > 0) || !double.IsFinite(sigma2)) throw new ArgumentOutOfRangeException(nameof(sigma2));
        double sd = Math.Sqrt(sigma2);
        double[] mean = ConditionalMean(_geometry.Project(residual));
        double[] z = [.. Enumerable.Range(0, _geometry.Dimension).Select(_ => Normal.SampleStandard(rng))];
        double[] noise = ConditionalNoise(z);
        for (int j = 0; j < noise.Length; j++) mean[j] += sd * noise[j];
        return _geometry.Effects(mean);
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
