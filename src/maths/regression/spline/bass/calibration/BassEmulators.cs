using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// A scalar or direct functional <see cref="BassFit"/> as a calibration emulator (C-X2): each retained draw's latent
/// surface at θ, with any fixed inputs the observations carry, and σ²_e·I as its noise (D09's amendment).
/// <list type="bullet">
/// <item><b>Inputs.</b> R's calibration sets every emulator input from θ (<c>calib.R</c>'s <c>func</c> predicts at the
/// transformed θ alone). Here each predictor is either a calibration parameter, matched by name, or given by
/// <c>fixedInputs</c>, whose rows are the observations' cases (one case when none is fixed).</item>
/// <item><b>Outputs.</b> m = cases × locations, case-major: a scalar fit's locations are 1, a direct functional fit's are
/// the declared output locations, its training locations by default.</item>
/// <item><b>Prediction</b> equals <see cref="BassFit.Predict"/>'s latent draw bit for bit. Each term's continuous part is
/// evaluated once per θ and its categorical part per level combination (Francom et al. 2019).</item>
/// </list>
/// </summary>
public sealed class BassFitEmulator : ICalibrationEmulator
{
    private readonly List<BassDraw> _draws;
    private readonly ThetaSurface _surface;
    private readonly double[]? _weights;

    /// <param name="fit">The fit.</param>
    /// <param name="data">The data the fit was run on; other data are refused by fingerprint.</param>
    /// <param name="parameters">The calibration parameters, each naming one of the fit's predictors.</param>
    /// <param name="fixedInputs">The other predictors' values, one row per observation case; null when every predictor
    /// is a parameter.</param>
    /// <param name="outputLocations">A direct functional fit's output locations; null for the training locations.</param>
    public BassFitEmulator(BassFit fit, BassData data, IReadOnlyList<CalibrationParameter> parameters, PredictionCases? fixedInputs = null,
                           IReadOnlyList<ContinuousColumn>? outputLocations = null)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(data);
        string fingerprint = Bass.Fingerprint(fit.Config, data);
        if (fingerprint != fit.Identity.Fingerprint)
            throw new ArgumentException($"The fit was run on other data: its fingerprint is {fit.Identity.Fingerprint}, theirs {fingerprint}.", nameof(data));
        BassModel model = BassModel.Create(fit.Config, data);
        (_draws, _) = BassDrawFit.Pool(fit.Draws, fit.Terms, false, nameof(fit));
        if (fit.DrawWeights is { } w)
        {
            if (w.Length != _draws.Count) throw new ArgumentException($"The fit's {w.Length} weights do not match its {_draws.Count} draws.", nameof(fit));
            _weights = (double[])w.Clone();
        }

        OutputInputs? outputs = null;
        if (data.Response is FunctionalResponse functional)
        {
            IReadOnlyList<ContinuousColumn> given = outputLocations
                ?? Array.ConvertAll([.. functional.OutputVariables], v => new ContinuousColumn(v.Name, v.Values.ToArray()));
            (outputs, _) = PredictionCases.PrepareLocations(model.Transform, given);
        }
        else if (outputLocations is not null) throw new ArgumentException("A scalar fit has no output locations.", nameof(outputLocations));

        _surface = new ThetaSurface(model, fit.Terms, parameters, fixedInputs, outputs);
        Parameters = _surface.Parameters;
        Outputs = _surface.Cases * _surface.Locations;
    }

    /// <inheritdoc/>
    public IReadOnlyList<CalibrationParameter> Parameters { get; }

    /// <inheritdoc/>
    public int DrawCount => _draws.Count;

    /// <inheritdoc/>
    public IReadOnlyList<double>? DrawWeights => _weights;

    /// <inheritdoc/>
    public int Outputs { get; }

    /// <summary>The observations' cases.</summary>
    public int Cases => _surface.Cases;

    /// <summary>The output locations per case: 1 for a scalar fit.</summary>
    public int Locations => _surface.Locations;

    /// <summary>The continuous parameters whose declared bounds leave the fit's training range, where the hinges
    /// extrapolate linearly.</summary>
    public IReadOnlyList<string> Extrapolating => _surface.Extrapolating;

    /// <inheritdoc/>
    public EmulatorNoise Noise(int draw) => new IndependentNoise(Draw(draw).Sigma2);

    /// <inheritdoc/>
    public double[,]? TruncationResiduals() => null;

    /// <inheritdoc/>
    public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output)
    {
        if (output.Length != Outputs) throw new ArgumentException($"The output needs {Outputs} values.", nameof(output));
        BassDraw d = Draw(draw);
        _surface.Evaluate(d.Structure.TermIds, d.Coefficients, continuous, [levels.ToArray()], output);
    }

    /// <inheritdoc/>
    public void PredictLevels(int draw, ReadOnlySpan<double> continuous, IReadOnlyList<int[]> combinations, double[,] output)
    {
        ArgumentNullException.ThrowIfNull(combinations);
        ArgumentNullException.ThrowIfNull(output);
        if (output.GetLength(0) != combinations.Count || output.GetLength(1) != Outputs)
            throw new ArgumentException($"The output needs one row per combination and {Outputs} columns.", nameof(output));
        BassDraw d = Draw(draw);
        var flat = new double[combinations.Count * Outputs];
        _surface.Evaluate(d.Structure.TermIds, d.Coefficients, continuous, combinations, flat);
        Buffer.BlockCopy(flat, 0, output, 0, flat.Length * sizeof(double));
    }

    private BassDraw Draw(int draw)
    {
        if ((uint)draw >= (uint)_draws.Count) throw new ArgumentOutOfRangeException(nameof(draw), $"The fit holds {_draws.Count} draws.");
        return _draws[draw];
    }
}

/// <summary>
/// A <see cref="BassBasisFit"/> as a calibration emulator (C-X2): joint draw s's components at θ reconstructed,
/// y_m + y_s ∘ K w_s(θ) (D07), with the component noise (y_s ∘ K)·diag(σ²_{s,k})·(y_s ∘ K)ᵀ and the decomposition's
/// truncation residuals y_s ∘ eᵢ. The basis exists at its training locations, so the observations are one case; a
/// predictor that is not a parameter is fixed by <c>fixedInputs</c>'s one row. Prediction equals
/// <see cref="BassBasisFit.Predict"/>'s latent draw bit for bit.
/// </summary>
public sealed class BassBasisEmulator : ICalibrationEmulator
{
    private readonly BassBasisFit _fit;
    private readonly List<BassDraw>[] _pooled;
    private readonly ThetaSurface[] _surfaces;
    private readonly double[,] _factor;

    /// <param name="fit">The basis fit.</param>
    /// <param name="parameters">The calibration parameters, each naming one of the fit's predictors.</param>
    /// <param name="fixedInputs">The other predictors' values, one row; null when every predictor is a parameter.</param>
    public BassBasisEmulator(BassBasisFit fit, IReadOnlyList<CalibrationParameter> parameters, PredictionCases? fixedInputs = null)
    {
        ArgumentNullException.ThrowIfNull(fit);
        if (fixedInputs is not null && fixedInputs.Cases != 1)
            throw new ArgumentException("A basis emulator's observations are one case: its fixed inputs take one row.", nameof(fixedInputs));
        _fit = fit;
        int r = fit.Decomposition.Rank;
        _pooled = new List<BassDraw>[r];
        _surfaces = new ThetaSurface[r];
        for (int k = 0; k < r; k++)
        {
            BassFit component = fit.Components[k];
            (_pooled[k], _) = BassDrawFit.Pool(component.Draws, component.Terms, false, nameof(fit));
            _surfaces[k] = new ThetaSurface(fit.Model(k), component.Terms, parameters, fixedInputs, null);
        }
        Parameters = _surfaces[0].Parameters;
        LinearReconstruction map = fit.Decomposition.Map;
        int m = map.Locations;
        _factor = new double[m, r];
        for (int l = 0; l < m; l++)
            for (int k = 0; k < r; k++) _factor[l, k] = map.Scales[l] * map.Basis(l, k);
    }

    /// <inheritdoc/>
    public IReadOnlyList<CalibrationParameter> Parameters { get; }

    /// <inheritdoc/>
    public int DrawCount => _fit.JointDraws.Count;

    /// <inheritdoc/>
    public IReadOnlyList<double>? DrawWeights => null;

    /// <inheritdoc/>
    public int Outputs => _fit.Decomposition.Locations;

    /// <summary>The continuous parameters whose declared bounds leave the fit's training range.</summary>
    public IReadOnlyList<string> Extrapolating => _surfaces[0].Extrapolating;

    /// <inheritdoc/>
    public EmulatorNoise Noise(int draw)
    {
        Check(draw);
        var variances = new double[_pooled.Length];
        for (int k = 0; k < variances.Length; k++) variances[k] = _pooled[k][_fit.JointDraws.Draw(k, draw)].Sigma2;
        return new ComponentNoise(_factor, variances);
    }

    /// <inheritdoc/>
    public double[,]? TruncationResiduals()
    {
        LinearReconstruction map = _fit.Decomposition.Map;
        double[,] e = _fit.Decomposition.Truncation.Values;
        int n = e.GetLength(0), m = e.GetLength(1);
        var physical = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) physical[i, l] = map.Scales[l] * e[i, l];
        return physical;
    }

    /// <inheritdoc/>
    public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output)
    {
        if (output.Length != Outputs) throw new ArgumentException($"The output needs {Outputs} values.", nameof(output));
        int[][] one = [levels.ToArray()];
        double[] w = Coefficients(draw, continuous, one);
        _fit.Decomposition.Map.Reconstruct(w, [], output);
    }

    /// <inheritdoc/>
    public void PredictLevels(int draw, ReadOnlySpan<double> continuous, IReadOnlyList<int[]> combinations, double[,] output)
    {
        ArgumentNullException.ThrowIfNull(combinations);
        ArgumentNullException.ThrowIfNull(output);
        int m = Outputs, r = _pooled.Length, count = combinations.Count;
        if (output.GetLength(0) != count || output.GetLength(1) != m)
            throw new ArgumentException($"The output needs one row per combination and {m} columns.", nameof(output));
        double[] w = Coefficients(draw, continuous, combinations);   // combination-major, r per combination
        var y = new double[m];
        for (int c = 0; c < count; c++)
        {
            _fit.Decomposition.Map.Reconstruct(w.AsSpan(c * r, r), [], y);
            for (int l = 0; l < m; l++) output[c, l] = y[l];
        }
    }

    // Each combination's coefficients w_k(θ), combination-major.
    private double[] Coefficients(int draw, ReadOnlySpan<double> continuous, IReadOnlyList<int[]> combinations)
    {
        Check(draw);
        int r = _pooled.Length, count = combinations.Count;
        var w = new double[count * r];
        var f = new double[count];
        for (int k = 0; k < r; k++)
        {
            BassDraw d = _pooled[k][_fit.JointDraws.Draw(k, draw)];
            _surfaces[k].Evaluate(d.Structure.TermIds, d.Coefficients, continuous, combinations, f);
            for (int c = 0; c < count; c++) w[c * r + k] = f[c];
        }
        return w;
    }

    private void Check(int draw)
    {
        if ((uint)draw >= (uint)_fit.JointDraws.Count) throw new ArgumentOutOfRangeException(nameof(draw), $"The fit holds {_fit.JointDraws.Count} joint draws.");
    }
}

/// <summary>
/// One fit's surface over θ: the cases the observations carry, with θ's columns set per call, evaluated as
/// <see cref="BassSurface"/> evaluates them (its accumulation, so the values are bit for bit the same). Each term is split
/// once into its part on θ's categorical variables and the rest, so that a call evaluates the rest once per θ and the
/// categorical part, an indicator, per level combination. Safe for concurrent calls.
/// </summary>
internal sealed class ThetaSurface
{
    private readonly BassBasis _basis;
    private readonly BassTermTable _terms;
    private readonly OutputInputs? _outputs;
    private readonly double[][] _continuous;          // template columns; θ's are replaced per call
    private readonly int[][] _categorical;
    private readonly (int Column, RangeScale Range)[] _thetaContinuous;
    private readonly int[] _thetaCategoricalVariable; // per categorical parameter, its variable among the fit's categoricals
    private readonly int[][] _thetaLevels;            // per categorical parameter, each declared level's training index
    private readonly ConcurrentDictionary<int, Split> _splits = new();

    private sealed record Split(BassTerm? Rest, (int Parameter, CategoricalFactor Factor)[] OnTheta, double[]? Output);

    public ThetaSurface(BassModel model, BassTermTable terms, IReadOnlyList<CalibrationParameter> parameters, PredictionCases? fixedInputs,
                        OutputInputs? outputs)
    {
        ArgumentNullException.ThrowIfNull(parameters);
        var layout = new ParameterLayout(parameters);
        Parameters = [.. parameters];
        _basis = model.Basis;
        _terms = terms;
        _outputs = outputs;
        BassInputTransform transform = model.Transform;

        var byName = new Dictionary<string, CalibrationParameter>(StringComparer.Ordinal);
        foreach (CalibrationParameter p in parameters) byName.Add(p.Name, p);
        var fixedNames = new HashSet<string>(StringComparer.Ordinal);
        if (fixedInputs is not null)
            foreach (PredictionColumn c in fixedInputs.Columns)
            {
                if (byName.ContainsKey(c.Name)) throw new ArgumentException($"'{c.Name}' is a calibration parameter and a fixed input.", nameof(fixedInputs));
                fixedNames.Add(c.Name);
            }
        var predictorNames = new HashSet<string>(StringComparer.Ordinal);
        foreach (PredictorScale scale in transform.Predictors)
        {
            predictorNames.Add(scale.Name);
            if (!byName.ContainsKey(scale.Name) && !fixedNames.Contains(scale.Name))
                throw new ArgumentException($"The predictor '{scale.Name}' is neither a calibration parameter nor a fixed input.", nameof(parameters));
        }
        foreach (string name in byName.Keys)
            if (!predictorNames.Contains(name)) throw new ArgumentException($"The parameter '{name}' names no predictor of the fit.", nameof(parameters));
        Cases = fixedInputs?.Cases ?? 1;

        // The template: the fixed columns and θ's at placeholder values, prepared as prediction prepares cases.
        var columns = new List<PredictionColumn>();
        if (fixedInputs is not null) columns.AddRange(fixedInputs.Columns);
        var extrapolating = new List<string>();
        var thetaContinuous = new (int, RangeScale)[layout.Continuous.Count];
        var thetaVariable = new int[layout.Categorical.Count];
        var thetaLevels = new int[layout.Categorical.Count][];
        for (int j = 0; j < transform.Predictors.Count; j++)
        {
            PredictorScale scale = transform.Predictors[j];
            if (!byName.TryGetValue(scale.Name, out CalibrationParameter? p)) continue;
            (bool isContinuous, int index) = model.Inputs.PredictorPosition(j);
            switch (p, scale)
            {
                case (ContinuousParameter c, RangeScale range) when isContinuous:
                    thetaContinuous[IndexOf(layout.Continuous, c)] = (index, range);
                    if (!(range.Contains(c.Lower) && range.Contains(c.Upper))) extrapolating.Add(c.Name);
                    columns.Add(new ContinuousColumn(c.Name, [.. Enumerable.Repeat(c.Lower, Cases)]));
                    break;
                case (CategoricalParameter f, TrainingLevels training) when !isContinuous:
                    int q = IndexOf(layout.Categorical, f);
                    thetaVariable[q] = index;
                    thetaLevels[q] = new int[f.Levels.Length];
                    for (int l = 0; l < f.Levels.Length; l++)
                    {
                        thetaLevels[q][l] = training.IndexOf(f.Levels[l]);
                        if (thetaLevels[q][l] < 0)
                            throw new ArgumentException($"The level '{f.Levels[l]}' of '{f.Name}' was taken by no training case, so the emulator has no value there.", nameof(parameters));
                    }
                    columns.Add(new CategoricalColumn(f.Name, [.. Enumerable.Repeat(f.Levels[0], Cases)]));
                    break;
                default:
                    throw new ArgumentException($"The parameter '{p.Name}' is not of its predictor's kind.", nameof(parameters));
            }
        }
        (CaseInputs template, _) = new PredictionCases(columns).Prepare(transform, model.Inputs, UnseenLevel.Refuse);
        _continuous = new double[template.ContinuousCount][];
        for (int v = 0; v < _continuous.Length; v++) _continuous[v] = template.Continuous(v).ToArray();
        _categorical = new int[template.CategoricalCount][];
        for (int v = 0; v < _categorical.Length; v++) _categorical[v] = template.Categorical(v).ToArray();
        _thetaContinuous = thetaContinuous;
        _thetaCategoricalVariable = thetaVariable;
        _thetaLevels = thetaLevels;
        Extrapolating = extrapolating;
    }

    public IReadOnlyList<CalibrationParameter> Parameters { get; }

    public int Cases { get; }

    public int Locations => _outputs?.Locations ?? 1;

    public IReadOnlyList<string> Extrapolating { get; }

    /// <summary>
    /// The function of <paramref name="coefficients"/> (the intercept, then one per id of <paramref name="ids"/>) at θ's
    /// continuous values and each of <paramref name="combinations"/>' levels, into <paramref name="f"/>:
    /// combination-major, then case by case, then location.
    /// </summary>
    public void Evaluate(ReadOnlySpan<int> ids, ReadOnlySpan<double> coefficients, ReadOnlySpan<double> continuous, IReadOnlyList<int[]> combinations,
                         Span<double> f)
    {
        if (continuous.Length != _thetaContinuous.Length) throw new ArgumentException($"θ has {_thetaContinuous.Length} continuous values.", nameof(continuous));
        int cases = Cases, locations = Locations, block = cases * locations;
        if (f.Length != combinations.Count * block) throw new ArgumentException("The output needs one block per combination.", nameof(f));
        foreach (int[] levels in combinations) CheckLevels(levels);

        var columns = (double[][])_continuous.Clone();
        for (int i = 0; i < _thetaContinuous.Length; i++)
        {
            if (!double.IsFinite(continuous[i])) throw new ArgumentException("θ's continuous values must be finite.", nameof(continuous));
            (int column, RangeScale range) = _thetaContinuous[i];
            var scaled = new double[cases];
            Array.Fill(scaled, range.Scale(continuous[i]));
            columns[column] = scaled;
        }
        var at = new CaseInputs(columns, _categorical);

        var rest = new double[ids.Length][];
        var splits = new Split[ids.Length];
        for (int q = 0; q < ids.Length; q++)
        {
            splits[q] = SplitOf(ids[q]);
            rest[q] = new double[cases];
            if (splits[q].Rest is { } term) _basis.CaseColumn(term, at, rest[q]);
            else Array.Fill(rest[q], 1.0);
        }

        var a = new double[cases];
        for (int c = 0; c < combinations.Count; c++)
        {
            int[] levels = combinations[c];
            Span<double> g = f.Slice(c * block, block);
            g.Fill(coefficients[0]);
            for (int q = 0; q < ids.Length; q++)
            {
                bool inside = true;
                foreach ((int parameter, CategoricalFactor factor) in splits[q].OnTheta)
                    if (!factor.Contains(_thetaLevels[parameter][levels[parameter]])) { inside = false; break; }
                for (int i = 0; i < cases; i++) a[i] = inside ? rest[q][i] : 0.0;
                double coefficient = coefficients[q + 1];
                double[]? h = splits[q].Output;
                if (h is null)
                {
                    for (int i = 0; i < cases; i++) g[i] += coefficient * a[i];
                    continue;
                }
                for (int i = 0; i < cases; i++)
                {
                    double ca = coefficient * a[i];
                    if (ca == 0.0) continue;
                    Span<double> row = g.Slice(i * locations, locations);
                    for (int l = 0; l < locations; l++) row[l] += ca * h[l];
                }
            }
        }
    }

    private void CheckLevels(int[] levels)
    {
        if (levels is null || levels.Length != _thetaLevels.Length)
            throw new ArgumentException($"θ has {_thetaLevels.Length} categorical levels.", nameof(levels));
        for (int q = 0; q < levels.Length; q++)
            if ((uint)levels[q] >= (uint)_thetaLevels[q].Length)
                throw new ArgumentOutOfRangeException(nameof(levels), $"Categorical parameter {q} has {_thetaLevels[q].Length} levels.");
    }

    // Term `id` without its factors on θ's categorical variables, those factors by parameter, and its output column.
    private Split SplitOf(int id) => _splits.GetOrAdd(id, key =>
    {
        BassTerm term = _terms[key];
        try
        {
            _basis.Validate(term);
        }
        catch (ArgumentException e)
        {
            throw new ArgumentException($"Term {key} ({term}) does not fit the model's inputs: {e.Message}", "terms", e);
        }
        var kept = new List<CategoricalFactor>();
        var onTheta = new List<(int, CategoricalFactor)>();
        foreach (CategoricalFactor factor in term.Categorical)
        {
            int parameter = Array.IndexOf(_thetaCategoricalVariable, factor.Variable);
            if (parameter >= 0) onTheta.Add((parameter, factor));
            else kept.Add(factor);
        }
        BassTerm? rest = onTheta.Count == 0 ? term
                       : term.Continuous.Count + kept.Count + term.Functional.Count == 0 ? null
                       : new BassTerm(term.Continuous, kept, term.Functional);
        double[]? h = null;
        if (_outputs is { } outputs)
        {
            h = new double[outputs.Locations];
            _basis.OutputColumn(term, outputs, h);
        }
        return new Split(rest, [.. onTheta], h);
    });

    private static int IndexOf<T>(IReadOnlyList<T> list, T item) where T : class
    {
        for (int i = 0; i < list.Count; i++)
            if (ReferenceEquals(list[i], item)) return i;
        throw new InvalidOperationException("The parameter is not in the layout.");
    }
}
