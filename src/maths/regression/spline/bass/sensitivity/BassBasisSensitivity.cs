using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// Basis sensitivity (D08, O7; R's <c>sobolBasis</c>): the variance of a basis fit's reconstructed latent output,
/// y_m(ℓ) + y_s(ℓ)·Σₖ K_ℓk·wₖ(x), under a declared product measure over the inputs, decomposed for each joint posterior
/// draw into Sobol main and total effects, per output location and integrated over the locations. Each component's surface
/// wₖ = aₖ₀ + Σⱼ aₖⱼBₖⱼ is a sum of term products, so the components' covariances on every subset u,
/// Cov(wₖ,u, wₖ′,u) = Σⱼₗ aₖⱼaₖ′ₗ·Cov_u(Bₖⱼ, Bₖ′ₗ), come exactly from <see cref="BassMoments"/>' term pairs, across
/// components as within one: fitting the components independently does not remove the cross terms. Those r × r matrices are
/// free of ℓ and are kept per draw; each form weighs them by g_ℓk = y_s(ℓ)·K_ℓk, so variances are in output units. The
/// truncation residual and the noise are outside the decomposition, as observation noise is for direct fits.
/// </summary>
public static class BassBasisSensitivity
{
    /// <summary>
    /// The decomposition of <paramref name="fit"/>'s joint draws under <paramref name="spec"/>. The spec's measures, order,
    /// budget and tolerance apply as for direct fits; the order counts variables, the location as one in the integrated form,
    /// whose default order is therefore one more. <see cref="SensitivitySpec.Functional"/> is ignored, since both forms are
    /// views of one analysis. Refused, with its reason: output coordinates or output-variable measures (the basis exists
    /// only at its training locations, which are integrated equally weighted), components fitted at different hinge degrees,
    /// and the spec's own refusals.
    /// </summary>
    public static BasisSobolAnalysis Of(BassBasisFit fit, SensitivitySpec? spec = null)
    {
        ArgumentNullException.ThrowIfNull(fit);
        spec ??= new SensitivitySpec();
        if (!(spec.RoundingTolerance >= 0.0) || double.IsPositiveInfinity(spec.RoundingTolerance))
            throw new ArgumentOutOfRangeException(nameof(spec), "The rounding tolerance must be nonnegative and finite.");
        if (spec.Surface != DecomposedSurface.Mean) throw new ArgumentException($"The decomposed surface {spec.Surface} is not built.", nameof(spec));
        if (spec.Coordinates is not null)
            throw new ArgumentException("A basis exists only at its training locations, which have no output coordinates: read a location's decomposition by its index.", nameof(spec));
        if (spec.OutputMeasures is { Count: > 0 })
            throw new ArgumentException("A basis fit's locations are integrated equally weighted, its training grid; output-variable measures apply to direct functional fits.", nameof(spec));

        int r = fit.Components.Count;
        int degree = fit.Components[0].Config.Structure.Degree;
        for (int k = 1; k < r; k++)
            if (fit.Components[k].Config.Structure.Degree != degree)
                throw new ArgumentException($"Component {k} was fitted at hinge degree {fit.Components[k].Config.Structure.Degree} and component 0 at {degree}: " +
                                            "the cross-component moments need one degree.", nameof(fit));
        BassFit first = fit.Components[0];
        BassInputs inputs = BassInputs.Prepare(fit.ComponentData(0), first.Transform);
        ProductMeasure measure = ProductMeasure.Resolve(first.Transform, inputs, spec.Measures);
        BassMoments moments = BassMoments.On(inputs, degree, measure);

        // Each component's draws pooled in chain order; joint draw d takes pooled[k][JointDraws.Draw(k, d)].
        BassJointDraws joint = fit.JointDraws;
        int n = joint.Count;
        var pooled = new List<BassDraw>[r];
        for (int k = 0; k < r; k++) pooled[k] = BassDrawFit.Pool(fit.Components[k].Draws, fit.Components[k].Terms, false, nameof(fit)).Pooled;
        var held = new HashSet<BassTerm>();
        for (int k = 0; k < r; k++)
            for (int d = 0; d < n; d++)
                foreach (int id in pooled[k][joint.Draw(k, d)].Structure.TermIds) held.Add(fit.Components[k].Terms[id]);
        int largest = held.Select(t => moments.Variables(t).Count).DefaultIfEmpty(0).Max();
        int maxOrder = spec.MaxOrder ?? Math.Max(1, largest);
        int integratedOrder = spec.MaxOrder ?? largest + 1;
        int[][] subsets = BassSensitivity.Subsets(held.Select(t => moments.Variables(t)), maxOrder, spec.SubsetBudget);
        var index = new Dictionary<int[], int>(subsets.Length, BassSensitivity.SubsetComparer.Instance);
        for (int u = 0; u < subsets.Length; u++) index[subsets[u]] = u;

        var plans = new Dictionary<(BassTerm, BassTerm), PairPlan>();
        PairPlan Plan(BassTerm a, BassTerm b)
        {
            if (plans.TryGetValue((a, b), out PairPlan? plan)) return plan;
            TermPairMoments pair = moments.Pair(a, b);
            var mains = pair.Components.Where(c => c.Variables.Length <= maxOrder).ToArray();
            plan = new PairPlan(pair.Covariance, Math.Abs(pair.Cross) + Math.Abs(pair.Product),
                                [.. mains.Select(c => index[c.Variables])], [.. mains.Select(c => c.Value)],
                                [.. pair.Totals.Select(t => t.Variable)], [.. pair.Totals.Select(t => t.Value)]);
            plans[(a, b)] = plan;
            return plan;
        }

        // Per joint draw, each component pair k ≤ k′ (packed): Cov(wₖ, wₖ′) on every subset, each input's total, and the
        // whole; the moment magnitude behind them; and each component's mean E[wₖ] with its magnitude.
        var state = new BasisState(n, r, subsets.Length, measure.Count);
        var draws = new BassDraw[r];
        for (int d = 0; d < n; d++)
        {
            for (int k = 0; k < r; k++)
            {
                draws[k] = pooled[k][joint.Draw(k, d)];
                double[] a = draws[k].Coefficients;
                int[] ids = draws[k].Structure.TermIds;
                double mean = a[0], magnitude = Math.Abs(a[0]);
                for (int j = 0; j < ids.Length; j++)
                {
                    double e = a[1 + j] * moments.Mean(fit.Components[k].Terms[ids[j]]);
                    mean += e;
                    magnitude += Math.Abs(e);
                }
                state.Means[d, k] = mean;
                state.MeanSize[d, k] = magnitude;
            }
            int pk = 0;
            for (int k = 0; k < r; k++)
                for (int k2 = k; k2 < r; k2++, pk++)
                {
                    double[] a = draws[k].Coefficients, b = draws[k2].Coefficients;
                    int[] ia = draws[k].Structure.TermIds, ib = draws[k2].Structure.TermIds;
                    BassTermTable ta = fit.Components[k].Terms, tb = fit.Components[k2].Terms;
                    for (int j = 0; j < ia.Length; j++)
                        for (int l = k == k2 ? j : 0; l < ib.Length; l++)
                        {
                            // Within a component each unordered term pair once, doubled off the diagonal (the direct
                            // fit's sum); across components every ordered pair.
                            double w = a[1 + j] * b[1 + l] * (k == k2 && j != l ? 2.0 : 1.0);
                            if (w == 0.0) continue;
                            PairPlan plan = Plan(ta[ia[j]], tb[ib[l]]);
                            state.Covariance[d, pk] += w * plan.Covariance;
                            state.Size[d, pk] += Math.Abs(w) * plan.Magnitude;
                            for (int q = 0; q < plan.Mains.Length; q++) state.Main[d, plan.Mains[q], pk] += w * plan.MainValues[q];
                            for (int q = 0; q < plan.Totals.Length; q++) state.Total[d, plan.Totals[q], pk] += w * plan.TotalValues[q];
                        }
                }
        }

        var chains = new int[n];
        var within = new int[n];
        if (joint.ChainLengths is { } lengths)
        {
            int at = 0;
            for (int c = 0; c < lengths.Count; c++)
                for (int i = 0; i < lengths[c]; i++, at++)
                {
                    chains[at] = c;
                    within[at] = i;
                }
        }
        else
            for (int d = 0; d < n; d++) within[d] = d;
        return new BasisSobolAnalysis(fit.Decomposition.Map, measure, maxOrder, integratedOrder, subsets, state, chains, within,
                                      spec.RoundingTolerance, moments.LargestError);
    }

    // A term pair's cached contributions: its covariance, its moment magnitude, and its shares by subset index and input.
    private sealed record PairPlan(double Covariance, double Magnitude, int[] Mains, double[] MainValues, int[] Totals, double[] TotalValues);
}

/// <summary>Each joint draw's component covariances, packed over the pairs k ≤ k′, in the coefficient responses' units.</summary>
internal sealed class BasisState
{
    public BasisState(int draws, int components, int subsets, int inputs)
    {
        int pairs = components * (components + 1) / 2;
        Components = components;
        Main = new double[draws, subsets, pairs];
        Total = new double[draws, inputs, pairs];
        Covariance = new double[draws, pairs];
        Size = new double[draws, pairs];
        Means = new double[draws, components];
        MeanSize = new double[draws, components];
    }

    public int Components { get; }

    /// <summary>Draws × subsets × pairs: Cov(wₖ,u, wₖ′,u).</summary>
    public double[,,] Main { get; }

    /// <summary>Draws × inputs × pairs: E[Cov(wₖ, wₖ′ | every other input)].</summary>
    public double[,,] Total { get; }

    /// <summary>Draws × pairs: Cov(wₖ, wₖ′).</summary>
    public double[,] Covariance { get; }

    /// <summary>Draws × pairs: Σ|aₖⱼaₖ′ₗ|(|E[BB′]| + |E B·E B′|), the magnitude the covariance's rounding scales with.</summary>
    public double[,] Size { get; }

    /// <summary>Draws × components: E[wₖ].</summary>
    public double[,] Means { get; }

    /// <summary>Draws × components: |aₖ₀| + Σⱼ|aₖⱼ·E Bₖⱼ|, the magnitude of E[wₖ]'s sum.</summary>
    public double[,] MeanSize { get; }

    /// <summary>Pair (k, k′)'s packed index, either order.</summary>
    public int Pair(int k, int k2)
    {
        if (k > k2) (k, k2) = (k2, k);
        return k * Components - k * (k - 1) / 2 + (k2 - k);
    }
}

/// <summary>
/// A basis fit's sensitivity analysis (<see cref="BassBasisSensitivity.Of"/>, D08): each joint draw's component covariances
/// on every subset, the ℓ-free state, and two views of it as <see cref="SobolDecomposition"/>s in output units:
/// <list type="bullet">
/// <item><see cref="AtLocation"/>, R's form: the decomposition of the output at one location, its variance
/// Σₖₖ′ g_ℓk g_ℓk′·Cov(wₖ, wₖ′) with g_ℓk = y_s(ℓ)·K_ℓk, its mains and totals alike.</item>
/// <item><see cref="Integrated"/>: the location as one more input, the training locations equally weighted. A subset u of
/// the inputs has main variance ḡᵀC_uḡ (the location-averaged surface's), u with the location tr(Σ_G·C_u), with Σ_G the
/// locations' covariance of g, and the location alone the variance over ℓ of the location means y_m(ℓ) + g_ℓᵀE[w]. An
/// input's total is the locations' mean of its totals; the law of total variance links the forms.</item>
/// </list>
/// Zero total variance, at a location or integrated, is <see cref="DecompositionStatus.ZeroVariance"/> with NaN indices,
/// and negative contributions are clamped or flagged under the spec's tolerance, as for direct fits. In the integrated form
/// the location means' part, a sum of squared deviations, counts as zero when its standard deviation is at most the
/// tolerance times the means' magnitude. The joint draws weigh equally.
/// </summary>
public sealed class BasisSobolAnalysis
{
    /// <summary>The integrated form's name for the location input.</summary>
    public const string LocationInput = "location";

    private readonly LinearReconstruction _map;
    private readonly BasisState _state;
    private readonly int[][] _subsets;
    private readonly int[] _chains;
    private readonly int[] _draws;
    private readonly double _tolerance;

    internal BasisSobolAnalysis(LinearReconstruction map, ProductMeasure measure, int maxOrder, int integratedOrder, int[][] subsets,
                                BasisState state, int[] chains, int[] draws, double tolerance, double momentError)
    {
        _map = map;
        _state = state;
        _subsets = subsets;
        _chains = chains;
        _draws = draws;
        _tolerance = tolerance;
        Measure = measure;
        MaxOrder = maxOrder;
        IntegratedMaxOrder = integratedOrder;
        Subsets = Array.AsReadOnly(subsets);
        MomentError = momentError;
        Integrated = BuildIntegrated();
    }

    /// <summary>The product measure, as declared with defaults filled in and resolved.</summary>
    public ProductMeasure Measure { get; }

    /// <summary>The order a location's main indices are reported to.</summary>
    public int MaxOrder { get; }

    /// <summary>The order of <see cref="Integrated"/>'s main indices, the location counted as a variable.</summary>
    public int IntegratedMaxOrder { get; }

    /// <summary>A location's reported subsets of the inputs, as <see cref="SobolDecomposition.Subsets"/>.</summary>
    public IReadOnlyList<int[]> Subsets { get; }

    /// <summary>The output locations m.</summary>
    public int Locations => _map.Locations;

    /// <summary>The components r.</summary>
    public int Components => _map.Coefficients;

    /// <summary>The joint draws.</summary>
    public int Count => _state.Covariance.GetLength(0);

    /// <summary>The largest quadrature error estimate among the one-variable moments (0 unless a Student component was
    /// integrated).</summary>
    public double MomentError { get; }

    /// <summary>The integrated form: inputs the predictors then <see cref="LocationInput"/>.</summary>
    public SobolDecomposition Integrated { get; }

    /// <summary>Location <paramref name="location"/>'s decomposition, computed on each call.</summary>
    public SobolDecomposition AtLocation(int location)
    {
        if ((uint)location >= (uint)Locations) throw new ArgumentOutOfRangeException(nameof(location), $"The basis has {Locations} locations.");
        int r = Components;
        var g = new double[r];
        for (int k = 0; k < r; k++) g[k] = _map.Scales[location] * _map.Basis(location, k);
        double[] weight = Weights((k, k2) => g[k] * g[k2]);
        double[] magnitude = Weights((k, k2) => Math.Abs(g[k] * g[k2]));
        int n = Count, inputs = Measure.Count;
        var main = new double[n, _subsets.Length];
        var total = new double[n, inputs];
        var variance = new double[n];
        var status = new DecompositionStatus[n];
        for (int d = 0; d < n; d++)
        {
            variance[d] = Dot(weight, _state.Covariance, d);
            for (int u = 0; u < _subsets.Length; u++) main[d, u] = Dot(weight, _state.Main, d, u);
            for (int v = 0; v < inputs; v++) total[d, v] = Dot(weight, _state.Total, d, v);
            status[d] = BassSensitivity.Status(variance[d], Dot(magnitude, _state.Size, d), _tolerance, main, total, d);
        }
        return new SobolDecomposition(Measure.Names, _subsets, main, total, variance, status, _chains, _draws, null);
    }

    /// <summary>Every location's decomposition, in location order.</summary>
    public IReadOnlyList<SobolDecomposition> PerLocation() => [.. Enumerable.Range(0, Locations).Select(AtLocation)];

    /// <summary>Draw <paramref name="draw"/>'s Cov(wₖ, wₖ′), r × r, in the coefficient responses' units.</summary>
    public double[,] ComponentCovariance(int draw) => Matrix(p => _state.Covariance[draw, p]);

    /// <summary>Draw <paramref name="draw"/>'s Cov(wₖ,u, wₖ′,u) on <see cref="Subsets"/>[<paramref name="subset"/>], r × r.</summary>
    public double[,] ComponentMainCovariance(int draw, int subset) => Matrix(p => _state.Main[draw, subset, p]);

    /// <summary>Draw <paramref name="draw"/>'s total-effect covariances for input <paramref name="input"/>, r × r.</summary>
    public double[,] ComponentTotalCovariance(int draw, int input) => Matrix(p => _state.Total[draw, input, p]);

    /// <summary>Draw <paramref name="draw"/>'s E[wₖ] under the measure.</summary>
    public double[] ComponentMeans(int draw) => [.. Enumerable.Range(0, Components).Select(k => _state.Means[draw, k])];

    // The integrated form: the inputs' subsets, the location alone, and the inputs' subsets with the location, up to the
    // integrated order, smaller first, then lexicographic with the location last.
    private SobolDecomposition BuildIntegrated()
    {
        int r = Components, m = Locations, p = Measure.Count, n = Count;
        var g = new double[m, r];
        var mean = new double[r];
        for (int l = 0; l < m; l++)
            for (int k = 0; k < r; k++)
            {
                g[l, k] = _map.Scales[l] * _map.Basis(l, k);
                mean[k] += g[l, k] / m;
            }
        double Over(Func<int, double> f)
        {
            double s = 0.0;
            for (int l = 0; l < m; l++) s += f(l);
            return s / m;
        }
        double[] averaged = Weights((k, k2) => mean[k] * mean[k2]);
        double[] spread = Weights((k, k2) => Over(l => (g[l, k] - mean[k]) * (g[l, k2] - mean[k2])));
        double[] second = Weights((k, k2) => Over(l => g[l, k] * g[l, k2]));
        double[] magnitude = Weights((k, k2) => Over(l => Math.Abs(g[l, k] * g[l, k2])));
        double offset = Over(l => _map.Offsets[l]);

        var kinds = new List<(int[] Subset, int Source, bool WithLocation)>();
        for (int u = 0; u < _subsets.Length; u++)
        {
            if (_subsets[u].Length <= IntegratedMaxOrder) kinds.Add((_subsets[u], u, false));
            if (_subsets[u].Length + 1 <= IntegratedMaxOrder) kinds.Add((_subsets[u].Append(p).ToArray(), u, true));
        }
        kinds.Add((new[] { p }, -1, true));
        kinds.Sort((x, y) => Compare(x.Subset, y.Subset));

        var main = new double[n, kinds.Count];
        var total = new double[n, p + 1];
        var variance = new double[n];
        var status = new DecompositionStatus[n];
        for (int d = 0; d < n; d++)
        {
            // The location means' variance, from their deviations, and their magnitude.
            double location = 0.0, size2 = 0.0;
            for (int l = 0; l < m; l++)
            {
                double deviation = _map.Offsets[l] - offset, size = Math.Abs(_map.Offsets[l]);
                for (int k = 0; k < r; k++)
                {
                    deviation += (g[l, k] - mean[k]) * _state.Means[d, k];
                    size += Math.Abs(g[l, k]) * _state.MeanSize[d, k];
                }
                location += deviation * deviation / m;
                size2 += size * size / m;
            }
            for (int s = 0; s < kinds.Count; s++)
            {
                (_, int u, bool withLocation) = kinds[s];
                main[d, s] = u < 0 ? location : Dot(withLocation ? spread : averaged, _state.Main, d, u);
            }
            for (int v = 0; v < p; v++) total[d, v] = Dot(second, _state.Total, d, v);
            double across = Dot(spread, _state.Covariance, d);
            total[d, p] = location + across;
            variance[d] = Dot(averaged, _state.Covariance, d) + across + location;
            status[d] = BassSensitivity.Status(variance[d], Dot(magnitude, _state.Size, d) + _tolerance * size2, _tolerance, main, total, d);
        }
        string[] names = [.. Measure.Names, LocationInput];
        return new SobolDecomposition(names, [.. kinds.Select(x => x.Subset)], main, total, variance, status, _chains, _draws, null);
    }

    // Per packed pair k ≤ k′: f(k, k′), doubled off the diagonal, so a quadratic form is one sum over the pairs.
    private double[] Weights(Func<int, int, double> f)
    {
        int r = Components;
        var w = new double[r * (r + 1) / 2];
        for (int k = 0; k < r; k++)
            for (int k2 = k; k2 < r; k2++) w[_state.Pair(k, k2)] = f(k, k2) * (k == k2 ? 1.0 : 2.0);
        return w;
    }

    private static double Dot(double[] w, double[,] values, int d)
    {
        double s = 0.0;
        for (int q = 0; q < w.Length; q++) s += w[q] * values[d, q];
        return s;
    }

    private static double Dot(double[] w, double[,,] values, int d, int i)
    {
        double s = 0.0;
        for (int q = 0; q < w.Length; q++) s += w[q] * values[d, i, q];
        return s;
    }

    private double[,] Matrix(Func<int, double> packed)
    {
        int r = Components;
        var c = new double[r, r];
        for (int k = 0; k < r; k++)
            for (int k2 = 0; k2 < r; k2++) c[k, k2] = packed(_state.Pair(k, k2));
        return c;
    }

    private static int Compare(int[] x, int[] y)
    {
        if (x.Length != y.Length) return x.Length.CompareTo(y.Length);
        for (int k = 0; k < x.Length; k++)
            if (x[k] != y[k]) return x[k].CompareTo(y[k]);
        return 0;
    }
}
