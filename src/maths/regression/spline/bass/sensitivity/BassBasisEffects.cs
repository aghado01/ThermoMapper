using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// Effect functions of a basis fit's reconstructed latent output (O8 for basis fits; R's <c>sobolBasis</c> has none). An
/// effect is linear in the surface, so it has no cross-component terms: at location ℓ the effect on u is
/// Σₖ g_ℓk·eₖ,u, with g_ℓk = y_s(ℓ)·K_ℓk and eₖ,u component k's effect (<see cref="BassSensitivity.Effects"/> on its
/// coefficient response), and each draw's mean is y_m(ℓ) + Σₖ g_ℓk·E[wₖ]. The components' effects are read per joint draw
/// and kept, free of ℓ, as basis sensitivity keeps its covariances (<see cref="BasisSobolAnalysis"/>).
/// </summary>
public static class BassBasisEffects
{
    /// <summary>
    /// The effects of <paramref name="fit"/>'s joint draws under <paramref name="spec"/>'s measure and each component's
    /// degree. With no subsets named, every input's main effect and the two-way effect of every pair a term of some joint
    /// draw holds, in any component. Refused, with its reason: output coordinates and output-variable measures (as basis
    /// sensitivity refuses them), more values than <see cref="EffectSpec.CellBudget"/> over the components and joint draws,
    /// and the direct path's own refusals.
    /// </summary>
    public static BasisEffectAnalysis Of(BassBasisFit fit, SensitivitySpec? spec = null, EffectSpec? effects = null)
    {
        ArgumentNullException.ThrowIfNull(fit);
        spec ??= new SensitivitySpec();
        effects ??= new EffectSpec();
        if (spec.Coordinates is not null)
            throw new ArgumentException("A basis exists only at its training locations, which have no output coordinates: read a location's effects by its index.", nameof(spec));
        if (spec.OutputMeasures is { Count: > 0 })
            throw new ArgumentException("A basis fit's locations are integrated equally weighted, its training grid; output-variable measures apply to direct functional fits.", nameof(spec));
        if (effects.GridPoints < 2) throw new ArgumentOutOfRangeException(nameof(effects), "An effect grid needs two points.");

        int r = fit.Components.Count;
        BassJointDraws joint = fit.JointDraws;
        int n = joint.Count;
        BassFit first = fit.Components[0];
        BassInputs inputs = BassInputs.Prepare(fit.ComponentData(0), first.Transform);
        ProductMeasure measure = ProductMeasure.Resolve(first.Transform, inputs, spec.Measures);

        // The effects, by name, so that every component computes the same set.
        string[][] named;
        if (effects.Subsets is { } given) named = [.. given];
        else
        {
            BassMoments moments = BassMoments.On(inputs, first.Config.Structure.Degree, measure);
            var pairs = new SortedSet<(int, int)>();
            for (int k = 0; k < r; k++)
            {
                List<BassDraw> pooled = BassDrawFit.Pool(fit.Components[k].Draws, fit.Components[k].Terms, false, nameof(fit)).Pooled;
                var seen = new HashSet<int>();
                for (int d = 0; d < n; d++)
                    foreach (int id in pooled[joint.Draw(k, d)].Structure.TermIds)
                    {
                        if (!seen.Add(id)) continue;
                        IReadOnlyList<int> v = moments.Variables(fit.Components[k].Terms[id]);
                        for (int i = 0; i < v.Count; i++)
                            for (int j = i + 1; j < v.Count; j++) pairs.Add((v[i], v[j]));
                    }
            }
            named = [.. Enumerable.Range(0, measure.Count).Select(v => new[] { measure.Names[v] }),
                     .. pairs.Select(q => new[] { measure.Names[q.Item1], measure.Names[q.Item2] })];
        }

        // The cells each effect holds, from its inputs' grids, checked against the budget before any component runs.
        int Size(string name)
        {
            for (int v = 0; v < measure.Count; v++)
                if (measure.Names[v] == name)
                    return measure.Continuous(v) is null ? ((TrainingLevels)first.Transform.Predictors[v]).Levels.Count : effects.GridPoints;
            throw new ArgumentException($"No input is named '{name}'.", nameof(effects));
        }
        long cells = 0;
        for (int e = 0; e < named.Length; e++)
        {
            string[] names = named[e] ?? throw new ArgumentException($"Effect {e} is null.", nameof(effects));
            cells += names.Aggregate(1L, (acc, name) => acc * Size(name));
        }
        long values = cells * n * r;
        if (values > effects.CellBudget)
            throw new ArgumentException($"The effects hold {values} values ({cells} cells × {n} joint draws × {r} components), above the budget of {effects.CellBudget}: name fewer effects, a coarser grid, or raise the budget.", nameof(effects));

        // Each component's effects over its pooled draws, read at the joint draws.
        EffectSpec each = effects with { Subsets = named };
        var state = new double[r][][,];
        var means = new double[n, r];
        IReadOnlyList<Effect>? layout = null;
        for (int k = 0; k < r; k++)
        {
            EffectFunctions own = BassSensitivity.Effects(fit.Components[k], fit.ComponentData(k), spec, each).At[0];
            layout ??= own.Effects;
            state[k] = new double[own.Effects.Count][,];
            for (int e = 0; e < own.Effects.Count; e++)
            {
                double[,] source = own.Effects[e].Values;
                int width = source.GetLength(1);
                var target = new double[n, width];
                for (int d = 0; d < n; d++)
                {
                    int from = joint.Draw(k, d);
                    for (int c = 0; c < width; c++) target[d, c] = source[from, c];
                }
                state[k][e] = target;
            }
            for (int d = 0; d < n; d++) means[d, k] = own.Mean[joint.Draw(k, d)];
        }
        return new BasisEffectAnalysis(fit.Decomposition.Map, measure, [.. layout!.Select(x => (x.Inputs, x.Grids))], state, means, effects.CellBudget);
    }
}

/// <summary>
/// A basis fit's effect functions (<see cref="BassBasisEffects.Of"/>): each component's effects per joint draw, in its
/// coefficient response's units and free of ℓ, and views of them as <see cref="EffectFunctions"/> in output units:
/// <list type="bullet">
/// <item><see cref="AtLocation"/>: the output at one location, effects Σₖ g_ℓk·eₖ,u and mean y_m(ℓ) + Σₖ g_ℓk·E[wₖ].</item>
/// <item><see cref="Integrated"/>: the location as one more input, the training locations equally weighted, as basis
/// sensitivity's integrated form. An effect on inputs u is weighted by ḡ, the locations' mean of g; the location's main
/// effect is y_m(ℓ) + g_ℓᵀE[w] less its mean over ℓ, y_m included; and each one-input effect has a two-way effect with the
/// location, Σₖ (g_ℓk − ḡₖ)·eₖ,v.</item>
/// </list>
/// The joint draws weigh equally.
/// </summary>
public sealed class BasisEffectAnalysis
{
    private readonly LinearReconstruction _map;
    private readonly (int[] Inputs, IReadOnlyList<EffectGrid> Grids)[] _layout;
    private readonly double[][][,] _state;   // [component][effect] → joint draws × cells
    private readonly double[,] _means;       // joint draws × components: E[wₖ]
    private readonly long _budget;

    internal BasisEffectAnalysis(LinearReconstruction map, ProductMeasure measure, (int[] Inputs, IReadOnlyList<EffectGrid> Grids)[] layout,
                                 double[][][,] state, double[,] means, long budget)
    {
        _map = map;
        _layout = layout;
        _state = state;
        _means = means;
        _budget = budget;
        Measure = measure;
    }

    /// <summary>The integrated form's name for the location input.</summary>
    public const string LocationInput = BasisSobolAnalysis.LocationInput;

    /// <summary>The product measure, as declared with defaults filled in and resolved.</summary>
    public ProductMeasure Measure { get; }

    /// <summary>The output locations m.</summary>
    public int Locations => _map.Locations;

    /// <summary>The components r.</summary>
    public int Components => _map.Coefficients;

    /// <summary>The joint draws.</summary>
    public int Count => _means.GetLength(0);

    /// <summary>Component <paramref name="component"/>'s effects and means per joint draw, in its coefficient response's
    /// units.</summary>
    public EffectFunctions ComponentEffects(int component)
    {
        if ((uint)component >= (uint)Components) throw new ArgumentOutOfRangeException(nameof(component), $"The basis has {Components} components.");
        var g = new double[Components];
        g[component] = 1.0;
        return Combine(g, 0.0);
    }

    /// <summary>Location <paramref name="location"/>'s effects in output units, computed on each call.</summary>
    public EffectFunctions AtLocation(int location)
    {
        if ((uint)location >= (uint)Locations) throw new ArgumentOutOfRangeException(nameof(location), $"The basis has {Locations} locations.");
        return Combine(G(location), _map.Offsets[location]);
    }

    /// <summary>Every location's effects, in location order.</summary>
    public IReadOnlyList<EffectFunctions> PerLocation() => [.. Enumerable.Range(0, Locations).Select(AtLocation)];

    /// <summary>
    /// The integrated form's effects, inputs the predictors then <see cref="LocationInput"/>, computed on each call: the
    /// chosen effects weighted by ḡ, in their order, then the location's main effect, then each one-input effect with the
    /// location, the location's index fastest. Refused when they hold more values than the analysis's cell budget.
    /// </summary>
    public EffectFunctions Integrated()
    {
        int r = Components, m = Locations, n = Count, p = Measure.Count;
        var mean = new double[r];
        for (int l = 0; l < m; l++)
            for (int k = 0; k < r; k++) mean[k] += _map.Scales[l] * _map.Basis(l, k) / m;
        double offset = 0.0;
        for (int l = 0; l < m; l++) offset += _map.Offsets[l] / m;

        long cells = m;
        for (int e = 0; e < _layout.Length; e++)
            cells += Width(e) * (_layout[e].Inputs.Length == 1 ? 1L + m : 1L);
        if (cells * n > _budget)
            throw new ArgumentException($"The integrated effects hold {cells * n} values ({cells} cells × {n} joint draws), above the budget of {_budget}.");

        EffectFunctions averaged = Combine(mean, offset);
        var effects = new List<Effect>(averaged.Effects);
        var location = new EffectGrid([.. Enumerable.Range(0, m).Select(l => (double)l)], null);
        var main = new double[n, m];
        for (int l = 0; l < m; l++)
        {
            double[] g = G(l);
            for (int d = 0; d < n; d++)
            {
                double f = _map.Offsets[l];
                for (int k = 0; k < r; k++) f += g[k] * _means[d, k];
                main[d, l] = f - averaged.Mean[d];
            }
        }
        effects.Add(new Effect([p], [location], main));
        for (int e = 0; e < _layout.Length; e++)
        {
            if (_layout[e].Inputs.Length != 1) continue;
            int width = Width(e);
            var values = new double[n, width * m];
            for (int l = 0; l < m; l++)
            {
                double[] g = G(l);
                for (int k = 0; k < r; k++)
                {
                    double weight = g[k] - mean[k];
                    if (weight == 0.0) continue;
                    double[,] source = _state[k][e];
                    for (int d = 0; d < n; d++)
                        for (int c = 0; c < width; c++) values[d, c * m + l] += weight * source[d, c];
                }
            }
            effects.Add(new Effect([_layout[e].Inputs[0], p], [_layout[e].Grids[0], location], values));
        }
        return new EffectFunctions([.. Measure.Names, LocationInput], effects, averaged.Mean);
    }

    // g_ℓ = y_s(ℓ)·K_ℓ.
    private double[] G(int location)
    {
        var g = new double[Components];
        for (int k = 0; k < g.Length; k++) g[k] = _map.Scales[location] * _map.Basis(location, k);
        return g;
    }

    private int Width(int effect) => _state[0][effect].GetLength(1);

    // The effects Σₖ gₖ·eₖ,u and means offset + Σₖ gₖ·E[wₖ], the components summed in order.
    private EffectFunctions Combine(double[] g, double offset)
    {
        int n = Count, r = Components;
        var effects = new Effect[_layout.Length];
        for (int e = 0; e < _layout.Length; e++)
        {
            int width = Width(e);
            var values = new double[n, width];
            for (int k = 0; k < r; k++)
            {
                if (g[k] == 0.0) continue;
                double[,] source = _state[k][e];
                for (int d = 0; d < n; d++)
                    for (int c = 0; c < width; c++) values[d, c] += g[k] * source[d, c];
            }
            effects[e] = new Effect(_layout[e].Inputs, _layout[e].Grids, values);
        }
        var mean = new double[n];
        for (int d = 0; d < n; d++)
        {
            double f = offset;
            for (int k = 0; k < r; k++) f += g[k] * _means[d, k];
            mean[d] = f;
        }
        return new EffectFunctions(Measure.Names, effects, mean);
    }
}
