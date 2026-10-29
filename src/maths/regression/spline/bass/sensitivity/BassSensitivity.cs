using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// Direct sensitivity of a BASS fit (D08; O6, O8): the variance of the latent surface under a declared product measure over
/// the inputs, decomposed for each posterior draw into Sobol main and total effects (Francom et al. 2018; JSS §3), and the
/// main and two-way effect functions. Each draw's surface f = a₀ + Σⱼ aⱼBⱼ is a sum of products of one-variable factors,
/// so its decomposition is exact from <see cref="BassMoments"/>' term-pair ANOVA structure: the variance is Σⱼₗ aⱼaₗ·Cov(Bⱼ,
/// Bₗ), and a subset's main variance and an input's total sum the pairs' shares on it. Observation noise is not in the
/// variance.
/// <para>A direct functional fit's output variables are integrated as more inputs (one decomposition, main and total
/// indices over inputs and outputs), or held at declared coordinates, where the surface's coefficients become aⱼ·hⱼ(t) and
/// each coordinate has its own decomposition.</para>
/// </summary>
public static class BassSensitivity
{
    /// <summary>
    /// The Sobol decomposition of <paramref name="fit"/>'s draws under <paramref name="spec"/>, on <paramref name="data"/>,
    /// which must be the data the fit was run on: the terms' knots and levels index its inputs, and the empirical measure
    /// reads its values.
    /// </summary>
    public static SobolAnalysis Of(BassFit fit, BassData data, SensitivitySpec? spec = null)
    {
        spec ??= new SensitivitySpec();
        Setup s = Setup.Of(fit, data, spec);
        int maxOrder = spec.MaxOrder ?? Math.Max(1, s.TermIds.Select(id => s.Moments.Variables(fit.Terms[id]).Count).DefaultIfEmpty(0).Max());
        int[][] subsets = Subsets(s.TermIds.Select(id => s.Moments.Variables(fit.Terms[id])), maxOrder, spec.SubsetBudget);
        var index = new Dictionary<int[], int>(subsets.Length, SubsetComparer.Instance);
        for (int k = 0; k < subsets.Length; k++) index[subsets[k]] = k;

        var plans = new Dictionary<(int, int), PairPlan>();
        PairPlan Plan(int a, int b)
        {
            if (a > b) (a, b) = (b, a);
            if (plans.TryGetValue((a, b), out PairPlan? plan)) return plan;
            TermPairMoments pair = s.Moments.Pair(fit.Terms[a], fit.Terms[b]);
            var mains = pair.Components.Where(c => c.Variables.Length <= maxOrder).ToArray();
            plan = new PairPlan(pair.Covariance, Math.Abs(pair.Cross) + Math.Abs(pair.Product),
                                [.. mains.Select(c => index[c.Variables])], [.. mains.Select(c => c.Value)],
                                [.. pair.Totals.Select(t => t.Variable)], [.. pair.Totals.Select(t => t.Value)]);
            plans[(a, b)] = plan;
            return plan;
        }

        int inputs = s.Moments.Measure.Count, n = s.Draws.Length;
        var decompositions = new SobolDecomposition[s.Surfaces];
        double[] b = new double[s.MaxTerms + 1];
        for (int surface = 0; surface < s.Surfaces; surface++)
        {
            var main = new double[n, subsets.Length];
            var total = new double[n, inputs];
            var variance = new double[n];
            var status = new DecompositionStatus[n];
            for (int d = 0; d < n; d++)
            {
                BassDraw draw = s.Draws[d];
                int[] ids = draw.Structure.TermIds;
                s.Coefficients(surface, draw, b);
                double v = 0.0, size = 0.0;
                for (int j = 0; j < ids.Length; j++)
                    for (int l = j; l < ids.Length; l++)
                    {
                        double w = b[1 + j] * b[1 + l] * (j == l ? 1.0 : 2.0);
                        if (w == 0.0) continue;
                        PairPlan plan = Plan(ids[j], ids[l]);
                        v += w * plan.Covariance;
                        size += Math.Abs(w) * plan.Magnitude;
                        for (int k = 0; k < plan.Mains.Length; k++) main[d, plan.Mains[k]] += w * plan.MainValues[k];
                        for (int k = 0; k < plan.Totals.Length; k++) total[d, plan.Totals[k]] += w * plan.TotalValues[k];
                    }
                variance[d] = v;
                status[d] = Status(v, size, spec.RoundingTolerance, main, total, d);
            }
            decompositions[surface] = new SobolDecomposition(s.Moments.Measure.Names, subsets, main, total, variance, status,
                                                             s.Chains, s.DrawIndex, s.Weights);
        }
        return new SobolAnalysis(s.Form, s.Moments.Measure, maxOrder, decompositions, s.Coordinates, s.CoordinatesOutside,
                                 s.Moments.LargestError);
    }

    /// <summary>
    /// The main and two-way effect functions of <paramref name="fit"/>'s draws under <paramref name="spec"/>'s measure and the
    /// fit's degree, on <paramref name="data"/>, the data the fit was run on. Refused, with its reason, for an effect of any
    /// other size, an input no variable has, a grid under two points, or more values than <see cref="EffectSpec.CellBudget"/>.
    /// </summary>
    public static EffectAnalysis Effects(BassFit fit, BassData data, SensitivitySpec? spec = null, EffectSpec? effects = null)
    {
        spec ??= new SensitivitySpec();
        effects ??= new EffectSpec();
        if (effects.GridPoints < 2) throw new ArgumentOutOfRangeException(nameof(effects), "An effect grid needs two points.");
        Setup s = Setup.Of(fit, data, spec);
        ProductMeasure measure = s.Moments.Measure;

        int[][] chosen;
        if (effects.Subsets is { } named)
        {
            chosen = new int[named.Count][];
            for (int e = 0; e < named.Count; e++)
            {
                string[] names = named[e] ?? throw new ArgumentException($"Effect {e} is null.", nameof(effects));
                if (names.Length is not (1 or 2)) throw new ArgumentException($"Effects are main or two-way; effect {e} names {names.Length} inputs.", nameof(effects));
                int[] vars = [.. names.Select(name => IndexOf(measure, name))];
                if (vars.Length == 2 && vars[0] == vars[1]) throw new ArgumentException($"Effect {e} names '{names[0]}' twice.", nameof(effects));
                Array.Sort(vars);
                chosen[e] = vars;
            }
        }
        else
        {
            var pairs = new SortedSet<(int, int)>();
            foreach (int id in s.TermIds)
            {
                IReadOnlyList<int> v = s.Moments.Variables(fit.Terms[id]);
                for (int i = 0; i < v.Count; i++)
                    for (int k = i + 1; k < v.Count; k++) pairs.Add((v[i], v[k]));
            }
            chosen = [.. Enumerable.Range(0, measure.Count).Select(v => new[] { v }), .. pairs.Select(p => new[] { p.Item1, p.Item2 })];
        }

        // Each variable's grid: scaled points (or level indices) and their presentation.
        var scaled = new double[measure.Count][];
        var grids = new EffectGrid[measure.Count];
        for (int v = 0; v < measure.Count; v++)
        {
            if (measure.Continuous(v) is { } c)
            {
                scaled[v] = [.. Enumerable.Range(0, effects.GridPoints).Select(i => c.Lower + (c.Upper - c.Lower) * i / (effects.GridPoints - 1))];
                RangeScale range = s.RangeOf(v);
                grids[v] = new EffectGrid([.. scaled[v].Select(x => range.Min + x * (range.Max - range.Min))], null);
            }
            else
            {
                var levels = ((TrainingLevels)fit.Transform.Predictors[v]).Levels;
                scaled[v] = [.. Enumerable.Range(0, levels.Count).Select(l => (double)l)];
                grids[v] = new EffectGrid([.. scaled[v]], levels);
            }
        }
        long cells = chosen.Sum(u => u.Aggregate(1L, (acc, v) => acc * scaled[v].Length));
        long values = cells * s.Draws.Length * s.Surfaces;
        if (values > effects.CellBudget)
            throw new ArgumentException($"The effects hold {values} values ({cells} cells × {s.Draws.Length} draws × {s.Surfaces} surfaces), above the budget of {effects.CellBudget}: name fewer effects, a coarser grid, or raise the budget.", nameof(effects));

        // Per term: its centred factor values on each of its variables' grids, and κ·Π E φ over the variables outside an effect.
        var centred = new Dictionary<int, double[][]>();
        foreach (int id in s.TermIds)
        {
            BassTerm term = fit.Terms[id];
            IReadOnlyList<int> v = s.Moments.Variables(term);
            IReadOnlyList<double> means = s.Moments.FactorMeans(term);
            var rows = new double[v.Count][];
            for (int k = 0; k < v.Count; k++)
                rows[k] = [.. scaled[v[k]].Select(x => s.Moments.FactorValue(term, k, x) - means[k])];
            centred[id] = rows;
        }

        var at = new EffectFunctions[s.Surfaces];
        double[] b = new double[s.MaxTerms + 1];
        int n = s.Draws.Length;
        for (int surface = 0; surface < s.Surfaces; surface++)
        {
            var mean = new double[n];
            var tables = chosen.Select(u => new double[n, u.Aggregate(1, (acc, v) => acc * scaled[v].Length)]).ToArray();
            for (int d = 0; d < n; d++)
            {
                BassDraw draw = s.Draws[d];
                int[] ids = draw.Structure.TermIds;
                s.Coefficients(surface, draw, b);
                mean[d] = b[0];
                for (int j = 0; j < ids.Length; j++)
                {
                    if (b[1 + j] == 0.0) continue;
                    BassTerm term = fit.Terms[ids[j]];
                    IReadOnlyList<int> v = s.Moments.Variables(term);
                    IReadOnlyList<double> means = s.Moments.FactorMeans(term);
                    double kappa = s.Moments.Scale(term);
                    mean[d] += b[1 + j] * s.Moments.Mean(term);
                    for (int e = 0; e < chosen.Length; e++)
                    {
                        int[] u = chosen[e];
                        int p0 = Position(v, u[0]), p1 = u.Length == 2 ? Position(v, u[1]) : -1;
                        if (p0 < 0 || (u.Length == 2 && p1 < 0)) continue;
                        double rho = b[1 + j] * kappa;
                        for (int k = 0; k < v.Count; k++)
                            if (k != p0 && k != p1) rho *= means[k];
                        if (rho == 0.0) continue;
                        double[] g0 = centred[ids[j]][p0];
                        double[,] table = tables[e];
                        if (p1 < 0)
                            for (int c = 0; c < g0.Length; c++) table[d, c] += rho * g0[c];
                        else
                        {
                            double[] g1 = centred[ids[j]][p1];
                            for (int c0 = 0; c0 < g0.Length; c0++)
                                for (int c1 = 0; c1 < g1.Length; c1++) table[d, c0 * g1.Length + c1] += rho * g0[c0] * g1[c1];
                        }
                    }
                }
            }
            at[surface] = new EffectFunctions(measure.Names,
                [.. chosen.Select((u, e) => new Effect(u, [.. u.Select(v => grids[v])], tables[e]))], mean);
        }
        return new EffectAnalysis(s.Form, measure, at, s.Coordinates);
    }

    private static int Position(IReadOnlyList<int> variables, int v)
    {
        for (int k = 0; k < variables.Count; k++)
            if (variables[k] == v) return k;
        return -1;
    }

    private static int IndexOf(ProductMeasure measure, string name)
    {
        for (int v = 0; v < measure.Count; v++)
            if (measure.Names[v] == name) return v;
        throw new ArgumentException($"No input is named '{name}'.");
    }

    // A draw's status: zero variance within rounding, else its negative contributions clamped within the tolerance.
    internal static DecompositionStatus Status(double variance, double size, double tolerance, double[,] main, double[,] total, int d)
    {
        if (!(variance > tolerance * size)) return DecompositionStatus.ZeroVariance;
        var status = DecompositionStatus.Defined;
        double floor = -tolerance * variance;
        foreach (double[,] table in new[] { main, total })
            for (int k = 0; k < table.GetLength(1); k++)
            {
                if (table[d, k] >= 0.0) continue;
                if (table[d, k] >= floor)
                {
                    table[d, k] = 0.0;
                    status |= DecompositionStatus.Clamped;
                }
                else status |= DecompositionStatus.NegativeBeyondTolerance;
            }
        return status;
    }

    // The nonempty subsets, up to maxOrder, of each posterior term's variables: smaller first, then lexicographic.
    internal static int[][] Subsets(IEnumerable<IReadOnlyList<int>> terms, int maxOrder, int budget)
    {
        if (maxOrder < 1) throw new ArgumentOutOfRangeException(nameof(maxOrder), "The order is a positive integer.");
        if (budget < 1) throw new ArgumentOutOfRangeException(nameof(budget), "The subset budget is a positive integer.");
        var found = new HashSet<int[]>(SubsetComparer.Instance);
        foreach (IReadOnlyList<int> v in terms)
        {
            for (int mask = 1; mask < 1 << v.Count; mask++)
            {
                if (System.Numerics.BitOperations.PopCount((uint)mask) > maxOrder) continue;
                int[] u = [.. Enumerable.Range(0, v.Count).Where(k => (mask >> k & 1) != 0).Select(k => v[k])];
                if (found.Add(u) && found.Count > budget)
                    throw new ArgumentException($"The posterior's terms span more than {budget} subsets up to order {maxOrder}: lower the order or raise the subset budget.");
            }
        }
        var list = found.ToList();
        list.Sort((x, y) =>
        {
            if (x.Length != y.Length) return x.Length.CompareTo(y.Length);
            for (int k = 0; k < x.Length; k++)
                if (x[k] != y[k]) return x[k].CompareTo(y[k]);
            return 0;
        });
        return [.. list];
    }

    // A term pair's cached contributions: its covariance, its moment magnitude, and its shares by subset index and input.
    private sealed record PairPlan(double Covariance, double Magnitude, int[] Mains, double[] MainValues, int[] Totals, double[] TotalValues);

    internal sealed class SubsetComparer : IEqualityComparer<int[]>
    {
        public static readonly SubsetComparer Instance = new();

        public bool Equals(int[]? x, int[]? y) => x is not null && y is not null && x.AsSpan().SequenceEqual(y);

        public int GetHashCode(int[] obj)
        {
            var hash = new HashCode();
            foreach (int v in obj) hash.Add(v);
            return hash.ToHashCode();
        }
    }

    // What both analyses share: the inputs, measure and moment service, the draws in chain order, and each surface's
    // coefficients (one surface, or one per conditioned coordinate).
    private sealed class Setup
    {
        private double[][]? _outputSide;   // conditioned form: per term id, its output side at each coordinate

        public required BassFit Fit { get; init; }
        public required BassMoments Moments { get; init; }
        public required BassDraw[] Draws { get; init; }
        public required int[] Chains { get; init; }
        public required int[] DrawIndex { get; init; }
        public required double[]? Weights { get; init; }
        public required int[] TermIds { get; init; }
        public required int MaxTerms { get; init; }
        public required FunctionalForm? Form { get; init; }
        public required IReadOnlyList<double[]>? Coordinates { get; init; }
        public required IReadOnlyList<bool>? CoordinatesOutside { get; init; }

        public int Surfaces => Coordinates?.Count ?? 1;

        // Variable v's range scale: a continuous predictor's or an integrated output variable's.
        public RangeScale RangeOf(int v)
        {
            int predictors = Fit.Transform.Predictors.Count;
            return v < predictors ? (RangeScale)Fit.Transform.Predictors[v] : Fit.Transform.OutputVariables[v - predictors];
        }

        // b₀ and one coefficient per term of the draw's structure: the draw's own, or aⱼ·hⱼ(t) at a conditioned coordinate,
        // the terms with no case factor folded into b₀.
        public void Coefficients(int surface, BassDraw draw, double[] b)
        {
            double[] a = draw.Coefficients;
            int[] ids = draw.Structure.TermIds;
            if (_outputSide is null)
            {
                Array.Copy(a, b, ids.Length + 1);
                return;
            }
            b[0] = a[0];
            for (int j = 0; j < ids.Length; j++)
            {
                double h = _outputSide[ids[j]][surface];
                if (Fit.Terms[ids[j]].HasCaseFactor) b[1 + j] = a[1 + j] * h;
                else
                {
                    b[0] += a[1 + j] * h;
                    b[1 + j] = 0.0;
                }
            }
        }

        public static Setup Of(BassFit fit, BassData data, SensitivitySpec spec)
        {
            ArgumentNullException.ThrowIfNull(fit);
            ArgumentNullException.ThrowIfNull(data);
            ArgumentNullException.ThrowIfNull(spec);
            string fingerprint = Bass.Fingerprint(fit.Config, data);
            if (fingerprint != fit.Identity.Fingerprint)
                throw new ArgumentException($"The fit was run on other data: its fingerprint is {fit.Identity.Fingerprint}, theirs {fingerprint}.", nameof(data));
            if (!(spec.RoundingTolerance >= 0.0) || double.IsPositiveInfinity(spec.RoundingTolerance))
                throw new ArgumentOutOfRangeException(nameof(spec), "The rounding tolerance must be nonnegative and finite.");
            if (spec.Surface != DecomposedSurface.Mean) throw new ArgumentException($"The decomposed surface {spec.Surface} is not built.", nameof(spec));

            BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
            FunctionalForm? form = inputs.Outputs is null ? null : spec.Functional;
            if (spec.Coordinates is not null && form != FunctionalForm.Condition)
                throw new ArgumentException("Output coordinates apply only to a functional fit's conditioned form.", nameof(spec));
            var measure = ProductMeasure.Resolve(fit.Transform, inputs, spec.Measures, spec.OutputMeasures, form == FunctionalForm.Integrate);
            var moments = BassMoments.On(inputs, fit.Config.Structure.Degree, measure);

            int total = fit.Draws.Sum(c => c.Count);
            var draws = new BassDraw[total];
            var chains = new int[total];
            var index = new int[total];
            int at = 0;
            for (int c = 0; c < fit.Draws.Count; c++)
                for (int i = 0; i < fit.Draws[c].Count; i++, at++)
                {
                    draws[at] = fit.Draws[c][i];
                    chains[at] = c;
                    index[at] = i;
                }
            if (fit.DrawWeights is { } w && w.Length != total)
                throw new ArgumentException($"The fit has {total} draws and {w.Length} weights.", nameof(fit));
            int[] ids = [.. draws.SelectMany(d => d.Structure.TermIds).Distinct().Order()];

            var setup = new Setup
            {
                Fit = fit, Moments = moments, Draws = draws, Chains = chains, DrawIndex = index, Weights = fit.DrawWeights,
                TermIds = ids, MaxTerms = draws.Length == 0 ? 0 : draws.Max(d => d.Structure.TermIds.Length), Form = form,
                Coordinates = null, CoordinatesOutside = null,
            };
            if (form != FunctionalForm.Condition) return setup;

            // The conditioned form's coordinates, in their own units, scaled; each term's output side at them.
            var functional = (FunctionalResponse)data.Response;
            int dims = fit.Transform.OutputVariables.Count;
            IReadOnlyList<double[]> coordinates = spec.Coordinates
                ?? [.. Enumerable.Range(0, functional.Locations).Select(r => Enumerable.Range(0, dims).Select(k => functional.OutputVariables[k].Values.Span[r]).ToArray())];
            var columns = new double[dims][];
            for (int k = 0; k < dims; k++) columns[k] = new double[coordinates.Count];
            var outside = new bool[coordinates.Count];
            for (int r = 0; r < coordinates.Count; r++)
            {
                double[] t = coordinates[r] ?? throw new ArgumentException($"Coordinate {r} is null.", nameof(spec));
                if (t.Length != dims) throw new ArgumentException($"Coordinate {r} has {t.Length} values for {dims} output variables.", nameof(spec));
                for (int k = 0; k < dims; k++)
                {
                    if (!double.IsFinite(t[k])) throw new ArgumentException($"Coordinate {r} is not finite.", nameof(spec));
                    columns[k][r] = fit.Transform.OutputVariables[k].Scale(t[k]);
                    outside[r] |= !fit.Transform.OutputVariables[k].Contains(t[k]);
                }
            }
            var basis = new BassBasis(inputs, fit.Config.Structure.Degree);
            var at_ = new OutputInputs(columns);
            var outputSide = new double[fit.Terms.Count][];
            foreach (int id in ids)
            {
                outputSide[id] = new double[coordinates.Count];
                basis.OutputColumn(fit.Terms[id], at_, outputSide[id]);
            }
            return new Setup
            {
                Fit = fit, Moments = moments, Draws = draws, Chains = chains, DrawIndex = index, Weights = fit.DrawWeights,
                TermIds = ids, MaxTerms = setup.MaxTerms, Form = form,
                Coordinates = [.. coordinates.Select(c => (double[])c.Clone())], CoordinatesOutside = outside,
                _outputSide = outputSide,
            };
        }
    }
}
