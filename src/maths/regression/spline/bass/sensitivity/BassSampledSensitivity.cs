using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Sensitivity;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// The sampled route of a BASS fit's Sobol decomposition (the sampled-sensitivity brief, slice 1): pick-freeze
/// (<see cref="PickFreeze"/>) on two base sets of N rows drawn once from the product law (E8's sampler, from the spec's seed) and shared
/// by every draw, so the draws' spread is the posterior's and not the Monte Carlo error's. It accumulates each draw's surface on A, B
/// and each block's hybrid, every distinct term's column read once per batch of draws over the stacked rows, and reduces them to a
/// <see cref="SobolDecomposition"/> whose every piece carries its standard error. Exact values per draw; only the integral over x is
/// approximate.
/// </summary>
public static partial class BassSensitivity
{
    internal const long ColumnBudget = 1L << 23;   // term-column values one batch of draws keeps, about 64 MB

    /// <summary><see cref="Of"/> with the main variances read in <paramref name="form"/> on the sampled route: the gates read both forms
    /// while the form's default waits on the user (the slice 0 report's follow-up).</summary>
    internal static SobolAnalysis OfInForm(BassFit fit, BassData data, SensitivitySpec? spec, MainVarianceForm form) => OfCore(fit, data, spec, null, form);

    /// <summary>A fit of several outputs read on <paramref name="output"/>'s functional, its main variances in <paramref name="form"/>.</summary>
    internal static SobolAnalysis Projected(BassFit fit, BassData data, double[] output, SensitivitySpec? spec, MainVarianceForm form) => OfCore(fit, data, spec, output, form);

    // Decision 1: the exact route where every term the draws hold has one-variable factors, the sampled route otherwise; Exact refuses a
    // term without them up front, with its reason.
    private static SensitivityRoute Resolve(BassFit fit, SensitivitySpec spec, Setup s)
    {
        var kernels = new DesignKernels(s.Moments.Inputs, fit.Config.Structure.Degree, fit.Config.Structure.ResolvedHingeNormalization);
        BassTerm? without = s.TermIds.Select(id => fit.Terms[id]).FirstOrDefault(t => !kernels.Separable(t));
        return spec.Route switch
        {
            SensitivityRoute.Auto => without is null ? SensitivityRoute.Exact : SensitivityRoute.Sampled,
            SensitivityRoute.Exact when without is not null => throw new ArgumentException(
                $"Term {without} is not a product of one-variable factors, so the exact route cannot integrate it; take the sampled route "
                + "(Route = Sampled, or Auto).", nameof(spec)),
            SensitivityRoute.Exact or SensitivityRoute.Sampled => spec.Route,
            _ => throw new ArgumentOutOfRangeException(nameof(spec), $"No sensitivity route {spec.Route}."),
        };
    }

    // A fit of several named outputs (E34) reads its declared output functional's surface: Setup's coefficients project each draw's
    // output coefficients, so the route reads it as it reads a scalar fit.
    private static SobolAnalysis Sampled(BassFit fit, SensitivitySpec spec, Setup s, MainVarianceForm? form)
    {
        if (s.Form is not null) return SampledFunctional(fit, spec, s, form);
        if (spec.Points < 2) throw new ArgumentOutOfRangeException(nameof(spec), "The sampled route needs two rows in each base set.");
        BassInputs inputs = s.Moments.Inputs;
        ProductMeasure measure = s.Moments.Measure;
        int p = measure.Count, n = spec.Points;
        int degree = fit.Config.Structure.Degree;
        HingeNormalization hinges = fit.Config.Structure.ResolvedHingeNormalization;
        var kernels = new DesignKernels(inputs, degree, hinges);

        Func<int, int[]> Variables = TermVariables(fit, s, kernels);
        int maxOrder = spec.MaxOrder ?? Math.Max(1, s.TermIds.Select(id => Variables(id).Length).DefaultIfEmpty(0).Max());
        int[][] subsets = Subsets(s.TermIds.Select(id => (IReadOnlyList<int>)Variables(id)), maxOrder, spec.SubsetBudget);
        PickFreezeDesign design = form is { } chosen ? new PickFreezeDesign(p, subsets, chosen, p + spec.SubsetBudget)
                                                     : new PickFreezeDesign(p, subsets, budget: p + spec.SubsetBudget);
        var read = new bool[p];
        foreach (int id in s.TermIds)
            foreach (int v in Variables(id)) read[v] = true;

        // The stacked rows: A, B, then each block reading only inputs some term reads, hybrid then reverse; a block holding an input no
        // term reads is a singleton whose values are A's (its reverse B's), bit for bit, so it is not evaluated.
        var segments = new List<(bool FromB, int[] Other)> { (false, []), (true, []) };
        var hybridAt = new int[design.Blocks.Count];
        var reverseAt = new int[design.Blocks.Count];
        for (int k = 0; k < design.Blocks.Count; k++)
        {
            int[] u = design.Blocks[k];
            if (!u.All(v => read[v]))
            {
                hybridAt[k] = 0;
                reverseAt[k] = 1;
                continue;
            }
            hybridAt[k] = segments.Count;
            segments.Add((false, u));
            if (!design.Reverses) continue;
            reverseAt[k] = segments.Count;
            segments.Add((true, u));
        }
        PredictionCases drawn = BassActiveSubspaces.Points(fit.Transform, measure, 2 * n, spec.Seed);
        (CaseInputs cases, _) = Stacked(drawn, measure, n, segments).Prepare(fit.Transform, inputs, UnseenLevel.Refuse);
        var basis = new BassBasis(inputs, degree, hinges);

        int draws = s.Draws.Length, reported = subsets.Length, rows = n * segments.Count;
        int[] blockOf = [.. subsets.Select(design.IndexOf)];
        var main = new double[draws, reported];
        var mainError = new double[draws, reported];
        var mainIndexError = new double[draws, reported];
        var total = new double[draws, p];
        var totalError = new double[draws, p];
        var totalIndexError = new double[draws, p];
        var variance = new double[draws];
        var varianceError = new double[draws];
        var status = new DecompositionStatus[draws];
        var f = new double[rows];
        double[] b = new double[s.MaxTerms + 1];
        double[][] values = [.. segments.Select(_ => new double[n])];
        double[][] hybrids = [.. hybridAt.Select(at => values[at])];
        double[][]? reverses = design.Reverses ? [.. reverseAt.Select(at => values[at])] : null;
        foreach ((int first, int count) in Batches(s.Draws, rows))
        {
            var surface = new BassSurface(basis, fit.Terms, cases, null);
            for (int d = first; d < first + count; d++)
            {
                BassDraw draw = s.Draws[d];
                int[] ids = draw.Structure.TermIds;
                s.Coefficients(0, draw, b);
                b[0] = 0.0;   // every piece is shift-invariant, and the intercept would inflate the rounding scale (the exact route leaves it out)
                surface.Evaluate(ids, b.AsSpan(0, ids.Length + 1), f);
                for (int g = 0; g < segments.Count; g++) f.AsSpan(g * n, n).CopyTo(values[g]);
                PickFreezeEstimate e = PickFreeze.Estimate(design, values[0], values[1], hybrids, reverses, spec.RoundingTolerance);
                variance[d] = e.Variance;
                varianceError[d] = e.VarianceError;
                for (int k = 0; k < reported; k++)
                {
                    main[d, k] = e.Main[blockOf[k]];
                    mainError[d, k] = e.MainError[blockOf[k]];
                    mainIndexError[d, k] = e.MainIndexError[blockOf[k]];
                }
                for (int v = 0; v < p; v++)
                {
                    total[d, v] = e.Total[v];
                    totalError[d, v] = e.TotalError[v];
                    totalIndexError[d, v] = e.TotalIndexError[v];
                }
                status[d] = ((e.Status & PickFreezeStatus.ZeroVariance) != 0 ? DecompositionStatus.ZeroVariance : DecompositionStatus.Defined)
                            | ((e.Status & PickFreezeStatus.BelowZeroBeyondError) != 0 ? DecompositionStatus.BelowZeroBeyondError : DecompositionStatus.Defined);
            }
        }
        var decomposition = new SobolDecomposition(measure.Names, subsets, main, total, variance, status, spec.RoundingTolerance, s.Chains, s.DrawIndex, s.Weights)
        {
            VarianceError = varianceError, MainVarianceError = mainError, TotalVarianceError = totalError,
            MainIndexError = mainIndexError, TotalIndexError = totalIndexError,
        };
        return new SobolAnalysis(null, measure, maxOrder, [decomposition], null, null, null, 0.0) { Route = SensitivityRoute.Sampled, Points = n };
    }

    // Each term's variables in the measure's numbering, by id (TermVariablesOf).
    private static Func<int, int[]> TermVariables(BassFit fit, Setup s, DesignKernels kernels)
    {
        Func<BassTerm, int[]> of = TermVariablesOf(s.Moments, kernels);
        return id => of(fit.Terms[id]);
    }

    // A term's variables in the measure's numbering: the moment service's for a separable term; otherwise the continuous variables a
    // case kernel reading the term says it reads (a projection's active coordinates), through the design-kernel contract, and the variables
    // of its categorical factors, which gate a direction (E36's slice 5) and would otherwise go unvaried.
    internal static Func<BassTerm, int[]> TermVariablesOf(BassMoments moments, DesignKernels kernels)
    {
        var variables = new Dictionary<BassTerm, int[]>();
        int continuous = moments.Inputs.Cases.ContinuousCount;
        ProductMeasure measure = moments.Measure;
        return term =>
        {
            if (variables.TryGetValue(term, out int[]? known)) return known;
            int[] v = kernels.Separable(term) ? [.. moments.Variables(term)]
                : [.. Enumerable.Range(0, continuous).Where(k => kernels.Case.Any(c => c.Reads(term) && c.ReadsVariable(term, k))).Select(measure.OfContinuous)
                     .Concat(term.Categorical.Select(c => measure.OfCategorical(c.Variable))).Order()];
            variables[term] = v;
            return v;
        };
    }

    // The effects asked for, or every input's main effect and the two-way effect of each pair some held term reads; each variable's grid
    // in scaled units (level indices for a categorical one) and its presentation; and the cell budget over `draws` × `surfaces` sets of
    // values, checked before anything is computed.
    internal static (int[][] Chosen, double[][] Scaled, EffectGrid[] Grids) EffectRequest(ProductMeasure measure, BassInputTransform transform,
        IEnumerable<IReadOnlyList<int>> held, int draws, int surfaces, EffectSpec effects)
    {
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
            foreach (IReadOnlyList<int> v in held)
                for (int i = 0; i < v.Count; i++)
                    for (int k = i + 1; k < v.Count; k++) pairs.Add((v[i], v[k]));
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
                RangeScale range = measure.OutputOf(v) is int o and >= 0 ? transform.OutputVariables[o] : (RangeScale)transform.Predictors[v];
                grids[v] = new EffectGrid([.. scaled[v].Select(x => range.Min + x * (range.Max - range.Min))], null);
            }
            else
            {
                var levels = ((TrainingLevels)transform.Predictors[v]).Levels;
                scaled[v] = [.. Enumerable.Range(0, levels.Count).Select(l => (double)l)];
                grids[v] = new EffectGrid([.. scaled[v]], levels);
            }
        }
        long cells = chosen.Sum(u => u.Aggregate(1L, (acc, v) => acc * scaled[v].Length));
        long values = cells * draws * surfaces;
        if (values > effects.CellBudget)
            throw new ArgumentException($"The effects hold {values} values ({cells} cells × {draws} draws × {surfaces} surfaces), above the budget of {effects.CellBudget}: name fewer effects, a coarser grid, or raise the budget.", nameof(effects));
        return (chosen, scaled, grids);
    }

    private const long StoreBudget = 1L << 23;   // the base and single-input blocks' values one batch of draws keeps, about 64 MB

    // Decision 5 for a direct fit: each cell's ANOVA component on the effect's inputs, by pinning (EffectCore with one component).
    private static EffectAnalysis SampledEffects(BassFit fit, SensitivitySpec spec, EffectSpec effects, Setup s)
    {
        if (spec.Points < 2) throw new ArgumentOutOfRangeException(nameof(spec), "The sampled route needs two rows in each base set.");
        if (s.Form is not null) return SampledFunctionalEffects(fit, spec, effects, s);
        BassInputs inputs = s.Moments.Inputs;
        ProductMeasure measure = s.Moments.Measure;
        int degree = fit.Config.Structure.Degree;
        HingeNormalization hinges = fit.Config.Structure.ResolvedHingeNormalization;
        var kernels = new DesignKernels(inputs, degree, hinges);
        Func<int, int[]> variables = TermVariables(fit, s, kernels);
        (int[][] chosen, double[][] scaled, EffectGrid[] grids) = EffectRequest(measure, fit.Transform, s.TermIds.Select(id => (IReadOnlyList<int>)variables(id)),
                                                                                 s.Draws.Length, s.Surfaces, effects);
        var b = new double[s.MaxTerms + 1];
        EffectState state = EffectCore(new BassBasis(inputs, degree, hinges), fit.Transform, inputs, measure, [fit.Terms], s.Draws.Length,
            (_, d) => { s.Coefficients(0, s.Draws[d], b); return (s.Draws[d].Structure.TermIds, b); },
            chosen, scaled, grids, spec, effects);
        var at0 = new EffectFunctions(measure.Names,
            [.. chosen.Select((u, k) => new Effect(u, [.. u.Select(v => grids[v])], state.Values[0][k]) { Errors = state.Error(k, [1.0]) })], state.Mean([1.0], 0.0))
        { MeanError = state.MeanError([1.0]) };
        return new EffectAnalysis(null, measure, [at0], null, null) { Route = SensitivityRoute.Sampled, Points = spec.Points };
    }

    /// <summary>
    /// The sampled effects of r components sharing their inputs and points (a direct fit is one): per effect and cell, each component's
    /// estimate and the packed covariance of the components' estimates; per draw, each component's mean, its intercept included, and the
    /// means' covariance. A combination Σₖ gₖ·eₖ is then read with its exact error √(gᵀΣg).
    /// </summary>
    internal sealed class EffectState(double[][][,] values, double[][,,] covariance, double[,] means, double[,] meanCovariance)
    {
        /// <summary>[component][effect] → draws × cells.</summary>
        public double[][][,] Values { get; } = values;

        /// <summary>[effect] → draws × cells × packed pairs k ≤ k′: the covariance of the components' cell estimates.</summary>
        public double[][,,] Covariance { get; } = covariance;

        /// <summary>Draws × components: each component's mean.</summary>
        public double[,] Means { get; } = means;

        /// <summary>Draws × packed pairs: the covariance of the components' mean estimates.</summary>
        public double[,] MeanCovariance { get; } = meanCovariance;

        public int Components => Means.GetLength(1);

        public static int Pair(int k, int k2, int r)
        {
            if (k > k2) (k, k2) = (k2, k);
            return k * r - k * (k - 1) / 2 + (k2 - k);
        }

        private double Quadratic(double[] g, Func<int, double> packed)
        {
            int r = Components;
            double q = 0.0;
            for (int k = 0; k < r; k++)
                for (int k2 = k; k2 < r; k2++) q += (k == k2 ? 1.0 : 2.0) * g[k] * g[k2] * packed(Pair(k, k2, r));
            return Math.Sqrt(Math.Max(q, 0.0));
        }

        /// <summary>The errors of Σₖ gₖ·eₖ on effect <paramref name="effect"/>'s cells, draws × cells.</summary>
        public double[,] Error(int effect, double[] g)
        {
            double[,,] c = Covariance[effect];
            int draws = c.GetLength(0), cells = c.GetLength(1);
            var e = new double[draws, cells];
            for (int d = 0; d < draws; d++)
                for (int x = 0; x < cells; x++) e[d, x] = Quadratic(g, q => c[d, x, q]);
            return e;
        }

        /// <summary>Each draw's offset + Σₖ gₖ·E[wₖ].</summary>
        public double[] Mean(double[] g, double offset)
        {
            var mean = new double[Means.GetLength(0)];
            for (int d = 0; d < mean.Length; d++)
            {
                double f = offset;
                for (int k = 0; k < Components; k++) f += g[k] * Means[d, k];
                mean[d] = f;
            }
            return mean;
        }

        /// <summary>Each draw's error of Σₖ gₖ·E[wₖ].</summary>
        public double[] MeanError(double[] g) => [.. Enumerable.Range(0, Means.GetLength(0)).Select(d => Quadratic(g, q => MeanCovariance[d, q]))];
    }

    // Decision 5's estimator over r components: each cell's ANOVA component on the effect's inputs by pinning them over the 2N shared base
    // rows, A and B as the Sobol route draws them from the same seed. A row's alternating sum Σ_{v⊆u} (−1)^{|u|−|v|} f(g_v, x_{−v}) has
    // expectation f_u(g) under the product law, so the cell is its mean over the rows (PinnedEffects' estimator), and the components' sums
    // at a row give the estimates' covariance, their sample covariance over 2N. The base block and every input's blocks at its grid values
    // are kept per batch of draws; a pair's cells are evaluated in chunks and reduced against them at once.
    internal static EffectState EffectCore(BassBasis basis, BassInputTransform transform, BassInputs inputs, ProductMeasure measure,
        IReadOnlyList<BassTermTable> terms, int draws, Func<int, int, (int[] Ids, double[] B)> coefficients,
        int[][] chosen, double[][] scaled, EffectGrid[] grids, SensitivitySpec spec, EffectSpec effects,
        Func<List<(int[] Inputs, int[] Cells)>, OutputInputs>? outputsOf = null)
    {
        int n = 2 * spec.Points, r = terms.Count, pairCount = r * (r + 1) / 2;
        object Pinned(int v, int g) => grids[v].Levels is { } levels ? levels[g] : grids[v].Values[g];
        int[] inputsUsed = [.. chosen.SelectMany(u => u).Distinct().Order()];
        var singleAt = new Dictionary<int, int>();   // an input's first block among the stored ones, the base being block 0
        var stored = new List<(int[] Inputs, int[] Cells)> { ([], []) };
        foreach (int v in inputsUsed)
        {
            singleAt[v] = stored.Count;
            for (int g = 0; g < scaled[v].Length; g++) stored.Add(([v], [g]));
        }
        long pinned = stored.Count + chosen.Where(u => u.Length == 2).Sum(u => (long)scaled[u[0]].Length * scaled[u[1]].Length);

        // The draws in batches: their stored blocks within the store budget, and the pinned blocks' term values within the evaluation
        // budget summed over the batches.
        int batch = (int)Math.Clamp(StoreBudget / ((long)stored.Count * n * r), 1, Math.Max(1, draws));
        var batches = new List<(int First, int Count, int Terms)>();
        for (int first = 0; first < draws; first += batch)
        {
            int count = Math.Min(batch, draws - first);
            var held = new HashSet<(int, int)>();
            for (int d = first; d < first + count; d++)
                for (int k = 0; k < r; k++)
                    foreach (int id in coefficients(k, d).Ids) held.Add((k, id));
            batches.Add((first, count, held.Count));
        }
        long evaluations = batches.Sum(b => pinned * n * (long)b.Terms);
        if (evaluations > effects.EvaluationBudget)
            throw new ArgumentException($"The sampled effects evaluate about {evaluations:E2} term values ({pinned} pinned blocks of {n} rows, for the "
                + $"terms of {batches.Count} batches of draws), above the budget of {effects.EvaluationBudget:E2}: name fewer effects, a coarser "
                + "grid, fewer points, or raise the budget.", nameof(effects));

        PredictionCases drawn = BassActiveSubspaces.Points(transform, measure, n, spec.Seed);
        int[] widths = [.. chosen.Select(u => u.Aggregate(1, (acc, v) => acc * scaled[v].Length))];
        double[][][,] values = [.. Enumerable.Range(0, r).Select(_ => widths.Select(w => new double[draws, w]).ToArray())];
        double[][,,] covariance = [.. widths.Select(w => new double[draws, w, pairCount])];
        var means = new double[draws, r];
        var meanCovariance = new double[draws, pairCount];
        var intercept = new double[draws, r];
        double[][] e = [.. Enumerable.Range(0, r).Select(_ => new double[n])];

        // The components' alternating sums at one cell (e), reduced to their means and their estimates' covariance.
        void Reduce(int d, int effect, int cell)
        {
            var mean = new double[r];
            for (int k = 0; k < r; k++)
            {
                mean[k] = PickFreeze.Mean(e[k]);
                values[k][effect][d, cell] = mean[k];
            }
            for (int k = 0; k < r; k++)
                for (int k2 = k; k2 < r; k2++)
                {
                    double s = 0.0;
                    for (int j = 0; j < n; j++) s += (e[k][j] - mean[k]) * (e[k2][j] - mean[k2]);
                    covariance[effect][d, cell, EffectState.Pair(k, k2, r)] = s / ((n - 1.0) * n);
                }
        }

        // Evaluates `blocks` for the batch's draws in chunks whose term columns fit the column budget, handing each draw's components' values on.
        void Evaluate(List<(int[] Inputs, int[] Cells)> blocks, (int First, int Count, int Terms) at, Action<int, int, double[][], int> take)
        {
            int perChunk = (int)Math.Clamp(ColumnBudget / ((long)n * Math.Max(1, at.Terms)), 1, Math.Max(1, blocks.Count));
            if ((long)perChunk * n > int.MaxValue) throw new InvalidOperationException("A chunk outgrows an array.");
            double[][] f = [.. Enumerable.Range(0, r).Select(_ => new double[perChunk * n])];
            for (int start = 0; start < blocks.Count; start += perChunk)
            {
                var chunk = blocks.GetRange(start, Math.Min(perChunk, blocks.Count - start));
                (CaseInputs cases, _) = Pinning(drawn, measure, n, chunk, Pinned).Prepare(transform, inputs, UnseenLevel.Refuse);
                OutputInputs? locations = outputsOf?.Invoke(chunk);   // a direct functional fit's rows' output locations
                BassSurface[] surfaces = [.. terms.Select(t => new BassSurface(basis, t, cases, locations))];
                for (int d = at.First; d < at.First + at.Count; d++)
                {
                    for (int k = 0; k < r; k++)
                    {
                        (int[] ids, double[] b) = coefficients(k, d);
                        intercept[d, k] = b[0];
                        b[0] = 0.0;   // the alternating sums' signs add to zero; the mean takes the intercept back
                        if (locations is null) surfaces[k].Evaluate(ids, b.AsSpan(0, ids.Length + 1), f[k].AsSpan(0, chunk.Count * n));
                        else RowValues(surfaces[k], ids, b, f[k].AsSpan(0, chunk.Count * n));
                    }
                    for (int c = 0; c < chunk.Count; c++) take(d, start + c, f, c * n);
                }
            }
        }

        foreach ((int First, int Count, int Terms) at in batches)
        {
            var store = new double[at.Count, r][];
            for (int d = 0; d < at.Count; d++)
                for (int k = 0; k < r; k++) store[d, k] = new double[stored.Count * n];
            Evaluate(stored, at, (d, block, f, offset) =>
            {
                for (int k = 0; k < r; k++) f[k].AsSpan(offset, n).CopyTo(store[d - at.First, k].AsSpan(block * n, n));
            });

            for (int d = at.First; d < at.First + at.Count; d++)
            {
                // Each component's mean over the base rows, its intercept taken back, and the means' covariance.
                for (int k = 0; k < r; k++)
                {
                    store[d - at.First, k].AsSpan(0, n).CopyTo(e[k]);
                    means[d, k] = intercept[d, k] + PickFreeze.Mean(e[k]);
                }
                for (int k = 0; k < r; k++)
                    for (int k2 = k; k2 < r; k2++)
                    {
                        double mk = PickFreeze.Mean(e[k]), mk2 = PickFreeze.Mean(e[k2]), s = 0.0;
                        for (int j = 0; j < n; j++) s += (e[k][j] - mk) * (e[k2][j] - mk2);
                        meanCovariance[d, EffectState.Pair(k, k2, r)] = s / ((n - 1.0) * n);
                    }
                for (int x = 0; x < chosen.Length; x++)
                {
                    if (chosen[x].Length != 1) continue;
                    int v = chosen[x][0];
                    for (int g = 0; g < scaled[v].Length; g++)
                    {
                        for (int k = 0; k < r; k++)
                        {
                            ReadOnlySpan<double> own = store[d - at.First, k];
                            ReadOnlySpan<double> pinnedRows = own.Slice((singleAt[v] + g) * n, n), baseRows = own[..n];
                            for (int j = 0; j < n; j++) e[k][j] = pinnedRows[j] - baseRows[j];
                        }
                        Reduce(d, x, g);
                    }
                }
            }

            for (int x = 0; x < chosen.Length; x++)
            {
                if (chosen[x].Length != 2) continue;
                int v = chosen[x][0], w = chosen[x][1], mw = scaled[w].Length;
                var cells = new List<(int[] Inputs, int[] Cells)>();
                for (int g = 0; g < scaled[v].Length; g++)
                    for (int h = 0; h < mw; h++) cells.Add(([v, w], [g, h]));
                int effect = x;
                Evaluate(cells, at, (d, cell, f, offset) =>
                {
                    int g = cell / mw, h = cell % mw;
                    for (int k = 0; k < r; k++)
                    {
                        ReadOnlySpan<double> own = store[d - at.First, k];
                        ReadOnlySpan<double> baseRows = own[..n], atV = own.Slice((singleAt[v] + g) * n, n), atW = own.Slice((singleAt[w] + h) * n, n);
                        ReadOnlySpan<double> both = f[k].AsSpan(offset, n);
                        for (int j = 0; j < n; j++) e[k][j] = both[j] - atV[j] - atW[j] + baseRows[j];
                    }
                    Reduce(d, effect, cell);
                });
            }
        }
        return new EffectState(values, covariance, means, meanCovariance);
    }

    // ── A direct functional fit's forms (slice 3, decision 8) ─────────────────────────────────────────────────────────────────

    // Each integrated output variable's 2N values in box units, drawn from its law as E8's sampler draws a predictor's, from a stream apart
    // from the predictors'; null for an output variable the form holds.
    private static double[]?[] OutputDraws(BassInputTransform transform, ProductMeasure measure, int count, int seed)
    {
        var draws = new double[]?[transform.OutputVariables.Count];
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(seed, 2)[1]);
        for (int v = 0; v < measure.Count; v++)
        {
            if (measure.OutputOf(v) is not (int o and >= 0)) continue;
            ContinuousMeasure law = measure.Continuous(v)!;
            var values = new double[count];
            for (int r = 0; r < count; r++)
                values[r] = law.Measure is Empirical e ? e.Values[rng.NextInt(e.Count)] : ((IContinuousDistribution)law.Measure).Sample(rng);
            draws[o] = values;
        }
        return draws;
    }

    // A surface's output locations at stacked rows: a held output variable at the coordinate's scaled value, an integrated one at each
    // row's own value, `rowOf` giving the drawn row (or a pinned value through `pinned`).
    private static OutputInputs Locations(BassInputTransform transform, ProductMeasure measure, double[]?[] draws, double[]? held, IReadOnlyList<string>? heldNames,
                                          int rows, Func<int, int, double> integrated)
    {
        int count = transform.OutputVariables.Count;
        var columns = new double[count][];
        for (int o = 0; o < count; o++)
        {
            columns[o] = new double[rows];
            if (draws[o] is not null)
            {
                int v = Enumerable.Range(0, measure.Count).First(x => measure.OutputOf(x) == o);
                for (int r = 0; r < rows; r++) columns[o][r] = integrated(v, r);
            }
            else
            {
                RangeScale range = transform.OutputVariables[o];
                int at = heldNames!.ToList().IndexOf(range.Name);
                Array.Fill(columns[o], range.Scale(held![at]));
            }
        }
        return new OutputInputs(columns);
    }

    // f at each stacked row r of a direct functional fit: Σⱼ bⱼ·aⱼ(x_r)·hⱼ(t_r), each term's case and output factors read once per
    // surface (BassSurface's columns), the intercept left out.
    private static void RowValues(BassSurface surface, int[] ids, ReadOnlySpan<double> b, Span<double> f)
    {
        f.Clear();
        for (int j = 0; j < ids.Length; j++)
        {
            double bj = b[1 + j];
            if (bj == 0.0) continue;
            (double[] a, double[]? h) = surface.Columns(ids[j]);
            for (int r = 0; r < f.Length; r++) f[r] += bj * a[r] * h![r];
        }
    }

    // The integrated form reads the output variables as more inputs, drawn with the points; the conditioned form holds them (or those
    // ConditionOn names) at each coordinate, a decomposition per coordinate. Each row's value is Σⱼ aⱼ·Bⱼ(x_r)·hⱼ(t_r) at its own output
    // location, so the pick-freeze blocks swap integrated output variables as they swap inputs.
    private static SobolAnalysis SampledFunctional(BassFit fit, SensitivitySpec spec, Setup s, MainVarianceForm? form)
    {
        if (spec.Points < 2) throw new ArgumentOutOfRangeException(nameof(spec), "The sampled route needs two rows in each base set.");
        BassInputs inputs = s.Moments.Inputs;
        ProductMeasure measure = s.Moments.Measure;
        int p = measure.Count, n = spec.Points, draws = s.Draws.Length, surfaces = s.Surfaces;
        int degree = fit.Config.Structure.Degree;
        HingeNormalization hinges = fit.Config.Structure.ResolvedHingeNormalization;
        var kernels = new DesignKernels(inputs, degree, hinges);
        Func<int, int[]> variables = TermVariables(fit, s, kernels);
        int maxOrder = spec.MaxOrder ?? Math.Max(1, s.TermIds.Select(id => variables(id).Length).DefaultIfEmpty(0).Max());
        int[][] subsets = Subsets(s.TermIds.Select(id => (IReadOnlyList<int>)variables(id)), maxOrder, spec.SubsetBudget);
        PickFreezeDesign design = form is { } chosen ? new PickFreezeDesign(p, subsets, chosen, p + spec.SubsetBudget)
                                                     : new PickFreezeDesign(p, subsets, budget: p + spec.SubsetBudget);
        var segments = new List<(bool FromB, int[] Other)> { (false, []), (true, []) };
        var hybridAt = new int[design.Blocks.Count];
        var reverseAt = new int[design.Blocks.Count];
        for (int k = 0; k < design.Blocks.Count; k++)
        {
            hybridAt[k] = segments.Count;
            segments.Add((false, design.Blocks[k]));
            if (!design.Reverses) continue;
            reverseAt[k] = segments.Count;
            segments.Add((true, design.Blocks[k]));
        }
        int rows = n * segments.Count;
        PredictionCases drawn = BassActiveSubspaces.Points(fit.Transform, measure, 2 * n, spec.Seed);
        (CaseInputs cases, _) = Stacked(drawn, measure, n, segments).Prepare(fit.Transform, inputs, UnseenLevel.Refuse);
        double[]?[] outputDraws = OutputDraws(fit.Transform, measure, 2 * n, spec.Seed);
        OutputInputs[] locations = [.. Enumerable.Range(0, surfaces).Select(c => Locations(fit.Transform, measure, outputDraws, s.Coordinates?[c], s.CoordinateVariables, rows,
            (v, r) =>
            {
                (bool fromB, int[] other) = segments[r / n];
                return outputDraws[measure.OutputOf(v)]![(fromB != (Array.IndexOf(other, v) >= 0) ? n : 0) + r % n];
            }))];
        var basis = new BassBasis(inputs, degree, hinges);

        int[] reported = [.. subsets.Select(design.IndexOf)];
        SampledTables[] tables = [.. Enumerable.Range(0, surfaces).Select(_ => new SampledTables(draws, subsets.Length, p))];
        var f = new double[rows];
        double[][] values = [.. segments.Select(_ => new double[n])];
        double[][] hybrids = [.. hybridAt.Select(at => values[at])];
        double[][]? reverses = design.Reverses ? [.. reverseAt.Select(at => values[at])] : null;
        foreach ((int first, int count) in Batches(s.Draws, rows * 2 * surfaces))
        {
            BassSurface[] surface = [.. locations.Select(l => new BassSurface(basis, fit.Terms, cases, l))];
            for (int d = first; d < first + count; d++)
            {
                BassDraw draw = s.Draws[d];
                for (int c = 0; c < surfaces; c++)
                {
                    RowValues(surface[c], draw.Structure.TermIds, draw.Coefficients, f);
                    for (int g = 0; g < segments.Count; g++) f.AsSpan(g * n, n).CopyTo(values[g]);
                    tables[c].Take(d, PickFreeze.Estimate(design, values[0], values[1], hybrids, reverses, spec.RoundingTolerance), reported);
                }
            }
        }
        SobolDecomposition[] decompositions = [.. tables.Select(t => t.Decomposition(measure.Names, subsets, spec.RoundingTolerance, s.Chains, s.DrawIndex, s.Weights))];
        return new SobolAnalysis(s.Form, measure, maxOrder, decompositions, s.Coordinates, s.CoordinatesOutside, s.CoordinateVariables, 0.0)
        {
            Route = SensitivityRoute.Sampled, Points = n,
        };
    }

    // A direct functional fit's effects: the integrated form's output variables pinned and drawn as inputs, or the conditioned form's
    // held at each coordinate, one set of effects per coordinate; each row's value at its own output location.
    private static EffectAnalysis SampledFunctionalEffects(BassFit fit, SensitivitySpec spec, EffectSpec effects, Setup s)
    {
        BassInputs inputs = s.Moments.Inputs;
        ProductMeasure measure = s.Moments.Measure;
        int degree = fit.Config.Structure.Degree, n = 2 * spec.Points;
        HingeNormalization hinges = fit.Config.Structure.ResolvedHingeNormalization;
        var kernels = new DesignKernels(inputs, degree, hinges);
        Func<int, int[]> variables = TermVariables(fit, s, kernels);
        (int[][] chosen, double[][] scaled, EffectGrid[] grids) = EffectRequest(measure, fit.Transform, s.TermIds.Select(id => (IReadOnlyList<int>)variables(id)),
                                                                                 s.Draws.Length, s.Surfaces, effects);
        double[]?[] outputDraws = OutputDraws(fit.Transform, measure, n, spec.Seed);
        var basis = new BassBasis(inputs, degree, hinges);
        var at = new EffectFunctions[s.Surfaces];
        for (int c = 0; c < s.Surfaces; c++)
        {
            int surface = c;
            OutputInputs Outputs(List<(int[] Inputs, int[] Cells)> blocks) => Locations(fit.Transform, measure, outputDraws, s.Coordinates?[surface], s.CoordinateVariables,
                blocks.Count * n, (v, r) =>
                {
                    (int[] pinnedInputs, int[] cells) = blocks[r / n];
                    int position = Array.IndexOf(pinnedInputs, v);
                    return position >= 0 ? scaled[v][cells[position]] : outputDraws[measure.OutputOf(v)]![r % n];
                });
            EffectState state = EffectCore(basis, fit.Transform, inputs, measure, [fit.Terms], s.Draws.Length,
                (_, d) => (s.Draws[d].Structure.TermIds, (double[])s.Draws[d].Coefficients.Clone()), chosen, scaled, grids, spec, effects, Outputs);
            at[c] = new EffectFunctions(measure.Names,
                [.. chosen.Select((u, k) => new Effect(u, [.. u.Select(v => grids[v])], state.Values[0][k]) { Errors = state.Error(k, [1.0]) })], state.Mean([1.0], 0.0))
            { MeanError = state.MeanError([1.0]) };
        }
        return new EffectAnalysis(s.Form, measure, at, s.Coordinates, s.CoordinateVariables) { Route = SensitivityRoute.Sampled, Points = spec.Points };
    }

    // The base rows stacked once per block, each block's inputs pinned to its cells' values.
    private static PredictionCases Pinning(PredictionCases drawn, ProductMeasure measure, int n, List<(int[] Inputs, int[] Cells)> blocks,
                                           Func<int, int, object> pinned)
    {
        var columns = new List<PredictionColumn>();
        foreach (PredictionColumn column in drawn.Columns)
        {
            int v = Array.IndexOf([.. measure.Names], column.Name);
            int length = n * blocks.Count;
            switch (column)
            {
                case ContinuousColumn c:
                {
                    ReadOnlySpan<double> source = c.Values.Span;
                    var stacked = new double[length];
                    for (int k = 0; k < blocks.Count; k++)
                    {
                        int at = Array.IndexOf(blocks[k].Inputs, v);
                        if (at >= 0) stacked.AsSpan(k * n, n).Fill((double)pinned(v, blocks[k].Cells[at]));
                        else source.CopyTo(stacked.AsSpan(k * n, n));
                    }
                    columns.Add(new ContinuousColumn(c.Name, stacked));
                    break;
                }
                case CategoricalColumn c:
                {
                    var stacked = new string[length];
                    for (int k = 0; k < blocks.Count; k++)
                    {
                        int at = Array.IndexOf(blocks[k].Inputs, v);
                        for (int j = 0; j < n; j++) stacked[k * n + j] = at >= 0 ? (string)pinned(v, blocks[k].Cells[at]) : c.Levels[j];
                    }
                    columns.Add(new CategoricalColumn(c.Name, stacked));
                    break;
                }
                default:
                    throw new InvalidOperationException($"Column '{column.Name}' is of no known kind.");
            }
        }
        return new PredictionCases(columns);
    }

    // The draws in batches whose distinct terms' columns over the stacked rows fit the column budget, one draw at least.
    private static IEnumerable<(int First, int Count)> Batches(BassDraw[] draws, int rows)
    {
        int first = 0;
        var ids = new HashSet<int>();
        for (int d = 0; d < draws.Length; d++)
        {
            foreach (int id in draws[d].Structure.TermIds) ids.Add(id);
            if (d > first && (long)ids.Count * rows > ColumnBudget)
            {
                yield return (first, d - first);
                first = d;
                ids = [.. draws[d].Structure.TermIds];
            }
        }
        if (draws.Length > first) yield return (first, draws.Length - first);
    }

    // The base sets' rows stacked by segment: a segment takes its inputs from A (B when FromB), and those in Other from the other set.
    // A's rows are the drawn points' first N, B's the next N.
    internal static PredictionCases Stacked(PredictionCases drawn, ProductMeasure measure, int n, List<(bool FromB, int[] Other)> segments)
    {
        var columns = new List<PredictionColumn>();
        foreach (PredictionColumn column in drawn.Columns)
        {
            int v = Array.IndexOf([.. measure.Names], column.Name);
            int Source(int segment, int row)
            {
                (bool fromB, int[] other) = segments[segment];
                return (fromB != (Array.IndexOf(other, v) >= 0) ? n : 0) + row;
            }
            int length = n * segments.Count;
            switch (column)
            {
                case ContinuousColumn c:
                {
                    ReadOnlySpan<double> source = c.Values.Span;
                    var stacked = new double[length];
                    for (int g = 0; g < segments.Count; g++)
                        for (int j = 0; j < n; j++) stacked[g * n + j] = source[Source(g, j)];
                    columns.Add(new ContinuousColumn(c.Name, stacked));
                    break;
                }
                case CategoricalColumn c:
                {
                    var stacked = new string[length];
                    for (int g = 0; g < segments.Count; g++)
                        for (int j = 0; j < n; j++) stacked[g * n + j] = c.Levels[Source(g, j)];
                    columns.Add(new CategoricalColumn(c.Name, stacked));
                    break;
                }
                default:
                    throw new InvalidOperationException($"Column '{column.Name}' is of no known kind.");
            }
        }
        return new PredictionCases(columns);
    }
}
