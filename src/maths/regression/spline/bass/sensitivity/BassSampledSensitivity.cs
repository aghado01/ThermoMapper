using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Sensitivity;

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
    private const long ColumnBudget = 1L << 23;   // term-column values one batch of draws keeps, about 64 MB

    /// <summary><see cref="Of"/> with the main variances read in <paramref name="form"/> on the sampled route: the gates read both forms
    /// while the form's default waits on the user (the slice 0 report's follow-up).</summary>
    internal static SobolAnalysis OfInForm(BassFit fit, BassData data, SensitivitySpec? spec, MainVarianceForm form) => OfCore(fit, data, spec, null, form);

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

    private static SobolAnalysis Sampled(BassFit fit, SensitivitySpec spec, Setup s, bool projected, MainVarianceForm? form)
    {
        if (projected) throw new ArgumentException("The sampled route reads a scalar fit for now: several outputs are the sampled-sensitivity brief's slice 3.", nameof(fit));
        if (s.Form is not null) throw new ArgumentException("The sampled route reads a scalar fit for now: a functional fit's forms are the sampled-sensitivity brief's slice 3.", nameof(fit));
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

    // Each term's variables in the measure's numbering, by id: the moment service's for a separable term; otherwise the continuous
    // variables a case kernel reading the term says it reads (a projection's active coordinates), through the design-kernel contract.
    private static Func<int, int[]> TermVariables(BassFit fit, Setup s, DesignKernels kernels)
    {
        var variables = new Dictionary<int, int[]>();
        int continuous = s.Moments.Inputs.Cases.ContinuousCount;
        ProductMeasure measure = s.Moments.Measure;
        return id =>
        {
            if (variables.TryGetValue(id, out int[]? known)) return known;
            BassTerm term = fit.Terms[id];
            int[] v = kernels.Separable(term) ? [.. s.Moments.Variables(term)]
                : [.. Enumerable.Range(0, continuous).Where(k => kernels.Case.Any(c => c.Reads(term) && c.ReadsVariable(term, k))).Select(measure.OfContinuous).Order()];
            variables[id] = v;
            return v;
        };
    }

    // The effects asked for, or every input's main effect and the two-way effect of each pair some term holds; each variable's grid in
    // scaled units (level indices for a categorical one) and its presentation; and the cell budget, checked before anything is computed.
    private static (int[][] Chosen, double[][] Scaled, EffectGrid[] Grids) EffectRequest(BassFit fit, Setup s, EffectSpec effects,
                                                                                      Func<int, IReadOnlyList<int>> variablesOf)
    {
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
                IReadOnlyList<int> v = variablesOf(id);
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
        return (chosen, scaled, grids);
    }

    private const long StoreBudget = 1L << 23;   // the base and single-input blocks' values one batch of draws keeps, about 64 MB

    // Decision 5: each cell's ANOVA component on the effect's inputs by pinning them over the 2N shared base rows, A and B as the Sobol
    // route draws them from the same seed. A row's alternating sum Σ_{v⊆u} (−1)^{|u|−|v|} f(g_v, x_{−v}) has expectation f_u(g) under the
    // product law, so the cell is its mean over the rows, with error sd/√2N (PinnedEffects' estimator). The base block and every input's
    // blocks at its grid values are kept per draw; a pair's cells are evaluated in chunks and reduced against them at once.
    private static EffectAnalysis SampledEffects(BassFit fit, SensitivitySpec spec, EffectSpec effects, Setup s)
    {
        if (s.Form is not null) throw new ArgumentException("The sampled route reads a scalar fit for now: a functional fit's forms are the sampled-sensitivity brief's slice 3.", nameof(fit));
        if (spec.Points < 2) throw new ArgumentOutOfRangeException(nameof(spec), "The sampled route needs two rows in each base set.");
        BassInputs inputs = s.Moments.Inputs;
        ProductMeasure measure = s.Moments.Measure;
        int degree = fit.Config.Structure.Degree;
        HingeNormalization hinges = fit.Config.Structure.ResolvedHingeNormalization;
        var kernels = new DesignKernels(inputs, degree, hinges);
        Func<int, int[]> variables = TermVariables(fit, s, kernels);
        (int[][] chosen, double[][] scaled, EffectGrid[] grids) = EffectRequest(fit, s, effects, id => variables(id));
        int n = 2 * spec.Points;

        // The pinned values each block holds, in the inputs' own units (a level's name for a categorical input).
        object Pinned(int v, int g) => grids[v].Levels is { } levels ? levels[g] : grids[v].Values[g];
        int[] inputsUsed = [.. chosen.SelectMany(u => u).Distinct().Order()];
        var singleAt = new Dictionary<int, int>();   // an input's first block among the stored ones, the base being block 0
        var stored = new List<(int[] Inputs, int[] Cells)> { ([], []) };
        foreach (int v in inputsUsed)
        {
            singleAt[v] = stored.Count;
            for (int g = 0; g < scaled[v].Length; g++) stored.Add(([v], [g]));
        }
        int[][] pairs = [.. chosen.Where(u => u.Length == 2)];
        long pinned = stored.Count + pairs.Sum(u => (long)scaled[u[0]].Length * scaled[u[1]].Length);

        // The draws in batches: their stored blocks within the store budget, and the pinned blocks' term values within the evaluation
        // budget summed over the batches.
        int draws = s.Draws.Length;
        int batch = (int)Math.Clamp(StoreBudget / ((long)stored.Count * n), 1, Math.Max(1, draws));
        var batches = new List<(int First, int Count, int Terms)>();
        for (int first = 0; first < draws; first += batch)
        {
            int count = Math.Min(batch, draws - first);
            batches.Add((first, count, s.Draws.Skip(first).Take(count).SelectMany(d => d.Structure.TermIds).Distinct().Count()));
        }
        long evaluations = batches.Sum(b => pinned * n * (long)b.Terms);
        if (evaluations > effects.EvaluationBudget)
            throw new ArgumentException($"The sampled effects evaluate about {evaluations:E2} term values ({pinned} pinned blocks of {n} rows, for the "
                + $"terms of {batches.Count} batches of draws), above the budget of {effects.EvaluationBudget:E2}: name fewer effects, a coarser "
                + "grid, fewer points, or raise the budget.", nameof(effects));

        PredictionCases drawn = BassActiveSubspaces.Points(fit.Transform, measure, n, spec.Seed);
        var basis = new BassBasis(inputs, degree, hinges);
        double[][,] values = [.. chosen.Select(u => new double[draws, u.Aggregate(1, (acc, v) => acc * scaled[v].Length)])];
        double[][,] errors = [.. chosen.Select(u => new double[draws, u.Aggregate(1, (acc, v) => acc * scaled[v].Length)])];
        var mean = new double[draws];
        var meanError = new double[draws];
        var e = new double[n];
        double[] b = new double[s.MaxTerms + 1];
        double[] intercept = new double[draws];

        // Evaluates `blocks` for the batch's draws in chunks whose term columns fit the column budget, handing each draw's values on.
        void Evaluate(List<(int[] Inputs, int[] Cells)> blocks, (int First, int Count, int Terms) at, Action<int, int, ReadOnlySpan<double>> take)
        {
            int perChunk = (int)Math.Clamp(ColumnBudget / ((long)n * Math.Max(1, at.Terms)), 1, Math.Max(1, blocks.Count));
            var f = new double[(long)perChunk * n > int.MaxValue ? throw new InvalidOperationException("A chunk outgrows an array.") : perChunk * n];
            for (int start = 0; start < blocks.Count; start += perChunk)
            {
                var chunk = blocks.GetRange(start, Math.Min(perChunk, blocks.Count - start));
                (CaseInputs cases, _) = Pinning(drawn, measure, n, chunk, Pinned).Prepare(fit.Transform, inputs, UnseenLevel.Refuse);
                var surface = new BassSurface(basis, fit.Terms, cases, null);
                for (int d = at.First; d < at.First + at.Count; d++)
                {
                    BassDraw draw = s.Draws[d];
                    int[] ids = draw.Structure.TermIds;
                    s.Coefficients(0, draw, b);
                    intercept[d] = b[0];
                    b[0] = 0.0;   // the alternating sums' signs add to zero; the mean takes the intercept back
                    surface.Evaluate(ids, b.AsSpan(0, ids.Length + 1), f.AsSpan(0, chunk.Count * n));
                    for (int k = 0; k < chunk.Count; k++) take(d, start + k, f.AsSpan(k * n, n));
                }
            }
        }

        foreach ((int First, int Count, int Terms) at in batches)
        {
            var store = new double[at.Count][];
            for (int d = 0; d < at.Count; d++) store[d] = new double[stored.Count * n];
            Evaluate(stored, at, (d, block, f) => f.CopyTo(store[d - at.First].AsSpan(block * n, n)));

            for (int d = at.First; d < at.First + at.Count; d++)
            {
                double[] own = store[d - at.First];
                ReadOnlySpan<double> baseRows = own.AsSpan(0, n);
                mean[d] = intercept[d] + PickFreeze.Mean(baseRows);
                meanError[d] = PickFreeze.Deviation(baseRows) / Math.Sqrt(n);
                for (int k = 0; k < chosen.Length; k++)
                {
                    if (chosen[k].Length != 1) continue;
                    int v = chosen[k][0];
                    for (int g = 0; g < scaled[v].Length; g++)
                    {
                        ReadOnlySpan<double> pinnedRows = own.AsSpan((singleAt[v] + g) * n, n);
                        for (int j = 0; j < n; j++) e[j] = pinnedRows[j] - baseRows[j];
                        values[k][d, g] = PickFreeze.Mean(e);
                        errors[k][d, g] = PickFreeze.Deviation(e) / Math.Sqrt(n);
                    }
                }
            }

            for (int k = 0; k < chosen.Length; k++)
            {
                if (chosen[k].Length != 2) continue;
                int v = chosen[k][0], w = chosen[k][1], mw = scaled[w].Length;
                var cells = new List<(int[] Inputs, int[] Cells)>();
                for (int g = 0; g < scaled[v].Length; g++)
                    for (int h = 0; h < mw; h++) cells.Add(([v, w], [g, h]));
                int effect = k;
                Evaluate(cells, at, (d, cell, f) =>
                {
                    double[] own = store[d - at.First];
                    int g = cell / mw, h = cell % mw;
                    ReadOnlySpan<double> baseRows = own.AsSpan(0, n), atV = own.AsSpan((singleAt[v] + g) * n, n), atW = own.AsSpan((singleAt[w] + h) * n, n);
                    for (int j = 0; j < n; j++) e[j] = f[j] - atV[j] - atW[j] + baseRows[j];
                    values[effect][d, cell] = PickFreeze.Mean(e);
                    errors[effect][d, cell] = PickFreeze.Deviation(e) / Math.Sqrt(n);
                });
            }
        }

        var at0 = new EffectFunctions(measure.Names,
            [.. chosen.Select((u, k) => new Effect(u, [.. u.Select(v => grids[v])], values[k]) { Errors = errors[k] })], mean) { MeanError = meanError };
        return new EffectAnalysis(null, measure, [at0], null, null) { Route = SensitivityRoute.Sampled, Points = spec.Points };
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
    private static PredictionCases Stacked(PredictionCases drawn, ProductMeasure measure, int n, List<(bool FromB, int[] Other)> segments)
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
