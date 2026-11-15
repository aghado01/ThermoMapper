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

        // Each term's variables in the measure's numbering: the moment service's for a separable term; otherwise the continuous variables
        // a case kernel reading the term says it reads (a projection's active coordinates), through the design-kernel contract.
        var variables = new Dictionary<int, int[]>();
        int[] Variables(int id)
        {
            if (variables.TryGetValue(id, out int[]? known)) return known;
            BassTerm term = fit.Terms[id];
            int[] v = kernels.Separable(term) ? [.. s.Moments.Variables(term)]
                : [.. Enumerable.Range(0, inputs.Cases.ContinuousCount)
                       .Where(k => kernels.Case.Any(c => c.Reads(term) && c.ReadsVariable(term, k))).Select(measure.OfContinuous).Order()];
            variables[id] = v;
            return v;
        }
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
