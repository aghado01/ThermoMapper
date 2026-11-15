using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Sensitivity;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// A basis fit's sensitivity on the sampled route (the sampled-sensitivity brief, slice 3, decision 8), for joint draws whose
/// components hold a term without one-variable factors. It accumulates each component's surface wₖ at shared pick-freeze points,
/// and reduces each location's reconstructed surface Σₖ g_ℓk·wₖ, and the integrated form's, by <see cref="PickFreeze"/>, every piece
/// with its error. The pieces are bilinear in the components, so each location's estimate is the exact route's form read at the
/// points; the integrated form reads the location as one more input, drawn uniformly over the training locations beside each row.
/// </summary>
public static partial class BassBasisSensitivity
{
    // Decision 1 for a basis fit: exact when every term the components' selected draws hold has one-variable factors, sampled
    // otherwise; Exact refuses the term it meets up front.
    private static SensitivityRoute Resolve(SensitivitySpec spec, HashSet<BassTerm> held, DesignKernels kernels)
    {
        BassTerm? without = held.FirstOrDefault(t => !kernels.Separable(t));
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

    /// <summary><see cref="Of"/> with the main variances read in <paramref name="form"/> on the sampled route: the gates read both forms
    /// while the form's default waits on the user.</summary>
    internal static BasisSobolAnalysis OfInForm(BassBasisFit fit, SensitivitySpec spec, MainVarianceForm form) => OfCore(fit, spec, form);

    private static BasisSobolAnalysis Sampled(BassBasisFit fit, SensitivitySpec spec, BassInputs inputs, ProductMeasure measure, BassMoments moments,
                                              DesignKernels kernels, BassJointDraws joint, int[] selected, List<BassDraw>[] pooled,
                                              HashSet<BassTerm> held, int[] chains, int[] within, MainVarianceForm? form)
    {
        if (spec.Points < 2) throw new ArgumentOutOfRangeException(nameof(spec), "The sampled route needs two rows in each base set.");
        int r = fit.Components.Count, p = measure.Count, n = spec.Points, draws = selected.Length;
        BassFit first = fit.Components[0];
        LinearReconstruction map = fit.Decomposition.Map;
        int m = map.Locations;

        // The subsets as the exact route reports them, at a location and in the integrated form (the location last, index p).
        Func<BassTerm, int[]> variables = BassSensitivity.TermVariablesOf(moments, kernels);
        int largest = held.Select(t => variables(t).Length).DefaultIfEmpty(0).Max();
        int maxOrder = spec.MaxOrder ?? Math.Max(1, largest);
        int integratedOrder = spec.MaxOrder ?? largest + 1;
        int[][] subsets = BassSensitivity.Subsets(held.Select(t => (IReadOnlyList<int>)variables(t)), maxOrder, spec.SubsetBudget);
        var kinds = new List<int[]>();
        foreach (int[] u in subsets)
        {
            if (u.Length <= integratedOrder) kinds.Add(u);
            if (u.Length + 1 <= integratedOrder) kinds.Add([.. u, p]);
        }
        kinds.Add([p]);
        kinds.Sort((x, y) => x.Length != y.Length ? x.Length.CompareTo(y.Length) : x.AsSpan().SequenceCompareTo(y));
        PickFreezeDesign Design(int inputs, IReadOnlyList<int[]> asked) =>
            form is { } chosen ? new PickFreezeDesign(inputs, asked, chosen, inputs + spec.SubsetBudget) : new PickFreezeDesign(inputs, asked, budget: inputs + spec.SubsetBudget);
        PickFreezeDesign atLocation = Design(p, subsets), integrated = Design(p + 1, kinds);

        // The stacked x-rows: A, B, then each location block reading only inputs some term reads (hybrid, then reverse); a block holding
        // an unread input is a singleton whose values are A's (its reverse B's). The integrated blocks' inputs are among them.
        var read = new bool[p];
        foreach (BassTerm t in held)
            foreach (int v in variables(t)) read[v] = true;
        var segments = new List<(bool FromB, int[] Other)> { (false, []), (true, []) };
        var hybridAt = new int[atLocation.Blocks.Count];
        var reverseAt = new int[atLocation.Blocks.Count];
        for (int k = 0; k < atLocation.Blocks.Count; k++)
        {
            int[] u = atLocation.Blocks[k];
            if (!u.All(v => read[v]))
            {
                (hybridAt[k], reverseAt[k]) = (0, 1);
                continue;
            }
            hybridAt[k] = segments.Count;
            segments.Add((false, u));
            if (!atLocation.Reverses) continue;
            reverseAt[k] = segments.Count;
            segments.Add((true, u));
        }
        int XSegment(int[] u, bool reverse)
        {
            int[] x = [.. u.Where(v => v < p)];
            if (x.Length == 0) return reverse ? 1 : 0;
            int k = atLocation.IndexOf(x);
            return reverse ? reverseAt[k] : hybridAt[k];
        }
        PredictionCases drawn = BassActiveSubspaces.Points(first.Transform, measure, 2 * n, spec.Seed);
        (CaseInputs cases, _) = BassSensitivity.Stacked(drawn, measure, n, segments).Prepare(first.Transform, inputs, UnseenLevel.Refuse);
        var basis = new BassBasis(inputs, first.Config.Structure.Degree, BassBasisEffects.OneNormalization(fit));

        // Each row's location for the integrated form, uniform over the training locations, from a stream apart from the points'.
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(spec.Seed, 1)[0]);
        int[] locationA = [.. Enumerable.Range(0, n).Select(_ => rng.NextInt(m))], locationB = [.. Enumerable.Range(0, n).Select(_ => rng.NextInt(m))];
        var g = new double[m, r];
        for (int l = 0; l < m; l++)
            for (int k = 0; k < r; k++) g[l, k] = map.Scales[l] * map.Basis(l, k);

        SampledTables[] at = [.. Enumerable.Range(0, m).Select(_ => new SampledTables(draws, subsets.Length, p))];
        var whole = new SampledTables(draws, kinds.Count, p + 1);
        int[] reported = [.. subsets.Select(atLocation.IndexOf)], reportedWhole = [.. kinds.Select(integrated.IndexOf)];
        int rows = n * segments.Count;
        double[][] w = [.. Enumerable.Range(0, r).Select(_ => new double[rows])];
        var a0 = new double[r];
        double[][] values = [.. segments.Select(_ => new double[n])];
        double[][] hybrids = [.. hybridAt.Select(s => values[s])];
        double[][]? reverses = atLocation.Reverses ? [.. reverseAt.Select(s => values[s])] : null;
        double[] wholeA = new double[n], wholeB = new double[n];
        double[][] wholeHybrids = [.. integrated.Blocks.Select(_ => new double[n])];
        double[][]? wholeReverses = integrated.Reverses ? [.. integrated.Blocks.Select(_ => new double[n])] : null;
        int maxTerms = Enumerable.Range(0, r).Max(k => selected.Select(j => pooled[k][joint.Draw(k, j)].Structure.TermIds.Length).DefaultIfEmpty(0).Max());
        var b = new double[maxTerms + 1];

        // Batches of joint draws whose components' distinct terms' columns over the stacked rows fit the column budget.
        foreach ((int start, int count) in Batches(selected, joint, pooled, rows))
        {
            BassSurface[] surfaces = [.. Enumerable.Range(0, r).Select(k => new BassSurface(basis, fit.Components[k].Terms, cases, null))];
            for (int d = start; d < start + count; d++)
            {
                for (int k = 0; k < r; k++)
                {
                    BassDraw draw = pooled[k][joint.Draw(k, selected[d])];
                    int[] ids = draw.Structure.TermIds;
                    Array.Copy(draw.Coefficients, b, ids.Length + 1);
                    a0[k] = b[0];
                    b[0] = 0.0;
                    surfaces[k].Evaluate(ids, b.AsSpan(0, ids.Length + 1), w[k]);
                }

                // Each location's surface Σₖ g_ℓk·wₖ, the intercepts left out as on the direct route.
                for (int l = 0; l < m; l++)
                {
                    for (int s = 0; s < segments.Count; s++)
                    {
                        double[] y = values[s];
                        Array.Clear(y);
                        for (int k = 0; k < r; k++)
                        {
                            double gk = g[l, k];
                            if (gk == 0.0) continue;
                            ReadOnlySpan<double> wk = w[k].AsSpan(s * n, n);
                            for (int j = 0; j < n; j++) y[j] += gk * wk[j];
                        }
                    }
                    at[l].Take(d, PickFreeze.Estimate(atLocation, values[0], values[1], hybrids, reverses, spec.RoundingTolerance), reported);
                }

                // The integrated form: y_m(ℓ) + Σₖ g_ℓk·(wₖ + aₖ₀) at each row's own location, B's where the block holds the location.
                void Fill(double[] y, int xSegment, int[] locations)
                {
                    for (int j = 0; j < n; j++)
                    {
                        int l = locations[j];
                        double v = map.Offsets[l];
                        for (int k = 0; k < r; k++) v += g[l, k] * (w[k][xSegment * n + j] + a0[k]);
                        y[j] = v;
                    }
                }
                Fill(wholeA, 0, locationA);
                Fill(wholeB, 1, locationB);
                for (int k = 0; k < integrated.Blocks.Count; k++)
                {
                    int[] u = integrated.Blocks[k];
                    bool location = u[^1] == p;
                    Fill(wholeHybrids[k], XSegment(u, false), location ? locationB : locationA);
                    if (wholeReverses is not null) Fill(wholeReverses[k], XSegment(u, true), location ? locationA : locationB);
                }
                whole.Take(d, PickFreeze.Estimate(integrated, wholeA, wholeB, wholeHybrids, wholeReverses, spec.RoundingTolerance), reportedWhole);
            }
        }

        int[] chainOf = [.. selected.Select(j => chains[j])], drawOf = [.. selected.Select(j => within[j])];
        SobolDecomposition[] locations = [.. at.Select(t => t.Decomposition(measure.Names, subsets, spec.RoundingTolerance, chainOf, drawOf))];
        SobolDecomposition wholeForm = whole.Decomposition([.. measure.Names, BasisSobolAnalysis.LocationInput], [.. kinds], spec.RoundingTolerance, chainOf, drawOf);
        return new BasisSobolAnalysis(map, measure, maxOrder, integratedOrder, subsets, locations, wholeForm, chainOf, drawOf, spec.RoundingTolerance, n);
    }

    // The joint draws in batches whose components' distinct terms over the stacked rows fit the column budget, one draw at least.
    private static IEnumerable<(int Start, int Count)> Batches(int[] selected, BassJointDraws joint, List<BassDraw>[] pooled, int rows)
    {
        int start = 0;
        var ids = new HashSet<(int, int)>();
        for (int d = 0; d < selected.Length; d++)
        {
            for (int k = 0; k < pooled.Length; k++)
                foreach (int id in pooled[k][joint.Draw(k, selected[d])].Structure.TermIds) ids.Add((k, id));
            if (d > start && (long)ids.Count * rows > BassSensitivity.ColumnBudget)
            {
                yield return (start, d - start);
                start = d;
                ids.Clear();
                for (int k = 0; k < pooled.Length; k++)
                    foreach (int id in pooled[k][joint.Draw(k, selected[d])].Structure.TermIds) ids.Add((k, id));
            }
        }
        if (selected.Length > start) yield return (start, selected.Length - start);
    }
}

    // One decomposition's tables on the sampled route, filled draw by draw from pick-freeze estimates.
internal sealed class SampledTables(int draws, int subsets, int inputs)
{
    private readonly double[,] _main = new double[draws, subsets], _mainError = new double[draws, subsets], _mainIndexError = new double[draws, subsets];
    private readonly double[,] _total = new double[draws, inputs], _totalError = new double[draws, inputs], _totalIndexError = new double[draws, inputs];
    private readonly double[] _variance = new double[draws], _varianceError = new double[draws];
    private readonly DecompositionStatus[] _status = new DecompositionStatus[draws];

    public void Take(int d, PickFreezeEstimate e, int[] reported)
    {
        _variance[d] = e.Variance;
        _varianceError[d] = e.VarianceError;
        for (int k = 0; k < reported.Length; k++)
        {
            _main[d, k] = e.Main[reported[k]];
            _mainError[d, k] = e.MainError[reported[k]];
            _mainIndexError[d, k] = e.MainIndexError[reported[k]];
        }
        for (int v = 0; v < e.Total.Length; v++)
        {
            _total[d, v] = e.Total[v];
            _totalError[d, v] = e.TotalError[v];
            _totalIndexError[d, v] = e.TotalIndexError[v];
        }
        _status[d] = ((e.Status & PickFreezeStatus.ZeroVariance) != 0 ? DecompositionStatus.ZeroVariance : DecompositionStatus.Defined)
                     | ((e.Status & PickFreezeStatus.BelowZeroBeyondError) != 0 ? DecompositionStatus.BelowZeroBeyondError : DecompositionStatus.Defined);
    }

    public SobolDecomposition Decomposition(IReadOnlyList<string> names, int[][] subsets, double tolerance, int[] chains, int[] draws, double[]? weights = null) =>
        new(names, subsets, _main, _total, _variance, _status, tolerance, chains, draws, weights)
        {
            VarianceError = _varianceError, MainVarianceError = _mainError, TotalVarianceError = _totalError,
            MainIndexError = _mainIndexError, TotalIndexError = _totalIndexError,
        };
}
