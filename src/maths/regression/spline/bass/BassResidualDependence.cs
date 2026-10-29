using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass;

/// <summary>What a residual-dependence check reads per case.</summary>
public enum ResidualStatistic
{
    /// <summary>The lag-one autocorrelation along an order of the locations.</summary>
    LagOne,

    /// <summary>Moran's I over each location's nearest neighbours in output coordinates.</summary>
    MoransI,
}

/// <summary>
/// The residual-dependence check across output locations (F11): a statistic per case of its residual, and a summary
/// across cases. A model whose noise is independent across locations (D07) leaves residuals reading about
/// <see cref="Expected"/>; a misfit smooth across the locations, as a deterministic simulator's is, reads high, and then
/// the fit counts n·m independent observations where it has fewer, so its intervals come out too narrow.
/// <list type="bullet">
/// <item><see cref="LagOne"/>: along an order of the locations, r₁ = Σₖ(eₖ − ē)(eₖ₊₁ − ē) / Σₖ(eₖ − ē)², with ē the case's
/// mean residual; about −1/m under independence, with a standard deviation about 1/√m over m locations.</item>
/// <item><see cref="MoransI"/>: over each location's K nearest neighbours in output coordinates (binary weights),
/// I = Σₖ Σ_{l ∈ N(k)} (eₖ − ē)(e_l − ē) / (K·Σₖ(eₖ − ē)²); −1/(m − 1) in expectation under independence (Moran 1950;
/// Cliff and Ord 1981).</item>
/// </list>
/// </summary>
/// <param name="Statistic">The statistic read per case.</param>
/// <param name="PerCase">Each case's value; NaN for a case whose residual does not vary.</param>
/// <param name="Mean">The mean of <paramref name="PerCase"/> over the cases where it is defined; NaN when none is.</param>
/// <param name="Expected">The statistic's mean under independent noise: −1/m for r₁, to first order, and −1/(m − 1) for
/// Moran's I.</param>
/// <param name="Probabilities">The probabilities of <paramref name="Quantiles"/>.</param>
/// <param name="Quantiles">The type-7 quantiles of <paramref name="PerCase"/> over the cases where it is defined; NaN when
/// none is.</param>
public sealed record ResidualDependence(ResidualStatistic Statistic, double[] PerCase, double Mean, double Expected, double[] Probabilities,
                                        double[] Quantiles)
{
    /// <summary>The summary's default probabilities: 5%, 25%, 50%, 75% and 95%.</summary>
    public static IReadOnlyList<double> DefaultProbabilities { get; } = Array.AsReadOnly(new[] { 0.05, 0.25, 0.5, 0.75, 0.95 });

    /// <summary>
    /// r₁ of the cases' residuals <paramref name="residuals"/> (cases × locations) along the locations ordered by
    /// <paramref name="index"/>, ties kept in location order, summarized at <paramref name="probabilities"/> (null:
    /// <see cref="DefaultProbabilities"/>). Refused: no case, fewer than two locations, a non-finite residual, an index of
    /// another length or with a non-finite value, and a probability outside [0, 1].
    /// </summary>
    public static ResidualDependence LagOne(double[,] residuals, ReadOnlySpan<double> index, IReadOnlyList<double>? probabilities = null)
    {
        (int n, int m) = CheckedResiduals(residuals);
        if (index.Length != m) throw new ArgumentException("The index needs one value per location.", nameof(index));
        foreach (double t in index)
            if (!double.IsFinite(t)) throw new ArgumentException("The index has a non-finite value.", nameof(index));
        double[] p = CheckedProbabilities(probabilities);

        double[] keys = index.ToArray();
        var order = new int[m];
        for (int r = 0; r < m; r++) order[r] = r;
        Array.Sort(order, (a, b) => keys[a] != keys[b] ? keys[a].CompareTo(keys[b]) : a.CompareTo(b));

        var lagOne = new double[n];
        var e = new double[m];
        for (int i = 0; i < n; i++)
        {
            double mean = 0.0;
            for (int k = 0; k < m; k++) mean += e[k] = residuals[i, order[k]];
            mean /= m;
            double cross = 0.0, squares = 0.0;
            for (int k = 0; k < m; k++)
            {
                double d = e[k] - mean;
                squares += d * d;
                if (k + 1 < m) cross += d * (e[k + 1] - mean);
            }
            lagOne[i] = squares > 0.0 ? cross / squares : double.NaN;
        }
        return Summarize(ResidualStatistic.LagOne, lagOne, -1.0 / m, p);
    }

    /// <summary>
    /// Moran's I of the cases' residuals <paramref name="residuals"/> (cases × locations) over each location's
    /// <paramref name="neighbours"/> nearest others by Euclidean distance in <paramref name="coordinates"/> (one array of
    /// the locations' values per coordinate), ties to the lower location, summarized at <paramref name="probabilities"/>
    /// (null: <see cref="DefaultProbabilities"/>). Refused: no case, fewer than two locations, a non-finite residual or
    /// coordinate, no coordinate or one of another length, a neighbour count outside [1, m − 1], and a probability
    /// outside [0, 1].
    /// </summary>
    public static ResidualDependence MoransI(double[,] residuals, IReadOnlyList<double[]> coordinates, int neighbours = 4,
                                             IReadOnlyList<double>? probabilities = null)
    {
        (int n, int m) = CheckedResiduals(residuals);
        ArgumentNullException.ThrowIfNull(coordinates);
        if (coordinates.Count == 0) throw new ArgumentException("Moran's I needs a coordinate.", nameof(coordinates));
        foreach (double[] c in coordinates)
        {
            if (c is null || c.Length != m) throw new ArgumentException("Every coordinate needs one value per location.", nameof(coordinates));
            foreach (double v in c)
                if (!double.IsFinite(v)) throw new ArgumentException("A coordinate has a non-finite value.", nameof(coordinates));
        }
        if (neighbours < 1 || neighbours > m - 1)
            throw new ArgumentOutOfRangeException(nameof(neighbours), $"The neighbour count lies in [1, m − 1] = [1, {m - 1}].");
        double[] p = CheckedProbabilities(probabilities);

        // Each location's nearest others, once for every case.
        var near = new int[m][];
        var distance = new double[m];
        var others = new int[m - 1];
        for (int k = 0; k < m; k++)
        {
            for (int l = 0; l < m; l++)
            {
                double s = 0.0;
                foreach (double[] c in coordinates)
                {
                    double d = c[k] - c[l];
                    s += d * d;
                }
                distance[l] = s;
            }
            for (int l = 0, j = 0; l < m; l++)
                if (l != k) others[j++] = l;
            Array.Sort(others, (a, b) => distance[a] != distance[b] ? distance[a].CompareTo(distance[b]) : a.CompareTo(b));
            near[k] = others[..neighbours];
        }

        var values = new double[n];
        var d0 = new double[m];
        for (int i = 0; i < n; i++)
        {
            double mean = 0.0;
            for (int k = 0; k < m; k++) mean += residuals[i, k];
            mean /= m;
            double squares = 0.0;
            for (int k = 0; k < m; k++)
            {
                d0[k] = residuals[i, k] - mean;
                squares += d0[k] * d0[k];
            }
            double cross = 0.0;
            for (int k = 0; k < m; k++)
            {
                double around = 0.0;
                foreach (int l in near[k]) around += d0[l];
                cross += d0[k] * around;
            }
            values[i] = squares > 0.0 ? cross / (neighbours * squares) : double.NaN;
        }
        return Summarize(ResidualStatistic.MoransI, values, -1.0 / (m - 1), p);
    }

    /// <summary>A copy of <paramref name="probabilities"/> (null: <see cref="DefaultProbabilities"/>), each checked to lie
    /// in [0, 1].</summary>
    internal static double[] CheckedProbabilities(IReadOnlyList<double>? probabilities)
    {
        double[] p = [.. probabilities ?? DefaultProbabilities];
        foreach (double q in p)
            if (!(q >= 0.0 && q <= 1.0)) throw new ArgumentOutOfRangeException(nameof(probabilities), "A probability lies in [0, 1].");
        return p;
    }

    private static (int Cases, int Locations) CheckedResiduals(double[,] residuals)
    {
        ArgumentNullException.ThrowIfNull(residuals);
        int n = residuals.GetLength(0), m = residuals.GetLength(1);
        if (n < 1) throw new ArgumentException("The check needs a case.", nameof(residuals));
        if (m < 2) throw new ArgumentException("The check needs at least two locations.", nameof(residuals));
        foreach (double v in residuals)
            if (!double.IsFinite(v)) throw new ArgumentException("A residual is not finite.", nameof(residuals));
        return (n, m);
    }

    // The record over the per-case values, its summary over the cases where a value is defined.
    private static ResidualDependence Summarize(ResidualStatistic statistic, double[] values, double expected, double[] p)
    {
        var defined = new List<double>(values.Length);
        foreach (double v in values)
            if (!double.IsNaN(v)) defined.Add(v);
        double[] sorted = [.. defined];
        Array.Sort(sorted);
        var quantiles = new double[p.Length];
        double average = double.NaN;
        if (sorted.Length > 0)
        {
            average = 0.0;
            foreach (double v in sorted) average += v;
            average /= sorted.Length;
            for (int q = 0; q < p.Length; q++) quantiles[q] = PosteriorInterval.Quantile(sorted, p[q]);
        }
        else Array.Fill(quantiles, double.NaN);
        return new ResidualDependence(statistic, values, average, expected, p, quantiles);
    }
}

/// <summary>
/// The output index a residual-dependence check reads a direct functional fit's residuals along, declared on the call
/// (decided 2026-10-28). A fit with one output variable reads along it by default; with several, no order is decided
/// until one is declared, since an output variable may be a covariate known at each location and not a coordinate.
/// </summary>
public abstract record ResidualIndex
{
    private protected ResidualIndex()
    {
    }

    /// <summary>r₁ along output variable <paramref name="Variable"/>'s order of the locations, ties in location order.</summary>
    public sealed record AlongVariable(int Variable) : ResidualIndex;

    /// <summary>r₁ along an explicit order: <paramref name="Locations"/>[k] is the k-th location read, a permutation of
    /// the locations.</summary>
    public sealed record AlongOrder(IReadOnlyList<int> Locations) : ResidualIndex;

    /// <summary>Moran's I over each location's <paramref name="Neighbours"/> nearest others in the scaled coordinates of
    /// output variables <paramref name="Variables"/>, for d ≥ 2 coordinates.</summary>
    public sealed record NearestNeighbours(IReadOnlyList<int> Variables, int Neighbours = 4) : ResidualIndex;
}

/// <summary>
/// BASS's residual-dependence check (F11; the direct functional algebra's proposal) over retained draws: each case's
/// residual y − f̄ under the posterior mean f̄ read along an output index by <see cref="ResidualDependence"/>. The fit's
/// <c>ResidualDependence</c> field is its default, a view over the draws; <c>BassFit.ResidualDependenceFor</c> takes a
/// declared index. The posterior mean goes by term: each distinct term's column is evaluated once, against its
/// coefficients summed over the draws, in O(S·M + T·N) for S draws of M terms over T distinct terms and N observations.
/// </summary>
public static class BassResidualDependence
{
    /// <summary>
    /// The check for a direct functional fit along <paramref name="index"/> (null: along its one output variable): the
    /// retained draws <paramref name="draws"/> of <paramref name="model"/>, per chain, whose structures index
    /// <paramref name="terms"/>, their fits averaged with <paramref name="weights"/> over the draws pooled in chain order
    /// (null: equal, as for MCMC), summarized at <paramref name="probabilities"/>. Null where no index is declared and none
    /// is decided: a scalar response, and several output variables. Refused: a declared index on a scalar response, and an
    /// index the outputs do not have.
    /// </summary>
    public static ResidualDependence? Of(BassModel model, BassTermTable terms, IReadOnlyList<IReadOnlyList<BassDraw>> draws,
                                         ResidualIndex? index = null, IReadOnlyList<double>? weights = null,
                                         IReadOnlyList<double>? probabilities = null)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(terms);
        (List<BassDraw> pooled, _) = BassDrawFit.Pool(draws, terms, false, nameof(draws));
        ResidualDependence.CheckedProbabilities(probabilities);
        if (weights is not null)
        {
            if (weights.Count != pooled.Count) throw new ArgumentException("One weight per draw is required.", nameof(weights));
            double total = 0.0;
            foreach (double w in weights)
            {
                if (!(w >= 0.0) || double.IsPositiveInfinity(w)) throw new ArgumentException("A weight must be non-negative and finite.", nameof(weights));
                total += w;
            }
            if (!(total > 0.0)) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
        }
        OutputInputs? outputs = model.Inputs.Outputs;
        if (outputs is null)
        {
            if (index is not null) throw new ArgumentException("A scalar response has no output index.", nameof(index));
            return null;
        }
        if (index is null)
        {
            if (outputs.Count != 1) return null;
            index = new ResidualIndex.AlongVariable(0);
        }
        Func<double[,], ResidualDependence> read = Reader(index, outputs, probabilities);

        BassResponse response = model.Data.Response;
        ReadOnlySpan<double> y = response.Values;
        int n = response.Cases, m = response.Locations;
        double[] mean = PosteriorMean(new BassDrawFit(model, terms), pooled, weights, n, m);
        var residuals = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int r = 0; r < m; r++) residuals[i, r] = y[i * m + r] - mean[i * m + r];
        return read(residuals);
    }

    // The declared index, checked against the outputs, as the check it makes of the residuals.
    private static Func<double[,], ResidualDependence> Reader(ResidualIndex index, OutputInputs outputs, IReadOnlyList<double>? probabilities)
    {
        int m = outputs.Locations;
        switch (index)
        {
            case ResidualIndex.AlongVariable along:
                if (along.Variable < 0 || along.Variable >= outputs.Count)
                    throw new ArgumentException($"The outputs have {outputs.Count} variable(s), not variable {along.Variable}.", nameof(index));
                double[] values = outputs.Variable(along.Variable).ToArray();
                return e => ResidualDependence.LagOne(e, values, probabilities);
            case ResidualIndex.AlongOrder order:
                if (order.Locations is not { } locations || locations.Count != m)
                    throw new ArgumentException($"The order must name each of the {m} locations once.", nameof(index));
                var keys = new double[m];
                var seen = new bool[m];
                for (int k = 0; k < m; k++)
                {
                    int r = locations[k];
                    if (r < 0 || r >= m || seen[r]) throw new ArgumentException($"The order must name each of the {m} locations once.", nameof(index));
                    seen[r] = true;
                    keys[r] = k;
                }
                return e => ResidualDependence.LagOne(e, keys, probabilities);
            case ResidualIndex.NearestNeighbours near:
                if (near.Variables is not { Count: > 0 } variables)
                    throw new ArgumentException("Moran's I needs an output variable.", nameof(index));
                var coordinates = new double[variables.Count][];
                var used = new HashSet<int>();
                for (int v = 0; v < variables.Count; v++)
                {
                    int variable = variables[v];
                    if (variable < 0 || variable >= outputs.Count || !used.Add(variable))
                        throw new ArgumentException($"The outputs have {outputs.Count} variable(s); each is named once.", nameof(index));
                    coordinates[v] = outputs.Variable(variable).ToArray();
                }
                if (near.Neighbours < 1 || near.Neighbours > m - 1)
                    throw new ArgumentException($"The neighbour count lies in [1, m − 1] = [1, {m - 1}].", nameof(index));
                return e => ResidualDependence.MoransI(e, coordinates, near.Neighbours, probabilities);
            default:
                throw new ArgumentException("Not a known output index.", nameof(index));
        }
    }

    // f̄ at every observation by term: the intercepts' and each term's coefficients summed over the draws with their
    // weights, then each distinct term's column once.
    private static double[] PosteriorMean(BassDrawFit fit, List<BassDraw> pooled, IReadOnlyList<double>? weights, int n, int m)
    {
        double intercept = 0.0, total = 0.0;
        var coefficients = new Dictionary<int, double>();
        for (int s = 0; s < pooled.Count; s++)
        {
            double w = weights?[s] ?? 1.0;
            if (w == 0.0) continue;
            BassDraw draw = pooled[s];
            intercept += w * draw.Coefficients[0];
            int[] ids = draw.Structure.TermIds;
            for (int q = 0; q < ids.Length; q++)
                coefficients[ids[q]] = coefficients.GetValueOrDefault(ids[q]) + w * draw.Coefficients[q + 1];
            total += w;
        }
        var mean = new double[n * m];
        Array.Fill(mean, intercept / total);
        foreach ((int id, double sum) in coefficients)
        {
            (double[] a, double[]? h) = fit.Columns(id);
            double c = sum / total;
            for (int i = 0; i < n; i++)
            {
                double ca = c * a[i];
                if (ca == 0.0) continue;
                if (h is null)
                {
                    mean[i] += ca;
                    continue;
                }
                for (int r = 0; r < m; r++) mean[i * m + r] += ca * h[r];
            }
        }
        return mean;
    }
}
