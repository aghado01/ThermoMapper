using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// The residual-dependence check along an output index (F11): each case's lag-one autocorrelation of its residual along
/// the locations ordered by the index, r₁ = Σₖ(eₖ − ē)(eₖ₊₁ − ē) / Σₖ(eₖ − ē)², with ē the case's mean residual, and a
/// summary across cases. A model whose noise is independent across locations (D07) leaves residuals with r₁ about −1/m
/// and a standard deviation about 1/√m over m locations; a misfit smooth along the index, as a deterministic simulator's
/// is, reads high, and then the fit counts n·m independent observations where it has fewer, so its intervals come out
/// too narrow.
/// </summary>
/// <param name="LagOne">Each case's r₁; NaN for a case whose residual does not vary.</param>
/// <param name="Mean">The mean of <paramref name="LagOne"/> over the cases where it is defined; NaN when none is.</param>
/// <param name="Probabilities">The probabilities of <paramref name="Quantiles"/>.</param>
/// <param name="Quantiles">The type-7 quantiles of <paramref name="LagOne"/> over the cases where it is defined; NaN when
/// none is.</param>
public sealed record ResidualDependence(double[] LagOne, double Mean, double[] Probabilities, double[] Quantiles)
{
    /// <summary>The summary's default probabilities: 5%, 25%, 50%, 75% and 95%.</summary>
    public static IReadOnlyList<double> DefaultProbabilities { get; } = Array.AsReadOnly(new[] { 0.05, 0.25, 0.5, 0.75, 0.95 });

    /// <summary>
    /// The check of the cases' residuals <paramref name="residuals"/> (cases × locations) along the locations ordered by
    /// <paramref name="index"/>, ties kept in location order, summarized at <paramref name="probabilities"/> (null:
    /// <see cref="DefaultProbabilities"/>). Refused: no case, fewer than two locations, a non-finite residual, an index of
    /// another length or with a non-finite value, and a probability outside [0, 1].
    /// </summary>
    public static ResidualDependence Of(double[,] residuals, ReadOnlySpan<double> index, IReadOnlyList<double>? probabilities = null)
    {
        ArgumentNullException.ThrowIfNull(residuals);
        int n = residuals.GetLength(0), m = residuals.GetLength(1);
        if (n < 1) throw new ArgumentException("The check needs a case.", nameof(residuals));
        if (m < 2) throw new ArgumentException("A lag needs at least two locations.", nameof(residuals));
        foreach (double v in residuals)
            if (!double.IsFinite(v)) throw new ArgumentException("A residual is not finite.", nameof(residuals));
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
        var defined = new List<double>(n);
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
            if (squares > 0.0) defined.Add(lagOne[i]);
        }

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
        return new ResidualDependence(lagOne, average, p, quantiles);
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
}

/// <summary>
/// BASS's residual-dependence check (F11; the direct functional algebra's proposal) over retained draws: each case's
/// residual y − f̄ under the posterior mean f̄, the draws' fits averaged, read along the output index by
/// <see cref="ResidualDependence"/>. It is the fit's <c>ResidualDependence</c> field, a view over the draws.
/// </summary>
public static class BassResidualDependence
{
    /// <summary>
    /// The check for a direct functional fit with one output variable, along its locations in that variable's order: the
    /// retained draws <paramref name="draws"/> of <paramref name="model"/>, per chain, whose structures index
    /// <paramref name="terms"/>, their fits averaged with <paramref name="weights"/> over the draws pooled in chain order
    /// (null: equal, as for MCMC), summarized at <paramref name="probabilities"/>. Null where no order is decided: a
    /// scalar response, and several output variables.
    /// </summary>
    public static ResidualDependence? Of(BassModel model, BassTermTable terms, IReadOnlyList<IReadOnlyList<BassDraw>> draws,
                                         IReadOnlyList<double>? weights = null, IReadOnlyList<double>? probabilities = null)
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
        if (model.Inputs.Outputs is not { Count: 1 } outputs) return null;

        BassResponse response = model.Data.Response;
        ReadOnlySpan<double> y = response.Values;
        int n = response.Cases, m = response.Locations;
        var fit = new BassDrawFit(model, terms);
        var mean = new double[y.Length];
        var f = new double[y.Length];
        double sum = 0.0;
        for (int s = 0; s < pooled.Count; s++)
        {
            double w = weights?[s] ?? 1.0;
            if (w == 0.0) continue;
            fit.Evaluate(pooled[s], f);
            for (int k = 0; k < f.Length; k++) mean[k] += w * f[k];
            sum += w;
        }
        var residuals = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int r = 0; r < m; r++) residuals[i, r] = y[i * m + r] - mean[i * m + r] / sum;
        return ResidualDependence.Of(residuals, outputs.Variable(0), probabilities);
    }
}
