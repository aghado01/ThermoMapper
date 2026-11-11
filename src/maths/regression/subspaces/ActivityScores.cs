using System;
using System.Collections.Generic;

namespace Maths.Regression.Subspaces;

/// <summary>
/// Activity scores (E8, slice 6): each input's share of the gradient energy, in the metric C was read in, per draw and under the
/// draws' weights.
/// <list type="bullet">
/// <item><b>The diagonal,</b> α_i = C_ii = Σ_{j=1}^{p} λ_j w_ij², the sum over every eigenpair: 2307.11241's score (slice 0 read it).
/// Always reported. Its weighted mean is E[C]'s diagonal, and the scores sum to tr C.</item>
/// <item><b>The truncation at k,</b> α_i(k) = Σ_{j≤k} λ_j w_ij², which reads only the directions kept (Constantine and Diaz, cited from
/// memory): at the readout's k per draw and for E[C]'s own eigensystem, and under a rank rule at each draw's k_s. At k = p it is the
/// diagonal.</item>
/// </list>
/// </summary>
/// <param name="Draws">Each draw's diagonal, draws × inputs.</param>
/// <param name="Mean">The diagonal's weighted mean, E[C]'s diagonal.</param>
/// <param name="Dimension">k.</param>
/// <param name="Truncated">Each draw's scores truncated at k, draws × inputs.</param>
/// <param name="TruncatedMean">Their weighted mean.</param>
/// <param name="MeanCTruncated">E[C]'s own scores truncated at k, from its eigensystem.</param>
/// <param name="RankTruncated">Each draw's scores truncated at its own k_s; null without a rank rule.</param>
/// <param name="RankTruncatedMean">Their weighted mean; null without a rank rule.</param>
public sealed record ActivityScores(double[,] Draws, double[] Mean, int Dimension, double[,] Truncated, double[] TruncatedMean,
                                    double[] MeanCTruncated, double[,]? RankTruncated, double[]? RankTruncatedMean)
{
    /// <summary>The scores of <paramref name="posterior"/>'s draws, and at each draw's own rank when <paramref name="ranks"/> is given.</summary>
    public static ActivityScores Of(SubspacePosterior posterior, RankPosterior? ranks = null)
    {
        ArgumentNullException.ThrowIfNull(posterior);
        IReadOnlyList<DrawSubspace> draws = posterior.Draws;
        IReadOnlyList<double> w = posterior.Weights;
        int n = draws.Count, p = posterior.MeanC.GetLength(0), k = posterior.Dimension;
        if (ranks is not null && ranks.Ranks.Count != n) throw new ArgumentException("One rank per draw.", nameof(ranks));
        var diagonal = new double[n, p];
        var truncated = new double[n, p];
        double[,]? atRank = ranks is null ? null : new double[n, p];
        var mean = new double[p];
        var truncatedMean = new double[p];
        double[]? rankMean = ranks is null ? null : new double[p];
        for (int s = 0; s < n; s++)
        {
            DrawSubspace d = draws[s];
            for (int i = 0; i < p; i++)
            {
                diagonal[s, i] = d.C[i, i];
                truncated[s, i] = Truncate(d.Eigenvalues, d.Vectors, i, k);
                mean[i] += w[s] * diagonal[s, i];
                truncatedMean[i] += w[s] * truncated[s, i];
                if (atRank is null) continue;
                atRank[s, i] = Truncate(d.Eigenvalues, d.Vectors, i, ranks!.Ranks[s]);
                rankMean![i] += w[s] * atRank[s, i];
            }
        }
        var meanC = new double[p];
        for (int i = 0; i < p; i++) meanC[i] = Truncate(posterior.MeanCEigenvalues, posterior.MeanCVectors, i, k);
        return new ActivityScores(diagonal, mean, k, truncated, truncatedMean, meanC, atRank, rankMean);
    }

    // Σ_{j<k} λ_j w_ij² over the leading eigenpairs; vectors beyond the factor's rank carry λ = 0 and are absent.
    private static double Truncate(IReadOnlyList<double> eigenvalues, double[,] vectors, int i, int k)
    {
        int kept = Math.Min(k, vectors.GetLength(1));
        double s = 0.0;
        for (int j = 0; j < kept; j++) s += eigenvalues[j] * vectors[i, j] * vectors[i, j];
        return s;
    }
}
