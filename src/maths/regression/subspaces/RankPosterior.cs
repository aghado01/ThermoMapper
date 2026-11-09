using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;

namespace Maths.Regression.Subspaces;

/// <summary>Which rule reads a draw's active dimension k_s from its eigenvalues (E8, slice 5).</summary>
public enum RankRuleKind
{
    /// <summary>The smallest k whose leading eigenvalues carry at least <see cref="RankRule.Share"/> of tr C.</summary>
    TraceShare,

    /// <summary>The k up to <see cref="RankRule.MaxRank"/> with the largest gap relative to the leading eigenvalue, (λ_k − λ_{k+1})/λ₁.</summary>
    Eigengap,
}

/// <summary>
/// A declared rule for each draw's active dimension k_s (E8, slice 5; the Grassmann discussion's item 1, settled by the user on
/// 2026-11-08). There is no default: a readout reads per-draw ranks only when given a rule. A draw whose C is zero has k_s = 0, the
/// empty subspace. The gap rule measures each gap against λ₁, with λ_{p+1} = 0, over k from 1 to min(k_max, p): a gap against λ_k
/// itself, or on a log scale, would read any trailing zero eigenvalue (an input no term holds) as the largest gap. Against λ₁ it favours
/// the first gap when the spectrum decays geometrically; the trace share is the other rule. A declarative record, JSON-serializable.
/// </summary>
public sealed record RankRule
{
    /// <summary>The rule.</summary>
    public RankRuleKind Kind { get; init; }

    /// <summary>The trace share, in (0, 1], for <see cref="RankRuleKind.TraceShare"/>.</summary>
    public double Share { get; init; }

    /// <summary>k_max, at least 1, for <see cref="RankRuleKind.Eigengap"/>.</summary>
    public int MaxRank { get; init; }

    /// <summary>The smallest k whose leading eigenvalues carry <paramref name="share"/> of tr C.</summary>
    public static RankRule TraceShare(double share) => new() { Kind = RankRuleKind.TraceShare, Share = share };

    /// <summary>The k up to <paramref name="maxRank"/> with the largest relative eigengap.</summary>
    public static RankRule Eigengap(int maxRank) => new() { Kind = RankRuleKind.Eigengap, MaxRank = maxRank };

    /// <summary>Refused: a share outside (0, 1], a k_max below 1, or an unknown kind.</summary>
    public void Validate()
    {
        switch (Kind)
        {
            case RankRuleKind.TraceShare when !(Share > 0.0 && Share <= 1.0):
                throw new ArgumentOutOfRangeException(nameof(Share), "A trace share is in (0, 1].");
            case RankRuleKind.Eigengap when MaxRank < 1:
                throw new ArgumentOutOfRangeException(nameof(MaxRank), "k_max is at least 1.");
            case RankRuleKind.TraceShare or RankRuleKind.Eigengap:
                return;
            default:
                throw new ArgumentOutOfRangeException(nameof(Kind), $"No rank rule {Kind}.");
        }
    }

    /// <summary>k_s for the descending <paramref name="eigenvalues"/> of one draw's C, negative ones read as 0.</summary>
    public int Apply(IReadOnlyList<double> eigenvalues)
    {
        ArgumentNullException.ThrowIfNull(eigenvalues);
        Validate();
        int p = eigenvalues.Count;
        double[] l = [.. eigenvalues.Select(x => Math.Max(0.0, x))];
        double trace = l.Sum();
        if (!(trace > 0.0) || !(l[0] > 0.0)) return 0;
        if (Kind == RankRuleKind.TraceShare)
        {
            double cumulative = 0.0;
            for (int k = 1; k <= p; k++)
            {
                cumulative += l[k - 1];
                if (cumulative >= Share * trace) return k;
            }
            return p;
        }
        int best = 1;
        double bestGap = double.NegativeInfinity;
        for (int k = 1; k <= Math.Min(MaxRank, p); k++)
        {
            double gap = (l[k - 1] - (k < p ? l[k] : 0.0)) / l[0];
            if (gap > bestGap)
            {
                bestGap = gap;
                best = k;
            }
        }
        return best;
    }
}

/// <summary>
/// A posterior over active subspaces of mixed rank (E8, slice 5): each draw at its own k_s by a declared <see cref="RankRule"/>, beside
/// the fixed-k readout.
/// <list type="bullet">
/// <item><b>The ranks:</b> each k_s and their weighted distribution over 0 to p.</item>
/// <item><b>P̄ = Σ w_s P_s</b> keeps its formula over mixed ranks; tr P̄ is the mean rank, and its ordered eigenvectors are the flag mean
/// for subspaces of mixed dimensions (Draper et al. 2014, cited from memory).</item>
/// <item><b>The rank summary.</b> Over projections P of every rank, Σ_s w_s‖P − P_s‖²_F = tr P + E[k_s] − 2 tr(P P̄), which for P = UUᵀ
/// is E[k_s] + Σ_i (1 − 2 u_iᵀP̄u_i). By Ky Fan the best P of each rank takes P̄'s leading eigenvectors, and the best rank keeps those
/// whose eigenvalue exceeds ½. On coordinate subspaces that is the median probability model (Barbieri and Berger 2004, cited from
/// memory). It reads the ranks' posterior as one subspace without choosing k, and is reported, not used as a default.</item>
/// <item><b>Each draw's distance to the summary,</b> the projection distance across ranks, d² = ½‖P_s − P*‖²_F = |k_s − k*|/2 +
/// Σ sin²θ_i over the min(k_s, k*) principal angles (Ye and Lim 2016, cited from memory): a rank term and a containment defect,
/// reported apart.</item>
/// </list>
/// </summary>
public sealed class RankPosterior
{
    private RankPosterior(RankRule rule, int[] ranks, double[] weights, double[] distribution, double[,] projection, double[] values,
                          double[,] vectors, int summaryRank, double summaryLoss, double[] rankTerms, double[] containment)
    {
        Rule = rule;
        Ranks = ranks;
        Weights = weights;
        Distribution = distribution;
        Projection = projection;
        ProjectionEigenvalues = values;
        ProjectionVectors = vectors;
        SummaryRank = summaryRank;
        SummaryLoss = summaryLoss;
        RankTerms = rankTerms;
        Containment = containment;
    }

    /// <summary>The rule the ranks were read by.</summary>
    public RankRule Rule { get; }

    /// <summary>Each draw's k_s.</summary>
    public IReadOnlyList<int> Ranks { get; }

    /// <summary>The draws' weights, normalized.</summary>
    public IReadOnlyList<double> Weights { get; }

    /// <summary>The weighted share of the draws at each rank from 0 to p.</summary>
    public IReadOnlyList<double> Distribution { get; }

    /// <summary>E[k_s], which is tr P̄.</summary>
    public double MeanRank => Ranks.Select((k, s) => k * Weights[s]).Sum();

    /// <summary>P̄ = Σ w_s P_s over the draws' own ranks.</summary>
    public double[,] Projection { get; }

    /// <summary>P̄'s eigenvalues, descending, in [0, 1].</summary>
    public double[] ProjectionEigenvalues { get; }

    /// <summary>P̄'s eigenvectors as columns, in order.</summary>
    public double[,] ProjectionVectors { get; }

    /// <summary>k*, the number of P̄'s eigenvalues above ½: the rank summary's dimension.</summary>
    public int SummaryRank { get; }

    /// <summary>The rank summary: P̄'s leading k* eigenvectors, p × k*.</summary>
    public double[,] Summary => Leading(ProjectionVectors, SummaryRank);

    /// <summary>Σ_s w_s‖P* − P_s‖²_F at the summary, the least over projections of every rank.</summary>
    public double SummaryLoss { get; }

    /// <summary>Each draw's rank term, |k_s − k*|/2.</summary>
    public IReadOnlyList<double> RankTerms { get; }

    /// <summary>Each draw's containment defect against the summary, Σ sin²θ_i over the min(k_s, k*) principal angles.</summary>
    public IReadOnlyList<double> Containment { get; }

    /// <summary>Σ_s w_s‖P − P_s‖²_F = tr P + E[k_s] − 2 tr(P P̄) for any orthonormal basis <paramref name="basis"/> (p × k, k ≥ 0).</summary>
    public double Loss(double[,] basis)
    {
        ArgumentNullException.ThrowIfNull(basis);
        int p = Projection.GetLength(0), k = basis.GetLength(1);
        if (basis.GetLength(0) != p) throw new ArgumentException("The basis has one row per input.", nameof(basis));
        double cross = 0.0;
        for (int j = 0; j < k; j++)
            for (int a = 0; a < p; a++)
                for (int b = 0; b < p; b++) cross += basis[a, j] * Projection[a, b] * basis[b, j];
        return k + MeanRank - 2.0 * cross;
    }

    /// <summary>
    /// The mixed-rank posterior of <paramref name="draws"/> under <paramref name="weights"/> (nonnegative, normalized here), each draw at
    /// the rank <paramref name="rule"/> reads from its eigenvalues, its subspace its leading k_s eigenvectors.
    /// </summary>
    public static RankPosterior Of(IReadOnlyList<DrawSubspace> draws, IReadOnlyList<double> weights, RankRule rule)
    {
        ArgumentNullException.ThrowIfNull(draws);
        ArgumentNullException.ThrowIfNull(weights);
        ArgumentNullException.ThrowIfNull(rule);
        rule.Validate();
        if (draws.Count == 0) throw new ArgumentException("A posterior needs a draw.", nameof(draws));
        if (weights.Count != draws.Count) throw new ArgumentException("One weight per draw.", nameof(weights));
        double[] w = Whitening.Normalized(weights, draws.Count, nameof(weights));
        int p = draws[0].C.GetLength(0);
        var ranks = new int[draws.Count];
        var bases = new double[draws.Count][,];
        var distribution = new double[p + 1];
        var projection = new double[p, p];
        for (int s = 0; s < draws.Count; s++)
        {
            DrawSubspace d = draws[s];
            if (d.C.GetLength(0) != p) throw new ArgumentException("The draws share their inputs.", nameof(draws));
            int k = ranks[s] = Math.Min(rule.Apply(d.Eigenvalues), d.Vectors.GetLength(1));
            bases[s] = Leading(d.Vectors, k);
            distribution[k] += w[s];
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                {
                    double pij = 0.0;
                    for (int r = 0; r < k; r++) pij += bases[s][i, r] * bases[s][j, r];
                    projection[i, j] += w[s] * pij;
                }
        }
        EigenResult e = ScaledEigen.Of(projection);
        var vectors = new double[p, p];
        for (int j = 0; j < p; j++)
            for (int i = 0; i < p; i++) vectors[i, j] = e.Eigenvectors[j][i];
        double[] values = [.. e.Eigenvalues];
        int summaryRank = values.Count(x => x > 0.5);
        double[,] summary = Leading(vectors, summaryRank);
        var rankTerms = new double[draws.Count];
        var containment = new double[draws.Count];
        for (int s = 0; s < draws.Count; s++)
        {
            rankTerms[s] = Math.Abs(ranks[s] - summaryRank) / 2.0;
            if (ranks[s] == 0 || summaryRank == 0) continue;
            foreach (double theta in SubspacePosterior.PrincipalAngles(bases[s], summary)) containment[s] += Math.Sin(theta) * Math.Sin(theta);
        }
        var posterior = new RankPosterior(rule, ranks, w, distribution, projection, values, vectors, summaryRank, 0.0, rankTerms, containment);
        return new RankPosterior(rule, ranks, w, distribution, projection, values, vectors, summaryRank, posterior.Loss(summary), rankTerms, containment);
    }

    private static double[,] Leading(double[,] vectors, int k)
    {
        int p = vectors.GetLength(0);
        var basis = new double[p, k];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < k; j++) basis[i, j] = vectors[i, j];
        return basis;
    }
}
