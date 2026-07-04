#nullable enable
using System;
using System.Collections.Generic;

namespace TDA.Ph;

/// <summary>How the two orientations of an undirected edge's prediction residual are folded.</summary>
public enum ResidualSymmetry
{
    /// <summary>Admit as soon as <em>either</em> orientation's prediction is satisfied (the earliest support).</summary>
    Min,
    /// <summary>Require <em>both</em> orientations to agree (the later support).</summary>
    Max,
    /// <summary>Average the two orientations' deviations.</summary>
    Mean,
}

/// <summary>
/// P1a of the Conditioned-Persistence synthesis — the prior generalizes. Turns an observation
/// field <c>{t_i}</c> plus an edge prior <c>τ</c> into residual-weighted content edges for
/// <see cref="ConditionedFiltration"/>: an edge <c>(i,j)</c> predicting <c>t_j ≈ t_i + τ</c> carries
/// filtration value <c>r = |t_j − (t_i + τ)|</c> — the deviation of the observed gap from the
/// prediction. Sweeping the slack <c>δ</c> admits edges in residual order; that sweep <em>is</em> the
/// δ-filtration. The residual is only a weight, so everything downstream is the P0 engine
/// (union → Rips → involuted persistence → barcode).
/// <para><c>τ≡0</c> gives <c>r = |t_j − t_i|</c> — P0's raw-distance similarity, recovered exactly. A
/// nonzero prior shifts a <em>return's birth</em>: a predicted return is born early (unsurprising), an
/// unpredicted one late (a persistent surprise). Directedness (Li's native <c>i→j</c> residual), the
/// <c>Δ</c> reach axis, and the zigzag reader are P1b/P2 — this is the undirected, monotone-δ core.</para>
/// </summary>
public static class ResidualPrior
{
    /// <summary>
    /// Residual-weighted content edges from an observation field and an edge prior. Each prior entry
    /// <c>(i, j, tau)</c> is a candidate edge predicting <c>t_j ≈ t_i + tau</c>; its filtration value is
    /// the prediction residual, folded across both orientations — forward <c>|t_j − (t_i + tau)|</c> and
    /// reverse <c>|t_i − (t_j + tau)|</c> — by <paramref name="symmetry"/>. With <c>tau = 0</c> the
    /// residual is the raw gap <c>|t_j − t_i|</c> and the fold is inert, so P0's similarity is recovered.
    /// </summary>
    /// <param name="observations">The scalar field <c>{t_i}</c>, indexed by vertex.</param>
    /// <param name="prior">Candidate content edges with their predictions <c>(i, j, tau)</c>.</param>
    /// <param name="symmetry">How to fold the two orientations into one undirected residual.</param>
    public static IReadOnlyList<(int i, int j, double r)> ResidualEdges(
        double[] observations,
        IReadOnlyList<(int i, int j, double tau)> prior,
        ResidualSymmetry symmetry = ResidualSymmetry.Min)
    {
        ArgumentNullException.ThrowIfNull(observations);
        ArgumentNullException.ThrowIfNull(prior);

        var edges = new List<(int, int, double)>(prior.Count);
        foreach (var (i, j, tau) in prior)
        {
            if ((uint)i >= (uint)observations.Length || (uint)j >= (uint)observations.Length)
                throw new ArgumentOutOfRangeException(nameof(prior),
                    $"Edge ({i},{j}) out of range for {observations.Length} observations.");

            double forward = Math.Abs(observations[j] - (observations[i] + tau));
            double reverse = Math.Abs(observations[i] - (observations[j] + tau));
            double r = symmetry switch
            {
                ResidualSymmetry.Min  => Math.Min(forward, reverse),
                ResidualSymmetry.Max  => Math.Max(forward, reverse),
                ResidualSymmetry.Mean => 0.5 * (forward + reverse),
                _ => throw new ArgumentOutOfRangeException(nameof(symmetry)),
            };
            edges.Add((i, j, r));
        }

        return edges;
    }
}
