using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;

namespace Maths.Regression.Subspaces;

/// <summary>
/// One draw pair's matrices in one metric: C_f = E[∇f∇fᵀ], C_g = E[∇g∇gᵀ] and C_fg = E[∇f∇gᵀ], with the expected inner products
/// t = tr C and the concordance t_fg/√(t_f t_g) (Rumsey, Hardy, Ahrens and Vander Wiel, arXiv 2311.18146).
/// </summary>
public sealed record CoActivePair(double[,] Cf, double[,] Cg, double[,] Cfg)
{
    /// <summary>t_f = tr C_f.</summary>
    public double Tf => Trace(Cf);

    /// <summary>t_g = tr C_g.</summary>
    public double Tg => Trace(Cg);

    /// <summary>t_fg = tr C_fg = E[∇f·∇g].</summary>
    public double Tfg => Trace(Cfg);

    /// <summary>conc(f, g) = t_fg/√(t_f t_g) in [−1, 1], 1 when g = a + |b|f; NaN where either function is constant.</summary>
    public double Concordance => Tf > 0.0 && Tg > 0.0 ? Tfg / Math.Sqrt(Tf * Tg) : double.NaN;

    /// <summary>discord(f, g) = √((1 − conc)/2), a pseudo-metric.</summary>
    public double Discordance => Math.Sqrt(Math.Max(0.0, (1.0 - Concordance) / 2.0));

    internal static double Trace(double[,] m)
    {
        double t = 0.0;
        for (int i = 0; i < m.GetLength(0); i++) t += m[i, i];
        return t;
    }
}

/// <summary>
/// A posterior over co-active pairs (E8, slice 3), fit-agnostic: E[C_f], E[C_g] and E[C_fg] under the pairs' weights; the symmetrized
/// V = (E[C_fg] + E[C_fg]ᵀ)/2, whose eigenvectors are the co-active directions, each maximizing the expected product of the two
/// directional derivatives subject to orthogonality to the earlier ones. V may be indefinite, so its eigenpairs are ordered by |λ|
/// (2311.18146, §3.2). The contributions π_i = λ_i/√(E[t_f]E[t_g]) sum to the concordance of the means, E[t_fg]/√(E[t_f]E[t_g]); each
/// pair's own concordance is reported beside it, since the concordance of the means is not the mean of the concordances.
/// </summary>
public sealed class CoActivePosterior
{
    private CoActivePosterior(IReadOnlyList<CoActivePair> pairs, double[] weights, double[,] meanCf, double[,] meanCg, double[,] meanCfg,
                              double[,] v, double[] values, double[,] directions)
    {
        Pairs = pairs;
        Weights = weights;
        MeanCf = meanCf;
        MeanCg = meanCg;
        MeanCfg = meanCfg;
        V = v;
        Values = values;
        Directions = directions;
    }

    /// <summary>Each draw pair's matrices.</summary>
    public IReadOnlyList<CoActivePair> Pairs { get; }

    /// <summary>The pairs' weights, normalized.</summary>
    public IReadOnlyList<double> Weights { get; }

    /// <summary>E[C_f].</summary>
    public double[,] MeanCf { get; }

    /// <summary>E[C_g].</summary>
    public double[,] MeanCg { get; }

    /// <summary>E[C_fg].</summary>
    public double[,] MeanCfg { get; }

    /// <summary>V = (E[C_fg] + E[C_fg]ᵀ)/2.</summary>
    public double[,] V { get; }

    /// <summary>V's eigenvalues, ordered by magnitude, descending.</summary>
    public double[] Values { get; }

    /// <summary>The co-active directions, V's eigenvectors as columns in <see cref="Values"/>' order.</summary>
    public double[,] Directions { get; }

    /// <summary>E[t_f], E[t_g] and E[t_fg].</summary>
    public (double Tf, double Tg, double Tfg) MeanTraces =>
        (CoActivePair.Trace(MeanCf), CoActivePair.Trace(MeanCg), CoActivePair.Trace(MeanCfg));

    /// <summary>The concordance of the means, E[t_fg]/√(E[t_f]E[t_g]).</summary>
    public double Concordance
    {
        get
        {
            (double tf, double tg, double tfg) = MeanTraces;
            return tf > 0.0 && tg > 0.0 ? tfg / Math.Sqrt(tf * tg) : double.NaN;
        }
    }

    /// <summary>The contributions π_i = λ_i/√(E[t_f]E[t_g]), in <see cref="Values"/>' order; they sum to <see cref="Concordance"/>.</summary>
    public double[] Contributions
    {
        get
        {
            (double tf, double tg, _) = MeanTraces;
            double scale = Math.Sqrt(tf * tg);
            var pi = new double[Values.Length];
            for (int i = 0; i < pi.Length; i++) pi[i] = scale > 0.0 ? Values[i] / scale : double.NaN;
            return pi;
        }
    }

    /// <summary>Each pair's concordance, in the pairs' order.</summary>
    public double[] Concordances
    {
        get
        {
            var c = new double[Pairs.Count];
            for (int s = 0; s < c.Length; s++) c[s] = Pairs[s].Concordance;
            return c;
        }
    }

    /// <summary>
    /// The co-activity scores α_i(q) = Σ_{j≤q} λ_j w_ij² over the q leading co-active directions (2311.18146's eq. signed), or
    /// with |λ_j| when <paramref name="signed"/> is false (its unsigned score).
    /// </summary>
    public double[] CoActivityScores(int q, bool signed = true)
    {
        int p = Values.Length;
        if (q < 1 || q > p) throw new ArgumentOutOfRangeException(nameof(q), "q is from 1 to the inputs.");
        var scores = new double[p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < q; j++) scores[i] += (signed ? Values[j] : Math.Abs(Values[j])) * Directions[i, j] * Directions[i, j];
        return scores;
    }

    /// <summary>
    /// The co-active projection mean over the pairs' top-<paramref name="q"/> co-active subspaces (<see cref="CoActiveProjection"/>).
    /// </summary>
    public CoActiveProjection ProjectionMean(int q)
    {
        int p = Values.Length;
        if (q < 1 || q > p) throw new ArgumentOutOfRangeException(nameof(q), "q is from 1 to the inputs.");
        var bases = new double[Pairs.Count][,];
        var selected = new double[Pairs.Count, q];
        var projection = new double[p, p];
        for (int s = 0; s < Pairs.Count; s++)
        {
            (double[] values, double[,] directions) = ByMagnitude(Symmetrized(Pairs[s].Cfg));
            var basis = new double[p, q];
            for (int j = 0; j < q; j++)
            {
                selected[s, j] = values[j];
                for (int i = 0; i < p; i++) basis[i, j] = directions[i, j];
            }
            bases[s] = basis;
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                {
                    double pij = 0.0;
                    for (int r = 0; r < q; r++) pij += basis[i, r] * basis[j, r];
                    projection[i, j] += Weights[s] * pij;
                }
        }
        var eigenvalues = new double[p];
        var vectors = new double[p, p];
        ScaledEigen.Into(projection, eigenvalues, vectors);
        var mean = new double[p, q];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < q; j++) mean[i, j] = vectors[i, j];
        var angles = new double[Pairs.Count][];
        for (int s = 0; s < Pairs.Count; s++) angles[s] = SubspacePosterior.PrincipalAngles(bases[s], mean);
        // Along each mean direction u: E[(u·∇f)(u·∇g)]/√(E[(u·∇f)²]E[(u·∇g)²]), in [−1, 1] by Cauchy–Schwarz.
        var concordances = new double[q];
        for (int j = 0; j < q; j++)
        {
            double fg = Quadratic(V, vectors, j), f = Quadratic(MeanCf, vectors, j), g = Quadratic(MeanCg, vectors, j);
            concordances[j] = f > 0.0 && g > 0.0 ? fg / Math.Sqrt(f * g) : double.NaN;
        }
        return new CoActiveProjection(q, projection, eigenvalues, vectors, bases, selected, angles, concordances);
    }

    private static double[,] Symmetrized(double[,] m)
    {
        int p = m.GetLength(0);
        var v = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) v[i, j] = 0.5 * (m[i, j] + m[j, i]);
        return v;
    }

    // The symmetric eigensystem ordered by |λ|, descending, ties by index; eigenvectors as columns.
    private static (double[] Values, double[,] Directions) ByMagnitude(double[,] v)
    {
        int p = v.GetLength(0);
        var signed = new double[p];
        var columns = new double[p, p];
        ScaledEigen.Into(v, signed, columns);
        int[] order = [.. Enumerable.Range(0, p).OrderByDescending(i => Math.Abs(signed[i])).ThenBy(i => i)];
        var values = new double[p];
        var directions = new double[p, p];
        for (int j = 0; j < p; j++)
        {
            values[j] = signed[order[j]];
            for (int i = 0; i < p; i++) directions[i, j] = columns[i, order[j]];
        }
        return (values, directions);
    }

    // uᵀMu for column j of the vectors.
    private static double Quadratic(double[,] m, double[,] vectors, int j)
    {
        int p = m.GetLength(0);
        double s = 0.0;
        for (int a = 0; a < p; a++)
            for (int b = 0; b < p; b++) s += vectors[a, j] * m[a, b] * vectors[b, j];
        return s;
    }

    /// <summary>The posterior of <paramref name="pairs"/> under <paramref name="weights"/> (nonnegative, normalized here).</summary>
    public static CoActivePosterior Of(IReadOnlyList<CoActivePair> pairs, IReadOnlyList<double> weights)
    {
        ArgumentNullException.ThrowIfNull(pairs);
        ArgumentNullException.ThrowIfNull(weights);
        if (pairs.Count == 0) throw new ArgumentException("A posterior needs a pair.", nameof(pairs));
        if (weights.Count != pairs.Count) throw new ArgumentException("One weight per pair.", nameof(weights));
        int p = pairs[0].Cfg.GetLength(0);
        double total = 0.0;
        foreach (double w in weights)
        {
            if (!(w >= 0.0) || double.IsPositiveInfinity(w)) throw new ArgumentException("Weights are nonnegative and finite.", nameof(weights));
            total += w;
        }
        if (!(total > 0.0)) throw new ArgumentException("The weights need a positive sum.", nameof(weights));
        var normalized = new double[pairs.Count];
        double[,] cf = new double[p, p], cg = new double[p, p], cfg = new double[p, p];
        for (int s = 0; s < pairs.Count; s++)
        {
            CoActivePair pair = pairs[s];
            if (pair.Cf.GetLength(0) != p || pair.Cg.GetLength(0) != p || pair.Cfg.GetLength(0) != p || pair.Cfg.GetLength(1) != p)
                throw new ArgumentException("The pairs share their inputs.", nameof(pairs));
            double w = normalized[s] = weights[s] / total;
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                {
                    cf[i, j] += w * pair.Cf[i, j];
                    cg[i, j] += w * pair.Cg[i, j];
                    cfg[i, j] += w * pair.Cfg[i, j];
                }
        }
        double[,] v = Symmetrized(cfg);
        (double[] values, double[,] directions) = ByMagnitude(v);
        return new CoActivePosterior(pairs, normalized, cf, cg, cfg, v, values, directions);
    }
}

/// <summary>
/// The co-active projection mean (E8, slice 6, by slice 3's rule): decision 4's second summary carried over to pairs. Each pair's
/// top-q co-active subspace is spanned by its V_s = (C_fg + C_fgᵀ)/2's eigenvectors of the q largest |λ|, and the mean subspace by the
/// top-q eigenvectors of P̄ = Σ w_s W_s W_sᵀ, whose eigenvalues are mean squared cosines as <see cref="SubspacePosterior"/>'s are. A
/// subspace has no sign, so the signs are reported beside it: each pair's q selected eigenvalues, signed, and along each mean direction
/// u the directional concordance uᵀE[V]u/√(uᵀE[C_f]u·uᵀE[C_g]u), in [−1, 1], positive where f and g change together along u and
/// negative where they oppose.
/// </summary>
/// <param name="Dimension">q.</param>
/// <param name="Projection">P̄.</param>
/// <param name="Eigenvalues">P̄'s eigenvalues, descending.</param>
/// <param name="Vectors">P̄'s eigenvectors as columns, in order: the leading q span the mean co-active subspace.</param>
/// <param name="Bases">Each pair's top-q co-active subspace, p × q.</param>
/// <param name="Selected">Each pair's q selected eigenvalues of V_s, signed, by |λ| descending: pairs × q.</param>
/// <param name="AnglesToMean">Each pair's q principal angles to the mean subspace, ascending, in radians.</param>
/// <param name="Concordances">The directional concordance along each of the q leading mean directions; NaN where either function's
/// gradient has no energy along it.</param>
public sealed record CoActiveProjection(int Dimension, double[,] Projection, double[] Eigenvalues, double[,] Vectors, IReadOnlyList<double[,]> Bases,
                                        double[,] Selected, double[][] AnglesToMean, double[] Concordances);
