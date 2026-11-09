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
        var v = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) v[i, j] = 0.5 * (cfg[i, j] + cfg[j, i]);
        EigenResult e = ScaledEigen.Of(v);
        int[] order = [.. Enumerable.Range(0, p).OrderByDescending(i => Math.Abs(e.Eigenvalues[i])).ThenBy(i => i)];
        var values = new double[p];
        var directions = new double[p, p];
        for (int j = 0; j < p; j++)
        {
            values[j] = e.Eigenvalues[order[j]];
            for (int i = 0; i < p; i++) directions[i, j] = e.Eigenvectors[order[j]][i];
        }
        return new CoActivePosterior(pairs, normalized, cf, cg, cfg, v, values, directions);
    }
}
