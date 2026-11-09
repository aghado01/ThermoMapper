using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;

namespace Maths.Regression.Subspaces;

/// <summary>
/// The output side of an active-subspace readout (E8, slice 5; the Grassmann discussion's item 3), fit-agnostic: per draw
/// G_s = E[J W_x Jᵀ], the m × m matrix of the output patterns that move with the inputs, given as a factor F_s with G_s = F_sF_sᵀ
/// (m × c, the stacked metric-weighted Jacobians or their moment-route factor), so G is never formed. Each draw's spectrum and leading
/// patterns come from F_s's thin SVD; E[G] = Σ w_s G_s's from the thin SVD of the draws' compressed factors side by side,
/// [√w_s U_sΣ_s]. tr G_s is C_s's trace under the same metrics. A pattern's sign is set so that its largest entry is positive.
/// </summary>
public sealed class OutputPosterior
{
    private OutputPosterior(double[] weights, double[] traces, double[][] drawValues, double[][,] drawPatterns, double[] meanValues,
                            double[,] meanPatterns)
    {
        Weights = weights;
        Traces = traces;
        DrawEigenvalues = drawValues;
        DrawPatterns = drawPatterns;
        MeanEigenvalues = meanValues;
        MeanPatterns = meanPatterns;
    }

    /// <summary>The draws' weights, normalized.</summary>
    public IReadOnlyList<double> Weights { get; }

    /// <summary>Each draw's tr G_s.</summary>
    public IReadOnlyList<double> Traces { get; }

    /// <summary>Each draw's leading eigenvalues of G_s, descending, as many as patterns were asked (fewer if G_s's rank is lower).</summary>
    public IReadOnlyList<double[]> DrawEigenvalues { get; }

    /// <summary>Each draw's leading patterns, m × the eigenvalues' count, as columns.</summary>
    public IReadOnlyList<double[,]> DrawPatterns { get; }

    /// <summary>E[G]'s nonzero eigenvalues, descending.</summary>
    public double[] MeanEigenvalues { get; }

    /// <summary>E[G]'s leading patterns as columns, m × min(patterns asked, rank).</summary>
    public double[,] MeanPatterns { get; }

    /// <summary>tr E[G].</summary>
    public double MeanTrace => MeanEigenvalues.Sum();

    /// <summary>The posterior of the factors <paramref name="factors"/> (each m × c_s, G_s = F_sF_sᵀ) under <paramref name="weights"/>,
    /// with <paramref name="patterns"/> leading patterns per draw and for E[G].</summary>
    public static OutputPosterior Of(IReadOnlyList<double[,]> factors, IReadOnlyList<double> weights, int patterns)
    {
        ArgumentNullException.ThrowIfNull(factors);
        if (factors.Count == 0) throw new ArgumentException("A posterior needs a draw.", nameof(factors));
        int m = factors[0].GetLength(0);
        if (patterns < 1 || patterns > m) throw new ArgumentOutOfRangeException(nameof(patterns), $"Patterns run from 1 to the {m} output locations.");
        double[] w = Whitening.Normalized(weights, factors.Count, nameof(weights));
        var traces = new double[factors.Count];
        var values = new double[factors.Count][];
        var drawPatterns = new double[factors.Count][,];
        var compressed = new List<double[]>();   // columns of [√w_s U_sΣ_s]
        for (int s = 0; s < factors.Count; s++)
        {
            double[,] f = factors[s];
            if (f.GetLength(0) != m) throw new ArgumentException("The factors share their output locations.", nameof(factors));
            (double[] sigma, double[][] left) = Left(f);
            traces[s] = sigma.Sum(x => x * x);
            int k = Math.Min(patterns, sigma.Length);
            values[s] = [.. sigma.Take(k).Select(x => x * x)];
            drawPatterns[s] = Columns(left, k, m);
            for (int j = 0; j < sigma.Length; j++)
            {
                double scale = Math.Sqrt(w[s]) * sigma[j];
                compressed.Add([.. left[j].Select(x => x * scale)]);
            }
        }
        double[] meanValues;
        double[,] meanPatterns;
        if (compressed.Count == 0)
        {
            meanValues = [];
            meanPatterns = new double[m, 0];
        }
        else
        {
            var all = new double[m, compressed.Count];
            for (int j = 0; j < compressed.Count; j++)
                for (int l = 0; l < m; l++) all[l, j] = compressed[j][l];
            (double[] sigma, double[][] left) = Left(all);
            meanValues = [.. sigma.Select(x => x * x)];
            meanPatterns = Columns(left, Math.Min(patterns, sigma.Length), m);
        }
        return new OutputPosterior(w, traces, values, drawPatterns, meanValues, meanPatterns);
    }

    // The nonzero singular values of an m × c factor and their left vectors, signs set.
    internal static (double[] Sigma, double[][] Left) Left(double[,] f)
    {
        ThinSvdResult svd = ThinSvd.Decompose(f);
        int r = svd.NumericalRank;
        var sigma = new double[r];
        var left = new double[r][];
        for (int j = 0; j < r; j++)
        {
            sigma[j] = svd.SingularValues[j];
            double[] u = [.. svd.LeftVectors[j]];
            int big = 0;
            for (int l = 1; l < u.Length; l++)
                if (Math.Abs(u[l]) > Math.Abs(u[big])) big = l;
            if (u[big] < 0.0)
                for (int l = 0; l < u.Length; l++) u[l] = -u[l];
            left[j] = u;
        }
        return (sigma, left);
    }

    private static double[,] Columns(double[][] vectors, int k, int m)
    {
        var c = new double[m, k];
        for (int j = 0; j < k; j++)
            for (int l = 0; l < m; l++) c[l, j] = vectors[j][l];
        return c;
    }
}

/// <summary>
/// The matrix Poincaré check (E8, slice 5; decision 6 in matrix form), fit-agnostic. Under a product law, Efron–Stein and each input's
/// one-variable Poincaré inequality give Cov_ρ(f) ≼ E[J D Jᵀ], with D the diagonal of the inputs' Poincaré constants ((b − a)²/π² for a
/// uniform input on [a, b]). The ratio is the largest generalized eigenvalue of the covariance against the bound, at most 1 when the
/// inequality holds; the equality case is cos(πx) on [0, 1]. It is infinite where the covariance has mass outside the bound's range.
/// </summary>
public static class MatrixPoincare
{
    /// <summary>max_v vᵀAv / vᵀBv for symmetric positive semidefinite <paramref name="covariance"/> A and <paramref name="bound"/> B.</summary>
    public static double Ratio(double[,] covariance, double[,] bound)
    {
        ArgumentNullException.ThrowIfNull(covariance);
        ArgumentNullException.ThrowIfNull(bound);
        int m = bound.GetLength(0);
        if (bound.GetLength(1) != m || covariance.GetLength(0) != m || covariance.GetLength(1) != m)
            throw new ArgumentException("The covariance and the bound are m × m alike.", nameof(bound));
        double scale = 0.0;
        foreach (double x in bound) scale = Math.Max(scale, Math.Abs(x));
        foreach (double x in covariance) scale = Math.Max(scale, Math.Abs(x));
        if (!(scale > 0.0)) return 0.0;
        EigenResult e = ScaledEigen.Of(bound);
        double top = Math.Max(e.Eigenvalues.Max(), 0.0);
        var kept = new List<int>();
        for (int q = 0; q < m; q++)
            if (e.Eigenvalues[q] > 1e-12 * top) kept.Add(q);
        // The covariance outside the bound's range: any mass there makes the ratio infinite.
        double outside = 0.0;
        for (int q = 0; q < m; q++)
        {
            if (kept.Contains(q)) continue;
            outside = Math.Max(outside, Quadratic(covariance, e.Eigenvectors[q], e.Eigenvectors[q]));
        }
        if (outside > 1e-12 * scale) return double.PositiveInfinity;
        int r = kept.Count;
        if (r == 0) return 0.0;
        var reduced = new double[r, r];   // diag(β)^{-1/2} Vᵀ A V diag(β)^{-1/2} on the bound's range
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++)
                reduced[i, j] = Quadratic(covariance, e.Eigenvectors[kept[i]], e.Eigenvectors[kept[j]])
                                / Math.Sqrt(e.Eigenvalues[kept[i]] * e.Eigenvalues[kept[j]]);
        for (int i = 0; i < r; i++)
            for (int j = i + 1; j < r; j++) reduced[i, j] = reduced[j, i] = 0.5 * (reduced[i, j] + reduced[j, i]);
        return Math.Max(0.0, ScaledEigen.Of(reduced).Eigenvalues.Max());
    }

    private static double Quadratic(double[,] a, double[] u, double[] v)
    {
        int m = u.Length;
        double s = 0.0;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++) s += u[i] * a[i, j] * v[j];
        return s;
    }
}

/// <summary>
/// Basis adequacy (E8, slice 5): for a basis fit beside a direct-mode fit of the same cases, the share of the direct fit's G outside the
/// basis's span, tr((I − P_K)G)/tr G, and the containment defect of G's leading k-dimensional eigenspace in span(K),
/// k − ‖QᵀU_k‖²_F = Σ sin²θ_i, per draw and for E[G]. Both are 0 when the structure that moves with the inputs lies inside the span. A
/// basis fit's own G lies there by construction, which is why the direct fit is needed. It asks D07's question: did the truncation drop
/// structure that moves with the inputs?
/// </summary>
/// <param name="Shares">Each draw's share of tr G_s outside the span.</param>
/// <param name="Containment">Each draw's containment defect of its top-k patterns in the span.</param>
/// <param name="MeanShare">E[G]'s share outside the span.</param>
/// <param name="MeanContainment">E[G]'s top-k patterns' containment defect.</param>
/// <param name="Dimension">k, the patterns whose containment is read.</param>
/// <param name="SpanRank">The span's dimension.</param>
/// <param name="Outputs">The direct fit's output side the shares were read from.</param>
public sealed record BasisAdequacy(IReadOnlyList<double> Shares, IReadOnlyList<double> Containment, double MeanShare, double MeanContainment,
                                   int Dimension, int SpanRank, OutputPosterior Outputs)
{
    /// <summary>
    /// Adequacy of the factors <paramref name="factors"/> (G_s = F_sF_sᵀ, m × c_s) against the column span of <paramref name="span"/>
    /// (m × K, any basis of it), with <paramref name="dimension"/> leading patterns, under <paramref name="weights"/>.
    /// </summary>
    public static BasisAdequacy Of(IReadOnlyList<double[,]> factors, IReadOnlyList<double> weights, double[,] span, int dimension)
    {
        ArgumentNullException.ThrowIfNull(factors);
        ArgumentNullException.ThrowIfNull(span);
        OutputPosterior g = OutputPosterior.Of(factors, weights, dimension);
        int m = span.GetLength(0);
        (_, double[][] q) = OutputPosterior.Left(span);   // an orthonormal basis of the span
        var shares = new double[factors.Count];
        var containment = new double[factors.Count];
        double inside = 0.0;
        for (int s = 0; s < factors.Count; s++)
        {
            double[,] f = factors[s];
            if (f.GetLength(0) != m) throw new ArgumentException("The factors and the span share their output locations.", nameof(span));
            double total = 0.0, held = 0.0;
            for (int c = 0; c < f.GetLength(1); c++)
            {
                for (int l = 0; l < m; l++) total += f[l, c] * f[l, c];
                foreach (double[] qj in q)
                {
                    double d = 0.0;
                    for (int l = 0; l < m; l++) d += qj[l] * f[l, c];
                    held += d * d;
                }
            }
            shares[s] = total > 0.0 ? Math.Max(0.0, total - held) / total : double.NaN;
            inside += g.Weights[s] * held;
            containment[s] = Defect(g.DrawPatterns[s], q);
        }
        double meanTrace = g.Traces.Select((t, s) => t * g.Weights[s]).Sum();
        double meanShare = meanTrace > 0.0 ? Math.Max(0.0, meanTrace - inside) / meanTrace : double.NaN;
        return new BasisAdequacy(shares, containment, meanShare, Defect(g.MeanPatterns, q), dimension, q.Length, g);
    }

    // k − ‖QᵀU‖²_F for orthonormal patterns U (m × k) and an orthonormal basis Q of the span.
    private static double Defect(double[,] patterns, double[][] q)
    {
        int m = patterns.GetLength(0), k = patterns.GetLength(1);
        double held = 0.0;
        for (int j = 0; j < k; j++)
            foreach (double[] qi in q)
            {
                double d = 0.0;
                for (int l = 0; l < m; l++) d += qi[l] * patterns[l, j];
                held += d * d;
            }
        return Math.Max(0.0, k - held);
    }
}
