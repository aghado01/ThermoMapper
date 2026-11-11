using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// Several outputs' resolution on their data (E34): the outputs' names, means and centre levels (the means under level centring and the
/// flat intercept, 0 under the zero centre), S_Y, their mean squares and cross-products about the levels with scalar BASS's divisor n,
/// and Σ's prior, decision 2's IW(K − 1 + 2c, 2c·S_Y). The arrays are the resolution's own, not to be written.
/// </summary>
public sealed record BassResolvedOutputs(IReadOnlyList<string> Names, double[] Means, double[] Levels, double[,] MeanSquares, InverseWishart Prior,
                                        ResolvedOutputRates? Rates = null)
{
    /// <summary>K, the outputs.</summary>
    public int Count => Means.Length;
}

/// <summary>
/// A structure's least-squares fit of K outputs in the centred basis (E34): the K × K residual form R = Σᵢ rᵢrᵢᵀ, from the explicit
/// residuals or their bordered updates, and the fitted form F = B_cᵀG_c⁻¹B_c, the scalar ‖r‖² and Q_c carried to K outputs. The arrays
/// are not to be written: one structure's forms are shared by every state that holds it.
/// </summary>
public sealed record OutputForms(double[,] Residual, double[,] Fitted, OutputPencil? Prior = null)
{
    /// <summary>K, the outputs.</summary>
    public int Outputs => Residual.GetLength(0);
}

/// <summary>A design's K-output statistics as held, for a restart bit for bit: each term's cross products, the coefficients, the
/// residuals and the two forms.</summary>
internal sealed record DesignOutputsSnapshot(double[][] Cross, double[][] Coefficients, double[][] Residuals, double[,] Residual, double[,] Fitted);

/// <summary>
/// The K-output side of a design (E34): the outputs centred on their means, ỹ_k = y_k − ȳ_k (read under the reading's weights and a
/// random-effect block's whitening, E34's slice 4), each term's cross products z_jᵀỹ_k, the
/// least-squares coefficients β̂ (one row per term, one column per output), the explicit residuals and the K × K forms. The design's
/// structure, factor and admission are the scalar design's; this holds what the response adds, in place of its scalar statistics.
/// Not thread-safe: one per design.
/// </summary>
internal sealed class DesignOutputs
{
    private readonly double[][] _y;          // ỹ_k, one per output
    private readonly List<double[]> _cross = [];   // per term: z_jᵀỹ_k over the outputs
    private double[][] _coefficients = [];   // per term: β̂_jk over the outputs
    private double[][] _residuals;           // r_k, one per output
    private double[,] _residual, _fitted;    // R and F

    /// <summary>The K-output side on <paramref name="reads"/>, each output as the design reads it (W^½(y_k − ȳ_k) on the reading's values),
    /// which it takes over.</summary>
    internal DesignOutputs(double[][] reads)
    {
        int k = reads.Length;
        _y = reads;
        _residuals = new double[k][];
        for (int l = 0; l < k; l++) _residuals[l] = (double[])_y[l].Clone();
        _residual = Gram(_residuals);
        _fitted = new double[k, k];
    }

    /// <summary>K.</summary>
    public int Count => _y.Length;

    /// <summary>The centred outputs ỹ_k.</summary>
    public IReadOnlyList<double[]> Centred => _y;

    /// <summary>The current structure's forms.</summary>
    public OutputForms Forms => new((double[,])_residual.Clone(), (double[,])_fitted.Clone());

    /// <summary>β̂, one row per term in the factor's order, one column per output; copied.</summary>
    public double[,] Coefficients
    {
        get
        {
            int m = _coefficients.Length, k = Count;
            var beta = new double[m, k];
            for (int j = 0; j < m; j++)
                for (int l = 0; l < k; l++) beta[j, l] = _coefficients[j][l];
            return beta;
        }
    }

    /// <summary>The current residual form R.</summary>
    public double[,] Residual => _residual;

    /// <summary>Term <paramref name="index"/>'s coefficients over the outputs.</summary>
    public double[] CoefficientsOf(int index) => _coefficients[index];

    /// <summary>The current terms' cross products, with term <paramref name="skip"/> left out (−1 keeps all), and
    /// <paramref name="added"/> last when given: a candidate's B, one column per output.</summary>
    public double[][] CandidateCross(int skip, double[]? added)
    {
        int k = Count, m = _cross.Count - (skip >= 0 ? 1 : 0) + (added is null ? 0 : 1);
        var b = new double[k][];
        for (int l = 0; l < k; l++)
        {
            var column = new double[m];
            int at = 0;
            for (int j = 0; j < _cross.Count; j++)
                if (j != skip) column[at++] = _cross[j][l];
            if (added is not null) column[at] = added[l];
            b[l] = column;
        }
        return b;
    }

    /// <summary>zᵀv_k over the outputs, for <paramref name="v"/> the centred outputs or the residuals.</summary>
    public static double[] Dot(double[] column, IReadOnlyList<double[]> v)
    {
        var d = new double[v.Count];
        for (int l = 0; l < d.Length; l++) d[l] = PairwiseSum.Dot(column, v[l]);
        return d;
    }

    /// <summary>The residuals r_k.</summary>
    public IReadOnlyList<double[]> Residuals => _residuals;

    /// <summary>Takes the accepted move's cross products: term <paramref name="removed"/>'s for a death or a change, and
    /// <paramref name="added"/> last for a birth or a change.</summary>
    public void Accept(int removed, double[]? added)
    {
        if (removed >= 0) _cross.RemoveAt(removed);
        if (added is not null) _cross.Add(added);
    }

    /// <summary>Replaces every term's cross products with <paramref name="columns"/>' against the centred outputs.</summary>
    public void Load(IReadOnlyList<double[]> columns)
    {
        _cross.Clear();
        foreach (double[] z in columns) _cross.Add(Dot(z, _y));
    }

    /// <summary>
    /// Solves β̂_k = G⁻¹b_k with <paramref name="factor"/> for each output, recomputes each residual explicitly, r_k = ỹ_k − Z_cβ̂_k, so it
    /// never drifts, and the forms: R_kl = r_kᵀr_l and F_kl = b_kᵀβ̂_l, symmetrized, which at K = 1 is the scalar design's b·β̂.
    /// </summary>
    public void Refit(ISpdFactor factor, List<double[]> columns)
    {
        int k = Count, m = _cross.Count;
        double[][] b = CandidateCross(-1, null);
        var beta = new double[k][];
        for (int l = 0; l < k; l++) beta[l] = m == 0 ? [] : factor.Solve(b[l]);
        _coefficients = new double[m][];
        for (int j = 0; j < m; j++)
        {
            var row = new double[k];
            for (int l = 0; l < k; l++) row[l] = beta[l][j];
            _coefficients[j] = row;
        }
        for (int l = 0; l < k; l++)
        {
            var r = (double[])_y[l].Clone();
            if (m > 0) ColumnOps.SubtractScaledColumns(r, System.Runtime.InteropServices.CollectionsMarshal.AsSpan(columns)[..m], beta[l]);
            _residuals[l] = r;
        }
        _residual = Gram(_residuals);
        var fitted = new double[k, k];
        if (m > 0)
            for (int a = 0; a < k; a++)
                for (int d = 0; d <= a; d++)
                    fitted[a, d] = fitted[d, a] = a == d ? PairwiseSum.Dot(b[a], beta[a]) : 0.5 * (PairwiseSum.Dot(b[a], beta[d]) + PairwiseSum.Dot(b[d], beta[a]));
        _fitted = fitted;
    }

    /// <summary>The statistics as held.</summary>
    public DesignOutputsSnapshot Snapshot()
        => new([.. _cross.ConvertAll(c => (double[])c.Clone())], [.. Array.ConvertAll(_coefficients, c => (double[])c.Clone())],
               [.. Array.ConvertAll(_residuals, r => (double[])r.Clone())], (double[,])_residual.Clone(), (double[,])_fitted.Clone());

    /// <summary>Takes <paramref name="snapshot"/>'s statistics as they stand.</summary>
    public void Restore(DesignOutputsSnapshot snapshot)
    {
        if (snapshot.Residuals.Length != Count) throw new ArgumentException("The snapshot holds another count of outputs.", nameof(snapshot));
        _cross.Clear();
        foreach (double[] c in snapshot.Cross) _cross.Add((double[])c.Clone());
        _coefficients = Array.ConvertAll(snapshot.Coefficients, c => (double[])c.Clone());
        _residuals = Array.ConvertAll(snapshot.Residuals, r => (double[])r.Clone());
        _residual = (double[,])snapshot.Residual.Clone();
        _fitted = (double[,])snapshot.Fitted.Clone();
    }

    /// <summary>The Gram of <paramref name="vectors"/>, Σᵢ v_a,i v_d,i by pairwise sums.</summary>
    public static double[,] Gram(IReadOnlyList<double[]> vectors)
    {
        int k = vectors.Count;
        var g = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int d = 0; d <= a; d++) g[a, d] = g[d, a] = PairwiseSum.Dot(vectors[a], vectors[d]);
        return g;
    }
}
