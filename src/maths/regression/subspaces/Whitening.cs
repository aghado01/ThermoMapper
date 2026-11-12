using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;

namespace Maths.Regression.Subspaces;

/// <summary>
/// Whitening for active subspaces under a dependent law (E8, slice 4; decision 2), fit-agnostic: a law's mean and covariance Σ from a
/// weighted joint sample, and the symmetric root L = Σ^{1/2}, which maps whitened coordinates to the inputs, x = μ + Lu. C in whitened
/// coordinates is LᵀCL, the generalized problem C v = λ Σ⁻¹ v: its eigenvalues are invariant under every invertible affine map of the
/// inputs, and its subspaces, read as functions of the inputs through L⁻ᵀ, map to the same functions. The symmetric root is the
/// whitening nearest the inputs' own coordinates (Kessy, Lewin and Strimmer, 2018, cited from memory); a Cholesky factor would depend on
/// the inputs' order.
/// </summary>
public static class Whitening
{
    /// <summary>
    /// The weighted mean and covariance of a sample read as its own law (divided by the weights' sum, not by n − 1): coordinate j of point
    /// r is <paramref name="columns"/>[j][r]; <paramref name="weights"/> are nonnegative and finite with a positive sum, normalized here,
    /// and uniform when null. Two passes, the covariance about the mean.
    /// </summary>
    public static (double[] Mean, double[,] Covariance) Moments(IReadOnlyList<double[]> columns, IReadOnlyList<double>? weights = null)
    {
        ArgumentNullException.ThrowIfNull(columns);
        int p = columns.Count;
        if (p == 0) throw new ArgumentException("A sample needs a coordinate.", nameof(columns));
        int n = columns[0]?.Length ?? throw new ArgumentException("Column 0 is null.", nameof(columns));
        if (n == 0) throw new ArgumentException("A sample needs a point.", nameof(columns));
        for (int j = 0; j < p; j++)
        {
            if (columns[j] is null || columns[j].Length != n) throw new ArgumentException($"Every column holds the {n} points.", nameof(columns));
            foreach (double x in columns[j])
                if (!double.IsFinite(x)) throw new ArgumentException("A sample's coordinates are finite.", nameof(columns));
        }
        double[] w = Normalized(weights, n, nameof(weights));
        var mean = new double[p];
        for (int j = 0; j < p; j++)
        {
            double s = 0.0;
            for (int r = 0; r < n; r++) s += w[r] * columns[j][r];
            mean[j] = s;
        }
        var covariance = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = i; j < p; j++)
            {
                double s = 0.0;
                for (int r = 0; r < n; r++) s += w[r] * (columns[i][r] - mean[i]) * (columns[j][r] - mean[j]);
                covariance[i, j] = covariance[j, i] = s;
            }
        return (mean, covariance);
    }

    /// <summary>
    /// Σ^{1/2}, symmetric, from Σ's eigendecomposition. Refused unless <paramref name="covariance"/> is square, finite, symmetric to
    /// rounding and positive definite, its smallest eigenvalue above 10⁻¹² of its largest: a degenerate law has no whitened coordinates.
    /// </summary>
    public static double[,] SymmetricRoot(double[,] covariance) => SymmetricRoots(covariance).Root;

    /// <summary>
    /// Σ^{1/2} and Σ^{-1/2}, both symmetric, from one eigendecomposition of Σ (spectral-solver brief, C5): the whitening map and the
    /// inverse a projection reads, with the eigenvectors shared rather than recomputed from the root. Refused as
    /// <see cref="SymmetricRoot"/> refuses.
    /// </summary>
    public static (double[,] Root, double[,] InverseRoot) SymmetricRoots(double[,] covariance)
    {
        ArgumentNullException.ThrowIfNull(covariance);
        int p = covariance.GetLength(0);
        if (p == 0 || covariance.GetLength(1) != p) throw new ArgumentException("A covariance is square.", nameof(covariance));
        double scale = 0.0;
        foreach (double x in covariance)
        {
            if (!double.IsFinite(x)) throw new ArgumentException("A covariance is finite.", nameof(covariance));
            scale = Math.Max(scale, Math.Abs(x));
        }
        var sym = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
            {
                if (Math.Abs(covariance[i, j] - covariance[j, i]) > 1e-12 * scale)
                    throw new ArgumentException("A covariance is symmetric.", nameof(covariance));
                sym[i, j] = 0.5 * (covariance[i, j] + covariance[j, i]);
            }
        EigenResult e = ScaledEigen.Of(sym);
        double top = e.Eigenvalues.Max(), bottom = e.Eigenvalues.Min();
        if (!(top > 0.0) || !(bottom > 1e-12 * top))
            throw new ArgumentException($"Whitening needs a positive definite covariance; its eigenvalues run from {bottom:G3} to {top:G3}.", nameof(covariance));
        return (MatrixFunctions.Symmetric(e, Math.Sqrt), MatrixFunctions.Symmetric(e, x => 1.0 / Math.Sqrt(x)));
    }

    // Weights normalized to sum to 1, uniform when null; refused unless nonnegative, finite, one per point and of positive sum.
    internal static double[] Normalized(IReadOnlyList<double>? weights, int n, string name)
    {
        if (weights is null) return [.. Enumerable.Repeat(1.0 / n, n)];
        if (weights.Count != n) throw new ArgumentException($"One weight per point: {weights.Count} for {n}.", name);
        double total = 0.0;
        foreach (double x in weights)
        {
            if (!(x >= 0.0) || double.IsPositiveInfinity(x)) throw new ArgumentException("Weights are nonnegative and finite.", name);
            total += x;
        }
        if (!(total > 0.0)) throw new ArgumentException("The weights need a positive sum.", name);
        return [.. weights.Select(x => x / total)];
    }
}
