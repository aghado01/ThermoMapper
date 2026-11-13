using System;
using System.Collections.Generic;
using Maths.LinAlg;

namespace Maths.Regression.Subspaces;

/// <summary>One draw's C with its eigenvalues, descending, and its top-k eigenspace as a p × k orthonormal basis.</summary>
public sealed record DrawSubspace(double[,] C, double[] Eigenvalues, double[,] Basis)
{
    private readonly double[,]? _vectors;

    /// <summary>All the eigenvectors the draw's subspace was read from, as columns in the eigenvalues' order (p × p from C, or as many
    /// as the factor's rank); the basis when none were kept. A per-draw rank (slice 5) reads its leading k_s.</summary>
    public double[,] Vectors { get => _vectors ?? Basis; init => _vectors = value; }

    /// <summary>The draw's subspace from C itself, by the symmetric eigendecomposition.</summary>
    public static DrawSubspace FromC(double[,] c, int dimension)
    {
        ArgumentNullException.ThrowIfNull(c);
        int p = Check(c, dimension);
        var values = new double[p];
        var vectors = new double[p, p];
        DenseEigen.DecomposeSymmetricInto(c, values, vectors, tol: 1e-14, stoppingRule: JacobiStoppingRule.ElementWiseRelative);
        var basis = new double[p, dimension];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < dimension; j++) basis[i, j] = vectors[i, j];
        return new DrawSubspace(c, values, basis) { Vectors = vectors };
    }

    /// <summary>
    /// The draw's subspace from a factor F with C = FᵀF, N × p (a sampled route's square-root-weighted gradient matrix): the right
    /// singular vectors of F's thin SVD, whose squared singular values are C's eigenvalues, so the gaps are not squared as forming C
    /// squares them (decision 5). <paramref name="c"/> is C as the caller accumulated it, kept as the draw's C rather than formed again
    /// from F (spectral-solver brief, C4); with fewer factor rows than <paramref name="dimension"/>, the subspace is C's own eigenspace.
    /// </summary>
    public static DrawSubspace FromFactor(double[,] factor, double[,] c, int dimension)
    {
        ArgumentNullException.ThrowIfNull(factor);
        ArgumentNullException.ThrowIfNull(c);
        int p = factor.GetLength(1);
        if (dimension < 1 || dimension > p) throw new ArgumentOutOfRangeException(nameof(dimension), "The dimension is from 1 to the inputs.");
        if (c.GetLength(0) != p || c.GetLength(1) != p) throw new ArgumentException("C is p × p, one row and column per factor column.", nameof(c));
        ThinSvdResult svd = ThinSvd.Decompose(factor);
        var eigenvalues = new double[p];
        for (int i = 0; i < svd.SingularValues.Length; i++) eigenvalues[i] = svd.SingularValues[i] * svd.SingularValues[i];
        if (svd.RightVectors.Length < dimension) return FromC(c, dimension);   // fewer points than the dimension: C's own eigenspace
        var basis = new double[p, dimension];
        var vectors = new double[p, svd.RightVectors.Length];
        for (int j = 0; j < svd.RightVectors.Length; j++)
            for (int i = 0; i < p; i++)
            {
                vectors[i, j] = svd.RightVectors[j][i];
                if (j < dimension) basis[i, j] = vectors[i, j];
            }
        return new DrawSubspace(c, eigenvalues, basis) { Vectors = vectors };
    }

    /// <summary>The gap λ_k − λ_{k+1} at the basis's dimension (λ_{p+1} = 0).</summary>
    public double Gap => Eigenvalues[Basis.GetLength(1) - 1] - (Basis.GetLength(1) < Eigenvalues.Length ? Eigenvalues[Basis.GetLength(1)] : 0.0);

    private static int Check(double[,] c, int dimension)
    {
        int p = c.GetLength(0);
        if (c.GetLength(1) != p) throw new ArgumentException("C is square.", nameof(c));
        if (dimension < 1 || dimension > p) throw new ArgumentOutOfRangeException(nameof(dimension), "The dimension is from 1 to the inputs.");
        return p;
    }
}

/// <summary>
/// A posterior over active subspaces from per-draw C (E8, decision 4; G3, SG-5), fit-agnostic. Two summaries, reported apart:
/// <list type="bullet">
/// <item><b>E[C]</b> = Σ w_s C_s and its eigensystem: the directions that carry the average squared gradient, a draw weighted by its
/// gradient's size.</item>
/// <item><b>The projection mean:</b> the top-k eigenvectors of P̄ = Σ w_s W_s W_sᵀ, the Bayes action under chordal loss, each draw
/// counted once. P̄'s eigenvalues are mean squared cosines, Σ w_s‖W_sᵀu_j‖², summing to k: 0.5 arises both when every draw sits at 45°
/// to u_j and when half the draws contain it, so each draw's principal angles to the mean are reported beside them. Its eigenvectors
/// in order give the mean subspace for every k (the flag mean, cited from memory).</item>
/// </list>
/// The angles between E[C]'s top-k eigenspace and the mean subspace say how far apart the two summaries sit; they differ when the
/// draws disagree (SG-5). Angles come from LinAlg's stable routine (slice 5), accurate to about √(p·k)·ε rad; slice 2 took the arccos
/// of the cosines, which resolves about 1.5·10⁻⁸ rad (decision 5).
/// Outlying draws are reported through their angles, not down-weighted.
/// </summary>
public sealed class SubspacePosterior
{
    private SubspacePosterior(IReadOnlyList<DrawSubspace> draws, double[] weights, int dimension, double[,] meanC, double[] meanCValues,
                              double[,] meanCVectors, double[,] projection, double[] projectionValues, double[,] projectionVectors,
                              double[][] anglesToMean, double[] meanCToMean)
    {
        Draws = draws;
        Weights = weights;
        Dimension = dimension;
        MeanC = meanC;
        MeanCEigenvalues = meanCValues;
        MeanCVectors = meanCVectors;
        Projection = projection;
        ProjectionEigenvalues = projectionValues;
        ProjectionVectors = projectionVectors;
        AnglesToMean = anglesToMean;
        MeanCToMean = meanCToMean;
    }

    /// <summary>Each draw's C, eigenvalues and top-k basis, in the metric the caller read them in.</summary>
    public IReadOnlyList<DrawSubspace> Draws { get; }

    /// <summary>The draws' weights, normalized.</summary>
    public IReadOnlyList<double> Weights { get; }

    /// <summary>k.</summary>
    public int Dimension { get; }

    /// <summary>E[C] = Σ w_s C_s.</summary>
    public double[,] MeanC { get; }

    /// <summary>E[C]'s eigenvalues, descending.</summary>
    public double[] MeanCEigenvalues { get; }

    /// <summary>E[C]'s eigenvectors as columns, in the eigenvalues' order.</summary>
    public double[,] MeanCVectors { get; }

    /// <summary>P̄ = Σ w_s W_s W_sᵀ.</summary>
    public double[,] Projection { get; }

    /// <summary>P̄'s eigenvalues, descending: mean squared cosines in [0, 1], summing to k.</summary>
    public double[] ProjectionEigenvalues { get; }

    /// <summary>P̄'s eigenvectors as columns, in order: the leading k span the mean subspace.</summary>
    public double[,] ProjectionVectors { get; }

    /// <summary>Each draw's k principal angles to the mean subspace, ascending, in radians.</summary>
    public double[][] AnglesToMean { get; }

    /// <summary>The k principal angles between E[C]'s top-k eigenspace and the mean subspace, ascending, in radians.</summary>
    public double[] MeanCToMean { get; }

    /// <summary>Each draw's gap λ_k − λ_{k+1}.</summary>
    public double[] Gaps
    {
        get
        {
            var gaps = new double[Draws.Count];
            for (int s = 0; s < gaps.Length; s++) gaps[s] = Draws[s].Gap;
            return gaps;
        }
    }

    /// <summary>The posterior of per-draw <paramref name="c"/>, each draw's subspace from C's eigendecomposition.</summary>
    public static SubspacePosterior Of(IReadOnlyList<double[,]> c, IReadOnlyList<double> weights, int dimension)
    {
        ArgumentNullException.ThrowIfNull(c);
        var draws = new DrawSubspace[c.Count];
        for (int s = 0; s < draws.Length; s++) draws[s] = DrawSubspace.FromC(c[s], dimension);
        return Of(draws, weights);
    }

    /// <summary>The posterior of <paramref name="draws"/>, each with its subspace already read, under <paramref name="weights"/>
    /// (nonnegative, normalized here).</summary>
    public static SubspacePosterior Of(IReadOnlyList<DrawSubspace> draws, IReadOnlyList<double> weights)
    {
        ArgumentNullException.ThrowIfNull(draws);
        ArgumentNullException.ThrowIfNull(weights);
        if (draws.Count == 0) throw new ArgumentException("A posterior needs a draw.", nameof(draws));
        if (weights.Count != draws.Count) throw new ArgumentException("One weight per draw.", nameof(weights));
        int p = draws[0].C.GetLength(0), k = draws[0].Basis.GetLength(1);
        double total = 0.0;
        foreach (double w in weights)
        {
            if (!(w >= 0.0) || double.IsPositiveInfinity(w)) throw new ArgumentException("Weights are nonnegative and finite.", nameof(weights));
            total += w;
        }
        if (!(total > 0.0)) throw new ArgumentException("The weights need a positive sum.", nameof(weights));
        var normalized = new double[draws.Count];
        var meanC = new double[p, p];
        var projection = new double[p, p];
        for (int s = 0; s < draws.Count; s++)
        {
            DrawSubspace d = draws[s];
            if (d.C.GetLength(0) != p || d.Basis.GetLength(0) != p || d.Basis.GetLength(1) != k)
                throw new ArgumentException("The draws share their inputs and dimension.", nameof(draws));
            double w = normalized[s] = weights[s] / total;
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                {
                    meanC[i, j] += w * d.C[i, j];
                    double pij = 0.0;
                    for (int r = 0; r < k; r++) pij += d.Basis[i, r] * d.Basis[j, r];
                    projection[i, j] += w * pij;
                }
        }
        (double[] meanCValues, double[,] meanCVectors) = Eigen(meanC);
        (double[] projectionValues, double[,] projectionVectors) = Eigen(projection);
        double[,] mean = Leading(projectionVectors, k), meanCSubspace = Leading(meanCVectors, k);
        var angles = new double[draws.Count][];
        for (int s = 0; s < draws.Count; s++) angles[s] = PrincipalAngles(draws[s].Basis, mean);
        return new SubspacePosterior(draws, normalized, k, meanC, meanCValues, meanCVectors, projection, projectionValues, projectionVectors,
                                     angles, PrincipalAngles(meanCSubspace, mean));
    }

    /// <summary>The principal angles between the spans of the orthonormal bases <paramref name="a"/> (p × k_a) and <paramref name="b"/>
    /// (p × k_b), min(k_a, k_b) of them, ascending, none when either is empty: LinAlg's stable routine
    /// (<see cref="global::Maths.LinAlg.PrincipalAngles.Between(double[,], double[,])"/>), which pairs the cosines of aᵀb with the sines
    /// of the residual off the other span through atan2, so small angles keep their accuracy where an arccos of the cosines alone
    /// resolves only about 2·10⁻⁸ rad (E8, slice 5, the Grassmann discussion's item 2).</summary>
    public static double[] PrincipalAngles(double[,] a, double[,] b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        if (b.GetLength(0) != a.GetLength(0)) throw new ArgumentException("The bases have one row per input.", nameof(b));
        if (a.GetLength(1) == 0 || b.GetLength(1) == 0) return [];
        return global::Maths.LinAlg.PrincipalAngles.Between(a, b);
    }

    // The symmetric eigendecomposition, eigenvalues descending, eigenvectors as columns.
    private static (double[] Values, double[,] Vectors) Eigen(double[,] m)
    {
        int p = m.GetLength(0);
        var values = new double[p];
        var vectors = new double[p, p];
        DenseEigen.DecomposeSymmetricInto(m, values, vectors, tol: 1e-14, stoppingRule: JacobiStoppingRule.ElementWiseRelative);
        return (values, vectors);
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

