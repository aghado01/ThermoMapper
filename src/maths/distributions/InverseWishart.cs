using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The inverse Wishart IW(ν, Ψ) on K × K positive definite matrices, in the conventional degrees of freedom: density
/// <c>|Ψ|^(ν/2) / (2^(νK/2) Γ_K(ν/2)) · |Σ|^(−(ν + K + 1)/2) · exp(−½ tr ΨΣ⁻¹)</c>, proper for ν &gt; K − 1, with mean Ψ/(ν − K − 1) for
/// ν &gt; K + 1. At K = 1 it is IG(ν/2, Ψ/2). Brown, Vannucci and Fearn's (1998) IW(δ; Q) is Dawid's shape, δ = ν − K + 1. The
/// conjugate draw of a multivariate normal's covariance: E34's Σ | Y, with A integrated, is IW(ν₀ + ν, Ψ₀ + R).
/// <para>A draw is Bartlett's, through Ψ's own factor: with Ψ = LLᵀ and A lower triangular, A_ii² ~ χ²(ν − i) and A_ij ~ N(0, 1)
/// below the diagonal, Σ⁻¹ = L⁻ᵀAAᵀL⁻¹ ~ W(ν, Ψ⁻¹), so Σ = (LA⁻ᵀ)(LA⁻ᵀ)ᵀ. Neither Ψ nor the Wishart is inverted, so a tiny χ²
/// pivot at small ν leaves the draw exact where W itself would be singular in doubles. At K = 1 it reads the inverse gamma's gamma
/// variate on the same stream, so the two draws agree to rounding.</para>
/// Immutable.
/// </summary>
public sealed class InverseWishart
{
    private readonly double[,] _scale;    // Ψ
    private readonly double[,] _root;     // L, Ψ = LLᵀ
    private readonly double _logNormalizer;   // (ν/2) ln|Ψ| − (νK/2) ln 2 − ln Γ_K(ν/2)

    private InverseWishart(double nu, double[,] scale, double[,] root, double logDetScale)
    {
        DegreesOfFreedom = nu;
        _scale = scale;
        _root = root;
        LogDeterminantScale = logDetScale;
        int k = scale.GetLength(0);
        _logNormalizer = 0.5 * nu * logDetScale - 0.5 * nu * k * SpecialFunctions.Ln2 - SpecialFunctions.LogMultivariateGamma(k, 0.5 * nu);
    }

    /// <summary>IW(<paramref name="degreesOfFreedom"/>, <paramref name="scale"/>): ν &gt; K − 1 and finite, Ψ symmetric positive definite
    /// (its lower triangle read, its symmetry checked); copied.</summary>
    public static InverseWishart Of(double degreesOfFreedom, double[,] scale)
    {
        ArgumentNullException.ThrowIfNull(scale);
        int k = scale.GetLength(0);
        if (k < 1 || scale.GetLength(1) != k) throw new ArgumentException("The scale is a nonempty square matrix.", nameof(scale));
        if (!(degreesOfFreedom > k - 1) || double.IsPositiveInfinity(degreesOfFreedom))
            throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom), $"IW(ν, Ψ) on {k} × {k} matrices is proper for finite ν > {k - 1}.");
        var copy = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j < k; j++)
            {
                double v = scale[i, j];
                if (!double.IsFinite(v)) throw new ArgumentException("The scale has a non-finite entry.", nameof(scale));
                if (v != scale[j, i]) throw new ArgumentException("The scale is not symmetric.", nameof(scale));
                copy[i, j] = v;
            }
        double[,] root = Cholesky(copy, out double logDet) ?? throw new ArgumentException("The scale is not positive definite.", nameof(scale));
        return new InverseWishart(degreesOfFreedom, copy, root, logDet);
    }

    /// <summary>K, the matrices' order.</summary>
    public int Dimension => _scale.GetLength(0);

    /// <summary>ν, the conventional degrees of freedom.</summary>
    public double DegreesOfFreedom { get; }

    /// <summary>Ψ, copied.</summary>
    public double[,] Scale => (double[,])_scale.Clone();

    /// <summary>ln|Ψ|.</summary>
    public double LogDeterminantScale { get; }

    /// <summary>Ψ/(ν − K − 1) for ν &gt; K + 1; null where the mean diverges.</summary>
    public double[,]? Mean
    {
        get
        {
            int k = Dimension;
            if (!(DegreesOfFreedom > k + 1)) return null;
            var mean = new double[k, k];
            double d = DegreesOfFreedom - k - 1;
            for (int i = 0; i < k; i++)
                for (int j = 0; j < k; j++) mean[i, j] = _scale[i, j] / d;
            return mean;
        }
    }

    /// <summary>ln p(<paramref name="sigma"/>); −∞ for a matrix that is not positive definite.</summary>
    public double LogPdf(double[,] sigma)
    {
        ArgumentNullException.ThrowIfNull(sigma);
        int k = Dimension;
        if (sigma.GetLength(0) != k || sigma.GetLength(1) != k) throw new ArgumentException($"Σ is {k} × {k}.", nameof(sigma));
        if (Cholesky(sigma, out double logDet) is not { } factor) return double.NegativeInfinity;
        return _logNormalizer - 0.5 * (DegreesOfFreedom + k + 1) * logDet - 0.5 * TraceSolve(factor, _scale);
    }

    /// <summary>One draw, Bartlett's through Ψ's factor; consumes, for i = 0 … K − 1, one gamma variate and then i normals.</summary>
    public double[,] Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        int k = Dimension;
        var a = new double[k, k];
        for (int i = 0; i < k; i++)
        {
            a[i, i] = Math.Sqrt(2.0 * Gamma.SampleStandard(rng, 0.5 * (DegreesOfFreedom - i)));
            for (int j = 0; j < i; j++) a[i, j] = Normal.SampleStandard(rng);
        }
        // T = A⁻¹, lower triangular; M = L·Tᵀ, so Σ = MMᵀ.
        var t = new double[k, k];
        for (int i = 0; i < k; i++)
        {
            t[i, i] = 1.0 / a[i, i];
            for (int j = 0; j < i; j++)
            {
                double sum = 0.0;
                for (int r = j; r < i; r++) sum -= a[i, r] * t[r, j];
                t[i, j] = sum / a[i, i];
            }
        }
        var m = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j < k; j++)
            {
                double sum = 0.0;
                for (int r = 0; r <= i; r++) sum += _root[i, r] * t[j, r];   // (L Tᵀ)_ij = Σ_r L_ir T_jr, L lower
                m[i, j] = sum;
            }
        var sigma = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = 0.0;
                for (int r = 0; r < k; r++) sum += m[i, r] * m[j, r];
                sigma[i, j] = sigma[j, i] = sum;
            }
        return sigma;
    }

    // The lower Cholesky factor and ln|A|, the lower triangle read; null when a pivot is not positive and finite.
    internal static double[,]? Cholesky(double[,] a, out double logDeterminant)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        logDeterminant = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = a[i, j];
                for (int r = 0; r < j; r++) sum -= l[i, r] * l[j, r];
                if (i != j)
                {
                    l[i, j] = sum / l[j, j];
                    continue;
                }
                if (!(sum > 0.0) || double.IsPositiveInfinity(sum)) return null;
                l[i, i] = Math.Sqrt(sum);
                logDeterminant += Math.Log(sum);
            }
        return l;
    }

    // tr(A⁻¹B) for A = LLᵀ: the diagonal of A⁻¹B, column by column.
    internal static double TraceSolve(double[,] factor, double[,] b)
    {
        int n = factor.GetLength(0);
        var column = new double[n];
        double trace = 0.0;
        for (int c = 0; c < n; c++)
        {
            for (int i = 0; i < n; i++)
            {
                double sum = b[i, c];
                for (int r = 0; r < i; r++) sum -= factor[i, r] * column[r];
                column[i] = sum / factor[i, i];
            }
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = column[i];
                for (int r = i + 1; r < n; r++) sum -= factor[r, i] * column[r];
                column[i] = sum / factor[i, i];
            }
            trace += column[c];
        }
        return trace;
    }
}
