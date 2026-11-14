using System;
using System.Collections.Generic;

namespace Maths.LinAlg;

/// <summary>
/// Functions of a symmetric matrix read from its eigensystem A = VΛVᵀ (spectral-solver brief, B6): f(A) = V·f(Λ)·Vᵀ, and the factor
/// forms V·f(Λ) over chosen eigenpairs, as columns or as rows, so a root F with FFᵀ = A or a row root R with RᵀR = A comes from one
/// decomposition. Which eigenpairs count, and what is refused, are the caller's rules: the helpers take the selection, not a cutoff.
/// The exponential of a skew matrix and the logarithm of a rotation are read the same way, from symmetric matrices that commute
/// with them.
/// </summary>
public static class MatrixFunctions
{
    /// <summary>
    /// f(A) = Σ_q f(λ_q) v_q v_qᵀ, symmetrized to (B + Bᵀ)/2, from <paramref name="eigen"/> (eigenvectors as rows of
    /// <see cref="EigenResult.Eigenvectors"/>). Terms with f(λ_q) = 0 and an eigenvector's zero components are skipped: they would add a
    /// signed zero to a sum that starts at +0, which leaves it unchanged, so the result is the same bits and a diagonal A costs O(n²).
    /// </summary>
    public static double[,] Symmetric(EigenResult eigen, Func<double, double> f)
    {
        ArgumentNullException.ThrowIfNull(eigen);
        ArgumentNullException.ThrowIfNull(f);
        int n = eigen.Eigenvalues.Length;
        var result = new double[n, n];
        var nonzero = new int[n];
        for (int q = 0; q < n; q++)
        {
            double weight = f(eigen.Eigenvalues[q]);
            if (weight == 0.0) continue;
            double[] v = eigen.Eigenvectors[q];
            int count = 0;
            for (int i = 0; i < n; i++)
                if (v[i] != 0.0) nonzero[count++] = i;
            for (int a = 0; a < count; a++)
            {
                int i = nonzero[a];
                for (int b = 0; b < count; b++)
                {
                    int j = nonzero[b];
                    result[i, j] += weight * v[i] * v[j];
                }
            }
        }
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) result[i, j] = result[j, i] = 0.5 * (result[i, j] + result[j, i]);
        return result;
    }

    /// <summary>The factor with columns f(λ_q)·v_q for the eigenpairs <paramref name="kept"/> selects, in its order: n × kept.Count,
    /// so with f = √ and the positive eigenpairs F·Fᵀ = A.</summary>
    public static double[,] FactorColumns(EigenResult eigen, IReadOnlyList<int> kept, Func<double, double> f)
    {
        ArgumentNullException.ThrowIfNull(eigen);
        ArgumentNullException.ThrowIfNull(kept);
        ArgumentNullException.ThrowIfNull(f);
        int n = eigen.Eigenvalues.Length;
        var factor = new double[n, kept.Count];
        for (int c = 0; c < kept.Count; c++)
        {
            double weight = f(eigen.Eigenvalues[kept[c]]);
            double[] v = eigen.Eigenvectors[kept[c]];
            for (int a = 0; a < n; a++) factor[a, c] = weight * v[a];
        }
        return factor;
    }

    /// <summary>The factor with rows f(λ_q)·v_qᵀ for the eigenpairs <paramref name="kept"/> selects, in its order: kept.Count × n, so
    /// with f = √ and the positive eigenpairs Rᵀ·R = A.</summary>
    public static double[,] FactorRows(EigenResult eigen, IReadOnlyList<int> kept, Func<double, double> f)
    {
        ArgumentNullException.ThrowIfNull(eigen);
        ArgumentNullException.ThrowIfNull(kept);
        ArgumentNullException.ThrowIfNull(f);
        int n = eigen.Eigenvalues.Length;
        var factor = new double[kept.Count, n];
        for (int r = 0; r < kept.Count; r++)
        {
            double weight = f(eigen.Eigenvalues[kept[r]]);
            double[] v = eigen.Eigenvectors[kept[r]];
            for (int a = 0; a < n; a++) factor[r, a] = weight * v[a];
        }
        return factor;
    }

    /// <summary>
    /// exp(Ω) of a real skew-symmetric Ω, read from the eigensystem of the symmetric X = ΩᵀΩ = −Ω². The even and odd terms of the
    /// series give exp(Ω) = cos(√X) + (sin √X / √X)·Ω, both functions of X, which commutes with Ω. The result is orthogonal with
    /// determinant 1, to rounding. Ω is read through its skew part (Ω − Ωᵀ)/2.
    /// </summary>
    public static double[,] SkewExponential(double[,] skew)
    {
        double[,] omega = SkewPart(skew, nameof(skew));
        int m = omega.GetLength(0);
        var x = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = i; j < m; j++)
            {
                double s = 0.0;
                for (int k = 0; k < m; k++) s += omega[k, i] * omega[k, j];
                x[i, j] = x[j, i] = s;
            }
        EigenResult eigen = Eigen.DecomposeSymmetric(x, tol: EigenTolerance);
        double[,] even = Symmetric(eigen, lambda => Math.Cos(Math.Sqrt(Math.Max(lambda, 0.0))));
        double[,] odd = Symmetric(eigen, SinOverRoot);
        var result = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
            {
                double s = even[i, j];
                for (int k = 0; k < m; k++) s += odd[i, k] * omega[k, j];
                result[i, j] = s;
            }
        return result;
    }

    /// <summary>
    /// The principal logarithm of a real orthogonal V with no eigenvalue −1, which also makes det V = 1: the skew L with exp(L) = V
    /// whose eigenvalues ±iθ have θ in [0, π). V is normal, so S = (V + Vᵀ)/2 and K = (V − Vᵀ)/2 commute. On each invariant plane S is
    /// cos θ and K is sin θ times a unit rotation, so L = G·K with G = Σ (θ_q/sin θ_q)·v_q v_qᵀ over S's eigenvectors v_q. Each angle
    /// is θ_q = atan2(‖K v_q‖, λ_q), the sine read from K and the cosine from S's eigenvalue, never an arccos: a small rotation keeps
    /// its digits through K, and one near a half-turn through its sine, which an arccos of a cosine near −1 would amplify by 1/sin θ.
    /// Throws where an angle lies within <paramref name="cutLocus"/> radians of π, where the logarithm stops being unique;
    /// <see cref="TryRotationLogarithm"/> reports it.
    /// </summary>
    public static double[,] RotationLogarithm(double[,] rotation, double cutLocus = 1e-8)
    {
        if (!TryRotationLogarithm(rotation, out double[,] log, cutLocus))
            throw new ArgumentException(
                "The matrix has an eigenvalue at −1 (a half-turn, or determinant −1), so it has no unique real logarithm.",
                nameof(rotation));
        return log;
    }

    /// <summary><see cref="RotationLogarithm"/>, returning false instead of throwing at an angle within
    /// <paramref name="cutLocus"/> radians of π.</summary>
    public static bool TryRotationLogarithm(double[,] rotation, out double[,] log, double cutLocus = 1e-8)
    {
        ArgumentNullException.ThrowIfNull(rotation);
        if (!(cutLocus >= 0.0 && cutLocus < Math.PI))
            throw new ArgumentOutOfRangeException(nameof(cutLocus), "The cut-locus margin must be in [0, π).");
        int m = rotation.GetLength(0);
        if (rotation.GetLength(1) != m) throw new ArgumentException("The matrix must be square.", nameof(rotation));
        var sym = new double[m, m];
        var skew = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
            {
                sym[i, j] = 0.5 * (rotation[i, j] + rotation[j, i]);
                skew[i, j] = 0.5 * (rotation[i, j] - rotation[j, i]);
            }
        EigenResult eigen = Eigen.DecomposeSymmetric(sym, tol: EigenTolerance);
        var g = new double[m, m];
        for (int q = 0; q < m; q++)
        {
            double[] v = eigen.Eigenvectors[q];
            double sineSquared = 0.0;
            for (int i = 0; i < m; i++)
            {
                double s = 0.0;
                for (int j = 0; j < m; j++) s += skew[i, j] * v[j];
                sineSquared += s * s;
            }
            double sine = Math.Sqrt(sineSquared), angle = Math.Atan2(sine, eigen.Eigenvalues[q]);
            if (angle >= Math.PI - cutLocus)
            {
                log = new double[0, 0];
                return false;
            }
            double weight = AngleOverSine(angle, sine);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++) g[i, j] += weight * v[i] * v[j];
        }
        log = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
            {
                double s = 0.0;
                for (int k = 0; k < m; k++) s += g[i, k] * skew[k, j];
                log[i, j] = s;
            }
        for (int i = 0; i < m; i++)
        {
            log[i, i] = 0.0;
            for (int j = i + 1; j < m; j++)
            {
                double a = 0.5 * (log[i, j] - log[j, i]);
                log[i, j] = a;
                log[j, i] = -a;
            }
        }
        return true;
    }

    // Element-wise relative Jacobi tolerance for the skew and rotation functions, about five ulps: a pivot left below
    // 1e-15·√|a_pp a_qq| mixes eigenvectors by about 1e-15 over their gap. The logarithm near a half-turn amplifies that mixing by
    // about 1/(π − θ), so 1e-14 left errors of 3e-11 at π − θ = 1e-3, and this tolerance 4e-12.
    private const double EigenTolerance = 1e-15;

    private static double[,] SkewPart(double[,] a, string name)
    {
        ArgumentNullException.ThrowIfNull(a, name);
        int m = a.GetLength(0);
        if (a.GetLength(1) != m) throw new ArgumentException("The matrix must be square.", name);
        var skew = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = i + 1; j < m; j++)
            {
                double s = 0.5 * (a[i, j] - a[j, i]);
                skew[i, j] = s;
                skew[j, i] = -s;
            }
        return skew;
    }

    // sin(√λ)/√λ, with its series where √λ is small.
    private static double SinOverRoot(double lambda)
    {
        if (lambda <= 0.0) return 1.0;
        double t = Math.Sqrt(lambda);
        return t < 1e-4 ? 1.0 - lambda / 6.0 + lambda * lambda / 120.0 : Math.Sin(t) / t;
    }

    // θ/sin θ from an angle in [0, π) and its sine, with the series 1 + θ²/6 + 7θ⁴/360 where θ is small (and the 0/0 at θ = 0).
    private static double AngleOverSine(double angle, double sine)
    {
        if (angle < 1e-4) return 1.0 + angle * angle / 6.0 + 7.0 * angle * angle * angle * angle / 360.0;
        return angle / sine;
    }
}
