using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using static Maths.LinAlg.Tests.FactorReferenceTests;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The strict dense reference factor: no pivot floor, no inverse. Solves and log-determinants are checked against
/// an independent LU, the innovation map against the exact inverse, strictness against matrices whose failing pivot
/// is exact, and the relative pivot against its two meanings (the squared sine of a column's angle to the earlier
/// columns, and a lower bound on the condition number of the equilibrated matrix).
/// </summary>
public sealed class DenseCholeskyTests
{
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(7)]
    [InlineData(40)]
    public void Factor_ReconstructsTheMatrix_AndSolvesAgainstLu(int n)
    {
        var rng = new Xoshiro256PlusPlus(seed: 500 + n);
        double[,] a = DenseSpd(n, rng);
        double[] b = Vector(n, rng);

        var f = new DenseCholesky(n);
        f.Decompose(a);
        Assert.True(f.Status.PositiveDefinite);
        Assert.Equal(0, f.Status.FlooredPivots);
        Assert.Equal(-1, f.FailedAt);

        double scale = 0.0;
        for (int i = 0; i < n; i++) scale = Math.Max(scale, Math.Abs(a[i, i]));
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = 0.0;
                for (int k = 0; k <= j; k++) s += f.Factor(i, k) * f.Factor(j, k);
                Assert.True(Math.Abs(s - a[i, j]) <= 4.0 * n * 2.2e-16 * scale, $"(LLᵀ)[{i},{j}] = {s:R} vs {a[i, j]:R}");
            }
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) Assert.Equal(0.0, f.Factor(i, j));

        (double[] x, double logAbsDet) = LuSolve(a, b);
        AssertClose(x, f.Solve(b), 1e-11);
        Assert.True(Math.Abs(f.LogDet - logAbsDet) <= 1e-11 * Math.Max(1.0, Math.Abs(logAbsDet)), $"ln|A| {f.LogDet:R} vs {logAbsDet:R}");

        double bx = 0.0;
        for (int i = 0; i < n; i++) bx += b[i] * x[i];
        Assert.True(Math.Abs(f.InverseQuadraticForm(b) - bx) <= 1e-11 * Math.Abs(bx), $"bᵀA⁻¹b {f.InverseQuadraticForm(b):R} vs {bx:R}");
    }

    [Fact]
    public void Innovation_HasExactlyTheInverseCovariance()
    {
        // v_k = L⁻ᵀ e_k, so Σ_k v_k v_kᵀ = L⁻ᵀ L⁻¹ = A⁻¹ and A·Σ v_k v_kᵀ = I.
        const int n = 9;
        double[,] a = DenseSpd(n, new Xoshiro256PlusPlus(seed: 61));
        var f = new DenseCholesky(n);
        f.Decompose(a);

        var cov = new double[n, n];
        var e = new double[n];
        var v = new double[n];
        for (int k = 0; k < n; k++)
        {
            Array.Clear(e);
            e[k] = 1.0;
            f.SampleInnovation(e, v);
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cov[i, j] += v[i] * v[j];
        }
        AssertIdentity(Multiply(a, cov), 1e-12);
    }

    [Fact]
    public void Factor_IsTheFlooredFactor_WhereNoPivotIsFloored()
    {
        // Same Cholesky–Banachiewicz order as CholeskyDecomposition, so with no pivot near the floor they agree exactly.
        const int n = 12;
        double[,] a = DenseSpd(n, new Xoshiro256PlusPlus(seed: 62));
        var strict = new DenseCholesky(n);
        strict.Decompose(a);
        var floored = new CholeskyDecomposition(n);
        floored.Decompose(a);
        var l = new double[n * n];
        floored.WriteLTo(l);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Assert.Equal(l[i * n + j], strict.Factor(i, j));
        Assert.Equal(floored.LogDet, strict.LogDet);
        Assert.Equal(floored.Status, strict.Status);
    }

    [Theory]
    [InlineData(1.0, 1.0, 1.0, 1, 0.0)]      // singular: the second pivot is exactly 0
    [InlineData(1.0, 2.0, 1.0, 1, -3.0)]     // indefinite: the second pivot is 1 − 4
    [InlineData(-1.0, 0.0, 1.0, 0, -1.0)]    // negative first diagonal
    public void NonPositiveDefinite_StopsAtTheFailingPivot_WithoutFlooring(double a00, double a10, double a11, int failedAt, double pivot)
    {
        var f = new DenseCholesky(2);
        f.Decompose(new double[,] { { a00, a10 }, { a10, a11 } });
        Assert.True(f.Status.Finite);
        Assert.False(f.Status.PositiveDefinite);
        Assert.Equal(0, f.Status.FlooredPivots);
        Assert.Equal(failedAt, f.FailedAt);
        Assert.Equal(pivot, f.Status.MinPivot);
        Assert.True(double.IsNaN(f.LogDet));
        Assert.Throws<InvalidOperationException>(() => f.Solve([1.0, 1.0]));
        Assert.Throws<InvalidOperationException>(() => f.InverseQuadraticForm([1.0, 1.0]));
        Assert.Throws<InvalidOperationException>(() => f.SampleInnovation([1.0, 1.0], new double[2]));
    }

    [Fact]
    public void NonFiniteEntry_FailsWhereItEntersAPivot()
    {
        var f = new DenseCholesky(3);
        f.Decompose(new double[,] { { 2, 0, 0 }, { 0, 2, 0 }, { 0, double.NaN, 2 } });
        Assert.False(f.Status.Finite);
        Assert.False(f.Status.PositiveDefinite);
        Assert.Equal(2, f.FailedAt);
    }

    [Fact]
    public void TinyScaleIdentity_IsFactoredExactly()
    {
        // The floored factors clamp every pivot of 1e-20·I; the strict factor has no floor to hit.
        var f = new DenseCholesky(2);
        f.Decompose(new double[,] { { 1e-20, 0 }, { 0, 1e-20 } });
        Assert.True(f.Status.PositiveDefinite);
        Assert.Equal(0, f.Status.FlooredPivots);
        Assert.Equal(2.0 * Math.Log(1e-20), f.LogDet, 12);
        AssertClose([1e20, 2e20], f.Solve([1.0, 2.0]), 1e-15);
        Assert.Equal(1.0, f.Status.MinRelativePivot);
    }

    [Fact]
    public void IllConditionedHilbert_IsFactored_AndSolvedBackwardStably()
    {
        // H_ij = 1/(i + j + 1), κ(H₁₀) ≈ 1.6e13: positive definite, so the strict factor proceeds. The solve's
        // residual is small relative to ‖H‖‖x‖ even though x itself carries the conditioning.
        const int n = 10;
        var h = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) h[i, j] = 1.0 / (i + j + 1);
        var f = new DenseCholesky(n);
        f.Decompose(h);
        Assert.True(f.Status.PositiveDefinite);
        // The last column is nearly in the span of the others: its relative pivot is 4.2e-10.
        Assert.True(f.Status.MinRelativePivot < 1e-8, $"relative pivot {f.Status.MinRelativePivot:E3}");

        double[] b = Vector(n, new Xoshiro256PlusPlus(seed: 63));
        double[] x = f.Solve(b);
        double rMax = 0.0, xMax = 0.0;
        for (int i = 0; i < n; i++)
        {
            double r = -b[i];
            for (int j = 0; j < n; j++) r += h[i, j] * x[j];
            rMax = Math.Max(rMax, Math.Abs(r));
            xMax = Math.Max(xMax, Math.Abs(x[i]));
        }
        Assert.True(rMax <= 1e-12 * xMax * 2.0, $"residual {rMax:E3} for ‖x‖∞ {xMax:E3}");
    }

    [Fact]
    public void RelativePivot_OfTheLastColumn_IsTheSquaredSineOfItsAngleToTheOthers()
    {
        // Gram of [x₁, x₂, x₃] with x₃ = cos θ·u + sin θ·w, u in span(x₁, x₂) and w orthogonal to it: the last pivot
        // over its diagonal is sin²θ, the D06 admission statistic for a birth.
        const int m = 50;
        const double theta = 1e-3;
        var rng = new Xoshiro256PlusPlus(seed: 64);
        double[] x1 = Vector(m, rng), x2 = Vector(m, rng), w = Vector(m, rng);
        double[] u = new double[m];
        for (int i = 0; i < m; i++) u[i] = 0.6 * x1[i] - 1.3 * x2[i];
        double[] q1 = (double[])x1.Clone(), q2 = (double[])x2.Clone();   // orthonormal basis of span(x₁, x₂)
        Scale(q1, 1.0 / Norm(q1));
        Orthogonalize(q2, q1);
        Scale(q2, 1.0 / Norm(q2));
        Orthogonalize(w, q1);
        Orthogonalize(w, q2);
        Orthogonalize(w, q1);
        Scale(w, Norm(u) / Norm(w));
        var x3 = new double[m];
        for (int i = 0; i < m; i++) x3[i] = Math.Cos(theta) * u[i] + Math.Sin(theta) * w[i];

        double[][] cols = [x1, x2, x3];
        var g = new double[3, 3];
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) g[i, j] = Dot(cols[i], cols[j]);
        var f = new DenseCholesky(3);
        f.Decompose(g);

        double sin2 = Math.Sin(theta) * Math.Sin(theta);
        Assert.True(Math.Abs(f.Status.MinRelativePivot - sin2) <= 1e-6 * sin2, $"relative pivot {f.Status.MinRelativePivot:R} vs sin²θ {sin2:R}");
        Assert.Equal(f.Status.MinRelativePivot, f.RelativePivot(2));
    }

    public static TheoryData<string> AllFactorizations => new() { "strict", "dense", "band-chol", "band-ldlt" };

    [Theory]
    [MemberData(nameof(AllFactorizations))]
    public void RelativePivot_DoesNotChangeWhenColumnsAreRescaled(string kind)
    {
        // S·A·S has the pivots s_j²·d_j and the diagonal s_j²·A_jj, so d_j/A_jj is unchanged; the pivot ratio is not.
        const int n = 6;
        double[,] a = DenseSpd(n, new Xoshiro256PlusPlus(seed: 65));
        double[] s = [1e-4, 1.0, 30.0, 1e3, 0.2, 1e-2];
        var sas = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) sas[i, j] = s[i] * a[i, j] * s[j];

        FactorStatus plain = Status(kind, a), scaled = Status(kind, sas);
        Assert.True(Math.Abs(plain.MinRelativePivot - scaled.MinRelativePivot) <= 1e-12 * plain.MinRelativePivot,
                    $"{kind}: {plain.MinRelativePivot:R} vs {scaled.MinRelativePivot:R}");
        Assert.True(scaled.PivotRatio < 1e-6 * plain.PivotRatio, $"{kind}: the pivot ratio should move with the scaling");
    }

    [Theory]
    [MemberData(nameof(AllFactorizations))]
    public void RelativePivot_BoundsTheEquilibratedConditionNumber(string kind)
    {
        // Ã = D^(−½) A D^(−½) has unit diagonal, so λ_max(Ã) ≥ 1, and its pivots are d_j/A_jj ≥ λ_min(Ã):
        // κ(Ã) ≥ 1/MinRelativePivot.
        double[] lambda = [1e-7, 0.01, 0.3, 1.0, 2.5, 40.0];
        int n = lambda.Length;
        double[,] q = Householder(Vector(n, new Xoshiro256PlusPlus(seed: 66)));
        var d = new double[n, n];
        for (int i = 0; i < n; i++) d[i, i] = lambda[i];
        double[,] a = Multiply(q, Multiply(d, Transpose(q)));

        var eq = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) eq[i, j] = a[i, j] / Math.Sqrt(a[i, i] * a[j, j]);
        double[] ev = DenseEigen.DecomposeSymmetric(eq).Eigenvalues;
        double kappa = Max(ev) / Min(ev);

        FactorStatus st = Status(kind, a);
        Assert.True(1.0 / st.MinRelativePivot <= kappa * (1.0 + 1e-6), $"{kind}: 1/relative pivot {1.0 / st.MinRelativePivot:E3} > κ(Ã) {kappa:E3}");
        Assert.True(st.MinRelativePivot <= 1.0);
    }

    [Fact]
    public void Capacity_IsReused_AndBounded()
    {
        var f = new DenseCholesky(5);
        Assert.Equal(5, f.Capacity);
        Assert.Throws<InvalidOperationException>(() => f.Solve([1.0]));   // nothing factored yet
        f.Decompose(new double[,] { { 4.0 } });
        Assert.Equal(1, f.Dimension);
        AssertClose([0.25], f.Solve([1.0]), 1e-16);
        f.Decompose(DenseSpd(5, new Xoshiro256PlusPlus(seed: 67)));
        Assert.Equal(5, f.Dimension);
        Assert.Throws<ArgumentException>(() => f.Decompose(new double[6, 6]));
        Assert.Throws<ArgumentException>(() => f.Decompose(new double[2, 3]));
        Assert.Throws<ArgumentException>(() => f.Solve([1.0, 2.0]));
    }

    private static FactorStatus Status(string kind, double[,] a)
    {
        int n = a.GetLength(0);
        switch (kind)
        {
            case "strict":
                var strict = new DenseCholesky(n);
                strict.Decompose(a);
                return strict.Status;
            case "dense":
                var chol = new CholeskyDecomposition(n);
                chol.Decompose(a);
                return chol.Status;
            default:
                var band = new BandCholesky(n, n - 1, kind == "band-ldlt" ? BandFactorization.Ldlt : BandFactorization.Cholesky);
                band.Decompose(a);
                return band.Status;
        }
    }

    private static double Dot(double[] a, double[] b) { double s = 0.0; for (int i = 0; i < a.Length; i++) s += a[i] * b[i]; return s; }
    private static double Norm(double[] a) => Math.Sqrt(Dot(a, a));
    private static void Scale(double[] a, double c) { for (int i = 0; i < a.Length; i++) a[i] *= c; }
    private static double Max(double[] a) { double m = double.NegativeInfinity; foreach (double v in a) m = Math.Max(m, v); return m; }
    private static double Min(double[] a) { double m = double.PositiveInfinity; foreach (double v in a) m = Math.Min(m, v); return m; }

    // Gram–Schmidt twice for a clean orthogonal complement.
    private static void Orthogonalize(double[] w, double[] x)
    {
        for (int pass = 0; pass < 2; pass++)
        {
            double c = Dot(w, x) / Dot(x, x);
            for (int i = 0; i < w.Length; i++) w[i] -= c * x[i];
        }
    }
}
