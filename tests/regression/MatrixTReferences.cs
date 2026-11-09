using System;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;

namespace Maths.Regression.Tests;

/// <summary>
/// One case of E34's model for its references (the multi-output brief, slice 0): Y (n × K) = XA + E, each row of E ~ N(0, Σ),
/// the likelihood raised to the inverse temperature b, and A | Σ ~ MN(A₀, Λ₀⁻¹, Σ) on X's p columns. vec stacks A's columns, one
/// per output, so vec(A) ~ N(vec A₀, Σ ⊗ Λ₀⁻¹). The intercept is one of X's columns, inside the prior, unless
/// <see cref="FlatIntercept"/> sets a level μ ∈ ℝᴷ apart from X under a flat prior, integrated out (Liang et al.'s convention);
/// <see cref="PriorMean"/> null is zero.
/// </summary>
internal sealed record MultiOutputCase(double[,] Y, double[,] X, double[,] PriorPrecision, double[,]? PriorMean, double InverseTemperature,
                                       bool FlatIntercept = false)
{
    public int Rows => Y.GetLength(0);

    public int Outputs => Y.GetLength(1);

    public int Columns => X.GetLength(1);

    /// <summary>
    /// BASS's g-form on centred columns X_c (D03) carried to K outputs: the intercept inside the prior with Gram entry n, centred on ȳ
    /// (level centring), on <paramref name="levels"/> (a frozen level) or at 0 (the zero centre), or flat; Λ₀ = τ·blockdiag(n, X_cᵀX_c).
    /// </summary>
    public static MultiOutputCase GForm(double[,] y, double[,] xc, GPriorCentre centre, double tau, double b, double[]? levels = null)
    {
        int n = y.GetLength(0), k = y.GetLength(1), m = xc.GetLength(1);
        double[,] gram = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xc), xc);
        if (centre == GPriorCentre.Flat)
        {
            var flat = new double[m, m];
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++) flat[i, j] = tau * gram[i, j];
            return new MultiOutputCase(y, xc, flat, null, b, true);
        }
        var x = new double[n, m + 1];
        var lambda = new double[m + 1, m + 1];
        lambda[0, 0] = tau * n;
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = 1.0;
            for (int j = 0; j < m; j++) x[i, j + 1] = xc[i, j];
        }
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++) lambda[i + 1, j + 1] = tau * gram[i, j];
        if (centre == GPriorCentre.Zero) return new MultiOutputCase(y, x, lambda, null, b);
        var a0 = new double[m + 1, k];
        for (int l = 0; l < k; l++)
        {
            double mean = 0.0;
            for (int i = 0; i < n; i++) mean += y[i, l];
            a0[0, l] = levels?[l] ?? mean / n;
        }
        return new MultiOutputCase(y, x, lambda, a0, b);
    }
}

/// <summary>Σ's prior: IW(ν, Ψ), or with <see cref="Scale"/> null the improper |Σ|^(−(ν + K + 1)/2), Jeffreys' at ν = 0.</summary>
internal sealed record ReferenceSigmaPrior(double Nu, double[,]? Scale)
{
    /// <summary>
    /// The multi-output brief's decision 2 on the centre levels <paramref name="levels"/>: ν₀ = K − 1 + 2c and Ψ₀ = 2c·S_Y, S_Y the
    /// outputs' mean cross products about the levels (divisor n, as the scalar s² reads), so that K = 1 is IG(c, c·s²).
    /// </summary>
    public static ReferenceSigmaPrior ScaledToResponse(double[,] y, double[] levels, double c = 1e-3)
    {
        int n = y.GetLength(0), k = y.GetLength(1);
        var s = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++)
            {
                double sum = 0.0;
                for (int i = 0; i < n; i++) sum += (y[i, a] - levels[a]) * (y[i, b] - levels[b]);
                s[a, b] = 2.0 * c * sum / n;
            }
        return new ReferenceSigmaPrior(k - 1 + 2.0 * c, s);
    }
}

/// <summary>
/// A case with A integrated at Σ given: the count ν = nb − f (f = 1 under the flat intercept), the K × K form R and the coefficient
/// integral's determinant term, so that ln p_b(Y | Σ) = −(νK/2) ln 2π − (ν/2) ln|Σ| + LogDeterminantTerm − ½ tr(Σ⁻¹R).
/// </summary>
internal sealed record MatrixTCollapse(double[,] R, double LogDeterminantTerm, double Count, double InverseTemperature, bool FlatIntercept)
{
    public int Outputs => R.GetLength(0);
}

/// <summary>
/// E34's collapsed marginal by dense matrix algebra, the closed forms the multi-output brief's slice 0 gates: A integrated against
/// the matrix-normal prior, then Σ against its inverse Wishart. With P = Λ₀ + bXᵀX and Ā = P⁻¹(bXᵀY + Λ₀A₀),
/// R = b(Y − XĀ)ᵀ(Y − XĀ) + (Ā − A₀)ᵀΛ₀(Ā − A₀), a sum of positive semidefinite terms, and the determinant term is
/// (K/2)(ln|Λ₀| − ln|P|). The flat intercept centres Y and X, spends one observation and adds −(K/2) ln(nb). Then
/// Σ | Y ~ IW(ν₀ + ν, Ψ₀ + R), proper where ν₀ + ν &gt; K − 1, and
/// ln p_b(Y) = −(νK/2) ln π + LogDeterminantTerm + (ν₀/2) ln|Ψ₀| − ln Γ_K(ν₀/2) + ln Γ_K((ν₀ + ν)/2) − ((ν₀ + ν)/2) ln|Ψ₀ + R|.
/// At K = 1 it is the scalar collapsed g-prior's marginal under IG(ν₀/2, Ψ₀/2). Untempered, it is Brown, Vannucci and Fearn's
/// (1998) eqs. (19)–(21), whose Σ ∼ IW(δ; Q) is in Dawid's shape, δ = ν₀ − K + 1, and whose vague intercept α ∼ N(α₀, hΣ), h → ∞,
/// counts n observations as level centring does.
/// </summary>
internal static class MatrixTReference
{
    private static readonly double LogTwoPi = Math.Log(2.0 * Math.PI);

    /// <summary>ln Γ_K(a) = K(K − 1)/4 · ln π + Σⱼ ln Γ(a − j/2), j = 0 … K − 1.</summary>
    public static double LogMultivariateGamma(int k, double a)
    {
        double sum = 0.25 * k * (k - 1) * Math.Log(Math.PI);
        for (int j = 0; j < k; j++) sum += SpecialFunctions.LogGamma(a - 0.5 * j);
        return sum;
    }

    /// <summary>A integrated at Σ given (<see cref="MatrixTCollapse"/>).</summary>
    public static MatrixTCollapse Collapse(MultiOutputCase c)
    {
        int n = c.Rows, k = c.Outputs, p = c.Columns;
        double b = c.InverseTemperature;
        if (!(b >= 0.0) || double.IsPositiveInfinity(b)) throw new ArgumentOutOfRangeException(nameof(c), "b must be non-negative and finite.");
        if (c.FlatIntercept && !(b > 0.0)) throw new ArgumentOutOfRangeException(nameof(c), "The flat intercept's integral diverges at b = 0.");
        double[,] y = c.FlatIntercept ? ReferenceAlgebra.Centred(c.Y) : c.Y, x = c.FlatIntercept ? ReferenceAlgebra.Centred(c.X) : c.X;
        double[,] lambda = c.PriorPrecision, a0 = c.PriorMean ?? new double[p, k];

        var precision = new double[p, p];
        var right = new double[p, k];
        for (int i = 0; i < p; i++)
        {
            for (int j = 0; j < p; j++)
            {
                double sum = 0.0;
                for (int r = 0; r < n; r++) sum += x[r, i] * x[r, j];
                precision[i, j] = lambda[i, j] + b * sum;
            }
            for (int l = 0; l < k; l++)
            {
                double sum = 0.0, prior = 0.0;
                for (int r = 0; r < n; r++) sum += x[r, i] * y[r, l];
                for (int j = 0; j < p; j++) prior += lambda[i, j] * a0[j, l];
                right[i, l] = b * sum + prior;
            }
        }
        double[,] factor = ReferenceAlgebra.Cholesky(precision) ?? throw new ArgumentException("Λ₀ + bXᵀX is not positive definite.", nameof(c));
        double[,] lambdaFactor = ReferenceAlgebra.Cholesky(lambda) ?? throw new ArgumentException("Λ₀ is not positive definite.", nameof(c));
        double[,] mean = ReferenceAlgebra.Solve(factor, right);

        var residual = new double[n, k];
        for (int r = 0; r < n; r++)
            for (int l = 0; l < k; l++)
            {
                double fit = 0.0;
                for (int j = 0; j < p; j++) fit += x[r, j] * mean[j, l];
                residual[r, l] = y[r, l] - fit;
            }
        var deviation = new double[p, k];
        for (int j = 0; j < p; j++)
            for (int l = 0; l < k; l++) deviation[j, l] = mean[j, l] - a0[j, l];
        var form = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int d = 0; d <= a; d++)
            {
                double sum = 0.0, prior = 0.0;
                for (int r = 0; r < n; r++) sum += residual[r, a] * residual[r, d];
                for (int i = 0; i < p; i++)
                    for (int j = 0; j < p; j++) prior += deviation[i, a] * lambda[i, j] * deviation[j, d];
                form[a, d] = form[d, a] = b * sum + prior;
            }

        double logDet = 0.5 * k * (ReferenceAlgebra.LogDeterminant(lambdaFactor) - ReferenceAlgebra.LogDeterminant(factor));
        if (c.FlatIntercept) logDet -= 0.5 * k * Math.Log(n * b);
        return new MatrixTCollapse(form, logDet, n * b - (c.FlatIntercept ? 1.0 : 0.0), b, c.FlatIntercept);
    }

    /// <summary>ln p_b(Y | Σ) with A integrated; exactly 0 at b = 0, where the likelihood is gone.</summary>
    public static double LogLikelihoodGiven(MatrixTCollapse collapse, double[,] sigma)
    {
        int k = collapse.Outputs;
        if (sigma.GetLength(0) != k || sigma.GetLength(1) != k) throw new ArgumentException("Σ must be K × K.", nameof(sigma));
        if (collapse.InverseTemperature == 0.0) return 0.0;
        double[,] factor = ReferenceAlgebra.Cholesky(sigma) ?? throw new ArgumentException("Σ is not positive definite.", nameof(sigma));
        double nu = collapse.Count;
        return -0.5 * nu * k * LogTwoPi - 0.5 * nu * ReferenceAlgebra.LogDeterminant(factor) + collapse.LogDeterminantTerm
               - 0.5 * ReferenceAlgebra.TraceSolve(factor, collapse.R);
    }

    /// <summary>
    /// ln p_b(Y) with A and Σ integrated. Proper prior: absolute, and exactly 0 at b = 0. Improper prior: up to its constant, with
    /// ln p = −(νK/2) ln 2π + LogDeterminantTerm + ((ν₀ + ν)K/2) ln 2 + ln Γ_K((ν₀ + ν)/2) − ((ν₀ + ν)/2) ln|R|. Refused where Σ's
    /// posterior is improper, ν₀ + ν ≤ K − 1, and for a proper prior with ν₀ ≤ K − 1.
    /// </summary>
    public static double LogMarginal(MatrixTCollapse collapse, ReferenceSigmaPrior prior)
    {
        int k = collapse.Outputs;
        double nu = collapse.Count, nu0 = prior.Nu, posterior = nu0 + nu;
        if (!(posterior > k - 1)) throw new InvalidOperationException($"Σ's posterior is improper: ν₀ + ν = {posterior:R} ≤ K − 1 = {k - 1}.");
        if (prior.Scale is not { } psi)
        {
            double[,] rFactor = ReferenceAlgebra.Cholesky(collapse.R) ?? throw new InvalidOperationException("The improper prior needs R positive definite.");
            return -0.5 * nu * k * LogTwoPi + collapse.LogDeterminantTerm + 0.5 * posterior * k * Math.Log(2.0)
                   + LogMultivariateGamma(k, 0.5 * posterior) - 0.5 * posterior * ReferenceAlgebra.LogDeterminant(rFactor);
        }
        if (!(nu0 > k - 1)) throw new ArgumentException($"IW(ν₀, Ψ₀) is improper at ν₀ = {nu0:R} ≤ K − 1.", nameof(prior));
        double[,] psiFactor = ReferenceAlgebra.Cholesky(psi) ?? throw new ArgumentException("Ψ₀ is not positive definite.", nameof(prior));
        if (collapse.InverseTemperature == 0.0) return 0.0;   // no likelihood, and a proper prior integrates to one
        double[,] total = ReferenceAlgebra.Sum(psi, collapse.R);
        double[,] totalFactor = ReferenceAlgebra.Cholesky(total) ?? throw new InvalidOperationException("Ψ₀ + R is not positive definite.");
        return -0.5 * nu * k * Math.Log(Math.PI) + collapse.LogDeterminantTerm + 0.5 * nu0 * ReferenceAlgebra.LogDeterminant(psiFactor)
               - LogMultivariateGamma(k, 0.5 * nu0) + LogMultivariateGamma(k, 0.5 * posterior)
               - 0.5 * posterior * ReferenceAlgebra.LogDeterminant(totalFactor);
    }

    /// <summary>Σ | Y ~ IW(ν₀ + ν, Ψ₀ + R): its degrees of freedom and scale (R alone under an improper prior).</summary>
    public static (double Nu, double[,] Scale) SigmaPosterior(MatrixTCollapse collapse, ReferenceSigmaPrior prior)
        => (prior.Nu + collapse.Count, prior.Scale is { } psi ? ReferenceAlgebra.Sum(psi, collapse.R) : (double[,])collapse.R.Clone());

    /// <summary>
    /// The matrix-t density of Y (n × K) that Y | Σ ~ MN(M, V, Σ) and Σ ~ IW(ν₀, Ψ) give, stated in the inverse Wishart's ν₀: the
    /// matrix-t as Brown, Vannucci and Fearn (1998, §2) define it, this mixture in Dawid's notation. The density is derived here by
    /// the integral; its textbook forms, Dickey 1967 and Gupta and Nagar 2000, ch. 4, are cited from memory:
    /// ln Γ_K((ν₀ + n)/2) − ln Γ_K(ν₀/2) − (nK/2) ln π − (K/2) ln|V| + (ν₀/2) ln|Ψ| − ((ν₀ + n)/2) ln|Ψ + (Y − M)ᵀV⁻¹(Y − M)|.
    /// E34's untempered marginal with the intercept inside the prior is the case M = XA₀, V = I + XΛ₀⁻¹Xᵀ, which reads n × n
    /// matrices where the collapse reads p × p and K × K ones.
    /// </summary>
    public static double MatrixTLogDensity(double[,] y, double[,] location, double[,] rowScale, double[,] columnScale, double nu0)
    {
        int n = y.GetLength(0), k = y.GetLength(1);
        double[,] vFactor = ReferenceAlgebra.Cholesky(rowScale) ?? throw new ArgumentException("V is not positive definite.", nameof(rowScale));
        double[,] psiFactor = ReferenceAlgebra.Cholesky(columnScale) ?? throw new ArgumentException("Ψ is not positive definite.", nameof(columnScale));
        var e = new double[n, k];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < k; l++) e[i, l] = y[i, l] - location[i, l];
        double[,] solved = ReferenceAlgebra.Solve(vFactor, e);
        var total = (double[,])columnScale.Clone();
        for (int a = 0; a < k; a++)
            for (int d = 0; d < k; d++)
            {
                double sum = 0.0;
                for (int i = 0; i < n; i++) sum += e[i, a] * solved[i, d];
                total[a, d] += sum;
            }
        ReferenceAlgebra.Symmetrize(total);
        double[,] totalFactor = ReferenceAlgebra.Cholesky(total) ?? throw new InvalidOperationException("Ψ + EᵀV⁻¹E is not positive definite.");
        return LogMultivariateGamma(k, 0.5 * (nu0 + n)) - LogMultivariateGamma(k, 0.5 * nu0) - 0.5 * n * k * Math.Log(Math.PI)
               - 0.5 * k * ReferenceAlgebra.LogDeterminant(vFactor) + 0.5 * nu0 * ReferenceAlgebra.LogDeterminant(psiFactor)
               - 0.5 * (nu0 + n) * ReferenceAlgebra.LogDeterminant(totalFactor);
    }
}

/// <summary>Which Kronecker product the dense vec form places as vec(A)'s prior covariance under vec's column stacking.</summary>
internal enum KroneckerOrder
{
    /// <summary>Σ ⊗ Λ₀⁻¹, the model's: block (k, l) is Σ_kl·Λ₀⁻¹.</summary>
    SigmaFirst,

    /// <summary>Λ₀⁻¹ ⊗ Σ in the same slot, the control a swapped convention would build.</summary>
    Swapped,
}

/// <summary>
/// ln p_b(Y | Σ) with A integrated as one dense Gaussian on vec(Y), stacking Y's columns, built entry by entry: vec(Y) = (I_K ⊗ X) vec(A)
/// + vec(E), with Cov(vec A) = Σ ⊗ Λ₀⁻¹ and Cov(vec E) = Σ ⊗ I/b, so Cov(vec Y) = (I_K ⊗ X)(Σ ⊗ Λ₀⁻¹)(I_K ⊗ X)ᵀ + Σ ⊗ I/b, formed by
/// explicit products rather than the Kronecker identities. The tempering is p(Y | A, Σ)^b = N(vec Y; vec XA, Σ ⊗ I/b)·
/// (2π)^(−nK(b − 1)/2)·b^(−nK/2)·|Σ|^(−n(b − 1)/2). The flat intercept rotates the rows by Helmert's orthonormal basis of 1's
/// complement, whose first row carries √n·μ; integrating μ leaves n^(−K/2) and the other n − 1 rows.
/// </summary>
internal static class DenseVecReference
{
    public static double LogLikelihoodGiven(MultiOutputCase c, double[,] sigma, KroneckerOrder order = KroneckerOrder.SigmaFirst)
    {
        int n = c.Rows, k = c.Outputs, p = c.Columns;
        double b = c.InverseTemperature;
        if (!(b > 0.0)) throw new ArgumentOutOfRangeException(nameof(c), "The dense form reads Σ/b, so it needs b > 0.");
        double[,] y = c.Y, x = c.X;
        int m = n;
        if (c.FlatIntercept)
        {
            double[,] helmert = Helmert(n);
            y = ReferenceAlgebra.Multiply(helmert, y);
            x = ReferenceAlgebra.Multiply(helmert, x);
            m = n - 1;
        }
        double[,] a0 = c.PriorMean ?? new double[p, k];
        double[,] priorCovariance = ReferenceAlgebra.Inverse(c.PriorPrecision);
        double[,] coefficientCovariance = order == KroneckerOrder.SigmaFirst
            ? ReferenceAlgebra.Kronecker(sigma, priorCovariance)
            : ReferenceAlgebra.Kronecker(priorCovariance, sigma);

        // I_K ⊗ X, row (l, r) at l·m + r, column (l, j) at l·p + j.
        var lift = new double[m * k, p * k];
        for (int l = 0; l < k; l++)
            for (int r = 0; r < m; r++)
                for (int j = 0; j < p; j++) lift[l * m + r, l * p + j] = x[r, j];
        double[,] covariance = ReferenceAlgebra.Multiply(ReferenceAlgebra.Multiply(lift, coefficientCovariance), ReferenceAlgebra.Transpose(lift));
        var noise = new double[m, m];
        for (int r = 0; r < m; r++) noise[r, r] = 1.0 / b;
        double[,] noiseCovariance = ReferenceAlgebra.Kronecker(sigma, noise);
        var deviation = new double[m * k];
        for (int l = 0; l < k; l++)
            for (int r = 0; r < m; r++)
            {
                double fit = 0.0;
                for (int j = 0; j < p; j++) fit += x[r, j] * a0[j, l];
                deviation[l * m + r] = y[r, l] - fit;
                for (int s = 0; s < m * k; s++) covariance[l * m + r, s] += noiseCovariance[l * m + r, s];
            }
        ReferenceAlgebra.Symmetrize(covariance);
        double[,] factor = ReferenceAlgebra.Cholesky(covariance) ?? throw new ArgumentException("Cov(vec Y) is not positive definite.");
        double[] whitened = ReferenceAlgebra.ForwardSolve(factor, deviation);
        double quadratic = 0.0;
        foreach (double w in whitened) quadratic += w * w;
        double gaussian = -0.5 * m * k * Math.Log(2.0 * Math.PI) - 0.5 * ReferenceAlgebra.LogDeterminant(factor) - 0.5 * quadratic;

        double[,] sigmaFactor = ReferenceAlgebra.Cholesky(sigma) ?? throw new ArgumentException("Σ is not positive definite.", nameof(sigma));
        double tempering = -0.5 * n * k * (b - 1.0) * Math.Log(2.0 * Math.PI) - 0.5 * n * k * Math.Log(b)
                           - 0.5 * n * (b - 1.0) * ReferenceAlgebra.LogDeterminant(sigmaFactor);
        return gaussian + tempering - (c.FlatIntercept ? 0.5 * k * Math.Log(n) : 0.0);
    }

    // The n − 1 rows of Helmert's matrix after its first: row j (1-based) is (1, …, 1, −j, 0, …)/√(j(j + 1)) with j ones.
    private static double[,] Helmert(int n)
    {
        var h = new double[n - 1, n];
        for (int j = 1; j < n; j++)
        {
            double scale = 1.0 / Math.Sqrt(j * (j + 1.0));
            for (int i = 0; i < j; i++) h[j - 1, i] = scale;
            h[j - 1, j] = -j * scale;
        }
        return h;
    }
}

/// <summary>
/// The inverse Wishart IW(ν, Ψ) on K × K matrices: density |Ψ|^(ν/2) / (2^(νK/2) Γ_K(ν/2)) · |Σ|^(−(ν + K + 1)/2) · exp(−½ tr ΨΣ⁻¹), mean
/// Ψ/(ν − K − 1) for ν &gt; K + 1, and draws by Bartlett's decomposition of the Wishart W ~ W(ν, Ψ⁻¹) inverted, Σ = W⁻¹. At K = 1 it is
/// IG(ν/2, Ψ/2).
/// </summary>
internal static class ReferenceInverseWishart
{
    public static double LogDensity(double[,] sigma, double nu, double[,] psi)
    {
        int k = sigma.GetLength(0);
        double[,] sigmaFactor = ReferenceAlgebra.Cholesky(sigma) ?? throw new ArgumentException("Σ is not positive definite.", nameof(sigma));
        double[,] psiFactor = ReferenceAlgebra.Cholesky(psi) ?? throw new ArgumentException("Ψ is not positive definite.", nameof(psi));
        return 0.5 * nu * ReferenceAlgebra.LogDeterminant(psiFactor) - 0.5 * nu * k * Math.Log(2.0) - MatrixTReference.LogMultivariateGamma(k, 0.5 * nu)
               - 0.5 * (nu + k + 1) * ReferenceAlgebra.LogDeterminant(sigmaFactor) - 0.5 * ReferenceAlgebra.TraceSolve(sigmaFactor, psi);
    }

    public static double[,] Sample(double nu, double[,] psi, Xoshiro256PlusPlus rng) => SampleFromRoot(nu, SamplingRoot(psi), rng);

    /// <summary>The lower Cholesky factor of Ψ⁻¹, which every draw reads.</summary>
    public static double[,] SamplingRoot(double[,] psi)
        => ReferenceAlgebra.Cholesky(ReferenceAlgebra.Inverse(psi)) ?? throw new ArgumentException("Ψ is not positive definite.", nameof(psi));

    /// <summary>A draw with Ψ⁻¹'s factor <paramref name="root"/> taken once.</summary>
    public static double[,] SampleFromRoot(double nu, double[,] root, Xoshiro256PlusPlus rng)
    {
        int k = root.GetLength(0);
        var bartlett = new double[k, k];
        for (int i = 0; i < k; i++)
        {
            bartlett[i, i] = Math.Sqrt(2.0 * Gamma.SampleStandard(rng, 0.5 * (nu - i)));
            for (int j = 0; j < i; j++) bartlett[i, j] = Normal.SampleStandard(rng);
        }
        // Σ = W⁻¹ = NᵀN with N = (RA)⁻¹, triangular: a tiny χ² pivot (small ν) leaves W singular in doubles, but not RA.
        double[,] inverse = ReferenceAlgebra.LowerInverse(ReferenceAlgebra.Multiply(root, bartlett));
        double[,] sigma = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(inverse), inverse);
        ReferenceAlgebra.Symmetrize(sigma);
        return sigma;
    }
}

/// <summary>
/// ∫ f(Σ) dΣ over 2 × 2 positive definite Σ by the tensor trapezoid rule in Bartlett's coordinates of the precision: Σ⁻¹ = (RA)(RA)ᵀ,
/// with R the lower Cholesky factor of a reference matrix S and A lower triangular, A₁₁ = e^θ₀, A₂₁ = θ₁, A₂₂ = e^θ₂. The Jacobian
/// is |Σ⁻¹|^(−3)·|R|³·4A₁₁³A₂₂² = 4|R|^(−3)A₁₁^(−3)A₂₂^(−4): inversion, the congruence by R and the triangular square, with the
/// log diagonal. Under an inverse Wishart IW(ν, S⁻¹) the coordinates are independent, A₁₁² ~ χ²(ν) and A₂₂² ~ χ²(ν − 1) and A₂₁
/// standard normal; in ln A_ii a χ²(k) tail falls as e^(k·θ) below the peak and faster than any exponential above it, and A₂₁'s is
/// Gaussian, so a box holds them. <see cref="Adaptive"/> reads a first pass's E[Σ] and E[Σ⁻¹] as an inverse Wishart's, ν̂ = 3r/(r − 1)
/// with r = tr(E[Σ⁻¹]E[Σ])/2 and Ψ̂ = ν̂·E[Σ⁻¹]⁻¹, and integrates again on S = Ψ̂⁻¹; nothing of a closed form enters, and an integrand
/// that is not an inverse Wishart only gets coordinates less independent. On each line the rule is spectrally accurate for an
/// integrand analytic in a strip and negligible at the box's faces: <see cref="Result.BoundaryGap"/> reports how far below the peak
/// the faces sit, and a coarser grid is the check on the rest.
/// </summary>
internal static class SigmaQuadrature
{
    /// <summary>ln f at Σ = [[s₁₁, s₂₁], [s₂₁, s₂₂]].</summary>
    internal delegate double LogIntegrand(double s11, double s21, double s22);

    internal sealed record Result(double LogIntegral, double BoundaryGap, double[,] Mean, double[,] MeanInverse);

    /// <summary>−(power/2)·ln|Σ| − ½ tr(Σ⁻¹M) at 2 × 2 Σ and M: an inverse Wishart's kernel, and the collapse's at Σ given.</summary>
    public static double Kernel(double s11, double s21, double s22, double[,] m, double power)
    {
        double det = s11 * s22 - s21 * s21;
        return -0.5 * power * Math.Log(det) - 0.5 * (s22 * m[0, 0] - 2.0 * s21 * m[1, 0] + s11 * m[1, 1]) / det;
    }

    /// <summary>
    /// The second pass's rule, its coarser check and the first pass, from the reference <paramref name="first"/>: the first over
    /// θ ∈ peak ± (8, 12, 8), the second over ln A_ii from 44/(ν̂ − 1) + ½ below the peak to 2.5 above it and A₂₁ within ±9.5.
    /// </summary>
    public static (Result Fine, Result Check, Result First) Adaptive(LogIntegrand logF, double[,] first)
    {
        double[] centre = Peak(logF, first, [1.0, 0.0, 1.0], [3.0, 4.0, 3.0], 25);
        Result pass = Integrate(logF, first, [centre[0] - 8.0, centre[1] - 12.0, centre[2] - 8.0], [centre[0] + 8.0, centre[1] + 12.0, centre[2] + 8.0], [0.15, 0.5, 0.15]);
        double[,] product = ReferenceAlgebra.Multiply(pass.MeanInverse, pass.Mean);
        double r = 0.5 * (product[0, 0] + product[1, 1]), nu = 3.0 * r / (r - 1.0);
        double[,] reference = ReferenceAlgebra.Multiply(pass.MeanInverse, new[,] { { 1.0 / nu, 0.0 }, { 0.0, 1.0 / nu } });
        ReferenceAlgebra.Symmetrize(reference);
        double[] peak = Peak(logF, reference, [0.5 * Math.Log(nu), 0.0, 0.5 * Math.Log(Math.Max(nu - 1.0, 0.5))], [1.5, 1.5, 1.5], 31);
        // ln A₂₂'s χ²(k = ν̂ − 1) density falls k·d − (k/2)(1 − e^(−2d)) at d below its peak: 44 at d = 44/k + ½.
        double below = Math.Max(2.5, 44.0 / Math.Max(nu - 1.0, 0.5) + 0.5);
        double[] lower = [peak[0] - below, peak[1] - 9.5, peak[2] - below], upper = [peak[0] + 2.5, peak[1] + 9.5, peak[2] + 2.5];
        return (Integrate(logF, reference, lower, upper, [0.06, 0.3, 0.06]), Integrate(logF, reference, lower, upper, [0.08, 0.4, 0.08]), pass);
    }

    /// <summary>The argmax of <paramref name="logF"/>, with its Jacobian, on a grid of <paramref name="points"/>³ over centre ± halfWidths.</summary>
    public static double[] Peak(LogIntegrand logF, double[,] reference, double[] centre, double[] halfWidths, int points)
    {
        double[,] root = ReferenceAlgebra.Cholesky(reference) ?? throw new ArgumentException("The reference is not positive definite.", nameof(reference));
        double best = double.NegativeInfinity;
        double[] arg = centre;
        for (int a = 0; a < points; a++)
            for (int t = 0; t < points; t++)
                for (int c = 0; c < points; c++)
                {
                    double[] theta = [Grid(centre, halfWidths, 0, a, points), Grid(centre, halfWidths, 1, t, points), Grid(centre, halfWidths, 2, c, points)];
                    double v = Evaluate(logF, root, theta, out _, out _, out _);
                    if (v > best) { best = v; arg = theta; }
                }
        return arg;

        static double Grid(double[] centre, double[] widths, int d, int i, int points) => centre[d] - widths[d] + 2.0 * widths[d] * i / (points - 1);
    }

    /// <summary>ln ∫ f dΣ over θ in [lower, upper] at the steps <paramref name="steps"/>, with E[Σ] and E[Σ⁻¹] under f.</summary>
    public static Result Integrate(LogIntegrand logF, double[,] reference, double[] lower, double[] upper, double[] steps)
    {
        double[,] root = ReferenceAlgebra.Cholesky(reference) ?? throw new ArgumentException("The reference is not positive definite.", nameof(reference));
        var count = new int[3];
        for (int d = 0; d < 3; d++) count[d] = (int)Math.Ceiling((upper[d] - lower[d]) / steps[d]) + 1;
        double max = double.NegativeInfinity, boundary = double.NegativeInfinity, sum = 0.0;
        var moment = new double[6];   // Σ₁₁, Σ₂₁, Σ₂₂ and Σ⁻¹'s, weighted, on the running scale of sum
        var theta = new double[3];
        for (int a = 0; a < count[0]; a++)
            for (int t = 0; t < count[1]; t++)
                for (int c = 0; c < count[2]; c++)
                {
                    theta[0] = lower[0] + a * steps[0];
                    theta[1] = lower[1] + t * steps[1];
                    theta[2] = lower[2] + c * steps[2];
                    double v = Evaluate(logF, root, theta, out double s11, out double s21, out double s22);
                    if (a == 0 || t == 0 || c == 0 || a == count[0] - 1 || t == count[1] - 1 || c == count[2] - 1) boundary = Math.Max(boundary, v);
                    // Unit weights at the faces too: the trapezoid's halves there sit below its resolution by the gap.
                    if (v > max)
                    {
                        double scale = Math.Exp(max - v);
                        sum *= scale;
                        for (int i = 0; i < moment.Length; i++) moment[i] *= scale;
                        max = v;
                    }
                    double w = Math.Exp(v - max), det = s11 * s22 - s21 * s21;
                    sum += w;
                    moment[0] += w * s11;
                    moment[1] += w * s21;
                    moment[2] += w * s22;
                    moment[3] += w * s22 / det;
                    moment[4] -= w * s21 / det;
                    moment[5] += w * s11 / det;
                }
        var mean = new double[2, 2];
        var inverse = new double[2, 2];
        mean[0, 0] = moment[0] / sum;
        mean[1, 0] = mean[0, 1] = moment[1] / sum;
        mean[1, 1] = moment[2] / sum;
        inverse[0, 0] = moment[3] / sum;
        inverse[1, 0] = inverse[0, 1] = moment[4] / sum;
        inverse[1, 1] = moment[5] / sum;
        return new Result(max + Math.Log(sum) + Math.Log(steps[0] * steps[1] * steps[2]), max - boundary, mean, inverse);
    }

    // ln f + ln Jacobian at θ, with Σ = NᵀN for N = (RA)⁻¹.
    private static double Evaluate(LogIntegrand logF, double[,] root, double[] theta, out double s11, out double s21, out double s22)
    {
        double a11 = Math.Exp(theta[0]), a21 = theta[1], a22 = Math.Exp(theta[2]);
        double m11 = root[0, 0] * a11, m21 = root[1, 0] * a11 + root[1, 1] * a21, m22 = root[1, 1] * a22;
        double n11 = 1.0 / m11, n21 = -m21 / (m11 * m22), n22 = 1.0 / m22;
        s11 = n11 * n11 + n21 * n21;
        s21 = n21 * n22;
        s22 = n22 * n22;
        return logF(s11, s21, s22) + Math.Log(4.0) - 3.0 * Math.Log(root[0, 0] * root[1, 1]) - 3.0 * theta[0] - 4.0 * theta[2];
    }
}

/// <summary>Dense helpers for the references, apart from the engine's factors: Cholesky, its solves and log-determinant, products.</summary>
internal static class ReferenceAlgebra
{
    /// <summary>The lower Cholesky factor of a symmetric matrix, its lower triangle read; null when a pivot is not positive.</summary>
    public static double[,]? Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = a[i, j];
                for (int r = 0; r < j; r++) sum -= l[i, r] * l[j, r];
                if (i == j)
                {
                    if (!(sum > 0.0) || double.IsPositiveInfinity(sum)) return null;
                    l[i, i] = Math.Sqrt(sum);
                }
                else l[i, j] = sum / l[j, j];
            }
        return l;
    }

    public static double LogDeterminant(double[,] factor)
    {
        double sum = 0.0;
        for (int i = 0; i < factor.GetLength(0); i++) sum += Math.Log(factor[i, i]);
        return 2.0 * sum;
    }

    public static double[] ForwardSolve(double[,] factor, double[] b)
    {
        int n = b.Length;
        var z = new double[n];
        for (int i = 0; i < n; i++)
        {
            double sum = b[i];
            for (int r = 0; r < i; r++) sum -= factor[i, r] * z[r];
            z[i] = sum / factor[i, i];
        }
        return z;
    }

    /// <summary>A⁻¹B for A = LLᵀ, column by column.</summary>
    public static double[,] Solve(double[,] factor, double[,] b)
    {
        int n = factor.GetLength(0), m = b.GetLength(1);
        var x = new double[n, m];
        var column = new double[n];
        for (int c = 0; c < m; c++)
        {
            for (int i = 0; i < n; i++) column[i] = b[i, c];
            double[] z = ForwardSolve(factor, column);
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = z[i];
                for (int r = i + 1; r < n; r++) sum -= factor[r, i] * x[r, c];
                x[i, c] = sum / factor[i, i];
            }
        }
        return x;
    }

    /// <summary>tr(A⁻¹B) for A = LLᵀ: ‖L⁻¹B^½‖² read as Σᵢ (A⁻¹B)ᵢᵢ.</summary>
    public static double TraceSolve(double[,] factor, double[,] b)
    {
        double[,] x = Solve(factor, b);
        double trace = 0.0;
        for (int i = 0; i < x.GetLength(0); i++) trace += x[i, i];
        return trace;
    }

    /// <summary>L⁻¹ for a lower triangular L with a nonzero diagonal, by forward substitution.</summary>
    public static double[,] LowerInverse(double[,] l)
    {
        int n = l.GetLength(0);
        var t = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            t[i, i] = 1.0 / l[i, i];
            for (int j = 0; j < i; j++)
            {
                double sum = 0.0;
                for (int k = j; k < i; k++) sum -= l[i, k] * t[k, j];
                t[i, j] = sum / l[i, i];
            }
        }
        return t;
    }

    public static double[,] Inverse(double[,] a)
    {
        int n = a.GetLength(0);
        double[,] factor = Cholesky(a) ?? throw new ArgumentException("The matrix is not positive definite.", nameof(a));
        var identity = new double[n, n];
        for (int i = 0; i < n; i++) identity[i, i] = 1.0;
        double[,] inverse = Solve(factor, identity);
        Symmetrize(inverse);
        return inverse;
    }

    public static double[,] Multiply(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), m = a.GetLength(1), q = b.GetLength(1);
        var c = new double[n, q];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < q; j++)
            {
                double sum = 0.0;
                for (int r = 0; r < m; r++) sum += a[i, r] * b[r, j];
                c[i, j] = sum;
            }
        return c;
    }

    public static double[,] Transpose(double[,] a)
    {
        int n = a.GetLength(0), m = a.GetLength(1);
        var t = new double[m, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) t[j, i] = a[i, j];
        return t;
    }

    public static double[,] Sum(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), m = a.GetLength(1);
        var c = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) c[i, j] = a[i, j] + b[i, j];
        return c;
    }

    /// <summary>A ⊗ B: block (i, j) is A_ij·B.</summary>
    public static double[,] Kronecker(double[,] a, double[,] b)
    {
        int ar = a.GetLength(0), ac = a.GetLength(1), br = b.GetLength(0), bc = b.GetLength(1);
        var k = new double[ar * br, ac * bc];
        for (int i = 0; i < ar; i++)
            for (int j = 0; j < ac; j++)
                for (int r = 0; r < br; r++)
                    for (int s = 0; s < bc; s++) k[i * br + r, j * bc + s] = a[i, j] * b[r, s];
        return k;
    }

    /// <summary>The columns less their means.</summary>
    public static double[,] Centred(double[,] a)
    {
        int n = a.GetLength(0), m = a.GetLength(1);
        var c = new double[n, m];
        for (int j = 0; j < m; j++)
        {
            double mean = 0.0;
            for (int i = 0; i < n; i++) mean += a[i, j];
            mean /= n;
            for (int i = 0; i < n; i++) c[i, j] = a[i, j] - mean;
        }
        return c;
    }

    public static void Symmetrize(double[,] a)
    {
        int n = a.GetLength(0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < i; j++) a[i, j] = a[j, i] = 0.5 * (a[i, j] + a[j, i]);
    }
}
