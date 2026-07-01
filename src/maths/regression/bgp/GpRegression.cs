using System;
using Maths.LinAlg;

namespace Maths.Regression.Bgp;

/// <summary>A GP fit at one bandwidth: the marginal evidence and the posterior weights α = (Σ_t+σ²I)⁻¹Y that
/// produce the conditional-mean prediction. Returned by <see cref="GpRegression.Fit"/>.</summary>
/// <param name="Bandwidth">The kernel bandwidth t this fit was computed at.</param>
/// <param name="LogMarginal">log p(Y | t), the GP marginal log-likelihood (the t-posterior's likelihood factor).</param>
/// <param name="Alpha">Posterior weights α = (Σ_t+σ²I)⁻¹Y; the conditional mean is k_t(·,X)·α.</param>
public sealed record GpFit(double Bandwidth, double LogMarginal, double[] Alpha);

/// <summary>
/// Bayesian GP regression with an ambient radial kernel (Tang, Wu, Cheng &amp; Dunson 2025): the conjugate core
/// behind the dimension-adaptive method. With the regression function marginalized out by GP conjugacy, all that
/// remains is the bandwidth t — and this class supplies the two quantities a t-sampler needs: the marginal
/// evidence log p(Y|t) (reduce data+model → the t-likelihood) and the posterior weights α for the conditional
/// mean. It also exposes the kernel-affinity observable v̂_n(t) (the dimension-implicit statistic the empirical
/// Bayes prior is built from). The linear algebra is dense (an n×n SPD Cholesky), unlike the banded free-knot
/// world; pairwise squared distances are t-independent, so they are computed once and reused across every fit.
/// </summary>
public sealed class GpRegression
{
    private readonly double[,] _x;        // n × D train predictors (ambient coordinates)
    private readonly double[] _y;         // n responses
    private readonly double _sigma2;      // known noise variance σ²
    private readonly IGpKernel _kernel;
    private readonly int _n;
    private readonly int _d;              // ambient dimension D
    private readonly double[,] _sqDist;   // n × n precomputed ‖X_i − X_j‖² (independent of t)

    /// <param name="x">n × D matrix of ambient predictor coordinates.</param>
    /// <param name="y">Length-n responses.</param>
    /// <param name="sigma2">Known observation-noise variance σ² &gt; 0.</param>
    /// <param name="kernel">Radial covariance kernel (squared-exponential is the canonical choice).</param>
    public GpRegression(double[,] x, double[] y, double sigma2, IGpKernel kernel)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(kernel);
        if (!(sigma2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(sigma2), "σ² must be positive.");
        _n = x.GetLength(0);
        _d = x.GetLength(1);
        if (y.Length != _n) throw new ArgumentException("y length must match the number of rows of x.", nameof(y));
        _x = x; _y = y; _sigma2 = sigma2; _kernel = kernel;

        _sqDist = new double[_n, _n];
        for (int i = 0; i < _n; i++)
            for (int j = i + 1; j < _n; j++)
            {
                double s = 0.0;
                for (int k = 0; k < _d; k++) { double diff = x[i, k] - x[j, k]; s += diff * diff; }
                _sqDist[i, j] = s;
                _sqDist[j, i] = s;
            }
    }

    public int Count => _n;

    /// <summary>
    /// v̂_n(t) = mean off-diagonal kernel affinity (TWCD2025 eq. 14): (1/n(n−1)) Σ_{i≠j} k_t(X_i,X_j). The
    /// dimension-implicit statistic — it scales as t^{d/2} in the adaptive range, so the empirical Bayes prior reads
    /// d off it without ever estimating d. Cheaper than a full fit (no factorization).
    /// </summary>
    public double KernelAffinity(double t)
    {
        double sum = 0.0;
        for (int i = 0; i < _n; i++)
            for (int j = i + 1; j < _n; j++)
                sum += _kernel.Evaluate(_sqDist[i, j], t);
        return 2.0 * sum / (_n * (double)(_n - 1));   // off-diagonal is symmetric ⇒ count each pair twice
    }

    /// <summary>
    /// Factor A = Σ_t + σ²I and reduce to the marginal evidence log p(Y|t) = −½ YᵀA⁻¹Y − ½ log|A| − (n/2)log2π,
    /// together with α = A⁻¹Y for the conditional mean. The scratch matrix + Cholesky are local, so the call is
    /// pure: an ensemble can drive many bandwidth chains concurrently over one shared <see cref="GpRegression"/>
    /// (only the immutable <see cref="_sqDist"/> precompute is reused across fits).
    /// </summary>
    public GpFit Fit(double t)
    {
        var a = new double[_n, _n];   // working A = Σ_t + σ²I (thread-local — Fit must stay reentrant)
        for (int i = 0; i < _n; i++)
        {
            a[i, i] = _kernel.Evaluate(0.0, t) + _sigma2;
            for (int j = 0; j < i; j++)
            {
                double v = _kernel.Evaluate(_sqDist[i, j], t);
                a[i, j] = v;
                a[j, i] = v;
            }
        }

        var chol = new CholeskyDecomposition(_n);
        chol.Decompose(a);

        // α = A⁻¹Y by triangular solve (O(n²)); the quadratic form YᵀA⁻¹Y is then just Y·α.
        double[] alpha = chol.Solve(_y);
        double quad = 0.0;
        for (int i = 0; i < _n; i++) quad += _y[i] * alpha[i];

        double logMarginal = -0.5 * quad - 0.5 * chol.LogDet - 0.5 * _n * Math.Log(2.0 * Math.PI);
        return new GpFit(t, logMarginal, alpha);
    }

    /// <summary>
    /// Posterior conditional mean f̂(x*) = Σ_j k_t(x*,X_j) α_j at each row of <paramref name="xTest"/>, using the
    /// weights of a <see cref="Fit"/>. Pass the training inputs as <paramref name="xTest"/> for the in-sample fit.
    /// </summary>
    public double[] PredictMean(GpFit fit, double[,] xTest)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(xTest);
        if (xTest.GetLength(1) != _d) throw new ArgumentException("Test inputs must have the training dimension.", nameof(xTest));
        int m = xTest.GetLength(0);
        double t = fit.Bandwidth;
        double[] alpha = fit.Alpha;

        var mean = new double[m];
        for (int a = 0; a < m; a++)
        {
            double f = 0.0;
            for (int j = 0; j < _n; j++)
            {
                double s = 0.0;
                for (int k = 0; k < _d; k++) { double diff = xTest[a, k] - _x[j, k]; s += diff * diff; }
                f += _kernel.Evaluate(s, t) * alpha[j];
            }
            mean[a] = f;
        }
        return mean;
    }
}
