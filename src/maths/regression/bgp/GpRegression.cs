using System;
using Maths.LinAlg;

namespace Maths.Regression.Bgp;

/// <summary>A GP fit at one bandwidth: the marginal evidence with the constant mean integrated out, the GLS mean m̂ and
/// the posterior weights α = K⁻¹(Y − m̂1) that produce the conditional-mean prediction. Returned by
/// <see cref="GpRegression.Fit"/>.</summary>
/// <param name="Bandwidth">The kernel bandwidth t this fit was computed at, in the GP's mapped frame.</param>
/// <param name="LogMarginal">log ∫ p(Y | m, t) dm under a flat prior on the constant mean m (the t-posterior's likelihood factor).</param>
/// <param name="Alpha">Posterior weights α = K⁻¹(Y − m̂1) with K = Σ_t + σ²I.</param>
/// <param name="Mean">The GLS mean m̂ = 1ᵀK⁻¹Y / 1ᵀK⁻¹1, the posterior mean of m; the conditional mean is m̂ + k_t(·,X)·α.</param>
public sealed record GpFit(double Bandwidth, double LogMarginal, double[] Alpha, double Mean);

/// <summary>A <see cref="GpFit"/> with what its conditional variance needs: the factor of K = Σ_t + σ²I, K⁻¹1 and
/// 1ᵀK⁻¹1. It costs n² numbers, so a sampler holds one per chain while it is current and never retains it.</summary>
internal sealed record GpFactoredFit(GpFit Fit, CholeskyDecomposition Factor, double[] KInvOne, double OneKOne);

/// <summary>
/// Bayesian GP regression with an ambient radial kernel (Tang, Wu, Cheng &amp; Dunson 2025): the conjugate core
/// behind the dimension-adaptive method. With the regression function marginalized out by GP conjugacy, all that
/// remains is the bandwidth t — and this class supplies the two quantities a t-sampler needs: the marginal
/// evidence log p(Y|t) (reduce data+model → the t-likelihood) and the posterior weights α for the conditional
/// mean. It also exposes the kernel-affinity observable v̂_n(t) (the dimension-implicit statistic the empirical
/// Bayes prior is built from). The linear algebra is dense (an n×n SPD Cholesky), unlike the banded free-knot
/// world; pairwise squared distances are t-independent, so they are computed once and reused across every fit.
/// One fitted <see cref="GpTransform"/> is shared by training and prediction: inputs are mapped into its frame (by
/// default an isotropic similarity into the unit cube), and bandwidths are in that frame. The GP's constant mean is
/// integrated out under a flat prior, so the evidence carries the t-dependent term −½ log 1ᵀK⁻¹1 and predictions
/// revert to the GLS mean away from the data.
/// </summary>
public sealed class GpRegression
{
    private readonly double[,] _x;        // n × D train predictors in the mapped frame
    private readonly double[] _y;         // n responses less _yShift
    private readonly double _yShift;      // training mean, subtracted for conditioning only (the mean is integrated out)
    private readonly double _sigma2;      // known noise variance σ²
    private readonly IGpKernel _kernel;
    private readonly int _n;
    private readonly int _d;              // ambient dimension D
    private readonly double[,] _sqDist;   // n × n precomputed ‖X_i − X_j‖² in the mapped frame (independent of t)

    /// <param name="x">n × D matrix of ambient predictor coordinates, in any units.</param>
    /// <param name="y">Length-n responses.</param>
    /// <param name="sigma2">Known observation-noise variance σ² &gt; 0, in response units.</param>
    /// <param name="kernel">Radial covariance kernel (squared-exponential is the canonical choice).</param>
    /// <param name="scaling">How the inputs are mapped before fitting (default: isotropic similarity into the unit cube).</param>
    public GpRegression(double[,] x, double[] y, double sigma2, IGpKernel kernel,
                        GpInputScaling scaling = GpInputScaling.Similarity)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(kernel);
        if (!(sigma2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(sigma2), "σ² must be positive.");
        Transform = GpTransform.Fit(x, scaling);
        _n = x.GetLength(0);
        _d = x.GetLength(1);
        if (y.Length != _n) throw new ArgumentException("y length must match the number of rows of x.", nameof(y));
        double mean = 0.0;
        foreach (double v in y)
        {
            if (!double.IsFinite(v)) throw new ArgumentException("Responses must be finite.", nameof(y));
            mean += v;
        }
        _yShift = mean / _n;
        _x = Transform.MapInputs(x);
        _y = new double[_n];
        for (int i = 0; i < _n; i++) _y[i] = y[i] - _yShift;
        _sigma2 = sigma2; _kernel = kernel;

        _sqDist = new double[_n, _n];
        for (int i = 0; i < _n; i++)
            for (int j = i + 1; j < _n; j++)
            {
                double s = 0.0;
                for (int k = 0; k < _d; k++) { double diff = _x[i, k] - _x[j, k]; s += diff * diff; }
                _sqDist[i, j] = s;
                _sqDist[j, i] = s;
            }
    }

    public int Count => _n;

    /// <summary>The fitted input map and response offset shared by training and prediction.</summary>
    public GpTransform Transform { get; }

    /// <summary>Pairwise squared distances in the mapped frame — the frame the bandwidth prior must read.</summary>
    internal double[,] SquaredDistances => _sqDist;

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
    /// Factor K = Σ_t + σ²I and reduce to the marginal evidence with the constant mean m integrated out under a flat
    /// prior: log ∫ p(Y|m,t) dm = −½[YᵀK⁻¹Y − (1ᵀK⁻¹Y)²/(1ᵀK⁻¹1)] − ½ log|K| − ½ log 1ᵀK⁻¹1 − ((n−1)/2)log2π,
    /// together with the GLS mean m̂ = 1ᵀK⁻¹Y / 1ᵀK⁻¹1 and α = K⁻¹(Y − m̂1); one extra solve (K⁻¹1) on the same
    /// factor. The evidence is invariant to a constant shift of Y. The scratch matrix + Cholesky are local, so the
    /// call is pure: an ensemble can drive many bandwidth chains concurrently over one shared
    /// <see cref="GpRegression"/> (only the immutable <see cref="_sqDist"/> precompute is reused across fits).
    /// </summary>
    public GpFit Fit(double t) => FitFactored(t).Fit;

    /// <summary>As <see cref="Fit"/>, keeping the factor, K⁻¹1 and 1ᵀK⁻¹1 that the conditional variance needs.</summary>
    internal GpFactoredFit FitFactored(double t)
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
        if (chol.Status is not { PositiveDefinite: true, FlooredPivots: 0 })   // the nugget makes it positive definite; clamping means bad input
            throw new InvalidOperationException($"Σ_t + σ²I is not numerically positive definite at t = {t} ({chol.Status}).");

        // K⁻¹Y and K⁻¹1 by triangular solves (O(n²) each); the GLS mean and the residual quadratic follow.
        double[] kInvY = chol.Solve(_y);
        var ones = new double[_n];
        Array.Fill(ones, 1.0);
        double[] kInvOne = chol.Solve(ones);
        double oneKOne = 0.0, oneKY = 0.0;
        for (int i = 0; i < _n; i++) { oneKOne += kInvOne[i]; oneKY += kInvY[i]; }
        double mHat = oneKY / oneKOne;   // GLS mean of the shifted response

        var alpha = new double[_n];      // K⁻¹(Y − m̂1); its quadratic (Y − m̂1)ᵀK⁻¹(Y − m̂1) = YᵀK⁻¹Y − (1ᵀK⁻¹Y)²/1ᵀK⁻¹1
        double quad = 0.0;
        for (int i = 0; i < _n; i++)
        {
            alpha[i] = kInvY[i] - mHat * kInvOne[i];
            quad += (_y[i] - mHat) * alpha[i];
        }

        double logMarginal = -0.5 * quad - 0.5 * chol.LogDet - 0.5 * Math.Log(oneKOne)
                             - 0.5 * (_n - 1) * Math.Log(2.0 * Math.PI);
        return new GpFactoredFit(new GpFit(t, logMarginal, alpha, _yShift + mHat), chol, kInvOne, oneKOne);
    }

    /// <summary>
    /// The conditional variance of the regression function f(x*) given the bandwidth, at each row of
    /// <paramref name="xTest"/> (original units), with the constant mean integrated out under its flat prior:
    /// Var = k(x*,x*) − k*ᵀK⁻¹k* + (1 − 1ᵀK⁻¹k*)² / 1ᵀK⁻¹1 (universal kriging; Rasmussen &amp; Williams 2006, §2.7).
    /// The last term is the mean's own uncertainty. It is the variance of f, not of a new observation, which adds σ².
    /// This refactors K at the fit's bandwidth (O(n³)); <see cref="BgpSampler"/> computes it during the run instead,
    /// while each fit's factor still exists.
    /// </summary>
    public double[] PredictVariance(GpFit fit, double[,] xTest)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(xTest);
        if (xTest.GetLength(1) != _d) throw new ArgumentException("Test inputs must have the training dimension.", nameof(xTest));
        return ConditionalVarianceInFrame(FitFactored(fit.Bandwidth), Transform.MapInputs(xTest));
    }

    /// <summary>As <see cref="PredictVariance"/>, from a factored fit and test inputs already in the fitted frame: one
    /// O(n²) solve per test point. Rounding can leave a variance a hair below zero where the data pin f; it is floored
    /// at 0.</summary>
    internal double[] ConditionalVarianceInFrame(GpFactoredFit factored, double[,] mappedTest)
    {
        int m = mappedTest.GetLength(0);
        double t = factored.Fit.Bandwidth;
        double prior = _kernel.Evaluate(0.0, t);
        var k = new double[_n];
        var variance = new double[m];
        for (int a = 0; a < m; a++)
        {
            for (int j = 0; j < _n; j++)
            {
                double s = 0.0;
                for (int c = 0; c < _d; c++) { double diff = mappedTest[a, c] - _x[j, c]; s += diff * diff; }
                k[j] = _kernel.Evaluate(s, t);
            }
            double[] kInvK = factored.Factor.Solve(k);
            double quad = 0.0, oneKInvK = 0.0;
            for (int j = 0; j < _n; j++) { quad += k[j] * kInvK[j]; oneKInvK += kInvK[j]; }
            double u = 1.0 - oneKInvK;
            variance[a] = Math.Max(0.0, prior - quad + u * u / factored.OneKOne);
        }
        return variance;
    }

    /// <summary>
    /// Posterior conditional mean f̂(x*) = m̂ + Σ_j k_t(x*,X_j) α_j at each row of <paramref name="xTest"/> (original
    /// units, mapped by the fitted <see cref="Transform"/>), using the weights of a <see cref="Fit"/>. Pass the
    /// training inputs as <paramref name="xTest"/> for the in-sample fit.
    /// </summary>
    public double[] PredictMean(GpFit fit, double[,] xTest)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(xTest);
        if (xTest.GetLength(1) != _d) throw new ArgumentException("Test inputs must have the training dimension.", nameof(xTest));
        return PredictMeanInFrame(fit, Transform.MapInputs(xTest));
    }

    /// <summary>As <see cref="PredictMean"/>, for test inputs already mapped into the fitted frame.</summary>
    internal double[] PredictMeanInFrame(GpFit fit, double[,] mappedTest)
    {
        int m = mappedTest.GetLength(0);
        double t = fit.Bandwidth;
        double[] alpha = fit.Alpha;

        var mean = new double[m];
        for (int a = 0; a < m; a++)
        {
            double f = fit.Mean;
            for (int j = 0; j < _n; j++)
            {
                double s = 0.0;
                for (int k = 0; k < _d; k++) { double diff = mappedTest[a, k] - _x[j, k]; s += diff * diff; }
                f += _kernel.Evaluate(s, t) * alpha[j];
            }
            mean[a] = f;
        }
        return mean;
    }
}
