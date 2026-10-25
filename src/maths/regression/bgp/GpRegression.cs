using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Bgp;

/// <summary>How a <see cref="GpRegression"/> treats the noise and signal variances.</summary>
public enum GpNoiseModel
{
    /// <summary>
    /// σ² is known and given, and the kernel's signal variance is 1: K = Σ_t + σ²I, as in TWCD2025's algorithm. The
    /// response must then be on a scale where a unit signal variance is sensible.
    /// </summary>
    Known,

    /// <summary>
    /// The signal variance s² and the noise variance σ² = s²η are both unknown: K = s²(Σ_t + ηI), where η is the
    /// noise-to-signal ratio. s² is integrated out under the 1/s² prior, leaving a marginal in (t, η), and drawn exactly
    /// given them, s² | t, η, y ~ IG((n − 1)/2, Q/2). Rescaling the response rescales s² and σ² only.
    /// </summary>
    UnknownScale,
}

/// <summary>A GP fit at one bandwidth (and, for <see cref="GpNoiseModel.UnknownScale"/>, one noise ratio): the marginal
/// evidence with the constant mean integrated out, the GLS mean m̂ and the posterior weights α = R⁻¹(Y − m̂1) that
/// produce the conditional-mean prediction. Returned by <see cref="GpRegression.Fit(double)"/> and
/// <see cref="GpRegression.Fit(double, double)"/>.</summary>
/// <param name="Bandwidth">The kernel bandwidth t this fit was computed at, in the GP's mapped frame.</param>
/// <param name="LogMarginal">The log evidence with the constant mean integrated out under its flat prior: log p(Y | t)
/// for a known σ², or log p(Y | t, η) with s² also integrated out under 1/s² for an unknown scale.</param>
/// <param name="Alpha">Posterior weights α = R⁻¹(Y − m̂1) with R = Σ_t + ηI (η = σ² for a known σ²).</param>
/// <param name="Mean">The GLS mean m̂ = 1ᵀR⁻¹Y / 1ᵀR⁻¹1, the posterior mean of m; the conditional mean is m̂ + k_t(·,X)·α,
/// which does not depend on s².</param>
/// <param name="NoiseRatio">η, the nugget added to the unit-amplitude kernel: σ² itself for a known σ² (s² = 1).</param>
/// <param name="Quadratic">Q = (Y − m̂1)ᵀR⁻¹(Y − m̂1), the residual quadratic that s²'s conditional posterior reads.</param>
public sealed record GpFit(double Bandwidth, double LogMarginal, double[] Alpha, double Mean, double NoiseRatio, double Quadratic);

/// <summary>A <see cref="GpFit"/> with what its conditional variance needs: the factor of R = Σ_t + ηI, R⁻¹1 and
/// 1ᵀR⁻¹1. It costs n² numbers, so a sampler holds one per chain while it is current and never retains it.</summary>
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
    private readonly double _sigma2;      // known noise variance σ² (NaN for an unknown scale)
    private readonly IGpKernel _kernel;
    private readonly int _n;
    private readonly int _d;              // ambient dimension D
    private readonly double[,] _sqDist;   // n × n precomputed ‖X_i − X_j‖² in the mapped frame (independent of t)

    /// <summary>A GP with a known noise variance and unit signal variance (<see cref="GpNoiseModel.Known"/>).</summary>
    /// <param name="x">n × D matrix of ambient predictor coordinates, in any units.</param>
    /// <param name="y">Length-n responses.</param>
    /// <param name="sigma2">Known observation-noise variance σ² &gt; 0, in response units.</param>
    /// <param name="kernel">Radial covariance kernel (squared-exponential is the canonical choice).</param>
    /// <param name="scaling">How the inputs are mapped before fitting (default: isotropic similarity into the unit cube).</param>
    public GpRegression(double[,] x, double[] y, double sigma2, IGpKernel kernel,
                        GpInputScaling scaling = GpInputScaling.Similarity)
        : this(x, y, kernel, scaling, sigma2 > 0.0
            ? sigma2 : throw new ArgumentOutOfRangeException(nameof(sigma2), "σ² must be positive."))
    {
    }

    /// <summary>A GP whose signal and noise variances are both unknown (<see cref="GpNoiseModel.UnknownScale"/>): fit it
    /// with <see cref="Fit(double, double)"/>. A constant response has no scale and is refused.</summary>
    public static GpRegression WithUnknownScale(double[,] x, double[] y, IGpKernel kernel,
                                                GpInputScaling scaling = GpInputScaling.Similarity)
    {
        var gp = new GpRegression(x, y, kernel, scaling, double.NaN);
        bool constant = true;
        foreach (double v in gp._y) constant &= v == 0.0;
        if (constant) throw new ArgumentException("A constant response has no scale to infer.", nameof(y));
        return gp;
    }

    private GpRegression(double[,] x, double[] y, IGpKernel kernel, GpInputScaling scaling, double sigma2)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(kernel);
        NoiseModel = double.IsNaN(sigma2) ? GpNoiseModel.UnknownScale : GpNoiseModel.Known;
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

    /// <summary>How this GP treats the noise and signal variances.</summary>
    public GpNoiseModel NoiseModel { get; }

    /// <summary>The fitted input map and response offset shared by training and prediction.</summary>
    public GpTransform Transform { get; }

    /// <summary>Pairwise squared distances in the mapped frame — the frame the bandwidth prior must read.</summary>
    internal double[,] SquaredDistances => _sqDist;

    /// <summary>Add what a fit depends on to a run fingerprint: the mapped inputs and their scale, the shifted responses
    /// and the shift, the noise variance and the kernel's type.</summary>
    internal void AddTo(RunFingerprint fingerprint)
    {
        fingerprint.Add(_n).Add(_d);
        foreach (double v in _x) fingerprint.Add(v);
        fingerprint.Add(_y).Add(_yShift).Add((int)NoiseModel).Add(_sigma2).Add(Transform.InputScale).Add(_kernel.GetType().FullName!);
    }

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
    /// For a known σ² only.
    /// </summary>
    public GpFit Fit(double t) => FitFactored(t).Fit;

    /// <summary>
    /// For an unknown scale (<see cref="GpNoiseModel.UnknownScale"/>): factor R = Σ_t + ηI and reduce to the evidence
    /// with the constant mean integrated out under its flat prior and s² under 1/s²:
    /// log p(Y | t, η) = −½ log|R| − ½ log 1ᵀR⁻¹1 − ((n − 1)/2) log Q + log Γ((n − 1)/2) − ((n − 1)/2) log π, with
    /// Q = (Y − m̂1)ᵀR⁻¹(Y − m̂1). The GLS mean and α are as for a known σ², with η in its place; neither depends on s².
    /// </summary>
    public GpFit Fit(double t, double eta) => FitFactored(t, eta).Fit;

    /// <summary>As <see cref="Fit(double)"/>, keeping the factor, K⁻¹1 and 1ᵀK⁻¹1 that the conditional variance needs.</summary>
    internal GpFactoredFit FitFactored(double t)
    {
        if (NoiseModel != GpNoiseModel.Known)
            throw new InvalidOperationException("This GP has an unknown scale; fit it with Fit(t, η).");
        return FitCore(t, _sigma2);
    }

    /// <summary>As <see cref="Fit(double, double)"/>, keeping the factor, R⁻¹1 and 1ᵀR⁻¹1.</summary>
    internal GpFactoredFit FitFactored(double t, double eta)
    {
        if (NoiseModel != GpNoiseModel.UnknownScale)
            throw new InvalidOperationException("This GP has a known σ²; fit it with Fit(t).");
        if (!(eta > 0.0 && double.IsFinite(eta))) throw new ArgumentOutOfRangeException(nameof(eta), "η must be positive and finite.");
        return FitCore(t, eta);
    }

    // A fit at its own bandwidth and nugget, whichever noise model made it: what refitting needs.
    private GpFactoredFit Refit(GpFit fit) => FitCore(fit.Bandwidth, fit.NoiseRatio);

    // The fit of R = Σ_t + nugget·I: σ² for a known σ², η for an unknown scale. Only the log evidence differs.
    private GpFactoredFit FitCore(double t, double nugget)
    {
        var a = new double[_n, _n];   // working A = Σ_t + nugget·I (thread-local — Fit must stay reentrant)
        for (int i = 0; i < _n; i++)
        {
            a[i, i] = _kernel.Evaluate(0.0, t) + nugget;
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
            throw new InvalidOperationException($"Σ_t + {nugget:G4}·I is not numerically positive definite at t = {t} ({chol.Status}).");

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

        double logMarginal = NoiseModel == GpNoiseModel.Known
            ? -0.5 * quad - 0.5 * chol.LogDet - 0.5 * Math.Log(oneKOne) - 0.5 * (_n - 1) * Math.Log(2.0 * Math.PI)
            : -0.5 * chol.LogDet - 0.5 * Math.Log(oneKOne) - 0.5 * (_n - 1) * Math.Log(quad)
              + SpecialFunctions.LogGamma(0.5 * (_n - 1)) - 0.5 * (_n - 1) * Math.Log(Math.PI);
        return new GpFactoredFit(new GpFit(t, logMarginal, alpha, _yShift + mHat, nugget, quad), chol, kInvOne, oneKOne);
    }

    /// <summary>
    /// Each unit's leave-unit-out predictive log p(y_U | y_rest, t, η) at one fit (OQ17's collapsed content): f, the
    /// constant mean and, for an unknown scale, s² are integrated out against their posterior given the observations
    /// outside the unit. The flat mean leaves the projected precision P = R⁻¹ − R⁻¹11ᵀR⁻¹/1ᵀR⁻¹1, whose PY is the fit's α
    /// and YᵀPY its Q. A known σ² gives a Gaussian. An unknown scale under 1/s² gives a Student t over n − 1 effective
    /// observations. Either way each term is the evidence difference log p(Y | t, η) − log p(Y_rest | t, η). One refit and
    /// one inverse, O(n³), serve every unit. A unit must leave one observation outside it, two for an unknown scale.
    /// </summary>
    public double[] LeaveOutLogDensities(GpFit fit, UnitPartition units)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(units);
        if (units.Observations != _n) throw new ArgumentException("The partition must cover the fitted observations.", nameof(units));
        if (fit.Alpha.Length != _n) throw new ArgumentException("The fit is not of this GP's data.", nameof(fit));
        GpFactoredFit factored = Refit(fit);
        double[] alpha = factored.Fit.Alpha, kInvOne = factored.KInvOne;
        double oneKOne = factored.OneKOne;
        var rInv = new double[_n, _n];
        factored.Factor.WriteInverseTo(rInv);
        var terms = new double[units.Units];
        for (int u = 0; u < units.Units; u++)
        {
            int[] members = units.Members[u];
            int m = members.Length;
            if (m >= _n) throw new ArgumentException("A unit must leave some observations outside it.", nameof(units));
            var block = new double[m, m];   // P_UU
            var pr = new double[m];         // (PY)_U = α_U
            for (int a = 0; a < m; a++)
            {
                int i = members[a];
                pr[a] = alpha[i];
                for (int b = 0; b < m; b++)
                {
                    int j = members[b];
                    block[a, b] = rInv[i, j] - kInvOne[i] * kInvOne[j] / oneKOne;
                }
            }
            terms[u] = NoiseModel == GpNoiseModel.Known
                ? LeaveOutPredictive.Gaussian(block, pr)
                : LeaveOutPredictive.StudentT(block, pr, factored.Fit.Quadratic, _n - 1);
        }
        return terms;
    }

    /// <summary>
    /// The conditional variance of the regression function f(x*) given the bandwidth, at each row of
    /// <paramref name="xTest"/> (original units), with the constant mean integrated out under its flat prior:
    /// Var = k(x*,x*) − k*ᵀK⁻¹k* + (1 − 1ᵀK⁻¹k*)² / 1ᵀK⁻¹1 (universal kriging; Rasmussen &amp; Williams 2006, §2.7).
    /// The last term is the mean's own uncertainty. It is the variance of f, not of a new observation, which adds σ².
    /// For an unknown scale it is in units of the signal variance: multiply by s². This refactors K at the fit's
    /// bandwidth and nugget (O(n³)); <see cref="BgpSampler"/> computes it during the run instead, while each fit's
    /// factor still exists.
    /// </summary>
    public double[] PredictVariance(GpFit fit, double[,] xTest)
    {
        ArgumentNullException.ThrowIfNull(fit);
        ArgumentNullException.ThrowIfNull(xTest);
        if (xTest.GetLength(1) != _d) throw new ArgumentException("Test inputs must have the training dimension.", nameof(xTest));
        return ConditionalVarianceInFrame(Refit(fit), Transform.MapInputs(xTest));
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
