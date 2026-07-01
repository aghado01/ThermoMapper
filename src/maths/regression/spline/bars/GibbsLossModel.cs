using System;
using Maths.LinAlg;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Gibbs bounded-loss observation model — the robustness-INSIDE-the-waist route (BM2021). The marginal is a
/// pseudo-likelihood from Tukey's bisquare loss, <c>log p ≈ −Σ ρ_c(r̂/σ̂) − (ν/2) log m</c>, with β̂ the Tukey
/// M-estimate (IRLS, LS warm start) and σ̂ a <i>config-free</i> robust noise scale (1.4826·MAD of the scaled
/// first differences) so the comparison across knot configs is consistent. ρ_c is normalized so
/// <c>ρ_c(u) ≈ u²/2</c> for small <c>u</c> (Gaussian-comparable to the Normal deviance). Near-redescending
/// rejection of gross outliers, but a <b>generalized / Gibbs posterior</b>: its k-selection consistency and
/// credible-interval interpretation are weaker than the scale-mixture route (which keeps a true likelihood).
/// Prefer the scale-mixture resamplers unless rejection (not downweighting) is the priority.
/// </summary>
public sealed class GibbsLossModel : IObservationModel
{
    private const int MaxIterations = 50;
    private const double Tolerance = 1e-9;
    private readonly double _c;

    public GibbsLossModel(double tuning = 4.685)
    {
        if (!(tuning > 0.0)) throw new ArgumentOutOfRangeException(nameof(tuning), "Tuning constant must be positive.");
        _c = tuning;
    }

    /// <inheritdoc/>
    public int EffectiveDimension(double[,] design, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        return design.GetLength(1);
    }

    /// <inheritdoc/>
    public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        int m = design.GetLength(0);
        int nu = design.GetLength(1);
        if (y.Length != m)
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        if (m <= nu)
            return double.NegativeInfinity;

        double sigma = RobustNoise(y);
        double[] beta = FitM(design, y, sigma);
        double loss = 0.0;
        for (int i = 0; i < m; i++)
            loss += Rho((y[i] - Fit(design, beta, i)) / sigma);
        return -loss - 0.5 * nu * Math.Log(m);
    }

    /// <inheritdoc/>
    public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        return FitM(design, y, RobustNoise(y));
    }

    // Tukey bisquare M-estimation by banded IRLS at fixed scale σ, warm-started with one LS step (iter 0, w = 1).
    private double[] FitM(double[,] design, double[] y, double sigma)
    {
        int nu = design.GetLength(1);
        var bd = new BandedDesign(design);
        int m = bd.Rows;
        var beta = new double[nu];
        var w = new double[m];
        var band = new double[bd.Bandwidth + 1, nu];
        var b = new double[nu];
        var chol = new BandCholesky(nu, bd.Bandwidth, BandFactorization.Ldlt);

        for (int iter = 0; iter < MaxIterations; iter++)
        {
            for (int i = 0; i < m; i++)
                w[i] = iter == 0 ? 1.0 : TukeyWeight((y[i] - Fit(design, beta, i)) / sigma);
            Array.Clear(band);
            Array.Clear(b);
            bd.Accumulate(w, band, y, b);

            chol.DecomposeBanded(band);
            double[] sol = chol.Solve(b);

            double shift = 0.0;
            for (int p = 0; p < nu; p++)
            {
                double d = sol[p] - beta[p];
                shift += d * d;
                beta[p] = sol[p];
            }
            if (iter > 0 && shift < Tolerance) break;
        }
        return beta;
    }

    // Tukey bisquare, normalized so ρ_c(u) ≈ u²/2 for small u.
    private double Rho(double u)
    {
        double a = u / _c;
        if (Math.Abs(a) >= 1.0) return _c * _c / 6.0;
        double t = 1.0 - a * a;
        return (_c * _c / 6.0) * (1.0 - t * t * t);
    }

    private double TukeyWeight(double u)
    {
        double a = u / _c;
        if (Math.Abs(a) >= 1.0) return 0.0;
        double t = 1.0 - a * a;
        return t * t;
    }

    private static double Fit(double[,] design, double[] beta, int row)
    {
        double f = 0.0;
        for (int j = 0; j < beta.Length; j++) f += design[row, j] * beta[j];
        return f;
    }

    // Config-free robust noise scale: 1.4826 · MAD of the √2-scaled first differences of y.
    private static double RobustNoise(double[] y)
    {
        int m = y.Length;
        if (m < 2) return 1.0;
        var d = new double[m - 1];
        for (int i = 0; i < m - 1; i++) d[i] = (y[i + 1] - y[i]) / Math.Sqrt(2.0);
        double sigma = 1.4826 * Mad(d);
        return sigma > 0.0 ? sigma : 1e-6;
    }

    private static double Mad(double[] r)
    {
        double med = Median(r);
        var dev = new double[r.Length];
        for (int i = 0; i < r.Length; i++) dev[i] = Math.Abs(r[i] - med);
        return Median(dev);
    }

    private static double Median(double[] values)
    {
        var sorted = (double[])values.Clone();
        Array.Sort(sorted);
        int n = sorted.Length;
        return n % 2 == 1 ? sorted[n / 2] : 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
    }
}
