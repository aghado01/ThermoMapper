using System;
using Maths.Distributions;
using Maths.Rng;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The weight block of scale-mixture robust BARS: draw the per-point weights <c>w</c> that the Normal model's marginal
/// conditions on, from their full conditional given the residuals <c>r = y − Zβ</c> of a coefficient draw and the
/// noise variance σ² drawn with it. One sweep is: knots given w (β and σ² integrated out), then (σ², β) given knots
/// and w, then w given (β, σ²), a valid collapse-then-redraw scan. A null resampler on an ensemble is the
/// non-robust case.
/// </summary>
public interface IWeightResampler
{
    /// <summary>Draw new per-point weights from their full conditional given <paramref name="residuals"/> and
    /// <paramref name="sigma2"/>.</summary>
    double[] Resample(ReadOnlySpan<double> residuals, double sigma2, Xoshiro256PlusPlus rng);
}

/// <summary>
/// Student-t error augmentation: <c>εᵢ ~ N(0, σ²/wᵢ)</c> with <c>wᵢ ~ Gamma(ν/2, rate ν/2)</c>, so
/// <c>wᵢ | rᵢ, σ² ~ Gamma((ν+1)/2, scale 2/(ν + rᵢ²/σ²))</c>, the standard t-as-scale-mixture-of-normals Gibbs step.
/// A true likelihood, so BIC/EBIC-over-k stays valid. <b>Bounded influence</b> — outliers are downweighted
/// (<c>E[wᵢ] → 0</c> as the residual grows), not rejected.
/// </summary>
public sealed class StudentTWeights : IWeightResampler
{
    private readonly double _nu;

    public StudentTWeights(double degreesOfFreedom = 4.0)
    {
        if (!(degreesOfFreedom > 0.0))
            throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "ν must be positive.");
        _nu = degreesOfFreedom;
    }

    public double[] Resample(ReadOnlySpan<double> residuals, double sigma2, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!(sigma2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(sigma2), "σ² must be positive.");

        double shape = (_nu + 1.0) / 2.0;
        var w = new double[residuals.Length];
        for (int i = 0; i < w.Length; i++)
        {
            double z2 = residuals[i] * residuals[i] / sigma2;
            w[i] = Gamma.Sample(rng, shape, 2.0 / (_nu + z2));
        }
        return w;
    }

    private static double[] Residuals(KnotConfig config, double[]? currentWeights, double[] x, double[] y,
                                      IBasis basis, IObservationModel model, out int m)
    {
        double[,] z = basis.Design(config, x);
        double[] coef = model.PosteriorMeanCoefficients(z, y, currentWeights);
        m = y.Length;
        var r = new double[m];
        for (int i = 0; i < m; i++)
        {
            double f = 0.0;
            for (int j = 0; j < coef.Length; j++) f += z[i, j] * coef[j];
            r[i] = y[i] - f;
        }
        return r;
    }

    internal static double[] Residuals(KnotConfig config, double[]? currentWeights, double[] x, double[] y,
                                       IBasis basis, IObservationModel model)
        => Residuals(config, currentWeights, x, y, basis, model, out _);
}

/// <summary>
/// Contaminated-normal augmentation: <c>εᵢ ~ (1−π) N(0,σ²) + π N(0, κ²σ²)</c> — a two-component scale mixture
/// giving near-<b>rejection</b> (not merely downweighting) of gross outliers. Each point's weight is 1 (good
/// component) or 1/κ² ≈ 0 (outlier), drawn from the posterior responsibility — which is itself the per-point
/// outlier probability, a free diagnostic. Still a true likelihood, so BIC/EBIC-over-k stays valid.
/// </summary>
public sealed class ContaminatedNormalWeights : IWeightResampler
{
    private readonly double _pi;
    private readonly double _kappa;

    public ContaminatedNormalWeights(double outlierFraction = 0.05, double inflation = 10.0)
    {
        if (!(outlierFraction > 0.0 && outlierFraction < 1.0))
            throw new ArgumentOutOfRangeException(nameof(outlierFraction), "Outlier fraction must be in (0,1).");
        if (!(inflation > 1.0))
            throw new ArgumentOutOfRangeException(nameof(inflation), "Inflation κ must exceed 1.");
        _pi = outlierFraction;
        _kappa = inflation;
    }

    public double[] Resample(ReadOnlySpan<double> residuals, double sigma2, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!(sigma2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(sigma2), "σ² must be positive.");

        double sigma = Math.Sqrt(sigma2);
        double invKappaSq = 1.0 / (_kappa * _kappa);
        var w = new double[residuals.Length];
        for (int i = 0; i < w.Length; i++)
        {
            double zg = residuals[i] / sigma;
            double zo = residuals[i] / (_kappa * sigma);
            double good = (1.0 - _pi) * Math.Exp(-0.5 * zg * zg) / sigma;
            double outlier = _pi * Math.Exp(-0.5 * zo * zo) / (_kappa * sigma);
            double pOutlier = outlier / (good + outlier);
            w[i] = rng.NextDouble() < pOutlier ? invKappaSq : 1.0;
        }
        return w;
    }
}

/// <summary>Robust dispersion of residuals — the MAD-based σ̂ = 1.4826·median|r − median r|, floored.</summary>
internal static class RobustScale
{
    public static double MadSigma(double[] r)
    {
        double sigma = 1.4826 * Mad(r);
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
