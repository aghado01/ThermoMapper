using System;
using Maths.Distributions;
using Maths.Rng;

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
/// <c>wᵢ | rᵢ, σ² ~ Gamma((ν+1)/2, rate (ν + rᵢ²/σ²)/2)</c>, the standard t-as-scale-mixture-of-normals Gibbs step.
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
            w[i] = Gamma.ShapeRate(shape, 0.5 * (_nu + z2)).Sample(rng);
        }
        return w;
    }
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

        // The responsibility from its log odds, log(π/((1−π)κ)) + ½·(r²/σ²)·(1 − 1/κ²): the two densities underflow
        // together past about 386σ at κ = 10, where their ratio would be 0/0.
        double invKappaSq = 1.0 / (_kappa * _kappa);
        double logPriorOdds = Math.Log(_pi / ((1.0 - _pi) * _kappa));
        double gain = 0.5 * (1.0 - invKappaSq);
        var w = new double[residuals.Length];
        for (int i = 0; i < w.Length; i++)
        {
            double logOdds = logPriorOdds + gain * residuals[i] * residuals[i] / sigma2;
            double pOutlier = 1.0 / (1.0 + Math.Exp(-logOdds));
            w[i] = rng.NextDouble() < pOutlier ? invKappaSq : 1.0;
        }
        return w;
    }
}
