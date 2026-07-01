using System;
using Maths.Distributions;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Prior on the number of interior knots — the dimension penalty the (k,ξ) chain pays for complexity, the
/// declared component ranging from flat (<see cref="UniformPrior"/>) through <see cref="PoissonPrior"/> to the
/// model-space-cardinality <see cref="EbicPrior"/>. Returned up to an additive constant (cancels in ratios).
/// </summary>
public interface IComplexityPrior
{
    /// <summary>Log prior mass on a configuration with <paramref name="knotCount"/> interior knots.</summary>
    double LogPrior(int knotCount);
}

/// <summary>Flat prior on 0..<c>maxKnots</c> interior knots; −∞ outside (a hard cap on dimension).</summary>
public sealed class UniformPrior : IComplexityPrior
{
    private readonly int _maxKnots;

    public UniformPrior(int maxKnots)
    {
        if (maxKnots < 0) throw new ArgumentOutOfRangeException(nameof(maxKnots));
        _maxKnots = maxKnots;
    }

    public double LogPrior(int knotCount)
        => knotCount >= 0 && knotCount <= _maxKnots ? 0.0 : double.NegativeInfinity;
}

/// <summary>Poisson(mean) prior on the knot count — DMGK's default smoothing knob.</summary>
public sealed class PoissonPrior : IComplexityPrior
{
    private readonly double _mean;
    private readonly double _logMean;

    public PoissonPrior(double mean)
    {
        if (!(mean > 0.0)) throw new ArgumentOutOfRangeException(nameof(mean), "Mean must be positive.");
        _mean = mean;
        _logMean = Math.Log(mean);
    }

    public double LogPrior(int knotCount)
    {
        if (knotCount < 0) return double.NegativeInfinity;
        // log p(k) = k·logλ − λ − log(k!)
        double logFactorial = 0.0;
        for (int i = 2; i <= knotCount; i++)
            logFactorial += Math.Log(i);
        return knotCount * _logMean - _mean - logFactorial;
    }
}

/// <summary>
/// EBIC-style complexity prior (He, Yang &amp; Kang 2024) adapted to continuous knots: the joint (k,ξ) prior is
/// <c>τ(M_k)^(−γ)</c> with <c>τ(M_k)=C(n,k)</c> the size of the k-knot model space over <c>n</c> effective
/// candidates, so the k-marginal penalty is <c>log p(k) = −γ·log C(n,k)</c> (the uniform-over-positions factor
/// is supplied by the continuous proposal). γ = 0 is flat; γ → 1 strongly penalizes large model spaces — the
/// guard against over-knotting a dense candidate grid.
/// </summary>
public sealed class EbicPrior : IComplexityPrior
{
    private readonly double _gamma;
    private readonly int _candidates;

    public EbicPrior(double gamma, int candidateCount)
    {
        if (gamma < 0.0 || gamma > 1.0) throw new ArgumentOutOfRangeException(nameof(gamma), "γ must be in [0,1].");
        if (candidateCount < 1) throw new ArgumentOutOfRangeException(nameof(candidateCount));
        _gamma = gamma;
        _candidates = candidateCount;
    }

    public double LogPrior(int knotCount)
    {
        if (knotCount < 0 || knotCount > _candidates) return double.NegativeInfinity;
        return -_gamma * LogBinomial(_candidates, knotCount);
    }

    private static double LogBinomial(int n, int k)
        => SpecialFunctions.LogGamma(n + 1) - SpecialFunctions.LogGamma(k + 1) - SpecialFunctions.LogGamma(n - k + 1);
}
