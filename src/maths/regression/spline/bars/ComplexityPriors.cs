using System;
using Maths.Distributions;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Prior on the number of interior knots — the dimension penalty the (k,ξ) chain pays for complexity, the
/// declared component ranging from flat (<see cref="UniformPrior"/>) through <see cref="PoissonPrior"/> to the
/// model-space-cardinality <see cref="EbicPrior"/>. Returned up to an additive constant (cancels in ratios). The
/// target supplies the ordered-uniform location density, so <c>p(k)</c> is the knot-count marginal.
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
/// Conway–Maxwell–Poisson prior on the knot count (Conway and Maxwell 1962; Shmueli et al. 2005):
/// <c>p(k) ∝ λᵏ/(k!)^ν</c>. ν sets how fast the penalty grows with k. ν = 1 is Poisson(λ). ν = 2 is the law BARS
/// ran under before the target carried the ordered-uniform location density (R01), since that left
/// <c>Poisson(λ)/k!</c>. ν &lt; 1 has heavier tails. Unnormalized: the sampler uses only the ratio
/// <c>p(k+1)/p(k) = λ/(k+1)^ν</c>.
/// </summary>
public sealed class ConwayMaxwellPoissonPrior : IComplexityPrior
{
    private readonly double _logLambda;
    private readonly double _nu;

    public ConwayMaxwellPoissonPrior(double lambda, double nu)
    {
        if (!(lambda > 0.0)) throw new ArgumentOutOfRangeException(nameof(lambda), "λ must be positive.");
        if (!(nu >= 0.0)) throw new ArgumentOutOfRangeException(nameof(nu), "ν must be non-negative.");
        if (nu == 0.0 && !(lambda < 1.0))
            throw new ArgumentOutOfRangeException(nameof(lambda), "With ν = 0 the prior is geometric and needs λ < 1.");
        _logLambda = Math.Log(lambda);
        _nu = nu;
    }

    public double LogPrior(int knotCount)
        => knotCount < 0 ? double.NegativeInfinity : knotCount * _logLambda - _nu * SpecialFunctions.LogGamma(knotCount + 1.0);
}

/// <summary>
/// EBIC-style complexity prior (He, Yang &amp; Kang 2024) adapted to continuous knots: the joint (k,ξ) prior is
/// <c>τ(M_k)^(−γ)</c> with <c>τ(M_k)=C(n,k)</c> the size of the k-knot model space over <c>n</c> effective
/// candidates, so the k-marginal is <c>log p(k) = −γ·log C(n,k)</c>. γ = 0 is flat; γ → 1 strongly penalizes large
/// model spaces, the guard against over-knotting a dense candidate grid. The law is symmetric in k ↔ n − k, so it
/// penalizes only below n/2. It is not He, Yang &amp; Kang's count law: their prior puts C(n,k)^(−γ) on each of the
/// C(n,k) models, so theirs is C(n,k)^(1−γ).
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
