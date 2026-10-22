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
/// He, Yang and Kang's EBIC count law (ABMRS: <i>Biometrics</i> 82(3), 2026, ujag121; preprint 2024). Over C discrete
/// candidates their joint prior is <c>τ(M_k)^(−γ)</c>, with <c>τ(M_k) = C(C,k)</c> the number of k-knot models, so
/// the knot count's law is <c>log p(k) = (1 − γ)·log C(C,k)</c> on [0, C]. γ = 0 is BIC's uniform prior over models,
/// which favors k near C/2; γ = 1 is flat in k; γ &gt; 1 penalizes. The authors recommend γ ∈ [0.5, 2] (default 1);
/// any γ ≥ 0 is accepted. BARS keeps its continuous knot measure, so only the count law is theirs, and C is a
/// complexity parameter rather than a grid. For γ &gt; 1 the law falls until C/2 and rises after it, so choose C at
/// least twice the largest knot count the data can support.
/// </summary>
public sealed class EbicPrior : IComplexityPrior
{
    private readonly double _gamma;
    private readonly int _candidates;

    /// <param name="gamma">He, Yang and Kang's γ: 1 is flat, above 1 penalizes, below 1 favors more knots.</param>
    /// <param name="candidateCount">C, the candidate count; the support is 0 ≤ k ≤ C.</param>
    public EbicPrior(double gamma, int candidateCount)
    {
        if (!(gamma >= 0.0 && double.IsFinite(gamma)))
            throw new ArgumentOutOfRangeException(nameof(gamma), "γ must be finite and non-negative.");
        if (candidateCount < 1) throw new ArgumentOutOfRangeException(nameof(candidateCount));
        _gamma = gamma;
        _candidates = candidateCount;
    }

    public double LogPrior(int knotCount)
    {
        if (knotCount < 0 || knotCount > _candidates) return double.NegativeInfinity;
        return (1.0 - _gamma) * LogBinomial(_candidates, knotCount);
    }

    private static double LogBinomial(int n, int k)
        => SpecialFunctions.LogGamma(n + 1) - SpecialFunctions.LogGamma(k + 1) - SpecialFunctions.LogGamma(n - k + 1);
}
