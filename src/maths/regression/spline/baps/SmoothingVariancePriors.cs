using System;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// A prior on a smoothing variance τ², in the form the samplers of the log variance <c>ρ = ln τ²</c> use: the log density
/// of ρ, which is the prior's log density at <c>e^ρ</c> plus the Jacobian ρ, up to a constant, with its first and second
/// derivatives in ρ. The anisotropic schemes put it on each component's τⱼ² (<see cref="BapsScheme.AnisotropicGibbs"/>,
/// <see cref="HierarchicalPSpline"/> with anisotropic smoothing).
/// </summary>
public interface ISmoothingVariancePrior
{
    /// <summary><c>ln q(ρ)</c> up to a constant, with <c>d ln q/dρ</c> and <c>d² ln q/dρ²</c>.</summary>
    double LogDensityOfLog(double rho, out double first, out double second);
}

/// <summary>
/// τ² ~ IG(shape a, scale b), the prior the conjugate schemes give τ² and the anisotropic schemes' default:
/// <c>ln q(ρ) = −aρ − b e^(−ρ)</c>. With small a and b it is nearly flat in ln τ² over a wide range, and its upper tail
/// in ρ decays only like <c>e^(−aρ)</c>, so the likelihood alone bounds large smoothing variances.
/// </summary>
public sealed record InverseGammaVariancePrior : ISmoothingVariancePrior
{
    public InverseGammaVariancePrior(double shape = 1e-3, double scale = 1e-3)
    {
        if (!(shape > 0.0) || double.IsInfinity(shape)) throw new ArgumentOutOfRangeException(nameof(shape), "Shape must be positive and finite.");
        if (!(scale > 0.0) || double.IsInfinity(scale)) throw new ArgumentOutOfRangeException(nameof(scale), "Scale must be positive and finite.");
        Shape = shape;
        Scale = scale;
    }

    /// <summary>Shape a.</summary>
    public double Shape { get; }

    /// <summary>Scale b.</summary>
    public double Scale { get; }

    public double LogDensityOfLog(double rho, out double first, out double second)
    {
        double pull = Scale * Math.Exp(-rho);
        first = -Shape + pull;
        second = -pull;
        return -Shape * rho - pull;
    }
}

/// <summary>
/// τ² ~ W(shape k, scale s) (<see cref="Weibull"/>): <c>ln q(ρ) = kρ − e^(k(ρ − ln s))</c>. With k = ½, the default, it is the
/// law of τ² for an exponential τ with rate 1/√s, the scale-dependent prior of Klein &amp; Kneib (2016) that Bach &amp;
/// Klein (2025) prefer for smoothing variances. Its upper tail in ρ decays doubly exponentially, so it holds the variances
/// within range where the inverse gamma lets them wander, at the price of a scale that means something:
/// <see cref="ScaledTo"/> sets s from a bound on the function's size, and <see cref="ScaledToResponse"/> takes the bound
/// from the data.
/// <para>
/// Parameterizations: Bach &amp; Klein write the prior as Weibull(½, λ) with λ a rate, density ∝ x^(−½) e^(−√(λx)), so their
/// λ is <see cref="Rate"/> = 1/s (their kernel −√λ e^(ρ/2) + ρ/2 is this one's); Klein &amp; Kneib's θ is the rate of the
/// exponential τ, θ² = 1/s.
/// </para>
/// </summary>
public sealed record WeibullVariancePrior : ISmoothingVariancePrior
{
    private readonly double _logScale;

    public WeibullVariancePrior(double scale, double shape = 0.5)
    {
        if (!(scale > 0.0) || double.IsInfinity(scale)) throw new ArgumentOutOfRangeException(nameof(scale), "Scale must be positive and finite.");
        if (!(shape > 0.0) || double.IsInfinity(shape)) throw new ArgumentOutOfRangeException(nameof(shape), "Shape must be positive and finite.");
        Scale = scale;
        Shape = shape;
        _logScale = Math.Log(scale);
    }

    /// <summary>Scale s.</summary>
    public double Scale { get; }

    /// <summary>Shape k.</summary>
    public double Shape { get; }

    /// <summary>The rate 1/s: Bach &amp; Klein's λ.</summary>
    public double Rate => 1.0 / Scale;

    public double LogDensityOfLog(double rho, out double first, out double second)
    {
        double t = Math.Exp(Shape * (rho - _logScale));   // (τ²/s)^k
        first = Shape - Shape * t;
        second = -Shape * Shape * t;
        return Shape * rho - t;
    }

    /// <summary>
    /// Prior scaling, Bach &amp; Klein's (2025, Supplement Appendix D; their R function <c>prior_scaling_WB</c>) version of
    /// Klein &amp; Kneib's (2016): the common scale at which the prior's function draws exceed ±<paramref name="bound"/>
    /// somewhere on the rows of <paramref name="design"/> with probability <paramref name="tailProbability"/>,
    /// <c>P(‖f‖∞ ≥ bound) = α</c>. Every component's τⱼ² is sWⱼ with Wⱼ ~ W(k, 1), so a draw of the function's penalized
    /// part, <c>f = Zβ</c> with <c>β ~ N(0, K(τ)⁺)</c>, is √s times one drawn at s = 1, and s is <c>(bound/q)²</c> with q the
    /// (1 − α)-quantile of <c>max|f|</c> at s = 1, estimated from <paramref name="draws"/> draws (their λ = (q/c)² is its
    /// rate). The penalty's null space has a flat prior and is left out: the bound speaks to the penalized part only.
    /// <para>
    /// α = 0.05 is their default. They draw 1,000 times; 2,000 here, which only sharpens the quantile. Klein &amp; Kneib's
    /// original criterion is pointwise, <c>P(|z′β| ≤ c) ≥ 1 − α</c> at one point (c = 3, α = 0.01), computed by integration
    /// (their R package sdPrior); the simultaneous criterion here is the stricter statement over the design.
    /// </para>
    /// </summary>
    public static WeibullVariancePrior ScaledTo(double[,] design, IAnisotropicPenalty penalty, double bound,
                                                double tailProbability = 0.05, double shape = 0.5, int draws = 2000, int seed = 1)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(penalty);
        if (!(bound > 0.0) || double.IsInfinity(bound)) throw new ArgumentOutOfRangeException(nameof(bound), "The bound must be positive and finite.");
        if (!(tailProbability > 0.0 && tailProbability < 1.0)) throw new ArgumentOutOfRangeException(nameof(tailProbability));
        if (draws < 10) throw new ArgumentOutOfRangeException(nameof(draws));
        int n = design.GetLength(0), nu = design.GetLength(1), p = penalty.Components;
        Weibull unit = Weibull.ShapeScale(shape, 1.0);
        var rng = new Xoshiro256PlusPlus(seed);
        var weights = new double[p];
        var beta = new double[nu];
        var maxima = new double[draws];
        for (int d = 0; d < draws; d++)
        {
            for (int j = 0; j < p; j++) weights[j] = 1.0 / unit.Sample(rng);   // precisions 1/Wⱼ
            penalty.SamplePrior(weights, rng, beta);
            double max = 0.0;
            for (int i = 0; i < n; i++)
            {
                double f = 0.0;
                for (int k = 0; k < nu; k++) f += design[i, k] * beta[k];
                max = Math.Max(max, Math.Abs(f));
            }
            maxima[d] = max;
        }
        Array.Sort(maxima);
        double q = PosteriorInterval.Quantile(maxima, 1.0 - tailProbability);
        return new WeibullVariancePrior(bound * bound / (q * q), shape);
    }

    /// <summary>
    /// <see cref="ScaledTo"/> with Bach &amp; Klein's default bound, the observations' largest absolute value. They
    /// standardize y first, so in the response's own units the bound is <c>max |yᵢ − ȳ|</c>, which gives the same prior on
    /// the unstandardized smoothing variances (τ² scales with the response's variance).
    /// </summary>
    public static WeibullVariancePrior ScaledToResponse(double[,] design, IAnisotropicPenalty penalty, double[] response,
                                                        double tailProbability = 0.05, double shape = 0.5, int draws = 2000,
                                                        int seed = 1)
    {
        ArgumentNullException.ThrowIfNull(response);
        if (response.Length == 0) throw new ArgumentException("Need at least one observation.", nameof(response));
        double mean = 0.0;
        foreach (double v in response) mean += v / response.Length;
        double bound = 0.0;
        foreach (double v in response) bound = Math.Max(bound, Math.Abs(v - mean));
        if (!(bound > 0.0)) throw new ArgumentException("A constant response gives no bound to scale to.", nameof(response));
        return ScaledTo(design, penalty, bound, tailProbability, shape, draws, seed);
    }
}
