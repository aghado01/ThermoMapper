using System;
using System.Globalization;
using Maths.Distributions;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// BARS's noise beyond the Gaussian (OQ20, E1's slice 6): a normal variance-mean mixture, εᵢ = σ(β·vᵢ + √(c·vᵢ)·ζᵢ) with each vᵢ from
/// its mixing law (<see cref="IMixingLaw"/>, the law BASS's noises share), so the noise is Gaussian given its latent. The asymmetric
/// Laplace with its skew fixed (<see cref="Quantile"/>) makes the fitted curve the q-quantile; a drawn skew (<see cref="Skewed"/>) lets the
/// data choose it.
/// <list type="bullet">
/// <item><b>A pure scale mixture</b> (β = 0, not drawn): σ² stays integrated in the knot step and is drawn with the coefficients, robust
/// BARS's sweep with the weights 1/(c·vᵢ) (<see cref="StudentTWeights"/> is <see cref="StudentT"/>'s weights alone).</item>
/// <item><b>A skew</b>, held away from 0 or drawn: σ² sits in the offset σβ·v, so it is held on the chain and drawn by its own block, and
/// the knot step reads the marginal given σ² on the working response y − σβ·v.</item>
/// </list>
/// The noise block, on the weight block's cadence, draws the coefficients given the knots, then each vᵢ by its law's step, then σ²
/// under a skew from its modified half-normal and a drawn β from its normal (<see cref="VarianceMeanConditionals"/>), GBASS's conditionals
/// as BASS draws them. Its priors read the data once, not the working response the noise moves.
/// </summary>
public sealed class BarsNoise
{
    private BarsNoise(IMixingLaw mixing, double varianceScale, double skew, Normal? skewPrior, string name)
    {
        ArgumentNullException.ThrowIfNull(mixing);
        if (!(varianceScale > 0.0 && double.IsFinite(varianceScale)))
            throw new ArgumentOutOfRangeException(nameof(varianceScale), "c must be positive and finite.");
        if (!double.IsFinite(skew)) throw new ArgumentOutOfRangeException(nameof(skew), "The skew must be finite.");
        Mixing = mixing;
        VarianceScale = varianceScale;
        Skew = skew;
        SkewPrior = skewPrior;
        Name = name;
    }

    /// <summary>v's mixing law.</summary>
    public IMixingLaw Mixing { get; }

    /// <summary>c, the variance per unit v and σ².</summary>
    public double VarianceScale { get; }

    /// <summary>β when held; a drawn skew's start, its prior's mean.</summary>
    public double Skew { get; }

    /// <summary>β's normal prior when the skew is drawn; null when it is held.</summary>
    public Normal? SkewPrior { get; }

    /// <summary>What the noise is, for a run's fingerprint and its reports.</summary>
    public string Name { get; }

    /// <summary>Whether σ² sits in the offset: a skew held away from 0, or drawn.</summary>
    public bool HasSkew => Skew != 0.0 || SkewPrior is not null;

    /// <summary>Student t on ν = <paramref name="degreesOfFreedom"/>, v ~ IG(ν/2, ν/2), robust BARS with v kept; c = 1 by default, so σ is t's
    /// scale, or (ν − 2)/ν, under which σ² is the noise's variance.</summary>
    public static BarsNoise StudentT(double degreesOfFreedom, double varianceScale = 1.0)
        => new(new StudentTMixing(degreesOfFreedom), varianceScale, 0.0, null, $"t({Text(degreesOfFreedom)})");

    /// <summary>
    /// The asymmetric Laplace at quantile level <paramref name="q"/> in (0, 1) (GBASS's <c>qbass</c>; Kozumi and Kobayashi 2011): v ~ Exp(1),
    /// c = 2/(q(1 − q)) and β = (1 − 2q)/(q(1 − q)) fixed, so the noise's q-quantile is 0 and the fitted curve is the q-quantile curve.
    /// </summary>
    public static BarsNoise Quantile(double q)
    {
        if (!(q > 0.0 && q < 1.0)) throw new ArgumentOutOfRangeException(nameof(q), "The quantile level lies in (0, 1).");
        double spread = q * (1.0 - q);
        return new(ExponentialMixing.Instance, 2.0 / spread, (1.0 - 2.0 * q) / spread, null, $"AL({Text(q)})");
    }

    /// <summary>A pure scale mixture on <paramref name="mixing"/>: β = 0, σ² integrated in the knot step.</summary>
    public static BarsNoise ScaleMixture(IMixingLaw mixing, double varianceScale = 1.0)
        => new(mixing, varianceScale, 0.0, null, Describe(mixing));

    /// <summary>A skew held at <paramref name="skew"/> on <paramref name="mixing"/>.</summary>
    public static BarsNoise HeldSkew(IMixingLaw mixing, double skew, double varianceScale = 1.0)
        => new(mixing, varianceScale, skew, null, $"{Describe(mixing)}, β {Text(skew)}");

    /// <summary>A skew drawn under the normal prior <paramref name="skewPrior"/>, starting at its mean, on <paramref name="mixing"/>.</summary>
    public static BarsNoise Skewed(IMixingLaw mixing, Normal skewPrior, double varianceScale = 1.0)
    {
        ArgumentNullException.ThrowIfNull(skewPrior);
        if (!(skewPrior.Sd > 0.0 && double.IsFinite(skewPrior.Sd) && double.IsFinite(skewPrior.Mean)))
            throw new ArgumentOutOfRangeException(nameof(skewPrior), "The skew's prior needs a finite mean and a positive finite sd.");
        return new(mixing, varianceScale, skewPrior.Mean, skewPrior, $"{Describe(mixing)}, β ~ N({Text(skewPrior.Mean)}, {Text(skewPrior.Sd)}²)");
    }

    /// <summary>
    /// The noise's log CDF at a draw, e ↦ ln F(ε = e | σ², β), v integrated out and each tail kept (E16, slice 5): Student t's closed form
    /// without a skew, the asymmetric Laplace's under exponential mixing, and otherwise the mixture's own law, the generalized hyperbolic or
    /// the normal–generalized-beta-prime, built once per draw. It is the law <see cref="IMixingLaw.MarginalLogDensity"/> reads at
    /// s² = σ²c and m = σβ.
    /// </summary>
    internal Func<double, double> LogCdfAt(double sigma2, double skew)
    {
        double unitVariance = sigma2 * VarianceScale, drift = Math.Sqrt(sigma2) * skew;
        IContinuousDistribution law = Mixing switch
        {
            StudentTMixing t when drift == 0.0 => Maths.Distributions.StudentT.LocationScale(t.DegreesOfFreedom, 0.0, Math.Sqrt(unitVariance)),
            StudentTMixing t => new GigMixing(t.Prior).Noise(unitVariance, drift),
            ExponentialMixing => AsymmetricLaplaceOf(unitVariance, drift),
            GigMixing mixing => mixing.Noise(unitVariance, drift),
            GbpMixing mixing => mixing.Noise(unitVariance, drift),
            _ => throw new NotSupportedException($"{Mixing.GetType().Name} gives no noise distribution, so no CDF."),
        };
        return e => NoiseLogCdf.Of(law, e);
    }

    // Exponential mixing's marginal, (m² + 2s²)^(−½)·exp(m·r/s² − |r|·√(m² + 2s²)/s²), as ALD(0, 1/(a + b), b/(a + b)) with left rate
    // a = (√(m² + 2s²) + m)/s² and right rate b = (√(m² + 2s²) − m)/s²; the smaller rate is taken as 2/(√(m² + 2s²) + |m|), with no cancellation.
    private static AsymmetricLaplace AsymmetricLaplaceOf(double unitVariance, double drift)
    {
        double root = Math.Sqrt(drift * drift + 2.0 * unitVariance);
        double small = 2.0 / (root + Math.Abs(drift)), large = (root + Math.Abs(drift)) / unitVariance;
        (double left, double right) = drift >= 0.0 ? (large, small) : (small, large);
        return AsymmetricLaplace.Of(0.0, 1.0 / (left + right), right / (left + right));
    }

    /// <summary>The noise's description with c: what a run's fingerprint takes.</summary>
    public string Fingerprint => $"{Name}; c {Text(VarianceScale)}";

    /// <summary>A start for σ² under a skew: the response's variance over the noise's per unit σ², β²Var[v] + c·E[v], or c where v's law
    /// has no finite moments (BASS's start).</summary>
    internal double StartSigma2(double[] y)
    {
        double mean = 0.0, squares = 0.0;
        foreach (double v in y) mean += v;
        mean /= y.Length;
        foreach (double v in y) squares += (v - mean) * (v - mean);
        double unit = Skew * Skew * Mixing.Variance + VarianceScale * Mixing.Mean;
        if (!(unit > 0.0 && double.IsFinite(unit))) unit = VarianceScale;
        double start = squares / y.Length / unit;
        return start > 0.0 && double.IsFinite(start) ? start : 1.0;
    }

    private static string Describe(IMixingLaw mixing) => mixing switch
    {
        StudentTMixing t => $"t({Text(t.DegreesOfFreedom)})",
        ExponentialMixing => "Exp(1)",
        GigMixing g => $"GIG({Text(g.Prior.P)}, {Text(g.Prior.A)}, {Text(g.Prior.B)})",
        GbpMixing b => $"GBP({Text(b.Prior.P)}, {Text(b.Prior.A)}, {Text(b.Prior.B)})",
        _ => mixing.GetType().FullName!,
    };

    private static string Text(double value) => value.ToString("R", CultureInfo.InvariantCulture);
}
