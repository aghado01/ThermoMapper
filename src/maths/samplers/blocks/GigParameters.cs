using System;
using System.Text.Json.Serialization;
using Maths.Distributions;

namespace Maths.Samplers.Blocks;

/// <summary>
/// The parameters of a generalized inverse Gaussian, GIG(p, a, b): density ∝ x^(p−1)·exp(−(a·x + b/x)/2) on x &gt; 0, with
/// a ≥ 0 and b ≥ 0, a &gt; 0 when p &gt; 0 and b &gt; 0 when p &lt; 0 (Barndorff-Nielsen's family; Jørgensen 1982). It is the
/// one form that the conjugate updates of a normal scale (<see cref="ConjugateUpdates"/>) are limits of, and the one that
/// generalized-hyperbolic noise's latent-scale block updates:
/// <list type="bullet">
/// <item>b → 0 is the <see cref="Gamma"/>, shape p, rate a/2 (<see cref="OfGamma"/>);</item>
/// <item>a → 0 is the <see cref="InverseGamma"/>, shape −p, scale b/2 (<see cref="OfInverseGamma"/>).</item>
/// </list>
/// A parameterization: its distribution and sampler are <see cref="GeneralizedInverseGaussian"/>'s (<see cref="ToDistribution"/>), and no
/// normalizing constant lives in either.
/// </summary>
public readonly record struct GigParameters(double P, double A, double B)
{
    /// <summary>The Gamma(<paramref name="shape"/>, <paramref name="rate"/>) as the limit b = 0: GIG(shape, 2·rate, 0).</summary>
    public static GigParameters OfGamma(double shape, double rate) => new(shape, 2.0 * rate, 0.0);

    /// <summary>The InverseGamma(<paramref name="shape"/>, <paramref name="scale"/>) as the limit a = 0: GIG(−shape, 0, 2·scale).</summary>
    public static GigParameters OfInverseGamma(double shape, double scale) => new(-shape, 0.0, 2.0 * scale);

    /// <summary>Whether this is the gamma limit: b = 0 with p, a &gt; 0. Read from the parameters, so not serialized.</summary>
    [JsonIgnore]
    public bool IsGamma => B == 0.0 && A > 0.0 && P > 0.0;

    /// <summary>Whether this is the inverse-gamma limit: a = 0 with p &lt; 0, b &gt; 0. Read from the parameters, so not serialized.</summary>
    [JsonIgnore]
    public bool IsInverseGamma => A == 0.0 && B > 0.0 && P < 0.0;

    /// <summary>The gamma of the limit b = 0.</summary>
    public Gamma ToGamma()
        => IsGamma ? Gamma.ShapeRate(P, 0.5 * A) : throw new InvalidOperationException("Only the limit b = 0 (p, a > 0) is a gamma.");

    /// <summary>The inverse gamma of the limit a = 0.</summary>
    public InverseGamma ToInverseGamma()
        => IsInverseGamma ? InverseGamma.ShapeScale(-P, 0.5 * B) : throw new InvalidOperationException("Only the limit a = 0 (p < 0, b > 0) is an inverse gamma.");

    /// <summary>The distribution these parameters name (<see cref="GeneralizedInverseGaussian"/>), whose sampler draws it.</summary>
    public GeneralizedInverseGaussian ToDistribution() => GeneralizedInverseGaussian.Of(P, A, B);

    /// <summary>ln of the unnormalized density, (p − 1)·ln x − (a·x + b/x)/2, for x &gt; 0.</summary>
    public double LogKernel(double x)
    {
        if (!(x > 0.0)) throw new ArgumentOutOfRangeException(nameof(x), "A GIG lives on x > 0.");
        return (P - 1.0) * Math.Log(x) - 0.5 * (A * x + B / x);
    }
}
