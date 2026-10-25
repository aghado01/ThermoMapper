using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Samplers.Blocks;

/// <summary>
/// A Gaussian conditional draw from its precision factor: for <c>x | · ~ N(Q⁻¹b, s²Q⁻¹)</c> with Q factored, the mean
/// solves <c>Qμ = b</c> and the draw adds <c>s·v</c>, v ~ N(0, Q⁻¹), from the factor's innovation. With a banded Q this
/// is the block update of a Gaussian Markov random field (Rue 2001), as in Bayesian P-splines (Lang &amp; Brezger 2004);
/// a dense factor serves the same call. The mean comes back beside the draw, since a Rao–Blackwellized estimate
/// averages the means. An instance holds one dimension's scratch, so each chain owns its own.
/// </summary>
public sealed class GaussianConditional
{
    private readonly double[] _z;
    private readonly double[] _innovation;

    public GaussianConditional(int dimension)
    {
        if (dimension < 1) throw new ArgumentOutOfRangeException(nameof(dimension));
        _z = new double[dimension];
        _innovation = new double[dimension];
    }

    /// <summary>The dimension the scratch serves.</summary>
    public int Dimension => _z.Length;

    /// <summary>Solve <c>Qμ = rhs</c>, then write <c>μ + s·v</c> into <paramref name="result"/>
    /// (<see cref="Perturb"/>). Returns μ, a fresh array.</summary>
    public double[] Draw(ISpdFactor precision, ReadOnlySpan<double> rhs, Xoshiro256PlusPlus rng, Span<double> result,
                         double scale = 1.0)
    {
        ArgumentNullException.ThrowIfNull(precision);
        double[] mean = precision.Solve(rhs);
        Perturb(precision, mean, rng, result, scale);
        return mean;
    }

    /// <summary>Write <c>mean + s·v</c> into <paramref name="result"/>, with v ~ N(0, Q⁻¹): <see cref="Dimension"/>
    /// standard normals in order, then the factor's innovation.</summary>
    public void Perturb(ISpdFactor precision, ReadOnlySpan<double> mean, Xoshiro256PlusPlus rng, Span<double> result,
                        double scale = 1.0)
    {
        int p = _z.Length;
        if (mean.Length != p || result.Length != p)
            throw new ArgumentException($"The mean and the result need {p} entries.", nameof(mean));
        Perturb(precision, mean, rng, result, scale, _z, _innovation);
    }

    /// <summary>
    /// <see cref="Perturb(ISpdFactor, ReadOnlySpan{double}, Xoshiro256PlusPlus, Span{double}, double)"/> on the caller's
    /// scratch: <paramref name="z"/> receives the standard normals and <paramref name="innovation"/> v. For a caller that
    /// keeps its scratch on the stack, or reads the normals back, as an independence proposal's density −½zᵀz does.
    /// </summary>
    public static void Perturb(ISpdFactor precision, ReadOnlySpan<double> mean, Xoshiro256PlusPlus rng, Span<double> result,
                               double scale, Span<double> z, Span<double> innovation)
    {
        ArgumentNullException.ThrowIfNull(precision);
        ArgumentNullException.ThrowIfNull(rng);
        int p = mean.Length;
        if (result.Length != p || z.Length != p || innovation.Length != p)
            throw new ArgumentException($"The result and the scratch need {p} entries, as the mean has.", nameof(result));
        for (int j = 0; j < p; j++) z[j] = Normal.SampleStandard(rng);
        precision.SampleInnovation(z, innovation);
        for (int j = 0; j < p; j++) result[j] = mean[j] + scale * innovation[j];
    }
}
