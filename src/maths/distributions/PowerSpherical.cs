using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The power spherical law PS(μ, κ) on the unit sphere S^(d−1) (De Cao and Aziz 2020, arXiv 2006.04437): density (1 + μᵀx)^κ / N
/// against surface measure, N = 2^(α+β) π^β Γ(α)/Γ(α + β) with α = (d − 1)/2 + κ and β = (d − 1)/2. Its cosine t = μᵀx has
/// (1 + t)/2 ~ Beta(α, β) (<see cref="CosineLaw"/>), and given t the rest of x is uniform on the (d − 2)-sphere orthogonal to μ, so a
/// draw is a beta variate, a uniform tangent and an orthogonal map carrying e₁ to μ, with no rejection. The proposal for a ridge
/// direction's active coordinates in BPPR's change move (Collins, Francom and Rumsey, arXiv 2210.09181, κ = 1000): the density is
/// symmetric in μ and x, so the proposal cancels from the move's ratio. On S⁰ (d = 1) it is the point mass at μ for κ &gt; 0, and
/// uniform on ±1 at κ = 0.
/// <para>Numerically: the beta variate is read as two gamma draws, so 1 − t = 2G_β/(G_α + G_β) and √(1 − t²) keep their digits
/// where κ pushes t toward 1; the map is the reflection along e₁ − μ, or minus the one along e₁ + μ when μ₁ &gt; 0, so its vector never
/// subtracts nearly equal numbers; the density reads 1 + μᵀx as ‖x + μ‖²/2, exact near the antipode and symmetric in μ and x bit for
/// bit; and the normalizer takes Γ(α + β)/Γ(α) from <see cref="SpecialFunctions.LogGammaRatio"/>, which holds its digits at κ = 10⁵.</para>
/// Immutable.
/// </summary>
public sealed class PowerSpherical
{
    // |‖v‖² − 1| allowed of a mean or a point said to lie on the sphere.
    private const double UnitTolerance = 1e-10;

    private readonly double[] _mean;
    private readonly double _alpha;
    private readonly double _beta;

    private PowerSpherical(double[] mean, double concentration)
    {
        _mean = mean;
        Concentration = concentration;
        _beta = 0.5 * (mean.Length - 1);
        _alpha = _beta + concentration;
        LogNormalizer = _alpha > 0.0
            ? (_alpha + _beta) * SpecialFunctions.Ln2 + _beta * Math.Log(Math.PI) - SpecialFunctions.LogGammaRatio(_alpha, _beta)
            : SpecialFunctions.Ln2;   // S⁰ at κ = 0: two points of mass ½
    }

    /// <summary>PS(<paramref name="mean"/>, <paramref name="concentration"/>): a unit mean (to 10⁻¹⁰ in its squared norm, kept as
    /// given) and a finite κ ≥ 0, κ = 0 being the uniform law; copied.</summary>
    public static PowerSpherical Of(ReadOnlySpan<double> mean, double concentration)
    {
        if (mean.Length == 0) throw new ArgumentException("The mean is a point of a sphere in at least one dimension.", nameof(mean));
        if (!(concentration >= 0.0) || double.IsPositiveInfinity(concentration))
            throw new ArgumentOutOfRangeException(nameof(concentration), "The concentration is finite and at least 0.");
        CheckUnit(mean, nameof(mean));
        return new PowerSpherical(mean.ToArray(), concentration);
    }

    /// <summary>d, the ambient dimension: the law lives on S^(d−1).</summary>
    public int Dimension => _mean.Length;

    /// <summary>κ.</summary>
    public double Concentration { get; }

    /// <summary>μ, the mode for κ &gt; 0.</summary>
    public ReadOnlySpan<double> Mean => _mean;

    /// <summary>ln N, the log of the normalizer against surface measure (counting measure on S⁰).</summary>
    public double LogNormalizer { get; }

    /// <summary>The law of (1 + μᵀx)/2, Beta((d − 1)/2 + κ, (d − 1)/2); refused on S⁰, where the cosine is ±1.</summary>
    public Beta CosineLaw => Dimension > 1
        ? Beta.Shapes(_alpha, _beta)
        : throw new InvalidOperationException("On S⁰ the cosine is ±1 and has no beta law.");

    /// <summary>ln of the density at a unit <paramref name="x"/> (to 10⁻¹⁰ in its squared norm): κ ln(‖x + μ‖²/2) − ln N, −∞ at the
    /// antipode for κ &gt; 0.</summary>
    public double LogPdf(ReadOnlySpan<double> x)
    {
        if (x.Length != Dimension) throw new ArgumentException($"The point has {x.Length} coordinates, the law {Dimension}.", nameof(x));
        CheckUnit(x, nameof(x));
        if (Concentration == 0.0) return -LogNormalizer;
        double sum = 0.0;
        for (int j = 0; j < x.Length; j++)
        {
            double v = x[j] + _mean[j];
            sum += v * v;
        }
        return Concentration * Math.Log(0.5 * sum) - LogNormalizer;
    }

    /// <summary>One draw, written to <paramref name="destination"/> (length d).</summary>
    public void Sample(Xoshiro256PlusPlus rng, Span<double> destination)
    {
        ArgumentNullException.ThrowIfNull(rng);
        int d = Dimension;
        if (destination.Length != d) throw new ArgumentException($"The destination has {destination.Length} coordinates, the law {d}.", nameof(destination));
        if (d == 1)
        {
            destination[0] = Concentration > 0.0 ? _mean[0] : rng.NextDouble() < 0.5 ? -1.0 : 1.0;
            return;
        }

        // t = z − w and √(1 − t²) = 2√(zw), z = G_α/(G_α + G_β) and w its complement, from the gammas scaled by the larger.
        double ga = Gamma.SampleStandard(rng, _alpha), gb = Gamma.SampleStandard(rng, _beta);
        double top = Math.Max(ga, gb), t = 0.0, radius = 1.0;
        if (top > 0.0)
        {
            double a = ga / top, b = gb / top, sum = a + b, z = a / sum, w = b / sum;
            t = z - w;
            radius = 2.0 * Math.Sqrt(z * w);
        }

        Span<double> y = d <= 64 ? stackalloc double[d] : new double[d];
        double norm;
        do
        {
            norm = 0.0;
            for (int j = 1; j < d; j++)
            {
                y[j] = Normal.SampleStandard(rng);
                norm += y[j] * y[j];
            }
        }
        while (!(norm > 0.0));
        double scale = radius / Math.Sqrt(norm);
        y[0] = t;
        for (int j = 1; j < d; j++) y[j] *= scale;

        // x = ρ(I − 2uuᵀ/‖u‖²)y: u = e₁ − μ and ρ = 1 when μ₁ ≤ 0, u = e₁ + μ and ρ = −1 when μ₁ > 0, so x = μ at y = e₁ and
        // ‖u‖² ≥ 2. Written as u₁ = 1 + |μ₁| and u_j = σμ_j, σ = −ρ.
        double sigma = _mean[0] > 0.0 ? 1.0 : -1.0, rho = -sigma;
        double u0 = 1.0 + Math.Abs(_mean[0]), uu = u0 * u0, uy = u0 * y[0];
        for (int j = 1; j < d; j++)
        {
            double uj = sigma * _mean[j];
            uu += uj * uj;
            uy += uj * y[j];
        }
        double c = 2.0 * uy / uu;
        destination[0] = rho * (y[0] - c * u0);
        for (int j = 1; j < d; j++) destination[j] = rho * (y[j] - c * sigma * _mean[j]);
    }

    /// <summary>One draw.</summary>
    public double[] Sample(Xoshiro256PlusPlus rng)
    {
        var x = new double[Dimension];
        Sample(rng, x);
        return x;
    }

    private static void CheckUnit(ReadOnlySpan<double> v, string name)
    {
        double norm = 0.0;
        foreach (double c in v)
        {
            if (!double.IsFinite(c)) throw new ArgumentException("A point of the sphere has finite coordinates.", name);
            norm += c * c;
        }
        if (!(Math.Abs(norm - 1.0) <= UnitTolerance)) throw new ArgumentException($"A point of the sphere has unit norm; its squared norm is {norm:R}.", name);
    }
}
