using System;
using Maths.Rng;

namespace Maths.Regression.Spline.Bars;

/// <summary>Which base measure a BARS birth's exploratory component draws from (<see cref="KnotBirthDensity"/>).</summary>
public enum BirthBase
{
    /// <summary>Uniform on (0, 1): every location equally (<see cref="UniformKnotBase"/>).</summary>
    Uniform,

    /// <summary>The inputs' own density, a histogram with a uniform floor (<see cref="DesignHistogramBase"/>).</summary>
    Design,
}

/// <summary>A base measure on (0, 1) for new knots: the exploratory births draw from it, and so does the first knot.</summary>
public interface IKnotBase
{
    /// <summary>Draw a location.</summary>
    double Sample(Xoshiro256PlusPlus rng);

    /// <summary>Log density at <paramref name="x"/>.</summary>
    double LogDensity(double x);
}

/// <summary>The uniform base on (0, 1): the birth BARS has always made from an empty configuration.</summary>
public sealed class UniformKnotBase : IKnotBase
{
    /// <summary>The one instance; it holds no state.</summary>
    public static readonly UniformKnotBase Instance = new();

    private UniformKnotBase() { }

    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        return rng.NextDouble();
    }

    public double LogDensity(double x) => x >= 0.0 && x <= 1.0 ? 0.0 : double.NegativeInfinity;
}

/// <summary>
/// The inputs' own density as a base: a histogram of the design points over equal bins of (0, 1), mixed with the
/// uniform at weight <see cref="Floor"/> so that every location keeps a positive density. It is the continuous
/// analogue of BASS's knots at observed values: exploratory births land where the data are, which is where the
/// support gate admits them.
/// </summary>
public sealed class DesignHistogramBase : IKnotBase
{
    private readonly double[] _cumulative;   // cumulative bin probabilities of the histogram part
    private readonly double[] _logDensity;   // log of the mixed density in each bin

    /// <param name="x">The design points, each in [0, 1].</param>
    /// <param name="bins">Equal bins on (0, 1); 0 chooses ⌈√n⌉, between 5 and 100.</param>
    /// <param name="floor">The uniform's weight in the mixture, in (0, 1].</param>
    public DesignHistogramBase(double[] x, int bins = 0, double floor = 0.1)
    {
        ArgumentNullException.ThrowIfNull(x);
        if (x.Length == 0) throw new ArgumentException("The design needs at least one point.", nameof(x));
        if (bins < 0) throw new ArgumentOutOfRangeException(nameof(bins));
        if (!(floor > 0.0 && floor <= 1.0)) throw new ArgumentOutOfRangeException(nameof(floor), "The floor must lie in (0, 1].");
        int b = bins > 0 ? bins : Math.Clamp((int)Math.Ceiling(Math.Sqrt(x.Length)), 5, 100);
        var counts = new double[b];
        foreach (double v in x)
        {
            if (!(v >= 0.0 && v <= 1.0)) throw new ArgumentOutOfRangeException(nameof(x), "Design points must lie in [0, 1].");
            counts[Math.Min((int)(v * b), b - 1)] += 1.0;
        }
        Bins = b;
        Floor = floor;
        _cumulative = new double[b];
        _logDensity = new double[b];
        double running = 0.0;
        for (int i = 0; i < b; i++)
        {
            running += counts[i] / x.Length;
            _cumulative[i] = running;
            _logDensity[i] = Math.Log(floor + (1.0 - floor) * counts[i] / x.Length * b);   // a bin of width 1/b
        }
        _cumulative[b - 1] = 1.0;
    }

    /// <summary>The number of equal bins.</summary>
    public int Bins { get; }

    /// <summary>The uniform's weight in the mixture.</summary>
    public double Floor { get; }

    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (rng.NextDouble() < Floor) return rng.NextDouble();
        double u = rng.NextDouble();
        int bin = Array.BinarySearch(_cumulative, u);
        bin = bin >= 0 ? bin : ~bin;
        if (bin >= Bins) bin = Bins - 1;
        return (bin + rng.NextDouble()) / Bins;
    }

    public double LogDensity(double x)
        => x >= 0.0 && x <= 1.0 ? _logDensity[Math.Min((int)(x * Bins), Bins - 1)] : double.NegativeInfinity;
}

/// <summary>
/// The birth proposal's density for a new knot given the current ones (preparation step 3.9, OQ11). With probability
/// ε the knot comes from a base measure; otherwise it perturbs an existing knot chosen uniformly (DMGK's locality
/// proposal); with no knots it comes from the base. So
/// <c>q(t | ξ) = ε·b(t) + (1 − ε)·(1/k)·Σᵢ h(t | ξᵢ)</c>. A death's acceptance carries this density at the removed knot,
/// so ε &gt; 0 bounds it below by ε·b(t), and a knot far from the others can still be removed. At ε = 0 births are
/// local only, and such a knot is stranded. Birth and death share one instance, so the two always agree.
/// </summary>
public sealed class KnotBirthDensity
{
    private readonly double _logExploration, _logLocal;

    /// <param name="kernel">The local kernel h, shared with the relocation move.</param>
    /// <param name="exploration">ε ∈ [0, 1], the share of births drawn from the base (0: local births only).</param>
    /// <param name="knotBase">The base measure b; uniform on (0, 1) when null.</param>
    public KnotBirthDensity(IKnotKernel kernel, double exploration = 0.0, IKnotBase? knotBase = null)
    {
        ArgumentNullException.ThrowIfNull(kernel);
        if (!(exploration >= 0.0 && exploration <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(exploration), "The exploratory share must lie in [0, 1].");
        Kernel = kernel;
        Exploration = exploration;
        Base = knotBase ?? UniformKnotBase.Instance;
        _logExploration = Math.Log(exploration);
        _logLocal = Math.Log(1.0 - exploration);
    }

    /// <summary>The local kernel h.</summary>
    public IKnotKernel Kernel { get; }

    /// <summary>ε, the share of births drawn from the base.</summary>
    public double Exploration { get; }

    /// <summary>The base measure b.</summary>
    public IKnotBase Base { get; }

    /// <summary>Draw a new knot given the current ones.</summary>
    public double Sample(double[] knots, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(knots);
        ArgumentNullException.ThrowIfNull(rng);
        int k = knots.Length;
        if (k == 0) return Base.Sample(rng);
        if (Exploration > 0.0 && rng.NextDouble() < Exploration) return Base.Sample(rng);
        return Kernel.Sample(knots[rng.NextInt(k)], rng);
    }

    /// <summary>Log density of proposing <paramref name="t"/> given <paramref name="knots"/>.</summary>
    public double LogDensity(double t, double[] knots)
    {
        ArgumentNullException.ThrowIfNull(knots);
        if (knots.Length == 0) return Base.LogDensity(t);
        if (Exploration == 0.0) return ProposalMath.LogMixtureDensity(t, knots, Kernel);
        double exploratory = _logExploration + Base.LogDensity(t);
        if (Exploration == 1.0) return exploratory;
        double local = _logLocal + ProposalMath.LogMixtureDensity(t, knots, Kernel);
        double max = Math.Max(exploratory, local);
        if (double.IsNegativeInfinity(max)) return double.NegativeInfinity;
        return max + Math.Log(Math.Exp(exploratory - max) + Math.Exp(local - max));
    }
}
