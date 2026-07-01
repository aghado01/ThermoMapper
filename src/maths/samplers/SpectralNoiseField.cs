using System;
using System.Numerics;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Samplers;

/// <summary>
/// Gaussian random field on a cubic N-D grid with a power-law spatial spectrum
/// (PSD ∝ |k|^−β): white (β=0), pink (β=1), or brown/red (β=2) "colored"
/// noise. Samples white noise, shapes its spectrum through <see cref="Fft"/>,
/// and exposes multilinear point queries — a sampler that materializes a
/// correlated field downstream generators reduce against (e.g. a background
/// density blanket, in 3-D or higher). Color presets: <see cref="White"/>,
/// <see cref="Pink"/>, <see cref="Brown"/>.
/// </summary>
/// <remarks>
/// Build: fill the grid with i.i.d. N(0,1) (Box–Muller off
/// <see cref="Xoshiro256PlusPlus"/>), forward FFT, multiply each bin by
/// |k|^(−β/2) so the power spectrum |F|² scales as |k|^−β, inverse FFT, keep
/// the real part, then standardize to zero mean / unit variance. The DC bin is
/// zeroed (the β-pole and the field mean). Larger β piles power into low
/// frequencies → smoother, lumpier fields. The dimension is the length of the
/// box passed to <see cref="Generate"/>; the synthesized grid is periodic and
/// <see cref="Sample(ReadOnlySpan{double})"/> clamps queries to the box faces.
/// </remarks>
public sealed class SpectralNoiseField
{
    public const double White = 0.0;
    public const double Pink = 1.0;
    public const double Brown = 2.0;

    private readonly int[] _dims;
    private readonly int[] _stride;
    private readonly double[] _min;   // box lower corner per axis
    private readonly double[] _span;  // box extent per axis
    private readonly double[] _field; // standardized real field, length ∏ dims

    /// <summary>Ambient dimension of the field.</summary>
    public int Dimension => _dims.Length;

    /// <summary>Maximum field value over the grid (for rejection envelopes).</summary>
    public double Max { get; }

    /// <summary>Minimum field value over the grid.</summary>
    public double Min { get; }

    private SpectralNoiseField(
        int[] dims, int[] stride, double[] min, double[] span, double[] field, double lo, double hi)
    {
        _dims = dims; _stride = stride; _min = min; _span = span; _field = field;
        Min = lo; Max = hi;
    }

    /// <summary>
    /// Synthesize a colored-noise field of side <paramref name="gridSize"/>
    /// (power of two) per axis over the axis-aligned box [<paramref name="boxMin"/>,
    /// <paramref name="boxMax"/>] (whose length sets the dimension), with spectral
    /// exponent <paramref name="beta"/>.
    /// </summary>
    public static SpectralNoiseField Generate(
        Xoshiro256PlusPlus rng, int gridSize, double beta, double[] boxMin, double[] boxMax)
    {
        ArgumentNullException.ThrowIfNull(rng);
        ArgumentNullException.ThrowIfNull(boxMin);
        ArgumentNullException.ThrowIfNull(boxMax);
        if (!Fft.IsPowerOfTwo(gridSize))
            throw new ArgumentException($"gridSize must be a power of two; got {gridSize}.", nameof(gridSize));
        if (boxMin.Length != boxMax.Length || boxMin.Length == 0)
            throw new ArgumentException("boxMin and boxMax must be non-empty and equal length.");

        int d = boxMin.Length;
        var dims = new int[d];
        for (int i = 0; i < d; i++) dims[i] = gridSize;

        var stride = new int[d];
        stride[d - 1] = 1;
        for (int i = d - 2; i >= 0; i--) stride[i] = stride[i + 1] * dims[i + 1];

        int total = 1;
        for (int i = 0; i < d; i++) total *= gridSize;

        var grid = new Complex[total];
        for (int i = 0; i < total; i++)
            grid[i] = new Complex(NextGaussian(rng), 0.0);

        Fft.ForwardND(grid, dims);

        // Power-law shaping: bin k scaled by |k|^(−β/2) = (k²)^(−β/4); DC → 0.
        int half = gridSize / 2;
        for (int idx = 0; idx < total; idx++)
        {
            double k2 = 0.0;
            int rem = idx;
            for (int ax = 0; ax < d; ax++)
            {
                int digit = rem / stride[ax];
                rem -= digit * stride[ax];
                double f = digit <= half ? digit : digit - gridSize;
                k2 += f * f;
            }
            double filter = k2 == 0.0 ? 0.0 : Math.Pow(k2, -beta / 4.0);
            grid[idx] *= filter;
        }

        Fft.InverseND(grid, dims);

        var field = new double[total];
        double mean = 0.0;
        for (int i = 0; i < total; i++) { field[i] = grid[i].Real; mean += field[i]; }
        mean /= total;

        double variance = 0.0;
        for (int i = 0; i < total; i++) { field[i] -= mean; variance += field[i] * field[i]; }
        variance /= total;
        double sd = Math.Sqrt(variance);
        if (sd < 1e-300) sd = 1.0;

        double lo = double.PositiveInfinity, hi = double.NegativeInfinity;
        for (int i = 0; i < total; i++)
        {
            field[i] /= sd;
            if (field[i] < lo) lo = field[i];
            if (field[i] > hi) hi = field[i];
        }

        var min = new double[d];
        var span = new double[d];
        for (int i = 0; i < d; i++) { min[i] = boxMin[i]; span[i] = boxMax[i] - boxMin[i]; }
        return new SpectralNoiseField(dims, stride, min, span, field, lo, hi);
    }

    /// <summary>
    /// Multilinearly interpolated field value at a world point; coordinates
    /// outside the box are clamped to its faces. The point must have at least
    /// <see cref="Dimension"/> components.
    /// </summary>
    public double Sample(ReadOnlySpan<double> point)
    {
        int d = _dims.Length;
        if (point.Length < d)
            throw new ArgumentException($"point must have at least {d} components.", nameof(point));

        Span<int> lo = d <= 16 ? stackalloc int[d] : new int[d];
        Span<double> frac = d <= 16 ? stackalloc double[d] : new double[d];
        for (int ax = 0; ax < d; ax++)
        {
            double t = _span[ax] <= 0.0 ? 0.0 : (point[ax] - _min[ax]) / _span[ax];
            t = Math.Clamp(t, 0.0, 1.0) * (_dims[ax] - 1);
            int l = (int)Math.Floor(t);
            if (l >= _dims[ax] - 1) { l = _dims[ax] - 1; frac[ax] = 0.0; }
            else frac[ax] = t - l;
            lo[ax] = l;
        }

        // Multilinear interpolation over the 2^d corners of the cell.
        double acc = 0.0;
        int corners = 1 << d;
        for (int mask = 0; mask < corners; mask++)
        {
            double w = 1.0;
            int flat = 0;
            for (int ax = 0; ax < d; ax++)
            {
                int bit = (mask >> ax) & 1;
                int coord = lo[ax] + bit;
                if (coord >= _dims[ax]) coord = _dims[ax] - 1;
                w *= bit == 1 ? frac[ax] : 1.0 - frac[ax];
                flat += coord * _stride[ax];
            }
            if (w != 0.0) acc += w * _field[flat];
        }
        return acc;
    }

    /// <summary>3-D convenience overload of <see cref="Sample(ReadOnlySpan{double})"/>.</summary>
    public double Sample(double x, double y, double z)
    {
        Span<double> p = stackalloc double[3] { x, y, z };
        return Sample(p);
    }

    private static double NextGaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
