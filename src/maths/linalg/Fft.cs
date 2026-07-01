// src/maths/linalg/Fft.cs
#nullable enable
using System;
using System.Numerics;

namespace Maths.LinAlg;

/// <summary>
/// Radix-2 Cooley–Tukey fast Fourier transform over power-of-two lengths,
/// in place on <see cref="Complex"/> buffers, with a separable N-D driver.
/// A field→field transform carrying space/time samples to frequency
/// coefficients and back — the spectral primitive behind colored-noise
/// synthesis and convolution. Reference: Cooley &amp; Tukey 1965.
/// </summary>
/// <remarks>
/// Iterative bit-reversal + butterfly formulation (no recursion). The N-D
/// transform is separable — one 1-D pass along each axis in turn — for total
/// cost O(N log N) in the cell count N = ∏ dims. Every axis length must be a
/// power of two; grids are row-major (dims[0] outermost, dims[^1] contiguous).
/// The inverse applies the 1/length normalization per pass, so the round trip
/// is 1/N overall.
/// </remarks>
public static class Fft
{
    public static bool IsPowerOfTwo(int n) => n > 0 && (n & (n - 1)) == 0;

    /// <summary>Forward DFT, in place. Length must be a power of two.</summary>
    public static void Forward(Complex[] data) => Transform1D(data, inverse: false);

    /// <summary>Inverse DFT with 1/length normalization, in place.</summary>
    public static void Inverse(Complex[] data) => Transform1D(data, inverse: true);

    /// <summary>Forward 3-D DFT over a flat grid indexed (x·ny + y)·nz + z.</summary>
    public static void Forward3D(Complex[] grid, int nx, int ny, int nz)
        => TransformND(grid, new[] { nx, ny, nz }, inverse: false);

    /// <summary>Inverse 3-D DFT with 1/(nx·ny·nz) normalization.</summary>
    public static void Inverse3D(Complex[] grid, int nx, int ny, int nz)
        => TransformND(grid, new[] { nx, ny, nz }, inverse: true);

    /// <summary>Forward separable N-D DFT over a row-major grid (dims[0] outermost).</summary>
    public static void ForwardND(Complex[] grid, int[] dims) => TransformND(grid, dims, inverse: false);

    /// <summary>Inverse separable N-D DFT with 1/∏dims normalization.</summary>
    public static void InverseND(Complex[] grid, int[] dims) => TransformND(grid, dims, inverse: true);

    private static void Transform1D(Complex[] a, bool inverse)
    {
        int n = a.Length;
        if (!IsPowerOfTwo(n))
            throw new ArgumentException($"FFT length must be a power of two; got {n}.", nameof(a));

        // Bit-reversal permutation.
        for (int i = 1, j = 0; i < n; i++)
        {
            int bit = n >> 1;
            for (; (j & bit) != 0; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                (a[i], a[j]) = (a[j], a[i]);
        }

        // Butterflies.
        for (int len = 2; len <= n; len <<= 1)
        {
            double ang = 2.0 * Math.PI / len * (inverse ? 1.0 : -1.0);
            var wLen = new Complex(Math.Cos(ang), Math.Sin(ang));
            int half = len >> 1;
            for (int i = 0; i < n; i += len)
            {
                var w = Complex.One;
                for (int k = 0; k < half; k++)
                {
                    Complex u = a[i + k];
                    Complex v = a[i + k + half] * w;
                    a[i + k] = u + v;
                    a[i + k + half] = u - v;
                    w *= wLen;
                }
            }
        }

        if (inverse)
            for (int i = 0; i < n; i++)
                a[i] /= n;
    }

    private static void TransformND(Complex[] data, int[] dims, bool inverse)
    {
        if (dims.Length == 0)
            throw new ArgumentException("dims must be non-empty.", nameof(dims));

        int total = 1;
        foreach (int d in dims)
        {
            if (!IsPowerOfTwo(d))
                throw new ArgumentException($"Every FFT axis length must be a power of two; got {d}.", nameof(dims));
            total *= d;
        }
        if (data.Length != total)
            throw new ArgumentException($"Grid length {data.Length} != product of dims {total}.", nameof(data));

        var stride = new int[dims.Length];
        stride[^1] = 1;
        for (int i = dims.Length - 2; i >= 0; i--) stride[i] = stride[i + 1] * dims[i + 1];

        // One 1-D pass per axis, innermost (contiguous) → outermost. Each line
        // along axis `ax` has elements at b + k·stride[ax], k = 0..dims[ax]-1.
        for (int ax = dims.Length - 1; ax >= 0; ax--)
        {
            int s = stride[ax];
            int m = dims[ax];
            int block = m * s;
            int outer = total / block;
            var line = new Complex[m];
            for (int o = 0; o < outer; o++)
            {
                int blockBase = o * block;
                for (int inner = 0; inner < s; inner++)
                {
                    int b = blockBase + inner;
                    for (int k = 0; k < m; k++) line[k] = data[b + k * s];
                    Transform1D(line, inverse);
                    for (int k = 0; k < m; k++) data[b + k * s] = line[k];
                }
            }
        }
    }
}
