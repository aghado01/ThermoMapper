using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace Maths.LinAlg;

/// <summary>
/// The dense product C = A·B of row-major matrices, A rows × inner and B inner × columns. Each entry is its inner sum taken in
/// order by fused multiply-adds, c ← fma(a_ik, b_kj, c) for k = 0, 1, …, from c = +0, so every entry is the scalar loop's,
/// bit for bit, on any hardware (<c>Vector256</c>'s software fallback fuses too), and two equal columns of B give two equal
/// columns of C. <see cref="FusedDot"/> is the same sum for one entry. B is packed once (<see cref="PackedColumns"/>, panels of
/// eight columns) and serves any number of left operands. The kernel holds a 4 × 8 tile of C in registers and runs the
/// inner index in chunks of <see cref="Chunk"/>, so a chunk of B's panels stays in cache across the tiles that read it; the
/// chunks continue each entry's sum in order, so they change no bit. Compiled fully optimized from the first call
/// (<see cref="MethodImplOptions.AggressiveOptimization"/>): tiered compilation would start the vector loops in unoptimized
/// code. Single-threaded; the caller splits rows or columns across threads, which changes no entry.
/// </summary>
public static class DenseProduct
{
    /// <summary>The columns of one packed panel, two 256-bit vectors.</summary>
    public const int Width = 8;

    /// <summary>The inner indices one pass of the kernel runs before moving to the next tile.</summary>
    public const int Chunk = 256;

    private const int Rows = 4;

    /// <summary>B, inner × columns row-major, packed into panels of <see cref="Width"/> columns.</summary>
    public static PackedColumns Pack(ReadOnlySpan<double> b, int inner, int columns)
    {
        var packed = new PackedColumns(inner, columns);
        packed.Fill(b);
        return packed;
    }

    /// <summary>C = A·B into <paramref name="c"/>, rows × B's columns row-major, with A rows × B's inner row-major.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static void Multiply(ReadOnlySpan<double> a, int rows, PackedColumns b, Span<double> c)
    {
        ArgumentNullException.ThrowIfNull(b);
        if (rows < 0) throw new ArgumentOutOfRangeException(nameof(rows));
        int inner = b.Inner, columns = b.Columns;
        if (a.Length < (long)rows * inner) throw new ArgumentException("A must hold rows × inner entries.", nameof(a));
        if (c.Length < (long)rows * columns) throw new ArgumentException("C must hold rows × columns entries.", nameof(c));
        if (rows == 0 || columns == 0) return;
        Span<double> into = c[..(rows * columns)];
        if (inner == 0)
        {
            into.Clear();
            return;
        }
        int panels = b.Panels, stride = columns;
        double[]? padded = null;
        if (columns % Width != 0)
        {
            // The last panel's tile is written whole, so C is taken padded to whole panels and copied out at the end.
            stride = panels * Width;
            padded = new double[rows * stride];
        }
        ref double ar = ref MemoryMarshal.GetReference(a);
        ref double cr = ref (padded is null ? ref MemoryMarshal.GetReference(into) : ref padded[0]);
        ref double br = ref b.Data[0];
        for (int k0 = 0; k0 < inner; k0 += Chunk)
        {
            int length = Math.Min(Chunk, inner - k0);
            bool first = k0 == 0;
            int i = 0;
            for (; i + Rows <= rows; i += Rows)
            {
                ref double ai = ref Unsafe.Add(ref ar, (nint)i * inner + k0);
                for (int p = 0; p < panels; p++)
                    Tile4(ref ai, inner, ref Unsafe.Add(ref br, ((nint)p * inner + k0) * Width), length,
                          ref Unsafe.Add(ref cr, (nint)i * stride + p * Width), stride, first);
            }
            for (; i < rows; i++)
            {
                ref double ai = ref Unsafe.Add(ref ar, (nint)i * inner + k0);
                for (int p = 0; p < panels; p++)
                    Tile1(ref ai, ref Unsafe.Add(ref br, ((nint)p * inner + k0) * Width), length,
                          ref Unsafe.Add(ref cr, (nint)i * stride + p * Width), first);
            }
        }
        if (padded is not null)
            for (int i = 0; i < rows; i++) padded.AsSpan(i * stride, columns).CopyTo(into.Slice(i * columns, columns));
    }

    /// <summary>C = A·B of two row-major matrices.</summary>
    public static double[,] Multiply(double[,] a, double[,] b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        int rows = a.GetLength(0), inner = a.GetLength(1), columns = b.GetLength(1);
        if (b.GetLength(0) != inner) throw new ArgumentException("B must have as many rows as A has columns.", nameof(b));
        var c = new double[rows, columns];
        if (rows == 0 || columns == 0) return c;
        ReadOnlySpan<double> left = inner == 0 ? default : MemoryMarshal.CreateReadOnlySpan(ref a[0, 0], a.Length);
        ReadOnlySpan<double> right = inner == 0 ? default : MemoryMarshal.CreateReadOnlySpan(ref b[0, 0], b.Length);
        Multiply(left, rows, Pack(right, inner, columns), MemoryMarshal.CreateSpan(ref c[0, 0], c.Length));
        return c;
    }

    /// <summary>Σ_k a_k·b_k by fused multiply-adds in order from +0: one entry of <see cref="Multiply(ReadOnlySpan{double}, int, PackedColumns, Span{double})"/>, bit for bit.</summary>
    public static double FusedDot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Vector dimensions must match.", nameof(b));
        double sum = 0.0;
        for (int k = 0; k < a.Length; k++) sum = Math.FusedMultiplyAdd(a[k], b[k], sum);
        return sum;
    }

    // Four rows of A against one panel over one chunk of the inner index: eight accumulators, two per row.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static void Tile4(ref double a, int inner, ref double b, int length, ref double c, int stride, bool first)
    {
        ref double a1 = ref Unsafe.Add(ref a, inner);
        ref double a2 = ref Unsafe.Add(ref a, 2 * inner);
        ref double a3 = ref Unsafe.Add(ref a, 3 * inner);
        ref double c1 = ref Unsafe.Add(ref c, stride);
        ref double c2 = ref Unsafe.Add(ref c, 2 * stride);
        ref double c3 = ref Unsafe.Add(ref c, 3 * stride);
        Vector256<double> c00, c01, c10, c11, c20, c21, c30, c31;
        if (first)
        {
            c00 = c01 = c10 = c11 = c20 = c21 = c30 = c31 = Vector256<double>.Zero;
        }
        else
        {
            c00 = Vector256.LoadUnsafe(ref c); c01 = Vector256.LoadUnsafe(ref c, 4);
            c10 = Vector256.LoadUnsafe(ref c1); c11 = Vector256.LoadUnsafe(ref c1, 4);
            c20 = Vector256.LoadUnsafe(ref c2); c21 = Vector256.LoadUnsafe(ref c2, 4);
            c30 = Vector256.LoadUnsafe(ref c3); c31 = Vector256.LoadUnsafe(ref c3, 4);
        }
        for (int k = 0; k < length; k++)
        {
            Vector256<double> b0 = Vector256.LoadUnsafe(ref b, (nuint)(k * Width));
            Vector256<double> b1 = Vector256.LoadUnsafe(ref b, (nuint)(k * Width + 4));
            Vector256<double> x = Vector256.Create(Unsafe.Add(ref a, k));
            c00 = Vector256.FusedMultiplyAdd(x, b0, c00);
            c01 = Vector256.FusedMultiplyAdd(x, b1, c01);
            x = Vector256.Create(Unsafe.Add(ref a1, k));
            c10 = Vector256.FusedMultiplyAdd(x, b0, c10);
            c11 = Vector256.FusedMultiplyAdd(x, b1, c11);
            x = Vector256.Create(Unsafe.Add(ref a2, k));
            c20 = Vector256.FusedMultiplyAdd(x, b0, c20);
            c21 = Vector256.FusedMultiplyAdd(x, b1, c21);
            x = Vector256.Create(Unsafe.Add(ref a3, k));
            c30 = Vector256.FusedMultiplyAdd(x, b0, c30);
            c31 = Vector256.FusedMultiplyAdd(x, b1, c31);
        }
        c00.StoreUnsafe(ref c); c01.StoreUnsafe(ref c, 4);
        c10.StoreUnsafe(ref c1); c11.StoreUnsafe(ref c1, 4);
        c20.StoreUnsafe(ref c2); c21.StoreUnsafe(ref c2, 4);
        c30.StoreUnsafe(ref c3); c31.StoreUnsafe(ref c3, 4);
    }

    // One row of A against one panel over one chunk: the rows left after the tiles of four.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static void Tile1(ref double a, ref double b, int length, ref double c, bool first)
    {
        Vector256<double> c0 = first ? Vector256<double>.Zero : Vector256.LoadUnsafe(ref c);
        Vector256<double> c1 = first ? Vector256<double>.Zero : Vector256.LoadUnsafe(ref c, 4);
        for (int k = 0; k < length; k++)
        {
            Vector256<double> x = Vector256.Create(Unsafe.Add(ref a, k));
            c0 = Vector256.FusedMultiplyAdd(x, Vector256.LoadUnsafe(ref b, (nuint)(k * Width)), c0);
            c1 = Vector256.FusedMultiplyAdd(x, Vector256.LoadUnsafe(ref b, (nuint)(k * Width + 4)), c1);
        }
        c0.StoreUnsafe(ref c);
        c1.StoreUnsafe(ref c, 4);
    }
}

/// <summary>
/// The right operand of <see cref="DenseProduct"/>, inner × columns, in panels of <see cref="DenseProduct.Width"/> columns:
/// panel p holds columns 8p … 8p + 7 for every inner index in turn, the last panel padded with zeros. Refilled in place
/// (<see cref="Fill"/>), it serves a sequence of products of one shape without new storage.
/// </summary>
public sealed class PackedColumns
{
    /// <summary>Storage for an inner × columns operand, filled with zeros.</summary>
    public PackedColumns(int inner, int columns)
    {
        if (inner < 0) throw new ArgumentOutOfRangeException(nameof(inner));
        if (columns < 0) throw new ArgumentOutOfRangeException(nameof(columns));
        Inner = inner;
        Columns = columns;
        Panels = (columns + DenseProduct.Width - 1) / DenseProduct.Width;
        Data = new double[Math.Max(1, (long)Panels * inner * DenseProduct.Width)];
    }

    /// <summary>The inner dimension: B's rows.</summary>
    public int Inner { get; }

    /// <summary>B's columns.</summary>
    public int Columns { get; }

    /// <summary>The panels, ⌈columns / 8⌉.</summary>
    public int Panels { get; }

    internal double[] Data { get; }

    /// <summary>Packs <paramref name="b"/>, inner × columns row-major, over the current contents.</summary>
    public void Fill(ReadOnlySpan<double> b)
    {
        if (b.Length < (long)Inner * Columns) throw new ArgumentException("B must hold inner × columns entries.", nameof(b));
        const int width = DenseProduct.Width;
        for (int p = 0; p < Panels; p++)
        {
            int j0 = p * width, count = Math.Min(width, Columns - j0);
            Span<double> panel = Data.AsSpan(p * Inner * width, Inner * width);
            for (int k = 0; k < Inner; k++)
            {
                b.Slice(k * Columns + j0, count).CopyTo(panel.Slice(k * width, count));
                if (count < width) panel.Slice(k * width + count, width - count).Clear();
            }
        }
    }
}
