using System;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// The element-wise passes over a design's columns that a term's column evaluation and the commit's explicit residual make, four
/// lanes at a time. Every operation is a per-lane add, subtract, multiply, divide or compare in the order the scalar loop
/// makes them, with no fused multiply-add, so each result is the scalar loop's, bit for bit, on any hardware (the software
/// fallback of <c>Vector256</c> runs the same operations). The reductions (<see cref="NonzeroCount"/> aside, which is integer) stay in
/// <see cref="Maths.LinAlg.PairwiseSum"/>.
/// </summary>
internal static class ColumnOps
{
    private const int Lanes = 4;

    /// <summary>The entries of <paramref name="values"/> that are not zero (a NaN counts, as <c>v != 0.0</c> does).</summary>
    public static int NonzeroCount(ReadOnlySpan<double> values)
    {
        int n = values.Length, nonzero = 0, i = 0;
        if (n >= Lanes)
        {
            ref double v = ref MemoryMarshal.GetReference(values);
            for (; i + Lanes <= n; i += Lanes)
            {
                uint zeros = Vector256.Equals(Vector256.LoadUnsafe(ref v, (nuint)i), Vector256<double>.Zero).ExtractMostSignificantBits();
                nonzero += Lanes - System.Numerics.BitOperations.PopCount(zeros);
            }
        }
        for (; i < n; i++)
            if (values[i] != 0.0) nonzero++;
        return nonzero;
    }

    /// <summary><c>values[i] -= amount</c>.</summary>
    public static void Subtract(Span<double> values, double amount)
    {
        int n = values.Length, i = 0;
        if (n >= Lanes)
        {
            ref double v = ref MemoryMarshal.GetReference(values);
            Vector256<double> a = Vector256.Create(amount);
            for (; i + Lanes <= n; i += Lanes)
                Vector256.Subtract(Vector256.LoadUnsafe(ref v, (nuint)i), a).StoreUnsafe(ref v, (nuint)i);
        }
        for (; i < n; i++) values[i] -= amount;
    }

    /// <summary><c>values[i] = roots[i] * (values[i] - amount)</c>: the whitened, centred column.</summary>
    public static void CentreAndScale(Span<double> values, double amount, ReadOnlySpan<double> roots)
    {
        int n = values.Length, i = 0;
        if (roots.Length != n) throw new ArgumentException("One root per row.", nameof(roots));
        if (n >= Lanes)
        {
            ref double v = ref MemoryMarshal.GetReference(values);
            ref double r = ref MemoryMarshal.GetReference(roots);
            Vector256<double> a = Vector256.Create(amount);
            for (; i + Lanes <= n; i += Lanes)
                Vector256.Multiply(Vector256.LoadUnsafe(ref r, (nuint)i), Vector256.Subtract(Vector256.LoadUnsafe(ref v, (nuint)i), a)).StoreUnsafe(ref v, (nuint)i);
        }
        for (; i < n; i++) values[i] = roots[i] * (values[i] - amount);
    }

    /// <summary><c>into[i] = (values[i] - amount) * roots[i]</c>, or without the roots when they are null: the prior side's centred column.</summary>
    public static void CentreInto(Span<double> into, ReadOnlySpan<double> values, double amount, ReadOnlySpan<double> roots)
    {
        int n = values.Length, i = 0;
        if (into.Length != n) throw new ArgumentException("One entry per value.", nameof(into));
        bool scaled = !roots.IsEmpty;
        if (scaled && roots.Length != n) throw new ArgumentException("One root per row.", nameof(roots));
        if (n >= Lanes)
        {
            ref double d = ref MemoryMarshal.GetReference(into);
            ref double v = ref MemoryMarshal.GetReference(values);
            ref double r = ref MemoryMarshal.GetReference(roots);
            Vector256<double> a = Vector256.Create(amount);
            for (; i + Lanes <= n; i += Lanes)
            {
                Vector256<double> c = Vector256.Subtract(Vector256.LoadUnsafe(ref v, (nuint)i), a);
                if (scaled) c = Vector256.Multiply(c, Vector256.LoadUnsafe(ref r, (nuint)i));
                c.StoreUnsafe(ref d, (nuint)i);
            }
        }
        for (; i < n; i++)
        {
            double c = values[i] - amount;
            into[i] = scaled ? c * roots[i] : c;
        }
    }

    /// <summary><c>into[l] = factor * values[l]</c>: one row of an outer product.</summary>
    public static void Scale(Span<double> into, double factor, ReadOnlySpan<double> values)
    {
        int n = values.Length, i = 0;
        if (into.Length != n) throw new ArgumentException("One entry per value.", nameof(into));
        if (n >= Lanes)
        {
            ref double d = ref MemoryMarshal.GetReference(into);
            ref double v = ref MemoryMarshal.GetReference(values);
            Vector256<double> f = Vector256.Create(factor);
            for (; i + Lanes <= n; i += Lanes)
                Vector256.Multiply(f, Vector256.LoadUnsafe(ref v, (nuint)i)).StoreUnsafe(ref d, (nuint)i);
        }
        for (; i < n; i++) into[i] = factor * values[i];
    }

    /// <summary>
    /// <c>residual[i] -= betas[j] * columns[j][i]</c> for every column in order, the columns taken in blocks so that one pass over the
    /// residual subtracts a block: each entry still meets its columns one at a time, in column order, so it is the column-by-column loop's.
    /// </summary>
    public static void SubtractScaledColumns(Span<double> residual, ReadOnlySpan<double[]> columns, ReadOnlySpan<double> betas)
    {
        int count = columns.Length, j = 0;
        if (betas.Length != count) throw new ArgumentException("One coefficient per column.", nameof(betas));
        for (; j + 4 <= count; j += 4) Subtract4(residual, columns[j], columns[j + 1], columns[j + 2], columns[j + 3], betas[j], betas[j + 1], betas[j + 2], betas[j + 3]);
        for (; j < count; j++) Subtract1(residual, columns[j], betas[j]);
    }

    /// <summary><c>residual[i] -= beta * column[i]</c>.</summary>
    public static void SubtractScaled(Span<double> residual, ReadOnlySpan<double> column, double beta) => Subtract1(residual, column, beta);

    private static void Subtract1(Span<double> r, ReadOnlySpan<double> z0, double b0)
    {
        int n = r.Length, i = 0;
        if (z0.Length != n) throw new ArgumentException("One entry per row.");
        if (n >= Lanes)
        {
            ref double rr = ref MemoryMarshal.GetReference(r);
            ref double a0 = ref MemoryMarshal.GetReference(z0);
            Vector256<double> v0 = Vector256.Create(b0);
            for (; i + Lanes <= n; i += Lanes)
                Vector256.Subtract(Vector256.LoadUnsafe(ref rr, (nuint)i), Vector256.Multiply(v0, Vector256.LoadUnsafe(ref a0, (nuint)i))).StoreUnsafe(ref rr, (nuint)i);
        }
        for (; i < n; i++) r[i] -= b0 * z0[i];
    }

    private static void Subtract4(Span<double> r, double[] z0, double[] z1, double[] z2, double[] z3, double b0, double b1, double b2, double b3)
    {
        int n = r.Length, i = 0;
        if (z0.Length != n || z1.Length != n || z2.Length != n || z3.Length != n) throw new ArgumentException("One entry per row.");
        if (n >= Lanes)
        {
            ref double rr = ref MemoryMarshal.GetReference(r);
            ref double a0 = ref MemoryMarshal.GetArrayDataReference(z0);
            ref double a1 = ref MemoryMarshal.GetArrayDataReference(z1);
            ref double a2 = ref MemoryMarshal.GetArrayDataReference(z2);
            ref double a3 = ref MemoryMarshal.GetArrayDataReference(z3);
            Vector256<double> v0 = Vector256.Create(b0), v1 = Vector256.Create(b1), v2 = Vector256.Create(b2), v3 = Vector256.Create(b3);
            for (; i + Lanes <= n; i += Lanes)
            {
                Vector256<double> v = Vector256.LoadUnsafe(ref rr, (nuint)i);
                v = Vector256.Subtract(v, Vector256.Multiply(v0, Vector256.LoadUnsafe(ref a0, (nuint)i)));
                v = Vector256.Subtract(v, Vector256.Multiply(v1, Vector256.LoadUnsafe(ref a1, (nuint)i)));
                v = Vector256.Subtract(v, Vector256.Multiply(v2, Vector256.LoadUnsafe(ref a2, (nuint)i)));
                v = Vector256.Subtract(v, Vector256.Multiply(v3, Vector256.LoadUnsafe(ref a3, (nuint)i)));
                v.StoreUnsafe(ref rr, (nuint)i);
            }
        }
        for (; i < n; i++)
        {
            double v = r[i];
            v -= b0 * z0[i];
            v -= b1 * z1[i];
            v -= b2 * z2[i];
            v -= b3 * z3[i];
            r[i] = v;
        }
    }

    /// <summary><c>residual[l] -= beta * (a * h[l] - level)</c>: a direct functional column's update over one case's locations.</summary>
    public static void SubtractScaledProduct(Span<double> residual, double beta, double a, ReadOnlySpan<double> h, double level)
    {
        int n = h.Length, i = 0;
        if (residual.Length != n) throw new ArgumentException("One entry per location.", nameof(residual));
        if (n >= Lanes)
        {
            ref double rr = ref MemoryMarshal.GetReference(residual);
            ref double hh = ref MemoryMarshal.GetReference(h);
            Vector256<double> vb = Vector256.Create(beta), va = Vector256.Create(a), vl = Vector256.Create(level);
            for (; i + Lanes <= n; i += Lanes)
            {
                Vector256<double> t = Vector256.Multiply(vb, Vector256.Subtract(Vector256.Multiply(va, Vector256.LoadUnsafe(ref hh, (nuint)i)), vl));
                Vector256.Subtract(Vector256.LoadUnsafe(ref rr, (nuint)i), t).StoreUnsafe(ref rr, (nuint)i);
            }
        }
        for (; i < n; i++) residual[i] -= beta * (a * h[i] - level);
    }
}
