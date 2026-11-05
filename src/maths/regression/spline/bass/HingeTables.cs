using System;
using System.Linq;
using System.Runtime.CompilerServices;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// One continuous variable's cases in ascending order, and where each knot's cases start in it: the order a hinge's prefix table
/// is scanned in. Knots are the variable's distinct values (D01), so each case sits at exactly one. Data only: it reads no response,
/// so a design builds it once per variable.
/// </summary>
internal sealed class HingeOrder
{
    private HingeOrder(int[] order, int[] start, double[] knots)
    {
        Order = order;
        Start = start;
        Knots = knots;
    }

    /// <summary>The case indices in ascending order of the variable (ties in case order).</summary>
    public int[] Order { get; }

    /// <summary>The sorted position of knot k's first case; the last entry is the case count.</summary>
    public int[] Start { get; }

    /// <summary>The variable's knots, ascending.</summary>
    public double[] Knots { get; }

    /// <summary>The order of <paramref name="values"/> against <paramref name="knots"/>, or null when the knots are not exactly the
    /// values' distinct entries (a table would then index the wrong knot).</summary>
    public static HingeOrder? Of(ReadOnlySpan<double> values, ReadOnlySpan<double> knots)
    {
        int n = values.Length;
        double[] keys = values.ToArray();
        int[] order = new int[n];
        for (int i = 0; i < n; i++) order[i] = i;
        Array.Sort(keys, order);
        var start = new int[knots.Length + 1];
        int k = -1;
        for (int p = 0; p < n; p++)
        {
            if (p == 0 || keys[p] != keys[p - 1])
            {
                k++;
                if (k >= knots.Length || knots[k] != keys[p]) return null;
                start[k] = p;
            }
        }
        if (k != knots.Length - 1) return null;
        start[knots.Length] = n;
        return new HingeOrder(order, start, knots.ToArray());
    }
}

/// <summary>
/// The prefix table of a hinge, lever 1 of the efficiency pass: for a held column z and a variable, <c>T[k]</c> = Σ zᵢ (s(xᵢ − t_k))₊^q
/// at every knot k, one pass over the cases in the variable's order, so that the cross product of any one-way hinge at that variable with
/// z is a lookup. Each sum about a knot is formed from its neighbour's by the anchored update (A₀ ← A₀ + G, A₂ ← A₂ + δ(2A₁ + δA₀),
/// A₁ ← A₁ + δA₀, G the neighbour group's total of z, δ the gap), so every term is the size of the sum itself and no difference of raw
/// moments is taken. The result is a different order of the same sum than the dot's, so it differs from it by rounding, within
/// 4·K·ε·‖z‖‖u‖ (K the knots), not in bits. Degrees one and two.
/// </summary>
internal static class HingeScan
{
    /// <summary>T⁺[k] = Σ_{xᵢ &gt; t_k} gᵢ (xᵢ − t_k)^q, for g = z, or roots∘z when <paramref name="roots"/> is given.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double[] Above(ReadOnlySpan<double> z, double[]? roots, HingeOrder order, int degree)
    {
        int count = order.Knots.Length;
        var table = new double[count];
        double[] knots = order.Knots;
        int[] idx = order.Order, start = order.Start;
        double a0 = 0.0, a1 = 0.0, a2 = 0.0;
        for (int k = count - 2; k >= 0; k--)
        {
            double delta = knots[k + 1] - knots[k], group = 0.0;
            for (int p = start[k + 1]; p < start[k + 2]; p++) group += roots is null ? z[idx[p]] : z[idx[p]] * roots[idx[p]];
            a0 += group;
            if (degree == 1)
            {
                a1 += delta * a0;
                table[k] = a1;
            }
            else
            {
                a2 += delta * (2.0 * a1 + delta * a0);
                a1 += delta * a0;
                table[k] = a2;
            }
        }
        return table;
    }

    /// <summary>T⁻[k] = Σ_{xᵢ &lt; t_k} gᵢ (t_k − xᵢ)^q.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double[] Below(ReadOnlySpan<double> z, double[]? roots, HingeOrder order, int degree)
    {
        int count = order.Knots.Length;
        var table = new double[count];
        double[] knots = order.Knots;
        int[] idx = order.Order, start = order.Start;
        double b0 = 0.0, b1 = 0.0, b2 = 0.0;
        for (int k = 1; k < count; k++)
        {
            double delta = knots[k] - knots[k - 1], group = 0.0;
            for (int p = start[k - 1]; p < start[k]; p++) group += roots is null ? z[idx[p]] : z[idx[p]] * roots[idx[p]];
            b0 += group;
            if (degree == 1)
            {
                b1 += delta * b0;
                table[k] = b1;
            }
            else
            {
                b2 += delta * (2.0 * b1 + delta * b0);
                b1 += delta * b0;
                table[k] = b2;
            }
        }
        return table;
    }
}
