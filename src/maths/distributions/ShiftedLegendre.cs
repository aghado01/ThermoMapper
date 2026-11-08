using System;
using System.Buffers;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace Maths.Distributions;

/// <summary>
/// The standardized shifted Legendre polynomials ψ_α(x) = √(2α + 1)·P_α(2x − 1), orthonormal on [0, 1] under the uniform
/// measure: KHAOS's univariate basis (E19). Values come from Bonnet's three-term recurrence, written
/// P_{k+1}(u) = aₖ·u·P_k(u) − bₖ·P_{k−1}(u) with aₖ = (2k + 1)/(k + 1) and bₖ = k/(k + 1) tabled once, so a step multiplies and
/// never divides; slopes come from P′_{k+1} = P′_{k−1} + (2k + 1)P_k, which holds at u = ±1 where the closed form
/// k(uP_k − P_{k−1})/(u² − 1) divides by zero. Both are O(α) per point and keep their digits to degree 30 and beyond on [0, 1]. A
/// column's values come in lanes, four points at a time, with each lane's operations the scalar recurrence's in its order, so the
/// bits match the scalar form on any hardware.
/// </summary>
public static class ShiftedLegendre
{
    private const int Tabled = 512;
    private static readonly double[] A = Table(k => (2.0 * k + 1.0) / (k + 1.0)), B = Table(k => k / (k + 1.0));

    private static double[] Table(Func<int, double> f)
    {
        var t = new double[Tabled];
        for (int k = 0; k < Tabled; k++) t[k] = f(k);
        return t;
    }

    // aₖ and bₖ, from the table below its end.
    private static double StepA(int k) => k < Tabled ? A[k] : (2.0 * k + 1.0) / (k + 1.0);

    private static double StepB(int k) => k < Tabled ? B[k] : k / (k + 1.0);

    /// <summary>ψ_<paramref name="degree"/>(<paramref name="x"/>).</summary>
    public static double Value(int degree, double x)
    {
        if (degree < 0) throw new ArgumentOutOfRangeException(nameof(degree), "A degree is not negative.");
        return Math.Sqrt(2.0 * degree + 1.0) * Polynomial(degree, 2.0 * x - 1.0);
    }

    /// <summary>ψ′_<paramref name="degree"/>(<paramref name="x"/>) = 2√(2α + 1)·P′_α(2x − 1).</summary>
    public static double Derivative(int degree, double x)
    {
        if (degree < 0) throw new ArgumentOutOfRangeException(nameof(degree), "A degree is not negative.");
        if (degree == 0) return 0.0;
        double u = 2.0 * x - 1.0;
        double p0 = 1.0, p1 = u;        // P_{k−1}, P_k
        double d0 = 0.0, d1 = 1.0;      // P′_{k−1}, P′_k
        for (int k = 1; k < degree; k++)
        {
            double p2 = StepA(k) * u * p1 - StepB(k) * p0;
            double d2 = d0 + (2 * k + 1) * p1;
            p0 = p1; p1 = p2;
            d0 = d1; d1 = d2;
        }
        return 2.0 * Math.Sqrt(2.0 * degree + 1.0) * d1;
    }

    /// <summary>Fills <paramref name="values"/>[α] with ψ_α(<paramref name="x"/>) for α below its length.</summary>
    public static void Values(double x, Span<double> values)
    {
        if (values.Length == 0) return;
        double u = 2.0 * x - 1.0;
        double p0 = 1.0, p1 = u;
        values[0] = 1.0;
        if (values.Length > 1) values[1] = Math.Sqrt(3.0) * u;
        for (int k = 1; k + 1 < values.Length; k++)
        {
            double p2 = StepA(k) * u * p1 - StepB(k) * p0;
            p0 = p1;
            p1 = p2;
            values[k + 1] = Math.Sqrt(2.0 * (k + 1) + 1.0) * p2;
        }
    }

    /// <summary>
    /// Multiplies <paramref name="column"/>[i] by ψ_<paramref name="degree"/>(<paramref name="x"/>[i]). The recurrence runs one step
    /// at a time over the whole column, four lanes at once, multiplications and subtractions only, so the steps of different points
    /// overlap where one point's chain of α steps would wait on each; then the column is scaled by √(2α + 1). Each element sees
    /// the same operations in the same order as <see cref="Value"/>, so the bits are its own.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static void Multiply(int degree, ReadOnlySpan<double> x, Span<double> column)
    {
        if (degree < 0) throw new ArgumentOutOfRangeException(nameof(degree), "A degree is not negative.");
        int n = column.Length;
        if (x.Length != n) throw new ArgumentException("One value per row.", nameof(x));
        if (degree == 0 || n == 0) return;
        double root = Math.Sqrt(2.0 * degree + 1.0);
        double[] scratch = ArrayPool<double>.Shared.Rent(3 * n);
        try
        {
            Span<double> u = scratch.AsSpan(0, n), p0 = scratch.AsSpan(n, n), p1 = scratch.AsSpan(2 * n, n);
            for (int i = 0; i < n; i++)
            {
                u[i] = 2.0 * x[i] - 1.0;
                p0[i] = 1.0;
                p1[i] = u[i];
            }
            for (int k = 1; k < degree; k++)
            {
                Step(StepA(k), StepB(k), u, p0, p1);   // P_{k+1} into p0
                Span<double> t = p0;
                p0 = p1;
                p1 = t;
            }
            int j = 0;
            if (n >= Vector256<double>.Count)
            {
                ref double c = ref MemoryMarshal.GetReference(column);
                ref double pr = ref MemoryMarshal.GetReference(p1);
                Vector256<double> vroot = Vector256.Create(root);
                for (; j + Vector256<double>.Count <= n; j += Vector256<double>.Count)
                    Vector256.Multiply(Vector256.LoadUnsafe(ref c, (nuint)j), Vector256.Multiply(vroot, Vector256.LoadUnsafe(ref pr, (nuint)j))).StoreUnsafe(ref c, (nuint)j);
            }
            for (; j < n; j++) column[j] *= root * p1[j];
        }
        finally
        {
            ArrayPool<double>.Shared.Return(scratch);
        }
    }

    // previous[i] ← a·u[i]·current[i] − b·previous[i], the scalar form's operations in its order, four lanes at a time.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static void Step(double a, double b, ReadOnlySpan<double> u, Span<double> previous, ReadOnlySpan<double> current)
    {
        int n = u.Length, i = 0;
        if (n >= Vector256<double>.Count)
        {
            ref double ur = ref MemoryMarshal.GetReference(u);
            ref double pr = ref MemoryMarshal.GetReference(previous);
            ref double cr = ref MemoryMarshal.GetReference(current);
            Vector256<double> va = Vector256.Create(a), vb = Vector256.Create(b);
            for (; i + Vector256<double>.Count <= n; i += Vector256<double>.Count)
            {
                Vector256<double> next = Vector256.Subtract(Vector256.Multiply(Vector256.Multiply(va, Vector256.LoadUnsafe(ref ur, (nuint)i)), Vector256.LoadUnsafe(ref cr, (nuint)i)),
                                                            Vector256.Multiply(vb, Vector256.LoadUnsafe(ref pr, (nuint)i)));
                next.StoreUnsafe(ref pr, (nuint)i);
            }
        }
        for (; i < n; i++) previous[i] = a * u[i] * current[i] - b * previous[i];
    }

    /// <summary>
    /// Writes ψ_<paramref name="degree"/>'s coefficients about <paramref name="centre"/>, cₖ of (x − c)ᵏ, into
    /// <paramref name="coefficients"/> (one per power up to the degree, the rest zeroed): the recurrence run on coefficient
    /// arrays in y = x − c, with u = 2y + (2c − 1). The expansion a one-piece polynomial hands the moment service
    /// (<c>PiecewisePolynomial</c>). Its coefficients grow with the degree, about 2^(2α)/√α at the centre ½, so the moments
    /// it serves lose digits past degree 15 or so; the uniform measure's are exact by orthonormality without it.
    /// </summary>
    public static void Coefficients(int degree, double centre, Span<double> coefficients)
    {
        if (degree < 0) throw new ArgumentOutOfRangeException(nameof(degree), "A degree is not negative.");
        if (coefficients.Length < degree + 1) throw new ArgumentException("One coefficient per power up to the degree.", nameof(coefficients));
        coefficients.Clear();
        double shift = 2.0 * centre - 1.0;
        Span<double> p0 = degree < 128 ? stackalloc double[degree + 1] : new double[degree + 1];
        Span<double> p1 = degree < 128 ? stackalloc double[degree + 1] : new double[degree + 1];
        Span<double> p2 = degree < 128 ? stackalloc double[degree + 1] : new double[degree + 1];
        p0[0] = 1.0;
        if (degree == 0)
        {
            coefficients[0] = 1.0;
            return;
        }
        p1[0] = shift;
        p1[1] = 2.0;
        for (int k = 1; k < degree; k++)
        {
            // (k + 1)P_{k+1} = (2k + 1)(2y + s)P_k − kP_{k−1}, on coefficients in y.
            p2.Clear();
            for (int r = 0; r <= k; r++)
            {
                p2[r] += (2 * k + 1) * shift * p1[r];
                p2[r + 1] += (2 * k + 1) * 2.0 * p1[r];
            }
            for (int r = 0; r < k; r++) p2[r] -= k * p0[r];
            for (int r = 0; r <= k + 1; r++) p2[r] /= k + 1;
            Span<double> t = p0;
            p0 = p1;
            p1 = p2;
            p2 = t;
        }
        double root = Math.Sqrt(2.0 * degree + 1.0);
        for (int r = 0; r <= degree; r++) coefficients[r] = root * p1[r];
    }

    // P_α(u) by the recurrence, the scalar form every lane matches.
    private static double Polynomial(int degree, double u)
    {
        if (degree == 0) return 1.0;
        double p0 = 1.0, p1 = u;
        for (int k = 1; k < degree; k++)
        {
            double p2 = StepA(k) * u * p1 - StepB(k) * p0;
            p0 = p1;
            p1 = p2;
        }
        return p1;
    }
}
