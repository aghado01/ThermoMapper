using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Maths.Regression.Tests;

/// <summary>
/// A reference for the sampled-sensitivity brief (slice 0): a function from outside BASS's model with its ANOVA decomposition in
/// closed form under the product uniform law on its box. Each family gives the mean, the closed variance τ_u = Var(E[f | x_u]) of
/// every subset and the conditional mean E[f | x_u]; the main variance of u (its ANOVA component's), each total and each effect follow
/// by Möbius inversion, τ_∅ = 0, unless the family states them directly to keep their digits.
/// </summary>
internal abstract class SobolReference(string name, ReferenceBox law)
{
    public string Name { get; } = name;

    public ReferenceBox Law { get; } = law;

    public int Dimension => Law.Dimension;

    public abstract double Value(ReadOnlySpan<double> x);

    /// <summary>E[f].</summary>
    public abstract double Mean { get; }

    /// <summary>τ_u = Var(E[f | x_u]) for u nonempty and ascending.</summary>
    public abstract double ClosedVariance(int[] u);

    /// <summary>E[f | x_u = <paramref name="values"/>], the values in u's order.</summary>
    public abstract double ConditionalMean(int[] u, ReadOnlySpan<double> values);

    /// <summary>Whether f varies with input <paramref name="v"/>; quadrature gives an inert input one node.</summary>
    public virtual bool Varies(int v) => true;

    /// <summary>The nodes a quadrature of f needs per panel on input <paramref name="v"/>, and the panels.</summary>
    public virtual (int Nodes, int Panels) Quadrature(int v) => Varies(v) ? (24, 2) : (1, 1);

    public double Variance => ClosedVariance([.. Enumerable.Range(0, Dimension)]);

    /// <summary>The variance of f's ANOVA component on u: Σ_{v⊆u} (−1)^{|u|−|v|} τ_v.</summary>
    public virtual double MainVariance(int[] u) => Mobius(u, v => v.Length == 0 ? 0.0 : ClosedVariance(v));

    /// <summary>Input i's total variance, V − τ_{−i}.</summary>
    public virtual double TotalVariance(int i) =>
        Dimension == 1 ? Variance : Variance - ClosedVariance([.. Enumerable.Range(0, Dimension).Where(v => v != i)]);

    /// <summary>The ANOVA component on u at <paramref name="values"/>: Σ_{v⊆u} (−1)^{|u|−|v|} E[f | x_v], E[f] at v = ∅.</summary>
    public double Effect(int[] u, ReadOnlySpan<double> values)
    {
        double[] at = values.ToArray();
        return Mobius(u, v => v.Length == 0 ? Mean : ConditionalMean(v, [.. v.Select(i => at[Array.IndexOf(u, i)])]));
    }

    public override string ToString() => Name;

    protected static double Mobius(int[] u, Func<int[], double> closed)
    {
        double s = 0.0;
        for (int mask = 0; mask < 1 << u.Length; mask++)
        {
            int[] v = [.. Enumerable.Range(0, u.Length).Where(k => (mask >> k & 1) != 0).Select(k => u[k])];
            s += ((u.Length - v.Length) % 2 == 0 ? 1.0 : -1.0) * closed(v);
        }
        return s;
    }
}

/// <summary>
/// Ishigami's function sin x₁ + a sin² x₂ + b x₃⁴ sin x₁ on [−π, π]³ (Ishigami and Homma 1990, cited from memory). With
/// E sin = 0, E sin² = ½, E sin⁴ = 3/8, E x₃⁴ = π⁴/5 and E x₃⁸ = π⁸/9: E[f] = a/2, τ₂ = a²/8, τ₁ = ½(1 + bπ⁴/5)²,
/// τ₁₃ = ½(1 + 2bπ⁴/5 + b²π⁸/9), τ₃ = 0, so V₁₃ = 8b²π⁸/225 and x₃ has no main effect but a total one.
/// </summary>
internal sealed class IshigamiReference(double a, double b) : SobolReference($"Ishigami (a = {a}, b = {b})", ReferenceBox.Cube(3, -Math.PI, Math.PI))
{
    private static readonly double Pi4 = Math.Pow(Math.PI, 4), Pi8 = Math.Pow(Math.PI, 8);

    public double A { get; } = a;

    public double B { get; } = b;

    public override double Value(ReadOnlySpan<double> x)
    {
        double s = Math.Sin(x[0]), x3 = x[2] * x[2];
        return s + A * Math.Sin(x[1]) * Math.Sin(x[1]) + B * x3 * x3 * s;
    }

    public override double Mean => 0.5 * A;

    public override double ClosedVariance(int[] u)
    {
        double tau = u.Contains(1) ? A * A / 8.0 : 0.0;
        if (u.Contains(0))
            tau += u.Contains(2) ? 0.5 * (1.0 + 2.0 * B * Pi4 / 5.0 + B * B * Pi8 / 9.0) : 0.5 * Math.Pow(1.0 + B * Pi4 / 5.0, 2);
        return tau;
    }

    public override double MainVariance(int[] u) => u switch
    {
        [0] => 0.5 * Math.Pow(1.0 + B * Pi4 / 5.0, 2),
        [1] => A * A / 8.0,
        [0, 2] => 8.0 * B * B * Pi8 / 225.0,
        _ => 0.0,
    };

    public override double TotalVariance(int i) => i switch
    {
        0 => MainVariance([0]) + MainVariance([0, 2]),
        1 => MainVariance([1]),
        _ => MainVariance([0, 2]),
    };

    public override double ConditionalMean(int[] u, ReadOnlySpan<double> values)
    {
        double[] given = values.ToArray();
        double At(int v) => given[Array.IndexOf(u, v)];
        double second = u.Contains(1) ? A * Math.Pow(Math.Sin(At(1)), 2) : 0.5 * A;
        if (!u.Contains(0)) return second;
        double x3 = u.Contains(2) ? Math.Pow(At(2), 4) : Pi4 / 5.0;
        return Math.Sin(At(0)) * (1.0 + B * x3) + second;
    }

    public override bool Varies(int v) => v != 2 || B != 0.0;
}

/// <summary>
/// Sobol's g-function Π_i (|4x_i − 2| + a_i)/(1 + a_i) on [0, 1]^p (Saltelli and Sobol' 1995, cited from memory). Each factor has
/// mean 1 and variance v_i = 1/(3(1 + a_i)²), so E[f | x_u] = Π_{i∈u} g_i(x_i), τ_u = Π_{i∈u}(1 + v_i) − 1, the main variance of u
/// is Π_{i∈u} v_i and the total of i is v_i Π_{j≠i}(1 + v_j): a large a_i makes a small index.
/// </summary>
internal sealed class SobolGReference(params double[] a) : SobolReference($"g-function (a = {string.Join(", ", a)})", ReferenceBox.Cube(a.Length, 0.0, 1.0))
{
    public IReadOnlyList<double> A { get; } = a;

    private double Part(int i) => 1.0 / (3.0 * (1.0 + A[i]) * (1.0 + A[i]));

    public override double Value(ReadOnlySpan<double> x)
    {
        double f = 1.0;
        for (int i = 0; i < x.Length; i++) f *= (Math.Abs(4.0 * x[i] - 2.0) + A[i]) / (1.0 + A[i]);
        return f;
    }

    public override double Mean => 1.0;

    // Π(1 + v_i) − 1 by t ← t + v_i(1 + t), which adds positive terms only.
    public override double ClosedVariance(int[] u) => u.Aggregate(0.0, (t, i) => t + Part(i) * (1.0 + t));

    public override double MainVariance(int[] u) => u.Aggregate(1.0, (s, i) => s * Part(i));

    public override double TotalVariance(int i) => Part(i) * Enumerable.Range(0, Dimension).Where(j => j != i).Aggregate(1.0, (s, j) => s * (1.0 + Part(j)));

    public override double ConditionalMean(int[] u, ReadOnlySpan<double> values)
    {
        double f = 1.0;
        for (int k = 0; k < u.Length; k++) f *= (Math.Abs(4.0 * values[k] - 2.0) + A[u[k]]) / (1.0 + A[u[k]]);
        return f;
    }

    // The kink at ½: two panels make each piece linear, and two nodes integrate the squares exactly.
    public override (int Nodes, int Panels) Quadrature(int v) => (2, 2);
}

/// <summary>
/// Friedman's function 10 sin(πx₁x₂) + 20(x₃ − ½)² + 10x₄ + 5x₅ on [0, 1]^p, p ≥ 5 (<see cref="FriedmanReference"/>), its inputs past
/// the fifth inert. It is additive across {x₁, x₂}, x₃, x₄ and x₅, so τ_u sums each group's share. With g = 10 sin(πx₁x₂):
/// E[g | x₁] = 10(1 − cos πx₁)/(πx₁), E[g] = 10·Cin(π)/π, E[g²] = 50(1 − Si(2π)/(2π)), and E[E[g | x₁]²] = (100/π²)(2π Si(π) − π Si(2π) − 4)
/// by parts; τ₃ = 20/9, τ₄ = 25/3 and τ₅ = 25/12.
/// </summary>
internal sealed class FriedmanSobolReference : SobolReference
{
    private static readonly double SiPi = SineIntegral(PiDecimal), SiTwoPi = SineIntegral(2m * PiDecimal), CinPi = CosineIntegral(PiDecimal);
    private const decimal PiDecimal = 3.1415926535897932384626433833m;
    private readonly FriedmanReference _f;

    public FriedmanSobolReference(int dimension) : base("Friedman's function", ReferenceBox.Cube(dimension, 0.0, 1.0)) => _f = new FriedmanReference(dimension);

    public static double GMean => 10.0 * CinPi / Math.PI;

    public static double GVariance => 50.0 * (1.0 - SiTwoPi / (2.0 * Math.PI)) - GMean * GMean;

    public static double GOneVariance => 100.0 / (Math.PI * Math.PI) * (2.0 * Math.PI * SiPi - Math.PI * SiTwoPi - 4.0) - GMean * GMean;

    public override double Value(ReadOnlySpan<double> x) => _f.Value(x);

    public override double Mean => GMean + 20.0 / 12.0 + 5.0 + 2.5;

    public override double ClosedVariance(int[] u)
    {
        bool one = u.Contains(0), two = u.Contains(1);
        double tau = one && two ? GVariance : one || two ? GOneVariance : 0.0;
        if (u.Contains(2)) tau += 20.0 / 9.0;
        if (u.Contains(3)) tau += 25.0 / 3.0;
        if (u.Contains(4)) tau += 25.0 / 12.0;
        return tau;
    }

    public override double ConditionalMean(int[] u, ReadOnlySpan<double> values)
    {
        double[] given = values.ToArray();
        double At(int v) => given[Array.IndexOf(u, v)];
        bool one = u.Contains(0), two = u.Contains(1);
        double g = one && two ? 10.0 * Math.Sin(Math.PI * At(0) * At(1)) : one ? GivenOne(At(0)) : two ? GivenOne(At(1)) : GMean;
        double x3 = u.Contains(2) ? 20.0 * (At(2) - 0.5) * (At(2) - 0.5) : 20.0 / 12.0;
        double x4 = u.Contains(3) ? 10.0 * At(3) : 5.0;
        double x5 = u.Contains(4) ? 5.0 * At(4) : 2.5;
        return g + x3 + x4 + x5;
    }

    public override bool Varies(int v) => v < 5;

    // x₃'s square needs three nodes; x₄ and x₅ are linear.
    public override (int Nodes, int Panels) Quadrature(int v) => v switch { 0 or 1 => (24, 2), 2 => (3, 1), 3 or 4 => (2, 1), _ => (1, 1) };

    /// <summary>E[g | x₁ = x] = 10(1 − cos πx)/(πx) = 20 sin²(πx/2)/(πx), 0 at x = 0.</summary>
    public static double GivenOne(double x)
    {
        double z = Math.PI * x;
        if (Math.Abs(z) < 1e-4) return 10.0 * (z / 2.0 - z * z * z / 24.0);
        double s = Math.Sin(z / 2.0);
        return 20.0 * s * s / z;
    }

    /// <summary>Si(z) = Σ_k (−1)^k z^(2k+1)/((2k + 1)(2k + 1)!), summed in decimal.</summary>
    public static double SineIntegral(decimal z)
    {
        decimal z2 = z * z, term = z, sum = 0m;   // term = z^(2k+1)/(2k+1)!
        for (int k = 0; term / (2 * k + 1) > 1e-27m; k++)
        {
            sum += (k % 2 == 0 ? term : -term) / (2 * k + 1);
            term *= z2 / ((2 * k + 2) * (2 * k + 3));
        }
        return (double)sum;
    }

    /// <summary>Cin(z) = ∫₀^z (1 − cos t)/t dt = Σ_k (−1)^(k+1) z^(2k)/(2k·(2k)!), summed in decimal.</summary>
    public static double CosineIntegral(decimal z)
    {
        decimal z2 = z * z, term = z2 / 2m, sum = 0m;   // term = z^(2k)/(2k)!
        for (int k = 1; term / (2 * k) > 1e-27m; k++)
        {
            sum += (k % 2 == 1 ? term : -term) / (2 * k);
            term *= z2 / ((2 * k + 1) * (2 * k + 2));
        }
        return (double)sum;
    }
}

/// <summary>
/// A sum of exponential ridges f = Σ_k Re(A_k e^{λ_k(w_kᵀx − c_k)}) (<see cref="RidgeSum"/> with <see cref="ExponentialProfile"/>s) under
/// the product uniform law on a box with centre m and half-widths h: a product of one-variable exponentials about m. With
/// δ_k = w_kᵀm − c_k and the centred moment E[e^{μ(x_j − m_j)}] = sinhc(μh_j), E[f | x_u] = Σ_k Re(α_k e^{λ_k w_{k,u}ᵀ(x_u − m_u)}),
/// α_k = A_k e^{λ_kδ_k} Π_{j∉u} sinhc(λ_k w_kj h_j), and E[E[f | x_u]²] = ½ Σ_{k,l} Re(α_kα_l Π_{j∈u} sinhc((λ_k w_kj + λ_l w_lj)h_j)
/// + α_kᾱ_l Π_{j∈u} sinhc((λ_k w_kj + λ̄_l w_lj)h_j)). Every exponent is the projection's own range, so a box far off the origin
/// keeps its digits; inputs no ridge reads leave τ unchanged, so they are dropped from u, and τ of them alone is 0.
/// </summary>
internal sealed class RidgeSobolReference : SobolReference
{
    private readonly RidgeSum _f;
    private readonly (Complex Amplitude, Complex Rate, double[] W, double Shift)[] _ridges;

    public RidgeSobolReference(string name, RidgeSum f, ReferenceBox law) : base(name, law)
    {
        if (f.Dimension != law.Dimension) throw new ArgumentException("The law has other inputs than the function.");
        if (!f.Ridges.All(r => r.Profile is ExponentialProfile)) throw new ArgumentException("Every ridge's profile is exponential.");
        _f = f;
        _ridges = [.. f.Ridges.Select(r =>
        {
            var e = (ExponentialProfile)r.Profile;
            return (e.Amplitude, e.Rate, r.Direction, r.Direction.Select((w, j) => w * law.Center[j]).Sum() - r.Offset);
        })];
    }

    /// <summary>E36's ridge references (<see cref="ProjectionPursuitReferences"/>) under their laws.</summary>
    public static IReadOnlyList<RidgeSobolReference> FromProjectionPursuit() =>
        [.. ProjectionPursuitReferences.All().Where(r => r.Function is RidgeSum).Select(r => new RidgeSobolReference(r.Name, (RidgeSum)r.Function, r.Law))];

    public override double Value(ReadOnlySpan<double> x) => _f.Value(x);

    public override double Mean => Conditional([], []);

    public override double ConditionalMean(int[] u, ReadOnlySpan<double> values) => Conditional(u, values);

    // τ_u = E[(E[f | x_u] − f₀)²], centred so a small τ keeps its digits: with h_k = α_k(e^{λ_k w_{k,u}ᵀ(x_u − m_u)} − Π_{j∈u} S_kj),
    // E[h_k h_l] = α_kα_l(Π_j S(a_j + b_j) − Π_j S(a_j)S(b_j)), a_j = λ_k w_kj, b_j = λ_l w_lj, S(μ) = sinhc(μh_j), the difference telescoped.
    public override double ClosedVariance(int[] u)
    {
        u = [.. u.Where(Varies)];
        if (u.Length == 0) return 0.0;
        Complex[] alpha = [.. Enumerable.Range(0, _ridges.Length).Select(k => Alpha(k, u))];
        double square = 0.0;
        for (int k = 0; k < _ridges.Length; k++)
            for (int l = 0; l < _ridges.Length; l++)
            {
                Complex rate = _ridges[l].Rate;
                Complex same = alpha[k] * alpha[l] * Telescoped(u, k, l, rate);
                Complex mixed = alpha[k] * Complex.Conjugate(alpha[l]) * Telescoped(u, k, l, Complex.Conjugate(rate));
                square += 0.5 * (same + mixed).Real;
            }
        return square;
    }

    // T_i = ½E[(f(x) − f(x'_i, x_{−i}))²] = ¼ Σ_{k,l} Re(c_kc_l Π_{j≠i} S(a_j + b_j)·2E_i(a_i, b_i) + the same at λ̄_l, c̄_l), c_k = A_k e^{λ_kδ_k}.
    public override double TotalVariance(int i)
    {
        if (!Varies(i)) return 0.0;
        int[] none = [.. Enumerable.Range(0, Dimension)];
        double total = 0.0;
        for (int k = 0; k < _ridges.Length; k++)
            for (int l = 0; l < _ridges.Length; l++)
            {
                Complex ck = Alpha(k, none), cl = Alpha(l, none);
                Complex Part(Complex rate, Complex c)
                {
                    Complex product = ck * c * 2.0 * Excess(i, _ridges[k].Rate * _ridges[k].W[i], rate * _ridges[l].W[i]);
                    for (int j = 0; j < Dimension; j++)
                        if (j != i) product *= Centred(j, _ridges[k].Rate * _ridges[k].W[j] + rate * _ridges[l].W[j]);
                    return product;
                }
                total += 0.25 * (Part(_ridges[l].Rate, cl) + Part(Complex.Conjugate(_ridges[l].Rate), Complex.Conjugate(cl))).Real;
            }
        return total;
    }

    // Π_{j∈u} S(a_j + b_j) − Π_{j∈u} S(a_j)S(b_j) = Σ_t Π_{s<t} S(a_s + b_s)·E_t·Π_{s>t} S(a_s)S(b_s), E_t = S(a_t + b_t) − S(a_t)S(b_t).
    private Complex Telescoped(int[] u, int k, int l, Complex rate)
    {
        Complex sum = Complex.Zero;
        for (int t = 0; t < u.Length; t++)
        {
            Complex term = Excess(u[t], _ridges[k].Rate * _ridges[k].W[u[t]], rate * _ridges[l].W[u[t]]);
            for (int s = 0; s < u.Length; s++)
            {
                if (s == t) continue;
                Complex a = _ridges[k].Rate * _ridges[k].W[u[s]], b = rate * _ridges[l].W[u[s]];
                term *= s < t ? Centred(u[s], a + b) : Centred(u[s], a) * Centred(u[s], b);
            }
            sum += term;
        }
        return sum;
    }

    // S(a + b) − S(a)S(b) on side j, through the parts less one: S − 1 = Σ_{n≥1} z^{2n}/(2n + 1)! where |z| < 1.
    private Complex Excess(int j, Complex a, Complex b)
    {
        double h = Law.HalfWidths[j];
        Complex ea = SinhcLessOne(a * h), eb = SinhcLessOne(b * h);
        return SinhcLessOne((a + b) * h) - ea - eb - ea * eb;
    }

    private static Complex SinhcLessOne(Complex z)
    {
        if (Complex.Abs(z) >= 1.0) return ReferenceBox.Sinhc(z) - 1.0;
        Complex z2 = z * z, term = Complex.One, sum = Complex.Zero;
        for (int n = 1; n <= 14; n++)
        {
            term *= z2 / ((2.0 * n) * (2.0 * n + 1.0));
            sum += term;
        }
        return sum;
    }

    public override bool Varies(int v) => _ridges.Any(r => r.W[v] != 0.0);

    // The profiles are entire and slow on these boxes: one panel of 16 nodes.
    public override (int Nodes, int Panels) Quadrature(int v) => Varies(v) ? (16, 1) : (1, 1);

    // α_k = A_k e^{λ_kδ_k} Π_{j∉u} sinhc(λ_k w_kj h_j).
    private Complex Alpha(int k, int[] u)
    {
        var (amplitude, rate, w, shift) = _ridges[k];
        Complex alpha = amplitude * Complex.Exp(rate * shift);
        for (int j = 0; j < Dimension; j++)
            if (Array.IndexOf(u, j) < 0) alpha *= Centred(j, rate * w[j]);
        return alpha;
    }

    private double Conditional(int[] u, ReadOnlySpan<double> values)
    {
        double f = 0.0;
        for (int k = 0; k < _ridges.Length; k++)
        {
            Complex s = Complex.Zero;
            for (int t = 0; t < u.Length; t++) s += _ridges[k].W[u[t]] * (values[t] - Law.Center[u[t]]);
            f += (Alpha(k, u) * Complex.Exp(_ridges[k].Rate * s)).Real;
        }
        return f;
    }

    // E[e^{μ(x_j − m_j)}] under the uniform law on side j.
    private Complex Centred(int j, Complex mu) => ReferenceBox.Sinhc(mu * Law.HalfWidths[j]);
}

/// <summary>A constant on a box: zero variance, every piece and effect 0.</summary>
internal sealed class ConstantReference(double value, int dimension) : SobolReference($"constant {value}", ReferenceBox.Cube(dimension, 0.0, 1.0))
{
    public override double Value(ReadOnlySpan<double> x) => value;

    public override double Mean => value;

    public override double ClosedVariance(int[] u) => 0.0;

    public override double ConditionalMean(int[] u, ReadOnlySpan<double> values) => value;

    public override bool Varies(int v) => false;
}

/// <summary>The slice's chosen references.</summary>
internal static class SobolReferences
{
    public static IshigamiReference Ishigami() => new(7.0, 0.1);

    /// <summary>Two large, two middling and two small indices: a = 99 makes S about 7·10⁻⁵ of the variance.</summary>
    public static SobolGReference G() => new(0.0, 1.0, 4.5, 9.0, 99.0, 99.0);

    public static IReadOnlyList<SobolReference> All() =>
        [Ishigami(), G(), new FriedmanSobolReference(6), .. RidgeSobolReference.FromProjectionPursuit()];
}

/// <summary>
/// The truth from a method apart from the closed forms: E[f | x_u] at a point, by tensor Gauss–Legendre quadrature over the other inputs,
/// and τ_u = E[(E[f | x_u] − E[f])²], by quadrature over u of those inner integrals, centred so a small τ keeps its digits.
/// </summary>
internal static class SobolQuadrature
{
    public static IReadOnlyList<QuadratureRule> Rules(SobolReference f) =>
        [.. Enumerable.Range(0, f.Dimension).Select(v => { var (nodes, panels) = f.Quadrature(v); return QuadratureRule.Uniform(f.Law.Lower[v], f.Law.Upper[v], nodes, panels); })];

    public static double Mean(SobolReference f, IReadOnlyList<QuadratureRule> rules) => ConditionalMean(f, rules, [], []);

    public static double ConditionalMean(SobolReference f, IReadOnlyList<QuadratureRule> rules, int[] u, ReadOnlySpan<double> values)
    {
        int p = f.Dimension;
        var x = new double[p];
        for (int t = 0; t < u.Length; t++) x[u[t]] = values[t];
        int[] rest = [.. Enumerable.Range(0, p).Where(v => Array.IndexOf(u, v) < 0)];
        double sum = 0.0;
        foreach (double w in Tensor(rules, rest, x)) sum += w * f.Value(x);
        return sum;
    }

    public static double ClosedVariance(SobolReference f, IReadOnlyList<QuadratureRule> rules, int[] u)
    {
        double mean = Mean(f, rules), sum = 0.0;
        var x = new double[f.Dimension];
        foreach (double w in Tensor(rules, u, x))
        {
            double[] at = [.. u.Select(v => x[v])];
            double c = ConditionalMean(f, rules, u, at) - mean;
            sum += w * c * c;
        }
        return sum;
    }

    // Walks the tensor grid over the inputs in `over`, writing each node into x and yielding its weight.
    private static IEnumerable<double> Tensor(IReadOnlyList<QuadratureRule> rules, int[] over, double[] x)
    {
        var index = new int[over.Length];
        while (true)
        {
            double w = 1.0;
            for (int t = 0; t < over.Length; t++)
            {
                x[over[t]] = rules[over[t]].Nodes[index[t]];
                w *= rules[over[t]].Weights[index[t]];
            }
            yield return w;
            int k = over.Length - 1;
            while (k >= 0 && ++index[k] == rules[over[k]].Nodes.Length) index[k--] = 0;
            if (k < 0) yield break;
        }
    }
}
