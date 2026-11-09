using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Maths.Distributions;

namespace Maths.Regression.Tests;

/// <summary>
/// The law an E8 reference function's closed-form C = E[∇f ∇fᵀ] is under (the active-subspace brief, slice 0): its centre m and
/// its centred moment generating function K(μ) = E[e^{μᵀ(x − m)}] at a complex μ, which is all a sum of exponential ridges needs.
/// </summary>
internal abstract class ReferenceMeasure
{
    private protected ReferenceMeasure(double[] center) => Center = center;

    /// <summary>The centre m: the box's midpoint, the Gaussian's mean.</summary>
    public double[] Center { get; }

    /// <summary>The inputs.</summary>
    public int Dimension => Center.Length;

    /// <summary>K(μ) = E[e^{μᵀ(x − m)}].</summary>
    public abstract Complex CenteredMgf(ReadOnlySpan<Complex> mu);
}

/// <summary>
/// The product uniform law on the box Π [Lower_j, Upper_j], the law BASS's measures default to: K(μ) = Π sinh(μ_j h_j)/(μ_j h_j),
/// h_j the half-width, read through its series near 0 so that narrow boxes keep their digits.
/// </summary>
internal sealed class ReferenceBox : ReferenceMeasure
{
    public ReferenceBox(double[] lower, double[] upper) : base(Midpoints(lower, upper))
    {
        Lower = lower;
        Upper = upper;
        HalfWidths = [.. lower.Select((a, j) => 0.5 * (upper[j] - a))];
    }

    /// <summary>The box [<paramref name="lower"/>, <paramref name="upper"/>]^<paramref name="dimension"/>.</summary>
    public static ReferenceBox Cube(int dimension, double lower, double upper) =>
        new([.. Enumerable.Repeat(lower, dimension)], [.. Enumerable.Repeat(upper, dimension)]);

    public double[] Lower { get; }

    public double[] Upper { get; }

    public double[] HalfWidths { get; }

    /// <summary>Whether <paramref name="other"/> is the same box.</summary>
    public bool SameAs(ReferenceBox other) => Lower.SequenceEqual(other.Lower) && Upper.SequenceEqual(other.Upper);

    public override Complex CenteredMgf(ReadOnlySpan<Complex> mu)
    {
        Complex k = Complex.One;
        for (int j = 0; j < mu.Length; j++) k *= Sinhc(mu[j] * HalfWidths[j]);
        return k;
    }

    /// <summary>sinh(z)/z, by its series where |z| &lt; 1.</summary>
    public static Complex Sinhc(Complex z)
    {
        if (Complex.Abs(z) < 1.0)
        {
            Complex z2 = z * z, term = Complex.One, sum = Complex.One;
            for (int n = 1; n <= 12; n++)
            {
                term *= z2 / ((2.0 * n) * (2.0 * n + 1.0));
                sum += term;
            }
            return sum;
        }
        return new Complex(Math.Sinh(z.Real) * Math.Cos(z.Imaginary), Math.Cosh(z.Real) * Math.Sin(z.Imaginary)) / z;
    }

    private static double[] Midpoints(double[] lower, double[] upper)
    {
        ArgumentNullException.ThrowIfNull(lower);
        ArgumentNullException.ThrowIfNull(upper);
        if (lower.Length != upper.Length || lower.Length == 0) throw new ArgumentException("A box needs as many upper ends as lower ends, and an input.");
        for (int j = 0; j < lower.Length; j++)
            if (!double.IsFinite(lower[j]) || !double.IsFinite(upper[j]) || !(upper[j] > lower[j]))
                throw new ArgumentException("Each side of a box is finite, lower below upper.");
        return [.. lower.Select((a, j) => 0.5 * (a + upper[j]))];
    }
}

/// <summary>
/// The Gaussian N(m, Σ), a dependent law: K(μ) = exp(½ μᵀΣμ), the bilinear form at a complex μ. Its Cholesky factor L
/// (Σ = LLᵀ) carries standard normal z to x = m + Lz for quadrature.
/// </summary>
internal sealed class ReferenceGaussian : ReferenceMeasure
{
    public ReferenceGaussian(double[] mean, double[,] covariance) : base(mean)
    {
        int p = mean.Length;
        if (covariance.GetLength(0) != p || covariance.GetLength(1) != p) throw new ArgumentException("The covariance is p × p.");
        Covariance = covariance;
        Factor = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = covariance[i, j];
                for (int k = 0; k < j; k++) s -= Factor[i, k] * Factor[j, k];
                if (i == j)
                {
                    if (!(s > 0.0)) throw new ArgumentException("The covariance is positive definite.");
                    Factor[i, i] = Math.Sqrt(s);
                }
                else Factor[i, j] = s / Factor[j, j];
            }
    }

    public double[,] Covariance { get; }

    /// <summary>The lower-triangular L with Σ = LLᵀ.</summary>
    public double[,] Factor { get; }

    public override Complex CenteredMgf(ReadOnlySpan<Complex> mu)
    {
        Complex q = Complex.Zero;
        for (int i = 0; i < mu.Length; i++)
            for (int j = 0; j < mu.Length; j++) q += mu[i] * Covariance[i, j] * mu[j];
        return Complex.Exp(0.5 * q);
    }
}

/// <summary>
/// A reference function for E8 (the active-subspace brief, slice 0): a function from outside BASS's model, with its value, its
/// gradient in the inputs' own units, and its gradient outer-product matrix C = E[∇f ∇fᵀ] in closed form under the laws its
/// family has one for. <see cref="GradientQuadrature"/> reads C under any law by quadrature of the gradient.
/// </summary>
internal abstract class ReferenceFunction(string name, int dimension)
{
    public string Name { get; } = name;

    public int Dimension { get; } = dimension;

    public abstract double Value(ReadOnlySpan<double> x);

    public abstract void Gradient(ReadOnlySpan<double> x, Span<double> gradient);

    /// <summary>C under <paramref name="measure"/> in closed form, in the inputs' own units; null where the family has none
    /// under that law.</summary>
    public abstract double[,]? C(ReferenceMeasure measure);

    public override string ToString() => Name;
}

/// <summary>A ridge's profile g, with its derivative.</summary>
internal abstract record RidgeProfile
{
    public abstract double Value(double t);

    public abstract double Slope(double t);
}

/// <summary>
/// g(t) = Re(A·e^{λt}): an exponential at a real rate, a wave at an imaginary one (cos at A = 1, sin at A = −i), a damped wave
/// between. Products of two are exponentials again, so a sum of such ridges has C in closed form under any law with a moment
/// generating function.
/// </summary>
internal sealed record ExponentialProfile(Complex Amplitude, Complex Rate) : RidgeProfile
{
    public static ExponentialProfile Exp(double rate, double amplitude = 1.0) => new(amplitude, rate);

    public static ExponentialProfile Sin(double frequency, double amplitude = 1.0) => new(new Complex(0.0, -amplitude), new Complex(0.0, frequency));

    public static ExponentialProfile Cos(double frequency, double amplitude = 1.0) => new(amplitude, new Complex(0.0, frequency));

    public override double Value(double t) => (Amplitude * Complex.Exp(Rate * t)).Real;

    public override double Slope(double t) => (Amplitude * Rate * Complex.Exp(Rate * t)).Real;
}

/// <summary>A profile given by its function and derivative, with no closed form for C: its ridge's direction is still known.</summary>
internal sealed record FunctionProfile(string Name, Func<double, double> G, Func<double, double> DG) : RidgeProfile
{
    public override double Value(double t) => G(t);

    public override double Slope(double t) => DG(t);
}

/// <summary>One ridge, g(wᵀx − c).</summary>
internal sealed record RidgeTerm(RidgeProfile Profile, double[] Direction, double Offset = 0.0);

/// <summary>
/// A sum of ridges f(x) = Σ_k g_k(w_kᵀx − c_k), whose gradient Σ_k g_k′ w_k lies in the directions' span, so C's range does too.
/// With exponential profiles, g_k′ = Re(a_k e^{λ_k s_k}) with a_k = A_kλ_k, and
/// E[g_k′g_l′] = ½ Re(a_k a_l E[e^{λ_k s_k + λ_l s_l}] + a_k ā_l E[e^{λ_k s_k + λ̄_l s_l}]), each expectation
/// e^{α(w_kᵀm − c_k) + β(w_lᵀm − c_l)}·K(αw_k + βw_l) about the law's centre m; C = Σ_{k,l} E[g_k′g_l′] w_k w_lᵀ.
/// </summary>
internal sealed class RidgeSum : ReferenceFunction
{
    private readonly RidgeTerm[] _ridges;

    public RidgeSum(string name, params RidgeTerm[] ridges) : base(name, ridges[0].Direction.Length)
    {
        if (ridges.Any(r => r.Direction.Length != Dimension)) throw new ArgumentException("Every ridge's direction has the same inputs.");
        _ridges = ridges;
    }

    public IReadOnlyList<RidgeTerm> Ridges => _ridges;

    public override double Value(ReadOnlySpan<double> x)
    {
        double f = 0.0;
        foreach (RidgeTerm r in _ridges) f += r.Profile.Value(Argument(r, x));
        return f;
    }

    public override void Gradient(ReadOnlySpan<double> x, Span<double> gradient)
    {
        gradient.Clear();
        foreach (RidgeTerm r in _ridges)
        {
            double slope = r.Profile.Slope(Argument(r, x));
            for (int j = 0; j < Dimension; j++) gradient[j] += slope * r.Direction[j];
        }
    }

    public override double[,]? C(ReferenceMeasure measure)
    {
        if (measure.Dimension != Dimension) throw new ArgumentException("The law has other inputs than the function.");
        if (!_ridges.All(r => r.Profile is ExponentialProfile)) return null;
        int p = Dimension, count = _ridges.Length;
        var shift = new double[count];   // w_kᵀm − c_k
        for (int k = 0; k < count; k++)
        {
            double s = 0.0;
            for (int j = 0; j < p; j++) s += _ridges[k].Direction[j] * measure.Center[j];
            shift[k] = s - _ridges[k].Offset;
        }
        var mu = new Complex[p];
        Complex Moment(Complex alpha, int k, Complex beta, int l)
        {
            for (int j = 0; j < p; j++) mu[j] = alpha * _ridges[k].Direction[j] + beta * _ridges[l].Direction[j];
            return Complex.Exp(alpha * shift[k] + beta * shift[l]) * measure.CenteredMgf(mu);
        }

        var e = new double[count, count];
        for (int k = 0; k < count; k++)
            for (int l = 0; l <= k; l++)
            {
                var (gk, gl) = ((ExponentialProfile)_ridges[k].Profile, (ExponentialProfile)_ridges[l].Profile);
                Complex ak = gk.Amplitude * gk.Rate, al = gl.Amplitude * gl.Rate;
                Complex same = Moment(gk.Rate, k, gl.Rate, l), mixed = Moment(gk.Rate, k, Complex.Conjugate(gl.Rate), l);
                e[k, l] = e[l, k] = 0.5 * (ak * al * same + ak * Complex.Conjugate(al) * mixed).Real;
            }
        var c = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
            {
                double s = 0.0;
                for (int k = 0; k < count; k++)
                    for (int l = 0; l < count; l++) s += e[k, l] * _ridges[k].Direction[i] * _ridges[l].Direction[j];
                c[i, j] = s;
            }
        return c;
    }

    private double Argument(RidgeTerm r, ReadOnlySpan<double> x)
    {
        double s = 0.0;
        for (int j = 0; j < Dimension; j++) s += r.Direction[j] * x[j];
        return s - r.Offset;
    }
}

/// <summary>One monomial c·Π_j y_j^{e_j}.</summary>
internal sealed record Monomial(double Coefficient, int[] Exponents);

/// <summary>
/// A polynomial in y = x − o about an origin o, Σ_t c_t Π_j y_j^{e_tj}. Under the product uniform law on a box each entry of C is
/// a sum of products of one-variable moments: C_ik = Σ_{t,s} c_t c_s e_ti e_sk Π_j E[y_j^{γ_j}], γ = e_t + e_s − 1_i − 1_k, and
/// E[(d + h u)^n] = Σ_{r even} C(n, r) d^{n−r} h^r/(r + 1) for u uniform on [−1, 1], whose terms share one sign.
/// </summary>
internal sealed class PolynomialReference : ReferenceFunction
{
    private readonly Monomial[] _terms;
    private readonly double[] _origin;

    public PolynomialReference(string name, IReadOnlyList<Monomial> terms, double[]? origin = null) : base(name, terms[0].Exponents.Length)
    {
        if (terms.Any(t => t.Exponents.Length != Dimension || t.Exponents.Any(e => e < 0)))
            throw new ArgumentException("Every monomial has a nonnegative exponent per input.");
        _terms = [.. terms];
        _origin = origin ?? new double[Dimension];
    }

    public IReadOnlyList<Monomial> Terms => _terms;

    /// <summary>The ridge g(wᵀx − c), g(t) = Σ_k a_k t^k with <paramref name="profile"/> the a_k: <see cref="OfLinearForms"/>
    /// with one form.</summary>
    public static PolynomialReference Ridge(string name, double[] direction, double offset, double[] profile, double[]? origin = null) =>
        OfLinearForms(name, [.. profile.Select((a, k) => new Monomial(a, [k]))], [direction], [offset], origin);

    /// <summary>
    /// h(W x − c) for a polynomial h in r variables, given by its monomials <paramref name="outer"/>, and r linear forms
    /// w_kᵀx − c_k (<paramref name="forms"/>, <paramref name="offsets"/>), expanded about <paramref name="origin"/>:
    /// w_kᵀx − c_k = w_kᵀy + (w_kᵀo − c_k). Its gradient Wᵀ∇h lies in the forms' span, so C's range does too.
    /// </summary>
    public static PolynomialReference OfLinearForms(string name, IReadOnlyList<Monomial> outer, IReadOnlyList<double[]> forms,
                                                    IReadOnlyList<double> offsets, double[]? origin = null)
    {
        int p = forms[0].Length, r = forms.Count;
        origin ??= new double[p];
        var one = new Dictionary<string, Monomial> { [Key(new int[p])] = new Monomial(1.0, new int[p]) };
        var linear = new Dictionary<string, Monomial>[r];
        for (int k = 0; k < r; k++)
        {
            double shift = -offsets[k];
            for (int j = 0; j < p; j++) shift += forms[k][j] * origin[j];
            linear[k] = [];
            Accumulate(linear[k], new Monomial(shift, new int[p]));
            for (int j = 0; j < p; j++)
            {
                var e = new int[p];
                e[j] = 1;
                Accumulate(linear[k], new Monomial(forms[k][j], e));
            }
        }
        var result = new Dictionary<string, Monomial>();
        foreach (Monomial m in outer)
        {
            if (m.Exponents.Length != r) throw new ArgumentException("Each outer monomial has an exponent per form.");
            Dictionary<string, Monomial> term = one;
            for (int k = 0; k < r; k++)
                for (int e = 0; e < m.Exponents[k]; e++) term = Multiply(term, linear[k]);
            foreach (Monomial t in term.Values) Accumulate(result, t with { Coefficient = m.Coefficient * t.Coefficient });
        }
        Monomial[] kept = [.. result.Values.Where(m => m.Coefficient != 0.0)];
        return new PolynomialReference(name, kept.Length > 0 ? kept : [new Monomial(0.0, new int[p])], origin);
    }

    /// <summary>½ yᵀAy, whose gradient is Ay for a symmetric A.</summary>
    public static PolynomialReference Quadratic(string name, double[,] a, double[]? origin = null)
    {
        int p = a.GetLength(0);
        var terms = new List<Monomial>();
        for (int i = 0; i < p; i++)
            for (int j = i; j < p; j++)
            {
                var e = new int[p];
                e[i]++;
                e[j]++;
                terms.Add(new Monomial(i == j ? 0.5 * a[i, i] : 0.5 * (a[i, j] + a[j, i]), e));
            }
        return new PolynomialReference(name, terms, origin);
    }

    public override double Value(ReadOnlySpan<double> x)
    {
        double f = 0.0;
        foreach (Monomial t in _terms)
        {
            double v = t.Coefficient;
            for (int j = 0; j < Dimension; j++) v *= Power(x[j] - _origin[j], t.Exponents[j]);
            f += v;
        }
        return f;
    }

    public override void Gradient(ReadOnlySpan<double> x, Span<double> gradient)
    {
        gradient.Clear();
        foreach (Monomial t in _terms)
            for (int i = 0; i < Dimension; i++)
            {
                if (t.Exponents[i] == 0) continue;
                double v = t.Coefficient * t.Exponents[i];
                for (int j = 0; j < Dimension; j++) v *= Power(x[j] - _origin[j], t.Exponents[j] - (j == i ? 1 : 0));
                gradient[i] += v;
            }
    }

    public override double[,]? C(ReferenceMeasure measure)
    {
        if (measure.Dimension != Dimension) throw new ArgumentException("The law has other inputs than the function.");
        if (measure is not ReferenceBox box) return null;
        int p = Dimension;
        var moments = new double[p][];
        for (int j = 0; j < p; j++)
            moments[j] = UniformMoments(box.Center[j] - _origin[j], box.HalfWidths[j], 2 * _terms.Max(t => t.Exponents[j]));
        var c = new double[p, p];
        var gamma = new int[p];
        for (int i = 0; i < p; i++)
            for (int k = i; k < p; k++)
            {
                double s = 0.0;
                foreach (Monomial t in _terms)
                {
                    if (t.Exponents[i] == 0) continue;
                    foreach (Monomial u in _terms)
                    {
                        if (u.Exponents[k] == 0) continue;
                        for (int j = 0; j < p; j++) gamma[j] = t.Exponents[j] + u.Exponents[j];
                        gamma[i]--;
                        gamma[k]--;
                        double v = t.Coefficient * u.Coefficient * t.Exponents[i] * u.Exponents[k];
                        for (int j = 0; j < p; j++) v *= moments[j][gamma[j]];
                        s += v;
                    }
                }
                c[i, k] = c[k, i] = s;
            }
        return c;
    }

    // E[(d + h·u)^n] for n ≤ top, u uniform on [−1, 1].
    private static double[] UniformMoments(double d, double h, int top)
    {
        var m = new double[top + 1];
        for (int n = 0; n <= top; n++)
        {
            double s = 0.0, binomial = 1.0;
            for (int r = 0; r <= n; r++)
            {
                if (r % 2 == 0) s += binomial * Power(d, n - r) * Power(h, r) / (r + 1);
                binomial = binomial * (n - r) / (r + 1);
            }
            m[n] = s;
        }
        return m;
    }

    private static double Power(double y, int n)
    {
        double v = 1.0;
        for (int k = 0; k < n; k++) v *= y;
        return v;
    }

    private static string Key(int[] exponents) => string.Join(",", exponents);

    private static void Accumulate(Dictionary<string, Monomial> into, Monomial m)
    {
        string key = Key(m.Exponents);
        into[key] = into.TryGetValue(key, out Monomial? held) ? held with { Coefficient = held.Coefficient + m.Coefficient } : m;
    }

    private static Dictionary<string, Monomial> Multiply(Dictionary<string, Monomial> a, Dictionary<string, Monomial> b)
    {
        var product = new Dictionary<string, Monomial>();
        foreach (Monomial x in a.Values)
            foreach (Monomial y in b.Values)
                Accumulate(product, new Monomial(x.Coefficient * y.Coefficient, [.. x.Exponents.Select((e, j) => e + y.Exponents[j])]));
        return product;
    }
}

/// <summary>One term b·Π_j ψ_{α_j} of a Legendre expansion.</summary>
internal sealed record LegendreTerm(double Coefficient, int[] Degrees);

/// <summary>
/// A Legendre expansion on a box, f(x) = Σ_α b_α Π_j ψ_{α_j}(u_j) with u_j = (x_j − a_j)/(b_j − a_j) and ψ_n = √(2n + 1)·P_n(2u − 1)
/// orthonormal on [0, 1], evaluated by its own recurrence. Under the uniform law on its own box, ψ′_m = Σ_{n &lt; m, m − n odd}
/// 2√((2m + 1)(2n + 1)) ψ_n gives ∫ψ′_mψ_n = 2√((2m + 1)(2n + 1)) for n &lt; m with m − n odd, and ∫ψ′_mψ′_n =
/// 2√((2m + 1)(2n + 1))·r(r + 1), r = min(m, n), for m − n even; so C_ik = Σ_{α,β} b_α b_β ∫∂_iψ_α ∂_kψ_β / ((b_i − a_i)(b_k − a_k))
/// with the other inputs' degrees matched.
/// </summary>
internal sealed class LegendreExpansion : ReferenceFunction
{
    private readonly LegendreTerm[] _terms;

    public LegendreExpansion(string name, ReferenceBox box, IReadOnlyList<LegendreTerm> terms) : base(name, box.Dimension)
    {
        if (terms.Any(t => t.Degrees.Length != Dimension || t.Degrees.Any(d => d < 0)))
            throw new ArgumentException("Every term has a nonnegative degree per input.");
        Box = box;
        _terms = [.. terms];
    }

    public ReferenceBox Box { get; }

    public IReadOnlyList<LegendreTerm> Terms => _terms;

    public override double Value(ReadOnlySpan<double> x)
    {
        double f = 0.0;
        foreach (LegendreTerm t in _terms)
        {
            double v = t.Coefficient;
            for (int j = 0; j < Dimension; j++) v *= Psi(t.Degrees[j], Unit(x, j)).Value;
            f += v;
        }
        return f;
    }

    public override void Gradient(ReadOnlySpan<double> x, Span<double> gradient)
    {
        gradient.Clear();
        Span<double> value = stackalloc double[Dimension], slope = stackalloc double[Dimension];
        foreach (LegendreTerm t in _terms)
        {
            for (int j = 0; j < Dimension; j++) (value[j], slope[j]) = Psi(t.Degrees[j], Unit(x, j));
            for (int i = 0; i < Dimension; i++)
            {
                double v = t.Coefficient * slope[i] / (Box.Upper[i] - Box.Lower[i]);
                for (int j = 0; j < Dimension; j++)
                    if (j != i) v *= value[j];
                gradient[i] += v;
            }
        }
    }

    public override double[,]? C(ReferenceMeasure measure)
    {
        if (measure.Dimension != Dimension) throw new ArgumentException("The law has other inputs than the function.");
        if (measure is not ReferenceBox box || !box.SameAs(Box)) return null;
        int p = Dimension;
        var c = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int k = i; k < p; k++)
            {
                double s = 0.0;
                foreach (LegendreTerm a in _terms)
                    foreach (LegendreTerm b in _terms)
                    {
                        bool matched = true;
                        for (int j = 0; j < p && matched; j++)
                            if (j != i && j != k && a.Degrees[j] != b.Degrees[j]) matched = false;
                        if (!matched) continue;
                        double integral = i == k
                            ? SlopeGram(a.Degrees[i], b.Degrees[i])
                            : SlopeProduct(a.Degrees[i], b.Degrees[i]) * SlopeProduct(b.Degrees[k], a.Degrees[k]);
                        s += a.Coefficient * b.Coefficient * integral;
                    }
                c[i, k] = c[k, i] = s / ((Box.Upper[i] - Box.Lower[i]) * (Box.Upper[k] - Box.Lower[k]));
            }
        return c;
    }

    /// <summary>(ψ_n(u), ψ′_n(u)) by Bonnet's recurrence and P′_{k+1} = P′_{k−1} + (2k + 1)P_k.</summary>
    public static (double Value, double Slope) Psi(int n, double u)
    {
        if (n == 0) return (1.0, 0.0);
        double t = 2.0 * u - 1.0;
        double p0 = 1.0, p1 = t, d0 = 0.0, d1 = 1.0;
        for (int k = 1; k < n; k++)
        {
            double p2 = ((2 * k + 1) * t * p1 - k * p0) / (k + 1);
            double d2 = d0 + (2 * k + 1) * p1;
            (p0, p1, d0, d1) = (p1, p2, d1, d2);
        }
        double s = Math.Sqrt(2.0 * n + 1.0);
        return (s * p1, 2.0 * s * d1);
    }

    /// <summary>∫₀¹ ψ′_m ψ_n.</summary>
    public static double SlopeProduct(int m, int n) =>
        n < m && ((m - n) & 1) == 1 ? 2.0 * Math.Sqrt((2.0 * m + 1.0) * (2.0 * n + 1.0)) : 0.0;

    /// <summary>∫₀¹ ψ′_m ψ′_n.</summary>
    public static double SlopeGram(int m, int n)
    {
        if (((m - n) & 1) == 1) return 0.0;
        int r = Math.Min(m, n);
        return 2.0 * Math.Sqrt((2.0 * m + 1.0) * (2.0 * n + 1.0)) * r * (r + 1);
    }

    private double Unit(ReadOnlySpan<double> x, int j) => (x[j] - Box.Lower[j]) / (Box.Upper[j] - Box.Lower[j]);
}

/// <summary>One input's quadrature rule for a probability law: nodes, and weights summing to 1.</summary>
internal sealed record QuadratureRule(double[] Nodes, double[] Weights)
{
    /// <summary>The uniform law on [<paramref name="lower"/>, <paramref name="upper"/>]: <paramref name="nodes"/>-point
    /// Gauss–Legendre on each of <paramref name="panels"/> equal panels.</summary>
    public static QuadratureRule Uniform(double lower, double upper, int nodes, int panels = 1) => Density(null, lower, upper, nodes, panels);

    /// <summary>The law with density proportional to <paramref name="density"/> on [<paramref name="lower"/>,
    /// <paramref name="upper"/>], by the same panels, renormalized by the rule's own mass.</summary>
    public static QuadratureRule Density(Func<double, double>? density, double lower, double upper, int nodes, int panels = 1)
    {
        if (!(upper > lower) || panels < 1) throw new ArgumentException("A rule needs an interval and a panel.");
        GaussLegendre rule = GaussLegendre.Of(nodes);
        var x = new double[nodes * panels];
        var w = new double[nodes * panels];
        double width = (upper - lower) / panels, mass = 0.0;
        for (int q = 0; q < panels; q++)
        {
            double a = lower + q * width, b = q == panels - 1 ? upper : lower + (q + 1) * width;
            double mid = 0.5 * (a + b), half = 0.5 * (b - a);
            for (int i = 0; i < nodes; i++)
            {
                int at = q * nodes + i;
                x[at] = mid + half * rule.Nodes[i];
                w[at] = half * rule.Weights[i] * (density?.Invoke(x[at]) ?? 1.0);
                mass += w[at];
            }
        }
        for (int i = 0; i < w.Length; i++) w[i] /= mass;
        return new QuadratureRule(x, w);
    }
}

/// <summary>A gradient field ∇f, written into <c>gradient</c> at <c>x</c>.</summary>
internal delegate void GradientField(ReadOnlySpan<double> x, Span<double> gradient);

/// <summary>
/// C = E[∇f ∇fᵀ] by tensor-product quadrature: the sum of w·∇f∇fᵀ over the product of one rule per input, each entry a
/// compensated (Neumaier) sum. It reads any gradient field, so it serves as brute force for any C on few inputs.
/// </summary>
internal static class GradientQuadrature
{
    public static double[,] C(GradientField gradient, IReadOnlyList<QuadratureRule> rules)
    {
        int p = rules.Count;
        var index = new int[p];
        var x = new double[p];
        var g = new double[p];
        var sum = new double[p * p];
        var carry = new double[p * p];
        while (true)
        {
            double w = 1.0;
            for (int j = 0; j < p; j++)
            {
                x[j] = rules[j].Nodes[index[j]];
                w *= rules[j].Weights[index[j]];
            }
            gradient(x, g);
            for (int a = 0; a < p; a++)
            {
                double wa = w * g[a];
                for (int b = a; b < p; b++) Add(ref sum[a * p + b], ref carry[a * p + b], wa * g[b]);
            }
            int k = p - 1;
            while (k >= 0 && ++index[k] == rules[k].Nodes.Length) index[k--] = 0;
            if (k < 0) break;
        }
        var c = new double[p, p];
        for (int a = 0; a < p; a++)
            for (int b = a; b < p; b++) c[a, b] = c[b, a] = sum[a * p + b] + carry[a * p + b];
        return c;
    }

    /// <summary>Under the product uniform law on <paramref name="box"/>: <paramref name="nodes"/> Gauss–Legendre nodes on each of
    /// <paramref name="panels"/>[j] panels of input j.</summary>
    public static double[,] Uniform(ReferenceFunction f, ReferenceBox box, int nodes, IReadOnlyList<int> panels) =>
        C(f.Gradient, [.. Enumerable.Range(0, f.Dimension).Select(j => QuadratureRule.Uniform(box.Lower[j], box.Upper[j], nodes, panels[j]))]);

    /// <summary>Under N(m, Σ): over z standard normal on [−<paramref name="radius"/>, <paramref name="radius"/>] per input, in
    /// <paramref name="panels"/> panels of <paramref name="nodes"/> nodes, at x = m + Lz. The mass beyond the radius is left out.</summary>
    public static double[,] Gaussian(ReferenceFunction f, ReferenceGaussian law, double radius, int nodes, int panels)
    {
        int p = f.Dimension;
        QuadratureRule rule = QuadratureRule.Density(z => Math.Exp(-0.5 * z * z), -radius, radius, nodes, panels);
        var x = new double[p];
        return C((z, g) =>
        {
            for (int i = 0; i < p; i++)
            {
                double v = law.Center[i];
                for (int j = 0; j <= i; j++) v += law.Factor[i, j] * z[j];
                x[i] = v;
            }
            f.Gradient(x, g);
        }, [.. Enumerable.Repeat(rule, p)]);
    }

    private static void Add(ref double sum, ref double carry, double value)
    {
        double t = sum + value;
        carry += Math.Abs(sum) >= Math.Abs(value) ? sum - t + value : value - t + sum;
        sum = t;
    }
}
