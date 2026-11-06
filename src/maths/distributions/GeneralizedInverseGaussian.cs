using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The generalized inverse Gaussian GIG(p, a, b), density ∝ x^(p−1)·exp(−(a·x + b/x)/2) on x &gt; 0 (Barndorff-Nielsen's family;
/// Jørgensen 1982), with a ≥ 0, b ≥ 0, a &gt; 0 when p ≥ 0 and b &gt; 0 when p ≤ 0. The limits are the gamma, b = 0: Gamma(p, rate a/2);
/// and the inverse gamma, a = 0: IG(−p, scale b/2).
/// <para><b>Draws.</b> In u = ln x the log density is p·u − (a·eᵘ + b·e⁻ᵘ)/2, whose second derivative −(a·eᵘ + b·e⁻ᵘ)/2 is negative:
/// the GIG is log-concave in ln x for every parameter. So one exact method serves every case, far ranges and intervals included: the
/// tangent hat of <see cref="LogConcaveEnvelope"/> about the mode in u, with the log ratio written about it as a sum of two non-positive terms,
/// so it keeps its digits (below); a draw on an interval is the same hat cut to it. An untruncated draw is Hörmann and Leydold's (2014),
/// GIGrvg's, where <see cref="DrawsByRatio"/> takes it: ratio-of-uniforms (<see cref="RatioOfUniforms"/>) with and without a shift by the
/// mode, and a dominating density of three pieces where the GIG is not T₋½-concave. It needs no search for a design point, and on the
/// parameters fits produce it is 2.5 to 5.6 times faster a draw. The limits draw their own families.</para>
/// <para><b>Normalizer and moments.</b> The normalizer is Bessel K of real order (<see cref="SpecialFunctions.LogBesselKScaled"/>), closed at the
/// limits; the mean and variance are its ratios, taken in the expon-scaled logs so the exponentials cancel exactly and the variance's
/// K_{p+2}K_p/K²_{p+1} − 1 keeps its digits at large √(ab).</para>
/// A readonly value type.
/// </summary>
public readonly struct GeneralizedInverseGaussian
{
    private GeneralizedInverseGaussian(double p, double a, double b)
    {
        P = p;
        A = a;
        B = b;
    }

    /// <summary>GIG(<paramref name="p"/>, <paramref name="a"/>, <paramref name="b"/>), each finite, a and b non-negative, a &gt; 0 when
    /// p ≥ 0 and b &gt; 0 when p ≤ 0.</summary>
    public static GeneralizedInverseGaussian Of(double p, double a, double b)
    {
        if (!double.IsFinite(p)) throw new ArgumentOutOfRangeException(nameof(p), "p must be finite.");
        if (!(a >= 0.0 && double.IsFinite(a))) throw new ArgumentOutOfRangeException(nameof(a), "a must be non-negative and finite.");
        if (!(b >= 0.0 && double.IsFinite(b))) throw new ArgumentOutOfRangeException(nameof(b), "b must be non-negative and finite.");
        if (p >= 0.0 && !(a > 0.0)) throw new ArgumentException("A GIG with p ≥ 0 needs a > 0.", nameof(a));
        if (p <= 0.0 && !(b > 0.0)) throw new ArgumentException("A GIG with p ≤ 0 needs b > 0.", nameof(b));
        return new GeneralizedInverseGaussian(p, a, b);
    }

    /// <summary>p.</summary>
    public double P { get; }

    /// <summary>a, the coefficient of x.</summary>
    public double A { get; }

    /// <summary>b, the coefficient of 1/x.</summary>
    public double B { get; }

    /// <summary>ln of the unnormalized density, (p − 1)·ln x − (a·x + b/x)/2; −∞ for x ≤ 0.</summary>
    public double LogKernel(double x) => x > 0.0 ? (P - 1.0) * Math.Log(x) - 0.5 * (A * x + B / x) : double.NegativeInfinity;

    /// <summary>The density's mode, the positive root of a·x² − 2(p − 1)x − b = 0; 0 at b = 0 with p ≤ 1.</summary>
    public double Mode
    {
        get
        {
            double q = P - 1.0;
            if (B == 0.0) return q > 0.0 ? 2.0 * q / A : 0.0;
            if (A == 0.0) return B / (2.0 * (1.0 - P));
            // The root without cancellation: (q + √(q² + ab))/a, or b/(√(q² + ab) − q) when q < 0.
            double root = Math.Sqrt(q * q + A * B);
            return q >= 0.0 ? (q + root) / A : B / (root - q);
        }
    }

    /// <summary>
    /// ln ∫₀^∞ x^(p−1)·exp(−(a·x + b/x)/2) dx, the normalizer of GIG(<paramref name="p"/>, <paramref name="a"/>, <paramref name="b"/>) for any
    /// p and finite a, b ≥ 0: ln 2 + (p/2)·ln(b/a) + ln K_p(√(ab)) where both are positive; ln Γ(p) + p·ln(2/a) at b = 0 with p &gt; 0, the
    /// gamma's; ln Γ(−p) + p·ln(b/2) at a = 0 with p &lt; 0, the inverse gamma's; +∞ where the integral diverges.
    /// </summary>
    public static double LogNormalizer(double p, double a, double b)
    {
        if (double.IsNaN(p)) return double.NaN;
        if (!(a >= 0.0 && double.IsFinite(a))) throw new ArgumentOutOfRangeException(nameof(a), "a must be non-negative and finite.");
        if (!(b >= 0.0 && double.IsFinite(b))) throw new ArgumentOutOfRangeException(nameof(b), "b must be non-negative and finite.");
        if (a > 0.0 && b > 0.0)
        {
            double omega = Math.Sqrt(a) * Math.Sqrt(b);
            return SpecialFunctions.Ln2 + 0.5 * p * (Math.Log(b) - Math.Log(a)) + SpecialFunctions.LogBesselKScaled(p, omega) - omega;
        }
        if (a > 0.0) return p > 0.0 ? SpecialFunctions.LogGamma(p) + p * (SpecialFunctions.Ln2 - Math.Log(a)) : double.PositiveInfinity;
        if (b > 0.0) return p < 0.0 ? SpecialFunctions.LogGamma(-p) + p * (Math.Log(b) - SpecialFunctions.Ln2) : double.PositiveInfinity;
        return double.PositiveInfinity;
    }

    /// <summary>ln of the density at <paramref name="x"/>; −∞ for x ≤ 0.</summary>
    public double LogPdf(double x) => LogKernel(x) - LogNormalizer(P, A, B);

    /// <summary>E[x]: √(b/a)·K_{p+1}(√(ab))/K_p(√(ab)); 2p/a at b = 0; b/(2(−p − 1)) at a = 0 for p &lt; −1, +∞ otherwise.</summary>
    public double Mean
    {
        get
        {
            if (B == 0.0) return 2.0 * P / A;
            if (A == 0.0) return P < -1.0 ? 0.5 * B / (-P - 1.0) : double.PositiveInfinity;
            double omega = Math.Sqrt(A) * Math.Sqrt(B);
            return Math.Exp(0.5 * (Math.Log(B) - Math.Log(A)) + SpecialFunctions.LogBesselKScaled(P + 1.0, omega) - SpecialFunctions.LogBesselKScaled(P, omega));
        }
    }

    /// <summary>Var[x]: E[x]²·(K_{p+2}K_p/K²_{p+1} − 1); 4p/a² at b = 0; (b/2)²/((−p − 1)²(−p − 2)) at a = 0 for p &lt; −2, +∞ otherwise.</summary>
    public double Variance
    {
        get
        {
            if (B == 0.0) return 4.0 * P / (A * A);
            if (A == 0.0) return P < -2.0 ? 0.25 * B * B / ((-P - 1.0) * (-P - 1.0) * (-P - 2.0)) : double.PositiveInfinity;
            double omega = Math.Sqrt(A) * Math.Sqrt(B);
            double k0 = SpecialFunctions.LogBesselKScaled(P, omega), k1 = SpecialFunctions.LogBesselKScaled(P + 1.0, omega);
            double k2 = SpecialFunctions.LogBesselKScaled(P + 2.0, omega), mean = Mean;
            return mean * mean * SpecialFunctions.Expm1(k2 + k0 - 2.0 * k1);
        }
    }

    /// <summary>One draw: the gamma or inverse gamma at a limit; otherwise Hörmann and Leydold's ratio-of-uniforms where
    /// <see cref="DrawsByRatio"/> takes it, the tangent hat in ln x where not.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        return DrawsByRatio(P, A, B) ? SampleByRatio(rng) : SampleOn(rng, 0.0, double.PositiveInfinity);
    }

    /// <summary>
    /// Whether an untruncated draw of GIG(<paramref name="p"/>, <paramref name="a"/>, <paramref name="b"/>) is Hörmann and Leydold's: both
    /// a and b positive, |p| ≤ 10³ and 10⁻⁸ ≤ ω = √(ab) ≤ 10⁸, the range its law has been checked over, where it was the faster draw at every
    /// point measured, 2.5 to 5.6 times on the parameters quantile and Normal-Wald fits produce. A function of the parameters alone, so a
    /// chain's draws do not depend on its history and a resumed run draws as the whole one did.
    /// </summary>
    public static bool DrawsByRatio(double p, double a, double b)
    {
        if (!(a > 0.0 && b > 0.0 && Math.Abs(p) <= 1e3)) return false;
        double omega = Math.Sqrt(a) * Math.Sqrt(b);
        return omega >= 1e-8 && omega <= 1e8;
    }

    /// <summary>
    /// Hörmann and Leydold's (2014) draw, GIGrvg's: in the standard form y ~ GIG(λ, ω, ω) ∝ y^(λ−1)·e^(−ω(y + 1/y)/2), λ = |p| by y ↦ 1/y,
    /// ω = √(ab), x = √(b/a)·y. Ratio-of-uniforms with a shift by the mode where λ &gt; 2 or ω &gt; 3, its rectangle's sides at the roots of
    /// a cubic by Cardano; without the shift where the density is T₋½-concave, λ ≥ 1 − 9ω²/4, or ω &gt; 0.2, its sides at 0 and the
    /// maximum of y·√f(y); and below both, where it is not T₋½-concave, a hat of three pieces, constant to ω/(1 − λ), then a power, then
    /// an exponential. Each is exact; untruncated.
    /// </summary>
    internal double SampleByRatio(Xoshiro256PlusPlus rng)
    {
        double lambda = Math.Abs(P), omega = Math.Sqrt(A) * Math.Sqrt(B), scale = Math.Sqrt(B) / Math.Sqrt(A);
        double y = lambda > 2.0 || omega > 3.0 ? ShiftedRatio(rng, lambda, omega)
                 : lambda >= 1.0 - 2.25 * omega * omega || omega > 0.2 ? Ratio(rng, lambda, omega)
                 : Dominated(rng, lambda, omega);
        return P < 0.0 ? scale / y : scale * y;
    }


    // The standard form's mode: ((λ − 1) + √((λ − 1)² + ω²))/ω, or for λ < 1 the same root written without cancellation.
    private static double StandardMode(double lambda, double omega)
        => lambda >= 1.0 ? (Math.Sqrt((lambda - 1.0) * (lambda - 1.0) + omega * omega) + (lambda - 1.0)) / omega
                         : omega / (Math.Sqrt((1.0 - lambda) * (1.0 - lambda) + omega * omega) + (1.0 - lambda));

    // Without the shift: u on [0, y_m·√(f(y_m)/f(mode))], y_m = ((λ + 1) + √((λ + 1)² + ω²))/ω the maximum of y·√f(y).
    private static double Ratio(Xoshiro256PlusPlus rng, double lambda, double omega)
    {
        double mode = StandardMode(lambda, omega);
        var h = new StandardRatio(lambda, omega, mode);
        double top = ((lambda + 1.0) + Math.Sqrt((lambda + 1.0) * (lambda + 1.0) + omega * omega)) / omega;
        return RatioOfUniforms.Draw(rng, h, 0.0, 0.0, top * Math.Exp(0.5 * h.Value(top)));
    }

    // With the shift by the mode m: the extremes of (y − m)·√f(y) are the roots in (0, m) and (m, ∞) of
    // y³ − (2(λ + 1)/ω + m)y² + (2(λ − 1)m/ω − 1)y + m = 0, by Cardano's rule in its trigonometric form (three real roots).
    private static double ShiftedRatio(Xoshiro256PlusPlus rng, double lambda, double omega)
    {
        double mode = StandardMode(lambda, omega);
        var h = new StandardRatio(lambda, omega, mode);
        double a = -(2.0 * (lambda + 1.0) / omega + mode), b = 2.0 * (lambda - 1.0) * mode / omega - 1.0, c = mode;
        double p = b - a * a / 3.0, q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
        double phi = Math.Acos(Math.Clamp(-q / (2.0 * Math.Sqrt(-(p * p * p) / 27.0)), -1.0, 1.0)), factor = 2.0 * Math.Sqrt(-p / 3.0);
        double above = factor * Math.Cos(phi / 3.0) - a / 3.0, below = factor * Math.Cos(phi / 3.0 + 4.0 / 3.0 * Math.PI) - a / 3.0;
        double plus = (above - mode) * Math.Exp(0.5 * h.Value(above)), minus = (below - mode) * Math.Exp(0.5 * h.Value(below));
        return RatioOfUniforms.Draw(rng, h, mode, minus, plus);
    }

    // λ < 1 and ω small, where the density is not T₋½-concave: rejection from a hat of three pieces, the density's maximum on [0, x₀],
    // x₀ = ω/(1 − λ); e^(−ω)·y^(λ−1) on [x₀, 2/ω]; (2/ω)^(λ−1)·e^(−ωy/2) beyond, or x₀^(λ−1)·e^(−ωy/2) from x₀ where x₀ ≥ 2/ω.
    private static double Dominated(Xoshiro256PlusPlus rng, double lambda, double omega)
    {
        double mode = StandardMode(lambda, omega), x0 = omega / (1.0 - lambda), edge = 2.0 / omega;
        double k0 = Math.Exp((lambda - 1.0) * Math.Log(mode) - 0.5 * omega * (mode + 1.0 / mode)), k1, k2;
        double area0 = k0 * x0, area1, area2;
        if (x0 >= edge)
        {
            k1 = 0.0;
            area1 = 0.0;
            k2 = Math.Pow(x0, lambda - 1.0);
            area2 = k2 * 2.0 * Math.Exp(-omega * x0 / 2.0) / omega;
        }
        else
        {
            k1 = Math.Exp(-omega);
            area1 = lambda == 0.0 ? k1 * Math.Log(2.0 / (omega * omega)) : k1 / lambda * (Math.Pow(edge, lambda) - Math.Pow(x0, lambda));
            k2 = Math.Pow(edge, lambda - 1.0);
            area2 = k2 * 2.0 * Math.Exp(-1.0) / omega;
        }
        double total = area0 + area1 + area2, start = Math.Max(x0, edge);
        for (int attempt = 0; attempt < 100_000; attempt++)
        {
            double pick = total * rng.NextDouble(), x, hat;
            if (pick <= area0)
            {
                x = x0 * pick / area0;
                hat = k0;
            }
            else if ((pick -= area0) <= area1)
            {
                x = lambda == 0.0 ? omega * Math.Exp(Math.Exp(omega) * pick) : Math.Pow(Math.Pow(x0, lambda) + lambda / k1 * pick, 1.0 / lambda);
                hat = k1 * Math.Pow(x, lambda - 1.0);
            }
            else
            {
                pick -= area1;
                x = -edge * Math.Log(Math.Exp(-omega / 2.0 * start) - omega / (2.0 * k2) * pick);
                hat = k2 * Math.Exp(-omega / 2.0 * x);
            }
            double u = (1.0 - rng.NextDouble()) * hat;
            if (Math.Log(u) <= (lambda - 1.0) * Math.Log(x) - 0.5 * omega * (x + 1.0 / x)) return x;
        }
        throw new ArithmeticException($"The GIG's dominated draw at λ = {lambda}, ω = {omega} accepted nothing.");
    }

    // The standard form's log density about its mode: (λ − 1)·ln(y/m) − ω/2·(y + 1/y − m − 1/m), −∞ off y > 0.
    private readonly struct StandardRatio(double lambda, double omega, double mode) : IRelativeLogDensity
    {
        public double Value(double y) => y > 0.0 ? (lambda - 1.0) * Math.Log(y / mode) - 0.5 * omega * (y - mode + (1.0 / y - 1.0 / mode)) : double.NegativeInfinity;
    }

    /// <summary>One draw conditioned on [<paramref name="lower"/>, <paramref name="upper"/>]: the hat in ln x cut to the interval, or the
    /// limit's own truncated draw.</summary>
    public double SampleOn(Xoshiro256PlusPlus rng, double lower, double upper)
    {
        ArgumentNullException.ThrowIfNull(rng);
        (lower, upper) = Tails.Clip(lower, upper, 0.0, double.PositiveInfinity);
        if (B == 0.0) return Gamma.ShapeRate(P, 0.5 * A).SampleOn(rng, lower, upper);
        if (A == 0.0) return InverseGamma.ShapeScale(-P, 0.5 * B).SampleOn(rng, lower, upper);
        // In u = ln x: p·u − (a·eᵘ + b·e⁻ᵘ)/2 = p·u − ω·cosh(u − ln s) with ω = √(ab), s = √(b/a); its peak is at u* = ln s + asinh(p/ω).
        // The hat's anchor is the peak, or the interval's end nearest it, d from it, as the gamma's is.
        double sa = Math.Sqrt(A), sb = Math.Sqrt(B), omega = sa * sb;
        double peak = Math.Log(sb) - Math.Log(sa) + Math.Asinh(P / omega);
        double lo = lower > 0.0 ? Math.Log(lower) : double.NegativeInfinity, hi = Math.Log(upper);
        double anchor = Math.Clamp(peak, lo, hi), d = anchor - peak;
        // R = √(p² + ω²) = ω·cosh(asinh(p/ω)); A = (R + p)/2 and B = (R − p)/2, the smaller one as ω²/(2(R + |p|)), without cancellation.
        double radius = Math.Sqrt(P * P + omega * omega), small = omega * omega / (2.0 * (radius + Math.Abs(P)));
        double up = P >= 0.0 ? 0.5 * (radius + P) : small, down = P >= 0.0 ? small : 0.5 * (radius - P);
        var h = new LogRatio(up, down, d);
        double leftReach, rightReach;
        if (d == 0.0)
        {
            // From the peak h falls by 1 no later than either term alone does: −A·g(δ) grows as eᵟ to the right and as |δ| to the left, and
            // −B·g(−δ) the other way. The envelope draws from any reach; a reach near that point saves it a search: the curvature's,
            // far past it where one term is tiny, made a draw four to seven times slower there.
            leftReach = Math.Min(ExponentialReach(down), LinearReach(up));
            rightReach = Math.Min(ExponentialReach(up), LinearReach(down));
        }
        else
        {
            double slope = h.Slope(0.0), curvature = -h.Curvature(0.0);
            leftReach = LogConcaveEnvelope.Reach(Math.Max(0.0, slope), curvature);
            rightReach = LogConcaveEnvelope.Reach(Math.Max(0.0, -slope), curvature);
        }
        double delta = LogConcaveEnvelope.Draw(rng, h, anchor - lo, hi - anchor, leftReach, rightReach);
        return Math.Clamp(Math.Exp(anchor + delta), lower, upper);
    }

    // The t > 0 where c·(eᵗ − 1 − t) = 1: √(2/c) for c ≥ 1, where the quadratic rules, else t = ln(1 + 1/c + t) iterated from ln(1/c).
    private static double ExponentialReach(double c)
    {
        if (c >= 1.0) return Math.Sqrt(2.0 / c);
        double t = -Math.Log(c);
        for (int i = 0; i < 4; i++) t = Math.Log(1.0 + 1.0 / c + t);
        return t;
    }

    // The t > 0 where c·(e⁻ᵗ − 1 + t) = 1: √(2/c) for c ≥ 1, else t = 1 + 1/c − e⁻ᵗ iterated from 1 + 1/c.
    private static double LinearReach(double c)
    {
        if (c >= 1.0) return Math.Sqrt(2.0 / c);
        double t = 1.0 + 1.0 / c;
        for (int i = 0; i < 4; i++) t = 1.0 + 1.0 / c - Math.Exp(-t);
        return t;
    }

    // About the peak, L(δ) = p·δ − ω·(cosh(u* + δ) − cosh u*) = −A·g(δ) − B·g(−δ), g(x) = eˣ − 1 − x ≥ 0, A = (R + p)/2, B = (R − p)/2,
    // R = √(p² + ω²): both terms non-positive, so nothing cancels however far the parameters range, and concave. The ratio about the anchor
    // d from the peak is h(δ) = L(d + δ) − L(d) = −A·e^d·(e^δ − 1) − B·e^(−d)·(e^(−δ) − 1) + (A − B)·δ.
    internal readonly struct LogRatio(double up, double down, double d) : IConcaveLogRatio
    {
        public double Value(double delta)
        {
            if (d == 0.0) return -up * ExpLessLinear(delta) - down * ExpLessLinear(-delta);
            return -up * Math.Exp(d) * SpecialFunctions.Expm1(delta) - down * Math.Exp(-d) * SpecialFunctions.Expm1(-delta) + (up - down) * delta;
        }

        public double Slope(double delta) => -up * SpecialFunctions.Expm1(d + delta) + down * SpecialFunctions.Expm1(-(d + delta));

        public double Curvature(double delta) => -up * Math.Exp(d + delta) - down * Math.Exp(-(d + delta));

        // eˣ − 1 − x, by its series where the two would cancel.
        private static double ExpLessLinear(double x)
        {
            if (Math.Abs(x) > 0.1) return SpecialFunctions.Expm1(x) - x;
            double x2 = x * x;
            return x2 * (0.5 + x * (1.0 / 6.0 + x * (1.0 / 24.0 + x * (1.0 / 120.0 + x * (1.0 / 720.0 + x * (1.0 / 5040.0 + x / 40320.0))))));
        }
    }
}
