using System;

namespace Maths.Distributions;

public static partial class SpecialFunctions
{
    // At and above this order the uniform expansion's first four terms hold to about 10⁻¹⁴; below it the recurrence's steps cost less.
    private const double DebyeOrder = 300.0;

    /// <summary>ln K_ν(x), the modified Bessel function of the second kind, for real ν and x ≥ 0: +∞ at x = 0, −∞ at x = +∞. In logs, so
    /// it holds far past where K itself over- or underflows.</summary>
    public static double LogBesselK(double nu, double x) => LogBesselKScaled(nu, x) - x;

    /// <summary>
    /// ln(eˣ·K_ν(x)), R's <c>besselK(x, ν, expon.scaled = TRUE)</c> in logs, for real ν and x ≥ 0. K is even in ν, so |ν| is read. For the
    /// fractional order μ = ν − n, |μ| ≤ ½: Temme's series for x ≤ 2 and Steed's continued fraction CF2 beyond (Temme 1975; Thompson and
    /// Barnett 1986; Numerical Recipes §6.7), each giving K_μ and K_{μ+1}; then the upward recurrence to ν, which is stable for K, carried in
    /// the scaled ratio sₖ = x·K_{μ+k+1}/K_{μ+k} = 2(μ + k) + x²/sₖ₋₁, which neither overflows nor underflows at any x. At ν ≥ 300 the
    /// uniform expansion in the order (Debye; DLMF 10.41.4) to four terms, ν·η written as ν(asinh(1/z) − 1/(z + √(1 + z²))) so neither
    /// end of z = x/ν cancels.
    /// </summary>
    public static double LogBesselKScaled(double nu, double x)
    {
        if (double.IsNaN(nu) || double.IsNaN(x)) return double.NaN;
        if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x), "K_ν(x) is real for x ≥ 0.");
        nu = Math.Abs(nu);
        if (x == 0.0 || double.IsPositiveInfinity(nu)) return double.PositiveInfinity;
        if (double.IsPositiveInfinity(x)) return double.NegativeInfinity;
        if (nu >= DebyeOrder) return DebyeScaled(nu, x);
        int n = (int)Math.Floor(nu + 0.5);
        double mu = nu - n;
        (double logK, double s) = x <= 2.0 ? Temme(mu, x) : Steed(mu, x);
        double logX = Math.Log(x);
        for (int k = 1; k <= n; k++)
        {
            logK += Math.Log(s) - logX;
            s = 2.0 * (mu + k) + x * (x / s);
        }
        return logK;
    }

    // Temme's series, |μ| ≤ ½ and x ≤ 2: ln(eˣK_μ) and s₀ = x·K_{μ+1}/K_μ, from K_μ = Σ cₖfₖ and K_{μ+1} = (2/x)Σ cₖ(pₖ − k·fₖ).
    private static (double LogScaled, double Ratio) Temme(double mu, double x)
    {
        double half = 0.5 * x, pimu = Math.PI * mu;
        double fact = Math.Abs(pimu) < Eps ? 1.0 : pimu / Math.Sin(pimu);
        double d = -Math.Log(half), e = mu * d;
        double fact2 = Math.Abs(e) < Eps ? 1.0 : Math.Sinh(e) / e;
        (double gam1, double gam2, double plus, double minus) = TemmeGammas(mu);
        double ff = fact * (gam1 * Math.Cosh(e) + gam2 * fact2 * d);
        double sum = ff, expE = Math.Exp(e);
        double p = 0.5 * expE / plus;      // ½(x/2)^(−μ)Γ(1 + μ)
        double q = 0.5 / (expE * minus);   // ½(x/2)^μ Γ(1 − μ)
        double c = 1.0, square = half * half, sum1 = p;
        for (int i = 1; i < MaxIterations; i++)
        {
            ff = (i * ff + p + q) / ((double)i * i - mu * mu);
            c *= square / i;
            p /= i - mu;
            q /= i + mu;
            double term = c * ff;
            sum += term;
            sum1 += c * (p - i * ff);
            if (Math.Abs(term) < Math.Abs(sum) * Eps) break;
        }
        return (Math.Log(sum) + x, 2.0 * sum1 / sum);
    }

    // 1/Γ(1 + μ), 1/Γ(1 − μ), their half-sum and their half-difference over μ, which near μ = 0 is its series −Σ d₂ⱼ₊₁μ²ʲ in the
    // coefficients dₖ of 1/Γ(1 + z) = Σ dₖzᵏ (Abramowitz and Stegun 6.1.34), where the difference would cancel.
    private static (double Gam1, double Gam2, double Plus, double Minus) TemmeGammas(double mu)
    {
        double plus = Math.Exp(-LogGamma(1.0 + mu)), minus = Math.Exp(-LogGamma(1.0 - mu));
        double gam1;
        if (Math.Abs(mu) > 0.05) gam1 = (minus - plus) / (2.0 * mu);
        else
        {
            double m2 = mu * mu;
            gam1 = -(0.5772156649015329 + m2 * (-0.0420026350340952 + m2 * (-0.0421977345555443 + m2 * (0.0072189432466630
                     + m2 * -0.0002152416741149))));
        }
        return (gam1, 0.5 * (minus + plus), plus, minus);
    }

    // Steed's method for CF2, |μ| ≤ ½ and x > 2: ln(eˣK_μ) = ½ln(π/(2x)) − ln s and s₀ = x·K_{μ+1}/K_μ = μ + x + ½ − h.
    private static (double LogScaled, double Ratio) Steed(double mu, double x)
    {
        double b = 2.0 * (1.0 + x), d = 1.0 / b, h = d, delh = d, q1 = 0.0, q2 = 1.0;
        double a1 = 0.25 - mu * mu, q = a1, c = a1, a = -a1, s = 1.0 + q * delh;
        for (int i = 2; i < MaxIterations; i++)
        {
            a -= 2 * (i - 1);
            c = -a * c / i;
            double next = (q1 - b * q2) / a;
            q1 = q2;
            q2 = next;
            q += c * next;
            b += 2.0;
            d = 1.0 / (b + a * d);
            delh = (b * d - 1.0) * delh;
            h += delh;
            double step = q * delh;
            s += step;
            if (Math.Abs(step / s) < Eps) break;
        }
        return (0.5 * Math.Log(Math.PI / (2.0 * x)) - Math.Log(s), mu + x + 0.5 - a1 * h);
    }

    // ln(eˣK_ν(x)) by the uniform expansion at z = x/ν, t = 1/√(1 + z²): x − νη − ½ln(2ν/π) − ¼ln(1 + z²) + ln Σₖ(−1)ᵏuₖ(t)/νᵏ.
    private static double DebyeScaled(double nu, double x)
    {
        double z = x / nu;
        double root = z > 1e8 ? z * Math.Sqrt(1.0 + 1.0 / (z * z)) : Math.Sqrt(1.0 + z * z);
        double asinhInverse = z >= 1.0 ? Math.Asinh(1.0 / z) : Math.Log(1.0 + root) - Math.Log(z);
        double exponent = nu * (asinhInverse - 1.0 / (z + root));
        double t = 1.0 / root, t2 = t * t;
        double u1 = t * (3.0 - 5.0 * t2) / 24.0;
        double u2 = t2 * (81.0 + t2 * (-462.0 + t2 * 385.0)) / 1152.0;
        double u3 = t * t2 * (30375.0 + t2 * (-369603.0 + t2 * (765765.0 - t2 * 425425.0))) / 414720.0;
        double u4 = t2 * t2 * (4465125.0 + t2 * (-94121676.0 + t2 * (349922430.0 + t2 * (-446185740.0 + t2 * 185910725.0)))) / 39813120.0;
        double inverse = 1.0 / nu;
        double series = 1.0 + inverse * (-u1 + inverse * (u2 + inverse * (-u3 + inverse * u4)));
        return exponent + 0.5 * Math.Log(Math.PI / (2.0 * nu)) - 0.5 * Math.Log(root) + Math.Log(series);
    }
}
