using System;

namespace Maths.Distributions;

/// <summary>
/// The Poincaré constant of a law with density ρ on a bounded interval [a, b]: the least C with Var_ρ(f) ≤ C·E_ρ[f′²] for every f
/// of finite energy, 1/λ₁ for λ₁ the first nonzero eigenvalue of the weighted Neumann problem −(ρu′)′ = λρu, u′(a) = u′(b) = 0
/// (E8, slice 8).
/// <para>Read by piecewise-linear finite elements on a uniform grid with the consistent mass matrix, each element's integrals of ρ by
/// four-point Gauss–Legendre: a Rayleigh–Ritz method on nested spaces, so the grid's λ₁ʰ lies above λ₁ and falls as the grid is
/// halved (the min–max principle, cited from memory), and the constant read lies under the law's by O(h²). The grid doubles from
/// <see cref="FirstGrid"/> elements until Richardson's estimate of λ₁ʰ's relative error, (λ_h − λ_{h/2})/(3λ_{h/2}), is within the
/// tolerance, or the grid reaches its cap.</para>
/// <para>λ₁ʰ is found by bisection on the inertia of K − λM (Sylvester's law): the count of negative pivots of its LDLᵀ is the count of
/// the pencil's eigenvalues below λ. Each pivot is carried as its departure from the stiffness alone, so the constants' zero eigenvalue
/// does not cancel the pivots away and a small λ₁ keeps its digits.</para>
/// </summary>
/// <param name="Constant">1/λ₁ʰ on the finest grid, under the law's constant by about <paramref name="RelativeError"/> of it; +∞ where
/// the density underflows across a whole element, which splits the interval into pieces a function can hold at different
/// constants.</param>
/// <param name="RelativeError">Richardson's estimate of λ₁ʰ's relative error from the last two grids; 0 for an infinite constant.</param>
/// <param name="Elements">The finest grid's elements.</param>
public sealed record IntervalPoincare(double Constant, double RelativeError, int Elements)
{
    /// <summary>The first grid's elements.</summary>
    public const int FirstGrid = 1024;

    /// <summary>The default cap on the grid's elements.</summary>
    public const int MaxElements = 1 << 20;

    private const double Bracket = 1e-13;   // bisection stops when the bracket's ends are this close, relatively

    /// <summary>
    /// The constant of the law with log density <paramref name="logDensity"/> on [<paramref name="lower"/>, <paramref name="upper"/>],
    /// normalized or not (λ₁ is unchanged by a constant factor), its grid doubled until the estimated relative error is within
    /// <paramref name="tolerance"/> or the grid reaches <paramref name="maxElements"/>.
    /// </summary>
    public static IntervalPoincare Of(Func<double, double> logDensity, double lower, double upper, double tolerance = 1e-7,
                                      int maxElements = MaxElements)
    {
        ArgumentNullException.ThrowIfNull(logDensity);
        Interval(lower, upper);
        if (!(tolerance > 0.0)) throw new ArgumentOutOfRangeException(nameof(tolerance), "The tolerance is positive.");
        if (maxElements < 2 * FirstGrid) throw new ArgumentOutOfRangeException(nameof(maxElements), $"The cap is at least {2 * FirstGrid} elements.");
        int n = FirstGrid;
        double coarse = Gap(logDensity, lower, upper, n, null);
        if (coarse == 0.0) return new IntervalPoincare(double.PositiveInfinity, 0.0, n);
        while (true)
        {
            n *= 2;
            double fine = Gap(logDensity, lower, upper, n, coarse);
            if (fine == 0.0) return new IntervalPoincare(double.PositiveInfinity, 0.0, n);
            double error = Math.Abs(coarse - fine) / (3.0 * fine);
            if (error <= tolerance || n > maxElements / 2) return new IntervalPoincare(1.0 / fine, error, n);
            coarse = fine;
        }
    }

    /// <summary>λ₁ʰ on a uniform grid of <paramref name="elements"/> elements: above λ₁, and 0 where the density underflows across a
    /// whole element.</summary>
    public static double FirstGap(Func<double, double> logDensity, double lower, double upper, int elements)
    {
        ArgumentNullException.ThrowIfNull(logDensity);
        Interval(lower, upper);
        if (elements < 2) throw new ArgumentOutOfRangeException(nameof(elements), "The grid needs at least two elements.");
        return Gap(logDensity, lower, upper, elements, null);
    }

    private static void Interval(double lower, double upper)
    {
        if (!double.IsFinite(lower) || !double.IsFinite(upper) || !(upper > lower))
            throw new ArgumentException("The interval needs finite ends, lower below upper.");
    }

    // λ₁ʰ by bisection, bracketed from a coarser grid's λ₁ʰ when given (the spaces are nested, so the finer one is no larger), else from
    // the linear function's Rayleigh quotient 1/Var_ρ(x), which the grid represents exactly.
    private static double Gap(Func<double, double> logDensity, double lower, double upper, int n, double? coarser)
    {
        Pencil pencil = Pencil.Of(logDensity, lower, upper, n);
        double hi = coarser is { } c ? c * (1.0 + 1e-9) : pencil.Variance > 0.0 ? 1.0 / pencil.Variance : 1.0;
        for (int i = 0; pencil.Below(hi) < 2; i++)
        {
            if (i == 64 || !double.IsFinite(hi)) throw new InvalidOperationException("The pencil's second eigenvalue was not bracketed.");
            hi *= 2.0;
        }
        double lo = 0.5 * hi;
        while (pencil.Below(lo) >= 2)
        {
            lo *= 0.5;
            if (lo < 1e-300) return 0.0;
        }
        while (hi > lo * (1.0 + Bracket))
        {
            double mid = Math.Sqrt(lo * hi);
            if (pencil.Below(mid) >= 2) hi = mid;
            else lo = mid;
        }
        return 0.5 * (lo + hi);
    }

    // The stiffness K and consistent mass M on n uniform elements, scaled by the density's largest value at the quadrature points: per
    // element e, k_e = ∫_e ρ / h², and the mass entries ∫_e ρφ_lφ_l (α), ∫_e ρφ_rφ_r (β) and ∫_e ρφ_lφ_r (γ), φ the hat functions.
    private sealed class Pencil
    {
        private readonly double[] _k, _alpha, _beta, _gamma;

        private Pencil(double[] k, double[] alpha, double[] beta, double[] gamma, double variance)
        {
            _k = k;
            _alpha = alpha;
            _beta = beta;
            _gamma = gamma;
            Variance = variance;
        }

        // Var_ρ(x), from the same element integrals.
        public double Variance { get; }

        public static Pencil Of(Func<double, double> logDensity, double lower, double upper, int n)
        {
            GaussLegendre rule = GaussLegendre.Of(4);
            int q = rule.Count;
            double h = (upper - lower) / n, half = 0.5 * h;
            var top = new double[n];
            var s0 = new double[n];
            var s1 = new double[n];
            var s2 = new double[n];
            Span<double> logs = stackalloc double[q];
            double global = double.NegativeInfinity;
            for (int e = 0; e < n; e++)
            {
                double centre = lower + (e + 0.5) * h, local = double.NegativeInfinity;
                for (int r = 0; r < q; r++)
                {
                    double l = logDensity(centre + half * rule.Nodes[r]);
                    if (double.IsNaN(l) || double.IsPositiveInfinity(l)) throw new ArgumentException("The log density is finite or −∞ on the interval.", nameof(logDensity));
                    logs[r] = l;
                    local = Math.Max(local, l);
                }
                top[e] = local;
                global = Math.Max(global, local);
                if (double.IsNegativeInfinity(local)) continue;
                double m0 = 0.0, m1 = 0.0, m2 = 0.0;
                for (int r = 0; r < q; r++)
                {
                    double t = half * rule.Nodes[r], w = half * rule.Weights[r] * Math.Exp(logs[r] - local);
                    m0 += w;
                    m1 += w * t;
                    m2 += w * t * t;
                }
                s0[e] = m0;
                s1[e] = m1;
                s2[e] = m2;
            }
            if (double.IsNegativeInfinity(global)) throw new ArgumentException("The density is zero across the interval.", nameof(logDensity));
            // The element integrals scaled to one largest density, written over the sums: k into top, α, β, γ into s0, s1, s2. The
            // moments about the interval's midpoint give Var_ρ(x), which only brackets the bisection.
            double mass = 0.0, first = 0.0, second = 0.0, middle = 0.5 * (lower + upper);
            for (int e = 0; e < n; e++)
            {
                double scale = double.IsNegativeInfinity(top[e]) ? 0.0 : Math.Exp(top[e] - global);
                double m0 = s0[e] * scale, m1 = s1[e] * scale, m2 = s2[e] * scale, d = lower + (e + 0.5) * h - middle;
                mass += m0;
                first += m1 + d * m0;
                second += m2 + 2.0 * d * m1 + d * d * m0;
                top[e] = m0 / (h * h);
                s0[e] = 0.25 * m0 - m1 / h + m2 / (h * h);
                s1[e] = 0.25 * m0 + m1 / h + m2 / (h * h);
                s2[e] = 0.25 * m0 - m2 / (h * h);
            }
            double mean = first / mass;
            return new Pencil(top, s0, s1, s2, second / mass - mean * mean);
        }

        // The count of the pencil's eigenvalues below λ: the negative pivots of K − λM's LDLᵀ. With q_i = k_i + r_i (k_n = 0 at the last
        // node), r_0 = −λα_0 and r_i = −λ(β_{i−1} + α_i) + [k_{i−1}(r_{i−1} − 2λγ_{i−1}) − λ²γ_{i−1}²]/q_{i−1}; a zero pivot is taken as the
        // smallest positive one.
        public int Below(double lambda)
        {
            int n = _k.Length, count = 0;
            double r = -lambda * _alpha[0], q = _k[0] + r;
            if (q < 0.0) count++;
            for (int i = 1; i <= n; i++)
            {
                double kPrev = _k[i - 1], g = _gamma[i - 1];
                double pivot = q != 0.0 ? q : double.Epsilon;
                r = -lambda * (_beta[i - 1] + (i < n ? _alpha[i] : 0.0)) + (kPrev * (r - 2.0 * lambda * g) - lambda * lambda * g * g) / pivot;
                q = (i < n ? _k[i] : 0.0) + r;
                if (q < 0.0) count++;
            }
            return count;
        }
    }
}
