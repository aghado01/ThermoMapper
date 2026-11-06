using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The generalized hyperbolic as a normal variance-mean mixture (Barndorff-Nielsen 1977): x = m·v + s·√v·ζ with ζ ~ N(0, 1) and
/// v ~ GIG(p, a, b) (<see cref="GeneralizedInverseGaussian"/>). The Normal-Wald, or normal inverse Gaussian, is p = −½, b = 1; the skewed t
/// a = 0; the variance gamma b = 0.
/// <para><b>Density.</b> Closed: with L(p, a, b) the GIG's log normalizer (<see cref="GeneralizedInverseGaussian.LogNormalizer"/>),
/// ln f(x) = m·x/s² − ½ln(2πs²) + L(p − ½, a + m²/s², b + x²/s²) − L(p, a, b).</para>
/// <para><b>Tails.</b> The CDF has no closed form: F(x) = ∫ Φ(z(u))·g(u) du over u = ln v, with z(u) = (x − m·eᵘ)/(s·e^(u/2)) and g the GIG's
/// density in u, log-concave there. The integrand is taken in logs: a scan out from the GIG's mode finds its largest value, refined by
/// golden section, and the span beyond which the GIG's density alone bounds it 60 below that; panels widen geometrically from the peak,
/// each integrated by adaptive Gauss–Kronrod (7–15), so the step of Φ, sharp where the drift is strong and v large, is resolved. Each tail
/// is integrated on its own side, Φ(−z) for the survival, so a far tail keeps its digits as a log (<see cref="IntegratedLogCdf"/>,
/// <see cref="IntegratedLogSurvival"/>). Each read costs hundreds of normal tails, so the tails a reader takes come from a table built on the
/// first read: the normal score z(x) = Φ⁻¹(F(x)) = −Φ⁻¹(S(x)), near linear in the body and smooth in every tail, at points walked out
/// from the median to the 10⁻¹² quantiles on each side, each step moving z by at most ¼ by its exact slope f(x)/φ(z), so a power tail
/// takes steps growing with x; joined by cubic Hermite, the tails read back as ln Φ(±z); beyond the table, the integral itself. Against the
/// integral the table holds the log tails to 10⁻⁵ relative in the body and 3·10⁻⁴ out at a power tail's 10⁻⁹, its quantiles' levels to
/// 4·10⁻⁵, within what a predictive's Monte Carlo over draws resolves.</para>
/// <para><b>Quantiles and draws.</b> Quantiles by Newton on the log tails, kept in a bracket; a draw is v then ζ, a draw on an interval the
/// quantile of a uniform between its ends' masses.</para>
/// </summary>
public sealed class GeneralizedHyperbolic : IContinuousDistribution
{
    private const double Drop = 60.0;
    private const int MaxDepth = 40;
    private const double TableStep = 0.25, TableEnd = -27.631021115928547;   // ln 10⁻¹²
    private readonly double _logNormalizer, _mode;
    private readonly GeneralizedInverseGaussian _mixing;
    private TailTable? _table;   // built on the first read of a tail

    private GeneralizedHyperbolic(GeneralizedInverseGaussian mixing, double drift, double scale)
    {
        _mixing = mixing;
        Drift = drift;
        Scale = scale;
        _logNormalizer = GeneralizedInverseGaussian.LogNormalizer(mixing.P, mixing.A, mixing.B);
        _mode = ModeInLog(mixing.P, mixing.A, mixing.B);
    }

    /// <summary>x = <paramref name="drift"/>·v + <paramref name="scale"/>·√v·ζ, v ~ GIG(<paramref name="p"/>, <paramref name="a"/>,
    /// <paramref name="b"/>), a proper GIG; the drift finite, the scale positive and finite.</summary>
    public static GeneralizedHyperbolic Of(double p, double a, double b, double drift, double scale)
    {
        GeneralizedInverseGaussian mixing = GeneralizedInverseGaussian.Of(p, a, b);
        if (!double.IsFinite(drift)) throw new ArgumentOutOfRangeException(nameof(drift), "The drift must be finite.");
        if (!(scale > 0.0 && double.IsFinite(scale))) throw new ArgumentOutOfRangeException(nameof(scale), "The scale must be positive and finite.");
        return new GeneralizedHyperbolic(mixing, drift, scale);
    }

    /// <summary>v's law.</summary>
    public GeneralizedInverseGaussian Mixing => _mixing;

    /// <summary>m, the drift per unit v.</summary>
    public double Drift { get; }

    /// <summary>s, the scale per unit √v.</summary>
    public double Scale { get; }

    /// <inheritdoc/>
    public double LogPdf(double x)
    {
        double s2 = Scale * Scale, m = Drift;
        return m * x / s2 - 0.5 * Math.Log(2.0 * Math.PI * s2)
               + GeneralizedInverseGaussian.LogNormalizer(_mixing.P - 0.5, _mixing.A + m * m / s2, _mixing.B + x * x / s2) - _logNormalizer;
    }

    /// <inheritdoc/>
    public double Cdf(double x) => Math.Exp(LogCdf(x));

    /// <inheritdoc/>
    public double Survival(double x) => Math.Exp(LogSurvival(x));

    /// <inheritdoc/>
    public double LogCdf(double x) => Table.Read(x, below: true);

    /// <inheritdoc/>
    public double LogSurvival(double x) => Table.Read(x, below: false);

    /// <summary>ln P(X ≤ x) by the integral over v itself, not the table.</summary>
    public double IntegratedLogCdf(double x) => LogTail(x, below: true);

    /// <summary>ln P(X &gt; x) by the integral over v itself, not the table.</summary>
    public double IntegratedLogSurvival(double x) => LogTail(x, below: false);

    private TailTable Table => _table ??= new TailTable(this);

    /// <inheritdoc/>
    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    /// <inheritdoc/>
    public double LogIntervalProbability(double lower, double upper)
    {
        if (!(upper > lower)) return double.NegativeInfinity;
        // The difference on the side whose tails are small, so neither is a complement near 1.
        double fl = LogCdf(lower), fu = LogCdf(upper);
        if (fu <= -SpecialFunctions.Ln2) return fu + SpecialFunctions.Log1mExp(fl - fu);
        double sl = LogSurvival(lower), su = LogSurvival(upper);
        if (sl <= -SpecialFunctions.Ln2) return sl + SpecialFunctions.Log1mExp(su - sl);
        return SpecialFunctions.Log1p(-Math.Exp(fl) - Math.Exp(su));
    }

    /// <inheritdoc/>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double v = _mixing.Sample(rng);
        return Drift * v + Scale * Math.Sqrt(v) * Normal.SampleStandard(rng);
    }

    /// <inheritdoc/>
    public double SampleOn(Xoshiro256PlusPlus rng, double lower, double upper)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!(upper > lower)) throw new ArgumentException("The interval must be non-empty.", nameof(upper));
        double fl = LogCdf(lower), fu = LogCdf(upper), w = rng.NextDouble();
        if (fu <= -SpecialFunctions.Ln2)
        {
            // ln(F(l) + w·(F(u) − F(l))) = ln F(u) + ln(e^(ln F(l) − ln F(u)) + w·(1 − that))
            double ratio = Math.Exp(fl - fu);
            return Math.Clamp(InverseLogCdf(fu + Math.Log(ratio + w * (1.0 - ratio))), lower, upper);
        }
        double sl = LogSurvival(lower), su = LogSurvival(upper);
        double back = Math.Exp(su - sl);
        return Math.Clamp(InverseLogSurvival(sl + Math.Log(back + w * (1.0 - back))), lower, upper);
    }

    /// <inheritdoc/>
    public double Quantile(double p)
    {
        Tails.CheckProbability(p);
        if (p == 0.0) return double.NegativeInfinity;
        if (p == 1.0) return double.PositiveInfinity;
        return p <= 0.5 ? InverseLogCdf(Math.Log(p)) : InverseLogSurvival(SpecialFunctions.Log1p(-p));
    }

    /// <inheritdoc/>
    public double InverseLogCdf(double logP)
    {
        if (!(logP <= 0.0)) throw new ArgumentOutOfRangeException(nameof(logP), "A log probability is not positive.");
        return Invert(logP, below: true, tabled: true);
    }

    /// <inheritdoc/>
    public double InverseLogSurvival(double logS)
    {
        if (!(logS <= 0.0)) throw new ArgumentOutOfRangeException(nameof(logS), "A log probability is not positive.");
        return Invert(logS, below: false, tabled: true);
    }

    // The x whose tail on that side has log mass `target`, the table's or the integral's: a bracket grown from the bulk, then Newton on
    // the log tail, whose slope is ±f/F from the closed density, kept inside the bracket by bisection, to 10⁻¹⁴ relative.
    private double Invert(double target, bool below, bool tabled)
    {
        if (target == 0.0) return below ? double.PositiveInfinity : double.NegativeInfinity;
        double v = Math.Exp(_mode), centre = Drift * v, width = Scale * Math.Sqrt(v) + Math.Abs(Drift) * v;
        double Tail(double x) => tabled ? Table.Read(x, below) : LogTail(x, below);
        // The tail below rises with x, the tail above falls: step out until the target is bracketed.
        double lo = centre - width, hi = centre + width;
        for (int k = 0; k < 2000 && (below ? Tail(lo) > target : Tail(lo) < target); k++) lo = centre - (centre - lo) * 2.0;
        for (int k = 0; k < 2000 && (below ? Tail(hi) < target : Tail(hi) > target); k++) hi = centre + (hi - centre) * 2.0;
        double x = Math.Clamp(centre, lo, hi);
        for (int k = 0; k < 200; k++)
        {
            double tail = Tail(x);
            if (tail == target) return x;
            if (below ? tail < target : tail > target) lo = x;
            else hi = x;
            double slope = Math.Exp(LogPdf(x) - tail) * (below ? 1.0 : -1.0);
            double next = x - (tail - target) / slope;
            if (!(next > lo && next < hi)) next = lo + 0.5 * (hi - lo);
            if (Math.Abs(next - x) <= 1e-14 * Math.Max(1.0, Math.Abs(x)) || next <= lo || next >= hi) return next;
            x = next;
        }
        return x;
    }

    // ln P(X ≤ x) (below) or ln P(X > x): ln ∫ exp(h(u)) du with h(u) = ln Φ(±z(u)) + ln g(u).
    private double LogTail(double x, bool below)
    {
        if (double.IsNaN(x)) return double.NaN;
        if (double.IsPositiveInfinity(x)) return below ? 0.0 : double.NegativeInfinity;
        if (double.IsNegativeInfinity(x)) return below ? double.NegativeInfinity : 0.0;
        double p = _mixing.P, a = _mixing.A, b = _mixing.B, s = Scale, m = Drift, normalizer = _logNormalizer;
        double LogMixing(double u) => p * u - 0.5 * ((a > 0.0 ? a * Math.Exp(u) : 0.0) + (b > 0.0 ? b * Math.Exp(-u) : 0.0)) - normalizer;
        double H(double u)
        {
            // Each term only where it is not zero, so an infinite eᵘ or e⁻ᵘ far out does not meet a zero drift or x as 0·∞.
            double z = ((x != 0.0 ? x * Math.Exp(-0.5 * u) : 0.0) - (m != 0.0 ? m * Math.Exp(0.5 * u) : 0.0)) / s;
            return Normal.StandardLogCdf(below ? z : -z) + LogMixing(u);
        }

        // The scan out from the GIG's mode in u: the best point seen, and on each side the span until g alone is Drop below it.
        double best = _mode, top = H(best);
        double Scan(double direction)
        {
            double step = 0.25, at = _mode;
            for (int k = 0; k < 10_000; k++)
            {
                at += direction * step;
                double h = H(at);
                if (h > top) { top = h; best = at; }
                if (LogMixing(at) < top - Drop && LogMixing(at) < LogMixing(at - direction * step)) return at;
                if (k >= 8) step *= 1.25;
            }
            return at;
        }
        double right = Scan(1.0), left = Scan(-1.0);
        if (double.IsNegativeInfinity(top)) return double.NegativeInfinity;

        // Golden section about the best point for the peak, then its width from the curvature.
        double lo = Math.Max(left, best - 0.5), hi = Math.Min(right, best + 0.5);
        const double golden = 0.6180339887498949;
        for (int k = 0; k < 25; k++)
        {
            double c = hi - golden * (hi - lo), d = lo + golden * (hi - lo);
            if (H(c) > H(d)) hi = d;
            else lo = c;
        }
        double peak = 0.5 * (lo + hi);
        top = Math.Max(top, H(peak));
        double e = 1e-4, curvature = -(H(peak + e) - 2.0 * H(peak) + H(peak - e)) / (e * e);
        double width = curvature > 1.0 ? 1.0 / Math.Sqrt(curvature) : 1.0;
        width = Math.Clamp(width, 1e-6, 1.0);

        double sum = 0.0;
        foreach (double direction in new[] { 1.0, -1.0 })
        {
            double end = direction > 0 ? right : left, at = peak, step = width;
            while (direction > 0 ? at < end : at > end)
            {
                double next = direction > 0 ? Math.Min(at + step, end) : Math.Max(at - step, end);
                double panel = Kronrod(H, Math.Min(at, next), Math.Max(at, next), top, 0);
                sum += panel;
                at = next;
                step *= 2.0;
                // Past the peak a panel that adds nothing, ending where h is far below it, ends the walk: what is beyond is smaller still.
                if (panel <= 1e-17 * sum && H(at) < top - 45.0 && LogMixing(at) < LogMixing(at - direction * 1e-3)) break;
            }
        }
        return top + Math.Log(sum);
    }

    // ∫ exp(h − top) over [lo, hi] by Gauss–Kronrod 7–15, halving a panel until the two agree to 10⁻¹¹ of the panel at the peak's height.
    private static double Kronrod(Func<double, double> h, double lo, double hi, double top, int depth)
    {
        double centre = 0.5 * (lo + hi), half = 0.5 * (hi - lo), kronrod = 0.0, gauss = 0.0;
        for (int i = 0; i < KronrodNodes.Length; i++)
        {
            double f = KronrodNodes[i] == 0.0
                ? Math.Exp(h(centre) - top)
                : Math.Exp(h(centre - half * KronrodNodes[i]) - top) + Math.Exp(h(centre + half * KronrodNodes[i]) - top);
            kronrod += KronrodWeights[i] * f;
            gauss += GaussWeights[i] * f;
        }
        kronrod *= half;
        gauss *= half;
        if (!(Math.Abs(kronrod - gauss) > 1e-11 * 2.0 * half) || depth >= MaxDepth) return kronrod;
        return Kronrod(h, lo, centre, top, depth + 1) + Kronrod(h, centre, hi, top, depth + 1);
    }

    // The normal score z(x) = Φ⁻¹(F(x)) at points walked out from the median, with its slopes f/φ(z); cubic Hermite between them, ln Φ(±z)
    // for the tails, and the integral beyond the ends.
    private sealed class TailTable
    {
        private static readonly Normal Standard = Normal.MeanSd(0.0, 1.0);
        private readonly GeneralizedHyperbolic _owner;
        private readonly double[] _x, _z, _slope;

        public TailTable(GeneralizedHyperbolic owner)
        {
            _owner = owner;
            double median = owner.Invert(-SpecialFunctions.Ln2, below: true, tabled: false);
            var x = new System.Collections.Generic.List<double>();
            var z = new System.Collections.Generic.List<double>();
            var slope = new System.Collections.Generic.List<double>();
            foreach (bool below in new[] { true, false })
            {
                var side = new System.Collections.Generic.List<(double X, double Z, double Slope)>();
                double at = median;
                for (int k = 0; k < 10_000; k++)
                {
                    double tail = owner.LogTail(at, below);
                    double score = below ? Standard.InverseLogCdf(tail) : -Standard.InverseLogCdf(tail);
                    double rate = Math.Exp(owner.LogPdf(at) + 0.5 * score * score + 0.5 * Math.Log(2.0 * Math.PI));
                    side.Add((at, score, rate));
                    if (tail < TableEnd || !(rate > 0.0) || !double.IsFinite(rate)) break;
                    at += (below ? -1.0 : 1.0) * TableStep / rate;
                }
                if (below) side.Reverse();
                foreach ((double X, double Z, double Slope) point in below ? side : side.GetRange(1, side.Count - 1))
                {
                    x.Add(point.X);
                    z.Add(point.Z);
                    slope.Add(point.Slope);
                }
            }
            (_x, _z, _slope) = ([.. x], [.. z], [.. slope]);
        }

        public double Read(double x, bool below)
        {
            if (double.IsNaN(x)) return double.NaN;
            if (x < _x[0] || x > _x[^1]) return _owner.LogTail(x, below);
            int k = Array.BinarySearch(_x, x);
            double score;
            if (k >= 0) score = _z[k];
            else
            {
                k = Math.Clamp(~k - 1, 0, _x.Length - 2);
                double h = _x[k + 1] - _x[k], t = (x - _x[k]) / h, t2 = t * t, t3 = t2 * t;
                score = (2.0 * t3 - 3.0 * t2 + 1.0) * _z[k] + (t3 - 2.0 * t2 + t) * h * _slope[k]
                        + (-2.0 * t3 + 3.0 * t2) * _z[k + 1] + (t3 - t2) * h * _slope[k + 1];
            }
            return Normal.StandardLogCdf(below ? score : -score);
        }
    }

    // The GIG's mode in u = ln v, where p = (a·eᵘ − b·e⁻ᵘ)/2.
    private static double ModeInLog(double p, double a, double b)
    {
        if (a == 0.0) return Math.Log(-0.5 * b / p);
        if (b == 0.0) return Math.Log(2.0 * p / a);
        double root = Math.Sqrt(p * p + a * b);
        return Math.Log(p >= 0.0 ? (p + root) / a : b / (root - p));
    }

    // The 15 Kronrod nodes on [0, 1] (the centre last) and their weights; the 7-point Gauss weights at the shared nodes, 0 elsewhere.
    private static readonly double[] KronrodNodes =
    [
        0.991455371120812639, 0.949107912342758525, 0.864864423359769073, 0.741531185599394440,
        0.586087235467691130, 0.405845151377397167, 0.207784955007898468, 0.0,
    ];

    private static readonly double[] KronrodWeights =
    [
        0.022935322010529225, 0.063092092629978553, 0.104790010322250184, 0.140653259715525919,
        0.169004726639267903, 0.190350578064785410, 0.204432940075298892, 0.209482141084727828,
    ];

    private static readonly double[] GaussWeights =
    [
        0.0, 0.129484966168869693, 0.0, 0.279705391489276668, 0.0, 0.381830050505118945, 0.0, 0.417959183673469388,
    ];
}
