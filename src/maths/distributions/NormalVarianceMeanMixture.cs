using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// A normal variance-mean mixture (Barndorff-Nielsen 1977): x = m·v + s·√v·ζ with ζ ~ N(0, 1) and v from a mixing law whose density g in
/// u = ln v is log-concave, the generalized hyperbolic's GIG (<see cref="GeneralizedHyperbolic"/>) and the generalized beta prime
/// (<see cref="NormalGeneralizedBetaPrime"/>). A family states g, its mode in u, its power at v → 0 and its draw of v; the integrals are here.
/// <para><b>Density.</b> f(x) = ∫ φ(z(u))·g(u)/(s·e^(u/2)) du with z(u) = (x − m·eᵘ)/(s·e^(u/2)), integrated as the tails are (below) where a
/// family has no closed form. Its log integrand is concave in u, since z(u)² = (x²e⁻ᵘ − 2xm + m²eᵘ)/s² is convex, so a scan may stop where
/// the integrand itself has fallen far enough and a golden section over the scan's span finds its one peak; where that peak is under 0.1 wide,
/// far out on the drift's side, Newton finds it. At x = 0 the density is +∞ where v's density near 0 falls no faster than v^(−½)
/// (<see cref="LeftExponent"/> ≤ ½), the mixture's pole, which the tails and their table step around.</para>
/// <para><b>Tails.</b> The CDF has no closed form: F(x) = ∫ Φ(z(u))·g(u) du. The integrand is taken in logs: a scan out from g's mode finds
/// its largest value, refined by golden section, and the span beyond which g alone bounds it 60 below that; panels widen geometrically from
/// the peak, each integrated by adaptive Gauss–Kronrod (7–15), so the step of Φ, sharp where the drift is strong and v large, is resolved.
/// Each tail is integrated on its own side, Φ(−z) for the survival, so a far tail keeps its digits as a log (<see cref="IntegratedLogCdf"/>,
/// <see cref="IntegratedLogSurvival"/>). Where x and m share a sign the integrals are taken in t = u − ln(x/m), where z = −2·sign(x)·√(xm)·
/// sinh(t/2)/s: z's two terms cancel near that point, and far out its step narrows below u's rounding.</para>
/// <para><b>The table.</b> Each read costs hundreds of normal tails, so the tails a reader takes come from a table built on the first read:
/// the normal score z(x) = Φ⁻¹(F(x)) = −Φ⁻¹(S(x)), near linear in the body and smooth in every tail, at points walked out from the median to
/// the 10⁻¹² quantiles on each side, each step moving z by at most ¼ by its exact slope f(x)/φ(z), so a power tail takes steps growing with
/// x; joined by cubic Hermite in asinh(x/s), linear in x near 0 and ln|x| far out, where a power tail's score goes as √ln|x| (in x a cell
/// there spans a multiple of x, which a cubic in x follows poorly); the tails read back as ln Φ(±z); beyond the table, the integral itself.
/// Near a pole the score goes as x·ln|x|, its slope infinite at 0: the walk's steps there are at most ¼ of |x| (½ within 10⁻²·s) down to
/// 10⁻⁶·s, then the pole itself, and a cell with the pole at an end is read by the integral. Against the integral the table holds the log
/// tails to 10⁻⁵ relative in the body and 3·10⁻⁴ out at a power tail's 10⁻⁹, its quantiles' levels to 4·10⁻⁵, within what a predictive's
/// Monte Carlo over draws resolves.</para>
/// <para><b>Quantiles and draws.</b> Quantiles by Newton on the log tails, kept in a bracket; a draw is v then ζ, a draw on an interval the
/// quantile of a uniform between its ends' masses.</para>
/// </summary>
public abstract class NormalVarianceMeanMixture : IContinuousDistribution
{
    private const double Drop = 60.0;
    private const int MaxDepth = 40;
    private const double TableStep = 0.25, TableEnd = -27.631021115928547;   // ln 10⁻¹²
    // The table's cells near a pole: at most ¼ of |x| long, ½ within 10⁻²·s, down to 10⁻⁶·s.
    private const double PoleRatio = 0.25, PoleNearRatio = 0.5, PoleNear = 1e-2, PoleFloor = 1e-6;
    private readonly double _mode;
    private TailTable? _table;   // built on the first read of a tail

    private protected NormalVarianceMeanMixture(double drift, double scale, double modeInLog)
    {
        Drift = drift;
        Scale = scale;
        _mode = modeInLog;
    }

    /// <summary>m, the drift per unit v.</summary>
    public double Drift { get; }

    /// <summary>s, the scale per unit √v.</summary>
    public double Scale { get; }

    /// <summary>ln g(u), the mixing law's normalized log density in u = ln v, log-concave.</summary>
    private protected abstract double LogMixing(double u);

    /// <summary>A draw of v.</summary>
    private protected abstract double SampleMixing(Xoshiro256PlusPlus rng);

    /// <summary>The rate at which ln g(u) rises from u = −∞, v's density near 0 going as v to that power less 1; +∞ where it falls
    /// faster than any power.</summary>
    private protected abstract double LeftExponent { get; }

    // Whether the density is +∞ at x = 0.
    private bool HasPole => LeftExponent <= 0.5;

    /// <summary>The density by the integral over v (<see cref="IntegratedLogPdf"/>), where a family has no closed form.</summary>
    public virtual double LogPdf(double x) => IntegratedLogPdf(x);

    /// <summary>ln f(x) by the integral over v itself; +∞ at a pole.</summary>
    public double IntegratedLogPdf(double x)
    {
        if (double.IsNaN(x)) return double.NaN;
        if (double.IsInfinity(x)) return double.NegativeInfinity;
        if (x == 0.0 && HasPole) return double.PositiveInfinity;
        return LogIntegral(x, Integrand.Density);
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
        double v = SampleMixing(rng);
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
    // the log tail, whose slope is ±f/F from the density, kept inside the bracket by bisection, to 10⁻¹⁴ relative. At a pole the slope is
    // infinite and the step stays put, on the bracket's end, so bisection takes it.
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

    private enum Integrand { Below, Above, Density }

    // ln P(X ≤ x) (below) or ln P(X > x).
    private double LogTail(double x, bool below)
    {
        if (double.IsNaN(x)) return double.NaN;
        if (double.IsPositiveInfinity(x)) return below ? 0.0 : double.NegativeInfinity;
        if (double.IsNegativeInfinity(x)) return below ? double.NegativeInfinity : 0.0;
        return LogIntegral(x, below ? Integrand.Below : Integrand.Above);
    }

    // ln ∫ exp(h(u)) du with h(u) = ln Φ(±z(u)) + ln g(u) for a tail, ln φ(z(u)) − ln s − u/2 + ln g(u) for the density, taken over
    // t = u − u₀. Where x and m share a sign, z vanishes at m·e^(u₀) = x, u₀ = ln(x/m), and z = −2·sign(x)·√(xm)·sinh(t/2)/s there: the
    // two terms of z cancel near u₀, and far out on the drift's side its peak is narrower than u's own rounding, so the integral is taken
    // about it, without the cancellation. Otherwise u₀ = 0 and z is the two terms.
    private double LogIntegral(double x, Integrand integrand)
    {
        double s = Scale, m = Drift, logRoot = 0.5 * Math.Log(2.0 * Math.PI) + Math.Log(s);
        bool density = integrand == Integrand.Density, shared = x != 0.0 && m != 0.0 && (x > 0.0) == (m > 0.0);
        double logX = shared ? Math.Log(Math.Abs(x)) : 0.0, logM = shared ? Math.Log(Math.Abs(m)) : 0.0;
        double origin = logX - logM, front = shared ? -2.0 * Math.Sign(x) * Math.Exp(0.5 * (logX + logM)) / s : 0.0;
        double H(double t)
        {
            // Each term only where it is not zero, so an infinite eᵘ or e⁻ᵘ far out does not meet a zero drift or x as 0·∞.
            double u = t + origin;
            double z = shared ? front * Math.Sinh(0.5 * t)
                              : ((x != 0.0 ? x * Math.Exp(-0.5 * u) : 0.0) - (m != 0.0 ? m * Math.Exp(0.5 * u) : 0.0)) / s;
            return integrand switch
            {
                Integrand.Below => Normal.StandardLogCdf(z),
                Integrand.Above => Normal.StandardLogCdf(-z),
                _ => -0.5 * z * z - 0.5 * u - logRoot,
            } + LogMixing(u);
        }
        double G(double t) => LogMixing(t + origin);
        // What ends a scan and a walk: under a tail g alone, which bounds h since Φ ≤ 1; for the density h itself, concave, so once it
        // falls far below the peak and still falls outward all beyond is smaller. Both −∞ counts as falling.
        double Level(double t) => density ? H(t) : G(t);
        bool Falling(double at, double inward) => density ? !(H(at) > H(inward)) : G(at) < G(inward);

        // The scan out from g's mode: the best point seen, and on each side the span until the bound is Drop below it.
        double start = _mode - origin, best = start, top = H(best), narrow = 0.0;
        if (density && shared && front * front >= 400.0)
        {
            // Far out on the drift's side the density's peak is a width 2/|front| ≤ 0.1 about t = 0, too narrow for a scan or a golden
            // section to find: Newton on h′(t) = −(front²/4)·sinh t + r′(t), r = −u/2 + ln g(u) the rest, smooth on g's scale, its width
            // from h″.
            double R(double t) => -0.5 * (t + origin) + G(t);
            const double d = 1e-4;
            double at = 0.0, curvature = 0.0;
            for (int k = 0; k < 20; k++)
            {
                double r0 = R(at), r1 = (R(at + d) - R(at - d)) / (2.0 * d), r2 = (R(at + d) - 2.0 * r0 + R(at - d)) / (d * d);
                curvature = 0.25 * front * front * Math.Cosh(at) - r2;
                double next = at + (-0.25 * front * front * Math.Sinh(at) + r1) / curvature;
                if (next == at) break;
                at = next;
            }
            (best, top, narrow) = (at, H(at), 1.0 / Math.Sqrt(curvature));
        }
        double Scan(double direction)
        {
            double step = 0.25, at = start;
            for (int k = 0; k < 10_000; k++)
            {
                at += direction * step;
                double h = H(at);
                if (h > top) { top = h; best = at; }
                if (Level(at) < top - Drop && Falling(at, at - direction * step)) return at;
                if (k >= 8) step *= 1.25;
            }
            return at;
        }
        double right = Scan(1.0), left = Scan(-1.0);
        if (double.IsNegativeInfinity(top)) return double.NegativeInfinity;

        // Golden section for the peak, then its width from the curvature: about the best point under a tail; over the whole span for the
        // density, whose concave integrand has one peak there, so a peak the scan stepped over is found.
        double peak, width;
        if (narrow > 0.0) (peak, width) = (best, Math.Min(narrow, 1.0));
        else
        {
            double lo = density ? left : Math.Max(left, best - 0.5), hi = density ? right : Math.Min(right, best + 0.5);
            const double golden = 0.6180339887498949;
            for (int k = 0, sections = density ? 64 : 25; k < sections; k++)
            {
                double c = hi - golden * (hi - lo), d = lo + golden * (hi - lo);
                if (H(c) > H(d)) hi = d;
                else lo = c;
            }
            peak = 0.5 * (lo + hi);
            top = Math.Max(top, H(peak));
            double e = 1e-4, curvature = -(H(peak + e) - 2.0 * H(peak) + H(peak - e)) / (e * e);
            width = curvature > 1.0 ? 1.0 / Math.Sqrt(curvature) : 1.0;
            width = Math.Clamp(width, 1e-6, 1.0);
        }

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
                if (panel <= 1e-17 * sum && H(at) < top - 45.0 && Falling(at, at - direction * 1e-3)) break;
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

    // The normal score z(x) = Φ⁻¹(F(x)) at points walked out from the median, with its slopes f/φ(z); cubic Hermite between them in
    // asinh(x/s), ln Φ(±z) for the tails, and the integral beyond the ends and in a cell with the pole at an end.
    private sealed class TailTable
    {
        private static readonly Normal Standard = Normal.MeanSd(0.0, 1.0);
        private readonly NormalVarianceMeanMixture _owner;
        private readonly double[] _at, _z, _slope;

        public TailTable(NormalVarianceMeanMixture owner)
        {
            _owner = owner;
            double median = owner.Invert(-SpecialFunctions.Ln2, below: true, tabled: false);
            var x = new List<double>();
            var z = new List<double>();
            var slope = new List<double>();
            foreach (bool below in new[] { true, false })
            {
                var side = new List<(double X, double Z, double Slope)>();
                double at = median, toward = below ? -1.0 : 1.0, floor = PoleFloor * owner.Scale;
                for (int k = 0; k < 10_000; k++)
                {
                    double tail = owner.LogTail(at, below);
                    double score = below ? Standard.InverseLogCdf(tail) : -Standard.InverseLogCdf(tail);
                    double rate = Math.Exp(owner.LogPdf(at) + 0.5 * score * score + 0.5 * Math.Log(2.0 * Math.PI));
                    side.Add((at, score, rate));
                    if (tail < TableEnd || !(rate > 0.0)) break;
                    if (!owner.HasPole)
                    {
                        if (!double.IsFinite(rate)) break;
                        at += toward * TableStep / rate;
                    }
                    // Near a pole the score goes as x·ln|x|, which a cubic follows only on cells short against |x|, its error going as the
                    // step's fourth power times |x|: towards it, and within s of it, steps of at most ¼ of |x|, ½ within 10⁻²·s, down to
                    // the floor, then onto the pole, and out from it the same way.
                    else if (at == 0.0) at = toward * floor;
                    else
                    {
                        bool towards = Math.Sign(at) != Math.Sign(toward);
                        if (towards && Math.Abs(at) <= floor) { at = 0.0; continue; }
                        double step = TableStep / rate;
                        double ratio = Math.Abs(at) < PoleNear * owner.Scale ? PoleNearRatio : PoleRatio;
                        if (towards || Math.Abs(at) < owner.Scale) step = Math.Min(step, ratio * Math.Max(Math.Abs(at), floor));
                        at += toward * step;
                    }
                }
                if (below) side.Reverse();
                foreach ((double X, double Z, double Slope) point in below ? side : side.GetRange(1, side.Count - 1))
                {
                    x.Add(point.X);
                    z.Add(point.Z);
                    slope.Add(point.Slope);
                }
            }
            // In asinh(x/s) the slope is (f/φ(z))·√(x² + s²), s·cosh of that coordinate.
            double s = owner.Scale;
            _at = [.. x.Select(v => Math.Asinh(v / s))];
            _slope = [.. slope.Select((rate, k) => rate * s * Math.Cosh(_at[k]))];
            _z = [.. z];
        }

        public double Read(double x, bool below)
        {
            if (double.IsNaN(x)) return double.NaN;
            double at = Math.Asinh(x / _owner.Scale);
            if (at < _at[0] || at > _at[^1])
            {
                // Past an end the side whose mass is near 1 there is read as ln(1 − e^(far tail)), the far tail integrated on its own side:
                // its own integral would lose 1 − F's digits and could read above 0.
                bool nearOne = at > _at[^1] ? below : !below;
                return nearOne ? SpecialFunctions.Log1p(-Math.Exp(_owner.LogTail(x, !below))) : _owner.LogTail(x, below);
            }
            int k = Array.BinarySearch(_at, at);
            double score;
            if (k >= 0) score = _z[k];
            else
            {
                k = Math.Clamp(~k - 1, 0, _at.Length - 2);
                if (!double.IsFinite(_slope[k]) || !double.IsFinite(_slope[k + 1])) return _owner.LogTail(x, below);
                double h = _at[k + 1] - _at[k], t = (at - _at[k]) / h, t2 = t * t, t3 = t2 * t;
                score = (2.0 * t3 - 3.0 * t2 + 1.0) * _z[k] + (t3 - 2.0 * t2 + t) * h * _slope[k]
                        + (-2.0 * t3 + 3.0 * t2) * _z[k + 1] + (t3 - t2) * h * _slope[k + 1];
            }
            return Normal.StandardLogCdf(below ? score : -score);
        }
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
