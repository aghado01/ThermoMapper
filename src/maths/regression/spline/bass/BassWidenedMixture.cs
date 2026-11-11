using System;
using Maths.Distributions;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A new group's noise under GBASS's mixtures (E21's slice 5, the candidate seam): e | v ~ N(m·v, s²·v + k), v from the mixing law, with
/// m = σβ, s² = σ²c, and k = σ²·Σ s_k²·z̃_k² the new groups' effects integrated from their prior. Its variance is affine in v, so it is no
/// normal variance-mean mixture, and its density and both tails are integrals in u = ln v of exp(h(u)): h the log normal term, or ln Φ of
/// either tail's score, plus ln p(eᵘ) + u. A scan of h on a grid in u, widened until it holds the peak, and a golden section find h's
/// largest value. Panels walk out from it, each integrated by adaptive Gauss–Kronrod (7–15) as <see cref="NormalVarianceMeanMixture"/>'s
/// tails are, and doubling in width, until a panel adds nothing where h has fallen 45 below the peak. Each tail is integrated on its own side,
/// so a far tail keeps its digits as a log. At k = 0 it is the mixture's own marginal (<see cref="IMixingLaw.MarginalLogDensity"/>), which
/// the tests read beside a brute-force sum over u.
/// </summary>
internal static class BassWidenedMixture
{
    private const int MaxDepth = 40;
    private const double ScanStep = 0.5, ScanHalf = 40.0, ScanLimit = 600.0, Drop = 45.0;

    /// <summary>ln ∫ N(e | m·v, s²·v + k) p(v) dv.</summary>
    public static double LogDensity(IMixingLaw mixing, double residual, double unitVariance, double unitDrift, double widened)
    {
        Check(unitVariance, widened);
        return LogIntegral(u =>
        {
            double v = Math.Exp(u), variance = unitVariance * v + widened, d = residual - unitDrift * v;
            return -0.5 * Math.Log(2.0 * Math.PI * variance) - d * d / (2.0 * variance) + mixing.LogPrior(v) + u;
        });
    }

    /// <summary>ln F(e) = ln ∫ Φ((e − m·v)/√(s²·v + k)) p(v) dv, each half read on its own tail (<see cref="NoiseLogCdf.FromTails"/>).</summary>
    public static double LogCdf(IMixingLaw mixing, double residual, double unitVariance, double unitDrift, double widened)
    {
        Check(unitVariance, widened);
        double Tail(double u, double sign)
        {
            double v = Math.Exp(u), z = sign * (residual - unitDrift * v) / Math.Sqrt(unitVariance * v + widened);
            return Normal.StandardLogCdf(z) + mixing.LogPrior(v) + u;
        }
        return NoiseLogCdf.FromTails(LogIntegral(u => Tail(u, 1.0)), LogIntegral(u => Tail(u, -1.0)));
    }

    private static void Check(double unitVariance, double widened)
    {
        if (!(unitVariance > 0.0 && double.IsFinite(unitVariance))) throw new ArgumentOutOfRangeException(nameof(unitVariance), "s² must be positive and finite.");
        if (!(widened >= 0.0 && double.IsFinite(widened))) throw new ArgumentOutOfRangeException(nameof(widened), "k must be finite and not negative.");
    }

    // ln ∫ exp(h(u)) du over the line: the peak by a scan and a golden section, then panels walked out from it on both sides.
    private static double LogIntegral(Func<double, double> h)
    {
        double H(double u)
        {
            double value = h(u);
            return double.IsNaN(value) ? double.NegativeInfinity : value;
        }
        // The scan, widened on the side its best point sits on until that point is interior.
        double lo = -ScanHalf, hi = ScanHalf, best = 0.0, top = double.NegativeInfinity;
        while (true)
        {
            for (double u = lo; u <= hi; u += ScanStep)
            {
                double value = H(u);
                if (value > top) (top, best) = (value, u);
            }
            if (!double.IsFinite(top)) return double.NegativeInfinity;
            if (best - lo < ScanStep && lo > -ScanLimit) { hi = lo; lo -= 2.0 * ScanHalf; continue; }
            if (hi - best < ScanStep && hi < ScanLimit) { lo = hi; hi += 2.0 * ScanHalf; continue; }
            break;
        }
        // A golden section in the scan's bracket about its best point.
        double a = best - ScanStep, b = best + ScanStep;
        const double Ratio = 0.6180339887498949;
        double c = b - Ratio * (b - a), d = a + Ratio * (b - a), hc = H(c), hd = H(d);
        for (int i = 0; i < 60 && b - a > 1e-9; i++)
        {
            if (hc > hd) { b = d; d = c; hd = hc; c = b - Ratio * (b - a); hc = H(c); }
            else { a = c; c = d; hc = hd; d = a + Ratio * (b - a); hd = H(d); }
        }
        double peak = hc > hd ? c : d;
        top = Math.Max(top, Math.Max(hc, hd));
        double sum = 0.0;
        foreach (double direction in new[] { 1.0, -1.0 })
        {
            double at = peak, step = 0.25;
            while (Math.Abs(at - peak) < 2.0 * ScanLimit)
            {
                double next = at + direction * step;
                double panel = Kronrod(H, Math.Min(at, next), Math.Max(at, next), top, 0);
                sum += panel;
                at = next;
                step *= 2.0;
                if (panel <= 1e-17 * sum && H(at) < top - Drop) break;
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
