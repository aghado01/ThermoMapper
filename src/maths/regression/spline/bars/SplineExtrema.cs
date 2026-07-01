using System;
using System.Collections.Generic;
using System.Numerics;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>A significant peak with its level-crossing span: the peak <paramref name="Location"/>/<paramref name="Height"/>,
/// the <paramref name="Left"/>/<paramref name="Right"/> crossings of the drop level, and whether either edge was
/// clipped by the domain boundary (the curve never descended to the level before the edge — a truncated span).</summary>
public readonly record struct PeakSpan(
    double Location, double Height, double Left, double Right, bool LeftClipped, bool RightClipped);

/// <summary>
/// Exact critical-point analysis of a fitted spline by the closed-form candidate set (Fermat): the extrema of
/// a piecewise polynomial lie at span boundaries or at interior derivative roots. For cubic (degree ≤ 3) spans —
/// the default, hot path — the span polynomial is recovered exactly from four evaluations and its derivative (a
/// quadratic) solved in closed form. Higher-degree spans are reconstructed from <c>degree+1</c> evaluations,
/// differentiated, and root-found generally (Durand–Kerner), so the analysis is exact for any spline degree while
/// the cubic regime stays on the closed form unchanged. No grid scan, no optimizer, so the per-draw peak the
/// ensemble pools carries zero optimizer slop (the <c>argmax_in_closed_form_set</c> spec). Domain is [0,1].
/// </summary>
public static class SplineExtrema
{
    /// <summary>Location and height of the global maximum of <c>Σ coef_j B_j(x)</c> on [0,1].</summary>
    public static (double Location, double Height) Argmax(KnotConfig config, double[] coef, SplineBasis basis)
    {
        (List<double> xs, List<double> fs) = CriticalPoints(config, coef, basis);
        double bx = xs[0], bf = fs[0];
        for (int i = 1; i < fs.Count; i++)
            if (fs[i] > bf) { bf = fs[i]; bx = xs[i]; }
        return (bx, bf);
    }

    /// <summary>
    /// The significant peaks — every local maximum whose topographic prominence is at least
    /// <paramref name="relativeProminence"/> times the curve's total range, as (location, height), tiny
    /// wiggles filtered out. The matching-free per-draw peak set the intensity readout pools into λ(T).
    /// </summary>
    /// <remarks>
    /// <b>DOMAIN-PREMISE:</b> The prominence gate and the closed-form critical enumeration are method-intrinsic
    /// (owned here). Whether a <i>boundary</i> maximum — included here — is a real "transition" is a domain
    /// premise the consumer owns, not a property of the curve: a bounded sweep that rejects edge transitions
    /// filters them downstream; BARS does not decide it.
    /// </remarks>
    public static List<(double Location, double Height)> SignificantPeaks(
        KnotConfig config, double[] coef, SplineBasis basis, double relativeProminence)
    {
        var peaks = new List<(double Location, double Height)>();
        (List<double> xs, List<double> fs) = CriticalPoints(config, coef, basis);
        int n = fs.Count;
        if (n < 2) return peaks;

        double fmax = fs[0], fmin = fs[0];
        for (int i = 1; i < n; i++) { if (fs[i] > fmax) fmax = fs[i]; if (fs[i] < fmin) fmin = fs[i]; }
        double range = fmax - fmin;
        if (range <= 0.0) return peaks;
        double threshold = relativeProminence * range;

        for (int i = 0; i < n; i++)
        {
            bool isMax = i == 0 ? fs[0] > fs[1]
                       : i == n - 1 ? fs[n - 1] > fs[n - 2]
                       : fs[i] > fs[i - 1] && fs[i] > fs[i + 1];
            if (!isMax) continue;

            double leftCol = i == 0 ? double.NegativeInfinity : ValleyFloor(fs, n, i, -1);
            double rightCol = i == n - 1 ? double.NegativeInfinity : ValleyFloor(fs, n, i, +1);
            double prominence = fs[i] - Math.Max(leftCol, rightCol);
            if (prominence >= threshold) peaks.Add((xs[i], fs[i]));
        }
        return peaks;
    }

    /// <summary>
    /// Number of <see cref="SignificantPeaks"/> — the "how many transitions" count, tiny wiggles filtered out.
    /// </summary>
    public static int SignificantPeakCount(KnotConfig config, double[] coef, SplineBasis basis, double relativeProminence)
        => SignificantPeaks(config, coef, basis, relativeProminence).Count;

    /// <summary>
    /// The significant peaks with their level-crossing spans — for each peak, the interval where the curve sits
    /// within <paramref name="dropFraction"/> of its prominence below the apex (FWHM is the half-drop instance).
    /// Crossings are exact closed-form roots of the per-span polynomial (cubic closed form for degree ≤ 3, general
    /// reconstruction for higher degree), no scan; clipped flags mark a span the curve
    /// never descended out of before the domain edge. The matching-free per-draw span set the coverage readout
    /// pools into π(T).
    /// </summary>
    /// <remarks>
    /// <b>DOMAIN-PREMISE:</b> The crossing arithmetic is method-intrinsic (owned here). The structural span is the
    /// FWHM of the curve itself, computed per draw — distinct from the credible interval on the peak location. The
    /// <paramref name="dropFraction"/> (½ = FWHM-analogue) and the prominence-relative baseline are the consumer's
    /// premise — what "the span around the peak" means for their application. A <b>clipped</b> span means the curve
    /// never descends to the level within the [0,1] domain. Under a global uniform first-pass fit (BARS) this is
    /// <i>not</i> a "recon window too narrow" signal — there is no separate recon window; the whole curve is fit
    /// over [0,1] (endpoints inclusive, clamped knots pinning the boundary values). It is a domain premise the
    /// consumer owns: either a genuine boundary/edge transition, or a temperature bracket set too narrow
    /// <i>upstream</i> of the fit — BARS surfaces it, the consumer decides which.
    /// </remarks>
    public static List<PeakSpan> SignificantPeakSpans(
        KnotConfig config, double[] coef, SplineBasis basis, double relativeProminence, double dropFraction)
    {
        var spans = new List<PeakSpan>();
        int degree = basis.Degree;
        double[] clamped = BSpline.MakeClampedKnots(config.InteriorKnots, degree);
        (List<double> xs, List<double> fs) = CriticalPoints(config, coef, basis);
        int n = fs.Count;
        if (n < 2) return spans;

        double fmax = fs[0], fmin = fs[0];
        for (int i = 1; i < n; i++) { if (fs[i] > fmax) fmax = fs[i]; if (fs[i] < fmin) fmin = fs[i]; }
        double range = fmax - fmin;
        if (range <= 0.0) return spans;
        double threshold = relativeProminence * range;

        for (int i = 0; i < n; i++)
        {
            bool isMax = i == 0 ? fs[0] > fs[1]
                       : i == n - 1 ? fs[n - 1] > fs[n - 2]
                       : fs[i] > fs[i - 1] && fs[i] > fs[i + 1];
            if (!isMax) continue;

            double leftCol = i == 0 ? double.NegativeInfinity : ValleyFloor(fs, n, i, -1);
            double rightCol = i == n - 1 ? double.NegativeInfinity : ValleyFloor(fs, n, i, +1);
            double prominence = fs[i] - Math.Max(leftCol, rightCol);
            if (prominence < threshold) continue;

            double level = fs[i] - dropFraction * prominence;

            // Left edge: first critical/break point (scanning toward 0) at or below the level; the crossing lies in
            // the monotone bracket between it and its inward neighbour. If none before the boundary, the span clips.
            double left = xs[0];
            bool leftClip = true;
            for (int j = i - 1; j >= 0; j--)
                if (fs[j] <= level)
                {
                    left = CrossingInBracket(clamped, degree, coef, xs[j], xs[j + 1], level);
                    leftClip = false;
                    break;
                }

            // Right edge: symmetric, scanning toward 1.
            double right = xs[n - 1];
            bool rightClip = true;
            for (int j = i + 1; j < n; j++)
                if (fs[j] <= level)
                {
                    right = CrossingInBracket(clamped, degree, coef, xs[j - 1], xs[j], level);
                    rightClip = false;
                    break;
                }

            spans.Add(new PeakSpan(xs[i], fs[i], left, right, leftClip, rightClip));
        }
        return spans;
    }

    // Lowest point reached scanning in `dir` from peak i, until a strictly higher point (exclusive) or the boundary.
    private static double ValleyFloor(List<double> fs, int n, int i, int dir)
    {
        double h = fs[i];
        double floor = double.PositiveInfinity;
        for (int j = i + dir; j >= 0 && j < n; j += dir)
        {
            if (fs[j] > h) return floor;
            if (fs[j] < floor) floor = fs[j];
        }
        return floor;
    }

    // Critical points (boundaries + per-span derivative roots) with heights, sorted ascending by location.
    private static (List<double> Xs, List<double> Fs) CriticalPoints(KnotConfig config, double[] coef, SplineBasis basis)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(coef);
        ArgumentNullException.ThrowIfNull(basis);

        int degree = basis.Degree;
        double[] clamped = BSpline.MakeClampedKnots(config.InteriorKnots, degree);

        double Eval(double x) => EvalSpline(clamped, degree, coef, x);

        var breaks = new List<double> { 0.0 };
        foreach (double t in config.InteriorKnots)
            if (t > 1e-12 && t < 1.0 - 1e-12) breaks.Add(t);
        breaks.Add(1.0);
        breaks.Sort();

        var xs = new List<double>();
        var fs = new List<double>();
        void Record(double x) { xs.Add(x); fs.Add(Eval(x)); }

        Record(0.0);
        for (int s = 0; s + 1 < breaks.Count; s++)
        {
            double a = breaks[s], b = breaks[s + 1], h = b - a;
            if (h < 1e-12) continue;

            var roots = new List<double>();
            foreach (double u in InteriorCriticalParameters(Eval, a, b, degree)) roots.Add(a + h * u);
            roots.Sort();
            foreach (double rx in roots) Record(rx);
            Record(b);
        }
        return (xs, fs);
    }

    // Derivative roots of one span polynomial, as parameters u ∈ (0,1) with x = a + (b−a)·u.
    // Degree ≤ 3 (the default): the polynomial is recovered exactly from four evaluations and its derivative —
    // a quadratic — solved in closed form (the unchanged, allocation-light fast path). Degree ≥ 4: reconstruct
    // the degree-d polynomial from d+1 evaluations, differentiate, and root-find generally.
    private static IEnumerable<double> InteriorCriticalParameters(Func<double, double> eval, double a, double b, int degree)
    {
        double h = b - a;
        if (degree <= 3)
        {
            double y0 = eval(a), y1 = eval(a + h / 3.0), y2 = eval(a + 2.0 * h / 3.0), y3 = eval(b);
            double d1 = y1 - y0, d2 = y2 - y0, d3 = y3 - y0;
            double q1 = 9.0 * d1 - 4.5 * d2 + d3;
            double q2 = (27.0 * d1 - d3) / 2.0 - 4.0 * q1;
            double q3 = d3 - q1 - q2;
            return RootsInUnitInterval(3.0 * q3, 2.0 * q2, q1);
        }

        // Reconstruct p(u) on degree+1 equally-spaced nodes, differentiate, return p'(u) roots in (0,1).
        int n = degree + 1;
        var u = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            u[i] = i / (double)degree;
            y[i] = eval(a + h * u[i]);
        }
        double[] c = InterpolatingMonomialCoeffs(u, y);     // p(u) = Σ c[k] u^k
        var deriv = new double[degree];                     // p'(u) = Σ deriv[k] u^k
        for (int k = 1; k <= degree; k++) deriv[k - 1] = k * c[k];
        return RealRootsInUnit(deriv);
    }

    // Monomial coefficients of the polynomial interpolating (u_i, y_i) — a small Vandermonde solve.
    private static double[] InterpolatingMonomialCoeffs(double[] u, double[] y)
    {
        int n = u.Length;
        var v = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            double p = 1.0;
            for (int j = 0; j < n; j++) { v[i, j] = p; p *= u[i]; }
        }
        return SolveLinear(v, (double[])y.Clone());
    }

    // Dense linear solve by Gaussian elimination with partial pivoting (n is degree+1, small).
    private static double[] SolveLinear(double[,] a, double[] b)
    {
        int n = b.Length;
        for (int col = 0; col < n; col++)
        {
            int piv = col;
            double best = Math.Abs(a[col, col]);
            for (int r = col + 1; r < n; r++) { double m = Math.Abs(a[r, col]); if (m > best) { best = m; piv = r; } }
            if (piv != col)
            {
                for (int j = 0; j < n; j++) (a[col, j], a[piv, j]) = (a[piv, j], a[col, j]);
                (b[col], b[piv]) = (b[piv], b[col]);
            }
            double d = a[col, col];
            if (Math.Abs(d) < 1e-300) continue;
            for (int r = col + 1; r < n; r++)
            {
                double f = a[r, col] / d;
                for (int j = col; j < n; j++) a[r, j] -= f * a[col, j];
                b[r] -= f * b[col];
            }
        }
        var x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double sgi = b[i];
            for (int j = i + 1; j < n; j++) sgi -= a[i, j] * x[j];
            x[i] = Math.Abs(a[i, i]) < 1e-300 ? 0.0 : sgi / a[i, i];
        }
        return x;
    }

    // Real roots in (0,1) of Σ coeffs[k] u^k. Low degrees use the closed forms; degree ≥ 3 uses Durand–Kerner.
    private static IEnumerable<double> RealRootsInUnit(double[] coeffs)
    {
        int deg = coeffs.Length - 1;
        while (deg > 0 && Math.Abs(coeffs[deg]) < 1e-12) deg--;
        if (deg <= 0) yield break;
        if (deg == 1)
        {
            double u = -coeffs[0] / coeffs[1];
            if (u > 0.0 && u < 1.0) yield return u;
            yield break;
        }
        if (deg == 2)
        {
            foreach (double u in RootsInUnitInterval(coeffs[2], coeffs[1], coeffs[0])) yield return u;
            yield break;
        }
        foreach (double u in DurandKernerRealRootsInUnit(coeffs, deg)) yield return u;
    }

    // All complex roots of a real polynomial via Durand–Kerner (Weierstrass), keeping the real ones in (0,1).
    private static IEnumerable<double> DurandKernerRealRootsInUnit(double[] coeffs, int deg)
    {
        var a = new Complex[deg + 1];
        double lead = coeffs[deg];
        for (int k = 0; k <= deg; k++) a[k] = new Complex(coeffs[k] / lead, 0.0);   // monic

        var z = new Complex[deg];
        var seed = new Complex(0.4, 0.9);
        Complex pw = Complex.One;
        for (int k = 0; k < deg; k++) { z[k] = pw; pw *= seed; }

        for (int iter = 0; iter < 200; iter++)
        {
            double maxDelta = 0.0;
            for (int k = 0; k < deg; k++)
            {
                Complex den = Complex.One;
                for (int j = 0; j < deg; j++) if (j != k) den *= z[k] - z[j];
                if (den == Complex.Zero) continue;
                Complex delta = EvalComplex(a, z[k]) / den;
                z[k] -= delta;
                if (delta.Magnitude > maxDelta) maxDelta = delta.Magnitude;
            }
            if (maxDelta < 1e-14) break;
        }

        var found = new List<double>();
        foreach (Complex root in z)
            if (Math.Abs(root.Imaginary) < 1e-7 && root.Real > 1e-9 && root.Real < 1.0 - 1e-9)
                found.Add(root.Real);
        found.Sort();
        return found;
    }

    private static Complex EvalComplex(Complex[] a, Complex x)
    {
        Complex r = Complex.Zero;
        for (int k = a.Length - 1; k >= 0; k--) r = r * x + a[k];
        return r;
    }

    /// <summary>Roots of <c>A u² + B u + C</c> in the open interval (0,1).</summary>
    private static IEnumerable<double> RootsInUnitInterval(double a, double b, double c)
    {
        const double eps = 1e-12;
        if (Math.Abs(a) < eps)
        {
            if (Math.Abs(b) >= eps)
            {
                double u = -c / b;
                if (u > 0.0 && u < 1.0) yield return u;
            }
            yield break;
        }

        double disc = b * b - 4.0 * a * c;
        if (disc < 0.0) yield break;
        double sq = Math.Sqrt(disc);
        double u1 = (-b + sq) / (2.0 * a);
        double u2 = (-b - sq) / (2.0 * a);
        if (u1 > 0.0 && u1 < 1.0) yield return u1;
        if (u2 > 0.0 && u2 < 1.0 && Math.Abs(u2 - u1) > eps) yield return u2;
    }

    private static double EvalSpline(double[] clamped, int degree, double[] coef, double x)
    {
        double[] row = BSpline.EvaluateBasis(x, clamped, degree);
        double f = 0.0;
        for (int j = 0; j < row.Length; j++) f += row[j] * coef[j];
        return f;
    }

    // Exact level crossing f(x) = L inside a single-span monotone bracket [xa, xb], straddling L. Recovers that
    // span's cubic from four evaluations (exact, since the bracket lies in one polynomial piece) and roots it.
    private static double CrossingInBracket(double[] clamped, int degree, double[] coef, double xa, double xb, double L)
    {
        double w = xb - xa;
        if (w < 1e-15) return xa;

        if (degree <= 3)
        {
            // Cubic (or lower): the span level p(v) = L is recovered exactly from four evaluations. (Unchanged.)
            double y0 = EvalSpline(clamped, degree, coef, xa);
            double y1 = EvalSpline(clamped, degree, coef, xa + w / 3.0);
            double y2 = EvalSpline(clamped, degree, coef, xa + 2.0 * w / 3.0);
            double y3 = EvalSpline(clamped, degree, coef, xb);
            double d1 = y1 - y0, d2 = y2 - y0, d3 = y3 - y0;
            double q1 = 9.0 * d1 - 4.5 * d2 + d3;
            double q2 = (27.0 * d1 - d3) / 2.0 - 4.0 * q1;
            double q3 = d3 - q1 - q2;
            // f(v) = y0 + q1 v + q2 v² + q3 v³ = L  ⇒  q3 v³ + q2 v² + q1 v + (y0 − L) = 0, v ∈ (0,1).
            foreach (double v in RootsCubicInUnitInterval(q3, q2, q1, y0 - L))
                return xa + w * v;
            double t = Math.Abs(y3 - y0) < 1e-15 ? 0.5 : (L - y0) / (y3 - y0);   // linear degrade if no interior root
            return xa + w * Math.Clamp(t, 0.0, 1.0);
        }

        // Degree ≥ 4: reconstruct the degree-d span polynomial from d+1 evaluations and solve p(v) = L exactly —
        // the level-crossing twin of the critical-point generalization, reusing the same helpers with L folded
        // into the constant term.
        int m = degree + 1;
        var u = new double[m];
        var ys = new double[m];
        for (int i = 0; i < m; i++)
        {
            u[i] = i / (double)degree;
            ys[i] = EvalSpline(clamped, degree, coef, xa + w * u[i]);
        }
        double[] c = InterpolatingMonomialCoeffs(u, ys);
        c[0] -= L;                                          // roots of p(v) − L in (0,1)
        foreach (double v in RealRootsInUnit(c))
            return xa + w * v;
        double tl = Math.Abs(ys[m - 1] - ys[0]) < 1e-15 ? 0.5 : (L - ys[0]) / (ys[m - 1] - ys[0]);
        return xa + w * Math.Clamp(tl, 0.0, 1.0);
    }

    /// <summary>Real roots of <c>c3·v³ + c2·v² + c1·v + c0</c> in the open interval (0,1), closed-form (Cardano/trig).</summary>
    private static IEnumerable<double> RootsCubicInUnitInterval(double c3, double c2, double c1, double c0)
    {
        const double eps = 1e-12;
        if (Math.Abs(c3) < eps)
        {
            foreach (double u in RootsInUnitInterval(c2, c1, c0)) yield return u;
            yield break;
        }

        double b = c2 / c3, c = c1 / c3, d = c0 / c3;
        double p = c - b * b / 3.0;
        double q = 2.0 * b * b * b / 27.0 - b * c / 3.0 + d;
        double shift = b / 3.0;
        double disc = q * q / 4.0 + p * p * p / 27.0;

        var roots = new List<double>(3);
        if (disc > eps)                          // one real root (Cardano)
        {
            double sq = Math.Sqrt(disc);
            roots.Add(Math.Cbrt(-q / 2.0 + sq) + Math.Cbrt(-q / 2.0 - sq) - shift);
        }
        else if (disc < -eps)                    // three real roots (trig; p < 0 here)
        {
            double m = 2.0 * Math.Sqrt(-p / 3.0);
            double theta = Math.Acos(Math.Clamp(3.0 * q / (m * p), -1.0, 1.0)) / 3.0;
            for (int k = 0; k < 3; k++)
                roots.Add(m * Math.Cos(theta - 2.0 * Math.PI * k / 3.0) - shift);
        }
        else                                     // multiple root
        {
            roots.Add((Math.Abs(p) < eps ? 0.0 : 3.0 * q / p) - shift);
            if (Math.Abs(p) >= eps) roots.Add(-3.0 * q / (2.0 * p) - shift);
        }

        roots.Sort();
        foreach (double u in roots)
            if (u > 0.0 && u < 1.0) yield return u;
    }
}
