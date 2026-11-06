using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using Maths.Distributions;

namespace Maths.Oracle.Tests;

/// <summary>
/// What the GBASS parity classes share (E1's build brief, slice 0): running a <c>gbass_*_oracle.R</c> mode on a fixture written here,
/// reading its JSON, and a quadrature in logs for the laws the oracles integrate, so each reference is computed on both sides.
/// The modes run against CRAN's GBASS 2.0.1, pinned in <c>r/renv.lock</c>; a class returns only when R itself is absent, so a missing
/// package fails the class rather than passing it.
/// </summary>
internal static class GbassOracle
{
    /// <summary>Runs <paramref name="oracle"/> on <paramref name="rows"/>, written as a headerless CSV, with <paramref name="tail"/> after
    /// the output path; null when R is absent.</summary>
    public static JsonElement? Run(string oracle, IEnumerable<double[]> rows, params string[] tail)
    {
        if (!ROracle.IsAvailable) return null;
        string csv = WriteCsv(rows);
        try { return ROracle.Run($"oracles/{oracle}", csv, tail); }
        finally { File.Delete(csv); }
    }

    /// <summary>Writes <paramref name="rows"/> to a temporary headerless CSV, round-trip digits, and returns its path.</summary>
    public static string WriteCsv(IEnumerable<double[]> rows)
    {
        var text = new StringBuilder();
        foreach (double[] row in rows)
            text.AppendLine(string.Join(",", row.Select(v => v.ToString("R", CultureInfo.InvariantCulture))));
        string path = Path.Combine(Path.GetTempPath(), $"gbass_{Guid.NewGuid():N}.csv");
        File.WriteAllText(path, text.ToString());
        return path;
    }

    /// <summary>A number as the oracles write it: R's non-finite values come as strings ("NaN", "Inf", "-Inf", "NA").</summary>
    public static double Number(JsonElement e) => e.ValueKind switch
    {
        JsonValueKind.Number => e.GetDouble(),
        JsonValueKind.String => e.GetString() switch
        {
            "Inf" => double.PositiveInfinity,
            "-Inf" => double.NegativeInfinity,
            _ => double.NaN,
        },
        JsonValueKind.Null => double.NaN,
        _ => throw new InvalidOperationException($"Not a number: {e}"),
    };

    /// <summary>An array of numbers; a scalar reads as one.</summary>
    public static double[] Numbers(JsonElement e) => e.ValueKind == JsonValueKind.Array ? [.. e.EnumerateArray().Select(Number)] : [Number(e)];

    /// <summary>
    /// ln ∫ exp(h(u)) du over (<paramref name="lower"/>, <paramref name="upper"/>) for a unimodal h, by 40-point Gauss–Legendre panels
    /// walked out from h's peak, found from <paramref name="guess"/>, on each side until a panel adds less than e⁻⁴⁰ of the total and h is
    /// falling. A panel spans at most one standard deviation of h's curvature at the peak, and each is summed about its own largest term,
    /// so far ranges keep their digits.
    /// </summary>
    public static double LogIntegral(Func<double, double> h, double guess, double lower = double.NegativeInfinity, double upper = double.PositiveInfinity)
    {
        double peak = Peak(h, guess, lower, upper);
        double e = 1e-4 * Math.Max(1.0, Math.Abs(peak));
        double curvature = -(h(peak + e) - 2.0 * h(peak) + h(peak - e)) / (e * e);
        double width = curvature > 4.0 ? 1.0 / Math.Sqrt(curvature) : 0.5;
        var panels = new List<double>();
        foreach (int side in new[] { 1, -1 })
        {
            double best = double.NegativeInfinity, step = width, a = peak;
            for (int k = 0; k < 200000; k++)
            {
                if (side > 0 ? a >= upper : a <= lower) break;
                double b = side > 0 ? Math.Min(a + step, upper) : Math.Max(a - step, lower);
                double panel = LogPanel(h, Math.Min(a, b), Math.Max(a, b));
                panels.Add(panel);
                best = Math.Max(best, panel);
                if (panel < best - 40.0 && h(b) < h(a)) break;
                a = b;
                if (k > 20) step = Math.Min(0.5, step * 1.05);   // the tails are smooth: widen as they fall
            }
        }
        double top = panels.Max();
        return top + Math.Log(panels.Sum(p => Math.Exp(p - top)));
    }

    // The peak of a unimodal h: bracket by doubling steps from the guess, then golden section.
    private static double Peak(Func<double, double> h, double guess, double lower, double upper)
    {
        double Clamp(double u) => Math.Min(Math.Max(u, lower + 1e-12), upper - 1e-12);
        double x = Clamp(guess), step = 0.25;
        int dir = h(Clamp(x + 1e-6)) >= h(x) ? 1 : -1;
        double a = x, b = x;
        for (int k = 0; k < 200; k++)
        {
            double next = Clamp(b + dir * step);
            if (h(next) < h(b) || next == b) break;
            a = b; b = next; step *= 2.0;
        }
        double lo = Math.Min(a, Clamp(b + dir * step)), hi = Math.Max(a, Clamp(b + dir * step));
        double g = (Math.Sqrt(5.0) - 1.0) / 2.0;
        for (int k = 0; k < 200 && hi - lo > 1e-12 * Math.Max(1.0, Math.Abs(lo)); k++)
        {
            double c = hi - g * (hi - lo), d = lo + g * (hi - lo);
            if (h(c) >= h(d)) hi = d; else lo = c;
        }
        return 0.5 * (lo + hi);
    }

    private static double LogPanel(Func<double, double> h, double a, double b)
    {
        GaussLegendre rule = GaussLegendre.Of(40);
        double mid = 0.5 * (a + b), half = 0.5 * (b - a);
        var terms = new double[rule.Count];
        for (int i = 0; i < rule.Count; i++) terms[i] = Math.Log(rule.Weights[i] * half) + h(mid + half * rule.Nodes[i]);
        double top = terms.Max();
        if (double.IsNegativeInfinity(top)) return top;
        return top + Math.Log(terms.Sum(t => Math.Exp(t - top)));
    }

    /// <summary>The binomial z-scores of draws' sample deciles against a law's CDF there: (F(q̂ₖ) − pₖ)/√(pₖ(1 − pₖ)/n), pₖ = k/10.</summary>
    public static double[] DecileScores(IReadOnlyList<double> cdf, int n)
        => [.. cdf.Select((f, k) => { double p = 0.1 * (k + 1); return (f - p) / Math.Sqrt(p * (1.0 - p) / n); })];
}
