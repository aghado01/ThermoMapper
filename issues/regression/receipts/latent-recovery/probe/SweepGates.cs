using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace LatentRecoveryProbe;

/// <summary>Reads the sweeps' CSVs against the gates pre-stated in the brief's report and writes a summary.</summary>
public static class SweepGates
{
    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static List<Dictionary<string, string>> Csv(string path)
    {
        string[] lines = File.ReadAllLines(path);
        string[] head = lines[0].Split(',');
        return lines.Skip(1).Select(l => l.Split(',')).Select(f => head.Select((h, i) => (h, f[i])).ToDictionary(t => t.h, t => t.Item2)).ToList();
    }

    public static string B1(string dir)
    {
        var sb = new StringBuilder();
        var pairs = Csv(Path.Combine(dir, "b1-pairs.csv"));
        double maxZ = pairs.Max(r => Math.Abs(D(r["z"])));
        int marginal = pairs.Count(r => Math.Abs(D(r["z"])) is > 4 and <= 5), over = pairs.Count(r => Math.Abs(D(r["z"])) > 5);
        sb.AppendLine($"B1 pairs: {pairs.Count} kappas, max |z| = {maxZ:F2}, marginal {marginal}, over 5: {over}; " +
                      $"min |z| against the linear alternative = {pairs.Where(r => D(r["kappa"]) is > 0 and < 1).Min(r => Math.Abs(D(r["z_linear"]))):F0}");
        string pairsVerdict = over > 0 || marginal > 1 ? "FAIL" : marginal == 1 ? "MARGINAL (rerun)" : "PASS";
        sb.AppendLine($"  verdict: {pairsVerdict}");
        foreach (string noise in new[] { "complex", "real" })
        {
            var cells = Csv(Path.Combine(dir, $"b1-cells-{noise}.csv"));
            foreach (var family in cells.GroupBy(r => r["family"]))
            {
                var z = family.Select(r => Math.Abs(D(r["z"]))).ToList();
                var diff = family.Select(r => D(r["estimate"]) - D(r["predicted"])).ToList();
                sb.AppendLine($"B1 cells {noise}/{family.Key}: {z.Count} pairs, max |z| = {z.Max():F2}, in (4,5]: {z.Count(v => v is > 4 and <= 5)}, " +
                              $"over 5: {z.Count(v => v > 5)}, max |est − pred| = {diff.Max(Math.Abs):F4}, chi2/df = {family.Sum(r => D(r["z"]) * D(r["z"])) / z.Count:F2}");
            }
            var gated = cells.Where(r => noise == "complex" || r["family"] != "near-dc").ToList();
            int m = gated.Count(r => Math.Abs(D(r["z"])) is > 4 and <= 5), o = gated.Count(r => Math.Abs(D(r["z"])) > 5);
            sb.AppendLine($"  {noise} gated families verdict: {(o > 0 || m > 2 ? "FAIL" : m > 0 ? "MARGINAL (rerun)" : "PASS")}");
        }
        var nearDc = Csv(Path.Combine(dir, "b1-cells-real.csv")).Where(r => r["family"] == "near-dc").ToList();
        sb.AppendLine("  near-DC reading, |est − pred| against the pseudo entries (|pseudo(ref,other)|, |pseudo(ref,ref)|):");
        foreach (var r in nearDc.OrderByDescending(r => Math.Abs(D(r["estimate"]) - D(r["predicted"]))).Take(8))
            sb.AppendLine($"    w {r["w_ref"]} dt {r["dt"]} df·w {D(r["df"]) * D(r["w_ref"]):F1}: est {D(r["estimate"]):F4} pred {D(r["predicted"]):F4} z {D(r["z"]):F1} " +
                          $"pseudo {D(r["pseudo_ref_other"]):F3} self {D(r["pseudo_ref_self"]):F3}");
        int smallPseudoBad = nearDc.Count(r => D(r["pseudo_ref_other"]) < 0.05 && D(r["pseudo_ref_self"]) < 0.05 && Math.Abs(D(r["z"])) > 4);
        sb.AppendLine($"    near-DC pairs with both pseudo entries < 0.05 and |z| > 4: {smallPseudoBad} of {nearDc.Count(r => D(r["pseudo_ref_other"]) < 0.05 && D(r["pseudo_ref_self"]) < 0.05)}");
        return sb.ToString();
    }

    public static string B2(string dir)
    {
        var sb = new StringBuilder();
        var rows = Csv(Path.Combine(dir, "b2.csv"));
        bool all = true;
        void Check(bool ok, string what)
        {
            all &= ok;
            sb.AppendLine($"  [{(ok ? "pass" : "FAIL")}] {what}");
        }
        double GridStep = Math.Pow(2, 0.25);
        foreach (string wn in new[] { "rect", "hann" })
        {
            var tone = rows.Where(r => r["signal"] == "tone" && r["window"] == wn).OrderBy(r => int.Parse(r["w"])).ToList();
            double rel = tone.Max(r => Math.Abs(D(r["power"]) / D(r["formula"]) - 1));
            double over = tone.Max(r => D(r["grid_peak"]) / D(r["power"]) - 1);
            double slope = Stats.Slope(tone.Select(r => Math.Log(D(r["w"]))).ToList(), tone.Select(r => Math.Log(D(r["power"]))).ToList());
            Check(rel <= 1e-9 && over <= 1e-9 && Math.Abs(slope - 1) <= 0.01, $"{wn} tone: max rel dev from |Σg|² {rel:E1}, grid peak over f0 by {over:E1}, slope {slope:F4}");
            var imp = rows.Where(r => r["signal"] == "impulse" && r["window"] == wn).OrderBy(r => int.Parse(r["w"])).ToList();
            double islope = Stats.Slope(imp.Select(r => Math.Log(D(r["w"]))).ToList(), imp.Select(r => Math.Log(D(r["power"]))).ToList());
            double w4096 = imp.Where(r => r["w"] == "4096").Select(r => D(r["power"]) * 4096).First();
            Check(Math.Abs(islope + 1) <= 0.01, $"{wn} impulse: slope {islope:F4}, w·P at 4096 = {w4096:F4}");
            foreach (string kind in new[] { "complex", "real" })
            {
                var noise = rows.Where(r => r["signal"] == "noise-" + kind && r["window"] == wn).OrderBy(r => int.Parse(r["w"])).ToList();
                var z = noise.Select(r => (D(r["power"]) - 1) / D(r["extra"])).ToList();
                double nslope = Stats.Slope(noise.Select(r => Math.Log(D(r["w"]))).ToList(), noise.Select(r => Math.Log(D(r["power"]))).ToList());
                double pslope = Stats.Slope(noise.Select(r => Math.Log(D(r["w"]))).ToList(), noise.Select(r => Math.Log(D(r["grid_peak"]))).ToList());
                int marginal = z.Count(v => Math.Abs(v) is > 4 and <= 5), bad = z.Count(v => Math.Abs(v) > 5);
                Check(bad == 0 && marginal <= 1 && Math.Abs(nslope) <= 0.01,
                      $"{wn} noise {kind}: max |z| {z.Max(Math.Abs):F2} (marginal {marginal}), slope {nslope:F4}; frame-peak reading slope {pslope:F3}, " +
                      $"peak/mean at w=4096 {noise.Last()["grid_peak"]}");
            }
        }
        foreach (string wn in new[] { "rect", "hann" })
        {
            var ratios = new List<double>();
            foreach (var group in rows.Where(r => r["signal"] == "resonance" && r["window"] == wn).GroupBy(r => D(r["param"])).OrderBy(g => g.Key))
            {
                double tau = group.Key;
                var pts = group.OrderBy(r => int.Parse(r["w"])).ToList();
                var best = pts.MaxBy(r => D(r["power"]))!;
                double wStar = D(best["w"]);
                ratios.Add(wStar / tau);
                var small = pts.Where(r => D(r["w"]) <= tau / 8).ToList();
                var large = pts.Where(r => D(r["w"]) >= 8 * tau).ToList();
                double ss = Stats.Slope(small.Select(r => Math.Log(D(r["w"]))).ToList(), small.Select(r => Math.Log(D(r["power"]))).ToList());
                double ls = Stats.Slope(large.Select(r => Math.Log(D(r["w"]))).ToList(), large.Select(r => Math.Log(D(r["power"]))).ToList());
                double gridOverF0 = pts.Max(r => D(r["grid_peak"]) / D(r["power"]) - 1);
                string rectNote = wn == "rect" ? $", closed form max rel dev {pts.Max(r => Math.Abs(D(r["power"]) / D(r["formula"]) - 1)):E1}" : "";
                Check(wStar / tau is >= 0.5 and <= 4 && ss >= 0.85 && ls <= -0.9,
                      $"{wn} resonance τ={tau}: w* = {wStar} (w*/τ = {wStar / tau:F3}), slope small {ss:F3} ({small.Count} pts), large {ls:F3} ({large.Count} pts), grid peak over f0 by {gridOverF0:E1}{rectNote}");
            }
            Check(ratios.Max() / ratios.Min() <= 1.25, $"{wn} resonance: w*/τ spread {ratios.Max() / ratios.Min():F3}");
            var cratios = new List<double>();
            foreach (var group in rows.Where(r => r["signal"] == "chirp" && r["window"] == wn).GroupBy(r => D(r["param"])).OrderBy(g => g.Key))
            {
                double c = group.Key, scale = 1 / Math.Sqrt(c);
                var pts = group.OrderBy(r => int.Parse(r["w"])).ToList();
                var best = pts.MaxBy(r => D(r["power"]))!;
                double wStar = D(best["w"]);
                cratios.Add(wStar / scale);
                var small = pts.Where(r => D(r["w"]) <= scale / 4).ToList();
                var large = pts.Where(r => D(r["w"]) >= 8 * scale).ToList();
                double ss = Stats.Slope(small.Select(r => Math.Log(D(r["w"]))).ToList(), small.Select(r => Math.Log(D(r["power"]))).ToList());
                double ls = Stats.Slope(large.Select(r => Math.Log(D(r["w"]))).ToList(), large.Select(r => Math.Log(D(r["power"]))).ToList());
                Check(wStar / scale is >= 0.5 and <= 2 && ss >= 0.9 && ls <= -0.8,
                      $"{wn} chirp c={c}: w* = {wStar} (w*·√c = {wStar / scale:F3}), slope small {ss:F3} ({small.Count} pts), large {ls:F3} ({large.Count} pts)");
            }
            Check(cratios.Max() / cratios.Min() <= 1.25, $"{wn} chirp: w*·√c spread {cratios.Max() / cratios.Min():F3}");
        }
        sb.Insert(0, $"B2 verdict: {(all ? "PASS" : "FAIL (see rows)")}\n");
        return sb.ToString();
    }

    /// <summary>
    /// Post-gate readings on B2's own rows, not part of its gate: each turnover refined by a parabola in ln w through the
    /// grid's argmax and its neighbours, and each curve's local slope by octave of w/scale.
    /// </summary>
    public static string B2Supplement(string dir)
    {
        var sb = new StringBuilder("B2 supplement (post-gate readings on the gate's rows)\n");
        var rows = Csv(Path.Combine(dir, "b2.csv"));
        foreach (string signal in new[] { "resonance", "chirp" })
            foreach (string wn in new[] { "rect", "hann" })
                foreach (var group in rows.Where(r => r["signal"] == signal && r["window"] == wn).GroupBy(r => D(r["param"])).OrderBy(g => g.Key))
                {
                    double scale = signal == "resonance" ? group.Key : 1 / Math.Sqrt(group.Key);
                    var pts = group.OrderBy(r => int.Parse(r["w"])).Select(r => (X: Math.Log(D(r["w"])), Y: Math.Log(D(r["power"])))).ToList();
                    int i = pts.FindIndex(p => p.Y == pts.Max(q => q.Y));
                    double refined = Math.Exp(pts[i].X);
                    if (i > 0 && i < pts.Count - 1)
                    {
                        (double x0, double y0) = pts[i - 1]; (double x1, double y1) = pts[i]; (double x2, double y2) = pts[i + 1];
                        double denom = (x0 - x1) * (x0 - x2) * (x1 - x2);
                        double a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom;
                        double b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom;
                        refined = Math.Exp(-b / (2 * a));
                    }
                    var octaves = new List<string>();
                    for (double lo = 1.0 / 16; lo <= 64; lo *= 2)
                    {
                        var seg = pts.Where(p => Math.Exp(p.X) / scale >= lo && Math.Exp(p.X) / scale < 2 * lo).ToList();
                        if (seg.Count >= 2) octaves.Add($"[{lo:G3},{2 * lo:G3}) {Stats.Slope(seg.Select(p => p.X).ToList(), seg.Select(p => p.Y).ToList()):F3}");
                    }
                    sb.AppendLine($"{signal} {wn} {(signal == "resonance" ? "tau" : "c")}={group.Key}: refined w*/scale {refined / scale:F3}; local slopes by w/scale octave: {string.Join("; ", octaves)}");
                }
        return sb.ToString();
    }
}
