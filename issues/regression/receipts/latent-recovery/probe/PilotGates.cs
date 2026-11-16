using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;

namespace LatentRecoveryProbe;

/// <summary>Reads the pilot's replicate rows against the Stage 1 gate as pre-stated (R1–R5) and tabulates every cell.</summary>
public static class PilotGates
{
    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);
    private static double[] V(string s) => s.Split(';').Select(D).ToArray();

    public sealed record Row(int K, int Rep, string Emulator, string Mode, double Noise, int Ratio, int M, bool Disc, double[] Truth, double[] Mean,
                             double[] Q05, double[] Q95, double[] Sd, double Err, int Covered, int Converged, int Runs, double Sigma, double ReconTruth,
                             double ReconObs, double FitObs, double MassTruth, double MassDouble, double MassHalf, double Seconds, double Snr);

    public static List<Row> Load(string path)
    {
        string[] lines = File.ReadAllLines(path);
        string[] head = lines[0].Split(',');
        int I(string name) => Array.IndexOf(head, name);
        return lines.Skip(1).Select(l => l.Split(',')).Select(f => new Row(
            int.Parse(f[I("k")]), int.Parse(f[I("rep")]), f[I("emulator")], f[I("mode")], D(f[I("noise")]), int.Parse(f[I("ratio")]), int.Parse(f[I("m")]),
            f[I("disc")] == "1", V(f[I("truth_unit")]), V(f[I("mean_unit")]), V(f[I("q05_unit")]), V(f[I("q95_unit")]), V(f[I("sd_unit")]),
            D(f[I("err_rms")]), int.Parse(f[I("covered")]), int.Parse(f[I("converged")]), int.Parse(f[I("runs")]), D(f[I("sigma_post")]),
            D(f[I("recon_truth_held")]), D(f[I("recon_obs_held")]), D(f[I("fit_obs_kept")]), D(f[I("mass_truth")]), D(f[I("mass_double")]),
            D(f[I("mass_half")]), D(f[I("seconds")]), D(f[I("snr_db")]))).ToList();
    }

    public static double Median(IEnumerable<double> values)
    {
        double[] v = values.Where(x => !double.IsNaN(x)).OrderBy(x => x).ToArray();
        if (v.Length == 0) return double.NaN;
        return v.Length % 2 == 1 ? v[v.Length / 2] : 0.5 * (v[v.Length / 2 - 1] + v[v.Length / 2]);
    }

    /// <summary>The two-sided acceptance band [lo, hi] of Binomial(n, p) at level 1 − alpha: each tail at most alpha/2.</summary>
    public static (int Lo, int Hi) Band(int n, double p, double alpha)
    {
        var pmf = new double[n + 1];
        for (int x = 0; x <= n; x++)
            pmf[x] = Math.Exp(LogChoose(n, x) + x * Math.Log(p) + (n - x) * Math.Log(1 - p));
        int lo = 0;
        double tail = 0;
        while (lo <= n && tail + pmf[lo] <= alpha / 2) tail += pmf[lo++];
        int hi = n;
        tail = 0;
        while (hi >= 0 && tail + pmf[hi] <= alpha / 2) tail += pmf[hi--];
        return (lo, hi);
    }

    private static double LogChoose(int n, int k)
    {
        double s = 0;
        for (int i = 1; i <= k; i++) s += Math.Log(n - k + i) - Math.Log(i);
        return s;
    }

    public sealed record Cell(int K, string Emulator, double Noise, int Ratio, bool Disc, int N, double MedianErr, double MeanErr, int Covered, int Intervals,
                              double Coverage, string Band, bool InBand, bool InWideBand, double Unconverged, double MedianSd, double ReconTruth,
                              double ReconObs, double FitObs, double Sigma, double Seconds);

    public static List<Cell> Cells(IEnumerable<Row> rows) => rows.Where(r => r.Mode == "nested")
        .GroupBy(r => (r.K, r.Emulator, r.Noise, r.Ratio, r.Disc))
        .OrderBy(g => g.Key.K).ThenBy(g => g.Key.Emulator).ThenBy(g => g.Key.Disc).ThenBy(g => g.Key.Noise).ThenBy(g => g.Key.Ratio)
        .Select(g =>
        {
            int n = g.Count(), intervals = n * g.Key.K, covered = g.Sum(r => r.Covered);
            (int lo, int hi) = Band(intervals, 0.9, 0.01);
            (int wlo, int whi) = Band(intervals, 0.9, 0.001);
            return new Cell(g.Key.K, g.Key.Emulator, g.Key.Noise, g.Key.Ratio, g.Key.Disc, n, Median(g.Select(r => r.Err)), g.Average(r => r.Err), covered, intervals,
                            covered / (double)intervals, $"[{lo},{hi}]", covered >= lo && covered <= hi, covered >= wlo && covered <= whi,
                            g.Count(r => r.Converged < r.Runs) / (double)n, Median(g.Select(r => r.Sd.Average())), Median(g.Select(r => r.ReconTruth)),
                            Median(g.Select(r => r.ReconObs)), Median(g.Select(r => r.FitObs)), Median(g.Select(r => r.Sigma)), Median(g.Select(r => r.Seconds)));
        }).ToList();

    private static string F(double v, string f = "F3") => v.ToString(f, CultureInfo.InvariantCulture);

    public static string Table(IEnumerable<Cell> cells)
    {
        var sb = new StringBuilder("| k | emulator | disc | σ²ₙ | ratio | n | median err | mean err | coverage (band 99%) | unconverged | median sd | recon held (clean) | recon held (noisy) | fit kept | σ post | s/cal |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
        foreach (Cell c in cells)
            sb.AppendLine($"| {c.K} | {c.Emulator.Replace("bass-functional", "bass")} | {(c.Disc ? "on" : "off")} | {c.Noise:0e0} | {c.Ratio} | {c.N} | {F(c.MedianErr)} | {F(c.MeanErr)} | " +
                          $"{c.Covered}/{c.Intervals} = {F(c.Coverage, "F2")} {c.Band}{(c.InBand ? "" : c.InWideBand ? " (marginal)" : " OUT")} | {F(c.Unconverged, "F2")} | {F(c.MedianSd)} | " +
                          $"{F(c.ReconTruth)} | {F(c.ReconObs)} | {F(c.FitObs)} | {F(c.Sigma)} | {F(c.Seconds, "F1")} |");
        return sb.ToString();
    }

    public static string Stage1(string dir, Pilot.Settings s)
    {
        var sb = new StringBuilder();
        var rows = new List<Row>();
        var fits = new List<string>();
        for (int k = 1; k <= 5; k++)
        {
            string path = Path.Combine(dir, $"s1-k{k}-replicates.csv");
            if (File.Exists(path)) rows.AddRange(Load(path));
            string fit = Path.Combine(dir, $"s1-k{k}-fit.jsonl");
            if (File.Exists(fit)) fits.AddRange(File.ReadAllLines(fit));
        }
        sb.AppendLine("## Stage 1 gate readings\n\n### Fits\n");
        foreach (string f in fits) sb.AppendLine("- " + f);
        List<Cell> cells = Cells(rows);
        sb.AppendLine("\n### Cells (nested sampler)\n");
        sb.AppendLine(Table(cells));

        double[] lowNoise = s.NoiseVariances[..2];
        sb.AppendLine("### R1 and R2\n");
        bool r1 = true, r2 = true, marginal = false;
        foreach (int k in cells.Select(c => c.K).Distinct())
            foreach (double v in lowNoise)
            {
                var line = cells.Where(c => c.K == k && c.Emulator == "bass-functional" && !c.Disc && c.Noise == v).OrderBy(c => c.Ratio).ToList();
                int? first = null;
                foreach (Cell c in line)
                    if (line.Where(d => d.Ratio >= c.Ratio).All(d => d.MedianErr <= 0.05)) { first = c.Ratio; break; }
                bool ok = first is <= 4;
                bool near = !ok && line.Where(c => c.Ratio <= 4).Any(c => line.Where(d => d.Ratio >= c.Ratio).All(d => d.MedianErr <= 0.06));
                r1 &= ok;
                marginal |= near;
                var counted = first is { } f0 ? line.Where(c => c.Ratio >= f0).ToList() : [];
                bool cover = counted.All(c => c.InBand);
                bool coverNear = !cover && counted.All(c => c.InWideBand);
                if (ok) { r2 &= cover; marginal |= coverNear; }
                sb.AppendLine($"- k={k} σ²ₙ={v:0e0}: R1 {(ok ? "pass" : near ? "MARGINAL" : "FAIL")} (settles at ratio {first?.ToString() ?? "none"}); " +
                              $"R2 {(ok ? (cover ? "pass" : coverNear ? "MARGINAL" : "FAIL") : "not counted")} ({string.Join(", ", counted.Select(c => $"r{c.Ratio} {F(c.Coverage, "F2")}"))})");
            }
        sb.AppendLine($"\nR1 {(r1 ? "PASS" : "FAIL")}, R2 {(r2 ? "PASS" : "FAIL")}{(marginal ? " (a marginal criterion: rerun rule applies)" : "")}\n");

        sb.AppendLine("### R3, modality\n");
        foreach (Cell c in cells.Where(c => c.Unconverged > 0.2))
            sb.AppendLine($"- flagged by convergence: k={c.K} {c.Emulator} disc={(c.Disc ? 1 : 0)} σ²ₙ={c.Noise:0e0} ratio {c.Ratio}: {F(c.Unconverged, "F2")} of replicates with an unconverged inner run");
        var nested = rows.Where(r => r.Mode == "nested" && !r.Disc).ToDictionary(r => (r.K, r.Rep, r.Emulator, r.Noise, r.Ratio));
        foreach (var g in rows.Where(r => r.Mode != "nested").GroupBy(r => (r.K, r.Emulator, r.Mode, r.Ratio)).OrderBy(g => g.Key.K).ThenBy(g => g.Key.Emulator).ThenBy(g => g.Key.Ratio))
        {
            int moved = 0;
            var shifts = new List<double>();
            foreach (Row r in g)
            {
                Row b = nested[(r.K, r.Rep, r.Emulator, r.Noise, r.Ratio)];
                double worst = Enumerable.Range(0, r.K).Max(j => Math.Abs(r.Mean[j] - b.Mean[j]) / Math.Max(b.Sd[j], 1e-9));
                shifts.Add(worst);
                if (worst > 2) moved++;
            }
            sb.AppendLine($"- {g.Key.Mode} vs nested, k={g.Key.K} {g.Key.Emulator} ratio {g.Key.Ratio}: {moved}/{g.Count()} replicates moved > 2 sd " +
                          $"(max shift {F(shifts.Max(), "F2")} sd; median {F(Median(shifts), "F2")}); {(moved > 2 ? "FLAG" : "no flag")}; " +
                          $"median err {F(Median(g.Select(r => r.Err)))} vs nested {F(Median(g.Select(r => nested[(r.K, r.Rep, r.Emulator, r.Noise, r.Ratio)].Err)))}");
        }

        sb.AppendLine("\n### R4, emulator misfit\n");
        foreach (Cell b in cells.Where(c => c.Emulator == "bass-functional" && !c.Disc))
        {
            Cell? o = cells.FirstOrDefault(c => c.Emulator == "oracle" && c.K == b.K && c.Noise == b.Noise && c.Ratio == b.Ratio);
            if (o is null) continue;
            bool fires = b.MedianErr > 2 * o.MedianErr + 0.01 || (!b.InBand && b.Coverage < 0.9 && o.InBand);
            if (fires)
                sb.AppendLine($"- fires: k={b.K} σ²ₙ={b.Noise:0e0} ratio {b.Ratio}: bass err {F(b.MedianErr)} cover {F(b.Coverage, "F2")} vs reference err {F(o.MedianErr)} cover {F(o.Coverage, "F2")}" +
                              $"{(o.MedianErr <= 0.05 ? " (reference meets R1's threshold)" : "")}");
        }

        sb.AppendLine("\n### R5, cost\n");
        foreach (var g in rows.GroupBy(r => r.K))
            sb.AppendLine($"- k={g.Key}: median calibration {F(Median(g.Where(r => r.Mode == "nested").Select(r => r.Seconds)), "F1")} s (bass {F(Median(g.Where(r => r.Mode == "nested" && r.Emulator != "oracle").Select(r => r.Seconds)), "F1")} s, " +
                          $"reference {F(Median(g.Where(r => r.Mode == "nested" && r.Emulator == "oracle").Select(r => r.Seconds)), "F1")} s), total {F(g.Sum(r => r.Seconds) / 3600, "F2")} CPU-h over {g.Count()} calibrations");
        return sb.ToString();
    }
}
