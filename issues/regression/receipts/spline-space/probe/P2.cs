using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Maths.Distributions;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using static SplineSpaceProbe.SplineTools;

namespace SplineSpaceProbe;

/// <summary>The count held near K by a steep, finite penalty, log p(j) = −50·|j − K|, so a chain may start anywhere.</summary>
public sealed class HeldCountPrior(int target) : IComplexityPrior
{
    public double LogPrior(int knotCount) => knotCount < 0 ? double.NegativeInfinity : -50.0 * Math.Abs(knotCount - target);
}

/// <summary>P2: uniform-knot least squares against BARS on curves smooth but for a few events.</summary>
public static class P2
{
    private const int Master = 2026111702, N = 2000, GridPoints = 4001, Replicates = 5;
    private const double Noise = 1e-3, Window = 0.02;
    private static readonly int[] LsKnots = [2, 4, 8, 16, 32, 64], HeldKnots = [2, 4, 8, 16, 32];

    public sealed record Curve(string Family, int E, double[] Locations, double[] Amplitudes)
    {
        public double Background(double x) => Math.Sin(2 * Math.PI * x) + 0.5 * Math.Cos(5 * Math.PI * x);
        public double FourthDerivative(double x) => Math.Pow(2 * Math.PI, 4) * Math.Sin(2 * Math.PI * x) + 0.5 * Math.Pow(5 * Math.PI, 4) * Math.Cos(5 * Math.PI * x);
        public double Value(double x)
        {
            double v = Background(x);
            for (int e = 0; e < Locations.Length; e++)
            {
                double u = x - Locations[e];
                v += Family == "kink" ? Amplitudes[e] * Math.Abs(u) : Amplitudes[e] * (u > 0 ? u * u : 0.0);
            }
            return v;
        }
    }

    private static Curve Draw(string family, int e, Xoshiro256PlusPlus rng)
    {
        var loc = new List<double>();
        while (loc.Count < e)
        {
            double c = 0.1 + 0.8 * rng.NextDouble();
            if (loc.All(l => Math.Abs(l - c) >= 0.05)) loc.Add(c);
        }
        loc.Sort();
        double[] amp = loc.Select(_ => (rng.NextDouble() < 0.5 ? -1 : 1) * (0.5 + 0.5 * rng.NextDouble())).ToArray();
        return new Curve(family, e, [.. loc], amp);
    }

    public static void Run(string dir, bool rerun)
    {
        int[] seeds = SeedTree.Derive(Master, 8);
        (string Family, int E)[] families = [("kink", 0), ("kink", 1), ("kink", 2), ("kink", 4), ("kink", 8), ("f2jump", 2)];
        // Reruns draw fresh curves, noise and chains from the rerun child.
        int[] roots = rerun ? SeedTree.Derive(seeds[5], 4) : [seeds[0], seeds[1], seeds[3], seeds[4]];
        int[] curveSeeds = SeedTree.Derive(roots[0], families.Length * Replicates);
        int[] noiseSeeds = SeedTree.Derive(roots[1], families.Length * Replicates);
        int[] freeSeeds = SeedTree.Derive(roots[2], families.Length * Replicates);
        int[] heldSeeds = SeedTree.Derive(roots[3], families.Length * Replicates * HeldKnots.Length);
        double[] x = Enumerable.Range(0, N).Select(i => i / (N - 1.0)).ToArray();
        double[] grid = Enumerable.Range(0, GridPoints).Select(i => i / (GridPoints - 1.0)).ToArray();
        var jobs = new List<(int F, int Rep, string Method, int K)>();
        for (int f = 0; f < families.Length; f++)
            for (int r = 0; r < Replicates; r++)
            {
                jobs.Add((f, r, "ls", 0));
                jobs.Add((f, r, "bars-free", 0));
                foreach (int k in HeldKnots) jobs.Add((f, r, "bars-held", k));
            }
        var rows = new ConcurrentBag<(int, string)>();
        Parallel.For(0, jobs.Count, new ParallelOptions { MaxDegreeOfParallelism = Program.Threads }, ji =>
        {
            (int f, int r, string method, int kTarget) = jobs[ji];
            int idx = f * Replicates + r;
            Curve curve = Draw(families[f].Family, families[f].E, new Xoshiro256PlusPlus(curveSeeds[idx]));
            var noise = new Xoshiro256PlusPlus(noiseSeeds[idx]);
            double[] y = x.Select(v => curve.Value(v) + Noise * Normal.SampleStandard(noise)).ToArray();
            double[] truth = grid.Select(curve.Value).ToArray();
            string ev = string.Join(";", curve.Locations.Select(l => F(l, "F4")));
            var sw = Stopwatch.StartNew();
            if (method == "ls")
            {
                var basis = new SplineBasis(3);
                foreach (int k in LsKnots)
                {
                    var cfg = new KnotConfig(Enumerable.Range(1, k).Select(i => i / (k + 1.0)).ToArray());
                    double[] b = LeastSquares(basis.Design(cfg, x), y);
                    double err = Rms(Multiply(basis.Design(cfg, grid), b), truth);
                    rows.Add((ji * 100 + k, string.Join(",", curve.Family, curve.E, r, "ls-uniform", k, k, F(err, "E4"), "", "", "", F(sw.Elapsed.TotalSeconds, "F2"), ev)));
                }
                return;
            }
            var config = new BarsConfig
            {
                Degree = 3,
                Prior = method == "bars-free" ? new UniformPrior(80) : new HeldCountPrior(kTarget),
                Chains = 2,
                MasterSeed = method == "bars-free" ? freeSeeds[idx] : heldSeeds[idx * HeldKnots.Length + Array.IndexOf(HeldKnots, kTarget)],
                BurnIn = 2000,
                MaxSamples = 1000,
            };
            BarsResult result = Bars.Run(config, x, y, grid);
            double error = Rms(result.Fit, truth);
            // Knot placement: the share within ±Window of an event, knots per event, and the off-event histogram's
            // total-variation distance to the density ∝ |b⁽⁴⁾|^{2/9}.
            long inWindow = 0, total = 0;
            const int bins = 100;
            var hist = new double[bins];
            int draws = 0;
            foreach (var chain in result.Draws)
                foreach (BarsDraw d in chain)
                {
                    draws++;
                    foreach (double kn in d.Knots.InteriorKnots)
                    {
                        total++;
                        if (curve.Locations.Any(l => Math.Abs(kn - l) <= Window)) inWindow++;
                        else hist[Math.Min(bins - 1, (int)(kn * bins))]++;
                    }
                }
            var theory = new double[bins];
            for (int bI = 0; bI < bins; bI++)
            {
                double c = (bI + 0.5) / bins;
                bool masked = curve.Locations.Any(l => Math.Abs(c - l) <= Window);
                theory[bI] = masked ? 0 : Math.Pow(Math.Abs(curve.FourthDerivative(c)), 2.0 / 9.0);
                if (masked) hist[bI] = 0;
            }
            double hs = hist.Sum(), ts = theory.Sum();
            double tv = hs > 0 ? 0.5 * hist.Zip(theory, (a, b) => Math.Abs(a / hs - b / ts)).Sum() : double.NaN;
            double windowShare = Math.Min(1.0, curve.Locations.Length * 2 * Window);
            rows.Add((ji * 100, string.Join(",", curve.Family, curve.E, r, method, kTarget, F(result.MeanKnots, "F2"), F(error, "E4"),
                F(total > 0 ? inWindow / (double)total : 0, "F3"), F(windowShare, "F3"), F(curve.E > 0 ? inWindow / (double)draws / curve.E : 0, "F2") + ";" + F(tv, "F3"),
                F(sw.Elapsed.TotalSeconds, "F1"), ev)));
            Console.WriteLine($"p2 {curve.Family} E={curve.E} rep {r} {method} {kTarget}: k̄ {result.MeanKnots:F1}, err {error:E2}, {sw.Elapsed.TotalSeconds:F0} s");
        });
        File.WriteAllLines(Path.Combine(dir, rerun ? "p2-runs-rerun.csv" : "p2-runs.csv"),
            ["family,E,rep,method,k_target,k_mean,rms,share_in_event_windows,window_share_of_domain,knots_per_event;tv_off_event,seconds,events",
             .. rows.OrderBy(t => t.Item1).Select(t => t.Item2)]);
    }

    /// <summary>
    /// Post-gate: on replicate 0 of the E = 0 and E = 4 kink curves, a BARS run held at 32 knots (the gate's own seeds),
    /// with least squares refitted at each retained draw's knots, against the posterior mean and uniform least squares:
    /// whether the free knots are placed well and the error comes from averaging, or the knots themselves are the gap.
    /// </summary>
    public static void Diagnose(string dir)
    {
        int[] seeds = SeedTree.Derive(Master, 8);
        (string Family, int E)[] families = [("kink", 0), ("kink", 1), ("kink", 2), ("kink", 4), ("kink", 8), ("f2jump", 2)];
        int[] curveSeeds = SeedTree.Derive(seeds[0], families.Length * Replicates);
        int[] noiseSeeds = SeedTree.Derive(seeds[1], families.Length * Replicates);
        int[] heldSeeds = SeedTree.Derive(seeds[4], families.Length * Replicates * HeldKnots.Length);
        double[] x = Enumerable.Range(0, N).Select(i => i / (N - 1.0)).ToArray();
        double[] grid = Enumerable.Range(0, GridPoints).Select(i => i / (GridPoints - 1.0)).ToArray();
        var lines = new List<string> { "family,E,k,posterior_mean_rms,ls_at_draw_knots_median,ls_at_draw_knots_best,uniform_ls_rms,draw_knots_within_0.02_of_events_per_event" };
        foreach (int f in new[] { 0, 3 })
        {
            int idx = f * Replicates;
            Curve curve = Draw(families[f].Family, families[f].E, new Xoshiro256PlusPlus(curveSeeds[idx]));
            var noise = new Xoshiro256PlusPlus(noiseSeeds[idx]);
            double[] y = x.Select(v => curve.Value(v) + Noise * Normal.SampleStandard(noise)).ToArray();
            double[] truth = grid.Select(curve.Value).ToArray();
            const int k = 32;
            var config = new BarsConfig { Degree = 3, Prior = new HeldCountPrior(k), Chains = 2, MasterSeed = heldSeeds[idx * HeldKnots.Length + Array.IndexOf(HeldKnots, k)], BurnIn = 2000, MaxSamples = 1000 };
            BarsResult result = Bars.Run(config, x, y, grid);
            var basis = new SplineBasis(3);
            var errs = new List<double>();
            double perEvent = 0;
            int used = 0;
            foreach (var chain in result.Draws)
                for (int i = 0; i < chain.Count; i += 50)   // every 50th retained draw
                {
                    KnotConfig kc = chain[i].Knots;
                    double[] b = LeastSquares(basis.Design(kc, x), y);
                    errs.Add(Rms(Multiply(basis.Design(kc, grid), b), truth));
                    if (curve.E > 0) perEvent += kc.InteriorKnots.Count(kn => curve.Locations.Any(l => Math.Abs(kn - l) <= Window)) / (double)curve.E;
                    used++;
                }
            var uni = new KnotConfig(Enumerable.Range(1, k).Select(i => i / (k + 1.0)).ToArray());
            double uniform = Rms(Multiply(basis.Design(uni, grid), LeastSquares(basis.Design(uni, x), y)), truth);
            lines.Add(string.Join(",", curve.Family, curve.E, k, F(Rms(result.Fit, truth), "E3"), F(Median(errs), "E3"), F(errs.Min(), "E3"), F(uniform, "E3"), F(perEvent / used, "F2")));
            Console.WriteLine(lines[^1]);
        }
        File.WriteAllLines(Path.Combine(dir, "p2-diagnostic.csv"), lines);
    }

    public static string Summary(string dir, string file = "p2-runs.csv")
    {
        var d = File.ReadAllLines(Path.Combine(dir, file)).Skip(1).Select(l => l.Split(',')).ToList();
        double V(string s) => double.Parse(s, CultureInfo.InvariantCulture);
        var sb = new StringBuilder($"## P2 readings ({file})\n\n| family | E | method | k | median k̄ | median rms |\n|---|---|---|---|---|---|\n");
        foreach (var g in d.GroupBy(r => (r[0], r[1], r[3], r[4])).OrderBy(g => g.Key.Item1).ThenBy(g => int.Parse(g.Key.Item2)).ThenBy(g => g.Key.Item3).ThenBy(g => int.Parse(g.Key.Item4)))
            sb.AppendLine($"| {g.Key.Item1} | {g.Key.Item2} | {g.Key.Item3} | {g.Key.Item4} | {Median(g.Select(r => V(r[5]))):F1} | {Median(g.Select(r => V(r[6]))):E2} |");
        sb.AppendLine("\n### Slopes of ln rms on ln k over k ∈ [4, 32] (medians over replicates)\n");
        sb.AppendLine("| family | E | uniform LS | BARS held | gap (uniform − held) | verdict (E ≥ 1 kinks: gap ≥ 1.0) |\n|---|---|---|---|---|---|");
        foreach (var g in d.GroupBy(r => (r[0], r[1])).OrderBy(g => g.Key.Item1).ThenBy(g => int.Parse(g.Key.Item2)))
        {
            double SlopeOf(string method, bool useMean)
            {
                var pts = g.Where(r => r[3] == method && int.Parse(r[4]) >= 4 && int.Parse(r[4]) <= 32).GroupBy(r => int.Parse(r[4]))
                    .Select(h => (K: useMean ? Median(h.Select(r => V(r[5]))) : h.Key, E: Median(h.Select(r => V(r[6]))))).OrderBy(p => p.K).ToList();
                return pts.Count < 2 ? double.NaN : Slope(pts.Select(p => Math.Log(p.K)).ToList(), pts.Select(p => Math.Log(p.E)).ToList());
            }
            double ls = SlopeOf("ls-uniform", false), held = SlopeOf("bars-held", true), gap = ls - held;
            string verdict = g.Key.Item1 == "kink" && int.Parse(g.Key.Item2) >= 1 ? (gap >= 1.0 ? "pass" : gap >= 0.8 ? "MARGINAL" : "FAIL") : "reading";
            sb.AppendLine($"| {g.Key.Item1} | {g.Key.Item2} | {ls:F2} | {held:F2} | {gap:F2} | {verdict} |");
        }
        sb.AppendLine("\n### Free-count runs: knots against events\n");
        sb.AppendLine("| family | E | median k̄ | share of knots in event windows | windows' share of domain | knots per event | TV off-event vs |b⁽⁴⁾|^{2/9} |\n|---|---|---|---|---|---|---|");
        var free = d.Where(r => r[3] == "bars-free").ToList();
        foreach (var g in free.GroupBy(r => (r[0], r[1])).OrderBy(g => g.Key.Item1).ThenBy(g => int.Parse(g.Key.Item2)))
        {
            var kpe = g.Select(r => V(r[9].Split(';')[0])).ToList();
            var tv = g.Select(r => V(r[9].Split(';')[1])).ToList();
            sb.AppendLine($"| {g.Key.Item1} | {g.Key.Item2} | {Median(g.Select(r => V(r[5]))):F1} | {Median(g.Select(r => V(r[7]))):F3} | {Median(g.Select(r => V(r[8]))):F3} | {Median(kpe):F2} | {Median(tv):F3} |");
        }
        var kinks = free.Where(r => r[0] == "kink").ToList();
        var ex = kinks.Select(r => V(r[1])).ToList();
        var ky = kinks.Select(r => V(r[5])).ToList();
        double slope = Slope(ex, ky), icpt = ky.Average() - slope * ex.Average();
        double resid = Math.Sqrt(kinks.Select(r => Math.Pow(V(r[5]) - icpt - slope * V(r[1]), 2)).Sum() / (kinks.Count - 2));
        sb.AppendLine($"\nk̄ against E (kink families, {kinks.Count} runs): slope {slope:F2} knots per event, intercept {icpt:F1}, residual sd {resid:F2} knots, so E read from k̄ carries an sd of {resid / slope:F2} events.");
        return sb.ToString();
    }

    public static double Median(IEnumerable<double> v)
    {
        var a = v.Where(x => !double.IsNaN(x)).OrderBy(x => x).ToArray();
        return a.Length == 0 ? double.NaN : a.Length % 2 == 1 ? a[a.Length / 2] : 0.5 * (a[a.Length / 2 - 1] + a[a.Length / 2]);
    }
}
