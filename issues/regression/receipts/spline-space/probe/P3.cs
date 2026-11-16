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
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using static SplineSpaceProbe.SplineTools;

namespace SplineSpaceProbe;

/// <summary>
/// The damped oscillator's one-sided L-spline basis on [0, 1] (Schumaker § 10.3; the Lean chip's S8): e^{−γx}cos ωx,
/// e^{−γx}sin ωx, and ω·g(x − ξᵢ) at each interior knot, g(u) = e^{−γu}·sin(ωu)/ω for u ≥ 0 and 0 before. The columns
/// are scaled to unit amplitude (ω·g); scale does not change the g-prior fit. Probe code behind the engine's IBasis.
/// </summary>
public sealed class OscillatorBasis(double gamma, double omega) : IBasis
{
    public double Gamma { get; } = gamma;
    public double Omega { get; } = omega;

    public int Dimension(int knotCount) => knotCount + 2;

    public double Column(int j, double[] knots, double x)
    {
        if (j == 0) return Math.Exp(-Gamma * x) * Math.Cos(Omega * x);
        if (j == 1) return Math.Exp(-Gamma * x) * Math.Sin(Omega * x);
        double u = x - knots[j - 2];
        return u < 0 ? 0.0 : Math.Exp(-Gamma * u) * Math.Sin(Omega * u);
    }

    public double[,] Design(KnotConfig config, double[] x)
    {
        double[] knots = config.InteriorKnots;
        int nu = knots.Length + 2;
        var z = new double[x.Length, nu];
        for (int i = 0; i < x.Length; i++) for (int j = 0; j < nu; j++) z[i, j] = Column(j, knots, x[i]);
        return z;
    }

    public double Evaluate(KnotConfig config, double[] coef, double x)
    {
        double v = 0;
        for (int j = 0; j < coef.Length; j++) v += coef[j] * Column(j, config.InteriorKnots, x);
        return v;
    }
}

/// <summary>P3: the oscillator basis against cubic splines under BARS, and S8's spacing condition.</summary>
public static class P3
{
    private const int Master = 2026111703, N = 1000, Replicates = 5;
    private const double Gamma = 8.0, Omega = 2 * Math.PI * 20, Noise = 0.01;

    private static (double[] Onsets, double[] Amplitudes) Draw(int r, Xoshiro256PlusPlus rng)
    {
        var on = new List<double>();
        while (on.Count < r)
        {
            double c = 0.05 + 0.75 * rng.NextDouble();
            if (on.All(o => Math.Abs(o - c) >= 0.05)) on.Add(c);
        }
        on.Sort();
        return ([.. on], on.Select(_ => (rng.NextDouble() < 0.5 ? -1 : 1) * (0.5 + rng.NextDouble())).ToArray());
    }

    private static double Signal(double x, double[] onsets, double[] amps, double gamma, double omega)
    {
        double v = 0;
        for (int i = 0; i < onsets.Length; i++)
        {
            double u = x - onsets[i];
            if (u >= 0) v += amps[i] * Math.Exp(-gamma * u) * Math.Sin(omega * u);
        }
        return v;
    }

    public static void Run(string dir)
    {
        int[] seeds = SeedTree.Derive(Master, 4);
        var cases = new List<(string Variant, int R, int Rep)>();
        foreach (int r in new[] { 1, 2, 4 }) for (int rep = 0; rep < Replicates; rep++) cases.Add(("matched", r, rep));
        foreach (string v in new[] { "omega x1.05", "gamma x1.5" }) for (int rep = 0; rep < Replicates; rep++) cases.Add((v, 2, rep));
        int[] signalSeeds = SeedTree.Derive(seeds[0], cases.Count), noiseSeeds = SeedTree.Derive(seeds[1], cases.Count), barsSeeds = SeedTree.Derive(seeds[2], cases.Count * 2);
        double[] x = Enumerable.Range(0, N).Select(i => i / (N - 1.0)).ToArray();
        double[] grid = Enumerable.Range(0, 2001).Select(i => i / 2000.0).ToArray();
        var jobs = Enumerable.Range(0, cases.Count).SelectMany(c => new[] { (c, "cubic"), (c, "oscillator") }).ToList();
        var rows = new ConcurrentBag<(int, string)>();
        Parallel.For(0, jobs.Count, new ParallelOptions { MaxDegreeOfParallelism = Program.Threads }, ji =>
        {
            (int c, string method) = jobs[ji];
            (string variant, int r, int rep) = cases[c];
            (double[] onsets, double[] amps) = Draw(r, new Xoshiro256PlusPlus(signalSeeds[c]));
            double g = variant == "gamma x1.5" ? 1.5 * Gamma : Gamma, w = variant == "omega x1.05" ? 1.05 * Omega : Omega;
            var noise = new Xoshiro256PlusPlus(noiseSeeds[c]);
            double[] y = x.Select(v => Signal(v, onsets, amps, g, w) + Noise * Normal.SampleStandard(noise)).ToArray();
            double[] truth = grid.Select(v => Signal(v, onsets, amps, g, w)).ToArray();
            var sw = Stopwatch.StartNew();
            var basis = new OscillatorBasis(Gamma, Omega);
            string oracle = "";
            if (method == "oscillator")
            {
                // Least squares on the oscillator basis with knots at the true onsets.
                var cfg = new KnotConfig(onsets);
                double[] b = LeastSquares(basis.Design(cfg, x), y);
                oracle = F(Rms(Multiply(basis.Design(cfg, grid), b), truth), "E3");
            }
            var config = new BarsConfig
            {
                Basis = method == "cubic" ? null : basis,
                Degree = 3,
                Model = new WeightedNormalModel(GPriorCentre.Zero),
                Prior = method == "cubic" ? new UniformPrior(150) : new UniformPrior(20),
                MinColumnMass = method == "cubic" ? SplineTarget.DefaultMinColumnMass : 0.0,
                Chains = 2,
                MasterSeed = barsSeeds[2 * c + (method == "cubic" ? 0 : 1)],
                BurnIn = 2000,
                MaxSamples = 1000,
            };
            string line;
            try
            {
                BarsResult result = Bars.Run(config, x, y, grid);
                double error = Rms(result.Fit, truth);
                // Oscillator knots near onsets: the share of posterior knots within ±0.01 of a true onset.
                long near = 0, total = 0;
                foreach (var chain in result.Draws)
                    foreach (BarsDraw d in chain)
                        foreach (double kn in d.Knots.InteriorKnots) { total++; if (onsets.Any(o => Math.Abs(kn - o) <= 0.01)) near++; }
                line = string.Join(",", variant, r, rep, method, F(result.MeanKnots, "F2"), F(error, "E3"), oracle, F(total > 0 ? near / (double)total : 0, "F3"),
                                   string.Join(";", result.Evaluations.Entries.Where(e => e.Count > 0).Select(e => $"{e.Category}:{e.Count}")), F(sw.Elapsed.TotalSeconds, "F1"));
            }
            catch (Exception e)
            {
                line = string.Join(",", variant, r, rep, method, "", "", oracle, "", "threw:" + e.GetType().Name + ":" + e.Message.Replace(',', ';').Replace('\n', ' '), F(sw.Elapsed.TotalSeconds, "F1"));
            }
            rows.Add((ji, line));
            Console.WriteLine($"p3 {line}");
        });
        File.WriteAllLines(Path.Combine(dir, "p3-runs.csv"),
            ["variant,R,rep,method,k_mean,rms,oracle_ls_rms,share_near_onsets,evaluations,seconds", .. rows.OrderBy(t => t.Item1).Select(t => t.Item2)]);
    }

    /// <summary>The spacing sweep: one-sided basis conditioning, S8's local element, and its unnormalized form.</summary>
    public static void Sweep(string dir)
    {
        const int knots = 40, sites = 4000;
        double[] ratios = [.. Enumerable.Range(1, 24).Select(i => 0.05 * i), .. Enumerable.Range(0, 9).Select(i => 0.48 + 0.005 * i)];
        ratios = ratios.Distinct().OrderBy(v => v).ToArray();
        var lines = new List<string> { "gamma,h_omega_over_pi,onesided_condition,local_sup,local_negative,local_unnormalized_condition" };
        foreach (double gamma in new[] { 0.0, Gamma })
            foreach (double ratio in ratios)
            {
                double h = ratio * Math.PI / Omega, length = (knots + 1) * h;
                double[] x = Enumerable.Range(0, sites).Select(i => length * i / (sites - 1.0)).ToArray();
                double[] y = Enumerable.Range(0, knots + 2).Select(j => j * h).ToArray();   // knots y_0 = 0 … y_{K+1} = length
                // One-sided basis on [0, length]: the two null-space functions and g at the interior knots y_1 … y_K.
                var z = new double[sites, knots + 2];
                for (int i = 0; i < sites; i++)
                {
                    z[i, 0] = Math.Exp(-gamma * x[i]) * Math.Cos(Omega * x[i]);
                    z[i, 1] = Math.Exp(-gamma * x[i]) * Math.Sin(Omega * x[i]);
                    for (int j = 1; j <= knots; j++)
                    {
                        double u = x[i] - y[j];
                        z[i, j + 1] = u < 0 ? 0 : Math.Exp(-gamma * u) * Math.Sin(Omega * u);
                    }
                }
                double oneSided = EquilibratedCondition(z);
                // S8's local element Q on (y0, y1, y2) = (0, h, 2h), scaled by e^{−γ(x−y0)}; its sup and sign on a fine grid.
                double s1 = Math.Sin(Omega * h), s2 = Math.Sin(Omega * h), s02 = Math.Sin(Omega * 2 * h);
                double sup = 0;
                bool negative = false;
                for (int i = 0; i <= 2000; i++)
                {
                    double t = 2 * h * i / 2000.0;
                    double q = t < h ? Math.Sin(Omega * t) / (s1 * s02) : Math.Sin(Omega * (2 * h - t)) / (s2 * s02);
                    q *= Math.Exp(-gamma * t);
                    if (i > 0 && i < 2000 && q < 0) negative = true;
                    sup = Math.Max(sup, Math.Abs(q));
                }
                // The unnormalized local elements P_j on (y_j, y_{j+1}, y_{j+2}), j = 0 … K − 1, S8's numerators.
                var p = new double[sites, knots];
                for (int i = 0; i < sites; i++)
                    for (int j = 0; j < knots; j++)
                    {
                        double u = x[i] - y[j];
                        double v = u <= 0 || u >= 2 * h ? 0 : u < h ? Math.Sin(Omega * h) * Math.Sin(Omega * u) : Math.Sin(Omega * h) * Math.Sin(Omega * (2 * h - u));
                        p[i, j] = v * Math.Exp(-gamma * u);
                    }
                double local = EquilibratedCondition(p);
                lines.Add(string.Join(",", F(gamma), F(ratio, "F3"), F(oneSided, "E3"), F(sup, "E3"), negative ? 1 : 0, F(local, "E3")));
                Console.WriteLine(lines[^1]);
            }
        File.WriteAllLines(Path.Combine(dir, "p3-sweep.csv"), lines);
    }

    public static string Summary(string dir)
    {
        var sb = new StringBuilder("## P3 readings\n\n| variant | R | method | median k̄ | median rms | median oracle LS rms | share of knots near onsets | runs with failures |\n|---|---|---|---|---|---|---|---|\n");
        var d = File.ReadAllLines(Path.Combine(dir, "p3-runs.csv")).Skip(1).Select(l => l.Split(',')).ToList();
        double V(string s) => s == "" ? double.NaN : double.Parse(s, CultureInfo.InvariantCulture);
        foreach (var g in d.GroupBy(r => (r[0], r[1], r[3])).OrderBy(g => g.Key.Item1).ThenBy(g => g.Key.Item2).ThenBy(g => g.Key.Item3))
            sb.AppendLine($"| {g.Key.Item1} | {g.Key.Item2} | {g.Key.Item3} | {P2.Median(g.Select(r => V(r[4]))):F1} | {P2.Median(g.Select(r => V(r[5]))):E2} | " +
                          $"{P2.Median(g.Select(r => V(r[6]))):E2} | {P2.Median(g.Select(r => V(r[7]))):F2} | {g.Count(r => r[8].StartsWith("threw"))} |");
        sb.AppendLine("\nSpacing sweep: see p3-sweep.csv (hω/π, one-sided condition, Q's sup and sign, unnormalized local condition).");
        return sb.ToString();
    }
}
