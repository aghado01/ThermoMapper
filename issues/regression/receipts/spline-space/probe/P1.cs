using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using static SplineSpaceProbe.SplineTools;

namespace SplineSpaceProbe;

/// <summary>P1: the greedy ordered matching, the numerical rank and the engine's admission on seeded spline designs.</summary>
public static class P1
{
    private const int Master = 2026111701, Draws = 20, ProductDraws = 20;

    private static double[] Knots(Xoshiro256PlusPlus rng, int k, string law)
    {
        var xi = new double[k];
        if (law == "uniform" || k < 2)
            for (int i = 0; i < k; i++) xi[i] = rng.NextDouble();
        else
        {
            double centre = 0.01 + 0.98 * rng.NextDouble();
            int clustered = k / 2;
            for (int i = 0; i < k; i++)
                xi[i] = i < clustered ? Math.Clamp(centre + 0.02 * (rng.NextDouble() - 0.5), 1e-9, 1 - 1e-9) : rng.NextDouble();
        }
        Array.Sort(xi);
        for (int i = 1; i < k; i++) if (xi[i] <= xi[i - 1]) xi[i] = Math.BitIncrement(xi[i - 1]);   // strictly increasing
        return xi;
    }

    private static double[] Sites(Xoshiro256PlusPlus rng, int count, string law)
    {
        var x = new double[count];
        switch (law)
        {
            case "uniform":
                for (int i = 0; i < count; i++) x[i] = rng.NextDouble();
                break;
            case "clustered":
            {
                int clusters = 1 + rng.NextInt(3);
                double[] centres = Enumerable.Range(0, clusters).Select(_ => 0.025 + 0.95 * rng.NextDouble()).ToArray();
                for (int i = 0; i < count; i++) x[i] = Math.Clamp(centres[rng.NextInt(clusters)] + 0.05 * (rng.NextDouble() - 0.5), 0, 1);
                break;
            }
            default:   // gapped: uniform on [0, 1] without a random interval of width 0.2
            {
                double start = 0.8 * rng.NextDouble();
                for (int i = 0; i < count; i++)
                {
                    double u = 0.8 * rng.NextDouble();
                    x[i] = u < start ? u : u + 0.2;
                }
                break;
            }
        }
        Array.Sort(x);
        return x;
    }

    public static void Run(string dir)
    {
        int[] seeds = SeedTree.Derive(Master, 4);
        var cells = (from d in new[] { 1, 2, 3 } from k in Enumerable.Range(0, 13) from nOff in new[] { -1, 0, 1, 3, -2 }
                     from kl in new[] { "uniform", "clustered" } from sl in new[] { "uniform", "clustered", "gapped" }
                     select (d, k, nOff, kl, sl)).ToList();
        int[] cellSeeds = SeedTree.Derive(seeds[0], cells.Count);
        int[] responseSeeds = SeedTree.Derive(seeds[3], cells.Count * Draws);   // the spare child: response values only
        var rows = new ConcurrentBag<(int, string)>();
        double worstEngineDiff = 0;
        object gate = new();
        Parallel.For(0, cells.Count, new ParallelOptions { MaxDegreeOfParallelism = Program.Threads }, ci =>
        {
            (int d, int k, int nOff, string kl, string sl) = cells[ci];
            int order = d + 1, n = k + order;
            int count = nOff == -2 ? 2 * n : n + nOff;
            var rng = new Xoshiro256PlusPlus(cellSeeds[ci]);
            var engineBasis = new SplineBasis(d);
            for (int draw = 0; draw < Draws; draw++)
            {
                double[] xi = Knots(rng, k, kl);
                double[] x = Sites(rng, count, sl);
                double[] t = Clamped(xi, order);
                double[,] z = Design(t, order, x);
                double[,] ze = engineBasis.Design(new KnotConfig(xi), x);
                double diff = 0;
                for (int i = 0; i < count; i++) for (int j = 0; j < n; j++) diff = Math.Max(diff, Math.Abs(z[i, j] - ze[i, j]));
                lock (gate) worstEngineDiff = Math.Max(worstEngineDiff, diff);
                Support[] supports = Supports(t, order);
                bool matched = Matches(x, supports);
                double[] sv = SingularValues(z);
                double rho = sv[0] > 0 ? sv[^1] / sv[0] : 0;
                // The recheck design: every entry outside its support set to exactly zero.
                double rhoZeroed = double.NaN;
                if (!matched && rho > 1e-13)
                {
                    var zz = (double[,])z.Clone();
                    for (int i = 0; i < count; i++) for (int j = 0; j < n; j++) if (!supports[j].Contains(x[i])) zz[i, j] = 0;
                    double[] s2 = SingularValues(zz);
                    rhoZeroed = s2[0] > 0 ? s2[^1] / s2[0] : 0;
                }
                // The engine's verdict: SplineTarget's evaluation, its rejection counted into a fresh tally.
                var y = new double[count];
                var noise = new Xoshiro256PlusPlus(responseSeeds[ci * Draws + draw]);   // the Normal model's Gram ignores the response
                for (int i = 0; i < count; i++) y[i] = Normal.SampleStandard(noise);
                var tally = new Tally<DesignRejection>();
                string verdict;
                try
                {
                    var target = new SplineTarget(engineBasis, new WeightedNormalModel(), new UniformPrior(12), x, y, outcomes: tally);
                    target.Evaluate(new KnotConfig(xi));
                    verdict = tally.Entries.Where(e => e.Count > 0).Select(e => e.Category.ToString()).FirstOrDefault() ?? "none-counted";
                }
                catch (Exception e) { verdict = "threw:" + e.GetType().Name; }
                rows.Add((ci * Draws + draw, string.Join(",", d, k, n, count, kl, sl, draw, matched ? 1 : 0, F(rho, "E3"), F(rhoZeroed, "E3"), verdict, F(diff, "E2"))));
            }
        });
        File.WriteAllLines(Path.Combine(dir, "p1-designs.csv"),
            ["degree,k,n,N,knots,sites,draw,matched,rho,rho_zeroed,engine,engine_design_diff", .. rows.OrderBy(r => r.Item1).Select(r => r.Item2)]);
        Console.WriteLine($"p1: {rows.Count} designs; worst |probe − engine design| {worstEngineDiff:E2}");
        Products(dir, seeds[1]);
    }

    private static void Products(string dir, int seed)
    {
        var cells = (from kx in Enumerable.Range(0, 5) from ky in Enumerable.Range(0, 5) from mult in new[] { 1, 2, 4 }
                     from sl in new[] { "uniform", "diagonal", "near-diagonal", "clustered" } select (kx, ky, mult, sl)).ToList();
        int[] cellSeeds = SeedTree.Derive(seed, cells.Count);
        var rows = new ConcurrentBag<(int, string)>();
        Parallel.For(0, cells.Count, new ParallelOptions { MaxDegreeOfParallelism = Program.Threads }, ci =>
        {
            (int kx, int ky, int mult, string sl) = cells[ci];
            const int order = 4;
            var rng = new Xoshiro256PlusPlus(cellSeeds[ci]);
            for (int draw = 0; draw < ProductDraws; draw++)
            {
                double[] xiX = Knots(rng, kx, "uniform"), xiY = Knots(rng, ky, "uniform");
                double[] tx = Clamped(xiX, order), ty = Clamped(xiY, order);
                int nx = tx.Length - order, ny = ty.Length - order, n = nx * ny, count = mult * n;
                var px = new double[count];
                var py = new double[count];
                int clusters = 1 + rng.NextInt(3);
                var centres = Enumerable.Range(0, clusters).Select(_ => (0.05 + 0.9 * rng.NextDouble(), 0.05 + 0.9 * rng.NextDouble())).ToArray();
                for (int i = 0; i < count; i++)
                {
                    double u = rng.NextDouble(), v = rng.NextDouble();
                    (px[i], py[i]) = sl switch
                    {
                        "uniform" => (u, v),
                        "diagonal" => (u, u),
                        "near-diagonal" => (u, Math.Clamp(u + 0.04 * (v - 0.5), 0, 1)),
                        _ => (Math.Clamp(centres[i % clusters].Item1 + 0.1 * (u - 0.5), 0, 1), Math.Clamp(centres[i % clusters].Item2 + 0.1 * (v - 0.5), 0, 1)),
                    };
                }
                Support[] sx = Supports(tx, order), sy = Supports(ty, order);
                bool mx = Matches(px.Distinct().OrderBy(v => v).ToArray(), sx);
                bool my = Matches(py.Distinct().OrderBy(v => v).ToArray(), sy);
                var z = new double[count, n];
                var adjacent = new bool[count, n];
                for (int i = 0; i < count; i++)
                {
                    double[] bx = Basis(tx, order, px[i]), by = Basis(ty, order, py[i]);
                    for (int a = 0; a < nx; a++)
                        for (int b = 0; b < ny; b++)
                        {
                            z[i, a * ny + b] = bx[a] * by[b];
                            adjacent[i, a * ny + b] = sx[a].Contains(px[i]) && sy[b].Contains(py[i]);
                        }
                }
                bool bipartite = MaximumMatching(adjacent) == n;
                double[] sv = SingularValues(z);
                double rho = sv[0] > 0 ? sv[^1] / sv[0] : 0;
                rows.Add((ci * ProductDraws + draw, string.Join(",", kx, ky, n, count, sl, draw, mx ? 1 : 0, my ? 1 : 0, bipartite ? 1 : 0, F(rho, "E3"))));
            }
        });
        File.WriteAllLines(Path.Combine(dir, "p1-products.csv"),
            ["kx,ky,n,N,sites,draw,mx,my,bipartite,rho", .. rows.OrderBy(r => r.Item1).Select(r => r.Item2)]);
        Console.WriteLine($"p1 products: {rows.Count} designs");
    }

    public static string Summary(string dir)
    {
        var sb = new StringBuilder("## P1 readings\n\n");
        string[] lines = File.ReadAllLines(Path.Combine(dir, "p1-designs.csv"));
        var d = lines.Skip(1).Select(l => l.Split(',')).ToList();
        double Rho(string[] r) => double.Parse(r[8], System.Globalization.CultureInfo.InvariantCulture);
        double RhoZ(string[] r) => double.Parse(r[9], System.Globalization.CultureInfo.InvariantCulture);
        var unmatched = d.Where(r => r[7] == "0").ToList();
        var matched = d.Where(r => r[7] == "1").ToList();
        var violations = unmatched.Where(r => Rho(r) > 1e-13).ToList();
        sb.AppendLine($"- designs {d.Count}: matched {matched.Count}, unmatched {unmatched.Count}; worst probe−engine design difference {d.Max(r => double.Parse(r[11], System.Globalization.CultureInfo.InvariantCulture)):E2}");
        sb.AppendLine($"- (i) unmatched with rho > 1e-13: {violations.Count}; of them cleared by the zeroed recheck (rho ≤ 1e-13): {violations.Count(r => RhoZ(r) <= 1e-13)}; max rho among unmatched {unmatched.Select(Rho).DefaultIfEmpty(0).Max():E2}");
        foreach (var r in violations.Take(10)) sb.AppendLine($"    violation: {string.Join(",", r)}");
        var mr = matched.Select(Rho).OrderBy(v => v).ToList();
        sb.AppendLine($"- (ii) matched rho: min {mr[0]:E2}, 1% {mr[(int)(0.01 * mr.Count)]:E2}, median {mr[mr.Count / 2]:E2}; below 1e-12: {mr.Count(v => v < 1e-12)}, below 1e-10: {mr.Count(v => v < 1e-10)}");
        sb.AppendLine("\n| matching | engine verdict | designs | of them N > n |\n|---|---|---|---|");
        foreach (var g in d.GroupBy(r => (r[7], r[10])).OrderBy(g => g.Key))
            sb.AppendLine($"| {(g.Key.Item1 == "1" ? "true" : "false")} | {g.Key.Item2} | {g.Count()} | {g.Count(r => int.Parse(r[3]) > int.Parse(r[2]))} |");
        sb.AppendLine("\nMatched, N > n, refused by the engine, by reason and rho band:\n");
        foreach (var g in matched.Where(r => int.Parse(r[3]) > int.Parse(r[2]) && r[10] != "None").GroupBy(r => r[10]))
        {
            var rr = g.Select(Rho).ToList();
            sb.AppendLine($"- {g.Key}: {g.Count()}; rho < 1e-10: {rr.Count(v => v < 1e-10)}, 1e-10..1e-6: {rr.Count(v => v >= 1e-10 && v < 1e-6)}, ≥ 1e-6: {rr.Count(v => v >= 1e-6)}; by site law {string.Join(", ", g.GroupBy(r => r[5]).Select(h => $"{h.Key} {h.Count()}"))}");
        }
        string[] pl = File.ReadAllLines(Path.Combine(dir, "p1-products.csv"));
        var p = pl.Skip(1).Select(l => l.Split(',')).ToList();
        bool Full(string[] r) => double.Parse(r[9], System.Globalization.CultureInfo.InvariantCulture) > 1e-13;
        sb.AppendLine($"\n### Products ({p.Count} designs)\n");
        sb.AppendLine($"- full rank with Mx or My false (per-axis necessity): {p.Count(r => Full(r) && (r[6] == "0" || r[7] == "0"))}");
        sb.AppendLine($"- Mx and My true but rank-deficient (per-axis sufficiency fails): {p.Count(r => !Full(r) && r[6] == "1" && r[7] == "1")} of {p.Count(r => r[6] == "1" && r[7] == "1")} with both matched");
        sb.AppendLine($"- full rank without a bipartite matching (its necessity): {p.Count(r => Full(r) && r[8] == "0")}");
        sb.AppendLine($"- bipartite matching but rank-deficient (its sufficiency fails): {p.Count(r => !Full(r) && r[8] == "1")} of {p.Count(r => r[8] == "1")} with a matching");
        sb.AppendLine("\n| sites | designs | full rank | Mx∧My | bipartite | Mx∧My ∧ deficient | bipartite ∧ deficient |\n|---|---|---|---|---|---|---|");
        foreach (var g in p.GroupBy(r => r[4]))
            sb.AppendLine($"| {g.Key} | {g.Count()} | {g.Count(Full)} | {g.Count(r => r[6] == "1" && r[7] == "1")} | {g.Count(r => r[8] == "1")} | {g.Count(r => !Full(r) && r[6] == "1" && r[7] == "1")} | {g.Count(r => !Full(r) && r[8] == "1")} |");
        sb.AppendLine("\nKept counterexamples (bipartite matching, rank-deficient; first 8):\n");
        foreach (var r in p.Where(r => !Full(r) && r[8] == "1").Take(8)) sb.AppendLine($"- {string.Join(",", r)}");
        return sb.ToString();
    }
}
