using System;
using System.Collections.Generic;
using System.Diagnostics;
using Graphs;
using TDA.Ph;
using Xunit;
using Xunit.Abstractions;

namespace TDA.DimReduction.Tests;

/// <summary>
/// P0 profiler (measure-first): per-`Evaluate` stage breakdown for the SPRED objective across a few n,
/// so P1 optimization order follows the profile rather than a guess. Not an assertion test — a
/// benchmark that logs; re-run after each P1 change to measure the win. Reference barcode excluded (it
/// is built once per objective, not per eval); the projected pipeline is what runs `maxIters` times.
/// </summary>
public sealed class SpredProfileTests
{
    private readonly ITestOutputHelper _out;
    public SpredProfileTests(ITestOutputHelper output) => _out = output;

    private const int MaxDim = 2;
    private const int Reps = 3;

    private static double[][] Cylinder3D(int n, int seed)
    {
        var rng = new Random(seed);
        var pts = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double t = 2.0 * Math.PI * rng.NextDouble();
            double h = 4.0 * rng.NextDouble() - 2.0;
            pts[i] = new[] { Math.Cos(t), Math.Sin(t), h };
        }
        return pts;
    }

    private static GraphCompilerConfig Recipe() => new()
    {
        Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 10 },
        Filter = new FilterConfig { Kind = FilterKind.OrRule },
        Repair = new RepairConfig { Kind = RepairKind.NoRepair },
        Projection = new DistanceProjection(),
    };

    private static double[][] ProjectXY(double[][] data)
    {
        var y = new double[data.Length][];
        for (int i = 0; i < data.Length; i++) y[i] = new[] { data[i][0], data[i][1] };
        return y;
    }

    // Build a barcode, timing graph-construction / Rips / PH separately.
    private static (double graphMs, double ripsMs, double phMs, Barcode bc) BuildTimed(
        double[][] feats, GraphCompilerConfig recipe)
    {
        var sw = Stopwatch.StartNew();
        var metric = GraphMetric.FromFeatures(feats);
        var graph = GraphCompiler.Build(recipe, feats.Length, metric).Graph;
        double g = sw.Elapsed.TotalMilliseconds; sw.Restart();
        var filt = RipsFiltration.RipsFromGraph(graph, FiltrationWeights.RawDistance, MaxDim);
        double r = sw.Elapsed.TotalMilliseconds; sw.Restart();
        var bc = PersistentHomology.Compute(filt, MaxDim);
        double p = sw.Elapsed.TotalMilliseconds;
        return (g, r, p, bc);
    }

    [Fact]
    [Trait("Category", "Benchmark")]   // manual profiler; exclude from the default suite via --filter "Category!=Benchmark"
    public void Profile_PerEvaluateStageBreakdown()
    {
        var essential = DiagramMetrics.EssentialPolicy.FinitePenalty(1.0);

        // JIT warmup — the first Wasserstein/Hungarian call otherwise contaminates the smallest n.
        {
            var wd = Cylinder3D(40, 1);
            Barcode wref = BuildTimed(wd, Recipe()).bc;
            Barcode wb = BuildTimed(ProjectXY(wd), Recipe()).bc;
            DiagramMetrics.Wasserstein(wb, wref, 0, 2.0, essential);
            DiagramMetrics.Wasserstein(wb, wref, 1, 2.0, essential);
        }

        _out.WriteLine("n     #H0 #H1 |  graph   rips     ph    W(H0)   W(H1) |  total(ms)  dominant");
        foreach (int n in new[] { 60, 100, 150, 200 })
        {
            double[][] data = Cylinder3D(n, seed: 3);
            var recipe = Recipe();
            Barcode refBc = BuildTimed(data, recipe).bc;        // reference: once, excluded from per-eval

            double[][] proj = ProjectXY(data);
            double g = 0, r = 0, p = 0, w0 = 0, w1 = 0;
            Barcode projBc = refBc;
            for (int rep = 0; rep < Reps; rep++)
            {
                var b = BuildTimed(proj, recipe);
                projBc = b.bc;
                g += b.graphMs; r += b.ripsMs; p += b.phMs;
                var sw = Stopwatch.StartNew();
                DiagramMetrics.Wasserstein(b.bc, refBc, 0, 2.0, essential);
                w0 += sw.Elapsed.TotalMilliseconds; sw.Restart();
                DiagramMetrics.Wasserstein(b.bc, refBc, 1, 2.0, essential);
                w1 += sw.Elapsed.TotalMilliseconds;
            }
            g /= Reps; r /= Reps; p /= Reps; w0 /= Reps; w1 /= Reps;

            int nH0 = 0, nH1 = 0;
            foreach (Bar bar in projBc.Bars) { if (bar.Dimension == 0) nH0++; else if (bar.Dimension == 1) nH1++; }

            double total = g + r + p + w0 + w1;
            var stages = new (string nm, double v)[] { ("graph", g), ("rips", r), ("ph", p), ("W(H0)", w0), ("W(H1)", w1) };
            string dom = ""; double domv = -1;
            foreach (var (nm, v) in stages) if (v > domv) { domv = v; dom = nm; }
            _out.WriteLine($"{n,-5} {nH0,4} {nH1,3} | {g,6:F1} {r,6:F1} {p,6:F1} {w0,7:F1} {w1,7:F1} | {total,8:F1}  {dom} ({100 * domv / total:F0}%)");
        }
        Assert.True(true);
    }

    [Fact]
    [Trait("Category", "Benchmark")]
    public void Profile_PruningCollapsesH1Wasserstein()
    {
        var essential = DiagramMetrics.EssentialPolicy.FinitePenalty(1.0);
        const double Tau = 0.25;   // persistence threshold

        // warmup
        { var wd = Cylinder3D(40, 1); var a = BuildTimed(wd, Recipe()).bc; var b = BuildTimed(ProjectXY(wd), Recipe()).bc; DiagramMetrics.Wasserstein(b, a, 1, 2.0, essential); }

        _out.WriteLine($"prune τ={Tau}:  n     #H1  maxPersH1  #H1≥τ |  W(H1) ms   W(H1)@τ ms   speedup");
        foreach (int n in new[] { 60, 100, 150, 200 })
        {
            double[][] data = Cylinder3D(n, seed: 3);
            Barcode refBc = BuildTimed(data, Recipe()).bc;
            Barcode projBc = BuildTimed(ProjectXY(data), Recipe()).bc;

            int nH1 = 0; double maxPers = 0;
            foreach (Bar bar in projBc.Bars)
                if (bar.Dimension == 1 && !bar.IsInfinite) { nH1++; if (bar.Persistence > maxPers) maxPers = bar.Persistence; }

            Barcode refP = PruneLocal(refBc, Tau), projP = PruneLocal(projBc, Tau);
            int nH1p = 0;
            foreach (Bar bar in projP.Bars) if (bar.Dimension == 1) nH1p++;

            var sw = Stopwatch.StartNew();
            DiagramMetrics.Wasserstein(projBc, refBc, 1, 2.0, essential);
            double full = sw.Elapsed.TotalMilliseconds; sw.Restart();
            DiagramMetrics.Wasserstein(projP, refP, 1, 2.0, essential);
            double pruned = sw.Elapsed.TotalMilliseconds;

            _out.WriteLine($"           {n,-5} {nH1,4}   {maxPers,7:F3}   {nH1p,4}  | {full,8:F1}   {pruned,9:F2}    {full / Math.Max(pruned, 0.001),6:F0}x");
        }
        Assert.True(true);
    }

    private static Barcode PruneLocal(Barcode bc, double tau)
    {
        var kept = new List<Bar>(bc.Bars.Count);
        foreach (Bar b in bc.Bars) if (b.Persistence >= tau) kept.Add(b);
        return new Barcode(kept, bc.AxisLabel);
    }
}
