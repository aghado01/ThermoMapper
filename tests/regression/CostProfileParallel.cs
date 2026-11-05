using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using Maths.LinAlg;
using Maths.Regression.Bgp;
using Maths.Regression.Changepoint;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;

namespace Maths.Regression.Tests;

/// <summary>
/// Chains in parallel (the brief's <c>alloc:parallel</c>, <c>membound</c> and <c>settings</c> blocks): independent instances of
/// one workload on 1, 2, 4 and 8 threads, with the wall time, the bytes allocated, the GC pause and the collections of each
/// run, so that the speedup, the efficiency and the pause share are read from one set of runs. Parallel readings on a machine
/// that other work shares measure that machine too; the block's tier says how.
/// </summary>
internal static class CostParallel
{
    private static double Sink { get => CostKernels.Sink; set => CostKernels.Sink = value; }
    private static object? ObjectSink { set => CostKernels.ObjectSink = value; }

    internal sealed record Workload(string Name, Func<int, Action<int>> Make);

    /// <summary>Runs the workload on each thread count and writes the rows. Each thread has its own instance (<c>Make</c> is given
    /// 1-based worker indices); the iteration count is fixed from the single-thread rate so that one thread runs about
    /// <paramref name="seconds"/> (scaled by the trial).</summary>
    internal static void Run(string block, Workload w, int[] threads, double seconds = 0.4, long? memoryPerThreadBytes = null)
    {
        seconds *= CostProfile.Scale;
        Action<int> probe = w.Make(0);
        probe(1);
        int r = 2;
        double per;
        while (true)
        {
            long t0 = Stopwatch.GetTimestamp();
            probe(r);
            per = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency / r;
            if (per * r >= 0.05 * CostProfile.Scale || r >= 1 << 24) break;
            r *= 4;
        }
        int iterations = (int)Math.Clamp(seconds / per, 1, 1 << 24);
        double wall1 = double.NaN;
        foreach (int c in threads)
        {
            if (memoryPerThreadBytes is long m && m * c > CostProfile.FreeMemoryGb() * 1e9 / 3)
            {
                CostProfile.Row(block, w.Name, $"threads={c}", "skipped", "footprint above a third of the free memory");
                continue;
            }
            Action<int>[] workers = [.. Enumerable.Range(1, c).Select(w.Make)];
            foreach (Action<int> a in workers) a(Math.Max(1, iterations / 20));
            var walls = new List<double>();
            var bytes = new List<double>();
            var pauses = new List<double>();
            var g0 = new List<double>(); var g1 = new List<double>(); var g2 = new List<double>();
            for (int rep = 0; rep < 3; rep++)
            {
                GC.Collect();
                GC.WaitForPendingFinalizers();
                var barrier = new Barrier(c + 1);
                var threadsArray = workers.Select(a => new Thread(() => { barrier.SignalAndWait(); a(iterations); barrier.SignalAndWait(); }, 16 * 1024 * 1024)).ToArray();
                foreach (Thread t in threadsArray) t.Start();
                long b0 = GC.GetTotalAllocatedBytes(true);
                TimeSpan p0 = GC.GetTotalPauseDuration();
                int c0 = GC.CollectionCount(0), c1 = GC.CollectionCount(1), c2 = GC.CollectionCount(2);
                barrier.SignalAndWait();
                long t0 = Stopwatch.GetTimestamp();
                barrier.SignalAndWait();
                double wall = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
                foreach (Thread t in threadsArray) t.Join();
                walls.Add(wall);
                bytes.Add(GC.GetTotalAllocatedBytes(true) - b0);
                pauses.Add((GC.GetTotalPauseDuration() - p0).TotalSeconds);
                g0.Add(GC.CollectionCount(0) - c0); g1.Add(GC.CollectionCount(1) - c1); g2.Add(GC.CollectionCount(2) - c2);
            }
            double wallMedian = CostProfile.Median(walls);
            if (c == threads[0]) wall1 = wallMedian;
            string p = $"threads={c}";
            CostProfile.Row(block, w.Name, p, "iterations_per_thread", iterations);
            CostProfile.Row(block, w.Name, p, "wall_s", wallMedian);
            CostProfile.Row(block, w.Name, p, "wall_s_q1", CostProfile.Quantile(walls, 0.25));
            CostProfile.Row(block, w.Name, p, "wall_s_q3", CostProfile.Quantile(walls, 0.75));
            CostProfile.Row(block, w.Name, p, "ns_per_iteration_per_thread", wallMedian * 1e9 / iterations);
            CostProfile.Row(block, w.Name, p, "speedup", c * wall1 / wallMedian);
            CostProfile.Row(block, w.Name, p, "efficiency", wall1 / wallMedian);
            CostProfile.Row(block, w.Name, p, "bytes_per_iteration", CostProfile.Median(bytes) / ((double)c * iterations));
            CostProfile.Row(block, w.Name, p, "gc_pause_s", CostProfile.Median(pauses));
            CostProfile.Row(block, w.Name, p, "gc_pause_share_of_wall", CostProfile.Median(pauses) / wallMedian);
            CostProfile.Row(block, w.Name, p, "gen0_collections", CostProfile.Median(g0));
            CostProfile.Row(block, w.Name, p, "gen1_collections", CostProfile.Median(g1));
            CostProfile.Row(block, w.Name, p, "gen2_collections", CostProfile.Median(g2));
        }
    }

    // ---- The workloads ----

    internal static Workload BarsStep(int n)
    {
        var fx = new CostEngines.BarsFixture(n);
        return new($"BARS step n={n}", i =>
        {
            ReversibleJumpChain<KnotConfig> chain = fx.Chain(100 + i);
            for (int s = 0; s < 500; s++) chain.Step();
            return c => { for (int s = 0; s < c; s++) chain.Step(); };
        });
    }

    internal static Workload BapsGibbs(int n, int nu)
    {
        var rng = new Xoshiro256PlusPlus(91 + n);
        double[] x = Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();
        double[] y = x.Select(t => Math.Sin(2 * Math.PI * t) + 0.15 * (rng.NextDouble() - 0.5)).ToArray();
        double[,] z = new SplineBasis(3).Design(CostEngines.BarsFixture.Equispaced(nu - 4), x);
        var baps = new Baps(z, y, new DifferencePenalty(2));
        return new($"BAPS Gibbs n={n} nu={nu}", i => c => ObjectSink = baps.Run(BapsScheme.Gibbs, chains: 1, masterSeed: 100 + i, burn: 10, samples: Math.Max(2, c)));
    }

    internal static Workload BgpFit(int n)
    {
        var rng = new Xoshiro256PlusPlus(7 + n);
        var x = new double[n, 3];
        var y = new double[n];
        for (int i = 0; i < n; i++) { for (int c = 0; c < 3; c++) x[i, c] = rng.NextDouble(); y[i] = Math.Sin(6 * x[i, 0]) + 0.1 * rng.NextDouble(); }
        var gp = new GpRegression(x, y, 0.01, new SquaredExponentialKernel());
        return new($"BGP fit n={n}", i => c => { for (int s = 0; s < c; s++) Sink = gp.Fit(0.05 + 1e-4 * s + 1e-3 * i).LogMarginal; });
    }

    internal static Workload TemperingCheap(int levels) => new($"PT cheap {levels} rungs", i =>
    {
        var ladder = new ParallelTempering<int>(_ => new IRjMove<int>[] { new CountBirth(), new CountDeath(), new CountWalk() }, new CountTarget(), 4,
                                                ParallelTempering<int>.GeometricLadder(levels, 0.01), TemperingPath.Likelihood, 100 + i);
        return c => { for (int s = 0; s < c; s++) ladder.Step(); };
    });

    internal static Workload SmcAdvance(int particles) => new($"SMC advance N={particles} d=10", i =>
    {
        var target = new GaussTarget(0);
        var options = new SmcOptions { Particles = particles, MutationSteps = 1, ResampleBelow = 0.5, MaxStages = 1_000_000, MaxDegreeOfParallelism = 1 };
        SmcRun<double[]> run = SmcRun<double[]>.Start(new SmcStage<double[]>(target, new Temperature(0.05, TemperingPath.Likelihood)), new GaussParticles(10),
                                                      new RandomWalkKernel(), options, 100 + i);
        for (int s = 1; s <= 10; s++) run.Advance(new SmcStage<double[]>(target, new Temperature(0.05 + 0.95 * s / 10.0, TemperingPath.Likelihood)));
        var stage = new SmcStage<double[]>(target, new Temperature(1.0, TemperingPath.Likelihood));
        return c => { for (int s = 0; s < c; s++) run.Advance(stage); };
    });

    internal static Workload ChangepointInfer(int n, int k)
    {
        var rng = new Xoshiro256PlusPlus(3 + n);
        double[] levels = [0.0, 2.0, -1.0, 1.0];
        double[] y = Enumerable.Range(0, n).Select(i => levels[i * 4 / n] + 0.5 * (rng.NextDouble() - 0.5)).ToArray();
        (double sigma2, double nu, double rho2) = ExactChangepoint.DefaultHyperparameters(y);
        var model = new GaussianSegmentModel(y, sigma2, nu, rho2);
        return new($"Changepoint n={n} k={k}", _ => c => { for (int s = 0; s < c; s++) ObjectSink = ExactChangepoint.Infer(y, model, k); });
    }

    internal static Workload CountChain() => new("RJ chain, cheap count target", i =>
    {
        var chain = new ReversibleJumpChain<int>(new IRjMove<int>[] { new CountBirth(), new CountDeath(), new CountWalk() }, new CountTarget(), 4, new Xoshiro256PlusPlus(100 + i));
        return c => { for (int s = 0; s < c; s++) chain.Step(); };
    });

    internal static void AllocationBlock()
    {
        int[] threads = [1, 2, 4, 8];
        foreach (Workload w in new[] { BarsStep(500), BapsGibbs(1000, 104), BgpFit(200), BgpFit(800), TemperingCheap(10), SmcAdvance(4000), ChangepointInfer(300, 10), CountChain() })
            Run("alloc:parallel", w, threads);
    }

    // ---- Memory-bound regimes ----

    internal static void MemoryBoundBlock()
    {
        const string block = "membound";
        int[] threads = [1, 2, 4, 8];

        // BAPS's residual loop over the design: the dense n × nu scan and the banded apply, from inside the L2 to far past the L3.
        foreach (double megabytes in new[] { 0.1, 1.0, 8.0, 64.0, 512.0 })
        {
            const int nu = 100;
            int n = Math.Max(8, (int)(megabytes * 1e6 / (8.0 * nu)));
            if (CostProfile.FreeMemoryGb() * 1e9 < 3.0 * megabytes * 1e6) { CostProfile.Row(block, $"baps residual {megabytes} MB", "", "skipped", "free memory"); continue; }
            var rng = new Xoshiro256PlusPlus(5);
            var z = new double[n, nu];
            var y = new double[n];
            for (int i = 0; i < n; i++)
            {
                int lo = (int)((long)i * (nu - 4) / n);
                for (int j = lo; j < lo + 4; j++) z[i, j] = 0.25 + 0.01 * rng.NextDouble();
                y[i] = rng.NextDouble();
            }
            var beta = new double[nu];
            for (int j = 0; j < nu; j++) beta[j] = 0.1 * j;
            var bd = new BandedDesign(z);
            var fit = new double[n];
            string name = $"baps residual {megabytes} MB (n={n}, nu={nu})";
            CostProfile.Point(block, name, new List<CostProfile.Part>
            {
                new("residual_dense", c => { for (int i = 0; i < c; i++) { double ss = 0; for (int r = 0; r < n; r++) { double f = 0; for (int j = 0; j < nu; j++) f += z[r, j] * beta[j]; double d = y[r] - f; ss += d * d; } Sink = ss; } }, (double)n * nu, "madd"),
                new("residual_banded_apply", c => { for (int i = 0; i < c; i++) { bd.Apply(beta, fit); double ss = 0; for (int r = 0; r < n; r++) { double d = y[r] - fit[r]; ss += d * d; } Sink = ss; } }, (double)n * (bd.Bandwidth + 1), "madd"),
            });
            Run(block, new Workload(name + " dense, shared design", _ => c => { for (int i = 0; i < c; i++) { double ss = 0; for (int r = 0; r < n; r++) { double f = 0; for (int j = 0; j < nu; j++) f += z[r, j] * beta[j]; double d = y[r] - f; ss += d * d; } Sink = ss; } }), threads);
        }

        // BARS's design build and support gate: n × nu doubles written and scanned, per evaluation.
        foreach (int n in new[] { 1000, 10_000, 100_000, 500_000 })
        {
            var fx = new CostEngines.BarsFixture(n);
            KnotConfig cfg = CostEngines.BarsFixture.Equispaced(8);
            double footprint = 8.0 * n * 12;
            if (CostProfile.FreeMemoryGb() * 1e9 < 4.0 * footprint) { CostProfile.Row(block, $"bars evaluate n={n}", "", "skipped", "free memory"); continue; }
            CostProfile.Point(block, $"bars evaluate n={n} (design {footprint / 1e6:F1} MB)", new List<CostProfile.Part>
            {
                new("design_build", c => { for (int i = 0; i < c; i++) ObjectSink = fx.Basis.Design(cfg, fx.X); }, 4.0 * n, "basis value"),
                new("evaluate_whole", c => { for (int i = 0; i < c; i++) ObjectSink = fx.Target.Evaluate(cfg); }, 1, "evaluation"),
            });
            Run(block, new Workload($"bars evaluate n={n}", _ => c => { for (int i = 0; i < c; i++) ObjectSink = fx.Target.Evaluate(cfg); }), threads, memoryPerThreadBytes: (long)(3 * footprint));
        }

        // BGP's kernel matrix: n × n, one thread's matrix per chain.
        foreach (int n in new[] { 200, 400, 800, 1600, 3200 })
        {
            double footprint = 8.0 * n * n;
            if (CostProfile.FreeMemoryGb() * 1e9 < 4.0 * footprint) { CostProfile.Row(block, $"bgp kernel fill n={n}", "", "skipped", "free memory"); continue; }
            var kernel = new SquaredExponentialKernel();
            var rng = new Xoshiro256PlusPlus(n);
            var sq = new double[n, n];
            for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) sq[i, j] = rng.NextDouble();
            var a = new double[n, n];
            void Fill(double[,] target, double t)
            {
                for (int i = 0; i < n; i++)
                {
                    target[i, i] = 1.01;
                    for (int j = 0; j < i; j++) { double v = kernel.Evaluate(sq[i, j], t); target[i, j] = v; target[j, i] = v; }
                }
            }
            CostProfile.Point(block, $"bgp kernel fill n={n} ({footprint / 1e6:F1} MB per matrix)", new List<CostProfile.Part>
            {
                new("kernel_fill", c => { for (int i = 0; i < c; i++) Fill(a, 0.05); }, (double)n * (n + 1) / 2, "pair"),
            });
            Run(block, new Workload($"bgp kernel fill n={n}", _ => { var mine = new double[n, n]; return c => { for (int i = 0; i < c; i++) Fill(mine, 0.05); }; }), threads, memoryPerThreadBytes: (long)(2 * footprint));
        }

        // The changepoint's tables: four (n+1)² arrays.
        foreach (int n in new[] { 300, 1000, 2000, 3000 })
        {
            double footprint = 8.0 * 5 * (n + 1.0) * (n + 1.0);
            if (CostProfile.FreeMemoryGb() * 1e9 < 3.0 * footprint) { CostProfile.Row(block, $"changepoint n={n}", "", "skipped", "free memory"); continue; }
            Workload w = ChangepointInfer(n, 3);
            CostProfile.Row(block, w.Name, "", "table_bytes", footprint);
            Run(block, w, threads, seconds: 1.0, memoryPerThreadBytes: (long)footprint);
        }
    }

    // ---- Runtime settings ----

    internal static void SettingsBlock()
    {
        const string block = "settings";

        // The warm-up: BARS's step in 25 ms batches from the first call, no earlier warm-up in this process.
        {
            var fx = new CostEngines.BarsFixture(500);
            ReversibleJumpChain<KnotConfig> chain = fx.Chain(3);
            for (int k = 0; k < 24; k++)
            {
                long t0 = Stopwatch.GetTimestamp();
                int steps = 0;
                while ((Stopwatch.GetTimestamp() - t0) * 1000.0 / Stopwatch.Frequency < 25.0) { chain.Step(); steps++; }
                CostProfile.Row(block, "warm-up BARS n=500", $"batch{k:00}", "ns_per_step", (Stopwatch.GetTimestamp() - t0) * 1e9 / Stopwatch.Frequency / steps);
            }
        }

        {
            var fx = new CostEngines.BarsFixture(500);
            ReversibleJumpChain<KnotConfig> chain = fx.Chain(4);
            for (int i = 0; i < 1000; i++) chain.Step();
            var rngGp = new Xoshiro256PlusPlus(7);
            var x = new double[400, 3];
            var y = new double[400];
            for (int i = 0; i < 400; i++) { for (int c = 0; c < 3; c++) x[i, c] = rngGp.NextDouble(); y[i] = Math.Sin(6 * x[i, 0]) + 0.1 * rngGp.NextDouble(); }
            var gp = new GpRegression(x, y, 0.01, new SquaredExponentialKernel());
            var rngB = new Xoshiro256PlusPlus(91);
            double[] xs = Enumerable.Range(0, 1000).Select(i => (i + 0.5) / 1000).ToArray();
            double[] ys = xs.Select(t => Math.Sin(2 * Math.PI * t) + 0.15 * (rngB.NextDouble() - 0.5)).ToArray();
            double[,] z = new SplineBasis(3).Design(CostEngines.BarsFixture.Equispaced(100), xs);
            var baps = new Baps(z, ys, new DifferencePenalty(2));
            CostProfile.Point(block, "one thread", new List<CostProfile.Part>
            {
                new("bars_step_n500", c => { for (int i = 0; i < c; i++) chain.Step(); }, 1, "step"),
                new("bgp_fit_n400", c => { for (int i = 0; i < c; i++) CostKernels.Sink = gp.Fit(0.05 + 1e-5 * i).LogMarginal; }, 1, "fit"),
                new("baps_gibbs_nu104", c => { for (int i = 0; i < c; i++) CostKernels.ObjectSink = baps.Run(BapsScheme.Gibbs, chains: 1, masterSeed: 3, burn: 10, samples: 90); }, 100, "iteration"),
            });
        }
        foreach (Workload w in new[] { BarsStep(500), BgpFit(400), BapsGibbs(1000, 104) })
            Run(block, w, [1, 4]);
    }
}
