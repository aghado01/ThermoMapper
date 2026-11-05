using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Runtime.CompilerServices;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The efficiency pass's lever 1, step 1 (<c>efficiency-pass-brief.md</c>): the break-even of a prefix table for a main-effect
/// hinge, measured before any production change. A table T(z, v, s)[k] = Σ_{s(xᵢ − t_k) > 0} zᵢ (s(xᵢ − t_k))^q for a held
/// column z and a variable v costs one pass over the cases in sorted order to build, against one vector dot per proposal it
/// answers. Whether it pays depends on R, the build's cost over a dot's, and on H, how many one-way proposals reach one
/// (column, variable) before the column leaves.
/// <list type="bullet">
/// <item><see cref="Build_AgainstDot"/> reads R at n = 1,000, 10,000 and 50,000, the columns in the cache and out of it, and the
/// scan's error against the explicit dot.</item>
/// <item><see cref="Reach_PerColumnAndVariable"/> replays a burned chain at the cost grid's points and counts, for a gate that
/// builds a (column, variable)'s table at its (r + 1)th hit, the dots, builds and lookups a table path would make, over all the
/// one-way proposals that reached the design's border.</item>
/// </list>
/// Behind <c>BASS_QUALIFICATION</c> = 1; rows go to <c>BASS_COST_OUT</c>. Times wait for a quiet machine; counts do not.
/// </summary>
public sealed class BassPrefixBreakEvenTests
{
    private readonly ITestOutputHelper _out;
    public BassPrefixBreakEvenTests(ITestOutputHelper output) => _out = output;

    // ---- The scan ----

    /// <summary>The table of a held column over the cases of one variable sorted by x (the prototype of lever 1's production
    /// scan): <paramref name="plus"/>[k] = Σ_{x &gt; t_k} g (x − t_k)^q and <paramref name="minus"/>[k] = Σ_{x &lt; t_k} g (t_k − x)^q,
    /// each sum about its own knot from the sum about its neighbour (A₀ ← A₀ + G, A₂ ← A₂ + δ(2A₁ + δA₀), A₁ ← A₁ + δA₀, G the
    /// neighbour group's total), so every update is the size of the quantity itself.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal static void Scan(double[] g, int[] order, int[] start, double[] knots, int degree, double[] plus, double[] minus)
    {
        int count = knots.Length;
        double a0 = 0.0, a1 = 0.0, a2 = 0.0;
        plus[count - 1] = 0.0;
        for (int k = count - 2; k >= 0; k--)
        {
            double delta = knots[k + 1] - knots[k], group = 0.0;
            for (int p = start[k + 1]; p < start[k + 2]; p++) group += g[order[p]];
            a0 += group;
            if (degree == 1)
            {
                a1 += delta * a0;
                plus[k] = a1;
            }
            else
            {
                a2 += delta * (2.0 * a1 + delta * a0);
                a1 += delta * a0;
                plus[k] = a2;
            }
        }
        double b0 = 0.0, b1 = 0.0, b2 = 0.0;
        minus[0] = 0.0;
        for (int k = 1; k < count; k++)
        {
            double delta = knots[k] - knots[k - 1], group = 0.0;
            for (int p = start[k - 1]; p < start[k]; p++) group += g[order[p]];
            b0 += group;
            if (degree == 1)
            {
                b1 += delta * b0;
                minus[k] = b1;
            }
            else
            {
                b2 += delta * (2.0 * b1 + delta * b0);
                b1 += delta * b0;
                minus[k] = b2;
            }
        }
    }

    // A variable's cases sorted by x, its distinct knots and each knot's first sorted position (the last entry is n).
    private static (int[] Order, int[] Start, double[] Knots) Sorted(double[] x)
    {
        int n = x.Length;
        int[] order = [.. Enumerable.Range(0, n)];
        Array.Sort(order, (a, b) => x[a].CompareTo(x[b]));
        var start = new List<int>();
        var knots = new List<double>();
        for (int p = 0; p < n; p++)
            if (p == 0 || x[order[p]] != x[order[p - 1]])
            {
                start.Add(p);
                knots.Add(x[order[p]]);
            }
        start.Add(n);
        return (order, [.. start], [.. knots]);
    }

    // ---- R: the build against a dot ----

    [Fact]
    public void Build_AgainstDot()
    {
        if (!BassQualification.Enabled) return;
        var rows = new BassCost.Rows("breakeven", "build");
        foreach ((string q, string v) in BassCost.Machine()) rows.Add("machine", q, v);
        var load = new BassCost.Load();
        var rng = new Maths.Rng.Xoshiro256PlusPlus(seed: 331);
        foreach (int n in new[] { 1000, 10000, 50000 })
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            (int[] order, int[] start, double[] knots) = Sorted(x);
            int k = knots.Length;
            var plus = new double[k];
            var minus = new double[k];
            long cacheBytes = 4L << 20, farBytes = 256L << 20;
            int near = (int)Math.Clamp(cacheBytes / (8L * n) / 4, 2, 8), far = (int)Math.Clamp(farBytes / (8L * n), 8, 600);
            double[][] columns = [.. Enumerable.Range(0, far).Select(_ => Centred(rng, n))];
            double[] u = [.. x.Select(v => Math.Max(v - 0.4, 0.0))];

            foreach (int degree in new[] { 1, 2 })
            {
                // The error against the explicit dot at a few knots, relative to ‖z‖‖u‖.
                double worst = 0.0;
                foreach (int knot in new[] { 0, 1, k / 7, k / 2, k - k / 7, k - 2, k - 1 })
                {
                    Scan(columns[0], order, start, knots, degree, plus, minus);
                    double t = knots[knot], dot = 0.0, norm = 0.0, unorm = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        double h = Math.Pow(Math.Max(x[i] - t, 0.0), degree);
                        dot += columns[0][i] * h;
                        norm += columns[0][i] * columns[0][i];
                        unorm += h * h;
                    }
                    worst = Math.Max(worst, Math.Abs(plus[knot] - dot) / Math.Sqrt(norm * unorm + 1e-300));
                }
                string key = $"n={n} q={degree}";
                rows.Add(key, "scan error vs dot, relative to ‖z‖‖u‖ (plus sign)", worst);

                // Interleaved: the dot and the build, the columns near (a few, in the cache) and far (out of it).
                int i0 = 0, i1 = 0, i2 = 0, i3 = 0, repeats = Math.Max(1, 2_000_000 / n);
                Dictionary<string, double[]> times = BassCost.Interleave(
                [
                    ("dot near", Timer(repeats, () => PairwiseSum.Dot(columns[i0++ % near], u))),
                    ("dot far", Timer(repeats, () => PairwiseSum.Dot(columns[i1++ % far], u))),
                    ("build near", Timer(repeats, () => { Scan(columns[i2++ % near], order, start, knots, degree, plus, minus); return plus[k / 2] + minus[k / 2]; })),
                    ("build far", Timer(repeats, () => { Scan(columns[i3++ % far], order, start, knots, degree, plus, minus); return plus[k / 2] + minus[k / 2]; })),
                ]);
                double dotNear = BassCost.Median(times["dot near"]), dotFar = BassCost.Median(times["dot far"]);
                double buildNear = BassCost.Median(times["build near"]), buildFar = BassCost.Median(times["build far"]);
                rows.Add(key, "columns near, far", $"{near}; {far}");
                rows.Add(key, "dot near ns per case", dotNear * 1000 / n);
                rows.Add(key, "dot far ns per case", dotFar * 1000 / n);
                rows.Add(key, "build (both signs) near ns per knot", buildNear * 1000 / k);
                rows.Add(key, "build (both signs) far ns per knot", buildFar * 1000 / k);
                rows.Add(key, "R near (build over dot)", buildNear / dotNear);
                rows.Add(key, "R far (build over dot)", buildFar / dotFar);
                _out.WriteLine($"{key}: dot {dotNear * 1000 / n:F2}/{dotFar * 1000 / n:F2} ns/case (near/far), build {buildNear * 1000 / k:F2}/{buildFar * 1000 / k:F2} ns/knot, " +
                               $"R {buildNear / dotNear:F1}/{buildFar / dotFar:F1}, error {worst:E1}");
            }
        }
        (string[] end, double cpu, double seconds, double others) = load.End();
        rows.Add("load", "other test hosts at the start", string.Join(" ", load.AtStart));
        rows.Add("load", "other test hosts at the end", string.Join(" ", end));
        rows.Add("load", "own CPU over wall", cpu);
        rows.Add("load", "other processes' cores busy", others);
        rows.Add("load", "block seconds", seconds);
        rows.Add("load", "flagged", load.AtStart.Length > 0 || end.Length > 0 ? "yes" : "no");
        _out.WriteLine($"rows: {rows.Flush()}");
    }

    private static double[] Centred(Maths.Rng.Xoshiro256PlusPlus rng, int n)
    {
        double[] z = new double[n];
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += z[i] = BassQualification.Gaussian(rng);
        double mean = sum / n;
        for (int i = 0; i < n; i++) z[i] -= mean;
        return z;
    }

    private static Func<int, double> Timer(int repeats, Func<double> call)
    {
        double sink = 0.0;
        return [MethodImpl(MethodImplOptions.AggressiveOptimization)] (int _) =>
        {
            long t0 = Stopwatch.GetTimestamp();
            for (int r = 0; r < repeats; r++) sink += call();
            double us = (Stopwatch.GetTimestamp() - t0) * 1e6 / Stopwatch.Frequency / repeats;
            return us + sink * 0.0;
        };
    }

    // ---- H: the reach of one-way proposals ----

    private static readonly int[] Gates = [0, 1, 2, 4, 8, 16, 32, 64, 128, 1 << 20];

    [Fact]
    public void Reach_PerColumnAndVariable()
    {
        if (!BassQualification.Enabled) return;
        var rows = new BassCost.Rows("breakeven", "reach");
        foreach ((string q, string v) in BassCost.Machine()) rows.Add("machine", q, v);
        (BassCost.Point Point, int Sweeps)[] points =
        [
            (BassCost.Base, 20000),
            (BassCost.Base with { N = 200 }, 20000),
            (BassCost.Base with { N = 50000 }, 4000),
            (BassCost.Base with { N = 10000, Cap = 5 }, 20000),
            (BassCost.Base with { N = 10000, Cap = 200 }, 6000),
            (BassCost.Base with { P = 200 }, 20000),
            (BassCost.Base with { Order = 5 }, 20000),
            (new BassCost.Point(BassCost.Surface.Categorical, 4096, 7, BassCost.BaseCap, BassCost.BaseOrder, L: 1024), 10000),
            (BassCost.Base with { Surface = BassCost.Surface.S10, P = 5 }, 20000),
        ];
        foreach ((BassCost.Point x, int sweeps) in points)
        {
            (BassModel model, BassChain chain, int burn) = BassCost.Burned(x);
            var counting = new Counting(model.Target);
            BassChain run = BassChain.Restore(model, chain.Snapshot(), counting);
            var counters = Gates.Select(_ => new Dictionary<BassTerm, int[]>()).ToArray();
            var dots = new long[Gates.Length];
            var builds = new long[Gates.Length];
            var lookups = new long[Gates.Length];
            var liveMax = new long[Gates.Length];
            long proposals = 0, oneWay = 0, pairHits = 0, borders = 0;
            int p = model.Basis.Inputs.Cases.ContinuousCount, knotsMax = 0;
            for (int v = 0; v < p; v++) knotsMax = Math.Max(knotsMax, model.Basis.Inputs.ContinuousKnots(v).Length);
            counting.Seen = cand =>
            {
                BassState current = run.Current;
                if (ReferenceEquals(cand.TermArray, current.TermArray)) return;   // τ's evaluations
                proposals++;
                if (cand.Count < current.Count) return;                              // a death
                if (cand.Fit.Refusal is BassRefusal.Support or BassRefusal.Cap) return;   // refused before any border
                BassTerm added = cand.TermArray[^1];
                borders++;
                if (added.Continuous.Count != 1 || added.Categorical.Count != 0 || added.Functional.Count != 0) return;
                oneWay++;
                int variable = added.Continuous[0].Variable;
                pairHits += cand.Count - 1;
                for (int g = 0; g < Gates.Length; g++)
                    for (int j = 0; j < cand.Count - 1; j++)
                    {
                        BassTerm held = cand.TermArray[j];
                        if (!counters[g].TryGetValue(held, out int[]? c)) counters[g][held] = c = new int[p];
                        int hits = ++c[variable];
                        if (hits <= Gates[g]) dots[g]++;
                        else if (hits == Gates[g] + 1) builds[g]++;
                        else lookups[g]++;
                    }
            };
            double sinceLive = 0;
            for (int s = 0; s < sweeps; s++)
            {
                run.Step();
                // A column that left takes its counters with it.
                if (s % 16 == 0)
                {
                    var present = new HashSet<BassTerm>(run.Current.TermArray);
                    for (int g = 0; g < Gates.Length; g++)
                    {
                        foreach (BassTerm gone in counters[g].Keys.Where(t => !present.Contains(t)).ToArray()) counters[g].Remove(gone);
                        long live = counters[g].Values.Sum(c => c.Count(h => h > Gates[g]));
                        liveMax[g] = Math.Max(liveMax[g], live);
                    }
                    sinceLive++;
                }
            }
            string key = x.Key;
            rows.Add(key, "sweeps", sweeps);
            rows.Add(key, "burn-in sweeps", burn);
            rows.Add(key, "M at the end", run.Current.Count);
            rows.Add(key, "knots (largest variable)", knotsMax);
            rows.Add(key, "p", p);
            rows.Add(key, "proposals evaluated per sweep", (double)proposals / sweeps);
            rows.Add(key, "proposals with a border per sweep", (double)borders / sweeps);
            rows.Add(key, "one-way proposals per sweep", (double)oneWay / sweeps);
            rows.Add(key, "one-way share of the borders", borders > 0 ? (double)oneWay / borders : double.NaN);
            rows.Add(key, "held-column dots in one-way borders per sweep", (double)pairHits / sweeps);
            for (int g = 0; g < Gates.Length; g++)
            {
                string gate = Gates[g] == 1 << 20 ? "never" : Gates[g].ToString(CultureInfo.InvariantCulture);
                rows.Add(key, $"gate {gate}: dots", dots[g]);
                rows.Add(key, $"gate {gate}: builds", builds[g]);
                rows.Add(key, $"gate {gate}: lookups", lookups[g]);
                rows.Add(key, $"gate {gate}: most tables live at once", liveMax[g]);
            }
            _out.WriteLine($"{key}: M {run.Current.Count}, {(double)oneWay / sweeps:F3} one-way per sweep, {(double)pairHits / sweeps:F1} column dots per sweep; " +
                           string.Join(" ", Gates.Take(Gates.Length - 1).Select((r, g) => $"[{r}: {builds[g]} builds, {lookups[g]} lookups]")));
            Assert.True(sinceLive > 0);
            rows.Flush();
        }
    }

    private sealed class Counting(BassTarget target) : ISamplingTarget<BassState>
    {
        public Action<BassState>? Seen { get; set; }

        public ITargetEvaluation Evaluate(BassState state)
        {
            Seen?.Invoke(state);
            return target.Evaluate(state);
        }
    }
}
