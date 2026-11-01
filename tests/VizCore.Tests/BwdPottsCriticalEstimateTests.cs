using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Models.Potts;
using Graphs.Primitives;
using Maths.Rng;
using UserRepl;
using Xunit;
using Xunit.Abstractions;

namespace VizCore.Tests;

/// <summary>
/// The q-only temperature anchor (Blatt, Wiseman and Domany 1997, eq. 4.2) and its bracket. The physical references are
/// measured: Swendsen–Wang's own S₂, which knows nothing of the formula. Tolerances were fixed before the first run.
/// </summary>
public sealed class BwdPottsCriticalEstimateTests(ITestOutputHelper output)
{
    private const int Seed = 20261101;

    // ── The formula ─────────────────────────────────────────────────────────

    [Theory]
    [InlineData(2, 0.17204)]
    [InlineData(10, 0.10633)]
    [InlineData(20, 0.08922)]
    [InlineData(100, 0.06323)]
    public void TpsUpperBound_IsBwdEquation42(int q, double expected)
        => Assert.InRange(Math.Abs(BwdPottsCriticalEstimate.TpsUpperBound(q) - expected), 0.0, 5e-5);

    [Fact]
    public void TpsUpperBound_FallsWithQ()
    {
        // Connectivity falls as q rises at fixed p (Grimmett 2006, (3.22)), so no transition can move hotter with q.
        for (int q = 2; q < 1000; q++)
            Assert.True(BwdPottsCriticalEstimate.TpsUpperBound(q + 1) < BwdPottsCriticalEstimate.TpsUpperBound(q));
    }

    [Fact]
    public void Bracket_IsTheDeclaredMultiplesOfTheAnchor()
    {
        (double lo, double hi) = BwdPottsCriticalEstimate.Bracket(20);
        double tps = BwdPottsCriticalEstimate.TpsUpperBound(20);
        Assert.Equal(BwdPottsCriticalEstimate.BracketColdFactor * tps, lo);
        Assert.Equal(BwdPottsCriticalEstimate.BracketHotFactor * tps, hi);
    }

    // ── Against Swendsen–Wang ───────────────────────────────────────────────

    /// <summary>
    /// On a homogeneous cloud the giant's melt, read from where S₂ falls through half its cold value, sits within a factor
    /// of 2 of the anchor at each q, and moves colder as q rises. The old form, with the logarithm on top, misses by 3.3×
    /// at q = 20 and 6.6× at q = 100, and rises with q.
    /// </summary>
    [Fact]
    public void Melt_TracksTheAnchor_AndFallsWithQ()
    {
        var rng = new Xoshiro256PlusPlus(Seed);
        double[][] cloud = Enumerable.Range(0, 400).Select(_ => new[] { rng.NextDouble(), rng.NextDouble() }).ToArray();
        CsrGraph graph = SpcGraphBuilder.BuildResult(cloud, BwdParityTests.BuildGraphConfig(k: 10, bandwidth: 0)).Graph;

        var melts = new List<double>();
        foreach (int q in new[] { 5, 20, 100 })
        {
            double tps = BwdPottsCriticalEstimate.TpsUpperBound(q);
            double[] grid = SpcScheduleHelpers.LogSpaceGrid(tps / 4.0, tps * 4.0, 13);
            var s2 = new double[grid.Length];
            Parallel.For(0, grid.Length, i => s2[i] = S2(graph, grid[i], q, new RunBudget(2000, 500)));

            double melt = HalfPoint(grid, s2);
            output.WriteLine($"q={q}: T_ps={tps:F4} melt={melt:F4} melt/T_ps={melt / tps:F2} cold S2={s2[0]:F3}");
            Assert.InRange(melt / tps, 0.5, 2.0);
            melts.Add(melt);
        }
        Assert.True(melts[0] > melts[1] && melts[1] > melts[2], $"Melt should fall with q: {string.Join(", ", melts)}");
    }

    /// <summary>
    /// Past the bracket's hot end the parity data holds no cluster mass: S₂ is within a factor of 2 of its floor 1/N. That
    /// is what lets the hot multiplier stay at 4 on the corrected anchor.
    /// </summary>
    [Theory]
    [Trait("Category", "Parity")]
    [InlineData("toy", 10)]
    [InlineData("iris", 5)]
    [InlineData("iris", 10)]
    [InlineData("landsat", 5)]
    public void Bracket_HotEnd_LeavesNoClusterMass(string name, int k)
    {
        CsrGraph graph = ParityGraph(name, k);
        (_, double hi) = BwdPottsCriticalEstimate.Bracket(20);
        double s2 = S2(graph, hi, 20, new RunBudget(200, 400));
        double overFloor = s2 * graph.NodeCount;
        output.WriteLine($"{name} k={k}: T_hi={hi:F4} S2={s2:F5} S2*N={overFloor:F3}");
        Assert.True(overFloor <= 2.0, $"S2*N = {overFloor:F3} at the hot end");
    }

    /// <summary>
    /// The bracket's cold end holds the toy's stripes joined (S₂ above the plateau the split falls to, 0.25–0.64): the
    /// clipped cold edge the misplaced logarithm produced (stopping primer, C22 and C44).
    /// </summary>
    [Fact]
    [Trait("Category", "Parity")]
    public void Bracket_ColdEnd_HoldsTheToysStripesJoined()
    {
        CsrGraph graph = ParityGraph("toy", 10);
        (double lo, _) = BwdPottsCriticalEstimate.Bracket(20);
        double s2 = S2(graph, lo, 20, new RunBudget(1000, 1000));
        output.WriteLine($"toy k=10: T_lo={lo:F5} S2={s2:F3}");
        Assert.True(s2 >= 0.75, $"S2 = {s2:F3} at the cold end");
    }

    private static CsrGraph ParityGraph(string name, int k)
    {
        double[][] features = name switch
        {
            "toy" => SpcUserDataset.FromSyntheticDataset(Synthetic.Euclidean.Bwd1995Toy.Generate(seed: 42), null).Features,
            "iris" => SpcUserSession.FromCsv(Datasets.Path("iris.csv"), null, false, ',').Features,
            "landsat" => SpcUserSession.FromCsv(Datasets.Path("landsat.csv"), null, false, ',').Features,
            _ => throw new ArgumentOutOfRangeException(nameof(name)),
        };
        return SpcGraphBuilder.BuildResult(features, BwdParityTests.BuildGraphConfig(k, bandwidth: 0)).Graph;
    }

    private static double S2(CsrGraph graph, double temperature, int q, RunBudget budget)
    {
        Accumulator acc = SwRunner.Run(new SwRunSpec
        {
            Graph = graph, Temperature = temperature, Q = q, ReplicaIndex = 0,
            Seed = SpcSeedHelper.Derive(Seed, temperature, 0),
            Budget = budget,
        }).Accumulator;
        double n = graph.NodeCount;
        return acc.RunningSumSqClusterSizes / acc.DrawCount / (n * n);
    }

    /// <summary>The T where S₂ falls through half its coldest value, interpolated in ln T.</summary>
    private static double HalfPoint(double[] grid, double[] s2)
    {
        double half = 0.5 * s2[0];
        for (int i = 1; i < grid.Length; i++)
        {
            if (s2[i] > half) continue;
            double w = (s2[i - 1] - half) / (s2[i - 1] - s2[i]);
            return Math.Exp(Math.Log(grid[i - 1]) + w * (Math.Log(grid[i]) - Math.Log(grid[i - 1])));
        }
        return double.NaN;
    }
}
