using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Models.Potts.Observables;
using Graphs.Primitives;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace VizCore.Tests;

/// <summary>
/// The Binder cumulant U = 1 − ⟨m⁴⟩/(3⟨m²⟩²) (S2), a view over the per-draw magnetization series, assembled per
/// temperature from the replicas' pooled moments. Its validation is the crossing of U across lattice sizes on the
/// q = 2 Potts model, at the exact T_c = 1/ln(1 + √2) for H = −Σ J δ with J = 1.
/// </summary>
public sealed class BinderCumulantTests(ITestOutputHelper output)
{
    // ── The reduction ────────────────────────────────────────────────────────

    [Fact]
    public void Reduce_ConstantMagnetization_IsTwoThirds()
        => Assert.InRange(Math.Abs(BinderCumulant.Reduce([0.5, 0.5, 0.5]) - 2.0 / 3.0), 0.0, 1e-15);

    [Fact]
    public void Reduce_HalfOrderedHalfZero_IsOneThird()
        // ⟨m²⟩ = ⟨m⁴⟩ = 1/2, so U = 1 − (1/2)/(3/4).
        => Assert.InRange(Math.Abs(BinderCumulant.Reduce([1.0, 0.0, 1.0, 0.0]) - 1.0 / 3.0), 0.0, 1e-15);

    [Fact]
    public void Reduce_IsUndefinedWithoutPositiveSecondMoment()
    {
        Assert.True(double.IsNaN(BinderCumulant.Reduce(ReadOnlySpan<double>.Empty)));
        Assert.True(double.IsNaN(BinderCumulant.Reduce([0.0, 0.0])));
    }

    [Fact]
    public void Moments_AreTheSeriesMeans()
    {
        BinderMoments moments = BinderCumulant.Moments([0.5, 1.0]);

        Assert.Equal((0.25 + 1.0) / 2, moments.SecondMoment);
        Assert.Equal((0.0625 + 1.0) / 2, moments.FourthMoment);
    }

    // ── The assembly ─────────────────────────────────────────────────────────

    [Fact]
    public void SweepProfile_PoolsTheReplicasMoments_NotTheirCumulants()
    {
        // Each replica alone has constant m, so U = 2/3 for each; pooled, ⟨m²⟩ = 0.625 and ⟨m⁴⟩ = 0.53125.
        SweepProfile profile = SweepProfile.From([Run(1.0, [1.0, 1.0]), Run(1.0, [0.5, 0.5])]);

        double expected = 1.0 - 0.53125 / (3.0 * 0.625 * 0.625);
        Assert.InRange(Math.Abs(profile.AdditionalChannels["BinderCumulant"][0] - expected), 0.0, 1e-15);
    }

    [Fact]
    public void SweepProfile_IsNaN_WhereAFrameCarriesNoRecord()
    {
        SpcRunResult withRecord = Run(1.0, [1.0, 0.5]);
        SpcRunResult sumsAlone = withRecord with
        {
            Accumulator = withRecord.Accumulator with { Temperature = 2.0, MagPerDraw = [] },
        };

        IReadOnlyList<double> u = SweepProfile.From([withRecord, sumsAlone]).AdditionalChannels["BinderCumulant"];

        Assert.False(double.IsNaN(u[0]));
        Assert.True(double.IsNaN(u[1]));
    }

    // ── The check: the q = 2 lattice crossing ────────────────────────────────

    // Pre-registered before the first run (2026-10-29), from a timing pass alone: each pairwise crossing within 0.02 of
    // T_c, on this grid, budget and seed. The crossing is the linear interpolation at the one sign change of
    // U_larger − U_smaller, which runs from positive (ordered side) to negative.
    private const double CrossingTolerance = 0.02;
    private static readonly double CriticalTemperature = 1.0 / Math.Log(1.0 + Math.Sqrt(2.0));
    private static readonly double[] CrossingGrid = [1.09, 1.11, 1.13, 1.15, 1.17, 1.19];
    private static readonly (int L, int Cycles)[] Lattices = [(8, 25_000), (16, 12_500), (32, 6_250)];
    private const int CrossingReplicas = 4;
    private const int CrossingBurnIn = 1_000;
    private const int CrossingSeed = 20261029;

    [Fact]
    public void Crossings_OnTheIsingLattice_SitAtTheExactCriticalTemperature()
    {
        var clock = Stopwatch.StartNew();
        double[][] curves = BinderCurves();
        clock.Stop();

        output.WriteLine($"T_c = {CriticalTemperature:F5}; tolerance {CrossingTolerance}; runtime {clock.Elapsed.TotalSeconds:F1} s");
        for (int k = 0; k < Lattices.Length; k++)
            output.WriteLine($"L = {Lattices[k].L,2}: U = " + string.Join(", ", curves[k].Select(u => u.ToString("F4"))));

        for (int a = 0; a < Lattices.Length; a++)
            for (int b = a + 1; b < Lattices.Length; b++)
            {
                double crossing = Crossing(curves[a], curves[b]);
                output.WriteLine($"L = {Lattices[a].L} × {Lattices[b].L}: crossing {crossing:F4}, off T_c by {crossing - CriticalTemperature:+0.0000;-0.0000}");
                Assert.InRange(crossing, CriticalTemperature - CrossingTolerance, CriticalTemperature + CrossingTolerance);
            }
    }

    /// <summary>U over <see cref="CrossingGrid"/> for each lattice, through <see cref="SweepProfile.From"/>.</summary>
    private static double[][] BinderCurves()
    {
        int[] latticeSeeds = SeedTree.Derive(CrossingSeed, Lattices.Length);
        CsrGraph[] graphs = Lattices.Select(l => PeriodicLattice(l.L)).ToArray();
        var tasks = (from k in Enumerable.Range(0, Lattices.Length)
                     from t in CrossingGrid
                     from r in Enumerable.Range(0, CrossingReplicas)
                     select (k, t, r)).ToArray();

        var runs = new SpcRunResult[tasks.Length];
        Parallel.For(0, tasks.Length, i =>
        {
            (int k, double t, int r) = tasks[i];
            runs[i] = new SpcRunResult
            {
                Graph = graphs[k],
                Accumulator = SwRunner.Run(new SwRunSpec
                {
                    Graph = graphs[k], Temperature = t, Q = 2, ReplicaIndex = r,
                    Seed = SpcSeedHelper.Derive(latticeSeeds[k], t, r),
                    Budget = new RunBudget(CrossingBurnIn, Lattices[k].Cycles),
                }).Accumulator,
            };
        });

        return Enumerable.Range(0, Lattices.Length).Select(k =>
        {
            SweepProfile profile = SweepProfile.From(runs.Where((_, i) => tasks[i].k == k).ToList());
            Assert.Equal(CrossingGrid, profile.Temperatures);
            return profile.AdditionalChannels["BinderCumulant"].ToArray();
        }).ToArray();
    }

    private static double Crossing(double[] smaller, double[] larger)
    {
        double[] d = larger.Zip(smaller, (l, s) => l - s).ToArray();
        int[] changes = Enumerable.Range(0, d.Length - 1).Where(i => Math.Sign(d[i]) != Math.Sign(d[i + 1])).ToArray();

        Assert.True(changes.Length == 1, $"expected one sign change of U_larger − U_smaller, got {changes.Length}: [{string.Join(", ", d.Select(x => x.ToString("F4")))}]");
        int c = changes[0];
        Assert.True(d[c] > 0 && d[c + 1] < 0, $"the sign change at T = {CrossingGrid[c]} does not run from ordered to disordered");
        return CrossingGrid[c] + d[c] / (d[c] - d[c + 1]) * (CrossingGrid[c + 1] - CrossingGrid[c]);
    }

    /// <summary>The L × L square lattice with periodic boundaries and unit couplings (L ≥ 3, so no edge repeats).</summary>
    internal static CsrGraph PeriodicLattice(int L)
    {
        var edges = new Edge[2 * L * L];
        for (int y = 0; y < L; y++)
            for (int x = 0; x < L; x++)
            {
                int i = y * L + x;
                edges[2 * i]     = new Edge(i, y * L + (x + 1) % L, 1.0);
                edges[2 * i + 1] = new Edge(i, (y + 1) % L * L + x, 1.0);
            }
        return CsrGraph.FromEdges(edges, L * L);
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private static SpcRunResult Run(double temperature, double[] magPerDraw) => new()
    {
        Graph = CsrGraph.FromEdges([new Edge(0, 1, 1.0)], 2),
        Accumulator = new Accumulator
        {
            Temperature = temperature, Q = 2, DrawCount = magPerDraw.Length,
            Spins = new int[2], ClusterSizeHistogram = new long[] { 0, (long)magPerDraw.Length },
            RngState0 = 1, RngState1 = 2, RngState2 = 3, RngState3 = 4,
            RunningSumSqClusterSizes = 4.0 * magPerDraw.Length, RunningSumSqClusterSizesExcl = 0,
            RunningSumEnergy = 0, RunningSumEnergySq = 0,
            RunningSumMag = magPerDraw.Sum(), RunningSumMagSq = magPerDraw.Sum(m => m * m),
            MagPerDraw = magPerDraw,
        },
    };
}
