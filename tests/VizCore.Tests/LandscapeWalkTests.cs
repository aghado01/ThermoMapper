using System;
using System.Linq;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// R1 resolution layer: the Landscape carrier + the per-leaf walk
/// (cluster profiles → integrate / peak). Oracles are HAND-COMPUTED from the
/// walk math (landscape-carrier-and-walk.md), not derived from any in-repo
/// algorithm — validation independence. The L ≡ 1 fixture pins the
/// consistency identity M(C) = |C| × lifetime (classic HDBSCAN stability).
/// </summary>
public sealed class LandscapeWalkTests
{
    // Four leaves; (0,1) merge at h=1 → id4; (2,3) at h=2 → id5; root at h=4.
    private static Dendrogram FourLeafTree(string axis = "temperature") => new(
        new[]
        {
            new DendrogramNode(0, 1, 1.0, 2),
            new DendrogramNode(2, 3, 2.0, 2),
            new DendrogramNode(4, 5, 4.0, 4),
        },
        LeafCount: 4,
        CostAxis: axis);

    private static Landscape ConstantLandscape(double value, int nodes = 4) => Landscape.Create(
        "temperature",
        new[] { 1.0, 2.0, 3.0, 4.0 },
        Enumerable.Range(0, 4).Select(_ => Enumerable.Repeat(value, nodes).ToArray()).ToArray());

    [Fact]
    public void ClusterProfiles_ConstantLandscape_MassIsSizeTimesLifetime()
    {
        // Grid [1,2,3,4], widths [1,1,1,1(replicated)]. L ≡ 1:
        //   id4: birth 1, death 4 → cells at 1,2,3 → M = 2·3 = 6 = |C|·lifetime
        //   id5: birth 2, death 4 → cells at 2,3   → M = 2·2 = 4
        //   root: birth 4, death ∞ → cell at 4     → M = 4·1 = 4
        var report = LandscapeWalk.ClusterProfiles(FourLeafTree(), ConstantLandscape(1.0));

        Assert.Equal(6.0, report.Mass[0], precision: 12);
        Assert.Equal(4.0, report.Mass[1], precision: 12);
        Assert.Equal(4.0, report.Mass[2], precision: 12);
        Assert.Equal(1.0, report.Birth[0], precision: 12);
        Assert.Equal(4.0, report.Death[0], precision: 12);
        Assert.True(double.IsPositiveInfinity(report.Death[2]));
        // Constant landscape: density 1 everywhere; peak = first cell in the lifetime.
        Assert.Equal(0, report.PeakGridIndex[0]);
        Assert.Equal(1.0, report.PeakDensity[0], precision: 12);
    }

    [Fact]
    public void ClusterProfiles_WeightedLandscape_MassesAreHandComputed()
    {
        // L[g][p] = p+1 for every grid point.
        //   id4 = {0,1}: S = 3 → M = 3·3 = 9;  density 1.5
        //   id5 = {2,3}: S = 7 → M = 7·2 = 14; density 3.5
        //   root:        S = 10 → M = 10·1 = 10
        var values = Enumerable.Range(0, 4)
            .Select(_ => new[] { 1.0, 2.0, 3.0, 4.0 })
            .ToArray();
        var landscape = Landscape.Create("temperature", new[] { 1.0, 2.0, 3.0, 4.0 }, values);

        var report = LandscapeWalk.ClusterProfiles(FourLeafTree(), landscape);

        Assert.Equal(9.0,  report.Mass[0], precision: 12);
        Assert.Equal(14.0, report.Mass[1], precision: 12);
        Assert.Equal(10.0, report.Mass[2], precision: 12);
        Assert.Equal(1.5, report.PeakDensity[0], precision: 12);
        Assert.Equal(3.5, report.PeakDensity[1], precision: 12);
    }

    [Fact]
    public void ClusterProfiles_VaryingProfile_PeakFindsTheArgmaxCell()
    {
        // Two leaves, one root merged at h=1; grid [1,2,3]. Columns give
        // S(root,·) = 2, 10, 4 → peak at grid index 1, density 5.
        var tree = new Dendrogram(new[] { new DendrogramNode(0, 1, 1.0, 2) }, LeafCount: 2, CostAxis: "temperature");
        var landscape = Landscape.Create(
            "temperature",
            new[] { 1.0, 2.0, 3.0 },
            new[] { new[] { 1.0, 1.0 }, new[] { 5.0, 5.0 }, new[] { 2.0, 2.0 } });

        var report = LandscapeWalk.ClusterProfiles(tree, landscape);

        Assert.Equal(1, report.PeakGridIndex[0]);
        Assert.Equal(5.0, report.PeakDensity[0], precision: 12);
        Assert.Equal(16.0, report.Mass[0], precision: 12); // 2·1 + 10·1 + 4·1
    }

    [Fact]
    public void SelectByExcessOfMass_SelectsChildrenAndResolvesAssignment()
    {
        var tree = FourLeafTree();
        var report = LandscapeWalk.ClusterProfiles(tree, ConstantLandscape(1.0));

        bool[] selected = LandscapeWalk.SelectByExcessOfMass(tree, report.Mass);

        // Root ineligible by default; both children carry positive mass over empty subtrees.
        Assert.Equal(new[] { true, true, false }, selected);

        Assignment assignment = LandscapeWalk.ToAssignment(tree, selected);
        Assert.Equal(2, assignment.Count);
        Assert.Equal(assignment.Labels[0], assignment.Labels[1]);
        Assert.Equal(assignment.Labels[2], assignment.Labels[3]);
        Assert.NotEqual(assignment.Labels[0], assignment.Labels[2]);
        Assert.DoesNotContain(Assignment.Unassigned, assignment.Labels);
    }

    [Fact]
    public void SelectByExcessOfMass_ParentBeatsWeakChildren()
    {
        // Force the parent's mass above the child sum: select the root's
        // children's parent instead of the leaves' pairs. Tree of 4 with
        // masses overridden directly (selection reads masses only).
        var tree = FourLeafTree();
        var mass = new[] { 1.0, 1.0, 5.0 };

        bool[] selected = LandscapeWalk.SelectByExcessOfMass(tree, mass, allowRoot: true);

        Assert.Equal(new[] { false, false, true }, selected);

        Assignment assignment = LandscapeWalk.ToAssignment(tree, selected);
        Assert.Equal(1, assignment.Count);
        Assert.All(assignment.Labels, label => Assert.Equal(0, label));
    }

    [Fact]
    public void ClusterProfiles_AxisMismatch_Throws()
    {
        var tree = FourLeafTree(axis: "mutual_reachability_distance");
        var ex = Assert.Throws<InvalidOperationException>(
            () => LandscapeWalk.ClusterProfiles(tree, ConstantLandscape(1.0)));
        Assert.Contains("Axis-alignment", ex.Message);
    }

    [Fact]
    public void LandscapeCreate_ValidatesShape()
    {
        Assert.Throws<ArgumentException>(() => Landscape.Create(
            "temperature", new[] { 2.0, 1.0 }, new[] { new[] { 1.0 }, new[] { 1.0 } }));   // not ascending
        Assert.Throws<ArgumentException>(() => Landscape.Create(
            "temperature", new[] { 1.0, 2.0 }, new[] { new[] { 1.0 } }));                  // column count mismatch
        Assert.Throws<ArgumentException>(() => Landscape.Create(
            "temperature", new[] { 1.0, 2.0 }, new[] { new[] { 1.0 }, new[] { 1.0, 2.0 } })); // ragged
    }

    [Fact]
    public void SweepLandscapes_FromFrames_PoolsAcrossReplicasAndStacksByTemperature()
    {
        // Two replicas at T=1.0 concentrating the per-node sums on opposite
        // nodes: pooled column = (10+0)/20, (0+10)/20 = [0.5, 0.5] — linear
        // pooling is exact. One run at T=0.5 → column [1, 1]. Grid ascends.
        var frames = new[]
        {
            MakeFrame(1.0, draws: 10, perNode: new[] { 10.0, 0.0 }),
            MakeFrame(1.0, draws: 10, perNode: new[] { 0.0, 10.0 }),
            MakeFrame(0.5, draws: 4,  perNode: new[] { 4.0, 4.0 }),
        };

        Landscape landscape = SweepLandscapes.FromFrames(frames, SwLandscapeSink.MeanClusterSize, graphId: "test");

        Assert.Equal("temperature", landscape.Axis);
        Assert.Equal(new[] { 0.5, 1.0 }, landscape.Grid);
        Assert.Equal(new[] { 1.0, 1.0 }, landscape.ValuesByGridPoint[0]);
        Assert.Equal(new[] { 0.5, 0.5 }, landscape.ValuesByGridPoint[1]);
        Assert.Equal("MeanClusterSize", landscape.Provenance!.Sink);
        Assert.Equal("test", landscape.Provenance.GraphId);
    }

    private static Accumulator MakeFrame(double temperature, int draws, double[] perNode) => new()
    {
        Temperature = temperature,
        Q = 4,
        DrawCount = draws,
        Spins = new int[perNode.Length],
        ClusterSizeHistogram = new int[perNode.Length],
        RngState0 = 1,
        RngState1 = 2,
        RngState2 = 3,
        RngState3 = 4,
        RunningSumSqClusterSizes = 0.0,
        RunningSumSqClusterSizesExcl = 0.0,
        RunningSumEnergy = 0.0,
        RunningSumEnergySq = 0.0,
        RunningSumMag = 0.0,
        RunningSumMagSq = 0.0,
        SumClusterSizePerNode = perNode,
    };
}
