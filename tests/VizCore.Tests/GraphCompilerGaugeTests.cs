using System;
using System.Linq;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Models.Potts;
using Graphs.Pipeline;
using Graphs.Pipeline.Scalers;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The gauge is a global factor applied last, so it only rescales temperature
/// (graph-compiler-brief §4.1, decision 2): it changes no threshold's outcome, and
/// <c>(J, T)</c> and <c>(gJ, gT)</c> give the same bond probabilities.
/// </summary>
public sealed class GraphCompilerGaugeTests
{
    private static double[][] Points()
    {
        var points = new double[24][];
        for (int i = 0; i < 12; i++)
        {
            points[i] = new[] { 0.1 * i + 0.01 * i * i, 0.0 };
            points[12 + i] = new[] { 8.0 + 0.12 * i + 0.012 * i * i, 0.5 };
        }
        return points;
    }

    private static GraphBuildResult Build(GaugeKind gauge, BandwidthStrategy? strategy = null, IKernelDescriptor? kernel = null)
    {
        double[][] p = Points();
        var config = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 4 },
            Filter = new FilterConfig { Kind = FilterKind.OrRule },
            Repair = new RepairConfig { Kind = RepairKind.MstMin },
            Projection = new CouplingProjection
            {
                Kernel = kernel ?? new Cauchy(),
                LmpRescale = false,
                BandwidthOverride = strategy,
                Gauge = gauge,
            },
        };
        return GraphCompiler.Build(config, p.Length, new GraphMetric((i, j) =>
        {
            double dx = p[i][0] - p[j][0], dy = p[i][1] - p[j][1];
            return Math.Sqrt(dx * dx + dy * dy);
        }));
    }

    [Fact]
    public void Gauge_RescalesEveryCouplingByOneOverTheReturnedGraphsMeanDegree_AndNothingElse()
    {
        CsrGraph plain = Build(GaugeKind.None).Graph;
        CsrGraph gauged = Build(GaugeKind.InverseMeanDegree).Graph;

        Assert.Equal(plain.Targets, gauged.Targets);
        Assert.Equal(plain.RowPointers, gauged.RowPointers);
        double factor = plain.NodeCount / (double)plain.Targets.Length;
        for (int slot = 0; slot < plain.Weights.Length; slot++)
            Assert.InRange(gauged.Weights[slot] / (plain.Weights[slot] * factor), 1 - 1e-12, 1 + 1e-12);
    }

    [Fact]
    public void MeanEdgeDistanceBandwidth_NoLongerImpliesTheGauge()
    {
        CsrGraph plain = Build(GaugeKind.None, BandwidthStrategy.MeanEdgeDistance).Graph;
        CsrGraph gauged = Build(GaugeKind.InverseMeanDegree, BandwidthStrategy.MeanEdgeDistance).Graph;

        // The same bandwidth either way; only the gauge differs.
        Assert.All(Enumerable.Range(0, plain.Weights.Length), slot => Assert.True(gauged.Weights[slot] < plain.Weights[slot]));
    }

    [Fact]
    public void Gauge_IsExactlyATemperatureRescale_ForTheBondProbabilitiesSwSamples()
    {
        CsrGraph plain = Build(GaugeKind.None, BandwidthStrategy.MeanEdgeDistance, new Gaussian()).Graph;
        GraphBuildResult gaugedBuild = Build(GaugeKind.InverseMeanDegree, BandwidthStrategy.MeanEdgeDistance, new Gaussian());
        CsrGraph gauged = gaugedBuild.Graph;
        double g = plain.NodeCount / (double)plain.Targets.Length;

        foreach (double temperature in new[] { 0.05, 0.2, 1.0 })
            for (int slot = 0; slot < plain.Weights.Length; slot++)
            {
                double atPlain = FkKernel.BondProbability(plain.Weights[slot], temperature);
                double atGauged = FkKernel.BondProbability(gauged.Weights[slot], g * temperature);
                Assert.InRange(atGauged, atPlain - 1e-12, atPlain + 1e-12);
            }
    }

    [Fact]
    public void Gauge_DoesNotShiftWhichEdgesSurviveTheFloor()
    {
        // K4 with one edge at raw coupling 2e-10 (d² / 2 = ln(1/2e-10)). K-hat = 3, so the gauged value is
        // 6.7e-11: under the floor if the gauge came first, kept now because the floor sees the raw coupling.
        double d = Math.Sqrt(2.0 * Math.Log(1.0 / 2e-10));
        var rows = new[]
        {
            new[] { Near(1, 0.1), Near(2, 0.1), Near(3, d) },
            new[] { Near(0, 0.1), Near(2, 0.1), Near(3, 0.1) },
            new[] { Near(0, 0.1), Near(1, 0.1), Near(3, 0.1) },
            new[] { Near(0, d), Near(1, 0.1), Near(2, 0.1) },
        };
        var selection = new NeighborSelection(rows, new double[4], new[] { d, 0.1, 0.1, d });

        ScalerResult scaled = new GaugeScaler(
            new GlobalBandwidthScaler(KernelType.Gaussian, bandwidth: 1.0), GaugeKind.InverseMeanDegree).Scale(selection, 4);

        Assert.Equal(0, scaled.CutEdges);
        Assert.Equal(3.0, 1.0 / scaled.GaugeFactor, 12);
        double kept = double.NaN;
        for (int slot = scaled.Graph.RowPointers[0]; slot < scaled.Graph.RowPointers[1]; slot++)
            if (scaled.Graph.Targets[slot] == 3) kept = scaled.Graph.Weights[slot];
        Assert.InRange(kept / (2e-10 / 3.0), 1 - 1e-9, 1 + 1e-9);
    }

    [Fact]
    public void Gauge_FollowsLmp_AndCountsTheZeroEdgesLmpWritesInTheMeanDegree()
    {
        // Strict LMP writes zeros but keeps the topology, so K-hat is still the structural mean degree
        // and every weight, the zeros included, is the ungauged LMP weight times g.
        CsrGraph plain = BuildLmp(GaugeKind.None);
        CsrGraph gauged = BuildLmp(GaugeKind.InverseMeanDegree);

        Assert.Contains(plain.Weights, w => w == 0.0);
        double g = plain.NodeCount / (double)plain.Targets.Length;
        for (int slot = 0; slot < plain.Weights.Length; slot++)
            Assert.InRange(gauged.Weights[slot], plain.Weights[slot] * g * (1 - 1e-12), plain.Weights[slot] * g * (1 + 1e-12));
    }

    private static CsrGraph BuildLmp(GaugeKind gauge)
    {
        double[][] p = Points();
        return GraphCompiler.Build(
            new GraphCompilerConfig
            {
                Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 4 },
                Filter = new FilterConfig { Kind = FilterKind.OrRule },
                Repair = new RepairConfig { Kind = RepairKind.MstMin },
                Projection = new CouplingProjection { Kernel = new Cauchy(), LmpRescale = true, Gauge = gauge },
            },
            p.Length,
            new GraphMetric((i, j) =>
            {
                double dx = p[i][0] - p[j][0], dy = p[i][1] - p[j][1];
                return Math.Sqrt(dx * dx + dy * dy);
            })).Graph;
    }

    private static Neighbor Near(int index, double distance) => new() { Index = index, Distance = distance };
}
