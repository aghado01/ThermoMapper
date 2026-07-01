using System;
using System.Threading.Tasks;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.Geometry;
using Repo.TestHarness;
using Synthetic;
using Synthetic.Euclidean;
using Viz;
using Xunit;

namespace VizCore.Tests;

internal static class TestTimeouts
{
    public const int HeavyFixtureMs = 300_000;
}

public sealed class IntrinsicFlowTests
{
    [Fact]
    public void SpectralBridge_LineFieldFromEigenvector_OnPathGraph_AlignsWithAxis()
    {
        double[][] points =
        {
            new[] { 0.0, 0.0 },
            new[] { 1.0, 0.0 },
            new[] { 2.0, 0.0 },
            new[] { 3.0, 0.0 },
        };
        Edge[] edges =
        {
            new(0, 1, 1.0),
            new(1, 2, 1.0),
            new(2, 3, 1.0),
        };

        CsrGraph graph = CsrGraph.FromEdges(edges, points.Length);
        double[] eigenvector = { -1.5, -0.5, 0.5, 1.5 };

        double[] lineField = SpectralBridge.LineFieldFromEigenvector(graph, points, eigenvector);

        for (int i = 0; i < points.Length; i++)
        {
            int offset = i * 2;
            double x = lineField[offset];
            double y = lineField[offset + 1];
            double norm = Math.Sqrt((x * x) + (y * y));

            Assert.InRange(norm, 0.999999, 1.000001);
            Assert.True(Math.Abs(x) > 0.999999, $"Expected node {i} to align with the x-axis, got ({x}, {y}).");
            Assert.InRange(Math.Abs(y), 0.0, 1e-12);
        }
    }

    [Fact]
    public void LocalTangent_EuclideanManifoldMatchesFlatPath()
    {
        double[][] points =
        {
            new[] { 0.0, 0.0 },
            new[] { 1.0, 0.1 },
            new[] { 2.0, 0.0 },
            new[] { 3.0, -0.1 },
            new[] { 4.0, 0.0 },
        };
        int[][] adjacency =
        {
            new[] { 1, 2 },
            new[] { 0, 2, 3 },
            new[] { 0, 1, 3, 4 },
            new[] { 1, 2, 4 },
            new[] { 2, 3 },
        };

        double[] flat = LocalTangent.Compute(points, adjacency);
        double[] metricAware = LocalTangent.Compute(points, adjacency, new EuclideanVectorManifold(points[0].Length));

        for (int i = 0; i < points.Length; i++)
        {
            double agreement = UnorientedAgreement(flat, metricAware, i, points[0].Length);
            Assert.InRange(agreement, 0.999999999, 1.0);
        }
    }

    // Producer-vs-producer diagnostic: records SpectralBridge / LocalTangent agreement without asserting a threshold.
    // Two data-driven producers, neither knowing the generator, agreeing on orientation is evidence the bridge is correct.
    // Disagreement is also data. See POC-planning.md "producer-vs-producer agreement" framing.
    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Crescent_SpectralBridgeAndLocalTangent_AgreementDiagnostic()
    {
        await Task.Run(() =>
        {
            var dataset = CrescentAndEllipsoid.Generate(crescentPoints: 1500, ellipsoidPoints: 0, seed: 42);
            CsrGraph graph = FixtureGoldenHelpers.BuildDraftGraph(dataset.Features, k: FixtureGoldenHelpers.ConnectedFixtureK, ensureConnected: true);
            int[][] adjacency = BuildAdjacency(graph);

            double[] spectral = SpectralBridge.LineFieldFromFiedler(graph, dataset.Features, seed: FixtureGoldenHelpers.SpectralSeed);
            double[] tangent = LocalTangent.Compute(dataset.Features, adjacency);
            double meanAgreement = MeanUnorientedAgreement(spectral, tangent, dataset.Features[0].Length);

            ArtifactRun run = HarnessArtifacts.Create(
                runKind: "test-runs",
                suiteName: nameof(IntrinsicFlowTests),
                runName: nameof(Crescent_SpectralBridgeAndLocalTangent_AgreementDiagnostic));
            string analysisPath = run.WriteRunJson("analysis", new
            {
                MeanUnorientedAgreement = meanAgreement,
                N = dataset.Features.Length,
                ProducerA = "SpectralBridge.LineFieldFromFiedler",
                ProducerB = "LocalTangent.Compute (flat Euclidean)",
                Note = "diagnostic only — no threshold; disagreement is also informative",
            });
            Console.WriteLine($"RunRoot\t{run.RunDirectory}");
            Console.WriteLine($"Manifest\t{run.ManifestPath}");
            Console.WriteLine($"Analysis\t{analysisPath}");
            Console.WriteLine($"MeanUnorientedAgreement\t{meanAgreement:F6}");
        });
    }

    private static int[][] BuildAdjacency(CsrGraph graph)
    {
        var adjacency = new int[graph.NodeCount][];
        for (int i = 0; i < graph.NodeCount; i++)
        {
            int start = graph.RowPointers[i];
            int length = graph.RowPointers[i + 1] - start;
            var neighbors = new int[length];
            Array.Copy(graph.Targets, start, neighbors, 0, length);
            adjacency[i] = neighbors;
        }

        return adjacency;
    }

    private static double MeanUnorientedAgreement(double[] first, double[] second, int d)
    {
        int n = first.Length / d;
        double sum = 0.0;
        int count = 0;

        for (int i = 0; i < n; i++)
        {
            double agreement = UnorientedAgreement(first, second, i, d);
            if (double.IsNaN(agreement))
                continue;

            sum += agreement;
            count++;
        }

        return count == 0 ? 0.0 : sum / count;
    }

    private static double UnorientedAgreement(double[] first, double[] second, int index, int d)
    {
        int offset = index * d;
        double dot = 0.0;
        double normFirst = 0.0;
        double normSecond = 0.0;

        for (int dim = 0; dim < d; dim++)
        {
            double a = first[offset + dim];
            double b = second[offset + dim];
            dot += a * b;
            normFirst += a * a;
            normSecond += b * b;
        }

        if (normFirst <= 1e-20 || normSecond <= 1e-20)
            return double.NaN;

        return Math.Abs(dot / Math.Sqrt(normFirst * normSecond));
    }
}
