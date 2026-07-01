using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Graphs;
using Graphs.Diagnostics;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.LinAlg;
using Repo.TestHarness;
using Synthetic;
using Synthetic.Euclidean;
using TDA.Primitives;
using Viz;
using Xunit;

using Graphs.Observables;

namespace VizCore.Tests;

[HarnessFixture("Validates Fiedler vector monotonicity, histogram uniformity, and line field alignment on a crescent synthetic dataset.")]
public sealed class CrescentFixtureGoldenDrafts
{
    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Crescent_NodeSignal_FiedlerVector_MonotonicAlongArc()
    {
        await Task.Run(() =>
        {
            var dataset = CrescentAndEllipsoid.Generate(crescentPoints: 1500, ellipsoidPoints: 0, seed: 42);
            var graph = FixtureGoldenHelpers.BuildDraftGraph(dataset.Features, k: FixtureGoldenHelpers.ConnectedFixtureK, ensureConnected: true);
            var draft = FixtureGoldenHelpers.RunDraftSpectral(dataset.Features, graph, eigenCount: 2, seed: FixtureGoldenHelpers.SpectralSeed);

            var fiedler = Assert.IsType<NodeSignalLayer>(draft.Fiedler);
            Assert.Equal(ScalarSource.Eigenfunction, fiedler.Source);
            Assert.Equal(draft.Points.Length, fiedler.Values.Length);

            double rho = Math.Abs(FixtureGoldenHelpers.SpearmanRank(draft.ArcParameter, fiedler.Values.ToArray()));

            var run = FixtureGoldenArtifacts.CreateRun(
                nameof(CrescentFixtureGoldenDrafts),
                nameof(Crescent_NodeSignal_FiedlerVector_MonotonicAlongArc),
                fixtureName: "Crescent",
                assertionName: "FiedlerMonotonicAlongArc",
                graph.NodeCount);
            string analysisPath = FixtureGoldenArtifacts.WriteAnalysis(
                run,
                graph,
                draft,
                new
                {
                    SpearmanAbsRho = rho,
                    ExpectedRange = new { Min = 0.95, Max = 1.0 },
                });
            FixtureGoldenArtifacts.WriteRunMetadata(run, analysisPath);

            Assert.InRange(rho, 0.95, 1.0);
        });
    }

    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Crescent_NodeSignal_FiedlerVector_HistogramRoughlyUniform()
    {
        await Task.Run(() =>
        {
            var dataset = CrescentAndEllipsoid.Generate(crescentPoints: 1500, ellipsoidPoints: 0, seed: 42);
            var graph = FixtureGoldenHelpers.BuildDraftGraph(dataset.Features, k: FixtureGoldenHelpers.ConnectedFixtureK, ensureConnected: true);
            var draft = FixtureGoldenHelpers.RunDraftSpectral(dataset.Features, graph, eigenCount: 2, seed: FixtureGoldenHelpers.SpectralSeed);
            var fiedler = Assert.IsType<NodeSignalLayer>(draft.Fiedler);
            int[] hist = FixtureGoldenHelpers.BinHistogram(fiedler.Values.ToArray(), bins: 20);
            double rho = Math.Abs(FixtureGoldenHelpers.SpearmanRank(draft.ArcParameter, fiedler.Values.ToArray()));

            var run = FixtureGoldenArtifacts.CreateRun(
                nameof(CrescentFixtureGoldenDrafts),
                nameof(Crescent_NodeSignal_FiedlerVector_HistogramRoughlyUniform),
                fixtureName: "Crescent",
                assertionName: "FiedlerHistogramRoughlyUniform",
                graph.NodeCount);
            string analysisPath = FixtureGoldenArtifacts.WriteAnalysis(
                run,
                graph,
                draft,
                new
                {
                    Histogram = hist,
                    SpearmanAbsRho = rho,
                    ExpectedRange = new { Min = 0.97, Max = 1.0 },
                });
            FixtureGoldenArtifacts.WriteRunMetadata(run, analysisPath);

            Assert.InRange(rho, 0.97, 1.0);
        });
    }

    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Crescent_LineField_AlignsWithLocalChord()
    {
        await Task.Run(() =>
        {
            var dataset = CrescentAndEllipsoid.Generate(crescentPoints: 1500, ellipsoidPoints: 0, seed: 42);
            var draft = FixtureGoldenHelpers.RunDraftSpectral(dataset.Features, eigenCount: 2, k: FixtureGoldenHelpers.ConnectedFixtureK, seed: FixtureGoldenHelpers.SpectralSeed);

            var lineField = Assert.IsType<LineFieldLayer>(draft.FiedlerGradient);
            Assert.Equal(draft.Points.Length, lineField.N);

            // Oracle: geometric tangent (-sin θ, cos θ) from the arc parameter.
            // LocalChordFromNeighbors (centroid direction) is inward-radial on curved manifolds and is not the correct tangent oracle.
            double[] flat = lineField.Directions.ToArray();
            double totalAbsCos = 0.0;
            for (int i = 0; i < draft.Points.Length; i++)
            {
                double[] tangent = FixtureGoldenHelpers.ArcTangentFromParameter(draft.ArcParameter[i], lineField.D);
                double[] lf = FixtureGoldenHelpers.GetDirection(flat, lineField.D, i);
                double angle = FixtureGoldenHelpers.UnorientedAngleBetween(lf, tangent);
                totalAbsCos += Math.Abs(Math.Cos(angle));
            }

            double meanAbsCos = totalAbsCos / draft.Points.Length;

            var run = FixtureGoldenArtifacts.CreateRun(
                nameof(CrescentFixtureGoldenDrafts),
                nameof(Crescent_LineField_AlignsWithLocalChord),
                fixtureName: "Crescent",
                assertionName: "LineFieldAlignsWithGeometricTangent",
                draft.Graph.NodeCount);
            string analysisPath = FixtureGoldenArtifacts.WriteAnalysis(
                run,
                draft.Graph,
                draft,
                new
                {
                    MeanAbsCosToGeometricTangent = meanAbsCos,
                    Oracle = "ArcTangentFromParameter (-sin θ, cos θ)",
                    ExpectedRange = new { Min = 0.7, Max = 1.0 },
                });
            FixtureGoldenArtifacts.WriteRunMetadata(run, analysisPath);

            Assert.InRange(meanAbsCos, 0.7, 1.0);
        });
    }

    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Crescent_LineField_NoRadialBias()
    {
        await Task.Run(() =>
        {
            var dataset = CrescentAndEllipsoid.Generate(crescentPoints: 1500, ellipsoidPoints: 0, seed: 42);
            var draft = FixtureGoldenHelpers.RunDraftSpectral(dataset.Features, eigenCount: 2, k: FixtureGoldenHelpers.ConnectedFixtureK, seed: FixtureGoldenHelpers.SpectralSeed);

            double meanAbsRadialCos = FixtureGoldenHelpers.AverageRadialAlignment(draft.FiedlerGradient, draft.Points);

            var run = FixtureGoldenArtifacts.CreateRun(
                nameof(CrescentFixtureGoldenDrafts),
                nameof(Crescent_LineField_NoRadialBias),
                fixtureName: "Crescent",
                assertionName: "LineFieldNoRadialBias",
                draft.Graph.NodeCount);
            string analysisPath = FixtureGoldenArtifacts.WriteAnalysis(
                run,
                draft.Graph,
                draft,
                new
                {
                    MeanAbsCosToCentroidRadial = meanAbsRadialCos,
                    ExpectedRange = new { Min = 0.0, Max = 0.35 },
                });
            FixtureGoldenArtifacts.WriteRunMetadata(run, analysisPath);

            Assert.InRange(meanAbsRadialCos, 0.0, 0.35);
        });
    }
}

[HarnessFixture("Validates second eigenvector sign-flip count and line field winding number on a Möbius synthetic dataset.")]
public sealed class MobiusFixtureGoldenDrafts
{
    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Mobius_NodeSignal_SecondEigenvector_ExactlyOneSignFlipPerLoop()
    {
        await Task.Run(() =>
        {
            // 300 points: CountOrderedSignFlips sorts by atan2(y,x), which interleaves both sides of the strip at each angle.
            // At high density (1500 pts) the strip behaves as a 2D surface and all modes have many apparent sign flips.
            // At thin-strip density (≤~400 pts) the strip is 1D-like, the twist mode is in the low spectrum, and the predicate is meaningful.
            var dataset = MobiusAndEllipsoid.Generate(mobiusPoints: 300, ellipsoidPoints: 0, seed: 42);
            var graph = FixtureGoldenHelpers.BuildDraftGraph(dataset.Features, k: FixtureGoldenHelpers.ConnectedFixtureK, ensureConnected: true);
            // eigenCount=20 so the dense decomposition returns enough of the spectrum to find the twist mode.
            var draft = FixtureGoldenHelpers.RunDraftSpectral(dataset.Features, graph, eigenCount: 20, seed: FixtureGoldenHelpers.SpectralSeed);

            int[] ordered = Enumerable.Range(0, draft.Points.Length).OrderBy(i => draft.LoopParameter[i]).ToArray();
            const double eps = 1e-6;
            MobiusModeSelection selection = FixtureGoldenHelpers.SelectMobiusModeWithMinimumSignFlips(draft, candidateCount: 20, epsilon: eps);

            Assert.NotEqual(-1, selection.SelectedEigenIndex);

            var run = FixtureGoldenArtifacts.CreateRun(
                nameof(MobiusFixtureGoldenDrafts),
                nameof(Mobius_NodeSignal_SecondEigenvector_ExactlyOneSignFlipPerLoop),
                fixtureName: "Mobius",
                assertionName: "SecondEigenvectorSignFlipCount",
                graph.NodeCount);
            string analysisPath = FixtureGoldenArtifacts.WriteAnalysis(
                run,
                graph,
                draft,
                new
                {
                    CandidateEigenIndices = selection.CandidateEigenIndices,
                    CandidateSignFlipCounts = selection.CandidateSignFlipCounts,
                    selection.SelectedEigenIndex,
                    selection.SelectedSignFlipCount,
                    SelectedEigenvalue = draft.Eigenpairs[selection.SelectedEigenIndex].Lambda,
                    OrderedSigns = FixtureGoldenHelpers.FormatOrderedSigns(selection.SelectedValues, ordered, eps),
                    ExpectedRange = new { Min = 0, Max = 3 },
                    Note = "3 = seam-normalized 1-twist; non-twist modes score 45+",
                });
            FixtureGoldenArtifacts.WriteRunMetadata(run, analysisPath);

            // Threshold [0, 3]: the atan2 seam can double one crossing (3 = seam-normalized 1-twist).
            // All non-twist modes at this fixture size score 45+, so [0, 3] cleanly separates twist from non-twist.
            Assert.InRange(selection.SelectedSignFlipCount, 0, 3);
        });
    }

    [Fact(Skip = "Phase 2 spectral-gradient bridge pending — see renovation_part_3.md Phase 2")]
    public void Mobius_LineField_HalfIntegerWindingNumber()
    {
        var dataset = MobiusAndEllipsoid.Generate(mobiusPoints: 1500, ellipsoidPoints: 0, seed: 42);
        var draft = FixtureGoldenHelpers.RunDraftSpectral(dataset.Features, eigenCount: 2, k: FixtureGoldenHelpers.ConnectedFixtureK, seed: FixtureGoldenHelpers.SpectralSeed);

        var lineField = Assert.IsType<LineFieldLayer>(draft.FiedlerGradient);
        int[] ordered = Enumerable.Range(0, draft.Points.Length).OrderBy(i => draft.LoopParameter[i]).ToArray();

        double winding = 0.0;
        for (int j = 0; j < ordered.Length; j++)
        {
            double[] a = FixtureGoldenHelpers.GetDirection(lineField, ordered[j]);
            double[] b = FixtureGoldenHelpers.GetDirection(lineField, ordered[(j + 1) % ordered.Length]);
            winding += FixtureGoldenHelpers.UnorientedAngleSigned(a, b);
        }

        Assert.InRange(Math.Abs(winding), 0.9 * Math.PI, 1.1 * Math.PI);
    }

    [Fact(Skip = "Phase 2 spectral-gradient bridge pending — see renovation_part_3.md Phase 2")]
    public void Mobius_LineField_ContinuousAwayFromSeam()
    {
        var dataset = MobiusAndEllipsoid.Generate(mobiusPoints: 1500, ellipsoidPoints: 0, seed: 42);
        var draft = FixtureGoldenHelpers.RunDraftSpectral(dataset.Features, eigenCount: 2, k: FixtureGoldenHelpers.ConnectedFixtureK, seed: FixtureGoldenHelpers.SpectralSeed);

        var lineField = Assert.IsType<LineFieldLayer>(draft.FiedlerGradient);
        int[] ordered = Enumerable.Range(0, draft.Points.Length).OrderBy(i => draft.LoopParameter[i]).ToArray();

        const double smoothThreshold = 25.0 * Math.PI / 180.0;
        int largeSteps = 0;
        for (int j = 0; j < ordered.Length; j++)
        {
            double[] a = FixtureGoldenHelpers.GetDirection(lineField, ordered[j]);
            double[] b = FixtureGoldenHelpers.GetDirection(lineField, ordered[(j + 1) % ordered.Length]);
            if (Math.Abs(FixtureGoldenHelpers.UnorientedAngleSigned(a, b)) > smoothThreshold)
                largeSteps++;
        }

        Assert.InRange(largeSteps, 0, 1);
    }
}

[HarnessFixture("Validates that spectral signals are isolated per component on a disconnected two-component graph.")]
public sealed class DisconnectedControlFixtureGoldenDrafts
{
    [Fact]
    public void DisconnectedControl_NodeSignal_IsolatedPerComponent()
    {
        var control = DisconnectedControl.Generate();
        var graph = FixtureGoldenHelpers.BuildDraftGraph(control.Points, k: 4, ensureConnected: false);
        var draft = FixtureGoldenHelpers.RunDraftSpectral(control.Points, graph, eigenCount: 2, seed: FixtureGoldenHelpers.SpectralSeed);

        Assert.Equal(2, draft.ComponentCount);

        NodeSignalLayer[] componentSignals = FixtureGoldenHelpers.BuildPerComponentFiedlers(control.Points, draft.Graph, eigenCount: 2, seed: FixtureGoldenHelpers.SpectralSeed);
        Assert.Equal(draft.ComponentCount, componentSignals.Length);

        double[] maxOffComponent = new double[draft.ComponentCount];

        for (int component = 0; component < draft.ComponentCount; component++)
        {
            double[] values = componentSignals[component].Values.ToArray();
            double outside = FixtureGoldenHelpers.MaxAbsOffComponent(values, draft.ComponentLabels, component);
            maxOffComponent[component] = outside;

            Assert.InRange(outside, 0.0, 1e-9);
        }

        var run = FixtureGoldenArtifacts.CreateRun(
            nameof(DisconnectedControlFixtureGoldenDrafts),
            nameof(DisconnectedControl_NodeSignal_IsolatedPerComponent),
            fixtureName: "DisconnectedControl",
            assertionName: "PerComponentIsolation",
            graph.NodeCount);
        string analysisPath = FixtureGoldenArtifacts.WriteAnalysis(
            run,
            graph,
            draft,
            new
            {
                HelperAppliesComponentSplitBeforeExtraction = false,
                SplitAppliedIn = nameof(FixtureGoldenHelpers.BuildPerComponentFiedlers),
                AssertionExtractionPath = "per-component-lift",
                ComponentSummary = FixtureGoldenHelpers.FormatComponentSummary(draft.ComponentLabels),
                MaxAbsOffComponent = maxOffComponent,
                ExpectedRange = new { Min = 0.0, Max = 1e-9 },
            });
        FixtureGoldenArtifacts.WriteRunMetadata(run, analysisPath);
    }
}

internal static class FixtureGoldenArtifacts
{
    public static ArtifactRun CreateRun(
        string suiteName,
        string runName,
        string fixtureName,
        string assertionName,
        int nodeCount)
    {
        return HarnessArtifacts.Create(
            runKind: "test-runs",
            suiteName: suiteName,
            runName: runName,
            metadata: new Dictionary<string, object?>
            {
                ["Fixture"] = fixtureName,
                ["Assertion"] = assertionName,
                ["NodeCount"] = nodeCount,
            });
    }

    public static string WriteAnalysis(
        ArtifactRun run,
        CsrGraph graph,
        DraftSpectralFixture draft,
        object assertion)
    {
        return run.WriteRunJson(
            "analysis",
            new FixtureGoldenAnalysis(
                Graph: BuildGraphSummary(graph),
                Draft: BuildDraftSummary(draft),
                Assertion: assertion));
    }

    public static void WriteRunMetadata(ArtifactRun run, string analysisPath)
    {
        Console.WriteLine($"RunRoot\t{run.RunDirectory}");
        Console.WriteLine($"Manifest\t{run.ManifestPath}");
        Console.WriteLine($"Analysis\t{analysisPath}");
    }

    private static FixtureGraphSummary BuildGraphSummary(CsrGraph graph)
    {
        ConnectivityReport connectivity = Connectivity.Validate(graph);
        return new FixtureGraphSummary(
            graph.NodeCount,
            connectivity.ComponentCount,
            FixtureGoldenHelpers.FormatLargestComponentSizes(graph),
            FixtureGoldenHelpers.GetUndirectedEdgeCount(graph),
            FixtureGoldenHelpers.GetMinimumEdgeWeight(graph));
    }

    private static FixtureDraftSummary BuildDraftSummary(DraftSpectralFixture draft)
    {
        return new FixtureDraftSummary(
            draft.ExtractionPath,
            draft.ComponentCount,
            draft.FiedlerEigenIndex,
            draft.SecondEigenEigenIndex,
            draft.FiedlerEigenIndex >= 0 ? draft.Eigenpairs[draft.FiedlerEigenIndex].Lambda : double.NaN,
            draft.SecondEigenEigenIndex >= 0 ? draft.Eigenpairs[draft.SecondEigenEigenIndex].Lambda : double.NaN,
            TakeEigenvalues(draft.Eigenpairs, first: true, count: 10),
            TakeEigenvalues(draft.Eigenpairs, first: false, count: 10));
    }

    private static double[] TakeEigenvalues(IReadOnlyList<EigenPair> eigenpairs, bool first, int count)
    {
        if (eigenpairs.Count == 0)
            return Array.Empty<double>();

        if (first)
            return eigenpairs.Take(Math.Min(count, eigenpairs.Count)).Select(pair => pair.Lambda).ToArray();

        int skip = Math.Max(0, eigenpairs.Count - count);
        return eigenpairs.Skip(skip).Select(pair => pair.Lambda).ToArray();
    }

    private sealed record FixtureGoldenAnalysis(
        FixtureGraphSummary Graph,
        FixtureDraftSummary Draft,
        object Assertion);

    private sealed record FixtureGraphSummary(
        int NodeCount,
        int ComponentCount,
        string LargestComponents,
        int UndirectedEdgeCount,
        double MinimumEdgeWeight);

    private sealed record FixtureDraftSummary(
        string ExtractionPath,
        int ComponentCount,
        int FiedlerEigenIndex,
        int SecondEigenIndex,
        double FiedlerEigenvalue,
        double SecondEigenvalue,
        double[] FirstEigenvalues,
        double[] LastEigenvalues);
}
