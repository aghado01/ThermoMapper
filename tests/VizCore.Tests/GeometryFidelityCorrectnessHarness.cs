using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance.Geodesic;
using Graphs.Primitives;
using Maths.Geometry;
using Repo.TestHarness;
using Xunit;

namespace VizCore.Tests;

[HarnessFixture("Geometry fidelity correctness harness for manifold invariants, Tangent-vs-GeodesicLinear parity, and shared-base distortion diagnostics")]
public sealed class GeometryFidelityCorrectnessHarness
{
    private const double InvariantTolerance = 1e-9;
    private const double RoundTripTolerance = 1e-9;
    private const double WeightTolerance = 1e-12;
    private const double AmbientContrastThreshold = 1e-3;
    private const double SmallRingDistortionCeiling = 0.05;

    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Validate_GeometryFidelityCorrectness_WritesHarnessArtifacts()
    {
        await Task.Run(() =>
        {
            ArtifactRun run = HarnessArtifacts.Create(
                runKind: "test-runs",
                suiteName: nameof(GeometryFidelityCorrectnessHarness),
                runName: nameof(Validate_GeometryFidelityCorrectness_WritesHarnessArtifacts));

            IReadOnlyList<ManifoldInvariantRow> invariantRows = BuildInvariantRows();
            IReadOnlyList<TangentDistortionReport> distortionReports = BuildDistortionReports();

            string invariantPath = run.WriteRunJson("manifold-invariants", invariantRows);
            string distortionPath = run.WriteRunJson("shared-base-tangent-distortion", distortionReports);
            string summaryPath = run.WriteRunText("summary", BuildSummary(invariantRows, distortionReports));

            Console.WriteLine($"RunRoot\t{run.RunDirectory}");
            Console.WriteLine($"Manifest\t{run.ManifestPath}");
            Console.WriteLine($"ManifoldInvariants\t{invariantPath}");
            Console.WriteLine($"SharedBaseTangentDistortion\t{distortionPath}");
            Console.WriteLine($"Summary\t{summaryPath}");

            Assert.NotEmpty(invariantRows);
            Assert.All(
                invariantRows,
                row =>
                {
                    Assert.InRange(Math.Abs(row.Distance - row.ForwardRadius), 0.0, InvariantTolerance);
                    Assert.InRange(Math.Abs(row.Distance - row.BackwardRadius), 0.0, InvariantTolerance);
                    Assert.InRange(Math.Abs(row.ForwardRadius - row.BackwardRadius), 0.0, InvariantTolerance);
                    Assert.InRange(row.RoundTripMaxError, 0.0, RoundTripTolerance);
                    Assert.True(
                        Math.Abs(row.Distance - row.AmbientForwardRadius) > AmbientContrastThreshold,
                        $"Expected the ambient tangent norm to differ from the geodesic distance for d={row.Dimension}, fixture={row.FixtureIndex}. distance={row.Distance}, ambient={row.AmbientForwardRadius}");
                });

            Assert.Equal(3, distortionReports.Count);
            TangentDistortionReport[] ordered = distortionReports.OrderBy(report => report.RingRadius).ToArray();
            Assert.InRange(ordered[0].MeanRelativeError, 0.0, SmallRingDistortionCeiling);
            Assert.True(
                ordered[0].MeanRelativeError < ordered[1].MeanRelativeError && ordered[1].MeanRelativeError < ordered[2].MeanRelativeError,
                $"Expected shared-base tangent distortion to grow with ring radius. small={ordered[0].MeanRelativeError:G6}, medium={ordered[1].MeanRelativeError:G6}, large={ordered[2].MeanRelativeError:G6}");
        });
    }

    private static IReadOnlyList<ManifoldInvariantRow> BuildInvariantRows()
    {
        int[] dimensions = { 2, 3, 8 };
        var rows = new List<ManifoldInvariantRow>(dimensions.Length * 2);

        foreach (int dimension in dimensions)
        {
            for (int fixtureIndex = 0; fixtureIndex < 2; fixtureIndex++)
                rows.Add(BuildInvariantRow(dimension, fixtureIndex));
        }

        return rows;
    }

    private static ManifoldInvariantRow BuildInvariantRow(int dimension, int fixtureIndex)
    {
        var manifold = new PoincareBallManifold(dimension);
        (double[] basePoint, double[] targetPoint) = CreateInvariantFixture(dimension, fixtureIndex);
        var forward = new double[dimension];
        var backward = new double[dimension];
        var recovered = new double[dimension];

        manifold.LogMap(basePoint, targetPoint, forward);
        manifold.LogMap(targetPoint, basePoint, backward);
        manifold.ExpMap(basePoint, forward, recovered);

        double distance = manifold.Distance(basePoint, targetPoint);
        double roundTripMaxError = MaxAbsDifference(targetPoint, recovered);

        return new ManifoldInvariantRow(
            Dimension: dimension,
            FixtureIndex: fixtureIndex,
            Distance: distance,
            ForwardRadius: manifold.Norm(basePoint, forward),
            BackwardRadius: manifold.Norm(targetPoint, backward),
            AmbientForwardRadius: EuclideanNorm(forward),
            AmbientBackwardRadius: EuclideanNorm(backward),
            RoundTripMaxError: roundTripMaxError);
    }


    private static IReadOnlyList<TangentDistortionReport> BuildDistortionReports()
    {
        double[] ringRadii = { 0.05, 0.20, 0.50 };
        var reports = new List<TangentDistortionReport>(ringRadii.Length);

        foreach (double ringRadius in ringRadii)
            reports.Add(BuildDistortionReport(ringRadius));

        return reports;
    }

    private static TangentDistortionReport BuildDistortionReport(double ringRadius)
    {
        var manifold = new PoincareBallManifold(2);
        double[] basePoint = { 0.0, 0.0 };
        double[][] tangentVectors =
        {
            new[] { ringRadius, 0.0 },
            new[] { 0.0, ringRadius },
            new[] { -ringRadius, 0.0 },
            new[] { 0.0, -ringRadius },
        };

        var embeddedPoints = new double[tangentVectors.Length][];
        for (int index = 0; index < tangentVectors.Length; index++)
        {
            embeddedPoints[index] = new double[2];
            manifold.ExpMap(basePoint, tangentVectors[index], embeddedPoints[index]);
        }

        var rows = new List<TangentDistortionRow>();
        for (int left = 0; left < tangentVectors.Length; left++)
        {
            for (int right = left + 1; right < tangentVectors.Length; right++)
            {
                double tangentDistance = TangentDifferenceNorm(manifold, basePoint, tangentVectors[left], tangentVectors[right]);
                double geodesicDistance = manifold.Distance(embeddedPoints[left], embeddedPoints[right]);
                double absError = Math.Abs(tangentDistance - geodesicDistance);
                double relativeError = absError / Math.Max(1e-12, geodesicDistance);

                rows.Add(new TangentDistortionRow(
                    RingRadius: ringRadius,
                    Left: left,
                    Right: right,
                    TangentDistance: tangentDistance,
                    GeodesicDistance: geodesicDistance,
                    AbsError: absError,
                    RelativeError: relativeError));
            }
        }

        return new TangentDistortionReport(
            RingRadius: ringRadius,
            MeanRelativeError: rows.Average(row => row.RelativeError),
            MaxRelativeError: rows.Max(row => row.RelativeError),
            Rows: rows);
    }

    private static string BuildSummary(
        IReadOnlyList<ManifoldInvariantRow> invariantRows,
        IReadOnlyList<TangentDistortionReport> distortionReports)
    {
        var builder = new StringBuilder();
        builder.AppendLine("Geometry Fidelity Correctness Harness");
        builder.AppendLine();
        builder.AppendLine("Manifold invariants");
        builder.AppendLine("Dimension\tFixture\tDistance\tForwardRadius\tBackwardRadius\tAmbientForward\tRoundTripMaxError");
        foreach (ManifoldInvariantRow row in invariantRows.OrderBy(row => row.Dimension).ThenBy(row => row.FixtureIndex))
        {
            builder.AppendLine($"{row.Dimension}\t{row.FixtureIndex}\t{row.Distance:F9}\t{row.ForwardRadius:F9}\t{row.BackwardRadius:F9}\t{row.AmbientForwardRadius:F9}\t{row.RoundTripMaxError:E3}");
        }

        builder.AppendLine("Shared-base tangent distortion");
        builder.AppendLine("RingRadius\tMeanRelativeError\tMaxRelativeError");
        foreach (TangentDistortionReport report in distortionReports.OrderBy(report => report.RingRadius))
        {
            builder.AppendLine($"{report.RingRadius:F3}\t{report.MeanRelativeError:F9}\t{report.MaxRelativeError:F9}");
        }

        return builder.ToString();
    }

    private static (double[] BasePoint, double[] TargetPoint) CreateInvariantFixture(int dimension, int fixtureIndex)
    {
        return fixtureIndex switch
        {
            0 =>
            (
                CreatePoint(dimension, 0.22, -0.18, 0.11, 0.04),
                CreatePoint(dimension, -0.05, 0.08, -0.02, 0.03)
            ),
            1 =>
            (
                CreatePoint(dimension, 0.16, -0.12, 0.09, 0.02),
                CreatePoint(dimension, -0.09, 0.14, -0.04, -0.01)
            ),
            _ => throw new ArgumentOutOfRangeException(nameof(fixtureIndex), fixtureIndex, "Unknown invariant fixture index."),
        };
    }

    private static double[] CreatePoint(int dimension, params double[] seed)
    {
        var point = new double[dimension];
        for (int i = 0; i < dimension; i++)
            point[i] = i < seed.Length ? seed[i] : 0.01 * (i + 1);

        return point;
    }

    private static double TangentDifferenceNorm(
        PoincareBallManifold manifold,
        double[] basePoint,
        double[] left,
        double[] right)
    {
        Span<double> delta = stackalloc double[left.Length];
        for (int i = 0; i < left.Length; i++)
            delta[i] = left[i] - right[i];

        return manifold.Norm(basePoint, delta);
    }

    private static double[][] CreateLongRangePoincareFixture(int ambientDimension)
    {
        if (ambientDimension < 1)
            throw new ArgumentOutOfRangeException(nameof(ambientDimension), ambientDimension, "Ambient dimension must be positive.");

        double[] axisCoordinates = { -0.97, -0.90, -0.40, -0.05, 0.00, 0.05, 0.40, 0.90, 0.97 };
        var features = new double[axisCoordinates.Length][];

        for (int index = 0; index < axisCoordinates.Length; index++)
        {
            var point = new double[ambientDimension];
            point[0] = axisCoordinates[index];
            features[index] = point;
        }

        return features;
    }

    private static GraphCompilerConfig CreateHarnessConfig(CouplingFidelity fidelity) => new()
    {
        Topology = new TopologyConfig { Kind = TopologyKind.EpsilonBall, Epsilon = 10.0 },
        Filter = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
        Repair = new RepairConfig { Kind = RepairKind.NoRepair },
        Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
        Projection = new CouplingProjection
        {
            Kernel = new Gaussian(0.9),
            LmpRescale = false,
            Fidelity = fidelity,
        },
    };

    private static Dictionary<long, double> BuildUndirectedWeightMap(CsrGraph graph)
    {
        var map = new Dictionary<long, double>();

        for (int source = 0; source < graph.NodeCount; source++)
        {
            int rowStart = graph.RowPointers[source];
            int rowEnd = graph.RowPointers[source + 1];
            for (int edge = rowStart; edge < rowEnd; edge++)
            {
                int target = graph.Targets[edge];
                if (target <= source)
                    continue;

                map[EdgeKey(source, target)] = graph.Weights[edge];
            }
        }

        return map;
    }

    private static long EdgeKey(int left, int right)
    {
        int source = Math.Min(left, right);
        int target = Math.Max(left, right);
        return (((long)source) << 32) | (uint)target;
    }

    private static double EuclideanNorm(double[] vector)
    {
        double sumSq = 0.0;
        for (int i = 0; i < vector.Length; i++)
            sumSq += vector[i] * vector[i];

        return Math.Sqrt(sumSq);
    }

    private static double MaxAbsDifference(double[] left, double[] right)
    {
        double max = 0.0;
        for (int i = 0; i < left.Length; i++)
            max = Math.Max(max, Math.Abs(left[i] - right[i]));

        return max;
    }

    private sealed record ManifoldInvariantRow(
        int Dimension,
        int FixtureIndex,
        double Distance,
        double ForwardRadius,
        double BackwardRadius,
        double AmbientForwardRadius,
        double AmbientBackwardRadius,
        double RoundTripMaxError);

    private sealed record TangentDistortionReport(
        double RingRadius,
        double MeanRelativeError,
        double MaxRelativeError,
        IReadOnlyList<TangentDistortionRow> Rows);

    private sealed record TangentDistortionRow(
        double RingRadius,
        int Left,
        int Right,
        double TangentDistance,
        double GeodesicDistance,
        double AbsError,
        double RelativeError);
}
