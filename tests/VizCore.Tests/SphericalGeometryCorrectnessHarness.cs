using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Maths.Geometry;
using Repo.TestHarness;
using Xunit;

namespace VizCore.Tests;

[HarnessFixture("Spherical geometry correctness harness for manifold invariants, radial isometry, log/exp round-trip, symmetry, contrast, and antipodal zero guards")]
public sealed class SphericalGeometryCorrectnessHarness
{
    private const double InvariantTolerance = 1e-9;
    private const double RoundTripTolerance = 1e-9;
    private const double AmbientContrastThreshold = 1e-3;

    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Validate_SphericalGeometry_WritesHarnessArtifacts()
    {
        await Task.Run(() =>
        {
            ArtifactRun run = HarnessArtifacts.Create(
                runKind: "test-runs",
                suiteName: nameof(SphericalGeometryCorrectnessHarness),
                runName: nameof(Validate_SphericalGeometry_WritesHarnessArtifacts));

            IReadOnlyList<SphericalInvariantRow> invariantRows = BuildInvariantRows();
            IReadOnlyList<SphericalAntipodalRow> antipodalRows = BuildAntipodalRows();

            string invariantPath = run.WriteRunJson("manifold-invariants", invariantRows);
            string antipodalPath = run.WriteRunJson("antipodal-limits", antipodalRows);
            string summaryPath = run.WriteRunText("summary", BuildSummary(invariantRows, antipodalRows));

            Console.WriteLine($"RunRoot\t{run.RunDirectory}");
            Console.WriteLine($"Manifest\t{run.ManifestPath}");
            Console.WriteLine($"ManifoldInvariants\t{invariantPath}");
            Console.WriteLine($"AntipodalLimits\t{antipodalPath}");
            Console.WriteLine($"Summary\t{summaryPath}");

            Assert.NotEmpty(invariantRows);
            
            // I1 (radial isometry), I2 (round-trip), I3 (symmetry), I4 (ambient contrast) assertions
            foreach (var row in invariantRows)
            {
                // I1: ||log_p(q)||_2 = d(p,q)
                Assert.InRange(Math.Abs(row.Distance - row.ForwardRadius), 0.0, InvariantTolerance);
                Assert.InRange(Math.Abs(row.Distance - row.BackwardRadius), 0.0, InvariantTolerance);
                Assert.InRange(Math.Abs(row.ForwardRadius - row.BackwardRadius), 0.0, InvariantTolerance);

                // I2: exp_p(log_p(q)) = q
                Assert.InRange(row.RoundTripMaxError, 0.0, RoundTripTolerance);

                // I3: Directed symmetry (d(p,q) = d(q,p))
                Assert.InRange(Math.Abs(row.Distance - row.DistanceBackward), 0.0, InvariantTolerance);

                // I4: Ambient Contrast
                // Geodesic distance must differ from Euclidean chord length by > 10^-3, unless points are extremely close.
                // We'll assert that the geodesic distance is larger than the Euclidean distance by at least AmbientContrastThreshold.
                Assert.True(
                    row.Distance > row.EuclideanChordLength + AmbientContrastThreshold,
                    $"Expected geodesic distance ({row.Distance:F6}) to be strictly greater than Euclidean chord length ({row.EuclideanChordLength:F6}) by at least {AmbientContrastThreshold} due to curvature.");
            }

            // I5: Antipodal limit zero-guard
            Assert.NotEmpty(antipodalRows);
            foreach (var row in antipodalRows)
            {
                Assert.True(row.IsAntipodalGuardTriggered, $"Expected antipodal log map guard to fire for dimension {row.Dimension}");
                Assert.All(row.LogMapResult, val => Assert.Equal(0.0, val));
            }
        });
    }

    private static IReadOnlyList<SphericalInvariantRow> BuildInvariantRows()
    {
        // Ambient dimensions 2, 3, 8 (corresponds to intrinsic S^1, S^2, S^7)
        int[] dimensions = { 2, 3, 8 };
        var rows = new List<SphericalInvariantRow>();

        foreach (int dimension in dimensions)
        {
            // Sample angles to construct points. Let's do a few different angles: pi/6, pi/3, pi/2, 2*pi/3.
            // None of these are close to 0 or pi to avoid boundary issues.
            double[] angles = { Math.PI / 6.0, Math.PI / 3.0, Math.PI / 2.0, 2.0 * Math.PI / 3.0 };
            for (int i = 0; i < angles.Length; i++)
            {
                rows.Add(BuildInvariantRow(dimension, angles[i], i));
            }
        }

        return rows;
    }

    private static SphericalInvariantRow BuildInvariantRow(int dimension, double angle, int fixtureIndex)
    {
        var manifold = new SphericalManifold(dimension);
        
        // Let's create p and q deterministically using the angle
        double[] p = new double[dimension];
        double[] q = new double[dimension];
        
        p[0] = 1.0; // Rest are 0.0
        q[0] = Math.Cos(angle);
        q[1] = Math.Sin(angle); // Rest are 0.0

        var forward = new double[dimension];
        var backward = new double[dimension];
        var recovered = new double[dimension];

        manifold.LogMap(p, q, forward);
        manifold.LogMap(q, p, backward);
        manifold.ExpMap(p, forward, recovered);

        double distance = manifold.Distance(p, q);
        double distanceBackward = manifold.Distance(q, p);
        
        double roundTripMaxError = MaxAbsDifference(q, recovered);
        double euclideanChordLength = EuclideanDistance(p, q);

        return new SphericalInvariantRow(
            Dimension: dimension,
            Angle: angle,
            FixtureIndex: fixtureIndex,
            Distance: distance,
            DistanceBackward: distanceBackward,
            ForwardRadius: manifold.Norm(p, forward),
            BackwardRadius: manifold.Norm(q, backward),
            RoundTripMaxError: roundTripMaxError,
            EuclideanChordLength: euclideanChordLength);
    }

    private static IReadOnlyList<SphericalAntipodalRow> BuildAntipodalRows()
    {
        int[] dimensions = { 2, 3, 8 };
        var rows = new List<SphericalAntipodalRow>();

        foreach (int dimension in dimensions)
        {
            var manifold = new SphericalManifold(dimension);
            
            // Construct p and its exact antipode q = -p
            double[] p = new double[dimension];
            double[] q = new double[dimension];
            
            p[0] = 1.0;
            q[0] = -1.0; // Rest are 0.0

            var logMapResult = new double[dimension];
            manifold.LogMap(p, q, logMapResult);

            // Check if all are exactly 0.0
            bool isZero = logMapResult.All(val => val == 0.0);

            rows.Add(new SphericalAntipodalRow(
                Dimension: dimension,
                IsAntipodalGuardTriggered: isZero,
                LogMapResult: logMapResult));
        }

        return rows;
    }

    private static double MaxAbsDifference(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        double max = 0.0;
        for (int i = 0; i < a.Length; i++)
        {
            double diff = Math.Abs(a[i] - b[i]);
            if (diff > max) max = diff;
        }
        return max;
    }

    private static double EuclideanDistance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        double sum = 0.0;
        for (int i = 0; i < a.Length; i++)
        {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return Math.Sqrt(sum);
    }

    private static string BuildSummary(
        IReadOnlyList<SphericalInvariantRow> invariantRows,
        IReadOnlyList<SphericalAntipodalRow> antipodalRows)
    {
        var builder = new StringBuilder();
        builder.AppendLine("Spherical Geometry Correctness Harness Summary");
        builder.AppendLine();
        builder.AppendLine("Manifold invariants");
        builder.AppendLine("Dimension\tAngle\tFixture\tDistance\tDistanceBackward\tForwardRadius\tBackwardRadius\tRoundTripMaxError\tEuclideanChordLength");
        foreach (var row in invariantRows)
        {
            builder.AppendLine($"{row.Dimension}\t{row.Angle:F4}\t{row.FixtureIndex}\t{row.Distance:F9}\t{row.DistanceBackward:F9}\t{row.ForwardRadius:F9}\t{row.BackwardRadius:F9}\t{row.RoundTripMaxError:E3}\t{row.EuclideanChordLength:F9}");
        }

        builder.AppendLine();
        builder.AppendLine("Antipodal Limit Zero-Guards");
        builder.AppendLine("Dimension\tIsAntipodalGuardTriggered\tResultLogMap");
        foreach (var row in antipodalRows)
        {
            builder.AppendLine($"{row.Dimension}\t{row.IsAntipodalGuardTriggered}\t[{string.Join(",", row.LogMapResult.Select(v => v.ToString("F3")))}]");
        }

        return builder.ToString();
    }

    private sealed record SphericalInvariantRow(
        int Dimension,
        double Angle,
        int FixtureIndex,
        double Distance,
        double DistanceBackward,
        double ForwardRadius,
        double BackwardRadius,
        double RoundTripMaxError,
        double EuclideanChordLength);

    private sealed record SphericalAntipodalRow(
        int Dimension,
        bool IsAntipodalGuardTriggered,
        double[] LogMapResult);
}
