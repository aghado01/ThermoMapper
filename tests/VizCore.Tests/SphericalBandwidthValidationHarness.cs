using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Distance.Geodesic;
using Graphs.Pipeline.Scalers;
using Graphs.Primitives;
using Maths.Geometry;
using Maths.LinAlg;
using Repo.TestHarness;
using Xunit;

namespace VizCore.Tests;

[HarnessFixture("Spherical bandwidth validation and intrinsic-vs-linear A/B harness")]
public sealed class SphericalBandwidthValidationHarness
{
    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Validate_SphericalBandwidth_MonotonicityAndPositiveDefiniteness_WritesHarnessArtifacts()
    {
        await Task.Run(() =>
        {
            ArtifactRun run = HarnessArtifacts.Create(
                runKind: "test-runs",
                suiteName: nameof(SphericalBandwidthValidationHarness),
                runName: nameof(Validate_SphericalBandwidth_MonotonicityAndPositiveDefiniteness_WritesHarnessArtifacts));

            // Part 1: Generate samples and check monotonicity
            List<MonotonicityReportRow> monotonicityRows = RunMonotonicitySweeps();
            
            // Part 2: Generate dense points near antipode and check positive-definiteness under GlobalSchoenberg
            List<SchoenbergPdReportRow> pdRows = RunSchoenbergPdChecks();

            string monotonicityPath = run.WriteRunJson("monotonicity-recovery", monotonicityRows);
            string pdPath = run.WriteRunJson("schoenberg-pd-checks", pdRows);
            string summaryPath = run.WriteRunText("summary", BuildSummary(monotonicityRows, pdRows));

            Console.WriteLine($"RunRoot\t{run.RunDirectory}");
            Console.WriteLine($"Manifest\t{run.ManifestPath}");
            Console.WriteLine($"MonotonicityRecovery\t{monotonicityPath}");
            Console.WriteLine($"SchoenbergPdChecks\t{pdPath}");
            Console.WriteLine($"Summary\t{summaryPath}");

            // Assertions for Monotonicity
            Assert.NotEmpty(monotonicityRows);
            foreach (var group in monotonicityRows.GroupBy(r => r.K))
            {
                var ordered = group.OrderBy(r => r.TrueSigma).ToArray();
                Assert.True(ordered.Length >= 2, "Need at least two sigma levels to assert monotonicity.");
                for (int i = 0; i < ordered.Length - 1; i++)
                {
                    Assert.True(
                        ordered[i + 1].RecoveredBandwidth > ordered[i].RecoveredBandwidth,
                        $"Expected recovered bandwidth to be monotonic with true sigma at k={group.Key}. " +
                        $"Sigma={ordered[i].TrueSigma} -> Bandwidth={ordered[i].RecoveredBandwidth:F6}; " +
                        $"Sigma={ordered[i + 1].TrueSigma} -> Bandwidth={ordered[i + 1].RecoveredBandwidth:F6}");
                }
            }

            // Assertions for Schoenberg PD Checks
            Assert.NotEmpty(pdRows);
            foreach (var row in pdRows)
            {
                Assert.True(row.IsPositiveDefinite, $"Expected Gram matrix to be positive definite for dimension {row.Dimension}, delta {row.Delta}");
                Assert.True(row.MinEigenvalue >= -1e-11, $"Minimum eigenvalue {row.MinEigenvalue:E3} below tolerance -1e-11 for dimension {row.Dimension}, delta {row.Delta}");
            }
        });
    }

    private static List<MonotonicityReportRow> RunMonotonicitySweeps()
    {
        var rows = new List<MonotonicityReportRow>();
        int dimension = 3; // Ambient coordinate count n=3, S^2 sphere (m=2)
        int pointCount = 2048;
        int seed = 42;

        double[] sigmas = { 0.20, 0.35 };
        int[] kValues = { 16, 32, 48 };

        foreach (int k in kValues)
        {
            foreach (double sigma in sigmas)
            {
                double[][] points = GenerateTangentGaussianSpherePoints(dimension, sigma, pointCount, seed);
                var metric = new SphericalGeodesicMetric();

                // Set up the spherical intrinsic coupling configuration
                var config = new GraphCompilerConfig
                {
                    Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = k },
                    Filter = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
                    Repair = new RepairConfig { Kind = RepairKind.NoRepair },
                    Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
                    Projection = new CouplingProjection
                    {
                        Kernel = new Gaussian(0.0), // Auto-estimate bandwidth
                        LmpRescale = false,
                        Fidelity = CouplingFidelity.Intrinsic,
                        SphericalMode = SphericalIntrinsicMode.LocalParametrix
                    }
                };

                GraphBuildResult build = SpcGraphBuilder.BuildResult(points, config, metric);
                Assert.NotNull(build.SingleBandwidth);
                double recovered = build.SingleBandwidth.Value;

                rows.Add(new MonotonicityReportRow(
                    K: k,
                    TrueSigma: sigma,
                    RecoveredBandwidth: recovered));
            }
        }

        return rows;
    }

    private static List<SchoenbergPdReportRow> RunSchoenbergPdChecks()
    {
        var rows = new List<SchoenbergPdReportRow>();
        int[] dimensions = { 2, 3, 5 }; // S^1, S^2, S^4
        int countPerPole = 25; // 50 points total
        double sigma = 0.15; // points tightly clustered near poles (antipodal pairs)
        int seed = 12345;
        double[] deltas = { 0.2, 0.5, 1.0 };

        foreach (int dimension in dimensions)
        {
            double[][] points = GenerateAntipodalDensePoints(dimension, sigma, countPerPole, seed);
            int m = dimension - 1; // intrinsic dimension

            foreach (double delta in deltas)
            {
                // Instantiate the scaler to get access to EvaluateGlobalSchoenbergNormalized
                var scaler = new GlobalBandwidthScaler(
                    KernelType.Gaussian,
                    bandwidth: delta,
                    strategy: BandwidthStrategy.QuantileNormalized,
                    geometry: SpaceGeometry.Spherical,
                    fidelity: CouplingFidelity.Intrinsic,
                    ambientDimension: dimension,
                    sphericalMode: SphericalIntrinsicMode.GlobalSchoenberg);

                var method = typeof(GlobalBandwidthScaler).GetMethod(
                    "EvaluateGlobalSchoenbergNormalized",
                    System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.NonPublic);
                
                Assert.NotNull(method);

                int n = points.Length;
                var gram = new double[n, n];
                var manifold = new SphericalManifold(dimension);

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        double dist = manifold.Distance(points[i], points[j]);
                        gram[i, j] = (double)method.Invoke(scaler, new object[] { dist, delta, m })!;
                    }
                }

                // Decompose the Gram matrix to get eigenvalues
                EigenResult result = Eigen.DecomposeSymmetric(gram);
                double minEigen = result.Eigenvalues.Min();
                bool isPd = minEigen >= -1e-11;

                rows.Add(new SchoenbergPdReportRow(
                    Dimension: dimension,
                    Delta: delta,
                    MinEigenvalue: minEigen,
                    IsPositiveDefinite: isPd));
            }
        }

        return rows;
    }

    private static double[][] GenerateTangentGaussianSpherePoints(int dimension, double sigma, int count, int seed)
    {
        var rng = new Random(seed);
        var manifold = new SphericalManifold(dimension);
        double[] p = new double[dimension];
        p[0] = 1.0;

        double[][] points = new double[count][];
        for (int i = 0; i < count; i++)
        {
            double[] v = new double[dimension];
            for (int j = 1; j < dimension; j += 2)
            {
                double u1 = rng.NextDouble();
                double u2 = rng.NextDouble();
                double z0 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
                v[j] = z0 * sigma;
                if (j + 1 < dimension)
                {
                    double z1 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
                    v[j + 1] = z1 * sigma;
                }
            }
            
            points[i] = new double[dimension];
            manifold.ExpMap(p, v, points[i]);
        }
        return points;
    }

    private static double[][] GenerateAntipodalDensePoints(int dimension, double sigma, int countPerPole, int seed)
    {
        var rng = new Random(seed);
        var manifold = new SphericalManifold(dimension);
        double[] p1 = new double[dimension];
        p1[0] = 1.0;
        double[] p2 = new double[dimension];
        p2[0] = -1.0;

        double[][] points = new double[countPerPole * 2][];
        
        for (int i = 0; i < countPerPole; i++)
        {
            double[] v = new double[dimension];
            for (int j = 1; j < dimension; j += 2)
            {
                double u1 = rng.NextDouble();
                double u2 = rng.NextDouble();
                double z0 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
                v[j] = z0 * sigma;
                if (j + 1 < dimension)
                {
                    double z1 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
                    v[j + 1] = z1 * sigma;
                }
            }
            points[i] = new double[dimension];
            manifold.ExpMap(p1, v, points[i]);
        }

        for (int i = 0; i < countPerPole; i++)
        {
            double[] v = new double[dimension];
            for (int j = 1; j < dimension; j += 2)
            {
                double u1 = rng.NextDouble();
                double u2 = rng.NextDouble();
                double z0 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
                v[j] = z0 * sigma;
                if (j + 1 < dimension)
                {
                    double z1 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
                    v[j + 1] = z1 * sigma;
                }
            }
            points[countPerPole + i] = new double[dimension];
            manifold.ExpMap(p2, v, points[countPerPole + i]);
        }

        return points;
    }

    private static string BuildSummary(
        IReadOnlyList<MonotonicityReportRow> monotonicityRows,
        IReadOnlyList<SchoenbergPdReportRow> pdRows)
    {
        var builder = new StringBuilder();
        builder.AppendLine("Spherical Bandwidth Validation Harness Summary");
        builder.AppendLine();
        builder.AppendLine("Monotonicity recovery sweeps");
        builder.AppendLine("K\tTrueSigma\tRecoveredBandwidth");
        foreach (var row in monotonicityRows)
        {
            builder.AppendLine($"{row.K}\t{row.TrueSigma:F2}\t{row.RecoveredBandwidth:F6}");
        }

        builder.AppendLine();
        builder.AppendLine("Schoenberg positive-definiteness checks");
        builder.AppendLine("Dimension\tDelta\tMinEigenvalue\tIsPositiveDefinite");
        foreach (var row in pdRows)
        {
            builder.AppendLine($"{row.Dimension}\t{row.Delta:F1}\t{row.MinEigenvalue:E3}\t{row.IsPositiveDefinite}");
        }

        return builder.ToString();
    }

    private sealed record MonotonicityReportRow(
        int K,
        double TrueSigma,
        double RecoveredBandwidth);

    private sealed record SchoenbergPdReportRow(
        int Dimension,
        double Delta,
        double MinEigenvalue,
        bool IsPositiveDefinite);
}
