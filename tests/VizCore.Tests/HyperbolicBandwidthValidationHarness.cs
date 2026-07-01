using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance.Geodesic;
using Graphs.Primitives;
using Repo.TestHarness;
using Xunit;

namespace VizCore.Tests;

[HarnessFixture("Hyperbolic factor validation and intrinsic-vs-linear A/B harness")]
public sealed class HyperbolicBandwidthValidationHarness
{
    private const int SamplesPerCombo = 2048;
    private const int CandidatePoolSize = 32;
    private const int NeighborRank = 1;
    private const string ValidationModelsEnvVar = "PWSHSPC_HYPERBOLIC_FACTOR_MODELS";
    private const int LaplacianCalibrationDimension = 3;
    private const double LaplacianCalibrationMinVolumePressure = 0.10;
    private const double LaplacianCalibrationMaxVolumePressure = 0.40;
    private const double LaplacianCalibrationTolerance = 0.20;
    private const double IntrinsicCalibrationTolerance = 0.20;
    private const double KeepThreshold = 1e-10;

    [Fact(Timeout = TestTimeouts.HeavyFixtureMs)]
    public async Task Validate_HyperbolicFactors_AndIntrinsicPayoff_WritesHarnessArtifacts()
    {
        await Task.Run(() =>
        {
            ArtifactRun run = HarnessArtifacts.Create(
                runKind: "test-runs",
                suiteName: nameof(HyperbolicBandwidthValidationHarness),
                runName: nameof(Validate_HyperbolicFactors_AndIntrinsicPayoff_WritesHarnessArtifacts));

            List<FactorValidationRow> factorRows = BuildFactorValidationRows();
            IReadOnlyList<IntrinsicSelfConsistencyReport> intrinsicSelfConsistencyReports = BuildIntrinsicSelfConsistencyReports();
            IReadOnlyList<IntrinsicGroundTruthReport> intrinsicGroundTruthReports = BuildIntrinsicGroundTruthReports();
            IReadOnlyList<IntrinsicAbReport> intrinsicReports = BuildIntrinsicAbReports();

            string factorPath = run.WriteRunJson("factor-validation", factorRows);
            string intrinsicSelfConsistencyPath = run.WriteRunJson("intrinsic-self-consistency", intrinsicSelfConsistencyReports);
            string intrinsicGroundTruthPath = run.WriteRunJson("intrinsic-ground-truth", intrinsicGroundTruthReports);
            string abPath = run.WriteRunJson("intrinsic-ab", intrinsicReports);
            string summaryPath = run.WriteRunText("summary", BuildSummary(factorRows, intrinsicSelfConsistencyReports, intrinsicGroundTruthReports, intrinsicReports));

            Console.WriteLine($"RunRoot\t{run.RunDirectory}");
            Console.WriteLine($"Manifest\t{run.ManifestPath}");
            Console.WriteLine($"FactorValidation\t{factorPath}");
            Console.WriteLine($"IntrinsicSelfConsistency\t{intrinsicSelfConsistencyPath}");
            Console.WriteLine($"IntrinsicGroundTruth\t{intrinsicGroundTruthPath}");
            Console.WriteLine($"IntrinsicAB\t{abPath}");
            Console.WriteLine($"Summary\t{summaryPath}");

            FactorValidationModel[] models = ResolveValidationModels();
            Assert.NotEmpty(models);

            foreach (FactorValidationModel model in models)
            {
                FactorValidationRow[] modelRows = factorRows
                    .Where(row => row.Model == model)
                    .ToArray();
                Assert.NotEmpty(modelRows);

                FactorValidationRow[] gaussianRows = modelRows
                    .Where(row => row.Family == KernelFamily.Gaussian)
                    .ToArray();
                Assert.NotEmpty(gaussianRows);
                Assert.All(
                    gaussianRows,
                    row => Assert.Equal(row.RecoveredEuclidean, row.RecoveredConfigured, 12));

                FactorValidationRow[] laplacianRows = modelRows
                    .Where(row => row.Family == KernelFamily.Laplacian)
                    .ToArray();
                Assert.NotEmpty(laplacianRows);

                FactorValidationRow[] laplacianCalibrationRows = laplacianRows
                    .Where(IsLaplacianCalibrationRow)
                    .ToArray();
                Assert.NotEmpty(laplacianCalibrationRows);
                double meanCalibrationFactor = laplacianCalibrationRows.Average(row => row.ImpliedFactor);
                Assert.InRange(
                    meanCalibrationFactor,
                    BandwidthEstimation.LaplacianHyperbolicFactor - LaplacianCalibrationTolerance,
                    BandwidthEstimation.LaplacianHyperbolicFactor + LaplacianCalibrationTolerance);
                Assert.True(
                    laplacianCalibrationRows.Average(row => row.RelativeErrorConfigured)
                    < laplacianCalibrationRows.Average(row => row.RelativeErrorEuclidean),
                    $"Expected the configured Laplacian hyperbolic factor to outperform the Euclidean constant in the H3 bulk calibration window for {model}.");

                FactorValidationRow[] laplacianBulkRows = laplacianRows
                    .Where(row => row.VolumePressure >= 0.25 && row.VolumePressure <= 0.80)
                    .ToArray();
                Assert.NotEmpty(laplacianBulkRows);
                Assert.True(
                    laplacianBulkRows.Average(row => row.RecoveredConfigured / row.RecoveredEuclidean) > 1.02,
                    $"Expected the configured Laplacian hyperbolic factor to inflate recovered scales above the Euclidean constant in the bulk regime for {model}.");

                FactorValidationRow[] laplacianSmallRows = laplacianRows
                    .Where(row => row.VolumePressure <= 0.10)
                    .ToArray();
                Assert.NotEmpty(laplacianSmallRows);

                FactorValidationRow[] laplacianDegradedRows = laplacianRows
                    .Where(row => row.VolumePressure >= 1.05)
                    .ToArray();
                Assert.NotEmpty(laplacianDegradedRows);
                Assert.True(
                    laplacianSmallRows.Average(row => row.RelativeErrorConfigured)
                    < laplacianDegradedRows.Average(row => row.RelativeErrorConfigured),
                    $"Expected the small-scale Laplacian probe to stay better behaved than the degraded large-scale regime for {model}.");

                Assert.True(
                    laplacianDegradedRows.Average(row => row.RelativeErrorConfigured)
                    > laplacianBulkRows.Average(row => row.RelativeErrorConfigured),
                    $"Expected the Laplacian hyperbolic estimate to degrade once scale exceeds the bulk-valid regime for {model}.");
            }

            Assert.NotEmpty(intrinsicSelfConsistencyReports);
            Assert.All(
                intrinsicSelfConsistencyReports,
                report => Assert.InRange(report.RelativeErrorRecovered, 0.0, IntrinsicCalibrationTolerance));

            Assert.NotEmpty(intrinsicGroundTruthReports);
            Assert.All(intrinsicGroundTruthReports, report =>
            {
                Assert.True(report.RecoveredBandwidth > 0.0, "Expected a positive recovered bandwidth.");
            });

            foreach (var group in intrinsicGroundTruthReports.GroupBy(report => (report.AmbientDimension, report.ActualK)))
            {
                IntrinsicGroundTruthReport lowSigma = group.Single(report => Math.Abs(report.TrueSigma - 0.20) < 1e-12);
                IntrinsicGroundTruthReport highSigma = group.Single(report => Math.Abs(report.TrueSigma - 0.35) < 1e-12);
                Assert.True(
                    highSigma.RecoveredBandwidth > lowSigma.RecoveredBandwidth,
                    $"Expected the larger true sigma to recover the larger intrinsic bandwidth at d={group.Key.AmbientDimension}, k={group.Key.ActualK}. low={lowSigma.RecoveredBandwidth:G6}, high={highSigma.RecoveredBandwidth:G6}");
            }

            IntrinsicAbReport report2 = intrinsicReports.Single(report => report.AmbientDimension == 2);
            IntrinsicAbReport report3 = intrinsicReports.Single(report => report.AmbientDimension == 3);
            IntrinsicAbReport report8 = intrinsicReports.Single(report => report.AmbientDimension == 8);

            foreach (IntrinsicAbReport intrinsicReport in intrinsicReports)
            {
                Assert.NotNull(intrinsicReport.Bandwidth);
                Assert.True(
                    intrinsicReport.FarthestLinearKeptIntrinsicWeight < intrinsicReport.FarthestLinearKeptLinearWeight,
                    $"Expected intrinsic fidelity to suppress the farthest retained linear edge at d={intrinsicReport.AmbientDimension}. linear={intrinsicReport.FarthestLinearKeptLinearWeight:G6}, intrinsic={intrinsicReport.FarthestLinearKeptIntrinsicWeight:G6}");
                Assert.True(
                    intrinsicReport.DroppedByIntrinsicCount > 0,
                    $"Expected at least one long hyperbolic edge to fall below the keep threshold under Intrinsic fidelity at d={intrinsicReport.AmbientDimension}.");
            }

            Assert.True(
                report3.ShortEdgeDifference < 2e-3,
                $"Expected short-range intrinsic suppression to stay small in the d=3 baseline path. diff={report3.ShortEdgeDifference:G6}");

            long[] commonRetainedKeys = report2.Rows
                .Where(row => row.LinearKept && row.IntrinsicKept)
                .Select(row => EdgeKey(row.Left, row.Right))
                .Intersect(report3.Rows.Where(row => row.LinearKept && row.IntrinsicKept).Select(row => EdgeKey(row.Left, row.Right)))
                .Intersect(report8.Rows.Where(row => row.LinearKept && row.IntrinsicKept).Select(row => EdgeKey(row.Left, row.Right)))
                .ToArray();
            Assert.NotEmpty(commonRetainedKeys);

            long comparisonKey = report2.Rows
                .Where(row => commonRetainedKeys.Contains(EdgeKey(row.Left, row.Right)))
                .OrderBy(row => row.Distance)
                .Last()
                is IntrinsicEdgeRow comparisonRow
                    ? EdgeKey(comparisonRow.Left, comparisonRow.Right)
                    : throw new InvalidOperationException("Expected a retained comparison edge for the intrinsic dimension harness.");

            IntrinsicEdgeRow comparisonRow2 = report2.Rows.Single(row => EdgeKey(row.Left, row.Right) == comparisonKey);
            IntrinsicEdgeRow comparisonRow3 = report3.Rows.Single(row => EdgeKey(row.Left, row.Right) == comparisonKey);
            IntrinsicEdgeRow comparisonRow8 = report8.Rows.Single(row => EdgeKey(row.Left, row.Right) == comparisonKey);

            double suppressionRatio2 = comparisonRow2.IntrinsicWeight / comparisonRow2.LinearWeight;
            double suppressionRatio3 = comparisonRow3.IntrinsicWeight / comparisonRow3.LinearWeight;
            double suppressionRatio8 = comparisonRow8.IntrinsicWeight / comparisonRow8.LinearWeight;

            Assert.True(
                suppressionRatio8 < suppressionRatio3 && suppressionRatio3 < suppressionRatio2,
                $"Expected intrinsic far-edge suppression to strengthen with ambient dimension. d2={suppressionRatio2:G6}, d3={suppressionRatio3:G6}, d8={suppressionRatio8:G6}");
            Assert.True(
                comparisonRow2.IntrinsicWeight > comparisonRow2.LegacyPocWeight,
                $"Expected d=2 intrinsic correction to suppress less than the legacy H3 POC. d2 intrinsic={comparisonRow2.IntrinsicWeight:G6}, legacy={comparisonRow2.LegacyPocWeight:G6}");
            Assert.Equal(comparisonRow3.LegacyPocWeight, comparisonRow3.IntrinsicWeight);
        });
    }

    private static List<FactorValidationRow> BuildFactorValidationRows()
    {
        var rows = new List<FactorValidationRow>();
        int[] dimensions = { 2, 3, 8 };
        FactorValidationModel[] models = ResolveValidationModels();

        double[] gaussianScales = { 0.10, 0.30, 0.75 };
        double[] laplacianScales = { 0.05, 0.10, 0.20, 0.50, 1.00 };

        foreach (FactorValidationModel model in models)
        {
            foreach (int dimension in dimensions)
            {
                for (int index = 0; index < gaussianScales.Length; index++)
                    rows.Add(BuildFactorValidationRow(model, KernelFamily.Gaussian, dimension, gaussianScales[index], seed: 1000 + ((int)model + 1) * 10000 + dimension * 100 + index));

                for (int index = 0; index < laplacianScales.Length; index++)
                    rows.Add(BuildFactorValidationRow(model, KernelFamily.Laplacian, dimension, laplacianScales[index], seed: 2000 + ((int)model + 1) * 10000 + dimension * 100 + index));
            }
        }

        return rows;
    }

    private static IReadOnlyList<IntrinsicSelfConsistencyReport> BuildIntrinsicSelfConsistencyReports()
    {
        int[] dimensions = { 2, 3, 8 };
        double[] scales = { 0.10, 0.20, 0.35 };
        FactorValidationModel[] models = ResolveValidationModels();
        var reports = new List<IntrinsicSelfConsistencyReport>(dimensions.Length * scales.Length * models.Length);

        foreach (FactorValidationModel model in models)
        {
            foreach (int dimension in dimensions)
            {
                for (int index = 0; index < scales.Length; index++)
                {
                    reports.Add(BuildIntrinsicSelfConsistencyReport(
                        model,
                        dimension,
                        scales[index],
                        seed: 3000 + ((int)model + 1) * 10000 + dimension * 100 + index));
                }
            }
        }

        return reports;
    }

    private static IntrinsicSelfConsistencyReport BuildIntrinsicSelfConsistencyReport(
        FactorValidationModel model,
        int dimension,
        double scale,
        int seed)
    {
        NearestNeighborSample sample = model switch
        {
            FactorValidationModel.CandidatePoolMinimum => SimulateNearestNeighborSample(KernelFamily.IntrinsicGaussian, dimension, scale, SamplesPerCombo, seed),
            FactorValidationModel.InhomogeneousPointProcess => SimulatePointProcessNearestNeighborSample(KernelFamily.IntrinsicGaussian, dimension, scale, SamplesPerCombo, seed),
            _ => throw new ArgumentOutOfRangeException(nameof(model), model, "Unknown intrinsic calibration model."),
        };

        double meanSquaredRadius = sample.Values.Select(value => value * value).Average();
        double recoveredBandwidth = BandwidthEstimation.ForIntrinsicGaussianHyperbolic(sample.Values, dimension);

        return new IntrinsicSelfConsistencyReport(
            Model: model,
            Dimension: dimension,
            Scale: scale,
            CandidatePoolSize: sample.CandidatePoolSize,
            NeighborRank: sample.NeighborRank,
            RadiusMax: sample.RadiusMax,
            AcceptanceRate: sample.AcceptanceRate,
            MeanSquaredRadius: meanSquaredRadius,
            RecoveredBandwidth: recoveredBandwidth,
            RelativeErrorRecovered: RelativeError(recoveredBandwidth, scale));
    }

    private static IReadOnlyList<IntrinsicGroundTruthReport> BuildIntrinsicGroundTruthReports()
    {
        int[] dimensions = { 2, 3, 8 };
        int[] actualKs = { 16, 32, 48 };
        double[] trueSigmas = { 0.20, 0.35 };
        var reports = new List<IntrinsicGroundTruthReport>(dimensions.Length * actualKs.Length * trueSigmas.Length);

        foreach (int ambientDimension in dimensions)
        {
            foreach (int actualK in actualKs)
            {
                foreach (double trueSigma in trueSigmas)
                {
                    reports.Add(BuildIntrinsicGroundTruthReport(ambientDimension, actualK, trueSigma));
                }
            }
        }

        return reports;
    }

    private static IntrinsicGroundTruthReport BuildIntrinsicGroundTruthReport(int ambientDimension, int actualK, double trueSigma)
    {
        double[][] features = SamplePoincareGaussianFixture(
            ambientDimension,
            trueSigma,
            SamplesPerCombo,
            seed: 9000 + ambientDimension * 100 + actualK + (int)(trueSigma * 1000.0));

        GraphBuildResult build = SpcGraphBuilder.BuildResult(
            features,
            CreateIntrinsicGroundTruthHarnessConfig(actualK),
            new PoincareMetric());

        double recoveredBandwidth = build.SingleBandwidth
            ?? throw new InvalidOperationException($"Expected a single recovered bandwidth for the ground-truth intrinsic harness at d={ambientDimension}, k={actualK}.");

        return new IntrinsicGroundTruthReport(
            AmbientDimension: ambientDimension,
            ActualK: actualK,
            TrueSigma: trueSigma,
            MeanKthNeighborDistance: build.DirectedSelection.KthNeighborDistances.Average(),
            RecoveredBandwidth: recoveredBandwidth,
            RelativeErrorRecovered: RelativeError(recoveredBandwidth, trueSigma));
    }

    private static FactorValidationRow BuildFactorValidationRow(FactorValidationModel model, KernelFamily family, int dimension, double scale, int seed)
    {
        NearestNeighborSample sample = model switch
        {
            FactorValidationModel.CandidatePoolMinimum => SimulateNearestNeighborSample(family, dimension, scale, SamplesPerCombo, seed),
            FactorValidationModel.InhomogeneousPointProcess => SimulatePointProcessNearestNeighborSample(family, dimension, scale, SamplesPerCombo, seed),
            _ => throw new ArgumentOutOfRangeException(nameof(model), model, "Unknown factor validation model."),
        };
        double[] scratch = new double[sample.Values.Length];
        double configuredFactor = GetConfiguredFactor(family);
        double euclideanFactor = GetEuclideanFactor(family);
        double recoveredConfigured = BandwidthEstimation.LogScaleBandwidth(sample.Values, scratch, configuredFactor);
        double recoveredEuclidean = BandwidthEstimation.LogScaleBandwidth(sample.Values, scratch, euclideanFactor);
        (double logMedian, double logMad) = ComputeLogMedianAndMad(sample.Values);
        double impliedFactor = logMad <= 1e-12
            ? double.NaN
            : (Math.Log(scale) - logMedian) / logMad;

        return new FactorValidationRow(
            Model: model,
            Family: family,
            Dimension: dimension,
            Scale: scale,
            VolumePressure: scale * (dimension - 1),
            DecayExcess: (1.0 / scale) - (dimension - 1),
            CandidatePoolSize: sample.CandidatePoolSize,
            NeighborRank: sample.NeighborRank,
            TruncatedProbe: sample.TruncatedProbe,
            RadiusMax: sample.RadiusMax,
            AcceptanceRate: sample.AcceptanceRate,
            LogMedian: logMedian,
            LogMad: logMad,
            ImpliedFactor: impliedFactor,
            ConfiguredFactor: configuredFactor,
            EuclideanFactor: euclideanFactor,
            RecoveredConfigured: recoveredConfigured,
            RecoveredEuclidean: recoveredEuclidean,
            RelativeErrorConfigured: RelativeError(recoveredConfigured, scale),
            RelativeErrorEuclidean: RelativeError(recoveredEuclidean, scale));
    }

    private static IReadOnlyList<IntrinsicAbReport> BuildIntrinsicAbReports()
    {
        int[] dimensions = { 2, 3, 8 };
        var reports = new List<IntrinsicAbReport>(dimensions.Length);

        foreach (int ambientDimension in dimensions)
            reports.Add(BuildIntrinsicAbReport(ambientDimension));

        return reports;
    }

    private static IntrinsicAbReport BuildIntrinsicAbReport(int ambientDimension)
    {
        double[][] features = CreateLongRangePoincareFixture(ambientDimension);
        var metric = new PoincareMetric();

        GraphBuildResult linear = SpcGraphBuilder.BuildResult(
            features,
            CreateIntrinsicHarnessConfig(CouplingFidelity.GeodesicLinear),
            metric);
        GraphBuildResult intrinsic = SpcGraphBuilder.BuildResult(
            features,
            CreateIntrinsicHarnessConfig(CouplingFidelity.Intrinsic),
            metric);

        double bandwidth = linear.SingleBandwidth
            ?? throw new InvalidOperationException($"Expected a resolved bandwidth for the intrinsic A/B harness at d={ambientDimension}.");

        Dictionary<long, double> linearWeights = BuildUndirectedWeightMap(linear.Graph);
        Dictionary<long, double> intrinsicWeights = BuildUndirectedWeightMap(intrinsic.Graph);
        var rows = new List<IntrinsicEdgeRow>();

        for (int left = 0; left < features.Length; left++)
        {
            for (int right = left + 1; right < features.Length; right++)
            {
                long key = EdgeKey(left, right);
                double distance = metric.Distance(features[left], features[right]);
                bool linearKept = linearWeights.TryGetValue(key, out double linearWeight);
                bool intrinsicKept = intrinsicWeights.TryGetValue(key, out double intrinsicWeight);
                double legacyPocWeight = GaussianKernel.Evaluate(distance, bandwidth) * LegacyIntrinsicPocCorrection(distance);

                rows.Add(new IntrinsicEdgeRow(
                    AmbientDimension: ambientDimension,
                    Left: left,
                    Right: right,
                    Distance: distance,
                    LinearKept: linearKept,
                    LinearWeight: linearKept ? linearWeight : 0.0,
                    IntrinsicKept: intrinsicKept,
                    IntrinsicWeight: intrinsicKept ? intrinsicWeight : 0.0,
                    LegacyPocWeight: legacyPocWeight));
            }
        }

        rows.Sort((left, right) => left.Distance.CompareTo(right.Distance));

        IntrinsicEdgeRow shortEdge = rows.First(row => row.LinearKept && row.IntrinsicKept);
        IntrinsicEdgeRow farthestLinearKept = rows.Last(row => row.LinearKept);

        return new IntrinsicAbReport(
            AmbientDimension: ambientDimension,
            Bandwidth: linear.SingleBandwidth,
            LinearUndirectedEdgeCount: linearWeights.Count,
            IntrinsicUndirectedEdgeCount: intrinsicWeights.Count,
            DroppedByIntrinsicCount: rows.Count(row => row.LinearKept && !row.IntrinsicKept),
            KeepThreshold: KeepThreshold,
            ShortEdgeDistance: shortEdge.Distance,
            ShortEdgeDifference: Math.Abs(shortEdge.LinearWeight - shortEdge.IntrinsicWeight),
            FarthestLinearKeptDistance: farthestLinearKept.Distance,
            FarthestLinearKeptLinearWeight: farthestLinearKept.LinearWeight,
            FarthestLinearKeptIntrinsicWeight: farthestLinearKept.IntrinsicWeight,
            FarthestLinearKeptLegacyPocWeight: farthestLinearKept.LegacyPocWeight,
            Rows: rows);
    }

    private static string BuildSummary(
        IReadOnlyList<FactorValidationRow> factorRows,
        IReadOnlyList<IntrinsicSelfConsistencyReport> intrinsicSelfConsistencyReports,
        IReadOnlyList<IntrinsicGroundTruthReport> intrinsicGroundTruthReports,
        IReadOnlyList<IntrinsicAbReport> intrinsicReports)
    {
        FactorValidationModel[] models = ResolveValidationModels();

        var builder = new StringBuilder();
        builder.AppendLine("Hyperbolic factor validation");
        builder.AppendLine($"Gaussian configured factor = {BandwidthEstimation.GaussianHyperbolicFactor.ToString("F4", CultureInfo.InvariantCulture)}");
        builder.AppendLine($"Laplacian configured factor = {BandwidthEstimation.LaplacianHyperbolicFactor.ToString("F4", CultureInfo.InvariantCulture)}");
        builder.AppendLine($"Validation models = {string.Join(", ", models)}");
        builder.AppendLine($"NN surrogate pool size = {CandidatePoolSize.ToString(CultureInfo.InvariantCulture)}");
        builder.AppendLine($"NN surrogate rank = {NeighborRank.ToString(CultureInfo.InvariantCulture)}");
        builder.AppendLine(
            $"Laplacian calibration window = d={LaplacianCalibrationDimension}, volume pressure in [{LaplacianCalibrationMinVolumePressure.ToString("F2", CultureInfo.InvariantCulture)}, {LaplacianCalibrationMaxVolumePressure.ToString("F2", CultureInfo.InvariantCulture)}], tolerance ±{LaplacianCalibrationTolerance.ToString("F2", CultureInfo.InvariantCulture)}");
        builder.AppendLine();
        foreach (FactorValidationModel model in models)
        {
            FactorValidationRow[] modelRows = factorRows
                .Where(row => row.Model == model)
                .OrderBy(row => row.Family)
                .ThenBy(row => row.Dimension)
                .ThenBy(row => row.Scale)
                .ToArray();

            double laplacianBulkMeanFactor = modelRows
                .Where(row => row.Family == KernelFamily.Laplacian && row.VolumePressure >= 0.25 && row.VolumePressure <= 0.80)
                .Average(row => row.ImpliedFactor);
            double laplacianSmallMeanFactor = modelRows
                .Where(row => row.Family == KernelFamily.Laplacian && row.VolumePressure <= 0.10)
                .Average(row => row.ImpliedFactor);
            double laplacianCalibrationMeanFactor = modelRows
                .Where(IsLaplacianCalibrationRow)
                .Average(row => row.ImpliedFactor);

            builder.AppendLine($"Model = {model}");
            builder.AppendLine($"Laplacian small-scale mean implied factor = {laplacianSmallMeanFactor.ToString("F4", CultureInfo.InvariantCulture)}");
            builder.AppendLine($"Laplacian bulk mean implied factor = {laplacianBulkMeanFactor.ToString("F4", CultureInfo.InvariantCulture)}");
            builder.AppendLine($"Laplacian calibration mean implied factor = {laplacianCalibrationMeanFactor.ToString("F4", CultureInfo.InvariantCulture)}");
            builder.AppendLine("Model\tFamily\tDim\tScale\tVolPressure\tLocalPopulation\tRank\tImplied\tCfgRecovered\tCfgRelErr\tEucRecovered\tEucRelErr\tTruncated");
            foreach (FactorValidationRow row in modelRows)
            {
                builder.AppendLine(string.Join("\t", new object?[]
                {
                    row.Model,
                    row.Family,
                    row.Dimension.ToString(CultureInfo.InvariantCulture),
                    row.Scale.ToString("F3", CultureInfo.InvariantCulture),
                    row.VolumePressure.ToString("F3", CultureInfo.InvariantCulture),
                    row.CandidatePoolSize.ToString(CultureInfo.InvariantCulture),
                    row.NeighborRank.ToString(CultureInfo.InvariantCulture),
                    row.ImpliedFactor.ToString("F4", CultureInfo.InvariantCulture),
                    row.RecoveredConfigured.ToString("F4", CultureInfo.InvariantCulture),
                    row.RelativeErrorConfigured.ToString("F4", CultureInfo.InvariantCulture),
                    row.RecoveredEuclidean.ToString("F4", CultureInfo.InvariantCulture),
                    row.RelativeErrorEuclidean.ToString("F4", CultureInfo.InvariantCulture),
                    row.TruncatedProbe ? "truncated" : "full",
                }));
            }

            builder.AppendLine();
        }

        builder.AppendLine("Intrinsic self-consistency recovery");
        builder.AppendLine($"Tolerance = ±{IntrinsicCalibrationTolerance.ToString("F2", CultureInfo.InvariantCulture)} relative error");
        builder.AppendLine("Model\tDim\tScale\tRecovered\tRelErr\tMeanSq\tLocalPopulation\tRank");
        foreach (IntrinsicSelfConsistencyReport report in intrinsicSelfConsistencyReports
                     .OrderBy(report => report.Model)
                     .ThenBy(report => report.Dimension)
                     .ThenBy(report => report.Scale))
        {
            builder.AppendLine(string.Join("\t", new object?[]
            {
                report.Model,
                report.Dimension.ToString(CultureInfo.InvariantCulture),
                report.Scale.ToString("F3", CultureInfo.InvariantCulture),
                report.RecoveredBandwidth.ToString("F4", CultureInfo.InvariantCulture),
                report.RelativeErrorRecovered.ToString("F4", CultureInfo.InvariantCulture),
                report.MeanSquaredRadius.ToString("F4", CultureInfo.InvariantCulture),
                report.CandidatePoolSize.ToString(CultureInfo.InvariantCulture),
                report.NeighborRank.ToString(CultureInfo.InvariantCulture),
            }));
        }

        builder.AppendLine();

        builder.AppendLine("Intrinsic ground-truth validation");
        builder.AppendLine("Dim\tActualK\tSigma\tRecovered\tRelErr\tMeanKthNN");
        foreach (IntrinsicGroundTruthReport report in intrinsicGroundTruthReports
                     .OrderBy(report => report.AmbientDimension)
                     .ThenBy(report => report.ActualK)
                     .ThenBy(report => report.TrueSigma))
        {
            builder.AppendLine(string.Join("\t", new object?[]
            {
                report.AmbientDimension.ToString(CultureInfo.InvariantCulture),
                report.ActualK.ToString(CultureInfo.InvariantCulture),
                report.TrueSigma.ToString("F3", CultureInfo.InvariantCulture),
                report.RecoveredBandwidth.ToString("F4", CultureInfo.InvariantCulture),
                report.RelativeErrorRecovered.ToString("F4", CultureInfo.InvariantCulture),
                report.MeanKthNeighborDistance.ToString("F4", CultureInfo.InvariantCulture),
            }));
        }

        builder.AppendLine();

        builder.AppendLine("Intrinsic k-sweep summary");
        builder.AppendLine("ActualK\tCount\tMeanRecovered\tMeanRelErr\tMinRelErr\tMaxRelErr");
        foreach (var kGroup in intrinsicGroundTruthReports
                     .GroupBy(report => report.ActualK)
                     .OrderBy(group => group.Key))
        {
            builder.AppendLine(string.Join("\t", new object?[]
            {
                kGroup.Key.ToString(CultureInfo.InvariantCulture),
                kGroup.Count().ToString(CultureInfo.InvariantCulture),
                kGroup.Average(report => report.RecoveredBandwidth).ToString("F4", CultureInfo.InvariantCulture),
                kGroup.Average(report => report.RelativeErrorRecovered).ToString("F4", CultureInfo.InvariantCulture),
                kGroup.Min(report => report.RelativeErrorRecovered).ToString("F4", CultureInfo.InvariantCulture),
                kGroup.Max(report => report.RelativeErrorRecovered).ToString("F4", CultureInfo.InvariantCulture),
            }));
        }

        builder.AppendLine();

        builder.AppendLine("Intrinsic vs GeodesicLinear A/B");
        foreach (IntrinsicAbReport intrinsicReport in intrinsicReports.OrderBy(report => report.AmbientDimension))
        {
            builder.AppendLine($"Dimension\t{intrinsicReport.AmbientDimension}");
            builder.AppendLine($"Bandwidth\t{intrinsicReport.Bandwidth?.ToString("F4", CultureInfo.InvariantCulture) ?? "null"}");
            builder.AppendLine($"LinearEdges\t{intrinsicReport.LinearUndirectedEdgeCount}");
            builder.AppendLine($"IntrinsicEdges\t{intrinsicReport.IntrinsicUndirectedEdgeCount}");
            builder.AppendLine($"DroppedByIntrinsic\t{intrinsicReport.DroppedByIntrinsicCount}");
            builder.AppendLine($"ShortEdgeDiff\t{intrinsicReport.ShortEdgeDifference.ToString("G6", CultureInfo.InvariantCulture)}");
            builder.AppendLine($"FarthestIntrinsicOverLinear\t{(intrinsicReport.FarthestLinearKeptIntrinsicWeight / intrinsicReport.FarthestLinearKeptLinearWeight).ToString("G6", CultureInfo.InvariantCulture)}");
            builder.AppendLine($"FarthestLegacyPocWeight\t{intrinsicReport.FarthestLinearKeptLegacyPocWeight.ToString("G6", CultureInfo.InvariantCulture)}");
            builder.AppendLine("Dim\tLeft\tRight\tDistance\tLinearKept\tLinearWeight\tIntrinsicKept\tIntrinsicWeight\tLegacyPocWeight");
            foreach (IntrinsicEdgeRow row in intrinsicReport.Rows)
            {
                builder.AppendLine(string.Join("\t", new object?[]
                {
                    row.AmbientDimension.ToString(CultureInfo.InvariantCulture),
                    row.Left.ToString(CultureInfo.InvariantCulture),
                    row.Right.ToString(CultureInfo.InvariantCulture),
                    row.Distance.ToString("F4", CultureInfo.InvariantCulture),
                    row.LinearKept ? "yes" : "no",
                    row.LinearWeight.ToString("G6", CultureInfo.InvariantCulture),
                    row.IntrinsicKept ? "yes" : "no",
                    row.IntrinsicWeight.ToString("G6", CultureInfo.InvariantCulture),
                    row.LegacyPocWeight.ToString("G6", CultureInfo.InvariantCulture),
                }));
            }

            builder.AppendLine();
        }

        return builder.ToString();
    }

    private static NearestNeighborSample SimulateNearestNeighborSample(KernelFamily family, int dimension, double scale, int count, int seed)
    {
        double radiusMax = ResolveRadiusMax(family, dimension, scale);
        double logMax = EstimateLogDensityMax(family, dimension, scale, radiusMax);
        var rng = new Random(seed);
        var samples = new double[count];
        int accepted = 0;
        int attempts = 0;
        int maxAttempts = count * CandidatePoolSize * 20000;
        var pool = new double[CandidatePoolSize];

        while (accepted < count)
        {
            if (attempts >= maxAttempts)
            {
                throw new InvalidOperationException(
                    $"Rejection sampler stalled for {family} d={dimension} scale={scale:G4} rMax={radiusMax:G4}.");
            }

            int poolFill = 0;
            while (poolFill < CandidatePoolSize)
            {
                if (attempts >= maxAttempts)
                {
                    throw new InvalidOperationException(
                        $"Nearest-neighbor sampler stalled for {family} d={dimension} scale={scale:G4} rMax={radiusMax:G4}.");
                }

                attempts++;
                double radius = rng.NextDouble() * radiusMax;
                double logAccept = LogDensity(family, dimension, scale, radius) - logMax;
                if (Math.Log(rng.NextDouble()) <= logAccept)
                    pool[poolFill++] = radius;
            }

            Array.Sort(pool);
            samples[accepted++] = pool[NeighborRank - 1];
        }

        return new NearestNeighborSample(
            Values: samples,
            CandidatePoolSize: CandidatePoolSize,
            NeighborRank: NeighborRank,
            RadiusMax: radiusMax,
            AcceptanceRate: (count * CandidatePoolSize) / (double)attempts,
            TruncatedProbe: family == KernelFamily.Laplacian && scale * (dimension - 1) >= 1.0);
    }

    private static NearestNeighborSample SimulatePointProcessNearestNeighborSample(KernelFamily family, int dimension, double scale, int count, int seed)
    {
        double radiusMax = ResolveRadiusMax(family, dimension, scale);
        (double[] radiusGrid, double[] cumulativeGrid, double totalMass) = BuildCumulativeMassGrid(family, dimension, scale, radiusMax);
        double localPopulation = CandidatePoolSize;
        double truncatedProbability = Math.Exp(-localPopulation);
        double normalization = 1.0 - truncatedProbability;
        var rng = new Random(seed);
        var samples = new double[count];

        for (int index = 0; index < count; index++)
        {
            double u = rng.NextDouble();
            double targetMass = -Math.Log(1.0 - u * normalization) * (totalMass / localPopulation);
            samples[index] = InvertCumulativeMass(radiusGrid, cumulativeGrid, targetMass);
        }

        return new NearestNeighborSample(
            Values: samples,
            CandidatePoolSize: CandidatePoolSize,
            NeighborRank: NeighborRank,
            RadiusMax: radiusMax,
            AcceptanceRate: 1.0 - truncatedProbability,
            TruncatedProbe: family == KernelFamily.Laplacian && scale * (dimension - 1) >= 1.0);
    }

    private static (double[] RadiusGrid, double[] CumulativeGrid, double TotalMass) BuildCumulativeMassGrid(
        KernelFamily family,
        int dimension,
        double scale,
        double radiusMax)
    {
        const int GridCount = 4096;
        var radiusGrid = new double[GridCount + 1];
        var cumulativeGrid = new double[GridCount + 1];
        double previousRadius = 0.0;
        double previousWeight = Math.Exp(LogDensity(family, dimension, scale, 0.0));

        radiusGrid[0] = 0.0;
        cumulativeGrid[0] = 0.0;

        for (int index = 1; index <= GridCount; index++)
        {
            double radius = radiusMax * index / GridCount;
            double weight = Math.Exp(LogDensity(family, dimension, scale, radius));
            radiusGrid[index] = radius;
            cumulativeGrid[index] = cumulativeGrid[index - 1] + 0.5 * (weight + previousWeight) * (radius - previousRadius);
            previousRadius = radius;
            previousWeight = weight;
        }

        return (radiusGrid, cumulativeGrid, cumulativeGrid[GridCount]);
    }

    private static double InvertCumulativeMass(double[] radiusGrid, double[] cumulativeGrid, double targetMass)
    {
        int index = Array.BinarySearch(cumulativeGrid, targetMass);
        if (index >= 0)
            return radiusGrid[index];

        index = ~index;
        if (index <= 0)
            return radiusGrid[0];
        if (index >= cumulativeGrid.Length)
            return radiusGrid[^1];

        double leftMass = cumulativeGrid[index - 1];
        double rightMass = cumulativeGrid[index];
        double t = (targetMass - leftMass) / Math.Max(1e-12, rightMass - leftMass);
        return radiusGrid[index - 1] + (radiusGrid[index] - radiusGrid[index - 1]) * t;
    }

    private static double EstimateLogDensityMax(KernelFamily family, int dimension, double scale, double radiusMax)
    {
        const int GridCount = 4096;
        double logMax = double.NegativeInfinity;

        for (int index = 0; index <= GridCount; index++)
        {
            double radius = radiusMax * index / GridCount;
            double logDensity = LogDensity(family, dimension, scale, radius);
            if (logDensity > logMax)
                logMax = logDensity;
        }

        return logMax;
    }

    private static double LogDensity(KernelFamily family, int dimension, double scale, double radius)
    {
        double kernelTerm = family switch
        {
            KernelFamily.Gaussian => -(radius * radius) / (2.0 * scale * scale),
            KernelFamily.IntrinsicGaussian => -(radius * radius) / (2.0 * scale * scale),
            KernelFamily.Laplacian => -radius / scale,
            _ => throw new ArgumentOutOfRangeException(nameof(family), family, "Unknown kernel family."),
        };

        if (family == KernelFamily.IntrinsicGaussian)
        {
            if (dimension <= 1)
                return kernelTerm;
            if (radius <= 0.0)
                return double.NegativeInfinity;

            return kernelTerm + 0.5 * (dimension - 1) * (Math.Log(Math.Max(radius, 1e-12)) + LogSinh(radius));
        }

        return dimension <= 1
            ? kernelTerm
            : kernelTerm + (dimension - 1) * LogSinh(radius);
    }

    private static double ResolveRadiusMax(KernelFamily family, int dimension, double scale)
    {
        return family switch
        {
            KernelFamily.Gaussian => Math.Max(1.0, 12.0 * scale),
            KernelFamily.IntrinsicGaussian => Math.Max(8.0, 14.0 * scale + dimension),
            KernelFamily.Laplacian => scale * (dimension - 1) < 1.0
                ? Math.Max(8.0, 10.0 + (2.0 / Math.Max(1e-6, (1.0 / scale) - (dimension - 1))))
                : Math.Max(12.0, 16.0 * scale + dimension),
            _ => throw new ArgumentOutOfRangeException(nameof(family), family, "Unknown kernel family."),
        };
    }

    private static (double LogMedian, double LogMad) ComputeLogMedianAndMad(double[] sample)
    {
        double[] work = new double[sample.Length];
        for (int index = 0; index < sample.Length; index++)
            work[index] = Math.Log(Math.Max(sample[index], 1e-12));

        Array.Sort(work);
        double logMedian = MedianOfSorted(work);

        for (int index = 0; index < sample.Length; index++)
            work[index] = Math.Abs(Math.Log(Math.Max(sample[index], 1e-12)) - logMedian);

        Array.Sort(work);
        return (logMedian, MedianOfSorted(work));
    }

    private static double MedianOfSorted(double[] sorted)
    {
        int mid = sorted.Length / 2;
        return (sorted.Length & 1) == 0
            ? 0.5 * (sorted[mid - 1] + sorted[mid])
            : sorted[mid];
    }

    private static double LogSinh(double radius)
    {
        if (radius < 1e-8)
            return Math.Log(Math.Max(radius, 1e-12));

        if (radius < 20.0)
            return Math.Log(Math.Sinh(radius));

        return radius - Math.Log(2.0);
    }

    private static double RelativeError(double estimate, double truth)
        => Math.Abs(estimate - truth) / Math.Max(1e-12, truth);

    private static bool IsLaplacianCalibrationRow(FactorValidationRow row)
        => row.Family == KernelFamily.Laplacian
           && row.Dimension == LaplacianCalibrationDimension
           && !row.TruncatedProbe
           && row.VolumePressure >= LaplacianCalibrationMinVolumePressure
           && row.VolumePressure <= LaplacianCalibrationMaxVolumePressure;

    private static FactorValidationModel[] ResolveValidationModels()
    {
        string? raw = Environment.GetEnvironmentVariable(ValidationModelsEnvVar);
        if (string.IsNullOrWhiteSpace(raw))
            return new[] { FactorValidationModel.CandidatePoolMinimum, FactorValidationModel.InhomogeneousPointProcess };

        var models = new List<FactorValidationModel>();
        foreach (string token in raw.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
        {
            if (token.Equals("candidate-pool", StringComparison.OrdinalIgnoreCase)
                || token.Equals("pool", StringComparison.OrdinalIgnoreCase))
            {
                models.Add(FactorValidationModel.CandidatePoolMinimum);
            }
            else if (token.Equals("point-process", StringComparison.OrdinalIgnoreCase)
                     || token.Equals("poisson", StringComparison.OrdinalIgnoreCase)
                     || token.Equals("point", StringComparison.OrdinalIgnoreCase))
            {
                models.Add(FactorValidationModel.InhomogeneousPointProcess);
            }
            else if (token.Equals("both", StringComparison.OrdinalIgnoreCase))
            {
                models.Add(FactorValidationModel.CandidatePoolMinimum);
                models.Add(FactorValidationModel.InhomogeneousPointProcess);
            }
            else
            {
                throw new ArgumentException(
                    $"Unknown factor validation model token '{token}'. Supported: candidate-pool, point-process, both.");
            }
        }

        return models.Distinct().ToArray();
    }

    private static double GetConfiguredFactor(KernelFamily family)
    {
        return family switch
        {
            KernelFamily.Gaussian => BandwidthEstimation.GaussianHyperbolicFactor,
            KernelFamily.Laplacian => BandwidthEstimation.LaplacianHyperbolicFactor,
            _ => throw new ArgumentOutOfRangeException(nameof(family), family, "Unknown kernel family."),
        };
    }

    private static double GetEuclideanFactor(KernelFamily family)
    {
        return family switch
        {
            KernelFamily.Gaussian => BandwidthEstimation.GaussianFactor,
            KernelFamily.Laplacian => BandwidthEstimation.LaplacianFactor,
            _ => throw new ArgumentOutOfRangeException(nameof(family), family, "Unknown kernel family."),
        };
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

    private static double LegacyIntrinsicPocCorrection(double distance)
    {
        if (distance < 1e-12)
            return 1.0;

        return distance / Math.Sinh(distance);
    }

    private static GraphCompilerConfig CreateIntrinsicHarnessConfig(CouplingFidelity fidelity) => new()
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

    private static GraphCompilerConfig CreateIntrinsicGroundTruthHarnessConfig(int k) => new()
    {
        Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = k },
        Filter = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
        Repair = new RepairConfig { Kind = RepairKind.NoRepair },
        Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
        Projection = new CouplingProjection
        {
            Kernel = new Gaussian(0.0),
            LmpRescale = false,
            Fidelity = CouplingFidelity.Intrinsic,
        },
    };

    private static double[][] SamplePoincareGaussianFixture(int ambientDimension, double sigma, int count, int seed)
    {
        var manifold = new Maths.Geometry.PoincareBallManifold(ambientDimension);
        var rng = new Random(seed);
        var origin = new double[ambientDimension];
        var tangent = new double[ambientDimension];
        var features = new double[count][];

        for (int index = 0; index < count; index++)
        {
            for (int dimension = 0; dimension < ambientDimension; dimension++)
                tangent[dimension] = SampleStandardNormal(rng) * sigma;

            var point = new double[ambientDimension];
            manifold.ExpMap(origin, tangent, point);
            features[index] = point;
        }

        return features;
    }

    private static double SampleStandardNormal(Random rng)
    {
        double u1 = Math.Max(rng.NextDouble(), 1e-12);
        double u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

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

    private sealed record NearestNeighborSample(
        double[] Values,
        int CandidatePoolSize,
        int NeighborRank,
        double RadiusMax,
        double AcceptanceRate,
        bool TruncatedProbe);

    private sealed record FactorValidationRow(
        FactorValidationModel Model,
        KernelFamily Family,
        int Dimension,
        double Scale,
        double VolumePressure,
        double DecayExcess,
        int CandidatePoolSize,
        int NeighborRank,
        bool TruncatedProbe,
        double RadiusMax,
        double AcceptanceRate,
        double LogMedian,
        double LogMad,
        double ImpliedFactor,
        double ConfiguredFactor,
        double EuclideanFactor,
        double RecoveredConfigured,
        double RecoveredEuclidean,
        double RelativeErrorConfigured,
        double RelativeErrorEuclidean);

    private sealed record IntrinsicSelfConsistencyReport(
        FactorValidationModel Model,
        int Dimension,
        double Scale,
        int CandidatePoolSize,
        int NeighborRank,
        double RadiusMax,
        double AcceptanceRate,
        double MeanSquaredRadius,
        double RecoveredBandwidth,
        double RelativeErrorRecovered);

    private sealed record IntrinsicGroundTruthReport(
        int AmbientDimension,
        int ActualK,
        double TrueSigma,
        double MeanKthNeighborDistance,
        double RecoveredBandwidth,
        double RelativeErrorRecovered);

    private sealed record IntrinsicAbReport(
        int AmbientDimension,
        double? Bandwidth,
        int LinearUndirectedEdgeCount,
        int IntrinsicUndirectedEdgeCount,
        int DroppedByIntrinsicCount,
        double KeepThreshold,
        double ShortEdgeDistance,
        double ShortEdgeDifference,
        double FarthestLinearKeptDistance,
        double FarthestLinearKeptLinearWeight,
        double FarthestLinearKeptIntrinsicWeight,
        double FarthestLinearKeptLegacyPocWeight,
        IReadOnlyList<IntrinsicEdgeRow> Rows);

    private sealed record IntrinsicEdgeRow(
        int AmbientDimension,
        int Left,
        int Right,
        double Distance,
        bool LinearKept,
        double LinearWeight,
        bool IntrinsicKept,
        double IntrinsicWeight,
        double LegacyPocWeight);

    private enum KernelFamily
    {
        Gaussian,
        IntrinsicGaussian,
        Laplacian,
    }

    private enum FactorValidationModel
    {
        CandidatePoolMinimum,
        InhomogeneousPointProcess,
    }
}
