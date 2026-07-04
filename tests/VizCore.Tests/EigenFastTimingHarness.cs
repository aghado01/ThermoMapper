using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.LinAlg;
using Repo.TestHarness;
using Synthetic;
using Synthetic.Euclidean;
using TDA.Ph;
using Xunit;
using Xunit.Abstractions;

namespace VizCore.Tests;

public sealed class EigenFastTimingHarness
{
    private static int _spectralSink;
    private readonly ITestOutputHelper _output;

    public EigenFastTimingHarness(ITestOutputHelper output)
    {
        _output = output;
    }

    [Fact]
    [Trait("Category", "Benchmark")]
    public void CompareReferenceAndFastTimings()
    {
#if DEBUG
        _output.WriteLine("Benchmark harness only reports meaningful timings in Release builds.");
        return;
#else
        var run = HarnessArtifacts.Create(
            runKind: "test-runs",
            suiteName: nameof(EigenFastTimingHarness),
            runName: nameof(CompareReferenceAndFastTimings),
            metadata: new Dictionary<string, object?>
            {
                ["DenseEigenBackend"] = GetDenseEigenBackendLabel(),
                ["Configuration"] = "Release",
                ["EnableBenchmarks"] = true,
            });

        string directPath = run.WriteRunJson("direct-comparison", Array.Empty<DirectBenchmarkCase>());
        string variantPath = run.WriteRunJson("fast-variant-comparison", Array.Empty<FastVariantBenchmarkCase>());
        string spectralPath = run.WriteRunJson("spectral-dispatch-comparison", Array.Empty<SpectralBenchmarkCase>());
        string solverPath = run.WriteRunJson("spectral-solver-comparison", Array.Empty<SpectralSolverComparisonCase>());
        string denseBottomKPath = run.WriteRunJson("dense-bottomk-comparison", Array.Empty<DenseBottomKComparisonCase>());
        string summaryPath = run.WriteRunText("summary", "pending");

        var directCases = new List<DirectBenchmarkCase>();
        var variantCases = new List<FastVariantBenchmarkCase>();
        var spectralCases = new List<SpectralBenchmarkCase>();
        var solverCases = new List<SpectralSolverComparisonCase>();
        var denseBottomKCases = new List<DenseBottomKComparisonCase>();

        foreach (int size in new[] { 32, 64, 128, 256, 512 })
        {
            double[,] matrix = BuildRandomSymmetricMatrix(size, seed: 10_000 + size);
            directCases.Add(MeasureCase($"rand-{size}", matrix));
            if (size >= 64)
                denseBottomKCases.Add(MeasureDenseBottomKComparison($"rand-{size}", matrix, k: 8));
        }

        DisconnectedControl.Fixture control = DisconnectedControl.Generate(pointsPerComponent: 64, separation: 100.0);
        CsrGraph graph = FixtureGoldenHelpers.BuildDraftGraph(control.Points, k: 4, ensureConnected: true);
        directCases.Add(MeasureCase("lap-128", BuildCombinatorialLaplacian(graph)));

#if !EIGEN_REFERENCE
        foreach (int size in new[] { 128, 256, 512 })
        {
            double[,] matrix = BuildRandomSymmetricMatrix(size, seed: 20_000 + size);
            variantCases.Add(MeasureFastVariantCase($"rand-{size}", matrix));
        }
#endif

        double[][] crescent = CrescentAndEllipsoid.Generate(crescentPoints: 500, ellipsoidPoints: 0, seed: 42).Features;
        CsrGraph crescentGraph = FixtureGoldenHelpers.BuildDraftGraph(crescent, k: FixtureGoldenHelpers.ConnectedFixtureK, ensureConnected: true);
        spectralCases.Add(MeasureSpectralCase(
            "spectral-crescent-500",
            crescentGraph,
            seed: FixtureGoldenHelpers.SpectralSeed,
            k: 8,
            denseOptions: default,
            materialization: DenseLaplacianMaterialization.Rectangular));
        spectralCases.Add(MeasureSpectralCase(
            "spectral-crescent-500",
            crescentGraph,
            seed: FixtureGoldenHelpers.SpectralSeed,
            k: 8,
            denseOptions: default,
            materialization: DenseLaplacianMaterialization.FlatColumnMajor));

#if !EIGEN_REFERENCE
        spectralCases.Add(MeasureSpectralCase(
            "spectral-crescent-500",
            crescentGraph,
            seed: FixtureGoldenHelpers.SpectralSeed,
            k: 8,
            denseOptions: new DenseEigenOptions(DenseEigenFastVariant.Fma),
            materialization: DenseLaplacianMaterialization.FlatColumnMajor));
#endif

        double[][] mobius = MobiusAndEllipsoid.Generate(mobiusPoints: 500, ellipsoidPoints: 0, seed: 42).Features;
        CsrGraph mobiusGraph = FixtureGoldenHelpers.BuildDraftGraph(mobius, k: FixtureGoldenHelpers.ConnectedFixtureK, ensureConnected: true);
        spectralCases.Add(MeasureSpectralCase(
            "spectral-mobius-500",
            mobiusGraph,
            seed: FixtureGoldenHelpers.SpectralSeed,
            k: 8,
            denseOptions: default,
            materialization: DenseLaplacianMaterialization.Rectangular));
        spectralCases.Add(MeasureSpectralCase(
            "spectral-mobius-500",
            mobiusGraph,
            seed: FixtureGoldenHelpers.SpectralSeed,
            k: 8,
            denseOptions: default,
            materialization: DenseLaplacianMaterialization.FlatColumnMajor));

#if !EIGEN_REFERENCE
        spectralCases.Add(MeasureSpectralCase(
            "spectral-mobius-500",
            mobiusGraph,
            seed: FixtureGoldenHelpers.SpectralSeed,
            k: 8,
            denseOptions: new DenseEigenOptions(DenseEigenFastVariant.Fma),
            materialization: DenseLaplacianMaterialization.FlatColumnMajor));
#endif

        solverCases.Add(MeasureSpectralSolverComparison(
            "solver-crescent-500", crescentGraph, seed: FixtureGoldenHelpers.SpectralSeed, k: 8));
        solverCases.Add(MeasureSpectralSolverComparison(
            "solver-mobius-500", mobiusGraph, seed: FixtureGoldenHelpers.SpectralSeed, k: 8));

        directPath = run.WriteRunJson("direct-comparison", directCases);
        variantPath = run.WriteRunJson("fast-variant-comparison", variantCases);
        spectralPath = run.WriteRunJson("spectral-dispatch-comparison", spectralCases);
        solverPath = run.WriteRunJson("spectral-solver-comparison", solverCases);
        denseBottomKPath = run.WriteRunJson("dense-bottomk-comparison", denseBottomKCases);
        summaryPath = run.WriteRunText("summary", BuildSummary(directCases, variantCases, spectralCases, solverCases, denseBottomKCases));

        _output.WriteLine($"RunRoot\t{run.RunDirectory}");
        _output.WriteLine($"Manifest\t{run.ManifestPath}");
        _output.WriteLine($"DirectComparison\t{directPath}");
        _output.WriteLine($"FastVariantComparison\t{variantPath}");
        _output.WriteLine($"SpectralDispatchComparison\t{spectralPath}");
        _output.WriteLine($"SpectralSolverComparison\t{solverPath}");
        _output.WriteLine($"DenseBottomKComparison\t{denseBottomKPath}");
        _output.WriteLine($"Summary\t{summaryPath}");
#endif
    }

    private static string BuildSummary(
        IReadOnlyList<DirectBenchmarkCase> directCases,
        IReadOnlyList<FastVariantBenchmarkCase> variantCases,
        IReadOnlyList<SpectralBenchmarkCase> spectralCases,
        IReadOnlyList<SpectralSolverComparisonCase> solverCases,
        IReadOnlyList<DenseBottomKComparisonCase> denseBottomKCases)
    {
        var builder = new StringBuilder();
        builder.AppendLine("direct-comparison");
        builder.AppendLine("Label\tEigenMs\tEigenFastMs\tSpeedup\tResidual");
        foreach (DirectBenchmarkCase item in directCases)
            builder.AppendLine($"{item.Label}\t{item.EigenMilliseconds:F3}\t{item.EigenFastMilliseconds:F3}\t{item.Speedup:F2}x\t{item.Residual:E3}");

        if (variantCases.Count > 0)
        {
            builder.AppendLine();
            builder.AppendLine("fast-variant-comparison");
            builder.AppendLine("Label\tDefaultMs\tFmaMs\tSpeedup\tResidual");
            foreach (FastVariantBenchmarkCase item in variantCases)
                builder.AppendLine($"{item.Label}\t{item.DefaultMilliseconds:F3}\t{item.FmaMilliseconds:F3}\t{item.Speedup:F2}x\t{item.Residual:E3}");
        }

        builder.AppendLine();
        builder.AppendLine("spectral-dispatch-comparison");
        builder.AppendLine("Label\tMaterialization\tFastVariant\tComputeBottomKMs\tLowestLambda\tCount");
        foreach (SpectralBenchmarkCase item in spectralCases)
            builder.AppendLine($"{item.Label}\t{item.Materialization}\t{item.FastVariant}\t{item.ComputeBottomKMilliseconds:F3}\t{item.LowestLambda:G17}\t{item.Count}");

        if (solverCases.Count > 0)
        {
            builder.AppendLine();
            builder.AppendLine("spectral-solver-comparison (Dense vs LOBPCG)");
            builder.AppendLine("Label\tDenseMs\tLobpcgMs\tSpeedup\tMaxLambdaDelta\tMaxResidual\tConverged\tIters");
            foreach (SpectralSolverComparisonCase item in solverCases)
                builder.AppendLine($"{item.Label}\t{item.DenseMilliseconds:F3}\t{item.LobpcgMilliseconds:F3}\t{item.Speedup:F2}x\t{item.MaxLambdaDelta:E3}\t{item.MaxResidual:E3}\t{item.Converged}\t{item.Iterations}");
        }

        if (denseBottomKCases.Count > 0)
        {
            builder.AppendLine();
            builder.AppendLine("dense-bottomk-comparison (Direct vs LOBPCG)");
            builder.AppendLine("Label\tDirectMs\tLobpcgMs\tSpeedup\tMaxLambdaDelta\tMaxResidual");
            foreach (DenseBottomKComparisonCase item in denseBottomKCases)
                builder.AppendLine($"{item.Label}\t{item.DirectMilliseconds:F3}\t{item.LobpcgMilliseconds:F3}\t{item.Speedup:F2}x\t{item.MaxLambdaDelta:E3}\t{item.MaxResidual:E3}");
        }

        return builder.ToString();
    }

    // Correctness gate (runs in every configuration, unlike the Release-only
    // benchmark above): LOBPCG must recover the dense reference spectrum, and the
    // deflated path must skip the trivial null mode rather than collapse onto it.
    [Fact]
    public void LobpcgMatchesDenseBottomK()
    {
        double[][] crescent = CrescentAndEllipsoid.Generate(crescentPoints: 250, ellipsoidPoints: 0, seed: 42).Features;
        CsrGraph graph = FixtureGoldenHelpers.BuildDraftGraph(
            crescent, k: FixtureGoldenHelpers.ConnectedFixtureK, ensureConnected: true);

        const int k = 6;
        int seed = FixtureGoldenHelpers.SpectralSeed;
        double[,] laplacian = BuildCombinatorialLaplacian(graph);

        IReadOnlyList<EigenPair> dense = Spectral.ComputeBottomK(
            graph, seed: seed, k: k, lapType: LaplacianType.Combinatorial, solverKind: SolverKind.Dense);

        // (1) Trivial-inclusive LOBPCG must match the dense bottom-k eigenvalue for
        // eigenvalue, and each returned pair must be a genuine eigenpair of L.
        LOBPCG.Result iterative = GraphSpectral.ComputeBottomK(
            graph, k: k, lapType: LaplacianType.Combinatorial, seed: seed,
            maxIter: 500, deflateNullSpace: false);

        for (int i = 0; i < k; i++)
        {
            double delta = Math.Abs(dense[i].Lambda - iterative.Eigenpairs[i].Lambda);
            Assert.True(delta < 1e-3,
                $"bottom-k mismatch at {i}: dense={dense[i].Lambda:G6} lobpcg={iterative.Eigenpairs[i].Lambda:G6} (res {iterative.ResidualNorm:E2}, {iterative.Iterations} iters)");
        }

        double iterativeResidual = MaxResidualForPairs(laplacian, iterative.Eigenpairs);
        Assert.True(iterativeResidual < 1e-3, $"LOBPCG residual too large: {iterativeResidual:E3}");

        // (2) Regression guard for the deflation fix: with the null mode deflated,
        // the smallest returned eigenvalue must be the (non-zero) Fiedler value, not
        // the λ≈0 constant mode leaking back in through the Jacobi preconditioner.
        LOBPCG.Result deflated = GraphSpectral.ComputeBottomK(
            graph, k: k - 1, lapType: LaplacianType.Combinatorial, seed: seed,
            maxIter: 500, deflateNullSpace: true);

        double fiedler = dense[1].Lambda; // smallest non-trivial eigenvalue
        Assert.True(deflated.Eigenpairs[0].Lambda > 0.5 * fiedler,
            $"deflation failed: smallest mode {deflated.Eigenpairs[0].Lambda:E3} collapsed toward the null space (Fiedler {fiedler:E3})");

        for (int i = 0; i < k - 1; i++)
        {
            double delta = Math.Abs(dense[i + 1].Lambda - deflated.Eigenpairs[i].Lambda);
            Assert.True(delta < 1e-3,
                $"non-trivial mode mismatch at {i}: dense={dense[i + 1].Lambda:G6} deflated={deflated.Eigenpairs[i].Lambda:G6}");
        }
    }

    private SpectralSolverComparisonCase MeasureSpectralSolverComparison(
        string label, CsrGraph graph, int seed, int k)
    {
        _ = Spectral.ComputeBottomK(graph, seed: seed, k: k, lapType: LaplacianType.Combinatorial, solverKind: SolverKind.Dense);
        _ = Spectral.ComputeBottomK(graph, seed: seed, k: k, lapType: LaplacianType.Combinatorial, solverKind: SolverKind.Iterative);

        double denseMedian = MeasureMedianMilliseconds(() =>
        {
            IReadOnlyList<EigenPair> p = Spectral.ComputeBottomK(graph, seed: seed, k: k, lapType: LaplacianType.Combinatorial, solverKind: SolverKind.Dense);
            _spectralSink = p.Count;
        }, iterations: 5);

        double lobpcgMedian = MeasureMedianMilliseconds(() =>
        {
            IReadOnlyList<EigenPair> p = Spectral.ComputeBottomK(graph, seed: seed, k: k, lapType: LaplacianType.Combinatorial, solverKind: SolverKind.Iterative);
            _spectralSink = p.Count;
        }, iterations: 5);

        IReadOnlyList<EigenPair> dense = Spectral.ComputeBottomK(graph, seed: seed, k: k, lapType: LaplacianType.Combinatorial, solverKind: SolverKind.Dense);
        LOBPCG.Result lobpcg = GraphSpectral.ComputeBottomK(graph, k: k, lapType: LaplacianType.Combinatorial, seed: seed, deflateNullSpace: false);

        double maxLambdaDelta = MaxLambdaDelta(dense, lobpcg.Eigenpairs);
        double maxResidual = MaxResidualForPairs(BuildCombinatorialLaplacian(graph), lobpcg.Eigenpairs);
        double speedup = lobpcgMedian > 0.0 ? denseMedian / lobpcgMedian : double.NaN;

        return new SpectralSolverComparisonCase(
            label, denseMedian, lobpcgMedian, speedup, maxLambdaDelta, maxResidual, lobpcg.Converged, lobpcg.Iterations);
    }

    private DenseBottomKComparisonCase MeasureDenseBottomKComparison(string label, double[,] matrix, int k)
    {
        var options = new LOBPCG.Options { MaxIterations = 2000, Tolerance = 1e-11 };

        _ = SpectralMath.BottomK(matrix, k);
        _ = LOBPCG.BottomK(matrix, k, options);

        double directMedian = MeasureMedianMilliseconds(() =>
        {
            IReadOnlyList<EigenPair> p = SpectralMath.BottomK(matrix, k);
            _spectralSink = p.Count;
        }, iterations: 5);

        double lobpcgMedian = MeasureMedianMilliseconds(() =>
        {
            IReadOnlyList<EigenPair> p = LOBPCG.BottomK(matrix, k, options);
            _spectralSink = p.Count;
        }, iterations: 5);

        IReadOnlyList<EigenPair> direct = SpectralMath.BottomK(matrix, k);
        IReadOnlyList<EigenPair> lobpcg = LOBPCG.BottomK(matrix, k, options);

        double maxLambdaDelta = MaxLambdaDelta(direct, lobpcg);
        double maxResidual = MaxResidualForPairs(matrix, lobpcg);
        double speedup = lobpcgMedian > 0.0 ? directMedian / lobpcgMedian : double.NaN;

        return new DenseBottomKComparisonCase(label, directMedian, lobpcgMedian, speedup, maxLambdaDelta, maxResidual);
    }

    private static double MaxResidualForPairs(double[,] matrix, IReadOnlyList<EigenPair> pairs)
    {
        double max = 0.0;
        foreach (EigenPair pair in pairs)
        {
            double residual = ComputeResidualNorm(matrix, pair.Lambda, pair.Vector);
            if (residual > max) max = residual;
        }

        return max;
    }

    private static double MaxLambdaDelta(IReadOnlyList<EigenPair> a, IReadOnlyList<EigenPair> b)
    {
        int count = Math.Min(a.Count, b.Count);
        double max = 0.0;
        for (int i = 0; i < count; i++)
        {
            double delta = Math.Abs(a[i].Lambda - b[i].Lambda);
            if (delta > max) max = delta;
        }

        return max;
    }

    private DirectBenchmarkCase MeasureCase(string label, double[,] matrix)
    {
        WarmUp(matrix);

        double eigenMedian = MeasureMedianMilliseconds(() => Eigen.DecomposeSymmetric(matrix), iterations: 7);
        double fastMedian = MeasureMedianMilliseconds(() => EigenFast.DecomposeSymmetric(matrix), iterations: 7);
        EigenResult fast = EigenFast.DecomposeSymmetric(matrix);
        double residual = ComputeMaxResidualNorm(matrix, fast);
        double speedup = eigenMedian / fastMedian;

        return new DirectBenchmarkCase(label, eigenMedian, fastMedian, speedup, residual);
    }

    private FastVariantBenchmarkCase MeasureFastVariantCase(string label, double[,] matrix)
    {
        WarmUp(matrix);
        _ = EigenFast.DecomposeSymmetric(matrix, fastVariant: DenseEigenFastVariant.Fma);

        double defaultMedian = MeasureMedianMilliseconds(
            () => EigenFast.DecomposeSymmetric(matrix),
            iterations: 7);
        double fmaMedian = MeasureMedianMilliseconds(
            () => EigenFast.DecomposeSymmetric(matrix, fastVariant: DenseEigenFastVariant.Fma),
            iterations: 7);
        EigenResult fma = EigenFast.DecomposeSymmetric(matrix, fastVariant: DenseEigenFastVariant.Fma);
        double residual = ComputeMaxResidualNorm(matrix, fma);
        double speedup = defaultMedian / fmaMedian;

        return new FastVariantBenchmarkCase(label, defaultMedian, fmaMedian, speedup, residual);
    }

    private SpectralBenchmarkCase MeasureSpectralCase(
        string label,
        CsrGraph graph,
        int seed,
        int k,
        DenseEigenOptions denseOptions,
        DenseLaplacianMaterialization materialization)
    {
        WarmUpSpectral(graph, seed, k, denseOptions, materialization);

        double median = MeasureMedianMilliseconds(() =>
        {
            IReadOnlyList<EigenPair> pairs = Spectral.ComputeBottomK(
                graph,
                seed: seed,
                k: k,
                lapType: LaplacianType.Combinatorial,
                solverKind: SolverKind.Dense,
                denseOptions: denseOptions,
                denseMaterialization: materialization);
            _spectralSink = pairs.Count;
        }, iterations: 5);

        IReadOnlyList<EigenPair> result = Spectral.ComputeBottomK(
            graph,
            seed: seed,
            k: k,
            lapType: LaplacianType.Combinatorial,
            solverKind: SolverKind.Dense,
            denseOptions: denseOptions,
            denseMaterialization: materialization);

        double lowestLambda = result.Count > 0 ? result[0].Lambda : double.NaN;
        return new SpectralBenchmarkCase(
            label,
            materialization.ToString(),
            denseOptions.FastVariant.ToString(),
            median,
            lowestLambda,
            result.Count);
    }

    private static void WarmUp(double[,] matrix)
    {
        _ = Eigen.DecomposeSymmetric(matrix);
        _ = EigenFast.DecomposeSymmetric(matrix);
    }

    private static void WarmUpSpectral(
        CsrGraph graph,
        int seed,
        int k,
        DenseEigenOptions denseOptions,
        DenseLaplacianMaterialization materialization)
    {
        IReadOnlyList<EigenPair> pairs = Spectral.ComputeBottomK(
            graph,
            seed: seed,
            k: k,
            lapType: LaplacianType.Combinatorial,
            solverKind: SolverKind.Dense,
            denseOptions: denseOptions,
            denseMaterialization: materialization);
        _spectralSink = pairs.Count;
    }

    private static string GetDenseEigenBackendLabel()
    {
#if EIGEN_REFERENCE
        return "Eigen";
#else
        return "EigenFast";
#endif
    }

    private static double MeasureMedianMilliseconds(Action action, int iterations)
    {
        var samples = new List<double>(iterations);

        for (int i = 0; i < iterations; i++)
        {
            var stopwatch = Stopwatch.StartNew();
            action();
            stopwatch.Stop();
            samples.Add(stopwatch.Elapsed.TotalMilliseconds);
        }

        samples.Sort();
        return samples[samples.Count / 2];
    }

    private static double ComputeMaxResidualNorm(double[,] matrix, EigenResult result)
    {
        double max = 0.0;

        for (int i = 0; i < result.Eigenvalues.Length; i++)
        {
            double residual = ComputeResidualNorm(matrix, result.Eigenvalues[i], result.Eigenvectors[i]);
            if (residual > max)
                max = residual;
        }

        return max;
    }

    private static double ComputeResidualNorm(double[,] matrix, double eigenvalue, double[] eigenvector)
    {
        int n = eigenvector.Length;
        double sumSquares = 0.0;

        for (int row = 0; row < n; row++)
        {
            double projected = 0.0;
            for (int col = 0; col < n; col++)
                projected += matrix[row, col] * eigenvector[col];

            double residual = projected - eigenvalue * eigenvector[row];
            sumSquares += residual * residual;
        }

        return Math.Sqrt(sumSquares);
    }

    private static double[,] BuildRandomSymmetricMatrix(int size, int seed)
    {
        var random = new Random(seed);
        var matrix = new double[size, size];

        for (int i = 0; i < size; i++)
        {
            for (int j = i; j < size; j++)
            {
                double value = random.NextDouble() * 2.0 - 1.0;
                matrix[i, j] = value;
                matrix[j, i] = value;
            }
        }

        return matrix;
    }

    private static double[,] BuildCombinatorialLaplacian(CsrGraph graph)
    {
        int n = graph.NodeCount;
        var laplacian = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            int start = graph.RowPointers[i];
            int end = graph.RowPointers[i + 1];
            double degree = 0.0;

            for (int edge = start; edge < end; edge++)
            {
                degree += graph.Weights[edge];
                laplacian[i, graph.Targets[edge]] -= graph.Weights[edge];
            }

            laplacian[i, i] = degree;
        }

        return laplacian;
    }

    private readonly record struct DirectBenchmarkCase(
        string Label,
        double EigenMilliseconds,
        double EigenFastMilliseconds,
        double Speedup,
        double Residual);

    private readonly record struct FastVariantBenchmarkCase(
        string Label,
        double DefaultMilliseconds,
        double FmaMilliseconds,
        double Speedup,
        double Residual);

    private readonly record struct SpectralBenchmarkCase(
        string Label,
        string Materialization,
        string FastVariant,
        double ComputeBottomKMilliseconds,
        double LowestLambda,
        int Count);

    private readonly record struct SpectralSolverComparisonCase(
        string Label,
        double DenseMilliseconds,
        double LobpcgMilliseconds,
        double Speedup,
        double MaxLambdaDelta,
        double MaxResidual,
        bool Converged,
        int Iterations);

    private readonly record struct DenseBottomKComparisonCase(
        string Label,
        double DirectMilliseconds,
        double LobpcgMilliseconds,
        double Speedup,
        double MaxLambdaDelta,
        double MaxResidual);
}
