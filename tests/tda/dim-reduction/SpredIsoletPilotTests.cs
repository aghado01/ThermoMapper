using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using Graphs;
using Maths.Geometry;
using TDA.Ph;
using Xunit;
using Xunit.Abstractions;

namespace TDA.DimReduction.Tests;

/// <summary>
/// ISOLET Phase-2 pilot (S0, dimension 30, two seeds) — the standing pilot requirement of the
/// ISOLET brief's "H0 matching-cost gate": expose pathology and cost, and produce the full-screen
/// wall-clock extrapolation before any promotion decision. Not an assertion test — a recorded run
/// that logs and writes the immutable JSON artifact to <c>issues/spred/pilot/</c>.
/// <para>Protocol per the brief: raw617 (no preprocessing), deterministic label-free row shuffle
/// with recorded seed and permutation hash, 8 contiguous blocks after the shuffle, H0-only
/// objective with the adopted sliced-Wasserstein screening metric, parity-grade graph recipe
/// (mutual kNN K=10 + MST repair) on both ambient reference and projected clouds.</para>
/// </summary>
public sealed class SpredIsoletPilotTests
{
    private readonly ITestOutputHelper _out;
    public SpredIsoletPilotTests(ITestOutputHelper output) => _out = output;

    private const int TargetDim = 30;
    private const int BlockCount = 8;
    private const int PilotIters = 100;
    private const int ShuffleSeed = 41;
    private static readonly int[] PilotSeeds = [211, 223];

    [Fact]
    [Trait("Category", "Benchmark")]   // manual pilot; run alone with --filter "Category=Benchmark&FullyQualifiedName~SpredIsoletPilot"
    public void Pilot_S0_Dim30_TwoSeeds_WallClockExtrapolation()
    {
        string repoRoot = LocateRepoRoot();
        string gzPath = Path.Combine(repoRoot, "datasets", "isolet.csv.gz");
        (double[][] features, string datasetHash) = LoadIsoletFeatures(gzPath);
        int n = features.Length, d = features[0].Length;
        _out.WriteLine($"dataset: {n} x {d}, sha256 {datasetHash[..16]}…");

        // Deterministic label-free shuffle (brief: ISOLET source order is speaker-organized, so
        // contiguous blocks would encode cohorts). SplitMix64-driven Fisher-Yates — the algorithm
        // is fully specified here so seed + hash make the permutation reproducible and invertible
        // without storing 7,797 indices.
        int[] permutation = SplitMix64Permutation(n, ShuffleSeed);
        string permutationHash = Sha256Hex(string.Join(",", permutation));
        var shuffled = new double[n][];
        for (int i = 0; i < n; i++) shuffled[i] = features[permutation[i]];
        _out.WriteLine($"shuffle: seed {ShuffleSeed}, permutation sha256 {permutationHash[..16]}…");

        PersistenceObjectiveConfig config = PilotConfig();
        int parallelism = Math.Min(BlockCount, Environment.ProcessorCount);

        // Fixed overhead: maxIters = 0 returns the PCA warm start per block, so this run prices
        // everything except annealing iterations (block references, PCA, diagnostics evals, and
        // the full-data objective).
        var sw = Stopwatch.StartNew();
        DistributedSpredResult overhead = DistributedSpred.ComputeWithDiagnostics(
            shuffled, TargetDim, BlockCount, config, maxIters: 0, seed: PilotSeeds[0], parallelism);
        double fixedSeconds = sw.Elapsed.TotalSeconds;
        _out.WriteLine($"fixed overhead (maxIters=0): {fixedSeconds:F1} s, fullDataObjective {overhead.FullDataObjective:G6}");

        // The two pilot seeds.
        var results = new DistributedSpredResult[PilotSeeds.Length];
        var seedSeconds = new double[PilotSeeds.Length];
        for (int s = 0; s < PilotSeeds.Length; s++)
        {
            sw.Restart();
            results[s] = DistributedSpred.ComputeWithDiagnostics(
                shuffled, TargetDim, BlockCount, config, PilotIters, PilotSeeds[s], parallelism);
            seedSeconds[s] = sw.Elapsed.TotalSeconds;
            _out.WriteLine($"seed {PilotSeeds[s]}: {seedSeconds[s]:F1} s, fullDataObjective {results[s].FullDataObjective:G6}");
        }

        // Health diagnostics: local vs aggregate objectives and Grassmann angles.
        var grass = new GrassmannManifold(d, TargetDim);
        var seedReports = new List<object>();
        for (int s = 0; s < PilotSeeds.Length; s++)
        {
            DistributedSpredResult r = results[s];
            var blocks = new List<object>();
            _out.WriteLine($"seed {PilotSeeds[s]} blocks (local / aggregate objective / angle-to-aggregate):");
            foreach (DistributedSpredBlockResult b in r.Blocks)
            {
                double angle = grass.Distance(PackFrame(b.Projection, d), PackFrame(r.Projection, d));
                _out.WriteLine($"  block {b.Index}: {b.LocalObjective:G6} / {b.AggregateObjective:G6} / {angle:F4}");
                blocks.Add(new
                {
                    b.Index, b.Start, b.Count, b.Seed,
                    localObjective = b.LocalObjective,
                    aggregateObjective = b.AggregateObjective,
                    grassmannAngleToAggregate = angle,
                });
            }
            seedReports.Add(new
            {
                seed = PilotSeeds[s],
                wallClockSeconds = seedSeconds[s],
                fullDataObjective = r.FullDataObjective,
                blocks,
            });
        }
        double seedToSeedAngle = grass.Distance(
            PackFrame(results[0].Projection, d), PackFrame(results[1].Projection, d));
        _out.WriteLine($"aggregate-to-aggregate Grassmann distance across seeds: {seedToSeedAngle:F4}");

        // Wall-clock extrapolation. Marginal per-iteration cost is the annealing share of a run;
        // the fixed share (references, diagnostics, full-data objective) recurs per run. Target
        // dimension moves only the projected-kNN share, so per-run cost is treated as
        // dimension-independent across {20, 30, 50}.
        double meanSeed = 0.0;
        foreach (double t in seedSeconds) meanSeed += t;
        meanSeed /= seedSeconds.Length;
        double perIter = Math.Max(meanSeed - fixedSeconds, 0.0) / PilotIters;

        _out.WriteLine("");
        _out.WriteLine($"extrapolation: run(I) ≈ {fixedSeconds:F1} s + I · {perIter:F3} s   (8 blocks, parallelism {parallelism})");
        var screens = new List<object>();
        _out.WriteLine("full S0 screen (3 dims x 5 seeds) at iteration budgets:");
        foreach (int iters in new[] { 100, 500, 1000 })
        {
            double perRun = fixedSeconds + iters * perIter;
            double screen = perRun * 3 * 5;
            _out.WriteLine($"  I={iters,5}: {perRun,7:F1} s/run  ->  {screen / 60,7:F1} min for S0");
            screens.Add(new { iterations = iters, secondsPerRun = perRun, s0ScreenMinutes = screen / 60 });
        }

        WriteArtifact(repoRoot, new
        {
            pilot = "S0-dim30-two-seeds",
            date = DateTime.UtcNow.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture),
            dataset = new { path = "datasets/isolet.csv.gz", rows = n, features = d, sha256 = datasetHash },
            shuffle = new
            {
                algorithm = "SplitMix64-FisherYates (Lemire bounded)",
                seed = ShuffleSeed,
                permutationSha256 = permutationHash,
            },
            config = new
            {
                objective = "H0",
                dimensions = new[] { new { dim = 0, weight = 1.0 } },
                maxDimension = 1,
                diagramDistance = "SlicedWasserstein",
                slicedDirections = 50,
                wassersteinOrder = 2.0,
                minPersistence = 0.0,
                graph = new { topology = "Knn", k = 10, filter = "MutualKnn", repair = "MstMin" },
                targetDim = TargetDim,
                blockCount = BlockCount,
                pilotIterations = PilotIters,
                maxDegreeOfParallelism = parallelism,
                processorCount = Environment.ProcessorCount,
            },
            fixedOverheadSeconds = fixedSeconds,
            seeds = seedReports,
            aggregateToAggregateGrassmann = seedToSeedAngle,
            extrapolation = new
            {
                model = "run(I) = fixedOverheadSeconds + I * perIterationSeconds",
                perIterationSeconds = perIter,
                s0Screen = screens,
            },
        });

        Assert.True(true);
    }

    private static PersistenceObjectiveConfig PilotConfig() => new()
    {
        Graph = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 10 },
            Filter = new FilterConfig { Kind = FilterKind.MutualKnn },
            Repair = new RepairConfig { Kind = RepairKind.MstMin },
            Projection = new DistanceProjection(),
        },
        Dimensions = [(0, 1.0)],
        MaxDimension = 1,
        DiagramDistance = DiagramDistanceKind.SlicedWasserstein,
    };

    // SplitMix64 (Vigna) driving a Fisher-Yates shuffle with Lemire's multiply-shift bounding —
    // fully specified so the permutation is reproducible from the seed alone, cross-platform.
    private static int[] SplitMix64Permutation(int n, int seed)
    {
        ulong x = unchecked((ulong)seed);
        ulong Next()
        {
            x += 0x9e3779b97f4a7c15UL;
            ulong z = x;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9UL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebUL;
            return z ^ (z >> 31);
        }

        var perm = new int[n];
        for (int i = 0; i < n; i++) perm[i] = i;
        for (int i = n - 1; i > 0; i--)
        {
            int j = (int)(ulong)(((UInt128)Next() * (ulong)(i + 1)) >> 64);
            (perm[i], perm[j]) = (perm[j], perm[i]);
        }
        return perm;
    }

    private static (double[][] Features, string Sha256) LoadIsoletFeatures(string gzPath)
    {
        using var fs = File.OpenRead(gzPath);
        using var gz = new GZipStream(fs, CompressionMode.Decompress);
        using var ms = new MemoryStream();
        gz.CopyTo(ms);
        byte[] bytes = ms.ToArray();
        string hash = Convert.ToHexString(SHA256.HashData(bytes)).ToLowerInvariant();

        // Label column (last) is split off before any graph/SPRED API sees the data.
        var feats = new List<double[]>(8000);
        using var reader = new StreamReader(new MemoryStream(bytes));
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (string.IsNullOrWhiteSpace(line)) continue;
            string[] parts = line.Split(',');
            var row = new double[parts.Length - 1];
            for (int j = 0; j < row.Length; j++)
                row[j] = double.Parse(parts[j], CultureInfo.InvariantCulture);
            feats.Add(row);
        }
        return (feats.ToArray(), hash);
    }

    private static string Sha256Hex(string text) =>
        Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(text))).ToLowerInvariant();

    private static double[] PackFrame(double[][] projection, int ambientDim)
    {
        var frame = new double[ambientDim * projection.Length];
        for (int col = 0; col < projection.Length; col++)
            Array.Copy(projection[col], 0, frame, col * ambientDim, ambientDim);
        return frame;
    }

    private static void WriteArtifact(string repoRoot, object artifact)
    {
        string dir = Path.Combine(repoRoot, "issues", "spred", "pilot");
        Directory.CreateDirectory(dir);
        string path = Path.Combine(dir, "spred-pilot-s0-dim30.json");
        File.WriteAllText(path, JsonSerializer.Serialize(
            artifact, new JsonSerializerOptions { WriteIndented = true }));
    }

    private static string LocateRepoRoot()
    {
        for (DirectoryInfo? dir = new(AppContext.BaseDirectory); dir is not null; dir = dir.Parent)
        {
            if (File.Exists(Path.Combine(dir.FullName, "datasets", "isolet.csv.gz")))
                return dir.FullName;
        }
        throw new DirectoryNotFoundException("Could not locate the repo root (datasets/isolet.csv.gz) above " + AppContext.BaseDirectory);
    }
}
