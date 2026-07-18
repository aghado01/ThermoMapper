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
using Maths.Geometry.DimReduction;
using Maths.LinAlg;
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

    /// <summary>
    /// Iteration-budget probe (pilot follow-up): one seed at I = 1000 against the PCA warm-start
    /// baseline, to separate "the anneal needs more budget" from "the proposal/cooling scale is
    /// wrong at 617 -> 30". The pilot found I = 100 leaves five of eight blocks bit-identical to
    /// the warm start; if I = 1000 moves them roughly linearly the screen just needs budget, and
    /// if they stay pinned the annealer needs tuning first.
    /// </summary>
    [Fact]
    [Trait("Category", "Benchmark")]
    public void Probe_S0_Dim30_IterationBudget()
    {
        const int probeIters = 1000;
        int probeSeed = PilotSeeds[0];

        string repoRoot = LocateRepoRoot();
        (double[][] features, string datasetHash) = LoadIsoletFeatures(
            Path.Combine(repoRoot, "datasets", "isolet.csv.gz"));
        int n = features.Length, d = features[0].Length;
        int[] permutation = SplitMix64Permutation(n, ShuffleSeed);
        var shuffled = new double[n][];
        for (int i = 0; i < n; i++) shuffled[i] = features[permutation[i]];

        PersistenceObjectiveConfig config = PilotConfig();
        int parallelism = Math.Min(BlockCount, Environment.ProcessorCount);

        var sw = Stopwatch.StartNew();
        DistributedSpredResult warm = DistributedSpred.ComputeWithDiagnostics(
            shuffled, TargetDim, BlockCount, config, maxIters: 0, probeSeed, parallelism);
        double warmSeconds = sw.Elapsed.TotalSeconds;

        sw.Restart();
        DistributedSpredResult annealed = DistributedSpred.ComputeWithDiagnostics(
            shuffled, TargetDim, BlockCount, config, probeIters, probeSeed, parallelism);
        double annealedSeconds = sw.Elapsed.TotalSeconds;

        _out.WriteLine($"warm start (maxIters=0): {warmSeconds:F1} s, fullDataObjective {warm.FullDataObjective:G6}");
        _out.WriteLine($"I={probeIters}, seed {probeSeed}: {annealedSeconds:F1} s, fullDataObjective {annealed.FullDataObjective:G6}");
        _out.WriteLine($"extrapolation check: pilot model predicted {46.7 + probeIters * 0.025:F1} s");
        _out.WriteLine("");
        _out.WriteLine("block | warm-start local | annealed local | improvement %");

        var grass = new GrassmannManifold(d, TargetDim);
        var blocks = new List<object>();
        int moved = 0;
        for (int b = 0; b < BlockCount; b++)
        {
            double baseline = warm.Blocks[b].LocalObjective;
            double after = annealed.Blocks[b].LocalObjective;
            double improvementPct = 100.0 * (baseline - after) / baseline;
            double angleFromWarm = grass.Distance(
                PackFrame(warm.Blocks[b].Projection, d), PackFrame(annealed.Blocks[b].Projection, d));
            if (after < baseline) moved++;
            _out.WriteLine($"  {b}   | {baseline,14:G6}  | {after,12:G6}  | {improvementPct,6:F3}   (Grassmann from warm start: {angleFromWarm:F4})");
            blocks.Add(new
            {
                index = b,
                warmStartLocal = baseline,
                annealedLocal = after,
                improvementPercent = improvementPct,
                grassmannFromWarmStart = angleFromWarm,
            });
        }
        _out.WriteLine("");
        _out.WriteLine($"blocks improved: {moved}/{BlockCount}");

        File.WriteAllText(
            Path.Combine(repoRoot, "issues", "spred", "pilot", "spred-probe-s0-dim30-i1000.json"),
            JsonSerializer.Serialize(new
            {
                probe = "S0-dim30-iteration-budget",
                date = DateTime.UtcNow.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture),
                seed = probeSeed,
                iterations = probeIters,
                shuffleSeed = ShuffleSeed,
                datasetSha256 = datasetHash,
                warmStartSeconds = warmSeconds,
                annealedSeconds,
                warmStartFullDataObjective = warm.FullDataObjective,
                annealedFullDataObjective = annealed.FullDataObjective,
                blocksImproved = moved,
                blocks,
            }, new JsonSerializerOptions { WriteIndented = true }));

        Assert.True(true);
    }

    /// <summary>
    /// Mobility re-probe (dev-sequence P3; mobility-brief findings 1+2): the I = 1000 budget probe
    /// re-run through the two-plane annealer, with <c>InitialTemperature</c> calibrated to the
    /// recorded per-move objective increment — the budget probe's single accepted move improved
    /// block 5 by ~4e-3, so temperature 1e-3 puts a comparable worsening at exp(−4) ≈ 2%
    /// acceptance: refinement, not melt (the same increment-commensurate choice as the engine
    /// mobility fact). Runs the addendum's eigengap pre-check first: per-block PCA eigengaps
    /// λ25…λ35 around the k = 30 cut, read before spending anneal budget — small relative gaps
    /// mean the offset directions smear across near-degenerate retained columns and single-column
    /// Givens moves crawl (finding 1), which is how a still-pinned re-probe should then be read.
    /// </summary>
    [Fact]
    [Trait("Category", "Benchmark")]
    public void Probe_S0_Dim30_MobilityReprobe()
    {
        const int probeIters = 1000;
        int probeSeed = PilotSeeds[0];
        var annealerOptions = new SubspaceAnnealerOptions { InitialTemperature = 1e-3 };

        string repoRoot = LocateRepoRoot();
        (double[][] features, string datasetHash) = LoadIsoletFeatures(
            Path.Combine(repoRoot, "datasets", "isolet.csv.gz"));
        int n = features.Length, d = features[0].Length;
        int[] permutation = SplitMix64Permutation(n, ShuffleSeed);
        var shuffled = new double[n][];
        for (int i = 0; i < n; i++) shuffled[i] = features[permutation[i]];

        // Eigengap pre-check (finding 1's interpretation gate). Same contiguous slices the
        // distributed run uses; eigenvalue indices below are 1-based in the λ naming.
        _out.WriteLine("eigengap pre-check, per-block PCA λ25..λ35 (relative gap (λi−λi+1)/λi):");
        var eigenReports = new List<object>();
        for (int b = 0; b < BlockCount; b++)
        {
            int start = b * n / BlockCount, end = (b + 1) * n / BlockCount;
            var slice = new double[end - start][];
            Array.Copy(shuffled, start, slice, 0, slice.Length);
            double[] eig = Pca.Compute(slice, numComponents: 35).Eigenvalues;

            var window = new double[11];
            Array.Copy(eig, 24, window, 0, 11);           // λ25..λ35
            var relGaps = new double[10];
            double minGap = double.PositiveInfinity;
            for (int i = 0; i < 10; i++)
            {
                relGaps[i] = (window[i] - window[i + 1]) / window[i];
                if (relGaps[i] < minGap) minGap = relGaps[i];
            }
            double cutGap = (eig[29] - eig[30]) / eig[29];   // the k = 30 boundary: λ30 vs λ31
            _out.WriteLine($"  block {b}: cut gap {cutGap:P2}, min gap in window {minGap:P2}");
            eigenReports.Add(new
            {
                block = b,
                lambda25to35 = window,
                relativeGaps = relGaps,
                cutRelativeGap = cutGap,
                minRelativeGapInWindow = minGap,
            });
        }

        PersistenceObjectiveConfig config = PilotConfig();
        int parallelism = Math.Min(BlockCount, Environment.ProcessorCount);

        var sw = Stopwatch.StartNew();
        DistributedSpredResult warm = DistributedSpred.ComputeWithDiagnostics(
            shuffled, TargetDim, BlockCount, config, maxIters: 0, probeSeed, parallelism);
        double warmSeconds = sw.Elapsed.TotalSeconds;

        sw.Restart();
        DistributedSpredResult annealed = DistributedSpred.ComputeWithDiagnostics(
            shuffled, TargetDim, BlockCount, config, probeIters, probeSeed, parallelism,
            annealerOptions);
        double annealedSeconds = sw.Elapsed.TotalSeconds;

        _out.WriteLine("");
        _out.WriteLine($"warm start (maxIters=0): {warmSeconds:F1} s, fullDataObjective {warm.FullDataObjective:G6}");
        _out.WriteLine($"I={probeIters}, seed {probeSeed}, two-plane @ T0=1e-3: {annealedSeconds:F1} s, " +
            $"fullDataObjective {annealed.FullDataObjective:G6}");
        _out.WriteLine("reference: old-annealer I=1000 probe moved 1/8 blocks, best improvement 0.006% " +
            "(spred-probe-s0-dim30-i1000.json)");
        _out.WriteLine("");
        _out.WriteLine("block | warm-start local | annealed local | improvement % | Grassmann from warm");

        var grass = new GrassmannManifold(d, TargetDim);
        var blocks = new List<object>();
        int moved = 0;
        for (int b = 0; b < BlockCount; b++)
        {
            double baseline = warm.Blocks[b].LocalObjective;
            double after = annealed.Blocks[b].LocalObjective;
            double improvementPct = 100.0 * (baseline - after) / baseline;
            double angleFromWarm = grass.Distance(
                PackFrame(warm.Blocks[b].Projection, d), PackFrame(annealed.Blocks[b].Projection, d));
            if (after < baseline) moved++;
            _out.WriteLine($"  {b}   | {baseline,14:G6}  | {after,12:G6}  | {improvementPct,8:F4}  | {angleFromWarm:F4}");
            blocks.Add(new
            {
                index = b,
                warmStartLocal = baseline,
                annealedLocal = after,
                improvementPercent = improvementPct,
                grassmannFromWarmStart = angleFromWarm,
            });
        }
        _out.WriteLine("");
        _out.WriteLine($"blocks improved: {moved}/{BlockCount}");

        File.WriteAllText(
            Path.Combine(repoRoot, "issues", "spred", "pilot", "spred-probe-s0-dim30-mobility-reprobe.json"),
            JsonSerializer.Serialize(new
            {
                probe = "S0-dim30-mobility-reprobe",
                date = DateTime.UtcNow.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture),
                seed = probeSeed,
                iterations = probeIters,
                shuffleSeed = ShuffleSeed,
                datasetSha256 = datasetHash,
                annealer = new
                {
                    proposals = "two-plane Givens primary (shipped defaults)",
                    initialTemperature = annealerOptions.InitialTemperature,
                    temperatureCalibration =
                        "budget probe's recorded per-move increment ~4e-3 => T0 = 1e-3 (exp(-4) accept-worse)",
                    blockSeedDerivation = "SeedTree SplitMix64 children (post seed-aliasing audit)",
                },
                eigengapPreCheck = eigenReports,
                warmStartSeconds = warmSeconds,
                annealedSeconds,
                warmStartFullDataObjective = warm.FullDataObjective,
                annealedFullDataObjective = annealed.FullDataObjective,
                blocksImproved = moved,
                blocks,
            }, new JsonSerializerOptions { WriteIndented = true }));

        Assert.True(true);
    }

    /// <summary>
    /// Paired-move S0 probe (the mobility re-probe's successor): the eigengap pre-check confirmed
    /// a flat tail in every block (min λ25…λ35 gaps 0.65–1.33%), firing finding 1's gate — the
    /// two-plane re-probe's 0.05–0.1%/1000-iter descent is single-column crawl against smeared
    /// defects. This probe adds the paired two-column move at 50% mixture. StepFloor is raised to
    /// 0.05 because the paired kind's φ-window acceptance sits far below the 25% controller
    /// target, so its scale pins to the floor — the default 1e-3 floor would put excisions at
    /// diffusion scale. Reference numbers: two-plane re-probe moved 8/8 blocks by 0.048–0.096%
    /// (spred-probe-s0-dim30-mobility-reprobe.json); old annealer moved 1/8 by 0.006%.
    /// </summary>
    [Fact]
    [Trait("Category", "Benchmark")]
    public void Probe_S0_Dim30_PairedMoves()
    {
        const int probeIters = 1000;
        int probeSeed = PilotSeeds[0];
        var annealerOptions = new SubspaceAnnealerOptions
        {
            PairedFraction = 0.5,
            InitialTemperature = 1e-3,
            StepFloor = 0.05,
        };

        string repoRoot = LocateRepoRoot();
        (double[][] features, string datasetHash) = LoadIsoletFeatures(
            Path.Combine(repoRoot, "datasets", "isolet.csv.gz"));
        int n = features.Length, d = features[0].Length;
        int[] permutation = SplitMix64Permutation(n, ShuffleSeed);
        var shuffled = new double[n][];
        for (int i = 0; i < n; i++) shuffled[i] = features[permutation[i]];

        PersistenceObjectiveConfig config = PilotConfig();
        int parallelism = Math.Min(BlockCount, Environment.ProcessorCount);

        var sw = Stopwatch.StartNew();
        DistributedSpredResult warm = DistributedSpred.ComputeWithDiagnostics(
            shuffled, TargetDim, BlockCount, config, maxIters: 0, probeSeed, parallelism);
        double warmSeconds = sw.Elapsed.TotalSeconds;

        sw.Restart();
        DistributedSpredResult annealed = DistributedSpred.ComputeWithDiagnostics(
            shuffled, TargetDim, BlockCount, config, probeIters, probeSeed, parallelism,
            annealerOptions);
        double annealedSeconds = sw.Elapsed.TotalSeconds;

        _out.WriteLine($"warm start (maxIters=0): {warmSeconds:F1} s, fullDataObjective {warm.FullDataObjective:G6}");
        _out.WriteLine($"I={probeIters}, seed {probeSeed}, paired 0.5 @ T0=1e-3, floor 0.05: {annealedSeconds:F1} s, " +
            $"fullDataObjective {annealed.FullDataObjective:G6}");
        _out.WriteLine("reference: two-plane re-probe (same protocol, no pairing) improved 8/8 blocks by 0.048-0.096%");
        _out.WriteLine("");
        _out.WriteLine("block | warm-start local | annealed local | improvement % | Grassmann from warm");

        var grass = new GrassmannManifold(d, TargetDim);
        var blocks = new List<object>();
        int moved = 0;
        for (int b = 0; b < BlockCount; b++)
        {
            double baseline = warm.Blocks[b].LocalObjective;
            double after = annealed.Blocks[b].LocalObjective;
            double improvementPct = 100.0 * (baseline - after) / baseline;
            double angleFromWarm = grass.Distance(
                PackFrame(warm.Blocks[b].Projection, d), PackFrame(annealed.Blocks[b].Projection, d));
            if (after < baseline) moved++;
            _out.WriteLine($"  {b}   | {baseline,14:G6}  | {after,12:G6}  | {improvementPct,8:F4}  | {angleFromWarm:F4}");
            blocks.Add(new
            {
                index = b,
                warmStartLocal = baseline,
                annealedLocal = after,
                improvementPercent = improvementPct,
                grassmannFromWarmStart = angleFromWarm,
            });
        }
        _out.WriteLine("");
        _out.WriteLine($"blocks improved: {moved}/{BlockCount}");

        File.WriteAllText(
            Path.Combine(repoRoot, "issues", "spred", "pilot", "spred-probe-s0-dim30-paired.json"),
            JsonSerializer.Serialize(new
            {
                probe = "S0-dim30-paired-moves",
                date = DateTime.UtcNow.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture),
                seed = probeSeed,
                iterations = probeIters,
                shuffleSeed = ShuffleSeed,
                datasetSha256 = datasetHash,
                annealer = new
                {
                    proposals = "two-plane Givens + paired two-column at 0.5 mixture",
                    pairedFraction = annealerOptions.PairedFraction,
                    initialTemperature = annealerOptions.InitialTemperature,
                    stepFloor = annealerOptions.StepFloor,
                    stepFloorRationale =
                        "paired-kind phi-window acceptance sits far below the controller target; the default floor would pin excisions at diffusion scale",
                    blockSeedDerivation = "SeedTree SplitMix64 children (post seed-aliasing audit)",
                },
                warmStartSeconds = warmSeconds,
                annealedSeconds,
                warmStartFullDataObjective = warm.FullDataObjective,
                annealedFullDataObjective = annealed.FullDataObjective,
                blocksImproved = moved,
                blocks,
            }, new JsonSerializerOptions { WriteIndented = true }));

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
