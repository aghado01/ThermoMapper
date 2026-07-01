// ============================================================================
// Clustering.Graphical.HdbScanSmoke — Program.cs
// ============================================================================
// Three-scenario smoke for the HDBSCAN pipeline:
//   1. 4 well-separated SpatialBlobs via HdbscanSession           — checks the
//      algorithm's basic cluster recovery on a known ground truth.
//   2. Single dense blob + 3 outliers, allowSingleCluster=true    — checks the
//      mapper-friendly degenerate case (one cluster + low-prob outliers).
//   3. Same 4 blobs via HdbscanClusterer (IClusterer adapter)     — checks the
//      adapter path used by Mapper (which itself routes through HdbscanSession).
//
// Exit code = number of scenarios that failed expectation. Stdout-only, no
// artifacts. Run via: dotnet run --project Clustering.Graphical.HdbScanSmoke.csproj
// ============================================================================

using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.HdbScan;
using Synthetic.Euclidean;
using TDA.Mapper;
using TDA.Mapper.Clusterers;

namespace Clustering.Graphical.HdbScanSmoke;

internal static class Program
{
    private static int Main()
    {
        Console.WriteLine("HDBSCAN smoke driver");
        Console.WriteLine("====================\n");

        int failures = 0;
        failures += Scenario_FourBlobs();
        failures += Scenario_SingleBlobPlusOutliers();
        failures += Scenario_ClustererAdapter();

        Console.WriteLine();
        Console.WriteLine(failures == 0
            ? "All scenarios passed."
            : $"FAILED — {failures} scenario(s) did not meet expectations.");
        return failures;
    }

    // ── Scenario 1 ─────────────────────────────────────────────────────────
    private static int Scenario_FourBlobs()
    {
        Console.WriteLine("[1] 4 blobs × 50pts in 2D, via HdbscanSession");

        var ds = SpatialBlobs.Generate(
            clusterCount: 4, pointsPerCluster: 50, dimensions: 2,
            separation: 10.0, spread: 0.3, seed: 42);

        ReportGeometry(ds);

        var result = HdbscanSession.Run(ds.Features, new HdbscanSettings
        {
            MinPts             = 5,
            AllowSingleCluster = true,
        });

        int    noise  = result.Labels.Count(L => L < 0);
        double purity = PerClusterPurity(result.Labels, ds.Labels, result.ClusterCount);
        int[] sizes = Enumerable.Range(0, result.ClusterCount)
            .Select(c => result.Labels.Count(L => L == c))
            .OrderByDescending(s => s)
            .ToArray();
        Console.WriteLine($"    K={result.ClusterCount}, sizes=[{string.Join(",", sizes)}], noise={noise}, purity={purity:F3}");

        bool ok = result.ClusterCount == 4 && purity >= 0.95 && noise <= 5;
        Console.WriteLine(ok ? "    PASS\n" : "    FAIL\n");
        return ok ? 0 : 1;
    }

    private static void ReportGeometry(Synthetic.SyntheticDataset ds)
    {
        double maxWithin = 0;
        double minBetween = double.PositiveInfinity;
        for (int i = 0; i < ds.Features.Length; i++)
        for (int j = i + 1; j < ds.Features.Length; j++)
        {
            double d = EuclideanDist(ds.Features[i], ds.Features[j]);
            if (ds.Labels[i] == ds.Labels[j])
            {
                if (d > maxWithin) maxWithin = d;
            }
            else
            {
                if (d < minBetween) minBetween = d;
            }
        }
        Console.WriteLine($"    geometry: max-within={maxWithin:F3}, min-between={minBetween:F3}, gap-ratio={minBetween / maxWithin:F2}");
    }

    private static double EuclideanDist(double[] a, double[] b)
    {
        double sum = 0;
        for (int i = 0; i < a.Length; i++) sum += (a[i] - b[i]) * (a[i] - b[i]);
        return Math.Sqrt(sum);
    }

    // ── Scenario 2 ─────────────────────────────────────────────────────────
    private static int Scenario_SingleBlobPlusOutliers()
    {
        Console.WriteLine("[2] 50pt dense blob + 3 outliers, allowSingleCluster=true");

        var rng = new Random(7);
        const int dense = 50;
        var pts = new double[dense + 3][];
        for (int i = 0; i < dense; i++)
            pts[i] = new[] { rng.NextDouble() * 0.5, rng.NextDouble() * 0.5 };
        pts[dense]     = new[] { 10.0,  10.0 };
        pts[dense + 1] = new[] { -8.0,   9.0 };
        pts[dense + 2] = new[] {  9.0,  -8.0 };

        var result = HdbscanSession.Run(pts, new HdbscanSettings
        {
            MinPts             = 5,
            AllowSingleCluster = true,
        });

        double avgProbDense    = result.MembershipProbabilities.Take(dense).Average();
        double avgProbOutliers = result.MembershipProbabilities.Skip(dense).Average();
        int    noise           = result.Labels.Count(L => L < 0);

        Console.WriteLine($"    K={result.ClusterCount}, noise={noise}, prob(dense)={avgProbDense:F3}, prob(outliers)={avgProbOutliers:F3}");

        // With allowSingleCluster=true we expect one cluster covering both groups,
        // but outliers should have markedly lower membership probability.
        bool ok = result.ClusterCount == 1
                  && avgProbDense    > 0.5
                  && avgProbOutliers < avgProbDense;
        Console.WriteLine(ok ? "    PASS\n" : "    FAIL\n");
        return ok ? 0 : 1;
    }

    // ── Scenario 3 ─────────────────────────────────────────────────────────
    private static int Scenario_ClustererAdapter()
    {
        Console.WriteLine("[3] 4 blobs via HdbscanClusterer (IClusterer adapter)");

        var ds = SpatialBlobs.Generate(
            clusterCount: 4, pointsPerCluster: 50, dimensions: 2,
            separation: 10.0, spread: 0.3, seed: 42);

        IClusterer clusterer = new HdbscanClusterer(minPts: 5, allowSingleCluster: true);
        ClusterResult result = clusterer.Cluster(ds.Features);

        int    noise  = result.Labels.Count(L => L < 0);
        double purity = PerClusterPurity(result.Labels, ds.Labels, result.K);

        Console.WriteLine($"    Name='{clusterer.Name}'");
        Console.WriteLine($"    K={result.K}, noise={noise}, purity={purity:F3}");

        bool ok = result.K == 4 && purity >= 0.95 && noise <= 5;
        Console.WriteLine(ok ? "    PASS\n" : "    FAIL\n");
        return ok ? 0 : 1;
    }

    // ── helpers ────────────────────────────────────────────────────────────

    /// <summary>
    /// Per-cluster purity, weighted by cluster size. For each predicted cluster
    /// C, count the dominant ground-truth label's share; sum across clusters
    /// and normalise by total clustered points. Noise (label -1) is ignored.
    /// </summary>
    private static double PerClusterPurity(int[] predicted, int[] truth, int k)
    {
        if (k == 0) return 0.0;
        double totalPure = 0.0;
        int totalClustered = 0;
        for (int c = 0; c < k; c++)
        {
            var counts = new Dictionary<int, int>();
            int size = 0;
            for (int i = 0; i < predicted.Length; i++)
            {
                if (predicted[i] != c) continue;
                size++;
                counts.TryGetValue(truth[i], out int n);
                counts[truth[i]] = n + 1;
            }
            if (size == 0) continue;
            int dominant = counts.Values.Max();
            totalPure       += dominant;
            totalClustered  += size;
        }
        return totalClustered == 0 ? 0.0 : totalPure / totalClustered;
    }
}
