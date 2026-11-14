using System;
using System.Linq;
using Xunit;

namespace Clustering.Geometric.EMeans.Tests;

public class EMeansTests
{
    private static EMeansConfig Explicit(int k, int epochs, int batch, double peak) => new()
    {
        K = k,
        Epochs = epochs,
        Batch = new BatchConfig { Strategy = BatchStrategy.Explicit, Size = batch },
        Step = new StepConfig { Strategy = StepStrategy.Explicit, Peak = peak },
    };

    // Brute-force Lloyd: assign to the nearest centroid (ties to the lowest index), move each centroid to the mean
    // of its points, leave a centroid without points where it is.
    private static double[][] Lloyd(double[][] data, double[][] initial, int iterations)
    {
        int k = initial.Length;
        int d = data[0].Length;
        double[][] centroids = initial.Select(row => (double[])row.Clone()).ToArray();
        for (int iteration = 0; iteration < iterations; iteration++)
        {
            var sums = new double[k][];
            var counts = new int[k];
            for (int j = 0; j < k; j++) sums[j] = new double[d];
            foreach (double[] x in data)
            {
                int j = NearestIndex(x, centroids);
                counts[j]++;
                for (int c = 0; c < d; c++) sums[j][c] += x[c];
            }
            for (int j = 0; j < k; j++)
            {
                if (counts[j] == 0) continue;
                for (int c = 0; c < d; c++) centroids[j][c] = sums[j][c] / counts[j];
            }
        }
        return centroids;
    }

    private static int NearestIndex(double[] x, double[][] centroids)
    {
        int best = 0;
        double bestDistance = double.PositiveInfinity;
        for (int j = 0; j < centroids.Length; j++)
        {
            double distance = 0.0;
            for (int c = 0; c < x.Length; c++)
            {
                double diff = x[c] - centroids[j][c];
                distance += diff * diff;
            }
            if (distance < bestDistance)
            {
                bestDistance = distance;
                best = j;
            }
        }
        return best;
    }

    private static void AssertClose(double[][] expected, double[][] actual, double tolerance)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int j = 0; j < expected.Length; j++)
        {
            for (int c = 0; c < expected[j].Length; c++)
            {
                double scale = Math.Max(1.0, Math.Abs(expected[j][c]));
                Assert.True(
                    Math.Abs(expected[j][c] - actual[j][c]) <= tolerance * scale,
                    $"centroid {j}, coordinate {c}: expected {expected[j][c]:R}, got {actual[j][c]:R}");
            }
        }
    }

    // The symmetric integer grid {-5..5}^2: its mean is exactly the origin in floating point.
    private static double[][] SymmetricGrid()
    {
        var points = new double[121][];
        for (int x = -5, i = 0; x <= 5; x++)
            for (int y = -5; y <= 5; y++, i++)
                points[i] = new double[] { x, y };
        return points;
    }

    [Fact]
    public void FullBatchAtStepOne_IsLloyd()
    {
        (double[][] data, _) = Blobs.Draw(Blobs.Corners(4, 3, 2.0), perBlob: 75, sigma: 1.0, seed: 3);
        // Four data rows and one centroid no point is ever nearest to: the last exercises a cluster left with no
        // mass at step one, which must hold its position as Lloyd's empty cluster does.
        double[][] initial =
        {
            data[0], data[80], data[160], data[240], new[] { 1000.0, 1000.0, 1000.0 },
        };
        const int iterations = 6;

        EMeansConfig config = Explicit(k: 5, epochs: iterations, batch: data.Length, peak: 1.0) with
        {
            Step = new StepConfig
            {
                Strategy = StepStrategy.Explicit, Peak = 1.0, Schedule = StepSchedule.Constant,
            },
            Revival = new RevivalConfig { Enabled = false },
            Finalization = EMeansFinalization.Ema,
        };
        EMeansResult result = EMeans.Cluster(data, config, initial);

        double[][] reference = Lloyd(data, initial, iterations);
        AssertClose(reference, result.EmaCentroids, 1e-10);
        AssertClose(reference, result.LastEpochMeanCentroids, 1e-10);
        Assert.Equal(data.Select(x => NearestIndex(x, reference)).ToArray(), result.Labels);
        Assert.DoesNotContain(4, result.Labels);
    }

    [Fact]
    public void SlowEma_RunAtTheTrainingStep_ReproducesTheTrainingEma()
    {
        (double[][] data, _) = Blobs.Draw(Blobs.Corners(4, 5, 3.0), perBlob: 100, sigma: 1.0, seed: 5);
        double[][] initial = Enumerable.Range(0, 8).Select(j => (double[])data[j * 50].Clone()).ToArray();
        initial[7] = Enumerable.Repeat(500.0, 5).ToArray(); // dead from the start, so the run splits a donor

        // Constant schedule at 0.5 with floor fraction 0.5 and multiplier 2: the slow step is 2 * 0.5 * 0.5 = 0.5.
        EMeansConfig config = Explicit(k: 8, epochs: 3, batch: 32, peak: 0.5) with
        {
            Step = new StepConfig
            {
                Strategy = StepStrategy.Explicit, Peak = 0.5, Schedule = StepSchedule.Constant, FloorFraction = 0.5,
            },
            SlowStepMultiplier = 2.0,
        };
        EMeansResult result = EMeans.Cluster(data, config, initial);

        Assert.Equal(0.5, result.Manifest.SlowStep);
        Assert.True(result.Manifest.Epochs.Sum(epoch => epoch.Revivals) >= 1);
        Assert.Equal(result.Counts, result.SlowCounts);
        for (int j = 0; j < 8; j++)
        {
            Assert.Equal(result.Sums[j], result.SlowSums[j]);
            Assert.Equal(result.EmaCentroids[j], result.SlowEmaCentroids[j]);
        }
    }

    [Fact]
    public void Run_IsAFunctionOfDataConfigAndSeed()
    {
        (double[][] data, _) = Blobs.Draw(Blobs.Corners(5, 6, 2.0), perBlob: 120, sigma: 1.0, seed: 9);
        EMeansConfig config = Explicit(k: 12, epochs: 4, batch: 48, peak: 0.6) with { Seed = 17 };

        EMeansResult first = EMeans.Cluster(data, config);
        EMeansResult second = EMeans.Cluster(data, config);
        EMeansResult reseeded = EMeans.Cluster(data, config with { Seed = 18 });

        Assert.Equal(first.Labels, second.Labels);
        Assert.Equal(first.Distortion, second.Distortion);
        for (int j = 0; j < 12; j++)
        {
            Assert.Equal(first.Centroids[j], second.Centroids[j]);
            Assert.Equal(first.Sums[j], second.Sums[j]);
        }
        Assert.NotEqual(first.Manifest.StreamSeeds, reseeded.Manifest.StreamSeeds);
        Assert.False(first.Centroids.Zip(reseeded.Centroids).All(pair => pair.First.SequenceEqual(pair.Second)));
    }

    [Fact]
    public void FromOneSeedPerBlob_RecoversThePlantedPartition()
    {
        double[][] centres = Blobs.Corners(6, 8, 5.0);
        (double[][] data, int[] planted) = Blobs.Draw(centres, perBlob: 200, sigma: 1.0, seed: 21);
        double[][] initial = Enumerable.Range(0, 6).Select(g => (double[])data[g * 200].Clone()).ToArray();

        EMeansResult result = EMeans.Cluster(data, Explicit(k: 6, epochs: 5, batch: 64, peak: 0.5), initial);

        Assert.True(Blobs.SamePartition(planted, result.Labels));
        // Each centroid is a mean of 200 unit-variance points: within 5 standard errors of its centre per coordinate.
        for (int g = 0; g < 6; g++)
        {
            double[] centroid = result.Centroids[result.Labels[g * 200]];
            for (int c = 0; c < 8; c++)
                Assert.InRange(centroid[c] - centres[g][c], -5.0 / Math.Sqrt(200), 5.0 / Math.Sqrt(200));
        }
    }

    [Fact]
    public void Revival_SplitsADonorWhoseCentroidSitsAtTheOrigin()
    {
        // One centroid at the grid's exact mean, one that no point is nearest to. Scaling the donor's coordinates
        // by 1 ± δ (the paper's split) would leave both children at the origin and the revived cluster empty.
        double[][] data = SymmetricGrid();
        double[][] initial = { new[] { 0.0, 0.0 }, new[] { 1e6, 1e6 } };
        EMeansConfig config = Explicit(k: 2, epochs: 10, batch: data.Length, peak: 1.0) with
        {
            Step = new StepConfig
            {
                Strategy = StepStrategy.Explicit, Peak = 1.0, Schedule = StepSchedule.Constant,
            },
            Finalization = EMeansFinalization.Ema,
        };

        EMeansResult result = EMeans.Cluster(data, config, initial);

        Assert.Equal(1, result.Manifest.Epochs.Sum(epoch => epoch.Revivals));
        Assert.InRange(result.Labels.Count(label => label == 0), 30, 91);
        Assert.InRange(result.Labels.Count(label => label == 1), 30, 91);
        // One cluster holds the grid at its total variance of 20; two halves hold it well below.
        Assert.True(result.Distortion < 15.0, $"distortion {result.Distortion}");
    }

    [Fact]
    public void Revival_SplitsATightDonorFarFromTheOriginWithinItsOwnRadius()
    {
        // Two clumps of radius ~1e-3, 1.8e4 apart and far from the origin; two centroids share the first clump.
        // An offset sized by the data's overall spread, or by the centroid's own magnitude, would put both
        // children of a split outside their clump, where the clump's other centroid takes every point.
        double[][] centres = { new[] { 1e4, 1e4 }, new[] { -3e3, -3e3 } };
        (double[][] data, int[] planted) = Blobs.Draw(centres, perBlob: 100, sigma: 1e-3, seed: 43);
        double[][] initial =
        {
            (double[])data[0].Clone(), (double[])data[1].Clone(), (double[])data[100].Clone(), new[] { 1e7, 1e7 },
        };
        EMeansConfig config = Explicit(k: 4, epochs: 10, batch: data.Length, peak: 1.0) with
        {
            Step = new StepConfig
            {
                Strategy = StepStrategy.Explicit, Peak = 1.0, Schedule = StepSchedule.Constant,
            },
            Finalization = EMeansFinalization.Ema,
        };

        EMeansResult result = EMeans.Cluster(data, config, initial);

        Assert.Equal(1, result.Manifest.Epochs.Sum(epoch => epoch.Revivals));
        for (int j = 0; j < 4; j++)
        {
            int[] members = Enumerable.Range(0, data.Length).Where(i => result.Labels[i] == j).ToArray();
            Assert.True(members.Length >= 10, $"cluster {j} holds {members.Length} points");
            Assert.Single(members.Select(i => planted[i]).Distinct());
        }
    }

    [Fact]
    public void WithoutRevival_AClusterWithNoMassHoldsItsPosition()
    {
        double[][] data = SymmetricGrid();
        double[][] initial = { new[] { 0.0, 0.0 }, new[] { 1e6, 1e6 } };
        EMeansConfig config = Explicit(k: 2, epochs: 10, batch: data.Length, peak: 1.0) with
        {
            Step = new StepConfig
            {
                Strategy = StepStrategy.Explicit, Peak = 1.0, Schedule = StepSchedule.Constant,
            },
            Revival = new RevivalConfig { Enabled = false },
            Finalization = EMeansFinalization.Ema,
        };

        EMeansResult result = EMeans.Cluster(data, config, initial);

        Assert.All(result.Labels, label => Assert.Equal(0, label));
        Assert.True(result.Counts[1] > 0.0);
        Assert.Equal(1e6, result.Centroids[1][0], 3);
        Assert.Equal(1e6, result.Centroids[1][1], 3);
        Assert.Equal(new[] { 0.0, 0.0 }, result.Centroids[0]);
    }

    [Fact]
    public void Run_IsEquivariantUnderTranslationAndScale()
    {
        (double[][] data, _) = Blobs.Draw(Blobs.Corners(3, 4, 4.0), perBlob: 100, sigma: 1.0, seed: 13);
        double[][] initial =
        {
            (double[])data[0].Clone(), (double[])data[150].Clone(), new[] { 1e4, 1e4, 1e4, 1e4 },
        };
        const double scale = 3.0;
        double[] shift = { 100.0, -50.0, 25.0, 7.0 };
        double[] Map(double[] x) => x.Select((value, c) => scale * value + shift[c]).ToArray();

        EMeansConfig config = Explicit(k: 3, epochs: 4, batch: 50, peak: 0.6) with { Seed = 11 };
        EMeansResult original = EMeans.Cluster(data, config, initial);
        EMeansResult mapped = EMeans.Cluster(
            data.Select(Map).ToArray(), config, initial.Select(Map).ToArray());

        Assert.True(original.Manifest.Epochs.Sum(epoch => epoch.Revivals) >= 1);
        Assert.Equal(
            original.Manifest.Epochs.Select(epoch => epoch.Revivals),
            mapped.Manifest.Epochs.Select(epoch => epoch.Revivals));
        Assert.Equal(original.Labels, mapped.Labels);
        AssertClose(original.Centroids.Select(Map).ToArray(), mapped.Centroids, 1e-9);
    }

    [Fact]
    public void WithoutRevival_EveryCentroidViewStaysInsideTheDataBox()
    {
        (double[][] data, _) = Blobs.Draw(Blobs.Corners(4, 5, 2.0), perBlob: 125, sigma: 1.0, seed: 29);
        var config = new EMeansConfig
        {
            K = 10, Epochs = 3, Revival = new RevivalConfig { Enabled = false },
        };

        EMeansResult result = EMeans.Cluster(data, config);

        for (int c = 0; c < 5; c++)
        {
            double low = data.Min(x => x[c]) - 1e-9;
            double high = data.Max(x => x[c]) + 1e-9;
            foreach (double[][] view in new[]
                     {
                         result.EmaCentroids, result.SlowEmaCentroids, result.LastEpochMeanCentroids,
                     })
            {
                Assert.All(view, centroid => Assert.InRange(centroid[c], low, high));
            }
        }
    }

    [Fact]
    public void Manifest_RecordsWhatTheAutoPicksResolved()
    {
        (double[][] data, _) = Blobs.Draw(Blobs.Corners(5, 4, 3.0), perBlob: 200, sigma: 1.0, seed: 31);
        var config = new EMeansConfig { K = 5, Epochs = 5 };

        EMeansResult result = EMeans.Cluster(data, config);
        EMeansManifest manifest = result.Manifest;

        // n = 1000, k = 5, T = 5: the batch law gives 69.4 -> 64, the step law 0.620978 at that batch.
        Assert.Same(config, manifest.Config);
        Assert.Equal(1000, manifest.N);
        Assert.Equal(4, manifest.Dimension);
        Assert.Equal(64, manifest.BatchSize);
        Assert.Equal(0.620978, manifest.PeakStep, 5);
        Assert.Equal(0.1 * manifest.PeakStep, manifest.SlowStep, 12);
        Assert.Equal(5 * 16, manifest.Updates);
        Assert.Equal(3, manifest.StreamSeeds.Length);
        Assert.Equal(5, manifest.Epochs.Count);
        Assert.All(manifest.Epochs, epoch => Assert.True(epoch.OnlineDistortion > 0.0));
        Assert.Contains(manifest.AutomationLog, entry => entry.Contains("extrapolated"));
    }

    [Fact]
    public void Assign_ScoresTheSelectedCentroids()
    {
        (double[][] data, _) = Blobs.Draw(Blobs.Corners(4, 3, 3.0), perBlob: 80, sigma: 1.0, seed: 37);
        EMeansResult result = EMeans.Cluster(data, Explicit(k: 6, epochs: 3, batch: 40, peak: 0.5));

        (int[] labels, double distortion) = EMeans.Assign(data, result.Centroids);

        Assert.Equal(result.Labels, labels);
        Assert.Equal(result.Distortion, distortion);
    }

    [Fact]
    public void KMeansPlusPlusInit_SeedsFromTheData()
    {
        (double[][] data, _) = Blobs.Draw(Blobs.Corners(4, 3, 3.0), perBlob: 80, sigma: 1.0, seed: 41);
        EMeansConfig config = Explicit(k: 4, epochs: 3, batch: 40, peak: 0.5) with
        {
            Init = EMeansInit.KMeansPlusPlus,
        };

        EMeansResult result = EMeans.Cluster(data, config);

        Assert.Contains(result.Manifest.AutomationLog, entry => entry.Contains("k-means++"));
        Assert.All(result.Labels, label => Assert.InRange(label, 0, 3));
        Assert.All(result.Centroids, centroid => Assert.All(centroid, value => Assert.True(double.IsFinite(value))));
    }

    [Fact]
    public void StepAt_AnnealsFromThePeakTowardsTheFloor()
    {
        const int updates = 100;
        Assert.Equal(0.8, EMeans.StepAt(0, updates, 0.8, 0.1, StepSchedule.Cosine));
        for (int t = 1; t < updates; t++)
        {
            Assert.True(
                EMeans.StepAt(t, updates, 0.8, 0.1, StepSchedule.Cosine)
                < EMeans.StepAt(t - 1, updates, 0.8, 0.1, StepSchedule.Cosine));
        }
        Assert.InRange(EMeans.StepAt(updates - 1, updates, 0.8, 0.1, StepSchedule.Cosine), 0.08, 0.081);
        Assert.Equal(0.8, EMeans.StepAt(updates - 1, updates, 0.8, 0.1, StepSchedule.Constant));
    }

    [Fact]
    public void InvalidInputs_AreRefused()
    {
        double[][] data = SymmetricGrid();
        Assert.Throws<ArgumentException>(() => EMeans.Cluster(data, new EMeansConfig { K = 122, Epochs = 1 }));
        Assert.Throws<ArgumentException>(() => EMeans.Cluster(data, new EMeansConfig { K = 2, Epochs = 0 }));
        Assert.Throws<ArgumentException>(() => EMeans.Cluster(data, new EMeansConfig
        {
            K = 2, Epochs = 1, Batch = new BatchConfig { Strategy = BatchStrategy.Explicit },
        }));
        Assert.Throws<ArgumentException>(() => EMeans.Cluster(data, new EMeansConfig
        {
            K = 2, Epochs = 1, Step = new StepConfig { Strategy = StepStrategy.Explicit, Peak = 1.5 },
        }));
        Assert.Throws<ArgumentException>(() => EMeans.Cluster(
            new[] { new[] { 0.0, 1.0 }, new[] { double.NaN, 2.0 } }, new EMeansConfig { K = 1, Epochs = 1 }));
        Assert.Throws<ArgumentException>(() => EMeans.Cluster(
            data, new EMeansConfig { K = 2, Epochs = 1 }, new[] { new[] { 0.0, 0.0 } }));
    }
}
