using System;
using System.Linq;
using Maths.Rng;
using Xunit;

namespace Clustering.Geometric.EMeans.Tests;

/// <summary>
/// A seeded sweep over problem sizes, data shapes and settings, corners included (k = n, one cluster, batches of one
/// point, a step of exactly one, coincident points), holding every run to the invariants the method promises.
/// </summary>
public class EMeansSweepTests
{
    private const int Cases = 400;

    [Fact]
    public void EveryRunHoldsItsInvariants()
    {
        int[] seeds = SeedTree.Derive(20261114, Cases);
        for (int index = 0; index < Cases; index++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[index]);
            int d = Pick(rng, 1, 2, 7);
            int n = 2 + rng.NextInt(299);
            int k = Pick(rng, 1, n, 1 + rng.NextInt(n), 1 + rng.NextInt(Math.Min(n, 12)));
            double[][] data = Draw(rng, n, d, kind: rng.NextInt(4));

            var config = new EMeansConfig
            {
                K = k,
                Epochs = 1 + rng.NextInt(4),
                Seed = rng.NextInt(1 << 20),
                Init = rng.NextInt(2) == 0 ? EMeansInit.RandomPoints : EMeansInit.KMeansPlusPlus,
                Batch = rng.NextInt(3) == 0
                    ? new BatchConfig()
                    : new BatchConfig { Strategy = BatchStrategy.Explicit, Size = Pick(rng, 1, n, 1 + rng.NextInt(n)) },
                Step = rng.NextInt(3) == 0
                    ? new StepConfig()
                    : new StepConfig
                    {
                        Strategy = StepStrategy.Explicit,
                        Peak = Pick(rng, 1.0, 1e-3, 1.0 - rng.NextDouble()),
                        Schedule = rng.NextInt(2) == 0 ? StepSchedule.Cosine : StepSchedule.Constant,
                        FloorFraction = Pick(rng, 0.0, 0.1, 1.0),
                    },
                Revival = new RevivalConfig
                {
                    Enabled = rng.NextInt(4) != 0,
                    Period = 1 + rng.NextInt(5),
                    Threshold = Pick(rng, 0.0, 0.1, 0.9),
                },
                Finalization = (EMeansFinalization)rng.NextInt(3),
                SlowStepMultiplier = Pick(rng, 0.25, 1.0, 4.0),
            };

            string label = $"case {index} (seed {seeds[index]}): n={n} d={d} k={k}";
            EMeansResult result = EMeans.Cluster(data, config);
            EMeansManifest manifest = result.Manifest;

            int batches = (n + manifest.BatchSize - 1) / manifest.BatchSize;
            Assert.True(manifest.Updates == config.Epochs * batches, label);
            Assert.True(manifest.Epochs.Count == config.Epochs, label);
            Assert.True(manifest.PeakStep > 0.0 && manifest.PeakStep <= 1.0, label);

            // Statistics stay positive and finite, and the EMA views are their ratios.
            for (int j = 0; j < k; j++)
            {
                Assert.True(result.Counts[j] > 0.0 && double.IsFinite(result.Counts[j]), label);
                Assert.True(result.SlowCounts[j] > 0.0 && double.IsFinite(result.SlowCounts[j]), label);
                for (int c = 0; c < d; c++)
                {
                    Assert.True(result.EmaCentroids[j][c] == result.Sums[j][c] / result.Counts[j], label);
                    Assert.True(result.SlowEmaCentroids[j][c] == result.SlowSums[j][c] / result.SlowCounts[j], label);
                }
            }

            // A centroid is a weighted mean of data points and split children. A child sits one offset from its
            // donor, an offset is at most the perturbation times the box's diagonal, and a lineage split again
            // before it receives a point stacks offsets: without revival no view leaves the data's box, and with it
            // none leaves by more than a few offsets.
            double diagonal = 0.0;
            for (int c = 0; c < d; c++)
            {
                double side = data.Max(x => x[c]) - data.Min(x => x[c]);
                diagonal += side * side;
            }
            diagonal = Math.Sqrt(diagonal);
            double slack = config.Revival.Enabled ? 16.0 * config.Revival.Perturbation * diagonal : 0.0;
            for (int j = 0; j < k; j++)
            {
                Assert.True(result.Scatter[j] >= 0.0 && double.IsFinite(result.Scatter[j]), label);
                double radius = Math.Sqrt(result.Scatter[j] / result.Counts[j]);
                double scale = data.Max(x => x.Max(Math.Abs));
                Assert.True(
                    radius <= diagonal + slack + 1e-9 * Math.Max(1.0, scale),
                    $"{label}: cluster {j} radius {radius:G4} exceeds the box diagonal {diagonal:G4}");
            }
            string[] names = { "selected", "EMA", "slow EMA", "last-epoch mean" };
            double[][][] views =
            {
                result.Centroids, result.EmaCentroids, result.SlowEmaCentroids, result.LastEpochMeanCentroids,
            };
            for (int c = 0; c < d; c++)
            {
                double low = data.Min(x => x[c]);
                double high = data.Max(x => x[c]);
                double margin = slack + 1e-9 * Math.Max(1.0, Math.Max(Math.Abs(low), Math.Abs(high)));
                for (int v = 0; v < views.Length; v++)
                {
                    for (int j = 0; j < k; j++)
                    {
                        double value = views[v][j][c];
                        double excess = Math.Max(low - value, value - high);
                        Assert.True(
                            excess <= margin,
                            $"{label}: {names[v]} centroid {j} leaves the box by {excess:G4} (margin {margin:G4}, diagonal {diagonal:G4})");
                    }
                }
            }

            // Labels are the nearest selected centroid, and the distortion is theirs.
            double total = 0.0;
            for (int i = 0; i < n; i++)
            {
                int nearest = 0;
                double best = double.PositiveInfinity;
                for (int j = 0; j < k; j++)
                {
                    double distance = 0.0;
                    for (int c = 0; c < d; c++)
                    {
                        double diff = data[i][c] - result.Centroids[j][c];
                        distance += diff * diff;
                    }
                    if (distance < best)
                    {
                        best = distance;
                        nearest = j;
                    }
                }
                Assert.True(result.Labels[i] == nearest, $"{label}: point {i} is not at its nearest centroid");
                total += best;
            }
            Assert.True(Math.Abs(result.Distortion - total / n) <= 1e-12 * Math.Max(1.0, total / n), label);

            // A rerun reproduces the run.
            if (index % 8 == 0)
            {
                EMeansResult again = EMeans.Cluster(data, config);
                Assert.True(result.Labels.SequenceEqual(again.Labels), $"{label}: rerun differs");
                Assert.True(result.Distortion == again.Distortion, $"{label}: rerun differs");
            }
        }
    }

    private static T Pick<T>(Xoshiro256PlusPlus rng, params T[] options) => options[rng.NextInt(options.Length)];

    // kind 0: Gaussian cloud; 1: small-integer lattice (many coincident points); 2: every point the same;
    // 3: two tight clumps far from the origin and from each other.
    private static double[][] Draw(Xoshiro256PlusPlus rng, int n, int d, int kind)
    {
        var data = new double[n][];
        for (int i = 0; i < n; i++)
        {
            var x = new double[d];
            for (int c = 0; c < d; c++)
            {
                x[c] = kind switch
                {
                    0 => Normal(rng),
                    1 => rng.NextInt(3) - 1,
                    2 => 4.25,
                    _ => (i % 2 == 0 ? 1e4 : -3e3) + 1e-3 * Normal(rng),
                };
            }
            data[i] = x;
        }
        return data;
    }

    private static double Normal(Xoshiro256PlusPlus rng)
    {
        double u = 1.0 - rng.NextDouble();
        double v = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u)) * Math.Cos(2.0 * Math.PI * v);
    }
}
