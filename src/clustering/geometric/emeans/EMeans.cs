using System;
using System.Collections.Generic;
using Clustering.Geometric.KMeans;
using Maths.Rng;

namespace Clustering.Geometric.EMeans;

/// <summary>
/// e-means (Bessa, Mazaré, Szilvasy, Touvron, Picard, Douze &amp; Jégou 2026, arXiv:2611.09501): mini-batch k-means
/// that accumulates each batch's per-cluster counts and sums into exponentially decayed statistics <c>(N, S)</c> and
/// reads every centroid off as the ratio <c>S/N</c>, so a centroid is always a weighted mean of points assigned to
/// it. The step <c>α = 1 − γ</c> anneals on a cosine across the whole run, dead clusters are revived by splitting a
/// live donor, and the run closes on one of three centroid views. Lloyd's algorithm is the corner <c>B = n</c>,
/// <c>α = 1</c>, constant schedule.
/// </summary>
/// <remarks>
/// Four departures from the paper's text, each recorded with its reason in
/// <c>issues/clustering/emeans-design.md</c>: an epoch's trailing partial batch is trained on; a cluster left with
/// no mass holds its position at a floor mass; a split offsets its two children along a random direction by a
/// fraction of the donor's own radius, read off a decayed scatter kept beside <c>(N, S)</c>; a split in the last
/// epoch restarts the last-epoch records of the two slots it touches.
/// </remarks>
public static class EMeans
{
    // Mass a cluster is carried at once decay would take it lower, keeping S/N defined.
    private const double MassFloor = 1e-200;

    /// <summary>
    /// Runs e-means on <paramref name="data"/> (one row per point). <paramref name="initialCentroids"/>, when given,
    /// replaces the config's initialization. The run is a function of the data, the config and those centroids.
    /// </summary>
    public static EMeansResult Cluster(double[][] data, EMeansConfig config, double[][]? initialCentroids = null)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(config);
        int n = data.Length;
        if (n == 0) throw new ArgumentException("Data cannot be empty.", nameof(data));
        int d = ValidateRows(data, nameof(data));
        ValidateConfig(config, n);
        int k = config.K;
        int epochs = config.Epochs;

        var log = new List<string>();

        bool autoBatch = config.Batch.Strategy == BatchStrategy.Auto;
        int batch;
        if (autoBatch)
        {
            batch = EMeansScalingLaw.BatchSize(n, k, epochs);
            string divided = n > EMeansScalingLaw.UndividedBatchLimit ? ", the law divided by 8 above 2M points" : "";
            log.Add($"batch: scaling law gives B={batch}{divided}");
        }
        else
        {
            batch = config.Batch.Size
                ?? throw new ArgumentException("An explicit batch needs Batch.Size.", nameof(config));
            if (batch < 1 || batch > n)
                throw new ArgumentException("Batch.Size must be between 1 and n.", nameof(config));
            log.Add($"batch: explicit B={batch}");
        }

        bool autoStep = config.Step.Strategy == StepStrategy.Auto;
        double peak;
        if (autoStep)
        {
            double unclipped = EMeansScalingLaw.UnclippedPeakStep(n, k, epochs, batch);
            peak = Math.Min(unclipped, 1.0);
            log.Add(unclipped > 1.0
                ? $"step: scaling law gives {unclipped:G6}, capped to a peak of 1"
                : $"step: scaling law gives a peak of {peak:G6}");
        }
        else
        {
            peak = config.Step.Peak
                ?? throw new ArgumentException("An explicit step needs Step.Peak.", nameof(config));
            if (!(peak > 0.0 && peak <= 1.0))
                throw new ArgumentException("Step.Peak must lie in (0, 1].", nameof(config));
            log.Add($"step: explicit peak of {peak:G6}");
        }

        if (autoBatch || autoStep)
        {
            IReadOnlyList<string> crossed = EMeansScalingLaw.Extrapolations(n, k, epochs, batch);
            if (crossed.Count > 0)
                log.Add("scaling law extrapolated beyond its fitting grid: " + string.Join("; ", crossed));
        }

        int updatesPerEpoch = (int)(((long)n + batch - 1) / batch);
        long totalUpdates = (long)epochs * updatesPerEpoch;
        if (totalUpdates > int.MaxValue)
            throw new ArgumentException("Epochs times batches per epoch exceeds the update counter.", nameof(config));
        int updates = (int)totalUpdates;

        double floorFraction = config.Step.FloorFraction;
        double slowStep = Math.Min(1.0, config.SlowStepMultiplier * floorFraction * peak);

        int[] streams = SeedTree.Derive(config.Seed, 3);
        var initRng = new Xoshiro256PlusPlus(streams[0]);
        var shuffleRng = new Xoshiro256PlusPlus(streams[1]);
        var revivalRng = new Xoshiro256PlusPlus(streams[2]);

        var perm = new int[n];
        var centroids = new double[k * d];
        if (initialCentroids != null)
        {
            if (initialCentroids.Length != k)
                throw new ArgumentException("Initial centroids must number k.", nameof(initialCentroids));
            if (ValidateRows(initialCentroids, nameof(initialCentroids)) != d)
                throw new ArgumentException("Initial centroids must match the data's dimension.", nameof(initialCentroids));
            for (int j = 0; j < k; j++) Array.Copy(initialCentroids[j], 0, centroids, j * d, d);
            log.Add("init: centroids supplied by the caller");
        }
        else if (config.Init == EMeansInit.KMeansPlusPlus)
        {
            double[][] seeds = KMeansPlusPlusSeeding.Seed(data, k, initRng);
            for (int j = 0; j < k; j++) Array.Copy(seeds[j], 0, centroids, j * d, d);
            log.Add("init: k-means++ seeding");
        }
        else
        {
            // The first k entries of a partial Fisher-Yates shuffle are k distinct uniform draws.
            for (int i = 0; i < n; i++) perm[i] = i;
            for (int j = 0; j < k; j++)
            {
                int other = j + UniformIndex(initRng, n - j);
                (perm[j], perm[other]) = (perm[other], perm[j]);
                Array.Copy(data[perm[j]], 0, centroids, j * d, d);
            }
            log.Add("init: k distinct random points");
        }

        double initialCount = config.InitialCount;
        var counts = new double[k];
        var sums = new double[k * d];
        var scatter = new double[k];
        var slowCounts = new double[k];
        var slowSums = new double[k * d];
        for (int j = 0; j < k; j++)
        {
            counts[j] = initialCount;
            slowCounts[j] = initialCount;
            for (int c = j * d; c < (j + 1) * d; c++)
            {
                sums[c] = initialCount * centroids[c];
                slowSums[c] = sums[c];
            }
        }

        var batchCounts = new int[k];
        var batchSums = new double[k * d];
        var batchScatter = new double[k];
        var lastCounts = new int[k];
        var lastSums = new double[k * d];
        var donorWeights = new double[k];
        var direction = new double[d];
        var traces = new List<EMeansEpochTrace>(epochs);
        RevivalConfig revival = config.Revival;
        int update = 0;
        int totalRevivals = 0;

        for (int epoch = 0; epoch < epochs; epoch++)
        {
            bool lastEpoch = epoch == epochs - 1;
            Shuffle(perm, shuffleRng);
            double online = 0.0;
            int revivals = 0;

            int start = 0;
            while (start < n)
            {
                int end = start + Math.Min(batch, n - start);
                Array.Clear(batchCounts);
                Array.Clear(batchSums);
                Array.Clear(batchScatter);
                for (int p = start; p < end; p++)
                {
                    double[] x = data[perm[p]];
                    int j = Nearest(x, centroids, k, d, out double distance);
                    online += distance;
                    batchCounts[j]++;
                    batchScatter[j] += distance;
                    int offset = j * d;
                    for (int c = 0; c < d; c++) batchSums[offset + c] += x[c];
                    if (lastEpoch)
                    {
                        lastCounts[j]++;
                        for (int c = 0; c < d; c++) lastSums[offset + c] += x[c];
                    }
                }

                double alpha = StepAt(update, updates, peak, floorFraction, config.Step.Schedule);
                Decay(counts, sums, scatter, batchCounts, batchSums, batchScatter, alpha, k, d);
                Decay(slowCounts, slowSums, null, batchCounts, batchSums, batchScatter, slowStep, k, d);
                update++;

                if (revival.Enabled && update % revival.Period == 0)
                {
                    revivals += Revive(
                        counts, sums, scatter, slowCounts, slowSums, k, d, revival.Threshold, revival.Perturbation,
                        revivalRng, donorWeights, direction, lastEpoch, lastCounts, lastSums);
                }

                for (int j = 0, offset = 0; j < k; j++, offset += d)
                {
                    double count = counts[j];
                    for (int c = 0; c < d; c++) centroids[offset + c] = sums[offset + c] / count;
                }

                start = end;
            }

            traces.Add(new EMeansEpochTrace(epoch, online / n, revivals));
            totalRevivals += revivals;
        }

        double[][] ema = Ratios(sums, counts, k, d);
        double[][] slow = Ratios(slowSums, slowCounts, k, d);
        var lastEpochMean = new double[k][];
        int fallbacks = 0;
        for (int j = 0; j < k; j++)
        {
            if (lastCounts[j] > 0)
            {
                var mean = new double[d];
                for (int c = 0; c < d; c++) mean[c] = lastSums[j * d + c] / lastCounts[j];
                lastEpochMean[j] = mean;
            }
            else
            {
                lastEpochMean[j] = (double[])ema[j].Clone();
                fallbacks++;
            }
        }

        if (revival.Enabled) log.Add($"revival: {totalRevivals} dead clusters revived");
        if (fallbacks > 0)
            log.Add($"last-epoch mean: {fallbacks} clusters without last-epoch records carry their EMA centroid");

        double[][] selected = config.Finalization switch
        {
            EMeansFinalization.LastEpochMean => lastEpochMean,
            EMeansFinalization.Ema => ema,
            EMeansFinalization.SlowEma => slow,
            _ => throw new ArgumentOutOfRangeException(nameof(config), "Unknown finalization."),
        };
        double[][] final = Copy(selected);
        (int[] labels, double distortion) = Assign(data, final);

        var manifest = new EMeansManifest(
            config, n, d, batch, peak, slowStep, updates, streams, log, traces);
        return new EMeansResult(
            labels, final, k, distortion,
            counts, Rows(sums, k, d), scatter, slowCounts, Rows(slowSums, k, d),
            ema, slow, lastEpochMean, manifest);
    }

    /// <summary>
    /// Assigns every point to its nearest centroid (ties to the lowest index) and returns the labels with the mean
    /// squared distance. Scores any centroid view of a run, on the training points or on held-out ones.
    /// </summary>
    public static (int[] Labels, double Distortion) Assign(double[][] data, double[][] centroids)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(centroids);
        if (data.Length == 0) throw new ArgumentException("Data cannot be empty.", nameof(data));
        if (centroids.Length == 0) throw new ArgumentException("Centroids cannot be empty.", nameof(centroids));
        int d = ValidateRows(data, nameof(data));
        if (ValidateRows(centroids, nameof(centroids)) != d)
            throw new ArgumentException("Centroids must match the data's dimension.", nameof(centroids));

        int k = centroids.Length;
        var flat = new double[k * d];
        for (int j = 0; j < k; j++) Array.Copy(centroids[j], 0, flat, j * d, d);

        var labels = new int[data.Length];
        double total = 0.0;
        for (int i = 0; i < data.Length; i++)
        {
            labels[i] = Nearest(data[i], flat, k, d, out double distance);
            total += distance;
        }
        return (labels, total / data.Length);
    }

    /// <summary>
    /// The step <c>α</c> of update <paramref name="update"/> (zero-based) in a run of <paramref name="updates"/>:
    /// under the cosine schedule <c>α₀·(ε + (1 − ε)·½(1 + cos(π·t/M)))</c>, the peak on the first update and
    /// approaching <c>ε·α₀</c> on the last.
    /// </summary>
    public static double StepAt(int update, int updates, double peak, double floorFraction, StepSchedule schedule)
    {
        if (updates < 1) throw new ArgumentOutOfRangeException(nameof(updates));
        if (update < 0 || update >= updates) throw new ArgumentOutOfRangeException(nameof(update));
        if (schedule == StepSchedule.Constant) return peak;
        double cosine = 0.5 * (1.0 + Math.Cos(Math.PI * update / updates));
        return peak * (floorFraction + (1.0 - floorFraction) * cosine);
    }

    // One EMA update of (N, S), and of the scatter when one is kept, towards a batch's statistics. A cluster whose
    // mass would fall below the floor keeps the position and radius the update gives it (its previous ones when no
    // mass is left at all) at the floor mass.
    private static void Decay(
        double[] counts, double[] sums, double[]? scatter,
        int[] batchCounts, double[] batchSums, double[] batchScatter, double alpha, int k, int d)
    {
        double gamma = 1.0 - alpha;
        for (int j = 0, offset = 0; j < k; j++, offset += d)
        {
            double updated = gamma * counts[j] + alpha * batchCounts[j];
            if (updated >= MassFloor)
            {
                counts[j] = updated;
                for (int c = offset; c < offset + d; c++) sums[c] = gamma * sums[c] + alpha * batchSums[c];
                if (scatter != null) scatter[j] = gamma * scatter[j] + alpha * batchScatter[j];
            }
            else if (updated > 0.0)
            {
                for (int c = offset; c < offset + d; c++)
                    sums[c] = MassFloor * ((gamma * sums[c] + alpha * batchSums[c]) / updated);
                if (scatter != null)
                    scatter[j] = MassFloor * ((gamma * scatter[j] + alpha * batchScatter[j]) / updated);
                counts[j] = MassFloor;
            }
            else
            {
                double previous = counts[j];
                for (int c = offset; c < offset + d; c++) sums[c] = MassFloor * (sums[c] / previous);
                if (scatter != null) scatter[j] = MassFloor * (scatter[j] / previous);
                counts[j] = MassFloor;
            }
        }
    }

    // Revives each cluster whose count is below threshold times the mean count: a donor drawn without replacement
    // among the live clusters, with probability proportional to its count, splits into two children at plus and
    // minus one offset from its centroid, each taking half its count and scatter. The offset runs along a random
    // direction for perturbation times the donor's root-mean-square radius. The slow statistics split the same slots.
    private static int Revive(
        double[] counts, double[] sums, double[] scatter, double[] slowCounts, double[] slowSums, int k, int d,
        double threshold, double perturbation, Xoshiro256PlusPlus rng,
        double[] donorWeights, double[] direction, bool lastEpoch, int[] lastCounts, double[] lastSums)
    {
        double total = 0.0;
        for (int j = 0; j < k; j++) total += counts[j];
        double cut = threshold * total / k;

        int dead = 0;
        double pool = 0.0;
        for (int j = 0; j < k; j++)
        {
            if (counts[j] < cut)
            {
                donorWeights[j] = -1.0;
                dead++;
            }
            else
            {
                donorWeights[j] = counts[j];
                pool += counts[j];
            }
        }
        if (dead == 0) return 0;

        int revived = 0;
        for (int j = 0; j < k; j++)
        {
            if (donorWeights[j] >= 0.0) continue;

            double target = rng.NextDouble() * pool;
            double cumulative = 0.0;
            int donor = -1;
            int lastLive = -1;
            for (int i = 0; i < k; i++)
            {
                if (!(donorWeights[i] > 0.0)) continue;
                lastLive = i;
                cumulative += donorWeights[i];
                if (cumulative > target)
                {
                    donor = i;
                    break;
                }
            }
            if (donor < 0) donor = lastLive;
            if (donor < 0) break; // every live cluster has already given half this round

            pool -= donorWeights[donor];
            donorWeights[donor] = 0.0;

            double squaredNorm = 0.0;
            for (int c = 0; c < d; c++)
            {
                double u = 2.0 * rng.NextDouble() - 1.0;
                direction[c] = u;
                squaredNorm += u * u;
            }
            if (!(squaredNorm > 0.0))
            {
                Array.Clear(direction);
                direction[0] = 1.0;
                squaredNorm = 1.0;
            }
            double radius = Math.Sqrt(scatter[donor] / counts[donor]);
            double scale = perturbation * radius / Math.Sqrt(squaredNorm);
            for (int c = 0; c < d; c++) direction[c] *= scale;

            Split(counts, sums, donor, j, d, direction);
            Split(slowCounts, slowSums, donor, j, d, direction);
            scatter[donor] *= 0.5;
            scatter[j] = scatter[donor];

            if (lastEpoch)
            {
                // Records made before the split describe clusters that no longer exist.
                lastCounts[donor] = 0;
                lastCounts[j] = 0;
                Array.Clear(lastSums, donor * d, d);
                Array.Clear(lastSums, j * d, d);
            }

            donorWeights[j] = 0.0;
            revived++;
        }
        return revived;
    }

    private static void Split(double[] counts, double[] sums, int donor, int dead, int d, double[] offset)
    {
        double whole = counts[donor];
        double half = 0.5 * whole;
        int from = donor * d;
        int to = dead * d;
        for (int c = 0; c < d; c++)
        {
            double centre = sums[from + c] / whole;
            sums[to + c] = half * (centre + offset[c]);
            sums[from + c] = half * (centre - offset[c]);
        }
        counts[donor] = half;
        counts[dead] = half;
    }

    private static int Nearest(double[] x, double[] centroids, int k, int d, out double best)
    {
        int nearest = 0;
        best = double.PositiveInfinity;
        for (int j = 0, offset = 0; j < k; j++, offset += d)
        {
            double sum = 0.0;
            for (int c = 0; c < d; c++)
            {
                double diff = x[c] - centroids[offset + c];
                sum += diff * diff;
            }
            if (sum < best)
            {
                best = sum;
                nearest = j;
            }
        }
        return nearest;
    }

    private static void Shuffle(int[] perm, Xoshiro256PlusPlus rng)
    {
        for (int i = 0; i < perm.Length; i++) perm[i] = i;
        for (int i = perm.Length - 1; i > 0; i--)
        {
            int other = UniformIndex(rng, i + 1);
            (perm[i], perm[other]) = (perm[other], perm[i]);
        }
    }

    // Uniform index in [0, bound) from the high half of a 64-bit product; the bias, bound/2^64, is far below any
    // n a run can hold.
    private static int UniformIndex(Xoshiro256PlusPlus rng, int bound)
        => (int)Math.BigMul(rng.NextUInt64(), (ulong)bound, out _);

    private static double[][] Ratios(double[] sums, double[] counts, int k, int d)
    {
        var rows = new double[k][];
        for (int j = 0; j < k; j++)
        {
            var row = new double[d];
            for (int c = 0; c < d; c++) row[c] = sums[j * d + c] / counts[j];
            rows[j] = row;
        }
        return rows;
    }

    private static double[][] Rows(double[] flat, int k, int d)
    {
        var rows = new double[k][];
        for (int j = 0; j < k; j++)
        {
            rows[j] = new double[d];
            Array.Copy(flat, j * d, rows[j], 0, d);
        }
        return rows;
    }

    private static double[][] Copy(double[][] rows)
    {
        var copy = new double[rows.Length][];
        for (int j = 0; j < rows.Length; j++) copy[j] = (double[])rows[j].Clone();
        return copy;
    }

    private static int ValidateRows(double[][] rows, string name)
    {
        if (rows[0] == null || rows[0].Length == 0)
            throw new ArgumentException("Rows must have at least one coordinate.", name);
        int d = rows[0].Length;
        for (int i = 0; i < rows.Length; i++)
        {
            double[] row = rows[i];
            if (row == null || row.Length != d)
                throw new ArgumentException($"Row {i} does not have {d} coordinates.", name);
            for (int c = 0; c < d; c++)
            {
                if (!double.IsFinite(row[c]))
                    throw new ArgumentException($"Row {i} holds a non-finite value.", name);
            }
        }
        return d;
    }

    private static void ValidateConfig(EMeansConfig config, int n)
    {
        if (config.K < 1 || config.K > n)
            throw new ArgumentException("K must be between 1 and n.", nameof(config));
        if (config.Epochs < 1)
            throw new ArgumentException("Epochs must be at least 1.", nameof(config));
        if (!(config.InitialCount >= MassFloor) || !double.IsFinite(config.InitialCount))
            throw new ArgumentException("InitialCount must be finite and at least 1e-200.", nameof(config));
        if (!(config.Step.FloorFraction >= 0.0 && config.Step.FloorFraction <= 1.0))
            throw new ArgumentException("Step.FloorFraction must lie in [0, 1].", nameof(config));
        if (!(config.SlowStepMultiplier > 0.0) || !double.IsFinite(config.SlowStepMultiplier))
            throw new ArgumentException("SlowStepMultiplier must be positive and finite.", nameof(config));
        if (config.Revival.Period < 1)
            throw new ArgumentException("Revival.Period must be at least 1.", nameof(config));
        if (!(config.Revival.Threshold >= 0.0 && config.Revival.Threshold < 1.0))
            throw new ArgumentException("Revival.Threshold must lie in [0, 1).", nameof(config));
        if (!(config.Revival.Perturbation > 0.0) || !double.IsFinite(config.Revival.Perturbation))
            throw new ArgumentException("Revival.Perturbation must be positive and finite.", nameof(config));
    }
}
