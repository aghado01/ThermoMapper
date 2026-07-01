using System;

namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// Mode-tree merging: components whose means ascend to the same local mode of
    /// the mixture density are assigned to the same cluster. Geometry-only — does
    /// not consume the responsibility matrix. See docs/gmm.md.
    /// </summary>
    public sealed class ModalMergeStrategy : IComponentMergeStrategy
    {
        /// <summary>
        /// Euclidean distance within which two converged modes are treated as one.
        /// <see cref="double.NaN"/> (default) selects an adaptive threshold equal to
        /// 5% of the mean pairwise distance between component means.
        /// </summary>
        public double ModeTolerance { get; }

        /// <summary>Iteration cap forwarded to <c>ModeAscent.Ascend</c>.</summary>
        public int MaxAscentSteps { get; }

        /// <summary>Displacement tolerance forwarded to <c>ModeAscent.Ascend</c>.</summary>
        public double AscentTol { get; }

        public ModalMergeStrategy(
            double modeTolerance = double.NaN,
            int maxAscentSteps = 50,
            double ascentTol = 1e-7)
        {
            ModeTolerance = modeTolerance;
            MaxAscentSteps = maxAscentSteps;
            AscentTol = ascentTol;
        }

        /// <inheritdoc />
        public int[] Merge(GaussianComponent[] components, double[,]? responsibilities = null)
        {
            int k = components.Length;
            if (k == 0) return Array.Empty<int>();
            if (k == 1) return new[] { 0 };

            // ── Ascend from each component mean (single shared scratch) ──────────
            int d = components[0].Mean.Length;
            var scratch = new ModeAscentScratch(d, k);
            double[][] modes = new double[k][];
            for (int i = 0; i < k; i++)
            {
                modes[i] = new double[d];
                ModeAscent.Ascend(components[i].Mean, components, modes[i], scratch,
                                  MaxAscentSteps, AscentTol);
            }

            // ── Resolve tolerance ─────────────────────────────────────────────────
            double threshold = double.IsNaN(ModeTolerance)
                ? ComputeAdaptiveTolerance(components, k)
                : ModeTolerance;

            // ── Sequential grouping by mode proximity ──────────────────────────────
            // clusterModes[c] is the representative mode for cluster c.
            // Allocation: at most k clusters (all components distinct).
            int[] map = new int[k];
            double[][] clusterModes = new double[k][];
            int clusterCount = 0;

            for (int i = 0; i < k; i++)
            {
                int assigned = -1;
                for (int c = 0; c < clusterCount; c++)
                {
                    if (EuclideanDistance(modes[i], clusterModes[c]) < threshold)
                    {
                        assigned = c;
                        break;
                    }
                }

                if (assigned >= 0)
                {
                    map[i] = assigned;
                }
                else
                {
                    clusterModes[clusterCount] = modes[i];
                    map[i] = clusterCount++;
                }
            }

            return map;
        }

        // ── Internals ──────────────────────────────────────────────────────────────

        // 5% of the mean pairwise distance between component means — scale-free.
        private static double ComputeAdaptiveTolerance(GaussianComponent[] components, int k)
        {
            double sum = 0.0;
            int pairs = 0;
            for (int i = 0; i < k; i++)
                for (int j = i + 1; j < k; j++)
                {
                    sum += EuclideanDistance(components[i].Mean, components[j].Mean);
                    pairs++;
                }
            return sum / pairs * 0.05;
        }

        private static double EuclideanDistance(double[] a, double[] b)
        {
            double sum = 0.0;
            for (int i = 0; i < a.Length; i++)
            {
                double d = a[i] - b[i];
                sum += d * d;
            }
            return Math.Sqrt(sum);
        }
    }
}
