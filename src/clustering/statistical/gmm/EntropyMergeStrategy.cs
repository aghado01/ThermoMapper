using System;
using System.Collections.Generic;
using Clustering.Dendrograms;

namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// One level of the merge dendrogram produced by <see cref="EntropyMergeStrategy"/>:
    /// the cluster count at this level, the classification entropy
    /// H = −Σᵢ Σₖ r̂ᵢₖ log r̂ᵢₖ of the pooled responsibilities, and the dense
    /// component→cluster map.
    /// </summary>
    public sealed record MergeStep(
        int ClusterCount,
        double ClassificationEntropy,
        int[] ComponentToClusterMap);

    /// <summary>
    /// Agglomerative component merging by classification entropy (Baudry et al.,
    /// JCGS 19(2), 2010). Operates on the final responsibility matrix — no
    /// re-fitting. See docs/gmm.md.
    /// </summary>
    public sealed class EntropyMergeStrategy : IComponentMergeStrategy
    {
        /// <summary>
        /// Stops merging at this many clusters. Defaults to 1 (fully collapsed);
        /// set to a BIC-optimal K for use as a post-hoc assignment step.
        /// </summary>
        public int TargetClusters { get; }

        public EntropyMergeStrategy(int targetClusters = 1)
        {
            if (targetClusters < 1)
                throw new ArgumentOutOfRangeException(nameof(targetClusters), "Must be ≥ 1.");
            TargetClusters = targetClusters;
        }

        /// <inheritdoc />
        /// <exception cref="ArgumentNullException">
        /// Thrown when <paramref name="responsibilities"/> is <c>null</c> — this
        /// strategy requires the responsibility matrix.
        /// </exception>
        public int[] Merge(GaussianComponent[] components, double[,]? responsibilities = null)
        {
            if (responsibilities is null)
                throw new ArgumentNullException(nameof(responsibilities),
                    $"{nameof(EntropyMergeStrategy)} requires the responsibility matrix.");

            var steps = MergeSequence(components, responsibilities);
            MergeStep target = steps[steps.Length - 1];
            foreach (var s in steps)
            {
                if (s.ClusterCount <= TargetClusters) { target = s; break; }
            }
            return target.ComponentToClusterMap;
        }

        /// <summary>
        /// Full merge dendrogram from K clusters (entry 0, identity map) down to 1
        /// (entry K−1). Returns one <see cref="MergeStep"/> per level.
        /// </summary>
        public static MergeStep[] MergeSequence(GaussianComponent[] components, double[,] responsibilities)
        {
            int origK = components.Length;
            int n = responsibilities.GetLength(0);

            // Pooled responsibilities at full origK width; liveness in `alive[]`
            // avoids O(K·N) memmove on every merge.
            double[][] r = new double[n][];
            for (int i = 0; i < n; i++)
            {
                r[i] = new double[origK];
                for (int g = 0; g < origK; g++)
                    r[i][g] = responsibilities[i, g];
            }

            var groupIndices = new List<int>[origK];
            for (int g = 0; g < origK; g++) groupIndices[g] = new List<int> { g };

            bool[] alive = new bool[origK];
            for (int g = 0; g < origK; g++) alive[g] = true;

            // ΔH cache: built once at O(K²·N); after each merge only the row/column
            // touching the surviving primary is refreshed at O(K·N). Total: O(K²·N).
            double[,] delta = new double[origK, origK];
            for (int j = 0; j < origK; j++)
                for (int l = j + 1; l < origK; l++)
                    delta[j, l] = PairDelta(r, j, l, n);

            var steps = new MergeStep[origK];
            steps[0] = new MergeStep(
                ClusterCount: origK,
                ClassificationEntropy: ComputeEntropy(r, alive, n, origK),
                ComponentToClusterMap: BuildMap(groupIndices, alive, origK));

            int curK = origK;

            for (int step = 1; step < origK; step++)
            {
                // Find the alive pair with minimum ΔH (cache scan, no N-work).
                int bestJ = -1, bestL = -1;
                double bestDelta = double.MaxValue;
                for (int j = 0; j < origK; j++)
                {
                    if (!alive[j]) continue;
                    for (int l = j + 1; l < origK; l++)
                    {
                        if (!alive[l]) continue;
                        if (delta[j, l] < bestDelta) { bestDelta = delta[j, l]; bestJ = j; bestL = l; }
                    }
                }

                for (int i = 0; i < n; i++) r[i][bestJ] += r[i][bestL];
                groupIndices[bestJ].AddRange(groupIndices[bestL]);
                alive[bestL] = false;
                curK--;

                // Refresh only the ΔH entries touching the surviving primary.
                for (int other = 0; other < origK; other++)
                {
                    if (other == bestJ || !alive[other]) continue;
                    int a = Math.Min(bestJ, other);
                    int b = Math.Max(bestJ, other);
                    delta[a, b] = PairDelta(r, a, b, n);
                }

                steps[step] = new MergeStep(
                    ClusterCount: curK,
                    ClassificationEntropy: ComputeEntropy(r, alive, n, origK),
                    ComponentToClusterMap: BuildMap(groupIndices, alive, origK));
            }

            return steps;
        }

        /// <summary>
        /// The agglomeration as a shared <see cref="Dendrogram"/> — leaves are the
        /// K components, each merge node's height is the cumulative classification-
        /// entropy reduction H(K) − H(level) (monotone non-decreasing, so the
        /// dendrogram build invariant holds). Lets the merge sequence feed the
        /// shared resolution spine (condensation / cut / walk) instead of the
        /// flattening <see cref="Merge"/> → <c>int[]</c>. Cost axis
        /// <c>"entropy_reduction"</c>.
        /// </summary>
        /// <remarks>
        /// GMM is the third <see cref="Dendrogram"/> producer (HDBSCAN condensed,
        /// PKWang/SW thermal single-linkage, GMM agglomerative). Landed ALONGSIDE
        /// the flatten path; rewiring <see cref="Merge"/> onto a shared cut is
        /// fresh-look cleanup (post-unification ledger).
        /// </remarks>
        public static Dendrogram BuildDendrogram(GaussianComponent[] components, double[,] responsibilities)
        {
            ArgumentNullException.ThrowIfNull(components);
            ArgumentNullException.ThrowIfNull(responsibilities);
            int origK = components.Length;
            if (origK < 2)
                throw new ArgumentOutOfRangeException(nameof(components), "Need at least 2 components for a dendrogram.");
            int n = responsibilities.GetLength(0);

            double[][] r = new double[n][];
            for (int i = 0; i < n; i++)
            {
                r[i] = new double[origK];
                for (int g = 0; g < origK; g++) r[i][g] = responsibilities[i, g];
            }

            bool[] alive = new bool[origK];
            for (int g = 0; g < origK; g++) alive[g] = true;

            // currentNodeId[g] = the dendrogram node id currently representing
            // alive group g (a leaf id, then internal ids as it absorbs merges).
            int[] currentNodeId = new int[origK];
            for (int g = 0; g < origK; g++) currentNodeId[g] = g;

            double[,] delta = new double[origK, origK];
            for (int j = 0; j < origK; j++)
                for (int l = j + 1; l < origK; l++)
                    delta[j, l] = PairDelta(r, j, l, n);

            double baseEntropy = ComputeEntropy(r, alive, n, origK);
            var merges = new DendrogramNode[origK - 1];
            int nextId = origK;

            for (int step = 1; step < origK; step++)
            {
                int bestJ = -1, bestL = -1;
                double bestDelta = double.MaxValue;
                for (int j = 0; j < origK; j++)
                {
                    if (!alive[j]) continue;
                    for (int l = j + 1; l < origK; l++)
                    {
                        if (!alive[l]) continue;
                        if (delta[j, l] < bestDelta) { bestDelta = delta[j, l]; bestJ = j; bestL = l; }
                    }
                }

                for (int i = 0; i < n; i++) r[i][bestJ] += r[i][bestL];
                alive[bestL] = false;

                for (int other = 0; other < origK; other++)
                {
                    if (other == bestJ || !alive[other]) continue;
                    int a = Math.Min(bestJ, other);
                    int b = Math.Max(bestJ, other);
                    delta[a, b] = PairDelta(r, a, b, n);
                }

                int sizeJ = SubtreeSize(currentNodeId[bestJ], origK, merges);
                int sizeL = SubtreeSize(currentNodeId[bestL], origK, merges);
                double height = baseEntropy - ComputeEntropy(r, alive, n, origK);
                merges[step - 1] = new DendrogramNode(
                    currentNodeId[bestJ], currentNodeId[bestL], height, sizeJ + sizeL);

                currentNodeId[bestJ] = nextId++;
            }

            return new Dendrogram(merges, origK, CostAxis: "entropy_reduction");
        }

        private static int SubtreeSize(int nodeId, int leafCount, DendrogramNode[] merges)
            => nodeId < leafCount ? 1 : merges[nodeId - leafCount].Size;

        // ── Internals ──────────────────────────────────────────────────────────────

        private static double XlogX(double v) => v < 1e-300 ? 0.0 : v * Math.Log(v);

        /// <summary>
        /// ΔH(j, l) = Σᵢ [ XlogX(rᵢⱼ) + XlogX(rᵢₗ) − XlogX(rᵢⱼ + rᵢₗ) ] ≤ 0.
        /// O(N) — the per-step kernel.
        /// </summary>
        private static double PairDelta(double[][] r, int j, int l, int n)
        {
            double d = 0.0;
            for (int i = 0; i < n; i++)
            {
                double merged = r[i][j] + r[i][l];
                d -= XlogX(merged) - XlogX(r[i][j]) - XlogX(r[i][l]);
            }
            return d;
        }

        private static double ComputeEntropy(double[][] pooled, bool[] alive, int nPts, int origK)
        {
            double h = 0.0;
            for (int i = 0; i < nPts; i++)
                for (int g = 0; g < origK; g++)
                    if (alive[g]) h -= XlogX(pooled[i][g]);
            return h;
        }

        /// <summary>
        /// Builds the dense component→cluster map. Cluster indices are assigned in
        /// order of first appearance among the alive primaries, giving a stable
        /// dense range [0, curK).
        /// </summary>
        private static int[] BuildMap(List<int>[] groupIndices, bool[] alive, int origK)
        {
            int[] map = new int[origK];
            int clusterId = 0;
            for (int g = 0; g < origK; g++)
            {
                if (!alive[g]) continue;
                foreach (int idx in groupIndices[g]) map[idx] = clusterId;
                clusterId++;
            }
            return map;
        }
    }
}
