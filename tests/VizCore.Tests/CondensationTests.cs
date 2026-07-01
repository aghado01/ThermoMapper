using System;
using System.Linq;
using Clustering.Dendrograms;
using Clustering.Graphical.HdbScan;
using Clustering.Primitives;
using Graphs.Distance;
using Graphs.Distance.Euclidean;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// Unification cross-check: the shared producer-agnostic condensation + EOM
/// (Clustering.Dendrograms.Condensation) reproduces HdbscanRunner's own labels
/// EXACTLY when fed the runner's emitted dendrogram with the same minClusterSize.
/// Validation independence: the oracle is HDBSCAN's bespoke ExtractClusters
/// (a separate implementation), not the shared path's own assumptions — so a
/// match proves the shared spine subsumes the in-silo extraction. Landed
/// alongside the runner; the rewire/rip-out is fresh-look cleanup (ledger #7).
/// </summary>
public sealed class CondensationTests
{
    // Three Gaussian-ish blobs in 2-D + a couple of stragglers, deterministic.
    private static (double[] data, int dim, int n) ThreeBlobs()
    {
        var pts = new System.Collections.Generic.List<double>();
        void Blob(double cx, double cy, int count, int seedMix)
        {
            // Deterministic pseudo-jitter (no RNG dependency in the test).
            for (int k = 0; k < count; k++)
            {
                double a = ((k * 2654435761u + (uint)seedMix) % 1000) / 1000.0 - 0.5;
                double b = ((k * 40503u + (uint)(seedMix * 7)) % 1000) / 1000.0 - 0.5;
                pts.Add(cx + a); pts.Add(cy + b);
            }
        }
        Blob(0, 0, 20, 1);
        Blob(8, 0, 20, 2);
        Blob(4, 7, 20, 3);
        int n = pts.Count / 2;
        return (pts.ToArray(), 2, n);
    }

    [Theory]
    [InlineData(3, 3, true)]
    [InlineData(5, 5, true)]
    [InlineData(3, 5, false)]
    [InlineData(4, 4, false)]
    public void SharedCondensation_ReproducesHdbscanLabels(int minPts, int minClusterSize, bool allowSingle)
    {
        var (data, dim, n) = ThreeBlobs();

        var runner = new HdbscanRunner(n);
        HdbscanResult oracle = runner.Run(
            data, dim, minPts, new EuclideanMetric(),
            minClusterSize: minClusterSize, allowSingleCluster: allowSingle);

        // Shared path over the runner's OWN emitted dendrogram, same minClusterSize.
        CondensedTree condensed = Condensation.Condense(oracle.Dendrogram, minClusterSize);
        bool[] selected = condensed.SelectByExcessOfMass(allowSingle);
        Assignment shared = condensed.ToAssignment(selected);

        // Cluster COUNT matches.
        Assert.Equal(oracle.ClusterCount, shared.Count);

        // Labels match up to a relabeling (both dense, noise = -1 / Unassigned).
        AssertSameClustering(oracle.Labels, shared.Labels, n);
    }

    private static void AssertSameClustering(int[] expected, int[] actual, int n)
    {
        Assert.Equal(n, expected.Length);
        Assert.Equal(n, actual.Length);
        var fwd = new System.Collections.Generic.Dictionary<int, int>();
        var rev = new System.Collections.Generic.Dictionary<int, int>();
        for (int i = 0; i < n; i++)
        {
            bool en = expected[i] < 0, an = actual[i] < 0;
            Assert.Equal(en, an); // noise iff noise
            if (en) continue;
            if (fwd.TryGetValue(expected[i], out int mapped)) Assert.Equal(mapped, actual[i]);
            else fwd[expected[i]] = actual[i];
            if (rev.TryGetValue(actual[i], out int back)) Assert.Equal(back, expected[i]);
            else rev[actual[i]] = expected[i];
        }
    }

    [Fact]
    public void Condense_RejectsSubMinClusterSize()
    {
        var (data, dim, n) = ThreeBlobs();
        HdbscanResult r = new HdbscanRunner(n).Run(data, dim, 3, new EuclideanMetric());
        Assert.Throws<ArgumentOutOfRangeException>(() => Condensation.Condense(r.Dendrogram, minClusterSize: 1));
    }

    /// <summary>
    /// Leaf selection picks exactly the condensed-tree leaves (clusters with no
    /// condensed child), and is never coarser than EOM: leaf ≥ EOM on cluster
    /// count, because wherever EOM selects an internal cluster, leaf selects the
    /// ≥2 leaves beneath it. This is the finer half of the shared selector axis
    /// (the behavioural "recovers more multi-class structure" proof is the
    /// Landsat control fact).
    /// </summary>
    [Theory]
    [InlineData(3, 3, false)]
    [InlineData(5, 5, false)]
    [InlineData(3, 3, true)]
    public void SelectByLeaf_PicksCondensedLeaves_AndIsNeverCoarserThanEom(int minPts, int minClusterSize, bool allowSingle)
    {
        var (data, dim, n) = ThreeBlobs();
        HdbscanResult r = new HdbscanRunner(n).Run(
            data, dim, minPts, new EuclideanMetric(),
            minClusterSize: minClusterSize, allowSingleCluster: allowSingle);

        CondensedTree condensed = Condensation.Condense(r.Dendrogram, minClusterSize);
        bool[] eom  = condensed.SelectByExcessOfMass(allowSingle);
        bool[] leaf = condensed.SelectByLeaf(allowSingle);

        // Every leaf-selected cluster is a true condensed-tree leaf (no child
        // points to it as parent).
        var isParent = new bool[condensed.ClusterCount];
        for (int c = 1; c < condensed.ClusterCount; c++)
            if (condensed.Parent[c] >= 0) isParent[condensed.Parent[c]] = true;
        for (int c = 0; c < condensed.ClusterCount; c++)
            if (leaf[c]) Assert.False(isParent[c], $"leaf-selected cluster {c} has condensed children");

        // Leaf is never coarser than EOM.
        Assert.True(Count(leaf) >= Count(eom),
            $"leaf count {Count(leaf)} should be >= eom count {Count(eom)}");

        static int Count(bool[] s) { int k = 0; foreach (bool b in s) if (b) k++; return k; }
    }

    /// <summary>
    /// Membership probabilities from <see cref="CondensedTree.ResolveLabeled"/>
    /// honour the HDBSCAN contract: all in [0,1]; noise (label −1) ⇒ prob 0;
    /// assigned points ⇒ prob &gt; 0. Guards the λ-ratio formula ported onto the
    /// spine during the runner rewire.
    /// </summary>
    [Fact]
    public void ResolveLabeled_MembershipProbabilities_HonourContract()
    {
        var (data, dim, n) = ThreeBlobs();
        HdbscanResult r = new HdbscanRunner(n).Run(
            data, dim, 5, new EuclideanMetric(), minClusterSize: 5, allowSingleCluster: false);

        for (int i = 0; i < n; i++)
        {
            double p = r.MembershipProbabilities[i];
            Assert.InRange(p, 0.0, 1.0);
            if (r.Labels[i] < 0) Assert.Equal(0.0, p);
            else Assert.True(p > 0.0, $"assigned point {i} (label {r.Labels[i]}) has prob {p}");
        }
    }

    /// <summary>
    /// The sparse-kNN MST path (<see cref="MstAlgorithm.SparseKnn"/>): with full
    /// neighbours (graphK ≥ n−1, clamped) every pair is a candidate edge, so the
    /// kNN-restricted mutual-reachability MST IS the dense MST — labels match
    /// exactly. With small k it still recovers well-separated blobs (the kNN graph
    /// is intra-blob-connected; Boruvka bridges the blobs at the same closest pairs
    /// the dense MST uses). Validation-independent: the dense path is the oracle,
    /// and the k=n−1 case is exact equivalence, not parity to an external library.
    /// </summary>
    [Theory]
    [InlineData(5, 5, false, 1000)]   // full neighbours (clamped to n-1) ⇒ exact dense MST
    [InlineData(3, 3, true,  1000)]
    [InlineData(5, 5, false, 10)]     // small k ⇒ still recovers the separated blobs
    public void SparseKnn_MatchesDense_OnSeparatedBlobs(int minPts, int mcs, bool allowSingle, int graphK)
    {
        var (data, dim, n) = ThreeBlobs();

        HdbscanResult dense = new HdbscanRunner(n).Run(
            data, dim, minPts, new EuclideanMetric(), minClusterSize: mcs, allowSingleCluster: allowSingle);
        HdbscanResult sparse = new HdbscanRunner(n).RunSparse(
            data, dim, minPts, new EuclideanMetric(), graphK, minClusterSize: mcs, allowSingleCluster: allowSingle);

        Assert.Equal(dense.ClusterCount, sparse.ClusterCount);
        AssertSameClustering(dense.Labels, sparse.Labels, n);
    }

    /// <summary>
    /// <see cref="CondensedTree.SelectByEpsilon"/> (sklearn's
    /// <c>cluster_selection_epsilon</c>): a distance below every cluster's birth
    /// is a no-op (nothing merges); a distance above every birth merges the whole
    /// base selection up to the root (one cluster, allowSingleCluster). Validates
    /// the traverse-upwards merge against the tree's own birth distances — no
    /// external oracle.
    /// </summary>
    [Fact]
    public void SelectByEpsilon_HugeMergesToRoot_TinyIsNoOp()
    {
        var (data, dim, n) = ThreeBlobs();
        HdbscanResult r = new HdbscanRunner(n).Run(
            data, dim, 5, new EuclideanMetric(), minClusterSize: 5, allowSingleCluster: true);
        CondensedTree condensed = Condensation.Condense(r.Dendrogram, 5);
        bool[] baseLeaf = condensed.SelectByLeaf(allowSingleCluster: true);
        Assert.True(CountTrue(baseLeaf) >= 2, "leaf base selection should find multiple blobs");

        bool[] tiny = condensed.SelectByEpsilon(baseLeaf, 1e-9, allowSingleCluster: true);
        Assert.Equal(CountTrue(baseLeaf), CountTrue(tiny));   // nothing born below 1e-9 ⇒ no merge

        bool[] huge = condensed.SelectByEpsilon(baseLeaf, 1e9, allowSingleCluster: true);
        Assert.Equal(1, CountTrue(huge));                     // everything merges up to the root
        Assert.True(huge[0], "root (id 0) selected after the huge-epsilon merge");

        static int CountTrue(bool[] s) { int k = 0; foreach (bool b in s) if (b) k++; return k; }
    }

    /// <summary>
    /// End-to-end through the runner: epsilon coarsens leaf selection. A huge
    /// epsilon collapses the blobs to a single cluster (vs leaf's several),
    /// proving the knob threads settings → selection and merges as intended.
    /// </summary>
    [Fact]
    public void ClusterSelectionEpsilon_CoarsensLeaf_EndToEnd()
    {
        var (data, dim, n) = ThreeBlobs();
        int K(double eps) => new HdbscanRunner(n).Run(
            data, dim, 5, new EuclideanMetric(),
            minClusterSize: 5, allowSingleCluster: true,
            selectionMethod: ClusterSelectionMethod.Leaf, clusterSelectionEpsilon: eps).ClusterCount;

        int k0 = K(0.0);
        Assert.True(k0 >= 2, $"leaf should find multiple blobs, got {k0}");
        Assert.True(K(1e9) < k0, "a huge epsilon must merge clusters below the base count");
        Assert.Equal(1, K(1e9));
    }
}
