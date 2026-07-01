using System;
using System.Collections.Generic;
using Clustering.Primitives;

namespace Clustering.Dendrograms;

/// <summary>
/// Per-internal-node walk results over a (dendrogram, landscape) pair: the
/// lifetime-integrated landscape mass (the generalized excess-of-mass
/// quantity) and the peak of the cluster's density profile. Indexed by merge
/// index (internal node id − LeafCount).
/// </summary>
/// <param name="Mass">∑_g S(C,g)·Δ_g over the cluster's lifetime — at L ≡ 1
/// this is exactly |C| × lifetime, the classic HDBSCAN stability.</param>
/// <param name="PeakGridIndex">argmax_g S(C,g)/|C| within the lifetime;
/// −1 when no grid cell falls inside it.</param>
/// <param name="PeakDensity">The density at the peak; NaN when no cell.</param>
/// <param name="Birth">The merge height that created the cluster.</param>
/// <param name="Death">The parent merge height; +∞ for a root.</param>
public sealed record ClusterWalkReport(
    double[] Mass,
    int[]    PeakGridIndex,
    double[] PeakDensity,
    double[] Birth,
    double[] Death);

/// <summary>
/// The selector-side walk: accumulates each cluster's landscape profile
/// S(C,·) = Σ_{p∈C} L(p,·) incrementally up the merge tree (vector-summed
/// union of child profiles — additivity of the per-leaf primitive is the
/// licence), then reduces it two ways: <b>integrate</b> over the lifetime
/// (generalized excess-of-mass) and <b>peak</b> (wave_clus's select move).
/// </summary>
/// <remarks>
/// <para><b>Axis-alignment law.</b> A walk is only defined when structure and
/// height share the axis: <see cref="Dendrogram.CostAxis"/> must equal
/// <see cref="Landscape.Axis"/> (ordinal comparison); throws otherwise.</para>
///
/// <para><b>Discretization contract (v1).</b> Left-point Riemann cells:
/// Δ_g = Grid[g+1] − Grid[g], with the last width replicated. Leaves carry
/// no mass by convention (their resolution is periphery-policy territory;
/// condensation is a producer-side concern).</para>
///
/// <para><b>Orientation.</b> Agglomerative heights may ascend (distance-like
/// axes: leaves merge upward) or DESCEND (thermal: hot singletons couple as
/// T falls — <see cref="ThermalDendrogram"/>). Detected from the merge
/// sequence; lifetime windows follow it: ascending → [birth, death) with the
/// root open at the large end; descending → (death, birth] with the root
/// open at the small (cold) end. A single-merge tree defaults to ascending.</para>
///
/// <para><b>Cardinal consumer.</b> Masses depend on the landscape's gauge —
/// declared in <see cref="Landscape.Provenance"/>; rankings are
/// sink-relative by design.</para>
/// </remarks>
public static class LandscapeWalk
{
    /// <summary>
    /// Computes each internal node's lifetime-integrated mass and profile
    /// peak via one pass over the merges (vector-summed union-find;
    /// small-into-large is implicit in the fixed merge order).
    /// </summary>
    public static ClusterWalkReport ClusterProfiles(Dendrogram dendrogram, Landscape landscape)
    {
        ArgumentNullException.ThrowIfNull(dendrogram);
        ArgumentNullException.ThrowIfNull(landscape);
        if (!string.Equals(dendrogram.CostAxis, landscape.Axis, StringComparison.Ordinal))
            throw new InvalidOperationException(
                $"Axis-alignment law: dendrogram cost axis '{dendrogram.CostAxis}' != landscape axis " +
                $"'{landscape.Axis}'. Walks are only defined when structure and height share the axis.");
        if (landscape.NodeCount != dendrogram.LeafCount)
            throw new InvalidOperationException(
                $"Landscape node count ({landscape.NodeCount}) != dendrogram leaf count ({dendrogram.LeafCount}).");

        int n = dendrogram.LeafCount;
        int m = dendrogram.InternalNodeCount;
        int gridCount = landscape.Grid.Length;
        DendrogramNode[] merges = dendrogram.Merges;
        double[] grid = landscape.Grid;
        double[][] columns = landscape.ValuesByGridPoint;

        // Left-point Riemann widths; last width replicated.
        var width = new double[gridCount];
        for (int g = 0; g < gridCount - 1; g++) width[g] = grid[g + 1] - grid[g];
        width[gridCount - 1] = gridCount > 1 ? width[gridCount - 2] : 1.0;

        var mass    = new double[m];
        var peakIdx = new int[m];
        var peakDen = new double[m];
        var birth   = new double[m];
        var death   = new double[m];
        var profile = new double[m][];   // live until the parent consumes it

        // Orientation: descending iff the merge heights predominantly fall in
        // build order (thermal trees: hot couplings first).
        int up = 0, down = 0;
        for (int i = 1; i < m; i++)
        {
            if (merges[i].Distance > merges[i - 1].Distance) up++;
            else if (merges[i].Distance < merges[i - 1].Distance) down++;
        }
        bool descending = down > up;

        void Finalize(int k, double deathHeight)
        {
            double[] s = profile[k]!;
            double b = birth[k];
            int size = merges[k].Size;
            double sum = 0.0;
            int    pi  = -1;
            double pd  = double.NegativeInfinity;
            for (int g = 0; g < gridCount; g++)
            {
                double t = grid[g];
                bool inside = descending
                    ? t > deathHeight && t <= b    // (death, birth] — thermal: alive below creation T
                    : t >= b && t < deathHeight;   // [birth, death) — distance-like
                if (!inside) continue;
                sum += s[g] * width[g];
                double density = s[g] / size;
                if (density > pd) { pd = density; pi = g; }
            }
            mass[k]    = sum;
            peakIdx[k] = pi;
            peakDen[k] = pi < 0 ? double.NaN : pd;
            death[k]   = deathHeight;
        }

        double[] Acquire(int id, double parentHeight)
        {
            if (id < n)
            {
                var v = new double[gridCount];
                for (int g = 0; g < gridCount; g++) v[g] = columns[g][id];
                return v;
            }
            int k = id - n;
            Finalize(k, parentHeight);
            double[] s = profile[k]!;
            profile[k] = null!;
            return s;
        }

        for (int i = 0; i < m; i++)
        {
            DendrogramNode node = merges[i];
            double h = node.Distance;
            double[] left  = Acquire(node.LeftChild, h);
            double[] right = Acquire(node.RightChild, h);
            for (int g = 0; g < gridCount; g++) left[g] += right[g];
            profile[i] = left;
            birth[i] = h;
        }

        // Never-consumed internal nodes are roots: the open end of the lifetime
        // points away from the merge progression (hot end for distance-like
        // axes, cold end for thermal).
        double rootSentinel = descending ? double.NegativeInfinity : double.PositiveInfinity;
        for (int i = 0; i < m; i++)
            if (profile[i] is not null)
                Finalize(i, rootSentinel);

        return new ClusterWalkReport(mass, peakIdx, peakDen, birth, death);
    }

    /// <summary>
    /// Classic excess-of-mass selection over the generalized masses: a
    /// cluster is selected iff its own mass meets or exceeds the best mass
    /// attainable from its subtree; ancestors of a selected cluster suppress
    /// their descendants. Cost- and landscape-blind — it reads masses only.
    /// </summary>
    /// <param name="allowRoot">Roots (never-absorbed internal nodes) are
    /// ineligible by default — selecting a root is the no-structure answer.</param>
    /// <param name="minClusterSize">Eligibility floor: clusters below this
    /// size cannot be selected (their members fall to Unassigned / a
    /// periphery policy), and their mass does not block their ancestors.
    /// NOTE: an eligibility floor, NOT HDBSCAN tree condensation — sub-min
    /// splits still end their parent's lifetime; full producer-side
    /// condensation is the recorded follow-up.</param>
    public static bool[] SelectByExcessOfMass(
        Dendrogram dendrogram, double[] mass, bool allowRoot = false, int minClusterSize = 1)
    {
        ArgumentNullException.ThrowIfNull(dendrogram);
        ArgumentNullException.ThrowIfNull(mass);
        int n = dendrogram.LeafCount;
        int m = dendrogram.InternalNodeCount;
        if (mass.Length != m)
            throw new ArgumentException($"Mass length ({mass.Length}) != internal node count ({m}).", nameof(mass));
        DendrogramNode[] merges = dendrogram.Merges;

        var isChild = new bool[m];
        for (int i = 0; i < m; i++)
        {
            if (merges[i].LeftChild  >= n) isChild[merges[i].LeftChild  - n] = true;
            if (merges[i].RightChild >= n) isChild[merges[i].RightChild - n] = true;
        }

        var subtreeBest = new double[m];
        var selected    = new bool[m];
        for (int i = 0; i < m; i++) // children precede parents in build order
        {
            double childSum = 0.0;
            if (merges[i].LeftChild  >= n) childSum += subtreeBest[merges[i].LeftChild  - n];
            if (merges[i].RightChild >= n) childSum += subtreeBest[merges[i].RightChild - n];

            bool eligible = (allowRoot || isChild[i]) && merges[i].Size >= minClusterSize;
            if (eligible && mass[i] > 0.0 && mass[i] >= childSum)
            {
                selected[i]    = true;
                subtreeBest[i] = mass[i];
            }
            else
            {
                subtreeBest[i] = childSum;
            }
        }

        // Keep only the topmost selected cluster on every root-to-leaf path.
        var suppressed = new bool[m];
        for (int i = m - 1; i >= 0; i--)
        {
            if (suppressed[i]) selected[i] = false;
            bool blockBelow = suppressed[i] || selected[i];
            if (merges[i].LeftChild  >= n) suppressed[merges[i].LeftChild  - n] |= blockBelow;
            if (merges[i].RightChild >= n) suppressed[merges[i].RightChild - n] |= blockBelow;
        }
        return selected;
    }

    /// <summary>
    /// Resolves a selection to the <see cref="Assignment"/> currency: leaves
    /// of each selected cluster share a dense label; everything else is
    /// <see cref="Assignment.Unassigned"/> (the honest abstain — periphery
    /// policies may complete it downstream).
    /// </summary>
    public static Assignment ToAssignment(Dendrogram dendrogram, bool[] selected)
    {
        ArgumentNullException.ThrowIfNull(dendrogram);
        ArgumentNullException.ThrowIfNull(selected);
        int n = dendrogram.LeafCount;
        if (selected.Length != dendrogram.InternalNodeCount)
            throw new ArgumentException(
                $"Selection length ({selected.Length}) != internal node count ({dendrogram.InternalNodeCount}).",
                nameof(selected));

        var labels = new int[n];
        Array.Fill(labels, Assignment.Unassigned);
        int next = 0;
        DendrogramNode[] merges = dendrogram.Merges;
        var stack = new Stack<int>();

        for (int i = 0; i < selected.Length; i++)
        {
            if (!selected[i]) continue;
            int label = next++;
            stack.Push(n + i);
            while (stack.Count > 0)
            {
                int id = stack.Pop();
                if (id < n) { labels[id] = label; continue; }
                DendrogramNode node = merges[id - n];
                stack.Push(node.LeftChild);
                stack.Push(node.RightChild);
            }
        }

        return new Assignment { Labels = labels, Count = next };
    }
}
