#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;
using Maths.LinAlg;

namespace Graphs.Spectral;

/// <summary>
/// Effective resistance distances on a graph Laplacian eigenbasis (DBK2023 gate metric).
/// Conductance weights enter the weighted combinatorial Laplacian; output R_eff values
/// are the filtration metric — never raw couplings or ambient distances.
/// </summary>
public static class EffectiveResistance
{
    public const double DefaultTailEpsilon = 1e-2;
    public const int DefaultKMax = 256;
    private const double LambdaTolerance = 1e-12;
    private const int InitialK = 8;

    /// <summary>
    /// Bottom-K Laplacian eigenpairs with adaptive K: grow until
    /// <c>(1/λ_K)/(1/λ_1) ≤ tailEpsilon</c> on the first non-trivial mode, capped at
    /// <c>min(n−1, kMax)</c>.
    /// </summary>
    public static IReadOnlyList<EigenPair> ComputeEigenpairs(
        CsrGraph graph,
        double tailEpsilon = DefaultTailEpsilon,
        int kMax = DefaultKMax,
        int seed = 0,
        LaplacianType lapType = LaplacianType.Combinatorial,
        SolverKind solverKind = SolverKind.Auto)
    {
        if (tailEpsilon <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(tailEpsilon));
        if (kMax <= 0)
            throw new ArgumentOutOfRangeException(nameof(kMax));

        int n = graph.NodeCount;
        if (n <= 1)
            return Array.Empty<EigenPair>();

        int cap = Math.Min(n, kMax);
        int k = Math.Min(InitialK, cap);
        IReadOnlyList<EigenPair> pairs = Spectral.ComputeBottomK(
            graph, seed, k, lapType, solverKind);

        while (k < cap && !TailConverged(pairs, tailEpsilon, out _))
        {
            int next = Math.Min(Math.Max(k + InitialK, k * 2), cap);
            if (next == k)
                break;
            k = next;
            pairs = Spectral.ComputeBottomK(graph, seed, k, lapType, solverKind);
        }

        return pairs;
    }

    /// <summary>
    /// R_eff(i,j) = Σ_{λ_k &gt; 0} (1/λ_k)(φ_k(i) − φ_k(j))².
    /// Returns <see cref="double.PositiveInfinity"/> when <paramref name="i"/> and
    /// <paramref name="j"/> lie in different connected components.
    /// </summary>
    public static double Pair(
        int i,
        int j,
        IReadOnlyList<EigenPair> eigenpairs,
        ReadOnlySpan<int> componentId)
    {
        if (i < 0 || j < 0 || i >= componentId.Length || j >= componentId.Length)
            throw new ArgumentOutOfRangeException(i == j ? nameof(i) : nameof(j));

        if (componentId[i] != componentId[j])
            return double.PositiveInfinity;

        double sum = 0.0;
        foreach (EigenPair pair in eigenpairs)
        {
            if (pair.Lambda <= LambdaTolerance)
                continue;

            double diff = pair.Vector[i] - pair.Vector[j];
            sum += diff * diff / pair.Lambda;
        }

        return sum;
    }

    /// <summary>
    /// Effective resistance on each undirected CSR edge (canonical u &lt; v).
    /// Cross-component edges are +∞.
    /// </summary>
    public static Dictionary<(int Lo, int Hi), double> BuildEdgeWeights(
        CsrGraph graph,
        IReadOnlyList<EigenPair> eigenpairs)
    {
        ArgumentNullException.ThrowIfNull(eigenpairs);
        int[] componentId = BuildComponentIds(graph);
        var map = new Dictionary<(int Lo, int Hi), double>();

        for (int u = 0; u < graph.NodeCount; u++)
        {
            int rowEnd = graph.RowPointers[u + 1];
            for (int e = graph.RowPointers[u]; e < rowEnd; e++)
            {
                int v = graph.Targets[e];
                if (v <= u) continue;

                double r = Pair(u, v, eigenpairs, componentId);
                if (double.IsNaN(r) || r < 0.0)
                    throw new InvalidOperationException(
                        $"Effective resistance on edge ({u},{v}) is {r}; expected non-negative finite or +∞.");

                map[(u, v)] = r;
            }
        }

        return map;
    }

    internal static bool TailConverged(
        IReadOnlyList<EigenPair> pairs,
        double tailEpsilon,
        out double tailRatio)
    {
        tailRatio = double.PositiveInfinity;
        double lambda1 = 0.0;
        double lambdaK = 0.0;
        int seen = 0;

        foreach (EigenPair pair in pairs)
        {
            if (pair.Lambda <= LambdaTolerance)
                continue;

            if (seen == 0)
                lambda1 = pair.Lambda;
            lambdaK = pair.Lambda;
            seen++;
        }

        if (seen == 0)
            return true;

        tailRatio = (1.0 / lambdaK) / (1.0 / lambda1);
        return tailRatio <= tailEpsilon;
    }

    static int[] BuildComponentIds(CsrGraph graph)
    {
        var uf = new UnionFind(graph.NodeCount);
        for (int u = 0; u < graph.NodeCount; u++)
        {
            int rowEnd = graph.RowPointers[u + 1];
            for (int e = graph.RowPointers[u]; e < rowEnd; e++)
            {
                int v = graph.Targets[e];
                if (u != v)
                    uf.Union(u, v);
            }
        }

        var ids = new int[graph.NodeCount];
        var roots = new Dictionary<int, int>();
        int next = 0;
        for (int i = 0; i < graph.NodeCount; i++)
        {
            int root = uf.Find(i);
            if (!roots.TryGetValue(root, out int id))
            {
                id = next++;
                roots[root] = id;
            }

            ids[i] = id;
        }

        return ids;
    }
}
