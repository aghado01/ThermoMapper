#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.LinAlg;
using Maths.Topology;

namespace TDA.Ph;

/// <summary>
/// Selects edge filtration values for <see cref="RipsFiltration.GraphRips"/>.
/// Filtration values are never coupling affinities — use <see cref="RawDistanceWeights"/>
/// for CSR distance weights, or <see cref="EffectiveResistanceWeights"/> (gate metric, P0#3).
/// </summary>
public abstract record FiltrationWeights
{
    public static FiltrationWeights RawDistance { get; } = RawDistanceWeights.Instance;
}

/// <summary>CSR edge weights as Rips filtration values (SPRED SA path).</summary>
public sealed record RawDistanceWeights : FiltrationWeights
{
    internal static readonly RawDistanceWeights Instance = new();
    private RawDistanceWeights() { }
}

/// <summary>
/// Effective-resistance edge weights from Laplacian eigenpairs (faithfulness gate).
/// </summary>
public sealed record EffectiveResistanceWeights(IReadOnlyList<EigenPair> Eigenpairs) : FiltrationWeights
{
    /// <summary>Compute eigenpairs (adaptive K) and wrap for <see cref="RipsFiltration"/>.</summary>
    public static EffectiveResistanceWeights FromGraph(
        CsrGraph graph,
        double tailEpsilon = EffectiveResistance.DefaultTailEpsilon,
        int kMax = EffectiveResistance.DefaultKMax,
        int seed = 0,
        SolverKind solverKind = SolverKind.Auto) =>
        new(EffectiveResistance.ComputeEigenpairs(
            graph, tailEpsilon, kMax, seed, solverKind: solverKind));
}

/// <summary>
/// Graph-restricted (lazy) Vietoris–Rips filtration builder over a symmetric CSR skeleton.
/// </summary>
public static class RipsFiltration
{
    /// <summary>
    /// Materialize a <see cref="SimplicialFiltration"/> from a weighted kNN / conditioned graph.
    /// Triangle enumeration delegates to <see cref="FlagComplex.Triangles"/>; edge weights are
    /// resolved from a canonical <c>(lo, hi)</c> map built in O(|E|).
    /// </summary>
    /// <param name="g">Symmetric CSR graph. Use distance (not coupling) edge weights with
    /// <see cref="RawDistanceWeights"/>.</param>
    /// <param name="weights">Filtration-value source (raw distance or effective resistance).</param>
    /// <param name="maxDimension">2 for H0+H1 (loops need triangle fillers); 1 for H0-only.</param>
    /// <param name="label">Filtration label stored on the output complex.</param>
    public static SimplicialFiltration GraphRips(
        CsrGraph g,
        FiltrationWeights weights,
        int maxDimension = 2,
        string label = "Rips")
    {
        ArgumentNullException.ThrowIfNull(weights);
        if (maxDimension < 0)
            throw new ArgumentOutOfRangeException(nameof(maxDimension));

        Dictionary<(int Lo, int Hi), double> edgeWeights = ResolveEdgeWeights(g, weights);
        int n = g.NodeCount;
        int triCount = maxDimension >= 2 ? FlagComplex.Triangles(g).Length / 3 : 0;
        var simplices = new List<Simplex>(n + edgeWeights.Count + triCount);

        for (int i = 0; i < n; i++)
            simplices.Add(new Simplex(0.0, i));

        foreach (var ((u, v), w) in edgeWeights)
            simplices.Add(new Simplex(w, u, v));

        if (maxDimension >= 2)
        {
            int[] tris = FlagComplex.Triangles(g);
            for (int t = 0; t < tris.Length; t += 3)
            {
                int u = tris[t];
                int v = tris[t + 1];
                int w = tris[t + 2];
                double val = Math.Max(
                    edgeWeights[(u, v)],
                    Math.Max(edgeWeights[(u, w)], edgeWeights[(v, w)]));
                simplices.Add(new Simplex(val, u, v, w));
            }
        }

        return new SimplicialFiltration(simplices, label);
    }

    static Dictionary<(int Lo, int Hi), double> ResolveEdgeWeights(CsrGraph g, FiltrationWeights weights) =>
        weights switch
        {
            RawDistanceWeights => BuildRawDistanceWeights(g),
            EffectiveResistanceWeights er => EffectiveResistance.BuildEdgeWeights(g, er.Eigenpairs),
            _ => throw new ArgumentException($"Unknown {nameof(FiltrationWeights)}: {weights.GetType().Name}", nameof(weights)),
        };

    static Dictionary<(int Lo, int Hi), double> BuildRawDistanceWeights(CsrGraph g)
    {
        var map = new Dictionary<(int Lo, int Hi), double>();
        for (int u = 0; u < g.NodeCount; u++)
        {
            int rowEnd = g.RowPointers[u + 1];
            for (int e = g.RowPointers[u]; e < rowEnd; e++)
            {
                int v = g.Targets[e];
                if (v <= u) continue;

                double w = g.Weights[e];
                if (w < 0.0)
                    throw new ArgumentException("Rips filtration requires non-negative edge weights (distances).");

                map[(u, v)] = w;
            }
        }

        return map;
    }
}
