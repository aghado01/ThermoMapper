using System;
using System.Collections.Generic;
using Graphs.Primitives;
using TDA.Ph;
using Viz;
using Xunit;

namespace VizCore.Tests;

public sealed class TriangleLayerTests
{
    [Fact]
    public void FlagComplex_TetrahedronSkeleton_ProducesFourTriangles()
    {
        EdgeLayer edges = BuildEdgeLayer(
            id: "tetrahedron-edges",
            name: "Tetrahedron Edges",
            src: new[] { 0, 0, 0, 1, 1, 2 },
            dst: new[] { 1, 2, 3, 2, 3, 3 });

        int[] triangles = FlagComplex.FromEdges(edges.Src.Span, edges.Dst.Span);

        Assert.Equal(new[]
        {
            0, 1, 2,
            0, 1, 3,
            0, 2, 3,
            1, 2, 3,
        }, triangles);
    }

    [Fact]
    public void FlagComplex_CountTriangles_TetrahedronSkeleton_ReturnsFour()
    {
        EdgeLayer edges = BuildEdgeLayer(
            id: "tetrahedron-edges",
            name: "Tetrahedron Edges",
            src: new[] { 0, 0, 0, 1, 1, 2 },
            dst: new[] { 1, 2, 3, 2, 3, 3 });

        int triangleCount = FlagComplex.CountTriangles(BuildGraph(edges.Src.ToArray(), edges.Dst.ToArray()));

        Assert.Equal(4, triangleCount);
    }

    [Fact]
    public void FlagComplex_TriangleClique_ProducesOneTriangle()
    {
        EdgeLayer edges = BuildEdgeLayer(
            id: "triangle-edges",
            name: "Triangle Edges",
            src: new[] { 0, 0, 1 },
            dst: new[] { 1, 2, 2 });

        int[] triangles = FlagComplex.FromEdges(edges.Src.Span, edges.Dst.Span);

        Assert.Equal(new[] { 0, 1, 2 }, triangles);
    }

    [Fact]
    public void FlagComplex_CountTriangles_TriangleClique_ReturnsOne()
    {
        EdgeLayer edges = BuildEdgeLayer(
            id: "triangle-edges",
            name: "Triangle Edges",
            src: new[] { 0, 0, 1 },
            dst: new[] { 1, 2, 2 });

        int triangleCount = FlagComplex.CountTriangles(BuildGraph(edges.Src.ToArray(), edges.Dst.ToArray()));

        Assert.Equal(1, triangleCount);
    }

    [Fact]
    public void FlagComplex_SquareCycle_ProducesNoTriangles()
    {
        EdgeLayer edges = BuildEdgeLayer(
            id: "square-cycle-edges",
            name: "Square Cycle",
            src: new[] { 0, 1, 2, 3 },
            dst: new[] { 1, 2, 3, 0 });

        int[] triangles = FlagComplex.FromEdges(edges.Src.Span, edges.Dst.Span);

        Assert.Empty(triangles);
    }

    [Fact]
    public void FlagComplex_CountTriangles_SquareCycle_ReturnsZero()
    {
        EdgeLayer edges = BuildEdgeLayer(
            id: "square-cycle-edges",
            name: "Square Cycle",
            src: new[] { 0, 1, 2, 3 },
            dst: new[] { 1, 2, 3, 0 });

        int triangleCount = FlagComplex.CountTriangles(BuildGraph(edges.Src.ToArray(), edges.Dst.ToArray()));

        Assert.Equal(0, triangleCount);
    }

    [Fact]
    public void FlagComplex_CountTriangles_MatchesFromEdges_OnDeterministicRandomGraph()
    {
        const int nodeCount = 7;
        var random = new Random(12345);
        var edges = new List<(int Source, int Target)>();

        for (int source = 0; source < nodeCount; source++)
        {
            for (int target = source + 1; target < nodeCount; target++)
            {
                if (random.NextDouble() < 0.45)
                    edges.Add((source, target));
            }
        }

        int[] src = new int[edges.Count];
        int[] dst = new int[edges.Count];
        for (int i = 0; i < edges.Count; i++)
        {
            src[i] = edges[i].Source;
            dst[i] = edges[i].Target;
        }

        int[] triangles = FlagComplex.FromEdges(src, dst);
        int triangleCount = FlagComplex.CountTriangles(BuildGraph(src, dst));

        Assert.Equal(triangles.Length / 3, triangleCount);
    }

    [Fact]
    public void TriangleLayer_FromFlagComplex_ThreadsSourceEdgeLayerId()
    {
        EdgeLayer edges = BuildEdgeLayer(
            id: "edges-primary",
            name: "Primary Edge Layer",
            src: new[] { 0, 0, 1 },
            dst: new[] { 1, 2, 2 });

        int[] triangles = FlagComplex.FromEdges(edges.Src.Span, edges.Dst.Span);
        TriangleLayer layer = TriangleLayer.FromFlagComplex(
            edges,
            triangles,
            name: "Primary Triangles",
            id: "triangles-primary");

        Assert.Equal("triangles-primary", layer.Id);
        Assert.Equal("Primary Triangles", layer.Name);
        Assert.Equal(TriangleSource.FlagComplex, layer.Source);
        Assert.Equal(edges.Id, layer.SourceEdgeLayerId);
        Assert.Equal(triangles, layer.Vertices.ToArray());
    }

    private static EdgeLayer BuildEdgeLayer(string id, string name, int[] src, int[] dst)
    {
        if (src.Length != dst.Length)
            throw new ArgumentException("Edge arrays must have matching lengths.");

        var weights = new double[src.Length];
        Array.Fill(weights, 1.0);
        return new EdgeLayer(name, src, dst, weights, id: id);
    }

    private static CsrGraph BuildGraph(int[] src, int[] dst)
    {
        var edges = new Edge[src.Length];
        for (int i = 0; i < src.Length; i++)
            edges[i] = new Edge(src[i], dst[i], 1.0);

        int nodeCount = 0;
        for (int i = 0; i < src.Length; i++)
            nodeCount = Math.Max(nodeCount, Math.Max(src[i], dst[i]) + 1);

        return CsrGraph.FromEdges(edges, nodeCount);
    }
}
