#nullable enable
using System.Collections.Generic;
using Maths.Topology;

namespace TDA.Ph.Tests;

/// <summary>
/// Hand-built filtrations for persistent (co)homology parity tests.
/// </summary>
static class SimplicialFiltrationFixtures
{
    /// <summary>Tetrahedron boundary — triangulated 2-sphere (β₂ = 1).</summary>
    public static SimplicialFiltration TetrahedronBoundary() =>
        new(Build(
            vertices: 4,
            edges: new (int, int)[]
            {
                (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
            },
            triangles: new (int, int, int)[]
            {
                (0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3),
            }), "sphere");

    /// <summary>
    /// n×m grid torus with opposite edges identified — β₁ = 2, β₂ = 1 for n, m ≥ 2.
    /// </summary>
    public static SimplicialFiltration GridTorus(int n, int m) =>
        new(BuildGridTorus(n, m), "torus");

    /// <summary>Wedge of two triangular loops sharing vertex 0 — β₁ = 2.</summary>
    public static SimplicialFiltration TwoLoopWedge() =>
        new(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(0.0, 3),
            new Simplex(0.0, 4),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 1, 2),
            new Simplex(1.0, 0, 2),
            new Simplex(1.0, 0, 3),
            new Simplex(1.0, 3, 4),
            new Simplex(1.0, 0, 4),
        }, "two-loop");

    static IEnumerable<Simplex> Build(int vertices, (int, int)[] edges, (int, int, int)[] triangles)
    {
        for (int v = 0; v < vertices; v++)
            yield return new Simplex(0.0, v);

        foreach (var (a, b) in edges)
            yield return new Simplex(1.0, a, b);

        foreach (var (a, b, c) in triangles)
            yield return new Simplex(2.0, a, b, c);
    }

    static IEnumerable<Simplex> BuildGridTorus(int n, int m)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                yield return new Simplex(0.0, Vid(i, j, n, m));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                yield return new Simplex(1.0, Vid(i, j, n, m), Vid(i + 1, j, n, m));
                yield return new Simplex(1.0, Vid(i, j, n, m), Vid(i, j + 1, n, m));
            }

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                int a = Vid(i, j, n, m);
                int b = Vid(i + 1, j, n, m);
                int c = Vid(i + 1, j + 1, n, m);
                int d = Vid(i, j + 1, n, m);
                yield return new Simplex(1.0, a, c);
                yield return new Simplex(2.0, a, b, c);
                yield return new Simplex(2.0, a, c, d);
            }
    }

    static int Vid(int i, int j, int n, int m) => (i % n) * m + (j % m);
}
