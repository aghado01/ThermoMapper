#nullable enable
using System.Collections.Generic;
using System.Linq;
using Graphs.Primitives;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// P0 · <see cref="ConditionedFiltration"/> — the undirected content core. Proves a
/// similarity chord closing a path <em>relative to a backbone</em> is an H₁ generator, and
/// that SIFTS falls out as the <c>τ≡0</c> degenerate. The barcode routes through the existing
/// reducer unchanged — these tests exercise the union convention, not new engine code.
/// </summary>
public sealed class ConditionedFiltrationTests
{
    // ── Test 1 · a return relative to the backbone ──────────────────────────────

    [Fact]
    public void PathPlusChord_IsEssentialH1_BornAtChordDistance()
    {
        // Backbone: path v0–v1–…–v5 at ε₀ = 0.  Similarity: one chord (v0, v5) at distance d.
        const int n = 6;
        const double d = 2.5;
        List<(int, int)> backbone = PathEdges(n);
        var similarity = new List<(int, int, double)> { (0, 5, d) };

        Barcode barcode = ConditionedFiltration.ComputeBarcode(n, backbone, similarity);

        Bar loop = Assert.Single(SignificantH1(barcode));
        Assert.True(loop.IsInfinite, "the loop never fills (no triangle) → essential H₁");
        Assert.Equal(d, loop.Birth, 12); // born exactly when the chord enters
    }

    [Fact]
    public void PathPlusChord_ReconstructedCycle_IsPathPlusChord()
    {
        const int n = 6;
        const double d = 2.5;
        List<(int, int)> backbone = PathEdges(n);
        var similarity = new List<(int, int, double)> { (0, 5, d) };

        CsrGraph graph = ConditionedFiltration.BuildGraph(n, backbone, similarity);
        var filtration = RipsFiltration.RipsFromGraph(graph, FiltrationWeights.RawDistance, 2, "conditioned");
        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: true);

        (Bar Bar, IReadOnlyList<UndirectedEdge> Edges) loop = Assert.Single(BarCycleEdges.H1Loops(barcode, filtration));
        HashSet<(int, int)> reconstructed = loop.Edges.Select(e => (e.Lo, e.Hi)).ToHashSet();

        var expected = new HashSet<(int, int)>
        {
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), // backbone path
            (0, 5),                                 // the closing chord
        };
        Assert.Equal(expected, reconstructed);
    }

    [Fact]
    public void PathWithoutChord_HasNoH1()
    {
        const int n = 6;
        List<(int, int)> backbone = PathEdges(n);
        var similarity = new List<(int, int, double)>(); // no return

        Barcode barcode = ConditionedFiltration.ComputeBarcode(n, backbone, similarity);

        Assert.Empty(SignificantH1(barcode));
    }

    // ── Test 2 · SIFTS is the τ≡0 degenerate (reading-order backbone) ────────────
    // Backbone = reading order (the τ≡0 prior); similarity = token tie-backs. A long-range
    // tie-back ("spider" recurring across the rhyme) is an essential loop; a purely linear
    // rhyme has none. Synthetic analog per P0 brief §5.2.

    [Fact]
    public void Sifts_LinearRhyme_HasNoH1() // Row-Row-Row → β₁ = 0
    {
        const int n = 6;
        List<(int, int)> readingOrder = PathEdges(n);
        var similarity = new List<(int, int, double)>(); // no tie-back

        Barcode barcode = ConditionedFiltration.ComputeBarcode(n, readingOrder, similarity);

        Assert.Empty(SignificantH1(barcode));
    }

    [Fact]
    public void Sifts_TieBackRhyme_HasOneH1() // Itsy-Bitsy → β₁ = 1
    {
        // Token 1 recurs at token 5 → a long tie-back chord spanning the rhyme.
        const int n = 6;
        List<(int, int)> readingOrder = PathEdges(n);
        var similarity = new List<(int, int, double)> { (1, 5, 1.0) };

        Barcode barcode = ConditionedFiltration.ComputeBarcode(n, readingOrder, similarity);

        Assert.Single(SignificantH1(barcode));
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    static List<(int, int)> PathEdges(int n)
    {
        var edges = new List<(int, int)>(n - 1);
        for (int i = 0; i + 1 < n; i++)
            edges.Add((i, i + 1));
        return edges;
    }

    static IEnumerable<Bar> SignificantH1(Barcode barcode) =>
        barcode.Bars.Where(b => b.Dimension == 1 && (b.IsInfinite || b.Death > b.Birth));
}
