#nullable enable
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// P1a · <see cref="ResidualPrior"/> — the prior generalizes. Proves the residual band is a thin
/// weight producer feeding P0's <see cref="ConditionedFiltration"/>: <c>τ≡0</c> subsumes P0's raw-distance
/// similarity, and a nonzero prior shifts a return's <em>birth</em> (the payoff P0 cannot express).
/// Undirected, monotone δ — no Δ, directedness, or zigzag.
/// </summary>
public sealed class ResidualPriorTests
{
    // A linear field so gaps read off directly: |t_j − t_i| = 0.5·(j − i).
    static readonly double[] Field = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 };

    [Fact]
    public void TauZero_ResidualEqualsRawGap_AndFoldIsInert()
    {
        var prior = new List<(int, int, double)> { (0, 5, 0.0) };

        var min = ResidualPrior.ResidualEdges(Field, prior, ResidualSymmetry.Min);
        var max = ResidualPrior.ResidualEdges(Field, prior, ResidualSymmetry.Max);

        (int i, int j, double r) e = Assert.Single(min);
        Assert.Equal((0, 5), (e.i, e.j));
        Assert.Equal(2.5, e.r, 12);            // raw gap |t5 − t0|
        Assert.Equal(max[0].r, min[0].r, 12);  // τ=0 ⇒ forward == reverse ⇒ fold inert
    }

    [Fact]
    public void TauZero_SubsumesP0Similarity()
    {
        const int n = 6;
        List<(int, int)> backbone = PathEdges(n);
        var viaPrior = ResidualPrior.ResidualEdges(Field, new List<(int, int, double)> { (0, 5, 0.0) });

        Bar throughPrior = Assert.Single(SignificantH1(
            ConditionedFiltration.ComputeBarcode(n, backbone, viaPrior)));
        Bar throughP0 = Assert.Single(SignificantH1(
            ConditionedFiltration.ComputeBarcode(n, backbone, new List<(int, int, double)> { (0, 5, 2.5) })));

        Assert.True(throughPrior.IsInfinite);
        Assert.Equal(2.5, throughPrior.Birth, 12);
        Assert.Equal(throughP0.Birth, throughPrior.Birth, 12); // identical to P0's similarity path
    }

    [Fact]
    public void Prior_ShiftsReturnBirth() // the P1 payoff P0 cannot express
    {
        const int n = 6;
        List<(int, int)> backbone = PathEdges(n);

        var unpredicted = ResidualPrior.ResidualEdges(Field, new List<(int, int, double)> { (0, 5, 0.0) }); // τ=0
        var predicted = ResidualPrior.ResidualEdges(Field, new List<(int, int, double)> { (0, 5, 2.5) });   // τ = observed gap

        Bar late = Assert.Single(SignificantH1(ConditionedFiltration.ComputeBarcode(n, backbone, unpredicted)));
        Bar early = Assert.Single(SignificantH1(ConditionedFiltration.ComputeBarcode(n, backbone, predicted)));

        Assert.Equal(2.5, late.Birth, 12);  // unpredicted return: born at its raw residual
        Assert.Equal(0.0, early.Birth, 12); // predicted return: born at 0 — the prior expected it
        Assert.True(late.Birth > early.Birth);
    }

    [Fact]
    public void ResidualOrder_IsTheFiltration()
    {
        const int n = 6;
        List<(int, int)> backbone = PathEdges(n);
        // Two unpredicted chords at different raw gaps: (0,3) → 1.5, (0,5) → 2.5.
        var content = ResidualPrior.ResidualEdges(
            Field, new List<(int, int, double)> { (0, 3, 0.0), (0, 5, 0.0) });

        double[] births = SignificantH1(ConditionedFiltration.ComputeBarcode(n, backbone, content))
            .Select(b => b.Birth).OrderBy(x => x).ToArray();

        Assert.Equal(2, births.Length);
        Assert.Equal(1.5, births[0], 12);
        Assert.Equal(2.5, births[1], 12); // born in residual order
    }

    [Fact]
    public void Symmetry_MinPicksBetterOrientedPrediction()
    {
        var prior = new List<(int, int, double)> { (0, 5, 1.0) }; // τ ≠ 0 ⇒ orientations diverge

        double min = ResidualPrior.ResidualEdges(Field, prior, ResidualSymmetry.Min)[0].r;
        double max = ResidualPrior.ResidualEdges(Field, prior, ResidualSymmetry.Max)[0].r;
        double mean = ResidualPrior.ResidualEdges(Field, prior, ResidualSymmetry.Mean)[0].r;

        Assert.Equal(1.5, min, 12);  // |t5 − (t0 + 1)|
        Assert.Equal(3.5, max, 12);  // |t0 − (t5 + 1)|
        Assert.Equal(2.5, mean, 12);
        Assert.True(min < max);
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
