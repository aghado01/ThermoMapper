#nullable enable
using System;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// The two screening-scale diagram distances (sliced, Sinkhorn) validated against the in-repo
/// exact Hungarian oracle (<see cref="DiagramMetrics.Wasserstein"/>) and analytic anchors.
/// The residual T4transport numerical cross-check is tracked in the ISOLET brief's
/// "H0 matching-cost gate".
/// </summary>
public sealed class DiagramMetricsApproxTests
{
    // Small deterministic mixed-scale fixtures — no RNG, per repo hygiene.
    private static Barcode FixtureA() => new(new[]
    {
        new Bar(0.0, 2.0, 1),
        new Bar(1.0, 3.5, 1),
        new Bar(0.5, 0.9, 1),
        new Bar(2.0, 6.0, 1),
        new Bar(4.0, 4.4, 1),
    });

    private static Barcode FixtureB() => new(new[]
    {
        new Bar(0.2, 2.3, 1),
        new Bar(1.1, 3.0, 1),
        new Bar(2.5, 5.0, 1),
        new Bar(0.7, 1.0, 1),
    });

    // ── Sliced ────────────────────────────────────────────────────────────────

    [Fact]
    public void Sliced_IdenticalDiagrams_ReturnZero()
    {
        double sw = DiagramMetrics.SlicedWasserstein(FixtureA(), FixtureA(), dimension: 1);
        Assert.Equal(0.0, sw, precision: 12);
    }

    [Fact]
    public void Sliced_SinglePointVsEmpty_W1MatchesAnalyticIntegral()
    {
        // One bar (0, 2) vs empty: side A = {(0,2)}, side B = {(1,1)} (its diagonal projection).
        // Per slice θ: |2·sinθ − (cosθ + sinθ)| = √2·|sin(θ − π/4)|, and the mean of |sin| over a
        // half-period is 2/π — so SW_1 = 2√2/π. Midpoint slices integrate this to high accuracy.
        var a = new Barcode(new[] { new Bar(0.0, 2.0, 0) });
        var b = new Barcode(Array.Empty<Bar>());

        double sw = DiagramMetrics.SlicedWasserstein(a, b, dimension: 0, p: 1.0, directions: 720);
        Assert.Equal(2.0 * Math.Sqrt(2.0) / Math.PI, sw, precision: 4);
    }

    [Fact]
    public void Sliced_Symmetric()
    {
        double ab = DiagramMetrics.SlicedWasserstein(FixtureA(), FixtureB(), dimension: 1);
        double ba = DiagramMetrics.SlicedWasserstein(FixtureB(), FixtureA(), dimension: 1);
        Assert.Equal(ab, ba, precision: 12);
    }

    [Theory]
    [InlineData(1.0)]
    [InlineData(2.0)]
    public void Sliced_ScaleEquivariant(double p)
    {
        const double scale = 3.7;
        double sw = DiagramMetrics.SlicedWasserstein(FixtureA(), FixtureB(), dimension: 1, p);
        double swScaled = DiagramMetrics.SlicedWasserstein(
            Scale(FixtureA(), scale), Scale(FixtureB(), scale), dimension: 1, p);
        Assert.Equal(scale * sw, swScaled, precision: 9);
    }

    [Fact]
    public void Sliced_EssentialMismatch_InfiniteOnMismatch_ReturnsInfinity()
    {
        var a = new Barcode(new[] { new Bar(0.0, double.PositiveInfinity, 0) });
        var b = new Barcode(Array.Empty<Bar>());
        Assert.True(double.IsPositiveInfinity(
            DiagramMetrics.SlicedWasserstein(a, b, dimension: 0)));
    }

    [Fact]
    public void Sliced_EssentialMismatch_FinitePenalty_ChargesPerBar()
    {
        var a = new Barcode(new[] { new Bar(0.0, double.PositiveInfinity, 0) });
        var b = new Barcode(Array.Empty<Bar>());
        var policy = DiagramMetrics.EssentialPolicy.FinitePenalty(perBar: 2.0);

        double sw = DiagramMetrics.SlicedWasserstein(a, b, dimension: 0, p: 2.0, essential: policy);
        Assert.Equal(2.0, sw, precision: 12);
    }

    [Fact]
    public void Sliced_EssentialSameCountDifferentBirth_ChargesBirthDistance()
    {
        // The essential term is slice-independent and identical to the exact backend's:
        // |Δbirth| = 99 for (1, ∞) vs (100, ∞), and 0 once the births coincide.
        var a = new Barcode(new[] { new Bar(1.0, double.PositiveInfinity, 1) });
        var b = new Barcode(new[] { new Bar(100.0, double.PositiveInfinity, 1) });

        double sw = DiagramMetrics.SlicedWasserstein(a, b, dimension: 1, p: 2.0);
        Assert.Equal(99.0, sw, precision: 12);

        double self = DiagramMetrics.SlicedWasserstein(a, a, dimension: 1, p: 2.0);
        Assert.Equal(0.0, self, precision: 12);
    }

    [Fact]
    public void Sliced_RejectsNonPositiveDirections()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            DiagramMetrics.SlicedWasserstein(FixtureA(), FixtureB(), dimension: 1, directions: 0));
    }

    [Fact]
    public void Sliced_RejectsPBelowOne()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            DiagramMetrics.SlicedWasserstein(FixtureA(), FixtureB(), dimension: 1, p: 0.5));
    }

    // ── Sinkhorn ──────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(1.0)]
    [InlineData(2.0)]
    public void Sinkhorn_SmallEpsilon_ConvergesToExactHungarian(double p)
    {
        double exact = DiagramMetrics.Wasserstein(FixtureA(), FixtureB(), dimension: 1, p);
        double sink = DiagramMetrics.SinkhornWasserstein(
            FixtureA(), FixtureB(), dimension: 1, p, epsilon: 1e-3, maxIters: 5000);

        Assert.True(double.IsFinite(sink));
        Assert.InRange(Math.Abs(sink - exact) / exact, 0.0, 1e-2);
    }

    [Fact]
    public void Sinkhorn_UpperBoundsExact()
    {
        // The entropic plan is feasible for the assignment LP, so its cost cannot undercut the
        // exact optimum (modulo the residual marginal violation at the iteration cap).
        double exact = DiagramMetrics.Wasserstein(FixtureA(), FixtureB(), dimension: 1, p: 2.0);
        double sink = DiagramMetrics.SinkhornWasserstein(
            FixtureA(), FixtureB(), dimension: 1, p: 2.0, epsilon: 1e-3, maxIters: 5000);

        Assert.True(sink >= exact - 1e-6, $"Sinkhorn {sink} undercuts exact {exact}.");
    }

    [Fact]
    public void Sinkhorn_IdenticalDiagrams_SelfDistanceShrinksWithEpsilon()
    {
        // The entropic bias: self-distance is positive (mass smears onto the near-diagonal escape
        // cells whose cost is comparable to ε) but vanishes as ε → 0.
        double coarse = DiagramMetrics.SinkhornWasserstein(
            FixtureA(), FixtureA(), dimension: 1, p: 2.0, epsilon: 1e-3, maxIters: 5000);
        double fine = DiagramMetrics.SinkhornWasserstein(
            FixtureA(), FixtureA(), dimension: 1, p: 2.0, epsilon: 1e-4, maxIters: 20000);

        Assert.True(fine < coarse, $"Self-distance grew as ε shrank: {fine} !< {coarse}.");
        Assert.InRange(fine, 0.0, 1e-2);
    }

    [Fact]
    public void Sinkhorn_SinglePointVsEmpty_MatchesDiagonalDistance()
    {
        // Exact value 1 (see the exact-metric test); one admissible cell, so Sinkhorn must agree.
        var a = new Barcode(new[] { new Bar(0.0, 2.0, 0) });
        var b = new Barcode(Array.Empty<Bar>());

        double sink = DiagramMetrics.SinkhornWasserstein(
            a, b, dimension: 0, p: 2.0, epsilon: 1e-3, maxIters: 2000);
        Assert.Equal(1.0, sink, precision: 3);
    }

    [Fact]
    public void Sinkhorn_EssentialMismatch_FinitePenalty_ChargesPerBar()
    {
        var a = new Barcode(new[] { new Bar(0.0, double.PositiveInfinity, 0) });
        var b = new Barcode(Array.Empty<Bar>());
        var policy = DiagramMetrics.EssentialPolicy.FinitePenalty(perBar: 2.0);

        double sink = DiagramMetrics.SinkhornWasserstein(a, b, dimension: 0, p: 2.0, essential: policy);
        Assert.Equal(2.0, sink, precision: 12);
    }

    [Fact]
    public void Sinkhorn_EssentialSameCountDifferentBirth_ChargesBirthDistance()
    {
        // The essential term is the exact birth matching, never entropically smoothed —
        // essential-only diagrams give exactly |Δbirth| = 99 regardless of ε.
        var a = new Barcode(new[] { new Bar(1.0, double.PositiveInfinity, 1) });
        var b = new Barcode(new[] { new Bar(100.0, double.PositiveInfinity, 1) });

        double sink = DiagramMetrics.SinkhornWasserstein(a, b, dimension: 1, p: 2.0);
        Assert.Equal(99.0, sink, precision: 12);
    }

    [Theory]
    [InlineData(1.0)]
    [InlineData(10.0)]
    public void Sinkhorn_ForcedAdmissibleSupport_LargeEpsilon_MatchesExact(double epsilon)
    {
        // Two bars vs empty: each bar's only admissible column is its own diagonal cell — the
        // cross diagonal cells carry the forbidden sentinel, so the constrained plan is forced to
        // the identity and the entropic value must equal the exact distance at ANY ε. Before the
        // log-kernel −∞ mask, exp(−big′/ε) was no longer negligible at these ε and the sentinel
        // cells won real mass, drifting the value toward the sentinel scale (≈ 138 and ≈ 323 here
        // versus the exact 101).
        var a = new Barcode(new[] { new Bar(0.0, 2.0, 0), new Bar(0.0, 200.0, 0) });
        var b = new Barcode(Array.Empty<Bar>());

        double exact = DiagramMetrics.Wasserstein(a, b, dimension: 0, p: 1.0);
        double sink = DiagramMetrics.SinkhornWasserstein(a, b, dimension: 0, p: 1.0, epsilon: epsilon);

        Assert.Equal(101.0, exact, precision: 12);
        Assert.Equal(exact, sink, precision: 9);
    }

    [Fact]
    public void Sinkhorn_LargeEpsilon_RemainsFeasibleOverAdmissibleSupport()
    {
        // At large ε the value smears across admissible cells, but the plan must stay a feasible
        // point of the admissible-support assignment LP: finite, and never below the exact optimum.
        double exact = DiagramMetrics.Wasserstein(FixtureA(), FixtureB(), dimension: 1, p: 2.0);
        double sink = DiagramMetrics.SinkhornWasserstein(
            FixtureA(), FixtureB(), dimension: 1, p: 2.0, epsilon: 1.0, maxIters: 2000);

        Assert.True(double.IsFinite(sink));
        Assert.True(sink >= exact - 1e-6, $"Large-ε Sinkhorn {sink} undercuts exact {exact}.");
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-1.0)]
    [InlineData(double.NaN)]
    public void Sinkhorn_RejectsBadEpsilon(double epsilon)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            DiagramMetrics.SinkhornWasserstein(FixtureA(), FixtureB(), dimension: 1, epsilon: epsilon));
    }

    [Fact]
    public void Sinkhorn_RejectsNonPositiveIterations()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            DiagramMetrics.SinkhornWasserstein(FixtureA(), FixtureB(), dimension: 1, maxIters: 0));
    }

    [Fact]
    public void BothApproximations_EmptyVsEmpty_ReturnZero()
    {
        var empty = new Barcode(Array.Empty<Bar>());
        Assert.Equal(0.0, DiagramMetrics.SlicedWasserstein(empty, empty, dimension: 0));
        Assert.Equal(0.0, DiagramMetrics.SinkhornWasserstein(empty, empty, dimension: 0));
    }

    [Fact]
    public void AllBackends_EssentialOnly_AgreeOnMatchedBirthTerm()
    {
        // The matched-essential term is the same exact birth assignment in all three backends, so
        // essential-only diagrams pin the cross-backend semantics: identical values, no tolerance.
        var a = new Barcode(new[] { new Bar(1.0, double.PositiveInfinity, 1) });
        var b = new Barcode(new[] { new Bar(100.0, double.PositiveInfinity, 1) });

        double w = DiagramMetrics.Wasserstein(a, b, dimension: 1, p: 2.0);
        double sw = DiagramMetrics.SlicedWasserstein(a, b, dimension: 1, p: 2.0);
        double sink = DiagramMetrics.SinkhornWasserstein(a, b, dimension: 1, p: 2.0);

        Assert.Equal(99.0, w, precision: 12);
        Assert.Equal(w, sw, precision: 12);
        Assert.Equal(w, sink, precision: 12);
    }

    private static Barcode Scale(Barcode barcode, double factor)
    {
        var bars = new Bar[barcode.Bars.Count];
        for (int i = 0; i < bars.Length; i++)
        {
            Bar bar = barcode.Bars[i];
            bars[i] = bar with { Birth = bar.Birth * factor, Death = bar.Death * factor };
        }
        return new Barcode(bars);
    }
}
