using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Deterministic check of the closed-form span crossing (no ensemble): on a hand-built symmetric cubic-spline
/// hump, <see cref="SplineExtrema.SignificantPeakSpans"/> returns one span whose edges sit exactly at the
/// half-drop level — verified by independently evaluating the spline there — and are symmetric about the peak.
/// Pins the Cardano/trig cubic root and the bracket reconstruction the stochastic span test only exercises loosely.
/// </summary>
public sealed class SpanCrossingTests
{
    private readonly ITestOutputHelper _out;
    public SpanCrossingTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void SymmetricHump_SpanEdgesSitAtHalfDropLevel()
    {
        var basis = new SplineBasis(3);
        var config = new KnotConfig(new[] { 0.25, 0.5, 0.75 });   // symmetric interior knots
        double[] coef = { 0.0, 0.2, 0.6, 1.0, 0.6, 0.2, 0.0 };    // symmetric hump (7 = 3 knots + degree + 1)

        var spans = SplineExtrema.SignificantPeakSpans(config, coef, basis, relativeProminence: 0.1, dropFraction: 0.5);

        Assert.Single(spans);
        PeakSpan s = spans[0];

        // The clamped ends sit at 0, so the prominence baseline is 0 and the half-drop level is half the height.
        double level = 0.5 * s.Height;
        double fLeft = basis.Evaluate(config, coef, s.Left);
        double fRight = basis.Evaluate(config, coef, s.Right);
        _out.WriteLine($"peak={s.Location:F5} h={s.Height:F5} level={level:F5} L={s.Left:F5}(f={fLeft:F5}) R={s.Right:F5}(f={fRight:F5})");

        // Closed-form crossings: the spline evaluated at each returned edge equals the half-drop level.
        Assert.Equal(level, fLeft, 6);
        Assert.Equal(level, fRight, 6);
        // Symmetric hump ⇒ peak centred and the span symmetric about it.
        Assert.Equal(0.5, s.Location, 4);
        Assert.Equal(s.Location - s.Left, s.Right - s.Location, 4);
        Assert.False(s.LeftClipped);
        Assert.False(s.RightClipped);
    }

    [Fact]
    public void SymmetricHump_SpanEdgesSitAtHalfDropLevel_Quartic()
    {
        // Degree 4 exercises the general (Durand–Kerner) span reconstruction, not the cubic closed form: a
        // quartic span solved as a cubic would miss the half-drop level by more than the 6-digit tolerance.
        var basis = new SplineBasis(4);
        var config = new KnotConfig(new[] { 0.25, 0.5, 0.75 });
        double[] coef = { 0.0, 0.2, 0.6, 1.0, 1.0, 0.6, 0.2, 0.0 };   // symmetric quartic hump (8 = 3 knots + 4 + 1)

        var spans = SplineExtrema.SignificantPeakSpans(config, coef, basis, relativeProminence: 0.1, dropFraction: 0.5);

        Assert.Single(spans);
        PeakSpan s = spans[0];

        double level = 0.5 * s.Height;   // clamped ends sit at 0 ⇒ prominence baseline 0
        double fLeft = basis.Evaluate(config, coef, s.Left);
        double fRight = basis.Evaluate(config, coef, s.Right);
        _out.WriteLine($"[quartic] peak={s.Location:F5} h={s.Height:F5} level={level:F5} L={s.Left:F5}(f={fLeft:F5}) R={s.Right:F5}(f={fRight:F5})");

        Assert.Equal(level, fLeft, 6);
        Assert.Equal(level, fRight, 6);
        Assert.Equal(0.5, s.Location, 4);
        Assert.Equal(s.Location - s.Left, s.Right - s.Location, 4);
        Assert.False(s.LeftClipped);
        Assert.False(s.RightClipped);
    }
}
