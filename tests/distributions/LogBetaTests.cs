using System;
using Maths.Distributions;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Distributions.Tests;

/// <summary>
/// <see cref="SpecialFunctions.LogBeta"/> against ln B to 400 digits from outside the engine. Until 2026-11-15 it was ln Γ(a) + ln Γ(b) −
/// ln Γ(a + b), which cancels where the shapes are far apart (ln B(10¹⁶, 1) was 0, its value −36.84; 3.6·10⁻⁶ off at (10¹⁰, ½)) and loses
/// up to three digits where both are large. The qualification behind the change: 2428 points from 10⁻¹² to 10³⁰⁰ against Python's
/// decimal at 400 digits, cross-checked against mpmath 1.3.0 (the defects brief's slice 9).
/// </summary>
public sealed class LogBetaTests(ITestOutputHelper output)
{
    // (a, b, ln B) with ln B from mpmath 1.3.0 at 400 digits, rounded to the double: balanced and unequal shapes, small to 10³⁰⁰, and
    // the points where the three-term form or the ratio form alone was worst in the qualification.
    public static readonly TheoryData<double, double, double> Reference = new()
    {
        { 100000.0, 1.0, -11.512925464970229 },
        { 1e+16, 1.0, -36.841361487904734 },
        { 1e+300, 1.0, -690.7755278982137 },
        { 0.5, 0.5, 1.1447298858494002 },
        { 1.0, 1.0, 0.0 },
        { 3.0, 4.0, -4.0943445622221 },
        { 30.0, 40.0, -48.30174909591612 },
        { 1e-300, 1e-300, 691.4686750787737 },
        { 1e-08, 1.0, 18.420680743952367 },
        { 100000.0, 100000.0, -138633.92706134808 },
        { 1000000000000000.0, 1000000000000000.0, -1386294361119906.5 },
        { 1e+300, 1e+300, -1.3862943611198907e+300 },
        { 100000.0, 2.5, -28.497649541827652 },
        { 10000000000.0, 0.5, -10.940560522033028 },
        { 1e-10, 10000000000.0, 23.02585092758015 },
        { 16.0, 16.0, -22.293680785633526 },
        { 15.9, 0.1, 1.978934483062005 },
        { 100.0, 1e-12, 27.63102111592337 },
        { 250000.0, 250001.0, -346579.23252262815 },
        { 100000000.0, 3.0, -54.56889508129715 },
        { 7.0, 1000000000.0, -138.4836096676148 },
        { 0.3, 0.7, 1.3566652413497422 },
        { 0.001, 0.001, 7.600900817008347 },
        { 9.999, 10.001, -13.736229121870219 },
        { 10.0, 10.0, -13.736229227036555 },
        { 10.0, 1000000.0, -125.35332309941877 },
        { 1.1619262291079851e+199, 1.128459826262219e+199, -1.5873301276805681e+199 },
        { 3000000000000000.0, 2980000000000000.0, -4144986694870086.0 },
        { 41.76594050915185, 40.51039258481092, -57.60982819982397 },
        { 452473311.99775445, 507.53877243231705, -7463.3540619963305 },
        // a + b overflows the double in these two; ln B does not.
        { 1e+300, 1.7976931348623157e+308, -2.0007184997951635e+301 },
        { 9e+307, 9.5e+307, -1.2816465260771486e+308 },
    };

    [Theory]
    [MemberData(nameof(Reference))]
    public void MeetsTheReference_BalancedAndUnequal(double a, double b, double reference)
    {
        double got = SpecialFunctions.LogBeta(a, b), error = Math.Abs(got - reference) / Math.Max(1.0, Math.Abs(reference));
        output.WriteLine($"ln B({a:R}, {b:R}) = {got:R}, reference {reference:R}, relative {error:G3}");
        Assert.True(error <= 4e-15, $"relative {error:G3}");
    }

    [Fact]
    public void FarApart_ItIsMinusTheLogOfTheLarger_AndSymmetric()
    {
        // ln B(x, 1) = −ln x exactly; the three-term form returned 0 from x ≈ 10¹⁶.
        foreach (double x in new[] { 10.0, 1e5, 1e10, 1e16, 1e50, 1e200, 1e300, double.MaxValue })
        {
            Assert.True(Math.Abs(SpecialFunctions.LogBeta(x, 1.0) + Math.Log(x)) <= 4e-16 * Math.Log(x), $"ln B({x:R}, 1)");
            Assert.Equal(SpecialFunctions.LogBeta(x, 1.0), SpecialFunctions.LogBeta(1.0, x));
        }
        foreach ((double a, double b) in new[] { (0.3, 7.0), (12.0, 1e9), (2e15, 3e15), (1e300, 1.7e308), (9e307, 9.5e307) })
        {
            double ab = SpecialFunctions.LogBeta(a, b);
            Assert.True(double.IsFinite(ab), $"ln B({a:R}, {b:R}) = {ab:R}");
            Assert.Equal(ab, SpecialFunctions.LogBeta(b, a));
        }
        // At (M, M), M the double's largest, ln B is about −2M·ln 2, past the double's range: −∞ is its value.
        Assert.Equal(double.NegativeInfinity, SpecialFunctions.LogBeta(double.MaxValue, double.MaxValue));
    }

    [Fact]
    public void ShapesNotFiniteAndPositive_TakeTheThreeTermForm()
    {
        foreach ((double a, double b) in new[] { (double.NaN, 1.0), (1.0, double.NaN), (0.0, 1.0), (-0.5, 2.0), (double.PositiveInfinity, 1.0) })
            Assert.Equal(SpecialFunctions.LogGamma(a) + SpecialFunctions.LogGamma(b) - SpecialFunctions.LogGamma(a + b), SpecialFunctions.LogBeta(a, b));
    }
}
