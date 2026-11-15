using System;
using Maths.Distributions;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// <see cref="SpecialFunctions.LogGamma"/> below about 5.6·10⁻³⁰⁹, where the reflection's π/sin(πx) ≈ 1/x overflowed: it returned +∞
/// there until 2026-11-15, found by the ln B qualification of the defects brief's slice 9 (ln B(5·10⁻³²⁴, 1) was +∞, its value 744.44);
/// and above that range, unchanged.
/// </summary>
public sealed class LogGammaTests
{
    // ln Γ(x) from mpmath 1.3.0 at 60 digits; the first three below the overflow, 5.5·10⁻³⁰⁹ just below it, 6·10⁻³⁰⁹ just above.
    [Theory]
    [InlineData(4.9406564584124654e-324, 744.4400719213812)]
    [InlineData(1e-320, 736.8272408909739)]
    [InlineData(1e-310, 713.8013788281542)]
    [InlineData(5e-309, 709.889355822726)]
    [InlineData(5.5e-309, 709.7940456429217)]
    [InlineData(6e-309, 709.7070342659321)]
    [InlineData(2.2250738585072014e-308, 708.3964185322641)]
    public void AtSubnormals_ItIsFinite_AndTheReference(double x, double reference)
    {
        double got = SpecialFunctions.LogGamma(x);
        Assert.True(Math.Abs(got - reference) <= 2.3e-16 * reference, $"ln Γ({x:R}) = {got:R}, reference {reference:R}");
    }

    [Fact]
    public void AboveTheOverflow_ItIsTheReflectionAsBefore_AndZeroIsThePole()
    {
        foreach (double x in new[] { 6e-309, 1e-300, 1e-17, 0.25, 0.49 })
            Assert.Equal(Math.Log(Math.PI / Math.Sin(Math.PI * x)) - SpecialFunctions.LogGamma(1.0 - x), SpecialFunctions.LogGamma(x));
        Assert.Equal(double.PositiveInfinity, SpecialFunctions.LogGamma(0.0));
    }
}
