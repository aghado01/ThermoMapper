using System;
using Maths.Distributions;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// <see cref="SpecialFunctions.Log1p"/> where Goldberg's correction ln(u)·x/(u − 1) overflows its product, x above about 2.5·10³⁰⁵: it
/// returned +∞ there until 2026-11-14, found by the log-gamma ratio's additivity sweep; and below that range, unchanged.
/// </summary>
public sealed class Log1pTests
{
    [Theory]
    [InlineData(2.4e305)]
    [InlineData(2.6e305)]
    [InlineData(1e306)]
    [InlineData(1e308)]
    [InlineData(double.MaxValue)]
    public void HugeArguments_AreTheLogOfTheArgument(double x)
    {
        double got = SpecialFunctions.Log1p(x);
        Assert.True(double.IsFinite(got), $"Log1p({x:R}) = {got:R}");
        Assert.Equal(Math.Log(x), got, 12);
    }

    [Fact]
    public void BelowTheOverflow_ItIsGoldbergsForm()
    {
        foreach (double x in new[] { -0.5, -1e-300, 1e-300, 1e-17, 1e-8, 0.3, 7.0, 1e15, 2.0e305 })
        {
            double u = 1.0 + x;
            double goldberg = u == 1.0 ? x : Math.Log(u) * x / (u - 1.0);
            Assert.Equal(goldberg, SpecialFunctions.Log1p(x));
        }
        Assert.Equal(double.PositiveInfinity, SpecialFunctions.Log1p(double.PositiveInfinity));
        Assert.Equal(double.NegativeInfinity, SpecialFunctions.Log1p(-1.0));
    }
}
