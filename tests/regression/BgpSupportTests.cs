using System;
using Maths.Regression.Bgp;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The empirical-Bayes bandwidth prior's support <c>(γ₁T_n², 1]</c> must be nonempty and finite in the sampled
/// coordinate log t before any chain starts. Without the check, two inputs at 0 and 100 put the lower bound at
/// 2500: every warm-started chain sits at a bandwidth of zero prior density, never accepts a move, and reports a
/// constant R̂ that satisfies a stopping rule. Boundary and constant-coordinate inputs pin the edges of the rule.
/// </summary>
public sealed class BgpSupportTests
{
    private static EmpiricalBayesBandwidthPrior Prior(double[,] x)
    {
        int n = x.GetLength(0);
        var y = new double[n];
        for (int i = 0; i < n; i++) y[i] = 0.3 * Math.Cos(i);
        var gp = new GpRegression(x, y, sigma2: 0.01, new SquaredExponentialKernel());
        return new EmpiricalBayesBandwidthPrior(x, gp);
    }

    private static double[,] Column(params double[] values)
    {
        var x = new double[values.Length, 1];
        for (int i = 0; i < values.Length; i++) x[i, 0] = values[i];
        return x;
    }

    [Fact]
    public void EmptySupport_IsRejected_BeforeAnyChainStarts()
    {
        // T_n = 100 with the default γ₁ = 1/4 gives a lower bound of 2500 against the fixed upper bound 1.
        var ex = Assert.Throws<ArgumentException>(() => Prior(Column(0.0, 100.0)));
        Assert.Contains("2500", ex.Message);
        Assert.Contains("support", ex.Message);
    }

    [Fact]
    public void LowerBoundAtTheUpperBound_IsEmpty_AndJustBelowIsAccepted()
    {
        // T_n = 2 ⇒ γ₁T_n² = 1: the support (1, 1] is empty.
        Assert.Throws<ArgumentException>(() => Prior(Column(0.0, 2.0)));

        // T_n = 1.9 ⇒ support (0.9025, 1]: narrow but valid, with zero density on both edges.
        EmpiricalBayesBandwidthPrior prior = Prior(Column(0.0, 1.9));
        Assert.Equal(0.9025, prior.LowerBound, 12);
        Assert.True(double.IsNegativeInfinity(prior.LogDensity(prior.LowerBound)));
        Assert.True(double.IsNegativeInfinity(prior.LogDensity(1.0 + 1e-12)));
        Assert.True(double.IsFinite(prior.LogDensity(0.95)));
        Assert.True(double.IsFinite(prior.LogDensity(1.0)));
    }

    [Fact]
    public void NarrowSupportNearTheUpperBound_StartsTheChainsInside()
    {
        // T_n = 1.9995 ⇒ support (0.99950006…, 1]; a warm start clamped to 0.999 would sit below it.
        double[,] x = Column(0.0, 1.9995);
        var gp = new GpRegression(x, new[] { 0.3, -0.1 }, sigma2: 0.01, new SquaredExponentialKernel());
        var prior = new EmpiricalBayesBandwidthPrior(x, gp);
        BgpResult r = new BgpSampler(gp, prior).Run(x, draws: 200, burn: 50, proposalSd: 1e-4, seed: 3, chains: 2);

        Assert.True(r.BandwidthLo > prior.LowerBound && r.BandwidthHi <= 1.0,
            $"draws [{r.BandwidthLo}, {r.BandwidthHi}] must lie in ({prior.LowerBound}, 1]");
        Assert.True(r.AcceptanceRate > 0.0, "the chains must move inside the support");
    }

    [Fact]
    public void StartOutsideTheSupport_IsRejected()
    {
        double[,] x = Column(0.0, 0.3, 0.7, 1.0);
        var gp = new GpRegression(x, new[] { 0.1, 0.4, -0.2, 0.0 }, sigma2: 0.01, new SquaredExponentialKernel());
        var sampler = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(x, gp));
        Assert.Throws<ArgumentOutOfRangeException>(() => sampler.Run(x, tInit: 2.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => sampler.Run(x, tInit: 1e-6));
    }

    [Theory]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    [InlineData(1e200)]   // finite input whose squared distance overflows: T_n = ∞
    public void NonFiniteSupport_IsRejected(double bad)
    {
        Assert.Throws<ArgumentException>(() => Prior(Column(0.0, bad, 1.0)));
    }

    [Fact]
    public void DuplicatedInputs_GiveZeroTn_AndAreRejected()
    {
        // Every point's nearest other point is a duplicate, so T_n = 0 and log t is unbounded below.
        var ex = Assert.Throws<ArgumentException>(() => Prior(Column(0.0, 0.0, 1.0, 1.0)));
        Assert.Contains("support", ex.Message);
    }

    [Fact]
    public void SingleInput_IsRejected_WithAClearMessage()
    {
        var ex = Assert.Throws<ArgumentException>(() => Prior(Column(0.5)));
        Assert.Contains("two", ex.Message);
    }

    [Fact]
    public void ConstantCoordinate_LeavesPriorAndFitUnchanged()
    {
        const int n = 40;
        var x1 = new double[n, 1];
        var x2 = new double[n, 2];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double u = (i + 0.5) / n;
            x1[i, 0] = u;
            x2[i, 0] = u;
            x2[i, 1] = 5.0;   // constant coordinate
            y[i] = Math.Sin(2.0 * Math.PI * u);
        }
        var gp1 = new GpRegression(x1, y, sigma2: 0.01, new SquaredExponentialKernel());
        var gp2 = new GpRegression(x2, y, sigma2: 0.01, new SquaredExponentialKernel());
        var p1 = new EmpiricalBayesBandwidthPrior(x1, gp1);
        var p2 = new EmpiricalBayesBandwidthPrior(x2, gp2);

        Assert.Equal(p1.Tn, p2.Tn);
        Assert.Equal(p1.LowerBound, p2.LowerBound);
        foreach (double t in new[] { 0.01, 0.1, 0.5, 1.0 })
        {
            Assert.Equal(p1.LogDensity(t), p2.LogDensity(t));
            Assert.Equal(gp1.Fit(t).LogMarginal, gp2.Fit(t).LogMarginal);
        }

        var test1 = new double[,] { { 0.25 }, { 0.8 } };
        var test2 = new double[,] { { 0.25, 5.0 }, { 0.8, 5.0 } };
        Assert.Equal(gp1.PredictMean(gp1.Fit(0.05), test1), gp2.PredictMean(gp2.Fit(0.05), test2));
    }
}
