using System;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// The neutral target contract (step 3.1). An evaluation splits a state's log density into an untempered log prior
/// and a log likelihood at any inverse temperature b; a <see cref="Temperature"/> composes them along the powered
/// path <c>b·(ln prior + ℓ₁)</c> or the likelihood path <c>ln prior + ℓ_b</c>, which agree at b = 1. Off the support no
/// likelihood work is done. A chain samples the tempered density it is given: on a Gaussian prior and likelihood the
/// two paths have different, known laws.
/// </summary>
public sealed class SamplingTargetTests
{
    [Fact]
    public void Paths_AgreeAtOne_AndDifferBelowIt()
    {
        var e = new SeparableEvaluation(-1.25, -2.75);
        Assert.Equal(-4.0, new Temperature(1.0, TemperingPath.Powered).LogDensity(e));
        Assert.Equal(-4.0, new Temperature(1.0, TemperingPath.Likelihood).LogDensity(e));
        Assert.Equal(-2.0, new Temperature(0.5, TemperingPath.Powered).LogDensity(e));
        Assert.Equal(-1.25 - 1.375, new Temperature(0.5, TemperingPath.Likelihood).LogDensity(e));
        Assert.Equal(-4.0, Temperature.LogPosterior(e));
    }

    [Fact]
    public void OffTheSupport_NoLikelihoodWorkIsDone()
    {
        var e = new RefusingEvaluation();
        Assert.Equal(double.NegativeInfinity, new Temperature(0.3, TemperingPath.Powered).LogDensity(e));
        Assert.Equal(double.NegativeInfinity, new Temperature(0.3, TemperingPath.Likelihood).LogDensity(e));

        int likelihoodCalls = 0;
        ISamplingTarget<int> target = SamplingTarget.Separable<int>(
            s => s < 0 ? double.NegativeInfinity : 0.0, s => { likelihoodCalls++; return -s; });
        Assert.Equal(double.NegativeInfinity, target.LogPosterior(-1));
        Assert.Equal(0, likelihoodCalls);
        Assert.Equal(-2.0, target.LogPosterior(2));
        Assert.Equal(1, likelihoodCalls);
    }

    [Fact]
    public void ZeroLikelihood_StaysZeroAtZeroTemperature()
    {
        var e = new SeparableEvaluation(0.0, double.NegativeInfinity);
        Assert.Equal(double.NegativeInfinity, e.LogLikelihood(0.0));   // not 0·(−∞) = NaN
        Assert.Equal(double.NegativeInfinity, new Temperature(0.0, TemperingPath.Likelihood).LogDensity(e));
        Assert.Equal(-1.5, new SeparableEvaluation(0.0, -3.0).LogLikelihood(0.5));
    }

    [Theory]
    [InlineData(1.5)]
    [InlineData(-0.1)]
    [InlineData(double.NaN)]
    public void Temperature_OutsideTheUnitInterval_IsRejected(double b)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new Temperature(b, TemperingPath.Likelihood));
    }

    [Fact]
    public void Temperature_UndefinedPath_IsRejected()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new Temperature(0.5, (TemperingPath)7));
    }

    [Fact]
    public void SetState_RecomputesTheDensityAtTheChainsOwnTemperature()
    {
        ISamplingTarget<int> target = SamplingTarget.Separable<int>(_ => -1.0, s => -2.0 * s);
        IRjMove<int>[] moves = [new Stay()];
        var chain = new ReversibleJumpChain<int>(moves, target, 0, new Xoshiro256PlusPlus(1),
                                                 new Temperature(0.5, TemperingPath.Likelihood));
        chain.SetState(3, target.Evaluate(3));
        Assert.Equal(3, chain.Current);
        Assert.Equal(-1.0 + 0.5 * -6.0, chain.CurrentLogDensity);
    }

    [Theory]
    [InlineData(TemperingPath.Likelihood, 0.4, 0.8)]   // precision 1 + b, mean 2b/(1 + b)
    [InlineData(TemperingPath.Powered, 1.0, 2.0)]      // precision 2b, mean 1
    public void Chain_SamplesTheTemperedDensityItIsGiven(TemperingPath path, double mean, double variance)
    {
        // Prior N(0, 1), likelihood of one observation y = 2 under N(x, 1), at b = 0.25.
        ISamplingTarget<double> target = SamplingTarget.Separable<double>(x => -0.5 * x * x, x => -0.5 * (x - 2.0) * (x - 2.0));
        IRjMove<double>[] moves = [new RandomWalk(1.5)];
        var chain = new ReversibleJumpChain<double>(moves, target, 0.0, new Xoshiro256PlusPlus(17),
                                                    new Temperature(0.25, path));
        for (int i = 0; i < 5_000; i++) chain.Step();

        const int n = 400_000;
        double sum = 0.0, sumSq = 0.0;
        for (int i = 0; i < n; i++)
        {
            double x = chain.Step();
            sum += x;
            sumSq += x * x;
        }
        double m = sum / n;
        Assert.InRange(m, mean - 0.03, mean + 0.03);
        Assert.InRange(sumSq / n - m * m, variance - 0.06, variance + 0.06);
    }

    private sealed class RefusingEvaluation : ITargetEvaluation
    {
        public double LogPrior => double.NegativeInfinity;
        public double LogLikelihood(double inverseTemperature) => throw new InvalidOperationException("Likelihood work was done.");
    }

    private sealed class Stay : IRjMove<int>
    {
        public string Key => "stay";
        public string ReverseKey => "stay";
        public double Weight(int state) => 1.0;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => new(current, 0.0, 0.0);
    }

    private sealed class RandomWalk : IRjMove<double>
    {
        private readonly double _sd;
        public RandomWalk(double sd) => _sd = sd;
        public string Key => "rw";
        public string ReverseKey => "rw";
        public double Weight(double state) => 1.0;
        public Proposal<double>? Propose(double current, Xoshiro256PlusPlus rng)
        {
            double u1 = 1.0 - rng.NextDouble();
            double u2 = 1.0 - rng.NextDouble();
            return new Proposal<double>(current + _sd * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2), 0.0, 0.0);
        }
    }
}
