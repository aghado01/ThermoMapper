using System;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// Parallel tempering on its own terms. On a well-separated bimodal target a single chain started in one mode can't
/// cross the barrier and stays stuck; the tempered ladder's hot replicas traverse it freely and ferry the cold replica
/// between modes, so the cold chain visits both ≈ equally. Swaps exchange the replicas' cached evaluations: no state
/// is re-evaluated for a swap.
/// </summary>
public sealed class ParallelTemperingTests
{
    private readonly ITestOutputHelper _out;
    public ParallelTemperingTests(ITestOutputHelper output) => _out = output;

    // 0.5·N(−3, 0.5²) + 0.5·N(+3, 0.5²): two modes separated by a deep barrier.
    private static readonly ISamplingTarget<double> Bimodal = SamplingTarget.FromLogDensity<double>(x =>
    {
        double l1 = LogNormal(x, -3.0, 0.5);
        double l2 = LogNormal(x, 3.0, 0.5);
        double m = Math.Max(l1, l2);
        return m + Math.Log(0.5 * Math.Exp(l1 - m) + 0.5 * Math.Exp(l2 - m));
    });

    private static double LogNormal(double x, double mu, double sd)
    {
        double z = (x - mu) / sd;
        return -0.5 * z * z - Math.Log(sd * Math.Sqrt(2.0 * Math.PI));
    }

    [Fact]
    public void Tempering_VisitsBothModes_WhereSingleChainSticks()
    {
        ISamplingTarget<double> target = Bimodal;
        var moves = new IRjMove<double>[] { new RandomWalk(0.5) };
        const int draws = 100_000;

        // Single cold chain started in the left mode — it stays stuck.
        var single = new ReversibleJumpChain<double>(moves, target, -3.0, new Xoshiro256PlusPlus(1));
        for (int i = 0; i < 5_000; i++) single.Step();
        int singleRight = 0;
        for (int i = 0; i < draws; i++) if (single.Step() > 0.0) singleRight++;
        double singleFrac = (double)singleRight / draws;

        // Tempered ladder — the cold replica visits both modes.
        double[] betas = ParallelTempering<double>.GeometricLadder(levels: 6, betaMin: 0.02);
        var pt = new ParallelTempering<double>(_ => moves, target, -3.0, betas, masterSeed: 7);
        for (int i = 0; i < 5_000; i++) pt.Step();
        int tempRight = 0;
        for (int i = 0; i < draws; i++) if (pt.Step() > 0.0) tempRight++;
        double tempFrac = (double)tempRight / draws;

        _out.WriteLine($"[temper] single-chain right-mode frac={singleFrac:F3} (stuck); " +
                       $"tempering cold frac={tempFrac:F3}; swap acc={(double)pt.SwapAccepts / pt.SwapAttempts:F3}");

        Assert.True(singleFrac < 0.05, $"single chain should stay stuck in one mode (got {singleFrac:F3})");
        Assert.InRange(tempFrac, 0.35, 0.65);   // cold replica visits both modes ≈ equally
    }

    [Fact]
    public void Swaps_ExchangeEvaluations_WithoutReevaluating()
    {
        int evaluations = 0;
        ISamplingTarget<int> target = SamplingTarget.FromLogDensity<int>(s => { evaluations++; return -0.3 * Math.Abs(s); });
        IRjMove<int>[] MovesFor(int replica) => [new UnitStep("up", "down", +1), new UnitStep("down", "up", -1)];
        double[] betas = ParallelTempering<int>.GeometricLadder(levels: 4, betaMin: 0.05);
        var pt = new ParallelTempering<int>(MovesFor, target, 0, betas, masterSeed: 5);
        Assert.Equal(4, evaluations);   // each replica's start

        const int steps = 2000;
        for (int i = 0; i < steps; i++) pt.Step();

        Assert.True(pt.SwapAccepts > 0);
        Assert.Equal(4 + 4 * steps, evaluations);   // one candidate per replica per step; swaps evaluate nothing
        Assert.Equal(1.0, pt.ColdChain.Temperature.InverseTemperature);
        Assert.Equal(-0.3 * Math.Abs(pt.Cold), pt.ColdChain.CurrentLogDensity);   // recomputed at the cold temperature
    }

    [Theory]
    [InlineData(1.5)]
    [InlineData(0.0)]
    [InlineData(double.NaN)]
    public void Ladder_OutsideTheUnitInterval_IsRejected(double beta)
    {
        IRjMove<int>[] MovesFor(int replica) => [new UnitStep("up", "down", +1), new UnitStep("down", "up", -1)];
        Assert.Throws<ArgumentException>(() => new ParallelTempering<int>(
            MovesFor, SamplingTarget.FromLogDensity<int>(_ => 0.0), 0, [1.0, beta], masterSeed: 1));
    }

    // Symmetric random walk; too small to jump the barrier at β = 1.
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
            double g = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
            return new Proposal<double>(current + _sd * g, 0.0, 0.0);   // symmetric ⇒ no proposal/Jacobian correction
        }
    }

    private sealed class UnitStep : IRjMove<int>
    {
        private readonly int _delta;
        public UnitStep(string key, string reverseKey, int delta) { Key = key; ReverseKey = reverseKey; _delta = delta; }
        public string Key { get; }
        public string ReverseKey { get; }
        public double Weight(int state) => 1.0;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => new(current + _delta, 0.0, 0.0);
    }
}
