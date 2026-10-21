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
/// is re-evaluated for a swap. On an enumerable space with a likelihood that is not linear in β, every rung samples its
/// own path's law exactly, which the β-linear swap shortcut does not achieve on the likelihood path.
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
        var pt = new ParallelTempering<double>(_ => moves, target, -3.0, betas, TemperingPath.Powered, masterSeed: 7);
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
        var pt = new ParallelTempering<int>(MovesFor, target, 0, betas, TemperingPath.Powered, masterSeed: 5);
        Assert.Equal(4, evaluations);   // each replica's start

        const int steps = 2000;
        for (int i = 0; i < steps; i++) pt.Step();

        Assert.True(pt.SwapAccepts > 0);
        Assert.Equal(4 + 4 * steps, evaluations);   // one candidate per replica per step; swaps evaluate nothing
        Assert.Equal(1.0, pt.ColdChain.Temperature.InverseTemperature);
        Assert.Equal(-0.3 * Math.Abs(pt.Cold), pt.ColdChain.CurrentLogDensity);   // recomputed at the cold temperature
    }

    [Theory]
    [InlineData(TemperingPath.Likelihood)]
    [InlineData(TemperingPath.Powered)]
    public void EveryRung_SamplesItsPathsLaw(TemperingPath path)
    {
        // States 0..5 under a prior p and a likelihood ℓ_β = 3·ln(1 + β·a) − β·c, which is not linear in β (as for a
        // collapsed model). Rung β's law is p·e^{ℓ_β} on the likelihood path and (p·e^{ℓ₁})^β on the powered path.
        double[] betas = path == TemperingPath.Likelihood ? [1.0, 0.5, 0.2, 0.0] : [1.0, 0.5, 0.2, 0.05];
        IRjMove<int>[] MovesFor(int replica) => [new UnitStep("up", "down", +1), new UnitStep("down", "up", -1)];
        var pt = new ParallelTempering<int>(MovesFor, new SixStates(), 0, betas, path, masterSeed: 13);
        for (int i = 0; i < 5_000; i++) pt.Step();

        const int n = 400_000;
        var counts = new double[betas.Length, SixStates.States];
        for (int i = 0; i < n; i++)
        {
            pt.Step();
            for (int r = 0; r < betas.Length; r++) counts[r, pt.Replicas[r].Current] += 1.0;
        }

        for (int r = 0; r < betas.Length; r++)
        {
            double[] expected = SixStates.Law(betas[r], path);
            double tv = 0.0;
            for (int k = 0; k < SixStates.States; k++)
            {
                double observed = counts[r, k] / n;
                tv += 0.5 * Math.Abs(observed - expected[k]);
                Assert.True(Math.Abs(observed - expected[k]) < 0.01,
                    $"{path} rung β = {betas[r]}, state {k}: {observed:F4} vs {expected[k]:F4}");
            }
            _out.WriteLine($"[pt-law] {path} β={betas[r]}: TV={tv:F4}");
        }
    }

    public static TheoryData<double[], TemperingPath> InvalidLadders => new()
    {
        { new[] { 1.0 }, TemperingPath.Likelihood },               // one rung
        { new[] { 0.9, 0.5 }, TemperingPath.Likelihood },          // the cold rung is not 1
        { new[] { 1.0, 1.5 }, TemperingPath.Likelihood },          // above 1
        { new[] { 1.0, double.NaN }, TemperingPath.Likelihood },
        { new[] { 1.0, 0.3, 0.5 }, TemperingPath.Likelihood },     // not descending
        { new[] { 1.0, 0.5, 0.5 }, TemperingPath.Powered },        // not strictly descending
        { new[] { 1.0, -0.1 }, TemperingPath.Likelihood },
        { new[] { 1.0, 0.0 }, TemperingPath.Powered },             // β = 0 is the bare reference measure
    };

    [Theory]
    [MemberData(nameof(InvalidLadders))]
    public void Ladder_Invalid_IsRejected(double[] betas, TemperingPath path)
    {
        IRjMove<int>[] MovesFor(int replica) => [new UnitStep("up", "down", +1), new UnitStep("down", "up", -1)];
        Assert.Throws<ArgumentException>(() => new ParallelTempering<int>(
            MovesFor, SamplingTarget.FromLogDensity<int>(_ => 0.0), 0, betas, path, masterSeed: 1));
    }

    [Fact]
    public void Ladder_DeclaresItsPathOnEveryRung()
    {
        IRjMove<int>[] MovesFor(int replica) => [new UnitStep("up", "down", +1), new UnitStep("down", "up", -1)];
        double[] betas = [1.0, 0.5, 0.0];   // β = 0 is the prior on the likelihood path
        var pt = new ParallelTempering<int>(MovesFor, SamplingTarget.FromLogDensity<int>(_ => 0.0), 0, betas,
                                            TemperingPath.Likelihood, masterSeed: 1);
        Assert.Equal(TemperingPath.Likelihood, pt.Path);
        Assert.Equal(betas.Length, pt.Replicas.Count);
        for (int r = 0; r < betas.Length; r++)
            Assert.Equal(new Temperature(betas[r], TemperingPath.Likelihood), pt.Replicas[r].Temperature);
        Assert.Same(pt.ColdChain, pt.Replicas[0]);
    }

    // States 0..5; a likelihood that is not linear in β.
    private sealed class SixStates : ISamplingTarget<int>
    {
        public const int States = 6;
        private static readonly double[] Prior = [0.30, 0.10, 0.25, 0.05, 0.20, 0.10];
        private static readonly double[] A = [0.0, 5.0, 20.0, 1.0, 50.0, 2.0];
        private static readonly double[] C = [0.0, 2.0, 6.0, 0.5, 9.0, 1.0];

        private static double LogLikelihood(int k, double beta) => 3.0 * Math.Log(1.0 + beta * A[k]) - beta * C[k];

        public ITargetEvaluation Evaluate(int k) => new Evaluation(k);

        public static double[] Law(double beta, TemperingPath path)
        {
            var w = new double[States];
            double total = 0.0;
            for (int k = 0; k < States; k++)
            {
                double logPrior = Math.Log(Prior[k]);
                w[k] = Math.Exp(path == TemperingPath.Likelihood
                    ? logPrior + LogLikelihood(k, beta)
                    : beta * (logPrior + LogLikelihood(k, 1.0)));
                total += w[k];
            }
            for (int k = 0; k < States; k++) w[k] /= total;
            return w;
        }

        private sealed class Evaluation : ITargetEvaluation
        {
            private readonly int _k;
            public Evaluation(int k) => _k = k;
            public double LogPrior => _k is >= 0 and < States ? Math.Log(Prior[_k]) : double.NegativeInfinity;
            public double LogLikelihood(double inverseTemperature)
                => _k is >= 0 and < States ? SixStates.LogLikelihood(_k, inverseTemperature) : double.NegativeInfinity;
        }
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
