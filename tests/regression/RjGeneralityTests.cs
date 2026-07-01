using System;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Green-generality of the reversible-jump engine, on its own terms — no model consumer involved. A trans-
/// dimensional target (the dimension k itself, with a known Poisson stationary law) is sampled by birth/death
/// moves with deliberately <b>asymmetric</b> selection weights. Recovering the Poisson mean and variance is
/// possible only if the engine applies the move-selection ratio j_{m'}(x')/j_m(x); ignoring it (the special-case
/// "weights cancel" assumption) would bias the dimension upward toward the birth-favoured drift.
/// </summary>
public sealed class RjGeneralityTests
{
    private readonly ITestOutputHelper _out;
    public RjGeneralityTests(ITestOutputHelper output) => _out = output;

    // Target over the count k: Poisson(μ) up to the k-independent constant (mean = variance = μ).
    private sealed class PoissonDimension : IRjTarget<int>
    {
        private readonly double _logMu;
        public PoissonDimension(double mu) => _logMu = Math.Log(mu);
        public double LogPosterior(int k)
        {
            double logFactorial = 0.0;
            for (int i = 2; i <= k; i++) logFactorial += Math.Log(i);
            return k * _logMu - logFactorial;
        }
    }

    // Deterministic up/down moves (no auxiliary, no Jacobian): they isolate the selection-ratio factor.
    private sealed class Birth : IRjMove<int>
    {
        private readonly double _w;
        public Birth(double w) => _w = w;
        public string Key => "birth";
        public string ReverseKey => "death";
        public double Weight(int k) => _w;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => new Proposal<int>(current + 1, 0.0, 0.0);
    }

    private sealed class Death : IRjMove<int>
    {
        private readonly double _w;
        public Death(double w) => _w = w;
        public string Key => "death";
        public string ReverseKey => "birth";
        public double Weight(int k) => _w;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng)
            => current == 0 ? null : new Proposal<int>(current - 1, 0.0, 0.0);
    }

    [Fact]
    public void AsymmetricSelection_RecoversPoissonDimensionLaw()
    {
        const double mu = 4.0;
        var target = new PoissonDimension(mu);
        var moves = new IRjMove<int>[] { new Birth(0.7), new Death(0.3) };   // 0.7 ≠ 0.3 — the test's whole point
        var chain = new ReversibleJumpChain<int>(moves, target, 0, new Xoshiro256PlusPlus(seed: 11));

        for (int i = 0; i < 20_000; i++) chain.Step();   // burn-in

        const int n = 400_000;
        double sum = 0.0, sumSq = 0.0;
        for (int i = 0; i < n; i++)
        {
            int k = chain.Step();
            sum += k;
            sumSq += (double)k * k;
        }
        double mean = sum / n;
        double var = sumSq / n - mean * mean;
        _out.WriteLine($"[rj-general] asymmetric w=(0.7,0.3) → mean={mean:F3} var={var:F3} (Poisson μ={mu})");

        // Poisson(μ): mean = var = μ. Both must hold despite the asymmetric weights.
        Assert.InRange(mean, mu - 0.2, mu + 0.2);
        Assert.InRange(var, mu - 0.4, mu + 0.4);
    }

    [Fact]
    public void EqualSelection_AlsoRecoversLaw()
    {
        // Sanity twin: with symmetric weights the selection ratio is 1, so the same law must come out.
        const double mu = 4.0;
        var target = new PoissonDimension(mu);
        var moves = new IRjMove<int>[] { new Birth(0.5), new Death(0.5) };
        var chain = new ReversibleJumpChain<int>(moves, target, 0, new Xoshiro256PlusPlus(seed: 7));

        for (int i = 0; i < 20_000; i++) chain.Step();
        const int n = 400_000;
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += chain.Step();
        double mean = sum / n;
        _out.WriteLine($"[rj-general] symmetric w=(0.5,0.5) → mean={mean:F3} (Poisson μ={mu})");
        Assert.InRange(mean, mu - 0.2, mu + 0.2);
    }
}
