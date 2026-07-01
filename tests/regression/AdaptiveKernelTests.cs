using System;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Adaptive proposal-scale tuning: from a deliberately bad (too-tight) spread τ — where the relocate move barely
/// moves and over-accepts — the Robbins–Monro <see cref="AdaptiveLocalBetaKernel.Adapt"/> step, fed the
/// engine's per-move acceptance (<see cref="ReversibleJumpChain{TState}.MoveStats"/>), drives the relocate
/// acceptance to the 1-D target and loosens τ. Validates the mechanism a free-knot run will use in its burn-in.
/// </summary>
public sealed class AdaptiveKernelTests
{
    private readonly ITestOutputHelper _out;
    public AdaptiveKernelTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double RelocateAcceptance(ReversibleJumpChain<KnotConfig> chain)
    {
        MoveStat r = chain.MoveStats().First(m => m.Key == "relocate");
        return r.Attempts == 0 ? double.NaN : (double)r.Accepted / r.Attempts;
    }

    private static (long att, long acc) Relocate(ReversibleJumpChain<KnotConfig> chain)
    {
        MoveStat r = chain.MoveStats().First(m => m.Key == "relocate");
        return (r.Attempts, r.Accepted);
    }

    [Fact]
    public void AdaptiveKernel_DrivesRelocateAcceptance_FromBadStart()
    {
        var noise = new Xoshiro256PlusPlus(seed: 17);
        int n = 80;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(3.0 * Math.PI * x[i]) + 0.08 * Gaussian(noise);
        }

        const double target = 0.44;
        var kernel = new AdaptiveLocalBetaKernel(initialTau: 2000.0);   // far too tight ⇒ relocate over-accepts
        var spline = new SplineTarget(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0), x, y);
        var moves = new IRjMove<KnotConfig>[]
        {
            new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel),
        };
        var chain = new ReversibleJumpChain<KnotConfig>(
            moves, spline, new KnotConfig(new[] { 0.25, 0.5, 0.75 }), new Xoshiro256PlusPlus(seed: 99));

        for (int s = 0; s < 1000; s++) chain.Step();   // settle into a few-knot region
        double startTau = kernel.Tau;

        // Burn-in adaptation: per round, take the relocate acceptance over the round and nudge τ (vanishing step).
        (long att, long acc) prev = Relocate(chain);
        for (int round = 0; round < 60; round++)
        {
            for (int s = 0; s < 400; s++) chain.Step();
            (long att, long acc) now = Relocate(chain);
            long dA = now.att - prev.att, dC = now.acc - prev.acc;
            prev = now;
            if (dA > 0) kernel.Adapt((double)dC / dA, stepSize: 3.0 / Math.Sqrt(round + 1), target: target);
        }

        // Freeze; measure relocate acceptance over a fresh window at the tuned τ.
        (long att, long acc) a = Relocate(chain);
        for (int s = 0; s < 8000; s++) chain.Step();
        (long att, long acc) b = Relocate(chain);
        double finalAcc = (double)(b.acc - a.acc) / (b.att - a.att);

        _out.WriteLine($"[adaptive-τ] τ: {startTau:F0} → {kernel.Tau:F1}; relocate acc → {finalAcc:F3} (target {target})");

        Assert.True(kernel.Tau < startTau * 0.5, $"τ should loosen from the too-tight start ({kernel.Tau:F1} vs {startTau:F0})");
        Assert.InRange(finalAcc, 0.30, 0.60);   // converged near the 1-D optimum
    }

    [Fact]
    public void AdaptiveBurn_InEnsemble_RecoversMixing_FromBadTau()
    {
        var noise = new Xoshiro256PlusPlus(seed: 41);
        int n = 90;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = Math.Sin(3.0 * Math.PI * x[i]);
            y[i] = f[i] + 0.10 * Gaussian(noise);
        }

        var basis = new SplineBasis(3);
        var prior = new PoissonPrior(4.0);

        // A fixed, far-too-tight τ: relocate barely moves knots ⇒ the peak mixes slowly.
        BarsResult fixedBad = new BarsEnsemble(basis, new WeightedNormalModel(), prior, new LocalBetaKernel(2000.0))
            .Run(x, y, grid: x, chains: 4, masterSeed: 5, burn: 1500, samples: 3000);
        // Same bad start, but adapt τ during burn-in.
        BarsResult adaptive = new BarsEnsemble(basis, new WeightedNormalModel(), prior, new LocalBetaKernel(50.0), adaptKernelTau: 2000.0)
            .Run(x, y, grid: x, chains: 4, masterSeed: 5, burn: 1500, samples: 3000);

        double Mse(double[] fit) { double s = 0.0; for (int i = 0; i < n; i++) { double d = fit[i] - f[i]; s += d * d; } return s / n; }
        double mseFixed = Mse(fixedBad.Fit);
        double mseAdaptive = Mse(adaptive.Fit);
        _out.WriteLine($"[adaptive-ens] fixedBad: mse={mseFixed:F5} peakESS={fixedBad.Peak.LocationEss:F0} ESS(k)={fixedBad.EssKnots:F0}");
        _out.WriteLine($"[adaptive-ens] adaptive: mse={mseAdaptive:F5} peakESS={adaptive.Peak.LocationEss:F0} ESS(k)={adaptive.EssKnots:F0}");

        Assert.True(mseAdaptive < 0.01, $"adaptive fit MSE {mseAdaptive:F5} should recover the curve");
        Assert.True(adaptive.Peak.LocationEss > fixedBad.Peak.LocationEss,
            $"adaptive peak ESS {adaptive.Peak.LocationEss:F0} should beat the too-tight fixed kernel {fixedBad.Peak.LocationEss:F0}");
    }
}
