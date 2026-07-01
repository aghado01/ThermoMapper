using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Parallel tempering on its own terms — a well-separated bimodal target. A single chain started in one mode
/// can't cross the barrier and stays stuck; the tempered ladder's hot replicas traverse it freely and ferry the
/// cold replica between modes, so the cold chain visits both ≈ equally. The multimodal-escape property the
/// multi-peak T_c posteriors will need.
/// </summary>
public sealed class TemperingTests
{
    private readonly ITestOutputHelper _out;
    public TemperingTests(ITestOutputHelper output) => _out = output;

    // 0.5·N(−3, 0.5²) + 0.5·N(+3, 0.5²): two modes separated by a deep barrier.
    private sealed class BimodalTarget : IRjTarget<double>
    {
        public double LogPosterior(double x)
        {
            double l1 = LogNormal(x, -3.0, 0.5);
            double l2 = LogNormal(x, 3.0, 0.5);
            double m = Math.Max(l1, l2);
            return m + Math.Log(0.5 * Math.Exp(l1 - m) + 0.5 * Math.Exp(l2 - m));
        }

        private static double LogNormal(double x, double mu, double sd)
        {
            double z = (x - mu) / sd;
            return -0.5 * z * z - Math.Log(sd * Math.Sqrt(2.0 * Math.PI));
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

    [Fact]
    public void Tempering_VisitsBothModes_WhereSingleChainSticks()
    {
        var target = new BimodalTarget();
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
        var pt = new ParallelTempering<double>(moves, target, -3.0, betas, masterSeed: 7);
        for (int i = 0; i < 5_000; i++) pt.Step();
        int tempRight = 0;
        for (int i = 0; i < draws; i++) if (pt.Step() > 0.0) tempRight++;
        double tempFrac = (double)tempRight / draws;

        _out.WriteLine($"[temper] single-chain right-mode frac={singleFrac:F3} (stuck); " +
                       $"tempering cold frac={tempFrac:F3}; swap acc={(double)pt.SwapAccepts / pt.SwapAttempts:F3}");

        Assert.True(singleFrac < 0.05, $"single chain should stay stuck in one mode (got {singleFrac:F3})");
        Assert.InRange(tempFrac, 0.35, 0.65);   // cold replica visits both modes ≈ equally
    }

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void TemperedEnsemble_RecoversPeak_ViaColdReplica()
    {
        var rng = new Xoshiro256PlusPlus(seed: 11);
        int n = 150;
        const double truePeak = 0.4;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Exp(-120.0 * (x[i] - truePeak) * (x[i] - truePeak)) + 0.05 * Gaussian(rng);
        }

        // Each "chain" is a 5-replica tempering ladder; the readout reads its cold (β = 1) replica unchanged.
        BarsResult result = new BarsEnsemble(
                new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                temperLevels: 5, temperBetaMin: 0.05)
            .Run(x, y, grid: x, chains: 4, masterSeed: 3, burn: 1500, samples: 2000);
        PeakPosterior p = result.Peak;

        _out.WriteLine($"[tempered-ens] peak={p.LocationMean:F3} [{p.LocationLo:F3},{p.LocationHi:F3}] " +
                       $"R̂={p.LocationRHat:F3} ESS={p.LocationEss:F0} coldAcc={result.AcceptanceRate:F3}");

        Assert.True(Math.Abs(p.LocationMean - truePeak) < 0.06, $"tempered peak {p.LocationMean:F3} off {truePeak}");
        Assert.True(p.LocationLo <= truePeak && truePeak <= p.LocationHi, "95% interval should bracket the peak");
        Assert.InRange(result.AcceptanceRate, 1e-6, 1.0);   // cold replica genuinely mixed
    }

    [Fact]
    public void Tempering_IsRejected_AlongsideResampler()
    {
        Assert.Throws<ArgumentException>(() => new BarsEnsemble(
            new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
            resampler: new StudentTWeights(4.0), temperLevels: 5));
    }
}
