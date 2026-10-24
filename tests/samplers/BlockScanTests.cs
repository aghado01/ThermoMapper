using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// The shared block primitives and the block scan. The Gaussian conditional draw has the mean and covariance of an
/// explicitly inverted precision, through a banded and a dense factor; the conjugate updates match the posterior by
/// numerical integration; a two-block Gibbs scan samples a known bivariate normal; cadences and counters behave; and a
/// restored scan continues exactly, leaving its snapshot usable.
/// </summary>
public sealed class BlockScanTests
{
    private static readonly double[,] Q = { { 2.0, -0.6, 0.0 }, { -0.6, 2.5, 0.7 }, { 0.0, 0.7, 1.8 } };
    private static readonly double[] Rhs = { 1.0, -2.0, 0.5 };

    // The inverse of a 3 × 3 matrix by cofactors, independent of either factor.
    private static double[,] Inverse(double[,] a)
    {
        double c00 = a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1], c01 = a[1, 2] * a[2, 0] - a[1, 0] * a[2, 2];
        double c02 = a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0];
        double det = a[0, 0] * c00 + a[0, 1] * c01 + a[0, 2] * c02;
        var inv = new double[3, 3];
        inv[0, 0] = c00 / det; inv[1, 0] = c01 / det; inv[2, 0] = c02 / det;
        inv[0, 1] = (a[0, 2] * a[2, 1] - a[0, 1] * a[2, 2]) / det;
        inv[1, 1] = (a[0, 0] * a[2, 2] - a[0, 2] * a[2, 0]) / det;
        inv[2, 1] = (a[0, 1] * a[2, 0] - a[0, 0] * a[2, 1]) / det;
        inv[0, 2] = (a[0, 1] * a[1, 2] - a[0, 2] * a[1, 1]) / det;
        inv[1, 2] = (a[0, 2] * a[1, 0] - a[0, 0] * a[1, 2]) / det;
        inv[2, 2] = (a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]) / det;
        return inv;
    }

    public static TheoryData<string> Factors => new() { "banded", "dense" };

    [Theory]
    [MemberData(nameof(Factors))]
    public void GaussianConditional_HasTheConditionalsMeanAndCovariance(string kind)
    {
        ISpdFactor factor;
        if (kind == "banded")
        {
            var band = new BandCholesky(3, 1, BandFactorization.Ldlt);   // the form that samples, as the engines use
            band.Decompose(Q);
            factor = band;
        }
        else
        {
            var dense = new DenseCholesky(3);
            dense.Decompose(Q);
            factor = dense;
        }
        double[,] inv = Inverse(Q);
        var exact = new double[3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) exact[i] += inv[i, j] * Rhs[j];

        const double scale = 1.7;
        const int n = 200_000;
        var gauss = new GaussianConditional(3);
        var rng = new Xoshiro256PlusPlus(seed: 41);
        var draw = new double[3];
        var sum = new double[3];
        var cross = new double[3, 3];
        double[] mean = Array.Empty<double>();
        for (int s = 0; s < n; s++)
        {
            mean = gauss.Draw(factor, Rhs, rng, draw, scale);
            for (int i = 0; i < 3; i++)
            {
                sum[i] += draw[i];
                for (int j = 0; j < 3; j++) cross[i, j] += (draw[i] - exact[i]) * (draw[j] - exact[j]);
            }
        }
        for (int i = 0; i < 3; i++)
        {
            Assert.Equal(exact[i], mean[i], 12);
            double sd = scale * Math.Sqrt(inv[i, i]);
            Assert.True(Math.Abs(sum[i] / n - exact[i]) <= 5.0 * sd / Math.Sqrt(n), $"mean {i}");
            for (int j = 0; j < 3; j++)
            {
                double cov = scale * scale * inv[i, j];
                double se = scale * scale * Math.Sqrt((inv[i, i] * inv[j, j] + inv[i, j] * inv[i, j]) / n);
                Assert.True(Math.Abs(cross[i, j] / n - cov) <= 5.0 * se, $"cov {i},{j}: {cross[i, j] / n} vs {cov}");
            }
        }
    }

    [Fact]
    public void ConjugateUpdates_MatchThePosteriorByNumericalIntegration()
    {
        double[] x = { 0.3, -1.2, 0.8, 2.1, -0.4, 0.05, -1.7 };
        double ss = 0.0;
        foreach (double v in x) ss += v * v;
        const double a = 2.0, b = 1.5;

        // E[σ² | x] and E[1/σ² | x] by the trapezoid rule over u = log σ² on [−12, 12], against prior × likelihood. The
        // prior IG(a, b) on σ² is the prior Gamma(a, b) on 1/σ²; either way its density in u is ∝ exp(−a·u − b/σ²).
        double w0 = 0.0, wVariance = 0.0, wPrecision = 0.0;
        const int steps = 200_000;
        double h = 24.0 / steps;
        for (int k = 0; k <= steps; k++)
        {
            double u = -12.0 + k * h, sigma2 = Math.Exp(u);
            double logPrior = -a * u - b / sigma2;
            double logLikelihood = -0.5 * x.Length * u - 0.5 * ss / sigma2;
            double w = (k == 0 || k == steps ? 0.5 : 1.0) * Math.Exp(logPrior + logLikelihood);
            w0 += w;
            wVariance += w * sigma2;
            wPrecision += w / sigma2;
        }
        Assert.Equal(wVariance / w0, ConjugateUpdates.NormalVariance(a, b, x.Length, ss).Mean, 8);
        Assert.Equal(wPrecision / w0, ConjugateUpdates.NormalPrecision(a, b, x.Length, ss).Mean, 8);
    }

    private sealed class Pair
    {
        public double X, Y;
        public Pair Copy() => new() { X = X, Y = Y };
    }

    // x | y ~ N(ρy, 1 − ρ²), and the same for y | x: a two-block Gibbs scan on the standard bivariate normal.
    private sealed class Conditional(bool first, double rho) : IBlock<Pair>
    {
        public string Name => first ? "x" : "y";

        public bool Update(Pair state, Xoshiro256PlusPlus rng)
        {
            double sd = Math.Sqrt(1.0 - rho * rho);
            if (first) state.X = rho * state.Y + sd * Normal.SampleStandard(rng);
            else state.Y = rho * state.X + sd * Normal.SampleStandard(rng);
            return true;
        }
    }

    private static IBlock<Pair>[] Gibbs(double rho) => new IBlock<Pair>[] { new Conditional(true, rho), new Conditional(false, rho) };

    [Fact]
    public void TwoBlockGibbs_SamplesTheBivariateNormal()
    {
        const double rho = 0.8;
        var scan = new BlockScan<Pair>(Gibbs(rho), new Pair(), new Xoshiro256PlusPlus(seed: 3), p => p.Copy());
        const int n = 200_000;
        double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
        for (int i = 0; i < n; i++)
        {
            Pair p = scan.Step();
            sx += p.X; sy += p.Y; sxx += p.X * p.X; syy += p.Y * p.Y; sxy += p.X * p.Y;
        }
        // Each coordinate's lag-one autocorrelation is ρ² = 0.64, so about 22% of the draws are effective.
        Assert.InRange(sx / n, -0.03, 0.03);
        Assert.InRange(sy / n, -0.03, 0.03);
        Assert.InRange(sxx / n, 0.96, 1.04);
        Assert.InRange(syy / n, 0.96, 1.04);
        Assert.InRange(sxy / Math.Sqrt(sxx * syy), rho - 0.02, rho + 0.02);
        Assert.Equal(2L * n, scan.Attempts);
        Assert.Equal(2L * n, scan.Accepted);
    }

    private sealed class Counter(string name, int every, bool acceptEveryOther) : IBlock<Pair>
    {
        private int _calls;
        public string Name => name;
        public int Every => every;

        public bool Update(Pair state, Xoshiro256PlusPlus rng)
        {
            _calls++;
            state.X += 1.0;
            return !acceptEveryOther || _calls % 2 == 0;
        }
    }

    [Fact]
    public void Cadences_AndCounters_FollowTheSchedule()
    {
        var blocks = new IBlock<Pair>[] { new Counter("each", 1, false), new Counter("third", 3, true) };
        var scan = new BlockScan<Pair>(blocks, new Pair(), new Xoshiro256PlusPlus(seed: 1), p => p.Copy());
        for (int i = 0; i < 9; i++) scan.Step();

        IReadOnlyList<Rjmcmc.MoveStat> stats = scan.BlockStats();
        Assert.Equal(new Rjmcmc.MoveStat("each", 9, 9), stats[0]);
        Assert.Equal(new Rjmcmc.MoveStat("third", 3, 1), stats[1]);   // runs on steps 3, 6, 9; accepts its second call
        Assert.Equal(12, scan.Attempts);
        Assert.Equal(10, scan.Accepted);
        Assert.Equal(12.0, scan.State.X);
        Assert.Equal(9, scan.Steps);
    }

    [Fact]
    public void RestoredScan_ContinuesExactly_AndLeavesItsSnapshotUsable()
    {
        var live = new BlockScan<Pair>(Gibbs(0.6), new Pair(), new Xoshiro256PlusPlus(seed: 12), p => p.Copy());
        for (int i = 0; i < 50; i++) live.Step();
        BlockScanSnapshot<Pair> snapshot = live.Snapshot();

        var expected = new List<(double, double)>();
        for (int i = 0; i < 50; i++) { Pair p = live.Step(); expected.Add((p.X, p.Y)); }

        for (int attempt = 0; attempt < 2; attempt++)
        {
            BlockScan<Pair> resumed = BlockScan<Pair>.Restore(Gibbs(0.6), snapshot, p => p.Copy());
            Assert.Equal(50, resumed.Steps);
            Assert.Equal(100, resumed.Attempts);
            for (int i = 0; i < 50; i++)
            {
                Pair p = resumed.Step();
                Assert.Equal(expected[i], (p.X, p.Y));
            }
            Assert.Equal(live.BlockStats(), resumed.BlockStats());
        }
    }

    [Fact]
    public void InvalidScans_AreRefused()
    {
        Assert.Throws<ArgumentException>(() => new BlockScan<Pair>(Array.Empty<IBlock<Pair>>(), new Pair(), new Xoshiro256PlusPlus(1), p => p.Copy()));
        Assert.Throws<ArgumentException>(() => new BlockScan<Pair>(new IBlock<Pair>[] { new Counter("never", 0, false) }, new Pair(),
                                                                    new Xoshiro256PlusPlus(1), p => p.Copy()));
        BlockScanSnapshot<Pair> snapshot = new BlockScan<Pair>(Gibbs(0.5), new Pair(), new Xoshiro256PlusPlus(1), p => p.Copy()).Snapshot();
        Assert.Throws<ArgumentException>(() => BlockScan<Pair>.Restore(new IBlock<Pair>[] { new Conditional(true, 0.5) }, snapshot, p => p.Copy()));
    }
}
