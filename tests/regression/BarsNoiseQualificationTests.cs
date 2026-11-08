using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BARS under a variance-mean noise (OQ20, E1's slice 6), its long gates, behind <c>BASS_QUALIFICATION</c> = 1; without it each returns at
/// once.
/// <list type="bullet">
/// <item><b>Geweke's joint-distribution test</b> on a small design: marginal-conditional draws from the knot prior on its support, the
/// coefficients' g-prior, σ²'s inverse gamma, v's law and β's prior, against the successive simulator stepping one BARS chain (a knot move,
/// then the noise block) on each regenerated response.</item>
/// <item><b>Recovery</b> of known quantile curves under heteroscedastic skewed noise, and of a known skew.</item>
/// <item><b>The cost row:</b> a step's cost under each noise beside plain and robust BARS.</item>
/// </list>
/// </summary>
public sealed class BarsNoiseQualificationTests(ITestOutputHelper output)
{
    private const int GewekeN = 16;
    private const double Alpha0 = 4.0, Beta0 = 0.6;
    private static readonly double[] GewekeX = [.. Enumerable.Range(0, GewekeN).Select(i => (i + 0.5) / GewekeN)];
    private static readonly string[] Names = ["k", "ln sigma2", "ln v0", "mean ln v", "atan f0", "atan f8", "beta", "z0^2"];

    private static BarsNoise GewekeNoise(string form) => form switch
    {
        "t" => BarsNoise.StudentT(5.0),
        "quantile" => BarsNoise.Quantile(0.8),
        "skewed-gig" => BarsNoise.Skewed(new GigMixing(new GigParameters(0.7, 1.2, 0.8)), Normal.MeanSd(0.3, 0.5)),
        _ => BarsNoise.Skewed(new GbpMixing(new GbpParameters(1.0, 0.5, 1.5)), Normal.MeanSd(0.0, 0.5)),
    };

    private static WeightedNormalModel GewekeModel() => new(GPriorCentre.Zero, sigma2Prior: Sigma2Prior.InverseGamma(Alpha0, Beta0));

    [Theory]
    [InlineData("t", 20261108, 1)]
    [InlineData("quantile", 20261108, 1)]
    [InlineData("skewed-gig", 20261108, 1)]
    [InlineData("skewed-gbp", 20261108, 1)]
    [InlineData("skewed-gbp", 20261109, 1)]
    [InlineData("skewed-gbp", 20261108, 4)]
    public void Geweke_TheNoiseSweepLeavesTheJointInvariant(string form, int seed, int scale)
    {
        if (!BassQualification.Enabled) return;   // heavy
        int draws = 400_000 * scale, steps = 600_000 * scale;
        const int batches = 50;
        BarsNoise noise = GewekeNoise(form);
        double c = noise.VarianceScale;
        var basis = new SplineBasis(3);
        var prior = new PoissonPrior(1.0);
        WeightedNormalModel model = GewekeModel();
        IReferenceSampler<KnotConfig> knotPrior = new SplineTarget(basis, model, prior, GewekeX, GewekeX).Reference(TemperingPath.Likelihood);
        var sigmaPrior = InverseGamma.ShapeScale(Alpha0, Beta0);
        var rng = new Xoshiro256PlusPlus(seed);

        // θ from the prior: knots on the support, σ², the coefficients N(0, mσ²G⁻¹), v, β.
        (KnotConfig K, double Sigma2, double[] V, double[] F, double Beta) Draw()
        {
            KnotConfig k = knotPrior.Draw(rng).State;
            double sigma2 = sigmaPrior.Sample(rng);
            double[,] z = basis.Design(k, GewekeX);
            double[] f = PriorCurve(z, sigma2, rng);
            double[] v = [.. Enumerable.Range(0, GewekeN).Select(_ => noise.Mixing.SamplePrior(rng))];
            double beta = noise.SkewPrior is { } sp ? sp.Mean + sp.Sd * Normal.SampleStandard(rng) : noise.Skew;
            return (k, sigma2, v, f, beta);
        }
        double[] Respond(double[] f, double sigma2, double[] v, double beta)
            => [.. Enumerable.Range(0, GewekeN).Select(i => f[i] + Math.Sqrt(sigma2) * (beta * v[i] + Math.Sqrt(c * v[i]) * Normal.SampleStandard(rng)))];
        double[] Functions(KnotConfig k, double sigma2, double[] v, double[] f, double[] y, double beta)
            => [k.Count, Math.Log(sigma2), Math.Log(v[0]), v.Average(Math.Log), Math.Atan(f[0]), Math.Atan(f[GewekeN / 2]), beta,
                Math.Pow(y[0] - f[0] - Math.Sqrt(sigma2) * beta * v[0], 2) / (sigma2 * c * v[0])];

        var marginal = new double[draws][];
        for (int d = 0; d < draws; d++)
        {
            var t = Draw();
            marginal[d] = Functions(t.K, t.Sigma2, t.V, t.F, Respond(t.F, t.Sigma2, t.V, t.Beta), t.Beta);
        }

        var theta = Draw();
        var ensemble = new BarsEnsemble(basis, model, prior, new LocalBetaKernel(50.0), noise: noise);
        INoiseProbe probe = ensemble.Probe(GewekeX, Respond(theta.F, theta.Sigma2, theta.V, theta.Beta), theta.K, theta.V, theta.Sigma2, theta.Beta,
                                           rng.NextInt(int.MaxValue));
        var successive = new double[steps][];
        for (int t = 0; t < steps; t++)
        {
            double[] y = Respond(theta.F, theta.Sigma2, theta.V, theta.Beta);
            probe.Respond(y);
            probe.Step();
            double[] v = (double[])probe.V.Clone();
            theta = (probe.Knots, probe.Sigma2, v, probe.Fitted, probe.Skew);
            successive[t] = Functions(theta.K, theta.Sigma2, v, theta.F, y, theta.Beta);
        }

        double worst = 0.0;
        for (int g = 0; g < Names.Length; g++)
        {
            double[] a = [.. marginal.Select(r => r[g])], s = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1), meanB = s.Average();
            int size = s.Length / batches;
            double[] means = [.. Enumerable.Range(0, batches).Select(k => s.Skip(k * size).Take(size).Average())];
            double seB2 = means.Sum(m => (m - means.Average()) * (m - means.Average())) / (batches - 1) / batches;
            bool constant = a.Min() == a.Max() && s.Min() == s.Max();
            double z = constant ? (a[0] == s[0] ? 0.0 : double.PositiveInfinity) : (meanA - meanB) / Math.Sqrt(varA / a.Length + seB2);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{Names[g],-10} marginal-conditional {meanA,10:F5}  successive {meanB,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"{form}, seed {seed}, {scale}× the size: noise blocks refused {probe.BlockRefusals} of {steps}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // f = Zβ with β ~ N(0, mσ²G⁻¹), G = ZᵀZ: β = √(mσ²)·L⁻ᵀξ for G = LLᵀ.
    private static double[] PriorCurve(double[,] z, double sigma2, Xoshiro256PlusPlus rng)
    {
        int m = z.GetLength(0), p = z.GetLength(1);
        var g = new double[p, p];
        for (int a = 0; a < p; a++)
            for (int b = 0; b < p; b++)
                for (int i = 0; i < m; i++) g[a, b] += z[i, a] * z[i, b];
        var l = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = g[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = i == j ? Math.Sqrt(s) : s / l[j, j];
            }
        var beta = new double[p];
        double scale = Math.Sqrt(m * sigma2);
        for (int i = p - 1; i >= 0; i--)   // Lᵀβ = ξ
        {
            double s = Normal.SampleStandard(rng);
            for (int k = i + 1; k < p; k++) s -= l[k, i] * beta[k];
            beta[i] = s / l[i, i];
        }
        var f = new double[m];
        for (int i = 0; i < m; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < p; j++) sum += z[i, j] * beta[j];
            f[i] = scale * sum;
        }
        return f;
    }

    [Fact]
    public void Recovery_QuantileCurves_AndAKnownSkew()
    {
        if (!BassQualification.Enabled) return;   // heavy
        // Heteroscedastic skewed noise: s(x)·(E − 1), E ~ Exp(1), whose q-quantile curve is sin 2πx + s(x)(−ln(1 − q) − 1).
        const int n = 400;
        double Scale(double x) => 0.2 + 0.3 * x;
        var rng = new Xoshiro256PlusPlus(20261108);
        double[] x = [.. Enumerable.Range(0, n).Select(i => (i + 0.5) / n)];
        double[] grid = [.. Enumerable.Range(0, 101).Select(i => i / 100.0)];
        var y = new double[n];
        for (int i = 0; i < n; i++) y[i] = Math.Sin(2.0 * Math.PI * x[i]) + Scale(x[i]) * (-Math.Log(1.0 - rng.NextDouble()) - 1.0);
        foreach (double q in new[] { 0.1, 0.5, 0.9 })
        {
            var ensemble = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                                            noise: BarsNoise.Quantile(q));
            BarsResult result = ensemble.Run(x, y, grid, chains: 4, masterSeed: 5, burn: 4000, samples: 2000);
            double squares = 0.0;
            for (int g = 0; g < grid.Length; g++)
            {
                double truth = Math.Sin(2.0 * Math.PI * grid[g]) + Scale(grid[g]) * (-Math.Log(1.0 - q) - 1.0);
                squares += (result.Fit[g] - truth) * (result.Fit[g] - truth);
            }
            BarsResult atData = ensemble.Run(x, y, x, chains: 4, masterSeed: 5, burn: 4000, samples: 2000);
            double share = Enumerable.Range(0, n).Count(i => y[i] <= atData.Fit[i]) / (double)n;
            output.WriteLine($"AL({q}): RMSE against the true {q}-quantile curve {Math.Sqrt(squares / grid.Length):F4} over the grid; " +
                             $"{share:P1} of the data at or below the fit; mean knots {result.MeanKnots:F2}");
            Assert.InRange(share, q - 0.05, q + 0.05);
        }

        // A known skew over four data sets: y = sin 2πx + σ(βv + √v·ζ), v ~ GIG(0.7, 1.2, 0.8), β = 0.8, σ = 0.25. The noise's mean σβE[v]
        // and the curve's level trade off, so β is read from the residuals' asymmetry alone.
        var mixing = new GigMixing(new GigParameters(0.7, 1.2, 0.8));
        int covered = 0;
        double meanOfMeans = 0.0;
        for (int set = 0; set < 4; set++)
        {
            var skewed = new double[n];
            for (int i = 0; i < n; i++)
            {
                double v = mixing.SamplePrior(rng);
                skewed[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.25 * (0.8 * v + Math.Sqrt(v) * Normal.SampleStandard(rng));
            }
            var fit = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                                       noise: BarsNoise.Skewed(mixing, Normal.MeanSd(0.0, 2.0)));
            BarsResult r = fit.Run(x, skewed, grid, chains: 4, masterSeed: 7 + set, burn: 4000, samples: 2000);
            double[] betas = [.. r.Draws.SelectMany(d => d).Select(d => d.Skew).Order()];
            double[] sigmas = [.. r.Draws.SelectMany(d => d).Select(d => Math.Sqrt(d.Sigma2)).Order()];
            double lo = betas[(int)(0.05 * betas.Length)], hi = betas[(int)(0.95 * betas.Length) - 1];
            covered += lo <= 0.8 && 0.8 <= hi ? 1 : 0;
            meanOfMeans += betas.Average() / 4.0;
            output.WriteLine($"skewed GIG, set {set}: β mean {betas.Average():F3} [{lo:F3}, {hi:F3}] against 0.8; σ mean {sigmas.Average():F4} against 0.25");
        }
        output.WriteLine($"skewed GIG: 90% intervals cover β in {covered} of 4 sets; the posterior means average {meanOfMeans:F3}");
        Assert.True(covered >= 3, $"β covered in {covered} of 4 sets");
    }

    [Fact]
    public void CostRow_AStep_UnderEachNoise()
    {
        if (!BassQualification.Enabled) return;   // heavy
        // One chain stepped on this thread (BarsEnsemble.ChainOn): a run fans its chains out over the pool, which a thread's clock misses
        // and the process's CPU time inflates with the pool's spinning.
        string clock = BassCost.Cycles ? "thread cycles, indicative" : "wall clock";
        output.WriteLine($"Clock: {clock}. µs per BARS step (a knot move, the noise or weight block every 25), one chain on this thread, " +
                         "after 1,000 untimed steps.");
        foreach (int n in new[] { 200, 1000 })
        {
            var rng = new Xoshiro256PlusPlus(9);
            double[] x = [.. Enumerable.Range(0, n).Select(i => (i + 0.5) / n)];
            double[] y = [.. x.Select(t => Math.Sin(2.0 * Math.PI * t) + 0.3 * (-Math.Log(1.0 - rng.NextDouble()) - 1.0))];
            var cells = new List<string>();
            foreach ((string name, Func<BarsEnsemble> make) in new (string, Func<BarsEnsemble>)[]
            {
                ("plain", () => Ensemble(null, null)),
                ("StudentTWeights(4)", () => Ensemble(new StudentTWeights(4.0), null)),
                ("t(4)", () => Ensemble(null, BarsNoise.StudentT(4.0))),
                ("AL(0.9)", () => Ensemble(null, BarsNoise.Quantile(0.9))),
                ("skewed GIG", () => Ensemble(null, BarsNoise.Skewed(new GigMixing(new GigParameters(0.7, 1.2, 0.8)), Normal.MeanSd(0.0, 1.0)))),
            })
            {
                const int warm = 1000, timed = 2000;
                IChain<KnotConfig> chain = make().ChainOn(x, y, weightEvery: 25, seed: 1);
                for (int s = 0; s < warm; s++) chain.Step();
                long t0 = BassCost.Tick();
                KnotConfig last = chain.Step();
                for (int s = 1; s < timed; s++) last = chain.Step();
                double us = (BassCost.Tick() - t0) / (double)timed / BassCost.TicksPerSecond * 1e6;
                cells.Add($"{name} {us:F1} (k {last.Count})");
            }
            output.WriteLine($"n {n,5}: {string.Join("; ", cells)}");
        }
    }

    private static BarsEnsemble Ensemble(IWeightResampler? resampler, BarsNoise? noise)
        => new(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0), resampler: resampler, noise: noise);
}
