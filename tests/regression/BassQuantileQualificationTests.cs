using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The asymmetric Laplace's long gates (E1's slice 2), behind <c>BASS_QUALIFICATION</c> = 1; without it each returns at once.
/// <list type="bullet">
/// <item><b>Geweke's joint-distribution test</b> under the skew: σ² on the record, its modified half-normal block, v's GIG block with the
/// drift, the offset following both; and at q = ½, where the drift vanishes and σ² is integrated.</item>
/// <item><b>The skew's sign</b> (E1's fixture): at q = 0.9 and 0.1 on symmetric noise, the held-out share below each surface is its q, which a
/// flipped sign would swap.</item>
/// <item><b>The q-quantile recovered</b> on data whose noise is skewed and heteroscedastic, outside the model, with the true quantile known:
/// the surface's error, its 90% intervals' coverage of the true quantile, and the held-out share below it.</item>
/// <item><b>The cost row and retention</b>, before and after the draw keeps its latents alone.</item>
/// </list>
/// </summary>
public sealed class BassQuantileQualificationTests(ITestOutputHelper output)
{
    private static readonly BassRunOptions Capped = new() { MaxDegreeOfParallelism = 4 };

    // ---- Geweke --------------------------------------------------------------------------------------------------------------------

    private const int GewekeN = 12;
    private static readonly double[] GewekeX = [.. Enumerable.Range(0, GewekeN).Select(i => i / (GewekeN - 1.0))];
    private static readonly string[] GewekeNames = ["M", "ln tau", "ln sigma2", "ln v0", "mean ln v", "atan f0", "atan f6", "z0^2"];

    private static BassData GewekeData(double[] y) => new([new ContinuousPredictor("x", GewekeX)], new ScalarResponse(y));

    private static BassConfig GewekeConfig(double q, double tauShape) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 2, Support = new SupportRule(SupportKind.Count, 2) },
        Coefficients = BassCoefficientPrior.G(tauShape > 0 ? TauPrior.UnitInformationGamma(tauShape) : TauPrior.Held()),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.InverseGamma(3, 2),
        Noise = BassNoise.AsymmetricLaplace(q, BassNoiseReference.None),
        Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 0, Samples = 1 },
    };

    // The working response at (v, σ²), through the model's own reading; σ² rides the record under the skew.
    private static WorkingResponse At(BassModel model, double[] v, double sigma2)
    {
        var mixture = (VarianceMeanMixture)model.Observation;
        return mixture.At(model.Working, new ObservationLatents(v, mixture.Skew, mixture.Skew != 0.0 ? sigma2 : null));
    }

    [Theory]
    [InlineData(0.3, 5.0, 20261105)]
    [InlineData(0.3, 5.0, 20261106)]
    [InlineData(0.8, 0.0, 20261105)]
    [InlineData(0.5, 5.0, 20261105)]
    public void Geweke_TheSweepLeavesTheJointInvariant(double q, double tauShape, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int draws = 200_000, steps = 300_000, batches = 50;
        BassConfig config = GewekeConfig(q, tauShape);
        double c = config.Noise!.VarianceScale, beta = config.Noise.Skew;
        IMixingLaw law = config.Noise.Mixing;
        InverseGamma sigmaPrior = InverseGamma.ShapeScale(3, 2);
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(config, GewekeData([.. GewekeX.Select(x => Math.Sin(7 * x))]));
        var sampler = new BassPriorSampler(prior);

        double[] Functions(BassState s, double sigma2, double[] v, double[] f, double[] y)
            => [s.Count, s.LogTau, Math.Log(sigma2), Math.Log(v[0]), v.Average(x => Math.Log(x)), Math.Atan(f[0]), Math.Atan(f[GewekeN / 2]),
                Math.Pow(y[0] - f[0] - Math.Sqrt(sigma2) * beta * v[0], 2) / (sigma2 * c * v[0])];
        double[] Respond(double[] f, double sigma2, double[] v)
            => [.. Enumerable.Range(0, GewekeN).Select(i => f[i] + Math.Sqrt(sigma2) * (beta * v[i] + Math.Sqrt(c * v[i]) * Normal.SampleStandard(rng)))];

        // The marginal-conditional simulator. With the drift, σ² is drawn from its prior and held; at q = ½ the draw of a takes it with it.
        var marginal = new double[draws][];
        long redrawn = 0;
        (BassState S, double Sigma2, double[] V, double[] F) Draw()
        {
            while (true)
            {
                BassState state = sampler.DrawState(rng).State;
                double[] v = [.. Enumerable.Range(0, GewekeN).Select(_ => law.SamplePrior(rng))];
                double held = sigmaPrior.Sample(rng);
                var chain = new BassChain(prior, state.Terms, state.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
                if (chain.Follow(At(prior, v, held)) is not null) { redrawn++; continue; }
                BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, prior.Target, power: 0.0);
                return (chain.Current, fit.Sigma2, v, fit.Fitted);
            }
        }
        for (int k = 0; k < draws; k++)
        {
            (BassState s, double sigma2, double[] v, double[] f) = Draw();
            marginal[k] = Functions(s, sigma2, v, f, Respond(f, sigma2, v));
        }

        // The successive-conditional simulator: y given θ, then one sweep on that y. σ² after the sweep is the record's under the drift,
        // the coefficients' block's at q = ½.
        (BassState S, double Sigma2, double[] V, double[] F) theta = Draw();
        var successive = new double[steps][];
        long refused = 0, blockRefused = 0;
        for (int t = 0; t < steps; t++)
        {
            double[] y = Respond(theta.F, theta.Sigma2, theta.V);
            BassModel model = BassModel.Create(config, GewekeData(y));
            var chain = new BassChain(model, theta.S.Terms, theta.S.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            if (chain.Follow(At(model, theta.V, theta.Sigma2)) is not null) { refused++; continue; }
            chain.Step();
            foreach (MoveStat s in chain.Scan.BlockStats().Where(b => b.Key is "latent scales" or "scale")) blockRefused += s.Attempts - s.Accepted;
            BassSweepFit fit = chain.Structure.SweepFit!;
            ObservationLatents latents = chain.Current.Fit.Working!.Latents!;
            double sigma2 = latents.Sigma2 ?? fit.Sigma2;
            successive[t] = Functions(chain.Current, sigma2, latents.ValueArray, fit.Fitted, y);
            theta = (chain.Current, sigma2, latents.ValueArray, fit.Fitted);
        }
        Assert.Equal(0, refused);

        double worst = 0;
        for (int g = 0; g < GewekeNames.Length; g++)
        {
            double[] a = [.. marginal.Select(r => r[g])], b = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1), meanB = b.Average();
            int size = b.Length / batches;
            double[] means = [.. Enumerable.Range(0, batches).Select(k => b.Skip(k * size).Take(size).Average())];
            double seB2 = means.Sum(m => (m - means.Average()) * (m - means.Average())) / (batches - 1) / batches;
            bool constant = a.Min() == a.Max() && b.Min() == b.Max();
            double z = constant ? (a[0] == b[0] ? 0.0 : double.PositiveInfinity) : (meanA - meanB) / Math.Sqrt(varA / a.Length + seB2);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{GewekeNames[g],-10} marginal-conditional {meanA,10:F5}  successive {meanB,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"q {q}, tau shape {tauShape} (0 held), seed {seed}: redrawn {redrawn}; moves refused {refused}; block draws refused {blockRefused} of {steps}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- the sign and the quantile -------------------------------------------------------------------------------------------------

    private static double Truth(double x0, double x1) => Math.Sin(2 * Math.PI * x0) + x1;

    // f plus s(x)·(G − 2), G ~ Gamma(2, 1), s(x) = 0.3 + 0.4·x₁: skewed and heteroscedastic, so outside the model, with the q-quantile
    // f + s(x)·(Q_G(q) − 2) known. `symmetric` takes 0.3·N(0, 1) instead.
    private static BassData Noisy(int n, int seed, bool symmetric = false)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n];
        var gamma = Gamma.ShapeRate(2, 1);
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            y[i] = Truth(x0[i], x1[i]) + (symmetric ? 0.3 * Normal.SampleStandard(rng) : (0.3 + 0.4 * x1[i]) * (gamma.Sample(rng) - 2));
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
    }

    private static BassConfig QuantileConfig(double q, int seed) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 20 },
        Noise = BassNoise.AsymmetricLaplace(q),
        Sampling = new BassSampling { Chains = 4, MasterSeed = seed, BurnIn = 2000, Samples = 500, Thin = 2 },
    };

    // The share of `test`'s responses below the fit's posterior mean surface at its inputs.
    private static double ShareBelow(BassFit fit, BassData data, BassData test)
    {
        double[,] mean = fit.Predict(data, PredictionCases.Of(test), new BassPredictOptions { Predictive = false }).Mean();
        ReadOnlySpan<double> y = test.Response.Values;
        int below = 0;
        for (int i = 0; i < y.Length; i++) below += y[i] < mean[i, 0] ? 1 : 0;
        return (double)below / y.Length;
    }

    [Fact]
    public void TheSkewsSign_PutsQOfTheResponsesBelowTheSurface()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassData data = Noisy(400, 61, symmetric: true), test = Noisy(2000, 62, symmetric: true);
        foreach (double q in new[] { 0.9, 0.1 })
        {
            BassFit fit = Bass.Fit(QuantileConfig(q, 71), data, Capped);
            double held = ShareBelow(fit, data, test), inSample = ShareBelow(fit, data, data);
            output.WriteLine($"q {q}: held-out share below the surface {held:F3}, in sample {inSample:F3}");
            Assert.InRange(held, q - 0.05, q + 0.05);
        }
    }

    [Fact]
    public void TheQuantile_IsRecovered_OnNoiseOutsideTheModel()
    {
        if (!BassQualification.Enabled) return;   // heavy
        double[] grid = [.. Enumerable.Range(0, 7).Select(k => 0.05 + 0.15 * k)];
        var points = (from a in grid from b in grid select (a, b)).ToArray();
        var cases = new PredictionCases([new ContinuousColumn("x0", [.. points.Select(p => p.a)]), new ContinuousColumn("x1", [.. points.Select(p => p.b)])]);
        var gamma = Gamma.ShapeRate(2, 1);
        foreach (double q in new[] { 0.5, 0.9 })
        {
            int covered = 0, total = 0;
            double squares = 0, held = 0;
            const int sets = 4;
            for (int set = 0; set < sets; set++)
            {
                BassData data = Noisy(400, 100 + set), test = Noisy(2000, 200 + set);
                BassFit fit = Bass.Fit(QuantileConfig(q, 31 + set), data, Capped);
                BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Predictive = false });
                double[,] mean = p.Mean();
                for (int i = 0; i < points.Length; i++)
                {
                    double truth = Truth(points[i].a, points[i].b) + (0.3 + 0.4 * points[i].b) * (gamma.Quantile(q) - 2);
                    double[] f = [.. Enumerable.Range(0, p.DrawCount).Select(d => p.Latent[d, i, 0]).Order()];
                    covered += truth >= f[(int)(0.05 * f.Length)] && truth <= f[(int)(0.95 * f.Length) - 1] ? 1 : 0;
                    squares += (mean[i, 0] - truth) * (mean[i, 0] - truth);
                    total++;
                }
                held += ShareBelow(fit, data, test) / sets;
            }
            double coverage = (double)covered / total, rmse = Math.Sqrt(squares / total);
            output.WriteLine($"q {q}: the true quantile inside f's 90% intervals at {coverage:P1} of {total} grid points, RMSE {rmse:F4}; " +
                             $"held-out share below the surface {held:F3}");
            // The recovery is the held-out share below the surface: the surface is the q-quantile. The error and the intervals' coverage
            // are readings: the asymmetric Laplace is a working likelihood here, and its posterior intervals are not calibrated for the
            // true quantile under other noise (Yang, Wang and He 2016).
            Assert.InRange(held, q - 0.04, q + 0.04);
        }
    }

    // ---- cost and retention --------------------------------------------------------------------------------------------------------

    [Fact]
    public void CostRow_PerPartOfASweep_AndRetentionBeforeAndAfter()
    {
        if (!BassQualification.Enabled) return;   // heavy
        output.WriteLine($"Clock: {(BassCost.Cycles ? "thread cycles, indicative" : "wall clock")}. µs per sweep, each part's mean; ALD at q = 0.9.");
        output.WriteLine("point                               form   M     structure  tau    a       v draws  Follow v  sigma2  Follow s2  sweep  | Gaussian");
        var points = new (BassCost.Point Point, int Burn, int Timed)[]
        {
            (BassCost.Base, 2000, 400),
            (BassCost.Base with { N = 200 }, 2000, 400),
            (BassCost.Base with { N = 10000, Cap = 200 }, 0, 10),
        };
        foreach ((BassCost.Point x, int burn, int timed) in points)
            foreach (bool ridge in new[] { false, true })
            {
                if (burn == 0 && ridge) continue;   // the largest point once, under the g-form: its sweeps are seconds long
                BassCoefficientPrior form = ridge ? BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()) : BassCoefficientPrior.G();
                BassData data = BassCost.Data(x);
                BassConfig gaussian = BassCost.Config(x) with { Coefficients = form };
                IReadOnlyList<BassTerm> start = [];
                double logTau = 0.0;
                if (burn == 0)
                {
                    (_, BassChain burned, _) = BassCost.Burned(x);
                    start = burned.Current.Terms;
                    logTau = burned.Current.LogTau;
                }
                BassModel model = BassModel.Create(gaussian with { Noise = BassNoise.AsymmetricLaplace(0.9) }, data);
                var chain = new BassChain(model, start, logTau, new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(5));
                for (int s = 0; s < burn; s++) chain.Step();
                var parts = new double[7];
                double count = 0;
                var rng = new Xoshiro256PlusPlus(7);
                List<IBlock<ReversibleJumpChain<BassState>>> blocks = BassSweepSlots.Blocks(model);
                var noise = (VarianceMeanMixture)model.Observation;
                for (int s = 0; s < timed; s++)
                {
                    parts[0] += Timed(() => blocks[0].Update(chain.Moves, rng));
                    parts[1] += Timed(() => blocks[1].Update(chain.Moves, rng));
                    parts[2] += Timed(() => blocks[2].Update(chain.Moves, rng));
                    // The latent block in its parts: the N draws of v, then Follow; then σ²'s block in its: the modified half-normal, then Follow.
                    BassSweepFit fit = chain.Structure.SweepFit!;
                    WorkingResponse working = chain.Structure.DesignAt(chain.Current).Working;
                    ObservationLatents latents = working.Latents!;
                    double[] v = new double[latents.Values.Length];
                    double drift = Math.Sqrt(fit.Sigma2) * latents.Skew;
                    parts[3] += Timed(() =>
                    {
                        ReadOnlySpan<double> y = working.Response.Values;
                        for (int i = 0; i < v.Length; i++)
                            v[i] = noise.Mixing.SampleConditional(y[i] - fit.Fitted[i], fit.Sigma2 * noise.VarianceScale, drift, 1.0, rng);
                    });
                    parts[4] += Timed(() => BassStructureMoves.Follow(chain.Moves, noise.At(working, new ObservationLatents(v, latents.Skew, latents.Sigma2))));
                    // σ²'s block whole, then its draw alone read off by redoing the draw's arithmetic is not separable without the block's
                    // internals, so the scale block is timed whole and its Follow apart from a second, offset-only move to the same σ².
                    WorkingResponse before = chain.Structure.DesignAt(chain.Current).Working;
                    long scale = Timed(() => blocks[4].Update(chain.Moves, rng));
                    WorkingResponse after = chain.Structure.DesignAt(chain.Current).Working;
                    long follow = Timed(() => BassStructureMoves.Follow(chain.Moves, noise.At(after, after.Latents!.WithSigma2(before.Latents!.Sigma2!.Value), weightsMoved: false)));
                    BassStructureMoves.Follow(chain.Moves, after);   // back to the block's draw
                    parts[5] += Math.Max(0, scale - follow);
                    parts[6] += follow;
                    count += chain.Current.Count;
                }
                BassModel plain = BassModel.Create(gaussian, data);
                var reference = new BassChain(plain, chain.Current.Terms, chain.Current.LogTau, new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(5));
                double gaussianSweep = 0;
                int plainTimed = Math.Max(timed, 50);
                for (int s = 0; s < plainTimed; s++) gaussianSweep += Timed(() => reference.Step());
                double us(double ticks) => ticks / timed / BassCost.TicksPerSecond * 1e6;
                output.WriteLine($"{x.Key,-36}{(ridge ? "ridge" : "g"),-7}{count / timed,-6:F1}{us(parts[0]),-11:F1}{us(parts[1]),-7:F1}{us(parts[2]),-8:F1}" +
                                 $"{us(parts[3]),-9:F1}{us(parts[4]),-10:F1}{us(parts[5]),-8:F1}{us(parts[6]),-11:F1}{us(parts.Sum()),-7:F1}| {gaussianSweep / plainTimed / BassCost.TicksPerSecond * 1e6:F1}");
            }

    }

    /// <summary>
    /// Retention, before and after the draw keeps its latents alone: K retained draws of one chain held two ways, the working response
    /// each was composed on (the former draw) and the draw as it is now. Under the skew the working response holds the weights, the
    /// offset, the shifted response and the latents.
    /// </summary>
    [Fact]
    public void Retention_BeforeAndAfter()
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int k = 200;
        foreach (int n in new[] { 1_000, 10_000 })
        {
            long none = Settled();
            (List<WorkingResponse> workings, List<BassDraw> draws) = Retained(n, k);
            long both = Settled();
            int distinct = draws.Select(d => d.Latents).Distinct().Count();
            workings.Clear();
            long now = Settled();
            GC.KeepAlive(draws);
            output.WriteLine($"N {n:N0}: {distinct} distinct latents over {k} draws; a retained draw holds {(both - none) / (double)k / 1024:F1} KB with its " +
                             $"working response, {(now - none) / (double)k / 1024:F1} KB with its latents alone ({(both - none) / (double)Math.Max(1, now - none):F2}×)");
        }
    }

    // K draws of a quantile chain, each with the working response it was composed on; the chain itself left behind.
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.NoInlining)]
    private static (List<WorkingResponse> Workings, List<BassDraw> Draws) Retained(int n, int k)
    {
        BassData data = Noisy(n, 5);
        BassModel model = BassModel.Create(QuantileConfig(0.9, 3) with { Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 8 } }, data);
        var chain = new BassChain(model, [], 0.0, new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(5));
        var terms = new BassTermTable();
        var compose = new Xoshiro256PlusPlus(9);
        var workings = new List<WorkingResponse>();
        var draws = new List<BassDraw>();
        for (int s = 0; s < k; s++)
        {
            chain.Step();
            workings.Add(chain.Current.Fit.Working!);
            draws.Add(chain.Compose(compose, terms));
        }
        return (workings, draws);
    }

    private static long Settled()
    {
        GC.Collect();
        GC.WaitForPendingFinalizers();
        GC.Collect();
        return GC.GetTotalMemory(true);
    }

    private static long Timed(Action action)
    {
        long t0 = BassCost.Tick();
        action();
        return BassCost.Tick() - t0;
    }

    private static long Timed(Func<bool> action) => Timed(() => { action(); });

    private static long Timed(Func<BassRefusal?> action) => Timed(() => { action(); });

    private static long Timed(Func<BassState> action) => Timed(() => { action(); });
}
