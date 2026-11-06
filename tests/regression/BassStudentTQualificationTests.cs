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
/// Student t's long gates (E1's slice 1), behind <c>BASS_QUALIFICATION</c> = 1; without it each returns at once.
/// <list type="bullet">
/// <item><b>Geweke's joint-distribution test</b> on a design small enough to run long, with proper priors: the marginal-conditional
/// simulator (θ from the prior, y given θ) against the successive-conditional one (y given θ, then one sweep of the chain on that y),
/// eight functions of (θ, y) compared by their means.</item>
/// <item><b>Large ν reaches Gaussian BASS's readouts:</b> ν = 10⁷ against Gaussian BASS on one data set, sixteen chains each, the fit
/// at every training case, σ² and M by their chain means.</item>
/// <item><b>Recovery from noise outside the model:</b> a contaminated normal, t(3) BASS and Gaussian BASS beside it, the 90% intervals'
/// coverage of the true surface on a grid and of held-out responses, reported.</item>
/// <item><b>The cost row and the retention reading:</b> each part of a sweep per sweep at the cost report's base, n 200, and M 200 at
/// n = 10,000, under the g-form and the ridge, on the clock <c>BASS_COST_CLOCK</c> names; and the bytes a retained draw's working
/// response holds.</item>
/// </list>
/// </summary>
public sealed class BassStudentTQualificationTests(ITestOutputHelper output)
{
    // ---- Geweke ------------------------------------------------------------------------------------------------------------------

    private const int GewekeN = 12;
    private static readonly double[] GewekeX = [.. Enumerable.Range(0, GewekeN).Select(i => i / (GewekeN - 1.0))];
    private static readonly string[] GewekeNames = ["M", "ln tau", "ln sigma2", "ln v0", "mean ln v", "atan f0", "atan f6", "z0^2"];

    private static BassData GewekeData(double[] y) => new([new ContinuousPredictor("x", GewekeX)], new ScalarResponse(y));

    // Proper priors that read no response: the zero centre, a declared inverse gamma on σ², τ's unit-information Gamma of shape
    // `tauShape` (held at its auto value at 0), t(10) noise with c = `c`, so a c misplaced in a block would show.
    private static BassConfig GewekeConfig(double c, double tauShape) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 2, Support = new SupportRule(SupportKind.Count, 2) },
        Coefficients = BassCoefficientPrior.G(tauShape > 0 ? TauPrior.UnitInformationGamma(tauShape) : TauPrior.Held()),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.InverseGamma(3, 2),
        Noise = BassNoise.StudentT(10, c, BassNoiseReference.None),
        Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 0, Samples = 1 },
    };

    // Functions of (θ, y) with finite variance: f's prior is heavy-tailed through τ, so f enters through its arctangent.
    private static double[] Functions(BassState state, double sigma2, double[] v, double[] f, double[] y, double c)
        => [state.Count, state.LogTau, Math.Log(sigma2), Math.Log(v[0]), v.Average(x => Math.Log(x)), Math.Atan(f[0]), Math.Atan(f[GewekeN / 2]),
            (y[0] - f[0]) * (y[0] - f[0]) / (sigma2 * c * v[0])];

    private static WorkingResponse At(BassModel model, double[] v)
        => model.Working.Moved(default, [.. v.Select(x => 1.0 / (model.Config.Noise!.VarianceScale * x))], new ObservationLatents(v));

    // τ's prior is proper with tails the successive-conditional chain mixes through, or τ is held. Under Zellner–Siow's shape ½ the
    // chain sticks in ln τ's long left tail (a tiny τ makes the coefficients' prior wide, the simulated y wild and τ's posterior tiny
    // again): on two seeds of 300,000 steps its mean of ln τ sat 0.30 above and 0.58 below the prior's, beyond what batch means read,
    // with every other function in place. That is the test's reach, not the sweep's: the sweep's τ step is the same at every shape.
    [Theory]
    [InlineData(0.7, 5.0, 20261105)]
    [InlineData(0.7, 5.0, 20261106)]
    [InlineData(1.0, 5.0, 20261105)]
    [InlineData(0.7, 0.0, 20261105)]
    public void Geweke_TheSweepLeavesTheJointInvariant(double c, double tauShape, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int draws = 200_000, steps = 300_000, batches = 50;
        BassConfig config = GewekeConfig(c, tauShape);
        var law = (StudentTMixing)config.Noise!.Mixing;
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(config, GewekeData([.. GewekeX.Select(x => Math.Sin(7 * x))]));
        var sampler = new BassPriorSampler(prior);

        // The marginal-conditional simulator: (S, τ) from the prior on F, v from its prior, (σ², a) from theirs given (S, τ), y given all.
        // A structure outside F at the drawn v has no joint mass, so the pair is drawn again; none was, at this size.
        var marginal = new double[draws][];
        long redrawn = 0;
        (BassState S, double Sigma2, double[] V, double[] F, double[] Y) Draw()
        {
            while (true)
            {
                BassState state = sampler.DrawState(rng).State;
                double[] v = [.. Enumerable.Range(0, GewekeN).Select(_ => law.SamplePrior(rng))];
                var chain = new BassChain(prior, state.Terms, state.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
                if (chain.Follow(At(prior, v)) is not null) { redrawn++; continue; }
                BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, prior.Target, power: 0.0);
                double[] y = [.. Enumerable.Range(0, GewekeN).Select(i => fit.Fitted[i] + Math.Sqrt(fit.Sigma2 * c * v[i]) * Normal.SampleStandard(rng))];
                return (chain.Current, fit.Sigma2, v, fit.Fitted, y);
            }
        }
        for (int k = 0; k < draws; k++)
        {
            (BassState s, double sigma2, double[] v, double[] f, double[] y) = Draw();
            marginal[k] = Functions(s, sigma2, v, f, y, c);
        }

        // The successive-conditional simulator: y given θ, then one sweep of a chain on that y from θ.
        (BassState S, double Sigma2, double[] V, double[] F, double[] Y) theta = Draw();
        var successive = new double[steps][];
        long refused = 0, latentRefused = 0;
        for (int t = 0; t < steps; t++)
        {
            double[] y = [.. Enumerable.Range(0, GewekeN).Select(i => theta.F[i] + Math.Sqrt(theta.Sigma2 * c * theta.V[i]) * Normal.SampleStandard(rng))];
            BassModel model = BassModel.Create(config, GewekeData(y));
            var chain = new BassChain(model, theta.S.Terms, theta.S.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            if (chain.Follow(At(model, theta.V)) is not null) { refused++; continue; }
            chain.Step();
            MoveStat latent = chain.Scan.BlockStats().First(b => b.Key == "latent scales");
            latentRefused += latent.Attempts - latent.Accepted;
            BassSweepFit fit = chain.Structure.SweepFit!;
            double[] v = chain.Current.Fit.Working!.Latents!.ValueArray;
            successive[t] = Functions(chain.Current, fit.Sigma2, v, fit.Fitted, y, c);
            theta = (chain.Current, fit.Sigma2, v, fit.Fitted, y);
        }
        Assert.Equal(0, refused);

        double worst = 0;
        for (int g = 0; g < GewekeNames.Length; g++)
        {
            double[] a = [.. marginal.Select(r => r[g])], b = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1);
            double meanB = b.Average();
            int size = b.Length / batches;
            double[] means = [.. Enumerable.Range(0, batches).Select(k => b.Skip(k * size).Take(size).Average())];
            double seB2 = means.Sum(m => (m - means.Average()) * (m - means.Average())) / (batches - 1) / batches;
            double se = Math.Sqrt(varA / a.Length + seB2);
            bool constant = a.Min() == a.Max() && b.Min() == b.Max();   // a held τ: no spread to read
            double z = constant ? (a[0] == b[0] ? 0.0 : double.PositiveInfinity) : (meanA - meanB) / se;
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{GewekeNames[g],-10} marginal-conditional {meanA,10:F5}  successive {meanB,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"c {c}, tau shape {tauShape} (0 held), seed {seed}: structures redrawn off F {redrawn}; successive moves refused {refused}; latent draws refused {latentRefused} of {steps}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- large ν -------------------------------------------------------------------------------------------------------------------

    // The machine is shared: four chains at a time.
    private static readonly BassRunOptions Capped = new() { MaxDegreeOfParallelism = 4 };

    [Fact]
    public void LargeNu_ReachesGaussianBassReadouts()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassData data = BassStudentTTests.Data(n: 200, seed: 9);
        const int chains = 16;
        BassConfig gaussian = BassStudentTTests.Config(null, burnIn: 1000, samples: 1000, chains: chains) with
        {
            Sampling = new BassSampling { Chains = chains, MasterSeed = 23, BurnIn = 1000, Samples = 1000, Thin = 2 },
        };
        BassConfig t = gaussian with { Noise = BassNoise.StudentT(1e7), Sampling = gaussian.Sampling with { MasterSeed = 29 } };
        var a = Readouts(Bass.Fit(gaussian, data, Capped), data);
        var b = Readouts(Bass.Fit(t, data, Capped), data);
        double Z(double[] x, double[] y) => (x.Average() - y.Average()) / Math.Sqrt(Var(x) / x.Length + Var(y) / y.Length);
        double zSigma = Z(a.Sigma2, b.Sigma2), zCount = Z(a.Count, b.Count);
        double[] zFit = [.. Enumerable.Range(0, data.Response.Observations).Select(i => Z(a.Fit.Select(r => r[i]).ToArray(), b.Fit.Select(r => r[i]).ToArray()))];
        double rms = Math.Sqrt(zFit.Average(z => z * z)), max = zFit.Max(Math.Abs);
        output.WriteLine($"sigma2: Gaussian {a.Sigma2.Average():F5}, t(1e7) {b.Sigma2.Average():F5}, z {zSigma:F2}");
        output.WriteLine($"M: Gaussian {a.Count.Average():F3}, t(1e7) {b.Count.Average():F3}, z {zCount:F2}");
        output.WriteLine($"fit at {zFit.Length} cases: rms z {rms:F2}, max |z| {max:F2}; mean |difference| {Enumerable.Range(0, zFit.Length).Average(i => Math.Abs(a.Fit.Average(r => r[i]) - b.Fit.Average(r => r[i]))):E2}");
        Assert.True(Math.Abs(zSigma) < 4.0);
        Assert.True(Math.Abs(zCount) < 4.0);
        Assert.True(rms < 1.6, $"rms z {rms:F2}");
    }

    private static double Var(double[] x)
    {
        double m = x.Average();
        return x.Sum(v => (v - m) * (v - m)) / (x.Length - 1);
    }

    // Each chain's means: σ², M, and the fit at every training case.
    private static (double[] Sigma2, double[] Count, double[][] Fit) Readouts(BassFit fit, BassData data)
    {
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false });
        int chains = fit.Draws.Count, n = p.Cases, at = 0;
        var sigma2 = new double[chains];
        var count = new double[chains];
        var means = new double[chains][];
        for (int c = 0; c < chains; c++)
        {
            int s = fit.Draws[c].Count;
            sigma2[c] = fit.Draws[c].Average(d => d.Sigma2);
            count[c] = fit.Draws[c].Average(d => d.Structure.TermIds.Length);
            means[c] = new double[n];
            for (int d = 0; d < s; d++, at++)
                for (int i = 0; i < n; i++) means[c][i] += p.Latent[at, i, 0] / s;
        }
        return (sigma2, count, means);
    }

    // ---- recovery ------------------------------------------------------------------------------------------------------------------

    private static double Truth(double x0, double x1) => Math.Sin(2 * Math.PI * x0) + x1;

    // f = sin(2πx₀) + x₁ on [0, 1]², noise from 0.85·N(0, 0.25²) + 0.15·N(0, 2²): a contaminated normal, outside the t family.
    private static BassData Contaminated(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            y[i] = Truth(x0[i], x1[i]) + (rng.NextDouble() < 0.15 ? 2.0 : 0.25) * Normal.SampleStandard(rng);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
    }

    [Fact]
    public void Recovery_FromNoiseOutsideTheModel_WithCoverage()
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int sets = 6, n = 300;
        double[] grid = [.. Enumerable.Range(0, 7).Select(k => 0.05 + 0.15 * k)];
        var points = (from a in grid from b in grid select (a, b)).ToArray();
        var cases = new PredictionCases([new ContinuousColumn("x0", [.. points.Select(p => p.a)]), new ContinuousColumn("x1", [.. points.Select(p => p.b)])]);
        foreach ((string name, BassNoise? noise) in new[] { ("t(3)", BassNoise.StudentT(3)), ("Gaussian", (BassNoise?)null) })
        {
            int covered = 0, total = 0, heldCovered = 0, held = 0;
            double squares = 0;
            for (int set = 0; set < sets; set++)
            {
                BassData data = Contaminated(n, 100 + set), test = Contaminated(n, 200 + set);
                BassConfig config = BassStudentTTests.Config(noise, burnIn: 1500, samples: 500, chains: 4) with
                {
                    Sampling = new BassSampling { Chains = 4, MasterSeed = 31 + set, BurnIn = 1500, Samples = 500, Thin = 2 },
                };
                BassFit fit = Bass.Fit(config, data, Capped);
                BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Predictive = false });
                double[,] mean = p.Mean();
                for (int i = 0; i < points.Length; i++)
                {
                    double[] f = [.. Enumerable.Range(0, p.DrawCount).Select(d => p.Latent[d, i, 0]).Order()];
                    double truth = Truth(points[i].a, points[i].b), lo = f[(int)(0.05 * f.Length)], hi = f[(int)(0.95 * f.Length) - 1];
                    covered += truth >= lo && truth <= hi ? 1 : 0;
                    total++;
                    squares += (mean[i, 0] - truth) * (mean[i, 0] - truth);
                }
                var testCases = PredictionCases.Of(test);
                BassPrediction q = fit.Predict(data, testCases, new BassPredictOptions { Predictive = false });
                double[,] q05 = q.Quantile(0.05), q95 = q.Quantile(0.95);
                ReadOnlySpan<double> yTest = test.Response.Values;
                for (int i = 0; i < yTest.Length; i++)
                {
                    heldCovered += yTest[i] >= q05[i, 0] && yTest[i] <= q95[i, 0] ? 1 : 0;
                    held++;
                }
            }
            double coverage = (double)covered / total, heldCoverage = (double)heldCovered / held, rmse = Math.Sqrt(squares / total);
            output.WriteLine($"{name}: f's 90% intervals cover {coverage:P1} of {total} grid points over {sets} data sets, RMSE {rmse:F4}; " +
                             $"the 90% predictive intervals cover {heldCoverage:P1} of {held} held-out responses");
            if (noise is not null) Assert.True(coverage >= 0.7, $"t's coverage {coverage:P1}");
        }
    }

    // ---- the cost row and retention ------------------------------------------------------------------------------------------------

    [Fact]
    public void CostRow_PerPartOfASweep_AndRetention()
    {
        if (!BassQualification.Enabled) return;   // heavy
        string clock = BassCost.Cycles ? "thread cycles, indicative" : "wall clock";
        output.WriteLine($"Clock: {clock}. µs per sweep, each part's mean over the timed sweeps; M the mean count.");
        output.WriteLine("point                               form   M     structure  tau     a       v draws  Follow   sweep  | Gaussian sweep");
        var points = new (BassCost.Point Point, int Burn, int Timed)[]
        {
            (BassCost.Base, 2000, 400),
            (BassCost.Base with { N = 200 }, 2000, 400),
            (BassCost.Base with { N = 10000, Cap = 200 }, 0, 12),
        };
        foreach ((BassCost.Point x, int burn, int timed) in points)
            foreach (bool ridge in new[] { false, true })
            {
                BassCoefficientPrior form = ridge ? BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()) : BassCoefficientPrior.G();
                BassData data = BassCost.Data(x);
                BassConfig gaussian = BassCost.Config(x) with { Coefficients = form };
                IReadOnlyList<BassTerm> start;
                double logTau;
                if (burn == 0)
                {
                    // The largest point starts where a Gaussian chain burned to its cap stands: the t sweep's cost there, not its walk.
                    (_, BassChain burned, _) = BassCost.Burned(x);
                    start = burned.Current.Terms;
                    logTau = burned.Current.LogTau;
                }
                else
                {
                    start = [];
                    logTau = 0.0;
                }
                BassModel model = BassModel.Create(gaussian with { Noise = BassNoise.StudentT(4) }, data);
                var chain = new BassChain(model, start, logTau, new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(5));
                for (int s = 0; s < burn; s++) chain.Step();
                var parts = new double[5];
                double count = 0;
                var rng = new Xoshiro256PlusPlus(7);
                List<IBlock<ReversibleJumpChain<BassState>>> blocks = BassSweepSlots.Blocks(model);
                var noise = (VarianceMeanMixture)model.Observation;
                for (int s = 0; s < timed; s++)
                {
                    parts[0] += Timed(() => blocks[0].Update(chain.Moves, rng));
                    parts[1] += Timed(() => blocks[1].Update(chain.Moves, rng));
                    parts[2] += Timed(() => blocks[2].Update(chain.Moves, rng));
                    // The latent block in its two parts: the N draws of v, then Follow.
                    BassSweepFit fit = chain.Structure.SweepFit!;
                    WorkingResponse working = chain.Structure.DesignAt(chain.Current).Working;
                    double[] v = new double[working.Values.Length], w = new double[v.Length];
                    parts[3] += Timed(() =>
                    {
                        ReadOnlySpan<double> y = working.Values;
                        for (int i = 0; i < v.Length; i++)
                        {
                            v[i] = noise.Mixing.StepConditional(working.Latents!.Values[i], y[i] - fit.Fitted[i], fit.Sigma2 * noise.VarianceScale, 0.0, 1.0, rng);
                            w[i] = 1.0 / (noise.VarianceScale * v[i]);
                        }
                    });
                    WorkingResponse moved = working.Moved(working.Offset, w, new ObservationLatents(v));
                    parts[4] += Timed(() => BassStructureMoves.Follow(chain.Moves, moved));
                    count += chain.Current.Count;
                }
                // A Gaussian chain's sweep at the same structure, for the scale.
                BassModel plain = BassModel.Create(gaussian, data);
                var reference = new BassChain(plain, chain.Current.Terms, chain.Current.LogTau, new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(5));
                double gaussianSweep = 0;
                int plainTimed = Math.Max(timed, 50);
                for (int s = 0; s < plainTimed; s++) gaussianSweep += Timed(() => reference.Step());
                double us(double ticks) => ticks / timed / BassCost.TicksPerSecond * 1e6;
                output.WriteLine($"{x.Key,-36}{(ridge ? "ridge" : "g"),-7}{count / timed,-6:F1}{us(parts[0]),-11:F1}{us(parts[1]),-8:F1}{us(parts[2]),-8:F1}" +
                                 $"{us(parts[3]),-9:F1}{us(parts[4]),-9:F1}{us(parts.Sum()),-7:F1}| {gaussianSweep / plainTimed / BassCost.TicksPerSecond * 1e6:F1}");
            }

        // Retention: a retained draw shares the working response it was composed on; under t it holds the weights and v, 2N doubles.
        foreach (int n in new[] { 1_000, 100_000 })
        {
            var response = new ScalarResponse([.. Enumerable.Range(0, n).Select(i => Math.Sin(i))]);
            WorkingResponse start = new VarianceMeanMixture(BassNoise.StudentT(4)).Working(response);
            double[] v = [.. Enumerable.Range(0, n).Select(i => 1.0 + 0.5 * Math.Cos(i))];
            double[] w = [.. v.Select(x => 1.0 / x)];
            var latents = new ObservationLatents(v);
            GC.Collect();
            long before = GC.GetAllocatedBytesForCurrentThread();
            WorkingResponse moved = start.Moved(default, w, latents);
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
            long retained = (long)moved.Weights.Length * 8 + (long)moved.Latents!.Values.Length * 8;
            output.WriteLine($"N {n:N0}: a move allocates {allocated / 1024.0:F1} KB (its weights' copy, their logs, the object); a retained draw holds " +
                             $"{retained / 1024.0:F1} KB beyond its own, the weights and v, of which v is the generative state and the weights a function of it");
            GC.KeepAlive(moved);
        }
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
