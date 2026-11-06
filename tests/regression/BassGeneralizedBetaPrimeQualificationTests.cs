using System;
using System.Collections.Generic;
using System.Diagnostics;
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
/// The generalized beta prime's long gates (E1's slice 4), behind <c>BASS_QUALIFICATION</c> = 1; without it each returns at once.
/// <list type="bullet">
/// <item><b>Geweke's joint-distribution test</b> with v moved by its slice step: the horseshoe with σ² integrated, the horseshoe with a drawn
/// skew and σ² on the record, and the Strawderman–Berger.</item>
/// <item><b>Recovery from noise outside the model:</b> a contaminated normal, the horseshoe under <c>hbass</c>'s half-Cauchy σ² and Gaussian
/// BASS beside it, the 90% intervals' coverage of the true surface on a grid and of held-out responses, reported.</item>
/// <item><b>The cost row</b>, each part of a sweep at the cost report's points, v's slice steps with their evaluations; and the pointwise
/// density's price per read, the opt-in's reason.</item>
/// </list>
/// </summary>
public sealed class BassGeneralizedBetaPrimeQualificationTests(ITestOutputHelper output)
{
    // ---- Geweke --------------------------------------------------------------------------------------------------------------------

    private const int GewekeN = 12;
    private static readonly double[] GewekeX = [.. Enumerable.Range(0, GewekeN).Select(i => i / (GewekeN - 1.0))];
    private static readonly string[] GewekeNames = ["M", "ln tau", "ln sigma2", "ln v0", "mean ln v", "atan f0", "atan f6", "beta", "z0^2"];

    private static BassData GewekeData(double[] y) => new([new ContinuousPredictor("x", GewekeX)], new ScalarResponse(y));

    private static BassNoise GewekeNoise(string form) => form switch
    {
        "horseshoe" => BassNoise.Horseshoe(reference: BassNoiseReference.None),
        "horseshoe-skewed" => BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5), reference: BassNoiseReference.None),
        _ => BassNoise.StrawdermanBerger(reference: BassNoiseReference.None),
    };

    private static BassConfig GewekeConfig(string form) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 2, Support = new SupportRule(SupportKind.Count, 2) },
        Coefficients = BassCoefficientPrior.G(TauPrior.UnitInformationGamma(5)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.InverseGamma(3, 2),
        Noise = GewekeNoise(form),
        Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 0, Samples = 1 },
    };

    private static WorkingResponse At(BassModel model, double[] v, double sigma2, double beta)
    {
        var mixture = (VarianceMeanMixture)model.Observation;
        return mixture.At(model.Working, new ObservationLatents(v, beta, mixture.Skewed ? sigma2 : null));
    }

    [Theory]
    [InlineData("horseshoe", 20261106)]
    [InlineData("horseshoe", 20261107)]
    [InlineData("horseshoe-skewed", 20261106)]
    [InlineData("strawderman-berger", 20261106)]
    public void Geweke_TheSweepLeavesTheJointInvariant(string form, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int draws = 400_000, steps = 600_000, batches = 50;
        BassConfig config = GewekeConfig(form);
        BassNoise noise = config.Noise!;
        double c = noise.VarianceScale;
        InverseGamma sigmaPrior = InverseGamma.ShapeScale(3, 2);
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(config, GewekeData([.. GewekeX.Select(x => Math.Sin(7 * x))]));
        var sampler = new BassPriorSampler(prior);
        IMixingLaw law = noise.Mixing;

        double[] Functions(BassState s, double sigma2, double[] v, double[] f, double[] y, double beta)
            => [s.Count, s.LogTau, Math.Log(sigma2), Math.Log(v[0]), v.Average(x => Math.Log(x)), Math.Atan(f[0]), Math.Atan(f[GewekeN / 2]),
                beta, Math.Pow(y[0] - f[0] - Math.Sqrt(sigma2) * beta * v[0], 2) / (sigma2 * c * v[0])];
        double[] Respond(double[] f, double sigma2, double[] v, double beta)
            => [.. Enumerable.Range(0, GewekeN).Select(i => f[i] + Math.Sqrt(sigma2) * (beta * v[i] + Math.Sqrt(c * v[i]) * Normal.SampleStandard(rng)))];

        // The marginal-conditional simulator: β from its prior, v from the GBP, σ² from its prior (held under a skew; drawn with a otherwise).
        var marginal = new double[draws][];
        (BassState S, double Sigma2, double[] V, double[] F, double Beta) Draw()
        {
            while (true)
            {
                BassState state = sampler.DrawState(rng).State;
                double beta = noise.SkewPrior is { IsHeld: false } b ? b.Mean + b.Sd * Normal.SampleStandard(rng) : noise.Skew;
                double[] v = [.. Enumerable.Range(0, GewekeN).Select(_ => law.SamplePrior(rng))];
                double held = sigmaPrior.Sample(rng);
                var chain = new BassChain(prior, state.Terms, state.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
                if (chain.Follow(At(prior, v, held, beta)) is not null) continue;
                BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, prior.Target, power: 0.0);
                return (chain.Current, fit.Sigma2, v, fit.Fitted, beta);
            }
        }
        for (int k = 0; k < draws; k++)
        {
            var d = Draw();
            marginal[k] = Functions(d.S, d.Sigma2, d.V, d.F, Respond(d.F, d.Sigma2, d.V, d.Beta), d.Beta);
        }

        // The successive-conditional simulator: y given θ, then one sweep on that y, v stepped from θ's.
        var theta = Draw();
        var successive = new double[steps][];
        long refused = 0, blockRefused = 0;
        for (int t = 0; t < steps; t++)
        {
            double[] y = Respond(theta.F, theta.Sigma2, theta.V, theta.Beta);
            BassModel model = BassModel.Create(config, GewekeData(y));
            var chain = new BassChain(model, theta.S.Terms, theta.S.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            if (chain.Follow(At(model, theta.V, theta.Sigma2, theta.Beta)) is not null) { refused++; continue; }
            chain.Step();
            foreach (MoveStat s in chain.Scan.BlockStats().Where(b => b.Key is "latent scales" or "scale" or "noise shape")) blockRefused += s.Attempts - s.Accepted;
            BassSweepFit fit = chain.Structure.SweepFit!;
            ObservationLatents latents = chain.Current.Fit.Working!.Latents!;
            double sigma2 = latents.Sigma2 ?? fit.Sigma2;
            successive[t] = Functions(chain.Current, sigma2, latents.ValueArray, fit.Fitted, y, latents.Skew);
            theta = (chain.Current, sigma2, latents.ValueArray, fit.Fitted, latents.Skew);
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
        output.WriteLine($"{form}, seed {seed}: moves refused {refused}; block draws refused {blockRefused} of {steps}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- recovery ------------------------------------------------------------------------------------------------------------------

    private static readonly BassRunOptions Capped = new() { MaxDegreeOfParallelism = 4 };

    private static double Truth(double x0, double x1) => Math.Sin(2 * Math.PI * x0) + x1;

    // f = sin(2πx₀) + x₁ on [0, 1]², noise from 0.85·N(0, 0.25²) + 0.15·N(0, 2²): a contaminated normal, outside the GBP family.
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
        foreach ((string name, BassNoise? noise) in new[] { ("horseshoe", BassNoise.Horseshoe()), ("Gaussian", (BassNoise?)null) })
        {
            int covered = 0, total = 0, heldCovered = 0, held = 0;
            double squares = 0;
            for (int set = 0; set < sets; set++)
            {
                BassData data = Contaminated(n, 100 + set), test = Contaminated(n, 200 + set);
                BassConfig config = BassStudentTTests.Config(noise) with
                {
                    Sampling = new BassSampling { Chains = 4, MasterSeed = 31 + set, BurnIn = 1500, Samples = 500, Thin = 2 },
                };
                if (noise is not null) config = config with { Sigma2 = Sigma2Prior.HalfCauchy() };
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
                BassPrediction q = fit.Predict(data, PredictionCases.Of(test), new BassPredictOptions { Predictive = false });
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
            if (noise is not null) Assert.True(coverage >= 0.7, $"the horseshoe's coverage {coverage:P1}");
        }
    }

    // ---- the cost row --------------------------------------------------------------------------------------------------------------

    [Fact]
    public void CostRow_PerPartOfASweep_AndThePointwiseRead()
    {
        if (!BassQualification.Enabled) return;   // heavy
        string clock = BassCost.Cycles ? "thread cycles, indicative" : "wall clock";
        output.WriteLine($"Clock: {clock}. µs per sweep, each part's mean over the timed sweeps; M the mean count; evaluations per v step.");
        output.WriteLine("point                               form   M     structure  tau     a       v steps  evals  Follow   sweep  | t(4) v draws | Gaussian sweep");
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
                    (_, BassChain burned, _) = BassCost.Burned(x);
                    start = burned.Current.Terms;
                    logTau = burned.Current.LogTau;
                }
                else
                {
                    start = [];
                    logTau = 0.0;
                }
                BassModel model = BassModel.Create(gaussian with { Noise = BassNoise.Horseshoe() }, data);
                var chain = new BassChain(model, start, logTau, new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(5));
                for (int s = 0; s < burn; s++) chain.Step();
                var parts = new double[5];
                double count = 0, tDraws = 0;
                long evaluations = 0, stepsTaken = 0;
                var rng = new Xoshiro256PlusPlus(7);
                List<IBlock<ReversibleJumpChain<BassState>>> blocks = BassSweepSlots.Blocks(model);
                var noise = (VarianceMeanMixture)model.Observation;
                var mixing = (GbpMixing)noise.Mixing;
                var t4 = new StudentTMixing(4);
                for (int s = 0; s < timed; s++)
                {
                    parts[0] += Timed(() => blocks[0].Update(chain.Moves, rng));
                    parts[1] += Timed(() => blocks[1].Update(chain.Moves, rng));
                    parts[2] += Timed(() => blocks[2].Update(chain.Moves, rng));
                    // The latent block in its two parts: the N slice steps of v, then Follow; t(4)'s N exact draws at the same residuals beside.
                    BassSweepFit fit = chain.Structure.SweepFit!;
                    WorkingResponse working = chain.Structure.DesignAt(chain.Current).Working;
                    ObservationLatents latents = working.Latents!;
                    double[] v = new double[latents.Values.Length], w = new double[v.Length];
                    double unit = fit.Sigma2 * noise.VarianceScale;
                    parts[3] += Timed(() =>
                    {
                        ReadOnlySpan<double> y = working.Response.Values, current = latents.Values;
                        for (int i = 0; i < v.Length; i++)
                        {
                            v[i] = mixing.StepConditional(current[i], y[i] - fit.Fitted[i], unit, 0.0, 1.0, rng);
                            w[i] = 1.0 / (noise.VarianceScale * v[i]);
                        }
                    });
                    tDraws += Timed(() =>
                    {
                        ReadOnlySpan<double> y = working.Response.Values;
                        for (int i = 0; i < v.Length; i++) _ = t4.SampleConditional(y[i] - fit.Fitted[i], unit, 0.0, 1.0, rng);
                    });
                    // The evaluations a step spends, counted apart on a copy of the stream, outside the timing.
                    if (s % 20 == 0)
                    {
                        var counting = new Xoshiro256PlusPlus(11);
                        ReadOnlySpan<double> y = working.Response.Values, current = latents.Values;
                        for (int i = 0; i < v.Length; i++)
                        {
                            double u = Math.Log(current[i]), r = y[i] - fit.Fitted[i];
                            double logR = -Math.Log(2.0 * unit) + 2.0 * Math.Log(Math.Abs(r));
                            double Density(double z) => 0.0 * z - LogSoftplus(z) - Math.Exp(logR - z);
                            double target = Density(u);
                            evaluations += 1 + SliceSampler.Step(ref u, ref target, 2.0, 64, Density, counting);
                            stepsTaken++;
                        }
                    }
                    WorkingResponse moved = working.Moved(working.Offset, w, new ObservationLatents(v));
                    parts[4] += Timed(() => BassStructureMoves.Follow(chain.Moves, moved));
                    count += chain.Current.Count;
                }
                BassModel plain = BassModel.Create(gaussian, data);
                var reference = new BassChain(plain, chain.Current.Terms, chain.Current.LogTau, new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(5));
                double gaussianSweep = 0;
                int plainTimed = Math.Max(timed, 50);
                for (int s = 0; s < plainTimed; s++) gaussianSweep += Timed(() => reference.Step());
                double us(double ticks) => ticks / timed / BassCost.TicksPerSecond * 1e6;
                output.WriteLine($"{x.Key,-36}{(ridge ? "ridge" : "g"),-7}{count / timed,-6:F1}{us(parts[0]),-11:F1}{us(parts[1]),-8:F1}{us(parts[2]),-8:F1}" +
                                 $"{us(parts[3]),-9:F1}{(double)evaluations / stepsTaken,-7:F2}{us(parts[4]),-9:F1}{us(parts.Sum()),-7:F1}| {us(tDraws),-13:F1}| " +
                                 $"{gaussianSweep / plainTimed / BassCost.TicksPerSecond * 1e6:F1}");
            }

        // The pointwise density's price: one read is an integral over v, against the closed t's.
        var horseshoe = new GbpMixing(GbpParameters.Horseshoe);
        var t = new StudentTMixing(4);
        var reads = new Xoshiro256PlusPlus(13);
        double[] residuals = [.. Enumerable.Range(0, 2_000).Select(_ => 2.0 * Normal.SampleStandard(reads))];
        double sink = 0.0;
        foreach (double r in residuals.Take(200)) sink += horseshoe.MarginalLogDensity(r, 1.0, 0.0);   // warm-up
        var watch = Stopwatch.StartNew();
        foreach (double r in residuals) sink += horseshoe.MarginalLogDensity(r, 1.0, 0.0);
        double gbpRead = watch.Elapsed.TotalMilliseconds * 1e3 / residuals.Length;
        watch.Restart();
        for (int k = 0; k < 100; k++)
            foreach (double r in residuals) sink += t.MarginalLogDensity(r, 1.0, 0.0);
        double tRead = watch.Elapsed.TotalMilliseconds * 1e3 / (100.0 * residuals.Length);
        output.WriteLine($"The pointwise read (wall clock): the horseshoe's integral {gbpRead:F2} µs, t(4)'s closed form {tRead:F3} µs; " +
                         $"at N 1,000 and 2,000 draws the horseshoe's pointwise density costs {gbpRead * 2e6 / 1e6:F1} s ({sink:G3})");
    }

    // ln(1 + eᶻ), the horseshoe's (a + b)·ln(1 + e^(pu)) at p = 1, a + b = 1.
    private static double LogSoftplus(double z) => z > 0.0 ? z + Math.Log(1.0 + Math.Exp(-z)) : Math.Log(1.0 + Math.Exp(z));

    private static long Timed(Action action)
    {
        long t0 = BassCost.Tick();
        action();
        return BassCost.Tick() - t0;
    }

    private static long Timed(Func<bool> action) => Timed(() => { action(); });

    private static long Timed(Func<BassRefusal?> action) => Timed(() => { action(); });
}
