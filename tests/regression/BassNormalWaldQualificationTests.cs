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
/// The Normal-Wald's, the skewed t's and the generalized hyperbolic's long gates (E1's slice 3), behind <c>BASS_QUALIFICATION</c> = 1;
/// without it each returns at once.
/// <list type="bullet">
/// <item><b>Geweke's joint-distribution test</b> with the skew and γ drawn: σ² on the record, v's GIG block at the record's γ, the
/// noise-shape block; and with γ alone drawn, where σ² is integrated.</item>
/// <item><b>Recovery</b> of β, γ and σ² on NIG noise, and the surface.</item>
/// <item><b>The <c>lag_beta</c> read</b> (decision 3): chains with the skew free from the start against chains that hold it at 0 for a
/// burn-in, then free it, as GBASS's <c>lag_beta</c> does.</item>
/// <item><b>The cost row</b>, by block.</item>
/// </list>
/// </summary>
public sealed class BassNormalWaldQualificationTests(ITestOutputHelper output)
{
    // ---- Geweke --------------------------------------------------------------------------------------------------------------------

    private const int GewekeN = 12;
    private static readonly double[] GewekeX = [.. Enumerable.Range(0, GewekeN).Select(i => i / (GewekeN - 1.0))];
    private static readonly string[] GewekeNames = ["M", "ln tau", "ln sigma2", "ln v0", "mean ln v", "atan f0", "atan f6", "beta", "ln gamma", "z0^2"];

    private static BassData GewekeData(double[] y) => new([new ContinuousPredictor("x", GewekeX)], new ScalarResponse(y));

    private static BassNoise GewekeNoise(string form) => form switch
    {
        "nw" => BassNoise.NormalWald(new BassNoisePrior(1.5, 1.0), new BassNoisePrior(0.0, 0.7), reference: BassNoiseReference.None),
        "nw-symmetric" => BassNoise.NormalWald(new BassNoisePrior(1.5, 1.0), reference: BassNoiseReference.None),
        "skewed-t" => BassNoise.SkewedT(5, new BassNoisePrior(0.0, 0.5), reference: BassNoiseReference.None),
        _ => BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), new BassNoisePrior(0.3, 0.5), reference: BassNoiseReference.None),
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

    private static WorkingResponse At(BassModel model, double[] v, double sigma2, double beta, double? gamma)
    {
        var mixture = (VarianceMeanMixture)model.Observation;
        return mixture.At(model.Working, new ObservationLatents(v, beta, mixture.Skewed ? sigma2 : null, gamma));
    }

    [Theory]
    [InlineData("nw", 20261106)]
    [InlineData("nw", 20261107)]
    [InlineData("nw-symmetric", 20261106)]
    [InlineData("skewed-t", 20261106)]
    [InlineData("gh", 20261106)]
    public void Geweke_TheSweepLeavesTheJointInvariant(string form, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        // Twice slice 2's size: at 200,000 and 300,000 one function of one case reached |z| 4.08, which doubling the same seed took to 0.31.
        const int draws = 400_000, steps = 600_000, batches = 50;
        BassConfig config = GewekeConfig(form);
        BassNoise noise = config.Noise!;
        double c = noise.VarianceScale;
        InverseGamma sigmaPrior = InverseGamma.ShapeScale(3, 2);
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(config, GewekeData([.. GewekeX.Select(x => Math.Sin(7 * x))]));
        var sampler = new BassPriorSampler(prior);

        double[] Functions(BassState s, double sigma2, double[] v, double[] f, double[] y, double beta, double? gamma)
            => [s.Count, s.LogTau, Math.Log(sigma2), Math.Log(v[0]), v.Average(x => Math.Log(x)), Math.Atan(f[0]), Math.Atan(f[GewekeN / 2]),
                beta, Math.Log(gamma ?? 1.0), Math.Pow(y[0] - f[0] - Math.Sqrt(sigma2) * beta * v[0], 2) / (sigma2 * c * v[0])];
        double[] Respond(double[] f, double sigma2, double[] v, double beta)
            => [.. Enumerable.Range(0, GewekeN).Select(i => f[i] + Math.Sqrt(sigma2) * (beta * v[i] + Math.Sqrt(c * v[i]) * Normal.SampleStandard(rng)))];

        // The marginal-conditional simulator: β and γ from their priors, v at γ, σ² from its prior (held under a skew; drawn with a otherwise).
        var marginal = new double[draws][];
        (BassState S, double Sigma2, double[] V, double[] F, double Beta, double? Gamma) Draw()
        {
            while (true)
            {
                BassState state = sampler.DrawState(rng).State;
                double beta = noise.SkewPrior is { IsHeld: false } b ? b.Mean + b.Sd * Normal.SampleStandard(rng) : noise.Skew;
                double? gamma = noise.GammaPrior switch
                {
                    { IsHeld: false } g => Normal.MeanSd(g.Mean, g.Sd).SampleOn(rng, 0.0, double.PositiveInfinity),
                    { } g => g.Mean,
                    null => null,
                };
                IMixingLaw law = noise.MixingAt(gamma);
                double[] v = [.. Enumerable.Range(0, GewekeN).Select(_ => law.SamplePrior(rng))];
                double held = sigmaPrior.Sample(rng);
                var chain = new BassChain(prior, state.Terms, state.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
                if (chain.Follow(At(prior, v, held, beta, gamma)) is not null) continue;
                BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, prior.Target, power: 0.0);
                return (chain.Current, fit.Sigma2, v, fit.Fitted, beta, gamma);
            }
        }
        for (int k = 0; k < draws; k++)
        {
            var d = Draw();
            marginal[k] = Functions(d.S, d.Sigma2, d.V, d.F, Respond(d.F, d.Sigma2, d.V, d.Beta), d.Beta, d.Gamma);
        }

        // The successive-conditional simulator: y given θ, then one sweep on that y.
        var theta = Draw();
        var successive = new double[steps][];
        long refused = 0, blockRefused = 0;
        for (int t = 0; t < steps; t++)
        {
            double[] y = Respond(theta.F, theta.Sigma2, theta.V, theta.Beta);
            BassModel model = BassModel.Create(config, GewekeData(y));
            var chain = new BassChain(model, theta.S.Terms, theta.S.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            if (chain.Follow(At(model, theta.V, theta.Sigma2, theta.Beta, theta.Gamma)) is not null) { refused++; continue; }
            chain.Step();
            foreach (MoveStat s in chain.Scan.BlockStats().Where(b => b.Key is "latent scales" or "scale" or "noise shape")) blockRefused += s.Attempts - s.Accepted;
            BassSweepFit fit = chain.Structure.SweepFit!;
            ObservationLatents latents = chain.Current.Fit.Working!.Latents!;
            double sigma2 = latents.Sigma2 ?? fit.Sigma2;
            successive[t] = Functions(chain.Current, sigma2, latents.ValueArray, fit.Fitted, y, latents.Skew, latents.Gamma);
            theta = (chain.Current, sigma2, latents.ValueArray, fit.Fitted, latents.Skew, latents.Gamma);
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

    private const double TrueSigma = 0.3, TrueBeta = 1.0, TrueGamma = 1.5;

    private static double Truth(double x0, double x1) => Math.Sin(2 * Math.PI * x0) + x1;

    // NIG noise in the core's terms: σ(βv + √v·ζ), v ~ GIG(−½, γ², 1).
    private static BassData NigData(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var law = GigMixing.NormalWald(TrueGamma);
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            double v = law.SamplePrior(rng);
            y[i] = Truth(x0[i], x1[i]) + TrueSigma * (TrueBeta * v + Math.Sqrt(v) * Normal.SampleStandard(rng));
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
    }

    private static BassNoise RecoveryNoise()
    {
        (double m, double s) = NormalWaldCalibration.GammaPrior();
        return BassNoise.NormalWald(new BassNoisePrior(m, s), new BassNoisePrior(0.0, 2.0));
    }

    private static BassConfig RecoveryConfig(int seed, int burnIn = 3000, int samples = 1000) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 20 },
        Noise = RecoveryNoise(),
        Sampling = new BassSampling { Chains = 4, MasterSeed = seed, BurnIn = burnIn, Samples = samples, Thin = 2 },
    };

    [Fact]
    public void TheShape_IsRecovered_OnNigNoise()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassData data = NigData(600, 21);
        BassFit fit = Bass.Fit(RecoveryConfig(5), data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        List<ObservationLatents> latents = [.. fit.Draws.SelectMany(c => c).Select(d => d.Latents!)];
        (double Lo, double Mid, double Hi) Band(IEnumerable<double> xs)
        {
            double[] s = [.. xs.Order()];
            return (s[(int)(0.005 * s.Length)], s[s.Length / 2], s[(int)(0.995 * s.Length) - 1]);
        }
        var beta = Band(latents.Select(l => l.Skew));
        var gamma = Band(latents.Select(l => l.Gamma!.Value));
        var sigma2 = Band(latents.Select(l => l.Sigma2!.Value));
        output.WriteLine($"β {beta.Mid:F3} [{beta.Lo:F3}, {beta.Hi:F3}] truth {TrueBeta}; γ {gamma.Mid:F3} [{gamma.Lo:F3}, {gamma.Hi:F3}] truth {TrueGamma}; " +
                         $"σ² {sigma2.Mid:F4} [{sigma2.Lo:F4}, {sigma2.Hi:F4}] truth {TrueSigma * TrueSigma}");
        BassData grid = NigData(400, 22);
        BassPrediction p = fit.Predict(data, PredictionCases.Of(grid), new BassPredictOptions { Predictive = false });
        double[,] surface = p.Mean();
        ReadOnlySpan<double> g0 = ((ContinuousPredictor)grid.Predictors[0]).Values.Span, g1 = ((ContinuousPredictor)grid.Predictors[1]).Values.Span;
        int cases = grid.Response.Values.Length;
        double shift = TrueSigma * TrueBeta / TrueGamma, squares = 0.0, meanSquares = 0.0;
        double meanShift = Enumerable.Range(0, p.DrawCount).Sum(d => p.Weights[d] * p.MeanShifts[d]) / Enumerable.Range(0, p.DrawCount).Sum(d => p.Weights[d]);
        for (int i = 0; i < cases; i++)
        {
            double truth = Truth(g0[i], g1[i]);
            squares += Math.Pow(surface[i, 0] - truth, 2);
            meanSquares += Math.Pow(surface[i, 0] + meanShift - (truth + shift), 2);
        }
        output.WriteLine($"the surface f against the truth: RMSE {Math.Sqrt(squares / cases):F4}; the mean f + σβE[v] against the true mean: " +
                         $"RMSE {Math.Sqrt(meanSquares / cases):F4} (shift {meanShift:F4} against {shift:F4})");
        Assert.InRange(TrueBeta, beta.Lo, beta.Hi);
        Assert.InRange(TrueGamma, gamma.Lo, gamma.Hi);
        Assert.InRange(TrueSigma * TrueSigma, sigma2.Lo, sigma2.Hi);
    }

    // ---- lag_beta ------------------------------------------------------------------------------------------------------------------

    /// <summary>
    /// Decision 3's read of GBASS's <c>lag_beta</c>, which holds the skew at its prior's mean for the first iterations. Eight chains each way
    /// on NIG data from the empty structure: free from the start, or held at 0 for 500 sweeps (a chain of the same noise with β held, handed
    /// its state) and then free. Read per arm: the sweeps until M first reaches the median of a long reference run's M, β's mean over the
    /// first 200 sweeps, and M, β and σ² over sweeps 1,500 to 2,000. The rule stated beside it: adopt the hold only if the free arm takes
    /// half as long again to reach M, or ends elsewhere.
    /// </summary>
    [Fact]
    public void LagBeta_Read()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassData data = NigData(600, 21);
        BassConfig free = RecoveryConfig(5) with { Sampling = new BassSampling { Chains = 1, MasterSeed = 5, BurnIn = 0, Samples = 1 } };
        BassConfig held = free with { Noise = BassNoise.NormalWald(RecoveryNoise().GammaPrior) };
        BassModel freeModel = BassModel.Create(free, data), heldModel = BassModel.Create(held, data);
        BassFit reference = Bass.Fit(RecoveryConfig(5, 4000, 1000), data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        double target = reference.Draws.SelectMany(c => c).Select(d => (double)d.Structure.TermIds.Length).Order().ElementAt(reference.Draws.Sum(c => c.Count) / 2);
        const int chains = 8, sweeps = 2_000, hold = 500;
        foreach (bool lag in new[] { false, true })
        {
            var reach = new List<int>();
            var early = new List<double>();
            var lateM = new List<double>();
            var lateBeta = new List<double>();
            var lateSigma2 = new List<double>();
            for (int k = 0; k < chains; k++)
            {
                BassChain chain;
                int start = 0;
                var ms = new List<double>();
                var betas = new List<double>();
                if (lag)
                {
                    var holding = new BassChain(heldModel, [], 0.0, new Xoshiro256PlusPlus(100 + k), new Xoshiro256PlusPlus(200 + k));
                    for (int s = 0; s < hold; s++)
                    {
                        holding.Step();
                        ms.Add(holding.Current.Count);
                        betas.Add(0.0);
                    }
                    ObservationLatents h = holding.Current.Fit.Working!.Latents!;
                    chain = new BassChain(freeModel, holding.Current.Terms, holding.Current.LogTau, new Xoshiro256PlusPlus(300 + k), new Xoshiro256PlusPlus(400 + k));
                    Assert.Null(chain.Follow(At(freeModel, h.ValueArray, holding.Structure.SweepFit!.Sigma2, 0.0, h.Gamma)));
                    start = hold;
                }
                else chain = new BassChain(freeModel, [], 0.0, new Xoshiro256PlusPlus(100 + k), new Xoshiro256PlusPlus(200 + k));
                var sigma2s = new List<double>();
                for (int s = start; s < sweeps; s++)
                {
                    chain.Step();
                    ObservationLatents l = chain.Current.Fit.Working!.Latents!;
                    ms.Add(chain.Current.Count);
                    betas.Add(l.Skew);
                    sigma2s.Add(l.Sigma2!.Value);
                }
                int first = ms.FindIndex(m => m >= target);
                reach.Add(first < 0 ? sweeps : first);
                early.Add(betas.Take(200).Average());
                lateM.Add(ms.Skip(1500).Average());
                lateBeta.Add(betas.Skip(1500).Average());
                lateSigma2.Add(sigma2s.Skip(sigma2s.Count - 500).Average());
            }
            output.WriteLine($"{(lag ? $"held for {hold}" : "free"),-12} sweeps to M ≥ {target}: median {reach.Order().ElementAt(chains / 2)} [{string.Join(", ", reach)}]; " +
                             $"β over sweeps 0-200 {early.Average():F3}; over 1,500-2,000 M {lateM.Average():F2}, β {lateBeta.Average():F3}, σ² {lateSigma2.Average():F4}");
        }
        var referenceLatents = reference.Draws.SelectMany(c => c).Select(d => d.Latents!).ToList();
        output.WriteLine($"reference: M median {target}, β {referenceLatents.Average(l => l.Skew):F3}, σ² {referenceLatents.Average(l => l.Sigma2!.Value):F4}");
    }

    // ---- cost ----------------------------------------------------------------------------------------------------------------------

    [Fact]
    public void CostRow_PerBlockOfASweep()
    {
        if (!BassQualification.Enabled) return;   // heavy
        output.WriteLine($"Clock: {(BassCost.Cycles ? "thread cycles, indicative" : "wall clock")}. µs per sweep, each block's mean; Normal-Wald with β and γ drawn.");
        var points = new (BassCost.Point Point, int Burn, int Timed)[]
        {
            (BassCost.Base, 2000, 400),
            (BassCost.Base with { N = 200 }, 2000, 400),
            (BassCost.Base with { N = 10000, Cap = 200 }, 0, 10),
        };
        foreach ((BassCost.Point x, int burn, int timed) in points)
        {
            BassData data = BassCost.Data(x);
            BassConfig gaussian = BassCost.Config(x) with { Coefficients = BassCoefficientPrior.G() };
            IReadOnlyList<BassTerm> start = [];
            double logTau = 0.0;
            if (burn == 0)
            {
                (_, BassChain burned, _) = BassCost.Burned(x);
                start = burned.Current.Terms;
                logTau = burned.Current.LogTau;
            }
            BassModel model = BassModel.Create(gaussian with { Noise = RecoveryNoise() }, data);
            var chain = new BassChain(model, start, logTau, new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(5));
            for (int s = 0; s < burn; s++) chain.Step();
            List<IBlock<ReversibleJumpChain<BassState>>> blocks = BassSweepSlots.Blocks(model);
            var parts = new double[blocks.Count];
            double count = 0;
            var rng = new Xoshiro256PlusPlus(7);
            for (int s = 0; s < timed; s++)
            {
                for (int b = 0; b < blocks.Count; b++)
                {
                    IBlock<ReversibleJumpChain<BassState>> block = blocks[b];
                    long t0 = BassCost.Tick();
                    block.Update(chain.Moves, rng);
                    parts[b] += BassCost.Tick() - t0;
                }
                count += chain.Current.Count;
            }
            BassModel plain = BassModel.Create(gaussian, data);
            var reference = new BassChain(plain, chain.Current.Terms, chain.Current.LogTau, new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(5));
            double gaussianSweep = 0;
            int plainTimed = Math.Max(timed, 50);
            for (int s = 0; s < plainTimed; s++)
            {
                long t0 = BassCost.Tick();
                reference.Step();
                gaussianSweep += BassCost.Tick() - t0;
            }
            double us(double ticks, int n) => ticks / n / BassCost.TicksPerSecond * 1e6;
            output.WriteLine($"{x.Key,-36} M {count / timed:F1}: " + string.Join(", ", blocks.Select((b, i) => $"{b.Name} {us(parts[i], timed):F1}")) +
                             $"; sweep {us(parts.Sum(), timed):F1} | Gaussian {us(gaussianSweep, plainTimed):F1}");
        }
    }
}
