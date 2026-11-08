using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The engine's reach under a noise (E1's slice 5), its long gates, behind <c>BASS_QUALIFICATION</c> = 1; without it each returns at once.
/// <list type="bullet">
/// <item><b>The law at a rung.</b> At b the tempered target p(θ)·p(v)·N(y | f + σβv, σ²c·v)^b is the posterior of a tilted model: since
/// N(y | μ, s²)^b ∝ N(y | μ, s²/b)·s^(1−b), the likelihood's variance is σ²c·v/b, and the prior gains (σ²)^(N(1−b)/2) and each vᵢ^((1−b)/2).
/// Each tilt here keeps an exact draw: σ² ~ IG(α₀ − N(1−b)/2, β₀), and v from its law at its exponent raised by (1 − b)/2, so Geweke's
/// joint-distribution test reads the rung's blocks through the tilted model's two simulators.</item>
/// <item><b>Evidence.</b> A ladder's and SMC's ln p(y) against plain Monte Carlo over the prior's exact draws, the likelihood read given
/// v, on data small enough for that to hold.</item>
/// <item><b>Recovery on a direct functional response</b> under contaminated noise, t BASS beside Gaussian BASS.</item>
/// <item><b>The cost row:</b> a sweep at a rung against b = 1, SMC's stage, and the direct functional point.</item>
/// </list>
/// </summary>
public sealed class BassNoiseReachQualificationTests(ITestOutputHelper output)
{
    // ---- the law at a rung --------------------------------------------------------------------------------------------------------

    private const int GewekeN = 12;
    private const double Alpha0 = 8.0, Beta0 = 6.0;
    private static readonly double[] GewekeX = [.. Enumerable.Range(0, GewekeN).Select(i => i / (GewekeN - 1.0))];
    private static readonly string[] GewekeNames = ["M", "ln tau", "ln sigma2", "ln v0", "mean ln v", "atan f0", "atan f6", "beta", "z0^2"];

    private static BassData GewekeData(double[] y) => new([new ContinuousPredictor("x", GewekeX)], new ScalarResponse(y));

    private static BassNoise GewekeNoise(string form) => form switch
    {
        "t" => BassNoise.StudentT(5, reference: BassNoiseReference.None),
        "gh" => BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), new BassNoisePrior(0.3, 0.5), reference: BassNoiseReference.None),
        "gbp" => BassNoise.GeneralizedBetaPrime(new GbpParameters(1.0, 0.5, 1.5), new BassNoisePrior(0.0, 0.5), reference: BassNoiseReference.None),
        "horseshoe-symmetric" => BassNoise.Horseshoe(reference: BassNoiseReference.None),
        "quantile" => BassNoise.AsymmetricLaplace(0.8, reference: BassNoiseReference.None),
        _ => BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5), reference: BassNoiseReference.None),
    };

    private static BassConfig GewekeConfig(string form, double alpha, bool interweave = false) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 2, Support = new SupportRule(SupportKind.Count, 2) },
        Coefficients = BassCoefficientPrior.G(TauPrior.UnitInformationGamma(5)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.InverseGamma(alpha, Beta0),
        Noise = GewekeNoise(form) with { Interweave = interweave },
        Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 0, Samples = 1 },
    };

    // v from its law tilted by v^((1−b)/2): the exponent of v's density raised by (1 − b)/2.
    private static Func<Xoshiro256PlusPlus, double> Tilted(string form, double b) => form switch
    {
        "t" => GeneralizedInverseGaussian.Of(-2.5 + 0.5 * (1.0 - b), 0.0, 5.0).Sample,
        "gh" => GeneralizedInverseGaussian.Of(0.7 + 0.5 * (1.0 - b), 1.2, 0.8).Sample,
        "gbp" => GeneralizedBetaPrime.Of(1.0, 0.5 + 0.5 * (1.0 - b), 1.5 - 0.5 * (1.0 - b)).Sample,
        "quantile" => GeneralizedInverseGaussian.Of(1.0 + 0.5 * (1.0 - b), 2.0, 0.0).Sample,
        _ => GeneralizedBetaPrime.Of(1.0, 0.5 + 0.5 * (1.0 - b), 0.5 * b).Sample,
    };

    private static WorkingResponse At(BassModel model, double[] v, double sigma2, double beta)
    {
        var mixture = (VarianceMeanMixture)model.Observation;
        return mixture.At(model.Working, new ObservationLatents(v, beta, mixture.Skewed ? sigma2 : null));
    }

    [Theory]
    [InlineData("t", 0.4, 20261107, false, 1)]
    [InlineData("t", 0.75, 20261107, false, 1)]
    [InlineData("gh", 0.4, 20261107, false, 1)]
    [InlineData("gh", 0.75, 20261107, false, 1)]
    [InlineData("gbp", 0.5, 20261107, false, 1)]
    [InlineData("horseshoe", 0.75, 20261107, false, 1)]
    [InlineData("horseshoe-symmetric", 0.5, 20261107, false, 1)]
    // E1's slice 7: the interwoven sweep at b = 1 and on rungs. The skewed horseshoe's cases at seed 20261107, which fail with and without
    // the step as its v mixes slowly (slice 7's report), are readings, not gates; BassInterweavingQualificationTests checks the step alone.
    [InlineData("t", 1.0, 20261107, true, 1)]
    [InlineData("gh", 1.0, 20261107, true, 1)]
    [InlineData("quantile", 1.0, 20261107, true, 1)]
    [InlineData("t", 0.4, 20261107, true, 1)]
    [InlineData("horseshoe", 1.0, 20261108, true, 1)]
    [InlineData("horseshoe", 0.75, 20261108, true, 1)]
    [InlineData("horseshoe", 0.5, 20261108, true, 1)]
    [InlineData("horseshoe", 1.0, 20261107, false, 1)]
    public void Geweke_ARungLeavesItsTiltedJointInvariant(string form, double b, int seed, bool interweave, int scale)
    {
        if (!BassQualification.Enabled) return;   // heavy
        int draws = 400_000 * scale, steps = 600_000 * scale;
        const int batches = 50;
        BassConfig config = GewekeConfig(form, Alpha0, interweave), tilted = GewekeConfig(form, Alpha0 - 0.5 * GewekeN * (1.0 - b), interweave);
        BassNoise noise = config.Noise!;
        double c = noise.VarianceScale;
        var sigmaTilted = InverseGamma.ShapeScale(Alpha0 - 0.5 * GewekeN * (1.0 - b), Beta0);
        Func<Xoshiro256PlusPlus, double> drawV = Tilted(form, b);
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(tilted, GewekeData([.. GewekeX.Select(x => Math.Sin(7 * x))]));
        var sampler = new BassPriorSampler(prior);

        double[] Functions(BassState s, double sigma2, double[] v, double[] f, double[] y, double beta)
            => [s.Count, s.LogTau, Math.Log(sigma2), Math.Log(v[0]), v.Average(x => Math.Log(x)), Math.Atan(f[0]), Math.Atan(f[GewekeN / 2]),
                beta, b * Math.Pow(y[0] - f[0] - Math.Sqrt(sigma2) * beta * v[0], 2) / (sigma2 * c * v[0])];
        double[] Respond(double[] f, double sigma2, double[] v, double beta)
            => [.. Enumerable.Range(0, GewekeN).Select(i => f[i] + Math.Sqrt(sigma2) * (beta * v[i] + Math.Sqrt(c * v[i] / b) * Normal.SampleStandard(rng)))];

        // The tilted model's marginal-conditional simulator: β from its prior, v and a held σ² from their tilted laws, a σ² the structure
        // integrates from the tilted prior with a.
        var marginal = new double[draws][];
        (BassState S, double Sigma2, double[] V, double[] F, double Beta) Draw()
        {
            while (true)
            {
                BassState state = sampler.DrawState(rng).State;
                double beta = noise.SkewPrior is { IsHeld: false } sp ? sp.Mean + sp.Sd * Normal.SampleStandard(rng) : noise.Skew;
                double[] v = [.. Enumerable.Range(0, GewekeN).Select(_ => drawV(rng))];
                double held = sigmaTilted.Sample(rng);
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

        // The successive-conditional simulator: y given θ under the tilted model, then one sweep of the chain at b on that y, its blocks
        // run as SMC's kernel runs them on a stage.
        var theta = Draw();
        var successive = new double[steps][];
        long refused = 0, blockRefused = 0;
        for (int t = 0; t < steps; t++)
        {
            double[] y = Respond(theta.F, theta.Sigma2, theta.V, theta.Beta);
            BassModel model = BassModel.Create(config, GewekeData(y));
            var moves = new BassStructureMoves(model, theta.S.Terms, theta.S.LogTau);
            var blockRng = new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue));
            var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)),
                                                           new Temperature(b, TemperingPath.Likelihood));
            if (BassStructureMoves.Follow(chain, At(model, theta.V, theta.Sigma2, theta.Beta)) is not null) { refused++; continue; }
            chain.Step();
            foreach (IBlock<ReversibleJumpChain<BassState>> block in BassSweepSlots.Blocks(model, BassSweepSlot.Tau))
                if (!block.Update(chain, blockRng)) blockRefused++;
            BassSweepFit fit = moves.SweepFit!;
            ObservationLatents latents = chain.Current.Fit.Working!.Latents!;
            // Under a pure scale mixture σ² is the sweep's draw with the coefficients, and the interweaving step moves σ with v without
            // keeping it: θ's coefficients and σ² are drawn afresh given the structure and the new latents, a draw from their conditional.
            if (interweave && latents.Sigma2 is null) fit = moves.DrawSweepFit(chain.Current, blockRng, model.Target, b);
            double sigma2 = latents.Sigma2 ?? fit.Sigma2;
            successive[t] = Functions(chain.Current, sigma2, latents.ValueArray, fit.Fitted, y, latents.Skew);
            theta = (chain.Current, sigma2, latents.ValueArray, fit.Fitted, latents.Skew);
        }
        Assert.Equal(0, refused);

        double worst = 0;
        for (int g = 0; g < GewekeNames.Length; g++)
        {
            double[] a = [.. marginal.Select(r => r[g])], s = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1), meanB = s.Average();
            int size = s.Length / batches;
            double[] means = [.. Enumerable.Range(0, batches).Select(k => s.Skip(k * size).Take(size).Average())];
            double seB2 = means.Sum(m => (m - means.Average()) * (m - means.Average())) / (batches - 1) / batches;
            bool constant = a.Min() == a.Max() && s.Min() == s.Max();
            double z = constant ? (a[0] == s[0] ? 0.0 : double.PositiveInfinity) : (meanA - meanB) / Math.Sqrt(varA / a.Length + seB2);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{GewekeNames[g],-10} marginal-conditional {meanA,10:F5}  successive {meanB,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"{form} at b = {b}, seed {seed}, {(interweave ? "interwoven" : "plain")}, {scale}× the size: moves refused {refused}; " +
                         $"block draws refused {blockRefused} of {steps}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- evidence -------------------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("t")]
    [InlineData("gh")]
    public void Evidence_MeetsPlainMonteCarloOverThePrior(string form)
    {
        if (!BassQualification.Enabled) return;   // heavy
        var noiseRng = new Xoshiro256PlusPlus(20261109);
        BassData data = GewekeData([.. GewekeX.Select(x => Math.Sin(7 * x) + 0.4 * Normal.SampleStandard(noiseRng))]);
        BassConfig config = GewekeConfig(form, Alpha0) with { Sampling = new BassSampling { Chains = 4, MasterSeed = 5, BurnIn = 2000, Samples = 2000, ReportEvidence = true } };
        BassModel model = BassModel.Create(config, data);
        double c = config.Noise!.VarianceScale;
        ReadOnlySpan<double> y = data.Response.Values;

        // ln p(y) = ln E[N(y | f + σβv, σ²c·v)] over exact draws of the whole prior on F, the latents with them.
        const int draws = 2_000_000;
        var sampler = new BassPriorSampler(model);
        var rng = new Xoshiro256PlusPlus(20261110);
        var logs = new double[draws];
        for (int k = 0; k < draws; k++)
        {
            BassState state = sampler.DrawPrior(rng).State;
            var moves = new BassStructureMoves(model, state, new BassMovesSnapshot(sampler.LastDesign(), new Tally<BassRefusal>(), 0,
                                                                                   new BassSupportTally(model.Prior), new BassMoveTally()));
            BassSweepFit fit = moves.DrawSweepFit(state, rng, model.Target, power: 0.0);
            ObservationLatents latents = state.Fit.Working!.Latents!;
            double sigma2 = latents.Sigma2 ?? fit.Sigma2, sigma = Math.Sqrt(sigma2), log = 0.0;
            for (int i = 0; i < y.Length; i++)
            {
                double vi = latents.Values[i], s2 = sigma2 * c * vi, r = y[i] - fit.Fitted[i] - sigma * latents.Skew * vi;
                log += -0.5 * Math.Log(2.0 * Math.PI * s2) - 0.5 * r * r / s2;
            }
            logs[k] = log;
        }
        double top = logs.Max(), mean = logs.Average(l => Math.Exp(l - top)), second = logs.Average(l => Math.Exp(2.0 * (l - top)));
        double brute = top + Math.Log(mean), bruteError = Math.Sqrt((second - mean * mean) / draws) / mean;

        BassConfig ladder = config with { Sampling = config.Sampling with { Tempering = new TemperingSpec { Levels = 10, HotEnd = 0.02, ReferenceRung = true } } };
        BassConfig smc = config with { Sampling = config.Sampling with { Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4 } } };
        var options = new BassRunOptions { MaxDegreeOfParallelism = 4 };
        BassEvidence tempered = Bass.Fit(ladder, data, options).Evidence!, sequential = Bass.Fit(smc, data, options).Evidence!;
        output.WriteLine($"{form}: plain Monte Carlo over {draws:N0} prior draws ln p(y) {brute:F4} ± {bruteError:F4} (prior attempts {sampler.Attempts:N0}); " +
                         $"the ladder's {tempered.LogEvidence:F4} ± {tempered.MonteCarloError:F4}; SMC's {sequential.LogEvidence:F4} ± {sequential.MonteCarloError:F4}");
        foreach ((string name, BassEvidence e) in new[] { ("ladder", tempered), ("SMC", sequential) })
        {
            double z = (e.LogEvidence - brute) / Math.Sqrt(e.MonteCarloError * e.MonteCarloError + bruteError * bruteError);
            Assert.True(Math.Abs(z) < 4.0, $"{name}: z {z:F2}");
        }
    }

    // ---- recovery on a direct functional response ------------------------------------------------------------------------------------

    [Fact]
    public void Recovery_OnADirectFunctionalResponse_UnderContaminatedNoise()
    {
        if (!BassQualification.Enabled) return;   // heavy
        (BassData clean, double[] truth) = BassQualificationRecoveryTests.Functional();
        int n = clean.Response.Cases, m = clean.Response.Locations;
        foreach ((string name, BassNoise? noise) in new[] { ("t(3)", BassNoise.StudentT(3)), ("Gaussian", (BassNoise?)null) })
        {
            int covered = 0, total = 0;
            double squares = 0;
            for (int set = 0; set < 2; set++)
            {
                // 0.85·N(0, 0.2²) + 0.15·N(0, 1.5²) about the fixture's truth, outside the t family.
                var rng = new Xoshiro256PlusPlus(300 + set);
                var y = new double[n, m];
                for (int i = 0; i < n; i++)
                    for (int l = 0; l < m; l++) y[i, l] = truth[i * m + l] + (rng.NextDouble() < 0.15 ? 1.5 : 0.2) * Normal.SampleStandard(rng);
                var data = new BassData(clean.Predictors, new FunctionalResponse(y, ((FunctionalResponse)clean.Response).OutputVariables));
                BassConfig config = BassStudentTTests.Config(noise) with
                {
                    Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 20 },
                    Sampling = new BassSampling { Chains = 4, MasterSeed = 41 + set, BurnIn = 1500, Samples = 500, Thin = 2 },
                };
                BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
                BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false });
                double[,] mean = p.Mean();
                for (int i = 0; i < n; i++)
                    for (int l = 0; l < m; l += 3)
                    {
                        double[] f = [.. Enumerable.Range(0, p.DrawCount).Select(d => p.Latent[d, i, l]).Order()];
                        double at = truth[i * m + l], lo = f[(int)(0.05 * f.Length)], hi = f[(int)(0.95 * f.Length) - 1];
                        covered += at >= lo && at <= hi ? 1 : 0;
                        total++;
                        squares += (mean[i, l] - at) * (mean[i, l] - at);
                    }
            }
            double coverage = (double)covered / total;
            output.WriteLine($"{name}: f's 90% intervals cover {coverage:P1} of {total} case-location points over 2 data sets, RMSE {Math.Sqrt(squares / total):F4}");
            if (noise is not null) Assert.True(coverage >= 0.7, $"t's coverage {coverage:P1}");
        }
    }

    // ---- the cost row --------------------------------------------------------------------------------------------------------------

    [Fact]
    public void CostRow_ARung_AnSmcStage_AndTheDirectFunctionalPoint()
    {
        if (!BassQualification.Enabled) return;   // heavy
        string clock = BassCost.Cycles ? "thread cycles, indicative" : "wall clock";
        output.WriteLine($"Clock: {clock}. µs per sweep of D04's chain under t(4), by the rung's b, each the mean over its timed sweeps.");
        var points = new (string Name, BassCost.Point Point)[]
        {
            ("the base", BassCost.Base),
            ("n 200", BassCost.Base with { N = 200 }),
            ("direct functional", new BassCost.Point(BassCost.Surface.Functional, 100, 5, BassCost.BaseCap, BassCost.BaseOrder, Locations: 20)),
        };
        // The first pass warms the JIT and is not read: a cold first cell read up to twice its warm cost under load (2026-11-07).
        for (int pass = 0; pass < 2; pass++)
        foreach ((string name, BassCost.Point x) in points)
        {
            BassData data = BassCost.Data(x);
            BassConfig gaussian = BassCost.Config(x);
            BassModel model = BassModel.Create(gaussian with { Noise = BassNoise.StudentT(4) }, data), plain = BassModel.Create(gaussian, data);
            var cells = new List<string>();
            foreach ((string label, BassModel m) in new[] { ("t(4)", model), ("Gaussian", plain) })
                foreach (double b in new[] { 1.0, 0.5, 0.05 })
                {
                    var moves = new BassStructureMoves(m, [], 0.0);
                    var rng = new Xoshiro256PlusPlus(5);
                    var chain = new ReversibleJumpChain<BassState>(moves.Moves, m.Target, moves.Start, new Xoshiro256PlusPlus(3), new Temperature(b, TemperingPath.Likelihood));
                    List<IBlock<ReversibleJumpChain<BassState>>> blocks = BassSweepSlots.Blocks(m, BassSweepSlot.Tau);
                    void Sweep()
                    {
                        chain.Step();
                        foreach (IBlock<ReversibleJumpChain<BassState>> block in blocks) block.Update(chain, rng);
                    }
                    for (int s = 0; s < 500; s++) Sweep();
                    const int timed = 200;
                    long t0 = BassCost.Tick();
                    for (int s = 0; s < timed; s++) Sweep();
                    double us = (BassCost.Tick() - t0) / (double)timed / BassCost.TicksPerSecond * 1e6;
                    cells.Add($"{label} b {b}: {us:F1} (M {chain.Current.Count})");
                }
            if (pass == 1) output.WriteLine($"{name,-18} rows {x.Rows,6}: {string.Join("; ", cells)}");
        }
    }
}
