using System;
using System.Collections.Generic;
using System.Diagnostics;
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
/// E21's slice 4 behind <c>BASS_QUALIFICATION</c> = 1: random effects under GBASS's noises, u a latent on every rung and integrated by the
/// structure step at b = 1 (decision 7 as amended on 2026-11-10).
/// <list type="bullet">
/// <item><b>Geweke's joint-distribution test of one sweep:</b> at b = 1 a chain's sweep with the collapsed step, and without it (the
/// reference form B); below 1 a replica's sweep, its own move then the blocks, against the tilted joint whose b-rung that kernel leaves
/// invariant, as a ladder's rung and an SMC stage run it. Under t, the asymmetric Laplace, the skewed Normal-Wald, the skewed horseshoe with
/// its move along β·v, and t interwoven; s drawn, u and s among the functions.</item>
/// <item><b>A ladder and SMC end to end:</b> their evidence against each other and their posterior means of ln s and M against long
/// chains'.</item>
/// <item><b>Mixing:</b> a plain chain's effective sizes per second with the collapsed step against form B without it.</item>
/// <item><b>Recovery under heavy-tailed noise,</b> effects generated outside the model.</item>
/// </list>
/// </summary>
public sealed class RandomEffectsNoiseQualificationTests(ITestOutputHelper output)
{
    private const int GewekeN = 12;
    private const double Alpha0 = 8.0, Beta0 = 6.0;
    private static readonly double[] GewekeX = [.. Enumerable.Range(0, GewekeN).Select(i => i / (GewekeN - 1.0))];
    private static readonly CategoricalPredictor GewekeGrouping =
        new("g", [.. Enumerable.Range(0, GewekeN).Select(i => i % 3)], ["a", "b", "c"]);
    private static readonly string[] GewekeNames =
        ["M", "ln tau", "ln sigma2", "mean ln v", "beta", "ln s", "atan u0", "ln sum u^2", "u energy", "atan f0", "atan f6", "z0^2"];

    private static BassData GewekeData(double[] y) => new([new ContinuousPredictor("x", GewekeX)], new ScalarResponse(y), [GewekeGrouping]);

    private static string Base(string form) => form.Split('+')[0];

    private static BassNoise GewekeNoise(string form)
    {
        BassNoise noise = Base(form) switch
        {
            "nw" => BassNoise.NormalWald(skew: new BassNoisePrior(0.3, 0.5), reference: BassNoiseReference.None),
            "quantile" => BassNoise.AsymmetricLaplace(0.8, reference: BassNoiseReference.None),
            "horseshoe" => BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5), reference: BassNoiseReference.None),
            _ => BassNoise.StudentT(5, reference: BassNoiseReference.None),
        };
        if (Base(form) == "horseshoe") noise = noise with { SkewRescale = form.Contains("+rescale") };
        return noise with { Interweave = form.Contains("+interweave") };
    }

    // v from its law tilted by v^((1−b)/2), as the noise reach harness tilts it.
    private static Func<Xoshiro256PlusPlus, double> Tilted(string form, double b) => Base(form) switch
    {
        "nw" => GeneralizedInverseGaussian.Of(-0.5 + 0.5 * (1.0 - b), 1.0, 1.0).Sample,
        "quantile" => GeneralizedInverseGaussian.Of(1.0 + 0.5 * (1.0 - b), 2.0, 0.0).Sample,
        "horseshoe" => GeneralizedBetaPrime.Of(1.0, 0.5 + 0.5 * (1.0 - b), 0.5 * b).Sample,
        _ => GeneralizedInverseGaussian.Of(-2.5 + 0.5 * (1.0 - b), 0.0, 5.0).Sample,
    };

    private static BassConfig GewekeConfig(string form, double alpha) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 2, Support = new SupportRule(SupportKind.Count, 2) },
        Coefficients = BassCoefficientPrior.G(TauPrior.UnitInformationGamma(5)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.InverseGamma(alpha, Beta0),
        Noise = GewekeNoise(form),
        RandomIntercepts = [new BassRandomIntercept("g", 0.8, ScalePrior: new BassGroupScalePrior())],
        Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 0, Samples = 1 },
    };

    // A reading of `model` at latents (v, β, σ²) carrying u at its scales.
    private static WorkingResponse At(BassModel model, double[] v, double beta, double sigma2, double scale, double[] effects)
    {
        var mixture = (VarianceMeanMixture)model.Observation;
        WorkingResponse reading = mixture.At(model.Working, new ObservationLatents(v, beta, sigma2));
        return mixture.At(reading, new RandomEffectsLatent(model.RandomEffects!.At(0, scale), [effects]));
    }

    private sealed record Theta(BassState State, double Sigma2, double[] V, double Beta, double Scale, double[] Effects, double[] Fitted);

    [Theory]
    [InlineData("t", 1.0, "collapse", 20261110)]
    [InlineData("t", 1.0, "plain", 20261110)]
    [InlineData("t", 0.5, "replica", 20261110)]
    [InlineData("quantile", 1.0, "collapse", 20261110)]
    [InlineData("nw", 1.0, "collapse", 20261110)]
    [InlineData("nw", 0.5, "replica", 20261110)]
    [InlineData("horseshoe+rescale", 1.0, "collapse", 20261110)]
    [InlineData("horseshoe+rescale", 0.75, "replica", 20261110)]
    [InlineData("t+interweave", 1.0, "collapse", 20261110)]
    public void Geweke_TheSweepLeavesItsJointInvariant(string form, double b, string sweep, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int starts = 40_000, sweeps = 10;
        BassConfig config = GewekeConfig(form, Alpha0), tilted = GewekeConfig(form, Alpha0 - 0.5 * GewekeN * (1.0 - b));
        BassNoise noise = config.Noise!;
        double c = noise.VarianceScale;
        var sigmaTilted = InverseGamma.ShapeScale(Alpha0 - 0.5 * GewekeN * (1.0 - b), Beta0);
        var scalePrior = new BassGroupScalePrior();
        Func<Xoshiro256PlusPlus, double> drawV = Tilted(form, b);
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(tilted, GewekeData([.. GewekeX.Select(x => Math.Sin(7 * x))]));
        var sampler = new BassPriorSampler(prior);

        double[] Functions(Theta t, double[] y) =>
        [
            t.State.Count, t.State.LogTau, Math.Log(t.Sigma2), t.V.Average(Math.Log), t.Beta, Math.Log(t.Scale), Math.Atan(t.Effects[0]),
            Math.Log(t.Effects.Sum(u => u * u)), t.Effects.Sum(u => u * u) / (t.Sigma2 * t.Scale * t.Scale),
            Math.Atan(t.Fitted[0]), Math.Atan(t.Fitted[GewekeN / 2]),
            b * Math.Pow(y[0] - t.Fitted[0] - Math.Sqrt(t.Sigma2) * t.Beta * t.V[0], 2) / (t.Sigma2 * c * t.V[0]),
        ];
        double[] Respond(Theta t)
            => [.. Enumerable.Range(0, GewekeN).Select(i => t.Fitted[i] + Math.Sqrt(t.Sigma2) * (t.Beta * t.V[i] + Math.Sqrt(c * t.V[i] / b) * Normal.SampleStandard(rng)))];

        // The tilted model's marginal-conditional simulator: β, s and u from their priors, v and σ² from their tilted laws, a from its prior.
        Theta Draw()
        {
            while (true)
            {
                BassState state = sampler.DrawState(rng).State;
                double beta = noise.SkewPrior is { IsHeld: false } sp ? sp.Mean + sp.Sd * Normal.SampleStandard(rng) : noise.Skew;
                double[] v = [.. Enumerable.Range(0, GewekeN).Select(_ => drawV(rng))];
                double sigma2 = sigmaTilted.Sample(rng), scale = scalePrior.Sample(rng);
                double[] effects = RandomEffectsLatent.DrawPrior(prior.RandomEffects!.At(0, scale), sigma2, rng).Effects[0];
                var chain = new BassChain(prior, state.Terms, state.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
                if (chain.Follow(At(prior, v, beta, sigma2, scale, effects)) is not null) continue;
                BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, prior.Target, power: 0.0);
                return new Theta(chain.Current, sigma2, v, beta, scale, effects, fit.Fitted);
            }
        }
        // Exact starts (Geweke 2004's test with independent restarts): (θ, y) from the tilted joint, then `sweeps` sweeps at b on that y.
        // θ after them is again a draw of θ given y, so each statistic's change has mean 0, and independent starts make its z exact whatever
        // the chain's mixing. A successive-conditional chain read the half-Cauchy s's funnel as bias: with u a latent, s's excursions to
        // large s, where the data pin u, outlast any batch. u0 and Σu² have no variance under that s, so they are read bounded or in logs.
        var before = new double[starts][];
        var change = new double[starts][];
        for (int r = 0; r < starts; r++)
        {
            Theta theta = Draw();
            double[] y = Respond(theta);
            before[r] = Functions(theta, y);
            BassModel model = BassModel.Create(config, GewekeData(y));
            var moves = new BassStructureMoves(model, theta.State.Terms, theta.State.LogTau);
            var blockRng = new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue));
            var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)),
                                                           new Temperature(b, TemperingPath.Likelihood));
            Assert.Null(BassStructureMoves.Follow(chain, At(model, theta.V, theta.Beta, theta.Sigma2, theta.Scale, theta.Effects)));
            List<IBlock<ReversibleJumpChain<BassState>>> blocks = sweep == "replica"
                ? BassSweepSlots.Blocks(model, BassSweepSlot.Tau, collapse: true)
                : BassSweepSlots.Blocks(model, BassSweepSlot.Structure, collapse: sweep == "collapse");
            for (int k = 0; k < sweeps; k++)
            {
                if (sweep == "replica") chain.Step();
                foreach (IBlock<ReversibleJumpChain<BassState>> block in blocks) block.Update(chain, blockRng);
            }
            WorkingResponse reading = moves.DesignAt(chain.Current).Working;
            ObservationLatents latents = reading.Latents!;
            RandomEffectsLatent effects = reading.Effects!;
            var after = new Theta(chain.Current, latents.Sigma2!.Value, latents.ValueArray, latents.Skew, effects.Block.Scales.Span[0],
                                  effects.Effects[0], moves.SweepFit!.Fitted);
            double[] read = Functions(after, y);
            change[r] = [.. read.Select((x, g) => x - before[r][g])];
        }

        double worst = 0;
        for (int g = 0; g < GewekeNames.Length; g++)
        {
            double[] d = [.. change.Select(r => r[g])];
            double mean = d.Average(), variance = d.Sum(x => (x - mean) * (x - mean)) / (d.Length - 1);
            double z = variance == 0.0 ? (mean == 0.0 ? 0.0 : double.PositiveInfinity) : mean / Math.Sqrt(variance / d.Length);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{GewekeNames[g],-10} start {before.Average(r => r[g]),10:F5}  after {sweeps} sweeps {before.Average(r => r[g]) + mean,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"{form} at b = {b}, {sweep}, seed {seed}: worst |z| {worst:F2}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- a ladder and SMC end to end, the mixing measurement, recovery ----

    // f = sin 2πx₀ + x₁ plus group effects 0.8·(E − 1) sorted onto twelve unbalanced levels by size, plus scaled t₃ noise.
    private static (BassData Data, double[] Effects, double[] Truth) Grouped(int n, int seed, double noise = 0.3)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x0 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())], x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        CategoricalPredictor grouping = RandomInterceptTests.Grouping(n, 12, seed + 1, "g");
        int[] codes = grouping.Codes.ToArray();
        double[] counts = [.. Enumerable.Range(0, 12).Select(l => (double)codes.Count(c => c == l))];
        double[] raw = [.. Enumerable.Range(0, 12).Select(_ => 0.8 * (-Math.Log(1.0 - rng.NextDouble()) - 1.0)).Order()];
        int[] bySize = [.. Enumerable.Range(0, 12).OrderBy(l => counts[l])];
        var effects = new double[12];
        for (int k = 0; k < 12; k++) effects[bySize[k]] = raw[k];
        double centre = Enumerable.Range(0, 12).Sum(l => counts[l] * effects[l]) / n;
        for (int l = 0; l < 12; l++) effects[l] -= centre;   // the convention's, so the effects are the block's to recover
        var t3 = StudentT.Standard(3.0);
        double[] truth = [.. Enumerable.Range(0, n).Select(i => Math.Sin(2.0 * Math.PI * x0[i]) + x1[i])];
        double[] y = [.. Enumerable.Range(0, n).Select(i => truth[i] + effects[codes[i]] + noise * t3.Sample(rng))];
        return (new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y), [grouping]), effects, truth);
    }

    private static BassConfig Config(BassSampling sampling) => BassStudentTTests.Config(BassNoise.StudentT(4)) with
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 6 },
        Sigma2 = Sigma2Prior.InverseGamma(3.0, 0.5),
        RandomIntercepts = [new BassRandomIntercept("g", 1.0, ScalePrior: new BassGroupScalePrior())],
        Sampling = sampling,
    };

    [Fact]
    public void ALadderAndSmc_AgreeWithEachOtherAndWithChains()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassData data = Grouped(120, 31).Data;
        BassFit chains = Bass.Fit(Config(new BassSampling { Chains = 4, MasterSeed = 3, BurnIn = 3000, Samples = 8000 }), data,
                                  new BassRunOptions { MaxDegreeOfParallelism = 4 });
        BassFit ladder = Bass.Fit(Config(new BassSampling
        {
            Chains = 4, MasterSeed = 5, BurnIn = 2000, Samples = 4000, ReportEvidence = true,
            Tempering = new TemperingSpec { Levels = 6, HotEnd = 0.05, ReferenceRung = true },
        }), data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        BassFit smc = Bass.Fit(Config(new BassSampling
        {
            Chains = 1, MasterSeed = 7, ReportEvidence = true, Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4 },
        }), data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        (double Mean, double Error) Summary(BassFit fit, Func<BassDraw, double> read)
        {
            double[] means = [.. fit.Draws.Select(chain => chain.Average(read))];
            double[]? w = fit.DrawWeights;
            if (w is not null)
            {
                BassDraw[] all = [.. fit.Draws.SelectMany(c => c)];
                double mean = all.Select((d, i) => w[i] * read(d)).Sum() / w.Sum();
                return (mean, Math.Sqrt(means.Select(m => (m - means.Average()) * (m - means.Average())).Sum() / Math.Max(1, means.Length - 1) / means.Length));
            }
            return (means.Average(), Math.Sqrt(means.Select(m => (m - means.Average()) * (m - means.Average())).Sum() / (means.Length - 1) / means.Length));
        }
        double LogS(BassDraw d) => Math.Log(d.Groups![0].Scale);
        double M(BassDraw d) => d.Structure.TermIds.Length;
        BassEvidence a = ladder.Evidence!, s = smc.Evidence!;
        output.WriteLine($"ln p(y): ladder {a.LogEvidence:F3} ± {a.MonteCarloError:F3}, SMC {s.LogEvidence:F3} ± {s.MonteCarloError:F3}");
        Assert.True(Math.Abs(a.LogEvidence - s.LogEvidence) <= 4.0 * Math.Sqrt(a.MonteCarloError * a.MonteCarloError + s.MonteCarloError * s.MonteCarloError) + 0.05);
        if (smc.Diagnostics.Smc is { } readout) Assert.Empty(readout.Flags);
        foreach ((string name, Func<BassDraw, double> read) in new (string, Func<BassDraw, double>)[] { ("ln s", LogS), ("M", M) })
        {
            (double cm, double ce) = Summary(chains, read);
            foreach ((string sampler, BassFit fit) in new[] { ("ladder", ladder), ("SMC", smc) })
            {
                (double m, double e) = Summary(fit, read);
                output.WriteLine($"{name}: {sampler} {m:F4} ± {e:F4} against chains' {cm:F4} ± {ce:F4}");
                Assert.True(Math.Abs(m - cm) <= 5.0 * Math.Sqrt(e * e + ce * ce) + 0.03, $"{name}, {sampler}: {m:F4} against {cm:F4}");
            }
        }
    }

    [Theory]
    [InlineData(11)]
    [InlineData(23)]
    public void Mixing_TheCollapsedStep_AgainstTheLikelihoodGivenU(int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassData data = Grouped(300, 41).Data;
        BassModel model = BassModel.Create(Config(new BassSampling { Chains = 1, MasterSeed = seed, BurnIn = 0, Samples = 1 }), data);
        string[] names = ["ln sigma2", "ln s", "u0", "M"];
        var perSecond = new double[2][];
        var perSweep = new double[2][];
        foreach (bool collapse in new[] { true, false })
        {
            // M's integrated autocorrelation is about a thousand sweeps here, so 8000 kept per chain read its ESS to within a factor of 3.
            const int chains = 4, burn = 2000, kept = 40_000;
            var samples = new double[names.Length][][];
            for (int g = 0; g < names.Length; g++) samples[g] = new double[chains][];
            long t0 = BassCost.Tick();
            for (int c = 0; c < chains; c++)
            {
                var moves = new BassStructureMoves(model, [], 0.0);
                var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(seed * 10 + c));
                var rng = new Xoshiro256PlusPlus(seed * 100 + c);
                List<IBlock<ReversibleJumpChain<BassState>>> blocks = BassSweepSlots.Blocks(model, BassSweepSlot.Structure, collapse);
                for (int g = 0; g < names.Length; g++) samples[g][c] = new double[kept];
                for (int k = 0; k < burn + kept; k++)
                {
                    foreach (IBlock<ReversibleJumpChain<BassState>> block in blocks) block.Update(chain, rng);
                    if (k < burn) continue;
                    WorkingResponse reading = moves.DesignAt(chain.Current).Working;
                    samples[0][c][k - burn] = Math.Log(reading.Latents!.Sigma2!.Value);
                    samples[1][c][k - burn] = Math.Log(reading.Effects!.Block.Scales.Span[0]);
                    samples[2][c][k - burn] = reading.Effects.Effects[0][0];
                    samples[3][c][k - burn] = chain.Current.Count;
                }
            }
            double seconds = (BassCost.Tick() - t0) / BassCost.TicksPerSecond;
            perSecond[collapse ? 0 : 1] = new double[names.Length];
            perSweep[collapse ? 0 : 1] = new double[names.Length];
            var cells = new List<string>();
            for (int g = 0; g < names.Length; g++)
            {
                ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(samples[g]);
                perSecond[collapse ? 0 : 1][g] = d.BulkEss.Value / seconds;
                perSweep[collapse ? 0 : 1][g] = d.BulkEss.Value;
                double[] means = [.. samples[g].Select(chain => chain.Average())];
                double mean = means.Average(), se = Math.Sqrt(means.Sum(m => (m - mean) * (m - mean)) / (chains - 1) / chains);
                cells.Add($"{names[g]} {mean:F4} ± {se:F4} R̂ {d.RHat.Value:F3} ESS {d.BulkEss.Value:F0} ({d.BulkEss.Value / seconds:F1}/s)");
            }
            output.WriteLine($"seed {seed}, {(collapse ? "collapsed at b = 1" : "the likelihood given u")}: {seconds:F1} s; {string.Join("; ", cells)}");
        }
        output.WriteLine($"seed {seed}: bulk ESS, collapsed over not, per sweep: " +
                         string.Join(", ", names.Select((n, g) => $"{n} {perSweep[0][g] / perSweep[1][g]:F2}")) + "; per second: " +
                         string.Join(", ", names.Select((n, g) => $"{n} {perSecond[0][g] / perSecond[1][g]:F2}")));
    }

    [Fact]
    public void Recovery_UnderHeavyTailedNoise_EffectsGeneratedOutsideTheModel()
    {
        if (!BassQualification.Enabled) return;   // heavy
        int covered = 0, total = 0, held = 0, better = 0;
        double standardized = 0.0;
        for (int set = 0; set < 4; set++)
        {
            (BassData data, double[] effects, double[] truth) = Grouped(300, 101 + 7 * set);
            BassFit fit = Bass.Fit(Config(new BassSampling { Chains = 2, MasterSeed = 9 + set, BurnIn = 2000, Samples = 3000 }), data,
                                   new BassRunOptions { MaxDegreeOfParallelism = 2 });
            BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
            for (int l = 0; l < effects.Length; l++)
            {
                double[] u = [.. draws.Select(d => d.Groups![0].Effects[l]).Order()];
                double lo = u[(int)(0.05 * u.Length)], hi = u[(int)(0.95 * u.Length)], mean = u.Average();
                double sd = Math.Sqrt(u.Sum(x => (x - mean) * (x - mean)) / (u.Length - 1));
                covered += effects[l] >= lo && effects[l] <= hi ? 1 : 0;
                standardized += Math.Pow((mean - effects[l]) / sd, 2);
                total++;
            }
            // The fitted level and surface at the training cases: f + u_g against the truth, with and without the block.
            BassFit plain = Bass.Fit(Config(new BassSampling { Chains = 2, MasterSeed = 9 + set, BurnIn = 2000, Samples = 3000 }) with { RandomIntercepts = null },
                                     new BassData(data.Predictors, data.Response), new BassRunOptions { MaxDegreeOfParallelism = 2 });
            int[] codes = data.Groupings[0].Codes.ToArray();
            // The posterior mean of f + Σ u at the training cases (BassDrawFit), against the truth's f + u.
            double Rmse(BassFit f, BassData on)
            {
                var reader = new BassDrawFit(BassModel.Create(f.Config, on), f.Terms);
                var mean = new double[on.Response.Cases];
                var one = new double[on.Response.Cases];
                BassDraw[] all = [.. f.Draws.SelectMany(c => c)];
                foreach (BassDraw d in all)
                {
                    reader.Evaluate(d, one);
                    for (int i = 0; i < mean.Length; i++) mean[i] += one[i] / all.Length;
                }
                return Math.Sqrt(Enumerable.Range(0, mean.Length).Average(i => Math.Pow(mean[i] - truth[i] - effects[codes[i]], 2)));
            }
            double rmseBlock = Rmse(fit, data), rmsePlain = Rmse(plain, new BassData(data.Predictors, data.Response));
            output.WriteLine($"set {set}: f + u RMSE {rmseBlock:F3} with the block against {rmsePlain:F3} without");
            held++;
            better += rmseBlock < rmsePlain ? 1 : 0;
        }
        double msse = standardized / total;
        output.WriteLine($"effects' 90% intervals cover {covered} of {total}; standardized errors' mean square {msse:F2}; the block fits better in {better} of {held}");
        Assert.True(covered >= 0.78 * total, $"{covered} of {total}");
        Assert.InRange(msse, 0.5, 2.0);
        Assert.Equal(held, better);
    }
}
