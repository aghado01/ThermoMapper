using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's SMC particles and kernel (stage 8, C4's step 2): a fork is an exact, independent copy of its parent, the
/// parent's design and later moves untouched by the fork's; and a tempered run gives the same population, normalizer,
/// stage reports and account whether its particles move one at a time or in parallel. Each stage's law against the
/// enumeration oracle is <c>BassSmcParityTests</c>'.
/// </summary>
public sealed class BassSmcTests
{
    private readonly ITestOutputHelper _out;
    public BassSmcTests(ITestOutputHelper output) => _out = output;

    private static BassModel Model() => BassModel.Create(BassSweepTests.Config(maxBasis: 4), BassSweepTests.Data(large: true));

    private static SmcStage<BassParticle> Stage(BassModel model, double b) => new(new BassParticleTarget(model.Target), new Temperature(b, TemperingPath.Likelihood));

    // A particle's design as it stands: its terms, its statistics and its coefficients.
    private static (BassTerm[] Terms, double Residual, double Fitted, double[] Beta, double LogTau) Design(BassParticle p)
    {
        BassDesign d = p.Moves.DesignAt(p.State);
        return ([.. d.Terms], d.ResidualSquares, d.FittedSquares, d.Coefficients.ToArray(), p.State.LogTau);
    }

    private static List<(BassTerm[] Terms, double LogTau)> Walk(BassSmcKernel kernel, BassParticle p, SmcStage<BassParticle> stage, int seed, int sweeps)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var path = new List<(BassTerm[], double)>();
        ITargetEvaluation e = stage.Target.Evaluate(p);
        for (int s = 0; s < sweeps; s++)
        {
            e = kernel.Move(p, e, stage, rng).Evaluation;
            path.Add(([.. p.State.Terms], p.State.LogTau));
        }
        return path;
    }

    private static void SamePath(List<(BassTerm[] Terms, double LogTau)> a, List<(BassTerm[] Terms, double LogTau)> b)
    {
        Assert.Equal(a.Count, b.Count);
        for (int s = 0; s < a.Count; s++)
        {
            Assert.Equal(a[s].Terms, b[s].Terms);
            Assert.Equal(a[s].LogTau, b[s].LogTau);
        }
    }

    [Fact]
    public void Fork_IsAnExactCopy_AndIndependentOfItsParent()
    {
        BassModel model = Model();
        var particles = new BassParticles(model);
        var kernel = new BassSmcKernel(model);
        SmcStage<BassParticle> stage = Stage(model, 0.6);
        BassParticle Parent()
        {
            BassParticle p = particles.Draw(new Xoshiro256PlusPlus(41));
            Walk(kernel, p, stage, seed: 42, sweeps: 60);
            return p;
        }

        // A fork moves exactly as its parent would from the same stream: the design is copied bit for bit.
        BassParticle parent = Parent(), twin = Parent();
        BassParticle fork = particles.Fork(parent);
        Assert.NotSame(parent.State, fork.State);
        Assert.NotSame(parent.Moves, fork.Moves);
        Assert.False(ReferenceEquals(parent.State.Terms, fork.State.Terms));
        var before = Design(parent);
        List<(BassTerm[] Terms, double LogTau)> forked = Walk(kernel, fork, stage, seed: 43, sweeps: 200);
        Assert.True(forked.Select(s => s.Terms.Length).Distinct().Count() > 1, "the fork's walk never changed its structure");

        // The fork's walk left the parent as it was, and the parent then moves as its twin, which was never forked.
        var after = Design(parent);
        Assert.Equal(before.Terms, after.Terms);
        Assert.Equal(before.Residual, after.Residual);
        Assert.Equal(before.Fitted, after.Fitted);
        Assert.Equal(before.Beta, after.Beta);
        Assert.Equal(before.LogTau, after.LogTau);
        List<(BassTerm[] Terms, double LogTau)> parentPath = Walk(kernel, parent, stage, seed: 43, sweeps: 200);
        SamePath(parentPath, forked);
        SamePath(Walk(kernel, twin, stage, seed: 43, sweeps: 200), parentPath);
    }

    [Fact]
    public void Run_InParallel_MatchesTheSerialRun()
    {
        BassModel model = Model();
        (SmcResult<BassParticle> Result, BassSmcAccount Account, long Attempts) Run(int parallelism)
        {
            var particles = new BassParticles(model);
            var kernel = new BassSmcKernel(model);
            var options = new SmcOptions { Particles = 96, MutationSteps = 2, MaxDegreeOfParallelism = parallelism };
            SmcResult<BassParticle> result = SmcSampler.RunTempered(new BassParticleTarget(model.Target), TemperingPath.Likelihood, 0.9, particles, kernel,
                                                                     options, seed: 2026);
            return (result, kernel.Account, particles.PriorAttempts);
        }
        var serial = Run(1);
        var parallel = Run(6);
        _out.WriteLine($"{serial.Result.Stages.Count - 1} stages, ln Ẑ {serial.Result.LogNormalizerRatio:F4} ± {serial.Result.LogNormalizerError:F4}, " +
                       $"acceptance {serial.Account.AcceptanceRate:F3}, prior attempts {serial.Attempts}");
        Assert.Equal(1.0, serial.Result.Stages[^1].Temperature.InverseTemperature);
        Assert.Equal(serial.Result.LogNormalizerRatio, parallel.Result.LogNormalizerRatio);
        Assert.Equal(serial.Result.Weights, parallel.Result.Weights);
        Assert.Equal(serial.Result.Stages, parallel.Result.Stages);
        for (int i = 0; i < serial.Result.Particles.Count; i++)
        {
            Assert.Equal(serial.Result.Particles[i].State.Terms, parallel.Result.Particles[i].State.Terms);
            Assert.Equal(serial.Result.Particles[i].State.LogTau, parallel.Result.Particles[i].State.LogTau);
        }
        Assert.Equal(serial.Attempts, parallel.Attempts);
        Assert.Equal(serial.Account.Moves, parallel.Account.Moves);
        Assert.Equal(serial.Account.Kinds.Entries, parallel.Account.Kinds.Entries);
        Assert.Equal(serial.Account.Refusals.Entries, parallel.Account.Refusals.Entries);
        Assert.Equal(serial.Account.Support.Entries, parallel.Account.Support.Entries);
        Assert.Equal(serial.Account.EmptyBirths, parallel.Account.EmptyBirths);
        // Every sweep is counted: the stages' moves over the particles of positive weight.
        Assert.Equal((serial.Result.Stages.Count - 1) * 96L * 2, serial.Account.Sweeps);
        Assert.Equal(serial.Account.Moves.Sum(m => m.Attempts), serial.Account.Kinds.Entries.Sum(k => k.Attempts));
    }

    // ── The fit (step 3) ──────────────────────────────────────────────────────────────────────────────────────────

    private static readonly SmcSpec Small = new() { Particles = 150, Replicates = 2, MutationSteps = 2, TargetConditionalEss = 0.9 };

    private static BassConfig SmcConfig(SmcSpec? smc = null, BassConfig? config = null)
    {
        BassConfig c = config ?? BassSweepTests.Config(maxBasis: 4);
        return c with { Sampling = c.Sampling with { MasterSeed = 7, Smc = smc ?? Small } };
    }

    // Every draw's term ids by value, τ, σ², λ and coefficients, replicate by replicate.
    private static void SameDraws(BassFit a, BassFit b)
    {
        Assert.Equal(a.Draws.Count, b.Draws.Count);
        for (int c = 0; c < a.Draws.Count; c++)
        {
            Assert.Equal(a.Draws[c].Count, b.Draws[c].Count);
            for (int i = 0; i < a.Draws[c].Count; i++)
            {
                BassDraw x = a.Draws[c][i], y = b.Draws[c][i];
                Assert.Equal(x.Structure.TermIds.Select(id => a.Terms[id]), y.Structure.TermIds.Select(id => b.Terms[id]));
                Assert.Equal((x.Tau, x.Sigma2, x.Lambda), (y.Tau, y.Sigma2, y.Lambda));
                Assert.Equal(x.Coefficients, y.Coefficients);
            }
        }
        Assert.Equal(a.DrawWeights, b.DrawWeights);
    }

    [Fact]
    public void Presentation_ReportsTheSpec_AndRefusesTraces()
    {
        BassFit fit = Bass.Fit(SmcConfig(), BassSweepTests.Data(large: true));
        BassRunLengths lengths = fit.Metadata().Lengths;
        Assert.Equal(Small.Replicates, lengths.Chains);
        Assert.Null(lengths.BurnIn);
        Assert.Null(lengths.Samples);
        Assert.Null(lengths.Thin);
        Assert.Equal(Small, lengths.Smc);
        Assert.All(lengths.DrawsPerChain, n => Assert.Equal(Small.Particles, n));
        Assert.Equal(Small.Replicates * Small.Particles, lengths.Draws);
        Assert.Same(lengths.Smc, fit.Summary().Metadata.Lengths.Smc);
        Assert.Throws<InvalidOperationException>(() => fit.Traces());
    }

    [Fact]
    public void Fit_KeepsEachReplicatesPopulation_WithWeightsThatReachTheViews()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig config = SmcConfig();
        BassFit fit = Bass.Fit(config, data);
        int r = Small.Replicates, n = Small.Particles;
        Assert.Equal(new RunIdentity(7, r, Bass.Fingerprint(config, data)), fit.Identity);
        Assert.Equal(r, fit.Draws.Count);
        Assert.All(fit.Draws, store => Assert.Equal(n, store.Count));
        double[] w = fit.DrawWeights!;
        Assert.Equal(r * n, w.Length);
        Assert.Equal(1.0, w.Sum(), 12);
        for (int c = 0; c < r; c++) Assert.Equal(1.0 / r, w.Skip(c * n).Take(n).Sum(), 12);

        // The population's account; the chain readings do not apply.
        BassSmcReadout smc = fit.Diagnostics.Smc!;
        Assert.True(smc.ReachedPosterior);
        Assert.Equal(r, smc.Replicates.Count);
        foreach (BassSmcReplicate rep in smc.Replicates)
        {
            Assert.Equal(0.0, rep.Stages[0].Temperature.InverseTemperature);
            Assert.Equal(1.0, rep.Stages[^1].Temperature.InverseTemperature);
            Assert.True(rep.StageCount >= 1);
            Assert.InRange(rep.Families, 1, n);
            Assert.InRange(rep.EffectiveSampleSize, 1.0, n);
            Assert.True(rep.PriorAttempts >= n);
            Assert.Equal(rep.StageCount * (long)n * Small.MutationSteps, rep.Account.Sweeps);
            Assert.True(double.IsFinite(rep.LogNormalizer) && rep.LogNormalizerError > 0.0);
            _out.WriteLine($"replicate {rep.Seed}: {rep.StageCount} stages, ln Ẑ {rep.LogNormalizer:F3} ± {rep.LogNormalizerError:F3}, ESS {rep.EffectiveSampleSize:F1}, " +
                           $"families {rep.Families}, acceptance {rep.Account.AcceptanceRate:F3}");
        }
        Assert.Equal(["M", "tau", "log_posterior", "sigma2", .. fit.Resolved.MonitoredPoints.Select(p => $"fitted[{p.Case}]")], smc.Spread.Select(q => q.Name));
        Assert.All(smc.Spread, q => Assert.True(q.Spread >= 0.0 && q.PosteriorSd >= 0.0));
        BassDiagnostics d = fit.Diagnostics;
        Assert.Equal(DiagnosticState.Unavailable, d.Count.RHat.State);
        Assert.Equal(ConvergenceOutcome.Undetermined, d.Verdict.Outcome);
        Assert.Empty(d.Chains);
        Assert.Null(d.Tempering);
        Assert.Equal(n, d.SamplesUsed);
        Assert.Equal(d.Moves.Sum(m => m.Attempts), smc.Replicates.Sum(x => x.Account.Moves.Sum(m => m.Attempts)));
        double[] means = [.. smc.Quantity("log_posterior").ReplicateMeans];
        Assert.Equal(means.Max() - means.Min(), d.ModeSeparation.LogPosteriorGap, 12);

        // The views read the weights: the training fit's mean and σ²'s summary against the weighted draws by hand, and
        // sensitivity's decomposition carries them.
        BassModel model = BassModel.Create(config, data);
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        double[][] fitted = BassQualification.Fitted(model, fit);
        BassTrainingFit training = fit.TrainingFit(data);
        for (int i = 0; i < training.Cases; i++)
            Assert.Equal(Enumerable.Range(0, pooled.Length).Sum(s => w[s] * fitted[s][i]), training.Mean[i, 0], 9);
        Assert.Equal(Enumerable.Range(0, pooled.Length).Sum(s => w[s] * pooled[s].Sigma2), fit.Summary().Posterior.Sigma2.Mean, 12);
        Assert.Equal(w, fit.Sensitivity(data).Decompositions[0].Weights);
        Assert.NotNull(fit.Evidence);
        Assert.Equal(fit.Evidence!.LogEvidence, smc.Replicates.Average(x => x.LogNormalizer), 12);
    }

    [Fact]
    public void PointwiseLogDensity_TakesTheWeights_AgainstAnEqualWeightResample()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassFit fit = Bass.Fit(SmcConfig(Small with { Particles = 400 }), data);
        PointwiseLogDensity density = fit.PointwiseLogDensity(data);
        Assert.Equal(fit.DrawWeights, density.Weights);
        Assert.Null(density.ChainLengths);

        // Equal weights at any scale give the unweighted checks.
        var unweighted = density with { Weights = null };
        foreach (double scale in new[] { 1.0 / density.Draws, 3.0 })
        {
            var equal = new PointwiseLogDensity(density.Values, density.Units, density.Latents) { Weights = [.. Enumerable.Repeat(scale, density.Draws)] };
            Assert.Equal(PredictiveChecks.Waic(unweighted).Elpd, PredictiveChecks.Waic(equal).Elpd, 9);
            Assert.Equal(PredictiveChecks.PsisLoo(unweighted).Elpd, PredictiveChecks.PsisLoo(equal).Elpd, 9);
        }

        // Integer weights are a draw repeated: WAIC's lppd = elpd + p is exact.
        int s = density.Draws, units = density.Units.Units;
        var counts = Enumerable.Range(0, s).Select(i => 1 + i % 3).ToArray();
        var repeated = new double[counts.Sum(), units];
        for (int i = 0, row = 0; i < s; i++)
            for (int k = 0; k < counts[i]; k++, row++)
                for (int u = 0; u < units; u++) repeated[row, u] = density.Values[i, u];
        PredictiveCheck byWeight = PredictiveChecks.Waic(new PointwiseLogDensity(density.Values, density.Units, density.Latents) { Weights = [.. counts.Select(c => (double)c)] });
        PredictiveCheck byRepeat = PredictiveChecks.Waic(new PointwiseLogDensity(repeated, density.Units, density.Latents));
        Assert.Equal(byRepeat.Elpd + byRepeat.EffectiveParameters, byWeight.Elpd + byWeight.EffectiveParameters, 9);

        // The weighted checks against an equal-weight multinomial resample of 40,000 draws.
        var rng = new Xoshiro256PlusPlus(99);
        double[] w = fit.DrawWeights!;
        var cdf = new double[s];
        for (int i = 0; i < s; i++) cdf[i] = (i > 0 ? cdf[i - 1] : 0.0) + w[i];
        const int resampled = 40_000;
        var values = new double[resampled, units];
        for (int k = 0; k < resampled; k++)
        {
            int j = Array.BinarySearch(cdf, rng.NextDouble() * cdf[^1]);
            int i = Math.Min(s - 1, j < 0 ? ~j : j);
            for (int u = 0; u < units; u++) values[k, u] = density.Values[i, u];
        }
        var resample = new PointwiseLogDensity(values, density.Units, density.Latents);
        PredictiveCheck waic = PredictiveChecks.Waic(density), waicR = PredictiveChecks.Waic(resample);
        PredictiveCheck loo = PredictiveChecks.PsisLoo(density), looR = PredictiveChecks.PsisLoo(resample);
        _out.WriteLine($"WAIC {waic.Elpd:F3} (p {waic.EffectiveParameters:F3}) against the resample's {waicR.Elpd:F3} (p {waicR.EffectiveParameters:F3}); " +
                       $"PSIS-LOO {loo.Elpd:F3} ± {loo.MonteCarloError:F3} against {looR.Elpd:F3} ± {looR.MonteCarloError:F3}");
        Assert.Equal(waicR.Elpd, waic.Elpd, 0.1);
        Assert.Equal(looR.Elpd, loo.Elpd, 0.1);
        Assert.Throws<ArgumentException>(() => density with { ChainLengths = null, Weights = [1.0] });
        Assert.Throws<ArgumentException>(() => new PointwiseLogDensity(density.Values, density.Units, density.Latents, [s]) { Weights = w });
    }

    [Theory]
    [InlineData("fixed g")]
    [InlineData("functional")]
    [InlineData("mixed")]
    [InlineData("jeffreys, linear")]
    [InlineData("flat, linear")]
    [InlineData("linear, traced")]
    [InlineData("one replicate")]
    [InlineData("systematic")]
    [InlineData("importance sampling")]
    public void Fit_UnderEachOption_ReachesThePosterior(string option)
    {
        BassConfig baseConfig = BassSweepTests.Config(maxBasis: 4);
        BassData data = option switch
        {
            "functional" => BassKernelEnumeration.FunctionalData(),
            "mixed" => BassKernelEnumeration.CategoricalData(),
            _ => BassSweepTests.Data(large: true),
        };
        SmcSpec spec = Small;
        BassConfig config = baseConfig;
        switch (option)
        {
            case "fixed g": config = config with { Coefficients = BassCoefficientPrior.G(TauPrior.Held()) }; break;
            case "jeffreys, linear": config = config with { Sigma2 = Sigma2Prior.Jeffreys }; spec = spec with { Form = TemperedForm.Linear }; break;
            case "flat, linear": config = config with { Intercept = GPriorCentre.Flat }; spec = spec with { Form = TemperedForm.Linear }; break;
            case "linear, traced": spec = spec with { Form = TemperedForm.Linear, TraceStages = true }; break;
            case "one replicate": spec = spec with { Replicates = 1 }; break;
            case "systematic": spec = spec with { Resampling = ResamplingScheme.Systematic }; break;
            case "importance sampling": spec = spec with { MutationSteps = 0, Particles = 2000 }; break;
        }
        config = SmcConfig(spec, config);
        BassFit fit = Bass.Fit(config, data);
        BassSmcReadout smc = fit.Diagnostics.Smc!;
        Assert.True(smc.ReachedPosterior);
        Assert.Equal(1.0, fit.DrawWeights!.Sum(), 12);
        Assert.Contains(fit.Resolved.AutoPicks, a => a.Setting == "Sampling.Chains" && a.Value == "not read");
        bool proper = !option.StartsWith("jeffreys") && !option.StartsWith("flat");
        Assert.Equal(proper, fit.Evidence is not null);
        Assert.Equal(option == "fixed g", fit.Diagnostics.Tau is null);
        Assert.DoesNotContain(smc.Spread, q => option == "fixed g" && q.Name == "tau");
        if (option == "fixed g") Assert.All(fit.Draws.SelectMany(c => c), dr => Assert.Equal(fit.Resolved.Scale.HeldTau!.Value, dr.Tau));
        if (option == "one replicate")
        {
            Assert.All(smc.Spread, q => Assert.True(double.IsNaN(q.Spread)));
            Assert.True(double.IsNaN(fit.Diagnostics.ModeSeparation.LogPosteriorGap));
            Assert.Equal(fit.Evidence!.SteppingStone.WithinChainError, fit.Evidence.MonteCarloError);
        }
        if (option == "systematic")
        {
            // Only multinomial offspring are what the genealogy error describes: the error is NaN, and so is the evidence's.
            Assert.All(smc.Replicates, x => Assert.True(double.IsNaN(x.LogNormalizerError)));
            Assert.True(double.IsNaN(fit.Evidence!.MonteCarloError));
        }
        if (option == "importance sampling")
        {
            Assert.All(smc.Replicates, x => Assert.Equal(0, x.Account.Sweeps));
            Assert.True(double.IsNaN(fit.Diagnostics.AcceptanceRate));
        }
        if (option == "linear, traced")
        {
            BassEvidence e = fit.Evidence!;
            Assert.NotNull(e.PathSampling);
            Assert.All(e.SteppingStone.Chains, c => Assert.NotEmpty(c.Rungs));
            for (int c = 0; c < smc.Replicates.Count; c++) Assert.Equal(smc.Replicates[c].LogNormalizer, e.SteppingStone.Chains[c].LogRatio, 9);
            _out.WriteLine($"linear: normalizer {e.LogEvidence:F4} ± {e.MonteCarloError:F4}, path sampling {e.PathSampling!.LogEvidence:F4} ± {e.PathSampling.MonteCarloError:F4}");
        }
        fit.Predict(data, PredictionCases.Of(data));
    }

    [Fact]
    public void Fit_StopsOnCancellation_AndReportsProgress()
    {
        BassData data = BassSweepTests.Data(large: true);
        using var cancel = new System.Threading.CancellationTokenSource();
        var reports = new List<BassProgress>();
        var progress = new Recorder(p =>
        {
            lock (reports) reports.Add(p);
            if (p.Sweeps >= 3000) cancel.Cancel();
        });
        Assert.Throws<OperationCanceledException>(() => Bass.Fit(SmcConfig(), data, new BassRunOptions
        {
            Cancellation = cancel.Token, Progress = progress, ProgressInterval = 1000,
        }));
        Assert.NotEmpty(reports);
        Assert.All(reports, p => Assert.Equal((long)Small.Replicates * Small.Particles * Small.MutationSteps * Small.MaxStages, p.PlannedSweeps));
    }

    private sealed class Recorder(Action<BassProgress> report) : IProgress<BassProgress>
    {
        public void Report(BassProgress value) => report(value);
    }

    [Fact]
    public void Fit_RefusesTheSettingsSmcCannotRead()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig config = SmcConfig();
        void Refused(BassConfig c, string field) => Assert.Contains(field, Assert.Throws<ArgumentException>(() => Bass.Fit(c, data)).Message);
        Refused(config with { Sampling = config.Sampling with { Tempering = new TemperingSpec() } }, "Sampling.Tempering");
        Refused(config with { Sampling = config.Sampling with { RHatTarget = 1.01 } }, "Sampling.RHatTarget");
        Refused(SmcConfig(Small with { Particles = 1 }), "Sampling.Smc.Particles");
        Refused(SmcConfig(Small with { Replicates = 0 }), "Sampling.Smc.Replicates");
        Refused(SmcConfig(Small with { TargetConditionalEss = 1.0 }), "Sampling.Smc.TargetConditionalEss");
        Refused(SmcConfig(Small with { MutationSteps = -1 }), "Sampling.Smc.MutationSteps");
        Refused(SmcConfig(Small with { ResampleBelow = 1.5 }), "Sampling.Smc.ResampleBelow");
        Refused(SmcConfig(Small with { MaxStages = 0 }), "Sampling.Smc.MaxStages");
        Refused(SmcConfig(Small, BassSweepTests.Config(maxBasis: 4) with { Sigma2 = Sigma2Prior.Jeffreys }), "Sampling.Smc.Form");
        Refused(SmcConfig(Small, BassSweepTests.Config(maxBasis: 4) with { Intercept = GPriorCentre.Flat }), "Sampling.Smc.Form");
        // The config round-trips with its SMC, and the setting moves the fingerprint.
        Assert.Equal(config, BassConfig.FromJson(config.ToJson()));
        Assert.Contains("\"Smc\":{\"Particles\":150,\"Replicates\":2,", config.ToJson());
        Assert.NotEqual(Bass.Fingerprint(config with { Sampling = config.Sampling with { Smc = null } }, data), Bass.Fingerprint(config, data));
        Assert.Contains("\"Smc\":null", new BassConfig().ToJson());
    }

    [Fact]
    public void Fit_ThatSpendsItsStageBudget_IsFlagged_AndResumesBitForBit()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassFit whole = Bass.Fit(SmcConfig(), data);
        int stages = whole.Diagnostics.Smc!.Replicates.Max(x => x.StageCount);
        Assert.True(stages > 4, $"the run took {stages} stages; the budget test needs more");

        // Short of b = 1: flagged, no evidence or residual check, its posterior views refused, its checkpoint kept.
        BassConfig shortConfig = SmcConfig(Small with { MaxStages = 3 });
        var run = new BassRunOptions { Checkpoint = true };
        BassFit cut = Bass.Fit(shortConfig, data, run);
        BassSmcReadout readout = cut.Diagnostics.Smc!;
        Assert.False(readout.ReachedPosterior);
        Assert.All(readout.Replicates, x => Assert.Equal(3, x.StageCount));
        Assert.All(readout.Replicates, x => Assert.True(x.InverseTemperature < 1.0));
        Assert.Null(cut.Evidence);
        Assert.Throws<InvalidOperationException>(() => cut.Summary());
        Assert.Throws<InvalidOperationException>(() => cut.Predict(data, PredictionCases.Of(data)));
        Assert.Throws<InvalidOperationException>(() => cut.PointwiseLogDensity(data));
        Assert.Throws<InvalidOperationException>(() => cut.Sensitivity(data));
        cut.Metadata();
        Assert.NotNull(cut.Checkpoint);

        // Each resume runs up to three more stages, exactly as the uninterrupted run went: draws, weights, stage reports,
        // normalizers and the kernels' accounts.
        BassFit resumed = cut;
        int resumes = 0;
        while (!resumed.Diagnostics.Smc!.ReachedPosterior)
        {
            resumed = Bass.Resume(resumed.Checkpoint!, shortConfig, data, run);
            Assert.True(++resumes <= stages, "the resumes never reached b = 1");
        }
        SameDraws(whole, resumed);
        for (int c = 0; c < Small.Replicates; c++)
        {
            BassSmcReplicate a = whole.Diagnostics.Smc!.Replicates[c], b = resumed.Diagnostics.Smc!.Replicates[c];
            Assert.Equal(a.Stages, b.Stages);
            Assert.Equal((a.LogNormalizer, a.LogNormalizerError, a.PriorAttempts), (b.LogNormalizer, b.LogNormalizerError, b.PriorAttempts));
            Assert.Equal(a.Account.Moves, b.Account.Moves);
            Assert.Equal(a.Account.Kinds.Entries, b.Account.Kinds.Entries);
        }
        Assert.Equal(whole.Evidence!.LogEvidence, resumed.Evidence!.LogEvidence);

        // A complete run has nothing left: its resume gives the fit again. Another config's checkpoint, and an MCMC one, are refused.
        BassFit again = Bass.Resume(resumed.Checkpoint!, shortConfig, data, run);
        SameDraws(resumed, again);
        Assert.Throws<ArgumentException>(() => Bass.Resume(resumed.Checkpoint!, SmcConfig(), data));
        BassConfig mcmc = BassSweepTests.Config(maxBasis: 4) with { Sampling = new BassSampling { Chains = 2, BurnIn = 20, Samples = 10 } };
        BassFit chains = Bass.Fit(mcmc, data, run);
        Assert.Throws<ArgumentException>(() => Bass.Resume(chains.Checkpoint!, mcmc with { Sampling = mcmc.Sampling with { Smc = Small } }, data));
    }

    [Fact]
    public void Kernel_RefusesAnotherModelsStage()
    {
        BassModel model = Model(), other = Model();
        var kernel = new BassSmcKernel(model);
        BassParticle p = new BassParticles(model).Draw(new Xoshiro256PlusPlus(3));
        SmcStage<BassParticle> foreign = Stage(other, 0.5);
        Assert.Throws<ArgumentException>(() => kernel.Move(p, foreign.Target.Evaluate(p), foreign, new Xoshiro256PlusPlus(4)));
    }
}
