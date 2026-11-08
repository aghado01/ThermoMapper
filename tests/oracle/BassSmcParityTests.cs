using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.Json;
using Maths.Regression.Spline.Bars;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's SMC (stage 8, C4) against the stage-4 enumeration oracle (<c>bass_enumeration_oracle.R</c>), which computes every
/// law below independently of the C#.
/// <list type="bullet">
/// <item><b>Each stage's law, through the kernel alone:</b> the SMC kernel's sweeps at a fixed b, run as a chain from exact
/// prior draws, against the oracle's tempered law on fixture A at three b, on both forms and under fixed g: the count law,
/// E[τ] and every term's inclusion. A kernel that left π_b's prior tempered, or read another b, would miss them.</item>
/// <item><b>The evidence:</b> an SMC fit's normalizer under the proper defaults, each replicate's within four of its
/// genealogy errors and the pooled one within four of its error, against the oracle's exact ln Z on fixture B, on both
/// forms, and on fixture C with its categorical predictor; stepping-stone over the traced stages reproduces each
/// normalizer, and path sampling on the linear form stands within its error and diagnostic.</item>
/// <item><b>The posterior:</b> an SMC fit's weighted populations against the exact posterior on fixture A, the defaults at
/// M_max = 3 and fixed g at M_max = 2: the count law, E[τ] and every term's inclusion, each within five standard errors
/// of the replicates' spread.</item>
/// <item><b>Under the adaptive count of sweeps</b> (<see cref="AdaptiveMutation"/>): populations held at b = 0 keep the prior's law,
/// populations tempered to b = 1 and held there meet the exact posterior on fixture A, and the normalizer meets ln Z on fixture B.</item>
/// </list>
/// Opt-in: it runs where the R toolchain and the BASS checkout are present.
/// </summary>
public sealed class BassSmcParityTests
{
    private readonly ITestOutputHelper _out;
    public BassSmcParityTests(ITestOutputHelper output) => _out = output;

    private static string Rungs(IEnumerable<double> b) => "rungs=" + string.Join(",", b.Select(v => v.ToString("R", CultureInfo.InvariantCulture)));

    [Theory]
    [InlineData("integrated", false)]
    [InlineData("linear", false)]
    [InlineData("integrated", true)]
    public void Kernel_LeavesEachStagesLawInvariant(string form, bool fixedG)
    {
        if (!BassTargetParityTests.Available) return;
        double[][] rows = BassTargetParityTests.FixtureA();
        double[] stages = [0.05, 0.3, 1.0];
        int maxBasis = fixedG ? 2 : 3;
        JsonElement r = fixedG
            ? BassTargetParityTests.Posterior(rows, $"M_max={maxBasis}", "coef=fixed", Rungs(stages))
            : BassTargetParityTests.Posterior(rows, $"M_max={maxBasis}", Rungs(stages));
        BassTerm[] terms = BassTargetParityTests.Terms(r);
        int[][] states = BassTargetParityTests.States(r);
        var index = terms.Select((t, i) => (t, i)).ToDictionary(e => e.t, e => e.i);
        JsonElement[] rung = [.. r.GetProperty("posterior").GetProperty("tempered").GetProperty(form.Split(' ')[0]).EnumerateArray()];
        BassConfig config = BassTargetParityTests.Config(fixedG: fixedG);
        config = config with { Structure = config.Structure with { MaxBasis = maxBasis } };
        if (fixedG) form = $"{form} (fixed g)";
        BassModel model = BassModel.Create(config, BassTargetParityTests.Data(rows), form.StartsWith("linear") ? TemperedForm.Linear : TemperedForm.Integrated);
        var particles = new BassParticles(model);
        var target = new BassParticleTarget(model.Target);

        const int chains = 4, burn = 2_000, sweeps = 30_000, batches = 20, perBatch = sweeps / batches;
        for (int k = 0; k < stages.Length; k++)
        {
            JsonElement law = rung.Single(x => x.GetProperty("b").GetDouble() == stages[k]);
            var stage = new SmcStage<BassParticle>(target, new Temperature(stages[k], TemperingPath.Likelihood));
            var count = new double[model.MaxBasis + 1, chains * batches];
            var included = new double[terms.Length, chains * batches];
            var tau = new double[chains * batches];
            var kernel = new BassSmcKernel(model);
            for (int c = 0; c < chains; c++)
            {
                BassParticle p = particles.Draw(new Xoshiro256PlusPlus(700 + 10 * k + c));
                var rng = new Xoshiro256PlusPlus(800 + 10 * k + c);
                ITargetEvaluation e = target.Evaluate(p);
                for (int t = 0; t < burn + sweeps; t++)
                {
                    e = kernel.Move(p, e, stage, rng).Evaluation;
                    if (t < burn) continue;
                    int b = c * batches + (t - burn) / perBatch;
                    count[p.State.Count, b] += 1.0 / perBatch;
                    foreach (BassTerm term in p.State.Terms) included[index[term], b] += 1.0 / perBatch;
                    tau[b] += p.State.Tau / perBatch;
                }
            }
            double[] probabilities = Doubles(law.GetProperty("law"));
            double[] countLaw = Doubles(law.GetProperty("count"));
            var inclusion = new double[terms.Length];
            for (int s = 0; s < states.Length; s++)
                foreach (int j in states[s]) inclusion[j] += probabilities[s];
            for (int m = 0; m < countLaw.Length; m++) Check($"{form}, b = {stages[k]}, M = {m}", Row(count, m), countLaw[m], 0.004);
            if (fixedG) Assert.All(tau, t => Assert.Equal(model.Target.HeldTau!.Value, t, 12));   // held: no τ block runs
            else Check($"{form}, b = {stages[k]}, E[τ]", tau, law.GetProperty("tau_mean").GetDouble(), 0.0);
            for (int j = 0; j < terms.Length; j++) Check($"{form}, b = {stages[k]}, {terms[j]}", Row(included, j), inclusion[j], 0.005, quiet: true);
            _out.WriteLine($"{form}, b = {stages[k]}: acceptance {kernel.Account.AcceptanceRate:F3} over {kernel.Account.Sweeps} sweeps");
        }
    }

    private static BassConfig SmcConfig(int maxBasis, SmcSpec smc, int maxInteraction = 2, bool fixedG = false, BassCoefficientPrior? coefficients = null)
    {
        BassConfig config = BassTargetParityTests.Config(fixedG: fixedG, maxInteraction: maxInteraction, coefficients: coefficients);
        return config with
        {
            Structure = config.Structure with { MaxBasis = maxBasis },
            Sampling = config.Sampling with { MasterSeed = 2026, Smc = smc },
        };
    }

    [Theory]
    [InlineData("integrated")]
    [InlineData("linear")]
    [InlineData("categorical")]
    [InlineData("ridge")]
    [InlineData("ridge, categorical")]
    [InlineData("ridge, held")]
    public void Fit_NormalizerMatchesTheOraclesLogZ(string form)
    {
        if (!BassTargetParityTests.Available) return;
        bool categorical = form.EndsWith("categorical", StringComparison.Ordinal), ridge = form.StartsWith("ridge", StringComparison.Ordinal);
        bool held = form == "ridge, held";
        double[][] rows = categorical ? BassTargetParityTests.FixtureC() : BassTargetParityTests.FixtureB();
        string[] formArgs = ridge ? held ? ["form=ridge", "coef=fixed", "tau=0.4"] : ["form=ridge"] : [];
        JsonElement r = categorical
            ? BassTargetParityTests.Posterior(rows, ["M_max=2", "K_max=1", "categorical=1", .. formArgs])
            : BassTargetParityTests.Posterior(rows, ["M_max=2", .. formArgs]);
        double logZ = r.GetProperty("posterior").GetProperty("log_Z").GetDouble();
        var smc = new SmcSpec { Particles = 2000, Form = form == "linear" ? TemperedForm.Linear : TemperedForm.Integrated, TraceStages = true };
        BassCoefficientPrior? coefficients = ridge ? BassCoefficientPrior.Ridge(held ? TauPrior.Held(0.4) : TauPrior.UnitInformationGamma()) : null;
        BassFit fit = Bass.Fit(SmcConfig(2, smc, maxInteraction: categorical ? 1 : 2, coefficients: coefficients), BassTargetParityTests.Data(rows, categorical));
        BassEvidence evidence = fit.Evidence!;
        BassSmcReadout readout = fit.Diagnostics.Smc!;
        Assert.True(readout.ReachedPosterior);
        // SMC's flags (the KHAOS trap's decision, 2026-11-08) read clear on these fits, which meet R's enumeration.
        _out.WriteLine($"{form}: between-to-within {readout.NormalizerBetweenToWithin:F2}; last steps from b " +
                       string.Join(", ", readout.Replicates.Select(x => $"{x.LastStepFrom:G3} over {x.StructuresBeforeLastStep} of {x.Particles}")));
        Assert.Empty(readout.Flags);
        foreach (BassSmcReplicate x in readout.Replicates)
        {
            _out.WriteLine($"{form}: replicate {x.Seed}, {x.StageCount} stages, ln Ẑ {x.LogNormalizer:F4} ± {x.LogNormalizerError:F4} " +
                           $"(z {(x.LogNormalizer - logZ) / x.LogNormalizerError:F2}), families {x.Families}");
            Assert.True(Math.Abs(x.LogNormalizer - logZ) <= 4 * x.LogNormalizerError + 0.01, $"{x.LogNormalizer:F4} ± {x.LogNormalizerError:F4} against {logZ:F4}");
        }
        _out.WriteLine($"{form}: pooled {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} (within {evidence.SteppingStone!.WithinChainError:F4}, " +
                       $"between {evidence.SteppingStone.BetweenChainError:F4}) against the oracle's {logZ:F4}");
        Assert.True(Math.Abs(evidence.LogEvidence - logZ) <= 4 * evidence.MonteCarloError + 0.01,
                    $"{evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against {logZ:F4}");
        // Stepping-stone over the traced stages is the normalizer, rung by rung.
        for (int c = 0; c < readout.Replicates.Count; c++)
        {
            Assert.Equal(readout.Replicates[c].LogNormalizer, evidence.SteppingStone.Chains[c].LogRatio, 9);
            Assert.Equal(readout.Replicates[c].StageCount, evidence.SteppingStone.Chains[c].Rungs.Count);
        }
        if (evidence.PathSampling is { } ps)
        {
            _out.WriteLine($"{form}: path sampling {ps.LogEvidence:F4} ± {ps.MonteCarloError:F4}, diagnostic {ps.DiscretizationError:F4}");
            Assert.True(Math.Abs(ps.LogEvidence + ps.DiscretizationError - logZ) <= 4 * ps.MonteCarloError + 0.5 * Math.Abs(ps.DiscretizationError) + 0.01);
        }
        Assert.Equal(form == "linear", evidence.PathSampling is not null);
    }

    [Theory]
    [InlineData("defaults")]
    [InlineData("fixed g")]
    public void Fit_MatchesTheOraclesPosterior(string option)
    {
        if (!BassTargetParityTests.Available) return;
        bool fixedG = option == "fixed g";
        int maxBasis = fixedG ? 2 : 3;
        double[][] rows = BassTargetParityTests.FixtureA();
        JsonElement r = fixedG
            ? BassTargetParityTests.Posterior(rows, $"M_max={maxBasis}", "coef=fixed")
            : BassTargetParityTests.Posterior(rows, $"M_max={maxBasis}");
        JsonElement p = r.GetProperty("posterior");
        double[] countLaw = Doubles(p.GetProperty("count").GetProperty("posterior"));
        double[] inclusion = Doubles(p.GetProperty("inclusion"));
        BassTerm[] terms = BassTargetParityTests.Terms(r);
        var index = terms.Select((t, i) => (t, i)).ToDictionary(e => e.t, e => e.i);
        const int replicates = 8;
        BassFit fit = Bass.Fit(SmcConfig(maxBasis, new SmcSpec { Particles = 1000, Replicates = replicates }, fixedG: fixedG), BassTargetParityTests.Data(rows));
        Assert.True(fit.Diagnostics.Smc!.ReachedPosterior);
        BassSmcReadout flags = fit.Diagnostics.Smc;
        _out.WriteLine($"{option}: between-to-within {flags.NormalizerBetweenToWithin:F2}; last steps from b " +
                       string.Join(", ", flags.Replicates.Select(x => $"{x.LastStepFrom:G3} over {x.StructuresBeforeLastStep} of {x.Particles}")));
        Assert.Empty(flags.Flags);

        // Each replicate's weighted estimate; their spread gives each comparison's error.
        double[] w = fit.DrawWeights!;
        var count = new double[countLaw.Length, replicates];
        var included = new double[terms.Length, replicates];
        var tau = new double[replicates];
        for (int c = 0, s = 0; c < replicates; c++)
            foreach (BassDraw d in fit.Draws[c])
            {
                double weight = w[s++] * replicates;
                count[d.Structure.TermIds.Length, c] += weight;
                foreach (int id in d.Structure.TermIds) included[index[fit.Terms[id]], c] += weight;
                tau[c] += weight * d.Tau;
            }
        for (int m = 0; m < countLaw.Length; m++) Check($"{option}, M = {m}", Row(count, m), countLaw[m], 0.003);
        if (fixedG) Assert.All(fit.Draws.SelectMany(c => c), d => Assert.Equal(fit.Resolved.Scale.HeldTau!.Value, d.Tau));
        else Check($"{option}, E[τ]", tau, p.GetProperty("tau").GetProperty("mean").GetDouble(), 0.0);
        for (int j = 0; j < terms.Length; j++) Check($"{option}, {terms[j]}", Row(included, j), inclusion[j], 0.005, quiet: true);
        _out.WriteLine($"{option}: {terms.Length} inclusions within their allowance; stages " +
                       $"{string.Join(", ", fit.Diagnostics.Smc.Replicates.Select(x => x.StageCount))}");
    }

    [Fact]
    public void AdaptiveMutation_KeepsThePriorAtBZero_AndMeetsThePosteriorAndLogZAtBOne()
    {
        // Under the adaptive count every stage's sweeps keep its law: SMC populations through BASS's particles and kernel, held at b = 0
        // for twenty stages, and tempered to b = 1 and held there five more, against the oracle's laws on fixture A; and the normalizer
        // against its exact ln Z on fixture B. Each comparison's error is the replicates' spread.
        if (!BassTargetParityTests.Available) return;
        const int replicates = 8, particles = 1000, rate = 8;
        double[][] rows = BassTargetParityTests.FixtureA();
        BassConfig config = BassTargetParityTests.Config();
        config = config with { Structure = config.Structure with { MaxBasis = 3 } };
        BassModel model = BassModel.Create(config, BassTargetParityTests.Data(rows), TemperedForm.Integrated);
        JsonElement r = BassTargetParityTests.Posterior(rows, "M_max=3", Rungs([0.0]));
        BassTerm[] terms = BassTargetParityTests.Terms(r);
        int[][] states = BassTargetParityTests.States(r);
        var index = terms.Select((t, i) => (t, i)).ToDictionary(e => e.t, e => e.i);
        JsonElement prior = r.GetProperty("posterior").GetProperty("tempered").GetProperty("integrated").EnumerateArray().Single(x => x.GetProperty("b").GetDouble() == 0.0);
        JsonElement posterior = r.GetProperty("posterior");
        int[] seeds = SeedTree.Derive(2026_11_08, replicates);

        foreach (bool atOne in new[] { false, true })
        {
            var count = new double[model.MaxBasis + 1, replicates];
            var included = new double[terms.Length, replicates];
            var tau = new double[replicates];
            for (int c = 0; c < replicates; c++)
            {
                SmcRun<BassParticle> run = AdaptiveRun(model, particles, seeds[c], rate);
                SmcStage<BassParticle> stage = run.Stage;
                if (atOne)
                {
                    SmcSampler.TemperToPosterior(run, 0.99);
                    stage = run.Stage;
                }
                for (int s = 0; s < (atOne ? 5 : 20); s++) run.Advance(stage);
                SmcResult<BassParticle> result = run.Result();
                for (int i = 0; i < particles; i++)
                {
                    BassState state = result.Particles[i].State;
                    double w = result.Weights[i];
                    count[state.Count, c] += w;
                    foreach (BassTerm term in state.Terms) included[index[term], c] += w;
                    tau[c] += w * state.Tau;
                }
                int[] steps = SmcStageReadout.Of(result.Stages).Steps;
                if (c == 0) _out.WriteLine($"b = {(atOne ? 1 : 0)}: steps {string.Join(" ", steps)}");
                Assert.All(steps.Skip(1), k => Assert.InRange(k, 1, int.MaxValue));
                Assert.True(steps.Skip(1).Any(k => k != rate), "the count never moved off the budget, so the run did not exercise the rule");
            }
            string at = atOne ? "b = 1" : "b = 0";
            double[] countLaw = atOne ? Doubles(posterior.GetProperty("count").GetProperty("posterior")) : Doubles(prior.GetProperty("count"));
            double[] inclusion;
            if (atOne) inclusion = Doubles(posterior.GetProperty("inclusion"));
            else
            {
                double[] probabilities = Doubles(prior.GetProperty("law"));
                inclusion = new double[terms.Length];
                for (int s = 0; s < states.Length; s++)
                    foreach (int j in states[s]) inclusion[j] += probabilities[s];
            }
            for (int m = 0; m < countLaw.Length; m++) Check($"adaptive, {at}, M = {m}", Row(count, m), countLaw[m], 0.003);
            Check($"adaptive, {at}, E[τ]", tau, atOne ? posterior.GetProperty("tau").GetProperty("mean").GetDouble() : prior.GetProperty("tau_mean").GetDouble(), 0.0);
            for (int j = 0; j < terms.Length; j++) Check($"adaptive, {at}, {terms[j]}", Row(included, j), inclusion[j], 0.005, quiet: true);
        }

        // The normalizer on fixture B, M_max 2.
        double[][] rowsB = BassTargetParityTests.FixtureB();
        double logZ = BassTargetParityTests.Posterior(rowsB, "M_max=2").GetProperty("posterior").GetProperty("log_Z").GetDouble();
        BassConfig configB = BassTargetParityTests.Config();
        BassModel modelB = BassModel.Create(configB with { Structure = configB.Structure with { MaxBasis = 2 } }, BassTargetParityTests.Data(rowsB), TemperedForm.Integrated);
        var normalizers = new double[4];
        double within = 0.0;
        for (int c = 0; c < normalizers.Length; c++)
        {
            SmcResult<BassParticle> result = SmcSampler.TemperToPosterior(AdaptiveRun(modelB, 2000, seeds[c], rate), 0.99);
            normalizers[c] = result.LogNormalizerRatio;
            within += result.LogNormalizerError * result.LogNormalizerError;
            _out.WriteLine($"adaptive: ln Ẑ {result.LogNormalizerRatio:F4} ± {result.LogNormalizerError:F4} against {logZ:F4}");
            Assert.True(Math.Abs(result.LogNormalizerRatio - logZ) <= 4 * result.LogNormalizerError + 0.01, $"{result.LogNormalizerRatio:F4} ± {result.LogNormalizerError:F4} against {logZ:F4}");
        }
        // Pooled as the fit's evidence pools its replicates: the larger of the within-run (genealogy) and between-run errors.
        double pooled = normalizers.Average(), between = Math.Sqrt(normalizers.Sum(v => (v - pooled) * (v - pooled)) / (normalizers.Length - 1) / normalizers.Length);
        double error = Math.Max(Math.Sqrt(within) / normalizers.Length, between);
        _out.WriteLine($"adaptive, pooled ln Ẑ {pooled:F4} ± {error:F4} (within {Math.Sqrt(within) / normalizers.Length:F4}, between {between:F4}) against {logZ:F4}");
        Assert.True(Math.Abs(pooled - logZ) <= 4 * error + 0.01, $"pooled {pooled:F4} ± {error:F4} against {logZ:F4}");
    }

    // An SMC run through BASS's particles and kernel at b = 0, under the adaptive count on `rate` sweeps per stage's budget. A stay
    // probability of 0.3 lets fixture A's stages want fewer than the budget, so the count moves.
    private static SmcRun<BassParticle> AdaptiveRun(BassModel model, int particles, int seed, int rate)
    {
        var options = new SmcOptions
        {
            Particles = particles, MutationSteps = rate, AdaptiveMutation = new AdaptiveMutation { StayProbability = 0.3 }, MaxDegreeOfParallelism = 4,
        };
        return SmcRun<BassParticle>.Start(new SmcStage<BassParticle>(new BassParticleTarget(model.Target), new Temperature(0.0, TemperingPath.Likelihood)),
                                          new BassParticles(model), new BassSmcKernel(model), options, seed);
    }

    private static double[] Doubles(JsonElement e) => [.. e.EnumerateArray().Select(v => v.GetDouble())];

    private static double[] Row(double[,] batches, int row) => [.. Enumerable.Range(0, batches.GetLength(1)).Select(b => batches[row, b])];

    // A batch-means comparison: within five standard errors plus `slack`.
    private void Check(string what, double[] batchMeans, double expected, double slack, bool quiet = false)
    {
        double mean = batchMeans.Average();
        double error = Math.Sqrt(batchMeans.Sum(v => (v - mean) * (v - mean)) / (batchMeans.Length - 1) / batchMeans.Length);
        if (!quiet) _out.WriteLine($"{what}: {mean:G5} against {expected:G5} (± {error:G2})");
        Assert.True(Math.Abs(mean - expected) <= 5 * error + slack, $"{what}: {mean:G6} against {expected:G6} (± {error:G2})");
    }
}
