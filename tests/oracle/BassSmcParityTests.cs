using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.Json;
using Maths.Regression.Spline.Bars;
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

    private static BassConfig SmcConfig(int maxBasis, SmcSpec smc, int maxInteraction = 2, bool fixedG = false)
    {
        BassConfig config = BassTargetParityTests.Config(fixedG: fixedG, maxInteraction: maxInteraction);
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
    public void Fit_NormalizerMatchesTheOraclesLogZ(string form)
    {
        if (!BassTargetParityTests.Available) return;
        bool categorical = form == "categorical";
        double[][] rows = categorical ? BassTargetParityTests.FixtureC() : BassTargetParityTests.FixtureB();
        JsonElement r = categorical
            ? BassTargetParityTests.Posterior(rows, "M_max=2", "K_max=1", "categorical=1")
            : BassTargetParityTests.Posterior(rows, "M_max=2");
        double logZ = r.GetProperty("posterior").GetProperty("log_Z").GetDouble();
        var smc = new SmcSpec { Particles = 2000, Form = form == "linear" ? TemperedForm.Linear : TemperedForm.Integrated, TraceStages = true };
        BassFit fit = Bass.Fit(SmcConfig(2, smc, maxInteraction: categorical ? 1 : 2), BassTargetParityTests.Data(rows, categorical));
        BassEvidence evidence = fit.Evidence!;
        BassSmcReadout readout = fit.Diagnostics.Smc!;
        Assert.True(readout.ReachedPosterior);
        foreach (BassSmcReplicate x in readout.Replicates)
        {
            _out.WriteLine($"{form}: replicate {x.Seed}, {x.StageCount} stages, ln Ẑ {x.LogNormalizer:F4} ± {x.LogNormalizerError:F4} " +
                           $"(z {(x.LogNormalizer - logZ) / x.LogNormalizerError:F2}), families {x.Families}");
            Assert.True(Math.Abs(x.LogNormalizer - logZ) <= 4 * x.LogNormalizerError + 0.01, $"{x.LogNormalizer:F4} ± {x.LogNormalizerError:F4} against {logZ:F4}");
        }
        _out.WriteLine($"{form}: pooled {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} (within {evidence.SteppingStone.WithinChainError:F4}, " +
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
