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
            if (fixedG) Assert.All(tau, t => Assert.Equal(model.Target.FixedTau!.Value, t, 12));   // held: no τ block runs
            else Check($"{form}, b = {stages[k]}, E[τ]", tau, law.GetProperty("tau_mean").GetDouble(), 0.0);
            for (int j = 0; j < terms.Length; j++) Check($"{form}, b = {stages[k]}, {terms[j]}", Row(included, j), inclusion[j], 0.005, quiet: true);
            _out.WriteLine($"{form}, b = {stages[k]}: acceptance {kernel.Account.AcceptanceRate:F3} over {kernel.Account.Sweeps} sweeps");
        }
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
