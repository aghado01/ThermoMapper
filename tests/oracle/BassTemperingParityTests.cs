using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's tempering (stage 4's S8) against the stage-4 enumeration oracle (<c>bass_enumeration_oracle.R</c>), which
/// computes every law below independently of the C#.
/// <list type="bullet">
/// <item><b>The b = 1 end:</b> a tempered <see cref="Bass.Fit"/>'s draws, its b = 1 replica's, against the exact
/// posterior on fixture A, as <c>BassFitParityTests</c> checks the untempered run, on both forms.</item>
/// <item><b>Every rung:</b> a ladder's replicas against the oracle's tempered laws on fixture B (<c>posterior</c> with
/// <c>rungs</c>): the count law, E[τ] and every term's inclusion at each rung, on both forms. Swaps that read anything but
/// ℓ_b, or moves that tempered the prior, would shift the interior rungs.</item>
/// <item><b>The swap's detailed balance,</b> exactly: on fixture A under fixed g, for every adjacent pair of the oracle's
/// rungs, ln π_b₁(S) − ln π_b₂(S) − (ℓ_b₁(S) − ℓ_b₂(S)) is one constant over the enumerated F, with π_b the oracle's laws
/// and ℓ_b the C#'s, which is the swap's balance. R's swap density is a defect fixture, asserted corrected: on the
/// oracle's <c>swap</c> slice (F-O6), R's swap misses balance against its own rung laws and against every candidate
/// family, where D05's balances to rounding.</item>
/// </list>
/// Opt-in: it runs where the R toolchain and the BASS checkout are present.
/// </summary>
public sealed class BassTemperingParityTests
{
    private readonly ITestOutputHelper _out;
    public BassTemperingParityTests(ITestOutputHelper output) => _out = output;

    private static BassConfig Config(int maxBasis, TemperingSpec? tempering, int chains = 4, int burn = 2_000, int samples = 20_000)
    {
        BassConfig config = BassTargetParityTests.Config();
        return config with
        {
            Structure = config.Structure with { MaxBasis = maxBasis },
            Sampling = config.Sampling with { Chains = chains, MasterSeed = 2026, BurnIn = burn, Samples = samples, Tempering = tempering },
        };
    }

    private static string Rungs(IEnumerable<double> b) => "rungs=" + string.Join(",", b.Select(v => v.ToString("R", CultureInfo.InvariantCulture)));

    [Theory]
    [InlineData("integrated")]
    [InlineData("linear")]
    public void Fit_Tempered_MatchesTheOraclesPosterior(string form)
    {
        if (!BassTargetParityTests.Available) return;
        double[][] rows = BassTargetParityTests.FixtureA();
        JsonElement r = BassTargetParityTests.Posterior(rows, "M_max=3");
        var spec = form == "linear"
            ? new TemperingSpec { Form = TemperedForm.Linear, Levels = 4, HotEnd = 0.05 }
            : new TemperingSpec { Levels = 4, HotEnd = 0.1, ReferenceRung = true };
        BassFit fit = Bass.Fit(Config(3, spec), BassTargetParityTests.Data(rows));
        Assert.NotNull(fit.Diagnostics.Tempering);

        JsonElement p = r.GetProperty("posterior");
        double[] countLaw = Doubles(p.GetProperty("count").GetProperty("posterior"));
        double[] inclusion = Doubles(p.GetProperty("inclusion"));
        BassTerm[] terms = BassTargetParityTests.Terms(r);
        var index = terms.Select((t, i) => (t, i)).ToDictionary(e => e.t, e => e.i);
        const int batches = 20;
        int chains = fit.Draws.Count, perBatch = fit.Draws[0].Count / batches;
        var count = new double[countLaw.Length, chains * batches];
        var included = new double[terms.Length, chains * batches];
        var tau = new double[chains * batches];
        for (int c = 0; c < chains; c++)
            for (int i = 0; i < perBatch * batches; i++)
            {
                BassDraw d = fit.Draws[c][i];
                int b = c * batches + i / perBatch;
                count[d.Structure.TermIds.Length, b] += 1.0 / perBatch;
                foreach (int id in d.Structure.TermIds) included[index[fit.Terms[id]], b] += 1.0 / perBatch;
                tau[b] += d.Tau / perBatch;
            }
        for (int m = 0; m < countLaw.Length; m++) Check($"{form}, M = {m}", Row(count, m), countLaw[m], 0.003);
        Check($"{form}, E[τ]", tau, p.GetProperty("tau").GetProperty("mean").GetDouble(), 0.0);
        for (int j = 0; j < terms.Length; j++) Check($"{form}, {terms[j]}", Row(included, j), inclusion[j], 0.005, quiet: true);
        BassTemperingReadout readout = fit.Diagnostics.Tempering!;
        _out.WriteLine($"{form}: {terms.Length} inclusions within their allowance; round trips {readout.RoundTripsPerScan:F4}/scan, " +
                       $"swap acceptance [{string.Join(", ", readout.PairAcceptance.Select(a => a.ToString("F2")))}]");
    }

    [Theory]
    [InlineData("integrated")]
    [InlineData("linear")]
    public void Ladder_EachRungMatchesTheOraclesTemperedLaw(string form)
    {
        if (!BassTargetParityTests.Available) return;
        double[][] rows = BassTargetParityTests.FixtureB();
        var spec = new TemperingSpec
        {
            Form = form == "linear" ? TemperedForm.Linear : TemperedForm.Integrated, Levels = 4, HotEnd = 0.1, ReferenceRung = true, Adaptive = false,
        };
        double[] betas = ParallelTempering<BassState>.ReferenceLadder(spec.Levels, spec.HotEnd);
        JsonElement r = BassTargetParityTests.Posterior(rows, "M_max=2", Rungs(betas));
        BassTerm[] terms = BassTargetParityTests.Terms(r);
        int[][] states = BassTargetParityTests.States(r);
        var index = terms.Select((t, i) => (t, i)).ToDictionary(e => e.t, e => e.i);
        JsonElement[] rung = [.. r.GetProperty("posterior").GetProperty("tempered").GetProperty(form).EnumerateArray()];

        BassModel model = BassModel.Create(Config(2, spec), BassTargetParityTests.Data(rows));
        const int ladders = 4, burn = 2_000, scans = 30_000, batches = 20, perBatch = scans / batches;
        int levels = betas.Length;
        var count = new double[levels][,];
        var included = new double[levels][,];
        var tau = new double[levels][];
        for (int k = 0; k < levels; k++)
        {
            count[k] = new double[model.MaxBasis + 1, ladders * batches];
            included[k] = new double[terms.Length, ladders * batches];
            tau[k] = new double[ladders * batches];
        }
        var sampler = new BassPriorSampler(model);
        for (int l = 0; l < ladders; l++)
        {
            BassState start = sampler.DrawState(new Xoshiro256PlusPlus(300 + l)).State;
            var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, seed: 400 + l);
            for (int t = 0; t < burn; t++) ladder.Step();
            for (int t = 0; t < scans; t++)
            {
                ladder.Step();
                int b = l * batches + t / perBatch;
                for (int k = 0; k < levels; k++)
                {
                    BassState s = ladder.Tempering.Replicas[k].Current;
                    count[k][s.Count, b] += 1.0 / perBatch;
                    foreach (BassTerm term in s.Terms) included[k][index[term], b] += 1.0 / perBatch;
                    tau[k][b] += s.Tau / perBatch;
                }
            }
        }

        for (int k = 0; k < levels; k++)
        {
            Assert.Equal(betas[k], rung[k].GetProperty("b").GetDouble());
            double[] law = Doubles(rung[k].GetProperty("law"));
            double[] countLaw = Doubles(rung[k].GetProperty("count"));
            var inclusion = new double[terms.Length];
            for (int s = 0; s < states.Length; s++)
                foreach (int j in states[s]) inclusion[j] += law[s];
            for (int m = 0; m < countLaw.Length; m++) Check($"{form}, b = {betas[k]:G3}, M = {m}", Row(count[k], m), countLaw[m], 0.004);
            Check($"{form}, b = {betas[k]:G3}, E[τ]", tau[k], rung[k].GetProperty("tau_mean").GetDouble(), 0.0);
            for (int j = 0; j < terms.Length; j++) Check($"{form}, b = {betas[k]:G3}, {terms[j]}", Row(included[k], j), inclusion[j], 0.005, quiet: true);
        }
    }

    [Theory]
    [InlineData("integrated")]
    [InlineData("linear")]
    public void Swap_BalancesOnTheOraclesRungLaws(string form)
    {
        if (!BassTargetParityTests.Available) return;
        double[][] rows = BassTargetParityTests.FixtureA();
        double[] rungs = [0.0, 0.05, 0.2, 0.5, 1.0];
        JsonElement r = BassTargetParityTests.Posterior(rows, "M_max=2", "coef=fixed", Rungs(rungs));
        BassTerm[] terms = BassTargetParityTests.Terms(r);
        int[][] states = BassTargetParityTests.States(r);
        JsonElement[] rung = [.. r.GetProperty("posterior").GetProperty("tempered").GetProperty(form).EnumerateArray()];
        BassConfig config = BassTargetParityTests.Config(fixedG: true);
        config = config with { Structure = config.Structure with { MaxBasis = 2 } };
        BassModel model = BassModel.Create(config, BassTargetParityTests.Data(rows), form == "linear" ? TemperedForm.Linear : TemperedForm.Integrated);
        BassDesign design = model.NewDesign();
        var evaluations = new BassEvaluation[states.Length];
        for (int s = 0; s < states.Length; s++)
        {
            BassTerm[] structure = [.. states[s].Select(i => terms[i])];
            Assert.Null(design.Load(structure));
            evaluations[s] = model.Target.Evaluate(BassStructureFit.Of(design, structure.Sum(model.Prior.LogPrior)));
        }
        for (int k = 0; k + 1 < rungs.Length; k++)
        {
            double b1 = rungs[k + 1], b2 = rungs[k];
            double[] law1 = Doubles(rung[k + 1].GetProperty("law")), law2 = Doubles(rung[k].GetProperty("law"));
            double lo = double.PositiveInfinity, hi = double.NegativeInfinity;
            int used = 0;
            for (int s = 0; s < states.Length; s++)
            {
                if (law1[s] < 1e-250 || law2[s] < 1e-250) continue;
                // The swap of replica b₁ holding S with replica b₂ reads ℓ_b alone: d_b₁ − d_b₂ = ℓ_b₁ − ℓ_b₂, the prior cancelling.
                var t1 = new Temperature(b1, TemperingPath.Likelihood);
                var t2 = new Temperature(b2, TemperingPath.Likelihood);
                double d = Math.Log(law1[s]) - Math.Log(law2[s]) - (t1.LogDensity(evaluations[s]) - t2.LogDensity(evaluations[s]));
                lo = Math.Min(lo, d);
                hi = Math.Max(hi, d);
                used++;
            }
            _out.WriteLine($"{form}, b = {b1} ↔ {b2}: spread {hi - lo:G3} over {used} of {states.Length} states");
            Assert.True(used > states.Length / 2);
            Assert.True(hi - lo <= 1e-8, $"b = {b1} ↔ {b2}: the swap misses balance by {hi - lo:G3}");
        }
    }

    [Fact]
    public void Swap_RsDensity_IsADefectFixture()
    {
        if (!BassTargetParityTests.Available) return;
        string csv = Path.Combine(Path.GetTempPath(), $"bass_swap_{Guid.NewGuid():N}.csv");
        var text = new StringBuilder();
        foreach (double[] row in BassTargetParityTests.FixtureA()) text.AppendLine(string.Join(",", row.Select(v => v.ToString("R", CultureInfo.InvariantCulture))));
        File.WriteAllText(csv, text.ToString());
        JsonElement swap;
        try
        {
            swap = ROracle.Run("oracles/bass_enumeration_oracle.R", csv, "swap").GetProperty("swap");
        }
        finally
        {
            File.Delete(csv);
        }
        foreach (JsonElement pair in swap.GetProperty("pairs").EnumerateArray())
        {
            double b1 = pair.GetProperty("b1").GetDouble(), b2 = pair.GetProperty("b2").GetDouble();
            double own = pair.GetProperty("r_swap_vs_r_kernels").GetDouble();
            double corrected = pair.GetProperty("corrected_swap_vs_corrected_kernels").GetDouble();
            var families = pair.GetProperty("r_swap_vs_family").EnumerateObject().Select(f => (f.Name, Spread: f.Value.GetDouble())).ToList();
            _out.WriteLine($"b = {b1} ↔ {b2}: R's swap misses its own rung laws by {own:F3} and the families by " +
                           $"{string.Join(", ", families.Select(f => $"{f.Name} {f.Spread:F3}"))}; D05's misses its own by {corrected:G2}");
            // R's swap balances against no law, its own kernels' included: a defect, which no TemperingSpec reproduces.
            Assert.True(own > 0.3);
            Assert.Equal(4, families.Count);
            Assert.All(families, f => Assert.True(f.Spread > 0.3, $"{f.Name}: {f.Spread}"));
            Assert.True(corrected <= 1e-10);
        }
        foreach (JsonElement rung in swap.GetProperty("rung").EnumerateArray())
            Assert.True(rung.GetProperty("tv_corrected_target").GetDouble() <= 1e-10);
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
