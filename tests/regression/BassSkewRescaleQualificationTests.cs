using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Maths.Distributions;
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
/// The skew's rescaling move's worth (<see cref="BassNoise.SkewRescale"/>), behind <c>BASS_QUALIFICATION</c> = 1: R̂ and the bulk and tail
/// effective sizes per second of the chains' work with and without it, as slice 7's row reads them; a ladder's swap rates under the skewed
/// horseshoe with and without it; and its cost per sweep. Its law is the reach harness's Geweke cases
/// (<see cref="BassNoiseReachQualificationTests"/>, the "+rescale" forms) and its statement on a frozen state
/// (<see cref="BassSkewRescaleTests"/>).
/// </summary>
public sealed class BassSkewRescaleQualificationTests(ITestOutputHelper output)
{
    private static readonly string[] Names = ["ln sigma2", "beta", "mean ln v", "M"];

    private static readonly (string Name, BassNoise Noise)[] Noises =
    [
        ("skewed horseshoe", BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5))),
        ("skewed GH", BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), new BassNoisePrior(0.0, 1.0))),
        ("skewed Normal-Wald", BassNoise.NormalWald(skew: new BassNoisePrior(0.0, 1.0))),
    ];

    private static BassData Data(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x0 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())], x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        // Skewed noise outside every family here: 0.3·(E − 1), E ~ Exp(1).
        double[] y = [.. Enumerable.Range(0, n).Select(i => Math.Sin(2.0 * Math.PI * x0[i]) + x1[i] + 0.3 * (-Math.Log(1.0 - rng.NextDouble()) - 1.0))];
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
    }

    [Theory]
    [InlineData(11)]
    [InlineData(23)]
    public void EffectiveSizesPerUnitOfWork_WithAndWithoutTheMove(int masterSeed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassData data = Data(300, 20261108);
        output.WriteLine($"Seed {masterSeed}. Work: the fit's thread cycles on one thread (MaxDegreeOfParallelism 1), indicative under load; " +
                         "4 chains × 8,000 draws after 2,000, at most 6 additive terms.");
        foreach ((string name, BassNoise noise) in Noises)
        {
            var perSecond = new double[2][];
            string[] forms = ["plain", "with the move"];
            for (int w = 0; w < 2; w++)
            {
                bool move = w > 0;
                BassConfig config = BassStudentTTests.Config(noise with { SkewRescale = move }) with
                {
                    Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 6 },
                    Sampling = new BassSampling { Chains = 4, MasterSeed = masterSeed, BurnIn = 2000, Samples = 8000 },
                };
                var wall = Stopwatch.StartNew();
                long t0 = BassCost.Tick();
                BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
                double seconds = (BassCost.Tick() - t0) / BassCost.TicksPerSecond;
                wall.Stop();
                var cells = new List<string>();
                perSecond[w] = new double[Names.Length];
                for (int g = 0; g < Names.Length; g++)
                {
                    double[][] chains = [.. fit.Draws.Select(chain => chain.Select(d => Functional(d, g)).ToArray())];
                    ConvergenceDiagnostic diagnostic = ChainDiagnostics.Diagnose(chains);
                    double bulk = diagnostic.BulkEss.Value;
                    perSecond[w][g] = bulk / seconds;
                    cells.Add($"{Names[g]} R̂ {diagnostic.RHat.Value:F3} ESS {bulk:F0}/{diagnostic.TailEss.Value:F0} ({bulk / seconds:F1}/s)");
                }
                MoveStat[] block = [.. fit.Diagnostics.Chains.SelectMany(c => c.Blocks!).Where(k => k.Key == "skew rescale")];
                long accepted = block.Sum(k => k.Accepted), attempts = block.Sum(k => k.Attempts);
                output.WriteLine($"{name}, {forms[w]}: {seconds:F1} s ({wall.Elapsed.TotalSeconds:F1} s wall)" +
                                 (attempts == 0 ? "" : $", the move accepted {accepted} of {attempts}") + $"; {string.Join("; ", cells)}");
            }
            for (int w = 1; w < 2; w++)
                output.WriteLine($"{name}: bulk ESS per second, {forms[w]} over plain: " +
                                 string.Join(", ", Names.Select((label, g) => $"{label} {perSecond[w][g] / perSecond[0][g]:F2}")));
        }
    }

    [Fact]
    public void SwapRates_OnALadderUnderTheSkewedHorseshoe()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassData data = Data(300, 20261108);
        foreach (bool move in new[] { false, true })
        {
            BassConfig config = BassStudentTTests.Config(BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5)) with { SkewRescale = move }) with
            {
                Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 6 },
                Sampling = new BassSampling { Chains = 2, MasterSeed = 17, BurnIn = 2000, Samples = 4000,
                                              Tempering = new TemperingSpec { Levels = 6, HotEnd = 0.05 } },
            };
            BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 2 });
            BassTemperingReadout ladder = fit.Diagnostics.Tempering!;
            double[] acceptance = [.. Enumerable.Range(0, ladder.PairAcceptance[0].Length).Select(k => ladder.PairAcceptance.Average(chain => chain[k]))];
            double[][] sigma2 = [.. fit.Draws.Select(chain => chain.Select(d => Math.Log(d.Sigma2)).ToArray())];
            ConvergenceDiagnostic diagnostic = ChainDiagnostics.Diagnose(sigma2);
            output.WriteLine($"{(move ? "with the move" : "plain")}: swap acceptance by pair, cold first, {string.Join(", ", acceptance.Select(a => a.ToString("F3")))}; " +
                             $"round trips {ladder.RoundTripsPerScan:E2}/scan; the cold chain's ln σ² R̂ {diagnostic.RHat.Value:F3}, bulk ESS {diagnostic.BulkEss.Value:F0}");
        }
    }

    [Fact]
    public void CostRow_ASweep_WithAndWithoutTheMove()
    {
        if (!BassQualification.Enabled) return;   // heavy
        string clock = BassCost.Cycles ? "thread cycles, indicative" : "wall clock";
        output.WriteLine($"Clock: {clock}. µs per sweep of one chain on this thread under the skewed horseshoe, after 500 untimed sweeps; an untimed pass first.");
        var points = new (string Name, BassCost.Point Point)[] { ("the base", BassCost.Base), ("n 200", BassCost.Base with { N = 200 }) };
        for (int pass = 0; pass < 2; pass++)
            foreach ((string name, BassCost.Point x) in points)
            {
                BassData data = BassCost.Data(x);
                var cells = new List<string>();
                foreach (bool move in new[] { false, true })
                {
                    BassModel model = BassModel.Create(BassCost.Config(x) with { Noise = BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5)) with { SkewRescale = move } }, data);
                    var moves = new BassStructureMoves(model, [], 0.0);
                    var rng = new Xoshiro256PlusPlus(5);
                    var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(3));
                    List<IBlock<ReversibleJumpChain<BassState>>> blocks = BassSweepSlots.Blocks(model, BassSweepSlot.Tau);
                    long accepted = 0, attempts = 0;
                    void Sweep(bool count)
                    {
                        chain.Step();
                        foreach (IBlock<ReversibleJumpChain<BassState>> block in blocks)
                        {
                            bool ok = block.Update(chain, rng);
                            if (count && block.Name == "skew rescale") { attempts++; accepted += ok ? 1 : 0; }
                        }
                    }
                    for (int s = 0; s < 500; s++) Sweep(false);
                    const int timed = 200;
                    long t0 = BassCost.Tick();
                    for (int s = 0; s < timed; s++) Sweep(true);
                    double us = (BassCost.Tick() - t0) / (double)timed / BassCost.TicksPerSecond * 1e6;
                    cells.Add($"{(move ? "with the move" : "plain")} {us:F1} (M {chain.Current.Count}" + (move ? $", accepted {accepted} of {attempts})" : ")"));
                }
                if (pass == 1) output.WriteLine($"{name,-9} rows {x.Rows,5}: {string.Join("; ", cells)}");
            }
    }

    // ln σ², β, mean ln v and M of a retained draw.
    private static double Functional(BassDraw draw, int g) => g switch
    {
        0 => Math.Log(draw.Sigma2),
        1 => draw.Latents!.Skew,
        2 => draw.Latents!.ValueArray.Average(Math.Log),
        _ => draw.Structure.TermIds.Length,
    };
}
