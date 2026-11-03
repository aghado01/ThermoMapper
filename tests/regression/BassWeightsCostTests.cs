using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The unit-weight path's cost beside the weighted core: a timing block, the same sweeps on the same seeds, run before and after
/// the change it reads. Since the unit path's draws are bit-identical, block k does the same work in either binary, so blocks
/// pair across them. Each block's line records the load the cost harness reads (<see cref="BassCost.Load"/>): the other test
/// hosts at its ends and the cores other processes kept busy, against the statistical brief's quiet rule. It reads and judges
/// nothing. Run it alone, in Release, behind <c>BASS_WEIGHTS_TIMING</c> = 1. CPU time ticks every 15.6 ms, so a block of about 0.2 s
/// reads to 8%; <c>BASS_WEIGHTS_BLOCK</c> sets a longer one, to read a smaller change.
/// </summary>
public sealed class BassWeightsCostTests
{
    private readonly ITestOutputHelper _out;
    public BassWeightsCostTests(ITestOutputHelper output) => _out = output;

    private const int WarmUp = 2000, Blocks = 9;
    private static readonly int Block = int.TryParse(Environment.GetEnvironmentVariable("BASS_WEIGHTS_BLOCK"), out int sweeps) ? sweeps : 20_000;

    // The blocks of `chain`, each timed apart from its load's reading, one line each.
    private void Timed(string label, BassChain chain)
    {
        for (int s = 0; s < WarmUp; s++) chain.Step();
        for (int k = 0; k < Blocks; k++)
        {
            var load = new BassCost.Load();
            TimeSpan start = Environment.CpuUsage.UserTime;
            long ticks = Stopwatch.GetTimestamp();
            for (int s = 0; s < Block; s++) chain.Step();
            double wall = Stopwatch.GetElapsedTime(ticks).TotalMicroseconds / Block;
            double cpu = (Environment.CpuUsage.UserTime - start).TotalMicroseconds / Block;
            (string[] atEnd, _, _, double others) = load.End();
            bool quiet = load.AtStart.Length == 0 && atEnd.Length == 0 && others <= 1.5;
            _out.WriteLine(string.Create(CultureInfo.InvariantCulture,
                $"BLOCK\t{label}\t{k}\t{chain.Current.Count}\t{wall:F2}\t{cpu:F2}\t{others:F2}\t{Math.Max(load.AtStart.Length, atEnd.Length)}\t{(quiet ? 1 : 0)}"));
        }
    }

    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Functional)]
    public void UnitPath_TimingBlock(string name)
    {
        if (Environment.GetEnvironmentVariable("BASS_WEIGHTS_TIMING") != "1") return;
        BassMixing.Fixture f = BassMixing.Get(name, null);
        BassModel model = BassModel.Create(f.Config, f.Data);
        Timed(name, new BassChain(model, [], Math.Log(1.0 / f.Data.Response.Observations), new Xoshiro256PlusPlus(11), new Xoshiro256PlusPlus(12)));
    }

    // What each weighted reading costs on the same fixture and seeds: the closed forms on the weighted statistics, and the pencil's
    // decomposition per proposal. The stream differs from the unit path's, so M does too; each line states its own.
    [Theory]
    [InlineData(WeightedPriorGram.Weighted)]
    [InlineData(WeightedPriorGram.Unweighted)]
    public void WeightedPath_TimingBlock(WeightedPriorGram gram)
    {
        if (Environment.GetEnvironmentVariable("BASS_WEIGHTS_TIMING") != "1") return;
        BassMixing.Fixture f = BassMixing.Get(BassMixing.Friedman, null);
        var rng = new Xoshiro256PlusPlus(21);
        double[] weights = [.. Enumerable.Range(0, f.Data.Response.Observations).Select(_ => 0.5 + rng.NextDouble())];
        BassModel model = BassModel.Create(f.Config, f.Data, new Reading(weights, gram));
        Timed($"weighted {gram}", new BassChain(model, [], Math.Log(1.0 / f.Data.Response.Observations), new Xoshiro256PlusPlus(11), new Xoshiro256PlusPlus(12)));
    }

    private sealed class Reading(double[] weights, WeightedPriorGram gram) : IBassObservationModel
    {
        public WorkingResponse Working(BassResponse response) => new(response, weights: weights, priorGram: gram);

        public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => [];
    }
}
