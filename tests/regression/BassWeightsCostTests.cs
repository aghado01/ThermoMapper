using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The unit-weight path's cost beside the weighted core: a small timing block, the same sweeps on the same seeds, run
/// before and after the change it reads. It reads and judges nothing. Run it alone, in Release, behind
/// <c>BASS_WEIGHTS_TIMING</c> = 1.
/// </summary>
public sealed class BassWeightsCostTests
{
    private readonly ITestOutputHelper _out;
    public BassWeightsCostTests(ITestOutputHelper output) => _out = output;

    private const int WarmUp = 2000, Block = 2000, Blocks = 9;

    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Functional)]
    public void UnitPath_TimingBlock(string name)
    {
        if (Environment.GetEnvironmentVariable("BASS_WEIGHTS_TIMING") != "1") return;
        BassMixing.Fixture f = BassMixing.Get(name, null);
        BassModel model = BassModel.Create(f.Config, f.Data);
        var chain = new BassChain(model, [], Math.Log(1.0 / f.Data.Response.Observations), new Xoshiro256PlusPlus(11), new Xoshiro256PlusPlus(12));
        for (int s = 0; s < WarmUp; s++) chain.Step();
        var wall = new double[Blocks];
        var cpu = new double[Blocks];
        for (int k = 0; k < Blocks; k++)
        {
            TimeSpan start = Environment.CpuUsage.UserTime;
            long ticks = Stopwatch.GetTimestamp();
            for (int s = 0; s < Block; s++) chain.Step();
            wall[k] = Stopwatch.GetElapsedTime(ticks).TotalMicroseconds / Block;
            cpu[k] = (Environment.CpuUsage.UserTime - start).TotalMicroseconds / Block;
        }
        Array.Sort(wall);
        Array.Sort(cpu);
        _out.WriteLine($"{name}: M = {chain.Current.Count}, µs per sweep: wall min {wall[0]:F1}, median {wall[Blocks / 2]:F1}; user CPU min {cpu[0]:F1}, median {cpu[Blocks / 2]:F1}");
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
        var chain = new BassChain(model, [], Math.Log(1.0 / f.Data.Response.Observations), new Xoshiro256PlusPlus(11), new Xoshiro256PlusPlus(12));
        for (int s = 0; s < WarmUp; s++) chain.Step();
        var wall = new double[Blocks];
        for (int k = 0; k < Blocks; k++)
        {
            long ticks = Stopwatch.GetTimestamp();
            for (int s = 0; s < Block; s++) chain.Step();
            wall[k] = Stopwatch.GetElapsedTime(ticks).TotalMicroseconds / Block;
        }
        Array.Sort(wall);
        _out.WriteLine($"Friedman, weighted, prior on the {gram} columns: M = {chain.Current.Count}, µs per sweep: wall min {wall[0]:F1}, median {wall[Blocks / 2]:F1}");
    }

    private sealed class Reading(double[] weights, WeightedPriorGram gram) : IBassObservationModel
    {
        public WorkingResponse Working(BassResponse response) => new(response, weights: weights, priorGram: gram);

        public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => [];
    }
}
