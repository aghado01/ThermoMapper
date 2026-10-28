using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's checkpoints (stage 4's S7): a run resumed from a checkpoint is bit-identical to the uninterrupted run under every
/// option, resumed once or twice, and a checkpoint stays reusable; a checkpoint taken on another config or data is
/// refused. Beneath it, a chain restored from its snapshot continues bit for bit, with its design's statistics, which a
/// strict reload of the same structure would not reproduce.
/// </summary>
public sealed class BassRunResumeTests
{
    private readonly ITestOutputHelper _out;
    public BassRunResumeTests(ITestOutputHelper output) => _out = output;

    private static BassConfig Config(string option, int samples)
    {
        BassConfig config = BassRunTests.Config(chains: 3, burn: 800, samples: samples, maxBasis: 8);
        return option switch
        {
            "fixed g" => config with { Coefficients = new FixedG() },
            "cancelling" => config with { Sampling = config.Sampling with { Moves = new PriorCancelling(0.4) } },
            "jeffreys, flat" => config with { Sigma2 = Sigma2Prior.Jeffreys, Intercept = GPriorCentre.Flat },
            "thinned" => config with { Sampling = config.Sampling with { Thin = 3 } },
            "batches" => config with { Sampling = config.Sampling with { BatchSize = 70 } },
            _ => config,
        };
    }

    private static BassData Data(string option) => option switch
    {
        "functional" => BassRunTests.FunctionalData(),
        "categorical" => BassPriorSamplerTests.MixedData(),
        _ => BassSweepTests.Data(large: true),
    };

    [Theory]
    [InlineData("defaults")]
    [InlineData("functional")]
    [InlineData("categorical")]
    [InlineData("fixed g")]
    [InlineData("cancelling")]
    [InlineData("jeffreys, flat")]
    [InlineData("thinned")]
    [InlineData("batches")]
    public void Resume_IsBitIdentical(string option)
    {
        BassData data = Data(option);
        BassFit whole = Bass.Fit(Config(option, 600), data);
        BassConfig half = Config(option, 300);
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        Assert.NotNull(first.Checkpoint);
        Assert.Equal(300, first.Checkpoint!.DrawsPerChain);
        Assert.Equal(first.Identity, first.Checkpoint.Identity);
        BassFit resumed = Bass.Resume(first.Checkpoint, half, data);
        Assert.Equal(600, resumed.Diagnostics.SamplesUsed);
        Assert.False(resumed.Diagnostics.StoppedEarly);
        BassRunTests.AssertSameFit(whole, resumed);
        _out.WriteLine($"{option}: {whole.Terms.Count} terms, acceptance {whole.Diagnostics.AcceptanceRate:F3}, " +
                       $"{whole.Diagnostics.Moves.Sum(m => m.Accepted)} moves accepted");
    }

    [Fact]
    public void Resume_Twice_AndFromTheSameCheckpointAgain()
    {
        BassData data = Data("defaults");
        BassConfig third = Config("defaults", 200);
        BassFit whole = Bass.Fit(Config("defaults", 600), data);
        BassFit first = Bass.Fit(third, data, new BassRunOptions { Checkpoint = true });
        BassFit second = Bass.Resume(first.Checkpoint!, third, data, new BassRunOptions { Checkpoint = true });
        Assert.Equal(400, second.Checkpoint!.DrawsPerChain);
        BassRunTests.AssertSameFit(whole, Bass.Resume(second.Checkpoint, third, data));

        // The checkpoint is restart state beside the posterior: resuming does not consume it, and neither run changes
        // the fit it came with.
        double before = first.Draws[1][199].Sigma2;
        BassRunTests.AssertSameFit(second, Bass.Resume(first.Checkpoint!, third, data, new BassRunOptions { Checkpoint = true }));
        Assert.Equal(200, first.Draws[1].Count);
        Assert.Equal(before, first.Draws[1][199].Sigma2);
    }

    [Fact]
    public void Resume_RefusesAnotherRun()
    {
        BassData data = Data("defaults");
        BassConfig config = Config("defaults", 50) with { Sampling = Config("defaults", 50).Sampling with { BurnIn = 20 } };
        var checkpoint = Bass.Fit(config, data, new BassRunOptions { Checkpoint = true }).Checkpoint!;

        var seed = Assert.Throws<ArgumentException>(() => Bass.Resume(checkpoint, config with { Sampling = config.Sampling with { MasterSeed = 12 } }, data));
        Assert.Contains("fingerprint", seed.Message);
        Assert.Equal("checkpoint", seed.ParamName);
        Assert.Throws<ArgumentException>(() => Bass.Resume(checkpoint, config with { Count = new BassCountPrior(5, 5) }, data));
        double[] y = [.. ((ScalarResponse)data.Response).Y.ToArray()];
        y[0] += 1e-12;
        Assert.Throws<ArgumentException>(() => Bass.Resume(checkpoint, config, new BassData(data.Predictors, new ScalarResponse(y))));
        Assert.Throws<ArgumentNullException>(() => Bass.Resume(null!, config, data));
        Assert.Null(Bass.Fit(config, data).Checkpoint);
    }

    [Fact]
    public void Resume_ReportsItsOwnSweeps()
    {
        BassData data = Data("defaults");
        BassConfig config = Config("thinned", 40);
        var checkpoint = Bass.Fit(config, data, new BassRunOptions { Checkpoint = true }).Checkpoint!;
        var reports = new List<BassProgress>();
        Bass.Resume(checkpoint, config, data, new BassRunOptions { Progress = new Inline(p => { lock (reports) reports.Add(p); }), ProgressInterval = 100 });
        Assert.Equal(new BassProgress(3 * 40 * 3, 3 * 40 * 3), reports[^1]);   // no burn-in again
    }

    [Theory]
    [InlineData("defaults")]
    [InlineData("functional")]
    public void Chain_RestoredFromItsSnapshot_ContinuesBitForBit(string option)
    {
        BassConfig config = Config(option, 1);
        var model = BassModel.Create(config, Data(option));
        BassState start = new BassPriorSampler(model).Draw(new Xoshiro256PlusPlus(3)).State;
        var live = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(4), new Xoshiro256PlusPlus(5));
        for (int i = 0; i < 1500; i++) live.Step();
        BassChainSnapshot snapshot = live.Snapshot();

        // A strict reload of the same structure: the statistics agree to rounding, not to the bit.
        BassDesign reloaded = model.NewDesign();
        Assert.Null(reloaded.Load(live.Current.Terms));
        BassDesign design = live.Structure.DesignAt(live.Current);
        _out.WriteLine($"{option}: M = {live.Current.Count}; a strict reload's ‖r‖² differs by {reloaded.ResidualSquares - design.ResidualSquares:G3}, " +
                       $"its coefficients by up to {Enumerable.Range(0, design.Count).Max(j => Math.Abs(reloaded.Coefficients[j] - design.Coefficients[j])):G3}");

        BassChain restored = BassChain.Restore(model, snapshot), again = BassChain.Restore(model, snapshot);
        BassDesign restoredDesign = restored.Structure.DesignAt(restored.Current);
        Assert.Equal(design.ResidualSquares, restoredDesign.ResidualSquares);
        Assert.Equal(design.FittedSquares, restoredDesign.FittedSquares);
        Assert.Equal(design.Coefficients.ToArray(), restoredDesign.Coefficients.ToArray());
        Assert.Equal(design.Residual.ToArray(), restoredDesign.Residual.ToArray());

        var liveComposition = new Xoshiro256PlusPlus(6);
        var restoredComposition = new Xoshiro256PlusPlus(6);
        BassTermTable liveTable = new(), restoredTable = new();
        for (int i = 0; i < 1500; i++)
        {
            BassState a = live.Step(), b = restored.Step();
            Assert.Equal(a.Terms, b.Terms);
            Assert.Equal(a.LogTau, b.LogTau);
            Assert.Equal(a.Fit, b.Fit);
            if (i % 10 != 0) continue;
            BassDraw x = live.Compose(liveComposition, liveTable), y = restored.Compose(restoredComposition, restoredTable);
            Assert.Equal(x.Sigma2, y.Sigma2);
            Assert.Equal(x.Coefficients, y.Coefficients);
        }
        Assert.Equal(live.Moves.MoveStats(), restored.Moves.MoveStats());
        Assert.Equal(live.Scan.BlockStats(), restored.Scan.BlockStats());
        Assert.Equal(live.Structure.Refusals.Entries, restored.Structure.Refusals.Entries);
        Assert.Equal(live.Structure.Support.Entries, restored.Structure.Support.Entries);
        Assert.Equal(live.Structure.EmptyBirths, restored.Structure.EmptyBirths);
        // The snapshot is unchanged by the restored chain's run: a chain restored now runs as one restored before.
        BassChain later = BassChain.Restore(model, snapshot);
        for (int i = 0; i < 200; i++) Assert.Equal(again.Step().Fit, later.Step().Fit);
    }

    [Fact]
    public void Restore_RefusesAForeignState()
    {
        BassConfig config = Config("defaults", 1);
        var model = BassModel.Create(config, Data("defaults"));
        var ensemble = new BassEnsemble(model, Bass.Fingerprint(config, Data("defaults")), new BassRunOptions(), 0, 1);
        Assert.Throws<ArgumentException>(() => ensemble.RestoreChain(new Foreign()));

        var chain = new BassChain(model, [], Math.Log(1.0 / 75), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        for (int i = 0; i < 300; i++) chain.Step();
        BassChainSnapshot snapshot = chain.Snapshot();
        // A state the snapshot's design does not hold.
        BassState other = new BassPriorSampler(model).Draw(new Xoshiro256PlusPlus(9)).State;
        Assert.False(other.Terms.SequenceEqual(chain.Current.Terms));
        Assert.Throws<ArgumentException>(() => new BassStructureMoves(model, other, snapshot.Structure));
        // A design over other data.
        var small = BassModel.Create(BassSweepTests.Config(), BassSweepTests.Data());
        Assert.Throws<ArgumentException>(() => small.NewDesign().Restore(snapshot.Structure.Design));
    }

    private sealed class Foreign : Maths.Samplers.Ensemble.IChainState { }

    private sealed class Inline(Action<BassProgress> report) : IProgress<BassProgress>
    {
        public void Report(BassProgress value) => report(value);
    }
}
