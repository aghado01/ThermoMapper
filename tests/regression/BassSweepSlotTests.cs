using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The sweep's slots (F-X2): D04's robust order, structure, τ, a, latent scales, then σ² and the skew, is the sweep's
/// general form, and Gaussian BASS fills the first two.
/// <list type="bullet">
/// <item>The slots are written in one place and in D04's order, and the chain, the tempering replicas and the SMC kernel read it.</item>
/// <item>Gaussian BASS's sweep is its two blocks, three with nothing else, two with τ held; its other slots are empty, and an
/// observation model that fills them explicitly with nothing walks the same chain, stream and all, as Gaussian.</item>
/// <item>A model that fills the later slots has its blocks run after τ's, in slot order, on the scan's stream: the first sweep
/// agrees with Gaussian's, which the later blocks cannot reach, and later sweeps do not once a block draws.</item>
/// </list>
/// </summary>
public sealed class BassSweepSlotTests
{
    private sealed class Probe(string name, List<string> log, bool draw) : IBlock<ReversibleJumpChain<BassState>>
    {
        public string Name => name;

        public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
        {
            log.Add(name);
            if (draw) rng.NextDouble();
            return true;
        }
    }

    private sealed class Filling(List<string> log, bool draw) : IBassObservationModel
    {
        public WorkingResponse Working(BassResponse response) => new(response);

        public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot)
            => slot switch
            {
                BassSweepSlot.Coefficients => [new Probe("a", log, draw)],
                BassSweepSlot.Latent => [new Probe("latent-1", log, draw), new Probe("latent-2", log, draw)],
                BassSweepSlot.Scale => [new Probe("scale", log, draw)],
                _ => throw new ArgumentOutOfRangeException(nameof(slot)),
            };
    }

    private sealed class Empty : IBassObservationModel
    {
        public WorkingResponse Working(BassResponse response) => new(response);

        public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => [];
    }

    private static BassConfig Sampled() => SensitivityFixture.Config();

    private static BassConfig Held() => SensitivityFixture.Config() with { Coefficients = BassCoefficientPrior.G(TauPrior.Held(0.3)) };

    [Fact]
    public void TheSlots_AreD04sOrder()
    {
        Assert.Equal([BassSweepSlot.Structure, BassSweepSlot.Tau, BassSweepSlot.Coefficients, BassSweepSlot.Latent, BassSweepSlot.Scale],
                     BassSweepSlots.Order);
        Assert.Equal(BassSweepSlots.Order.OrderBy(s => s), BassSweepSlots.Order);
        foreach (BassSweepSlot slot in new[] { BassSweepSlot.Coefficients, BassSweepSlot.Latent, BassSweepSlot.Scale })
            Assert.Empty(GaussianObservation.Instance.Blocks(slot));
    }

    [Fact]
    public void GaussianBass_FillsTheStructureAndTauSlots_AndNothingElse()
    {
        BassData data = SensitivityFixture.Scalar();
        BassModel sampled = BassModel.Create(Sampled(), data);
        List<IBlock<ReversibleJumpChain<BassState>>> sweep = BassSweepSlots.Blocks(sampled);
        Assert.Equal(2, sweep.Count);
        Assert.IsType<MoveStep<BassState>>(sweep[0]);
        Assert.IsType<BassTauBlock>(sweep[1]);
        Assert.Single(BassSweepSlots.Blocks(sampled, BassSweepSlot.Tau));

        BassModel held = BassModel.Create(Held(), data);
        Assert.False(held.Target.SamplesTau);
        Assert.Single(BassSweepSlots.Blocks(held));
        Assert.Empty(BassSweepSlots.Blocks(held, BassSweepSlot.Tau));
    }

    [Fact]
    public void AModelThatFillsTheLaterSlots_RunsThemAfterTau_InSlotOrder()
    {
        BassData data = SensitivityFixture.Scalar();
        var log = new List<string>();
        BassModel model = BassModel.Create(Sampled(), data, new Filling(log, draw: false));
        List<IBlock<ReversibleJumpChain<BassState>>> sweep = BassSweepSlots.Blocks(model);
        Assert.IsType<MoveStep<BassState>>(sweep[0]);
        Assert.IsType<BassTauBlock>(sweep[1]);
        Assert.Equal(["a", "latent-1", "latent-2", "scale"], sweep.Skip(2).Select(b => b.Name));
        Assert.Equal(["tau", "a", "latent-1", "latent-2", "scale"], BassSweepSlots.Blocks(model, BassSweepSlot.Tau).Select(b => b.Name));

        var chain = new BassChain(model, [], Math.Log(0.05), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        for (int step = 0; step < 5; step++) chain.Step();
        Assert.Equal(Enumerable.Repeat(new[] { "a", "latent-1", "latent-2", "scale" }, 5).SelectMany(x => x), log);

        // Held τ: the later slots stay, behind the moves alone.
        BassModel heldModel = BassModel.Create(Held(), data, new Filling(new List<string>(), draw: false));
        Assert.Equal(["a", "latent-1", "latent-2", "scale"], BassSweepSlots.Blocks(heldModel).Skip(1).Select(b => b.Name));
    }

    [Fact]
    public void AnEmptySlot_DrawsNothing_AndConsumesNoRandomNumber()
    {
        BassData data = SensitivityFixture.Scalar();
        BassModel gaussian = BassModel.Create(Sampled(), data);
        BassModel empty = BassModel.Create(Sampled(), data, new Empty());
        var a = new BassChain(gaussian, [], Math.Log(0.05), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        var b = new BassChain(empty, [], Math.Log(0.05), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        for (int step = 0; step < 100; step++)
        {
            BassState sa = a.Step(), sb = b.Step();
            Assert.Equal(sa.LogTau, sb.LogTau);
            Assert.Equal(sa.Fit.ResidualSquares, sb.Fit.ResidualSquares);
        }
        BassChainSnapshot ca = a.Snapshot(), cb = b.Snapshot();
        Assert.Equal(ca.BlockRng, cb.BlockRng);                 // the scan's stream, word for word
        Assert.Equal((ca.Chain.Rng0, ca.Chain.Rng1, ca.Chain.Rng2, ca.Chain.Rng3), (cb.Chain.Rng0, cb.Chain.Rng1, cb.Chain.Rng2, cb.Chain.Rng3));   // and the moves'
        Assert.Equal(ca.BlockAttempts, cb.BlockAttempts);       // and no block ran that Gaussian lacks
    }

    [Fact]
    public void TheTemperingReplicas_RunTheLaterSlotsAtEveryRung()
    {
        var log = new List<string>();
        var spec = new TemperingSpec { Levels = 3, HotEnd = 0.2 };
        BassModel model = BassModel.Create(Sampled(), SensitivityFixture.Scalar(), new Filling(log, draw: false));
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, 2);
        ladder.Step();
        Assert.True(log.Count >= 3 * 4 && log.Count % 4 == 0, $"{log.Count} probe calls");   // each rung's replica: a, two latent blocks, scale
        Assert.All(log.Chunk(4), call => Assert.Equal(["a", "latent-1", "latent-2", "scale"], call));
    }

    [Fact]
    public void ABlockInALaterSlot_SharesTheScansStream_BehindTau()
    {
        BassData data = SensitivityFixture.Scalar();
        BassModel gaussian = BassModel.Create(Sampled(), data);
        BassModel filled = BassModel.Create(Sampled(), data, new Filling(new List<string>(), draw: true));
        var a = new BassChain(gaussian, [], Math.Log(0.05), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        var b = new BassChain(filled, [], Math.Log(0.05), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        BassState first = a.Step(), firstFilled = b.Step();
        Assert.Equal(first.LogTau, firstFilled.LogTau);         // the later blocks run after τ and change nothing here
        Assert.Equal(first.Count, firstFilled.Count);
        bool diverged = false;
        for (int step = 0; step < 20 && !diverged; step++) diverged = a.Step().LogTau != b.Step().LogTau;
        Assert.True(diverged, "the filled slots' draws should move the stream τ's next step reads");
    }
}
