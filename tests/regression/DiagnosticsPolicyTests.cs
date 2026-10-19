using System;
using System.Threading;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// B01 (R09) diagnostic failure policy: an undefined R̂ never certifies convergence, chains frozen at different
/// values stay distinct from a constant functional, and the early stop needs enough draws and at least two chains.
/// </summary>
public sealed class DiagnosticsPolicyTests
{
    [Fact]
    public void RHat_ConstantChains_NonRepresentableValue_IsConstant()
    {
        // Running Σ/Σ² read this as a finite estimate (0.9995); streaming moments keep it exactly constant.
        RHatDiagnostic r = ChainDiagnostics.RHat(Chains(4, 1000, (_, _) => 0.1));
        Assert.Equal(RHatState.Constant, r.State);
        Assert.True(r.MeetsTarget(1.01));
    }

    [Fact]
    public void RHat_ChainsFrozenAtDifferentValues_IsFrozen()
    {
        RHatDiagnostic r = ChainDiagnostics.RHat(Chains(4, 1000, (chain, _) => 0.1 * (chain + 1)));
        Assert.Equal(RHatState.Frozen, r.State);
        Assert.False(r.MeetsTarget(1e6));
    }

    [Fact]
    public void RHat_OneChainOrOneDraw_IsUnavailable()
    {
        Assert.Equal(RHatState.Unavailable, ChainDiagnostics.RHat(Chains(1, 100, (_, i) => i)).State);
        Assert.Equal(RHatState.Unavailable, ChainDiagnostics.RHat(Chains(4, 1, (chain, _) => chain)).State);
        Assert.False(RHatDiagnostic.Unavailable.MeetsTarget(1e6));
    }

    [Fact]
    public void RHat_NonFiniteDraw_IsUnavailable()
    {
        RHatDiagnostic r = ChainDiagnostics.RHat(Chains(4, 100, (chain, i) => chain == 2 && i == 50 ? double.NaN : i));
        Assert.Equal(RHatState.Unavailable, r.State);
    }

    [Fact]
    public void RHat_MomentsForm_MatchesHandComputedValue()
    {
        // Four identical chains of {1,2,3,4,5}: W = 2.5, B = 0, R̂ = √((4/5)·W / W) = √0.8.
        RHatDiagnostic r = ChainDiagnostics.RHat(Chains(4, 5, (_, i) => i + 1));
        Assert.Equal(RHatState.Estimated, r.State);
        Assert.Equal(Math.Sqrt(0.8), r.Value, 12);
    }

    [Fact]
    public void Worst_OrdersConstantEstimateFrozenUnavailable()
    {
        RHatDiagnostic[] set = [RHatDiagnostic.Constant, RHatDiagnostic.Estimate(1.02), RHatDiagnostic.Estimate(1.4)];
        Assert.Equal(RHatDiagnostic.Estimate(1.4), RHatDiagnostic.Worst(set));
        Assert.Equal(RHatState.Frozen, RHatDiagnostic.Worst([.. set, RHatDiagnostic.Frozen]).State);
        Assert.Equal(RHatState.Unavailable,
                     RHatDiagnostic.Worst([RHatDiagnostic.Unavailable, RHatDiagnostic.Frozen, .. set]).State);
        Assert.Equal(RHatState.Constant, RHatDiagnostic.Worst([RHatDiagnostic.Constant]).State);
    }

    [Fact]
    public void EarlyStop_ChainsFrozenAtDifferentValues_DoesNotStop()
    {
        var model = new ScriptedModel((index, _) => () => index);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 0, samples: 200,
                                                  rHatTarget: 1.1, batchSize: 10).Summary;
        Assert.Equal(200, run.SamplesUsed);
        Assert.Equal(RHatState.Frozen, run.FunctionalRHat[0].State);
    }

    [Fact]
    public void EarlyStop_TooFewDraws_DoesNotStop()
    {
        var model = new ScriptedModel(Noise);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 0, samples: 400,
                                                  rHatTarget: 1.1, batchSize: 1).Summary;
        Assert.True(run.SamplesUsed >= ChainEnsemble.MinDrawsForStop, $"stopped after {run.SamplesUsed} draw(s)");
    }

    [Fact]
    public void EarlyStop_NonFiniteDraw_DoesNotStop()
    {
        var model = new ScriptedModel((index, seed) =>
        {
            Func<double> noise = Noise(index, seed);
            return () => index == 0 ? double.NaN : noise();
        });
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 0, samples: 200,
                                                  rHatTarget: 1.1, batchSize: 10).Summary;
        Assert.Equal(200, run.SamplesUsed);
        Assert.Equal(RHatState.Unavailable, run.FunctionalRHat[0].State);
    }

    [Fact]
    public void EarlyStop_ConstantFunctional_DoesNotBlockStop()
    {
        var model = new ScriptedModel((_, _) => () => 0.1);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 0, samples: 400,
                                                  rHatTarget: 1.1, batchSize: 10).Summary;
        Assert.Equal(ChainEnsemble.MinDrawsForStop, run.SamplesUsed);
        Assert.Equal(RHatState.Constant, run.FunctionalRHat[0].State);
    }

    [Fact]
    public void EarlyStop_SingleChain_IsRejected()
    {
        var model = new ScriptedModel(Noise);
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new ChainEnsemble().Run(model, chains: 1, masterSeed: 1, burn: 0, samples: 200,
                                    rHatTarget: 1.1, batchSize: 10));
    }

    // c chains of n draws, value(chain, draw), folded into streaming moments.
    private static RunningMoments[] Chains(int c, int n, Func<int, int, double> value)
    {
        var chains = new RunningMoments[c];
        for (int ci = 0; ci < c; ci++)
            for (int i = 0; i < n; i++) chains[ci].Add(value(ci, i));
        return chains;
    }

    private static Func<double> Noise(int index, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        return () =>
        {
            double u1 = 1.0 - rng.NextDouble();
            double u2 = 1.0 - rng.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        };
    }

    // One scalar functional per chain; the source maps (chain index, chain seed) to that chain's draw generator.
    private sealed class ScriptedModel : IEnsembleModel<double>
    {
        private readonly Func<int, int, Func<double>> _source;
        private int _next = -1;

        public ScriptedModel(Func<int, int, Func<double>> source) => _source = source;

        public int FunctionalDim => 1;
        public int EssDim => 0;

        public IChainRun<double> StartChain(int chainSeed, int readoutSeed)
            => new ScriptedRun(new ScriptedChain(_source(Interlocked.Increment(ref _next), chainSeed)));
    }

    private sealed class ScriptedChain : IChain<double>
    {
        private readonly Func<double> _next;
        public ScriptedChain(Func<double> next) => _next = next;
        public double Step() => _next();
        public long Accepted => 0;
        public long Attempts => 0;
    }

    private sealed class ScriptedRun : IChainRun<double>
    {
        public ScriptedRun(ScriptedChain chain) => Chain = chain;
        public IChain<double> Chain { get; }
        public void Burn(int steps) { }
        public void Accumulate(in double draw, Span<double> functionals) => functionals[0] = draw;
    }
}
