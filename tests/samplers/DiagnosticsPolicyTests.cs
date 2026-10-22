using System;
using System.Threading;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Samplers.Tests;

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
        Assert.Equal(DiagnosticState.Constant, r.State);
        Assert.True(r.MeetsTarget(1.01));
    }

    [Fact]
    public void RHat_ChainsFrozenAtDifferentValues_IsFrozen()
    {
        RHatDiagnostic r = ChainDiagnostics.RHat(Chains(4, 1000, (chain, _) => 0.1 * (chain + 1)));
        Assert.Equal(DiagnosticState.Frozen, r.State);
        Assert.False(r.MeetsTarget(1e6));
    }

    [Fact]
    public void RHat_OneChainOrOneDraw_IsUnavailable()
    {
        Assert.Equal(DiagnosticState.Unavailable, ChainDiagnostics.RHat(Chains(1, 100, (_, i) => i)).State);
        Assert.Equal(DiagnosticState.Unavailable, ChainDiagnostics.RHat(Chains(4, 1, (chain, _) => chain)).State);
        Assert.False(RHatDiagnostic.Unavailable.MeetsTarget(1e6));
    }

    [Fact]
    public void RHat_NonFiniteDraw_IsUnavailable()
    {
        RHatDiagnostic r = ChainDiagnostics.RHat(Chains(4, 100, (chain, i) => chain == 2 && i == 50 ? double.NaN : i));
        Assert.Equal(DiagnosticState.Unavailable, r.State);
    }

    [Fact]
    public void RHat_MomentsForm_MatchesHandComputedValue()
    {
        // Four identical chains of {1,2,3,4,5}: W = 2.5, B = 0, R̂ = √((4/5)·W / W) = √0.8.
        RHatDiagnostic r = ChainDiagnostics.RHat(Chains(4, 5, (_, i) => i + 1));
        Assert.Equal(DiagnosticState.Estimated, r.State);
        Assert.Equal(Math.Sqrt(0.8), r.Value, 12);
    }

    [Fact]
    public void RunningMoments_Merge_MatchesFoldingTheConcatenation()
    {
        double[] a = [1.5, -2.0, 3.25, 0.5];
        double[] b = [10.0, 11.5, 9.0];
        RunningMoments merged = RunningMoments.Merge(RunningMoments.Of(a), RunningMoments.Of(b));
        RunningMoments whole = RunningMoments.Of([.. a, .. b]);

        Assert.Equal(whole.Count, merged.Count);
        Assert.Equal(whole.Mean, merged.Mean, 12);
        Assert.Equal(whole.SumSquaredDeviations, merged.SumSquaredDeviations, 10);
        Assert.Equal(whole.Mean, RunningMoments.Merge(RunningMoments.Of(b), RunningMoments.Of(a)).Mean, 12);
        Assert.Equal(RunningMoments.Of(a).Mean, RunningMoments.Merge(RunningMoments.Of(a), new RunningMoments()).Mean);
    }

    [Fact]
    public void Worst_OrdersConstantEstimateFrozenUnavailable()
    {
        RHatDiagnostic[] set = [RHatDiagnostic.Constant, RHatDiagnostic.Estimate(1.02), RHatDiagnostic.Estimate(1.4)];
        Assert.Equal(RHatDiagnostic.Estimate(1.4), RHatDiagnostic.Worst(set));
        Assert.Equal(DiagnosticState.Frozen, RHatDiagnostic.Worst([.. set, RHatDiagnostic.Frozen]).State);
        Assert.Equal(DiagnosticState.Unavailable,
                     RHatDiagnostic.Worst([RHatDiagnostic.Unavailable, RHatDiagnostic.Frozen, .. set]).State);
        Assert.Equal(DiagnosticState.Constant, RHatDiagnostic.Worst([RHatDiagnostic.Constant]).State);
    }

    [Fact]
    public void EarlyStop_ChainsFrozenAtDifferentValues_DoesNotStop()
    {
        var model = new ScriptedModel((index, _) => () => index);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 0, samples: 200,
                                                  rHatTarget: 1.1, batchSize: 10).Summary;
        Assert.Equal(200, run.SamplesUsed);
        Assert.Equal(DiagnosticState.Frozen, run.FunctionalRHat[0].State);
    }

    [Fact]
    public void EarlyStop_WaitsForTheEssGate()
    {
        // Independent draws meet an R̂ target almost at once; the stop still waits until the bulk and tail ESS reach
        // MinEssPerChain per chain, about 100 draws each here.
        var model = new ScriptedModel(Noise);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 0, samples: 400,
                                                  rHatTarget: 1.1, batchSize: 1).Summary;
        Assert.True(run.SamplesUsed < 400, "independent draws should stop early");
        Assert.True(run.SamplesUsed >= 60, $"stopped after {run.SamplesUsed} draw(s)");
        Assert.True(run.Monitored[0].Ess.MeetsMinimum(4 * ChainEnsemble.MinEssPerChain));
        Assert.Equal(4, run.Traces.Count);
        Assert.All(run.Traces, t => Assert.Equal(run.SamplesUsed, t.Count));
        Assert.Equal(run.Monitored[0], ChainDiagnostics.Diagnose(run.Traces, 0));
    }

    [Fact]
    public void EarlyStop_StronglyCorrelatedDraws_WaitLonger()
    {
        // AR(1) with φ = 0.95: about 1/39 of the draws are effective, so the gate needs about 3,900 draws per chain.
        var model = new ScriptedModel((_, seed) =>
        {
            Func<double> noise = Noise(0, seed);
            double x = 0.0;
            return () => x = 0.95 * x + Math.Sqrt(1.0 - 0.95 * 0.95) * noise();
        });
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 200, samples: 20_000,
                                                  rHatTarget: 1.1, batchSize: 100).Summary;
        Assert.InRange(run.SamplesUsed, 2_000, 8_000);
    }

    [Fact]
    public void EarlyStop_WithoutAMonitoredFunctional_IsRejected()
    {
        var model = new ScriptedModel(monitored: 0, Noise);
        Assert.Throws<ArgumentException>(() =>
            new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 0, samples: 200, rHatTarget: 1.1));
    }

    [Fact]
    public void AcceptanceRate_CountsTheSamplingPhaseOnly()
    {
        // The chain accepts every proposal while burning in and none while sampling.
        var model = new BurnAcceptingModel();
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 2, masterSeed: 1, burn: 50, samples: 50).Summary;
        Assert.Equal(0.0, run.AcceptanceRate);
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
        Assert.Equal(DiagnosticState.Unavailable, run.FunctionalRHat[0].State);
    }

    [Fact]
    public void EarlyStop_ConstantFunctional_DoesNotBlockStop()
    {
        var model = new ScriptedModel((_, _) => () => 0.1, Noise);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 0, samples: 400,
                                                  rHatTarget: 1.1, batchSize: 10).Summary;
        Assert.True(run.SamplesUsed < 400, "a constant functional beside a converged one should not block the stop");
        Assert.Equal(DiagnosticState.Constant, run.FunctionalRHat[0].State);
        Assert.Equal(DiagnosticState.Estimated, run.FunctionalRHat[1].State);
    }

    [Fact]
    public void EarlyStop_OnlyConstantFunctionals_RunsToTheEnd()
    {
        var model = new ScriptedModel((_, _) => () => 0.1, (_, _) => () => 2.5);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 1, burn: 0, samples: 400,
                                                  rHatTarget: 1.1, batchSize: 10).Summary;
        Assert.Equal(400, run.SamplesUsed);
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

    // One functional per source; each source maps (chain index, chain seed) to that chain's draw generator.
    private sealed class ScriptedModel : IEnsembleModel<double[]>
    {
        private readonly Func<int, int, Func<double>>[] _sources;
        private int _next = -1;

        public ScriptedModel(params Func<int, int, Func<double>>[] sources) : this(sources.Length, sources) { }

        public ScriptedModel(int monitored, params Func<int, int, Func<double>>[] sources)
        {
            _sources = sources;
            MonitoredDim = monitored;
        }

        public int FunctionalDim => _sources.Length;
        public int MonitoredDim { get; }

        public IChainRun<double[]> StartChain(int chainSeed, int readoutSeed)
        {
            int index = Interlocked.Increment(ref _next);
            var generators = new Func<double>[_sources.Length];
            for (int f = 0; f < _sources.Length; f++) generators[f] = _sources[f](index, chainSeed);
            return new ScriptedRun(new ScriptedChain(generators));
        }
    }

    private sealed class ScriptedChain : IChain<double[]>
    {
        private readonly Func<double>[] _generators;
        public ScriptedChain(Func<double>[] generators) => _generators = generators;
        public double[] Step()
        {
            var draw = new double[_generators.Length];
            for (int f = 0; f < draw.Length; f++) draw[f] = _generators[f]();
            return draw;
        }
        public long Accepted => 0;
        public long Attempts => 0;
    }

    private sealed class ScriptedRun : IChainRun<double[]>
    {
        public ScriptedRun(ScriptedChain chain) => Chain = chain;
        public IChain<double[]> Chain { get; }
        public void Burn(int steps) { }
        public void Accumulate(in double[] draw, Span<double> functionals) => draw.CopyTo(functionals);
    }

    private sealed class BurnAcceptingModel : IEnsembleModel<double>
    {
        public int FunctionalDim => 1;
        public int MonitoredDim => 1;
        public IChainRun<double> StartChain(int chainSeed, int readoutSeed) => new BurnAcceptingRun();
    }

    private sealed class BurnAcceptingRun : IChainRun<double>, IChain<double>
    {
        private bool _burning;
        private double _x;
        public IChain<double> Chain => this;
        public long Accepted { get; private set; }
        public long Attempts { get; private set; }

        public double Step()
        {
            Attempts++;
            if (_burning) Accepted++;
            return _x += 1.0;
        }

        public void Burn(int steps)
        {
            _burning = true;
            for (int i = 0; i < steps; i++) Step();
            _burning = false;
        }

        public void Accumulate(in double draw, Span<double> functionals) => functionals[0] = draw;
    }
}
