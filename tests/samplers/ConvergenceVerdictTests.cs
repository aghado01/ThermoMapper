using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// D10's convergence verdict (N14): read at the end of every ensemble run by B01's rules applied to each of the model's
/// groups, with a stopping scope naming the groups that gate an R̂ stop. Synthetic ensembles whose chains sit at two
/// levels, constant and frozen functionals, one chain, malformed declarations and scopes, and the rule read by the stop
/// against the verdict on the same traces.
/// </summary>
public sealed class ConvergenceVerdictTests
{
    [Fact]
    public void ChainsAtTwoLevels_FailTheirGroupAlone_AndNameItsFunctionals()
    {
        var model = new GroupedModel(masterSeed: 3, chains: 4, monitored: 3,
            [new FunctionalGroup("steady", [0]), new FunctionalGroup("split", [1, 2])],
            ["a", "b", "c"], Noise, Levels(3.0), Noise);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 3, burn: 0, samples: 1000).Summary;
        ConvergenceVerdict v = run.Verdict;

        Assert.Equal(ConvergenceOutcome.NotConverged, v.Outcome);
        Assert.Equal(ChainEnsemble.ReferenceRHatTarget, v.Target);
        Assert.Equal(ChainEnsemble.MinEssPerChain, v.MinEssPerChain);
        Assert.Equal(new[] { "steady", "split" }, v.Groups.Select(g => g.Name));
        Assert.All(v.Groups, g => Assert.True(g.GatesStop));

        GroupVerdict steady = v.Group("steady");
        Assert.Equal(ConvergenceOutcome.Converged, steady.Outcome);
        Assert.Empty(steady.Failing);

        GroupVerdict split = v.Group("split");
        Assert.Equal(ConvergenceOutcome.NotConverged, split.Outcome);
        FunctionalReading failing = Assert.Single(split.Failing);
        Assert.Equal(1, failing.Index);
        Assert.Equal("b", failing.Name);
        Assert.Equal(run.Monitored[1].RHat, failing.RHat);
        Assert.Equal(run.Monitored[1].Ess, failing.Ess);
        Assert.Equal(DiagnosticState.Estimated, failing.RHat.State);
        Assert.True(failing.RHat.Value > 1.5, $"R̂ {failing.RHat}");
        Assert.Throws<KeyNotFoundException>(() => v.Group("level"));
    }

    [Fact]
    public void ADefaultModel_ReadsOneGroup_NamedAll_ByIndexNames()
    {
        var model = new GroupedModel(masterSeed: 1, chains: 3, monitored: 2, groups: null, names: null, Noise, Levels(5.0));
        ConvergenceVerdict v = new ChainEnsemble().Run(model, chains: 3, masterSeed: 1, burn: 0, samples: 400).Summary.Verdict;
        GroupVerdict all = Assert.Single(v.Groups);
        Assert.Equal(FunctionalGroup.AllName, all.Name);
        Assert.Equal(ConvergenceOutcome.NotConverged, all.Outcome);
        Assert.Equal("f1", Assert.Single(all.Failing).Name);
    }

    [Fact]
    public void ConstantFunctionals_PassButCannotCertify()
    {
        // A constant beside an estimate converges; a group of constants alone is undetermined, and so is a run whose every
        // gating group is.
        var mixed = new GroupedModel(masterSeed: 2, chains: 4, monitored: 3,
            [new FunctionalGroup("with", [0, 1]), new FunctionalGroup("alone", [2])],
            ["k", "x", "k2"], Constant(0.1), Noise, Constant(2.5));
        ConvergenceVerdict v = new ChainEnsemble().Run(mixed, chains: 4, masterSeed: 2, burn: 0, samples: 400).Summary.Verdict;
        Assert.Equal(ConvergenceOutcome.Converged, v.Group("with").Outcome);
        Assert.Equal(ConvergenceOutcome.Undetermined, v.Group("alone").Outcome);
        Assert.Empty(v.Group("alone").Failing);
        Assert.Equal(ConvergenceOutcome.Converged, v.Outcome);

        var constants = new GroupedModel(masterSeed: 2, chains: 4, monitored: 2, groups: null, names: null, Constant(0.1), Constant(2.5));
        ConvergenceVerdict c = new ChainEnsemble().Run(constants, chains: 4, masterSeed: 2, burn: 0, samples: 400).Summary.Verdict;
        Assert.Equal(ConvergenceOutcome.Undetermined, c.Outcome);
        Assert.Equal(ConvergenceOutcome.Undetermined, Assert.Single(c.Groups).Outcome);

        // A gating group of constants beside a non-gating converged one leaves the run undetermined.
        ConvergenceVerdict scoped = new ChainEnsemble().Run(mixed, chains: 4, masterSeed: 2, burn: 0, samples: 400,
                                                            stoppingScope: ["alone"]).Summary.Verdict;
        Assert.Equal(ConvergenceOutcome.Undetermined, scoped.Outcome);
        Assert.False(scoped.Group("with").GatesStop);
    }

    [Fact]
    public void FrozenAndUnavailableReadings_Fail()
    {
        // Each chain constant at its own level, and a non-finite draw in one chain.
        var model = new GroupedModel(masterSeed: 4, chains: 4, monitored: 1,
            [new FunctionalGroup("frozen", [0]), new FunctionalGroup("broken", [1]), new FunctionalGroup("fine", [2])],
            ["frozen", "broken", "fine"],
            (index, _) => () => index, (index, rng) => { Func<double> n = Noise(index, rng); return () => index == 2 ? double.NaN : n(); },
            Noise);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 4, burn: 0, samples: 400).Summary;
        ConvergenceVerdict v = run.Verdict;

        FunctionalReading frozen = Assert.Single(v.Group("frozen").Failing);
        Assert.Equal(DiagnosticState.Frozen, frozen.RHat.State);
        Assert.Equal(DiagnosticState.Frozen, frozen.Ess!.Value.State);
        FunctionalReading broken = Assert.Single(v.Group("broken").Failing);
        Assert.Equal(DiagnosticState.Unavailable, broken.RHat.State);
        Assert.Null(broken.Ess);   // unmonitored: classic R̂ only
        Assert.Equal(ConvergenceOutcome.Converged, v.Group("fine").Outcome);   // unmonitored, passes on its classic R̂
        Assert.Equal(ConvergenceOutcome.NotConverged, v.Outcome);
    }

    [Fact]
    public void AnUnmonitoredFunctional_FailsOnItsClassicRHat()
    {
        var model = new GroupedModel(masterSeed: 5, chains: 4, monitored: 1, groups: null, names: ["m", "u"], Noise, Levels(2.0));
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 4, masterSeed: 5, burn: 0, samples: 500).Summary;
        FunctionalReading u = Assert.Single(run.Verdict.Groups[0].Failing);
        Assert.Equal("u", u.Name);
        Assert.Equal(run.FunctionalRHat[1], u.RHat);
        Assert.Null(u.Ess);
    }

    [Fact]
    public void OneChain_IsUndetermined()
    {
        var model = new GroupedModel(masterSeed: 6, chains: 1, monitored: 2,
            [new FunctionalGroup("a", [0]), new FunctionalGroup("b", [1])], ["x", "y"], Noise, Noise);
        EnsembleRun run = new ChainEnsemble().Run(model, chains: 1, masterSeed: 6, burn: 0, samples: 500).Summary;
        Assert.Equal(ConvergenceOutcome.Undetermined, run.Verdict.Outcome);
        Assert.All(run.Verdict.Groups, g =>
        {
            Assert.Equal(ConvergenceOutcome.Undetermined, g.Outcome);
            Assert.Empty(g.Failing);
        });
    }

    [Fact]
    public void Scopes_AndDeclarations_AreChecked()
    {
        var model = new GroupedModel(masterSeed: 1, chains: 2, monitored: 1,
            [new FunctionalGroup("m", [0]), new FunctionalGroup("u", [1])], ["x", "y"], Noise, Noise);
        EnsembleRun Execute(GroupedModel m, IReadOnlyCollection<string>? scope, double target = 0.0)
            => new ChainEnsemble().Run(m, chains: 2, masterSeed: 1, burn: 0, samples: 20, rHatTarget: target, stoppingScope: scope).Summary;

        Assert.Throws<ArgumentException>("stoppingScope", () => Execute(model, ["level"]));
        Assert.Throws<ArgumentException>("stoppingScope", () => Execute(model, []));
        Assert.Throws<ArgumentException>("stoppingScope", () => Execute(model, ["u"], target: 1.1));   // no monitored functional gates
        Assert.True(Execute(model, ["m", "m"]).Verdict.Group("m").GatesStop);
        Assert.False(Execute(model, ["m"]).Verdict.Group("u").GatesStop);

        GroupedModel Declared(FunctionalGroup[] groups, string[]? names = null)
            => new(masterSeed: 1, chains: 2, monitored: 1, groups, names ?? ["x", "y"], Noise, Noise);
        Assert.Throws<ArgumentException>("model", () => Execute(Declared([new("a", [0, 1]), new("b", [1])]), null));   // twice
        Assert.Throws<ArgumentException>("model", () => Execute(Declared([new("a", [0])]), null));                       // none
        Assert.Throws<ArgumentException>("model", () => Execute(Declared([new("a", [0]), new("a", [1])]), null));      // a name twice
        Assert.Throws<ArgumentException>("model", () => Execute(Declared([new("a", [0, 2]), new("b", [1])]), null));   // out of range
        Assert.Throws<ArgumentException>("model", () => Execute(Declared([new("a", [0, 1]), new("b", [])]), null));    // empty
        Assert.Throws<ArgumentException>("model", () => Execute(Declared([new(" ", [0, 1])]), null));                  // unnamed
        Assert.Throws<ArgumentException>("model", () => Execute(Declared([new("a", [0, 1])], ["x"]), null));           // a name short
    }

    [Fact]
    public void AScopeExcludingAFailingGroup_CanStop_AndStillReportsIt()
    {
        var model = new GroupedModel(masterSeed: 7, chains: 4, monitored: 2,
            [new FunctionalGroup("good", [0]), new FunctionalGroup("bad", [1])], ["x", "level"], Noise, Levels(3.0));

        EnsembleRun scoped = new ChainEnsemble().Run(model, chains: 4, masterSeed: 7, burn: 0, samples: 2000, rHatTarget: 1.1,
                                                     batchSize: 10, stoppingScope: ["good"]).Summary;
        Assert.True(scoped.SamplesUsed < 2000, "the good group alone should stop the run");
        Assert.Equal(ConvergenceOutcome.Converged, scoped.Verdict.Outcome);
        Assert.Equal(1.1, scoped.Verdict.Target);
        GroupVerdict bad = scoped.Verdict.Group("bad");
        Assert.False(bad.GatesStop);
        Assert.Equal(ConvergenceOutcome.NotConverged, bad.Outcome);
        Assert.Equal("level", Assert.Single(bad.Failing).Name);

        EnsembleRun every = new ChainEnsemble().Run(model, chains: 4, masterSeed: 7, burn: 0, samples: 2000, rHatTarget: 1.1,
                                                    batchSize: 10).Summary;
        Assert.Equal(2000, every.SamplesUsed);
        Assert.Equal(ConvergenceOutcome.NotConverged, every.Verdict.Outcome);
    }

    [Fact]
    public void AnEarlyStop_EndsConverged_AndAnExhaustedBudget_EndsNotConverged()
    {
        var noise = new GroupedModel(masterSeed: 1, chains: 4, monitored: 1, groups: null, names: null, Noise);
        EnsembleRun stopped = new ChainEnsemble().Run(noise, chains: 4, masterSeed: 1, burn: 0, samples: 400, rHatTarget: 1.1,
                                                      batchSize: 10).Summary;
        Assert.True(stopped.SamplesUsed < 400);
        Assert.Equal(ConvergenceOutcome.Converged, stopped.Verdict.Outcome);

        // AR(1) with φ = 0.95 needs about 3,900 draws per chain for the ESS gate; 300 run out of budget unconverged, which
        // a stop flag alone reads the same as no target at all.
        var slow = new GroupedModel(masterSeed: 1, chains: 4, monitored: 1, groups: null, names: null, Ar(0.95));
        EnsembleRun exhausted = new ChainEnsemble().Run(slow, chains: 4, masterSeed: 1, burn: 0, samples: 300, rHatTarget: 1.1,
                                                        batchSize: 50).Summary;
        Assert.Equal(300, exhausted.SamplesUsed);
        Assert.Equal(ConvergenceOutcome.NotConverged, exhausted.Verdict.Outcome);
        FunctionalReading f0 = Assert.Single(exhausted.Verdict.Groups[0].Failing);
        Assert.False(f0.Ess!.Value.MeetsMinimum(4 * ChainEnsemble.MinEssPerChain));
    }

    [Fact]
    public void TheStopsDecision_IsTheVerdictsOnTheSameTraces()
    {
        // Randomized ensembles: AR(1) noise, some with chains offset from each other, two groups, a random scope. The run
        // stops at round k; replayed for j rounds in one, never checked, its verdict against the same target is Converged
        // at j = k and not before, so the stop read exactly the verdict's rule on each round's traces.
        var rng = new Xoshiro256PlusPlus(20261029);
        int stops = 0, runsOut = 0;
        const int batch = 25, budget = 600;
        for (int trial = 0; trial < 24; trial++)
        {
            int chains = 2 + (int)(rng.NextDouble() * 3);
            int seed = (int)(rng.NextDouble() * int.MaxValue);
            double phiA = 0.9 * rng.NextDouble(), phiB = 0.9 * rng.NextDouble();
            double offset = rng.NextDouble() < 0.5 ? 0.0 : 0.6 * rng.NextDouble();
            string[]? scope = (int)(rng.NextDouble() * 3) switch { 0 => null, 1 => ["a"], _ => ["b"] };
            double target = 1.02 + 0.08 * rng.NextDouble();
            var model = new GroupedModel(seed, chains, monitored: 3,
                [new FunctionalGroup("a", [0, 2]), new FunctionalGroup("b", [1])], ["a0", "b0", "a1"],
                Ar(phiA), Shifted(Ar(phiB), offset), Ar(phiA));

            EnsembleRun run = new ChainEnsemble().Run(model, chains, seed, burn: 0, samples: budget, rHatTarget: target, batchSize: batch,
                                                      minEssPerChain: 20, stoppingScope: scope).Summary;
            int rounds = run.SamplesUsed / batch;
            bool stoppedEarly = run.SamplesUsed < budget;
            if (stoppedEarly) stops++; else runsOut++;
            for (int j = 1; j <= rounds; j++)
            {
                ConvergenceVerdict replay = new ChainEnsemble().Run(model, chains, seed, burn: 0, samples: j * batch, rHatTarget: target,
                                                                    minEssPerChain: 20, stoppingScope: scope).Summary.Verdict;
                bool converged = replay.Outcome == ConvergenceOutcome.Converged;
                if (j < rounds) Assert.False(converged, $"trial {trial}: the verdict converges at round {j}, the stop not until {rounds}");
                else if (stoppedEarly) Assert.True(converged, $"trial {trial}: the stop fired at round {j}, the verdict reads {replay.Outcome}");
                if (j == rounds) Assert.Equal(run.Verdict.Outcome, replay.Outcome);
            }
        }
        Assert.True(stops > 0 && runsOut > 0, $"{stops} runs stopped and {runsOut} ran out; the trials should hold both");
    }

    // Chain i's functional draws from rng, which the chain's seed starts, shared by its functionals in order.
    private delegate Func<double> Source(int chainIndex, Xoshiro256PlusPlus rng);

    private static Func<double> Noise(int index, Xoshiro256PlusPlus rng) => () =>
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    };

    private static Source Constant(double value) => (_, _) => () => value;

    // Noise about a level that alternates with the chain: 0 for even chains, `shift` for odd ones.
    private static Source Levels(double shift) => Shifted(Noise, shift);

    private static Source Shifted(Source source, double shift) => (index, rng) =>
    {
        Func<double> draw = source(index, rng);
        double level = index % 2 == 0 ? 0.0 : shift;
        return () => level + draw();
    };

    private static Source Ar(double phi) => (index, rng) =>
    {
        Func<double> noise = Noise(index, rng);
        double x = 0.0, scale = Math.Sqrt(1.0 - phi * phi);
        return () => x = phi * x + scale * noise();
    };

    // A model whose chain index comes from the shell's own seed derivation, so a shorter run replays a longer one's first
    // draws exactly; groups and names as declared, or the interface's defaults when null.
    private sealed class GroupedModel : IEnsembleModel<double[]>
    {
        private readonly Source[] _sources;
        private readonly int[] _seeds;
        private readonly IReadOnlyList<FunctionalGroup>? _groups;
        private readonly IReadOnlyList<string>? _names;

        public GroupedModel(int masterSeed, int chains, int monitored, IReadOnlyList<FunctionalGroup>? groups,
                            IReadOnlyList<string>? names, params Source[] sources)
        {
            _sources = sources;
            _seeds = SeedTree.Derive(masterSeed, 2 * chains)[..chains];
            _groups = groups;
            _names = names;
            MonitoredDim = monitored;
        }

        public int FunctionalDim => _sources.Length;
        public int MonitoredDim { get; }

        IReadOnlyList<string> IEnsembleModel<double[]>.FunctionalNames
            => _names ?? ((IEnsembleModel<double[]>)new DefaultsOf(FunctionalDim)).FunctionalNames;

        IReadOnlyList<FunctionalGroup> IEnsembleModel<double[]>.Groups
            => _groups ?? ((IEnsembleModel<double[]>)new DefaultsOf(FunctionalDim)).Groups;

        public IChainRun<double[]> StartChain(int chainSeed, int readoutSeed)
        {
            int index = Array.IndexOf(_seeds, chainSeed);
            var rng = new Xoshiro256PlusPlus(chainSeed);
            var generators = new Func<double>[_sources.Length];
            for (int f = 0; f < _sources.Length; f++) generators[f] = _sources[f](index, rng);
            return new ScriptedRun(generators);
        }
    }

    // The interface's defaults for a model of `dim` functionals that declares nothing.
    private sealed class DefaultsOf(int dim) : IEnsembleModel<double[]>
    {
        public int FunctionalDim => dim;
        public int MonitoredDim => dim;
        public IChainRun<double[]> StartChain(int chainSeed, int readoutSeed) => throw new NotSupportedException();
    }

    private sealed class ScriptedRun(Func<double>[] generators) : IChainRun<double[]>, IChain<double[]>
    {
        public IChain<double[]> Chain => this;
        public long Accepted => 0;
        public long Attempts => 0;
        public void Burn(int steps) { }
        public void Accumulate(in double[] draw, Span<double> functionals) => draw.CopyTo(functionals);

        public double[] Step()
        {
            var draw = new double[generators.Length];
            for (int f = 0; f < draw.Length; f++) draw[f] = generators[f]();
            return draw;
        }
    }
}
